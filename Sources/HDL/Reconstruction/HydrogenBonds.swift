//
//  HydrogenBonds.swift
//  HDL
//
//  Created by Philip Turner on 7/2/25.
//

extension Compilation {
  // A reduced form each hydrogen atom, with the 4th vector slot storing
  // the index of the carbon that spawned it.
  //
  // Inputs: material -> bond length
  //         topology.atoms
  //         topology.bonds -> orbitals
  func createHydrogenData() -> [SIMD4<Float>] {
    let orbitalLists = topology.nonbondingOrbitals()
    let bondLength = createBondLength()
    
    var output: [SIMD4<Float>] = []
    for i in topology.atoms.indices {
      let atom = topology.atoms[i]
      let orbitalList = orbitalLists[i]
      for orbital in orbitalList {
        let position = atom.position + bondLength * orbital
        let encodedID = Float(bitPattern: UInt32(i))
        output.append(SIMD4(position, encodedID))
      }
    }
    return output
  }
  
  // Inputs: topology.atoms
  private func hydrogenSiteCenter(_ atomList: [UInt32]) -> SIMD3<Float> {
    var center: SIMD3<Float> = .zero
    for atomID in atomList {
      let atom = topology.atoms[Int(atomID)]
      center += atom.position
    }
    center /= Float(atomList.count)
    return center
  }
  
  // Inputs:  topology.atoms.count
  // Outputs: atomsToHydrogensMap
  //          hydrogensToAtomsMap (length determined here)
  mutating func createHydrogenSites(
    hydrogenData: [SIMD4<Float>]
  ) {
    guard hydrogensToAtomsMap.count == 0,
          atomsToHydrogensMap.count == 0 else {
      fatalError("Maps were not empty.")
    }
    
    func createMatches() -> [Topology.MatchStorage] {
      let hydrogenAtoms = hydrogenData.map {
        var atom = $0
        atom.w = 1
        return atom
      }
      
      // Create a transient topology to de-duplicate the hydrogens and merge
      // references between them.
      var matcher = Topology()
      matcher.insert(atoms: hydrogenAtoms)
      return matcher.match(
        hydrogenAtoms, algorithm: .absoluteRadius(0.050))
    }
    
    func filter(
      matches: [Topology.MatchStorage]
    ) -> [Topology.MatchStorage] {
      var output: [Topology.MatchStorage] = []
      for i in hydrogenData.indices {
        let match = matches[i]
        
        func isCompatible() -> Bool {
          if match.count > 1 {
            for j in match where i != j {
              if i > j {
                return false
              }
            }
          }
          return true
        }
        
        if isCompatible() {
          output.append(match)
        }
      }
      return output
    }
    
    func createFilteredMatches() -> [Topology.MatchStorage] {
      let matches = createMatches()
      let filteredMatches = filter(matches: matches)
      return filteredMatches
    }
    
    // Initialize the hydrogen lists.
    atomsToHydrogensMap = Array(repeating: [], count: topology.atoms.count)
    
    // Fill each list of hydrogens.
    let filteredMatches = createFilteredMatches()
    for match in filteredMatches {
      func createAtomList() -> [UInt32] {
        var atomList: [UInt32] = []
        for j in match {
          let data = hydrogenData[Int(j)]
          let atomID = data.w.bitPattern
          atomList.append(atomID)
        }
        if atomList.count >= 4 {
          fatalError("4-way collisions are not handled yet.")
        }
        
        // Sort the atom list, in-place.
        atomList.sort()
        return atomList
      }
      
      // Integrate the atom list into the map.
      let atomList = createAtomList()
      let hydrogenID = UInt32(hydrogensToAtomsMap.count)
      hydrogensToAtomsMap.append(atomList)
      
      // Mutate the hydrogen list.
      for j in atomList {
        // Appending to the array in-place has a measurable performance
        // improvement, compared to extract + modify + insert.
        atomsToHydrogensMap[Int(j)].append(hydrogenID)
      }
    }
    
    // Sort each hydrogen list, in-place.
    for j in topology.atoms.indices {
      atomsToHydrogensMap[Int(j)].sort()
    }
  }
  
  // Inputs:  topology.atoms
  //          topology.bonds -> orbitals
  //          atomsToHydrogensMap
  //          hydrogensToAtomsMap
  // Outputs: topology.atoms (insert)
  //          topology.bonds (insert)
  mutating func createHydrogenBonds() {
    let orbitalLists = topology.nonbondingOrbitals()
    
    // TODO: Refactor to just return an atom and a bond.
    var insertedAtoms: [Atom] = []
    var insertedBonds: [SIMD2<UInt32>] = []
    func addBond(
      sourceAtomID: UInt32,
      orbital: SIMD3<Float>
    ) {
      let atom = topology.atoms[Int(sourceAtomID)]
      guard let element = Element(rawValue: atom.atomicNumber) else {
        fatalError("This should never happen.")
      }
      
      var bondLength = element.covalentRadius
      bondLength += Element.hydrogen.covalentRadius
      let position = atom.position + bondLength * orbital
      let hydrogenID = topology.atoms.count + insertedAtoms.count
      
      let hydrogen = Atom(position: position, element: .hydrogen)
      let bond = SIMD2(sourceAtomID, UInt32(hydrogenID))
      insertedAtoms.append(hydrogen)
      insertedBonds.append(bond)
    }
    
    // TODO: Refactor to just return an atom and a bond.
    func handleSingleAtom(
      hydrogenID: UInt32,
      atomID: UInt32
    ) {
      // Create a list of collision marks.
      let hydrogenList = atomsToHydrogensMap[Int(atomID)]
      let collisionMask: [Bool] = hydrogenList.map {
        let atomList = hydrogensToAtomsMap[Int($0)]
        switch atomList.count {
        case 0:
          fatalError("This should never happen.")
        case 1:
          return false
        default:
          return true
        }
      }
      
      // Retrieve the list of orbitals.
      let orbitalList = orbitalLists[Int(atomID)]
      guard orbitalList.count == hydrogenList.count else {
        fatalError("Unexpected orbital count.")
      }
      
      // All control flow paths result in just 1 call to 'addBond'.
      switch hydrogenList.count {
      case 1:
        guard collisionMask[0] == false else {
          fatalError("This should never happen.")
        }
        
        // 70 hydrogens from case 1
        addBond(
          sourceAtomID: atomID,
          orbital: orbitalList[0])
      case 2:
        if collisionMask[0] && collisionMask[1] {
          fatalError("This should never happen.")
        } else if collisionMask[0] || collisionMask[1] {
          let collisionID =
          collisionMask[0] ? hydrogenList[0] : hydrogenList[1]
          let nonCollisionID =
          collisionMask[0] ? hydrogenList[1] : hydrogenList[0]
          guard hydrogenID == nonCollisionID,
                hydrogenID != collisionID else {
            fatalError("Unexpected hydrogen IDs for collision.")
          }
          
          let atomList = hydrogensToAtomsMap[Int(collisionID)]
          let siteCenter = hydrogenSiteCenter(atomList)
          let atom = topology.atoms[Int(atomID)]
          
          let delta = siteCenter - atom.position
          let score0 = (orbitalList[0] * delta).sum()
          let score1 = (orbitalList[1] * delta).sum()
          
//          var isCulprit = false
//          do {
//            var matchedChains: [SIMD4<UInt32>] = []
//            for chain in rawDimerChains {
//              let matches0 = any(chain .== atomList[0])
//              let matches1 = any(chain .== atomList[1])
//              if matches0 && matches1 {
//                matchedChains.append(chain)
//              }
//            }
//            
//            if matchedChains.count > 0 {
//              let usedChain = matchedChains[0]
//              if usedChain[1] == usedChain[2] {
//                if usedChain[2] == atomList[0] ||
//                    usedChain[3] == atomList[1] {
//                  print()
//                  print("skipped2")
//                  
//                  // Mechanism for finding 6 of the culprit hydrogen sites.
//                  isCulprit = true
//                  print(usedChain)
//                  
//                  print(siteCenter)
//                  print(atomID, atom.position, delta)
//                  print("  \(orbitalList[0]) \(score0)")
//                  print("  \(orbitalList[1]) \(score1)")
//                }
//              }
//            }
//          }
          
          // orbitalList[1] is always causing the problem
          
          // Prefer the orbital that doesn't correspond to the collision site.
          if score0 > score1 {
            // 15 hydrogens from case 1, all 30 culprits
            addBond(
              sourceAtomID: atomID,
              orbital: orbitalList[1])
//            if isCulprit {
//              print("w \(orbitalList[1]) \(score1)")
//            }
          } else if score0 < score1 {
            // 15 hydrogens
            addBond(
              sourceAtomID: atomID,
              orbital: orbitalList[0])
//            if isCulprit {
//              print("r \(orbitalList[0]) \(score0)")
//            }
          } else {
            fatalError("Scores were equal.")
          }
        } else {
          // Assign the first hydrogen in the list to the first orbital.
          if hydrogenID == hydrogenList[0] {
            // 12 hydrogens
            addBond(
              sourceAtomID: atomID,
              orbital: orbitalList[0])
          } else {
            // 12 hydrogens
            addBond(
              sourceAtomID: atomID,
              orbital: orbitalList[1])
          }
        }
      default:
        fatalError("Unexpected size for hydrogen list.")
      }
    }
    
    for hydrogenID in hydrogensToAtomsMap.indices {
      let atomList = hydrogensToAtomsMap[hydrogenID]
      guard atomList.count > 0 else {
        // This collision was resolved.
        continue
      }
      
      switch atomList.count {
      case 1:
        // Move this out of the conditional statement for legibility.
        
        // 124 hydrogens
        handleSingleAtom(
          hydrogenID: UInt32(hydrogenID),
          atomID: atomList[0])
      case 2:
        let siteCenter = hydrogenSiteCenter(atomList)
        
//        var matchedChains: [SIMD4<UInt32>] = []
//        for chain in rawDimerChains {
//          let matches0 = any(chain .== atomList[0])
//          let matches1 = any(chain .== atomList[1])
//          if matches0 && matches1 {
//            matchedChains.append(chain)
//          }
//        }
//        
//        var isCulprit = false
//        if matchedChains.count > 0 {
//          let usedChain = matchedChains[0]
//          if usedChain[1] == usedChain[2] {
//            if usedChain[2] == atomList[0] ||
//                usedChain[3] == atomList[1] {
//              print()
//              print("skipped")
//              
//              // Mechanism for finding 6 of the culprit hydrogen sites.
//              isCulprit = true
//              print(usedChain)
//            }
//          }
//        }
//        
//        if isCulprit {
//          print(siteCenter)
//        }
        
        for atomID in atomList {
          let atom = topology.atoms[Int(atomID)]
          let delta = siteCenter - atom.position
//          if isCulprit {
//            print(atomID, atom.position, delta)
//          }
          
          var bestOrbital: SIMD3<Float>?
          var bestOrbitalID: Int?
          var bestScore: Float = -.greatestFiniteMagnitude
          let orbitalList = orbitalLists[Int(atomID)]
          
          for orbitalID in orbitalList.indices {
            let orbital = orbitalList[orbitalID]
            let score = (orbital * delta).sum()
//            if isCulprit {
//              print("  \(orbital) \(score)")
//            }
            
            if score > bestScore {
              bestOrbital = orbital
              bestOrbitalID = orbitalID
              bestScore = score // solution to the bug?
            }
          }
//          if isCulprit {
//            print("w \(bestOrbital!) \(bestScore)")
//          }
          guard let bestOrbital,
                let bestOrbitalID else {
            fatalError("Could not find best orbital.")
          }
          
          // 108 hydrogens
          //
          // 30 hydrogens from this group end up overlapping those from case 1
          // orbitalID = 0: 30 hydrogens, all of the culprits
          // orbitalID = 1: 78 hydrogens
          addBond(
            sourceAtomID: atomID, // change 'sourceAtomID' to 'atomID'
            orbital: bestOrbital)
        }
      default:
        fatalError("3/4-way collisions should have been caught.")
      }
    }
    topology.insert(atoms: insertedAtoms)
    topology.insert(bonds: insertedBonds)
  }
}

/*
private let rawDimerChains: [SIMD4<UInt32>] = [
  SIMD4<UInt32>(10, 136, 136, 132),
  SIMD4<UInt32>(12, 36, 36, 34),
  SIMD4<UInt32>(17, 144, 277, 273),
  SIMD4<UInt32>(19, 44, 68, 66),
  SIMD4<UInt32>(24, 152, 418, 414),
  SIMD4<UInt32>(26, 52, 100, 98),
  SIMD4<UInt32>(28, 193, 543, 544),
  SIMD4<UInt32>(30, 161, 161, 130),
  SIMD4<UInt32>(34, 36, 36, 12),
  SIMD4<UInt32>(58, 84, 108, 106),
  SIMD4<UInt32>(60, 227, 402, 403),
  SIMD4<UInt32>(62, 195, 302, 271),
  SIMD4<UInt32>(66, 68, 44, 19),
  SIMD4<UInt32>(90, 116, 116, 114),
  SIMD4<UInt32>(92, 261, 261, 262),
  SIMD4<UInt32>(94, 229, 443, 412),
  SIMD4<UInt32>(98, 100, 52, 26),
  SIMD4<UInt32>(106, 108, 84, 58),
  SIMD4<UInt32>(114, 116, 116, 90),
  SIMD4<UInt32>(125, 264, 550, 551),
  SIMD4<UInt32>(126, 266, 409, 410),
  SIMD4<UInt32>(127, 268, 268, 269),
  SIMD4<UInt32>(130, 161, 161, 30),
  SIMD4<UInt32>(132, 136, 136, 10),
  SIMD4<UInt32>(156, 293, 426, 422),
  SIMD4<UInt32>(160, 334, 509, 510),
  SIMD4<UInt32>(231, 370, 477, 445),
  SIMD4<UInt32>(262, 261, 261, 92),
  SIMD4<UInt32>(263, 405, 548, 549),
  SIMD4<UInt32>(269, 268, 268, 127),
  SIMD4<UInt32>(271, 302, 195, 62),
  SIMD4<UInt32>(273, 277, 144, 17),
  SIMD4<UInt32>(297, 434, 434, 430),
  SIMD4<UInt32>(301, 475, 475, 476),
  SIMD4<UInt32>(372, 511, 511, 479),
  SIMD4<UInt32>(403, 402, 227, 60),
  SIMD4<UInt32>(404, 546, 546, 547),
  SIMD4<UInt32>(410, 409, 266, 126),
  SIMD4<UInt32>(412, 443, 229, 94),
  SIMD4<UInt32>(414, 418, 152, 24),
  SIMD4<UInt32>(422, 426, 293, 156),
  SIMD4<UInt32>(430, 434, 434, 297),
  SIMD4<UInt32>(445, 477, 370, 231),
  SIMD4<UInt32>(476, 475, 475, 301),
  SIMD4<UInt32>(479, 511, 511, 372),
  SIMD4<UInt32>(510, 509, 334, 160),
  SIMD4<UInt32>(544, 543, 193, 28),
  SIMD4<UInt32>(547, 546, 546, 404),
  SIMD4<UInt32>(549, 548, 405, 263),
  SIMD4<UInt32>(551, 550, 264, 125),
  SIMD4<UInt32>(552, 557, 575, 576),
  SIMD4<UInt32>(553, 559, 568, 569),
  SIMD4<UInt32>(554, 561, 561, 562),
  SIMD4<UInt32>(556, 564, 573, 574),
  SIMD4<UInt32>(562, 561, 561, 554),
  SIMD4<UInt32>(563, 571, 571, 572),
  SIMD4<UInt32>(569, 568, 559, 553),
  SIMD4<UInt32>(572, 571, 571, 563),
  SIMD4<UInt32>(574, 573, 564, 556),
  SIMD4<UInt32>(576, 575, 557, 552),
]
*/
