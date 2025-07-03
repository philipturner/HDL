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
  // Inputs:  material -> bond length
  //          topology.atoms
  //          topology.bonds -> orbitals
  // Outputs: [SIMD4<Float>]
  func createHydrogenData() -> [SIMD4<Float>] {
    let orbitalLists = topology.nonbondingOrbitals(hybridization: .sp3)
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
  
  // Inputs:  topology.atoms
  //          list of 1-3 indices
  // Outputs: SIMD3<Float>
  private func hydrogenSiteCenter(_ atomList: [UInt32]) -> SIMD3<Float> {
    var center: SIMD3<Float> = .zero
    for atomID in atomList {
      let atom = topology.atoms[Int(atomID)]
      center += atom.position
    }
    center /= Float(atomList.count)
    return center
  }
  
  // Inputs:  hydrogenData
  //          topology.atoms.count
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
    
    // Initialize each list of hydrogens.
    atomsToHydrogensMap = Array(repeating: [], count: topology.atoms.count)
    
    // Fill each list of hydrogens.
    for match in createFilteredMatches() {
      // Create the sorted list of atoms.
      var atomList: [UInt32] = []
      for j in match {
        let data = hydrogenData[Int(j)]
        let atomID = data.w.bitPattern
        atomList.append(atomID)
      }
      guard atomList.count <= 3 else {
        fatalError("Edge case with 4 hydrogens in a site not handled yet.")
      }
      atomList.sort()
      
      // Integrate this list into the bidirectional map.
      let hydrogenID = UInt32(hydrogensToAtomsMap.count)
      hydrogensToAtomsMap.append(atomList)
      for j in atomList {
        // Appending to the array in-place has a measurable performance
        // improvement, compared to extract + modify + insert.
        atomsToHydrogensMap[Int(j)].append(hydrogenID)
      }
    }
    
    // Sort each list of hydrogens, in-place.
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
    let orbitalLists = topology.nonbondingOrbitals(hybridization: .sp3)
    
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
        addBond(
          sourceAtomID: atomID,
          orbital: orbitalList[orbitalList.startIndex])
      case 2:
        let orbital0 = orbitalList[orbitalList.startIndex]
        let orbital1 = orbitalList[orbitalList.endIndex - 1]
        
        if collisionMask[0] && collisionMask[1] {
          fatalError("This should never happen.")
        } else if collisionMask[0] || collisionMask[1] {
          let collisionID =
          (collisionMask[0]) ? hydrogenList[0] : hydrogenList[1]
          let nonCollisionID =
          (collisionMask[0]) ? hydrogenList[1] : hydrogenList[0]
          guard hydrogenID == nonCollisionID,
                hydrogenID != collisionID else {
            fatalError("Unexpected hydrogen IDs for collision.")
          }
          
          let atomList = hydrogensToAtomsMap[Int(collisionID)]
          let siteCenter = hydrogenSiteCenter(atomList)
          let atomCenter = topology.atoms[Int(atomID)].position
          
          let delta = siteCenter - atomCenter
          let score0 = (orbital0 * delta).sum()
          let score1 = (orbital1 * delta).sum()
          
          // Use a scoring function to match the collision to an orbital.
          if score0 > score1 {
            addBond(
              sourceAtomID: atomID,
              orbital: orbital1)
          } else if score0 < score1 {
            addBond(
              sourceAtomID: atomID,
              orbital: orbital0)
          } else {
            fatalError("Scores were equal.")
          }
        } else {
          // Assign the first hydrogen in the list to the first orbital.
          if hydrogenID == hydrogenList[0] {
            addBond(
              sourceAtomID: atomID,
              orbital: orbital0)
          } else {
            addBond(
              sourceAtomID: atomID,
              orbital: orbital1)
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
        handleSingleAtom(
          hydrogenID: UInt32(hydrogenID),
          atomID: atomList[0])
      case 2:
        let siteCenter = hydrogenSiteCenter(atomList)
        for atomID in atomList {
          let orbitalList = orbitalLists[Int(atomID)]
          
          // TODO: Refactor this statement. It's too dense.
          let delta = siteCenter - topology.atoms[Int(atomID)].position
          
          // TODO: Refactor this statement. It's too dense.
          var keyValuePairs = orbitalList.map { orbital -> (SIMD3<Float>, Float) in
            (orbital, (orbital * delta).sum())
          }
          keyValuePairs.sort(by: { $0.1 > $1.1 })
          
          let closestOrbital = keyValuePairs[0].0
          addBond(
            sourceAtomID: atomID, // change 'sourceAtomID' to 'atomID'
            orbital: closestOrbital)
        }
      default:
        fatalError("3/4-way collisions should have been caught.")
      }
    }
    topology.insert(atoms: insertedAtoms)
    topology.insert(bonds: insertedBonds)
  }
}
