//
//  Reconstruction.swift
//
//
//  Created by Philip Turner on 6/11/24.
//

public struct Reconstruction {
  public var _material: MaterialType?
  
  public var material: MaterialType {
    get {
      fatalError("There is no getter for 'material'.")
    }
    set {
      _material = newValue
    }
  }
  
  public var topology: Topology = Topology()
  
  var initialTypeRawValues: [UInt8] = []
  
  // These lists must always be sorted.
  var hydrogensToAtomsMap: [[UInt32]] = []
  var atomsToHydrogensMap: [[UInt32]] = []
  
  public init() {
    
  }
}

extension Reconstruction {
  func createBondLength() -> Float {
    var bondLength: Float
    switch _material {
    case .elemental(let element):
      bondLength = 2 * element.covalentRadius
    case .checkerboard(let element, let element2):
      bondLength = element.covalentRadius + element2.covalentRadius
    case nil:
      fatalError("Material not specified.")
    }
    return bondLength
  }
  
  // Fix the center atoms and add hydrogen termination.
  public mutating func compile() {
    // Invoke this function once.
    removePathologicalAtoms()
    
    // Loop over this several times.
    var converged = false
    for _ in 0..<100 {
      createBulkAtomBonds()
      createHydrogenSites()
      
      if hydrogensToAtomsMap.contains(where: { $0.count > 2 }) {
        // Add center atoms to problematic sites.
        resolveThreeWayCollisions()
        
        // Reverse the actions from the start of this iteration.
        topology.bonds = []
        initialTypeRawValues = []
        hydrogensToAtomsMap = []
        atomsToHydrogensMap = []
      } else {
        converged = true
        break
      }
    }
    guard converged else {
      fatalError("Could not resolve 3-way collisions.")
    }
    
    // Add hydrogens after the center atoms are fixed.
    resolveTwoWayCollisions()
    createHydrogenBonds()
  }
}

extension Reconstruction {
  // Remove atoms with less than two covalent bonds.
  mutating func removePathologicalAtoms() {
    var converged = false
    for _ in 0..<100 {
      let matches = topology.match(
        topology.atoms, algorithm: .absoluteRadius(createBondLength() * 1.1))
      var removedAtoms: [UInt32] = []
      
      for i in topology.atoms.indices {
        let match = matches[i]
        if match.count > 5 {
          fatalError("Unexpected situation: match count > 5")
        } else if match.count > 2 {
          
        } else {
          removedAtoms.append(UInt32(i))
        }
      }
      
      if removedAtoms.count == 0 {
        converged = true
        break
      } else {
        topology.remove(atoms: removedAtoms)
      }
    }
    
    guard converged else {
      fatalError("Could not remove pathological atoms.")
    }
  }
  
  // Form all center atom bonds in the lattice interior, assign center types.
  mutating func createBulkAtomBonds() {
    let matches = topology.match(
      topology.atoms, algorithm: .absoluteRadius(createBondLength() * 1.1))
    var insertedBonds: [SIMD2<UInt32>] = []
    
    for i in topology.atoms.indices {
      let match = matches[i]
      if match.count > 5 {
        fatalError("Unexpected situation: match count > 5")
      } else if match.count > 2 {
        initialTypeRawValues.append(UInt8(match.count - 1))
        
        for j in match where i < j {
          insertedBonds.append(SIMD2(UInt32(i), j))
        }
      } else {
        fatalError("Pathological atoms should be removed.")
      }
    }
    
    topology.insert(bonds: insertedBonds)
  }
  
  // Next, form the hydrogen bonds. Place hydrogens at the C-C bond length
  // instead of the C-H bond length.
  mutating func createHydrogenSites() {
    precondition(hydrogensToAtomsMap.count == 0, "Map not empty.")
    precondition(atomsToHydrogensMap.count == 0, "Map not empty.")
    atomsToHydrogensMap = Array(repeating: [], count: topology.atoms.count)
    
    let orbitals = topology.nonbondingOrbitals(hybridization: .sp3)
    let bondLength = createBondLength()
    var hydrogenData: [SIMD4<Float>] = []
    
    for i in topology.atoms.indices {
      let atom = topology.atoms[i]
      for orbital in orbitals[i] {
        let position = atom.position + bondLength * orbital
        let encodedID = Float(bitPattern: UInt32(i))
        hydrogenData.append(SIMD4(position, encodedID))
      }
    }
    
    // Create a transient topology to de-duplicate the hydrogens and merge
    // references between them.
    let hydrogenAtoms = hydrogenData.map {
      var atom = $0
      atom.w = 1
      return atom
    }
    var matcher = Topology()
    matcher.insert(atoms: hydrogenAtoms)
    let matches = matcher.match(
      hydrogenAtoms, algorithm: .absoluteRadius(0.050))
    
  outer:
    for i in hydrogenData.indices {
      let match = matches[i]
      if match.count > 1 {
        for j in match where i != j {
          if i > j {
            continue outer
          }
        }
      }
      
      let hydrogenID = UInt32(hydrogensToAtomsMap.count)
      var atomList: [UInt32] = []
      for j in match {
        let data = hydrogenData[Int(j)]
        let atomID = data.w.bitPattern
        atomList.append(atomID)
      }
      atomList.sort()
      hydrogensToAtomsMap.append(atomList)
      for j in atomList {
        atomsToHydrogensMap[Int(j)].append(hydrogenID)
      }
    }
    
    for j in topology.atoms.indices {
      atomsToHydrogensMap[Int(j)].sort()
    }
  }
  
  mutating func createHydrogenBonds() {
    var insertedAtoms: [Atom] = []
    var insertedBonds: [SIMD2<UInt32>] = []
    func createCenter(_ atomList: [UInt32]) -> SIMD3<Float>? {
      guard atomList.count > 1 else {
        return nil
      }
      var output: SIMD3<Float> = .zero
      for atomID in atomList {
        let atom = topology.atoms[Int(atomID)]
        output += atom.position
      }
      output /= Float(atomList.count)
      return output
    }
    func addBond(_ atomID: Int, orbital: SIMD3<Float>) {
      let atom = topology.atoms[atomID]
      guard let element = Element(rawValue: atom.atomicNumber) else {
        fatalError("This should never happen.")
      }
      
      var bondLength = element.covalentRadius
      bondLength += Element.hydrogen.covalentRadius
      let position = atom.position + bondLength * orbital
      let hydrogenID = topology.atoms.count + insertedAtoms.count
      
      let hydrogen = Atom(position: position, element: .hydrogen)
      let bond = SIMD2(UInt32(atomID), UInt32(hydrogenID))
      insertedAtoms.append(hydrogen)
      insertedBonds.append(bond)
    }
    func withClosestOrbitals(
      _ atomList: [UInt32],
      _ closure: (UInt32, SIMD3<Float>) -> Void
    ) {
      let siteCenter = createCenter(atomList)!
      for atomID in atomList {
        let orbital = orbitals[Int(atomID)]
        let delta = siteCenter - topology.atoms[Int(atomID)].position
        var keyValuePairs = orbital.map { orbital -> (SIMD3<Float>, Float) in
          (orbital, (orbital * delta).sum())
        }
        keyValuePairs.sort(by: { $0.1 > $1.1 })
        
        let closestOrbital = keyValuePairs[0].0
        closure(atomID, closestOrbital)
      }
    }
    let orbitals = topology.nonbondingOrbitals(hybridization: .sp3)
    
    for i in hydrogensToAtomsMap.indices {
      let atomList = hydrogensToAtomsMap[i]
      
      if atomList.count == 0 {
        // This collision was resolved.
        continue
      } else if atomList.count == 1 {
        let atomID = Int(atomList[0])
        let hydrogenList = atomsToHydrogensMap[atomID]
        let collisionMask = hydrogenList.map {
          let atomList = hydrogensToAtomsMap[Int($0)]
          precondition(atomList.count > 0)
          return atomList.count > 1
        }
        let orbital = orbitals[atomID]
        precondition(orbital.count > 0, "No orbitals.")
        
        // Switch over the different cases of the atom's hydrogen list.
        if hydrogenList.count == 1 {
          precondition(orbital.count == 1, "Unexpected orbital count.")
          
          // Easiest case:
          //
          // The list only has a single hydrogen.
          precondition(!collisionMask[0])
          addBond(atomID, orbital: orbital[orbital.startIndex])
        } else if hydrogenList.count == 2 {
          precondition(orbital.count == 2, "Unexpected orbital count.")
          let orbital0 = orbital[orbital.startIndex]
          let orbital1 = orbital[orbital.endIndex-1]
          
          if collisionMask[0] && collisionMask[1] {
            fatalError("This should never happen.")
          } else if collisionMask[0] || collisionMask[1] {
            // If 1 orbital has a collision:
            //
            // Use a scoring function to match collision(s) to orbitals.
            
            let collisionID =
            (collisionMask[0]) ? hydrogenList[0] : hydrogenList[1]
            let nonCollisionID =
            (collisionMask[0]) ? hydrogenList[1] : hydrogenList[0]
            precondition(collisionID != UInt32(i))
            precondition(nonCollisionID == UInt32(i))
            
            let atomList = hydrogensToAtomsMap[Int(collisionID)]
            let center = createCenter(atomList)!
            let delta = center - topology.atoms[atomID].position
            let score0 = (orbital0 * delta).sum()
            let score1 = (orbital1 * delta).sum()
            
            if score0 > score1 {
              addBond(atomID, orbital: orbital1)
            } else if score0 < score1 {
              addBond(atomID, orbital: orbital0)
            } else {
              fatalError("Scores were equal.")
            }
          } else {
            // If there are 2 orbitals and both are collision-free:
            //
            // The compiler uses a deterministic method to generate orbitals.
            // Plus, the orbitals are already generated once. Assign the first
            // hydrogen in the list to the first orbital.
            let isFirst = hydrogenList[0] == UInt32(i)
            let orbital = isFirst ? orbital0 : orbital1
            addBond(atomID, orbital: orbital)
          }
        } else {
          fatalError("Large hydrogen lists not handled yet.")
        }
      } else if atomList.count == 2 {
        withClosestOrbitals(atomList) { atomID, orbital in
          addBond(Int(atomID), orbital: orbital)
        }
      } else if atomList.count == 3 {
        withClosestOrbitals(atomList) { atomID, orbital in
          addBond(Int(atomID), orbital: orbital)
        }
      } else if atomList.count > 3 {
        fatalError("Edge case with >3 hydrogens in a site not handled yet.")
      }
    }
    topology.insert(atoms: insertedAtoms)
    topology.insert(bonds: insertedBonds)
  }
}
