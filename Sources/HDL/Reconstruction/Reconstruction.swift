//
//  Reconstruction.swift
//
//
//  Created by Philip Turner on 6/11/24.
//

// Only applies to crystals with sp3 hybridization. Makes heavy use of
// 'Topology.nonbondingOrbitals' with sp3 hybridization.

// Revised user-facing API:
//
// struct Reconstruction {
//   var material: MaterialType?
//   var atoms: [SIMD4<Float>]?
//
//   func compile() -> Topology
// }
//
// var reconstruction = Reconstruction()
// reconstruction.materialType = ...
// reconstruction.atoms = ...
// let topology = reconstruction.compile()

public struct Reconstruction {
  var _material: MaterialType?
  
  public var material: MaterialType {
    get {
      fatalError("There is no getter for 'material'.")
    }
    set {
      _material = newValue
    }
  }
  
  public var topology: Topology = Topology()
  
  // Each atom has a list corresponding to it. That list is guaranteed to
  // always be sorted in a deterministic order.
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
    
    // When there are no more 3-way collisions, the list of atoms will stop
    // growing. That means the center types will stabilize.
    var stableCenterTypes: [UInt8]?
    
    // Loop over this a few times (typically less than 10).
    for _ in 0..<100 {
      // Fill the list of center types.
      let centerTypes = createBulkAtomBonds()
      
      // Crash if 4-way collisions exist.
      createHydrogenSites()
      
      // Check whether there are still 3-way collisions.
      if hydrogensToAtomsMap.contains(where: { $0.count > 2 }) {
        // Add center atoms to problematic sites.
        resolveThreeWayCollisions()
        
        // Reverse the actions from the start of this iteration.
        topology.bonds = []
        hydrogensToAtomsMap = []
        atomsToHydrogensMap = []
      } else {
        stableCenterTypes = centerTypes
        break
      }
    }
    guard let stableCenterTypes else {
      fatalError("Could not resolve 3-way collisions.")
    }
    
    // Assert that there are no 3-way or 4-way collisions.
    for atomList in hydrogensToAtomsMap {
      if atomList.count > 2 {
        fatalError("3/4-way collisions should have been caught.")
      }
    }
    
    // Add hydrogens after the center atoms are fixed.
    validate(centerTypes: stableCenterTypes)
    resolveTwoWayCollisions(centerTypes: stableCenterTypes)
    createHydrogenBonds()
  }
}

// MARK: - Utilities

extension Reconstruction {
  // Remove atoms with less than two covalent bonds.
  private mutating func removePathologicalAtoms() {
    // Loop over this a few times (typically less than 10).
    var converged = false
    for _ in 0..<100 {
      let bondLength = createBondLength()
      let matches = topology.match(
        topology.atoms, algorithm: .absoluteRadius(bondLength * 1.1))
      
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
      
      if removedAtoms.count > 0 {
        topology.remove(atoms: removedAtoms)
      } else {
        converged = true
        break
      }
    }
    guard converged else {
      fatalError("Could not remove pathological atoms.")
    }
  }
  
  // Form all center atom bonds in the lattice interior.
  //
  // Returns the center type of each atom.
  private mutating func createBulkAtomBonds() -> [UInt8] {
    let bondLength = createBondLength()
    let matches = topology.match(
      topology.atoms, algorithm: .absoluteRadius(bondLength * 1.1))
    
    var insertedBonds: [SIMD2<UInt32>] = []
    var centerTypes: [UInt8] = []
    for i in topology.atoms.indices {
      let match = matches[i]
      if match.count > 5 {
        fatalError("Unexpected situation: match count > 5")
      } else if match.count > 2 {
        centerTypes.append(UInt8(match.count - 1))
        
        for j in match where i < j {
          insertedBonds.append(SIMD2(UInt32(i), j))
        }
      } else {
        fatalError("Pathological atoms should be removed.")
      }
    }
    topology.insert(bonds: insertedBonds)
    
    return centerTypes
  }
  
  // Next, form the hydrogen bonds. Place hydrogens at the C-C bond length
  // instead of the C-H bond length.
  private mutating func createHydrogenSites() {
    guard hydrogensToAtomsMap.count == 0,
          atomsToHydrogensMap.count == 0 else {
      fatalError("Maps were not empty.")
    }
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
    
    // TODO: Refactor this to remove incomprehensible control flow.
    
    var compatibleMatches: [Topology.MatchStorage] = []
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
      compatibleMatches.append(match)
    }
    
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
    
    for j in topology.atoms.indices {
      atomsToHydrogensMap[Int(j)].sort()
    }
  }
  
  private func validate(centerTypes: [UInt8]) {
    let orbitals = topology.nonbondingOrbitals(hybridization: .sp3)
    for i in orbitals.indices {
      let orbital = orbitals[i]
      var expectedRawValue: UInt8
      
      switch orbital.count {
      case 2:
        expectedRawValue = 2
      case 1:
        expectedRawValue = 3
      case 0:
        expectedRawValue = 4
      default:
        fatalError("This should never happen.")
      }
      
      guard centerTypes[i] == expectedRawValue else {
        fatalError("Incorrect raw value.")
      }
    }
  }
  
  private mutating func createHydrogenBonds() {
    // Utility function: 'withClosestOrbitals'
    let orbitals = topology.nonbondingOrbitals(hybridization: .sp3)
    func withClosestOrbitals(
      _ atomList: [UInt32],
      _ closure: (UInt32, SIMD3<Float>) -> Void
    ) {
      let siteCenter = hydrogenSiteCenter(atomList)
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
    
    // Utility function: 'addBond'
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
    
    // Utility function: 'handleAtomListCount1'
    func handleAtomListCount1(
      sourceAtomID: UInt32,
      neighborID: UInt32
    ) {
      // Create a list of collision marks.
      let hydrogenList = atomsToHydrogensMap[Int(neighborID)]
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
      let orbital = orbitals[Int(neighborID)]
      guard orbital.count == hydrogenList.count else {
        fatalError("Unexpected orbital count.")
      }
      
      // Switch over the different cases of the atom's hydrogen list.
      switch hydrogenList.count {
      case 1:
        // Easiest case:
        //
        // The list only has a single hydrogen.
        guard collisionMask[0] == false else {
          fatalError("This should never happen.")
        }
        addBond(
          sourceAtomID: neighborID,
          orbital: orbital[orbital.startIndex])
        
      case 2:
        let orbital0 = orbital[orbital.startIndex]
        let orbital1 = orbital[orbital.endIndex - 1]
        
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
          guard sourceAtomID == nonCollisionID,
                sourceAtomID != collisionID else {
            fatalError("Unexpected atom IDs for collision.")
          }
          
          let atomList = hydrogensToAtomsMap[Int(collisionID)]
          let siteCenter = hydrogenSiteCenter(atomList)
          let neighborCenter = topology.atoms[Int(neighborID)].position
          
          let delta = siteCenter - neighborCenter
          let score0 = (orbital0 * delta).sum()
          let score1 = (orbital1 * delta).sum()
          
          if score0 > score1 {
            addBond(
              sourceAtomID: neighborID,
              orbital: orbital1)
          } else if score0 < score1 {
            addBond(
              sourceAtomID: neighborID,
              orbital: orbital0)
          } else {
            fatalError("Scores were equal.")
          }
        } else {
          // If there are 2 orbitals and both are collision-free:
          //
          // The compiler uses a deterministic method to generate orbitals.
          // Plus, the orbitals are already generated once. Assign the first
          // hydrogen in the list to the first orbital.
          if sourceAtomID == hydrogenList[0] {
            addBond(
              sourceAtomID: neighborID,
              orbital: orbital0)
          } else {
            addBond(
              sourceAtomID: neighborID,
              orbital: orbital1)
          }
        }
        
      default:
        fatalError("Large hydrogen lists not handled yet.")
      }
    }
    
    // High-level specification of the algorithm structure.
    for i in hydrogensToAtomsMap.indices {
      let atomList = hydrogensToAtomsMap[i]
      guard atomList.count > 0 else {
        // This collision was resolved.
        continue
      }
      
      switch atomList.count {
      case 1:
        // Move this out of the conditional statement for legibility.
        handleAtomListCount1(
          sourceAtomID: UInt32(i),
          neighborID: atomList[0])
      case 2:
        withClosestOrbitals(atomList) { atomID, orbital in
          addBond(
            sourceAtomID: atomID,
            orbital: orbital)
        }
      case 3:
        withClosestOrbitals(atomList) { atomID, orbital in
          addBond(
            sourceAtomID: atomID,
            orbital: orbital)
        }
      default:
        fatalError("This should never happen.")
      }
    }
    topology.insert(atoms: insertedAtoms)
    topology.insert(bonds: insertedBonds)
  }
  
  private func hydrogenSiteCenter(_ atomList: [UInt32]) -> SIMD3<Float> {
    var center: SIMD3<Float> = .zero
    for atomID in atomList {
      let atom = topology.atoms[Int(atomID)]
      center += atom.position
    }
    center /= Float(atomList.count)
    return center
  }
}
