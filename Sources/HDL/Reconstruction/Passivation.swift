//
//  Passivation.swift
//  HDL
//
//  Created by Philip Turner on 7/2/25.
//

struct PassivationInput {
  var atomsToHydrogensMap: [[UInt32]] = []
  var hydrogensToAtomsMap: [[UInt32]] = []
}

struct PassivationResult {
  var insertedAtoms: [Atom] = []
  var insertedBonds: [SIMD2<UInt32>] = []
}

struct PassivationImpl {
  var atoms: [Atom] = []
  var bonds: [SIMD2<UInt32>] = []
  
  func createOrbitalLists() -> [Topology.OrbitalStorage] {
    var topology = Topology()
    topology.atoms = atoms
    topology.bonds = bonds
    return topology.nonbondingOrbitals()
  }
  
  func hydrogenSiteCenter(_ atomList: [UInt32]) -> SIMD3<Float> {
    var center: SIMD3<Float> = .zero
    for atomID in atomList {
      let atom = atoms[Int(atomID)]
      center += atom.position
    }
    center /= Float(atomList.count)
    return center
  }
  
  func createHydrogen(
    atomID: UInt32,
    orbital: SIMD3<Float>
  ) -> Atom {
    let atom = atoms[Int(atomID)]
    let atomElement = Element(rawValue: atom.atomicNumber)
    guard let atomElement else {
      fatalError("This should never happen.")
    }
    
    var bondLength = atomElement.covalentRadius
    bondLength += Element.hydrogen.covalentRadius
    
    let position = atom.position + bondLength * orbital
    return Atom(position: position, element: .hydrogen)
  }
  
  #if true
  func compile(input: PassivationInput) -> PassivationResult {
    let orbitalLists = createOrbitalLists()
    
    var output = PassivationResult()
    func appendBond(atomID: UInt32, hydrogen: Atom) {
      let hydrogenID = atoms.count + output.insertedAtoms.count
      output.insertedAtoms.append(hydrogen)
      
      let bond = SIMD2(atomID, UInt32(hydrogenID))
      output.insertedBonds.append(bond)
    }
    
    for atomID in atoms.indices {
      let orbitalList = orbitalLists[atomID]
      for orbital in orbitalList {
        let hydrogen = createHydrogen(
          atomID: UInt32(atomID),
          orbital: orbital)
        appendBond(
          atomID: UInt32(atomID),
          hydrogen: hydrogen)
      }
    }
    return output
  }
  #else
  func compile(input: PassivationInput) -> PassivationResult {
    let orbitalLists = createOrbitalLists()
    
    func createHydrogenOrbital(
      atomID: UInt32,
      hydrogenID: UInt32
    ) -> SIMD3<Float> {
      // Create a list of collision marks.
      let hydrogenList = input.atomsToHydrogensMap[Int(atomID)]
      let collisionMask: [Bool] = hydrogenList.map {
        let atomList = input.hydrogensToAtomsMap[Int($0)]
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
      
      switch hydrogenList.count {
      case 1:
        guard collisionMask[0] == false else {
          fatalError("This should never happen.")
        }
        return orbitalList[0]
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
          
          let atomList = input.hydrogensToAtomsMap[Int(collisionID)]
          let siteCenter = hydrogenSiteCenter(atomList)
          let atom = atoms[Int(atomID)]
          
          let delta = siteCenter - atom.position
          let score0 = (orbitalList[0] * delta).sum()
          let score1 = (orbitalList[1] * delta).sum()
          
          // Prefer the orbital that doesn't correspond to the collision site.
          if score0 > score1 {
            return orbitalList[1]
          } else if score0 < score1 {
            return orbitalList[0]
          } else {
            fatalError("Scores were equal.")
          }
        } else {
          // Assign the first hydrogen in the list to the first orbital.
          //
          // TODO: Add a unit test to check whether this order is respected.
          if hydrogenID == hydrogenList[0] {
            return orbitalList[0]
          } else {
            return orbitalList[1]
          }
        }
      default:
        fatalError("Unexpected size for hydrogen list.")
      }
    }
    
    var output = PassivationResult()
    func appendBond(atomID: UInt32, hydrogen: Atom) {
      let hydrogenID = atoms.count + output.insertedAtoms.count
      output.insertedAtoms.append(hydrogen)
      
      let bond = SIMD2(atomID, UInt32(hydrogenID))
      output.insertedBonds.append(bond)
    }
    
    for hydrogenID in input.hydrogensToAtomsMap.indices {
      let atomList = input.hydrogensToAtomsMap[hydrogenID]
      guard atomList.count > 0 else {
        // This collision was resolved.
        continue
      }
      
      switch atomList.count {
      case 1:
        let atomID = atomList[0]
        let orbital = createHydrogenOrbital(
          atomID: atomID,
          hydrogenID: UInt32(hydrogenID))
        
        let hydrogen = createHydrogen(
          atomID: atomID,
          orbital: orbital)
        appendBond(
          atomID: atomID,
          hydrogen: hydrogen)
      case 2:
        let siteCenter = hydrogenSiteCenter(atomList)
        for atomID in atomList {
          let atom = atoms[Int(atomID)]
          let delta = siteCenter - atom.position
          
          var bestOrbital: SIMD3<Float>?
          var bestScore: Float = -.greatestFiniteMagnitude
          let orbitalList = orbitalLists[Int(atomID)]
          for orbital in orbitalList {
            let score = (orbital * delta).sum()
            if score > bestScore {
              bestOrbital = orbital
              bestScore = score
            }
          }
          guard let bestOrbital else {
            fatalError("Could not find the best orbital.")
          }
          
          let hydrogen = createHydrogen(
            atomID: atomID,
            orbital: bestOrbital)
          appendBond(
            atomID: atomID,
            hydrogen: hydrogen)
        }
      default:
        fatalError("3/4-way collisions should have been caught.")
      }
    }
    return output
  }
  #endif
}
