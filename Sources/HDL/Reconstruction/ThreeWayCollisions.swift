//
//  ThreeWayCollisions.swift
//  
//
//  Created by Philip Turner on 6/11/24.
//

extension Compilation {
  // Adds carbon atoms to places where 3 hydrogens collide.
  //
  // Inputs:  material
  //          topology.atoms
  //          topology.bonds -> orbitals
  //          hydrogensToAtomsMap
  // Outputs: topology.atoms (remove)
  func plugThreeWayCollisions(
    siteMap: HydrogenSiteMap
  ) -> [Atom] {
    var insertedAtoms: [Atom] = []
    for hydrogenSiteID in siteMap.hydrogensToAtomsMap.indices {
      let atomList = siteMap.hydrogensToAtomsMap[hydrogenSiteID]
      guard atomList.count == 3 else {
        continue
      }
      
      // Iterate over all 3 atoms in the collision.
      var atomicNumbersDict: [UInt8: Int] = [:]
      for atomID in atomList {
        let atom = atoms[Int(atomID)]
        let atomicNumber = atom.atomicNumber
        if atomicNumbersDict[atomicNumber] == nil {
          atomicNumbersDict[atomicNumber] = 1
        } else {
          atomicNumbersDict[atomicNumber]! += 1
        }
      }
      
      // Find the atomic number with the greatest frequency in the neighbors
      // to the collision site.
      func createDominantAtomicNumber() -> UInt8 {
        var maxAtomicNumber: UInt8 = .zero
        var maxAtomicNumberCount: Int = .zero
        for atomicNumber in atomicNumbersDict.keys {
          let count = atomicNumbersDict[atomicNumber]!
          if count > maxAtomicNumberCount {
            maxAtomicNumberCount = count
            maxAtomicNumber = atomicNumber
          }
        }
        return maxAtomicNumber
      }
      let dominantAtomicNumber = createDominantAtomicNumber()
      
      // Fill in with atoms of the assigned material.
      var chosenAtomicNumber: UInt8
      switch material {
      case .elemental(let element):
        chosenAtomicNumber = element.rawValue
      case .checkerboard(let element1, let element2):
        if dominantAtomicNumber == element1.rawValue {
          chosenAtomicNumber = element2.rawValue
        } else if dominantAtomicNumber == element2.rawValue {
          chosenAtomicNumber = element1.rawValue
        } else {
          fatalError("Could not resolve identity of inserted atom in checkerboard structure.")
        }
      }
      
      let position = siteMap.hydrogenSiteCenters[hydrogenSiteID]
      let atom = Atom(
        position: position,
        atomicNumber: chosenAtomicNumber)
      insertedAtoms.append(atom)
    }
    return insertedAtoms
  }
}
