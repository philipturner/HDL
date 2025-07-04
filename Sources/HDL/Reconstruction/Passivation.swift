//
//  Passivation.swift
//  HDL
//
//  Created by Philip Turner on 7/2/25.
//

// TODO: Eliminate this from the library, and migrate it into the test suite.

struct PassivationImpl {
  static func compile(topology: inout Topology) {
    func createHydrogen(
      atomID: UInt32,
      orbital: SIMD3<Float>
    ) -> Atom {
      let atom = topology.atoms[Int(atomID)]
      let atomElement = Element(rawValue: atom.atomicNumber)
      guard let atomElement else {
        fatalError("This should never happen.")
      }
      
      var bondLength = atomElement.covalentRadius
      bondLength += Element.hydrogen.covalentRadius
      
      let position = atom.position + bondLength * orbital
      return Atom(position: position, element: .hydrogen)
    }
    
    let orbitalLists = topology.nonbondingOrbitals()
    
    var insertedAtoms: [Atom] = []
    var insertedBonds: [SIMD2<UInt32>] = []
    for atomID in topology.atoms.indices {
      let orbitalList = orbitalLists[atomID]
      for orbital in orbitalList {
        let hydrogen = createHydrogen(
          atomID: UInt32(atomID),
          orbital: orbital)
        let hydrogenID = topology.atoms.count + insertedAtoms.count
        insertedAtoms.append(hydrogen)
        
        let bond = SIMD2(
          UInt32(atomID),
          UInt32(hydrogenID))
        insertedBonds.append(bond)
      }
    }
    topology.insert(atoms: insertedAtoms)
    topology.insert(bonds: insertedBonds)
  }
}
