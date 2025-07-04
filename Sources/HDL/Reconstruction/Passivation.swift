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
}
