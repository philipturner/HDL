//
//  Compilation.swift
//  HDL
//
//  Created by Philip Turner on 7/2/25.
//

struct Compilation {
  var atoms: [SIMD4<Float>]
  let material: MaterialType
  
  init(
    atoms: [SIMD4<Float>],
    material: MaterialType
  ) {
    self.atoms = atoms
    self.material = material
  }
  
  // TODO: Refactor the data mutations / state variables throughout this
  // algorithm, making it easier to work with and eventually modify.
  // - Need to refactor this without introducing a performance regression.
  //
  // Implementation Plan
  //
  // Phase I:   Isolate 'resolveTwoWayCollisions', giving it sight of only
  //            the (mutable) atoms <-> hydrogens maps.
  // Phase II:  Isolate 'createHydrogenBonds', attempt to remove the dependency
  //            on the two maps. This is not yet a complete implementation of
  //            a 'Passivation' API.
  // Phase III: Minimize the persistence of the two maps across compilation
  //            steps.
  // Phase IV:  Eliminate passivation from the library code.
  // Phase V:   Fix up the handling of 'topology.bonds'.
  // Phase VI:  Remaining TODOs to the internal code for surface reconstruction.
  //            Clean up the annotations to the function signatures.
  mutating func compile() -> [SIMD2<UInt32>] {
    removePathologicalAtoms()
    
    // Loop over this a few times (typically less than 10).
    for _ in 0..<100 {
      let carbonSites = createCarbonSites()
      let orbitalLists = createOrbitals(
        bonds: carbonSites.bonds)
      let hydrogenSites = createHydrogenSites(
        orbitalLists: orbitalLists)
      
      // Check whether there are still 3-way collisions.
      if hydrogenSites.hydrogensToAtomsMap.contains(where: { $0.count > 2 }) {
        let plugs = plugThreeWayCollisions(
          hydrogensToAtomsMap: hydrogenSites.hydrogensToAtomsMap,
          orbitalLists: orbitalLists)
        self.atoms += plugs
      } else {
        var dimerProcessor = DimerProcessor()
        dimerProcessor.atomsToHydrogensMap = hydrogenSites.atomsToHydrogensMap
        dimerProcessor.hydrogensToAtomsMap = hydrogenSites.hydrogensToAtomsMap
        
        let hydrogenChains = dimerProcessor.createHydrogenChains(
          centerTypes: carbonSites.centerTypes)
        let dimerBonds = dimerProcessor.destroyCollisions(
          hydrogenChains: hydrogenChains)
        return carbonSites.bonds + dimerBonds
      }
    }
    
    fatalError("Could not resolve 3-way collisions.")
  }
}

// MARK: - Utilities

extension Compilation {
  // TODO: Modify the code to use the bulk lattice constant instead of the
  // covalent radii of the atoms. Try to eliminate 'createBondLength()' from
  // every place where it appears.
  
  // Place hydrogens at the C-C bond length instead of the C-H bond length.
  //
  // Inputs:  material
  // Outputs: Float
  func createBondLength() -> Float {
//    var bondLength: Float
//    switch material {
//    case .elemental(let element):
//      bondLength = 2 * element.covalentRadius
//    case .checkerboard(let element, let element2):
//      bondLength = element.covalentRadius + element2.covalentRadius
//    }
//    return bondLength
    
    var bulkBondLength = Constant(.square) { material }
    bulkBondLength *= Float(3).squareRoot() / 4
    return bulkBondLength
  }
  
  func createOrbitals(
    bonds: [SIMD2<UInt32>]
  ) -> [Topology.OrbitalStorage] {
    var output = Topology()
    output.atoms = atoms
    output.bonds = bonds
    return output.nonbondingOrbitals()
  }
  
  func createAtomMatches() -> [Topology.MatchStorage] {
    var topology = Topology()
    topology.atoms = atoms
    
    let bondLength = createBondLength()
    
    // Original problem, caused by inexact bond length:
    //
    // 1.037 - aluminum phosphide fails
    // 1.016 - surface reconstruction reproducers fail
    // 1.009 - hundreds of tests fail
    //
    // Switching to lattice-aligned covalent bond length:
    //
    // 1.001    - no tests fail
    // 1 + 3e-6 - no tests fail
    // 1 + 2e-6 - surface reconstruction reproducers fail
    // 1 + 1e-6 - a test crashes
    //
    // Rationale for new radius:
    //
    // surface reconstruction reproducer: 3.57 nm
    // world dimension: 256 nm
    // conservative overestimate: 1 μm
    //
    // (256 / 3.57) * 2e-6 = 1.4e-4
    // (1000 / 3.57) * 2e-6 = 5.6e-4
    // Both of these numbers are smaller than 1e-2.
    //
    // 0.01 * 0.154 nm = 1.5 pm
    return topology.match(
      atoms, algorithm: .absoluteRadius(bondLength * 1.008))
  }
  
  // Remove methyl groups and floating atoms from the list.
  //
  // Inputs:  material -> bond length
  //          topology.atoms
  // Outputs: topology.atoms (remove)
  private mutating func removePathologicalAtoms() {
    // Loop over this a few times (typically less than 10).
    for _ in 0..<100 {
      let matches = createAtomMatches()
      
      var removedAtoms: [UInt32] = []
      for i in atoms.indices {
        let match = matches[i]
        if match.count > 5 {
          fatalError("Unexpected situation: match count > 5")
        } else if match.count > 2 {
          // 5 matches -> quaternary
          // 4 matches -> tertiary
          // 3 matches -> secondary
        } else if match.count > 0 {
          // 2 matches -> primary
          // 1 match   -> methane
          removedAtoms.append(UInt32(i))
        } else {
          // 0 matches -> impossible
          fatalError("This should never happen.")
        }
      }
      
      if removedAtoms.count > 0 {
        var topology = Topology()
        topology.atoms = atoms
        topology.remove(atoms: removedAtoms)
        self.atoms = topology.atoms
      } else {
        return
      }
    }
    
    fatalError("Could not remove pathological atoms.")
  }
}
