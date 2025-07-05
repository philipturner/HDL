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
  mutating func compile() -> [SIMD2<UInt32>] {
    removePathologicalAtoms()
    
    // Loop over this a few times (typically less than 10).
    for _ in 0..<100 {
      // TODO: Isolate and clarify this state mutation. Perhaps make Topology
      // transient (which it can be, at no cost) and store/return the bonds
      // separately.
      let carbonSites = createCarbonSites()
      let hydrogenSites = createHydrogenSites(
        bonds: carbonSites.bonds)
      
      // Check whether there are still 3-way collisions.
      if hydrogenSites.hydrogensToAtomsMap.contains(where: { $0.count > 2 }) {
        resolveThreeWayCollisions(
          bonds: carbonSites.bonds,
          hydrogensToAtomsMap: hydrogenSites.hydrogensToAtomsMap)
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
  // Place hydrogens at the C-C bond length instead of the C-H bond length.
  //
  // Inputs:  material
  // Outputs: Float
  func createBondLength() -> Float {
    var bondLength: Float
    switch material {
    case .elemental(let element):
      bondLength = 2 * element.covalentRadius
    case .checkerboard(let element, let element2):
      bondLength = element.covalentRadius + element2.covalentRadius
    }
    return bondLength
  }
  
  // TODO: It might be possible to merge several calls to this. Don't make the
  // change until you've got the tests working again, because it will alter
  // performance.
  //
  // Another blocker to this change: the minor relocation of a couple lines
  // of code inside a loop.
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
    return topology.match(
      atoms, algorithm: .absoluteRadius(bondLength * 1.1))
  }
  
  // TODO: Perhaps change this to be non-mutating, so the state changes are
  // more explicit.
  //
  // Remove atoms with less than two covalent bonds.
  //
  // Inputs:  material -> bond length
  //          topology.atoms
  // Outputs: topology.atoms (remove)
  private mutating func removePathologicalAtoms() {
    // Loop over this a few times (typically less than 10).
    var converged = false
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
        converged = true
        break
      }
    }
    guard converged else {
      fatalError("Could not remove pathological atoms.")
    }
    
    // TODO: Use the 'nil' return convention to restructure unbounded loops.
  }
}
