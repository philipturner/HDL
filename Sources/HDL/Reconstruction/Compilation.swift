//
//  Compilation.swift
//  HDL
//
//  Created by Philip Turner on 7/2/25.
//

struct Compilation {
  var atoms: [SIMD4<Float>]
  var bonds: [SIMD2<UInt32>] // stop toward making this transient, not stored
  let material: MaterialType
  
  init(
    atoms: [SIMD4<Float>],
    material: MaterialType
  ) {
    self.atoms = atoms
    self.bonds = []
    self.material = material
  }
  
  func createOrbitals() -> [Topology.OrbitalStorage] {
    var output = Topology()
    output.atoms = atoms
    output.bonds = bonds
    return output.nonbondingOrbitals()
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
  mutating func compile() {
    removePathologicalAtoms()
    
    // Loop over this a few times (typically less than 10).
    var converged = false
    for _ in 0..<100 {
      // TODO: Isolate and clarify this state mutation. Perhaps make Topology
      // transient (which it can be, at no cost) and store/return the bonds
      // separately.
      let centerTypes = createBulkAtomBonds()
      let siteMap = createHydrogenSites()
      
      // Check whether there are still 3-way collisions.
      if siteMap.hydrogensToAtomsMap.contains(where: { $0.count > 2 }) {
        // Add center atoms to problematic sites.
        resolveThreeWayCollisions(
          hydrogensToAtomsMap: siteMap.hydrogensToAtomsMap)
        
        // Reverse the actions from the start of this iteration.
        //
        // TODO: Isolate and clarify this state mutation.
        self.bonds = []
      } else {
        var dimerProcessor = DimerProcessor()
        dimerProcessor.atomsToHydrogensMap = siteMap.atomsToHydrogensMap
        dimerProcessor.hydrogensToAtomsMap = siteMap.hydrogensToAtomsMap
        
        let hydrogenChains = dimerProcessor.createHydrogenChains(
          centerTypes: centerTypes)
        let insertedBonds = dimerProcessor.destroyCollisions(
          hydrogenChains: hydrogenChains)
        self.bonds += insertedBonds
        
        // We could return the bonds here.
        
        converged = true
        break
      }
    }
    guard converged else {
      fatalError("Could not resolve 3-way collisions.")
    }
    
    // We can just return 'nil' here to signal that bonds could not be
    // materialized. Or perhaps group together one of 3 internal convergence
    // errors, which are beautifully shown in an error message.
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
  
  private func createAtomMatches() -> [Topology.MatchStorage] {
    var topology = Topology()
    topology.atoms = atoms
    
    let bondLength = createBondLength()
    return topology.match(
      atoms, algorithm: .absoluteRadius(bondLength * 1.1))
  }
  
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
  }
  
  // Form all center atom bonds in the lattice interior.
  //
  // Returns the center type of each atom.
  //
  // Inputs:  material -> bond length
  //          topology.atoms
  // Outputs: topology.bonds (insert)
  //          center types
  private mutating func createBulkAtomBonds() -> [UInt8] {
    let matches = createAtomMatches()
    
    var insertedBonds: [SIMD2<UInt32>] = []
    var centerTypes: [UInt8] = []
    for i in atoms.indices {
      let match = matches[i]
      if match.count > 5 {
        fatalError("Unexpected situation: match count > 5")
      } else if match.count > 2 {
        for j in match where i < j {
          let bond = SIMD2(UInt32(i), j)
          insertedBonds.append(bond)
        }
        centerTypes.append(UInt8(match.count - 1))
      } else {
        fatalError("Pathological atoms should be removed.")
      }
    }
    self.bonds += insertedBonds
    
    return centerTypes
  }
}
