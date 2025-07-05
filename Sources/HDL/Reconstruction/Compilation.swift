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
      let hydrogenSites = createHydrogenSites(
        bonds: carbonSites.bonds)
      
      // Check whether there are still 3-way collisions.
      if hydrogenSites.hydrogensToAtomsMap.contains(where: { $0.count > 2 }) {
        let plugs = plugThreeWayCollisions(
          siteMap: hydrogenSites)
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
  
  // Place hydrogens at the C-C bond length instead of the C-H bond length.
  //
  // Inputs:  material
  // Outputs: Float
  func createBondLength() -> Float {
    var bulkBondLength = Constant(.square) { material }
    bulkBondLength *= Float(3).squareRoot() / 4
    return bulkBondLength
  }
  
  func createAtomMatches() -> [Topology.MatchStorage] {
    var topology = Topology()
    topology.atoms = atoms
    
    // Limit for a 10 micron shift: 1.008    -> 1.2 pm
    // Limit for a  2 micron shift: 1.0007   -> 0.11 pm
    // Limit for a    500 nm shift: 1.00032  -> 0.049 pm
    // Limit for a    100 nm shift: 1.000056 -> 0.008 pm
    //
    // Choice based on the data: 1.0021 (0.3 pm)
    // Most sensible choice from first principles: 1.03 (4.6 pm)
    let bondLength = createBondLength()
    return topology.match(
      atoms, algorithm: .absoluteRadius(bondLength * 1.03))
  }
  
  // Remove methyl groups and floating atoms from the list.
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
