//
//  Compilation.swift
//  HDL
//
//  Created by Philip Turner on 7/2/25.
//

import QuartzCore

private func display(_ start: Double, _ end: Double, _ name: String)  {
  let seconds = end - start
  let milliseconds = seconds * 1e3
  let repr = String(format: "%.3f", milliseconds)
  print("\(name):", repr, "ms")
}

struct Compilation {
  var atoms: [SIMD4<Float>]
  let material: MaterialType
  
  mutating func compile() -> [SIMD2<UInt32>] {
    let start = CACurrentMediaTime()
    removeMethylSites()
    let end = CACurrentMediaTime()
    display(start, end, "removeMethylSites")
    
    // Loop over this a few times (typically less than 10).
    for i in 0..<100 {
      let checkpoint0 = CACurrentMediaTime()
      let carbonSites = createCarbonSites()
      let checkpoint1 = CACurrentMediaTime()
      display(checkpoint0, checkpoint1, "\(i) - createCarbonSites")
      
      let hydrogenSites = createHydrogenSites(
        bonds: carbonSites.bonds)
      let checkpoint2 = CACurrentMediaTime()
      display(checkpoint1, checkpoint2, "\(i) - createHydrogenSites")
      
      // Check whether there are still 3-way collisions.
      if hydrogenSites.hydrogensToAtomsMap.contains(where: { $0.count > 2 }) {
        let plugs = plugThreeWayCollisions(
          siteMap: hydrogenSites)
        let checkpoint3 = CACurrentMediaTime()
        display(checkpoint2, checkpoint3, "\(i) - plugThreeWayCollisions")
        
        self.atoms += plugs
        let checkpoint4 = CACurrentMediaTime()
        display(checkpoint3, checkpoint4, "\(i) - self.atoms += plugs")
      } else {
        var dimerProcessor = DimerProcessor()
        dimerProcessor.atomsToHydrogensMap = hydrogenSites.atomsToHydrogensMap
        dimerProcessor.hydrogensToAtomsMap = hydrogenSites.hydrogensToAtomsMap
        
        let hydrogenChains = dimerProcessor.createHydrogenChains(
          centerTypes: carbonSites.centerTypes)
        let checkpoint5 = CACurrentMediaTime()
        display(checkpoint2, checkpoint5, "\(i) - createHydrogenChains")
        
        let dimerBonds = dimerProcessor.destroyCollisions(
          hydrogenChains: hydrogenChains)
        let checkpoint6 = CACurrentMediaTime()
        display(checkpoint5, checkpoint6, "\(i) - destroyCollisions")
        
        let output = carbonSites.bonds + dimerBonds
        let checkpoint7 = CACurrentMediaTime()
        display(checkpoint6, checkpoint7, "\(i) - carbonSites.bonds + dimerBonds")
        
        return output
      }
    }
    
    fatalError("Could not resolve 3-way collisions.")
  }
  
  // Use a bond length that perfectly agrees with the lattice constant.
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
}
