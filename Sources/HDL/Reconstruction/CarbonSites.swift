//
//  CarbonSites.swift
//  HDL
//
//  Created by Philip Turner on 7/5/25.
//

import QuartzCore

struct CarbonSiteMap {
  var bonds: [SIMD2<UInt32>] = []
  var centerTypes: [UInt8] = []
}

extension Compilation {
  // Remove methyl groups and floating atoms from the list.
  mutating func removeMethylSites() {
    // Loop over this a few times (typically less than 10).
    for i in 0..<100 {
      let checkpoint0 = CACurrentMediaTime()
      let matches = createAtomMatches()
      let checkpoint1 = CACurrentMediaTime()
      display(checkpoint0, checkpoint1, "\(i) - removeMethylSites/createAtomMatches")
      
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
      let checkpoint2 = CACurrentMediaTime()
      display(checkpoint1, checkpoint2, "\(i) - removeMethylSites/removedAtoms")
      
      if removedAtoms.count > 0 {
        var topology = Topology()
        topology.atoms = atoms
        topology.remove(atoms: removedAtoms)
        self.atoms = topology.atoms
        
        let checkpoint3 = CACurrentMediaTime()
        display(checkpoint2, checkpoint3, "\(i) - removeMethylSites/topology.remove(atoms:)")
      } else {
        return
      }
    }
    
    fatalError("Could not remove methyl sites.")
  }
  
  func createCarbonSites(i: Int) -> CarbonSiteMap {
    let checkpoint0 = CACurrentMediaTime()
    let matches = createAtomMatches()
    let checkpoint1 = CACurrentMediaTime()
    display(checkpoint0, checkpoint1, "\(i) - createCarbonSites/createAtomMatches")
    
    var output = CarbonSiteMap()
    for i in atoms.indices {
      let match = matches[i]
      if match.count > 5 {
        fatalError("Unexpected situation: match count > 5")
      } else if match.count > 2 {
        for j in match where i < j {
          let bond = SIMD2(UInt32(i), j)
          output.bonds.append(bond)
        }
      } else {
        fatalError("Pathological atoms should be removed.")
      }
      
      let centerType = UInt8(match.count - 1)
      output.centerTypes.append(centerType)
    }
    let checkpoint2 = CACurrentMediaTime()
    display(checkpoint1, checkpoint2, "\(i) - createCarbonSites/CarbonSiteMap")
    
    return output
  }
}
