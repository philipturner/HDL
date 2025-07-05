//
//  CarbonSites.swift
//  HDL
//
//  Created by Philip Turner on 7/5/25.
//

struct CarbonSiteMap {
  var bonds: [SIMD2<UInt32>] = []
  var centerTypes: [UInt8] = []
}

extension Compilation {
  // Remove methyl groups and floating atoms from the list.
  mutating func removeMethylSites() {
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
    
    fatalError("Could not remove methyl sites.")
  }
  
  func createCarbonSites() -> CarbonSiteMap {
    let matches = createAtomMatches()
    
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
    return output
  }
}
