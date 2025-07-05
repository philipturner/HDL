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
  // Form all center atom bonds in the lattice interior.
  //
  // Returns the center type of each atom.
  //
  // Inputs:  material -> bond length
  //          topology.atoms
  // Outputs: topology.bonds (insert)
  //          center types
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
        
        // TODO: Migrate this to the destination while checking for a
        // performance regression in Reconstruction tests.
        let centerType = UInt8(match.count - 1)
        output.centerTypes.append(centerType)
      } else {
        fatalError("Pathological atoms should be removed.")
      }
      
      // destination
    }
    return output
  }
}
