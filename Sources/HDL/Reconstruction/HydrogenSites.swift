//
//  HydrogenSites.swift
//  HDL
//
//  Created by Philip Turner on 7/4/25.
//

struct HydrogenSiteMap {
  var atomsToHydrogensMap: [[UInt32]] = []
  var hydrogensToAtomsMap: [[UInt32]] = []
}

extension Compilation {
  // A reduced form each hydrogen atom, with the 4th vector slot storing
  // the index of the carbon that spawned it.
  //
  // Inputs: material -> bond length
  //         topology.atoms
  //         topology.bonds -> orbitals
  private func createHydrogenData() -> [SIMD4<Float>] {
    let orbitalLists = createTopology().nonbondingOrbitals()
    let bondLength = createBondLength()
    
    var output: [SIMD4<Float>] = []
    for i in atoms.indices {
      let atom = atoms[i]
      let orbitalList = orbitalLists[i]
      for orbital in orbitalList {
        let position = atom.position + bondLength * orbital
        let encodedID = Float(bitPattern: UInt32(i))
        output.append(SIMD4(position, encodedID))
      }
    }
    
    return output
  }
  
  private static func createMatches(
    hydrogenData: [SIMD4<Float>]
  ) -> [Topology.MatchStorage] {
    let hydrogenAtoms = hydrogenData.map {
      var atom = $0
      atom.w = 1
      return atom
    }
    
    // Create a transient topology to de-duplicate the hydrogens and merge
    // references between them.
    var matcher = Topology()
    matcher.atoms = hydrogenAtoms
    return matcher.match(
      hydrogenAtoms, algorithm: .absoluteRadius(0.050))
  }
  
  // Inputs:  material -> hydrogen data -> bond length
  //          topology.atoms -> hydrogen data
  //          topology.bonds -> hydrogen data -> orbitals
  // Outputs: atomsToHydrogensMap
  //          hydrogensToAtomsMap (length determined here)
  func createHydrogenSites() -> HydrogenSiteMap {
    let hydrogenData = createHydrogenData()
    
    func filter(
      matches: [Topology.MatchStorage]
    ) -> [Topology.MatchStorage] {
      var output: [Topology.MatchStorage] = []
      for i in hydrogenData.indices {
        let match = matches[i]
        
        func isCompatible() -> Bool {
          if match.count > 1 {
            for j in match where i != j {
              if i > j {
                return false
              }
            }
          }
          return true
        }
        
        if isCompatible() {
          output.append(match)
        }
      }
      return output
    }
    
    var output = HydrogenSiteMap()
    output.atomsToHydrogensMap = Array(
      repeating: [],
      count: atoms.count)
    
    let matches = Self.createMatches(hydrogenData: hydrogenData)
    let filteredMatches = filter(matches: matches)
    for match in filteredMatches {
      func createAtomList() -> [UInt32] {
        var atomList: [UInt32] = []
        for j in match {
          let data = hydrogenData[Int(j)]
          let atomID = data.w.bitPattern
          atomList.append(atomID)
        }
        if atomList.count >= 4 {
          fatalError("4-way collisions are not handled yet.")
        }
        
        // Sort the atom list, in-place.
        atomList.sort()
        return atomList
      }
      
      // Integrate the atom list into the map.
      let atomList = createAtomList()
      let hydrogenID = UInt32(output.hydrogensToAtomsMap.count)
      output.hydrogensToAtomsMap.append(atomList)
      
      // Mutate the hydrogen list.
      for j in atomList {
        // Appending to the array in-place has a measurable performance
        // improvement, compared to extract + modify + insert.
        output.atomsToHydrogensMap[Int(j)].append(hydrogenID)
      }
    }
    
    // Sort each hydrogen list, in-place.
    for j in atoms.indices {
      output.atomsToHydrogensMap[Int(j)].sort()
    }
    return output
  }
}
