//
//  HydrogenSites.swift
//  HDL
//
//  Created by Philip Turner on 7/4/25.
//

struct HydrogenSiteMap {
  var atomsToHydrogensMap: [[UInt32]] = []
  var hydrogensToAtomsMap: [[UInt32]] = []
  var hydrogenSiteCenters: [SIMD3<Float>] = []
}

extension Compilation {
  // A reduced form each hydrogen atom, with the 4th vector slot storing
  // the index of the carbon that spawned it.
  //
  // Inputs: material -> bond length
  //         topology.atoms
  //         topology.bonds -> orbitals
  private func createHydrogenData(
    orbitalLists: [Topology.OrbitalStorage]
  ) -> [SIMD4<Float>] {
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
    
    // Limit for a 10 micron shift: 5.0 pm
    // Limit for a  2 micron shift: 0.9 pm
    // Limit for a    500 nm shift: 0.13 pm
    // Limit for a    100 nm shift: 0.026 pm
    //
    // Choice based on the data: 2.7 pm
    // Most sensible choice from first principles: 10 pm
    return matcher.match(
      hydrogenAtoms, algorithm: .absoluteRadius(0.010))
  }
  
  private static func filter(
    matches: [Topology.MatchStorage],
    hydrogenData: [SIMD4<Float>]
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
  
  private static func createAtomList(
    match: Topology.MatchStorage,
    hydrogenData: [SIMD4<Float>]
  ) -> [UInt32] {
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
  
  private static func createSiteCenter(
    match: Topology.MatchStorage,
    hydrogenData: [SIMD4<Float>]
  ) -> SIMD3<Float> {
    var sum: SIMD3<Float> = .zero
    for j in match {
      let data = hydrogenData[Int(j)]
      let position = unsafeBitCast(data, to: SIMD3<Float>.self)
      sum += position
    }
    guard match.count > 0 else {
      fatalError("Attempted to divide by zero.")
    }
    return sum / Float(match.count)
  }
  
  // Inputs:  material -> hydrogen data -> bond length
  //          topology.atoms -> hydrogen data
  //          topology.bonds -> hydrogen data -> orbitals
  // Outputs: atomsToHydrogensMap
  //          hydrogensToAtomsMap (length determined here)
  func createHydrogenSites(
    orbitalLists: [Topology.OrbitalStorage]
  ) -> HydrogenSiteMap {
    let hydrogenData = createHydrogenData(
      orbitalLists: orbitalLists)
    let rawMatches = Self.createMatches(
      hydrogenData: hydrogenData)
    let filteredMatches = Self.filter(
      matches: rawMatches,
      hydrogenData: hydrogenData)
    
    var output = HydrogenSiteMap()
    output.atomsToHydrogensMap = Array(
      repeating: [],
      count: atoms.count)
    for match in filteredMatches {
      let siteCenter = Self.createSiteCenter(
        match: match,
        hydrogenData: hydrogenData)
      output.hydrogenSiteCenters.append(siteCenter)
      
      // Integrate the atom list into the map.
      let atomList = Self.createAtomList(
        match: match,
        hydrogenData: hydrogenData)
      let hydrogenID = UInt32(
        output.hydrogensToAtomsMap.count)
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
