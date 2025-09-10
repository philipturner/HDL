//
//  HydrogenSites.swift
//  HDL
//
//  Created by Philip Turner on 7/4/25.
//

import QuartzCore

struct HydrogenSiteMap {
  var atomsToHydrogensMap: [[UInt32]] = []
  var hydrogensToAtomsMap: [[UInt32]] = []
  var hydrogenSiteCenters: [SIMD3<Float>] = []
}

extension Compilation {
  func createHydrogenSites(
    i: Int,
    bonds: [SIMD2<UInt32>]
  ) -> HydrogenSiteMap {
    let checkpoint0 = CACurrentMediaTime()
    let orbitalLists = createOrbitalLists(
      bonds: bonds)
    let checkpoint1 = CACurrentMediaTime()
    display(checkpoint0, checkpoint1, "\(i) - createHydrogenSites/createOrbitalLists")
    
    let hydrogenData = createHydrogenData(
      orbitalLists: orbitalLists)
    let checkpoint2 = CACurrentMediaTime()
    display(checkpoint1, checkpoint2, "\(i) - createHydrogenSites/createHydrogenData")
    
    let rawMatches = Self.createMatches(
      hydrogenData: hydrogenData)
    let checkpoint3 = CACurrentMediaTime()
    display(checkpoint2, checkpoint3, "\(i) - createHydrogenSites/createMatches")
    
    let filteredMatches = Self.filter(
      matches: rawMatches,
      hydrogenData: hydrogenData)
    let checkpoint4 = CACurrentMediaTime()
    display(checkpoint3, checkpoint4, "\(i) - createHydrogenSites/filter")
    
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
    let checkpoint5 = CACurrentMediaTime()
    display(checkpoint4, checkpoint5, "\(i) - createHydrogenSites/HydrogenSiteMap")
    
    // Sort each hydrogen list, in place.
    for j in atoms.indices {
      output.atomsToHydrogensMap[Int(j)].sort()
    }
    let checkpoint6 = CACurrentMediaTime()
    display(checkpoint5, checkpoint6, "\(i) - createHydrogenSites/sort")
    
    return output
  }
  
  private func createOrbitalLists(
    bonds: [SIMD2<UInt32>]
  ) -> [Topology.OrbitalStorage] {
    var output = Topology()
    output.atoms = atoms
    output.bonds = bonds
    return output.nonbondingOrbitals()
  }
  
  // A compact data structure for a hydrogen site.
  // slot 0: position.x
  // slot 1: position.y
  // slot 2: position.z
  // slot 3: index of the carbon that spawned it
  private func createHydrogenData(
    orbitalLists: [Topology.OrbitalStorage]
  ) -> [SIMD4<Float>] {
    let bondLength = createBondLength()
    
    var output: [SIMD4<Float>] = []
    for i in atoms.indices {
      let atom = atoms[i]
      let orbitalList = orbitalLists[i]
      for orbital in orbitalList {
        let atomPosition = unsafeBitCast(atom, to: SIMD3<Float>.self)
        let hydrogenPosition = atomPosition + bondLength * orbital
        let encodedID = Float(bitPattern: UInt32(i))
        output.append(SIMD4(hydrogenPosition, encodedID))
      }
    }
    
    return output
  }
  
  private static func createMatches(
    hydrogenData: [SIMD4<Float>]
  ) -> [Topology.MatchStorage] {
    // Give all the atoms a valid atomic number, for safety when invoking the
    // match procedure.
    let atoms = hydrogenData.map {
      var atom = $0
      atom.w = 1
      return atom
    }
    var topology = Topology()
    topology.atoms = atoms
    
    // Limit for a 10 micron shift: 5.0 pm
    // Limit for a  2 micron shift: 0.9 pm
    // Limit for a    500 nm shift: 0.13 pm
    // Limit for a    100 nm shift: 0.026 pm
    //
    // Choice based on the data: 2.7 pm
    // Most sensible choice from first principles: 10 pm
    return topology.match(
      atoms, algorithm: .absoluteRadius(0.010))
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
    
    // Sort the atom list, in place.
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
}
