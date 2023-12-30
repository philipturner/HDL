import XCTest
import HDL
import Numerics

final class TopologyTests: XCTestCase {
  func testInsertRemove() throws {
    let lonsdaleite = Lonsdaleite()
    var topology = Topology()
    
    topology.insert(atoms: lonsdaleite.atoms)
    XCTAssertEqual(topology.atoms, lonsdaleite.atoms)
    XCTAssertEqual(topology.bonds, [])
    
    topology.remove(atoms: lonsdaleite.atoms.indices.map(UInt32.init))
    XCTAssertEqual(topology.atoms, [])
    XCTAssertEqual(topology.bonds, [])
    
    topology.insert(atoms: lonsdaleite.atoms)
    topology.insert(bonds: lonsdaleite.bonds)
    XCTAssertEqual(topology.atoms, lonsdaleite.atoms)
    XCTAssertEqual(topology.bonds, lonsdaleite.bonds)
    
    topology.remove(bonds: lonsdaleite.bonds.indices.map(UInt32.init))
    XCTAssertEqual(topology.atoms, lonsdaleite.atoms)
    XCTAssertEqual(topology.bonds, [])
    
    topology.insert(bonds: lonsdaleite.bonds)
    XCTAssertEqual(topology.atoms, lonsdaleite.atoms)
    XCTAssertEqual(topology.bonds, lonsdaleite.bonds)
    
    topology.remove(atoms: lonsdaleite.atoms.indices.map(UInt32.init))
    XCTAssertEqual(topology.atoms, [])
    XCTAssertEqual(topology.bonds, [])
  }
  
  func testPartialRemove() {
    let lonsdaleite = Lonsdaleite()
    var topology = Topology()
    
    var compactedAtomCount = 0
    var mappedAtomIndices: [Int] = []
    var removedAtomIDs: [UInt32] = []
    var removedBondIDs: [UInt32] = []
    
    var carbonAtoms: [Entity] = []
    var originalCarbonBonds: [SIMD2<UInt32>] = []
    var compactedCarbonBonds: [SIMD2<UInt32>] = []
    
    for atomID in lonsdaleite.atoms.indices {
      let atom = lonsdaleite.atoms[atomID]
      if atom.atomicNumber == 1 {
        mappedAtomIndices.append(-1)
        removedAtomIDs.append(UInt32(atomID))
      } else {
        mappedAtomIndices.append(compactedAtomCount)
        carbonAtoms.append(atom)
        compactedAtomCount += 1
      }
    }
    
    for bondID in lonsdaleite.bonds.indices {
      let bond = lonsdaleite.bonds[bondID]
      let atom1 = lonsdaleite.atoms[Int(bond[0])]
      let atom2 = lonsdaleite.atoms[Int(bond[1])]
      if atom1.atomicNumber == 1 || atom2.atomicNumber == 1 {
        removedBondIDs.append(UInt32(bondID))
        continue
      }
      
      let mappedID1 = mappedAtomIndices[Int(bond[0])]
      let mappedID2 = mappedAtomIndices[Int(bond[1])]
      originalCarbonBonds.append(bond)
      compactedCarbonBonds.append(SIMD2(UInt32(mappedID1),
                                        UInt32(mappedID2)))
    }
    
    // Remove all carbon-hydrogen bonds and ensure the remaining topology
    // appears as expected.
    topology.insert(atoms: lonsdaleite.atoms)
    topology.insert(bonds: lonsdaleite.bonds)
    topology.remove(bonds: removedBondIDs)
    XCTAssertEqual(topology.atoms, lonsdaleite.atoms)
    XCTAssertEqual(topology.bonds, originalCarbonBonds)
    
    topology.remove(atoms: topology.atoms.indices.map(UInt32.init))
    XCTAssertEqual(topology.atoms, [])
    XCTAssertEqual(topology.bonds, [])
    
    // Remove all hydrogen atoms and ensure the remaining carbon-carbon bonds
    // appear as expected.
    topology.insert(atoms: lonsdaleite.atoms)
    topology.insert(bonds: lonsdaleite.bonds)
    topology.remove(atoms: removedAtomIDs)
    XCTAssertEqual(topology.atoms, carbonAtoms)
    XCTAssertEqual(topology.bonds, compactedCarbonBonds)
  }
  
  func testMap() throws {
    let lonsdaleite = Lonsdaleite()
    var topology = Topology()
    topology.insert(atoms: lonsdaleite.atoms)
    topology.insert(bonds: lonsdaleite.bonds)
    
    // bonds -> atoms
    do {
      let map = topology.map(.bonds, to: .atoms)
      XCTAssertEqual(map.count, topology.bonds.count)
      
      for bondID in map.indices {
        let bond = topology.bonds[bondID]
        let slice = map[bondID]
        XCTAssertEqual(slice.count, 2)
        XCTAssertEqual(slice[slice.startIndex + 0], bond[0])
        XCTAssertEqual(slice[slice.startIndex + 1], bond[1])
      }
    }
    
    // atoms -> bonds
    do {
      var bondCoverage = [Bool](
        repeating: false, count: lonsdaleite.bonds.count)
      let map = topology.map(.atoms, to: .bonds)
      XCTAssertEqual(map.count, topology.atoms.count)
      
      for atomID in map.indices {
        let atom = topology.atoms[atomID]
        let slice = map[atomID]
        if atom.atomicNumber == 1 {
          XCTAssertEqual(slice.count, 1)
        } else {
          XCTAssertEqual(slice.count, 4)
        }
        
        for i in slice.indices {
          let bondID = Int(slice[i])
          let bond = topology.bonds[bondID]
          bondCoverage[bondID] = true
          XCTAssert(bond[0] == UInt32(atomID) ||
                    bond[1] == UInt32(atomID))
        }
      }
      
      XCTAssert(bondCoverage.allSatisfy { $0 == true })
    }
    
    // atoms -> atoms
    do {
      var atomCoverage = [Bool](
        repeating: false, count: lonsdaleite.atoms.count)
      let map = topology.map(.atoms, to: .atoms)
      XCTAssertEqual(map.count, topology.atoms.count)
      
      for atomID in map.indices {
        let atom = topology.atoms[atomID]
        let slice = map[atomID]
        if atom.atomicNumber == 1 {
          XCTAssertEqual(slice.count, 1)
        } else {
          XCTAssertEqual(slice.count, 4)
        }
        
        for i in slice.indices {
          let neighborID = Int(slice[i])
          let neighbor = topology.atoms[neighborID]
          atomCoverage[neighborID] = true
          
          if atom.atomicNumber == 1 {
            XCTAssertEqual(neighbor.atomicNumber, 6)
          }
          XCTAssertNotEqual(atomID, neighborID)
        }
      }
      
      XCTAssert(atomCoverage.allSatisfy { $0 == true })
    }
  }
  
  func testSort() throws {
    let lonsdaleite = Lonsdaleite()
    let shuffledMapping = lonsdaleite.atoms.indices.map { $0 }.shuffled()
    var shuffledAtoms = lonsdaleite.atoms
    for originalID in shuffledMapping.indices {
      let reorderedID = shuffledMapping[originalID]
      shuffledAtoms[reorderedID] = lonsdaleite.atoms[originalID]
    }
    
    let shuffledBonds = lonsdaleite.bonds.map {
      let atomID1 = shuffledMapping[Int($0.x)]
      let atomID2 = shuffledMapping[Int($0.y)]
      return SIMD2<UInt32>(UInt32(atomID1),
                           UInt32(atomID2))
    }
    
    var topology = Topology()
    topology.insert(atoms: shuffledAtoms)
    topology.insert(bonds: shuffledBonds)
    
    func checkSelfConsistency() {
      let map = topology.map(.atoms, to: .atoms)
      for atomID in topology.atoms.indices {
        let atom = topology.atoms[atomID]
        let slice = map[atomID]
        if atom.atomicNumber == 1 {
          XCTAssertEqual(slice.count, 1)
        } else {
          XCTAssertEqual(slice.count, 4)
        }
        
        for i in slice.indices {
          XCTAssertNotEqual(atomID, Int(slice[i]))
        }
      }
    }
    
    checkSelfConsistency()
    let topologyReordering = topology.sort()
    checkSelfConsistency()
    
    // Check that the reordering doesn't just copy the input.
    let copyReordering = topology.atoms.indices.map(UInt32.init)
    XCTAssertNotEqual(copyReordering, topologyReordering)
    
    // Check that the atoms are exactly the same.
    let octreeSorter = OctreeSorter(atoms: shuffledAtoms)
    let octreeReordering = octreeSorter.mortonReordering()
    for (lhs, rhs) in zip(topologyReordering, octreeReordering) {
      XCTAssertEqual(lhs, rhs)
    }
    
    // Check that the correct bonds exist in some order.
    let correctBonds = shuffledBonds.map {
      let atomID1 = Int($0.x)
      let atomID2 = Int($0.y)
      return SIMD2<UInt32>(topologyReordering[atomID1],
                           topologyReordering[atomID2])
    }
    var bondsDictionary: [SIMD2<UInt32>: Bool] = [:]
    for reorderedBond in topology.bonds {
      bondsDictionary[reorderedBond] = true
    }
    for correctBond in correctBonds {
      let reversed = SIMD2(correctBond.y, correctBond.x)
      XCTAssert(bondsDictionary[correctBond] != nil ||
                bondsDictionary[reversed] != nil)
    }
    
    // Check that the bonds were rearranged to be sorted in order (>99% of the
    // time, this occurs).
    var different = false
    for (lhs, rhs) in zip(correctBonds, topology.bonds) {
      let reversed = SIMD2(lhs.y, lhs.x)
      if lhs != rhs && reversed != rhs {
        different = true
      }
    }
    XCTAssertTrue(different)
    
    // Sort the bonds, then check they are exactly the same.
    var sortedBonds = correctBonds.map {
      SIMD2($0.min(), $0.max())
    }
    sortedBonds.sort(by: {
      let hash1 = Int($0.y) + Int(1 << 32) * Int($0.x)
      let hash2 = Int($1.y) + Int(1 << 32) * Int($1.x)
      return hash1 < hash2
    })
    XCTAssertEqual(sortedBonds, topology.bonds)
  }
  
  // Test match() against phosphorus-doped silicon.
  func testMatchSilicon() throws {
    let lattice = Lattice<Cubic> { h, k, l in
      Bounds { 3 * (h + k + l) }
      Material { .elemental(.silicon) }
      
      Volume {
        Origin { 0.1 * h }
        Plane { -h }
        
        // The dopants replace sidewall silicons, so they only form 2 covalent
        // bonds with the bulk crystal.
        Replace { .atom(.phosphorus) }
      }
    }
    
    let expectedImpurity: Float = 0.25 / 3
    let siliconCount = lattice.atoms.filter { $0.atomicNumber == 14 }.count
    let phosphorusCount = lattice.atoms.filter { $0.atomicNumber == 15 }.count
    let dopantConcentration =
    Float(phosphorusCount) / Float(siliconCount + phosphorusCount)
    XCTAssertEqual(expectedImpurity, dopantConcentration, accuracy: 0.01)
    
    var topology = Topology()
    topology.insert(atoms: lattice.atoms)
    let matches = topology.match(
      lattice.atoms, algorithm: .covalentBondLength(2.1))
    XCTAssertEqual(matches.indices, lattice.atoms.indices)
    
    for i in matches.indices {
      let range = matches[i]
      XCTAssertGreaterThan(range.count, 0)
      XCTAssertGreaterThan(range.count, 4)
      XCTAssertLessThan(range.count, 30)
      XCTAssertEqual(range[range.startIndex], UInt32(i))
      
      let atom = lattice.atoms[i]
      var selfRadius: Float
      if atom.atomicNumber == 14 {
        selfRadius = Element.silicon.covalentRadius
      } else {
        selfRadius = Element.phosphorus.covalentRadius
      }
      
      var actualMatches: [Int: Bool] = [:]
      for match in range {
        actualMatches[Int(match)] = true
      }
      
      for j in lattice.atoms.indices {
        let neighbor = lattice.atoms[j]
        let delta = atom.position - neighbor.position
        let distance = (delta * delta).sum().squareRoot()
        
        var neighborRadius: Float
        if neighbor.atomicNumber == 14 {
          neighborRadius = Element.silicon.covalentRadius
        } else {
          neighborRadius = Element.phosphorus.covalentRadius
        }
        let searchRadius = 2.1 * (selfRadius + neighborRadius)
        
        if distance < searchRadius {
          XCTAssertNotNil(
            actualMatches[j], "\(i) -> \(j), \(distance / searchRadius)")
        } else {
          XCTAssertNil(
            actualMatches[j], "\(i) -> \(j), \(distance / searchRadius)")
        }
      }
      
      var previousDistance: Float = 0
      for match in range {
        let neighbor = lattice.atoms[Int(match)]
        let delta = atom.position - neighbor.position
        let distance = (delta * delta).sum().squareRoot()
        XCTAssertLessThanOrEqual(previousDistance, distance + 1e-3)
        previousDistance = distance
      }
    }
  }
  
  // Test that a raw carbon lattice produces the same results with
  // .covalentBondLength(1.1) and 2.2x the absolute radius of carbon. Then,
  // test that lonsdaleite with ~1.1-1.5x covalent bond length returns only
  // the C-C and C-H covalent bonds. Finally, perform an asymmetric search
  // using hydrogens detached from an unpassivated lonsdaleite lattice. Form a
  // topology and assert that, after sorting, it is the same as Lonsdaleite().
  func testMatchLonsdaleite() throws {
    let lonsdaleite = Lonsdaleite()
    let carbons = lonsdaleite.atoms.filter { $0.atomicNumber == 6 }
    do {
      var topology = Topology()
      topology.insert(atoms: carbons)
      let matchesCovalent = topology.match(
        topology.atoms, algorithm: .covalentBondLength(1.1))
      XCTAssertEqual(matchesCovalent.count, carbons.count)
      
      for i in matchesCovalent.indices {
        let range = matchesCovalent[i]
        XCTAssertGreaterThan(range.count, 0)
        XCTAssertGreaterThan(range.count, 2)
        XCTAssertLessThan(range.count, 6)
      }
      
      var absoluteRadius = Element.carbon.covalentRadius
      absoluteRadius += Element.carbon.covalentRadius
      absoluteRadius *= 1.1
      let matchesAbsolute = topology.match(
        topology.atoms, algorithm: .absoluteRadius(absoluteRadius))
      XCTAssertEqual(matchesAbsolute.count, carbons.count)
      XCTAssertEqual(matchesCovalent, matchesAbsolute)
    }
    
    do {
      var topology = Topology()
      topology.insert(atoms: lonsdaleite.atoms)
      let matches = topology.match(topology.atoms)
      XCTAssertEqual(matches.count, lonsdaleite.atoms.count)
      
      for i in matches.indices {
        let atom = lonsdaleite.atoms[i]
        let range = matches[i]
        if atom.atomicNumber == 1 {
          XCTAssertEqual(range.count, 2)
        } else {
          XCTAssertEqual(range.count, 5)
        }
        XCTAssertEqual(range[range.startIndex], UInt32(i))
        
        if atom.atomicNumber == 1 {
          let carbonIndex = Int(range[range.startIndex + 1])
          let carbon = lonsdaleite.atoms[carbonIndex]
          XCTAssertEqual(carbon.atomicNumber, 6)
        }
      }
    }
    
    var expectedAtoms: [Entity]
    var expectedBonds: [SIMD2<UInt32>]
    do {
      var topology = Topology()
      topology.insert(atoms: lonsdaleite.atoms)
      topology.insert(bonds: lonsdaleite.bonds)
      topology.sort()
      XCTAssertGreaterThan(topology.atoms.count, 0)
      XCTAssertGreaterThan(topology.bonds.count, 0)
      
      expectedAtoms = topology.atoms
      expectedBonds = topology.bonds
    }
    
    let hydrogens = lonsdaleite.atoms.filter { $0.atomicNumber == 1 }
    do {
      var topology = Topology()
      topology.insert(atoms: carbons)
      let carbonMatches = topology.match(topology.atoms)
      
      for i in carbons.indices {
        let range = carbonMatches[i]
        XCTAssertGreaterThanOrEqual(range.count, 3)
        XCTAssertLessThanOrEqual(range.count, 5)
        
        var bonds: [SIMD2<UInt32>] = []
        XCTAssertEqual(range[range.startIndex], UInt32(i))
        for j in (range.startIndex + 1)..<range.endIndex {
          XCTAssertNotEqual(range[j], UInt32(i))
          
          // Don't create duplicate bonds.
          let bond = SIMD2<UInt32>(range[j], UInt32(i))
          if any(bond % 3 .== 0) {
            if bond[0] < bond[1] {
              bonds.append(bond)
            }
          } else {
            if bond[1] < bond[0] {
              bonds.append(bond)
            }
          }
        }
        topology.insert(bonds: bonds)
      }
      
      let hydrogenMatches = topology.match(hydrogens)
      XCTAssertEqual(hydrogenMatches.count, hydrogens.count)
      XCTAssertNotEqual(hydrogenMatches.count, carbons.count)
      XCTAssertNotEqual(hydrogenMatches.count, topology.atoms.count)
      XCTAssertNotEqual(hydrogenMatches.count, lonsdaleite.atoms.count)
      
      let hydrogenStart = topology.atoms.count
      topology.insert(atoms: hydrogens)
      XCTAssertGreaterThan(hydrogenMatches.count, 0)
      
      for i in hydrogens.indices {
        let range = hydrogenMatches[i]
        XCTAssertEqual(range.count, 1)
        
        let carbon: UInt32 = range[range.startIndex]
        let hydrogen: UInt32 = UInt32(hydrogenStart + i)
        var bond: SIMD2<UInt32>
        if Bool.random() {
          bond = SIMD2(hydrogen, carbon)
        } else {
          bond = SIMD2(carbon, hydrogen)
        }
        topology.insert(bonds: [bond])
      }
      
      XCTAssertNotEqual(topology.atoms, expectedAtoms)
      XCTAssertNotEqual(topology.bonds, expectedBonds)
      XCTAssertEqual(topology.atoms.count, expectedAtoms.count)
      XCTAssertEqual(topology.bonds.count, expectedBonds.count)
      
      topology.sort()
      XCTAssertEqual(topology.atoms, expectedAtoms)
      XCTAssertEqual(topology.bonds, expectedBonds)
      
      XCTAssertGreaterThan(topology.atoms.count, 0)
      XCTAssertGreaterThan(topology.bonds.count, 0)
      XCTAssertGreaterThan(expectedAtoms.count, 0)
      XCTAssertGreaterThan(expectedBonds.count, 0)
    }
  }
  
  // Test both cubic diamond (only sidewall or malformed carbons) and
  // lonsdaleite (only bridgehead carbons). For the sidewall carbons, assert
  // that the majority of generated hydrogens physically overlap a neighbor.
  func testNonbondingOrbitals() throws {
    func carbonTypeSummary(_ stats: SIMD8<Int>) -> [String] {
      var passivatedCarbonCount = 0
      var totalCarbonCount = 0
      passivatedCarbonCount += stats[1] + stats[2] + stats[3]
      totalCarbonCount += passivatedCarbonCount
      totalCarbonCount += stats[0] + stats[4]
      passivatedCarbonCount = max(1, passivatedCarbonCount)
      
      var lines: [String] = []
      lines.append(
        "type      " + " | " + "count" + " | " + "total" + " | " + "passivated")
      lines.append(
        "----------" + " | " + "-----" + " | " + "-----" + " | " + "-----")
      
      for i in 0...4 {
        var output: String
        switch i {
        case 0: output = "malformed "
        case 1: output = "primary   "
        case 2: output = "sidewall  "
        case 3: output = "bridgehead"
        case 4: output = "bulk      "
        default:
          fatalError("This should never happen.")
        }
        output += " | "
        
        let total = Float(stats[i]) / Float(totalCarbonCount)
        let passivated = Float(stats[i]) / Float(passivatedCarbonCount)
        var countRepr = "\(stats[i])"
        var totalRepr = String(format: "%.1f", 100 * total) + "%"
        var passivatedRepr = String(format: "%.1f", 100 * passivated) + "%"
        while countRepr.count < 5 {
          countRepr = " \(countRepr)"
        }
        while totalRepr.count < 5 {
          totalRepr = " \(totalRepr)"
        }
        while passivatedRepr.count < 5 {
          passivatedRepr = " \(passivatedRepr)"
        }
        
        output += countRepr
        output += " | "
        output += totalRepr
        output += " | "
        if i >= 1 && i <= 3 {
          output += passivatedRepr
        }
        lines.append(output)
      }
      return lines
    }
    
    // Test lonsdaleite and bridgehead carbons. Assert that a small fraction of
    // the carbons are sidewall, but none are malformed.
    do {
      let lonsdaleite = Lonsdaleite()
      var topology = Topology()
      topology.insert(atoms: lonsdaleite.atoms)
      topology.insert(bonds: lonsdaleite.bonds)
      
      let hydrogenIndices = lonsdaleite.atoms.indices.filter {
        lonsdaleite.atoms[$0].atomicNumber == 1
      }
      topology.remove(atoms: hydrogenIndices.map(UInt32.init))
      
      let matches = topology.match(
        topology.atoms, algorithm: .covalentBondLength(1.1))
      let atomsToAtomsMap = topology.map(.atoms, to: .atoms)
      let orbitals = topology.nonbondingOrbitals()
      XCTAssertEqual(matches.count, atomsToAtomsMap.count)
      XCTAssertGreaterThan(matches.count, 0)
      
      var hydrogenTopology = Topology()
      hydrogenTopology.insert(atoms: lonsdaleite.atoms.filter {
        $0.atomicNumber == 1
      })
      let hydrogenMatches = hydrogenTopology.match(topology.atoms)
      XCTAssertGreaterThan(hydrogenTopology.atoms.count, 0)
      XCTAssertEqual(hydrogenMatches.count, topology.atoms.count)
      
      var stats: SIMD8<Int> = .zero
      for i in topology.atoms.indices {
        let matchesRange = matches[i]
        let neighborRange = atomsToAtomsMap[i]
        XCTAssertEqual(matchesRange.count - 1, neighborRange.count)
        
        for j in (matchesRange.startIndex+1)..<matchesRange.endIndex {
          var matchCount: Int = 0
          for k in neighborRange.indices {
            if matchesRange[j] == neighborRange[k] {
              matchCount += 1
            }
          }
          XCTAssertEqual(matchCount, 1)
        }
        
        let orbitalsRange = orbitals[i]
        let hydrogenRange = hydrogenMatches[i]
        XCTAssertEqual(orbitalsRange.count, hydrogenRange.count)
        XCTAssertEqual(orbitalsRange.count, 4 - neighborRange.count)
        XCTAssertEqual(orbitalsRange.indices, hydrogenRange.indices)
        
        var orbitalDirections: [SIMD3<Float>] = []
        var hydrogenDirections: [SIMD3<Float>] = []
        
        for j in orbitalsRange.indices {
          let orbitalDirection = orbitalsRange[j]
          let orbitalLengthSq = (orbitalDirection * orbitalDirection).sum()
          orbitalDirections.append(orbitalDirection)
          XCTAssertEqual(orbitalLengthSq, 1, accuracy: 1e-3)
        }
        for j in hydrogenRange.indices {
          let hydrogenID = hydrogenRange[j]
          let hydrogenAtom = hydrogenTopology.atoms[Int(hydrogenID)]
          let carbonAtom = topology.atoms[Int(i)]
          
          var delta = hydrogenAtom.position - carbonAtom.position
          let deltaLength = (delta * delta).sum().squareRoot()
          delta /= deltaLength
          hydrogenDirections.append(delta)
          XCTAssertGreaterThan(deltaLength, 1e-3)
        }
        
        for orbitalDirection in orbitalDirections {
          var matchCount: Int = 0
          for hydrogenDirection in hydrogenDirections {
            let dotProduct = (orbitalDirection * hydrogenDirection).sum()
            if abs(dotProduct - 1) < 1e-3 {
              matchCount += 1
            }
          }
          XCTAssertEqual(matchCount, 1)
        }
        stats[neighborRange.count] += 1
      }
      
#if false
      let lines = carbonTypeSummary(stats)
      print()
      print("lonsdaleite:")
      for line in lines {
        print(line)
      }
#endif
      XCTAssertEqual(stats[0], 0) // none are malformed
      XCTAssertEqual(stats[1], 0) // none are primary
      XCTAssertEqual(stats[2], 20) // a small fraction are sidewall
      XCTAssertEqual(stats[3], 134) // predominantly bridgehead
      XCTAssertEqual(stats[4], 170)
    }
    
    // Test cubic diamond, malformed, and sidewall carbons. Assert that a small
    // fraction are primary, 4 are malformed, and none are bridgehead.
    do {
      let lattice = Lattice<Cubic> { h, k, l in
        Bounds { 3 * h + 3 * k + 4 * l }
        Material { .elemental(.carbon) }
      }
      var topology = Topology()
      topology.insert(atoms: lattice.atoms)
      let matches = topology.match(topology.atoms)
      XCTAssertEqual(lattice.atoms.count, topology.atoms.count)
      XCTAssertEqual(matches.count, topology.atoms.count)
      
      var statsBefore: SIMD8<Int> = .zero
      for i in topology.atoms.indices {
        let matchesRange = matches[i]
        XCTAssertGreaterThanOrEqual(matchesRange.count, 1)
        XCTAssertLessThanOrEqual(matchesRange.count, 5)
        
        var bonds: [SIMD2<UInt32>] = []
        for j in matchesRange[(matchesRange.startIndex+1)...] {
          if j < UInt32(i) {
            bonds.append(SIMD2(j, UInt32(i)))
          }
        }
        topology.insert(bonds: bonds)
        
        statsBefore[matchesRange.count - 1] += 1
      }
      
      // The topology is not sorted when we create these orbitals. This is
      // different from the conditions for generating lonsdaleite in the
      // unit test 'testMatchLonsdaleite()'.
      let atomsToAtomsMap = topology.map(.atoms, to: .atoms)
      let orbitals = topology.nonbondingOrbitals()
      XCTAssertEqual(atomsToAtomsMap.count, topology.atoms.count)
      XCTAssertEqual(orbitals.count, topology.atoms.count)
      XCTAssertGreaterThan(orbitals.count, 0)
      
      var stats: SIMD8<Int> = .zero
      for i in topology.atoms.indices {
        let matchesRange = matches[i]
        let neighborRange = atomsToAtomsMap[i]
        let orbitalsRange = orbitals[i]
        XCTAssertEqual(matchesRange.count - 1, neighborRange.count)
        
        var expectedOrbitals: Int
        switch neighborRange.count {
        case 0: expectedOrbitals = 0
        case 1: expectedOrbitals = 0
        case 2: expectedOrbitals = 2
        case 3: expectedOrbitals = 1
        case 4: expectedOrbitals = 0
        default:
          fatalError("This should never happen.")
        }
        XCTAssertEqual(orbitalsRange.count, expectedOrbitals)
        
        stats[neighborRange.count] += 1
      }
      
      // Assert that the statistics before and after C-C bond generation are the
      // same.
      for i in 0...4 {
        XCTAssertEqual(statsBefore[i], stats[i])
      }
      
#if false
      let lines = carbonTypeSummary(stats)
      print()
      print("cubic diamond:")
      for line in lines {
        print(line)
      }
#endif
      XCTAssertEqual(stats[0], 4) // 4 are malformed
      XCTAssertEqual(stats[1], 32) // a small fraction are primary
      XCTAssertEqual(stats[2], 98) // predominantly sidewall
      XCTAssertEqual(stats[3], 0) // none are bridgehead
      XCTAssertEqual(stats[4], 231)
    }
  }
  
  // Implementation plan:
  // - 1) Visualizer for Morton order and bond topology in GitHub gist. ✅
  //   - 1.1) Test against Lattice -> Diamondoid reordering. ✅
  //   - 1.2) Test against Topology.sort(). ✅
  // - 2) Test simple diamond and lonsdaleite lattice formation. ✅
  //   - 2.1) Demonstrate (100) surface reconstruction, validate that it's
  //          accepted by MM4RigidBody. ✅
  //   - 2.2) Add reference code to HDLTests, for generating lonsdaleite, cubic
  //          diamond, and curved shell structures with the new compiler. ✅
  //   - 2.3) Add reference code to HDLTests, for generating graphene thiol and
  //          HAbst tripods with the new compiler. ✅
  // - 3) Optimize the compiler.
  //   - 3.1) Add performance test cases for some snippets of real-world code
  //          from these experiments.
  //   - 3.2) Enumerate the time spent on different stages of compilation.
  //   - 3.3) Create an optimized match(), tuned against both extremes of search
  //          radius and asymmetry.
  // - 4) Add performance test for the infamous back board from the nanofactory.
  
  // Elaboration on sections 2, 3, and 4:
  //
  // Creating some demos for adding to HDL unit tests. Due to the volume of the
  // code for (100) surface reconstruction, it cannot be added to the test suite.
  // Therefore, I have to design an alternative strained shell structure built
  // from hexagonal diamond. This alternative may not be manufacturable.
  // - small lonsdaleite shell structure ✅
  // - graphene thiol ✅
  // - HAbst tripod ✅
  //
  // For performance tests, profile 2 different sizes of cubic diamond crystal.
  // The following match operations are plausible during (100) reconstruction.
  // - match carbons against carbons
  // - match a large number of duplicated hydrogens against each other
  //   - uses a small absolute radius
  // - match colliding hydrogens against each other and nearby carbons
  //   - not actually used, but a good test of asymmetry
  //
  // Finally, a performance test of the nanofactory backboard that reports the
  // time spent in each stage. This includes sorting.
}
