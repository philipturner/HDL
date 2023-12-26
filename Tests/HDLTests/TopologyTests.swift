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
      lines.append("type      " + " | " + "total" + " | " + "passivated")
      lines.append("----------" + " | " + "-----" + " | " + "-----")
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
        var totalRepr = String(format: "%.1f", 100 * total) + "%"
        var passivatedRepr = String(format: "%.1f", 100 * passivated) + "%"
        while totalRepr.count < 5 {
          totalRepr = " \(totalRepr)"
        }
        while passivatedRepr.count < 5 {
          passivatedRepr = " \(passivatedRepr)"
        }
        
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
        let matchRange = matches[i]
        let neighborRange = atomsToAtomsMap[i]
        XCTAssertEqual(matchRange.count - 1, neighborRange.count)
        
        for j in (matchRange.startIndex+1)..<matchRange.endIndex {
          var matchCount: Int = 0
          for k in neighborRange.indices {
            if matchRange[j] == neighborRange[k] {
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
        
        var orbitalID = 0
        for j in orbitalsRange.indices {
          defer { orbitalID += 1 }
          let orbitalDirection = orbitalsRange[j]
          let hydrogenID = hydrogenRange[j]
          let hydrogenAtom = hydrogenTopology.atoms[Int(hydrogenID)]
          let carbonAtom = topology.atoms[Int(i)]
          
          var delta = hydrogenAtom.position - carbonAtom.position
          let orbitalLengthSq = (orbitalDirection * orbitalDirection).sum()
          let deltaLength = (delta * delta).sum().squareRoot()
          delta /= deltaLength
          XCTAssertEqual(orbitalLengthSq, 1, accuracy: 1e-3)
          XCTAssertGreaterThan(deltaLength, 1e-3)
          
          if orbitalsRange.count == 2 &&
              (orbitalDirection * delta).sum() < 0.99999 &&
             (orbitalsRange[j + 1 - 2 * orbitalID] * delta).sum() < 0.99999 {
            print(i, neighborRange.count, orbitalDirection, delta)
            if j == orbitalsRange.startIndex {
              continue
            }
            
            var deltas: [SIMD3<Float>] = []
            for neighborID in neighborRange {
              let neighbor = topology.atoms[Int(neighborID)]
              var delta = neighbor.position - carbonAtom.position
              delta /= (delta * delta).sum().squareRoot()
              deltas.append(delta)
              print("-", delta)
            }
            var normal = deltas[0] + deltas[1]
            var axis = deltas[1] - deltas[0]
            normal /= (normal * normal).sum().squareRoot()
            axis /= (axis * axis).sum().squareRoot()
            
            var midPoint = (deltas[0] + deltas[1]) / 2
            var axis2 = deltas[1] - midPoint
            midPoint /= (midPoint * midPoint).sum().squareRoot()
            axis2 /= (axis2 * axis2).sum().squareRoot()
            print("--", normal, midPoint)
            print("--", axis, axis2)
            
            let quaternion = Quaternion<Float>(
              angle: 109.5/2 * Float.pi / 180,
              axis: axis)
            let rotated = quaternion.act(on: -normal)
            let quaternion2 = Quaternion<Float>(
              angle: -109.5/2 * Float.pi / 180,
              axis: axis)
            let rotated2 = quaternion2.act(on: -normal)
            print("---", rotated)
            print("---", rotated2)
            
            var crossProduct: SIMD3<Float> = .zero
            crossProduct.x = axis.y * normal.z - axis.z * normal.y
            crossProduct.y = axis.z * normal.x - axis.x * normal.z
            crossProduct.z = axis.x * normal.y - axis.y * normal.x
            let crossProductSquared = (crossProduct * crossProduct).sum()
            crossProduct /= crossProductSquared.squareRoot()
            
            var altCrossProduct = rotated2 - rotated
            let altCrossProductSq = (altCrossProduct * altCrossProduct).sum()
            altCrossProduct /= altCrossProductSq.squareRoot()
            
            print("---", crossProduct)
            print("---", altCrossProduct)
            print("---", (rotated * crossProduct).sum() * (rotated * crossProduct).sum())
            print("---", (rotated * normal).sum() * (rotated * normal).sum())
            
            /*
             let normal = _cross_platform_normalize(thisAtom.origin - midPoint)
             let axis = _cross_platform_normalize(neighborCenters[1] - midPoint)
             for angle in [-sp3BondAngle / 2, sp3BondAngle / 2] {
               let rotation = Quaternion<Float>(angle: angle, axis: axis)
               let direction = rotation.act(on: normal)
               addHydrogen(direction: direction)
             }
             */
          }
        }
        stats[neighborRange.count] += 1
      }
#if true
      let lines = carbonTypeSummary(stats)
      print()
      print("lonsdaleite:")
      for line in lines {
        print(line)
      }
#endif
    }
    
    // Test cubic diamond, malformed, and sidewall carbons. Assert that a small
    // fraction are primary, 4 are malformed, and none are bridgehead.
  }
  
  // TODO: - Instead of a time-consuming, exhaustive test suite, debug this
  // visually in the renderer. The covered functionality is not as complex as
  // MM4RigidBody and doesn't need the same approach to testing. Getting (100)
  // reconstruction to work will likely trigger edge cases where the compiler is
  // broken; reproducers can be added to the test suite.
  //
  // TODO: - The second test case is an interesting method of forming strained
  // shell structures. Passivate a crystalline lattice, then warp and remove
  // hydrogens bonded to carbons that will merge. Validate that both this
  // and the reconstructed (100) structure are accepted by the simulator.
  //
  // TODO: - Tutorials demonstrating correct formation of the following. The
  // audience is my future self.
  // - Graphene thiol
  // - HAbst tripod w/ energy minimization in xTB
  //
  // The old 'Diamondoid' topology generator can be used to bootstrap testing
  // of how sort() treats bonds, before the remainder of the functionality is
  // working. In addition, enter into MM4Parameters as a final validation test.
  // Finally, run the structures through the old MM4 simulator for a few ps.
  //
  // Implementation plan:
  // - 1) Visualizer for Morton order and bond topology in GitHub gist. ✅
  //   - 1.1) Test against Lattice -> Diamondoid reordering. ✅
  //   - 1.2) Test against Topology.sort(). ✅
  // - 2) Test simple diamond and lonsdaleite lattice formation. ✅
  //   - 2.1) Visually debug hydrogen generation against Diamondoid.
  // - 3) Reproduce graphene thiol and HAbst tripods from the
  //      nanofactory demo.
  //   - 3.1) Create tutorials for these two molecules.
  // - 4) Debug performance of match() and nonbondedOrbitals().
  //   - 4.1) Use nonbondedOrbitals() to create large performance tests for
  //          asymmetric match().
  //   - 4.2) Ensure there's no unusually slow code in nonbondedOrbitals().
  //   - 4.3) Create an optimized match(), tuned against two radii/algorithms
  //          for symmetric and two variants of plausible asymmetric search.
  // - 5) Demonstrate (100) surface and strained shell structure formation.
}
