import XCTest
import HDL

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
  
  // Good ideas for stuff the test suite should eventually cover:
  //
  // Idea for testing correctness of Topology.match(): Ensure the matched
  // indices actually appear in ascending order of distance. Ensure none of them
  // are more than the algorithm's specified distance cutoff.
  //
  // Idea for ultimate litmus test: use the new, flexible topology compiler to
  // reconstruct (100) surfaces.
  
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
  // The old 'Diamondoid' topology generator can be used to bootstrap testing
  // of how sort() treats bonds, before the remainder of the functionality is
  // working. In addition, enter into MM4Parameters as a final validation test.
  //
  // Implementation plan:
  // - 1) Visualizer for Morton order and bond topology in GitHub gist.
  //   - 1.1) Test against Lattice -> Diamondoid reordering.
  //   - 1.2) Test against Topology.sort().
  // - 2) Debug correctness of match() and nonbondedOrbitals().
  // - 3) Test simple diamond and lonsdaleite lattice formation.
  // - 4) Reproduce graphene thiol and HAbst tripods from the
  //      nanofactory demo.
  // - 5) Demonstrate (100) surface and strained shell structure formation.
}
