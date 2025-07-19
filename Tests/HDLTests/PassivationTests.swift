import XCTest
import HDL

final class PassivationTests: XCTestCase {
  static func commonLattice() -> Lattice<Cubic> {
    Lattice<Cubic> { h, k, l in
      Bounds { 4 * h + 4 * k + 4 * l }
      Material { .checkerboard(.silicon, .carbon) }
    }
  }
  
  static func checkConnectivity(_ topology: Topology) {
    let atomsToBondsMap = topology.map(.atoms, to: .bonds)
    for atomID in topology.atoms.indices {
      let atom = topology.atoms[atomID]
      let bondsMap = atomsToBondsMap[atomID]
      if atom.atomicNumber == 1 {
        XCTAssertEqual(bondsMap.count, 1)
      } else {
        XCTAssertEqual(bondsMap.count, 4)
      }
    }
  }
  
  static func checkNoOverlaps(_ topology: Topology) {
    let matchRanges = topology.match(
      topology.atoms, algorithm: .absoluteRadius(0.010))
    for atomID in topology.atoms.indices {
      let matchRange = matchRanges[atomID]
      XCTAssertEqual(matchRange.count, 1)
    }
  }
  
  static func passivate(topology: inout Topology) {
    func createHydrogen(
      atomID: UInt32,
      orbital: SIMD3<Float>
    ) -> Atom {
      let atom = topology.atoms[Int(atomID)]
      
      var bondLength = atom.element.covalentRadius
      bondLength += Element.hydrogen.covalentRadius
      
      let position = atom.position + bondLength * orbital
      return Atom(position: position, element: .hydrogen)
    }
    
    let orbitalLists = topology.nonbondingOrbitals()
    
    var insertedAtoms: [Atom] = []
    var insertedBonds: [SIMD2<UInt32>] = []
    for atomID in topology.atoms.indices {
      let orbitalList = orbitalLists[atomID]
      for orbital in orbitalList {
        let hydrogen = createHydrogen(
          atomID: UInt32(atomID),
          orbital: orbital)
        let hydrogenID = topology.atoms.count + insertedAtoms.count
        insertedAtoms.append(hydrogen)
        
        let bond = SIMD2(
          UInt32(atomID),
          UInt32(hydrogenID))
        insertedBonds.append(bond)
      }
    }
    topology.atoms += insertedAtoms
    topology.bonds += insertedBonds
  }
  
  func testCommonLattice() throws {
    let lattice = Self.commonLattice()
    XCTAssertEqual(lattice.atoms.count, 621)
    
    var reconstruction = Reconstruction()
    reconstruction.atoms = lattice.atoms
    reconstruction.material = .checkerboard(.silicon, .carbon)
    var topology = reconstruction.compile()
    PassivationTests.passivate(topology: &topology)
    
    var groupIVAtomCount: Int = .zero
    for atom in topology.atoms {
      if atom.atomicNumber != 1 {
        groupIVAtomCount += 1
      }
    }
    XCTAssertEqual(groupIVAtomCount, 577)
    
    let hydrogenAtomCount = topology.atoms.count - groupIVAtomCount
    XCTAssertEqual(hydrogenAtomCount, 232)
    XCTAssertEqual(topology.bonds.count, 1270)
  }
  
#if RELEASE
  func testNoOverlaps() throws {
    let lattice = Self.commonLattice()
    
    var reconstruction = Reconstruction()
    reconstruction.atoms = lattice.atoms
    reconstruction.material = .checkerboard(.silicon, .carbon)
    var topology = reconstruction.compile()
    PassivationTests.passivate(topology: &topology)
    PassivationTests.checkConnectivity(topology)
    
    var hydrogenAtoms: [Atom] = []
    for atom in topology.atoms {
      if atom.atomicNumber == 1 {
        hydrogenAtoms.append(atom)
      }
    }
    
    var hydrogenTopology = Topology()
    hydrogenTopology.atoms = hydrogenAtoms
    let matchRanges = hydrogenTopology.match(
      hydrogenAtoms, algorithm: .absoluteRadius(0.010))
    
    var matchCountStats: SIMD8<Int> = .zero
    for hydrogenID in hydrogenAtoms.indices {
      let matchRange = matchRanges[hydrogenID]
      matchCountStats[matchRange.count] += 1
    }
    XCTAssertEqual(hydrogenAtoms.count, 232)
    XCTAssertEqual(matchCountStats[0], 0)
    XCTAssertEqual(matchCountStats[1], 232)
    XCTAssertEqual(matchCountStats[2], 0)
  }
  
  func testHydrogenPositions() throws {
    func createHydrogens() -> [SIMD4<Float>] {
      let lattice = Self.commonLattice()
      
      var reconstruction = Reconstruction()
      reconstruction.atoms = lattice.atoms
      reconstruction.material = .checkerboard(.silicon, .carbon)
      var oldTopology = reconstruction.compile()
      PassivationTests.passivate(topology: &oldTopology)
      
      var topology = Topology()
      topology.atoms = oldTopology.atoms.filter { $0[3] == 1 }
      topology.sort() // TODO: remove this
      return topology.atoms
    }
    
    let expectedHydrogens = Self.expectedSortedHydrogens()
    let hydrogens = createHydrogens()
    XCTAssertEqual(expectedHydrogens.count, hydrogens.count)
    
    for hydrogenID in 0..<232 {
      let expectedHydrogen = expectedHydrogens[hydrogenID]
      let hydrogen = hydrogens[hydrogenID]
      let delta = hydrogen.position - expectedHydrogen.position
      XCTAssertLessThan(delta[0].magnitude, 0.001)
      XCTAssertLessThan(delta[1].magnitude, 0.001)
      XCTAssertLessThan(delta[2].magnitude, 0.001)
    }
  }
#endif
  
  // A sorted list of hydrogens from one invocation of 'Reconstruction' on the
  // common lattice. The exact results may vary from machine to machine.
  //
  // 'Reconstruction' has since been revised to remove the built-in hydrogen
  // passivation feature, thus simplifying and generalizing the API. This
  // unit test was critical to correctly implementing the elision.
  //
  // TODO: Revise this test to account for the possible change in order from a
  // new sorting algorithm. Passivation and Sort tests should cover distinct,
  // decoupled pieces of code.
  private static func expectedSortedHydrogens() -> [SIMD4<Float>] {
    return [
      SIMD4<Float>(0.047, 0.047, 0.047, 1),
      SIMD4<Float>(0.483, 0.047, 0.047, 1),
      SIMD4<Float>(0.300, 0.136, -0.082, 1),
      SIMD4<Float>(0.300, -0.082, 0.136, 1),
      SIMD4<Float>(0.136, 0.300, -0.082, 1),
      SIMD4<Float>(0.047, 0.483, 0.047, 1),
      SIMD4<Float>(-0.082, 0.300, 0.136, 1),
      SIMD4<Float>(0.370, 0.502, -0.107, 1),
      SIMD4<Float>(0.136, -0.082, 0.300, 1),
      SIMD4<Float>(-0.082, 0.136, 0.300, 1),
      SIMD4<Float>(0.047, 0.047, 0.483, 1),
      SIMD4<Float>(0.370, -0.107, 0.502, 1),
      SIMD4<Float>(-0.107, 0.370, 0.502, 1),
      SIMD4<Float>(0.720, 0.152, -0.107, 1),
      SIMD4<Float>(0.720, -0.107, 0.152, 1),
      SIMD4<Float>(0.919, 0.047, 0.047, 1),
      SIMD4<Float>(1.156, 0.152, -0.107, 1),
      SIMD4<Float>(1.156, -0.107, 0.152, 1),
      SIMD4<Float>(0.806, 0.502, -0.107, 1),
      SIMD4<Float>(0.806, -0.107, 0.502, 1),
      SIMD4<Float>(0.136, 0.736, -0.082, 1),
      SIMD4<Float>(-0.107, 0.720, 0.152, 1),
      SIMD4<Float>(0.300, 0.572, -0.082, 1),
      SIMD4<Float>(0.047, 0.919, 0.047, 1),
      SIMD4<Float>(0.136, 1.172, -0.082, 1),
      SIMD4<Float>(-0.107, 1.156, 0.152, 1),
      SIMD4<Float>(0.370, 0.938, -0.107, 1),
      SIMD4<Float>(0.300, 1.008, -0.082, 1),
      SIMD4<Float>(-0.107, 0.806, 0.502, 1),
      SIMD4<Float>(0.720, 0.588, -0.107, 1),
      SIMD4<Float>(1.156, 0.588, -0.107, 1),
      SIMD4<Float>(0.806, 0.938, -0.107, 1),
      SIMD4<Float>(0.720, 1.024, -0.107, 1),
      SIMD4<Float>(1.156, 1.024, -0.107, 1),
      SIMD4<Float>(0.136, -0.082, 0.736, 1),
      SIMD4<Float>(-0.082, 0.136, 0.736, 1),
      SIMD4<Float>(0.300, -0.082, 0.572, 1),
      SIMD4<Float>(-0.082, 0.300, 0.572, 1),
      SIMD4<Float>(0.047, 0.047, 0.919, 1),
      SIMD4<Float>(0.136, -0.082, 1.172, 1),
      SIMD4<Float>(-0.082, 0.136, 1.172, 1),
      SIMD4<Float>(0.370, -0.107, 0.938, 1),
      SIMD4<Float>(0.300, -0.082, 1.008, 1),
      SIMD4<Float>(-0.107, 0.370, 0.938, 1),
      SIMD4<Float>(-0.082, 0.300, 1.008, 1),
      SIMD4<Float>(0.720, -0.107, 0.588, 1),
      SIMD4<Float>(1.156, -0.107, 0.588, 1),
      SIMD4<Float>(0.806, -0.107, 0.938, 1),
      SIMD4<Float>(0.720, -0.107, 1.024, 1),
      SIMD4<Float>(1.156, -0.107, 1.024, 1),
      SIMD4<Float>(-0.107, 0.720, 0.588, 1),
      SIMD4<Float>(-0.107, 1.156, 0.588, 1),
      SIMD4<Float>(-0.107, 0.806, 0.938, 1),
      SIMD4<Float>(-0.107, 0.720, 1.024, 1),
      SIMD4<Float>(-0.107, 1.156, 1.024, 1),
      SIMD4<Float>(1.355, 0.047, 0.047, 1),
      SIMD4<Float>(1.592, 0.152, -0.107, 1),
      SIMD4<Float>(1.592, -0.107, 0.152, 1),
      SIMD4<Float>(1.242, 0.502, -0.107, 1),
      SIMD4<Float>(1.697, 0.389, 0.047, 1),
      SIMD4<Float>(1.242, -0.107, 0.502, 1),
      SIMD4<Float>(1.697, 0.047, 0.389, 1),
      SIMD4<Float>(1.851, 0.152, 0.152, 1),
      SIMD4<Float>(1.851, 0.502, 0.502, 1),
      SIMD4<Float>(1.592, 0.588, -0.107, 1),
      SIMD4<Float>(1.697, 0.825, 0.047, 1),
      SIMD4<Float>(1.242, 0.938, -0.107, 1),
      SIMD4<Float>(1.592, 1.024, -0.107, 1),
      SIMD4<Float>(1.851, 0.588, 0.152, 1),
      SIMD4<Float>(1.851, 1.024, 0.152, 1),
      SIMD4<Float>(1.851, 0.938, 0.502, 1),
      SIMD4<Float>(1.592, -0.107, 0.588, 1),
      SIMD4<Float>(1.697, 0.047, 0.825, 1),
      SIMD4<Float>(1.242, -0.107, 0.938, 1),
      SIMD4<Float>(1.592, -0.107, 1.024, 1),
      SIMD4<Float>(1.851, 0.152, 0.588, 1),
      SIMD4<Float>(1.851, 0.152, 1.024, 1),
      SIMD4<Float>(1.851, 0.502, 0.938, 1),
      SIMD4<Float>(1.851, 0.588, 0.588, 1),
      SIMD4<Float>(1.851, 1.024, 0.588, 1),
      SIMD4<Float>(1.851, 0.588, 1.024, 1),
      SIMD4<Float>(1.851, 0.938, 0.938, 1),
      SIMD4<Float>(1.851, 1.024, 1.024, 1),
      SIMD4<Float>(0.047, 1.355, 0.047, 1),
      SIMD4<Float>(0.370, 1.374, -0.107, 1),
      SIMD4<Float>(0.300, 1.444, -0.082, 1),
      SIMD4<Float>(0.136, 1.608, -0.082, 1),
      SIMD4<Float>(-0.107, 1.592, 0.152, 1),
      SIMD4<Float>(0.389, 1.697, 0.047, 1),
      SIMD4<Float>(-0.107, 1.242, 0.502, 1),
      SIMD4<Float>(0.047, 1.697, 0.389, 1),
      SIMD4<Float>(0.806, 1.374, -0.107, 1),
      SIMD4<Float>(0.736, 1.444, -0.082, 1),
      SIMD4<Float>(1.172, 1.444, -0.082, 1),
      SIMD4<Float>(0.572, 1.608, -0.082, 1),
      SIMD4<Float>(0.825, 1.697, 0.047, 1),
      SIMD4<Float>(1.008, 1.608, -0.082, 1),
      SIMD4<Float>(0.152, 1.851, 0.152, 1),
      SIMD4<Float>(0.502, 1.851, 0.502, 1),
      SIMD4<Float>(0.588, 1.851, 0.152, 1),
      SIMD4<Float>(1.024, 1.851, 0.152, 1),
      SIMD4<Float>(0.938, 1.851, 0.502, 1),
      SIMD4<Float>(-0.107, 1.592, 0.588, 1),
      SIMD4<Float>(0.047, 1.697, 0.825, 1),
      SIMD4<Float>(-0.107, 1.242, 0.938, 1),
      SIMD4<Float>(-0.107, 1.592, 1.024, 1),
      SIMD4<Float>(0.152, 1.851, 0.588, 1),
      SIMD4<Float>(0.152, 1.851, 1.024, 1),
      SIMD4<Float>(0.502, 1.851, 0.938, 1),
      SIMD4<Float>(0.588, 1.851, 0.588, 1),
      SIMD4<Float>(1.024, 1.851, 0.588, 1),
      SIMD4<Float>(0.588, 1.851, 1.024, 1),
      SIMD4<Float>(0.938, 1.851, 0.938, 1),
      SIMD4<Float>(1.024, 1.851, 1.024, 1),
      SIMD4<Float>(1.242, 1.374, -0.107, 1),
      SIMD4<Float>(1.697, 1.261, 0.047, 1),
      SIMD4<Float>(1.608, 1.444, -0.082, 1),
      SIMD4<Float>(1.826, 1.444, 0.136, 1),
      SIMD4<Float>(1.444, 1.608, -0.082, 1),
      SIMD4<Float>(1.261, 1.697, 0.047, 1),
      SIMD4<Float>(1.444, 1.826, 0.136, 1),
      SIMD4<Float>(1.697, 1.697, 0.047, 1),
      SIMD4<Float>(1.826, 1.608, 0.300, 1),
      SIMD4<Float>(1.608, 1.826, 0.300, 1),
      SIMD4<Float>(1.697, 1.697, 0.483, 1),
      SIMD4<Float>(1.851, 1.374, 0.502, 1),
      SIMD4<Float>(1.374, 1.851, 0.502, 1),
      SIMD4<Float>(1.826, 1.444, 0.572, 1),
      SIMD4<Float>(1.444, 1.826, 0.572, 1),
      SIMD4<Float>(1.826, 1.608, 0.736, 1),
      SIMD4<Float>(1.608, 1.826, 0.736, 1),
      SIMD4<Float>(1.826, 1.444, 1.008, 1),
      SIMD4<Float>(1.444, 1.826, 1.008, 1),
      SIMD4<Float>(1.697, 1.697, 0.919, 1),
      SIMD4<Float>(1.826, 1.608, 1.172, 1),
      SIMD4<Float>(1.608, 1.826, 1.172, 1),
      SIMD4<Float>(1.851, 1.374, 0.938, 1),
      SIMD4<Float>(1.374, 1.851, 0.938, 1),
      SIMD4<Float>(0.047, 0.047, 1.355, 1),
      SIMD4<Float>(0.370, -0.107, 1.374, 1),
      SIMD4<Float>(0.300, -0.082, 1.444, 1),
      SIMD4<Float>(-0.107, 0.370, 1.374, 1),
      SIMD4<Float>(-0.082, 0.300, 1.444, 1),
      SIMD4<Float>(0.136, -0.082, 1.608, 1),
      SIMD4<Float>(-0.082, 0.136, 1.608, 1),
      SIMD4<Float>(0.389, 0.047, 1.697, 1),
      SIMD4<Float>(0.047, 0.389, 1.697, 1),
      SIMD4<Float>(0.806, -0.107, 1.374, 1),
      SIMD4<Float>(0.736, -0.082, 1.444, 1),
      SIMD4<Float>(1.172, -0.082, 1.444, 1),
      SIMD4<Float>(0.572, -0.082, 1.608, 1),
      SIMD4<Float>(0.825, 0.047, 1.697, 1),
      SIMD4<Float>(1.008, -0.082, 1.608, 1),
      SIMD4<Float>(-0.107, 0.806, 1.374, 1),
      SIMD4<Float>(-0.082, 0.736, 1.444, 1),
      SIMD4<Float>(-0.082, 1.172, 1.444, 1),
      SIMD4<Float>(-0.082, 0.572, 1.608, 1),
      SIMD4<Float>(0.047, 0.825, 1.697, 1),
      SIMD4<Float>(-0.082, 1.008, 1.608, 1),
      SIMD4<Float>(0.152, 0.152, 1.851, 1),
      SIMD4<Float>(0.502, 0.502, 1.851, 1),
      SIMD4<Float>(0.588, 0.152, 1.851, 1),
      SIMD4<Float>(1.024, 0.152, 1.851, 1),
      SIMD4<Float>(0.938, 0.502, 1.851, 1),
      SIMD4<Float>(0.152, 0.588, 1.851, 1),
      SIMD4<Float>(0.152, 1.024, 1.851, 1),
      SIMD4<Float>(0.502, 0.938, 1.851, 1),
      SIMD4<Float>(0.588, 0.588, 1.851, 1),
      SIMD4<Float>(1.024, 0.588, 1.851, 1),
      SIMD4<Float>(0.588, 1.024, 1.851, 1),
      SIMD4<Float>(0.938, 0.938, 1.851, 1),
      SIMD4<Float>(1.024, 1.024, 1.851, 1),
      SIMD4<Float>(1.242, -0.107, 1.374, 1),
      SIMD4<Float>(1.697, 0.047, 1.261, 1),
      SIMD4<Float>(1.608, -0.082, 1.444, 1),
      SIMD4<Float>(1.826, 0.136, 1.444, 1),
      SIMD4<Float>(1.444, -0.082, 1.608, 1),
      SIMD4<Float>(1.261, 0.047, 1.697, 1),
      SIMD4<Float>(1.444, 0.136, 1.826, 1),
      SIMD4<Float>(1.697, 0.047, 1.697, 1),
      SIMD4<Float>(1.826, 0.300, 1.608, 1),
      SIMD4<Float>(1.608, 0.300, 1.826, 1),
      SIMD4<Float>(1.697, 0.483, 1.697, 1),
      SIMD4<Float>(1.851, 0.502, 1.374, 1),
      SIMD4<Float>(1.826, 0.572, 1.444, 1),
      SIMD4<Float>(1.826, 1.008, 1.444, 1),
      SIMD4<Float>(1.444, 0.572, 1.826, 1),
      SIMD4<Float>(1.826, 0.736, 1.608, 1),
      SIMD4<Float>(1.608, 0.736, 1.826, 1),
      SIMD4<Float>(1.444, 1.008, 1.826, 1),
      SIMD4<Float>(1.826, 1.172, 1.608, 1),
      SIMD4<Float>(1.697, 0.919, 1.697, 1),
      SIMD4<Float>(1.608, 1.172, 1.826, 1),
      SIMD4<Float>(1.851, 0.938, 1.374, 1),
      SIMD4<Float>(1.374, 0.502, 1.851, 1),
      SIMD4<Float>(1.374, 0.938, 1.851, 1),
      SIMD4<Float>(-0.107, 1.242, 1.374, 1),
      SIMD4<Float>(0.047, 1.697, 1.261, 1),
      SIMD4<Float>(-0.082, 1.608, 1.444, 1),
      SIMD4<Float>(0.136, 1.826, 1.444, 1),
      SIMD4<Float>(-0.082, 1.444, 1.608, 1),
      SIMD4<Float>(0.047, 1.261, 1.697, 1),
      SIMD4<Float>(0.136, 1.444, 1.826, 1),
      SIMD4<Float>(0.047, 1.697, 1.697, 1),
      SIMD4<Float>(0.300, 1.826, 1.608, 1),
      SIMD4<Float>(0.300, 1.608, 1.826, 1),
      SIMD4<Float>(0.483, 1.697, 1.697, 1),
      SIMD4<Float>(0.572, 1.826, 1.444, 1),
      SIMD4<Float>(1.008, 1.826, 1.444, 1),
      SIMD4<Float>(0.572, 1.444, 1.826, 1),
      SIMD4<Float>(1.008, 1.444, 1.826, 1),
      SIMD4<Float>(0.736, 1.826, 1.608, 1),
      SIMD4<Float>(0.736, 1.608, 1.826, 1),
      SIMD4<Float>(1.172, 1.826, 1.608, 1),
      SIMD4<Float>(1.172, 1.608, 1.826, 1),
      SIMD4<Float>(0.919, 1.697, 1.697, 1),
      SIMD4<Float>(0.502, 1.851, 1.374, 1),
      SIMD4<Float>(0.938, 1.851, 1.374, 1),
      SIMD4<Float>(0.502, 1.374, 1.851, 1),
      SIMD4<Float>(0.938, 1.374, 1.851, 1),
      SIMD4<Float>(1.826, 1.444, 1.444, 1),
      SIMD4<Float>(1.444, 1.826, 1.444, 1),
      SIMD4<Float>(1.697, 1.697, 1.355, 1),
      SIMD4<Float>(1.444, 1.444, 1.826, 1),
      SIMD4<Float>(1.697, 1.355, 1.697, 1),
      SIMD4<Float>(1.355, 1.697, 1.697, 1),
      SIMD4<Float>(1.826, 1.608, 1.608, 1),
      SIMD4<Float>(1.608, 1.826, 1.608, 1),
      SIMD4<Float>(1.608, 1.608, 1.826, 1),
      SIMD4<Float>(1.851, 1.374, 1.374, 1),
      SIMD4<Float>(1.374, 1.851, 1.374, 1),
      SIMD4<Float>(1.374, 1.374, 1.851, 1),
    ]
  }
}
