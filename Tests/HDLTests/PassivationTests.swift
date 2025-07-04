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
  
  func testCommonLattice() throws {
    let lattice = Self.commonLattice()
    XCTAssertEqual(lattice.atoms.count, 621)
    
    var reconstruction = Reconstruction()
    reconstruction.atoms = lattice.atoms
    reconstruction.material = .checkerboard(.silicon, .carbon)
    let topology = reconstruction.compile()
    
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
    let topology = reconstruction.compile()
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
  
  // Measure the "extended volume" of the structure, whichever exact
  // terminology was used in Nanosystems. Measure according to the center
  // position of the atoms, not the edge of their atomic radius. Fluorine
  // should have a larger volume than hydrogen, which has a larger volume
  // than unpassivated. In all cases, the number of C-Si bonds is the same.
  
  // Test hydrogen passivation during vs. after the reconstruction. If done
  // manually afterward, the "extended volume" should increase. In addition,
  // the hydrogens within a dimer should grow closer. Implement this closeness
  // test by matching hydrogen passivators extracted from the two different
  // crystals. Find the closest 1-2 hydrogens near each abnormally large
  // C-Si bond.
  //
  // Wait...it looks like the new bonding structure already affects the
  // placement of the hydrogens? Calculate the theoretical hydrogen
  // position in the absence of surface reconstruction. Ensure it agrees with
  // actual placed hydrogens in a primitive passivation algorithm that
  // generates a C(100)-(1×1) surface.
  
  // Setting up these tests now, to ensure further code changes don't alter
  // the positions of hydrogens.
  //
  // Will adapt this test afterward, to account for the variability between
  // machines.
  func testHydrogenPositions() throws {
    let lattice = Self.commonLattice()
    
    var reconstruction = Reconstruction()
    reconstruction.atoms = lattice.atoms
    reconstruction.material = .checkerboard(.silicon, .carbon)
    var topology = reconstruction.compile()
    topology.sort()
    
    // Only store data for atoms & bonds that relate to hydrogens
    print(topology.atoms.count)
    print(topology.bonds.count)
  }
#endif
}
