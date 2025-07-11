import XCTest
import HDL

// Unit test how various APIs respond when they receive arrays of zero length.

// APIs:
//
// Topology -> match
//   test finite topology.atoms, zero input
//   test zero topology.atoms, finite input
//   test zero topology.atoms, zero input
// Topology -> nonbondingOrbitals
//   test all of sp, sp2, sp3 procedurally in a loop
//     test finite atoms, no bonds
//     test no atoms, no bonds
// Topology -> remove
//   test all of remove(atoms:), remove(bonds:) procedurally in a loop
//     test finite atoms, finite bonds, no indices
//     test finite atoms, no bonds, no indices
//     test no atoms, no bonds, no indices
// Topology -> sort
//   test no atoms, no bonds
//   check that returned index array has zero length
// Reconstruction
//   compile with zero atoms, and a standard material
//
// Use silicon carbide as the standard material
//
// Run ZeroTests both in Swift release mode and Swift debug mode. Must catch
// runtime validation checks that were triggered (out of bounds, null
// unwrapping, integer overflow). Easiest approach is to always use debug mode
// while developing the tests.

final class ZeroTests: XCTestCase {
  static func createDefaultAtoms() -> [Atom] {
    return [
      Atom(position: SIMD3(0.000, 0.000, 0.000), element: .carbon),
      Atom(position: SIMD3(1.000, 0.000, 0.000), element: .carbon),
      Atom(position: SIMD3(0.000, 1.000, 0.000), element: .carbon),
      Atom(position: SIMD3(0.000, 0.000, 1.000), element: .carbon),
    ]
  }
  
  func testLatticeBounds() throws {
    // Test cubic lattices.
    do {
      let lattice = Lattice<Cubic> { h, k, l in
        Bounds { 1 * h + 1 * k + 1 * l }
        Material { .checkerboard(.silicon, .carbon) }
      }
      XCTAssertEqual(lattice.atoms.count, 18)
    }
    do {
      let lattice = Lattice<Cubic> { h, k, l in
        Bounds { 1 * h + 1 * k + 0 * l }
        Material { .checkerboard(.silicon, .carbon) }
      }
      XCTAssertEqual(lattice.atoms.count, 5)
    }
    do {
      let lattice = Lattice<Cubic> { h, k, l in
        Bounds { 1 * h + 0 * k + 0 * l }
        Material { .checkerboard(.silicon, .carbon) }
      }
      XCTAssertEqual(lattice.atoms.count, 2)
    }
    do {
      let lattice = Lattice<Cubic> { h, k, l in
        Bounds { 0 * h + 0 * k + 0 * l }
        Material { .checkerboard(.silicon, .carbon) }
      }
      XCTAssertEqual(lattice.atoms.count, 1)
    }
    
    // Test hexagonal lattices.
    do {
      let lattice = Lattice<Hexagonal> { h, k, l in
        let h2k = h + 2 * k
        Bounds { 1 * h + 1 * h2k + 1 * l }
        Material { .checkerboard(.silicon, .carbon) }
      }
      XCTAssertEqual(lattice.atoms.count, 12)
    }
    do {
      let lattice = Lattice<Hexagonal> { h, k, l in
        let h2k = h + 2 * k
        Bounds { 0 * h + 1 * h2k + 1 * l }
        Material { .checkerboard(.silicon, .carbon) }
      }
      XCTAssertEqual(lattice.atoms.count, 4)
    }
    do {
      let lattice = Lattice<Hexagonal> { h, k, l in
        let h2k = h + 2 * k
        Bounds { 1 * h + 0 * h2k + 1 * l }
        Material { .checkerboard(.silicon, .carbon) }
      }
      XCTAssertEqual(lattice.atoms.count, 0)
    }
    do {
      let lattice = Lattice<Hexagonal> { h, k, l in
        let h2k = h + 2 * k
        Bounds { 1 * h + 1 * h2k + 0 * l }
        Material { .checkerboard(.silicon, .carbon) }
      }
      XCTAssertEqual(lattice.atoms.count, 0)
    }
    do {
      let lattice = Lattice<Hexagonal> { h, k, l in
        let h2k = h + 2 * k
        Bounds { 1 * h + 0 * h2k + 0 * l }
        Material { .checkerboard(.silicon, .carbon) }
      }
      XCTAssertEqual(lattice.atoms.count, 0)
    }
    do {
      let lattice = Lattice<Hexagonal> { h, k, l in
        let h2k = h + 2 * k
        Bounds { 0 * h + 1 * h2k + 0 * l }
        Material { .checkerboard(.silicon, .carbon) }
      }
      XCTAssertEqual(lattice.atoms.count, 0)
    }
    // crasher:
    //   0 * h + 0 * h2k + 1 * l
    // explicitly forbidden right now:
    //   0 * h + 0 * h2k + 0 * l
  }
  
  func testLatticeBoundsPlane() throws {
    // Test cubic lattices.
    do {
      // 1, 1, 1
      let lattice = Lattice<Cubic> { h, k, l in
        Bounds { 1 * h + 1 * k + 1 * l }
        Material { .checkerboard(.silicon, .carbon) }
        
        Volume {
          Origin { 0.5 * h + 0.5 * k + 0.5 * l }
          Plane { h + k + l }
          Replace { .empty }
        }
      }
      XCTAssertEqual(lattice.atoms.count, 8)
    }
    do {
      // 1, 1, 0
      let lattice = Lattice<Cubic> { h, k, l in
        Bounds { 1 * h + 1 * k + 0 * l }
        Material { .checkerboard(.silicon, .carbon) }
        
        Volume {
          Origin { 0.5 * h + 0.5 * k + 0.5 * l }
          Plane { h + k + l }
          Replace { .empty }
        }
      }
      XCTAssertEqual(lattice.atoms.count, 4)
    }
    do {
      // 1, 0, 0
      let lattice = Lattice<Cubic> { h, k, l in
        Bounds { 1 * h + 0 * k + 0 * l }
        Material { .checkerboard(.silicon, .carbon) }
        
        Volume {
          Origin { 0.5 * h + 0.5 * k + 0.5 * l }
          Plane { h + k + l }
          Replace { .empty }
        }
      }
      XCTAssertEqual(lattice.atoms.count, 2)
    }
    do {
      // 0, 0, 0
      let lattice = Lattice<Cubic> { h, k, l in
        Bounds { 0 * h + 0 * k + 0 * l }
        Material { .checkerboard(.silicon, .carbon) }
        
        Volume {
          Origin { 0.5 * h + 0.5 * k + 0.5 * l }
          Plane { h + k + l }
          Replace { .empty }
        }
      }
      XCTAssertEqual(lattice.atoms.count, 1)
    }
    
    // Test hexagonal lattices.
    do {
      // 1, 1, 1
      let lattice = Lattice<Hexagonal> { h, k, l in
        let h2k = h + 2 * k
        Bounds { 1 * h + 1 * h2k + 1 * l }
        Material { .checkerboard(.silicon, .carbon) }
        
        Volume {
          Origin { 0.5 * h + 0.5 * h2k + 0.5 * l }
          Plane { h + h2k + l }
          Replace { .empty }
        }
      }
      XCTAssertEqual(lattice.atoms.count, 7)
    }
    do {
      // 0, 1, 1
      let lattice = Lattice<Hexagonal> { h, k, l in
        let h2k = h + 2 * k
        Bounds { 0 * h + 1 * h2k + 1 * l }
        Material { .checkerboard(.silicon, .carbon) }
        
        Volume {
          Origin { 0.5 * h + 0.5 * h2k + 0.5 * l }
          Plane { h + h2k + l }
          Replace { .empty }
        }
      }
      XCTAssertEqual(lattice.atoms.count, 3)
    }
    do {
      // 1, 0, 1
      let lattice = Lattice<Hexagonal> { h, k, l in
        let h2k = h + 2 * k
        Bounds { 1 * h + 0 * h2k + 1 * l }
        Material { .checkerboard(.silicon, .carbon) }
        
        Volume {
          Origin { 0.5 * h + 0.5 * h2k + 0.5 * l }
          Plane { h + h2k + l }
          Replace { .empty }
        }
      }
      XCTAssertEqual(lattice.atoms.count, 0)
    }
    do {
      // 1, 1, 0
      let lattice = Lattice<Hexagonal> { h, k, l in
        let h2k = h + 2 * k
        Bounds { 1 * h + 1 * h2k + 0 * l }
        Material { .checkerboard(.silicon, .carbon) }
        
        Volume {
          Origin { 0.5 * h + 0.5 * h2k + 0.5 * l }
          Plane { h + h2k + l }
          Replace { .empty }
        }
      }
      XCTAssertEqual(lattice.atoms.count, 0)
    }
    do {
      // 1, 0, 0
      let lattice = Lattice<Hexagonal> { h, k, l in
        let h2k = h + 2 * k
        Bounds { 1 * h + 0 * h2k + 0 * l }
        Material { .checkerboard(.silicon, .carbon) }
        
        Volume {
          Origin { 0.5 * h + 0.5 * h2k + 0.5 * l }
          Plane { h + h2k + l }
          Replace { .empty }
        }
      }
      XCTAssertEqual(lattice.atoms.count, 0)
    }
    do {
      // 0, 1, 0
      let lattice = Lattice<Hexagonal> { h, k, l in
        let h2k = h + 2 * k
        Bounds { 0 * h + 1 * h2k + 0 * l }
        Material { .checkerboard(.silicon, .carbon) }
        
        Volume {
          Origin { 0.5 * h + 0.5 * h2k + 0.5 * l }
          Plane { h + h2k + l }
          Replace { .empty }
        }
      }
      XCTAssertEqual(lattice.atoms.count, 0)
    }
    // crasher:
    //   0 * h + 0 * h2k + 1 * l
    // explicitly forbidden right now:
    //   0 * h + 0 * h2k + 0 * l
  }
  
  func testTopologyMap() throws {
    // finite atoms, no bonds
    do {
      var topology = Topology()
      topology.atoms = Self.createDefaultAtoms()
      
      // .atoms, .atoms
      do {
        let map = topology.map(.atoms, to: .atoms)
        XCTAssertEqual(map.count, 4)
        for atomID in 0..<4 {
          XCTAssertEqual(map[atomID].count, 0)
        }
      }
      
      // .atoms, .bonds
      do {
        let map = topology.map(.atoms, to: .bonds)
        XCTAssertEqual(map.count, 4)
        for atomID in 0..<4 {
          XCTAssertEqual(map[atomID].count, 0)
        }
      }
      
      // .bonds, .atoms
      do {
        let map = topology.map(.bonds, to: .atoms)
        XCTAssertEqual(map.count, 0)
      }
    }
    
    // no atoms, no bonds
    do {
      var topology = Topology()
      
      // .atoms, .atoms
      do {
        let map = topology.map(.atoms, to: .atoms)
        XCTAssertEqual(map.count, 0)
      }
      
      // .atoms, .bonds
      do {
        let map = topology.map(.atoms, to: .bonds)
        XCTAssertEqual(map.count, 0)
      }
      
      // .bonds, .atoms
      do {
        let map = topology.map(.bonds, to: .atoms)
        XCTAssertEqual(map.count, 0)
      }
    }
  }
}
