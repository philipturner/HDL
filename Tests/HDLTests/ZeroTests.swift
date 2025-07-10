import XCTest
import HDL

// Unit test how various APIs respond when they receive arrays of zero length.

// APIs:
//
// test all of Cubic, Hexagonal procedurally in a loop
//   Lattice -> Bounds -> atoms
//     test (0, 0, finite) for each of 3 dimensions
//     test zero bounds
//   Lattice -> Bounds -> Volume/Plane -> atoms
//     use a plane that changes the output in the finite case
//     test finite bounds
//     test (0, 0, finite) for each of 3 dimensions
//     test zero bounds
// Topology -> map
//   test finite atoms, no bonds
//   test no atoms, no bonds
// Topology -> match
//   test finite topology.atoms, zero input
//   test zero topology.atoms, finite input
//   test zero topology.atoms, zero input
// Topology -> orbitals
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
  }
  
  func testLatticeBoundsPlane() throws {
    
  }
}
