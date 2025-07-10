import XCTest
import HDL

// Unit test how various APIs respond when they receive arrays of zero length.
//
// APIs:
//
// Lattice -> Bounds -> atoms
//   test zero bounds
// Lattice -> Bounds -> Volume/Plane -> atoms
//   use a plane that changes the output in the finite case
//   test both finite and zero bounds
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
// Run ZeroTests both in Swift release mode and Swift debug mode. Must catch
// runtime validation checks that were triggered (out of bounds, null
// unwrapping, integer overflow). Easiest approach is to always use debug mode
// while developing the tests.

final class ZeroTests: XCTestCase {
  
}
