import HDL
import XCTest

// Temporary file for treating the Swift tests like a scratch workspace.

final class WorkspaceTests: XCTestCase {
  func testWorkspace() throws {
    let atoms: [SIMD4<Float>] = [
      SIMD4(0.00, 0.00, 0.00, 6),
      SIMD4(0.50, 0.00, 0.00, 6),
      SIMD4(0.00, 0.50, 0.00, 6),
      SIMD4(0.00, 0.00, 0.50, 6),
    ]
    
    // Objective: render these atoms in a 2D command-line render.
  }
}
