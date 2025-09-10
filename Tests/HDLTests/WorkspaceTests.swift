import HDL
import XCTest

// Temporary file for treating the Swift tests like a scratch workspace.

final class WorkspaceTests: XCTestCase {
  func testWorkspace() throws {
    let atoms: [SIMD4<Float>] = [
      SIMD4(0.00, 0.00, 0.00, 6),
      SIMD4(0.00, 0.00, 0.00, 6),
      SIMD4(0.00, 0.50, 0.00, 6),
      SIMD4(0.00, 0.00, 0.50, 6),
    ]
    
    // Specify the screen parameters.
    let binSize: Float = 0.1
    let binCount: Int = 10
    let origin = SIMD2<Float>(0.0, 0.0)
    var pixelArray = [Bool](
      repeating: true,
      count: binCount * binCount)
    
    // Mark the atoms on the screen.
    for atom in atoms {
      let coordsXY = SIMD2(atom.x, atom.y) - origin
      let bin = (coordsXY / binSize).rounded(.down)
      let binInt = SIMD2<Int>(bin)
      if any(binInt .< 0) || any(binInt .>= binCount) {
        fatalError("Atom \(atom) was outside of bounds: \(bin)")
      }
      
      let binMemoryIndex = binInt.y * binCount + binInt.x
      pixelArray[binMemoryIndex] = true
    }
    
    // Display the screen.
    
  }
}
