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
    
    // Specify the screen parameters.
    let binSizes = SIMD2<Float>(0.05, 0.1)
    let binCounts = SIMD2<Int>(20, 10)
    let origin = SIMD2<Float>(0.0, 0.0)
    var pixelArray = [Bool](
      repeating: false,
      count: binCounts.x * binCounts.y)
    
    // Mark the atoms on the screen.
    for atom in atoms {
      let coordsXY = SIMD2(atom.x, atom.y) - origin
      let bin = (coordsXY / binSizes).rounded(.down)
      let binInt = SIMD2<Int>(bin)
      if any(binInt .< 0) || any(binInt .>= binCounts) {
        fatalError("Atom \(atom) was outside of bounds: \(bin)")
      }
      
      let binMemoryIndex = binInt.y * binCounts.x + binInt.x
      pixelArray[binMemoryIndex] = true
    }
    
    // Display the screen.
    for pixelPositionY in (0..<binCounts.y).reversed() {
      var line: String = ""
      for pixelPositionX in 0..<binCounts.x {
        let pixelMemoryIndex = pixelPositionY * binCounts.x + pixelPositionX
        let pixelValue = pixelArray[pixelMemoryIndex]
        if pixelValue {
          line.append("X")
        } else {
          line.append(" ")
        }
      }
      print(line)
    }
  }
}
