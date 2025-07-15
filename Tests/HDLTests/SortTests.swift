import XCTest
@testable import HDL

final class SortTests: XCTestCase {
  // Temporary test for gathering data.
  func testWorkspace() throws {
    let latticeScale: Float = 3
    let material: MaterialType = .elemental(.carbon)
    print(Int(latticeScale), terminator: ", ")
    
    let lattice = Lattice<Cubic> { h, k, l in
      Bounds { latticeScale * (h + k + l) }
      Material { material }
    }
    
    var reconstruction = Reconstruction()
    reconstruction.atoms = lattice.atoms
    reconstruction.material = material
    var topology = reconstruction.compile()
    PassivationTests.passivate(topology: &topology)
    print(topology.atoms.count, terminator: ", ")
    
    
    
    print()
  }
}
