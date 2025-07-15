import XCTest
import HDL

final class SortTests: XCTestCase {
  // Temporary test for gathering data.
  func testWorkspace() throws {
    let latticeScale: Float = 3
    let material: MaterialType = .elemental(.carbon)
    
    let lattice = Lattice<Cubic> { h, k, l in
      Bounds { latticeScale * (h + k + l) }
      Material { material }
    }
    print(lattice.atoms.count)
    
    var reconstruction = Reconstruction()
    reconstruction.atoms = lattice.atoms
    reconstruction.material = material
    var topology = reconstruction.compile()
    print(topology.atoms.count)
    
    PassivationTests.passivate(topology: &topology)
    print(topology.atoms.count)
  }
}
