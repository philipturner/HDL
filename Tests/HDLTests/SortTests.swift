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
    
    let sorter = OctreeSorter(atoms: topology.atoms)
    print(sorter.dimensions.max(), terminator: ", ")
    
    func numVoxels(voxelSize: Float) -> Int {
      var output: Int = 1
      for laneID in 0..<3 {
        var dimension = sorter.dimensions[laneID]
        dimension /= voxelSize
        dimension.round(.up)
        output *= Int(exactly: dimension)!
      }
      return output
    }
    print(numVoxels(voxelSize: 2.0), terminator: ", ")
    print(numVoxels(voxelSize: 4.0), terminator: ", ")
    print(numVoxels(voxelSize: 8.0), terminator: ", ")
    print()
  }
}
