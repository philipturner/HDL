import XCTest
@testable import HDL

final class SortTests: XCTestCase {
  func testCube() throws {
    let latticeScale: Float = 5
    let materials: [MaterialType] = [
      .elemental(.carbon),
      .checkerboard(.silicon, .carbon),
      .elemental(.silicon),
    ]
    
    // Iterate over the materials.
    for material in materials {
      // print(Int(latticeScale), terminator: ", ")
      
      let lattice = Lattice<Cubic> { h, k, l in
        Bounds { latticeScale * (h + k + l) }
        Material { material }
      }
      
      var reconstruction = Reconstruction()
      reconstruction.atoms = lattice.atoms
      reconstruction.material = material
      var topology = reconstruction.compile()
      PassivationTests.passivate(topology: &topology)
      // print(topology.atoms.count, terminator: ", ")
      
      let sorter = OctreeSorter(atoms: topology.atoms)
      // print(sorter.dimensions.max(), terminator: ", ")
      
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
      // print(numVoxels(voxelSize: 2.0), terminator: ", ")
      // print(numVoxels(voxelSize: 4.0), terminator: ", ")
      // print(numVoxels(voxelSize: 8.0), terminator: ", ")
      // print()
    }
  }
  
  // Next test:
  // - round the octree origin and dimensions to a multiple of 2 nm
  // - populate a grid of 2 nm, 4 nm, 8 nm voxels with alignment to 2 nm
  // - disallow very small lattice scales
}
