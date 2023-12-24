import XCTest
#if DEBUG
@testable import HDL
#else
import HDL
#endif

final class TopologyTests: XCTestCase {
  func testTopology() throws {
    //  func testTopologyInit() {
    //    let lattice = Lattice<Cubic> { h, k, l in
    //      Bounds { 4 * h + 4 * k + 4 * l }
    //      Material { .elemental(.carbon) }
    //    }
    //
    //    let topology1 = Topology(lattice.entities)
    //    _ = topology1
    //
    //    let topology2 = Topology(lattice.entities.map {
    //      var copy = $0
    //      copy.position -= 0.5
    //      return copy
    //    })
    //    _ = topology2
    //  }
  }
  
  // Idea for testing correctness of TopologyGrid.mortonReordering: render a
  // trail of interpolated points between each atom in the list. If the Morton
  // reordering is correct, it will look like a Z-order curve.
  
  // Idea for testing correctness of Topology.sort(): Repeat the same test
  // multiple times with the input randomly shuffled beforehand. Assert
  // that the output (atoms, bonds) is always the same, and that the output is
  // different from the input.
  
  // Idea for testing correctness of Topology.match(): Ensure the matched
  // indices actually appear in ascending order of distance. Ensure none of them
  // are more than the algorithm's specified distance cutoff.
}
