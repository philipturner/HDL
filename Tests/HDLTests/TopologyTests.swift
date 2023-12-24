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
}
