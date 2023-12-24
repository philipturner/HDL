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
  //
  // Idea for testing correctness of Topology.sort(): Repeat the same test
  // multiple times with the input randomly shuffled beforehand. Assert
  // that the output (atoms, bonds) is always the same, and that the output is
  // different from the input.
  //
  // Idea for testing correctness of Topology.match(): Ensure the matched
  // indices actually appear in ascending order of distance. Ensure none of them
  // are more than the algorithm's specified distance cutoff.
  //
  // Idea for ultimate litmus test: use the new, flexible topology compiler to
  // reconstruct (100) surfaces.
  
  // TODO: - Instead of a time-consuming, exhaustive test suite, debug this
  // visually in the renderer. The covered functionality is not as complex as
  // MM4RigidBody and doesn't need the same approach to testing. Getting (100)
  // reconstruction to work will likely trigger edge cases where the compiler is
  // broken; reproducers can be added to the test suite.
  //
  // sort() should still be debugged as explained above; this method is a form
  // of visual debugging. It may be helpful to have a utility function for
  // debugging bonds. For example, explode the crystal lattice and place marker
  // atoms on an interpolated line between actual atoms. Presence of malformed
  // bonds will be extremely obvious.
  //
  // The old 'Diamondoid' topology generator can be used to bootstrap testing
  // of how sort() treats bonds, before the remainder of the functionality is
  // working.
  //
  // Implementation plan:
  // - 1) Visualizer for Morton order and bond topology in GitHub gist.
  //   - 1.1) Test against Diamondoid and Lattice.
  //   - 1.2) Test against Topology.sort().
  // - 2) Test simple diamond and lonsdaleite lattice formation.
  // - 3) Demonstrate (100) surface and strained shell structure formation.
  // - 4) Test the generated structures in the old MM4 simulator.
  
  // TODO: - The second test case is an interesting method of forming strained
  // shell structures. Passivate a crystalline lattice, then warp and remove
  // hydrogens bonded to carbons that will merge. Validate that both this
  // and the reconstructed (100) structure are accepted by the simulator.
}
