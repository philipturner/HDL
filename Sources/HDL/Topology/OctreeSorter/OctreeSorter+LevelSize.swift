//
//  LevelSizes.swift
//  HDL
//
//  Created by Philip Turner on 7/14/25.
//

import func Foundation.pow

// TODO: Change the algorithm to make the level size a power of 2? And
// clip the lower corner to a good power of 2? The current algorithm
// provides too erratic of performance, and something similar to chaos or
// nondeterminism in the results. The root cause is that the origin can
// vary drastically, depending on where the user offsets the lattice. I
// designed an algorithm that's robust to this fact, and some of the tests
// now depend on its exact behavior. We must revise it during the current
// round of source-breaking API changes.
//
// Implement this fix once the entire 'Sort' file is refactored and easier
// to work with, without introducing new regressions.
//
// New ideas:
// - Defer this to until you have a functioning renderer, to debug glitches
//   where visualizing the atoms helps massively.
// - Reduce the task count when the number of cells is very large,
//   merging nearby cells in the list.
//   - This creates a more homogeneous atom count and should reduce the spike
//     in overhead.
//   - Before doing this, thoroughly de-mystify the latency as a function of
//     lattice size, and have specific benchmarks (in milliseconds) to refer
//     to.
//     - Try doing the optimization, seeing whether there's an instant
//       noticeable improvement, then commenting out the optimization and
//       investigating more thoroughly. This minimizes the risk of wasted
//       time, if the optimization is ineffective.
//
// Choice: defer changes to the underlying algorithm until we have a working
// renderer on Windows. Instead, include an attempt to improve workload
// distribution in this PR.



// Clarified plans for new algorithm:
// - Fall back to OctreeSorter for small systems
// - Larger systems tuned for assumed ~50-176 nm^3 atom density
// - Assume (2 nm)^3 granularity of dense regions, in the case of very sparse
//   structures. This is a good choice used in the renderer.
// - Benchmark small, sparse, and highly anisotropic shapes. Prove that the new
//   algorithm serves these better/equal to the old one. The new renderer isn't
//   needed to implement these tests.
//   - Three new, distinct performance tests for sorting:
//     - small, between 100 and 10000 atoms
//     - sparse + highly anisotropic L shape
//     - 1 micron offset from the world origin
//   - All tests should include hydrogen passivation, which changes the
//     distribution of atom density.
//   - Leave the old test as-is, and create new tests from scratch that make
//     benchmarking easier. Include a fresh version of the test for large
//     problem sizes.
//     - Create a file, 'SortTests'.
//     - Figure out how to expose internal APIs in the current compilation
//       mode, for benchmarking. This can be a temporary feature.
//
// molecular-renderer is able to efficiently parallelize the radix sort when
// constructing an octree. HDL cannot exploit the restriction of atom positions
// to a ~(256 nm)^3, origin-centered world volume. There are fundamental
// differences between how the two codes must achieve a similar task.
//
// # Overview
//
// Existing algorithm has a threshold at 1-2 nm. This level size is the offset
// of a child node from the parent node. The offset is always half the child
// node side length. So the algorithm terminates when nodes are 2-4 nm large.
//
// Deepest level of the octree terminates when level size is 1/64 to 1/32 nm.
// This means the nodes are 1/32 to 1/16 nm large.
//
// Analyze the existing algorithm, accounting for the two extremes of level
// size and atom density. Model the target use cases of latticeScale=5 to 40.
// Note that performance should be good slightly outside this range as well.
//
// lattice scale | atom count | size (C) | size (Si) |
// ------------- | ---------- | -------- | --------- |
//             5 |       2100 |   2.3 nm |    3.5 nm |
//            10 |      16400 |   4.5 nm |    6.9 nm |
//            20 |     129600 |   9.0 nm |   13.7 nm |
//            30 |     435600 |  13.5 nm |   20.6 nm |
//            40 |    1030400 |  18.0 nm |   27.4 nm |
//
// atoms per cell | C     | Si    |
// -------------- | ----- | ----- |
//         2.0 nm |  1408 |   400 |
//         4.0 nm | 11264 |  3200 |
//
// The latency per atom might be very high. To achieve 20 µs task time with
// 0.05 µs/atom, that evaluates to 400 atoms/task. If the data is misleading
// because of multicore parallelism, that is instead 50 atoms/task. Most light,
// O(n) algorithms use a conservative (large) task size of 5000 atoms.
//
// The vast majority of cells have a size close to the average of all cells in
// the grid. So, in most cases, the presence of small outlier tasks does not
// worsen the overhead of task distribution. I hypothesize that smaller task
// size for Si + 2.0 nm is the principal issue. Multithreading is needed to
// make the grid algorithm competitive with octree, so even very small grids
// (e.g. task size 64) must use it to remain competitive. This fact leaked into
// the decision for a variable grid threshold size.
//
// # Physical justification for lowest level of the octree
//
// The center of the current node is (node size / 2, ...). Child nodes are
// placed at (node size / 4, ...) and (node size * 3 / 4, ...). This explains
// the confusion in what "level size" means.
//
// Child offset | Level size | Node size | Node volume |
// ------------ | ---------- | --------- | ----------- |
//   1 / 128 nm |  1 / 64 nm | 1 / 32 nm | 3.1e-5 nm^3 |
//    1 / 64 nm |  1 / 32 nm | 1 / 16 nm | 2.4e-4 nm^3 |
//       0.5 nm |     1.0 nm |    2.0 nm |      8 nm^3 |
//       1.0 nm |     2.0 nm |    4.0 nm |     64 nm^3 |
//
// Atom | Radius | Volume      | 1 / Atom Density |
// ---- | ------ | ----------- | ---------------- |
//    H |  31 pm | 1.2e-4 nm^3 |                  |
//    C |  66 pm | 1.2e-3 nm^3 |      5.7e-3 nm^3 |
//   Si | 111 pm | 5.7e-3 nm^3 |      2.0e-2 nm^3 |
//   Au | 136 pm | 1.1e-2 nm^3 |      1.7e-2 nm^3 |
//
// For two atoms to occupy the same voxel at the lowest level, their bond
// length must be under 54 pm (1 / 32 nm) or 108 pm (1 / 16 nm). Since the
// lowest level will be highly parallelized, we can afford an extra level for
// security.
//
// # Division of O(nlogn) work between serial and parallel stages
//
// Schemes:
// - current algorithm, low end:        2.0 nm + 1 / 32 nm
// - current algorithm, high end:       4.0 nm + 1 / 16 nm
// - possible choice for new algorithm: 4.0 nm + 1 / 32 nm
//
// Analyzing latticeScale vs material vs size in a Google Sheet.

// TODO: Implementation Plan
//
// Phase I:   Merge OctreeSorter into the main code base. Expose it to the
//            test module with @testable.
// Phase II:  Provide a means to easily switch between old and new algorithms.
// Phase III: Test various optimizations with the existing 'testSort'. The
//            latency data is very accessible, but crude and with low coverage.
// Phase IV:  Expand the test suite to the cases outlined further above. This
//            will provide a means for validating correctness in future, when
//            it is possible to render.
//
// Might have to tweak the objectives for Phase IV, as the first phases are
// checked off and more unknowns are resolved.

// TODO: Adjust all level sizes to be 2x larger, in all function bodies.

extension OctreeSorter {
  var highestLevelSize: Float {
    func truncatedDimensions() -> SIMD3<Float> {
      SIMD3<Float>(SIMD3<Int>(dimensions))
    }
    guard all(dimensions .== truncatedDimensions()) else {
      fatalError("Dimensions must be integers.")
    }
    
    var output: Float = 4
    for _ in 0..<100 {
      if output < dimensions.max() {
        output = 2 * output
      } else {
        return output
      }
    }
    
    fatalError("Took too many iterations to find highest level size.")
  }
}
