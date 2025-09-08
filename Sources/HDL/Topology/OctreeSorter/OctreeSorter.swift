//
//  OctreeSorter.swift
//  HDL
//
//  Created by Philip Turner on 7/14/25.
//

func debugProfile(
  _ start: Double,
  _ end: Double,
  _ name: String
) {
  let elapsedTime = end - start
  let elapsedMilliseconds = elapsedTime * 1000
  let formatted = String(format: "%.3f", elapsedMilliseconds)
  print(name, formatted)
}

// Performance tests:
//
// --filter testFlatSheetScaling
//
//   - sort: 0.056 µs/atom
//   - sort: 0.047 µs/atom
//   - sort: 0.055 µs/atom
//
// --filter testSort
//
//   atoms: 16400
//   dataset    | octree |  grid
//   ---------- | ------ | ------
//   pre-sorted |    924 |    617
//   lattice    |    904 |    592
//   shuffled   |   1046 |    646
//   reversed   |    901 |    595
//
//   atoms: 54900
//   dataset    | octree |  grid
//   ---------- | ------ | ------
//   pre-sorted |   3112 |   1874
//   lattice    |   3022 |   1972
//   shuffled   |   3669 |   2162
//   reversed   |   3179 |   1925
//
//   atoms: 129600
//   dataset    | octree |  grid
//   ---------- | ------ | ------
//   pre-sorted |   7929 |   4573
//   lattice    |   8033 |   4626
//   shuffled   |   9862 |   5254
//   reversed   |   8147 |   4752
//
// The sort inside the library is not actually slower than OctreeSorter.
// Use latticeScale=20 as the go-to test for quickly checking for a regression.

// Why is one algorithm faster than the other? What's going on at the lowest
// level, and is there more room for improvement?
//
// Grid algorithm for reordering:
//
//               | Part 1    | Part 2, parallel | Part 2, serial
// ------------- | --------- | ---------------- | --------------
//   16400 atoms |  0.1 ms   |  0.3 ms          |  0.9 ms
//  129600 atoms |  2.2 ms   |  1.2 ms          |  7.9 ms
//  435600 atoms | 11.3 ms   |  3.2 ms          | 23.4 ms
// 1030400 atoms | 32.5 ms   |  9.2 ms          | 70.3 ms
//
// Octree algorithm for reordering:
//
//               | Combined Pass
// ------------- | -------------
//   16400 atoms |  0.9 ms
//  129600 atoms |  8.8 ms
//  435600 atoms | 30.2 ms
// 1030400 atoms | 82.7 ms

// New, revised metric for checking for regressions:
//
//  atoms: 129600
//  dataset    | octree |  grid
//  ---------- | ------ | ------
//  pre-sorted |   7784 |   3596
//  lattice    |   7889 |   3699
//  shuffled   |   9654 |   4510
//  reversed   |   8054 |   3702

struct OctreeSorter {
  var atoms: [Atom] = []
  
  var origin: SIMD3<Float>
  var dimensions: SIMD3<Float>
  
  init(atoms: [Atom]) {
    self.atoms = atoms
    
    if atoms.count == 0 {
      origin = .zero
      dimensions = SIMD3(repeating: 1)
    } else {
      var minimum = SIMD3<Float>(repeating: .greatestFiniteMagnitude)
      var maximum = -minimum
      for atom in atoms {
        // @_transparent attribute is ineffective.
        let position = unsafeBitCast(atom, to: SIMD3<Float>.self)
        minimum.replace(with: position, where: position .< minimum)
        maximum.replace(with: position, where: position .> maximum)
      }
      minimum.round(.down)
      maximum.round(.up)
      
      origin = minimum
      dimensions = maximum - minimum
      dimensions.replace(
        with: SIMD3(repeating: 1),
        where: dimensions .< 1)
    }
  }
  
  static func invertOrder(_ input: [UInt32]) -> [UInt32] {
    return [UInt32](unsafeUninitializedCapacity: input.count) {
      $1 = input.count
      let baseAddress = $0.baseAddress.unsafelyUnwrapped
      
      for reorderedID in input.indices {
        let originalID = Int(input[reorderedID])
        baseAddress[originalID] = UInt32(reorderedID)
      }
    }
  }
}

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
