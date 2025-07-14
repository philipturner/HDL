//
//  Sort.swift
//  HDL
//
//  Created by Philip Turner on 12/2/23.
//

import Dispatch

extension Topology {
  @discardableResult
  public mutating func sort() -> [UInt32] {
    let gridSorter = GridSorter(atoms: atoms)
    let grid = gridSorter.createGrid()
    let reordering = gridSorter.mortonReordering(grid: grid)
    let previousAtoms = atoms
    
    for i in reordering.indices {
      let mortonIndex = reordering[i]
      atoms[i] = previousAtoms[Int(mortonIndex)]
    }
    
    let inverted = GridSorter.invertOrder(reordering)
    for i in bonds.indices {
      let bond = SIMD2<Int>(truncatingIfNeeded: bonds[i])
      var newBond: SIMD2<UInt32> = .zero
      newBond[0] = inverted[bond[0]]
      newBond[1] = inverted[bond[1]]
      newBond = SIMD2(newBond.min(), newBond.max())
      bonds[i] = newBond
    }
    bonds.sort {
      if $0.x != $1.x {
        return $0.x < $1.x
      } else {
        return $0.y < $1.y
      }
    }
    return inverted
  }
}

// TODO: Rename to OctreeSorter
struct GridSorter {
  var atoms: [Atom] = []
  
  var origin: SIMD3<Float>
  var dimensions: SIMD3<Float>
  
  init(atoms: [Atom]) {
    self.atoms = atoms
    
    if atoms.count == 0 {
      origin = .zero
      dimensions = .init(repeating: 0.5)
    } else {
      var minimum = SIMD3<Float>(repeating: .greatestFiniteMagnitude)
      var maximum = -minimum
      for atom in atoms {
        // @_transparent attribute is ineffective.
        let position = unsafeBitCast(atom, to: SIMD3<Float>.self)
        minimum.replace(with: position, where: position .< minimum)
        maximum.replace(with: position, where: position .> maximum)
      }
      
      origin = minimum
      dimensions = maximum - minimum
      dimensions.replace(
        with: SIMD3(repeating: 0.5),
        where: dimensions .< 0.5)
    }
  }
}

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

struct LevelSizes {
  var highest: Float
  var octreeStart: Float
  var threshold: Float
  
  init(dimensions: SIMD3<Float>) {
    // Make an initial guess of 67% for the top-level binary divider.
    let volume = dimensions.x * dimensions.y * dimensions.z
    let chunkVolume = volume / 27
    self.highest = 2 * pow(chunkVolume, 1.0 / 3)
    self.octreeStart = highest
    
    // If the grid has dimensions that vary wildly, 'highestLevelSize' does not
    // provide an accurate center. Increase it by powers of 2 until it at least
    // reaches 51% of each dimension.
    var iterationCount = 0
    while true {
      iterationCount += 1
      if iterationCount > 100 {
        fatalError("Too many iterations.")
      }
      
      let targetDimension = dimensions.max() * 0.51
      if octreeStart < targetDimension {
        octreeStart *= 2
      } else {
        break
      }
    }
    
    // Very important: deciding the granularity with which to parallelize the
    // grid. 0.5 looks way too small for sparser lattices like SiC and Si.
    self.threshold = min(2.0, max(0.5, highest * 0.51))
  }
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

extension GridSorter {
  struct Cell {
    var range: Range<Int>
    var origin: SIMD3<Float>
    var size: Float
  }
  
  struct Grid {
    var data: [UInt32] = []
    var cells: [Cell] = []
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
  
  func createGrid() -> Grid {
    let levelSizes = LevelSizes(dimensions: dimensions)
    var grid = Grid()
    
    // Create the scratch pad.
    let scratchPad: UnsafeMutablePointer<UInt32> =
      .allocate(capacity: 8 * atoms.count)
    defer { scratchPad.deallocate() }
    
    func traverseGrid(
      atomIDs: UnsafeBufferPointer<UInt32>,
      levelOrigin: SIMD3<Float>,
      levelSize: Float
    ) {
      if levelSize < levelSizes.threshold {
        let rangeStart = grid.data.count
        grid.data += atomIDs
        let rangeEnd = grid.data.count
        
        let cell = Cell(
          range: rangeStart..<rangeEnd,
          origin: levelOrigin,
          size: levelSize)
        grid.cells.append(cell)
        return
      }
      
      // Use the scratch pad.
      var childNodeCounts: SIMD8<Int> = .zero
      for atomID in atomIDs {
        @inline(__always)
        func createAtomOffset() -> SIMD3<Float> {
          // @_transparent attribute is ineffective.
          let atom = atoms[Int(atomID)]
          let position = unsafeBitCast(atom, to: SIMD3<Float>.self)
          return position - self.origin
        }
        
        var index = SIMD3<UInt32>(repeating: 1)
        index.replace(
          with: SIMD3.zero,
          where: createAtomOffset() .< levelOrigin)
        
        let key = (index &<< SIMD3(0, 1, 2)).wrappedSum()
        let previousCount = childNodeCounts[Int(key)]
        childNodeCounts[Int(key)] += 1
        scratchPad[Int(key) * atomIDs.count + previousCount] = atomID
      }
      
      // Create the temporary allocation.
      withUnsafeTemporaryAllocation(
        of: UInt32.self,
        capacity: childNodeCounts.wrappedSum()
      ) { allocationBuffer in
        @inline(__always)
        func allocationPointer() -> UnsafeMutablePointer<UInt32> {
          allocationBuffer.baseAddress.unsafelyUnwrapped
        }
        
        // Transfer the scratch pad to the temporary allocation.
        do {
          var cursor = 0
          for key in 0..<UInt32(8) {
            let childNodeCount = childNodeCounts[Int(key)]
            guard childNodeCount > 0 else {
              continue
            }
            
            let newPointer = allocationPointer() + cursor
            cursor += childNodeCount
            
            newPointer.initialize(
              from: scratchPad + Int(key) * atomIDs.count,
              count: childNodeCount)
          }
        }
        
        // Invoke the traversal function recursively.
        var cursor = 0
        for key in 0..<UInt32(8) {
          let childNodeCount = childNodeCounts[Int(key)]
          guard childNodeCount > 0 else {
            continue
          }
          
          let newPointer = allocationPointer() + cursor
          cursor += childNodeCount
          
          @inline(__always)
          func createNewOrigin() -> SIMD3<Float> {
            let intOffset = (key &>> SIMD3(0, 1, 2)) & 1
            let floatOffset = SIMD3<Float>(intOffset) * 2 - 1
            return levelOrigin + floatOffset * levelSize / 2
          }
          let newBufferPointer = UnsafeBufferPointer(
            start: newPointer,
            count: childNodeCount)
          traverseGrid(
            atomIDs: newBufferPointer,
            levelOrigin: createNewOrigin(),
            levelSize: levelSize / 2)
        }
      }
    }
    
    // Invoke the traversal function the first time.
    let levelOrigin = SIMD3<Float>(repeating: levelSizes.octreeStart)
    let initialArray = atoms.indices.map(UInt32.init)
    initialArray.withUnsafeBufferPointer { bufferPointer in
      traverseGrid(
        atomIDs: bufferPointer,
        levelOrigin: levelOrigin,
        levelSize: levelSizes.octreeStart)
    }
    guard grid.data.count == atoms.count else {
      fatalError("This should never happen.")
    }
    
    return grid
  }
  
  // Multi-threaded algorithm.
  func mortonReordering(grid: Grid) -> [UInt32] {
    nonisolated(unsafe)
    var globalOutput = [UInt32](unsafeUninitializedCapacity: atoms.count) {
      $1 = atoms.count
    }
    
    @Sendable
    func execute(cell: Cell) {
      var localOutput: [UInt32] = []
      
      // Create the scratch pad.
      let scratchPad: UnsafeMutablePointer<UInt32> =
        .allocate(capacity: 8 * cell.range.count)
      defer { scratchPad.deallocate() }
      
      func traverseTree(
        atomIDs: UnsafeBufferPointer<UInt32>,
        levelOrigin: SIMD3<Float>,
        levelSize: Float
      ) {
        if levelSize < 1 / 32 {
          localOutput += atomIDs
          return
        }
        
        // Use the scratch pad.
        var childNodeCounts: SIMD8<Int> = .zero
        for atomID in atomIDs {
          @inline(__always)
          func createAtomOffset() -> SIMD3<Float> {
            // @_transparent attribute is ineffective.
            let atom = atoms[Int(atomID)]
            let position = unsafeBitCast(atom, to: SIMD3<Float>.self)
            return position - self.origin
          }
          
          var index = SIMD3<UInt32>(repeating: 1)
          index.replace(
            with: SIMD3.zero,
            where: createAtomOffset() .< levelOrigin)
          
          let key = (index &<< SIMD3(0, 1, 2)).wrappedSum()
          let previousCount = childNodeCounts[Int(key)]
          childNodeCounts[Int(key)] += 1
          scratchPad[Int(key) * atomIDs.count + previousCount] = atomID
        }
        
        // Create the temporary allocation.
        withUnsafeTemporaryAllocation(
          of: UInt32.self,
          capacity: childNodeCounts.wrappedSum()
        ) { allocationBuffer in
          @inline(__always)
          func allocationPointer() -> UnsafeMutablePointer<UInt32> {
            allocationBuffer.baseAddress.unsafelyUnwrapped
          }
          
          // Transfer the scratch pad to the temporary allocation.
          do {
            var cursor = 0
            for key in 0..<UInt32(8) {
              let childNodeCount = childNodeCounts[Int(key)]
              guard childNodeCount > 0 else {
                continue
              }
              
              let newPointer = allocationPointer() + cursor
              cursor += childNodeCount
              
              newPointer.initialize(
                from: scratchPad + Int(key) * atomIDs.count,
                count: childNodeCount)
            }
          }
          
          // Invoke the traversal function recursively.
          var cursor = 0
          for key in 0..<UInt32(8) {
            let childNodeCount = childNodeCounts[Int(key)]
            guard childNodeCount > 0 else {
              continue
            }
            
            let newPointer = allocationPointer() + cursor
            cursor += childNodeCount
            
            if childNodeCount == 1 {
              localOutput.append(newPointer[0])
              continue
            }
            
            @inline(__always)
            func createNewOrigin() -> SIMD3<Float> {
              let intOffset = (key &>> SIMD3(0, 1, 2)) & 1
              let floatOffset = SIMD3<Float>(intOffset) * 2 - 1
              return levelOrigin + floatOffset * levelSize / 2
            }
            let newBufferPointer = UnsafeBufferPointer(
              start: newPointer,
              count: childNodeCount)
            traverseTree(
              atomIDs: newBufferPointer,
              levelOrigin: createNewOrigin(),
              levelSize: levelSize / 2)
          }
        }
      }
      
      // Invoke the traversal function the first time.
      let initialArray = grid.data[cell.range]
      initialArray.withUnsafeBufferPointer { bufferPointer in
        traverseTree(
          atomIDs: bufferPointer,
          levelOrigin: cell.origin,
          levelSize: cell.size)
      }
      guard localOutput.count == cell.range.count else {
        fatalError("This should never happen.")
      }
      
      for localID in localOutput.indices {
        let globalID = cell.range.startIndex + localID
        globalOutput[globalID] = localOutput[localID]
      }
    }
    
    func createLargeCellCount() -> Int {
      var output = 0
      for cell in grid.cells {
        let atomCount = cell.range.count
        if atomCount > 64 {
          output += 1
        }
      }
      return output
    }
    
    let taskCount = grid.cells.count
    if createLargeCellCount() < 3 {
      for z in 0..<taskCount {
        let cell = grid.cells[z]
        execute(cell: cell)
      }
    } else {
      DispatchQueue.concurrentPerform(iterations: taskCount) { z in
        let cell = grid.cells[z]
        execute(cell: cell)
      }
    }
    
    return globalOutput
  }
  
  // Single threaded algorithm without a grid.
  func invertedMortonReordering() -> [UInt32] {
    var output: [UInt32] = []
    
    // Create the scratch pad.
    let scratchPad: UnsafeMutablePointer<UInt32> =
      .allocate(capacity: 8 * atoms.count)
    defer { scratchPad.deallocate() }
    
    func traverse(
      atomIDs: UnsafeBufferPointer<UInt32>,
      levelOrigin: SIMD3<Float>,
      levelSize: Float
    ) {
      if levelSize < 1 / 32 {
        output += atomIDs
        return
      }
      
      // Use the scratch pad.
      var childNodeCounts: SIMD8<Int> = .zero
      for atomID in atomIDs {
        @inline(__always)
        func createAtomOffset() -> SIMD3<Float> {
          // @_transparent attribute is ineffective.
          let atom = atoms[Int(atomID)]
          let position = unsafeBitCast(atom, to: SIMD3<Float>.self)
          return position - self.origin
        }
        var index = SIMD3<UInt32>(repeating: 1)
        index.replace(
          with: SIMD3.zero,
          where: createAtomOffset() .< levelOrigin)
        
        let key = (index &<< SIMD3(0, 1, 2)).wrappedSum()
        let previousCount = childNodeCounts[Int(key)]
        childNodeCounts[Int(key)] += 1
        scratchPad[Int(key) * atomIDs.count + previousCount] = atomID
      }
      
      // Create the temporary allocation.
      withUnsafeTemporaryAllocation(
        of: UInt32.self,
        capacity: childNodeCounts.wrappedSum()
      ) { allocationBuffer in
        @inline(__always)
        func allocationPointer() -> UnsafeMutablePointer<UInt32> {
          allocationBuffer.baseAddress.unsafelyUnwrapped
        }
        
        // Transfer the scratch pad to the temporary allocation.
        do {
          var cursor = 0
          for key in 0..<UInt32(8) {
            let childNodeCount = childNodeCounts[Int(key)]
            guard childNodeCount > 0 else {
              continue
            }
            
            let newPointer = allocationPointer() + cursor
            cursor += childNodeCount
            
            newPointer.initialize(
              from: scratchPad + Int(key) * atomIDs.count,
              count: childNodeCount)
          }
        }
        
        // Invoke the traversal function recursively.
        var cursor = 0
        for key in 0..<UInt32(8) {
          let childNodeCount = childNodeCounts[Int(key)]
          guard childNodeCount > 0 else {
            continue
          }
          
          let newPointer = allocationPointer() + cursor
          cursor += childNodeCount
          
          if childNodeCount == 1 {
            output.append(newPointer[0])
            continue
          }
          
          @inline(__always)
          func createNewOrigin() -> SIMD3<Float> {
            let intOffset = (key &>> SIMD3(0, 1, 2)) & 1
            let floatOffset = SIMD3<Float>(intOffset) * 2 - 1
            return levelOrigin + floatOffset * levelSize / 2
          }
          let newBufferPointer = UnsafeBufferPointer(
            start: newPointer,
            count: childNodeCount)
          traverse(
            atomIDs: newBufferPointer,
            levelOrigin: createNewOrigin(),
            levelSize: levelSize / 2)
        }
      }
    }
    
    // Invoke the traversal function the first time.
    let levelSizes = LevelSizes(dimensions: dimensions)
    let levelOrigin = SIMD3<Float>(repeating: levelSizes.octreeStart)
    let initialArray = atoms.indices.map(UInt32.init)
    initialArray.withUnsafeBufferPointer { bufferPointer in
      traverse(
        atomIDs: bufferPointer,
        levelOrigin: levelOrigin,
        levelSize: levelSizes.octreeStart)
    }
    guard output.count == atoms.count else {
      fatalError("This should never happen.")
    }
    
    // Note: Do not include this part when benchmarking performance of the
    // octree sorter. It unfairly biases performance in favor of the grid
    // sorter.
    
    // TODO: Clarify what is going on with what/where the reordering happens.
    // This function doesn't have the same semantic meaning as morton
    // reordering from the grid sort.
    var reordering = [UInt32](repeating: .max, count: atoms.count)
    for reorderedID in output.indices {
      let originalID32 = output[reorderedID]
      let originalID = Int(originalID32)
      let reorderedID32 = UInt32(reorderedID)
      reordering[originalID] = reorderedID32
    }
    return reordering
  }
}
