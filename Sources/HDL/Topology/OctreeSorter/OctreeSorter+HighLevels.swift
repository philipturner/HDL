//
//  OctreeSorter+HighLevels.swift
//  HDL
//
//  Created by Philip Turner on 9/7/25.
//

extension OctreeSorter {
  struct Cell {
    var range: Range<Int>
    var origin: SIMD3<Float>
  }
  
  // WARNING: This is not thread-safe. Need to prepare a different data
  // organization scheme to make it thread-safe.
  struct Thread {
    var cells: [Cell]
  }
  
  // Pseudocode of traversal algorithm:
  //
  // Each parent thread gets allocated n * 8 child threads, where n is the cell
  // count in the parent thread. Each possible child up to 8 cells allocated.
  //
  // Each parent thread gets allocated a destination region of n * 8 cells,
  // where n is the cell count in the parent thread.
  //
  // All space is allocated upfront in a large array, and written to in a
  // thread-safe manner. The Swift 'Array' data type is never used.
  
  // After a pass in the high levels, the source thread and source cell objects
  // must be deleted. In the low levels, the source thread doesn't perish
  // because we just directly recurse all remaining levels of the octree in
  // order.
  //
  // A function takes a single cell in, and returns an array of threads. If
  // the array has count one, it doesn't fork off the parent thread. The
  // overhead of array creation is relatively negligible at the high levels,
  // but vastly simplifies the implementation.
  //
  // The in-place atom buffer survives across passes / levels. It is
  // overwritten in-place by several threads, concurrently. While individual
  // threads may not access it contiguously, the memory operations should
  // theoretically be atomic / safe.
  //
  // A single task on the dispatch queue will be multiple function calls into
  // the tree traversal function. However, on first glance, the structure of
  // the calling loop might make it possible to fully inline. Don't worry
  // about ensuring it's fully inlined, as the overhead at high levels is low.
}
