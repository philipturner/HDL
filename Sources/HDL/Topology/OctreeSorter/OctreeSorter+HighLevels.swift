//
//  OctreeSorter+HighLevels.swift
//  HDL
//
//  Created by Philip Turner on 9/7/25.
//

import Dispatch

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

// Instead of passing an offset'd pointer into each child function call, as
// with previous algorithms (and likely the lower levels), we retain a
// reference to the global base address. Instead, we identify the location
// by a range. This range can recalculate the atom count too, simplifying the
// amount of data passed between functions.
//
// 'Cell' owns the range.

// ## Loop Over Levels of the Tree
//
// levelSize specified in this loop, which is the source of truth. No
// division by 2 for child functions. This is a point of divergence with the
// lower levels. The 'traverse' function should accept levelSize as an
// argument.
//
// High-level data structure containing [Thread]:
//
// Start with a single Thread, and a single Cell scoped to the entire octree.
// Invoke this thread with DispatchQueue.concurrentPerform, to simplify the
// code. All passes now use concurrent dispatch.
//
// Inside each dispatch task, the number of cells is unpacked. Create a new
// list of output cells. Sort them into ones belonging to the parent thread,
// and a separate array of child Thread objects. Write the results to a
// global buffer in a thread-safe way.
//
// In a post-processing part of the code, after the concurrent dispatch,
// scan the global buffer of threads and cells. Eliminate any parent threads
// that now have zero cells. Eliminate any allocated child threads that were
// never materialized.
//
// Thread-safe data structure for results:
//
// Order of cells in a cell buffer could be inconsequential. A cell's range
// always points to the same atoms, no matter where you store the cell in
// memory. What's important is the Thread, which must reference the cell
// objects.
//
// This process is especially tricky, because one thread could reference many
// non-contiguous patches of atoms in the list. Does this mean many
// non-contiguous patches of Cell objects? Not exactly.
//
// We can probably use a simple sequential scan-compact algorithm. Process
// each thread in order from first to last. Update the thread based on some
// running counters. Forget the current Thread object and move on to the next.

extension OctreeSorter {
  struct Cell {
    var range: Range<Int> = 0..<0
    var origin: SIMD3<Float> = .zero
  }
  
  struct Thread {
    var cells: [Cell] = []
  }
  
  struct WorkDistribution {
    // WARNING: Receiving function must explicitly deallocate this buffer,
    // otherwise there will be a memory leak. It appears that in all use cases,
    // we can theoretically get by without ever spawning a separate Array.
    // Keep this in mind and eventually implement this optimization.
    //
    // Could wrap this in a 'class' for automatic memory management, while
    // still having decent internal APIs. Perhaps wrap the output of the
    // function for low levels.
    var inPlaceBuffer: [UInt32] = []
    var threads: [Thread] = []
  }
  
  func mortonReorderingHighLevels() -> WorkDistribution {
    nonisolated(unsafe)
    let scratchBuffer: UnsafeMutablePointer<UInt32> =
      .allocate(capacity: 8 * atoms.count)
    defer { scratchBuffer.deallocate() }
    
    nonisolated(unsafe)
    var threads: [Thread] = [createFirstThread()]
    
    nonisolated(unsafe)
    var inPlaceBuffer = atoms.indices.map(UInt32.init)
    inPlaceBuffer.withUnsafeMutableBufferPointer { inPlaceBufferPointer in
      nonisolated(unsafe)
      var levelSize = highestLevelSize
      while levelSize > 1 {
        // Perform a prefix sum, to allocate memory for outputs.
        let threadCellOffsets = Self.cellOffsets(threads: threads)
        let inputCellCount = threads.map(\.cells.count).reduce(0, +)
        
        // Thread-safe buffers for storing results.
        nonisolated(unsafe)
        var results = LevelResults(
          threadCount: threads.count,
          inputCellCount: inputCellCount)
        
        DispatchQueue.concurrentPerform(
          iterations: threads.count
        ) { threadID in
          let thread = threads[threadID]
          let inputPrefixSum = Int(threadCellOffsets[threadID])
          
          var parentCells: [Cell] = []
          var children: [Thread] = []
          for cell in thread.cells {
            let output = Self.traverseHighLevel(
              inPlaceBuffer: inPlaceBufferPointer.baseAddress!,
              scratchBuffer: scratchBuffer,
              cell: cell,
              levelSize: levelSize)
            guard output.count > 0 else {
              fatalError("Unexpected output cell count.")
            }
            
            if output.count == 1 {
              parentCells += output[0].cells
            } else {
              children += output
            }
          }
          
          // Write the parent cells into the results.
          results.outputCellsPerParent[threadID] = UInt32(parentCells.count)
          for cellID in parentCells.indices {
            let cell = parentCells[cellID]
            let cellOffset = inputPrefixSum * 8 + cellID
            results.outputParentCells[cellOffset] = cell
          }
          
          // Write the children into the results.
          results.childrenPerParent[threadID] = UInt32(children.count)
          var childPrefixSum: Int = .zero
          for childID in children.indices {
            let child = children[childID]
            let childOffset = inputPrefixSum * 8 + childID
            let cellCount = child.cells.count
            results.outputCellsPerChild[childOffset] = UInt32(cellCount)
            
            for cellID in child.cells.indices {
              let cell = child.cells[cellID]
              let cellOffset = inputPrefixSum * 8 + childPrefixSum + cellID
              results.outputChildCells[cellOffset] = cell
            }
            childPrefixSum += cellCount
          }
        }
        
        // Reconstruct 'parentCells' and 'children', scan-compact the list.
        threads = results.nextLevel(
          threadCellOffsets: threadCellOffsets)
        levelSize /= 2
      }
    }
    
    fatalError("Not implemented.")
  }
}

extension OctreeSorter {
  private func createFirstThread() -> Thread {
    var cell = Cell()
    cell.range = atoms.indices
    cell.origin = SIMD3<Float>(
      repeating: highestLevelSize / 2)
    
    let thread = Thread(cells: [cell])
    return thread
  }
  
  private static func cellOffsets(threads: [Thread]) -> [UInt32] {
    var output: [UInt32] = []
    var counter: UInt32 = .zero
    for thread in threads {
      output.append(counter)
      counter += UInt32(thread.cells.count)
    }
    return output
  }
  
  private struct LevelResults {
    var outputCellsPerParent: [UInt32]
    var outputParentCells: [Cell]
    
    var childrenPerParent: [UInt32]
    var outputCellsPerChild: [UInt32]
    var outputChildCells: [Cell]
    
    init(
      threadCount: Int,
      inputCellCount: Int
    ) {
      outputCellsPerParent = [UInt32](
        repeating: 0, count: threadCount)
      outputParentCells = [Cell](
        repeating: Cell(), count: 8 * inputCellCount)
      
      childrenPerParent = [UInt32](
        repeating: 0, count: threadCount)
      outputCellsPerChild = [UInt32](
        repeating: 0, count: 8 * inputCellCount)
      outputChildCells = [Cell](
        repeating: Cell(), count: 8 * inputCellCount)
    }
    
    func nextLevel(threadCellOffsets: [UInt32]) -> [Thread] {
      var output: [Thread] = []
      for threadID in outputCellsPerParent.indices {
        let inputPrefixSum = Int(threadCellOffsets[threadID])
        
        // Reconstruct 'parentCells'.
        let parentCellCount = Int(outputCellsPerParent[threadID])
        var parentCells: [Cell] = []
        for cellID in 0..<parentCellCount {
          let cellOffset = inputPrefixSum * 8 + cellID
          let cell = outputParentCells[cellOffset]
          parentCells.append(cell)
        }
        
        // Reconstruct 'children'.
        let childCount = Int(childrenPerParent[threadID])
        var children: [Thread] = []
        var childPrefixSum: Int = .zero
        for childID in 0..<childCount {
          let childOffset = inputPrefixSum * 8 + childID
          let cellCount = Int(outputCellsPerChild[childOffset])
          
          var childCells: [Cell] = []
          for cellID in 0..<cellCount {
            let cellOffset = inputPrefixSum * 8 + childPrefixSum + cellID
            let cell = outputChildCells[cellOffset]
            childCells.append(cell)
          }
          childPrefixSum += cellCount
          
          let child = Thread(cells: childCells)
          children.append(child)
        }
        
        // Compact the vacant list entries.
        if parentCells.count > 0 {
          let parent = Thread(cells: parentCells)
          output.append(parent)
        }
        output += children
      }
      return output
    }
  }
  
  @Sendable
  private static func traverseHighLevel(
    inPlaceBuffer: UnsafeMutablePointer<UInt32>,
    scratchBuffer: UnsafeMutablePointer<UInt32>,
    cell: Cell,
    levelSize: Float
  ) -> [Thread] {
    fatalError("Not implemented.")
  }
}
