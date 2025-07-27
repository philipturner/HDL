//
//  OctreeSorter+DynamicThreaded.swift
//  HDL
//
//  Created by Philip Turner on 7/27/25.
//

import Dispatch

// compute ideal task count
//   retrieve total atom count
//   retrieve levels remaining (7 @ 4 nm)
//   compute total latency from 7.5 ns/atom/level
//   ideal task count = round_to_nearest(total latency / 20 μs)
//   restrict ideal task count to 1 to 8
//
// early returns to disallow work splitting
//   level size is 1.0 nm or smaller
//   task count is 1
//
// child count = 1 to 8
// task count ≥ child count
//   every child
//   likely at highest level of tree
//   likely leaving breadcrumbs
// task count < child count
//   continue with algorithm

// Tasks:
// - Implement work splitting, but make it single-threaded.
// - Start out with a condition: either 1 task or 8. Only split if there are
//   8 children.
// - Reorder, treating all children as nonzero. Although this is sub-optimal,
//   it simplifies the implementation for now.
// - Rearrange the recursive part into two nested loops.

extension OctreeSorter {
  // Algorithm that adaptively uses multi-threading, when a subset of the
  // octree has enough atoms.
  func mortonReorderingDynamic() -> [UInt32] {
    // Create the scratch pad.
    let scratchPad: UnsafeMutablePointer<UInt32> =
      .allocate(capacity: 8 * atoms.count)
    defer { scratchPad.deallocate() }
    
    func traverse(
      atomIDs: UnsafeMutableBufferPointer<UInt32>,
      levelOrigin: SIMD3<Float>,
      levelSize: Float
    ) {
      // Use the scratch pad.
      var childSizes: SIMD8<UInt32> = .zero
      for atomID in atomIDs {
        func createAtomOffset() -> SIMD3<Float> {
          let atom = atoms[Int(atomID)]
          let position = unsafeBitCast(atom, to: SIMD3<Float>.self)
          return position - self.origin
        }
        
        var index = SIMD3<UInt32>(repeating: 1)
        index.replace(
          with: SIMD3.zero,
          where: createAtomOffset() .< levelOrigin)
        
        let childID = (index &<< SIMD3(0, 1, 2)).wrappedSum()
        let previousSize = childSizes[Int(childID)]
        childSizes[Int(childID)] = previousSize + 1
        
        let scratchPadSlot = childID * UInt32(atomIDs.count) + previousSize
        scratchPad[Int(scratchPadSlot)] = atomID
      }
      
      // Retrieve the base pointer of the input buffer.
      func allocationPointer() -> UnsafeMutablePointer<UInt32> {
        atomIDs.baseAddress.unsafelyUnwrapped
      }
      
      // Transfer the scratch pad to the input buffer.
      // TODO: Define the cursor offset of each child before the second loop.
      do {
        var cursor: UInt32 = .zero
        for childID in 0..<UInt32(8) {
          let childSize = childSizes[Int(childID)]
          guard childSize > 0 else {
            continue
          }
          
          let newPointer = allocationPointer() + Int(cursor)
          newPointer.initialize(
            from: scratchPad + Int(childID) * atomIDs.count,
            count: Int(childSize))
          cursor += childSize
        }
      }
      if (levelSize / 2) <= Float(1 / 32) {
        return
      }
      
      // Organize the children into tasks.
      // var taskSizes: SIMD8<UInt8> = .zero
      // var taskChildren
      
      // Invoke the traversal function recursively.
      do {
        var cursor: UInt32 = .zero
        for childID in 0..<UInt32(8) {
          let newPointer = allocationPointer() + Int(cursor)
          
          let childSize = childSizes[Int(childID)]
          cursor += childSize
          if childSize <= 1 {
            continue
          }
          
          func createNewOrigin() -> SIMD3<Float> {
            let intOffset = (childID &>> SIMD3(0, 1, 2)) & 1
            let floatOffset = SIMD3<Float>(intOffset) * 2 - 1
            return levelOrigin + floatOffset * levelSize / 4
          }
          let newBufferPointer = UnsafeMutableBufferPointer(
            start: newPointer,
            count: Int(childSize))
          traverse(
            atomIDs: newBufferPointer,
            levelOrigin: createNewOrigin(),
            levelSize: levelSize / 2)
        }
      }
    }
    
    // Invoke the traversal function the first time.
    let levelOrigin = SIMD3<Float>(
      repeating: highestLevelSize / 2)
    var inPlaceBuffer = atoms.indices.map(UInt32.init)
    inPlaceBuffer.withUnsafeMutableBufferPointer { bufferPointer in
      traverse(
        atomIDs: bufferPointer,
        levelOrigin: levelOrigin,
        levelSize: highestLevelSize)
    }
    return inPlaceBuffer
  }
}
