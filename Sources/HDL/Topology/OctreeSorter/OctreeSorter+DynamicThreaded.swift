//
//  OctreeSorter+DynamicThreaded.swift
//  HDL
//
//  Created by Philip Turner on 7/27/25.
//

// compute ideal task count
//   retrieve total atom count
//   retrieve levels remaining (7 @ 4 nm)
//   compute total latency from 5.0 ns/atom/level
//     recent optimizations improved overall execution speed regardless of the
//     algorithm, so this dropped from 7.5 to 5.0
//   ideal task count = round_to_nearest(total latency / 20 μs)
//   restrict ideal task count to 1 to 8
//
// early returns to disallow work splitting
//   level size is 1.0 nm or smaller
//   task count is 1
//
// child count = 1 to 8
// task count ≥ child count
//   <-- 2.5 μs guard reduces the frequency that this branch is hit
//   every child gets a distinct task
//   likely at highest level of tree
//   likely leaving breadcrumbs
// task count < child count
//   continue with algorithm

// prevent breadcrumbs at the highest level, in similar situations that
// bypass the guard
// - find the number of children whose latency exceeds 2.5 μs
// - task count is limited to, at most, this number

// Tasks:
// - Implement work splitting, but make it single-threaded.
// - Migrate from either 1 or 8 tasks, to the correct number. Implement the
//   'restrictedTest' algorithm now.
// - Reorder, treating all children as nonzero. Although this is sub-optimal,
//   it simplifies the implementation for now.
//   - Situations with just 1 valid child will face horrible breadcrumbs.

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
      var childOffsets: SIMD8<UInt32> = .zero
      do {
        var offset: UInt32 = .zero
        for childID in 0..<8 {
          let childSize = childSizes[Int(childID)]
          guard childSize > 0 else {
            continue
          }
          
          let newPointer = allocationPointer() + Int(offset)
          newPointer.initialize(
            from: scratchPad + Int(childID) * atomIDs.count,
            count: Int(childSize))
          
          childOffsets[Int(childID)] = offset
          offset += childSize
        }
      }
      if (levelSize / 2) <= Float(1.0 / 32) {
        return
      }
      
      /*
       1 loop:
       
       atoms: 129600
       dataset    | octree |  grid
       ---------- | ------ | ------
       pre-sorted |   7534 |   6462
       lattice    |   7785 |   6691
       shuffled   |   9309 |   6682
       reversed   |   7974 |   6639
       
       2 nested loops:
       
       atoms: 129600
       dataset    | octree |  grid
       ---------- | ------ | ------
       pre-sorted |   7521 |   7039
       lattice    |   7860 |   7221
       shuffled   |   9436 |   7294
       reversed   |   7898 |   7167
       
       change:
       
       +9%
       +8%
       +9%
       +8%
       
       unexpected jump in cost, after making program output depend on result
       of work splitting:
       
       atoms: 129600
       dataset    | octree |  grid
       ---------- | ------ | ------
       pre-sorted |   7550 |   9793
       lattice    |   7808 |   9949
       shuffled   |   9311 |   9992
       reversed   |   7961 |   9897
       
       improvement after getting the compiler to inline a function:
       
       atoms: 129600
       dataset    | octree |  grid
       ---------- | ------ | ------
       pre-sorted |   7610 |   9189
       lattice    |   7790 |   9329
       shuffled   |   9427 |   9297
       reversed   |   8025 |   9212
       
       Keep this data around to track progress, as performance worsens with
       the inclusion of work splitting. Eventually, it may prove economical to
       provide an explicit 1-loop branch, earlier up in this function body.
       
       */
      
      func createLevelsRemaining() -> Int {
        let exponentFor4 = Float(4).exponentBitPattern
        let exponentForLevel = levelSize.exponentBitPattern
        return Int(exponentForLevel - exponentFor4) + 7
      }
      func createChildLatencies() -> SIMD8<Float> {
        // 5.0 ns/atom/level
        var output = SIMD8<Float>(childSizes)
        output *= Float(5.0e-9)
        output *= Float(createLevelsRemaining())
        return output
      }
      let childLatencies = createChildLatencies()
      
      func createMaximumTaskCount() -> Float {
        if levelSize <= 1 {
          return 1
        }
        
        // 2.5 μs = 20 μs / 8
        var marks: SIMD8<Float> = .zero
        marks.replace(
          with: SIMD8(repeating: 1),
          where: childLatencies .> 2.5e-6)
        
        var output = marks.sum()
        output = max(output, 1)
        return output
      }
      func createTaskCount(maximum: Float) -> Int {
        let totalLatency = childLatencies.sum()
        
        // 20 μs task size
        var output = totalLatency / Float(20e-6)
        output.round(.toNearestOrEven)
        output = max(output, 1)
        output = min(output, maximum)
        return Int(output)
      }
      
      // Child count is always 8, until we break through the barrier to entry
      // for implementing the full algorithm.
      struct WorkSplitting {
        var taskCount: Int = .zero
        var assignments: SIMD8<UInt8> = .zero
      }
      func createWorkSplitting() -> WorkSplitting {
        let maximumTaskCount = createMaximumTaskCount()
        let idealTaskCount = createTaskCount(maximum: maximumTaskCount)
        
        var output = WorkSplitting()
        if idealTaskCount > 1 {
          
          output.taskCount = 8
          output.assignments = SIMD8(0, 1, 2, 3, 4, 5, 6, 7)
        } else {
          output.taskCount = 1
          output.assignments = SIMD8.zero
        }
        return output
      }
      let workSplitting = createWorkSplitting()
      
      // Testing the impact when invoked on every call.
      do {
        var testInput = TestInput()
        testInput.taskCount = 3
        testInput.childCount = 8
        testInput.childLatencies = childLatencies
        _ = runRestrictedTest(testInput: testInput)
      }
      
      // Organize the children into tasks.
      var taskSizes: SIMD8<UInt8> = .zero
      var taskChildren: SIMD8<UInt64> = .zero
      for childID in 0..<8 {
        // taskID per child obtained from work splitting algorithm
        let taskID = workSplitting.assignments[childID]
        let offset = taskSizes[Int(taskID)]
        taskSizes[Int(taskID)] = offset + 1
        
        var children = unsafeBitCast(
          taskChildren[Int(taskID)], to: SIMD8<UInt8>.self)
        children[Int(offset)] = UInt8(childID)
        taskChildren[Int(taskID)] = unsafeBitCast(
          children, to: UInt64.self)
      }
      
      // Invoke the traversal function recursively.
      for taskID in 0..<workSplitting.taskCount {
        let size = taskSizes[taskID]
        let children = unsafeBitCast(
          taskChildren[taskID], to: SIMD8<UInt8>.self)
        
        for workItemID in 0..<size {
          let childID = children[Int(workItemID)]
          let childSize = childSizes[Int(childID)]
          if childSize <= 1 {
            continue
          }
          
          let offset = childOffsets[Int(childID)]
          let newBufferPointer = UnsafeMutableBufferPointer(
            start: allocationPointer() + Int(offset),
            count: Int(childSize))
          
          func createNewOrigin() -> SIMD3<Float> {
            let intOffset = (UInt8(childID) &>> SIMD3(0, 1, 2)) & 1
            let floatOffset = SIMD3<Float>(intOffset) * 2 - 1
            return levelOrigin + floatOffset * levelSize / 4
          }
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
