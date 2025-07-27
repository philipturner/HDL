//
//  OctreeSorter+WorkSplitting.swift
//  HDL
//
//  Created by Philip Turner on 7/27/25.
//

struct WorkSplitting {
  var taskCount: Int
  var taskSizes: SIMD8<UInt8> = .zero
  var taskChildren: SIMD8<UInt64> = .zero
  
  init(childLatencies: SIMD8<Float>) {
    // Utilities for finding the task count.
    func latencyThreshold() -> Float {
      Float(150e-6)
    }
    func createMaximumTaskCount() -> Float {
      // 2.5 μs = 20 μs / 8
      let reducedThreshold = latencyThreshold() / 8
      var marks: SIMD8<Float> = .zero
      marks.replace(
        with: SIMD8(repeating: 1),
        where: childLatencies .> reducedThreshold)
      
      var output = marks.sum()
      output = max(output, 1)
      return output
    }
    func createTaskCount(maximum: Float) -> Int {
      let totalLatency = childLatencies.sum()
      
      // 20 μs task size
      var output = totalLatency / latencyThreshold()
      output.round(.toNearestOrEven)
      output = max(output, 1)
      output = min(output, maximum)
      return Int(output)
    }
    
    let maximumTaskCount = createMaximumTaskCount()
    self.taskCount = createTaskCount(maximum: maximumTaskCount)
    
    func createAssignments() -> SIMD8<UInt8> {
      if taskCount == 1 {
        return SIMD8.zero
      } else if taskCount < 8 {
        var test = WorkSplittingTest()
        test.taskCount = taskCount
        test.childCount = 8
        test.childLatencies = childLatencies
        return test.run()
      } else {
        return SIMD8(0, 1, 2, 3, 4, 5, 6, 7)
      }
    }
    
    let assignments = createAssignments()
    for childID in 0..<8 {
      let taskID = assignments[childID]
      let workItemOffset = taskSizes[Int(taskID)]
      taskSizes[Int(taskID)] = workItemOffset + 1
      
      var children = unsafeBitCast(
        taskChildren[Int(taskID)], to: SIMD8<UInt8>.self)
      children[Int(workItemOffset)] = UInt8(childID)
      taskChildren[Int(taskID)] = unsafeBitCast(
        children, to: UInt64.self)
    }
  }
}

struct WorkSplittingTest {
  var taskCount: Int = .zero
  var childCount: Int = .zero
  var childLatencies: SIMD8<Float> = .zero
  
  // Input: which task each child is assigned to.
  func taskLatencies(
    assignments: SIMD8<UInt8>
  ) -> SIMD8<Float> {
    var output: SIMD8<Float> = .zero
    for childID in 0..<childCount {
      let latency = childLatencies[childID]
      let taskID = assignments[childID]
      output[Int(taskID)] += latency
    }
    return output
  }
    
  // (task count) to the power of (child count)
  func combinationCount(childCount: Int) -> Int {
    var output: Int = 1
    for _ in 0..<childCount {
      output = output * taskCount
    }
    return output
  }
  
  // The number of combinations in the restricted algorithm, without investing
  // further effort into shrinking the combinatorial space.
  var nativeCombinationCount: Int {
    var remainingChildCount = childCount - taskCount
    remainingChildCount -= 1
    return combinationCount(
      childCount: remainingChildCount)
  }
  
  func run() -> SIMD8<UInt8> {
    let preparationStage = PreparationStage(test: self)
    
    // Declare the state variables for the best assignment.
    var bestCounter: SIMD8<UInt8>?
    var bestCounterLatency: Float = .greatestFiniteMagnitude
    
    // Iterate over all combinations of variable children.
    var counter: SIMD8<UInt8> = .zero
    let combinationCount = combinationCount(
      childCount: preparationStage.sortedChildPairs.count)
    for _ in 0..<combinationCount {
      var taskLatencies = preparationStage.fixedTaskLatencies
      for sortedChildID in preparationStage.sortedChildPairs.indices {
        let pair = preparationStage.sortedChildPairs[sortedChildID]
        let latency = pair[1]
        let taskID = counter[sortedChildID]
        taskLatencies[Int(taskID)] += latency
      }
      
      let maxTaskLatency = taskLatencies.max()
      if maxTaskLatency < bestCounterLatency {
        bestCounter = counter
        bestCounterLatency = maxTaskLatency
      }
      
      for laneID in 0..<8 {
        counter[laneID] += 1
        if counter[laneID] >= taskCount {
          counter[laneID] = 0
        } else {
          break
        }
      }
    }
    
    // Merge the fixed and variable assignments.
    guard let bestCounter else {
      fatalError("This should never happen.")
    }
    var combinedAssignments = preparationStage.fixedChildAssignments
    for sortedChildID in preparationStage.sortedChildPairs.indices {
      let pair = preparationStage.sortedChildPairs[sortedChildID]
      let childID = Int(pair[0])
      let taskID = bestCounter[sortedChildID]
      combinedAssignments[childID] = UInt8(taskID)
    }
    return combinedAssignments
  }
}

private struct PreparationStage {
  var sortedChildPairs: [SIMD2<Float>]
  var fixedChildAssignments: SIMD8<UInt8>
  var fixedTaskLatencies: SIMD8<Float>
  
  init(test: WorkSplittingTest) {
    sortedChildPairs = Self.createChildPairs(test: test)
    sortedChildPairs.sort {
      $0[1] < $1[1]
    }
    fixedChildAssignments = Self.createFixedAssignments(
      test: test,
      pairs: sortedChildPairs)
    fixedTaskLatencies = Self.createFixedLatencies(
      test: test,
      fixedAssignments: fixedChildAssignments)
    
    let fixedChildCount = Self.createFixedChildCount(test: test)
    sortedChildPairs.removeLast(fixedChildCount)
  }
  
  static func createChildPairs(
    test: WorkSplittingTest
  ) -> [SIMD2<Float>] {
    var output: [SIMD2<Float>] = []
    for childID in 0..<test.childCount {
      let latency = test.childLatencies[childID]
      let pair = SIMD2(
        Float(childID),
        latency)
      output.append(pair)
    }
    return output
  }
  
  static func createFixedAssignments(
    test: WorkSplittingTest,
    pairs: [SIMD2<Float>]
  ) -> SIMD8<UInt8> {
    var output = SIMD8<UInt8>(repeating: .max)
    guard test.childCount > test.taskCount else {
      fatalError("Invalid conditions for the restricted algorithm.")
    }
    
    // Assign the highest-index children to the lowest-index tasks.
    for taskID in 0..<test.taskCount {
      let sortedChildID = test.childCount - 1 - taskID
      let pair = pairs[sortedChildID]
      let childID = Int(pair[0])
      output[childID] = UInt8(taskID)
    }
    
    // Assign the largest of the remaining children to the highest-index task.
    let remainingChildCount = test.childCount - test.taskCount
    if test.nativeCombinationCount < 20 {
      let pair = pairs[remainingChildCount - 1]
      let childID = Int(pair[0])
      let taskID = test.taskCount - 1
      output[childID] = UInt8(taskID)
    } else {
      let pair0 = pairs[remainingChildCount - 2]
      let pair1 = pairs[remainingChildCount - 1]
      let childID0 = Int(pair0[0])
      let childID1 = Int(pair1[0])
      output[childID0] = UInt8(test.taskCount - 2)
      output[childID1] = UInt8(test.taskCount - 1)
    }
    
    return output
  }
  
  static func createFixedLatencies(
    test: WorkSplittingTest,
    fixedAssignments: SIMD8<UInt8>
  ) -> SIMD8<Float> {
    var output: SIMD8<Float> = .zero
    for childID in 0..<test.childCount {
      let latency = test.childLatencies[childID]
      let taskID = fixedAssignments[childID]
      if taskID != UInt8.max {
        output[Int(taskID)] += latency
      }
    }
    return output
  }
  
  static func createFixedChildCount(
    test: WorkSplittingTest
  ) -> Int {
    if test.nativeCombinationCount < 20 {
      return test.taskCount + 1
    } else {
      return test.taskCount + 2
    }
  }
}
