//
//  OctreeSorter+WorkSplitting.swift
//  HDL
//
//  Created by Philip Turner on 7/27/25.
//

func runRestrictedTest(
  testInput: TestInput
) -> SIMD8<UInt8> {
  let preparationStage = PreparationStage(testInput: testInput)
  
  // Declare the state variables for the best assignment.
  var bestCounter: SIMD8<UInt8>?
  var bestCounterLatency: Float = .greatestFiniteMagnitude
  
  // Iterate over all combinations of variable children.
  var counter: SIMD8<UInt8> = .zero
  let combinationCount = testInput.combinationCount(
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
      if counter[laneID] >= testInput.taskCount {
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

struct TestInput {
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
}

private struct PreparationStage {
  var sortedChildPairs: [SIMD2<Float>]
  var fixedChildAssignments: SIMD8<UInt8>
  var fixedTaskLatencies: SIMD8<Float>
  
  init(testInput: TestInput) {
    sortedChildPairs = Self.createChildPairs(testInput: testInput)
    sortedChildPairs.sort {
      $0[1] < $1[1]
    }
    fixedChildAssignments = Self.createFixedAssignments(
      testInput: testInput,
      pairs: sortedChildPairs)
    fixedTaskLatencies = Self.createFixedLatencies(
      testInput: testInput,
      fixedAssignments: fixedChildAssignments)
    
    let fixedChildCount = Self.createFixedChildCount(testInput: testInput)
    sortedChildPairs.removeLast(fixedChildCount)
  }
  
  static func createChildPairs(
    testInput: TestInput
  ) -> [SIMD2<Float>] {
    var output: [SIMD2<Float>] = []
    for childID in 0..<testInput.childCount {
      let latency = testInput.childLatencies[childID]
      let pair = SIMD2(
        Float(childID),
        latency)
      output.append(pair)
    }
    return output
  }
  
  static func createFixedAssignments(
    testInput: TestInput,
    pairs: [SIMD2<Float>]
  ) -> SIMD8<UInt8> {
    var output = SIMD8<UInt8>(repeating: .max)
    guard testInput.childCount > testInput.taskCount else {
      fatalError("Invalid conditions for the restricted algorithm.")
    }
    
    // Assign the highest-index children to the lowest-index tasks.
    for taskID in 0..<testInput.taskCount {
      let sortedChildID = testInput.childCount - 1 - taskID
      let pair = pairs[sortedChildID]
      let childID = Int(pair[0])
      output[childID] = UInt8(taskID)
    }
    
    // Assign the largest of the remaining children to the highest-index task.
    let remainingChildCount = testInput.childCount - testInput.taskCount
    if testInput.nativeCombinationCount < 20 {
      let pair = pairs[remainingChildCount - 1]
      let childID = Int(pair[0])
      let taskID = testInput.taskCount - 1
      output[childID] = UInt8(taskID)
    } else {
      let pair0 = pairs[remainingChildCount - 2]
      let pair1 = pairs[remainingChildCount - 1]
      let childID0 = Int(pair0[0])
      let childID1 = Int(pair1[0])
      output[childID0] = UInt8(testInput.taskCount - 2)
      output[childID1] = UInt8(testInput.taskCount - 1)
    }
    
    return output
  }
  
  static func createFixedLatencies(
    testInput: TestInput,
    fixedAssignments: SIMD8<UInt8>
  ) -> SIMD8<Float> {
    var output: SIMD8<Float> = .zero
    for childID in 0..<testInput.childCount {
      let latency = testInput.childLatencies[childID]
      let taskID = fixedAssignments[childID]
      if taskID != UInt8.max {
        output[Int(taskID)] += latency
      }
    }
    return output
  }
  
  static func createFixedChildCount(
    testInput: TestInput
  ) -> Int {
    if testInput.nativeCombinationCount < 20 {
      return testInput.taskCount + 1
    } else {
      return testInput.taskCount + 2
    }
  }
}
