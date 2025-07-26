@testable import HDL
import XCTest

final class SortTests: XCTestCase {
  static let printPerformanceSummary = true
  
  // Expected performance on original benchmarked computer (M1 Max):
  //
  // bounds | atoms  | octree | serial | parallel | speedup
  // ------ | ------ | ------ | ------ | -------- | ----------
  // 10     |  16400 |   2370 |   1555 |      905 | 1.5 -> 2.6
  
  // With latest version of the code base (2 years later) + Swift 5.8:
  //
  // atoms: 16400
  // dataset    | octree | serial | parallel
  // ---------- | ------ | ------ | --------
  // pre-sorted |   2060 |   1597 |    872
  // lattice    |   2072 |   1571 |    907
  // shuffled   |   2341 |   1697 |    898
  // reversed   |   2069 |   1590 |    830
  
  // After removing inlining bottlenecks from Sort:
  //
  // atoms: 16400
  // dataset    | octree | serial | parallel
  // ---------- | ------ | ------ | --------
  // pre-sorted |   1864 |   1070 |    728
  // lattice    |   1808 |   1073 |    719
  // shuffled   |   2192 |   1234 |    676
  // reversed   |   1899 |   1017 |    711
  
  func testDiagonalOrder() throws {
    func createLattice() -> Lattice<Cubic> {
      Lattice<Cubic> { h, k, l in
        Bounds { 10 * (h + k + l) }
        Material { .checkerboard(.silicon, .carbon) }
      }
    }
    func transform(lattice: Lattice<Cubic>) -> [Atom] {
      var output: [Atom] = []
      for atomID in lattice.atoms.indices {
        var atom = lattice.atoms[atomID]
        atom.position = -atom.position
        output.append(atom)
      }
      return output
    }
    func sort(atoms: [Atom]) -> [Atom] {
      var topology = Topology()
      topology.atoms = atoms
      topology.sort()
      return topology.atoms
    }
    
    let lattice = createLattice()
    let transformedAtoms = transform(lattice: lattice)
    let sortedAtoms = sort(atoms: transformedAtoms)
    
    var diagonalAtoms: [SIMD4<Float>] = []
    for atom in sortedAtoms {
      guard (atom.x - atom.y).magnitude < 0.001,
            (atom.y - atom.z).magnitude < 0.001 else {
        continue
      }
      diagonalAtoms.append(atom)
    }
    XCTAssertEqual(diagonalAtoms.count, 21)
    
    for i in 1..<diagonalAtoms.count {
      let previousAtom = diagonalAtoms[i - 1]
      let currentAtom = diagonalAtoms[i]
      
      for laneID in 0..<3 {
        let previousCoordinate = previousAtom[laneID]
        let currentCoordinate = currentAtom[laneID]
        XCTAssertLessThan(previousCoordinate, currentCoordinate)
      }
    }
  }
  
  func testSortPerformance() throws {
    // Revert to 10 and true after any refactorings
    let latticeScale: Float = 20
    let testParallel = Bool.random() ? false : false
    let lattice = Lattice<Hexagonal> { h, k, l in
      let h2k = h + 2 * k
      Bounds { latticeScale * (2 * h + h2k + l) }
      Material { .elemental(.carbon) }
    }
    
    var output: [String] = []
    if testParallel {
      output.append("dataset    | octree | serial | parallel")
      output.append("---------- | ------ | ------ | --------")
    } else {
      output.append("dataset    | octree |  grid ")
      output.append("---------- | ------ | ------")
    }
    
    struct Trial {
      var atoms: [Atom] = []
      var name: String = ""
      
      init(lattice: Lattice<Hexagonal>, index: Int) {
        switch index {
        case 0:
          var topology = Topology()
          topology.atoms = lattice.atoms
          topology.sort()
          
          atoms = topology.atoms
          name = "pre-sorted"
        case 1:
          atoms = lattice.atoms
          name = "lattice   "
        case 2:
          atoms = lattice.atoms.shuffled()
          name = "shuffled  "
        case 3:
          atoms = lattice.atoms.reversed()
          name = "reversed  "
        default:
          fatalError("This should never happen.")
        }
      }
    }
    
    // Revert to 0..<4 after any refactorings
    for trialID in 0..<4 {
      let trial = Trial(lattice: lattice, index: trialID)
      
      let startParallel = Profiler.time()
      nonisolated(unsafe)
      var resultGrid1: [UInt32] = []
      nonisolated(unsafe)
      var resultGrid2: [UInt32] = []
      if testParallel {
        DispatchQueue.concurrentPerform(iterations: 2) { z in
          var topology = Topology()
          topology.atoms = trial.atoms
          let resultGrid = topology.sort()
          if z == 0 {
            resultGrid1 = resultGrid
          } else {
            resultGrid2 = resultGrid
          }
        }
      }
      let endParallel = Profiler.time()
      
      let startGrid = Profiler.time()
      var topology = Topology()
      topology.atoms = trial.atoms
      let resultGrid = topology.sort()
      if testParallel {
        var topology = Topology()
        topology.atoms = trial.atoms
        _ = topology.sort()
      }
      let endGrid = Profiler.time()
      
      let startOctree = Profiler.time()
      var resultOctree: [UInt32]
      do {
        let sorter = OctreeSorter(atoms: trial.atoms)
        let reordering = sorter.mortonReordering()
        resultOctree = OctreeSorter.invertOrder(reordering)
      }
      if testParallel {
        let sorter = OctreeSorter(atoms: trial.atoms)
        let reordering = sorter.mortonReordering()
        _ = OctreeSorter.invertOrder(reordering)
      }
      let endOctree = Profiler.time()
      
      XCTAssertEqual(resultGrid, resultOctree)
      if testParallel {
        XCTAssertEqual(resultGrid1, resultOctree)
        XCTAssertEqual(resultGrid2, resultOctree)
      }
      
      let usParallel = Int((endParallel - startParallel) * 1e6)
      let usGrid = Int((endGrid - startGrid) * 1e6)
      let usOctree = Int((endOctree - startOctree) * 1e6)
      
      var reprParallel = "\(usParallel)"
      var reprGrid = "\(usGrid)"
      var reprOctree = "\(usOctree)"
      
      while reprParallel.count < 6 {
        reprParallel = " \(reprParallel)"
      }
      while reprGrid.count < 6 {
        reprGrid = " \(reprGrid)"
      }
      while reprOctree.count < 6 {
        reprOctree = " \(reprOctree)"
      }
      if testParallel {
        output.append("\(trial.name) | \(reprOctree) | \(reprGrid) | \(reprParallel)")
      } else {
        output.append("\(trial.name) | \(reprOctree) | \(reprGrid)")
      }
    }
    
    if Self.printPerformanceSummary {
      print()
      print("atoms:", lattice.atoms.count)
      for line in output {
        print(line)
      }
    }
    
    // During Topology.match, there are some situations where two similarly
    // sized grids will be constructed. They might be the exact same, although
    // the program can't detect that fact in a generalizable/robust manner.
    // Parallelization offers a simpler alternative that, based on the data
    // below, provides about the same speedup as eliding the compute work.
    
    // 'lattice' configuration, serial
    //
    // bounds | atoms  | octree |  0.25 |  0.5 |    1 |    2 |    4 | optimized
    // ------ | ------ | ------ | ----- | ---- | ---- | ---- | ---- | ----------
    // 5      |   2100 |    136 |   286 |  142 |  146 |      |      |  174
    // 7      |   5684 |    411 |   472 |  299 |  315 |  508 |      |  293
    // 10     |  16400 |   1168 |  2345 |  866 |  698 |  686 | 1276 |  887
    // 14     |  44688 |   3333 |  3447 | 2122 | 1863 | 1775 | 3512 | 1891
    // 20     | 129600 |   9245 | 19899 | 6695 | 5882 | 5332 | 4959 | 5403
    
    // 'lattice' configuration, 2x duplicated
    //
    // bounds | atoms  | octree | serial | parallel | speedup
    // ------ | ------ | ------ | ------ | -------- | ----------
    // 5      |   2100 |    314 |    298 |      186 | 1.1 -> 1.7
    // 7      |   5684 |    750 |    562 |      344 | 1.3 -> 2.2
    // 10     |  16400 |   2370 |   1555 |      905 | 1.5 -> 2.6
    // 14     |  44688 |   6085 |   3789 |     2160 | 1.6 -> 2.8
    // 20     | 129600 |  19932 |  10811 |     6567 | 1.8 -> 3.0
  }
  
  // This test will prototype the work splitting algorithm, ensuring it
  // executes in a reasonable amount of time.
  //
  // Once the tests are in place, we can try optimizations without causing
  // correctness regressions.
  // - Current execution time: ~1.0-3.5 μs, depending on problem size
  //
  // Tasks:
  // - Perform small cleanups that coincidentally optimize performance, for the
  //   restricted algorithm.
  // - Remove the "restricted1" algorithm variant. The regular full algorithm
  //   will now serve as a proxy for the component dominated by combination
  //   count.
  // - Add profiler metrics to the main test ('testWorkSplittingMain')
  func testWorkSplittingMain() throws {
    var testCase = TestCase()
    testCase.taskCount = 2
    testCase.childCount = 3
    
    // Set the child latencies to random values.
    for childID in 0..<testCase.childCount {
      var latency = Float.random(in: 1...1000)
      latency.round(.toNearestOrEven)
      testCase.childLatencies[childID] = latency
    }
    
    // Generate assignments from the three algorithm variants.
    _ = runFullTest(
      testCase: testCase)
    _ = runRestrictedTest(
      testCase: testCase,
      restrictMaxCombinations: false)
    _ = runRestrictedTest(
      testCase: testCase,
      restrictMaxCombinations: true)
  }
  
  func testWorkSplittingUnit() throws {
    do {
      var test = CompleteTestCase()
      test.problemSize = (2, 8)
      test.childValues = [
        787.0,
        824.0,
        183.0,
        916.0,
        885.0,
        447.0,
        799.0,
        878.0,
      ]
      test.resultFull = 2862
      test.resultRestricted1 = 2949
      test.resultRestricted2 = 2974
      test.run()
    }
    
    do {
      var test = CompleteTestCase()
      test.problemSize = (3, 8)
      test.childValues = [
        459.0,
        713.0,
        657.0,
        672.0,
        358.0,
        345.0,
        845.0,
        202.0,
      ]
      test.resultFull = 1461
      test.resultRestricted1 = 1506
      test.resultRestricted2 = 1517
      test.run()
    }
    
    do {
      var test = CompleteTestCase()
      test.problemSize = (4, 8)
      test.childValues = [
        961.0,
        416.0,
        426.0,
        710.0,
        290.0,
        510.0,
        362.0,
        262.0,
      ]
      test.resultFull = 1040
      test.resultRestricted1 = 1072
      test.resultRestricted2 = 1104
      test.run()
    }
    
    do {
      var test = CompleteTestCase()
      test.problemSize = (5, 8)
      test.childValues = [
        770.0,
        82.0,
        267.0,
        329.0,
        770.0,
        205.0,
        720.0,
        636.0,
      ]
      test.resultFull = 801
      test.resultRestricted1 = 801
      test.resultRestricted2 = 841
      test.run()
    }
    
    do {
      var test = CompleteTestCase()
      test.problemSize = (6, 8)
      test.childValues = [
        608.0,
        866.0,
        311.0,
        289.0,
        737.0,
        921.0,
        935.0,
        596.0,
      ]
      test.resultFull = 935
      test.resultRestricted1 = 935
      test.resultRestricted2 = 935
      test.run()
    }
    
    do {
      var test = CompleteTestCase()
      test.problemSize = (7, 8)
      test.childValues = [
        453.0,
        708.0,
        120.0,
        358.0,
        15.0,
        62.0,
        45.0,
        905.0,
      ]
      test.resultFull = 905
      test.resultRestricted1 = 905
      test.resultRestricted2 = 905
      test.run()
    }
    
    do {
      var test = CompleteTestCase()
      test.problemSize = (3, 7)
      test.childValues = [
        836.0,
        860.0,
        544.0,
        411.0,
        493.0,
        765.0,
        151.0,
      ]
      test.resultFull = 1380
      test.resultRestricted1 = 1398
      test.resultRestricted2 = 1422
      test.run()
    }
    
    do {
      var test = CompleteTestCase()
      test.problemSize = (3, 6)
      test.childValues = [
        331.0,
        838.0,
        528.0,
        290.0,
        112.0,
        463.0,
      ]
      test.resultFull = 865
      test.resultRestricted1 = 906
      test.resultRestricted2 = 906
      test.run()
    }
    
    do {
      var test = CompleteTestCase()
      test.problemSize = (2, 5)
      test.childValues = [
        640.0,
        545.0,
        31.0,
        653.0,
        574.0,
      ]
      test.resultFull = 1227
      test.resultRestricted1 = 1229
      test.resultRestricted2 = 1229
      test.run()
    }
    
    do {
      var test = CompleteTestCase()
      test.problemSize = (4, 5)
      test.childValues = [
        946.0,
        642.0,
        527.0,
        464.0,
        413.0,
      ]
      test.resultFull = 946
      test.resultRestricted1 = 946
      test.resultRestricted2 = 946
      test.run()
    }
    
    do {
      var test = CompleteTestCase()
      test.problemSize = (2, 4)
      test.childValues = [
        691.0,
        380.0,
        124.0,
        382.0,
      ]
      test.resultFull = 815
      test.resultRestricted1 = 815
      test.resultRestricted2 = 815
      test.run()
    }
    
    do {
      var test = CompleteTestCase()
      test.problemSize = (2, 3)
      test.childValues = [
        475.0,
        479.0,
        60.0,
      ]
      test.resultFull = 535
      test.resultRestricted1 = 535
      test.resultRestricted2 = 535
      test.run()
    }
  }
}

// MARK: - Utilities

// All latencies and their sums will be 4 digits or less.
private func format(latency: Float) -> String {
  let rounded = latency.rounded(.toNearestOrEven)
  var repr = String(Int(rounded))
  guard repr.count > 0,
        repr.count <= 4 else {
    fatalError("Unexpected character count for latency.")
  }
  
  while repr.count < 4 {
    repr = " " + repr
  }
  return repr
}

struct TestCase {
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

private struct CompleteTestCase {
  var problemSize: (taskCount: Int, childCount: Int)?
  var childValues: [Float]?
  var resultFull: Float?
  var resultRestricted1: Float?
  var resultRestricted2: Float?
  
  func run() {
    guard let problemSize,
          let childValues,
          let resultFull,
          let resultRestricted1,
          let resultRestricted2 else {
      fatalError("Test was not fully specified.")
    }
    
    var testCase = TestCase()
    testCase.taskCount = problemSize.taskCount
    testCase.childCount = problemSize.childCount
    
    // Assign the child latencies.
    guard childValues.count == testCase.childCount else {
      fatalError("Incorrect number of child values.")
    }
    for childID in 0..<testCase.childCount {
      let latency = childValues[childID]
      testCase.childLatencies[childID] = latency
    }
    
    // Generate assignments from the three algorithm variants.
    let assignmentFull = runFullTest(
      testCase: testCase)
    let assignmentPartial1 = runRestrictedTest(
      testCase: testCase,
      restrictMaxCombinations: false)
    let assignmentPartial2 = runRestrictedTest(
      testCase: testCase,
      restrictMaxCombinations: true)
    
    // Check that restricted algorithms fill each task with ≥1 child.
    func validate(assignment: SIMD8<UInt8>) {
      let taskLatencies = testCase.taskLatencies(
        assignments: assignment)
      for taskID in 0..<testCase.taskCount {
        let latency = taskLatencies[taskID]
        guard latency > 0 else {
          fatalError("Unassigned task.")
        }
      }
    }
    validate(assignment: assignmentPartial1)
    validate(assignment: assignmentPartial2)
    
    // Check the exact value of the outputs.
    func latency(assignment: SIMD8<UInt8>) -> Float {
      let taskLatencies = testCase.taskLatencies(
        assignments: assignment)
      return taskLatencies.max()
    }
    func compare(assignment: SIMD8<UInt8>, expected: Float) {
      let latency = latency(assignment: assignment)
      XCTAssertEqual(latency, expected)
    }
    compare(assignment: assignmentFull, expected: resultFull)
    compare(assignment: assignmentPartial1, expected: resultRestricted1)
    compare(assignment: assignmentPartial2, expected: resultRestricted2)
  }
}

// MARK: - Algorithm Variants

private func runFullTest(
  testCase: TestCase
) -> SIMD8<UInt8> {
  var bestAssignment: SIMD8<UInt8>?
  var bestAssignmentLatency: Float = .greatestFiniteMagnitude
  
  // Iterate over all combinations.
  var counter: SIMD8<UInt8> = .zero
  let combinationCount = testCase.combinationCount(
    childCount: testCase.childCount)
  for combinationID in 0..<combinationCount {
    let taskLatencies = testCase.taskLatencies(assignments: counter)
    let maxTaskLatency = taskLatencies.max()
    if maxTaskLatency < bestAssignmentLatency {
      bestAssignment = counter
      bestAssignmentLatency = maxTaskLatency
    }
    
    for laneID in 0..<8 {
      counter[laneID] += 1
      if counter[laneID] >= testCase.taskCount {
        counter[laneID] = 0
      } else {
        break
      }
    }
  }
  
  guard let bestAssignment else {
    fatalError("This should never happen.")
  }
  return bestAssignment
}

private struct PreparationStage {
  var sortedChildPairs: [SIMD2<Float>]
  var fixedChildAssignments: SIMD8<UInt8>
  
  init(testCase: TestCase) {
    sortedChildPairs = createChildPairs(testCase: testCase)
    sortedChildPairs.sort {
      $0[1] < $1[1]
    }
    fixedChildAssignments = createFixedAssignments(
      pairs: sortedChildPairs)
    sortedChildPairs.removeLast(createFixedChildCount())
  }
  
  static func createChildPairs(
    testCase: TestCase
  ) -> [SIMD2<Float>] {
    var output: [SIMD2<Float>] = []
    for childID in 0..<testCase.childCount {
      let latency = testCase.childLatencies[childID]
      let pair = SIMD2(
        Float(childID),
        latency)
      output.append(pair)
    }
    return output
  }
  
  static func maxNativeCombinations(
    testCase: TestCase
  ) -> Int {
    restrictMaxCombinations ? 20 : 100
  }
  
  static func createFixedChildCount() -> Int {
    if testCase.nativeCombinationCount < maxNativeCombinations() {
      return testCase.taskCount + 1
    } else {
      return testCase.taskCount + 2
    }
  }
  
  static func createFixedAssignments(
    pairs: [SIMD2<Float>]
  ) -> SIMD8<UInt8> {
    var output = SIMD8<UInt8>(repeating: .max)
    guard testCase.childCount > testCase.taskCount else {
      fatalError("Invalid conditions for the restricted algorithm.")
    }
    
    // Assign the highest-index children to the lowest-index tasks.
    for taskID in 0..<testCase.taskCount {
      let sortedChildID = testCase.childCount - 1 - taskID
      let pair = pairs[sortedChildID]
      let childID = Int(pair[0])
      output[childID] = UInt8(taskID)
    }
    
    // Assign the largest of the remaining children to the highest-index task.
    let remainingChildCount = testCase.childCount - testCase.taskCount
    if testCase.nativeCombinationCount < maxNativeCombinations() {
      let pair = pairs[remainingChildCount - 1]
      let childID = Int(pair[0])
      let taskID = testCase.taskCount - 1
      output[childID] = UInt8(taskID)
    } else {
      let pair0 = pairs[remainingChildCount - 2]
      let pair1 = pairs[remainingChildCount - 1]
      let childID0 = Int(pair0[0])
      let childID1 = Int(pair1[0])
      output[childID0] = UInt8(testCase.taskCount - 2)
      output[childID1] = UInt8(testCase.taskCount - 1)
    }
    
    return output
  }
}

private func runRestrictedTest(
  testCase: TestCase,
  restrictMaxCombinations: Bool
) -> SIMD8<UInt8> {
  let preparationStage = PreparationStage()
  
  // Declare the state variables for the best assignment.
  var bestAssignment: SIMD8<UInt8>?
  var bestAssignmentLatency: Float = .greatestFiniteMagnitude
  
  // Iterate over all combinations of variable children.
  var counter: SIMD8<UInt8> = .zero
  let combinationCount = testCase.combinationCount(
    childCount: preparationStage.sortedChildPairs.count)
  for combinationID in 0..<combinationCount {
    // Merge the fixed and variable assignments.
    var combinedAssignments = preparationStage.fixedChildAssignments
    for sortedChildID in preparationStage.sortedChildPairs.indices {
      let pair = preparationStage.sortedChildPairs[sortedChildID]
      let childID = Int(pair[0])
      let taskID = counter[sortedChildID]
      combinedAssignments[childID] = UInt8(taskID)
    }
    
    let taskLatencies = testCase.taskLatencies(
      assignments: combinedAssignments)
    let maxTaskLatency = taskLatencies.max()
    if maxTaskLatency < bestAssignmentLatency {
      bestAssignment = combinedAssignments
      bestAssignmentLatency = maxTaskLatency
    }
    
    for laneID in 0..<8 {
      counter[laneID] += 1
      if counter[laneID] >= testCase.taskCount {
        counter[laneID] = 0
      } else {
        break
      }
    }
  }
  
  guard let bestAssignment else {
    fatalError("This should never happen.")
  }
  return bestAssignment
}
