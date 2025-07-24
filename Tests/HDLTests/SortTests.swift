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
  //
  // find the optimal grouping of children in tasks, which minimizes the
  // length of the longest task
  // - create an algorithm to explicitly iterate over all possible combinations
  // - study the characteristics with computationally tractable parameters
  // - study the effects of restricting the combinatorial space
  // - study the worst-case execution time of the refined algorithm
  func testWorkSplitting() throws {
    var testCase = TestCase()
    testCase.taskCount = 3
    testCase.childCount = 5
    
    // Set the child latencies to random values.
    for childID in 0..<testCase.childCount {
      let latency = Float.random(in: 0..<1000)
      testCase.childLatencies[childID] = latency
    }
    
    runFullTest(testCase: testCase)
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
  
  // combinations = tasks^children
  var combinationCount: Int {
    var output: Int = 1
    for _ in 0..<childCount {
      output = output * taskCount
    }
    return output
  }
  
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
  
  // Input: which task each child is assigned to.
  func combinationRepr(
    assignments: SIMD8<UInt8>
  ) -> String {
    var tasks = [[Float]](repeating: [], count: taskCount)
    for childID in 0..<childCount {
      let latency = childLatencies[childID]
      let taskID = assignments[childID]
      tasks[Int(taskID)].append(latency)
    }
    
    var maxChildCount: Int = .zero
    for task in tasks {
      let childCount = task.count
      if childCount > maxChildCount {
        maxChildCount = childCount
      }
    }
    
    var outputLines: [String] = []
    for lineID in 0..<maxChildCount {
      var lineEntries: [String] = []
      for taskID in 0..<taskCount {
        var entry: String
        if lineID < tasks[taskID].count {
          let latency = tasks[taskID][lineID]
          entry = format(latency: latency)
        } else {
          entry = "    "
        }
        lineEntries.append(entry)
      }
      
      let line = lineEntries.joined(separator: "  ")
      outputLines.append(line)
    }
    
    let output = outputLines.joined(separator: "\n")
    return output
  }
}

// A line of sorted combinations to display.
private struct CombinationLine {
  var entries: [SIMD2<Float>] = []
}

private func createCombinationLines(
  pairs: [SIMD2<Float>]
) -> [CombinationLine] {
  let maxEntriesPerLine: Int = 15
  var output: [CombinationLine] = []
  
  var currentLine: [SIMD2<Float>] = []
  for pair in pairs {
    currentLine.append(pair)
    if currentLine.count >= maxEntriesPerLine {
      let combinationLine = CombinationLine(entries: currentLine)
      output.append(combinationLine)
      currentLine = []
    }
  }
  
  if currentLine.count > 0 {
    let combinationLine = CombinationLine(entries: currentLine)
    output.append(combinationLine)
    currentLine = []
  }
  return output
}

private func display(combinationLines: [CombinationLine]) {
  for line in combinationLines {
    var entriesTaskID: [String] = []
    var entriesLatency: [String] = []
    for entry in line.entries {
      let taskID = entry[0]
      let reprTaskID = format(latency: taskID)
      entriesTaskID.append(reprTaskID)
      
      let latency = entry[1]
      let reprLatency = format(latency: latency)
      entriesLatency.append(reprLatency)
    }
    
    let lineTaskID = entriesTaskID.joined(separator: "  ")
    let lineLatency = entriesLatency.joined(separator: "  ")
    print(lineTaskID)
    print(lineLatency)
    print()
  }
}

// MARK: - Algorithm Variants

private func runFullTest(testCase: TestCase) {
  // Declare variables for finding the best combination.
  var combinationPairs: [SIMD2<Float>] = []
  
  // Iterate over all combinations.
  print()
  var counter: SIMD8<UInt8> = .zero
  for combinationID in 0..<testCase.combinationCount {
    print("#\(combinationID)")
    print(counter)
    let repr = testCase.combinationRepr(assignments: counter)
    print(repr)
    
    let taskLatencies = testCase.taskLatencies(assignments: counter)
    for taskID in 0..<testCase.taskCount {
      let latency = taskLatencies[taskID]
      let repr = format(latency: latency)
      print(repr, terminator: "  ")
    }
    print()
    print()
    
    let maxTaskLatency = taskLatencies.max()
    let pair = SIMD2(
      Float(combinationID),
      maxTaskLatency)
    combinationPairs.append(pair)
    
    for laneID in 0..<8 {
      counter[laneID] += 1
      if counter[laneID] >= testCase.taskCount {
        counter[laneID] = 0
      } else {
        break
      }
    }
  }
  combinationPairs.sort {
    $0[1] < $1[1]
  }
  
  let combinationLines = createCombinationLines(
    pairs: combinationPairs)
  display(combinationLines: combinationLines)
}
