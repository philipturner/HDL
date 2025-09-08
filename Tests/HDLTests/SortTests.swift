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
  
#if RELEASE
  func testSortPerformance() throws {
    // Revert to 10 and true after any refactorings
    // Default for this benchmarking period is 20 and false
    let latticeScale: Float = 40
    let testParallel = Bool.random() ? false : false
    let lattice = Lattice<Hexagonal> { h, k, l in
      let h2k = h + 2 * k
      Bounds { latticeScale * (2 * h + h2k + l) }
      Material { .elemental(.silicon) }
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
    for trialID in 2...2 {
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
#endif
  
  func testWorkSplittingUnit() throws {
    do {
      var test = TestCase()
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
      test.result = 2974
      test.run()
    }
    
    do {
      var test = TestCase()
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
      test.result = 1517
      test.run()
    }
    
    do {
      var test = TestCase()
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
      test.result = 1104
      test.run()
    }
    
    do {
      var test = TestCase()
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
      test.result = 841
      test.run()
    }
    
    do {
      var test = TestCase()
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
      test.result = 935
      test.run()
    }
    
    do {
      var test = TestCase()
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
      test.result = 905
      test.run()
    }
    
    do {
      var test = TestCase()
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
      test.result = 1422
      test.run()
    }
    
    do {
      var test = TestCase()
      test.problemSize = (3, 6)
      test.childValues = [
        331.0,
        838.0,
        528.0,
        290.0,
        112.0,
        463.0,
      ]
      test.result = 906
      test.run()
    }
    
    do {
      var test = TestCase()
      test.problemSize = (2, 5)
      test.childValues = [
        640.0,
        545.0,
        31.0,
        653.0,
        574.0,
      ]
      test.result = 1229
      test.run()
    }
    
    do {
      var test = TestCase()
      test.problemSize = (4, 5)
      test.childValues = [
        946.0,
        642.0,
        527.0,
        464.0,
        413.0,
      ]
      test.result = 946
      test.run()
    }
    
    do {
      var test = TestCase()
      test.problemSize = (2, 4)
      test.childValues = [
        691.0,
        380.0,
        124.0,
        382.0,
      ]
      test.result = 815
      test.run()
    }
    
    do {
      var test = TestCase()
      test.problemSize = (2, 3)
      test.childValues = [
        475.0,
        479.0,
        60.0,
      ]
      test.result = 535
      test.run()
    }
  }
}

private struct TestCase {
  var problemSize: (taskCount: Int, childCount: Int)?
  var childValues: [Float]?
  var result: Float?
  
  func run() {
    guard let problemSize,
          let childValues,
          let result else {
      fatalError("Test was not fully specified.")
    }
    
    var test = WorkSplittingTest()
    test.taskCount = problemSize.taskCount
    test.childCount = problemSize.childCount
    
    // Assign the child latencies.
    guard childValues.count == test.childCount else {
      fatalError("Incorrect number of child values.")
    }
    for childID in 0..<test.childCount {
      let latency = childValues[childID]
      test.childLatencies[childID] = latency
    }
    
    // Generate the output.
    let assignment = test.run()
    
    // Check that each task is filled with ≥1 child.
    func validate(assignment: SIMD8<UInt8>) {
      let taskLatencies = test.taskLatencies(
        assignments: assignment)
      for taskID in 0..<test.taskCount {
        let latency = taskLatencies[taskID]
        guard latency > 0 else {
          fatalError("Unassigned task.")
        }
      }
    }
    validate(assignment: assignment)
    
    // Check the exact value of the output.
    func latency(assignment: SIMD8<UInt8>) -> Float {
      let taskLatencies = test.taskLatencies(
        assignments: assignment)
      return taskLatencies.max()
    }
    func compare(assignment: SIMD8<UInt8>, expected: Float) {
      let latency = latency(assignment: assignment)
      XCTAssertEqual(latency, expected)
    }
    compare(assignment: assignment, expected: result)
  }
}
