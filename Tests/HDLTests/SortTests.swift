@testable import HDL
import XCTest

final class SortTests: XCTestCase {
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
    let latticeScale: Float = 10
    let lattice = Lattice<Hexagonal> { h, k, l in
      let h2k = h + 2 * k
      Bounds { latticeScale * (2 * h + h2k + l) }
      Material { .elemental(.carbon) }
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
    
    for trialID in 0..<4 {
      let trial = Trial(lattice: lattice, index: trialID)
      
      var resultGrid: [UInt32]
      do {
        var topology = Topology()
        topology.atoms = trial.atoms
        resultGrid = topology.sort()
      }
      
      var resultOctree: [UInt32]
      do {
        let sorter = OctreeSorter(atoms: trial.atoms)
        let reordering = sorter.mortonReordering()
        resultOctree = OctreeSorter.invertOrder(reordering)
      }
      XCTAssertEqual(resultGrid, resultOctree)
    }
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
