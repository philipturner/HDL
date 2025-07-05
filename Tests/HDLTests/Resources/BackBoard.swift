//
//  BackBoard.swift
//  HDLTests
//
//  Created by Philip Turner on 1/2/24.
//

import HDL
import XCTest

// A protocol that handles the compilation passes of all 3 back board pieces.
protocol BackBoardComponent {
  var topology: Topology { get set }
  
  // The expected atom and bond count at each stage of compilation.
  static var expectedTopologyState: [Int: SIMD2<Int>] { get }
  
  mutating func compilationPass0()
}

extension BackBoardComponent {
  // 'optimized' is whether to report out-of-the-box performance. One section
  // of the code is maximum possible performance available from the public API.
  // Another is performance of typical code segments.
  mutating func compile(reportingPerformance: Bool) {
    let checkpoint0 = Profiler.time()
    
    compilationPass0()
    XCTAssertEqual(topology.atoms.count, Self.expectedTopologyState[0]![0])
    XCTAssertEqual(topology.bonds.count, Self.expectedTopologyState[0]![1])
    
    let checkpoint1 = Profiler.time()
    
    compilationPass1()
    XCTAssertEqual(topology.atoms.count, Self.expectedTopologyState[1]![0])
    XCTAssertEqual(topology.bonds.count, Self.expectedTopologyState[1]![1])
    
    let checkpoint2 = Profiler.time()
    
    compilationPass2()
    XCTAssertEqual(topology.atoms.count, Self.expectedTopologyState[2]![0])
    XCTAssertEqual(topology.bonds.count, Self.expectedTopologyState[2]![1])
    
    let checkpoint3 = Profiler.time()
    
    guard reportingPerformance else {
      return
    }
    
    let checkpoints = [checkpoint0, checkpoint1, checkpoint2, checkpoint3]
    var elapsedTimes: [Double] = []
    var microseconds: [Int] = []
    for i in 0..<checkpoints.count - 1 {
      elapsedTimes.append(checkpoints[i + 1] - checkpoints[i])
      
      let thisMicroseconds = Int(elapsedTimes[i] * 1e6)
      microseconds.append(thisMicroseconds)
    }
    
    func milliseconds(_ microseconds: Int) -> String {
      let value = Double(microseconds) / 1000
      var repr = String(format: "%.1f", value)
      while repr.count < 4 {
        repr = " " + repr
      }
      return repr + " ms"
    }
    
    let atomCount = topology.atoms.count
    print()
    print("-    atoms: \(atomCount / 1000),\(atomCount % 1000)")
    print("-  lattice: \(milliseconds(microseconds[0]))")
    print("-    match: \(milliseconds(microseconds[1]))")
    print("- orbitals: \(milliseconds(microseconds[2]))")
    print("-    total: \(milliseconds(microseconds.reduce(0, +)))")
  }
  
  mutating func compilationPass1() {
    let matches = topology.match(topology.atoms)
    
    var removedAtoms: [UInt32] = []
    var insertedBonds: [SIMD2<UInt32>] = []
    for i in topology.atoms.indices {
      let matchRange = matches[i]
      if matchRange.count <= 2 {
        removedAtoms.append(UInt32(truncatingIfNeeded: i))
      } else if matchRange.count <= 5 {
        for j in matchRange where i < j {
          insertedBonds.append(SIMD2(UInt32(truncatingIfNeeded: i),
                                     UInt32(truncatingIfNeeded: j)))
        }
      } else {
        fatalError("More than 5 matches.")
      }
    }
    
    topology.bonds += insertedBonds
    topology.remove(atoms: removedAtoms)
  }
  
  mutating func compilationPass2() {
    let orbitalLists = topology.nonbondingOrbitals()
    let chBondLength = Element.carbon.covalentRadius +
    Element.hydrogen.covalentRadius
    
    var insertedAtoms: [Atom] = []
    var insertedBonds: [SIMD2<UInt32>] = []
    
    for i in topology.atoms.indices {
      let atom = topology.atoms[i]
      let orbitalList = orbitalLists[i]
      for orbital in orbitalList {
        let position = atom.position + orbital * chBondLength
        let hydrogen = Atom(position: position, element: .hydrogen)
        let hydrogenID = topology.atoms.count + insertedAtoms.count
        insertedAtoms.append(hydrogen)
        
        let bond = SIMD2(UInt32(i), UInt32(hydrogenID))
        insertedBonds.append(bond)
      }
    }
    topology.atoms += insertedAtoms
    topology.bonds += insertedBonds
  }
}

struct BackBoardSmallLeft: BackBoardComponent {
  var topology: Topology = .init()
  
  static let expectedTopologyState: [Int: SIMD2<Int>] = [
    0: [23492, 0],
    1: [23488, 40643],
    2: [36154, 53309],
  ]
  
  mutating func compilationPass0() {
    var archSpacing: Float = 18.2
    archSpacing /= Constant(.hexagon) { .elemental(.carbon) }
    archSpacing.round(.toNearestOrEven)
    
    let lattice = Lattice<Hexagonal> { h, k, l in
      let h2k = h + 2 * k
      Bounds { 120 * h + 100 * h2k + 1 * l }
      Material { .elemental(.carbon) }
      
      Volume {
        Origin { 9 * h }
        Origin { archSpacing/2 * h }
        
        Concave {
          Convex {
            Origin { -15 * h }
            Plane { h }
          }
          Convex {
            Origin { 15 * h }
            Plane { -h }
          }
          Origin { 35 * h2k }
          for direction in [k, h + k] {
            Convex {
              Origin { 12 * direction }
              Plane { -direction }
            }
          }
        }
        
        Convex {
          Origin { -0.75 * h }
          Plane { h }
        }
        Convex {
          Origin { (-archSpacing/2 + 0.5) * h }
          Plane { -h }
        }
        
        Replace { .empty }
      }
    }
    topology.atoms += lattice.atoms
  }
}

struct BackBoardSmallRight: BackBoardComponent {
  var topology: Topology = .init()
  
  static let expectedTopologyState: [Int: SIMD2<Int>] = [
    0: [13638, 0],
    1: [13634, 23423],
    2: [21324, 31113],
  ]
  
  mutating func compilationPass0() {
    var archSpacing: Float = 18.2
    archSpacing /= Constant(.hexagon) { .elemental(.carbon) }
    archSpacing.round(.toNearestOrEven)
    
    let lattice = Lattice<Hexagonal> { h, k, l in
      let h2k = h + 2 * k
      Bounds { 357 * h + 100 * h2k + 1 * l }
      Material { .elemental(.carbon) }
      
      func createArch(index: Int) {
        Origin { 9 * h }
        Origin { Float(index) * archSpacing * h }
        Origin { archSpacing/2 * h }
        
        Concave {
          Convex {
            Origin { -15 * h }
            Plane { h }
          }
          Convex {
            Origin { 15 * h }
            Plane { -h }
          }
          Origin { 35 * h2k }
          for direction in [k, h + k] {
            Convex {
              Origin { 12 * direction }
              Plane { -direction }
            }
          }
        }
        
        Convex {
          Origin { 0.75 * h }
          Plane { -h }
        }
      }
      
      Volume {
        Convex {
          createArch(index: 4)
        }
        
        for elevation in [Float(5), 59, 96] {
          Convex {
            Origin { 352 * h }
            Origin { elevation * h2k }
            Concave {
              Plane { h + k + h }
              Plane { -k + h }
            }
          }
        }
        
        Replace { .empty }
      }
    }
    topology.atoms += lattice.atoms
  }
}

struct BackBoardLarge: BackBoardComponent {
  var topology: Topology = .init()
  
  static let expectedTopologyState: [Int: SIMD2<Int>] = [
    0: [237706, 0],
    1: [237704, 414085],
    2: [360350, 536731],
  ]
  
  mutating func compilationPass0() {
    var archSpacing: Float = 18.2
    archSpacing /= Constant(.hexagon) { .elemental(.carbon) }
    archSpacing.round(.toNearestOrEven)
    
    let lattice = Lattice<Hexagonal> { h, k, l in
      let h2k = h + 2 * k
      Bounds { 357 * h + 100 * h2k + 1 * l }
      Material { .elemental(.carbon) }
      
      func createArch(index: Int) {
        Origin { 9 * h }
        Origin { Float(index) * archSpacing * h }
        Origin { archSpacing/2 * h }
        
        Concave {
          Convex {
            Origin { -15 * h }
            Plane { h }
          }
          Convex {
            Origin { 15 * h }
            Plane { -h }
          }
          Origin { 35 * h2k }
          for direction in [k, h + k] {
            Convex {
              Origin { 12 * direction }
              Plane { -direction }
            }
          }
        }
      }
      
      Volume {
        for index in 0...4 {
          Convex {
            createArch(index: index)
          }
        }
        
        for elevation in [Float(5), 59, 96] {
          Convex {
            Origin { 352 * h }
            Origin { elevation * h2k }
            Concave {
              Plane { h + k + h }
              Plane { -k + h }
            }
          }
        }
        
        Replace { .empty }
      }
    }
    topology.atoms += lattice.atoms
  }
}
