import XCTest
import HDL
import Numerics
import SystemPackage

private let startTime = ContinuousClock.now

func cross_platform_media_time() -> Double {
  let duration = ContinuousClock.now.duration(to: startTime)
  let seconds = duration.components.seconds
  let attoseconds = duration.components.attoseconds
  return -(Double(seconds) + Double(attoseconds) * 1e-18)
}

private func fmt(_ start: Double, _ end: Double) -> String {
  let seconds = end - start
  if seconds > 1 {
    return String(format: "%.3f", seconds) + " s"
  } else if seconds > 1e-3 {
    return String(format: "%.3f", seconds * 1e3) + " ms"
  } else {
    return String(format: "%.3f", seconds * 1e6) + " μs"
  }
}

final class PerformanceTests: XCTestCase {
  static let printPerformanceSummary = false
  
#if RELEASE
  func testGoldSurface() throws {
    let overallStart = cross_platform_media_time()
    var gridInitStart: Double = 0
    var gridInitEnd: Double = 0
    var intersectStart: Double = 0
    var intersectEnd: Double = 0
    var replaceStart: Double = 0
    var replaceEnd: Double = 0
    
    let scaleFactor: Float = 2
    let lattice = Lattice<Cubic> { h, k, l in
      gridInitStart = cross_platform_media_time()
      Bounds { scaleFactor * 40 * (h + k + l) }
      Material { .elemental(.gold) }
      
      Volume {
        gridInitEnd = cross_platform_media_time()
        intersectStart = cross_platform_media_time()
        
        Convex {
          Origin { scaleFactor * 20 * (h + k + l) }
          
          Convex {
            Origin { 0.5 * (h + k + l) }
            Plane { h + k + l }
          }
          Convex {
            Origin { -0.5 * (h + k + l) }
            Plane { -(h + k + l) }
          }
          
          // Change the chiseling on the 3 (110) sides.
          for direction in [h + k, h + l, k + l] {
            Convex {
              Origin { scaleFactor * 10 * direction }
              Plane { direction }
            }
          }
        }
        
        // Change the chiseling on the 3 (100) sides.
        for direction in [h - k - l, k - h - l, l - h - k] {
          Convex {
            Origin { scaleFactor * 6.75 * direction }
            Plane { direction }
          }
        }
        
        intersectEnd = cross_platform_media_time()
        replaceStart = cross_platform_media_time()
        Replace { .empty }
        replaceEnd = cross_platform_media_time()
      }
    }
    XCTAssertEqual(lattice.atoms.count, 57601)
    
    let overallEnd = cross_platform_media_time()
    if Self.printPerformanceSummary {
      print("gold surface:")
      print("- overall:   \(fmt(overallStart, overallEnd))")
      print("- grid init: \(fmt(gridInitStart, gridInitEnd))")
      print("- intersect: \(fmt(intersectStart, intersectEnd))")
      print("- replace:   \(fmt(replaceStart, replaceEnd))")
    }
    
    // Before optimizations: ~0.318 seconds
    // After optimization 1: ~0.066 seconds
    // After optimization 2: ~0.059 seconds
    // After optimization 3: ~0.044 seconds
    // After optimization 5: ~0.017 seconds
    // After optimization 6: ~0.009 seconds
    // 35.3x speedup
  }
  
  func testSiliconProbe() throws {
    let overallStart = cross_platform_media_time()
    var gridInitStart: Double = 0
    var gridInitEnd: Double = 0
    var intersectStart: Double = 0
    var intersectEnd: Double = 0
    var replaceStart: Double = 0
    var replaceEnd: Double = 0
    
    let lattice = Lattice<Cubic> { h, k, l in
      gridInitStart = cross_platform_media_time()
      Bounds { 50 * (h + k + l) }
      Material { .elemental(.silicon) }
      
      Volume {
        gridInitEnd = cross_platform_media_time()
        intersectStart = cross_platform_media_time()
        
        var directions: [SIMD3<Float>] = []
        directions.append([1, 1, 0])
        directions.append([1, 0, 1])
        directions.append([0, 1, 1])
        directions.append([1, 0, 0])
        directions.append([0, 1, 0])
        directions.append([0, 0, 1])
        let direction111 = -SIMD3<Float>(1, 1, 1) / Float(3).squareRoot()
        let rotation = Quaternion<Float>(angle: .pi / 6, axis: direction111)
        
        for directionID in directions.indices {
          var adjusted = directions[directionID]
          var adjustedLength = (adjusted * adjusted).sum().squareRoot()
          adjusted /= adjustedLength
          
          let dotProduct = (adjusted * direction111).sum()
          adjusted -= direction111 * dotProduct
          adjustedLength = (adjusted * adjusted).sum().squareRoot()
          adjusted /= adjustedLength
          adjusted = rotation.act(on: adjusted)
          
          directions[directionID] = adjusted
        }
        
        var passes: [SIMD2<Float>] = []
        passes.append([2.50, 6.50])
        passes.append([4.50, 6.83])
        passes.append([5.50, 7.17])
        passes.append([6.00, 7.50])
        passes.append([6.50, 8.17])
        passes.append([7.00, 8.83])
        passes.append([7.50, 10.50])
        
        Concave {
          for pass in passes {
            let distanceWidth: Float = pass[0]
            let distanceHeight: Float = pass[1]
            
            Convex {
              Convex {
                Origin { distanceHeight * (h + k + l) }
                Plane { -(h + k + l) }
              }
              Convex {
                Origin { 37 * (h + k + l) }
                Plane { h + k + l }
              }
              for direction in directions {
                Convex {
                  Origin { (distanceWidth + 0.01) * direction }
                  Plane { direction }
                }
              }
            }
          }
        }
        
        intersectEnd = cross_platform_media_time()
        replaceStart = cross_platform_media_time()
        Replace { .empty }
        replaceEnd = cross_platform_media_time()
      }
    }
    XCTAssertEqual(lattice.atoms.count, 81142)
    
    let overallEnd = cross_platform_media_time()
    if Self.printPerformanceSummary {
      print("silicon probe:")
      print("- overall:   \(fmt(overallStart, overallEnd))")
      print("- grid init: \(fmt(gridInitStart, gridInitEnd))")
      print("- intersect: \(fmt(intersectStart, intersectEnd))")
      print("- replace:   \(fmt(replaceStart, replaceEnd))")
    }
    
    // Before optimizations: ~0.300 seconds
    // After optimization 1: ~0.047 seconds
    // After optimization 3: ~0.031 seconds
    // After optimization 5: ~0.024 seconds
    // After optimization 6: ~0.012 seconds
    // 25.0x speedup
  }
  
  func testGoldSurface2() throws {
    let overallStart = cross_platform_media_time()
    var gridInitStart: Double = 0
    var gridInitEnd: Double = 0
    var intersectStart: Double = 0
    var intersectEnd: Double = 0
    var replaceStart: Double = 0
    var replaceEnd: Double = 0
    
    let scaleFactor: Float = 6
    let lattice = Lattice<Cubic> { h, k, l in
      gridInitStart = cross_platform_media_time()
      Bounds { scaleFactor * 40 * (h + k + l) }
      Material { .elemental(.gold) }
      
      Volume {
        gridInitEnd = cross_platform_media_time()
        intersectStart = cross_platform_media_time()
        
        Convex {
          Origin { scaleFactor * 20 * (h + k + l) }
          
          Convex {
            Origin { 0.5 * (h + k + l) }
            Plane { h + k + l }
          }
          Convex {
            Origin { 0.25 * (h + k + l) }
            Plane { -(h + k + l) }
          }
        }
        
        intersectEnd = cross_platform_media_time()
        replaceStart = cross_platform_media_time()
        Replace { .empty }
        replaceEnd = cross_platform_media_time()
      }
    }
    XCTAssertEqual(lattice.atoms.count, 173517)
    
    let overallEnd = cross_platform_media_time()
    if Self.printPerformanceSummary {
      print("gold surface 2:")
      print("- overall:   \(fmt(overallStart, overallEnd))")
      print("- grid init: \(fmt(gridInitStart, gridInitEnd))")
      print("- intersect: \(fmt(intersectStart, intersectEnd))")
      print("- replace:   \(fmt(replaceStart, replaceEnd))")
    }
    
    // Before optimizations: 4.474 seconds
    // After optimization 3: 0.862 seconds
    // After optimization 4: 0.822 seconds
    // After optimization 5: 0.148 seconds
    // After optimization 6: 0.093 seconds
    // 48.1x speedup
  }
  
  func testSiliconProbe2() throws {
    let overallStart = cross_platform_media_time()
    var gridInitStart: Double = 0
    var gridInitEnd: Double = 0
    var intersectStart: Double = 0
    var intersectEnd: Double = 0
    var replaceStart: Double = 0
    var replaceEnd: Double = 0
    var intersectStart2: Double = 0
    var intersectEnd2: Double = 0
    var replaceStart2: Double = 0
    var replaceEnd2: Double = 0
    
    let lattice = Lattice<Cubic> { h, k, l in
      gridInitStart = cross_platform_media_time()
      Bounds { 80 * (h + k + l) }
      Material { .elemental(.silicon) }
      let topCutoff: Float = 67
      
      Volume {
        gridInitEnd = cross_platform_media_time()
        intersectStart = cross_platform_media_time()
        
        var directions: [SIMD3<Float>] = []
        directions.append([1, 1, 0])
        directions.append([1, 0, 1])
        directions.append([0, 1, 1])
        directions.append([1, 0, 0])
        directions.append([0, 1, 0])
        directions.append([0, 0, 1])
        let direction111 = -SIMD3<Float>(1, 1, 1) / Float(3).squareRoot()
        let rotation = Quaternion<Float>(angle: .pi / 6, axis: direction111)
        
        for directionID in directions.indices {
          var adjusted = directions[directionID]
          var adjustedLength = (adjusted * adjusted).sum().squareRoot()
          adjusted /= adjustedLength
          
          let dotProduct = (adjusted * direction111).sum()
          adjusted -= direction111 * dotProduct
          adjustedLength = (adjusted * adjusted).sum().squareRoot()
          adjusted /= adjustedLength
          adjusted = rotation.act(on: adjusted)
          
          directions[directionID] = adjusted
        }
        
        var passes: [SIMD4<Float>] = []
        passes.append([2.50, 6.50, 1.0, 0.5])
        passes.append([4.50, 6.83, 1.0, 0.5])
        passes.append([5.50, 7.17, 1.0, 0.5])
        passes.append([6.00, 7.50, 1.0, 0.5])
        passes.append([6.50, 8.17, 1.0, 0.5])
        passes.append([7.00, 8.83, 1.0, 0.5])
        passes.append([7.50, 9.50, 1.0, 0.5])
        passes.append([8.00, 10.83, 1.0, 0.5])
        
        passes.append([8.25, 12.83, 1.0, 0.5])
        passes.append([8.50, 14.83, 1.0, 0.5])
        passes.append([9.00, 17.17, 1.0, 0.5])
        passes.append([9.50, 20.17, 1.0, 0.5])
        passes.append([9.75, 23.17, 1.0, 0.5])
        passes.append([10.00, 27.17, 0.5, 0.25])
        passes.append([10.25, 31.17, 0.5, 0.25])
        passes.append([10.75, 36.17, 0.75, 0.5])
        passes.append([11.00, 41.17, 0.5, 0.25])
        passes.append([11.50, 47.17, 0.75, 0.5])
        passes.append([11.75, 53.17, 0.5, 0.25])
        
        Concave {
          for pass in passes {
            let distanceWidth: Float = pass[0]
            let distanceHeight: Float = pass[1]
            
            Convex {
              Convex {
                Origin { distanceHeight * (h + k + l) }
                Plane { -(h + k + l) }
              }
              Convex {
                Origin { topCutoff * (h + k + l) }
                Plane { h + k + l }
              }
              for direction in directions {
                Convex {
                  Origin { (distanceWidth + 0.01) * direction }
                  Plane { direction }
                }
              }
            }
          }
        }
        
        Convex {
          for pass in passes {
            let distanceWidth: Float = pass[0]
            let distanceHeight: Float = pass[1]
            let wallThickness: Float = pass[2]
            
            Concave {
              Convex {
                Origin { (distanceHeight + wallThickness) * (h + k + l) }
                Plane { (h + k + l) }
              }
              Convex {
                Origin { 12 * (h + k + l) }
                Plane { h + k + l }
              }
              Convex {
                Origin { (topCutoff - 1) * (h + k + l) }
                Plane { -(h + k + l) }
              }
              for direction in directions {
                Convex {
                  Origin {
                    (distanceWidth + 0.01 - wallThickness) * direction
                  }
                  Plane { -direction }
                }
              }
            }
          }
        }
        
        intersectEnd = cross_platform_media_time()
        replaceStart = cross_platform_media_time()
        Replace { .empty }
        replaceEnd = cross_platform_media_time()
        intersectStart2 = cross_platform_media_time()
        
        Convex {
          for pass in passes {
            let distanceWidth: Float = pass[0]
            let distanceHeight: Float = pass[1]
            let wallThickness: Float = pass[2]
            let carbonDelta: Float = pass[3]
            
            Concave {
              Convex {
                Origin { (distanceHeight + wallThickness - carbonDelta) * (h + k + l) }
                Plane { (h + k + l) }
              }
              Convex {
                Origin { (12 - carbonDelta) * (h + k + l) }
                Plane { h + k + l }
              }
              Convex {
                Origin { (topCutoff - carbonDelta) * (h + k + l) }
                Plane { -(h + k + l) }
              }
              for direction in directions {
                Convex {
                  Origin {
                    (distanceWidth + 0.00 - wallThickness + carbonDelta) * direction
                  }
                  Plane { -direction }
                }
              }
            }
          }
        }
        
        intersectEnd2 = cross_platform_media_time()
        replaceStart2 = cross_platform_media_time()
        Replace { .atom(.carbon) }
        replaceEnd2 = cross_platform_media_time()
      }
    }
    XCTAssertEqual(lattice.atoms.count, 64226)
    
    let overallEnd = cross_platform_media_time()
    intersectEnd += intersectEnd2 - intersectStart2
    replaceEnd += replaceEnd2 - replaceStart2
    if Self.printPerformanceSummary {
      print("silicon probe 2:")
      print("- overall:   \(fmt(overallStart, overallEnd))")
      print("- grid init: \(fmt(gridInitStart, gridInitEnd))")
      print("- intersect: \(fmt(intersectStart, intersectEnd))")
      print("- replace:   \(fmt(replaceStart, replaceEnd))")
    }
    
    // Before optimizations: ~9.121 seconds
    // After optimization 3: ~0.525 seconds
    // After optimization 4: ~0.452 seconds
    // After optimization 5: ~0.426 seconds
    // After optimization 6: ~0.160 seconds
    // 57.0x speedup
  }
  
  func testSort() throws {
    let latticeScale: Float = 10
    let testParallel = Bool.random() ? true : true
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
    
    for trialID in 0..<4 {
      var trialAtoms: [Atom]
      var trialName: String
      
      switch trialID {
      case 0:
        var topology = Topology()
        topology.insert(atoms: lattice.atoms)
        topology.sort()
        
        trialAtoms = topology.atoms
        trialName = "pre-sorted"
      case 1:
        trialAtoms = lattice.atoms
        trialName = "lattice   "
      case 2:
        trialAtoms = lattice.atoms.shuffled()
        trialName = "shuffled  "
      case 3:
        trialAtoms = lattice.atoms.reversed()
        trialName = "reversed  "
      default:
        fatalError("This should never happen.")
      }
      
      let startParallel = cross_platform_media_time()
      var resultGrid1: [UInt32] = []
      var resultGrid2: [UInt32] = []
      if testParallel {
        DispatchQueue.concurrentPerform(iterations: 2) { z in
          var topology = Topology()
          topology.insert(atoms: trialAtoms)
          let resultGrid = topology.sort()
          if z == 0 {
            resultGrid1 = resultGrid
          } else {
            resultGrid2 = resultGrid
          }
        }
      }
      let endParallel = cross_platform_media_time()
      
      let startGrid = cross_platform_media_time()
      var topology = Topology()
      topology.insert(atoms: trialAtoms)
      let resultGrid = topology.sort()
      if testParallel {
        var topology = Topology()
        topology.insert(atoms: trialAtoms)
        _ = topology.sort()
      }
      let endGrid = cross_platform_media_time()
      
      let startOctree = cross_platform_media_time()
      let octree = OctreeSorter(atoms: trialAtoms)
      let resultOctree = octree.mortonReordering()
      if testParallel {
        let octree = OctreeSorter(atoms: trialAtoms)
        _ = octree.mortonReordering()
      }
      let endOctree = cross_platform_media_time()
      
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
        output.append("\(trialName) | \(reprOctree) | \(reprGrid) | \(reprParallel)")
      } else {
        output.append("\(trialName) | \(reprOctree) | \(reprGrid)")
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
  
  // The infamous nanofactory back board that took ~1000 ms to compile.
  func testBackBoard() throws {
    let reportingPerformance = Self.printPerformanceSummary
    
    var smallLeft = BackBoardSmallLeft()
    smallLeft.compile(reportingPerformance: reportingPerformance)
    XCTAssertEqual(smallLeft.topology.atoms.count, 36154)
    
    var smallRight = BackBoardSmallRight()
    smallRight.compile(reportingPerformance: reportingPerformance)
    XCTAssertEqual(smallRight.topology.atoms.count, 21324)
    
    var large = BackBoardLarge()
    large.compile(reportingPerformance: reportingPerformance)
    XCTAssertEqual(large.topology.atoms.count, 360_350)
  }
  
  // Old Diamondoid API:
  // lattice 1 - 62.4 ms
  // lattice 2 - 61.9 ms
  // lattice 3 - 537.2 ms
  //
  // New Topology API:
  // lattice 1 - 21.4 ms
  // lattice 2 - 33.3 ms
  // lattice 3 - 335.9 ms
  //
  // With bug fixes to sorting:
  // lattice 1 - 16.3 ms - 36,154 atoms
  // lattice 2 - 31.1 ms - 21,324 atoms
  // lattice 3 - 106.5 ms - 360,350 atoms
  
  // Before some optimizations to 'lattice' and 'orbitals':
  // -    atoms: 36,154
  // -  lattice:  7.0 ms
  // -    match:  3.9 ms
  // - orbitals:  2.7 ms
  // -    total: 13.6 ms
  //
  // -    atoms: 21,324
  // -  lattice: 18.8 ms
  // -    match:  2.4 ms
  // - orbitals:  1.5 ms
  // -    total: 22.7 ms
  //
  // -    atoms: 360,350
  // -  lattice: 34.8 ms
  // -    match: 42.0 ms
  // - orbitals: 30.0 ms
  // -    total: 106.9 ms
  
  // After some optimizations to 'lattice':
  //
  // -    atoms: 36,154
  // -  lattice:  4.3 ms
  // -    match:  4.0 ms
  // - orbitals:  2.8 ms
  // -    total: 11.1 ms
  //
  // -    atoms: 21,324
  // -  lattice:  8.5 ms
  // -    match:  2.5 ms
  // - orbitals:  1.5 ms
  // -    total: 12.6 ms
  //
  // -    atoms: 360,350
  // -  lattice: 15.0 ms
  // -    match: 39.2 ms
  // - orbitals: 28.3 ms
  // -    total: 82.5 ms
  
  // After the first round of optimizations to 'orbitals':
  // -    atoms: 36,154
  // -  lattice:  4.3 ms
  // -    match:  3.9 ms
  // - orbitals:  2.2 ms
  // -    total: 10.4 ms
  //
  // -    atoms: 21,324
  // -  lattice:  8.3 ms
  // -    match:  2.7 ms
  // - orbitals:  1.3 ms
  // -    total: 12.3 ms
  //
  // -    atoms: 360,350
  // -  lattice: 15.2 ms
  // -    match: 39.2 ms
  // - orbitals: 23.1 ms
  // -    total: 77.5 ms
  
  // After the second round of optimizations to 'orbitals':
  // -    atoms: 36,154
  // -  lattice:  4.3 ms
  // -    match:  3.8 ms
  // - orbitals:  1.9 ms
  // -    total: 10.0 ms
  //
  // -    atoms: 21,324
  // -  lattice:  8.4 ms
  // -    match:  2.2 ms
  // - orbitals:  1.1 ms
  // -    total: 11.7 ms
  //
  // -    atoms: 360,350
  // -  lattice: 15.1 ms
  // -    match: 40.2 ms
  // - orbitals: 19.2 ms
  // -    total: 74.5 ms
  
  // After the third round of optimizations to 'orbitals':
  // -    atoms: 36,154
  // -  lattice:  4.0 ms
  // -    match:  3.8 ms
  // - orbitals:  1.7 ms
  // -    total:  9.6 ms
  //
  // -    atoms: 21,324
  // -  lattice:  8.3 ms
  // -    match:  2.3 ms
  // - orbitals:  1.0 ms
  // -    total: 11.7 ms
  //
  // -    atoms: 360,350
  // -  lattice: 15.0 ms
  // -    match: 39.9 ms
  // - orbitals: 18.1 ms
  // -    total: 73.1 ms
  
  // After using a custom type for orbitals storage:
  // -    atoms: 36,154
  // -  lattice:  4.0 ms
  // -    match:  3.8 ms
  // - orbitals:  1.6 ms
  // -    total:  9.3 ms
  //
  // -    atoms: 21,324
  // -  lattice:  8.8 ms
  // -    match:  2.2 ms
  // - orbitals:  1.0 ms
  // -    total: 11.9 ms
  //
  // -    atoms: 360,350
  // -  lattice: 15.6 ms
  // -    match: 39.4 ms
  // - orbitals: 14.3 ms
  // -    total: 69.4 ms
  
  // After changing the storage of 'match':
  // -    atoms: 36,154
  // -  lattice:  3.9 ms
  // -    match:  3.8 ms
  // - orbitals:  1.7 ms
  // -    total:  9.4 ms
  //
  // -    atoms: 21,324
  // -  lattice:  7.9 ms
  // -    match:  2.1 ms
  // - orbitals:  0.9 ms
  // -    total: 10.9 ms
  //
  // -    atoms: 360,350
  // -  lattice: 15.1 ms
  // -    match: 38.0 ms
  // - orbitals: 15.5 ms
  // -    total: 68.7 ms
  
  // After applying a few more optimizations (optimized):
  // -    atoms: 36,154
  // -  lattice:  3.9 ms
  // -    match:  3.6 ms
  // - orbitals:  1.7 ms
  // -    total:  9.2 ms
  //
  // -    atoms: 21,324
  // -  lattice:  8.2 ms
  // -    match:  2.2 ms
  // - orbitals:  0.9 ms
  // -    total: 11.4 ms
  //
  // -    atoms: 360,350
  // -  lattice: 14.7 ms
  // -    match: 37.7 ms
  // - orbitals: 15.1 ms
  // -    total: 67.5 ms
  
  // Out of the box performance:
  //
  // -    atoms: 36,154
  // -  lattice:  4.1 ms
  // -    match:  4.1 ms
  // - orbitals:  1.6 ms
  // -    total:  9.8 ms
  //
  // -    atoms: 21,324
  // -  lattice:  8.2 ms
  // -    match:  2.3 ms
  // - orbitals:  1.0 ms
  // -    total: 11.5 ms
  //
  // -    atoms: 360,350
  // -  lattice: 14.4 ms
  // -    match: 40.6 ms
  // - orbitals: 15.6 ms
  // -    total: 70.7 ms
  
  // With the storage modification reverted:
  //
  // -    atoms: 36,154
  // -  lattice:  4.0 ms
  // -    match:  4.3 ms
  // - orbitals:  1.7 ms
  // -    total: 10.0 ms
  //
  // -    atoms: 21,324
  // -  lattice:  8.2 ms
  // -    match:  2.3 ms
  // - orbitals:  1.1 ms
  // -    total: 11.6 ms
  //
  // -    atoms: 360,350
  // -  lattice: 15.0 ms
  // -    match: 39.4 ms
  // - orbitals: 14.5 ms
  // -    total: 69.0 ms
  
  // Overall speedup:
  // lattice 1 -  62.4 ms -> 10.0 ms (6.24x)
  // lattice 2 -  61.9 ms -> 11.6 ms (5.33x)
  // lattice 3 - 537.2 ms -> 69.0 ms (7.79x)
  
#endif
}
