import XCTest
import HDL
import Numerics
import System

private let startTime = ContinuousClock.now

private func cross_platform_media_time() -> Double {
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
  
#if !DEBUG
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
  
  // This test compares the grid sorter to an alternative, recursive
  // implementation. It first measures cold-start speed, then speed once the
  // grid sorter can take advantage of the set already being sorted.
  //
  // The eventually chosen implementation might be a hybrid of these. It might
  // temporarily generate Morton indices to check whether segments of the atoms
  // are already sorted. It could also switch between different methods at
  // different levels of the hierarchy.
  func testSortingSpeed() throws {
    
  }
  
  // We need to run performance tests of Topology.match, to ensure the
  // acceleration algorithm is working properly. I could imagine subtle bugs
  // that make it incorrect, resulting in O(n^2) scaling.
#endif
}
