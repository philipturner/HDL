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
    XCTAssertEqual(lattice.entities.count, 57601)
    
    let overallEnd = cross_platform_media_time()
    print("gold surface:")
    print("- overall:   \(fmt(overallStart, overallEnd))")
    print("- grid init: \(fmt(gridInitStart, gridInitEnd))")
    print("- intersect: \(fmt(intersectStart, intersectEnd))")
    print("- replace:   \(fmt(replaceStart, replaceEnd))")
    
    // Before optimizations: ~0.318 seconds
    // After optimization 1: ~0.066 seconds
    // After optimization 2: ~0.059 seconds
    // After optimization 3: ~0.047 seconds
    // 6.8x speedup
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
    XCTAssertEqual(lattice.entities.count, 81142)
    
    let overallEnd = cross_platform_media_time()
    print("silicon probe:")
    print("- overall:   \(fmt(overallStart, overallEnd))")
    print("- grid init: \(fmt(gridInitStart, gridInitEnd))")
    print("- intersect: \(fmt(intersectStart, intersectEnd))")
    print("- replace:   \(fmt(replaceStart, replaceEnd))")
    
    // Before optimizations: ~0.300 seconds
    // After optimization 1: ~0.047 seconds
    // After optimization 2: ~0.031 seconds
    // 9.7x speedup
  }
#endif
}
