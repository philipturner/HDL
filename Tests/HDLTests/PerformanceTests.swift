import XCTest
import HDL
import Numerics

final class PerformanceTests: XCTestCase {
#if RELEASE
  func testGoldSurface() throws {
    let scaleFactor: Float = 2
    let lattice = Lattice<Cubic> { h, k, l in
      Bounds { scaleFactor * 40 * (h + k + l) }
      Material { .elemental(.gold) }
      
      Volume {
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
        
        Replace { .empty }
      }
    }
    XCTAssertEqual(lattice.atoms.count, 57601)
  }
  
  func testSiliconProbe() throws {
    let lattice = Lattice<Cubic> { h, k, l in
      Bounds { 50 * (h + k + l) }
      Material { .elemental(.silicon) }
      
      Volume {
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
        
        Replace { .empty }
      }
    }
    XCTAssertEqual(lattice.atoms.count, 81142)
  }
  
  func testGoldSurface2() throws {
    let scaleFactor: Float = 6
    let lattice = Lattice<Cubic> { h, k, l in
      Bounds { scaleFactor * 40 * (h + k + l) }
      Material { .elemental(.gold) }
      
      Volume {
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
        
        Replace { .empty }
      }
    }
    XCTAssertEqual(lattice.atoms.count, 173517)
  }
  
  func testSiliconProbe2() throws {
    let lattice = Lattice<Cubic> { h, k, l in
      Bounds { 80 * (h + k + l) }
      Material { .elemental(.silicon) }
      let topCutoff: Float = 67
      
      Volume {
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
        
        Replace { .empty }
        
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
        
        Replace { .atom(.carbon) }
      }
    }
    XCTAssertEqual(lattice.atoms.count, 64226)
  }
  
  func testSort() throws {
    let latticeScale: Float = 10
    let lattice = Lattice<Hexagonal> { h, k, l in
      let h2k = h + 2 * k
      Bounds { latticeScale * (2 * h + h2k + l) }
      Material { .elemental(.carbon) }
    }
    
    for trialID in 0..<4 {
      var trialAtoms: [Atom]
      
      switch trialID {
      case 0:
        var topology = Topology()
        topology.insert(atoms: lattice.atoms)
        topology.sort()
        
        trialAtoms = topology.atoms
      case 1:
        trialAtoms = lattice.atoms
      case 2:
        trialAtoms = lattice.atoms.shuffled()
      case 3:
        trialAtoms = lattice.atoms.reversed()
      default:
        fatalError("This should never happen.")
      }
      
      var resultGrid1: [UInt32] = []
      var resultGrid2: [UInt32] = []
      do {
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
      
      var topology = Topology()
      topology.insert(atoms: trialAtoms)
      let resultGrid = topology.sort()
      do {
        var topology = Topology()
        topology.insert(atoms: trialAtoms)
        _ = topology.sort()
      }
      
      let octree = OctreeSorter(atoms: trialAtoms)
      let resultOctree = octree.mortonReordering()
      do {
        let octree = OctreeSorter(atoms: trialAtoms)
        _ = octree.mortonReordering()
      }
      
      XCTAssertEqual(resultGrid, resultOctree)
      XCTAssertEqual(resultGrid1, resultOctree)
      XCTAssertEqual(resultGrid2, resultOctree)
    }
  }
  
  // The infamous nanofactory back board that took ~1000 ms to compile.
  func testBackBoard() throws {
    var smallLeft = BackBoardSmallLeft()
    var smallRight = BackBoardSmallRight()
    var large = BackBoardLarge()
    
    smallLeft.compile()
    smallRight.compile()
    large.compile()
    
    XCTAssertEqual(smallLeft.topology.atoms.count, 36154)
    XCTAssertEqual(smallRight.topology.atoms.count, 21324)
    XCTAssertEqual(large.topology.atoms.count, 360_350)
  }
#endif
}
