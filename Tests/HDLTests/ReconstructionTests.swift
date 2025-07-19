import HDL
import XCTest

final class ReconstructionTests: XCTestCase {
  static let printPerformanceSummary = false
  
  func testUnitTest() throws {
    let lattice = Lattice<Cubic> { h, k, l in
      Bounds { 5 * h + 5 * k + 5 * l }
      Material { .checkerboard(.carbon, .germanium) }
      
      Volume {
        Concave {
          Origin { 1.5 * k + 1.5 * l }
          
          // Create a groove for the rod.
          Concave {
            Plane { k }
            Plane { l }
          }
        }
        Replace { .empty }
      }
    }
    XCTAssertEqual(lattice.atoms.count, 651)
    
    var reconstruction = Reconstruction()
    reconstruction.atoms = lattice.atoms
    reconstruction.material = .checkerboard(.carbon, .germanium)
    var topology = reconstruction.compile()
    PassivationTests.passivate(topology: &topology)
    XCTAssertEqual(topology.atoms.count, 860)
    XCTAssertEqual(topology.bonds.count, 1318)
    PassivationTests.checkConnectivity(topology)
    PassivationTests.checkNoOverlaps(topology)
  }
  
  func testAluminumPhosphide() throws {
    let lattice = Lattice<Cubic> { h, k, l in
      Bounds { 3 * h + 3 * k + 3 * l }
      Material { .checkerboard(.phosphorus, .aluminum) }
    }
    XCTAssertEqual(lattice.atoms.count, 280)
    
    var reconstruction = Reconstruction()
    reconstruction.atoms = lattice.atoms
    reconstruction.material = .checkerboard(.phosphorus, .aluminum)
    var topology = reconstruction.compile()
    PassivationTests.passivate(topology: &topology)
    XCTAssertEqual(topology.atoms.count, 384)
    XCTAssertEqual(topology.bonds.count, 564)
    PassivationTests.checkConnectivity(topology)
    PassivationTests.checkNoOverlaps(topology)
    
    let orbitalLists = topology.nonbondingOrbitals()
    var hasUnfilledValences = false
    for atomID in topology.atoms.indices {
      let orbitalList = orbitalLists[atomID]
      if orbitalList.count > 0 {
        hasUnfilledValences = true
      }
    }
    XCTAssertFalse(hasUnfilledValences)
  }
  
#if RELEASE
  func testDoubleCompile() throws {
    let lattice = Lattice<Hexagonal> { h, k, l in
      let h2k = h + 2 * k
      Bounds { 4 * h + 4 * h2k + 4 * l }
      Material { .checkerboard(.silicon, .carbon) }
      
      Volume {
        Origin { 2 * h + 1.6 * l }
        
        // Create a corner cut.
        Concave {
          Plane { h }
          Plane { l }
        }
        Replace { .empty }
      }
    }
    XCTAssertEqual(lattice.atoms.count, 432)
    
    var reconstruction = Reconstruction()
    reconstruction.atoms = lattice.atoms
    reconstruction.material = .checkerboard(.silicon, .carbon)
    
    for _ in 0..<2 {
      var topology = reconstruction.compile()
      PassivationTests.passivate(topology: &topology)
      XCTAssertEqual(topology.atoms.count, 658)
      XCTAssertEqual(topology.bonds.count, 959)
      PassivationTests.checkConnectivity(topology)
      
      var groupIVAtomCount: Int = .zero
      for atom in topology.atoms {
        if atom.atomicNumber != 1 {
          groupIVAtomCount += 1
        }
      }
      XCTAssertEqual(groupIVAtomCount, 420)
    }
  }
  
  func testReproducerBefore() throws {
    let lattice = Lattice<Cubic> { h, k, l in
      Bounds { 10 * h + 9 * k + 7 * l }
      Material { .elemental(.carbon) }
      
      Volume {
        Concave {
          Origin { 1.5 * k + 1.5 * l }
          
          // Create a groove for the rod.
          Concave {
            Plane { k }
            Plane { l }
            Origin { 6.25 * k + 4 * l }
            Plane { -k }
            Plane { -l }
          }
          
          Concave {
            Origin { 2 * h }
            
            // Create a 45-degree inclined plane.
            Plane { h - k }
            
            // Correct the walls of the shaft.
            Convex {
              Origin { 4 * l }
              Origin { 0.25 * (k - l) }
              Plane { k - l }
            }
            
            // Correct the walls of the shaft.
            Convex {
              Origin { 0.25 * (k + l) }
              Plane { k + l }
            }
            
            // Correct the concave corner of the drive wall site.
            Convex {
              Origin { 1.75 * h }
              Plane { h }
            }
            
            // Correct the concave corner of the drive wall site.
            Convex {
              Origin { 1.75 * h + 0.5 * (h + k + l) }
              Plane { h + k + l }
            }
          }
        }
        
        Replace { .empty }
      }
    }
    
    var reconstruction = Reconstruction()
    reconstruction.atoms = lattice.atoms
    reconstruction.material = .elemental(.carbon)
    
    let start = Profiler.time()
    var topology = reconstruction.compile()
    let end = Profiler.time()
    if Self.printPerformanceSummary {
      // expected: 3.1 ms
      let elapsedSeconds = end - start
      let elapsedMilliseconds = 1000 * elapsedSeconds
      let formatted = String(format: "%.1f", elapsedMilliseconds)
      print("reproducerBefore: \(formatted) ms")
    }
    
    PassivationTests.passivate(topology: &topology)
    XCTAssertEqual(topology.atoms.count, 5735)
    XCTAssertEqual(topology.bonds.count, 9751)
    PassivationTests.checkConnectivity(topology)
    PassivationTests.checkNoOverlaps(topology)
  }
  
  func testReproducerAfter() throws {
    let lattice = Lattice<Cubic> { h, k, l in
      Bounds { 10 * h + 9 * k + 7 * l }
      Material { .elemental(.carbon) }
      
      Volume {
        Concave {
          Origin { 1.5 * k + 1.5 * l }
          
          // Create a groove for the rod.
          Concave {
            Plane { k }
            Plane { l }
            Origin { 6.25 * k + 4 * l }
            Plane { -k }
            Plane { -l }
          }
          
          Concave {
            Origin { 2 * h }
            
            // Create a 45-degree inclined plane.
            Plane { h - k }
          }
        }
        
        Replace { .empty }
      }
    }
    
    var reconstruction = Reconstruction()
    reconstruction.atoms = lattice.atoms
    reconstruction.material = .elemental(.carbon)
    
    let start = Profiler.time()
    var topology = reconstruction.compile()
    let end = Profiler.time()
    if Self.printPerformanceSummary {
      // expected: 6.6 ms
      let elapsedSeconds = end - start
      let elapsedMilliseconds = 1000 * elapsedSeconds
      let formatted = String(format: "%.1f", elapsedMilliseconds)
      print("reproducerAfter: \(formatted) ms")
    }
    
    PassivationTests.passivate(topology: &topology)
    XCTAssertEqual(topology.atoms.count, 5720)
    XCTAssertEqual(topology.bonds.count, 9697)
    PassivationTests.checkConnectivity(topology)
    PassivationTests.checkNoOverlaps(topology)
  }
  
  // The previous test, but shifted by 1 micron to exacerbate FP32 error.
  func testShiftedLattice() throws {
    let lattice = Lattice<Cubic> { h, k, l in
      Bounds { 10 * h + 9 * k + 7 * l }
      Material { .elemental(.carbon) }
      
      Volume {
        Concave {
          Origin { 1.5 * k + 1.5 * l }
          
          // Create a groove for the rod.
          Concave {
            Plane { k }
            Plane { l }
            Origin { 6.25 * k + 4 * l }
            Plane { -k }
            Plane { -l }
          }
          
          Concave {
            Origin { 2 * h }
            
            // Create a 45-degree inclined plane.
            Plane { h - k }
          }
        }
        
        Replace { .empty }
      }
    }
    
    var reconstruction = Reconstruction()
    reconstruction.atoms = lattice.atoms.map {
      var output = $0
      output.position += 1000 * SIMD3<Float>(-1, 1, -1)
      return output
    }
    reconstruction.material = .elemental(.carbon)
    
    var topology = reconstruction.compile()
    PassivationTests.passivate(topology: &topology)
    XCTAssertEqual(topology.atoms.count, 5720)
    XCTAssertEqual(topology.bonds.count, 9697)
    PassivationTests.checkConnectivity(topology)
    PassivationTests.checkNoOverlaps(topology)
  }
#endif
}
