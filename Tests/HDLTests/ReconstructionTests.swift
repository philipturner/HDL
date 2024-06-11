import XCTest
import HDL

final class ReconstructionTests: XCTestCase {
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
    
    var reconstruction = Reconstruction()
    reconstruction.material = .checkerboard(.carbon, .germanium)
    reconstruction.topology.insert(atoms: lattice.atoms)
    reconstruction.compile()
    
    let topology = reconstruction.topology
    XCTAssertEqual(topology.atoms.count, 860)
    XCTAssertEqual(topology.bonds.count, 1318)
  }
  
  #if RELEASE
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
    reconstruction.material = .elemental(.carbon)
    reconstruction.topology.insert(atoms: lattice.atoms)
    reconstruction.compile()
    
    let topology = reconstruction.topology
    XCTAssertEqual(topology.atoms.count, 5735)
    XCTAssertEqual(topology.bonds.count, 9751)
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
    reconstruction.material = .elemental(.carbon)
    reconstruction.topology.insert(atoms: lattice.atoms)
    reconstruction.compile()
    
    let topology = reconstruction.topology
    XCTAssertEqual(topology.atoms.count, 5720)
    XCTAssertEqual(topology.bonds.count, 9697)
  }
  #endif
}
