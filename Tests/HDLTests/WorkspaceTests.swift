import HDL
import XCTest

// Temporary file for treating the Swift tests like a scratch workspace.

final class WorkspaceTests: XCTestCase {
  func testWorkspace() throws {
    let lattice = Lattice<Cubic> { h, k, l in
      Bounds { 20 * (h + k + l) }
      Material { .checkerboard(.silicon, .carbon) }
      
      let beamWidth: Float = 2
      
      Volume {
        Concave {
          Convex {
            Origin { beamWidth * h }
            Plane { h }
          }
          Convex {
            Origin { beamWidth * k }
            Plane { k }
          }
        }
        
        Concave {
          Convex {
            Origin { beamWidth * h }
            Plane { h }
          }
          Convex {
            Origin { beamWidth * l }
            Plane { l }
          }
        }
        
        Concave {
          Convex {
            Origin { beamWidth * k }
            Plane { k }
          }
          Convex {
            Origin { beamWidth * l }
            Plane { l }
          }
        }
        
        Replace { .empty }
      }
    }
    print(lattice.atoms.count)
  }
}
