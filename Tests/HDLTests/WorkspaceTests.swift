import HDL
import XCTest

// Temporary file for treating the Swift tests like a scratch workspace.

final class WorkspaceTests: XCTestCase {
  func testWorkspace() throws {
    // 40 / 4 -  16281 atoms, ~3 ms
    // 80 / 8 - 122225 atoms, ~9 ms
    let lattice = Lattice<Cubic> { h, k, l in
      Bounds { 80 * (h + k + l) }
      Material { .checkerboard(.silicon, .carbon) }
      
      let beamWidth: Float = 8
      
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
    
    // 40 / 4 -  15806 atoms, ~14 ms
    // 80 / 8 - 121270 atoms, ~94 ms -> ~91 ms -> ~90 ms
    var reconstruction = Reconstruction()
    reconstruction.atoms = lattice.atoms
    reconstruction.material = .checkerboard(.silicon, .carbon)
    let topology = reconstruction.compile()
    print(topology.atoms.count)
  }
}
