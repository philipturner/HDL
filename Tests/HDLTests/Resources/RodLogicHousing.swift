//
//  RodLogicHousing.swift
//  HDLTests
//
//  Created by Philip Turner on 1/2/24.
//

import HDL

struct RodLogicHousing {
  var topology: Topology = .init()
  
  mutating func compilationPass0() {
    let lattice = createHousingLattice()
    topology.insert(atoms: lattice.atoms)
  }
}

private func createHousingLattice() -> Lattice<Hexagonal> {
  let housingLattice = Lattice<Hexagonal> { h, k, l in
    let h2k = h + 2 * k
    Bounds { 20 * h + 14 * h2k + 12 * l }
    Material { .elemental(.carbon) }
    
    Volume {
      Origin { 12 * h + 9 * h2k + 6 * l }
      
      // TODO: Always remember to comment your HDL code. Otherwise, it's
      // almost impossible to understand when looking back on it.
      
      // Cut the initial block into an L-shape.
      //
      // There's a compiler bug preventing me from wrapping
      // "Origin { 2.8 * l } - Plane { l }" neatly in a shared scope.
      Concave {
        Origin { 2.8 * l - 2.5 * h2k }
        Plane { l }
        Plane { -h2k }
      }
      
      // Cut a cool chiseled shape around the first rod's housing.
      for direction in [-2 * h - k, h - k] {
        Concave {
          Origin { 2.8 * l }
          if direction.x > 0 {
            Origin { 4.0 * direction }
          } else {
            Origin { 3.5 * direction }
          }
          Plane { l }
          Plane { direction }
        }
      }
      for direction in [h * 2 + k, -h + k] {
        Convex {
          Origin { 5 * direction }
          Plane { direction }
        }
      }
      Concave {
        Origin { 2.8 * l - 6.5 * h }
        Plane { l }
        Plane { -h }
      }
      
      // Chop off a slice of atoms that isn't needed for the second rod.
      Convex {
        Origin { -11 * h }
        Plane { -h }
      }
      
      // Create the hole for the first rod to go through.
      Concave {
        for direction in [h2k, -h2k] {
          Convex {
            if direction.y > 0 {
              Origin { 3 * direction }
            } else {
              Origin { 3.5 * direction }
            }
            Plane { -direction }
          }
        }
        for direction in [h, -h] {
          Convex {
            if direction.x > 0 {
              Origin { 4 * direction }
            } else {
              Origin { 3.5 * direction }
            }
            Plane { -direction }
          }
        }
        for direction in [h * 2 + k, -h * 2 - k] {
          Convex {
            if direction.y > 0 {
              if direction.x > 0 {
                Origin { 3 * direction }
              }
            } else {
              Origin { 3 * direction }
            }
            Plane { -direction }
          }
        }
        
        // Create the overhang that stops the first rod from falling out.
        //
        // It seems, the second rod needs to be placed far off-center. The
        // first rod doesn't move far enough for the second rod to run
        // through the middle.
        Volume {
          Concave {
            Origin { -0.5 * h2k }
            Plane { h2k }
            Origin { -2 * h }
            Plane { h + k }
          }
          Replace { .empty }
        }
        Volume {
          Concave {
            for direction in [l, -l] {
              Convex {
                if direction[2] > 0 {
                  Origin { 4 * direction }
                } else {
                  Origin { 4.2 * direction }
                }
                Plane { -direction }
              }
            }
          }
          Replace { .empty }
        }
        
        // Etch a protrusion that theoretically decreases the vdW binding
        // energy of the gate knob. This reduces the force necessary to
        // separate the rod from ~1600 pN to ~1200 pN, much less than
        // expected. The first experiment used 1600 pN in order to finish
        // faster than the break-even force, 1200 pN. Don't confuse that
        // 1600 pN with the ~1600 pN cited here.
        Volume {
          Concave {
            Plane { l }
            Origin { -1.5 * h2k + 4.75 * l }
            Plane { -l }
            Convex {
              Plane { h2k }
              for direction in [h, -h] {
                Convex {
                  Origin { 0.75 * direction }
                  Plane { direction }
                }
              }
            }
          }
          Replace { .empty }
        }
      }
    }
    
    createSecondHole(h, k, l)
  }
  return housingLattice
}

private func createSecondHole(
  _ h: SIMD3<Float>, _ k: SIMD3<Float>, _ l: SIMD3<Float>
) {
  Volume {
    let h2k = h + 2 * k
    Origin { 10 * h + 5 * h2k + 3 * l }
    
    Concave {
      for direction in [-l, l] {
        Convex {
          if direction.z > 0 {
            Origin { 2.8 * direction }
          } else {
            Origin { 1.2 * direction }
          }
          Plane { -direction }
        }
      }
      Convex {
        Origin { -2.5 * h2k }
        Plane { h2k }
        
        // Remove some atoms that are an intense source of energy dissipation.
        for direction in [k + h, k] {
          Convex {
            Origin { 4 * direction }
            Plane { direction }
          }
        }
      }
      Convex {
        Origin { h2k }
        Plane { -h2k }
      }
    }
    Replace { .empty }
  }
}
