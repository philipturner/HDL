import XCTest
@testable import HDL

final class HDLTests: XCTestCase {
    func testAdamantane() throws {
      for element in [Element.carbon, Element.silicon] {
        let lattice = Lattice<Cubic> { h, k, l in
          Bounds { 4 * h + 4 * k + 4 * l }
          Material { .elemental(element) }
          
          Volume {
            Origin { 2 * h + 2 * k + 2 * l }
            Origin { 0.25 * (h + k - l) }
            
            // Remove the front plane.
            Convex {
              Origin { 0.25 * (h + k + l) }
              Plane { h + k + l }
            }
            
            func triangleCut(sign: Float) {
              Convex {
                Origin { 0.25 * sign * (h - k - l) }
                Plane { sign * (h - k / 2 - l / 2) }
              }
              Convex {
                Origin { 0.25 * sign * (k - l - h) }
                Plane { sign * (k - l / 2 - h / 2) }
              }
              Convex {
                Origin { 0.25 * sign * (l - h - k) }
                Plane { sign * (l - h / 2 - k / 2) }
              }
            }
            
            // Remove three sides forming a triangle.
            triangleCut(sign: +1)
            
            // Remove their opposites.
            triangleCut(sign: -1)
            
            // Remove the back plane.
            Convex {
              Origin { -0.25 * (h + k + l) }
              Plane { -(h + k + l) }
            }
            
            Replace { .empty }
          }
        }
        
        XCTAssertTrue(lattice.entities.count > 0)
        XCTAssertTrue(lattice.entities.contains(where: {
          $0.type == .atom(element)
        }))
      }
    }
}
