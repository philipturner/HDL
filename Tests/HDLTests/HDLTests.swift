import XCTest
#if DEBUG
@testable import HDL
#else
import HDL
#endif

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
      
      XCTAssertGreaterThan(lattice.atoms.count, 0)
      XCTAssertTrue(lattice.atoms.contains(where: {
        $0.type == .atom(element)
      }))
    }
  }
  
  func testGold() throws {
    let carbonLattice = Lattice<Cubic> { h, k, l in
      Bounds { 4 * (h + k + l) }
      Material { .elemental(.carbon) }
      
      Volume {
        Origin { 2 * (h + k + l) }
        Plane { -(h + k + l) }
        Replace { .empty }
      }
    }
    let goldLattice = Lattice<Cubic> { h, k, l in
      Bounds { 4 * (h + k + l) }
      Material { .elemental(.gold) }
      
      Volume {
        Origin { 2 * (h + k + l) }
        Plane { -(h + k + l) }
        Replace { .empty }
      }
    }
    XCTAssertGreaterThanOrEqual(
      goldLattice.atoms.count, carbonLattice.atoms.count / 2)
    
    for atom in goldLattice.atoms {
      // Map the position to a fraction of the unit cell.
      var position = atom.position
      position /= Constant(.square) { .elemental(.gold) }
      position *= 2
      
      // Eliminate floating-point error.
      let rounded = position.rounded(.toNearestOrEven)
      var difference = position - rounded
      difference.replace(with: -difference, where: difference .< 0)
      
      // Ensure the gold atoms are snapped to a multiple of 0.5 cells.
      XCTAssertTrue(all(difference .< 1e-3))
    }
  }
  
  func testPlane111() throws {
    var parameters: [SIMD3<Int>] = [
      SIMD3(3, 613, 1018),
      SIMD3(4, 1337, 2313),
      SIMD3(5, 2481, 4406),
    ]
#if !DEBUG
    parameters += [
      SIMD3(6, 4141, 7489),
      SIMD3(7, 6413, 11754),
      SIMD3(8, 9393, 17393),
      SIMD3(9, 13177, 24598),
      SIMD3(10, 17861, 33561),
      SIMD3(11, 23541, 44474),
    ]
#endif
    
    for parameter in parameters {
      let scale = Float(parameter[0])
      let carbonLattice = Lattice<Cubic> { h, k, l in
        Bounds { scale * 2 * (h + k + l) }
        Material { .elemental(.carbon) }
        
        Volume {
          Origin { scale * (h + k + l) }
          Plane { -(h + k + l) }
          Replace { .empty }
        }
      }
#if !DEBUG
      let goldLattice = Lattice<Cubic> { h, k, l in
        Bounds { scale * 2 * (h + k + l) }
        Material { .elemental(.gold) }
        
        Volume {
          Origin { scale * (h + k + l) }
          Plane { -(h + k + l) }
          Replace { .empty }
        }
      }
      XCTAssertEqual(goldLattice.atoms.count, parameter[1])
#endif
      XCTAssertEqual(carbonLattice.atoms.count, parameter[2])
    }
  }
  
  // TODO: Enable this in debug mode once the matching algorithm has been
  // sped up.
  #if !DEBUG
  func testShellStructure() throws {
    var elements: [Element] = []
    elements.append(.silicon)
    #if !DEBUG
    elements.append(.carbon)
    #endif
    
    for element in elements {
      var structure = ShellStructure(element: element)
      structure.compilationPass0()
      XCTAssertEqual(structure.topology.atoms.count, 894)
      XCTAssertEqual(structure.topology.bonds.count, 0)
      
      structure.compilationPass1()
      XCTAssertEqual(structure.topology.atoms.count, 1514)
      XCTAssertEqual(structure.topology.bonds.count, 2098)
      
      structure.compilationPass2()
      XCTAssertEqual(structure.topology.atoms.count, 1514)
      XCTAssertEqual(structure.topology.bonds.count, 2098)
      
      if element == .carbon {
        structure.compilationPass3(onlyMergeHydrogens: true)
        XCTAssertEqual(structure.topology.atoms.count, 1358)
        XCTAssertEqual(structure.topology.bonds.count, 1942)
      } else {
        structure.compilationPass3(onlyMergeHydrogens: false)
        XCTAssertEqual(structure.topology.atoms.count, 1098)
        XCTAssertEqual(structure.topology.bonds.count, 1539)
      }
    }
  }
  #endif
}
