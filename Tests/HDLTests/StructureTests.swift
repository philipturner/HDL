import XCTest
import HDL

final class StructureTests: XCTestCase {
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
        $0.element == element
      }))
    }
  }
  
  func testBackpropagationFailure() throws {
    let lattice = Lattice<Cubic> { h, k, l in
      Bounds { 5 * h + 5 * k + 5 * l }
      Material { .elemental(.carbon) }
      
      Volume {
        Replace { .empty }
      }
    }
    
    XCTAssertEqual(
      lattice.atoms.count,
      (6 * 6 * 6) + (5 * 5 * 5) + 3 * (5 * 5 * 11))
  }
  
  func testElement() throws {
    var permittedAtomicNumbers: [UInt8] = []
    permittedAtomicNumbers += [1, 6, 7, 8, 9]
    permittedAtomicNumbers += [13, 14, 15, 16, 17]
    permittedAtomicNumbers += [31, 32, 33, 34, 35]
    permittedAtomicNumbers += [50, 79, 82]
    
    for atomicNumber in permittedAtomicNumbers {
      var description: String
      var covalentRadius: Float // in picometers
      
      switch atomicNumber {
      case 1:
        description = ".hydrogen"
        covalentRadius = 31
      case 6:
        description = ".carbon"
        covalentRadius = 76
      case 7:
        description = ".nitrogen"
        covalentRadius = 71
      case 8:
        description = ".oxygen"
        covalentRadius = 66
      case 9:
        description = ".fluorine"
        covalentRadius = 57
        
      case 13:
        description = ".aluminum"
        covalentRadius = 121
      case 14:
        description = ".silicon"
        covalentRadius = 111
      case 15:
        description = ".phosphorus"
        covalentRadius = 107
      case 16:
        description = ".sulfur"
        covalentRadius = 105
      case 17:
        description = ".chlorine"
        covalentRadius = 102
        
      case 31:
        description = ".gallium"
        covalentRadius = 122
      case 32:
        description = ".germanium"
        covalentRadius = 120
      case 33:
        description = ".arsenic"
        covalentRadius = 119
      case 34:
        description = ".selenium"
        covalentRadius = 120
      case 35:
        description = ".bromine"
        covalentRadius = 120
        
      case 50:
        description = ".tin"
        covalentRadius = 139
      case 79:
        description = ".gold"
        covalentRadius = 136
      case 82:
        description = ".lead"
        covalentRadius = 146
      default:
        fatalError("Unrecognized atomic number: \(atomicNumber)")
      }
      
      guard let element = Element(rawValue: UInt8(atomicNumber)) else {
        XCTAssert(false, "'Element' failed to initialize.")
        continue
      }
      XCTAssertEqual(description, element.description)
      XCTAssertEqual(Int(covalentRadius), Int(element.covalentRadius * 1e3))
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
#if RELEASE
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
#if RELEASE
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
  
  func testShellStructure() throws {
    var elements: [Element] = []
    elements.append(.silicon)
    #if RELEASE
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
  
  func testGrapheneThiol() throws {
    var structure = GrapheneThiol()
    func filter(_ atomicNumber: UInt8) -> [UInt32] {
      var output: [UInt32] = []
      for i in structure.topology.atoms.indices {
        let atom = structure.topology.atoms[i]
        if atom.atomicNumber == atomicNumber {
          output.append(UInt32(i))
        }
      }
      return output
    }
    func withOrbitalCount(_ closure: (Int) -> Void) {
      let orbitals = structure.topology.nonbondingOrbitals(hybridization: .sp2)
      let count = orbitals.reduce(0) { $0 + $1.count }
      closure(count)
    }
    
    structure.compilationPass0()
    XCTAssertEqual(structure.topology.atoms.count, 54)
    XCTAssertEqual(structure.topology.bonds.count, 0)
    XCTAssertEqual(filter(6).count, 54)
    withOrbitalCount { orbitalCount in
      XCTAssertEqual(orbitalCount, 0)
    }
    
    structure.compilationPass1()
    XCTAssertEqual(structure.topology.atoms.count, 54)
    XCTAssertEqual(structure.topology.bonds.count, 0)
    withOrbitalCount { orbitalCount in
      XCTAssertEqual(orbitalCount, 0)
    }
    
    structure.compilationPass2()
    XCTAssertEqual(structure.topology.atoms.count, 54)
    XCTAssertEqual(structure.topology.bonds.count, 71)
    withOrbitalCount { orbitalCount in
      XCTAssertGreaterThan(orbitalCount, 0)
    }
    
    structure.compilationPass3()
    XCTAssertEqual(structure.topology.atoms.count, 78)
    XCTAssertEqual(structure.topology.bonds.count, 95)
    XCTAssertEqual(filter(1).count, 20)
    XCTAssertEqual(filter(6).count, 54)
    XCTAssertEqual(filter(16).count, 4)
    withOrbitalCount { orbitalCount in
      XCTAssertEqual(orbitalCount, 0)
    }
  }
  
  func testCBNTripod() throws {
    var tripod = CBNTripod()
    var atomRecord: [UInt8: Int]
    
    func createAtomRecord(_ tripod: CBNTripod) -> [UInt8: Int] {
      var output: [UInt8: Int] = [:]
      for atom in tripod.createAtoms() {
        var previous = output[atom.atomicNumber] ?? 0
        previous += 1
        output[atom.atomicNumber] = previous
      }
      return output
    }
    
    tripod.rotateLegs(slantAngleDegrees: 62, swingAngleDegrees: 5)
    atomRecord = createAtomRecord(tripod)
    XCTAssertEqual(tripod.createAtoms().count, 63)
    XCTAssertEqual(atomRecord.keys.count, 6)
    XCTAssertEqual(atomRecord[1], 18)
    XCTAssertEqual(atomRecord[6], 29)
    XCTAssertEqual(atomRecord[7], 3)
    XCTAssertEqual(atomRecord[9], 9)
    XCTAssertEqual(atomRecord[14], 3)
    XCTAssertEqual(atomRecord[32], 1)
    
    tripod.passivateNHGroups(.hydrogen)
    atomRecord = createAtomRecord(tripod)
    XCTAssertEqual(tripod.createAtoms().count, 63)
    XCTAssertEqual(atomRecord.keys.count, 5)
    XCTAssertEqual(atomRecord[1], 21)
    XCTAssertEqual(atomRecord[6], 29)
    XCTAssertEqual(atomRecord[7], 3)
    XCTAssertEqual(atomRecord[9], 9)
    XCTAssertEqual(atomRecord[32], 1)
    
    tripod.passivateNHGroups(.silicon)
    atomRecord = createAtomRecord(tripod)
    XCTAssertEqual(tripod.createAtoms().count, 63)
    XCTAssertEqual(atomRecord.keys.count, 6)
    XCTAssertEqual(atomRecord[1], 18)
    XCTAssertEqual(atomRecord[6], 29)
    XCTAssertEqual(atomRecord[7], 3)
    XCTAssertEqual(atomRecord[9], 9)
    XCTAssertEqual(atomRecord[14], 3)
    XCTAssertEqual(atomRecord[32], 1)
  }
  
  // This is disabled in debug mode, due to the high computational prefactor
  // of the numerous lattice intersections.
#if RELEASE
  // This test is to ensure there are no bugs, if one attempts to optimize
  // intersections in 'Lattice<Hexagonal>' through sparsity.
  func testRodLogicHousing() throws {
    var housing = RodLogicHousing()
    XCTAssertEqual(housing.topology.atoms.count, 0)
    
    housing.compilationPass0()
    XCTAssertEqual(housing.topology.atoms.count, 14454)
  }
#endif
}
