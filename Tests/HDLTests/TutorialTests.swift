import XCTest
import HDL
import Numerics

final class TutorialTests: XCTestCase {
  // Swift implementation of the tutorial.
  func testTutorial() throws {
    let carbonLattice = Self.step1()
    let carbons = Self.step2(carbonLattice: carbonLattice)
    let siliconLattice = Self.step3()
    let initialSilicons = Self.step4(siliconLattice: siliconLattice)
    let finalSilicons = Self.step5(input: initialSilicons)
    
    let actualStep1 = Self.exportToXYZ(carbonLattice.atoms, comment: "Step 1")
    let actualStep2 = Self.exportToXYZ(carbons, comment: "Step 2")
    let actualStep3 = Self.exportToXYZ(siliconLattice.atoms, comment: "Step 3")
    let actualStep4 = Self.exportToXYZ(initialSilicons, comment: "Step 4")
    let actualStep5 = Self.exportToXYZ(carbons + finalSilicons, comment: "Step 5")
    
    XCTAssertEqual(actualStep1, expectedStep1)
    XCTAssertEqual(actualStep2, expectedStep2)
    XCTAssertEqual(actualStep3, expectedStep3)
    XCTAssertEqual(actualStep4, expectedStep4)
    XCTAssertEqual(actualStep5, expectedStep5)
  }
  
  // Utility function for exporting atoms.
  static func exportToXYZ(_ atoms: [Atom], comment: String = "") -> String {
    var output: String = ""
    output += String(describing: atoms.count)
    output += "\n"
    output += comment
    output += "\n"
    
    for atom in atoms {
      var elementSymbol: String
      switch atom.element {
      case .hydrogen: elementSymbol = "H "
      case .carbon:   elementSymbol = "C "
      case .silicon:  elementSymbol = "Si"
      default: fatalError("Unrecognized element symbol: \(atom.element)")
      }
      output += elementSymbol
      
      let position: SIMD3<Float> = atom.position
      for vectorLane in 0..<3 {
        // Convert the coordinate from nm -> angstrom.
        let coordinate = position[vectorLane] * 10
        output += " "
        output += String(format: "%.3f", coordinate)
      }
      output += "\n"
    }
    return output
  }
  
  // Code for step 1.
  static func step1() -> Lattice<Hexagonal> {
    let carbonLattice = Lattice<Hexagonal> { h, k, l in
      let h2k = h + 2 * k
      Bounds { 4 * h + 3 * h2k + 1 * l }
      Material { .elemental(.carbon) }
      
      Volume {
        // Move the player position from the origin to (0, 0, 0.25).
        Origin { 0.25 * l }
        
        // Create a plane pointing from the origin to positive `l`.
        Plane { l }
        
        // Remove all atoms on the positive side of the plane.
        Replace { .empty }
      }
    }
    
    return carbonLattice
  }
  
  // Code for step 2.
  static func step2(
    carbonLattice: Lattice<Hexagonal>
  ) -> [Atom] {
    var grapheneHexagonScale: Float
    do {
      // Convert graphene lattice constant from Å to nm.
      let grapheneConstant: Float = 2.45 / 10
      
      // Retrieve lonsdaleite lattice constant in nm.
      let lonsdaleiteConstant = Constant(.hexagon) { .elemental(.carbon) }
      
      // Each hexagon's current side length is the value of
      // `lonsdaleiteConstant`. Dividing by this constant, changes the hexagon
      // so its sides are all 1 nm.
      grapheneHexagonScale = 1 / lonsdaleiteConstant
      
      // Multiply by the graphene constant. This second transformation stretches
      // the hexagon, so its sides are all 0.245 nm.
      grapheneHexagonScale *= grapheneConstant
    }
    
    var carbons: [Atom] = carbonLattice.atoms
    for atomID in carbons.indices {
      // Flatten the sp3 sheet into an sp2 sheet.
      carbons[atomID].position.z = 0
      
      // Resize the hexagon side length, so it matches graphene.
      carbons[atomID].position.x *= grapheneHexagonScale
      carbons[atomID].position.y *= grapheneHexagonScale
    }
    
    return carbons
  }
  
  // Code for step 3.
  static func step3() -> Lattice<Hexagonal> {
    let siliconLattice = Lattice<Hexagonal> { h, k, l in
      let h2k = h + 2 * k
      Bounds { 3 * h + 2 * h2k + 1 * l }
      Material { .elemental(.silicon) }
      
      Volume {
        Origin { 0.25 * l }
        Plane { l }
        Replace { .empty }
      }
    }
    
    return siliconLattice
  }
  
  // Code for step 4.
  static func step4(
    siliconLattice: Lattice<Hexagonal>
  ) -> [Atom] {
    var siliceneHexagonScale: Float
    do {
      // Convert silicene lattice constant from Å to nm.
      let siliceneConstant: Float = 3.75 / 10
      
      // Retrieve the constant for 3D lonsdaleite-shaped silicon, in nm.
      let lonsdaleiteConstant = Constant(.hexagon) { .elemental(.silicon) }
      
      // Create a number that maps from 3D lattice spacing to silicene
      // lattice spacing.
      siliceneHexagonScale = siliceneConstant / lonsdaleiteConstant
    }
    
    var silicons: [Atom] = siliconLattice.atoms
    for atomID in silicons.indices {
      // Partially flatten the sp3 sheet, so the elevated atoms reach the
      // buckling distance from the literature.
      if silicons[atomID].position.z > 0 {
        silicons[atomID].position.z = 0.62 / 10
      }
      
      // Resize the hexagon side length, so it matches silicene.
      silicons[atomID].position.x *= siliceneHexagonScale
      silicons[atomID].position.y *= siliceneHexagonScale
    }
    
    return silicons
  }
  
  // Code for step 5.
  static func step5(
    input: [Atom]
  ) -> [Atom] {
    var silicons = input
    
    var rotation: Quaternion<Float>
    do {
      // Convert the angle from degrees to radians.
      let angle: Float = 10.9 * .pi / 180
      
      // Create a counterclockwise rotation around the Z axis.
      rotation = Quaternion(angle: angle, axis: [0, 0, 1])
    }
    
    for atomID in silicons.indices {
      // Elevate the silicon atom by 0.33 nm in the Z direction.
      silicons[atomID].position.z += 3.3 / 10
      
      // Rotate the silicon atom by 10.9° about the origin.
      silicons[atomID].position = rotation.act(on: silicons[atomID].position)
    }
    
    return silicons
  }
}

// MARK: - Expected XYZ Files

private let expectedStep1: String = """
54
Step 1
C  6.306 0.728 0.000
C  5.044 1.456 0.515
C  3.783 0.728 0.000
C  1.261 0.728 0.000
C  2.522 1.456 0.515
C  2.522 2.912 0.000
C  1.261 3.641 0.515
C  0.000 2.912 0.000
C  0.000 1.456 0.515
C  8.828 0.728 0.000
C  10.089 1.456 0.515
C  10.089 2.912 0.000
C  8.828 3.641 0.515
C  7.567 2.912 0.000
C  7.567 1.456 0.515
C  5.044 2.912 0.000
C  6.306 3.641 0.515
C  6.306 5.097 0.000
C  5.044 5.825 0.515
C  3.783 5.097 0.000
C  3.783 3.641 0.515
C  1.261 5.097 0.000
C  2.522 5.825 0.515
C  2.522 7.281 0.000
C  1.261 8.009 0.515
C  0.000 7.281 0.000
C  0.000 5.825 0.515
C  8.828 5.097 0.000
C  10.089 5.825 0.515
C  10.089 7.281 0.000
C  8.828 8.009 0.515
C  7.567 7.281 0.000
C  7.567 5.825 0.515
C  5.044 7.281 0.000
C  6.306 8.009 0.515
C  6.306 9.465 0.000
C  5.044 10.194 0.515
C  3.783 9.465 0.000
C  3.783 8.009 0.515
C  1.261 9.465 0.000
C  2.522 10.194 0.515
C  2.522 11.650 0.000
C  1.261 12.378 0.515
C  0.000 11.650 0.000
C  0.000 10.194 0.515
C  8.828 9.465 0.000
C  10.089 10.194 0.515
C  10.089 11.650 0.000
C  8.828 12.378 0.515
C  7.567 11.650 0.000
C  7.567 10.194 0.515
C  5.044 11.650 0.000
C  6.306 12.378 0.515
C  3.783 12.378 0.515

"""

private let expectedStep2: String = """
54
Step 2
C  6.125 0.707 0.000
C  4.900 1.415 0.000
C  3.675 0.707 0.000
C  1.225 0.707 0.000
C  2.450 1.415 0.000
C  2.450 2.829 0.000
C  1.225 3.536 0.000
C  0.000 2.829 0.000
C  0.000 1.415 0.000
C  8.575 0.707 0.000
C  9.800 1.415 0.000
C  9.800 2.829 0.000
C  8.575 3.536 0.000
C  7.350 2.829 0.000
C  7.350 1.415 0.000
C  4.900 2.829 0.000
C  6.125 3.536 0.000
C  6.125 4.951 0.000
C  4.900 5.658 0.000
C  3.675 4.951 0.000
C  3.675 3.536 0.000
C  1.225 4.951 0.000
C  2.450 5.658 0.000
C  2.450 7.073 0.000
C  1.225 7.780 0.000
C  0.000 7.073 0.000
C  0.000 5.658 0.000
C  8.575 4.951 0.000
C  9.800 5.658 0.000
C  9.800 7.073 0.000
C  8.575 7.780 0.000
C  7.350 7.073 0.000
C  7.350 5.658 0.000
C  4.900 7.073 0.000
C  6.125 7.780 0.000
C  6.125 9.194 0.000
C  4.900 9.902 0.000
C  3.675 9.194 0.000
C  3.675 7.780 0.000
C  1.225 9.194 0.000
C  2.450 9.902 0.000
C  2.450 11.316 0.000
C  1.225 12.023 0.000
C  0.000 11.316 0.000
C  0.000 9.902 0.000
C  8.575 9.194 0.000
C  9.800 9.902 0.000
C  9.800 11.316 0.000
C  8.575 12.023 0.000
C  7.350 11.316 0.000
C  7.350 9.902 0.000
C  4.900 11.316 0.000
C  6.125 12.023 0.000
C  3.675 12.023 0.000

"""

private let expectedStep3: String = """
28
Step 3
Si 9.601 1.109 0.000
Si 7.681 2.217 0.784
Si 5.760 1.109 0.000
Si 1.920 1.109 0.000
Si 3.840 2.217 0.784
Si 3.840 4.434 0.000
Si 1.920 5.543 0.784
Si 0.000 4.434 0.000
Si 0.000 2.217 0.784
Si 11.521 4.434 0.000
Si 11.521 2.217 0.784
Si 7.681 4.434 0.000
Si 9.601 5.543 0.784
Si 9.601 7.760 0.000
Si 7.681 8.869 0.784
Si 5.760 7.760 0.000
Si 5.760 5.543 0.784
Si 1.920 7.760 0.000
Si 3.840 8.869 0.784
Si 3.840 11.086 0.000
Si 1.920 12.195 0.784
Si 0.000 11.086 0.000
Si 0.000 8.869 0.784
Si 11.521 11.086 0.000
Si 11.521 8.869 0.784
Si 7.681 11.086 0.000
Si 9.601 12.195 0.784
Si 5.760 12.195 0.784

"""

private let expectedStep4: String = """
28
Step 4
Si 9.375 1.083 0.000
Si 7.500 2.165 0.620
Si 5.625 1.083 0.000
Si 1.875 1.083 0.000
Si 3.750 2.165 0.620
Si 3.750 4.330 0.000
Si 1.875 5.413 0.620
Si 0.000 4.330 0.000
Si 0.000 2.165 0.620
Si 11.250 4.330 0.000
Si 11.250 2.165 0.620
Si 7.500 4.330 0.000
Si 9.375 5.413 0.620
Si 9.375 7.578 0.000
Si 7.500 8.660 0.620
Si 5.625 7.578 0.000
Si 5.625 5.413 0.620
Si 1.875 7.578 0.000
Si 3.750 8.660 0.620
Si 3.750 10.825 0.000
Si 1.875 11.908 0.620
Si 0.000 10.825 0.000
Si 0.000 8.660 0.620
Si 11.250 10.825 0.000
Si 11.250 8.660 0.620
Si 7.500 10.825 0.000
Si 9.375 11.908 0.620
Si 5.625 11.908 0.620

"""

private let expectedStep5: String = """
82
Step 5
C  6.125 0.707 0.000
C  4.900 1.415 0.000
C  3.675 0.707 0.000
C  1.225 0.707 0.000
C  2.450 1.415 0.000
C  2.450 2.829 0.000
C  1.225 3.536 0.000
C  0.000 2.829 0.000
C  0.000 1.415 0.000
C  8.575 0.707 0.000
C  9.800 1.415 0.000
C  9.800 2.829 0.000
C  8.575 3.536 0.000
C  7.350 2.829 0.000
C  7.350 1.415 0.000
C  4.900 2.829 0.000
C  6.125 3.536 0.000
C  6.125 4.951 0.000
C  4.900 5.658 0.000
C  3.675 4.951 0.000
C  3.675 3.536 0.000
C  1.225 4.951 0.000
C  2.450 5.658 0.000
C  2.450 7.073 0.000
C  1.225 7.780 0.000
C  0.000 7.073 0.000
C  0.000 5.658 0.000
C  8.575 4.951 0.000
C  9.800 5.658 0.000
C  9.800 7.073 0.000
C  8.575 7.780 0.000
C  7.350 7.073 0.000
C  7.350 5.658 0.000
C  4.900 7.073 0.000
C  6.125 7.780 0.000
C  6.125 9.194 0.000
C  4.900 9.902 0.000
C  3.675 9.194 0.000
C  3.675 7.780 0.000
C  1.225 9.194 0.000
C  2.450 9.902 0.000
C  2.450 11.316 0.000
C  1.225 12.023 0.000
C  0.000 11.316 0.000
C  0.000 9.902 0.000
C  8.575 9.194 0.000
C  9.800 9.902 0.000
C  9.800 11.316 0.000
C  8.575 12.023 0.000
C  7.350 11.316 0.000
C  7.350 9.902 0.000
C  4.900 11.316 0.000
C  6.125 12.023 0.000
C  3.675 12.023 0.000
Si 9.001 2.836 3.300
Si 6.955 3.544 3.920
Si 5.319 2.127 3.300
Si 1.636 1.418 3.300
Si 3.273 2.835 3.920
Si 2.864 4.961 3.300
Si 0.818 5.670 3.920
Si -0.819 4.252 3.300
Si -0.409 2.126 3.920
Si 10.228 6.379 3.300
Si 10.638 4.253 3.920
Si 6.546 5.670 3.300
Si 8.182 7.088 3.920
Si 7.773 9.214 3.300
Si 5.727 9.922 3.920
Si 4.091 8.505 3.300
Si 4.500 6.379 3.920
Si 0.408 7.796 3.300
Si 2.045 9.213 3.920
Si 1.635 11.339 3.300
Si -0.411 12.048 3.920
Si -2.047 10.630 3.300
Si -1.638 8.504 3.920
Si 9.000 12.757 3.300
Si 9.409 10.631 3.920
Si 5.318 12.048 3.300
Si 6.954 13.466 3.920
Si 3.272 12.757 3.920

"""