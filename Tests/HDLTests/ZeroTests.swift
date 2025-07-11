import XCTest
import HDL

// Unit test how various APIs respond when they receive arrays of zero length.

final class ZeroTests: XCTestCase {
  static func createDefaultAtoms() -> [Atom] {
    return [
      Atom(position: SIMD3(0.000, 0.000, 0.000), element: .carbon),
      Atom(position: SIMD3(1.000, 0.000, 0.000), element: .carbon),
      Atom(position: SIMD3(0.000, 1.000, 0.000), element: .carbon),
      Atom(position: SIMD3(0.000, 0.000, 1.000), element: .carbon),
    ]
  }
  
  func testLatticeBounds() throws {
    // Test cubic lattices.
    do {
      let lattice = Lattice<Cubic> { h, k, l in
        Bounds { 1 * h + 1 * k + 1 * l }
        Material { .checkerboard(.silicon, .carbon) }
      }
      XCTAssertEqual(lattice.atoms.count, 18)
    }
    do {
      let lattice = Lattice<Cubic> { h, k, l in
        Bounds { 1 * h + 1 * k + 0 * l }
        Material { .checkerboard(.silicon, .carbon) }
      }
      XCTAssertEqual(lattice.atoms.count, 5)
    }
    do {
      let lattice = Lattice<Cubic> { h, k, l in
        Bounds { 1 * h + 0 * k + 0 * l }
        Material { .checkerboard(.silicon, .carbon) }
      }
      XCTAssertEqual(lattice.atoms.count, 2)
    }
    do {
      let lattice = Lattice<Cubic> { h, k, l in
        Bounds { 0 * h + 0 * k + 0 * l }
        Material { .checkerboard(.silicon, .carbon) }
      }
      XCTAssertEqual(lattice.atoms.count, 1)
    }
    
    // Test hexagonal lattices.
    do {
      let lattice = Lattice<Hexagonal> { h, k, l in
        let h2k = h + 2 * k
        Bounds { 1 * h + 1 * h2k + 1 * l }
        Material { .checkerboard(.silicon, .carbon) }
      }
      XCTAssertEqual(lattice.atoms.count, 12)
    }
    do {
      let lattice = Lattice<Hexagonal> { h, k, l in
        let h2k = h + 2 * k
        Bounds { 0 * h + 1 * h2k + 1 * l }
        Material { .checkerboard(.silicon, .carbon) }
      }
      XCTAssertEqual(lattice.atoms.count, 4)
    }
    do {
      let lattice = Lattice<Hexagonal> { h, k, l in
        let h2k = h + 2 * k
        Bounds { 1 * h + 0 * h2k + 1 * l }
        Material { .checkerboard(.silicon, .carbon) }
      }
      XCTAssertEqual(lattice.atoms.count, 0)
    }
    do {
      let lattice = Lattice<Hexagonal> { h, k, l in
        let h2k = h + 2 * k
        Bounds { 1 * h + 1 * h2k + 0 * l }
        Material { .checkerboard(.silicon, .carbon) }
      }
      XCTAssertEqual(lattice.atoms.count, 0)
    }
    do {
      let lattice = Lattice<Hexagonal> { h, k, l in
        let h2k = h + 2 * k
        Bounds { 1 * h + 0 * h2k + 0 * l }
        Material { .checkerboard(.silicon, .carbon) }
      }
      XCTAssertEqual(lattice.atoms.count, 0)
    }
    do {
      let lattice = Lattice<Hexagonal> { h, k, l in
        let h2k = h + 2 * k
        Bounds { 0 * h + 1 * h2k + 0 * l }
        Material { .checkerboard(.silicon, .carbon) }
      }
      XCTAssertEqual(lattice.atoms.count, 0)
    }
    // crasher:
    //   0 * h + 0 * h2k + 1 * l
    // explicitly forbidden right now:
    //   0 * h + 0 * h2k + 0 * l
  }
  
  func testLatticeBoundsPlane() throws {
    // Test cubic lattices.
    do {
      // 1, 1, 1
      let lattice = Lattice<Cubic> { h, k, l in
        Bounds { 1 * h + 1 * k + 1 * l }
        Material { .checkerboard(.silicon, .carbon) }
        
        Volume {
          Origin { 0.5 * h + 0.5 * k + 0.5 * l }
          Plane { h + k + l }
          Replace { .empty }
        }
      }
      XCTAssertEqual(lattice.atoms.count, 8)
    }
    do {
      // 1, 1, 0
      let lattice = Lattice<Cubic> { h, k, l in
        Bounds { 1 * h + 1 * k + 0 * l }
        Material { .checkerboard(.silicon, .carbon) }
        
        Volume {
          Origin { 0.5 * h + 0.5 * k + 0.5 * l }
          Plane { h + k + l }
          Replace { .empty }
        }
      }
      XCTAssertEqual(lattice.atoms.count, 4)
    }
    do {
      // 1, 0, 0
      let lattice = Lattice<Cubic> { h, k, l in
        Bounds { 1 * h + 0 * k + 0 * l }
        Material { .checkerboard(.silicon, .carbon) }
        
        Volume {
          Origin { 0.5 * h + 0.5 * k + 0.5 * l }
          Plane { h + k + l }
          Replace { .empty }
        }
      }
      XCTAssertEqual(lattice.atoms.count, 2)
    }
    do {
      // 0, 0, 0
      let lattice = Lattice<Cubic> { h, k, l in
        Bounds { 0 * h + 0 * k + 0 * l }
        Material { .checkerboard(.silicon, .carbon) }
        
        Volume {
          Origin { 0.5 * h + 0.5 * k + 0.5 * l }
          Plane { h + k + l }
          Replace { .empty }
        }
      }
      XCTAssertEqual(lattice.atoms.count, 1)
    }
    
    // Test hexagonal lattices.
    do {
      // 1, 1, 1
      let lattice = Lattice<Hexagonal> { h, k, l in
        let h2k = h + 2 * k
        Bounds { 1 * h + 1 * h2k + 1 * l }
        Material { .checkerboard(.silicon, .carbon) }
        
        Volume {
          Origin { 0.5 * h + 0.5 * h2k + 0.5 * l }
          Plane { h + h2k + l }
          Replace { .empty }
        }
      }
      XCTAssertEqual(lattice.atoms.count, 7)
    }
    do {
      // 0, 1, 1
      let lattice = Lattice<Hexagonal> { h, k, l in
        let h2k = h + 2 * k
        Bounds { 0 * h + 1 * h2k + 1 * l }
        Material { .checkerboard(.silicon, .carbon) }
        
        Volume {
          Origin { 0.5 * h + 0.5 * h2k + 0.5 * l }
          Plane { h + h2k + l }
          Replace { .empty }
        }
      }
      XCTAssertEqual(lattice.atoms.count, 3)
    }
    do {
      // 1, 0, 1
      let lattice = Lattice<Hexagonal> { h, k, l in
        let h2k = h + 2 * k
        Bounds { 1 * h + 0 * h2k + 1 * l }
        Material { .checkerboard(.silicon, .carbon) }
        
        Volume {
          Origin { 0.5 * h + 0.5 * h2k + 0.5 * l }
          Plane { h + h2k + l }
          Replace { .empty }
        }
      }
      XCTAssertEqual(lattice.atoms.count, 0)
    }
    do {
      // 1, 1, 0
      let lattice = Lattice<Hexagonal> { h, k, l in
        let h2k = h + 2 * k
        Bounds { 1 * h + 1 * h2k + 0 * l }
        Material { .checkerboard(.silicon, .carbon) }
        
        Volume {
          Origin { 0.5 * h + 0.5 * h2k + 0.5 * l }
          Plane { h + h2k + l }
          Replace { .empty }
        }
      }
      XCTAssertEqual(lattice.atoms.count, 0)
    }
    do {
      // 1, 0, 0
      let lattice = Lattice<Hexagonal> { h, k, l in
        let h2k = h + 2 * k
        Bounds { 1 * h + 0 * h2k + 0 * l }
        Material { .checkerboard(.silicon, .carbon) }
        
        Volume {
          Origin { 0.5 * h + 0.5 * h2k + 0.5 * l }
          Plane { h + h2k + l }
          Replace { .empty }
        }
      }
      XCTAssertEqual(lattice.atoms.count, 0)
    }
    do {
      // 0, 1, 0
      let lattice = Lattice<Hexagonal> { h, k, l in
        let h2k = h + 2 * k
        Bounds { 0 * h + 1 * h2k + 0 * l }
        Material { .checkerboard(.silicon, .carbon) }
        
        Volume {
          Origin { 0.5 * h + 0.5 * h2k + 0.5 * l }
          Plane { h + h2k + l }
          Replace { .empty }
        }
      }
      XCTAssertEqual(lattice.atoms.count, 0)
    }
    // crasher:
    //   0 * h + 0 * h2k + 1 * l
    // explicitly forbidden right now:
    //   0 * h + 0 * h2k + 0 * l
  }
  
  func testTopologyMap() throws {
    // finite atoms, no bonds
    do {
      var topology = Topology()
      topology.atoms = Self.createDefaultAtoms()
      
      // .atoms, .atoms
      do {
        let map = topology.map(.atoms, to: .atoms)
        XCTAssertEqual(map.count, 4)
        for atomID in 0..<4 {
          XCTAssertEqual(map[atomID].count, 0)
        }
      }
      
      // .atoms, .bonds
      do {
        let map = topology.map(.atoms, to: .bonds)
        XCTAssertEqual(map.count, 4)
        for atomID in 0..<4 {
          XCTAssertEqual(map[atomID].count, 0)
        }
      }
    }
    
    // no atoms, no bonds
    do {
      let topology = Topology()
      
      // .atoms, .atoms
      do {
        let map = topology.map(.atoms, to: .atoms)
        XCTAssertEqual(map.count, 0)
      }
      
      // .atoms, .bonds
      do {
        let map = topology.map(.atoms, to: .bonds)
        XCTAssertEqual(map.count, 0)
      }
    }
  }
  
  func testTopologyMatch() throws {
    // finite topology.atoms, zero input
    do {
      var topology = Topology()
      topology.atoms = Self.createDefaultAtoms()
      
      let matches = topology.match([])
      XCTAssertEqual(matches.count, 0)
    }
    
    // zero topology.atoms, finite input
    do {
      let topology = Topology()
      
      let matches = topology.match(Self.createDefaultAtoms())
      XCTAssertEqual(matches.count, 4)
      for atomID in 0..<4 {
        XCTAssertEqual(matches[atomID].count, 0)
      }
    }
    
    // zero topology.atoms, zero input
    do {
      let topology = Topology()
      
      let matches = topology.match([])
      XCTAssertEqual(matches.count, 0)
    }
  }
  
  func testTopologyNonbondingOrbitals() throws {
    let hybridizations: [Topology.OrbitalHybridization] = [
      .sp, .sp2, .sp3
    ]
    
    for hybridization in hybridizations {
      // finite atoms, no bonds
      do {
        var topology = Topology()
        topology.atoms = Self.createDefaultAtoms()
        
        let orbitalLists = topology.nonbondingOrbitals(
          hybridization: hybridization)
        XCTAssertEqual(orbitalLists.count, 4)
        for atomID in 0..<4 {
          XCTAssertEqual(orbitalLists[atomID].count, 0)
        }
      }
      
      // no atoms, no bonds
      do {
        let topology = Topology()
        
        let orbitalLists = topology.nonbondingOrbitals(
          hybridization: hybridization)
        XCTAssertEqual(orbitalLists.count, 0)
      }
    }
  }
  
  func testTopologyRemove() throws {
    // finite atoms, finite bonds, no indices
    do {
      var topology = Topology()
      topology.atoms = Self.createDefaultAtoms()
      topology.bonds = [
        SIMD2(0, 1),
        SIMD2(1, 2),
        SIMD2(2, 3),
      ]
      XCTAssertEqual(topology.atoms.count, 4)
      XCTAssertEqual(topology.bonds.count, 3)
      
      topology.remove(atoms: [])
      topology.remove(bonds: [])
      XCTAssertEqual(topology.atoms.count, 4)
      XCTAssertEqual(topology.bonds.count, 3)
    }
    
    // finite atoms, no bonds, no indices
    do {
      var topology = Topology()
      topology.atoms = Self.createDefaultAtoms()
      XCTAssertEqual(topology.atoms.count, 4)
      XCTAssertEqual(topology.bonds.count, 0)
      
      topology.remove(atoms: [])
      topology.remove(bonds: [])
      XCTAssertEqual(topology.atoms.count, 4)
      XCTAssertEqual(topology.bonds.count, 0)
    }
    
    // no atoms, no bonds, no indices
    do {
      var topology = Topology()
      XCTAssertEqual(topology.atoms.count, 0)
      XCTAssertEqual(topology.bonds.count, 0)
      
      topology.remove(atoms: [])
      topology.remove(bonds: [])
      XCTAssertEqual(topology.atoms.count, 0)
      XCTAssertEqual(topology.bonds.count, 0)
    }
    
    // finite atoms, no bonds, finite atom indices
    do {
      var topology = Topology()
      topology.atoms = Self.createDefaultAtoms()
      XCTAssertEqual(topology.atoms.count, 4)
      XCTAssertEqual(topology.bonds.count, 0)
      
      topology.remove(atoms: [1, 3])
      XCTAssertEqual(topology.atoms.count, 2)
      XCTAssertEqual(topology.bonds.count, 0)
    }
  }
  
  func testTopologySort() throws {
    // finite atoms, no bonds
    do {
      var topology = Topology()
      topology.atoms = Self.createDefaultAtoms()
      
      let sortedIndices = topology.sort()
      XCTAssertEqual(sortedIndices.count, 4)
    }
    
    // no atoms, no bonds
    do {
      var topology = Topology()
      
      let sortedIndices = topology.sort()
      XCTAssertEqual(sortedIndices.count, 0)
    }
  }
  
  func testReconstruction() throws {
    var reconstruction = Reconstruction()
    reconstruction.atoms = []
    reconstruction.material = .checkerboard(.silicon, .carbon)
    let topology = reconstruction.compile()
    XCTAssertEqual(topology.atoms.count, 0)
    XCTAssertEqual(topology.bonds.count, 0)
  }
}
