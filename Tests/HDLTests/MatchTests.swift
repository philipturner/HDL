import XCTest
import HDL
import Numerics

// This file contains a performance test for topology.match(), which differs
// significantly from the other performance tests. It was moved into a separate
// file for better organization.

final class MatchTests: XCTestCase {
  // We need to run performance tests of Topology.match, to ensure the
  // acceleration algorithm is working properly. One could imagine subtle bugs
  // that make it incorrect, resulting in O(n^2) scaling.
  
  func testMatch() throws {
    // Compute cost scales with the sixth power of lattice width.
    #if RELEASE
    let latticeSizes: [Float] = [1, 1, 2, 3, 4, 5, 6, 7, 8, 10, 16, 24]
    #else
    let latticeSizes: [Float] = [2, 3]
    #endif
    
    for latticeSize in latticeSizes {
      let lattice = Lattice<Cubic> { h, k, l in
        Bounds { latticeSize * (h + k + l) }
        Material { .elemental(.carbon) }
      }
      
      var topology = Topology()
      topology.atoms = lattice.atoms
      
      // Match carbons against carbons.
      do {
        let matches = topology.match(topology.atoms)
        
        var removedAtoms: [UInt32] = []
        var insertedBonds: [SIMD2<UInt32>] = []
        for i in topology.atoms.indices {
          let matchRange = matches[i]
          if matchRange.count <= 2 {
            removedAtoms.append(UInt32(i))
          }
          for j in matchRange where i < j {
            insertedBonds.append(SIMD2(UInt32(i), UInt32(j)))
          }
        }
        topology.insert(bonds: insertedBonds)
        topology.remove(atoms: removedAtoms)
      }
      
      // Define the hydrogen placement distance and the search radius for a
      // subsequent step.
      var ccBondLength = Constant(.square) { .elemental(.carbon) }
      ccBondLength *= Float(3).squareRoot() / 4
      
      // Generate hydrogens and de-duplicate them.
      var hydrogenTopology = Topology()
      do {
        let orbitals = topology.nonbondingOrbitals()
        
        var insertedAtoms: [Atom] = []
        for i in topology.atoms.indices {
          let atom = topology.atoms[i]
          for orbital in orbitals[i] {
            let position = atom.position + orbital * ccBondLength
            let hydrogen = Atom(position: position, element: .hydrogen)
            insertedAtoms.append(hydrogen)
          }
        }
        hydrogenTopology.insert(atoms: insertedAtoms)
        
        let matches = hydrogenTopology.match(
          hydrogenTopology.atoms, algorithm: .absoluteRadius(0.020))
        
        var removedAtoms: [UInt32] = []
        for i in hydrogenTopology.atoms.indices {
          let matchRange = matches[i]
          var survived = true
          for match in matchRange where match < i {
            survived = false
          }
          if !survived {
            removedAtoms.append(UInt32(i))
          }
        }
        hydrogenTopology.remove(atoms: removedAtoms)
      }
      
      // Locate the carbons attached to each hydrogen collision site.
      do {
        _ = topology.match(
          hydrogenTopology.atoms,
          algorithm: .absoluteRadius(1.1 * ccBondLength))
      }
    }
  }
  
  #if RELEASE
  func testFlatSheetScaling() {
    let sheetSizes: [Float] = [100, 200, 300]
    
    for sheetSize in sheetSizes {
      let lattice = Lattice<Hexagonal> { h, k, l in
        let h2k = h + 2 * k
        Bounds { sheetSize * h + 100 * h2k + 1 * l }
        Material { .elemental(.carbon) }
      }
      
      do {
        var topology = Topology()
        topology.insert(atoms: lattice.atoms)
        topology.sort()
      }
      
      do {
        var topology = Topology()
        topology.insert(atoms: lattice.atoms)
        let matches = topology.match(topology.atoms)
        
        var averageMatchCount: Double = 0
        var maxMatchCount: Double = 0
        for matchRange in matches {
          averageMatchCount += Double(matchRange.count)
          maxMatchCount = max(maxMatchCount, Double(matchRange.count))
        }
        averageMatchCount /= Double(topology.atoms.count)
        XCTAssertGreaterThan(averageMatchCount, 4)
        XCTAssertEqual(maxMatchCount, 5)
      }
    }
  }
  #endif
}
