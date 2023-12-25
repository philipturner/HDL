import XCTest
import HDL

final class TopologyTests: XCTestCase {
  func testInsertRemove() throws {
    let lonsdaleite = Lonsdaleite()
    var topology = Topology()
    
    topology.insert(atoms: lonsdaleite.atoms)
    XCTAssertEqual(topology.atoms, lonsdaleite.atoms)
    XCTAssertEqual(topology.bonds, [])
    
    topology.remove(atoms: lonsdaleite.atoms.indices.map(UInt32.init))
    XCTAssertEqual(topology.atoms, [])
    XCTAssertEqual(topology.bonds, [])
    
    topology.insert(atoms: lonsdaleite.atoms)
    topology.insert(bonds: lonsdaleite.bonds)
    XCTAssertEqual(topology.atoms, lonsdaleite.atoms)
    XCTAssertEqual(topology.bonds, lonsdaleite.bonds)
    
    topology.remove(bonds: lonsdaleite.bonds.indices.map(UInt32.init))
    XCTAssertEqual(topology.atoms, lonsdaleite.atoms)
    XCTAssertEqual(topology.bonds, [])
    
    topology.insert(bonds: lonsdaleite.bonds)
    XCTAssertEqual(topology.atoms, lonsdaleite.atoms)
    XCTAssertEqual(topology.bonds, lonsdaleite.bonds)
    
    topology.remove(atoms: lonsdaleite.atoms.indices.map(UInt32.init))
    XCTAssertEqual(topology.atoms, [])
    XCTAssertEqual(topology.bonds, [])
  }
  
  func testPartialRemove() {
    let lonsdaleite = Lonsdaleite()
    var topology = Topology()
    
    var compactedAtomCount = 0
    var mappedAtomIndices: [Int] = []
    var removedAtomIDs: [UInt32] = []
    var removedBondIDs: [UInt32] = []
    
    var carbonAtoms: [Entity] = []
    var originalCarbonBonds: [SIMD2<UInt32>] = []
    var compactedCarbonBonds: [SIMD2<UInt32>] = []
    
    for atomID in lonsdaleite.atoms.indices {
      let atom = lonsdaleite.atoms[atomID]
      if atom.atomicNumber == 1 {
        mappedAtomIndices.append(-1)
        removedAtomIDs.append(UInt32(atomID))
      } else {
        mappedAtomIndices.append(compactedAtomCount)
        carbonAtoms.append(atom)
        compactedAtomCount += 1
      }
    }
    
    for bondID in lonsdaleite.bonds.indices {
      let bond = lonsdaleite.bonds[bondID]
      let atom1 = lonsdaleite.atoms[Int(bond[0])]
      let atom2 = lonsdaleite.atoms[Int(bond[1])]
      if atom1.atomicNumber == 1 || atom2.atomicNumber == 1 {
        removedBondIDs.append(UInt32(bondID))
        continue
      }
      
      let mappedID1 = mappedAtomIndices[Int(bond[0])]
      let mappedID2 = mappedAtomIndices[Int(bond[1])]
      originalCarbonBonds.append(bond)
      compactedCarbonBonds.append(SIMD2(UInt32(mappedID1),
                                        UInt32(mappedID2)))
    }
    
    // Remove all carbon-hydrogen bonds and ensure the remaining topology
    // appears as expected.
    topology.insert(atoms: lonsdaleite.atoms)
    topology.insert(bonds: lonsdaleite.bonds)
    topology.remove(bonds: removedBondIDs)
    XCTAssertEqual(topology.atoms, lonsdaleite.atoms)
    XCTAssertEqual(topology.bonds, originalCarbonBonds)
    
    topology.remove(atoms: topology.atoms.indices.map(UInt32.init))
    XCTAssertEqual(topology.atoms, [])
    XCTAssertEqual(topology.bonds, [])
    
    // Remove all hydrogen atoms and ensure the remaining carbon-carbon bonds
    // appear as expected.
    topology.insert(atoms: lonsdaleite.atoms)
    topology.insert(bonds: lonsdaleite.bonds)
    topology.remove(atoms: removedAtomIDs)
    XCTAssertEqual(topology.atoms, carbonAtoms)
    XCTAssertEqual(topology.bonds, compactedCarbonBonds)
  }
  
  // This test currently covers Morton reordering, but not correctness of
  // bond remapping afterward. Before testing the latter, we need to visualize
  // results of the atom and bond reordering in the renderer.
  func testSorting() throws {
    /*
     // TODO: Shuffle the indices for atoms in a random order, then re-map the
     // bonds accordingly.
     
     // TODO: In a separate test case, invert the positions of the atoms and
     // ensure the bonds are reordered properly.
     */
  }
  
  // Good ideas for stuff the test suite should eventually cover:
  //
  // Idea for testing correctness of Topology.sort(): Repeat the same test
  // multiple times with the input randomly shuffled beforehand. Assert
  // that the output (atoms, bonds) is always the same, and that the output is
  // different from the input.
  //
  // Idea for testing correctness of Topology.match(): Ensure the matched
  // indices actually appear in ascending order of distance. Ensure none of them
  // are more than the algorithm's specified distance cutoff.
  //
  // Idea for ultimate litmus test: use the new, flexible topology compiler to
  // reconstruct (100) surfaces.
  
  // TODO: - Instead of a time-consuming, exhaustive test suite, debug this
  // visually in the renderer. The covered functionality is not as complex as
  // MM4RigidBody and doesn't need the same approach to testing. Getting (100)
  // reconstruction to work will likely trigger edge cases where the compiler is
  // broken; reproducers can be added to the test suite.
  //
  // TODO: - The second test case is an interesting method of forming strained
  // shell structures. Passivate a crystalline lattice, then warp and remove
  // hydrogens bonded to carbons that will merge. Validate that both this
  // and the reconstructed (100) structure are accepted by the simulator.
  //
  // sort() should still be debugged as explained above; this method is a form
  // of visual debugging. It may be helpful to have a utility function for
  // debugging bonds. For example, explode the crystal lattice and place marker
  // atoms on an interpolated line between actual atoms. Presence of malformed
  // bonds will be extremely obvious.
  //
  // The old 'Diamondoid' topology generator can be used to bootstrap testing
  // of how sort() treats bonds, before the remainder of the functionality is
  // working. In addition, enter into MM4Parameters as a final validation test.
  //
  // Implementation plan:
  // - 1) Visualizer for Morton order and bond topology in GitHub gist.
  //   - 1.1) Test against Lattice -> Diamondoid reordering.
  //   - 1.2) Test against Topology.sort().
  // - 2) Test simple diamond and lonsdaleite lattice formation.
  // - 3) Reproduce graphene thiol and HAbst tripods from the
  //      nanofactory demo.
  // - 4) Demonstrate (100) surface and strained shell structure formation.
}
