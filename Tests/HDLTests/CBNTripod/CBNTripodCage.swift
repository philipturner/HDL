//
//  CBNTripodCage.swift
//  HDLTests
//
//  Created by Philip Turner on 12/29/23.
//

import HDL
import Numerics
import XCTest

struct CBNTripodCage: CBNTripodComponent {
  var topology: Topology
  
  init() {
    self.topology = Topology()
    
    var bondRecord: [SIMD2<UInt8>: Int]
    var expectedRecord: [SIMD2<UInt8>: Int]
    
    compilationPass0()
    XCTAssertEqual(topology.atoms.count, 13)
    XCTAssertEqual(topology.bonds.count, 0)
    
    compilationPass1()
    bondRecord = createBondRecord()
    XCTAssertEqual(topology.atoms.count, 13)
    XCTAssertEqual(topology.bonds.count, 15)
    XCTAssertEqual(bondRecord.keys.count, 2)
    expectedRecord = [
      SIMD2<UInt8>(6, 6): 12,
      SIMD2<UInt8>(6, 32): 3,
    ]
    for (key, value) in expectedRecord {
      XCTAssertEqual(bondRecord[key], value)
    }
    
    compilationPass2()
    bondRecord = createBondRecord()
    XCTAssertEqual(topology.atoms.count, 19)
    XCTAssertEqual(topology.bonds.count, 21)
    XCTAssertEqual(bondRecord.keys.count, 4)
    expectedRecord = [
      SIMD2<UInt8>(1, 6): 3,
      SIMD2<UInt8>(6, 6): 12,
      SIMD2<UInt8>(6, 8): 3,
      SIMD2<UInt8>(6, 32): 3,
    ]
    for (key, value) in expectedRecord {
      XCTAssertEqual(bondRecord[key], value)
    }
    
    compilationPass3()
    bondRecord = createBondRecord()
    XCTAssertEqual(topology.atoms.count, 31)
    XCTAssertEqual(topology.bonds.count, 33)
    XCTAssertEqual(bondRecord.keys.count, 4)
    expectedRecord = [
      SIMD2<UInt8>(1, 6): 15,
      SIMD2<UInt8>(6, 6): 12,
      SIMD2<UInt8>(6, 8): 3,
      SIMD2<UInt8>(6, 32): 3,
    ]
    for (key, value) in expectedRecord {
      XCTAssertEqual(bondRecord[key], value)
    }
    
    compilationPass4()
    bondRecord = createBondRecord()
    XCTAssertEqual(topology.atoms.count, 33)
    XCTAssertEqual(topology.bonds.count, 35)
    XCTAssertEqual(bondRecord.keys.count, 4)
    expectedRecord = [
      SIMD2<UInt8>(1, 6): 15,
      SIMD2<UInt8>(6, 6): 13,
      SIMD2<UInt8>(6, 8): 3,
      SIMD2<UInt8>(6, 32): 4,
    ]
    for (key, value) in expectedRecord {
      XCTAssertEqual(bondRecord[key], value)
    }
    
    compilationPass5()
    XCTAssertEqual(topology.atoms.count, 33)
    XCTAssertEqual(topology.bonds.count, 35)
  }
  
  // Create the initial structure based on the lattice.
  mutating func compilationPass0() {
    let atoms = createLattice()
    topology.atoms += atoms
  }
  
  // Form C-C and Ge-C bonds, then make the structure a little closer to the
  // minimized one.
  mutating func compilationPass1() {
    // MM4 force field:
    // - C-C bond length is 1.5270 Å.
    // - Ge-C bond length is 1.949 Å.
    //
    // Change the Ge-C bond length to something else for more realistic results
    // before minimization.
    let ccBondLength: Float = 1.5270 / 10
    let geCBondLength: Float = 1.850 / 10
    
    let matches = topology.match(
      topology.atoms, algorithm: .absoluteRadius(ccBondLength * 1.1))
    var insertedBonds: [SIMD2<UInt32>] = []
    for i in topology.atoms.indices {
      for j in matches[i] where i < j {
        let bond = SIMD2(UInt32(i), UInt32(j))
        insertedBonds.append(bond)
      }
    }
    topology.bonds += insertedBonds
    
    // Move the germanium atom slightly upward.
    for i in topology.atoms.indices {
      var atom = topology.atoms[i]
      if atom.atomicNumber == 32 {
        atom.position.y += 20 / 1000
      }
      topology.atoms[i] = atom
    }
    
    // Adjust the nearby carbons.
    let atomsToAtomsMap = topology.map(.atoms, to: .atoms)
    for i in topology.atoms.indices {
      var atom = topology.atoms[i]
      let neighbors = atomsToAtomsMap[i]
      guard atom.atomicNumber == 6 && neighbors.count == 2 else {
        continue
      }
      
      var carbonID: Int = -1
      var germaniumID: Int = -1
      for neighbor in neighbors {
        let atom = topology.atoms[Int(neighbor)]
        if atom.atomicNumber == 32 {
          germaniumID = Int(neighbor)
        } else {
          carbonID = Int(neighbor)
        }
      }
      guard germaniumID >= 0 else {
        continue
      }
      XCTAssertGreaterThanOrEqual(carbonID, 0)
      let germanium = topology.atoms[germaniumID]
      let carbon = topology.atoms[carbonID]
      
      // Iterate until both bonds reach equilibrium length.
      for _ in 0..<5 {
        var geCDelta = atom.position - germanium.position
        geCDelta /= (geCDelta * geCDelta).sum().squareRoot()
        atom.position = germanium.position + geCDelta * geCBondLength
        
        var ccDelta = atom.position - carbon.position
        ccDelta /= (ccDelta * ccDelta).sum().squareRoot()
        atom.position = carbon.position + ccDelta * ccBondLength
      }
      topology.atoms[i] = atom
    }
  }
  
  // Create CHO groups at the methyl sites.
  // - replace carbons near the primary carbons with silicon markers
  // - remove the methyl groups
  // - restore once you know each orbital's direction
  mutating func compilationPass2() {
    // MM3 Tinker parameters:
    // - carbonyl sp2 C-H bond length is 1.1180 Å.
    // - carbonyl sp2 C - sp3 C bond ength is 1.5090 Å.
    // - carbonyl sp2 C=O bond length is 1.2080 Å.
    let chBondLength: Float = 1.1180 / 10
    let ccBondLength: Float = 1.5090 / 10
    let coBondLength: Float = 1.2080 / 10
    
    do {
      let atomsToAtomsMap = topology.map(.atoms, to: .atoms)
      
      var removedAtoms: [UInt32] = []
      for i in topology.atoms.indices {
        let neighbors = atomsToAtomsMap[i]
        if neighbors.count == 1 {
          let neighborID = Int(neighbors.first!)
          topology.atoms[neighborID].atomicNumber = 14
          removedAtoms.append(UInt32(i))
        }
      }
      topology.remove(atoms: removedAtoms)
    }
    
    do {
      let orbitalLists = topology.nonbondingOrbitals(hybridization: .sp3)
      
      var insertedAtoms: [Atom] = []
      var insertedBonds: [SIMD2<UInt32>] = []
      for i in topology.atoms.indices {
        let atom = topology.atoms[i]
        guard atom.atomicNumber == 14 else {
          continue
        }
        topology.atoms[i].atomicNumber = 6
        
        let orbital = orbitalLists[i].first!
        let carbonPosition = atom.position + orbital * ccBondLength
        let carbon = Atom(position: carbonPosition, element: .carbon)
        let carbonID = topology.atoms.count + insertedAtoms.count
        insertedAtoms.append(carbon)
        insertedBonds.append(SIMD2(UInt32(i), UInt32(carbonID)))
        
        let ground = SIMD3<Float>(0, -1, 0)
        let groundRotation = Quaternion<Float>(from: -orbital, to: ground)
        let axis = groundRotation.axis
        let rotation = Quaternion(angle: 2 * .pi / 3, axis: axis)
        
        let oxygenOrbital = rotation.act(on: -orbital)
        let oxygenPosition = carbon.position + oxygenOrbital * coBondLength
        let oxygen = Atom(position: oxygenPosition, element: .oxygen)
        let oxygenID = topology.atoms.count + insertedAtoms.count
        insertedAtoms.append(oxygen)
        insertedBonds.append(SIMD2(UInt32(carbonID), UInt32(oxygenID)))
        
        let hydrogenOrbital = rotation.act(on: oxygenOrbital)
        let hydrogenPosition = carbon.position + hydrogenOrbital * chBondLength
        let hydrogen = Atom(
          position: hydrogenPosition, element: .hydrogen)
        let hydrogenID = topology.atoms.count + insertedAtoms.count
        insertedAtoms.append(hydrogen)
        insertedBonds.append(SIMD2(UInt32(carbonID), UInt32(hydrogenID)))
      }
      topology.atoms += insertedAtoms
      topology.bonds += insertedBonds
    }
  }
  
  // Hydrogenate the remainder of the lattice, except the germanium.
  mutating func compilationPass3() {
    // MM4 force field:
    // - C-H bond length is 1.1120 Å.
    let chBondLength: Float = 1.1120 / 10
    let atomsToAtomsMap = topology.map(.atoms, to: .atoms)
    let orbitalLists = topology.nonbondingOrbitals(hybridization: .sp3)
    
    var insertedAtoms: [Atom] = []
    var insertedBonds: [SIMD2<UInt32>] = []
    for i in topology.atoms.indices {
      let atom = topology.atoms[i]
      if atom.atomicNumber == 32 {
        continue
      }
      
      let neighbors = atomsToAtomsMap[i]
      var oxygenID: Int = -1
      for neighbor in neighbors {
        let atom = topology.atoms[Int(neighbor)]
        if atom.atomicNumber == 8 {
          oxygenID = Int(neighbor)
        }
      }
      if oxygenID >= 0 {
        continue
      }
      
      let orbitalList = orbitalLists[i]
      for orbital in orbitalList {
        let position = atom.position + orbital * chBondLength
        let hydrogen = Atom(position: position, element: .hydrogen)
        let hydrogenID = topology.atoms.count + insertedAtoms.count
        insertedAtoms.append(hydrogen)
        insertedBonds.append(SIMD2(UInt32(i), UInt32(hydrogenID)))
      }
    }
    topology.atoms += insertedAtoms
    topology.bonds += insertedBonds
  }
  
  // Add the carbon dimer to the germanium.
  mutating func compilationPass4() {
    // MM3 Tinker parameters:
    // - sp1 C - sp1 C bond length is 1.2100 Å.
    //
    // For the C-Ge bond length, extrapolate from the sp1 C-C and sp2 C-Ge
    // lengths. Take the delta between sp1 C - sp3 C and sp2 C - sp3 C lengths,
    // then add to the sp2 C - sp3 Ge length.
    let ccBondLength: Float = 1.2100 / 10
    let cGeBondLength: Float = (1.4700 - 1.4990 + 1.9350) / 10
    
    do {
      let orbitalLists = topology.nonbondingOrbitals(hybridization: .sp3)
      
      var germaniumID: Int = -1
      for i in topology.atoms.indices {
        if topology.atoms[i].atomicNumber == 32 {
          germaniumID = i
        }
      }
      
      let germanium = topology.atoms[germaniumID]
      let orbital = orbitalLists[germaniumID].first!
      let carbonPosition1 = germanium.position + orbital * cGeBondLength
      let carbon1 = Atom(position: carbonPosition1, element: .carbon)
      let carbonID1 = topology.atoms.count
      topology.atoms.append(carbon1)
      
      let bond = SIMD2(
        UInt32(germaniumID),
        UInt32(carbonID1))
      topology.bonds.append(bond)
    }
    
    // Test out the new sp1 orbital generation functionality by adding the
    // second carbon in a separate pass.
    do {
      let orbitalLists = topology.nonbondingOrbitals(hybridization: .sp)
      
      var carbonID1: Int = -1
      for i in topology.atoms.indices {
        let atom = topology.atoms[i]
        let orbitalList = orbitalLists[i]
        if atom.atomicNumber == 6,
           orbitalList.count == 1 {
          carbonID1 = i
        }
      }
      
      let carbon1 = topology.atoms[carbonID1]
      let orbital = orbitalLists[carbonID1].first!
      let carbonPosition2 = carbon1.position + orbital * ccBondLength
      let carbon2 = Atom(position: carbonPosition2, element: .carbon)
      let carbonID2 = topology.atoms.count
      topology.atoms.append(carbon2)
      
      let bond = SIMD2(
        UInt32(carbonID1),
        UInt32(carbonID2))
      topology.bonds.append(bond)
    }
  }
  
  // Replace the atom positions with the energy-minimized ones from xTB.
  mutating func compilationPass5() {
    let xtbOptimizedAtoms: [Atom] = [
      Atom(position: SIMD3( 0.0000, -0.2471, -0.1449), element: .carbon),
      Atom(position: SIMD3(-0.1255, -0.2471,  0.0725), element: .carbon),
      Atom(position: SIMD3( 0.1255, -0.2471,  0.0725), element: .carbon),
      Atom(position: SIMD3( 0.0000, -0.2024,  0.1490), element: .carbon),
      Atom(position: SIMD3( 0.0000, -0.0523,  0.1782), element: .carbon),
      Atom(position: SIMD3( 0.1290, -0.2024, -0.0745), element: .carbon),
      Atom(position: SIMD3( 0.1543, -0.0523, -0.0891), element: .carbon),
      Atom(position: SIMD3(-0.1290, -0.2024, -0.0745), element: .carbon),
      Atom(position: SIMD3(-0.1543, -0.0523, -0.0891), element: .carbon),
      Atom(position: SIMD3( 0.0000,  0.0341,  0.0000), element: .germanium),
      Atom(position: SIMD3( 0.0000, -0.2795,  0.2808), element: .carbon),
      Atom(position: SIMD3(-0.0000, -0.3987,  0.2894), element: .oxygen),
      Atom(position: SIMD3( 0.0000, -0.2153,  0.3710), element: .hydrogen),
      Atom(position: SIMD3( 0.2431, -0.2795, -0.1404), element: .carbon),
      Atom(position: SIMD3( 0.2506, -0.3987, -0.1447), element: .oxygen),
      Atom(position: SIMD3( 0.3213, -0.2153, -0.1856), element: .hydrogen),
      Atom(position: SIMD3(-0.2431, -0.2795, -0.1404), element: .carbon),
      Atom(position: SIMD3(-0.2506, -0.3987, -0.1447), element: .oxygen),
      Atom(position: SIMD3(-0.3213, -0.2153, -0.1856), element: .hydrogen),
      Atom(position: SIMD3(-0.0000, -0.3563, -0.1495), element: .hydrogen),
      Atom(position: SIMD3( 0.0000, -0.2088, -0.2475), element: .hydrogen),
      Atom(position: SIMD3(-0.2143, -0.2088,  0.1237), element: .hydrogen),
      Atom(position: SIMD3(-0.1295, -0.3563,  0.0748), element: .hydrogen),
      Atom(position: SIMD3( 0.1295, -0.3563,  0.0748), element: .hydrogen),
      Atom(position: SIMD3( 0.2143, -0.2088,  0.1237), element: .hydrogen),
      Atom(position: SIMD3( 0.0880, -0.0254,  0.2368), element: .hydrogen),
      Atom(position: SIMD3(-0.0880, -0.0254,  0.2368), element: .hydrogen),
      Atom(position: SIMD3( 0.1610, -0.0254, -0.1946), element: .hydrogen),
      Atom(position: SIMD3( 0.2490, -0.0254, -0.0422), element: .hydrogen),
      Atom(position: SIMD3(-0.2490, -0.0254, -0.0422), element: .hydrogen),
      Atom(position: SIMD3(-0.1610, -0.0254, -0.1946), element: .hydrogen),
      Atom(position: SIMD3(-0.0000,  0.2273, -0.0000), element: .carbon),
      Atom(position: SIMD3(-0.0000,  0.3471, -0.0000), element: .carbon),
    ]
    
    topology.atoms = xtbOptimizedAtoms
  }
}

extension CBNTripodCage {
  // Create a lattice with a germanium dopant. This will be extremely far from
  // equilibrium, but the compiled structure doesn't need to be close. It only
  // needed to be close for the tripod leg, where replacing with the minimized
  // structure constituted unacceptable information loss.
  func createLattice() -> [Atom] {
    // Create the adamantane cage with Ge and with 3 methyl groups in place of
    // the leg structures.
    let lattice = Lattice<Cubic> { h, k, l in
      Bounds { 4 * h + 4 * k + 4 * l }
      Material { .elemental(.carbon) }
      
      Volume {
        Origin { 2 * h + 2 * k + 2 * l }
        Origin { 0.25 * (h + k - l) }
        
        // Remove the front plane.
        Convex {
          Origin { 0.25 * (h + k + l) }
          Plane { h + k + l }
        }
        
        Volume {
          Convex {
            Origin { 0.2 * (h + k + l) }
            Plane { h + k + l }
          }
          
          Replace { .atom(.germanium) }
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
        
        // Keep the 3 carbons representing the legs.
        // triangleCut(sign: +1)
        
        // Remove the remaining carbons.
        triangleCut(sign: -1)
        
        // Remove the back plane.
        Convex {
          Origin { -0.25 * (h + k + l) }
          Plane { -(h + k + l) }
        }
        
        Replace { .empty }
      }
    }
    var atoms = lattice.atoms
    
    // Rotate the cage so the germanium points straight up, and one of the
    // legs points toward +Z.
    let basisX = SIMD3<Float>(1, 0, -1) / Float(2).squareRoot()
    let basisY = SIMD3<Float>(1, 1, 1) / Float(3).squareRoot()
    XCTAssertLessThan((basisX * basisY).sum().magnitude, 1e-3)
    
    func cross<T: Real & SIMDScalar>(
      _ x: SIMD3<T>, _ y: SIMD3<T>
    ) -> SIMD3<T> {
      // Source: https://en.wikipedia.org/wiki/Cross_product#Computing
      let s1 = x[1] * y[2] - x[2] * y[1]
      let s2 = x[2] * y[0] - x[0] * y[2]
      let s3 = x[0] * y[1] - x[1] * y[0]
      return SIMD3(s1, s2, s3)
    }
    let basisZ = -cross(basisX, basisY)
    let basisZLength = (basisZ * basisZ).sum().squareRoot()
    XCTAssertLessThan((basisZLength - 1).magnitude, 1e-3)
    
    for i in atoms.indices {
      var atom = atoms[i]
      let componentX = (atom.position * basisX).sum()
      let componentY = (atom.position * basisY).sum()
      let componentZ = (atom.position * basisZ).sum()
      atom.position = SIMD3(componentX, componentY, componentZ)
      atoms[i] = atom
    }
    
    var germaniumID: Int = -1
    for i in atoms.indices {
      let atom = atoms[i]
      if atom.atomicNumber == 32 {
        germaniumID = i
      }
    }
    XCTAssertNotEqual(germaniumID, -1)
    
    // Center the cage so the germanium is at (0, 0, 0).
    let germaniumPosition = atoms[germaniumID].position
    let translation = SIMD3<Float>.zero - germaniumPosition
    for i in atoms.indices {
      atoms[i].position += translation
    }
    return atoms
  }
}
