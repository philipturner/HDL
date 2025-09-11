//
//  GrapheneThiol.swift
//  HDLTests
//
//  Created by Philip Turner on 12/28/23.
//

import HDL
import QuaternionModule

struct GrapheneThiol {
  var topology: Topology = .init()
  
  mutating func compilationPass0() {
    let lattice = Lattice<Hexagonal> { h, k, l in
      let h2k = h + 2 * k
      Bounds { 4 * h + 3 * h2k + 1 * l }
      Material { .elemental(.carbon) }
      
      Volume {
        Origin { 0.25 * l }
        Plane { l }
        Replace { .empty }
      }
    }
    topology.atoms += lattice.atoms
  }
  
  mutating func compilationPass1() {
    var grapheneHexagonScale: Float
    
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
    
    for atomID in topology.atoms.indices {
      // Flatten the sp3 sheet into an sp2 sheet.
      topology.atoms[atomID].position.z = 0
      
      // Resize the hexagon side length, so it matches graphene.
      topology.atoms[atomID].position.x *= grapheneHexagonScale
      topology.atoms[atomID].position.y *= grapheneHexagonScale
    }
  }
  
  mutating func compilationPass2() {
    // Graphene's covalent bond length is 1.42 Å.
    let covalentBondLength: Float = 1.42 / 10
    let matches = topology.match(
      topology.atoms, algorithm: .absoluteRadius(covalentBondLength * 1.01))
    
    var insertedBonds: [SIMD2<UInt32>] = []
    for i in topology.atoms.indices {
      for j in matches[i] where i < j {
        let bond = SIMD2(UInt32(i), UInt32(j))
        insertedBonds.append(bond)
      }
    }
    topology.bonds += insertedBonds
  }
  
  mutating func compilationPass3() {
    // From the MM4 alkene paper, the C-H bond length in graphene is 1.103 Å.
    // This is very close to the length of bonds to sp3 carbon.
    //
    // Assume the C-S bond length is the same as in sp3 carbon. From the MM4
    // parameters, the length is 1.814 Å. The S-H bond length is 1.332 Å. Of
    // course, there are electronegativity corrections to bond length, but we
    // just need a number to cite from the literature.
    let chBondLength: Float = 1.103 / 10
    let csBondLength: Float = 1.814 / 10
    let shBondLength: Float = 1.332 / 10
    let orbitalLists = topology.nonbondingOrbitals(hybridization: .sp2)
    
    func createCenterOfMass() -> SIMD3<Float> {
      var output: SIMD3<Double> = .zero
      for atom in topology.atoms {
        output += SIMD3(atom.position)
      }
      output /= Double(topology.atoms.count)
      return SIMD3<Float>(output)
    }
    let centerOfMass = createCenterOfMass()
    
    let thiolRotation = Quaternion<Float>(
      angle: 109.5 * .pi / 180, axis: [0, 0, 1])
    
    var insertedAtoms: [Atom] = []
    var insertedBonds: [SIMD2<UInt32>] = []
    for i in topology.atoms.indices {
      let atom = topology.atoms[i]
      let atomDelta = atom.position - centerOfMass
      
      let orbitalList = orbitalLists[i]
      for orbital in orbitalList {
        if abs(atomDelta.x) > 0.4 && abs(atomDelta.y) > 0.4 {
          let sulfurID = topology.atoms.count + insertedAtoms.count
          insertedBonds.append(SIMD2(UInt32(i), UInt32(sulfurID)))
          
          let sulfurPosition = atom.position + orbital * csBondLength
          let sulfur = Atom(
            position: sulfurPosition,
            element: .sulfur)
          insertedAtoms.append(sulfur)
          
          var thiolOrbital = thiolRotation.act(on: -orbital)
          if atomDelta.x * atomDelta.y > 0 {
            // Mess with the orientation to see whether energy minimization
            // forces them all to the same configuration.
            let rotation = Quaternion<Float>(angle: -.pi / 2, axis: orbital)
            thiolOrbital = rotation.act(on: thiolOrbital)
          }
          let hydrogenID = topology.atoms.count + insertedAtoms.count
          insertedBonds.append(SIMD2(UInt32(i), UInt32(hydrogenID)))
          
          let hydrogenPosition = sulfur.position + thiolOrbital * shBondLength
          let hydrogen = Atom(
            position: hydrogenPosition,
            element: .hydrogen)
          insertedAtoms.append(hydrogen)
        } else {
          let hydrogenID = topology.atoms.count + insertedAtoms.count
          let bond = SIMD2(UInt32(i), UInt32(hydrogenID))
          insertedBonds.append(bond)
          
          let hydrogenPosition = atom.position + orbital * chBondLength
          let hydrogen = Atom(
            position: hydrogenPosition,
            element: .hydrogen)
          insertedAtoms.append(hydrogen)
        }
      }
    }
    topology.atoms += insertedAtoms
    topology.bonds += insertedBonds
  }
}
