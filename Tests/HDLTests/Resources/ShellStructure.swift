//
//  ShellStructure.swift
//  HDLTests
//
//  Created by Philip Turner on 12/28/23.
//

import HDL
import QuaternionModule
import XCTest

struct ShellStructure {
  var element: Element
  var topology: Topology
  
  init(element: Element) {
    self.element = element
    self.topology = Topology()
  }
  
  mutating func compilationPass0() {
    let lattice = createLonsdaleiteLattice(element: element)
    topology.atoms += lattice.atoms
  }
  
  mutating func compilationPass1() {
    do {
      let matches = topology.match(topology.atoms)
      var insertedBonds: [SIMD2<UInt32>] = []
      for i in topology.atoms.indices {
        let match = matches[i]
        for j in match where i < j {
          insertedBonds.append(
            SIMD2(UInt32(i), UInt32(j)))
        }
      }
      topology.bonds += insertedBonds
    }
    
    let orbitalLists = topology.nonbondingOrbitals()
    
    var insertedAtoms: [Atom] = []
    var insertedBonds: [SIMD2<UInt32>] = []
    for i in topology.atoms.indices {
      let atom = topology.atoms[i]
      let bondLength = Element.hydrogen.covalentRadius +
      Element(rawValue: atom.atomicNumber)!.covalentRadius
      
      let orbitalList = orbitalLists[i]
      for orbital in orbitalList {
        let hydrogenID = topology.atoms.count + insertedAtoms.count
        let position = atom.position + bondLength * orbital
        let hydrogen = Atom(position: position, element: .hydrogen)
        insertedAtoms.append(hydrogen)
        insertedBonds.append(
          SIMD2(UInt32(i), UInt32(hydrogenID)))
      }
    }
    topology.atoms += insertedAtoms
    topology.bonds += insertedBonds
  }
  
  mutating func compilationPass2() {
    let latticeConstant = Constant(.hexagon) { .elemental(element) }
    let perimeter = 3 * 8 * latticeConstant
    let radius = perimeter / (2 * Float.pi)
    
    var averageZ: Double = 0
    for atom in topology.atoms {
      averageZ += Double(atom.position.z)
    }
    averageZ /= Double(topology.atoms.count)
    
    for i in topology.atoms.indices {
      let previous = topology.atoms[i].position
      let theta = previous.x / perimeter * (2 * Float.pi)
      let r = (Float(averageZ) + radius) - previous.z
      
      let vectorI = SIMD3<Float>(0, 0, -1)
      let vectorJ = SIMD3<Float>(1, 0, 0)
      
      var position: SIMD3<Float> = .zero
      position += vectorI * r * Float.cos(theta)
      position += vectorJ * r * Float.sin(theta)
      position.y = previous.y
      topology.atoms[i].position = position
    }
  }
  
  mutating func compilationPass3(onlyMergeHydrogens: Bool) {
    let matches = topology.match(
      topology.atoms, algorithm: .absoluteRadius(0.010))
    let atomsToAtomsMap = topology.map(.atoms, to: .atoms)
    var removedAtoms: Set<UInt32> = []
    var insertedBonds: Set<SIMD2<UInt32>> = []
    
    for i in topology.atoms.indices {
      let atomI = topology.atoms[i]
      guard matches[i].count > 1 else {
        continue
      }
      guard matches[i].count == 2 else {
        fatalError("Too many overlapping atoms.")
      }
      
      var j: Int = -1
      for match in matches[i] where i != match {
        j = Int(match)
      }
      let atomJ = topology.atoms[j]
      XCTAssertEqual(atomI.atomicNumber, atomJ.atomicNumber)
      
      // Choose the carbon with the lowest index, or the H duplicate associated
      // with that carbon.
      let neighborsI = atomsToAtomsMap[i]
      let neighborsJ = atomsToAtomsMap[j]
      XCTAssertEqual(neighborsI.count, neighborsJ.count)
      if atomI.atomicNumber == 1 {
        XCTAssertEqual(neighborsI.count, 1)
        guard neighborsI.first! < neighborsJ.first! else {
          continue
        }
      } else {
        XCTAssertEqual(neighborsI.count, 4)
        guard i < j else {
          continue
        }
      }
      
      if atomI.atomicNumber == 1 {
        removedAtoms.insert(UInt32(j))
        continue
      }
      
      if onlyMergeHydrogens {
        continue
      }
      
      struct Orbital {
        var neighborID: UInt32
        var neighborElement: UInt8
        var delta: SIMD3<Float>
      }
      func createOrbitals(_ index: Int) -> [Orbital] {
        let neighbors = atomsToAtomsMap[index]
        let selfAtom = topology.atoms[index]
        
        var output: [Orbital] = []
        for neighborID in neighbors {
          let otherAtom = topology.atoms[Int(neighborID)]
          var delta = otherAtom.position - selfAtom.position
          delta /= (delta * delta).sum().squareRoot()
          
          let orbital = Orbital(
            neighborID: neighborID,
            neighborElement: otherAtom.atomicNumber,
            delta: delta)
          output.append(orbital)
        }
        return output
      }
      let orbitalsI = createOrbitals(i)
      var orbitalsJ = createOrbitals(j)
      var orbitalJMatches: [Int] = []
      for orbitalJ in orbitalsJ {
        var maxScore: Float = -.greatestFiniteMagnitude
        var maxIndex: Int = -1
        for indexI in 0..<4 {
          let orbitalI = orbitalsI[indexI]
          let score = (orbitalI.delta * orbitalJ.delta).sum()
          if score > maxScore {
            maxScore = score
            maxIndex = indexI
          }
        }
        XCTAssertGreaterThanOrEqual(maxIndex, 0)
        XCTAssertFalse(orbitalJMatches.contains(maxIndex))
        orbitalJMatches.append(maxIndex)
      }
      let nullOrbital = Orbital(
        neighborID: 0, neighborElement: 0, delta: .zero)
      var newOrbitalsJ = Array(repeating: nullOrbital, count: 4)
      for indexJ in 0..<4 {
        let maxIndex = orbitalJMatches[indexJ]
        newOrbitalsJ[maxIndex] = orbitalsJ[indexJ]
      }
      orbitalsJ = newOrbitalsJ
      
      for (orbitalI, orbitalJ) in zip(orbitalsI, orbitalsJ) {
        switch (orbitalI.neighborElement, orbitalJ.neighborElement) {
        case (1, 1):
          // The overlapping hydrogens should already be removed.
          break
        case (6, 6), (14, 14):
          if orbitalI.neighborID < orbitalJ.neighborID {
            // The sigma bond to the other carbon was duplicated, and will be
            // automatically removed.
            break
          } else {
            fatalError("Edge case not handled.")
          }
        case (6, 1), (14, 1):
          // The overlapping hydrogen and carbon are not already removed by
          // other code.
          removedAtoms.insert(orbitalJ.neighborID)
        case (1, 6), (1, 14):
          // The hydrogen from the first atom must be superseded by the carbon
          // from the second atom. That carbon is not registered as overlapping
          // anything, because its position differs from the replaced hydrogen.
          XCTAssertFalse(removedAtoms.contains(orbitalJ.neighborID))
          removedAtoms.insert(orbitalI.neighborID)
          insertedBonds.insert(SIMD2(UInt32(i), orbitalJ.neighborID))
        default:
          fatalError("Unrecognized bond.")
        }
      }
      removedAtoms.insert(UInt32(j))
    }
    
    topology.bonds += Array(insertedBonds)
    topology.remove(atoms: Array(removedAtoms))
  }
}

private func createLonsdaleiteLattice(
  element: Element
) -> Lattice<Hexagonal> {
  Lattice<Hexagonal> { h, k, l in
    let h2k = h + 2 * k
    
    let thickness: Float = 3
    let spacing: Float = 8
    let boundsH: Float = (3+1)*spacing
    Bounds { boundsH * h + 4 * h2k + 1 * l }
    Material { .elemental(element) }
    
    Volume {
      let cuts = 1 + Int((boundsH/spacing).rounded(.up))
      for i in 0..<cuts {
        Concave {
          Origin { 3.5 * h2k }
          Origin { Float(i) * spacing * h }
          Convex {
            Plane { h + k }
          }
          Convex {
            Plane { h2k }
          }
          Convex {
            Origin { thickness * h }
            Plane { k }
          }
        }
      }
      
      Origin { (-spacing/2) * h }
      for i in 0..<(cuts+1) {
        Concave {
          Origin { 0.5 * h2k }
          Origin { Float(i) * spacing * h }
          Convex {
            Plane { -k }
          }
          Convex {
            Plane { -h2k }
          }
          Convex {
            Origin { thickness * h }
            Plane { -k - h }
          }
        }
      }
      
      Replace { .empty }
    }
  }
}
