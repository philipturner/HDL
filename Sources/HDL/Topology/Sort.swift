//
//  Sort.swift
//  MolecularRenderer
//
//  Created by Philip Turner on 12/2/23.
//

// Notes for when you get around to optimizing this:
//
// A test should compare the grid sorter to an alternative, recursive
// implementation. It first measures cold-start speed, then speed once the
// grid sorter can take advantage of the set already being sorted.
//
// The eventually chosen implementation might be a hybrid of these. It might
// temporarily generate Morton indices to check whether segments of the atoms
// are already sorted. It could also switch between different methods at
// different levels of the hierarchy.
//
// Note: random shuffling creates a different data distribution than the
// typical distribution from Lattice. It may unfairly increase execution time
// of the grid sorter.

import System

extension Topology {
  @discardableResult
  public mutating func sort() -> [UInt32] {
    let grid = GridSorter(atoms: atoms)
    let reordering = grid.mortonReordering()
    let previousAtoms = atoms
    for originalID in reordering.indices {
      let reorderedID32 = reordering[originalID]
      let reorderedID = Int(truncatingIfNeeded: reorderedID32)
      atoms[reorderedID] = previousAtoms[originalID]
    }
    
    for i in bonds.indices {
      let bond = SIMD2<Int>(truncatingIfNeeded: bonds[i])
      var newBond: SIMD2<UInt32> = .zero
      newBond[0] = reordering[bond[0]]
      newBond[1] = reordering[bond[1]]
      newBond = SIMD2(newBond.min(), newBond.max())
      bonds[i] = newBond
    }
    bonds.sort {
      if $0.x != $1.x {
        return $0.x < $1.x
      } else {
        return $0.y < $1.y
      }
    }
    return reordering
  }
}

private let startTime = ContinuousClock.now

func cross_platform_media_time() -> Double {
  let duration = ContinuousClock.now.duration(to: startTime)
  let seconds = duration.components.seconds
  let attoseconds = duration.components.attoseconds
  return -(Double(seconds) + Double(attoseconds) * 1e-18)
}

struct GridSorter {
  var atoms: [Entity] = []
  
  var origin: SIMD3<Float>
  var dimensions: SIMD3<Float>
  
  init(atoms: [Entity]) {
    self.atoms = atoms
    
    if atoms.count == 0 {
      origin = .zero
      dimensions = .init(repeating: 1)
    } else {
      var minimum = SIMD3<Float>(repeating: .greatestFiniteMagnitude)
      var maximum = -minimum
      for atom in atoms {
        let position = atom.position
        minimum.replace(with: position, where: position .< minimum)
        maximum.replace(with: position, where: position .> maximum)
      }
      for atom in atoms {
        let position = atom.position
        minimum.replace(with: position, where: position .< minimum)
      }
      
      origin = minimum
      dimensions = maximum - minimum
      dimensions.replace(with: .init(repeating: 0.5), where: dimensions .< 0.5)
    }
  }
}

// MARK: - Morton Reordering

extension GridSorter {
  func mortonReordering() -> [UInt32] {
    var output: [UInt32] = []
    let dictionary: UnsafeMutablePointer<UInt32> =
      .allocate(capacity: 8 * atoms.count)
    defer {
      dictionary.deallocate()
    }
    
    func traverse(
      atomIDs: UnsafeBufferPointer<UInt32>,
      levelOrigin: SIMD3<Float>,
      levelSize: Float
    ) {
      if levelSize < 1 / 32 {
        output += atomIDs
        return
      }
      
      var dictionaryCount: SIMD8<Int> = .zero
      
      for atomID32 in atomIDs {
        var index: SIMD3<UInt32> = .init(repeating: 1)
        let atomID = Int(truncatingIfNeeded: atomID32)
        let atomPosition = atoms[atomID].position - self.origin
        index.replace(
          with: .init(repeating: 0),
          where: atomPosition .< levelOrigin)
        
        let key = Int(
          truncatingIfNeeded: (index &<< SIMD3(0, 1, 2)).wrappedSum())
        let previousCount = dictionaryCount[key]
        let pointer = dictionary.advanced(
          by: key &* atoms.count &+ previousCount)
        
        dictionaryCount[key] &+= 1
        pointer.pointee = atomID32
      }
      
      var temporaryAllocationSize = 0
      for laneID in 0..<8 {
        let allocationSize = dictionaryCount[laneID]
        temporaryAllocationSize &+= allocationSize
      }
      
      withUnsafeTemporaryAllocation(
        of: UInt32.self,
        capacity: temporaryAllocationSize
      ) { bufferPointer in
        var start = 0
        for laneID in 0..<8 {
          let allocationSize = dictionaryCount[laneID]
          guard allocationSize > 0 else {
            continue
          }
          
          let oldPointer = dictionary.advanced(by: laneID &* atoms.count)
          let newPointer = bufferPointer.baseAddress.unsafelyUnwrapped + start
          newPointer.initialize(from: oldPointer, count: allocationSize)
          start &+= allocationSize
        }
        
        start = 0
        for laneID in 0..<8 {
          let allocationSize = dictionaryCount[laneID]
          guard allocationSize > 0 else {
            continue
          }
          
          let newPointer = bufferPointer.baseAddress.unsafelyUnwrapped + start
          start &+= allocationSize
          if allocationSize == 1 {
            output.append(newPointer.pointee)
            continue
          }
          
          let newBufferPointer = UnsafeBufferPointer(
            start: newPointer, count: allocationSize)
          let key32 = UInt32(truncatingIfNeeded: laneID)
          let intOffset = (key32 &>> SIMD3(0, 1, 2)) & 1
          let floatOffset = SIMD3<Float>(intOffset) * 2 - 1
          let newOrigin = levelOrigin + floatOffset * levelSize / 2
          traverse(
            atomIDs: newBufferPointer,
            levelOrigin: newOrigin,
            levelSize: levelSize / 2)
        }
      }
    }
    
    let levelSize = self.dimensions.max() / 2
    let levelOrigin = SIMD3<Float>(repeating: levelSize)
    let initialArray = atoms.indices.map(UInt32.init(truncatingIfNeeded:))
    initialArray.withUnsafeBufferPointer { bufferPointer in
      traverse(
        atomIDs: bufferPointer,
        levelOrigin: levelOrigin,
        levelSize: levelSize)
    }
    precondition(output.count == atoms.count)
    
    var reordering = [UInt32](repeating: .max, count: atoms.count)
    for reorderedID in output.indices {
      let originalID32 = output[reorderedID]
      let originalID = Int(truncatingIfNeeded: originalID32)
      let reorderedID32 = UInt32(truncatingIfNeeded: reorderedID)
      reordering[originalID] = reorderedID32
    }
    precondition(!reordering.contains(.max))
    return reordering
  }
}
