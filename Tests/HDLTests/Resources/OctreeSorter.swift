//
//  OctreeSorter.swift
//  HDLTests
//
//  Created by Philip Turner on 12/25/23.
//

import Foundation
import HDL

private struct LevelSizes {
  var highest: Float
  var octreeStart: Float
  
  init(dimensions: SIMD3<Float>) {
    let volume = dimensions.x * dimensions.y * dimensions.z
    let chunkVolume = volume / 27
    self.highest = 2 * pow(chunkVolume, 1.0 / 3)
    self.octreeStart = highest
    
    var iterationCount = 0
    while true {
      iterationCount += 1
      if iterationCount > 100 {
        fatalError("Too many iterations.")
      }
      
      let targetDimension = dimensions.max() * 0.51
      if octreeStart < targetDimension {
        octreeStart *= 2
      } else {
        break
      }
    }
  }
}

struct OctreeSorter {
  var atoms: [Atom] = []
  
  var origin: SIMD3<Float>
  var dimensions: SIMD3<Float>
  
  init(atoms: [Atom]) {
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
      
      origin = minimum
      dimensions = maximum - minimum
      dimensions.replace(
        with: SIMD3(repeating: 0.5),
        where: dimensions .< 0.5)
    }
  }
  
  func mortonReordering() -> [UInt32] {
    var output: [UInt32] = []
    let dictionary: UnsafeMutablePointer<UInt32> =
      .allocate(capacity: 8 * atoms.count)
    defer { dictionary.deallocate() }
    
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
        @inline(__always)
        func createAtomPosition() -> SIMD3<Float> {
          let atomID = Int(atomID32)
          return atoms[atomID].position - self.origin
        }
        
        var index = SIMD3<UInt32>(repeating: 1)
        index.replace(
          with: .init(repeating: 0),
          where: createAtomPosition() .< levelOrigin)
        
        let key = Int((index &<< SIMD3(0, 1, 2)).wrappedSum())
        let previousCount = dictionaryCount[key]
        let pointer = dictionary.advanced(
          by: key * atoms.count + previousCount)
        
        dictionaryCount[key] += 1
        pointer.pointee = atomID32
      }
      
      var temporaryAllocationSize = 0
      for laneID in 0..<8 {
        let allocationSize = dictionaryCount[laneID]
        temporaryAllocationSize += allocationSize
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
          
          let key32 = UInt32(truncatingIfNeeded: laneID)
          let intOffset = (key32 &>> SIMD3(0, 1, 2)) & 1
          let floatOffset = SIMD3<Float>(intOffset) * 2 - 1
          let newOrigin = levelOrigin + floatOffset * levelSize / 2
          
          let newBufferPointer = UnsafeBufferPointer(
            start: newPointer, count: allocationSize)
          traverse(
            atomIDs: newBufferPointer,
            levelOrigin: newOrigin,
            levelSize: levelSize / 2)
        }
      }
    }
    
    let levelSizes = LevelSizes(dimensions: dimensions)
    let levelOrigin = SIMD3<Float>(repeating: levelSizes.octreeStart)
    let initialArray = atoms.indices.map(UInt32.init(truncatingIfNeeded:))
    initialArray.withUnsafeBufferPointer { bufferPointer in
      traverse(
        atomIDs: bufferPointer,
        levelOrigin: levelOrigin,
        levelSize: levelSizes.octreeStart)
    }
    guard output.count == atoms.count else {
      fatalError("This should never happen.")
    }
    
    var reordering = [UInt32](repeating: .max, count: atoms.count)
    for reorderedID in output.indices {
      let originalID32 = output[reorderedID]
      let originalID = Int(truncatingIfNeeded: originalID32)
      let reorderedID32 = UInt32(truncatingIfNeeded: reorderedID)
      reordering[originalID] = reorderedID32
    }
    return reordering
  }
}
