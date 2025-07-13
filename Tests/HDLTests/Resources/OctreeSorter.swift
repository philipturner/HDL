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
    
    // Create the scratch pad.
    let scratchPad: UnsafeMutablePointer<UInt32> =
      .allocate(capacity: 8 * atoms.count)
    defer { scratchPad.deallocate() }
    
    func traverse(
      atomIDs: UnsafeBufferPointer<UInt32>,
      levelOrigin: SIMD3<Float>,
      levelSize: Float
    ) {
      if levelSize < 1 / 32 {
        output += atomIDs
        return
      }
      
      // Use the scratch pad.
      var childNodeCounts: SIMD8<Int> = .zero
      for atomID in atomIDs {
        @inline(__always)
        func createAtomOffset() -> SIMD3<Float> {
          let atom = atoms[Int(atomID)]
          return atom.position - self.origin
        }
        var index = SIMD3<UInt32>(repeating: 1)
        index.replace(
          with: SIMD3.zero,
          where: createAtomOffset() .< levelOrigin)
        
        let key = (index &<< SIMD3(0, 1, 2)).wrappedSum()
        let previousCount = childNodeCounts[Int(key)]
        childNodeCounts[Int(key)] += 1
        scratchPad[Int(key) * atomIDs.count + previousCount] = atomID
      }
      
      // Create the temporary allocation.
      withUnsafeTemporaryAllocation(
        of: UInt32.self,
        capacity: childNodeCounts.wrappedSum()
      ) { allocationBuffer in
        @inline(__always)
        func allocationPointer() -> UnsafeMutablePointer<UInt32> {
          allocationBuffer.baseAddress.unsafelyUnwrapped
        }
        
        // Transfer the scratch pad to the temporary allocation.
        do {
          var cursor = 0
          for key in 0..<UInt32(8) {
            let childNodeCount = childNodeCounts[Int(key)]
            guard childNodeCount > 0 else {
              continue
            }
            
            let newPointer = allocationPointer() + cursor
            cursor += childNodeCount
            
            newPointer.initialize(
              from: scratchPad + Int(key) * atomIDs.count,
              count: childNodeCount)
          }
        }
        
        // Invoke the traversal function recursively.
        var cursor = 0
        for key in 0..<UInt32(8) {
          let childNodeCount = childNodeCounts[Int(key)]
          guard childNodeCount > 0 else {
            continue
          }
          
          let newPointer = allocationPointer() + cursor
          cursor += childNodeCount
          
          if childNodeCount == 1 {
            output.append(newPointer[0])
            continue
          }
          
          @inline(__always)
          func createNewOrigin() -> SIMD3<Float> {
            let intOffset = (key &>> SIMD3(0, 1, 2)) & 1
            let floatOffset = SIMD3<Float>(intOffset) * 2 - 1
            return levelOrigin + floatOffset * levelSize / 2
          }
          let newBufferPointer = UnsafeBufferPointer(
            start: newPointer,
            count: childNodeCount)
          traverse(
            atomIDs: newBufferPointer,
            levelOrigin: createNewOrigin(),
            levelSize: levelSize / 2)
        }
      }
    }
    
    // Invoke the traversal function the first time.
    let levelSizes = LevelSizes(dimensions: dimensions)
    let levelOrigin = SIMD3<Float>(repeating: levelSizes.octreeStart)
    let initialArray = atoms.indices.map(UInt32.init)
    initialArray.withUnsafeBufferPointer { bufferPointer in
      traverse(
        atomIDs: bufferPointer,
        levelOrigin: levelOrigin,
        levelSize: levelSizes.octreeStart)
    }
    guard output.count == atoms.count else {
      fatalError("This should never happen.")
    }
    
    // Note: Do not include this part when benchmarking performance of the
    // octree sorter. It unfairly biases performance in favor of the grid
    // sorter.
    var reordering = [UInt32](repeating: .max, count: atoms.count)
    for reorderedID in output.indices {
      let originalID32 = output[reorderedID]
      let originalID = Int(originalID32)
      let reorderedID32 = UInt32(reorderedID)
      reordering[originalID] = reorderedID32
    }
    return reordering
  }
}
