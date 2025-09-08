//
//  OctreeSorter+LowLevels.swift
//  HDL
//
//  Created by Philip Turner on 9/8/25.
//

// Reminder: If a cell has just 1 atom at any point, leave the function
// immediately. It should not have too much overhead to defer this choice
// to 0.5 nm iterations, if a 1 nm cell has just 1 atom.

extension OctreeSorter {
  @Sendable
  private func traverseLowLevel(
    inPlaceBuffer: UnsafeMutablePointer<UInt32>,
    scratchBuffer: UnsafeMutablePointer<UInt32>,
    cellRange: Range<Int>,
    cellOrigin: SIMD3<Float>,
    levelSize: Float
  ) {
    // Use the scratch buffer.
    func createChildNodeSizes() -> SIMD8<UInt32> {
      var childNodeSizes: SIMD8<UInt32> = .zero
      let scratchStart = UInt32(cellRange.startIndex * 8)
      let scratchStride = UInt32(cellRange.count)
      
      for inPlaceOffset in cellRange {
        let atomID = inPlaceBuffer[inPlaceOffset]
        func createAtomOffset() -> SIMD3<Float> {
          let atom = atoms[Int(atomID)]
          let position = unsafeBitCast(atom, to: SIMD3<Float>.self)
          return position - self.origin
        }
        
        var index = SIMD3<UInt32>(repeating: 1)
        index.replace(
          with: SIMD3.zero,
          where: createAtomOffset() .< cellOrigin)
        
        let childNodeID = (index &<< SIMD3(0, 1, 2)).wrappedSum()
        let previousSize = childNodeSizes[Int(childNodeID)]
        childNodeSizes[Int(childNodeID)] = previousSize + 1
        
        var scratchOffset = scratchStart + childNodeID * scratchStride
        scratchOffset += previousSize
        scratchBuffer[Int(scratchOffset)] = atomID
      }
      return childNodeSizes
    }
    let childNodeSizes = createChildNodeSizes()
    
    // Transfer the scratch buffer to the input buffer.
    func createChildNodeOffsets() -> SIMD8<UInt32> {
      var childNodeOffsets: SIMD8<UInt32> = .zero
      let scratchStart = UInt32(cellRange.startIndex * 8)
      let scratchStride = UInt32(cellRange.count)
      
      var inPlaceOffset = cellRange.startIndex
      for childNodeID in 0..<UInt32(8) {
        let childNodeSize = childNodeSizes[Int(childNodeID)]
        guard childNodeSize > 0 else {
          continue
        }
        
        let scratchOffset = scratchStart + childNodeID * scratchStride
        (inPlaceBuffer + inPlaceOffset).initialize(
          from: scratchBuffer + Int(scratchOffset),
          count: Int(childNodeSize))
        
        childNodeOffsets[Int(childNodeID)] = UInt32(inPlaceOffset)
        inPlaceOffset += Int(childNodeSize)
      }
      return childNodeOffsets
    }
    let childNodeOffsets = createChildNodeOffsets()
    if (levelSize / 2) <= Float(1.0 / 32) {
      return
    }
    
    for childNodeID in 0..<UInt32(8) {
      let childNodeSize = Int(childNodeSizes[Int(childNodeID)])
      let inPlaceOffset = Int(childNodeOffsets[Int(childNodeID)])
      if childNodeSize <= 1 {
        continue
      }
      
      func createNewOrigin() -> SIMD3<Float> {
        let intOffset = (UInt8(childNodeID) &>> SIMD3(0, 1, 2)) & 1
        let floatOffset = SIMD3<Float>(intOffset) * 2 - 1
        return cellOrigin + floatOffset * levelSize / 4
      }
      
      traverseLowLevel(
        inPlaceBuffer: inPlaceBuffer,
        scratchBuffer: scratchBuffer,
        cellRange: inPlaceOffset..<(inPlaceOffset + childNodeSize),
        cellOrigin: createNewOrigin(),
        levelSize: levelSize / 2)
    }
  }
}
