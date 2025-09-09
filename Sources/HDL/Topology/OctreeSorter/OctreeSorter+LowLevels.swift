//
//  OctreeSorter+LowLevels.swift
//  HDL
//
//  Created by Philip Turner on 9/8/25.
//

import Dispatch

extension OctreeSorter {
  func traverseLowLevels(state: TraversalState) -> [UInt32] {
    guard let inPlaceBuffer_bypass = state.inPlaceBuffer,
          let levelSize = state.levelSize,
          let scratchBuffer_bypass = state.scratchBuffer,
          let threads = state.threads else {
      fatalError("State was not fully specified.")
    }
    
    // Bypass the annoying Swift concurrency warning.
    nonisolated(unsafe)
    let inPlaceBuffer = inPlaceBuffer_bypass
    nonisolated(unsafe)
    let scratchBuffer = scratchBuffer_bypass
    
    // Iterate over the threads (via concurrent dispatch).
    // Iterate over the cells within the threads.
    if threads.count == 0 {
      fatalError("This should never happen.")
    } else if threads.count == 1 {
      let thread = threads[0]
      for cell in thread.cells {
        traverseLowLevel(
          inPlaceBuffer: inPlaceBuffer,
          scratchBuffer: scratchBuffer,
          cellRange: cell.range,
          cellOrigin: cell.origin,
          levelSize: levelSize)
      }
    } else {
      DispatchQueue.concurrentPerform(
        iterations: threads.count
      ) { threadID in
        let thread = threads[threadID]
        for cell in thread.cells {
          traverseLowLevel(
            inPlaceBuffer: inPlaceBuffer,
            scratchBuffer: scratchBuffer,
            cellRange: cell.range,
            cellOrigin: cell.origin,
            levelSize: levelSize)
        }
      }
    }
    
    // Migrate the pointer's contents to an array, then deallocate.
    let output = [UInt32](unsafeUninitializedCapacity: atoms.count) {
      let baseAddress = $0.baseAddress!
      baseAddress.initialize(
        from: inPlaceBuffer, count: atoms.count)
      $1 = atoms.count
    }
    
    withExtendedLifetime(state) { }
    return output
  }
  
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
