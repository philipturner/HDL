//
//  OctreeSorter+MultiThreaded.swift
//  HDL
//
//  Created by Philip Turner on 7/14/25.
//

import Dispatch

extension OctreeSorter {
  // Multi-threaded algorithm.
  func mortonReordering(grid: Grid) -> [UInt32] {
    nonisolated(unsafe)
    var globalOutput = [UInt32](unsafeUninitializedCapacity: atoms.count) {
      $1 = atoms.count
    }
    
    @Sendable
    func execute(cell: Cell) {
      var localOutput: [UInt32] = []
      
      // Create the scratch pad.
      let scratchPad: UnsafeMutablePointer<UInt32> =
        .allocate(capacity: 8 * cell.range.count)
      defer { scratchPad.deallocate() }
      
      func traverseTree(
        atomIDs: UnsafeBufferPointer<UInt32>,
        levelOrigin: SIMD3<Float>,
        levelSize: Float
      ) {
        if levelSize < 1 / 32 {
          localOutput += atomIDs
          return
        }
        
        // Use the scratch pad.
        var childNodeCounts: SIMD8<Int> = .zero
        for atomID in atomIDs {
          @inline(__always)
          func createAtomOffset() -> SIMD3<Float> {
            // @_transparent attribute is ineffective.
            let atom = atoms[Int(atomID)]
            let position = unsafeBitCast(atom, to: SIMD3<Float>.self)
            return position - self.origin
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
//        let childNodeCounts = createChildNodeCounts(
//          atomIDs: atomIDs,
//          levelOrigin: levelOrigin,
//          scratchPad: scratchPad)
        
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
              localOutput.append(newPointer[0])
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
            traverseTree(
              atomIDs: newBufferPointer,
              levelOrigin: createNewOrigin(),
              levelSize: levelSize / 2)
          }
        }
      }
      
      // Invoke the traversal function the first time.
      let initialArray = grid.data[cell.range]
      initialArray.withUnsafeBufferPointer { bufferPointer in
        traverseTree(
          atomIDs: bufferPointer,
          levelOrigin: cell.origin,
          levelSize: cell.size)
      }
      guard localOutput.count == cell.range.count else {
        fatalError("This should never happen.")
      }
      
      for localID in localOutput.indices {
        let globalID = cell.range.startIndex + localID
        globalOutput[globalID] = localOutput[localID]
      }
    }
    
    func createLargeCellCount() -> Int {
      var output = 0
      for cell in grid.cells {
        let atomCount = cell.range.count
        if atomCount > 64 {
          output += 1
        }
      }
      return output
    }
    
    let taskCount = grid.cells.count
    if createLargeCellCount() < 3 {
      for z in 0..<taskCount {
        let cell = grid.cells[z]
        execute(cell: cell)
      }
    } else {
      DispatchQueue.concurrentPerform(iterations: taskCount) { z in
        let cell = grid.cells[z]
        execute(cell: cell)
      }
    }
    
    return globalOutput
  }
}
