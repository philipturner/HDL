//
//  OctreeSorter+MultiThreaded.swift
//  HDL
//
//  Created by Philip Turner on 7/14/25.
//

import Dispatch
import QuartzCore

extension OctreeSorter {
  // Multi-threaded algorithm.
  func mortonReordering(grid: Grid) -> [UInt32] {
//    let start = CACurrentMediaTime()
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
        // Return early.
        if levelSize == 1 / 32 {
          localOutput += atomIDs
          return
        } else if levelSize < 1 / 32 {
          fatalError("This should never happen.")
        }
        
        // Use the scratch pad.
        var childNodeCounts: SIMD8<Int> = .zero
        for atomID in atomIDs {
          func createAtomOffset() -> SIMD3<Float> {
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
        
        // Create the temporary allocation.
        withUnsafeTemporaryAllocation(
          of: UInt32.self,
          capacity: childNodeCounts.wrappedSum()
        ) { allocationBuffer in
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
            
            func createNewOrigin() -> SIMD3<Float> {
              let intOffset = (key &>> SIMD3(0, 1, 2)) & 1
              let floatOffset = SIMD3<Float>(intOffset) * 2 - 1
              return levelOrigin + floatOffset * levelSize / 4
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
          levelSize: 4)
      }
      guard localOutput.count == cell.range.count else {
        fatalError("This should never happen.")
      }
      
      for localID in localOutput.indices {
        let globalID = cell.range.startIndex + localID
        globalOutput[globalID] = localOutput[localID]
      }
    }
    
    let taskCount = grid.cells.count
    if taskCount == 0 {
      fatalError("This should never happen.")
    } else if taskCount == 1 {
      let cell = grid.cells[0]
      execute(cell: cell)
    } else {
      DispatchQueue.concurrentPerform(iterations: taskCount) { z in
        let cell = grid.cells[z]
        execute(cell: cell)
      }
    }
    
//    let end = CACurrentMediaTime()
//    debugProfile(start, end, "part 2")
    return globalOutput
  }
}
