//
//  OctreeSorter+SingleThreaded.swift
//  HDL
//
//  Created by Philip Turner on 7/14/25.
//

extension OctreeSorter {
  // Single threaded algorithm without a grid.
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
            output.append(newPointer[0])
            continue
          }
          
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
    
    return output
  }
}
