//
//  OctreeSorter+Grid.swift
//  HDL
//
//  Created by Philip Turner on 7/14/25.
//

extension OctreeSorter {
  struct Cell {
    var range: Range<Int>
    var origin: SIMD3<Float>
    var size: Float
  }
  
  struct Grid {
    var data: [UInt32] = []
    var cells: [Cell] = []
  }
  
  func createGrid() -> Grid {
    let levelSizes = LevelSizes(dimensions: dimensions)
    var grid = Grid()
    
    // Create the scratch pad.
    let scratchPad: UnsafeMutablePointer<UInt32> =
      .allocate(capacity: 8 * atoms.count)
    defer { scratchPad.deallocate() }
    
    func traverseGrid(
      atomIDs: UnsafeBufferPointer<UInt32>,
      levelOrigin: SIMD3<Float>,
      levelSize: Float
    ) {
      if levelSize < levelSizes.threshold {
        let rangeStart = grid.data.count
        grid.data += atomIDs
        let rangeEnd = grid.data.count
        
        let cell = Cell(
          range: rangeStart..<rangeEnd,
          origin: levelOrigin,
          size: levelSize)
        grid.cells.append(cell)
        return
      }
      
      // Use the scratch pad.
      var childNodeCounts: SIMD8<Int> = .zero
      for atomID in atomIDs {
        // Compiler may not inline a nested function, so use a do statement.
        var atomOffset: SIMD3<Float>
        do {
          // @_transparent attribute is ineffective.
          let atom = atoms[Int(atomID)]
          let position = unsafeBitCast(atom, to: SIMD3<Float>.self)
          atomOffset = position - self.origin
        }
        
        var index = SIMD3<UInt32>(repeating: 1)
        index.replace(
          with: SIMD3.zero,
          where: atomOffset .< levelOrigin)
        
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
        // Compiler may not inline a nested function, so use a variable.
        let allocationPointer = allocationBuffer.baseAddress.unsafelyUnwrapped
        
        // Transfer the scratch pad to the temporary allocation.
        do {
          var cursor = 0
          for key in 0..<UInt32(8) {
            let childNodeCount = childNodeCounts[Int(key)]
            guard childNodeCount > 0 else {
              continue
            }
            
            let newPointer = allocationPointer + cursor
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
          
          let newPointer = allocationPointer + cursor
          cursor += childNodeCount
          
          // Compiler may not inline a nested function, so use a do statement.
          var newOrigin: SIMD3<Float>
          do {
            let intOffset = (key &>> SIMD3(0, 1, 2)) & 1
            let floatOffset = SIMD3<Float>(intOffset) * 2 - 1
            newOrigin = floatOffset * levelSize / 2
          }
          let newBufferPointer = UnsafeBufferPointer(
            start: newPointer,
            count: childNodeCount)
          traverseGrid(
            atomIDs: newBufferPointer,
            levelOrigin: newOrigin,
            levelSize: levelSize / 2)
        }
      }
    }
    
    // Invoke the traversal function the first time.
    let levelOrigin = SIMD3<Float>(repeating: levelSizes.octreeStart)
    let initialArray = atoms.indices.map(UInt32.init)
    initialArray.withUnsafeBufferPointer { bufferPointer in
      traverseGrid(
        atomIDs: bufferPointer,
        levelOrigin: levelOrigin,
        levelSize: levelSizes.octreeStart)
    }
    guard grid.data.count == atoms.count else {
      fatalError("This should never happen.")
    }
    
    return grid
  }
}
