//
//  OctreeSorter+HighLevels.swift
//  HDL
//
//  Created by Philip Turner on 9/7/25.
//

import Dispatch

extension OctreeSorter {
  struct Cell {
    var range: Range<Int> = 0..<0
    var origin: SIMD3<Float> = .zero
  }
  
  struct Thread {
    var cells: [Cell]
  }
  
  struct TraversalState {
    var atomIDs: [UInt32] = []
    
    // Allows the lowest level size to be changed at will, without needing
    // to modify the function for low levels.
    var levelSize: Float = .zero
    
    var threads: [Thread] = []
  }
  
  func mortonReorderingHighLevels() -> TraversalState {
    nonisolated(unsafe)
    let inPlaceBuffer: UnsafeMutablePointer<UInt32> =
      .allocate(capacity: atoms.count)
    nonisolated(unsafe)
    let scratchBuffer: UnsafeMutablePointer<UInt32> =
      .allocate(capacity: 8 * atoms.count)
    defer { scratchBuffer.deallocate() }
    
    // Initialize the list of atom IDs.
    for i in 0..<atoms.count {
      inPlaceBuffer[i] = UInt32(i)
    }
    
    nonisolated(unsafe)
    var threads: [Thread] = [createFirstThread()]
    nonisolated(unsafe)
    var levelSize = highestLevelSize
    while levelSize > 1 {
      // Perform a prefix sum, to allocate memory for outputs.
      let threadCellOffsets = Self.cellOffsets(threads: threads)
      let inputCellCount = threads.map(\.cells.count).reduce(0, +)
      
      // Thread-safe buffers for storing results.
      nonisolated(unsafe)
      var results = LevelResults(
        threadCount: threads.count,
        inputCellCount: inputCellCount)
      
      DispatchQueue.concurrentPerform(
        iterations: threads.count
      ) { threadID in
        let thread = threads[threadID]
        let inputPrefixSum = Int(threadCellOffsets[threadID])
        
        var parentCells: [Cell] = []
        var children: [Thread] = []
        for cell in thread.cells {
          let output = traverseHighLevel(
            inPlaceBuffer: inPlaceBuffer,
            scratchBuffer: scratchBuffer,
            cell: cell,
            levelSize: levelSize)
          guard output.count > 0 else {
            fatalError("Unexpected output cell count.")
          }
          
          if output.count == 1 {
            parentCells += output[0].cells
          } else {
            children += output
          }
        }
        
        // Write the parent cells into the results.
        results.outputCellsPerParent[threadID] = UInt32(parentCells.count)
        for cellID in parentCells.indices {
          let cell = parentCells[cellID]
          let cellOffset = inputPrefixSum * 8 + cellID
          results.outputParentCells[cellOffset] = cell
        }
        
        // Write the children into the results.
        results.childrenPerParent[threadID] = UInt32(children.count)
        var childPrefixSum: Int = .zero
        for childID in children.indices {
          let child = children[childID]
          let childOffset = inputPrefixSum * 8 + childID
          let cellCount = child.cells.count
          results.outputCellsPerChild[childOffset] = UInt32(cellCount)
          
          for cellID in child.cells.indices {
            let cell = child.cells[cellID]
            let cellOffset = inputPrefixSum * 8 + childPrefixSum + cellID
            results.outputChildCells[cellOffset] = cell
          }
          childPrefixSum += cellCount
        }
      }
      
      // Reconstruct 'parentCells' and 'children', scan-compact the list.
      threads = results.nextLevel(
        threadCellOffsets: threadCellOffsets)
      levelSize /= 2
    }
    
    // Migrate the pointer's contents to an array, then deallocate.
    let atomIDs = [UInt32](unsafeUninitializedCapacity: atoms.count) {
      let baseAddress = $0.baseAddress!
      baseAddress.initialize(
        from: inPlaceBuffer, count: atoms.count)
      $1 = atoms.count
    }
    inPlaceBuffer.deallocate()
    
    var output = TraversalState()
    output.atomIDs = atomIDs
    output.levelSize = 1
    output.threads = threads
    return output
  }
}

extension OctreeSorter {
  private func createFirstThread() -> Thread {
    var cell = Cell()
    cell.range = atoms.indices
    cell.origin = SIMD3<Float>(
      repeating: highestLevelSize / 2)
    
    let thread = Thread(cells: [cell])
    return thread
  }
  
  private static func cellOffsets(threads: [Thread]) -> [UInt32] {
    var output: [UInt32] = []
    var counter: UInt32 = .zero
    for thread in threads {
      output.append(counter)
      counter += UInt32(thread.cells.count)
    }
    return output
  }
  
  private struct LevelResults {
    var outputCellsPerParent: [UInt32]
    var outputParentCells: [Cell]
    
    var childrenPerParent: [UInt32]
    var outputCellsPerChild: [UInt32]
    var outputChildCells: [Cell]
    
    init(
      threadCount: Int,
      inputCellCount: Int
    ) {
      outputCellsPerParent = [UInt32](
        repeating: 0, count: threadCount)
      outputParentCells = [Cell](
        repeating: Cell(), count: 8 * inputCellCount)
      
      childrenPerParent = [UInt32](
        repeating: 0, count: threadCount)
      outputCellsPerChild = [UInt32](
        repeating: 0, count: 8 * inputCellCount)
      outputChildCells = [Cell](
        repeating: Cell(), count: 8 * inputCellCount)
    }
    
    func nextLevel(threadCellOffsets: [UInt32]) -> [Thread] {
      var output: [Thread] = []
      for threadID in outputCellsPerParent.indices {
        let inputPrefixSum = Int(threadCellOffsets[threadID])
        
        // Reconstruct 'parentCells'.
        let parentCellCount = Int(outputCellsPerParent[threadID])
        var parentCells: [Cell] = []
        for cellID in 0..<parentCellCount {
          let cellOffset = inputPrefixSum * 8 + cellID
          let cell = outputParentCells[cellOffset]
          parentCells.append(cell)
        }
        
        // Reconstruct 'children'.
        let childCount = Int(childrenPerParent[threadID])
        var children: [Thread] = []
        var childPrefixSum: Int = .zero
        for childID in 0..<childCount {
          let childOffset = inputPrefixSum * 8 + childID
          let cellCount = Int(outputCellsPerChild[childOffset])
          
          var childCells: [Cell] = []
          for cellID in 0..<cellCount {
            let cellOffset = inputPrefixSum * 8 + childPrefixSum + cellID
            let cell = outputChildCells[cellOffset]
            childCells.append(cell)
          }
          childPrefixSum += cellCount
          
          let child = Thread(cells: childCells)
          children.append(child)
        }
        
        // Compact the vacant list entries.
        if parentCells.count > 0 {
          let parent = Thread(cells: parentCells)
          output.append(parent)
        }
        output += children
      }
      return output
    }
  }
  
  @Sendable
  private func traverseHighLevel(
    inPlaceBuffer: UnsafeMutablePointer<UInt32>,
    scratchBuffer: UnsafeMutablePointer<UInt32>,
    cell: Cell,
    levelSize: Float
  ) -> [Thread] {
    // Use the scratch buffer.
    func createChildNodeSizes() -> SIMD8<UInt32> {
      var childNodeSizes: SIMD8<UInt32> = .zero
      let scratchStart = UInt32(cell.range.startIndex * 8)
      let scratchStride = UInt32(cell.range.count)
      
      for inPlaceOffset in cell.range {
        let atomID = inPlaceBuffer[inPlaceOffset]
        func createAtomOffset() -> SIMD3<Float> {
          let atom = atoms[Int(atomID)]
          let position = unsafeBitCast(atom, to: SIMD3<Float>.self)
          return position - self.origin
        }
        
        var index = SIMD3<UInt32>(repeating: 1)
        index.replace(
          with: SIMD3.zero,
          where: createAtomOffset() .< cell.origin)
        
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
      let scratchStart = UInt32(cell.range.startIndex * 8)
      let scratchStride = UInt32(cell.range.count)
      
      var prefixSum: Int = .zero
      for childNodeID in 0..<UInt32(8) {
        let childNodeSize = childNodeSizes[Int(childNodeID)]
        guard childNodeSize > 0 else {
          continue
        }
        
        let scratchOffset = scratchStart + childNodeID * scratchStride
        (inPlaceBuffer + prefixSum).initialize(
          from: scratchBuffer + Int(scratchOffset),
          count: Int(childNodeSize))
        
        childNodeOffsets[Int(childNodeID)] = UInt32(prefixSum)
        prefixSum += Int(childNodeSize)
      }
      return childNodeOffsets
    }
    
    func createRemainingLevels() -> Int {
      let exponentFor4 = Int(Float(4).exponentBitPattern)
      let exponentForLevel = Int(levelSize.exponentBitPattern)
      return (exponentForLevel - exponentFor4) + 7
    }
    func createChildLatencies() -> SIMD8<Float> {
      // 5.0 ns/atom/level
      var output = SIMD8<Float>(childNodeSizes)
      output *= Float(5.0e-9)
      output *= Float(createRemainingLevels())
      return output
    }
    let childLatencies = createChildLatencies()
    let workSplitting = WorkSplitting(childLatencies: childLatencies)
    
    
    
    fatalError("Not implemented.")
  }
}
