//
//  Sort.swift
//  HDL
//
//  Created by Philip Turner on 12/2/23.
//

import Dispatch

extension Topology {
  @discardableResult
  public mutating func sort() -> [UInt32] {
    let grid = GridSorter(atoms: atoms)
    let reordering = grid.mortonReordering()
    let previousAtoms = atoms
    
    for i in reordering.indices {
      let mortonIndex = reordering[i]
      atoms[i] = previousAtoms[Int(mortonIndex)]
    }
    
    let inverted = GridSorter.invertOrder(reordering)
    for i in bonds.indices {
      let bond = SIMD2<Int>(truncatingIfNeeded: bonds[i])
      var newBond: SIMD2<UInt32> = .zero
      newBond[0] = inverted[bond[0]]
      newBond[1] = inverted[bond[1]]
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
    return inverted
  }
}

struct GridSorter {
  var atoms: [Atom] = []
  
  var origin: SIMD3<Float>
  var dimensions: SIMD3<Float>
  
  init(atoms: [Atom]) {
    self.atoms = atoms
    
    if atoms.count == 0 {
      origin = .zero
      dimensions = .init(repeating: 0.5)
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
}

// TODO: Change the algorithm to make the level size a power of 2? And
// clip the lower corner to a good power of 2? The current algorithm
// provides too erratic of performance, and something similar to chaos or
// nondeterminism in the results. The root cause is that the origin can
// vary drastically, depending on where the user offsets the lattice. I
// designed an algorithm that's robust to this fact, and some of the tests
// now depend on its exact behavior. We must revise it during the current
// round of source-breaking API changes.
//
// Implement this fix once the entire 'Sort' file is refactored and easier
// to work with, without introducing new regressions.
//
// New ideas:
// - Defer this to until you have a functioning renderer, to debug glitches
//   where visualizing the atoms helps massively.
// - Reduce the task count when the number of cells is very large,
//   merging nearby cells in the list.
//   - This creates a more homogeneous atom count and should reduce the spike
//     in overhead.
//   - Before doing this, thoroughly de-mystify the latency as a function of
//     lattice size, and have specific benchmarks (in milliseconds) to refer
//     to.
//
// Choice: defer changes to the underlying algorithm until we have a working
// renderer on Windows. Instead, include an attempt to improve workload
// distribution in this PR.
private struct LevelSizes {
  var highest: Float
  var octreeStart: Float
  var threshold: Float
  
  init(dimensions: SIMD3<Float>) {
    // Make an initial guess of 67% for the top-level binary divider.
    let volume = dimensions.x * dimensions.y * dimensions.z
    let chunkVolume = volume / 27
    self.highest = 2 * pow(chunkVolume, 1.0 / 3)
    self.octreeStart = highest
    
    // If the grid has dimensions that vary wildly, 'highestLevelSize' does not
    // provide an accurate center. Increase it by powers of 2 until it at least
    // reaches 51% of each dimension.
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
    
    // Very important: deciding the granularity with which to parallelize the
    // grid. 0.5 looks way too small for sparser lattices like SiC and Si.
    self.threshold = min(2.0, max(0.5, highest * 0.51))
  }
}

private struct GridCell {
  var range: Range<Int>
  var origin: SIMD3<Float>
  var size: Float
}

private struct Grid {
  var data: [UInt32] = []
  var cells: [GridCell] = []
}

// Performance tests:
//
// --filter testFlatSheetScaling
//
//   - sort: 0.056 µs/atom
//   - sort: 0.047 µs/atom
//   - sort: 0.055 µs/atom
//
// --filter testSort
//
//   atoms: 16400
//   dataset    | octree |  grid
//   ---------- | ------ | ------
//   pre-sorted |    924 |    617
//   lattice    |    904 |    592
//   shuffled   |   1046 |    646
//   reversed   |    901 |    595
//
//   atoms: 54900
//   dataset    | octree |  grid
//   ---------- | ------ | ------
//   pre-sorted |   3112 |   1874
//   lattice    |   3022 |   1972
//   shuffled   |   3669 |   2162
//   reversed   |   3179 |   1925
//
//   atoms: 129600
//   dataset    | octree |  grid
//   ---------- | ------ | ------
//   pre-sorted |   7929 |   4573
//   lattice    |   8033 |   4626
//   shuffled   |   9862 |   5254
//   reversed   |   8147 |   4752
//
// The sort inside the library is not actually slower than OctreeSorter.
// Use latticeScale=20 as the go-to test for quickly checking for a regression.

// Why is one algorithm faster than the other? What's going on at the lowest
// level, and is there more room for improvement?
//
// Grid algorithm for reordering:
//
//               | Part 1    | Part 2, parallel | Part 2, serial
// ------------- | --------- | ---------------- | --------------
//   16400 atoms |  0.1 ms   |  0.3 ms          |  0.9 ms
//  129600 atoms |  2.2 ms   |  1.2 ms          |  7.9 ms
//  435600 atoms | 11.3 ms   |  3.2 ms          | 23.4 ms
// 1030400 atoms | 32.5 ms   |  9.2 ms          | 70.3 ms
//
// Octree algorithm for reordering:
//
//               | Combined Pass
// ------------- | -------------
//   16400 atoms |  0.9 ms
//  129600 atoms |  8.8 ms
//  435600 atoms | 30.2 ms
// 1030400 atoms | 82.7 ms

extension GridSorter {
  static func invertOrder(_ input: [UInt32]) -> [UInt32] {
    return [UInt32](unsafeUninitializedCapacity: input.count) {
      $1 = input.count
      let baseAddress = $0.baseAddress.unsafelyUnwrapped
      
      for reorderedID in input.indices {
        let originalID = Int(input[reorderedID])
        baseAddress[originalID] = UInt32(reorderedID)
      }
    }
  }
  
  private func createGrid() -> Grid {
    let levelSizes = LevelSizes(dimensions: dimensions)
    var grid = Grid()
    
    let dictionary: UnsafeMutablePointer<UInt32> =
      .allocate(capacity: 8 * atoms.count)
    defer { dictionary.deallocate() }
    
    // TODO: Refactor this to move it outside of the enclosing function,
    // isolating the mutable context it sees. Do all of this without causing
    // a performance regression.
    func traverseGrid(
      atomIDs: UnsafeBufferPointer<UInt32>,
      levelOrigin: SIMD3<Float>,
      levelSize: Float
    ) {
      if levelSize < levelSizes.threshold {
        let rangeStart = grid.data.count
        grid.data += atomIDs
        let rangeEnd = grid.data.count
        
        let cell = GridCell(
          range: rangeStart..<rangeEnd,
          origin: levelOrigin,
          size: levelSize)
        grid.cells.append(cell)
        return
      }
      
      // Write to the dictionary.
      var dictionaryCount: SIMD8<Int> = .zero
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
        let previousCount = dictionaryCount[Int(key)]
        dictionaryCount[Int(key)] += 1
        dictionary[Int(key) * atomIDs.count + previousCount] = atomID
      }
      
      withUnsafeTemporaryAllocation(
        of: UInt32.self,
        capacity: dictionaryCount.wrappedSum()
      ) { allocationBuffer in
        @inline(__always)
        func allocationPointer() -> UnsafeMutablePointer<UInt32> {
          allocationBuffer.baseAddress.unsafelyUnwrapped
        }
        
        var start = 0
        for key in 0..<UInt32(8) {
          let allocationSize = dictionaryCount[Int(key)]
          guard allocationSize > 0 else {
            continue
          }
          
          let newPointer = allocationPointer() + start
          start += allocationSize
          
          newPointer.initialize(
            from: dictionary + Int(key) * atomIDs.count,
            count: allocationSize)
        }
        
        // TODO: Use a different variable, instead of letting a reference
        // to a mutable state variable survive across 2 different contexts.
        start = 0
        for key in 0..<UInt32(8) {
          let allocationSize = dictionaryCount[Int(key)]
          guard allocationSize > 0 else {
            continue
          }
          
          let newPointer = allocationPointer() + start
          start += allocationSize
          
          @inline(__always)
          func createNewOrigin() -> SIMD3<Float> {
            let intOffset = (key &>> SIMD3(0, 1, 2)) & 1
            let floatOffset = SIMD3<Float>(intOffset) * 2 - 1
            return levelOrigin + floatOffset * levelSize / 2
          }
          let newBufferPointer = UnsafeBufferPointer(
            start: newPointer,
            count: allocationSize)
          traverseGrid(
            atomIDs: newBufferPointer,
            levelOrigin: createNewOrigin(),
            levelSize: levelSize / 2)
        }
      }
    }
    
    let levelOrigin = SIMD3<Float>(repeating: levelSizes.octreeStart)
    let initialArray = atoms.indices.map(UInt32.init)
    initialArray.withUnsafeBufferPointer { bufferPointer in
      traverseGrid(
        atomIDs: bufferPointer,
        levelOrigin: levelOrigin,
        levelSize: levelSizes.octreeStart)
    }
    
    return grid
  }
  
  func mortonReordering() -> [UInt32] {
    let grid = createGrid()
    
    var finalOutput = [UInt32](unsafeUninitializedCapacity: atoms.count) {
      $1 = atoms.count
    }
    
    func createMaxCellAtomCount() -> Int {
      var output = 0
      for cell in grid.cells {
        let atomCount = cell.range.count
        if atomCount > output {
          output = atomCount
        }
      }
      return output
    }
    let maxCellAtomCount = createMaxCellAtomCount()
    
    // TODO: Fix the multiple errors that spawn when marking this function
    // as @Sendable.
    func execute(taskID: Int) {
      let dictionary: UnsafeMutablePointer<UInt32> =
        .allocate(capacity: 8 * maxCellAtomCount)
      defer { dictionary.deallocate() }
      
      var output: [UInt32] = []
      
      let cell = grid.cells[taskID]
      let initialArray = grid.data[cell.range]
      initialArray.withUnsafeBufferPointer { bufferPointer in
        traverseTree(
          atomIDs: bufferPointer,
          levelOrigin: cell.origin,
          levelSize: cell.size)
      }
      
      let finalStartI = cell.range.startIndex
      for outputI in output.indices {
        let finalI = finalStartI + outputI
        finalOutput[finalI] = output[outputI]
      }
      
      func traverseTree(
        atomIDs: UnsafeBufferPointer<UInt32>,
        levelOrigin: SIMD3<Float>,
        levelSize: Float
      ) {
        if levelSize < 1 / 32 {
          output += atomIDs
          return
        }
        
        // Write to the dictionary.
        var dictionaryCount: SIMD8<Int> = .zero
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
          let previousCount = dictionaryCount[Int(key)]
          dictionaryCount[Int(key)] += 1
          dictionary[Int(key) * atomIDs.count + previousCount] = atomID
        }
        
        withUnsafeTemporaryAllocation(
          of: UInt32.self,
          capacity: dictionaryCount.wrappedSum()
        ) { allocationBuffer in
          @inline(__always)
          func allocationPointer() -> UnsafeMutablePointer<UInt32> {
            allocationBuffer.baseAddress.unsafelyUnwrapped
          }
          
          var start = 0
          for key in 0..<UInt32(8) {
            let allocationSize = dictionaryCount[Int(key)]
            guard allocationSize > 0 else {
              continue
            }
            
            let newPointer = allocationPointer() + start
            start += allocationSize
            
            newPointer.initialize(
              from: dictionary + Int(key) * atomIDs.count,
              count: allocationSize)
          }
          
          // TODO: Use a different variable, instead of letting a reference
          // to a mutable state variable survive across 2 different contexts.
          start = 0
          for key in 0..<UInt32(8) {
            let allocationSize = dictionaryCount[Int(key)]
            guard allocationSize > 0 else {
              continue
            }
            
            let newPointer = allocationPointer() + start
            start += allocationSize
            
            if allocationSize == 1 {
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
              count: allocationSize)
            traverseTree(
              atomIDs: newBufferPointer,
              levelOrigin: createNewOrigin(),
              levelSize: levelSize / 2)
          }
        }
      }
    }
    
    do {
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
      let largeGridCellCount = createLargeCellCount()
      
      if largeGridCellCount >= 3 {
        let taskCount = grid.cells.count
//        DispatchQueue.concurrentPerform(
//          iterations: taskCount,
//          execute: execute(taskID:))
        
        // TODO: Check for a negative impact from extra function calls, when
        // latticeScale=5
        DispatchQueue.concurrentPerform(iterations: taskCount) { z in
          execute(taskID: z)
        }
      } else {
        for z in grid.cells.indices {
          execute(taskID: z)
        }
      }
    }
    
    return finalOutput
  }
}
