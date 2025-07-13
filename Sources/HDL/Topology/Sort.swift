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
    
    let inverted = grid.invertOrder(reordering)
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

// MARK: - Morton Reordering

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
  func invertOrder(_ input: [UInt32]) -> [UInt32] {
    return [UInt32](unsafeUninitializedCapacity: input.count) {
      $1 = input.count
      let baseAddress = $0.baseAddress.unsafelyUnwrapped
      
      for reorderedID in input.indices {
        let originalID = Int(input[reorderedID])
        baseAddress[originalID] = UInt32(reorderedID)
      }
    }
  }
  
  func mortonReordering() -> [UInt32] {
    var gridData: [UInt32] = []
    var gridCells: [(Range<Int>, SIMD3<Float>, Float)] = []
    let levelSizes = LevelSizes(dimensions: dimensions)
    
    do {
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
          let rangeStart = gridData.count
          gridData += atomIDs
          let rangeEnd = gridData.count
          gridCells.append((rangeStart..<rangeEnd, levelOrigin, levelSize))
          return
        }
        
        // Write to the dictionary.
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
          
          let key = (index &<< SIMD3(0, 1, 2)).wrappedSum()
          let previousCount = dictionaryCount[Int(key)]
          let pointer = dictionary.advanced(
            by: Int(key) * atoms.count + previousCount)
          
          dictionaryCount[Int(key)] += 1
          pointer.pointee = atomID32
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
            
            let oldPointer = dictionary + Int(key) * atoms.count
            let newPointer = allocationPointer() + start
            newPointer.initialize(from: oldPointer, count: allocationSize)
            start += allocationSize
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
            
            let newBufferPointer = UnsafeBufferPointer(
              start: newPointer, count: allocationSize)
            let intOffset = (key &>> SIMD3(0, 1, 2)) & 1
            let floatOffset = SIMD3<Float>(intOffset) * 2 - 1
            let newOrigin = levelOrigin + floatOffset * levelSize / 2
            traverseGrid(
              atomIDs: newBufferPointer,
              levelOrigin: newOrigin,
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
    }
    
    var finalOutput = [UInt32](unsafeUninitializedCapacity: atoms.count) {
      $1 = atoms.count
    }
    
    var maxCellSize = 0
    var largeGridCellCount = 0
    for gridCell in gridCells {
      maxCellSize = max(maxCellSize, gridCell.0.count)
      if gridCell.0.count > 64 {
        largeGridCellCount += 1
      }
    }
    
    if largeGridCellCount >= 3 {
      DispatchQueue.concurrentPerform(
        iterations: gridCells.count, execute: executeIteration(_:))
    } else {
      for z in gridCells.indices {
        executeIteration(z)
      }
    }
    
    // TODO: Fix the multiple errors that spawn when marking this function
    // as @Sendable.
    func executeIteration(_ z: Int) {
      let gridCell = gridCells[z]
      
      let dictionary: UnsafeMutablePointer<UInt32> =
        .allocate(capacity: 8 * maxCellSize)
      defer { dictionary.deallocate() }
      
      var output: [UInt32] = []
      
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
          
          let key = (index &<< SIMD3(0, 1, 2)).wrappedSum()
          let previousCount = dictionaryCount[Int(key)]
          let pointer = dictionary.advanced(
            by: Int(key) * maxCellSize + previousCount)
          
          dictionaryCount[Int(key)] += 1
          pointer.pointee = atomID32
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
            
            let oldPointer = dictionary + Int(key) * maxCellSize
            let newPointer = allocationPointer() + start
            newPointer.initialize(from: oldPointer, count: allocationSize)
            start += allocationSize
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
              output.append(newPointer.pointee)
              continue
            }
            
            let newBufferPointer = UnsafeBufferPointer(
              start: newPointer, count: allocationSize)
            let intOffset = (key &>> SIMD3(0, 1, 2)) & 1
            let floatOffset = SIMD3<Float>(intOffset) * 2 - 1
            let newOrigin = levelOrigin + floatOffset * levelSize / 2
            traverseTree(
              atomIDs: newBufferPointer,
              levelOrigin: newOrigin,
              levelSize: levelSize / 2)
          }
        }
      }
      
      let levelSize = gridCell.2
      let levelOrigin = gridCell.1
      let initialArray = gridData[gridCell.0]
      initialArray.withUnsafeBufferPointer { bufferPointer in
        traverseTree(
          atomIDs: bufferPointer,
          levelOrigin: levelOrigin,
          levelSize: levelSize)
      }
      
      let finalStartI = gridCell.0.startIndex
      for outputI in output.indices {
        let finalI = finalStartI + outputI
        finalOutput[finalI] = output[outputI]
      }
    }
    
    return finalOutput
  }
}
