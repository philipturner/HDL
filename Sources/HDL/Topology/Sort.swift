//
//  Sort.swift
//  HDL
//
//  Created by Philip Turner on 12/2/23.
//

import Dispatch

// TODO: Remove this temporary import
import QuartzCore

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

struct GridSorter {
  var atoms: [Atom] = []
  
  var origin: SIMD3<Float>
  var dimensions: SIMD3<Float>
  
  init(atoms: [Atom]) {
    self.atoms = atoms
    
    // TODO: Revisit the bin sizes, measure whether they're optimal for
    // sparser crystals like silicon carbide and silicon. Is there a net
    // slowdown for these crystals?
    //
    // Is a grid cell inherently holding 22 or 176 atoms for diamond? Profile
    // this before continuing the code cleanup.
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
//
// Why is one algorithm faster than the other? What's going on at the lowest
// level, and is there more room for improvement?

extension GridSorter {
  func invertOrder(_ input: [UInt32]) -> [UInt32] {
    return [UInt32](unsafeUninitializedCapacity: input.count) {
      $1 = input.count
      let baseAddress = $0.baseAddress.unsafelyUnwrapped
      
      for reorderedID in input.indices {
        let originalID = Int(truncatingIfNeeded: input[reorderedID])
        baseAddress[originalID] = UInt32(truncatingIfNeeded: reorderedID)
      }
    }
  }
  
  func mortonReordering() -> [UInt32] {
    var gridData: [UInt32] = []
    var gridCells: [(Range<Int>, SIMD3<Float>, Float)] = []
    
    // Make an initial guess of 67% for the top-level binary divider.
    let volume = dimensions.x * dimensions.y * dimensions.z
    let chunkVolume = volume / 27
    let highestLevelSize = 2 * pow(chunkVolume, 1.0 / 3)
    var octreeStartSize = highestLevelSize
    
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
      if octreeStartSize < targetDimension {
        octreeStartSize *= 2
      } else {
        break
      }
    }
    
    do {
      let dictionary: UnsafeMutablePointer<UInt32> =
        .allocate(capacity: 8 * atoms.count)
      defer { dictionary.deallocate() }
      
      let threshold = min(2.0, max(0.5, highestLevelSize * 0.51))
      
      // TODO: Refactor this to move it outside of the enclosing function,
      // isolating the mutable context it sees. Do all of this without causing
      // a performance regression.
      //
      // TODO: Observe the similarities between 'traverseGrid' and
      // 'traverseTree'.
      func traverseGrid(
        atomIDs: UnsafeBufferPointer<UInt32>,
        levelOrigin: SIMD3<Float>,
        levelSize: Float
      ) {
        if levelSize < threshold {
          let rangeStart = gridData.count
          gridData += atomIDs
          let rangeEnd = gridData.count
          gridCells.append((rangeStart..<rangeEnd, levelOrigin, levelSize))
          return
        }
        
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
          
          let key = Int((index &<< SIMD3(0, 1, 2)).wrappedSum())
          let previousCount = dictionaryCount[key]
          let pointer = dictionary.advanced(
            by: key * atoms.count + previousCount)
          
          dictionaryCount[key] += 1
          pointer.pointee = atomID32
        }
        
        var temporaryAllocationSize = 0
        for laneID in 0..<8 {
          let allocationSize = dictionaryCount[laneID]
          temporaryAllocationSize += allocationSize
        }
        
        withUnsafeTemporaryAllocation(
          of: UInt32.self,
          capacity: temporaryAllocationSize
        ) { bufferPointer in
          var start = 0
          for laneID in 0..<8 {
            let allocationSize = dictionaryCount[laneID]
            guard allocationSize > 0 else {
              continue
            }
            
            let oldPointer = dictionary.advanced(by: laneID &* atoms.count)
            let newPointer = bufferPointer.baseAddress.unsafelyUnwrapped + start
            newPointer.initialize(from: oldPointer, count: allocationSize)
            start &+= allocationSize
          }
          
          // TODO: Use a different variable, instead of letting a reference
          // to a mutable state variable survive across 2 different contexts.
          start = 0
          for laneID in 0..<8 {
            let allocationSize = dictionaryCount[laneID]
            guard allocationSize > 0 else {
              continue
            }
            
            let newPointer = bufferPointer.baseAddress.unsafelyUnwrapped + start
            start &+= allocationSize
            
            let newBufferPointer = UnsafeBufferPointer(
              start: newPointer, count: allocationSize)
            let key32 = UInt32(truncatingIfNeeded: laneID)
            let intOffset = (key32 &>> SIMD3(0, 1, 2)) & 1
            let floatOffset = SIMD3<Float>(intOffset) * 2 - 1
            let newOrigin = levelOrigin + floatOffset * levelSize / 2
            traverseGrid(
              atomIDs: newBufferPointer,
              levelOrigin: newOrigin,
              levelSize: levelSize / 2)
          }
        }
      }
      
      let levelOrigin = SIMD3<Float>(repeating: octreeStartSize)
      let initialArray = atoms.indices.map(UInt32.init(truncatingIfNeeded:))
      initialArray.withUnsafeBufferPointer { bufferPointer in
        traverseGrid(
          atomIDs: bufferPointer,
          levelOrigin: levelOrigin,
          levelSize: octreeStartSize)
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
          
          let key = Int((index &<< SIMD3(0, 1, 2)).wrappedSum())
          let previousCount = dictionaryCount[key]
          let pointer = dictionary.advanced(
            by: key * maxCellSize + previousCount)
          
          dictionaryCount[key] += 1
          pointer.pointee = atomID32
        }
        
        var temporaryAllocationSize = 0
        for laneID in 0..<8 {
          let allocationSize = dictionaryCount[laneID]
          temporaryAllocationSize += allocationSize
        }
        
        withUnsafeTemporaryAllocation(
          of: UInt32.self,
          capacity: temporaryAllocationSize
        ) { bufferPointer in
          var start = 0
          for laneID in 0..<8 {
            let allocationSize = dictionaryCount[laneID]
            guard allocationSize > 0 else {
              continue
            }
            
            let oldPointer = dictionary.advanced(by: laneID &* maxCellSize)
            let newPointer = bufferPointer.baseAddress.unsafelyUnwrapped + start
            newPointer.initialize(from: oldPointer, count: allocationSize)
            start &+= allocationSize
          }
          
          // TODO: Use a different variable, instead of letting a reference
          // to a mutable state variable survive across 2 different contexts.
          start = 0
          for laneID in 0..<8 {
            let allocationSize = dictionaryCount[laneID]
            guard allocationSize > 0 else {
              continue
            }
            
            let newPointer = bufferPointer.baseAddress.unsafelyUnwrapped + start
            start &+= allocationSize
            if allocationSize == 1 {
              output.append(newPointer.pointee)
              continue
            }
            
            let newBufferPointer = UnsafeBufferPointer(
              start: newPointer, count: allocationSize)
            let key32 = UInt32(truncatingIfNeeded: laneID)
            let intOffset = (key32 &>> SIMD3(0, 1, 2)) & 1
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
