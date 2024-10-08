//
//  Sort.swift
//  HDL
//
//  Created by Philip Turner on 12/2/23.
//

import Foundation

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
    
    if atoms.count == 0 {
      origin = .zero
      dimensions = .init(repeating: 1)
    } else {
      var minimum = SIMD3<Float>(repeating: .greatestFiniteMagnitude)
      var maximum = -minimum
      for atom in atoms {
        let position = atom.position
        minimum.replace(with: position, where: position .< minimum)
        maximum.replace(with: position, where: position .> maximum)
      }
      for atom in atoms {
        let position = atom.position
        minimum.replace(with: position, where: position .< minimum)
      }
      
      origin = minimum
      dimensions = maximum - minimum
      dimensions.replace(with: .init(repeating: 0.5), where: dimensions .< 0.5)
    }
  }
}

// MARK: - Morton Reordering

extension GridSorter {
  func invertOrder(_ input: [UInt32]) -> [UInt32] {
    if input.count == 0 {
      return []
    }
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
          var index: SIMD3<UInt32> = .init(repeating: 1)
          let atomID = Int(truncatingIfNeeded: atomID32)
          let atomPosition = atoms[atomID].position - self.origin
          index.replace(
            with: .init(repeating: 0),
            where: atomPosition .< levelOrigin)
          
          let key = Int(
            truncatingIfNeeded: (index &<< SIMD3(0, 1, 2)).wrappedSum())
          let previousCount = dictionaryCount[key]
          let pointer = dictionary.advanced(
            by: key &* atoms.count &+ previousCount)
          
          dictionaryCount[key] &+= 1
          pointer.pointee = atomID32
        }
        
        var temporaryAllocationSize = 0
        for laneID in 0..<8 {
          let allocationSize = dictionaryCount[laneID]
          temporaryAllocationSize &+= allocationSize
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
          var index: SIMD3<UInt32> = .init(repeating: 1)
          let atomID = Int(truncatingIfNeeded: atomID32)
          let atomPosition = atoms[atomID].position - self.origin
          index.replace(
            with: .init(repeating: 0),
            where: atomPosition .< levelOrigin)
          
          let key = Int(
            truncatingIfNeeded: (index &<< SIMD3(0, 1, 2)).wrappedSum())
          let previousCount = dictionaryCount[key]
          let pointer = dictionary.advanced(
            by: key &* maxCellSize &+ previousCount)
          
          dictionaryCount[key] &+= 1
          pointer.pointee = atomID32
        }
        
        var temporaryAllocationSize = 0
        for laneID in 0..<8 {
          let allocationSize = dictionaryCount[laneID]
          temporaryAllocationSize &+= allocationSize
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
