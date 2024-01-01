//
//  Match.swift
//
//
//  Created by Philip Turner on 12/23/23.
//

import Dispatch

extension Topology {
  public enum MatchAlgorithm {
    /// Search for neighbors within a fixed radius (in nanometers).
    case absoluteRadius(Float)
    
    /// Search for neighbors within a multiple of covalent bond length.
    case covalentBondLength(Float)
  }
  
  public func match(
    _ input: [Entity],
    algorithm: MatchAlgorithm = .covalentBondLength(1.5)
  ) -> [ArraySlice<UInt32>] {
    // Try creating a cutoff. For very small structures, sorting may
    // harm performance more than it helps. Multithreading also harms
    // performance for small-enough problem sizes and must be selectively
    // disabled. It would be simple to have just 2 ensemble members:
    // - single-threaded, no sorting
    // - multi-threaded, sorting
    #if false
    return matchImpl(lhs: input, rhs: atoms, algorithm: algorithm)
    
    #else
    var lhsReordering: [UInt32] = []
    var rhsReordering: [UInt32] = []
    DispatchQueue.concurrentPerform(iterations: 2) { z in
      if z == 0 {
        let lhsGrid = GridSorter(atoms: input)
        lhsReordering = lhsGrid.mortonReordering()
      } else {
        let rhsGrid = GridSorter(atoms: atoms)
        rhsReordering = rhsGrid.mortonReordering()
      }
    }
    
    let invalidEntity = Entity(storage: .init(repeating: -1))
    var lhs = Array(repeating: invalidEntity, count: input.count)
    var rhs = Array(repeating: invalidEntity, count: atoms.count)
    var rhsReorderedMap = [UInt32](repeating: .max, count: rhsReordering.count)
    for i in lhsReordering.indices {
      let reordered = Int(truncatingIfNeeded: lhsReordering[i])
      let entity = input[i]
      lhs[reordered] = entity
    }
    for i in rhsReordering.indices {
      let reordered = Int(truncatingIfNeeded: rhsReordering[i])
      let entity = atoms[i]
      rhs[reordered] = entity
      rhsReorderedMap[reordered] = UInt32(truncatingIfNeeded: i)
    }
    precondition(
      !rhsReorderedMap.contains(.max), "RHS indices were mapped incorrectly.")
    
    // Call the actual matching function.
    let slicesReordered = matchImpl(lhs: lhs, rhs: rhs, algorithm: algorithm)
    precondition(
      input.count == slicesReordered.count, "Incorrect number of match slices.")
    
    var outputArray: [UInt32] = []
    var outputRanges: [Range<Int>] = []
    var outputSlices: [ArraySlice<UInt32>] = []
    for lhsID in input.indices {
      let lhsReordered = Int(truncatingIfNeeded: lhsReordering[lhsID])
      let sliceReordered = slicesReordered[lhsReordered]
      
      let rangeStart = outputArray.count
      for rhsReorderedID32 in sliceReordered {
        let rhsReorderedID = Int(truncatingIfNeeded: rhsReorderedID32)
        let rhsID = rhsReorderedMap[rhsReorderedID]
        outputArray.append(rhsID)
      }
      let rangeEnd = outputArray.count
      outputRanges.append(rangeStart..<rangeEnd)
    }
    for range in outputRanges {
      let slice = outputArray[range]
      outputSlices.append(slice)
    }
    
    // WARNING: This code assumes the array slices are contiguous.
    guard input.count > 0 else {
      return []
    }
    precondition(slicesReordered.first!.startIndex == 0)
    precondition(
      outputArray.count == slicesReordered.last!.endIndex,
      "Output indices were mapped incorrectly.")
    return outputSlices
    #endif
  }
}

private func matchImpl(
  lhs lhsAtoms: [Entity],
  rhs rhsAtoms: [Entity],
  algorithm: Topology.MatchAlgorithm
) -> [ArraySlice<UInt32>] {
  // MARK: - Prepare LHS and RHS
  
  var lhs: Operand = .init()
  var rhs: Operand = .init()
  var outputMatches: [[SIMD2<Float>]] = []
  
  DispatchQueue.concurrentPerform(iterations: 3) { z in
    if z == 0 {
      lhs = transform8(lhsAtoms, algorithm: algorithm)
    } else if z == 1 {
      rhs = transform8(rhsAtoms, algorithm: algorithm)
    } else if z == 2 {
      for _ in 0..<UInt32(lhsAtoms.count + 7) / 8 * 8 {
        var array: [SIMD2<Float>] = []
        array.reserveCapacity(8)
        outputMatches.append(array)
      }
    }
  }
  
  let paddedCutoff = lhs.paddedCutoff + rhs.paddedCutoff
  
  var outputArray: [UInt32] = []
  var outputRanges: [Range<Int>] = []
  var outputSlices: [ArraySlice<UInt32>] = []
  
  // MARK: - Hierarchy
  
  let loopStartI: UInt32 = 0
  var loopEndI = loopStartI + UInt32(lhsAtoms.count * 2) + 256
  loopEndI = min(loopEndI, UInt32(lhsAtoms.count + 127) / 128)
  
  let loopStartJ: UInt32 = 0
  var loopEndJ = loopStartJ + UInt32(rhsAtoms.count * 2) + 256
  loopEndJ = min(loopEndJ, UInt32(rhsAtoms.count + 127) / 128)
  
  for vIDi3 in loopStartI..<loopEndI {
    for vIDj3 in loopStartJ..<loopEndJ {
      let lhsBlockBound = lhs.blockBounds128[Int(vIDi3)]
      let rhsBlockBound = rhs.blockBounds128[Int(vIDj3)]
      let mask = compareBlocks(
        lhsBlockBound, rhsBlockBound, paddedCutoff: paddedCutoff)
      if mask {
        innerLoop2(vIDi3, vIDj3)
      }
    }
  }
  
  func innerLoop2(_ vIDi3: UInt32, _ vIDj3: UInt32) {
    let loopStartI = vIDi3 * 4
    var loopEndI = loopStartI + 4
    loopEndI = min(loopEndI, UInt32(lhsAtoms.count + 31) / 32)
    
    let loopStartJ = vIDj3 * 4
    var loopEndJ = loopStartJ + 4
    loopEndJ = min(loopEndJ, UInt32(rhsAtoms.count + 31) / 32)
    
    for vIDi2 in loopStartI..<loopEndI {
      for vIDj2 in loopStartJ..<loopEndJ {
        let lhsBlockBound = lhs.blockBounds32[Int(vIDi2)]
        let rhsBlockBound = rhs.blockBounds32[Int(vIDj2)]
        let mask = compareBlocks(
          lhsBlockBound, rhsBlockBound, paddedCutoff: paddedCutoff)
        if mask {
          innerLoop1(vIDi2, vIDj2)
        }
      }
    }
  }
  
  // MARK: - Search
  
  func innerLoop1(_ vIDi2: UInt32, _ vIDj2: UInt32) {
    let loopStartI = vIDi2 * 4
    var loopEndI = loopStartI + 4
    loopEndI = min(loopEndI, UInt32(lhsAtoms.count + 7) / 8)
    
    let loopStartJ = vIDj2 * 4
    var loopEndJ = loopStartJ + 4
    loopEndJ = min(loopEndJ, UInt32(rhsAtoms.count + 7) / 8)
    
    for vIDi in loopStartI..<loopEndI {
      let pointer = Int(vIDi &* 4)
      let lhsX = lhs.transformed8[pointer &+ 0]
      let lhsY = lhs.transformed8[pointer &+ 1]
      let lhsZ = lhs.transformed8[pointer &+ 2]
      let lhsR = lhs.transformed8[pointer &+ 3]
      let lhsBlockBound = lhs.blockBounds8[Int(vIDi)]
      
      for vIDj in loopStartJ..<loopEndJ {
        let rhsBlockBound = rhs.blockBounds8[Int(vIDj)]
        guard compareBlocks(
          lhsBlockBound, rhsBlockBound, paddedCutoff: paddedCutoff) else {
          continue
        }
        
        for lane in 0..<UInt32(8) {
          let atomID = vIDj &* 8 &+ lane
          let atom = rhs.transformed4[Int(atomID)]
          let deltaX = lhsX - atom.x
          let deltaY = lhsY - atom.y
          let deltaZ = lhsZ - atom.z
          let deltaSquared = deltaX * deltaX + deltaY * deltaY + deltaZ * deltaZ
          
          let r = lhsR + atom.w
          let mask = deltaSquared .<= r * r
          if any(mask) {
            for laneID in 0..<UInt32(8) {
              let keyValuePair = SIMD2(
                Float(bitPattern: atomID), deltaSquared[Int(laneID)])
              if mask[Int(laneID)] {
                outputMatches[Int(vIDi &* 8 &+ laneID)].append(keyValuePair)
              }
            }
          }
        }
      }
    }
  }
  
  // MARK: - Sort
  
  for laneID in lhsAtoms.indices {
    // Remove bogus matches that result from padding the RHS.
    while let last = outputMatches[laneID].last,
          last[0].bitPattern >= rhsAtoms.count {
      outputMatches[laneID].removeLast()
    }
    
    // Sort the matches in ascending order of distance.
    outputMatches[laneID].sort { $0.y < $1.y }
    
    let rangeStart = outputArray.count
    for match in outputMatches[laneID] {
      let atomID = match[0].bitPattern
      outputArray.append(atomID)
    }
    let rangeEnd = outputArray.count
    outputRanges.append(rangeStart..<rangeEnd)
  }
  
  for range in outputRanges {
    let slice = outputArray[range]
    outputSlices.append(slice)
  }
  
  return outputSlices
}

// MARK: - Block Comparison

private struct Operand {
  var paddedCutoff: Float = .zero
  var transformed4: [SIMD4<Float>] = []
  var transformed8: [SIMD8<Float>] = []
  var blockBounds8: [SIMD4<Float>] = []
  var blockBounds32: [SIMD4<Float>] = []
  var blockBounds128: [SIMD4<Float>] = []
}

private func transform8(
  _ atoms: [Entity],
  algorithm: Topology.MatchAlgorithm
) -> Operand {
  var output = Operand()
  let paddedAtomCount = (atoms.count + 127) / 128 * 128
  output.transformed8.reserveCapacity(4 * paddedAtomCount / 8)
  output.transformed4.reserveCapacity(paddedAtomCount)
  output.blockBounds8.reserveCapacity(paddedAtomCount / 8)
  output.blockBounds32.reserveCapacity(paddedAtomCount / 32)
  output.blockBounds128.reserveCapacity(paddedAtomCount / 128)
  
  // Accumulate data at the 32-atom granularity. At the same time, generate most
  // of the other per-operand data.
  var currentBlockBox32: SIMD8<Float> = .zero
  var blockBoxes32: [SIMD8<Float>] = []
  blockBoxes32.reserveCapacity(paddedAtomCount / 32)
  
  let vIDEnd = (UInt32(atoms.count) + 7) / 8
  for vID in 0..<vIDEnd {
    var vecX: SIMD8<Float> = .zero
    var vecY: SIMD8<Float> = .zero
    var vecZ: SIMD8<Float> = .zero
    var vecR: SIMD8<Float> = .zero
    
    if vID &* 8 &+ 8 < UInt32(truncatingIfNeeded: atoms.count) {
      for laneID in 0..<UInt32(8) {
        let atom = atoms[Int(vID &* 8 &+ laneID)].storage
        vecX[Int(laneID)] = atom.x
        vecY[Int(laneID)] = atom.y
        vecZ[Int(laneID)] = atom.z
        
        let atomicNumber = Int(atom.w)
        vecR[Int(laneID)] = Element.covalentRadii[atomicNumber]
        output.transformed4.append(atom)
      }
    } else {
      let validLanes = UInt32(truncatingIfNeeded: atoms.count) &- vID &* 8
      for laneID in 0..<validLanes {
        let atom = atoms[Int(vID &* 8 &+ laneID)].storage
        vecX[Int(laneID)] = atom.x
        vecY[Int(laneID)] = atom.y
        vecZ[Int(laneID)] = atom.z
        
        let atomicNumber = Int(atom.w)
        vecR[Int(laneID)] = Element.covalentRadii[atomicNumber]
        output.transformed4.append(atom)
      }
      let last = output.transformed4.last ?? .zero
      for laneID in validLanes..<8 {
        vecX[Int(laneID)] = vecX[Int(validLanes &- 1)]
        vecY[Int(laneID)] = vecY[Int(validLanes &- 1)]
        vecZ[Int(laneID)] = vecZ[Int(validLanes &- 1)]
        vecR[Int(laneID)] = vecR[Int(validLanes &- 1)]
        output.transformed4.append(last)
      }
    }
    switch algorithm {
    case .absoluteRadius(let radius):
      vecR = SIMD8(repeating: radius / 2)
    case .covalentBondLength(let scale):
      vecR *= scale
    }
    for laneID in 0..<UInt32(8) {
      output.transformed4[Int(vID &* 8 &+ laneID)].w = vecR[Int(laneID)]
    }
    
    output.transformed8.append(vecX)
    output.transformed8.append(vecY)
    output.transformed8.append(vecZ)
    output.transformed8.append(vecR)
    output.paddedCutoff = max(output.paddedCutoff, vecR.max())
    
    let minimum = SIMD4(vecX.min(), vecY.min(), vecZ.min(), 0)
    let maximum = SIMD4(vecX.max(), vecY.max(), vecZ.max(), 0)
    var bounds = (minimum + maximum) / 2
    
    // Find the maximum deviation from the center.
    let deltaX = vecX - bounds.x
    let deltaY = vecY - bounds.y
    let deltaZ = vecZ - bounds.z
    let deltaSq = deltaX * deltaX + deltaY * deltaY + deltaZ * deltaZ
    bounds.w = deltaSq.max().squareRoot()
    output.blockBounds8.append(bounds)
    
    if (vID % 4 == 0) {
      currentBlockBox32 = SIMD8(lowHalf: minimum, highHalf: maximum)
    } else {
      currentBlockBox32.lowHalf.replace(
        with: minimum, where: minimum .< currentBlockBox32.lowHalf)
      currentBlockBox32.highHalf.replace(
        with: maximum, where: maximum .> currentBlockBox32.highHalf)
    }
    if (vID % 4 == 3) || (vID == vIDEnd &- 1) {
      let bounds = (currentBlockBox32.lowHalf + currentBlockBox32.highHalf) / 2
      output.blockBounds32.append(bounds)
      blockBoxes32.append(currentBlockBox32)
    }
  }
  
  // Pad the smaller arrays to the granularity of 128 atoms.
  do {
    let vID8End = (UInt32(atoms.count) + 7) / 8
    let vID32End = (UInt32(atoms.count) + 31) / 32
    let vID128End = (UInt32(atoms.count) + 127) / 128
    
    let last4 = output.transformed4.last ?? .zero
    for _ in vID8End * 8..<vID128End * 128 {
      output.transformed4.append(last4)
    }
    
    let vecX = output.transformed8[Int(vID8End &* 4 &- 4)]
    let vecY = output.transformed8[Int(vID8End &* 4 &- 3)]
    let vecZ = output.transformed8[Int(vID8End &* 4 &- 2)]
    let vecR = output.transformed8[Int(vID8End &* 4 &- 1)]
    let bound8 = output.blockBounds8[Int(vID8End &- 1)]
    for _ in vID8End..<vID128End * 16 {
      output.transformed8.append(vecX)
      output.transformed8.append(vecY)
      output.transformed8.append(vecZ)
      output.transformed8.append(vecR)
      output.blockBounds8.append(bound8)
    }
    
    let box32 = blockBoxes32.last ?? .zero
    let bound32 = output.blockBounds32.last ?? .zero
    for _ in vID32End..<vID128End * 4 {
      blockBoxes32.append(box32)
      output.blockBounds32.append(bound32)
    }
  }
  
  // Generate the centers at 128-atom granularity in a standalone pass. Then,
  // the array for 'blockBoxes32' should be released.
  for vID128 in 0..<UInt32(paddedAtomCount / 128) {
    var blockBox128: SIMD8<Float> = .zero
    blockBox128 = blockBoxes32[Int(vID128 &* 4)]
    
    for laneID in 1..<UInt32(4) {
      let vID32 = vID128 &* 4 &+ laneID
      let blockBox32 = blockBoxes32[Int(vID32)]
      let (minimum, maximum) = (blockBox32.lowHalf, blockBox32.highHalf)
      blockBox128.lowHalf.replace(
        with: minimum, where: minimum .< blockBox128.lowHalf)
      blockBox128.highHalf.replace(
        with: maximum, where: maximum .> blockBox128.highHalf)
    }
    
    let bounds = (blockBox128.lowHalf + blockBox128.highHalf) / 2
    output.blockBounds128.append(bounds)
  }
  
  // Pre-compute the distance of each atom from the center.
  for vID128 in 0..<UInt32(paddedAtomCount / 128) {
    let center128 = output.blockBounds128[Int(vID128)]
    var distance32: SIMD4<Float> = .zero
    var distance128: SIMD4<Float> = .zero
    
    for laneID in 0..<UInt32(4) {
      let vID32 = vID128 &* 4 &+ laneID
      let center32 = output.blockBounds32[Int(vID32)]
      var distance32Cache: SIMD4<Float> = .zero
      var distance128Cache: SIMD4<Float> = .zero
      
      for laneID in 0..<UInt32(4) {
        let vID8 = vID32 &* 4 &+ laneID
        let vecX = output.transformed8[Int(vID8 &* 4 &+ 0)]
        let vecY = output.transformed8[Int(vID8 &* 4 &+ 1)]
        let vecZ = output.transformed8[Int(vID8 &* 4 &+ 2)]
        
        // Find the maximum deviation from the center.
        @inline(__always)
        func createDistanceSq(_ center: SIMD4<Float>) -> Float {
          let deltaX = vecX - center.x
          let deltaY = vecY - center.y
          let deltaZ = vecZ - center.z
          let deltaSq = deltaX * deltaX + deltaY * deltaY + deltaZ * deltaZ
          return deltaSq.max()
        }
        distance32Cache[Int(laneID)] = createDistanceSq(center32)
        distance128Cache[Int(laneID)] = createDistanceSq(center128)
      }
      distance32[Int(laneID)] = distance32Cache.max()
      distance128[Int(laneID)] = distance128Cache.max()
    }
    
    distance32.formSquareRoot()
    distance128.formSquareRoot()
    for laneID in 0..<UInt32(4) {
      let vID32 = vID128 &* 4 &+ laneID
      let distance = distance32[Int(laneID)]
      output.blockBounds32[Int(vID32)].w = distance
    }
    output.blockBounds128[Int(vID128)].w = distance128.max()
  }
  
  return output
}

@inline(__always)
private func compareBlocks(
  _ lhsBlockBound: SIMD4<Float>,
  _ rhsBlockBound: SIMD4<Float>,
  paddedCutoff: Float
) -> Bool {
  // Code from the OpenMM GPU shader for finding interacting blocks.
  /*
   real4 blockCenterY = sortedBlockCenter[block2];
//                real3 blockSizeY = sortedBlockBoundingBox[block2];
   real3 blockSizeY = real3(vloada_half3(block2, sortedBlockBoundingBox));
   real4 blockDelta = blockCenterX-blockCenterY;
#ifdef USE_PERIODIC
   APPLY_PERIODIC_TO_DELTA(blockDelta)
#endif
   includeBlock2 &= (blockDelta.x*blockDelta.x+blockDelta.y*blockDelta.y+blockDelta.z*blockDelta.z < (PADDED_CUTOFF+blockCenterX.w+blockCenterY.w)*(PADDED_CUTOFF+blockCenterX.w+blockCenterY.w));
   blockDelta.x = max((real) 0, fabs(blockDelta.x)-blockSizeX.x-blockSizeY.x);
   blockDelta.y = max((real) 0, fabs(blockDelta.y)-blockSizeX.y-blockSizeY.y);
   blockDelta.z = max((real) 0, fabs(blockDelta.z)-blockSizeX.z-blockSizeY.z);
   includeBlock2 &= (blockDelta.x*blockDelta.x+blockDelta.y*blockDelta.y+blockDelta.z*blockDelta.z < PADDED_CUTOFF_SQUARED);
   */
  
  let blockCenterX = lhsBlockBound
  let blockCenterY = rhsBlockBound
  let largeCutoff = paddedCutoff + blockCenterX.w + blockCenterY.w
  let largeCutoffSquared = largeCutoff * largeCutoff
  
  var blockDelta = blockCenterX - blockCenterY
  blockDelta.w = 0
  return (blockDelta * blockDelta).sum() < largeCutoffSquared
}
