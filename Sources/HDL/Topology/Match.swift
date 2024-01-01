//
//  Match.swift
//
//
//  Created by Philip Turner on 12/23/23.
//

import Dispatch

// Goals for optimizing this:
//
// Higher levels of the hierarchy: recursive bounding box-bounding box test,
// also with search radius in the same manner as OpenMM PADDED_CUTOFF. Store
// block position as 3 x SIMD8<Float>, block bounds as 3 x SIMD8<Half>, block
// radii as SIMD8<Half>.
// - use FloatingPoint.nextUp to round toward positive infinity
// - make the two halves of each SIMD8 independent, so they can be paged into
//   ISA registers without spilling

extension Topology {
  public enum MatchAlgorithm {
    // Search for neighbors within a fixed radius (in nanometers).
    case absoluteRadius(Float)
    
    // Search for neighbors within a multiple of covalent bond length.
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

// The lhs and rhs must already be sorted in Morton order.
private func matchImpl(
  lhs lhsAtoms: [Entity],
  rhs rhsAtoms: [Entity],
  algorithm: Topology.MatchAlgorithm
) -> [ArraySlice<UInt32>] {
  // MARK: - Prepare LHS and RHS
  
  // A future optimization could change the number of hierarchy levels with the
  // problem size. Perhaps use a ratio over the smallest problem dimension.
  var lhs: Operand = .init()
  var rhs: Operand = .init()
  var outputMatches: [[SIMD2<Float>]] = []
  
//  var lhs8: [SIMD4<Float>] = []
//  var rhs8: [SIMD4<Float>] = []
  var lhs32: [SIMD4<Float>] = []
  var rhs32: [SIMD4<Float>] = []
  var lhs128: [SIMD4<Float>] = []
  var rhs128: [SIMD4<Float>] = []
  
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
    
    if z == 0 {
      lhs32 = blockBounds(lhs.transformed4, size: 32)
      lhs128 = blockBounds(lhs.transformed4, size: 128)
    } else if z == 1 {
      rhs32 = blockBounds(rhs.transformed4, size: 32)
      rhs128 = blockBounds(rhs.transformed4, size: 128)
    }
  }
  
  let paddedCutoff = lhs.paddedCutoff + rhs.paddedCutoff
  
  var outputArray: [UInt32] = []
  var outputRanges: [Range<Int>] = []
  var outputSlices: [ArraySlice<UInt32>] = []
  
  // MARK: - Hierarchy
  
  let loopStartI: UInt32 = 0
  var loopEndI = loopStartI + UInt32(lhsAtoms.count * 2)
  loopEndI = min(loopEndI, UInt32(lhsAtoms.count + 127) / 128)
  
  let loopStartJ: UInt32 = 0
  var loopEndJ = loopStartJ + UInt32(rhsAtoms.count * 2)
  loopEndJ = min(loopEndJ, UInt32(rhsAtoms.count + 127) / 128)
  
  for vIDi3 in loopStartI..<loopEndI {
    for vIDj3 in loopStartJ..<loopEndJ {
      let lhsBlockBound = lhs128[Int(vIDi3)]
      let rhsBlockBound = rhs128[Int(vIDj3)]
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
        let lhsBlockBound = lhs32[Int(vIDi2)]
        let rhsBlockBound = rhs32[Int(vIDj2)]
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
            // If this part is a bottleneck, there may be clever ways to improve
            // the parallelism of conditional writes to a compacted array. 8
            // arrays can be generated concurrently, using zero control flow.
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

private func blockBounds(
  _ transformed: [SIMD4<Float>],
  size: UInt32
) -> [SIMD4<Float>] {
  var blockBounds: [SIMD4<Float>] = []
  // some scratch arrays for storing 16 partial reductions in local memory
  
  for vID in 0..<UInt32(transformed.count) / size {
    var minimum = SIMD4<Float>(repeating: .greatestFiniteMagnitude)
    var maximum = -minimum
    for lane in 0..<size {
      let atom = transformed[Int(vID &* size &+ lane)]
      minimum.replace(with: atom, where: atom .< minimum)
      maximum.replace(with: atom, where: atom .> maximum)
    }
    var bounds = (minimum + maximum) / 2
    
    // Write the maximum deviation^2 from the center to bounds[3]
    var maxCenterDeviationSq: Float = .zero
    for lane in 0..<size {
      let atom = transformed[Int(vID &* size &+ lane)]
      var delta = atom - bounds
      delta.w = 0
      
      let deviationSq = (delta * delta).sum()
      maxCenterDeviationSq = max(maxCenterDeviationSq, deviationSq)
    }
    bounds[3] = maxCenterDeviationSq.squareRoot()
    blockBounds.append(bounds)
  }
  return blockBounds
}

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
  // One could even modify the function to generate block bounds on the fly.
  // However, that task should be addressed in a separate optimization.
  var output = Operand()
  let paddedAtomCount = (atoms.count + 127) / 128 * 128
  output.transformed8.reserveCapacity(4 * paddedAtomCount / 8)
  output.transformed4.reserveCapacity(paddedAtomCount)
  output.blockBounds8.reserveCapacity(paddedAtomCount / 8)
  
  for vID in 0..<(UInt32(atoms.count) + 7) / 8 {
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
    
    let minimum = SIMD3(vecX.min(), vecY.min(), vecZ.min())
    let maximum = SIMD3(vecX.max(), vecY.max(), vecZ.max())
    var bounds = SIMD4((minimum + maximum) / 2, .zero)
    
    // Find the maximum deviation from the center.
    let deltaX = vecX - bounds.x
    let deltaY = vecY - bounds.y
    let deltaZ = vecZ - bounds.z
    let deltaSq = deltaX * deltaX + deltaY * deltaY + deltaZ * deltaZ
    bounds.w = deltaSq.max().squareRoot()
    output.blockBounds8.append(bounds)
  }
  
  do {
    let vID8End = (UInt32(atoms.count) + 7) / 8
    let vID128End = (UInt32(atoms.count) + 127) / 128
    let vecX = output.transformed8[Int(vID8End &* 4 &- 4)]
    let vecY = output.transformed8[Int(vID8End &* 4 &- 3)]
    let vecZ = output.transformed8[Int(vID8End &* 4 &- 2)]
    let vecR = output.transformed8[Int(vID8End &* 4 &- 1)]
    let bounds = output.blockBounds8[Int(vID8End &- 1)]
    
    for _ in vID8End..<vID128End * 16 {
      output.transformed8.append(vecX)
      output.transformed8.append(vecY)
      output.transformed8.append(vecZ)
      output.transformed8.append(vecR)
      output.blockBounds8.append(bounds)
    }
    let last = output.transformed4.last ?? .zero
    for _ in vID8End * 8..<vID128End * 128 {
      output.transformed4.append(last)
    }
  }
  
  // Start by merging centers from this block bounds finder with the existing
  // one. Assert that both results are the exact same, and print out a few pairs
  // of old + new centers.
  //
  // Then, proceed with pre-computing the distance of each atom from the center.
  output.blockBounds32.reserveCapacity(paddedAtomCount / 32)
  output.blockBounds128.reserveCapacity(paddedAtomCount / 128)
  for vID128 in 0..<UInt32(paddedAtomCount / 8) {
    // We lose information about the actual block bounds, instead going with a
    // sphere around a block center. This provides a more conservative estimate
    // of the bounding box. However, the final number - maximum atom distance -
    // is computed using the information of every atom.
    //
    // Using SIMD8<Float>, we store the full bounding box at the 32-atom
    // granularity. Thus, the bounding box information is only destroyed and
    // recreated during the vec8 -> vec32 stage.
    //
    // TODO: - Better idea - keep a running sum of the current 32-element
    // bounding box while making the first pass through the atoms. The avoids
    // the destruction of information. It also allows for the idea of debugging
    // - checking that centers are exactly the same before and after the code
    // change.
//    withUnsafeTemporaryAllocation(of: SIMD8<Float>.self, capacity: 4) {
//      let blockBounds32 = $0.baseAddress.unsafelyUnwrapped
//      
//      let vID32Start = vID128 &* 4
//      for lane32 in 0..<UInt32(4) {
//        let vID32 = vID32Start &+ lane32
//        var min32 = SIMD4<Float>(repeating: .greatestFiniteMagnitude)
//        var max32 = -min32
//        
//        let vID8Start = vID32 &* 4
//        for lane8 in 0..<UInt32(4) {
//          let vID8 = vID8Start &+ lane8
//          let bound8 = output.blockBounds8[Int(vID8)]
//          let min8 = bound8
//        }
//      }
//    }
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

