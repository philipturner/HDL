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
  let size2: UInt32 = 32
  let size3: UInt32 = 128
  var lhsTransformed8: [SIMD8<Float>] = []
  var rhsTransformed8: [SIMD8<Float>] = []
  var lhsTransformed4: [SIMD4<Float>] = []
  var rhsTransformed4: [SIMD4<Float>] = []
  var outputMatches: [[SIMD2<Float>]] = []
  var paddedCutoffLHS: Float = .zero
  var paddedCutoffRHS: Float = .zero
  
  var lhs8: [SIMD4<Float>] = []
  var rhs8: [SIMD4<Float>] = []
  var lhs32: [SIMD4<Float>] = []
  var rhs32: [SIMD4<Float>] = []
  var lhs128: [SIMD4<Float>] = []
  var rhs128: [SIMD4<Float>] = []
  
  DispatchQueue.concurrentPerform(iterations: 3) { z in
    if z == 0 {
      (lhsTransformed8, lhsTransformed4, paddedCutoffLHS) =
      transform8(lhsAtoms, size: size3, algorithm: algorithm)
    } else if z == 1 {
      (rhsTransformed8, rhsTransformed4, paddedCutoffRHS) =
      transform8(rhsAtoms, size: size3, algorithm: algorithm)
    } else if z == 2 {
      for _ in 0..<UInt32(lhsAtoms.count + 7) / 8 * 8 {
        var array: [SIMD2<Float>] = []
        array.reserveCapacity(8)
        outputMatches.append(array)
      }
    }
    
    if z == 0 {
      lhs8 = blockBounds(lhsTransformed4, size: 8)
      lhs32 = blockBounds(lhsTransformed4, size: size2)
      lhs128 = blockBounds(lhsTransformed4, size: size3)
    } else if z == 1 {
      rhs8 = blockBounds(rhsTransformed4, size: 8)
      rhs32 = blockBounds(rhsTransformed4, size: size2)
      rhs128 = blockBounds(rhsTransformed4, size: size3)
    }
  }
  
  var outputArray: [UInt32] = []
  var outputRanges: [Range<Int>] = []
  var outputSlices: [ArraySlice<UInt32>] = []
  
  // MARK: - Hierarchy
  
  let paddedCutoff = paddedCutoffLHS + paddedCutoffRHS
  
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
      let lhsX = lhsTransformed8[pointer &+ 0]
      let lhsY = lhsTransformed8[pointer &+ 1]
      let lhsZ = lhsTransformed8[pointer &+ 2]
      let lhsR = lhsTransformed8[pointer &+ 3]
      let lhsBlockBound = lhs8[Int(vIDi)]
      
      for vIDj in loopStartJ..<loopEndJ {
        let rhsBlockBound = rhs8[Int(vIDj)]
        guard compareBlocks(
          lhsBlockBound, rhsBlockBound, paddedCutoff: paddedCutoff) else {
          continue
        }
        
        for lane in 0..<UInt32(8) {
          let atomID = vIDj &* 8 &+ lane
          let atom = rhsTransformed4[Int(atomID)]
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

func blockBounds(
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

func blockBounds3(
  _ transformed: [SIMD4<Float>],
  size: UInt32
) {
  // TODO: Swizzle into the format of 4 x SIMD8<Float> for the "transformed"
  // format before fixing up this function. It should reduce the overhead of
  // swizzling on-the-spot during the search function.
//  for blockID128 in 0..<UInt32(transformed.count) / 128 {
//    let blockID32Start = blockID128 &* 4
//    let blockID32End = blockID32Start &+ 4
//    for blockID32 in blockID32Start..<blockID32End {
//      let vIDStart = blockID32 &* 4
//      let vIDEnd = blockID32 &+ 4
//      for vID in vIDStart..<vIDEnd {
//        let atomIDStart = vID &* 8
//        let atomIDEnd = atomIDStart &+ 8
//      }
//    }
//  }
}

func transform(
  _ atoms: [Entity],
  size: UInt32,
  algorithm: Topology.MatchAlgorithm
) -> ([SIMD4<Float>], Float) {
  var transformed: [SIMD4<Float>] = []
  var paddedCutoff: Float = .zero
  for atomID in 0..<UInt32(atoms.count) {
    let atom = atoms[Int(atomID)]
    let atomicNumber = Int(atom.storage.w)
    var rhsR = Element.covalentRadii[atomicNumber]
    
    switch algorithm {
    case .absoluteRadius(let radius):
      rhsR = radius / 2
    case .covalentBondLength(let scale):
      rhsR *= scale
    }
    var rhs = atom.storage
    rhs.w = rhsR
    transformed.append(rhs)
    paddedCutoff = max(paddedCutoff, rhs.w)
  }
  
  // Pad to the granularity of blocks.
  let atomCount = UInt32(atoms.count)
  for _ in atomCount..<(atomCount &+ size &- 1)/size*size {
    let original = transformed[Int(atomCount &- 1)]
    transformed.append(original)
  }
  
  return (transformed, paddedCutoff)
}

// Start off by using the outputs of this function, just to substitute the
// swizzling part of the inner loop. See whether that is correct, and how it
// changes performance (2x overhead, but lower cost for the bulk of compute).
// Then, proceed with rewriting the functions for generating block bounds, and
// unrolling the inner loop of this function.
func transform8(
  _ atoms: [Entity],
  size: UInt32,
  algorithm: Topology.MatchAlgorithm
) -> ([SIMD8<Float>], [SIMD4<Float>], Float) {
  // This function has two parts: the vectorized part and the scalarized part.
  // The code for them could be de-duplicated by force-inlining and having
  // constant folding only for the first invocation. For debugging, start out
  // with a function that isn't modified to have force-inlined parts.
  //
  // One could even modify the function to generate block bounds on the fly.
  // However, that task should be addressed in a separate optimization.
  var output: [SIMD8<Float>] = []
  var expectedCapacity = (atoms.count + 127) / 128
  expectedCapacity *= 16 // 16 vectors per 128
  expectedCapacity *= 4 // 4 vectors interleaved in proximity
  output.reserveCapacity(Int(expectedCapacity))
  var paddedCutoff: Float = .zero
  
  var transformed4: [SIMD4<Float>] = []
  transformed4.reserveCapacity((atoms.count + 127) / 128 * 128)
  
  for vID in 0..<(UInt32(atoms.count) + 7) / 8 {
    var vecX: SIMD8<Float> = .zero
    var vecY: SIMD8<Float> = .zero
    var vecZ: SIMD8<Float> = .zero
    var vecR: SIMD8<Float> = .zero
    
    // Perhaps just make a simple branch here, instead of an inlined function.
    // If the branch is true, it activates the unrolled fast-path. Defer this
    // optimization until after the bulk of compute is multithreaded.
    let validLanes = min(8, UInt32(atoms.count) &- vID &* 8)
    for laneID in 0..<UInt32(validLanes) {
      let atom = atoms[Int(vID &* 8 &+ laneID)].storage
      vecX[Int(laneID)] = atom.x
      vecY[Int(laneID)] = atom.y
      vecZ[Int(laneID)] = atom.z
      
      let atomicNumber = Int(atom.w)
      vecR[Int(laneID)] = Element.covalentRadii[atomicNumber]
      transformed4.append(atom)
    }
    for laneID in UInt32(validLanes)..<8 {
      vecX[Int(laneID)] = vecX[Int(validLanes &- 1)]
      vecY[Int(laneID)] = vecY[Int(validLanes &- 1)]
      vecZ[Int(laneID)] = vecZ[Int(validLanes &- 1)]
      vecR[Int(laneID)] = vecR[Int(validLanes &- 1)]
      transformed4.append(transformed4.last!)
    }
    switch algorithm {
    case .absoluteRadius(let radius):
      vecR = SIMD8(repeating: radius / 2)
    case .covalentBondLength(let scale):
      vecR *= scale
    }
    for laneID in 0..<UInt32(8) {
      transformed4[Int(vID &* 8 &+ laneID)].w = vecR[Int(laneID)]
    }
    
    output.append(vecX)
    output.append(vecY)
    output.append(vecZ)
    output.append(vecR)
    paddedCutoff = max(paddedCutoff, vecR.max())
  }
  
  do {
    let vID8End = (UInt32(atoms.count) + 7) / 8
    let vID128End = (UInt32(atoms.count) + 127) / 128
    let vecX = output[Int(vID8End &* 4 &- 4)]
    let vecY = output[Int(vID8End &* 4 &- 3)]
    let vecZ = output[Int(vID8End &* 4 &- 2)]
    let vecR = output[Int(vID8End &* 4 &- 1)]
    for _ in vID8End..<vID128End * 16 {
      output.append(vecX)
      output.append(vecY)
      output.append(vecZ)
      output.append(vecR)
    }
    let last = transformed4.last ?? .zero
    for _ in vID8End * 8..<vID128End * 128 {
      transformed4.append(last)
    }
  }
  
  return (output, transformed4, paddedCutoff)
}

@inline(__always)
func compareBlocks(
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

