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
    // There is already some inherent order to atoms, since they originate from
    // spatially local and semi-Z-curve unit cells. Therefore, sorting should
    // be explored only as an optimization over an O(n)/O(nlogn) algorithm. We
    // should quantify how much it speeds up performance. Try different test
    // cases with truly random shuffling of data.
    //
    // You could also create a cutoff. For very small structures, sorting may
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
  let size2 = 32
  let lhsTransformed = transform(lhsAtoms, size: size2, algorithm: algorithm)
  let rhsTransformed = transform(rhsAtoms, size: size2, algorithm: algorithm)
  let lhsBlockBounds = blockBounds(lhsTransformed, size: 8)
  let rhsBlockBounds = blockBounds(rhsTransformed, size: 8)
  let lhs32 = blockBounds(lhsTransformed, size: size2)
  let rhs32 = blockBounds(rhsTransformed, size: size2)
  
  var mask32: [[Bool]] = []
  for i in 0..<(lhsAtoms.count + size2 - 1) / size2 {
    var row: [Bool] = []
    for j in 0..<(rhsAtoms.count + size2 - 1) / size2 {
      let lhsBlockBound = lhs32[i]
      let rhsBlockBound = rhs32[j]
      row.append(compareBlocks(
        lhsBlockBound, rhsBlockBound, executeCondition2: false))
    }
    mask32.append(row)
  }
  
  let lhs = (transformed: lhsTransformed, blockBounds: lhsBlockBounds)
  let rhs = (transformed: rhsTransformed, blockBounds: rhsBlockBounds)
  var outputArray: [UInt32] = []
  var outputRanges: [Range<Int>] = []
  var outputSlices: [ArraySlice<UInt32>] = []
  
  let vectorCount = (lhsAtoms.count + 7) / 8
  for vIDi in 0..<vectorCount {
    var lhsX: SIMD8<Float> = .zero
    var lhsY: SIMD8<Float> = .zero
    var lhsZ: SIMD8<Float> = .zero
    var lhsR: SIMD8<Float> = .zero
    for laneID in 0..<8 {
      let atom = lhs.transformed[vIDi &* 8 &+ laneID]
      lhsX[laneID] = atom.x
      lhsY[laneID] = atom.y
      lhsZ[laneID] = atom.z
      lhsR[laneID] = atom.w
    }
    let lhsBlockBound = lhs.blockBounds[vIDi]
    
    // MARK: - Search
    
    var matches: [[SIMD2<Float>]] = []
    for _ in 0..<8 {
      var array: [SIMD2<Float>] = []
      array.reserveCapacity(8)
      matches.append(array)
    }
    
    var vIDj = rhs.blockBounds.startIndex
    while vIDj < rhs.blockBounds.endIndex {
      if vIDj % (size2/8) == 0 {
        let mask8 = mask32[vIDi / (size2/8)][vIDj / (size2/8)]
        guard mask8 else {
          vIDj += size2/8
          continue
        }
      }
      defer { vIDj += 1 }
      
      let rhsBlockBound = rhs.blockBounds[vIDj]
      guard compareBlocks(lhsBlockBound, rhsBlockBound) else {
        continue
      }
      
      for lane in 0..<8 {
        let atomID = vIDj &* 8 &+ lane
        let atom = rhs.transformed[atomID]
        let deltaX = lhsX - atom.x
        let deltaY = lhsY - atom.y
        let deltaZ = lhsZ - atom.z
        let deltaSquared = deltaX * deltaX + deltaY * deltaY + deltaZ * deltaZ
        
        let r = lhsR + atom.w
        let mask = deltaSquared .<= r * r
        if any(mask) {
          let atomID32 = UInt32(truncatingIfNeeded: atomID)
          
          // If this part is a bottleneck, there may be clever ways to improve
          // the parallelism of conditional writes to a compacted array. 8
          // arrays can be generated concurrently, using zero control flow.
          for laneID in 0..<8 {
            let keyValuePair = SIMD2(
              Float(bitPattern: atomID32),
              Float(deltaSquared[laneID]))
            if mask[laneID] {
              matches[laneID].append(keyValuePair)
            }
          }
        }
      }
    }
    
    // MARK: - Sort
    
    let validLanes = min(8, lhsAtoms.count - vIDi * 8)
    for laneID in 0..<validLanes {
      // Remove bogus matches that result from padding the RHS.
      while let last = matches[laneID].last,
            last[0].bitPattern >= rhsAtoms.count {
        matches[laneID].removeLast()
      }
      
      // Sort the matches in ascending order of distance.
      matches[laneID].sort { $0.y < $1.y }
      
      let rangeStart = outputArray.count
      for match in matches[laneID] {
        let atomID = match[0].bitPattern
        outputArray.append(atomID)
      }
      let rangeEnd = outputArray.count
      outputRanges.append(rangeStart..<rangeEnd)
    }
  }
  
  for range in outputRanges {
    let slice = outputArray[range]
    outputSlices.append(slice)
  }
  return outputSlices
}

@inline(__always)
func createBoundingBox(
  vID: Int,
  size: Int,
  array: [SIMD4<Float>]
) -> SIMD8<Float> {
  var minimum = SIMD4<Float>(repeating: .greatestFiniteMagnitude)
  var maximum = -minimum
  for lane in 0..<size {
    let atom = array[vID &* size &+ lane]
    minimum.replace(with: atom, where: atom .< minimum)
    maximum.replace(with: atom, where: atom .> maximum)
  }
  var bounds = SIMD8(
    lowHalf: minimum + maximum,
    highHalf: maximum - minimum)
  bounds /= 2
  
  // Write the maximum deviation^2 from the center to bounds[3]
  var maxCenterDeviationSq: Float = .zero
  for lane in 0..<size {
    let atom = array[vID &* size &+ lane]
    var delta = atom - bounds.lowHalf
    delta.w = 0
    
    let deviationSq = (delta * delta).sum()
    maxCenterDeviationSq = max(maxCenterDeviationSq, deviationSq)
  }
  bounds[3] = maxCenterDeviationSq.squareRoot()
  
  // Write the maximum radius to bounds[7]
  bounds[7] = maximum.w
  return bounds
}

@inline(__always)
func blockBounds(
  _ transformed: [SIMD4<Float>],
  size: Int
) -> [SIMD8<Float>] {
  var blockBounds: [SIMD8<Float>] = []
  for vID in 0..<transformed.count/size {
    let bounds = createBoundingBox(vID: vID, size: size, array: transformed)
    blockBounds.append(bounds)
  }
  return blockBounds
}

@inline(__always)
func transform(
  _ atoms: [Entity],
  size: Int,
  algorithm: Topology.MatchAlgorithm
) -> [SIMD4<Float>] {
  var transformed: [SIMD4<Float>] = []
  for atomID in atoms.indices {
    let atom = atoms[atomID]
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
  }
  
  // Pad to the granularity of blocks.
  for _ in atoms.count..<(atoms.count + size - 1)/size*size {
    let original = transformed[atoms.count &- 1]
    transformed.append(original)
  }
  
  return transformed
}

@inline(__always)
func compareBlocks(
  _ lhsBlockBound: SIMD8<Float>,
  _ rhsBlockBound: SIMD8<Float>,
  executeCondition2: Bool = false
) -> Bool {
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
  
  let blockCenterX = lhsBlockBound.lowHalf
  let blockCenterY = rhsBlockBound.lowHalf
  let paddedCutoff = lhsBlockBound[7] + rhsBlockBound[7]
  let largeCutoff = paddedCutoff + blockCenterX.w + blockCenterY.w
  let largeCutoffSquared = largeCutoff * largeCutoff
  
  var blockDelta = blockCenterX - blockCenterY
  blockDelta.w = 0
  let condition1 = (blockDelta * blockDelta).sum() < largeCutoffSquared
  if !executeCondition2 {
    return condition1
  } else if condition1 {
    return true
  }
  
  let blockSizeX = lhsBlockBound.highHalf
  let blockSizeY = rhsBlockBound.highHalf
  let paddedCutoffSquared = paddedCutoff * paddedCutoff

  blockDelta.replace(with: -blockDelta, where: blockDelta .< 0)
  blockDelta = blockDelta - blockSizeX - blockSizeY
  blockDelta.replace(with: SIMD4.zero, where: blockDelta .< 0)
  blockDelta.w = 0
  let condition2 = (blockDelta * blockDelta).sum() < paddedCutoffSquared
  return condition2
}
