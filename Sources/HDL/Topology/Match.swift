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
  let (lhsTransformed, paddedCutoffLHS) =
  transform(lhsAtoms, size: size3, algorithm: algorithm)
  let (rhsTransformed, paddedCutoffRHS) =
  transform(rhsAtoms, size: size3, algorithm: algorithm)
  let paddedCutoff = paddedCutoffLHS + paddedCutoffRHS
  
  let lhsBlockBounds = blockBounds(lhsTransformed, size: 8)
  let rhsBlockBounds = blockBounds(rhsTransformed, size: 8)
  let lhs = (transformed: lhsTransformed, blockBounds: lhsBlockBounds)
  let rhs = (transformed: rhsTransformed, blockBounds: rhsBlockBounds)
  
  var outputMatches: [[SIMD2<Float>]] = []
  for _ in 0..<UInt32(lhsAtoms.count + 7) / 8 * 8 {
    var array: [SIMD2<Float>] = []
    array.reserveCapacity(8)
    outputMatches.append(array)
  }
  var outputArray: [UInt32] = []
  var outputRanges: [Range<Int>] = []
  var outputSlices: [ArraySlice<UInt32>] = []
  
  // MARK: - Hierarchy
  
  let lhs32 = blockBounds(lhsTransformed, size: size2)
  let rhs32 = blockBounds(rhsTransformed, size: size2)
  let lhs128 = blockBounds(lhsTransformed, size: size3)
  let rhs128 = blockBounds(rhsTransformed, size: size3)
  
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
      var lhsX: SIMD8<Float> = .zero
      var lhsY: SIMD8<Float> = .zero
      var lhsZ: SIMD8<Float> = .zero
      var lhsR: SIMD8<Float> = .zero
      for laneID in 0..<UInt32(8) {
        let atom = lhs.transformed[Int(vIDi &* 8 &+ laneID)]
        lhsX[Int(laneID)] = atom.x
        lhsY[Int(laneID)] = atom.y
        lhsZ[Int(laneID)] = atom.z
        lhsR[Int(laneID)] = atom.w
      }
      let lhsBlockBound = lhs.blockBounds[Int(vIDi)]
      
      for vIDj in loopStartJ..<loopEndJ {
        let rhsBlockBound = rhs.blockBounds[Int(vIDj)]
        guard compareBlocks(
          lhsBlockBound, rhsBlockBound, paddedCutoff: paddedCutoff) else {
          continue
        }
        
        for lane in 0..<UInt32(8) {
          let atomID = vIDj &* 8 &+ lane
          let atom = rhs.transformed[Int(atomID)]
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

@inline(__always)
func createBoundingBox(
  vID: UInt32,
  size: UInt32,
  array: [SIMD4<Float>]
) -> SIMD4<Float> {
  var minimum = SIMD4<Float>(repeating: .greatestFiniteMagnitude)
  var maximum = -minimum
  for lane in 0..<size {
    let atom = array[Int(vID &* size &+ lane)]
    minimum.replace(with: atom, where: atom .< minimum)
    maximum.replace(with: atom, where: atom .> maximum)
  }
  var bounds = (minimum + maximum) / 2
  
  // Write the maximum deviation^2 from the center to bounds[3]
  var maxCenterDeviationSq: Float = .zero
  for lane in 0..<size {
    let atom = array[Int(vID &* size &+ lane)]
    var delta = atom - bounds
    delta.w = 0
    
    let deviationSq = (delta * delta).sum()
    maxCenterDeviationSq = max(maxCenterDeviationSq, deviationSq)
  }
  bounds[3] = maxCenterDeviationSq.squareRoot()
  
  return bounds
}

@inline(__always)
func blockBounds(
  _ transformed: [SIMD4<Float>],
  size: UInt32
) -> [SIMD4<Float>] {
  var blockBounds: [SIMD4<Float>] = []
  for vID in 0..<UInt32(transformed.count) / size {
    let bounds = createBoundingBox(vID: vID, size: size, array: transformed)
    blockBounds.append(bounds)
  }
  return blockBounds
}

@inline(__always)
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

