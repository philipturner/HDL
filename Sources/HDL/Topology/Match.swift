//
//  Match.swift
//
//
//  Created by Philip Turner on 12/23/23.
//

import QuartzCore

#if arch(arm64)
typealias Half = Float16
#else
typealias Half = Float32
#endif

// Goals for optimizing this:
//
// MARK: - Near-Term Goal
//
// First step: a fully O(n^2) function that writes outputs into a format you
// want. After seeing this through, I can understand how to implement the upper
// levels of the hierarchy, and transfer data between them.
//
// Get the test working with this naive algorithm. Then, you can make informed
// optimizations based on experimental feedback, minimizing the time wasted on
// premature optimizations. Fix every other bottleneck before tackling O(n^2)
// scaling. Prove that the inner loop of the atom-by-atom comparison maximizes
// vALU% and avoids register spills (maybe even theorize about maximum possible
// performance using GFLOPS).
//
// Also, think about the overhead of the scalarized sorting code. Can you reduce
// the time spent in quicksort by exploiting properties of Morton order? The
// closest points might be in the middle of the list, while the farthest points
// are in the end. Perhaps a fixed halfway split or a poll to locate the
// smallest point, or a 50-50 mix between the two heuristics.
//
// MARK: - Far-Term Goal
//
// Higher levels of the hierarchy: recursive bounding box-bounding box test,
// also with search radius in the same manner as OpenMM PADDED_CUTOFF. Store
// block position as 3 x SIMD8<Float>, block bounds as 3 x SIMD8<Half>, block
// radii as SIMD8<Half>.
// - use FloatingPoint.nextUp to round toward positive infinity
// - make the two halves of each SIMD8 independent, so they can be paged into
//   ISA registers without spilling
//
// Lowest level of the hierarchy: 8 atoms in one block against 8 atoms in
// another block. Store atom data as SIMD4<Float> in memory, but unpack into
// 3 x SIMD8<Float> + SIMD8<Half> in registers.

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
    #if true
    return matchImpl(lhs: input, rhs: atoms, algorithm: algorithm)
    
    #else
    let lhsGrid = GridSorter(atoms: input)
    let lhsReordering = lhsGrid.mortonReordering()
    let rhsGrid = GridSorter(atoms: atoms)
    let rhsReordering = rhsGrid.mortonReordering()
    
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
  lhs: [Entity],
  rhs: [Entity],
  algorithm: Topology.MatchAlgorithm
) -> [ArraySlice<UInt32>] {
  var rhsTransformed: [SIMD4<Float>] = []
  for atomID in rhs.indices {
    let atom = rhs[atomID]
    let atomicNumber = Int(atom.storage.w)
    let covalentRadius = Element.covalentRadii[atomicNumber]
    var rhsR = Half(covalentRadius).nextUp
    
    switch algorithm {
    case .absoluteRadius(let radius):
      rhsR = Half(radius / 2).nextUp
    case .covalentBondLength(let scale):
      rhsR *= Half(scale).nextUp
    }
    var rhs = atom.storage
    rhs.w = Float(rhsR)
    rhsTransformed.append(rhs)
  }
  
  var outputArray: [UInt32] = []
  var outputRanges: [Range<Int>] = []
  var outputSlices: [ArraySlice<UInt32>] = []
  
  let vectorCount = (lhs.count + 7) / 8
#if PROFILE_MATCH
  var profilingOutput = "\(lhs.count) x \(rhs.count)\n"
  let divisor = vectorCount / 8
#endif
  
  for vID in 0..<vectorCount {
#if PROFILE_MATCH
    var checkpoint0: Double = .zero
    var checkpoint1: Double = .zero
    var checkpoint2: Double = .zero
    var checkpoint3: Double = .zero
    let sample = (vID % divisor == 0)
    
    if sample {
      checkpoint0 = CACurrentMediaTime()
    }
#endif
    
    var lhsX: SIMD8<Float> = .zero
    var lhsY: SIMD8<Float> = .zero
    var lhsZ: SIMD8<Float> = .zero
    var lhsR: SIMD8<Half> = .zero
    let validLanes = min(8, lhs.count - vID * 8)
    precondition(validLanes > 0, "Unexpected lane count.")
    
    for laneID in 0..<validLanes {
      let atom = lhs[vID &* 8 &+ laneID]
      lhsX[laneID] = atom.storage.x
      lhsY[laneID] = atom.storage.y
      lhsZ[laneID] = atom.storage.z
      
      let atomicNumber = Int(atom.storage.w)
      let covalentRadius = Element.covalentRadii[atomicNumber]
      lhsR[laneID] = Half(covalentRadius).nextUp
    }
    for laneID in validLanes..<8 {
      lhsX[laneID] = lhsX[validLanes &- 1]
      lhsY[laneID] = lhsY[validLanes &- 1]
      lhsZ[laneID] = lhsZ[validLanes &- 1]
      lhsR[laneID] = lhsR[validLanes &- 1]
    }
    
    switch algorithm {
    case .absoluteRadius(let radius):
      let radius16 = Half(radius / 2).nextUp
      lhsR = SIMD8(repeating: radius16)
    case .covalentBondLength(let scale):
      lhsR *= Half(scale).nextUp
    }
    
    // A future optimization may allocate this beforehand. However, doing so
    // complicates the code.
    var matches: [[SIMD2<Float>]] = []
    for _ in 0..<8 {
      matches.append([])
    }
    
#if PROFILE_MATCH
    if sample {
      checkpoint1 = CACurrentMediaTime()
    }
#endif
    
    for atomID in rhs.indices {
      // Maybe the loop will be faster if deltaX/deltaY/deltaZ are converted to
      // half-precision. Explore this possibility after getting the code to work
      // at all.
      let atom = rhsTransformed[atomID]
      let deltaX = lhsX - atom.x
      let deltaY = lhsY - atom.y
      let deltaZ = lhsZ - atom.z
      let deltaSquared32 = deltaX * deltaX + deltaY * deltaY + deltaZ * deltaZ
      let deltaSquared16 = SIMD8<Half>(deltaSquared32)
      
      let r = lhsR + Half(atom.w)
      let rSquared = r * r
      
      let mask = deltaSquared16 .<= rSquared
      if any(mask) {
        let atomID32 = UInt32(truncatingIfNeeded: atomID)
        let deltaSquared32 = SIMD8<Float>(deltaSquared16)
        
        // If this part is a bottleneck, there may be clever ways to improve the
        // parallelism of conditional writes to a compacted array. 8 arrays can
        // be simultaneously generated in parallel, using zero control flow.
        for laneID in 0..<8 {
          let keyValuePair = SIMD2(
            Float(bitPattern: atomID32),
            Float(deltaSquared32[laneID]))
          if mask[laneID] {
            matches[laneID].append(keyValuePair)
          }
        }
      }
    }
    
#if PROFILE_MATCH
    if sample {
      checkpoint2 = CACurrentMediaTime()
    }
#endif
    
    for laneID in 0..<validLanes {
      var sortedMatches = matches[laneID]
      
      // Sort the matches in ascending order of distance.
      sortedMatches.sort { $0.y < $1.y }
      
      let rangeStart = outputArray.count
      for match in sortedMatches {
        let atomID = match[0].bitPattern
        outputArray.append(atomID)
      }
      let rangeEnd = outputArray.count
      outputRanges.append(rangeStart..<rangeEnd)
    }
    
#if PROFILE_MATCH
    if sample {
      checkpoint3 = CACurrentMediaTime()
      
      let elapsedTimes = SIMD3<Double>(
        checkpoint1 - checkpoint0,
        checkpoint2 - checkpoint1,
        checkpoint3 - checkpoint2)
      let totalTime = elapsedTimes.sum()
      let proportions = elapsedTimes / totalTime
      let percents = SIMD3<Int>(proportions * 100)
      
      let output = """
      \(Int(totalTime * 1e6)) - \
      \(percents[0])% \(percents[1])% \(percents[2])%
      """
      profilingOutput += "vID=\(vID) time=\(output)\n"
    }
#endif
  }
  
  #if PROFILE_MATCH
  print(profilingOutput)
  #endif
  
  for range in outputRanges {
    let slice = outputArray[range]
    outputSlices.append(slice)
  }
  return outputSlices
}
