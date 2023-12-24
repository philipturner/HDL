//
//  Match.swift
//
//
//  Created by Philip Turner on 12/23/23.
//

#if arch(arm64)
typealias Half = Float16
#else
typealias Half = Float32
#endif

// Notes for when you get around to optimizing this:
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

extension Topology {
  public enum MatchAlgorithm {
    // Search for neighbors within a fixed radius (in nanometers).
    case absoluteRadius(Float)
    
    // Search for neighbors within a multiple of covalent bond length.
    case covalentBondScale(Float)
  }
  
  public func match(
    _ input: [Entity],
    algorithm: MatchAlgorithm = .covalentBondScale(1.5)
  ) -> [ArraySlice<UInt32>] {
    // TODO: Activate the more advanced implementation below, once the simpler
    // one without reordering is ironed out. There may be no need to activate
    // it, until you get around to optimization.
    #if true
    return matchImpl(lhs: input, rhs: atoms, algorithm: algorithm)
    
    #else
    let lhsGrid = TopologyGrid(atoms: input)
    let lhsReordering = lhsGrid.mortonReordering()
    let rhsGrid = TopologyGrid(atoms: atoms)
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
    case .covalentBondScale(let scale):
      rhsR *= Half(scale).nextUp
    }
    var rhs = atom.storage
    rhs.w = Float(rhsR)
    rhsTransformed.append(rhs)
  }
  
  var outputArray: [UInt32] = []
  var outputRanges: [Range<Int>] = []
  var outputSlices: [ArraySlice<UInt32>] = []
  
  for vID in 0..<(lhs.count + 7) / 8 {
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
    case .covalentBondScale(let scale):
      lhsR *= Half(scale).nextUp
    }
    
    // A future optimization may allocate this beforehand. However, doing so
    // complicates the code.
    var distancesSquared: [SIMD8<Half>] = []
    var matchMasks: [UInt8] = []
    distancesSquared.reserveCapacity((lhs.count + 7) / 8)
    matchMasks.reserveCapacity((lhs.count + 7) / 8)
    
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
      distancesSquared.append(deltaSquared16)
      
      let r = lhsR + Half(atom.w)
      let rSquared = r * r
      let flags = SIMD8<Half.RawSignificand>(
        1 << 0, 1 << 1, 1 << 2, 1 << 3,
        1 << 4, 1 << 5, 1 << 6, 1 << 7)
      
      var output: SIMD8<Half.RawSignificand> = .zero
      output.replace(with: flags, where: deltaSquared16 .<= rSquared)
      let mask16 = output.wrappedSum()
      let mask8 = UInt8(truncatingIfNeeded: mask16)
      matchMasks.append(mask8)
    }
    
    distancesSquared.withUnsafeBufferPointer {
      let opaque = OpaquePointer($0.baseAddress)
      let casted = UnsafePointer<Half>(opaque).unsafelyUnwrapped
      var sortedMatches: [SIMD2<Float>] = []
      
      for laneID in 0..<validLanes {
        sortedMatches.removeAll(keepingCapacity: true)
        
        // Eventually, we will loop over blocks that the previous hierarchy
        // level marked as candidates. The granularity is per-block instead of
        // per-atom, so we don't compact the mask list in the O(n^2)
        // implementation.
        //
        // If this part is a bottleneck, there may be clever ways to improve the
        // parallelism of conditional writes to a compacted array. 8 arrays can
        // be simultaneously generated in parallel, using zero control flow.
        let lhsMask = UInt8(1) << laneID
        for atomID in rhs.indices {
          let mask = matchMasks[atomID] & lhsMask
          if mask > 0 {
            let atomID32 = UInt32(truncatingIfNeeded: atomID)
            let address = atomID &* 8 &+ laneID
            let keyValuePair = SIMD2(
              Float(bitPattern: atomID32),
              Float(casted[address]))
            sortedMatches.append(keyValuePair)
          }
        }
        
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
    }
  }
  
  for range in outputRanges {
    let slice = outputArray[range]
    outputSlices.append(slice)
  }
  return outputSlices
}
