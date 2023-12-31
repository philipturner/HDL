//
//  Match.swift
//
//
//  Created by Philip Turner on 12/23/23.
//

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
    var rhsR = Element.covalentRadii[atomicNumber]
    
    switch algorithm {
    case .absoluteRadius(let radius):
      rhsR = radius / 2
    case .covalentBondLength(let scale):
      rhsR *= scale
    }
    var rhs = atom.storage
    rhs.w = rhsR
    rhsTransformed.append(rhs)
  }
  
  var outputArray: [UInt32] = []
  var outputRanges: [Range<Int>] = []
  var outputSlices: [ArraySlice<UInt32>] = []
  
  let vectorCount = (lhs.count + 7) / 8
  for vID in 0..<vectorCount {
    var lhsX: SIMD8<Float> = .zero
    var lhsY: SIMD8<Float> = .zero
    var lhsZ: SIMD8<Float> = .zero
    var lhsR: SIMD8<Float> = .zero
    let validLanes = min(8, lhs.count - vID * 8)
    precondition(validLanes > 0, "Unexpected lane count.")
    
    for laneID in 0..<validLanes {
      let atom = lhs[vID &* 8 &+ laneID]
      lhsX[laneID] = atom.storage.x
      lhsY[laneID] = atom.storage.y
      lhsZ[laneID] = atom.storage.z
      
      let atomicNumber = Int(atom.storage.w)
      let covalentRadius = Element.covalentRadii[atomicNumber]
      lhsR[laneID] = covalentRadius
    }
    for laneID in validLanes..<8 {
      lhsX[laneID] = lhsX[validLanes &- 1]
      lhsY[laneID] = lhsY[validLanes &- 1]
      lhsZ[laneID] = lhsZ[validLanes &- 1]
      lhsR[laneID] = lhsR[validLanes &- 1]
    }
    
    switch algorithm {
    case .absoluteRadius(let radius):
      lhsR = SIMD8(repeating: radius / 2)
    case .covalentBondLength(let scale):
      lhsR *= scale
    }
    
    var matches: [[SIMD2<Float>]] = []
    for _ in 0..<8 {
      var array: [SIMD2<Float>] = []
      array.reserveCapacity(8)
      matches.append(array)
    }
        
    for atomID in rhs.indices {
      let atom = rhsTransformed[atomID]
      let deltaX = lhsX - atom.x
      let deltaY = lhsY - atom.y
      let deltaZ = lhsZ - atom.z
      let deltaSquared = deltaX * deltaX + deltaY * deltaY + deltaZ * deltaZ
      
      let r = lhsR + atom.w
      let mask = deltaSquared .<= r * r
      if any(mask) {
        let atomID32 = UInt32(truncatingIfNeeded: atomID)
        
        // If this part is a bottleneck, there may be clever ways to improve the
        // parallelism of conditional writes to a compacted array. 8 arrays can
        // be generated concurrently, using zero control flow.
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
    
    for laneID in 0..<validLanes {
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
