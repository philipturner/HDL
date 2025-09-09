//
//  Match.swift
//  HDL
//
//  Created by Philip Turner on 12/23/23.
//

import Dispatch

extension Topology {
  public enum MatchAlgorithm: Sendable {
    /// Search for neighbors within a fixed radius (in nanometers).
    case absoluteRadius(Float)
    
    /// Search for neighbors within a multiple of covalent bond length.
    case covalentBondLength(Float)
  }
  
  public typealias MatchStorage = ArraySlice<UInt32>
  
  public func match(
    _ source: [Atom],
    algorithm: MatchAlgorithm = .covalentBondLength(1.5),
    maximumNeighborCount: Int = 8
  ) -> [MatchStorage] {
    // During Topology.match, there are some situations where two similarly
    // sized grids will be constructed. They might be the exact same, although
    // the program can't detect that fact in a generalizable/robust manner.
    // Parallelization offers a simpler alternative that, based on the data
    // below, provides about the same speedup as eliding the compute work.
    
    // 'lattice' configuration, serial
    //
    // bounds | atoms  | octree |  0.25 |  0.5 |    1 |    2 |    4 | optimized
    // ------ | ------ | ------ | ----- | ---- | ---- | ---- | ---- | ----------
    // 5      |   2100 |    136 |   286 |  142 |  146 |      |      |  174
    // 7      |   5684 |    411 |   472 |  299 |  315 |  508 |      |  293
    // 10     |  16400 |   1168 |  2345 |  866 |  698 |  686 | 1276 |  887
    // 14     |  44688 |   3333 |  3447 | 2122 | 1863 | 1775 | 3512 | 1891
    // 20     | 129600 |   9245 | 19899 | 6695 | 5882 | 5332 | 4959 | 5403
    
    // 'lattice' configuration, 2x duplicated
    //
    // bounds | atoms  | octree | serial | parallel | speedup
    // ------ | ------ | ------ | ------ | -------- | ----------
    // 5      |   2100 |    314 |    298 |      186 | 1.1 -> 1.7
    // 7      |   5684 |    750 |    562 |      344 | 1.3 -> 2.2
    // 10     |  16400 |   2370 |   1555 |      905 | 1.5 -> 2.6
    // 14     |  44688 |   6085 |   3789 |     2160 | 1.6 -> 2.8
    // 20     | 129600 |  19932 |  10811 |     6567 | 1.8 -> 3.0
    
    let rmsAtomCount = (Float(source.count) * Float(atoms.count)).squareRoot()
    @Sendable
    func reorder(_ atoms: [Atom]) -> [UInt32] {
      if rmsAtomCount < 10_000 {
        return atoms.indices.map(UInt32.init)
      } else {
        let sorter = OctreeSorter(atoms: atoms)
        let state = sorter.traverseHighLevels()
        return sorter.traverseLowLevels(state: state)
      }
    }
    
    let safeAtoms = self.atoms
    nonisolated(unsafe)
    var lhsOperand: Operand = .init(atomCount: 0)
    nonisolated(unsafe)
    var rhsOperand: Operand = .init(atomCount: 0)
    DispatchQueue.concurrentPerform(iterations: 2) { z in
      if z == 0 {
        let lhsReordering = reorder(source)
        lhsOperand = transform8(
          source, lhsReordering, include4: false, algorithm: algorithm)
        lhsOperand.reordering = lhsReordering
      } else if z == 1 {
        let rhsReordering = reorder(safeAtoms)
        rhsOperand = transform8(
          safeAtoms, rhsReordering, include4: true, algorithm: algorithm)
        rhsOperand.reordering = rhsReordering
      }
    }
    
    // Call the actual matching function.
    let outputSlices = matchImpl(
      lhs: lhsOperand,
      rhs: rhsOperand,
      algorithm: algorithm,
      matchLimit: UInt32(maximumNeighborCount))
    
    return outputSlices
  }
}

// MARK: - Prepare Acceleration Structure

private struct Operand {
  var atomCount: Int
  var paddedCutoff: Float = .zero
  var transformed4: [SIMD4<Float>] = []
  var transformed8: [SIMD8<Float>] = []
  var blockBounds8: [SIMD4<Float>] = []
  var blockBounds32: [SIMD4<Float>] = []
  var blockBounds128: [SIMD4<Float>] = []
  var reordering: [UInt32] = []
}

private func transform8(
  _ atoms: [Atom],
  _ reordering: [UInt32],
  include4: Bool,
  algorithm: Topology.MatchAlgorithm
) -> Operand {
  var output = Operand(atomCount: atoms.count)
  let paddedAtomCount = (atoms.count + 127) / 128 * 128
  output.transformed8.reserveCapacity(4 * paddedAtomCount / 8)
  if include4 {
    output.transformed4.reserveCapacity(paddedAtomCount)
  }
  output.blockBounds8.reserveCapacity(paddedAtomCount / 8)
  output.blockBounds32.reserveCapacity(paddedAtomCount / 32)
  output.blockBounds128.reserveCapacity(paddedAtomCount / 128)
  
  // Accumulate data at the 32-atom granularity. At the same time, generate
  // most of the other per-operand data.
  var currentBlockBox32: SIMD8<Float> = .zero
  var blockBoxes32: [SIMD8<Float>] = []
  blockBoxes32.reserveCapacity(paddedAtomCount / 32)
  
  do {
    let taskCount = paddedAtomCount / 128
    for vID128 in 0..<taskCount {
      block1(UInt32(truncatingIfNeeded: vID128))
    }
    if paddedAtomCount > 0 {
      block2(UInt32(truncatingIfNeeded: paddedAtomCount / 128) &- 1)
    }
    for vID128 in 0..<taskCount {
      block3(UInt32(truncatingIfNeeded: vID128))
      block4(UInt32(truncatingIfNeeded: vID128))
    }
  }
  
  @inline(__always)
  func block1(_ vID128: UInt32) {
    let vIDStart = vID128 &* 16
    let vIDEnd = min(vIDStart &+ 16, (UInt32(atoms.count) + 7) / 8)
    for vID in vIDStart..<vIDEnd {
      var vecX: SIMD8<Float> = .zero
      var vecY: SIMD8<Float> = .zero
      var vecZ: SIMD8<Float> = .zero
      var vecR: SIMD8<Float> = .zero
      
      if vID &* 8 &+ 8 < UInt32(truncatingIfNeeded: atoms.count) {
        for laneID in 0..<UInt32(8) {
          let mortonIndex = reordering[Int(vID &* 8 &+ laneID)]
          let atom = atoms[Int(mortonIndex)]
          vecX[Int(laneID)] = atom.x
          vecY[Int(laneID)] = atom.y
          vecZ[Int(laneID)] = atom.z
          
          let atomicNumber = Int(atom.w)
          vecR[Int(laneID)] = Element.covalentRadii[atomicNumber]
          if include4 {
            output.transformed4.append(atom)
          }
        }
      } else {
        let validLanes = UInt32(truncatingIfNeeded: atoms.count) &- vID &* 8
        for laneID in 0..<validLanes {
          let mortonIndex = reordering[Int(vID &* 8 &+ laneID)]
          let atom = atoms[Int(mortonIndex)]
          vecX[Int(laneID)] = atom.x
          vecY[Int(laneID)] = atom.y
          vecZ[Int(laneID)] = atom.z
          
          let atomicNumber = Int(atom.w)
          vecR[Int(laneID)] = Element.covalentRadii[atomicNumber]
          if include4 {
            output.transformed4.append(atom)
          }
        }
        let last = (paddedAtomCount > 0 && include4)
        ? output.transformed4.last! : .zero
        for laneID in validLanes..<8 {
          vecX[Int(laneID)] = vecX[Int(validLanes &- 1)]
          vecY[Int(laneID)] = vecY[Int(validLanes &- 1)]
          vecZ[Int(laneID)] = vecZ[Int(validLanes &- 1)]
          vecR[Int(laneID)] = vecR[Int(validLanes &- 1)]
          if include4 {
            output.transformed4.append(last)
          }
        }
      }
      switch algorithm {
      case .absoluteRadius(let radius):
        vecR = SIMD8(repeating: radius / 2)
      case .covalentBondLength(let scale):
        vecR *= scale
      }
      if include4 {
        for laneID in 0..<UInt32(8) {
          output.transformed4[Int(vID &* 8 &+ laneID)].w = vecR[Int(laneID)]
        }
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
        let bounds =
        (currentBlockBox32.lowHalf + currentBlockBox32.highHalf) / 2
        output.blockBounds32.append(bounds)
        blockBoxes32.append(currentBlockBox32)
      }
    }
  }
  
  @inline(__always)
  func block2(_ vID128: UInt32) {
    let vID8End = (UInt32(atoms.count) + 7) / 8
    let vID32End = (UInt32(atoms.count) + 31) / 32
    let vID128End = (UInt32(atoms.count) + 127) / 128
    
    if include4 {
      let last = (paddedAtomCount > 0)
      ? output.transformed4.last! : .zero
      for _ in vID8End * 8..<vID128End * 128 {
        output.transformed4.append(last)
      }
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
    
    let box32 = (paddedAtomCount > 0)
    ? blockBoxes32.last! : .zero
    
    let bound32 = (paddedAtomCount > 0)
    ? output.blockBounds32.last! : .zero
    for _ in vID32End..<vID128End * 4 {
      blockBoxes32.append(box32)
      output.blockBounds32.append(bound32)
    }
  }
  
  // Generate the centers at 128-atom granularity in a standalone pass.
  // Then, the array for 'blockBoxes32' should be released.
  @inline(__always)
  func block3(_ vID128: UInt32) {
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
  @inline(__always)
  func block4(_ vID128: UInt32) {
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

// MARK: - Compare Blocks

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

// MARK: - Compare Atoms

private func matchImpl(
  lhs: Operand,
  rhs: Operand,
  algorithm: Topology.MatchAlgorithm,
  matchLimit: UInt32
) -> [Topology.MatchStorage] {
  let lhsSize32 = UInt32(truncatingIfNeeded: lhs.atomCount &+ 31) / 32
  let rhsSize32 = UInt32(truncatingIfNeeded: rhs.atomCount &+ 31) / 32
  let lhsSize8 = UInt32(truncatingIfNeeded: lhs.atomCount &+ 7) / 8
  let rhsSize8 = UInt32(truncatingIfNeeded: rhs.atomCount &+ 7) / 8
  let paddedCutoff = lhs.paddedCutoff + rhs.paddedCutoff
  
  // The current implementation doesn't allow a maximum neighbor count that
  // exceeds 254.
  guard matchLimit < UInt8.max else {
    fatalError("Maximum neighbor count is too large.")
  }
  
  // Allocate the match buffer.
  nonisolated(unsafe)
  let matchBuffer: UnsafeMutablePointer<SIMD2<Float>> =
    .allocate(capacity: Int(lhsSize8) * 8 * Int(matchLimit + 1))
  
  // Allocate the match count buffer.
  nonisolated(unsafe)
  let matchCount: UnsafeMutablePointer<SIMD8<UInt8>> =
    .allocate(capacity: Int(lhsSize8))
  matchCount.initialize(repeating: .zero, count: Int(lhsSize8))
  let opaqueCount = OpaquePointer(matchCount)
  nonisolated(unsafe)
  let castedCount = UnsafeMutablePointer<UInt8>(opaqueCount)
  
  // Allocate the output range buffer.
  let outRangeCapacity = Int(lhs.atomCount)
  var outRangeBuffer = [ArraySlice<UInt32>?](
    unsafeUninitializedCapacity: outRangeCapacity
  ) { $1 = outRangeCapacity }
  nonisolated(unsafe)
  let outRangePointer = outRangeBuffer.withUnsafeMutableBufferPointer { $0 }
  
  struct LoopScope {
    var startI: UInt32 = .zero
    var endI: UInt32 = .zero
    var startJ: UInt32 = .zero
    var endJ: UInt32 = .zero
  }
  
  @Sendable
  func createLoopScope3() -> LoopScope {
    var scope = LoopScope()
    scope.startI = 0
    scope.endI = scope.startI + UInt32(lhs.atomCount * 2) + 256
    scope.endI = min(scope.endI, UInt32(lhs.atomCount + 127) / 128)
    
    scope.startJ = 0
    scope.endJ = scope.startJ + UInt32(rhs.atomCount * 2) + 256
    scope.endJ = min(scope.endJ, UInt32(rhs.atomCount + 127) / 128)
    return scope
  }
  
  do {
    let scope = createLoopScope3()
    let taskCount = scope.endI - scope.startI
    
    if taskCount == 0 {
      
    } else if taskCount == 1 {
      innerLoop3(0)
      finishInnerLoop3(0)
    } else {
      DispatchQueue.concurrentPerform(iterations: Int(taskCount)) { z in
        innerLoop3(UInt32(z))
        finishInnerLoop3(UInt32(z))
      }
    }
  }
  
  @Sendable
  func innerLoop3(_ vIDi3: UInt32) {
    let scope = createLoopScope3()
    for vIDj3 in scope.startJ..<scope.endJ {
      let lhsBlockBound = lhs.blockBounds128[Int(vIDi3)]
      let rhsBlockBound = rhs.blockBounds128[Int(vIDj3)]
      let mask = compareBlocks(
        lhsBlockBound, rhsBlockBound, paddedCutoff: paddedCutoff)
      if mask {
        innerLoop2(vIDi3, vIDj3)
      }
    }
  }
  
  @Sendable
  func finishInnerLoop3(_ vIDi3: UInt32) {
    let scalarStart = vIDi3 &* 128
    let scalarEnd = min(
      scalarStart &+ 128, UInt32(truncatingIfNeeded: lhs.atomCount))
    
    let outMatchBuffer = [UInt32](
      unsafeUninitializedCapacity: 128 &* Int(matchLimit &+ 1)
    ) { buffer, bufferCount in
      bufferCount = 128 &* Int(matchLimit &+ 1)
      let baseAddress = buffer.baseAddress.unsafelyUnwrapped
      
      for laneID in scalarStart..<scalarEnd {
        let address = laneID &* (matchLimit &+ 1)
        let outAddress = (laneID &- scalarStart) &* (matchLimit &+ 1)
        let count = Int(castedCount[Int(laneID)])
        
        // Remove bogus matches that result from padding the RHS.
        var localBuffer = UnsafeMutableBufferPointer<SIMD2<Float>>(
          start: matchBuffer.advanced(by: Int(address)), count: count)
        while let last = localBuffer.last, last[0].bitPattern >= rhs.atomCount {
          localBuffer = UnsafeMutableBufferPointer(
            start: localBuffer.baseAddress, count: localBuffer.count - 1)
        }
        castedCount[Int(laneID)] = UInt8(truncatingIfNeeded: localBuffer.count)
        
        // Sort the matches in order of ascending distance.
        localBuffer.sort { $0.y < $1.y }
        
        let compactedOpaque = OpaquePointer(
          localBuffer.baseAddress.unsafelyUnwrapped)
        let compactedCasted = UnsafeMutablePointer<UInt32>(compactedOpaque)
        for i in 0..<localBuffer.count {
          let index = localBuffer[i][0].bitPattern
          let mortonIndex = rhs.reordering[Int(index)]
          compactedCasted[i] = mortonIndex
          baseAddress[Int(outAddress) &+ i] = mortonIndex
        }
      }
    }
    
    for laneID in scalarStart..<scalarEnd {
      let outAddress = (laneID &- scalarStart) &* (matchLimit &+ 1)
      let count = Int(castedCount[Int(laneID)])
      let outID = lhs.reordering[Int(laneID)]
      
      let opaqueRange = OpaquePointer(outRangePointer.baseAddress!)
      let castedRange = UnsafeMutablePointer<UInt>(opaqueRange)
      let pointer = castedRange.advanced(by: Int(outID &* 4))
      for i in 0..<4 {
        // Initialize the array to 'nil', on multiple CPU cores at once,
        // without causing a crash.
        pointer[i] = 0
      }
      
      let rangeStart = Int(outAddress)
      let rangeEnd = rangeStart &+ count
      outRangePointer[Int(outID)] = outMatchBuffer[rangeStart..<rangeEnd]
    }
  }
  
  @Sendable
  @inline(__always)
  func innerLoop2(_ vIDi3: UInt32, _ vIDj3: UInt32) {
    var scope = LoopScope()
    scope.startI = vIDi3 * 4
    scope.endI = scope.startI + 4
    scope.endI = min(scope.endI, lhsSize32)
    
    scope.startJ = vIDj3 * 4
    scope.endJ = scope.startJ + 4
    scope.endJ = min(scope.endJ, rhsSize32)
    
    for vIDi2 in scope.startI..<scope.endI {
      for vIDj2 in scope.startJ..<scope.endJ {
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
  
  @Sendable
  @inline(__always)
  func innerLoop1(_ vIDi2: UInt32, _ vIDj2: UInt32) {
    var scope = LoopScope()
    scope.startI = vIDi2 * 4
    scope.endI = scope.startI + 4
    scope.endI = min(scope.endI, lhsSize8)
    
    scope.startJ = vIDj2 * 4
    scope.endJ = scope.startJ + 4
    scope.endJ = min(scope.endJ, rhsSize8)
    
    for vIDi in scope.startI..<scope.endI {
      let pointer = Int(vIDi &* 4)
      let lhsX = lhs.transformed8[pointer &+ 0]
      let lhsY = lhs.transformed8[pointer &+ 1]
      let lhsZ = lhs.transformed8[pointer &+ 2]
      let lhsR = lhs.transformed8[pointer &+ 3]
      let lhsBlockBound = lhs.blockBounds8[Int(vIDi)]
      var lhsMatchCount = matchCount[Int(vIDi)]
      defer {
        matchCount[Int(vIDi)] = lhsMatchCount
      }
      
      for vIDj in scope.startJ..<scope.endJ {
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
          let mask32 = deltaSquared .<= r * r
          
          // There is no official API for converting between masks with
          // different bitwidths.
          let bypass32 = unsafeBitCast(mask32, to: SIMD8<UInt32>.self)
          let bypass8 = SIMD8<UInt8>(truncatingIfNeeded: bypass32)
          let mask8 = unsafeBitCast(
            bypass8, to: SIMDMask<SIMD8<UInt8>.MaskStorage>.self)
          
          if unsafeBitCast(mask8, to: UInt64.self) > 0 {
            for laneID in 0..<UInt32(8) {
              let count = lhsMatchCount[Int(laneID)]
              let address = vIDi &* 8 &+ laneID
              let address2 = address &* (matchLimit &+ 1) &+ UInt32(count)
              let keyValuePair = SIMD2(Float(bitPattern: atomID),
                                       deltaSquared[Int(laneID)])
              matchBuffer[Int(address2)] = keyValuePair
            }
            
            var nextCount = lhsMatchCount &+ 1
            let limit = SIMD8(repeating: UInt8(truncatingIfNeeded: matchLimit))
            nextCount.replace(with: limit, where: nextCount .> limit)
            lhsMatchCount.replace(with: nextCount, where: mask8)
          }
        }
      }
    }
  }
  
  matchBuffer.deallocate()
  matchCount.deallocate()
  
  return unsafeBitCast(outRangeBuffer, to: [_].self)
}
