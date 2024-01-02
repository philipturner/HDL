//
//  Orbitals.swift
//  
//
//  Created by Philip Turner on 12/24/23.
//

import Dispatch
import QuartzCore

extension Topology {
  public enum OrbitalHybridization {
    case sp1
    case sp2
    case sp3
    
    // Private API for generating nonbonding orbitals.
    var piBondCount: Int {
      switch self {
      case .sp1: return 2
      case .sp2: return 1
      case .sp3: return 0
      }
    }
  }
  
  // After reducing the overhead of array slice generation, elide the creation
  // of the atoms-to-atoms map. Just use the connections buffer directly.
  public func nonbondingOrbitals(
    hybridization: OrbitalHybridization = .sp3
  ) -> [ArraySlice<SIMD3<Float>>] {
    /*
     Some performance data gathered before performing any optimizations:
     -   atoms: 237704
     -     map: 8060
     -   alloc: 12
     - compute: 4000
     -  object: 729
     */
    
    let connectionsMap = createConnnectionsMap(secondaryType: .atoms)
    
   
    
    
    
    var outRangeBuffer = [ArraySlice<SIMD3<Float>>?](
      unsafeUninitializedCapacity: atoms.count
    ) { $1 = atoms.count }
    let outRangePointer = outRangeBuffer.withUnsafeMutableBufferPointer { $0 }
    
    
    
    var orbitalCounts: [Int] = Array(repeating: 0, count: atoms.count)
    let outputArray = [SIMD3<Float>](
      unsafeUninitializedCapacity: atoms.count &* 2
    ) { buffer, bufferCount in
      bufferCount = atoms.count &* 2
      let baseAddress = buffer.baseAddress.unsafelyUnwrapped
      
      let taskSize: Int = 5_000
      let taskCount = (atoms.count + taskSize - 1) / taskSize
      
      func execute(taskID: Int) {
        let scalarStart = taskID &* taskSize
        let scalarEnd = min(scalarStart &+ taskSize, atoms.count)
        for atomID in scalarStart..<scalarEnd {
          let (orbitalCount, orbital1, orbital2) = addOrbitals(
            atoms: atoms,
            atomID: atomID,
            bonds: bonds,
            connectionsMap: connectionsMap,
            hybridization: hybridization)
          orbitalCounts[atomID &- 0] = orbitalCount
          
          if orbitalCount > 0 {
            baseAddress[(atomID &- 0) &* 2 &+ 0] = orbital1
            baseAddress[(atomID &- 0) &* 2 &+ 1] = orbital2
          }
        }
      }
      
      if taskCount == 0 {
        
      } else if taskCount == 1 {
        for taskID in 0..<taskCount {
          execute(taskID: taskID)
        }
      } else {
        DispatchQueue.concurrentPerform(iterations: taskCount) { z in
          execute(taskID: z)
        }
      }
    }
    
    for atomID in 0..<atoms.count {
      let rangeStart = (atomID &- 0) &* 2
      let rangeEnd = rangeStart &+ orbitalCounts[atomID &- 0]
      
      let opaqueRange = OpaquePointer(outRangePointer.baseAddress!)
      let castedRange = UnsafeMutablePointer<UInt>(opaqueRange)
      let pointer = castedRange.advanced(by: Int(atomID &* 4))
      for i in 0..<4 {
        // Initialize the array to 'nil', on multiple CPU cores at once,
        // without causing a crash.
        pointer[i] = 0
      }
      if rangeEnd > rangeStart {
        outRangePointer[atomID] = outputArray[rangeStart..<rangeEnd]
      }
    }
    
    return unsafeBitCast(outRangeBuffer, to: [_].self)
  }
}

@inline(__always)
private func addOrbitals(
  atoms: [Entity],
  atomID: Int,
  bonds: [SIMD2<UInt32>],
  connectionsMap: [SIMD8<Int32>],
  hybridization: Topology.OrbitalHybridization
) -> (Int, SIMD3<Float>, SIMD3<Float>) {
  let atom = atoms[atomID]
  let atomicNumber = UInt8(atom.storage.w)
  let valence: Int
  switch Element(rawValue: atomicNumber) {
  case .hydrogen: valence = 1
  case .carbon: valence = 4
  case .nitrogen: valence = 3
  case .oxygen: valence = 2
  case .fluorine: valence = 1
  case .silicon: valence = 4
  case .phosphorus: valence = 3
  case .sulfur: valence = 2
  case .germanium: valence = 4
  case .gold: valence = 0
  case nil:
    fatalError("Invalid atomic number.")
  }
  let sigmaBondCount = valence - hybridization.piBondCount
  guard sigmaBondCount >= 1 else {
    return (0, .zero, .zero)
  }
  
  let neighborIDs = connectionsMap[atomID]
  var returnEarly = true
  switch (hybridization) {
  case .sp1:
    if neighborIDs[0] != -1, neighborIDs[1] == -1 {
      returnEarly = false
    }
  case .sp2:
    if neighborIDs[1] != -1, neighborIDs[2] == -1 {
      returnEarly = false
    }
  case .sp3:
    if neighborIDs[1] != -1, neighborIDs[2] == -1 {
      returnEarly = false
    }
    if neighborIDs[2] != -1, neighborIDs[3] == -1 {
      returnEarly = false
    }
  }
  if returnEarly {
    return (0, .zero, .zero)
  }
  
  var neighborCount = 0
  for lane in 0..<8 {
    if neighborIDs[lane] == -1 {
      break
    }
    neighborCount &+= 1
  }
  
  return withUnsafeTemporaryAllocation(of: SIMD3<Float>.self, capacity: 3) {
    let deltas = $0
    var distances = SIMD3<Float>(repeating: .greatestFiniteMagnitude)
    var normal: SIMD3<Float> = .zero
    
    // Calculate deltas between the atom and its neighbors.
    for i in 0..<neighborCount {
//      let bondID = neighborIDs[i]
//      let bond = bonds[Int(bondID)]
//      let neighborID = (bond[0] == atomID) ? bond[1] : bond[0]
      let neighborID = neighborIDs[i]
      
      let neighbor = atoms[Int(neighborID)]
      let delta4 = neighbor.storage - atom.storage
      let delta = unsafeBitCast(delta4, to: SIMD3<Float>.self)
      
      let distance = (delta * delta).sum().squareRoot()
      deltas[i] = delta / distance
      distances[i] = distance
      normal += delta / distance
    }
    if any(distances .< 0.001) {
      // Reject any bonds smaller than 1 picometer.
      return (0, .zero, .zero)
    }
    
    let normalLength = (normal * normal).sum().squareRoot()
    normal /= normalLength
    if normalLength < 0.001 {
      // Reject a normal smaller than 1 picometer.
      return (0, .zero, .zero)
    }
    
    // Branch on whether the situation resembles a sidewall carbon.
    if hybridization == .sp3 && neighborCount == 2 {
      var crossProduct: SIMD3<Float> = .zero
      let axis = deltas[1] - deltas[0]
      crossProduct.x = axis.y * normal.z - axis.z * normal.y
      crossProduct.y = axis.z * normal.x - axis.x * normal.z
      crossProduct.z = axis.x * normal.y - axis.y * normal.x
      
      let crossProductSquared = (crossProduct * crossProduct).sum()
      if crossProductSquared < 0.001 * 0.001 {
        // Reject a situation that outputs NAN direction vectors.
        return (0, .zero, .zero)
      }
      crossProduct /= crossProductSquared.squareRoot()
      
      // The order of the returned bonds is ambiguous, but it will be
      // deterministic after calling 'sort()'.
      let normalWeight = -Float(1.0 / 3).squareRoot()
      let crossProductWeight = Float(2.0 / 3).squareRoot()
      return (2,
      normal * normalWeight - crossProduct * crossProductWeight,
      normal * normalWeight + crossProduct * crossProductWeight)
    } else {
      // In the remaining cases, simply return something pointing opposite
      // to the average of the deltas.
      return (1, -normal, .zero)
    }
  }
}
