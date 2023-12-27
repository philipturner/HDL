//
//  Orbitals.swift
//  
//
//  Created by Philip Turner on 12/24/23.
//

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
  
  public func nonbondingOrbitals(
    hybridization: OrbitalHybridization = .sp3
  ) -> [ArraySlice<SIMD3<Float>>] {
    let atomsToAtomsMap = map(.atoms, to: .atoms)
    var outputArray: [SIMD3<Float>] = []
    var outputRanges: [Range<UInt32>] = []
    var outputSlices: [ArraySlice<SIMD3<Float>>] = []
    
  outer:
    for atomID in atoms.indices {
      let rangeStart = UInt32(truncatingIfNeeded: outputArray.count)
      defer {
        let rangeEnd = UInt32(truncatingIfNeeded: outputArray.count)
        outputRanges.append(rangeStart..<rangeEnd)
      }
      
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
      case .none:
        fatalError("Invalid atomic number.")
      }
      let sigmaBondCount = valence - hybridization.piBondCount
      guard sigmaBondCount >= 1 else {
        continue outer
      }
      
      let neighborIDs = atomsToAtomsMap[atomID]
      switch (hybridization, neighborIDs.count) {
      case (.sp1, 1): break
      case (.sp2, 2): break
      case (.sp3, 2): break
      case (.sp3, 3): break
      default:
        continue outer
      }
      
      withUnsafeTemporaryAllocation(of: SIMD3<Float>.self, capacity: 3) {
        let deltas = $0
        var distances = SIMD3<Float>(repeating: .greatestFiniteMagnitude)
        var normal: SIMD3<Float> = .zero
        
        // Calculate deltas between the atom and its neighbors.
        for i in 0..<neighborIDs.count {
          let neighborID = neighborIDs[neighborIDs.startIndex &+ i]
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
          return
        }
        
        let normalLength = (normal * normal).sum().squareRoot()
        normal /= normalLength
        if normalLength < 0.001 {
          // Reject a normal smaller than 1 picometer.
          return
        }
        
        // Branch on whether the situation resembles a sidewall carbon.
        if hybridization == .sp3 && neighborIDs.count == 2 {
          var crossProduct: SIMD3<Float> = .zero
          let axis = deltas[1] - deltas[0]
          crossProduct.x = axis.y * normal.z - axis.z * normal.y
          crossProduct.y = axis.z * normal.x - axis.x * normal.z
          crossProduct.z = axis.x * normal.y - axis.y * normal.x
          
          let crossProductSquared = (crossProduct * crossProduct).sum()
          if crossProductSquared < 0.001 * 0.001 {
            // Reject a situation that outputs NAN direction vectors.
            return
          }
          crossProduct /= crossProductSquared.squareRoot()
          
          // The order of the returned bonds is ambiguous, but it will be
          // deterministic after calling 'sort()'.
          let normalWeight = -Float(1.0 / 3).squareRoot()
          let crossProductWeight = Float(2.0 / 3).squareRoot()
          outputArray.append(
            normal * normalWeight - crossProduct * crossProductWeight)
          outputArray.append(
            normal * normalWeight + crossProduct * crossProductWeight)
        } else {
          // In the remaining cases, simply return something pointing opposite
          // to the average of the deltas.
          outputArray.append(-normal)
        }
      }
    }
    
    for range in outputRanges {
      let rangeStart = Int(truncatingIfNeeded: range.lowerBound)
      let rangeEnd = Int(truncatingIfNeeded: range.upperBound)
      let slice = outputArray[rangeStart..<rangeEnd]
      outputSlices.append(slice)
    }
    return outputSlices
  }
}
