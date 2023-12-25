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
        continue
      }
      
      let neighbors = atomsToAtomsMap[atomID]
      switch (hybridization, neighbors.count) {
      case (.sp1, 1):
        break
      case (.sp2, 2):
        break
      case (.sp3, 2):
        break
      case (.sp3, 3):
        break
      default:
        continue
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
