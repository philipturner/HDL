//
//  Orbitals.swift
//  
//
//  Created by Philip Turner on 12/24/23.
//

import Dispatch

extension Topology {
  public enum OrbitalHybridization: Sendable {
    case sp
    case sp2
    case sp3
    
    // Private API for generating nonbonding orbitals.
    var piBondCount: Int {
      switch self {
      case .sp: return 2
      case .sp2: return 1
      case .sp3: return 0
      }
    }
  }
  
  public struct OrbitalStorage: Collection, Equatable {
    @usableFromInline var storage: SIMD16<Float>
    
    public typealias Index = Int
    
    public typealias Element = SIMD3<Float>
    
    @_transparent
    public var startIndex: Int { 0 }
    
    @_transparent
    public var endIndex: Int {
      Int(storage[15])
    }
    
    @_transparent
    public func index(after i: Int) -> Int {
      i &+ 1
    }
    
    public subscript(position: Int) -> SIMD3<Float> {
      @_transparent
      _read {
        var vec4: SIMD4<Float>
        switch position {
        case 0:
          vec4 = storage.lowHalf.lowHalf
        case 1:
          vec4 = storage.lowHalf.highHalf
        case 2:
          vec4 = storage.highHalf.lowHalf
        case 3:
          vec4 = storage.highHalf.highHalf
        default:
          fatalError("Invalid position: \(position)")
        }
        yield unsafeBitCast(vec4, to: SIMD3<Float>.self)
      }
    }
  }
  
  // After reducing the overhead of array slice generation, elide the creation
  // of the atoms-to-atoms map. Just use the connections buffer directly.
  public func nonbondingOrbitals(
    hybridization: OrbitalHybridization = .sp3
  ) -> [OrbitalStorage] {
    let connectionsMap = map(.atoms, to: .atoms)
    var storageBuffer = [OrbitalStorage](
      repeating: .init(storage: .zero), count: atoms.count)
    
    let taskSize: Int = 5_000
    let taskCount = (atoms.count + taskSize - 1) / taskSize
    
    // TODO: Fix the multiple errors that spawn when marking this function
    // as @Sendable.
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
        
        var storage = OrbitalStorage(storage: .zero)
        storage.storage.lowHalf.lowHalf = SIMD4(orbital1, 0)
        storage.storage.lowHalf.highHalf = SIMD4(orbital2, 0)
        storage.storage[15] = Float(orbitalCount)
        storageBuffer[atomID] = storage
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

    return storageBuffer
  }
}

@inline(__always)
private func addOrbitals(
  atoms: [Atom],
  atomID: Int,
  bonds: [SIMD2<UInt32>],
  connectionsMap: [Topology.MapStorage],
  hybridization: Topology.OrbitalHybridization
) -> (Int, SIMD3<Float>, SIMD3<Float>) {
  let atom = atoms[atomID]
  let atomicNumber = UInt8(atom.w)
  let valence: Int
  
  switch Element(rawValue: atomicNumber) {
  case .hydrogen: valence = 1
    
    // Always assume Group III elements are involved in a dative bond.
  case .boron: valence = 4
  case .carbon: valence = 4
  case .nitrogen: valence = 3
  case .oxygen: valence = 2
  case .fluorine: valence = 1
    
  case .aluminum: valence = 4
  case .silicon: valence = 4
  case .phosphorus: valence = 3
  case .sulfur: valence = 2
  case .chlorine: valence = 1
    
  case .gallium: valence = 4
  case .germanium: valence = 4
  case .arsenic: valence = 3
  case .selenium: valence = 2
  case .bromine: valence = 1
    
    // Always assume Group IV elements are tetravalent.
  case .tin: valence = 4
  case .gold: valence = 0
  case .lead: valence = 4
  case nil:
    fatalError("Invalid atomic number.")
  }
  let sigmaBondCount = valence - hybridization.piBondCount
  guard sigmaBondCount >= 1 else {
    return (0, .zero, .zero)
  }
  
  let neighborIDs = connectionsMap[atomID].storage
  switch (hybridization, neighborIDs[7]) {
  case (.sp, 1): break
  case (.sp2, 2): break
  case (.sp3, 2): break
  case (.sp3, 3): break
  default:
    return (0, .zero, .zero)
  }
  
  return withUnsafeTemporaryAllocation(of: SIMD3<Float>.self, capacity: 3) {
    let deltas = $0
    var distances = SIMD3<Float>(repeating: .greatestFiniteMagnitude)
    var normal: SIMD3<Float> = .zero
    
    // Calculate deltas between the atom and its neighbors.
    for i in 0..<Int(neighborIDs[7]) {
      let neighborID = neighborIDs[i]
      
      let neighbor = atoms[Int(neighborID)]
      let delta4 = neighbor - atom
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
    if hybridization == .sp3 && neighborIDs[7] == 2 {
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
