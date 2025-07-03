//
//  TwoWayCollisions.swift
//
//
//  Created by Philip Turner on 6/11/24.
//

extension Compilation {
  mutating func resolveTwoWayCollisions(centerTypes: [UInt8]) {
    // Extract all sets of connected hydrogen sites from the topology.
    var hydrogenChains: [[UInt32]] = []
    for hydrogenID in hydrogensToAtomsMap.indices {
      let atomList = hydrogensToAtomsMap[hydrogenID]
      guard atomList.count == 2 else {
        continue
      }
      
      // First elements: (C, H, C)
      let dimerGeometry = DimerGeometry(
        centerTypes: centerTypes,
        atomList: atomList)
      let initialDimerChain = startDimerChain(
        hydrogenID: UInt32(hydrogenID),
        dimerGeometry: dimerGeometry)
      guard let initialDimerChain else {
        continue
      }
      
      // Append several groups of (..., H, C)
      let dimerChain = growDimerChain(
        initialDimerChain: initialDimerChain)
      if dimerChain.count > 3 {
        print("SIMD4<UInt32>(\(dimerChain[0]), \(dimerChain[2]), \(dimerChain[dimerChain.count - 3]), \(dimerChain[dimerChain.count - 1])),")
      }
      
      // Even indices: carbon sites
      // Odd indices: hydrogen sites
      var hydrogenChain: [UInt32] = []
      for i in dimerChain.indices where i % 2 == 1 {
        let hydrogenID = dimerChain[i]
        hydrogenChain.append(hydrogenID)
      }
      hydrogenChains.append(hydrogenChain)
    }
    
    // With these chains available, decide which collision sites will be
    // replaced by a carbon-carbon bond.
    var collisionStates = [CollisionState](
      repeating: .noCollision,
      count: hydrogensToAtomsMap.count)
    for hydrogenChain in hydrogenChains {
      for i in hydrogenChain.indices {
        let hydrogenID = Int(hydrogenChain[i])
        guard collisionStates[hydrogenID] == .noCollision else {
          continue
        }
        
        if i % 2 == 0 {
          collisionStates[hydrogenID] = .mergeDimer
        } else {
          collisionStates[hydrogenID] = .keepCollision
        }
      }
    }
    mergeDimers(states: collisionStates)
  }
}

// MARK: - Utilities

private struct DimerGeometry {
  var bridgeheadID: UInt32?
  var sidewallID: UInt32?
  var bothBridgehead: Bool = true
  var bothSidewall: Bool = true
  
  init(
    centerTypes: [UInt8],
    atomList: [UInt32]
  ) {
    for atomID in atomList {
      switch centerTypes[Int(atomID)] {
      case 2:
        sidewallID = atomID
        bothBridgehead = false
      case 3:
        bridgeheadID = atomID
        bothSidewall = false
      default:
        fatalError("This should never happen.")
      }
    }
  }
}

private enum CollisionState {
  case noCollision
  case keepCollision
  case mergeDimer
}

extension Compilation {
  // The unsorted list must already be guaranteed to have 2 elements.
  private static func opposite(
    original: UInt32,
    unsortedList: [UInt32]
  ) -> UInt32 {
    if original == unsortedList[0] {
      return unsortedList[1]
    } else if original == unsortedList[1] {
      return unsortedList[0]
    } else {
      fatalError("Unexpected hydrogen list.")
    }
  }
  
  private func startDimerChain(
    hydrogenID: UInt32,
    dimerGeometry: DimerGeometry
  ) -> SIMD3<UInt32>? {
    // This list is guaranteed to have 2 elements.
    let atomList = hydrogensToAtomsMap[Int(hydrogenID)]
    
    if dimerGeometry.bothBridgehead {
      let hydrogenList = atomsToHydrogensMap[Int(atomList[0])]
      guard hydrogenList.count == 1 else {
        fatalError("Bridgehead site did not have 1 hydrogen.")
      }
      guard hydrogenID == hydrogenList[0] else {
        fatalError("Unexpected hydrogen list.")
      }
      
      // (C, H, C)
      return SIMD3(
        atomList[0],
        hydrogenID,
        atomList[1])
    } else if dimerGeometry.bothSidewall {
      // There is a guaranteed order for the dimers in a chain. Search for a
      // dimer at the start of this sequence, not the end.
      for atomID in atomList {
        let hydrogenList = atomsToHydrogensMap[Int(atomID)]
        guard hydrogenList.count == 2 else {
          fatalError("Sidewall site did not have 2 hydrogens.")
        }
        
        let oppositeHydrogenID = Self.opposite(
          original: hydrogenID,
          unsortedList: hydrogenList)
        let oppositeAtomList = hydrogensToAtomsMap[Int(oppositeHydrogenID)]
        
        // Cases for 'oppositeAtomList.count':
        // 0: guaranteed to be impossible
        // 1: (oppositeHydrogenID, atomID, hydrogenID) is a terminator
        //    (H, C, H, ...)
        // 2: (something else, oppositeHydrogenID, atomID, hydrogenID)
        //    (..., C, H, C, H, ...)
        // 3: supposedly impossible
        guard oppositeAtomList.count == 1 else {
          guard oppositeAtomList.count == 2 else {
            fatalError("This should never happen.")
          }
          continue
        }
        
        // (oppositeHydrogenID, atomID, hydrogenID, oppositeAtomID)
        // (H, C, H, C, ...)
        let oppositeAtomID = Self.opposite(
          original: atomID,
          unsortedList: atomList)
        return SIMD3(
          atomID,
          hydrogenID,
          oppositeAtomID)
      }
      
      // We found a dimer in the middle of a chain.
      //
      // If there is a long list of dimers, we don't try to travel the entire
      // row length when encountering each bond in the middle. That would be
      // an O(n^2) scaling algorithm. Instead, we only travel when the starting
      // point is one of the two ends.
      //
      // If the chain wraps back onto itself, forming a ring, the entire thing
      // fails to reconstruct.
      return nil
    } else if let bridgeheadID = dimerGeometry.bridgeheadID,
              let sidewallID = dimerGeometry.sidewallID {
      // (C, H, C, ...)
      return SIMD3(
        bridgeheadID,
        hydrogenID,
        sidewallID)
    } else {
      fatalError("This should never happen.")
    }
  }
  
  private func nextDimer(
    hydrogenID: UInt32,
    atomID: UInt32
  ) -> SIMD2<UInt32>? {
    let hydrogenList = atomsToHydrogensMap[Int(atomID)]
    switch hydrogenList.count {
    case 1:
      guard hydrogenID == hydrogenList[0] else {
        fatalError("Unexpected hydrogen list.")
      }
      
      // (..., H, C)
      return nil
    case 2:
      let oppositeHydrogenID = Self.opposite(
        original: hydrogenID,
        unsortedList: hydrogenList)
      let oppositeAtomList = hydrogensToAtomsMap[Int(oppositeHydrogenID)]
      
      // Cases for 'oppositeAtomList.count':
      // 0: guaranteed to be impossible
      // 1: (hydrogenID, atomID, oppositeHydrogenID) is a terminator
      //    (..., H, C, H)
      // 2: (hydrogenID, atomID, oppositeHydrogenID, something else)
      //    (..., H, C, H, C, ...)
      // 3: supposedly impossible
      guard oppositeAtomList.count == 2 else {
        guard oppositeAtomList.count == 1 else {
          fatalError("This should never happen.")
        }
        return nil
      }
      
      // (hydrogenID, atomID, oppositeHydrogenID, expandingAtomID)
      // (..., H, C, H, C, ...)
      let expandingAtomID = Self.opposite(
        original: atomID,
        unsortedList: oppositeAtomList)
      return SIMD2(
        oppositeHydrogenID,
        expandingAtomID)
    default:
      fatalError("This should never happen.")
    }
  }
  
  private func growDimerChain(
    initialDimerChain: SIMD3<UInt32>
  ) -> [UInt32] {
    var dimerChain = [
      initialDimerChain[0],
      initialDimerChain[1],
      initialDimerChain[2],
    ]
    
    // Iteratively search through the topology, seeing whether the chain
    // of linked carbon atoms finally ends.
    var converged = false
    for _ in 0..<4096 {
      let dimer = nextDimer(
        hydrogenID: dimerChain[dimerChain.count - 2],
        atomID: dimerChain[dimerChain.count - 1])
      if let dimer {
        dimerChain.append(dimer[0])
        dimerChain.append(dimer[1])
      } else {
        converged = true
        break
      }
    }
    guard converged else {
      fatalError("Took too many iterations to find length of dimer chain.")
    }
    
    return dimerChain
  }
  
  private mutating func mergeDimers(states: [CollisionState]) {
    var insertedBonds: [SIMD2<UInt32>] = []
    for hydrogenID in states.indices {
      let state = states[hydrogenID]
      guard state == .mergeDimer else {
        continue
      }
      
      // Fetch the references to the two carbons.
      let atomList = hydrogensToAtomsMap[hydrogenID]
      guard atomList.count == 2 else {
        fatalError("This should never happen.")
      }
      
      // Iterate over the two carbons.
      for atomID in atomList {
        let hydrogenList = atomsToHydrogensMap[Int(atomID)]
        
        func shrunkHydrogenList() -> [UInt32] {
          switch hydrogenList.count {
          case 1:
            return []
          case 2:
            let oppositeHydrogenID = Self.opposite(
              original: UInt32(hydrogenID),
              unsortedList: hydrogenList)
            return [oppositeHydrogenID]
          default:
            fatalError("This should never happen.")
          }
        }
        atomsToHydrogensMap[Int(atomID)] = shrunkHydrogenList()
      }
      
      // Erase this array slot because the conflict was resolved. The hydrogen
      // site no longer exists.
      hydrogensToAtomsMap[hydrogenID] = []
      
      // Register a new bond between the two carbons.
      let bond = SIMD2(
        atomList[0],
        atomList[1])
      insertedBonds.append(bond)
    }
    topology.insert(bonds: insertedBonds)
  }
}
