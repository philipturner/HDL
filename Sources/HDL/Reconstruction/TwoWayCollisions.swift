//
//  TwoWayCollisions.swift
//
//
//  Created by Philip Turner on 6/11/24.
//

extension Reconstruction {
  private enum CollisionState {
    case keep
    case bond
    case oneHydrogen(Int)
  }
  
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
  
  func validate(centerTypes: [UInt8]) {
    let orbitals = topology.nonbondingOrbitals(hybridization: .sp3)
    for i in orbitals.indices {
      let orbital = orbitals[i]
      var expectedRawValue: UInt8
      
      switch orbital.count {
      case 2:
        expectedRawValue = 2
      case 1:
        expectedRawValue = 3
      case 0:
        expectedRawValue = 4
      default:
        fatalError("This should never happen.")
      }
      
      guard centerTypes[i] == expectedRawValue else {
        fatalError("Incorrect raw value.")
      }
    }
  }
  
  mutating func resolveTwoWayCollisions(centerTypes: [UInt8]) {
    var updates = [CollisionState?](
      repeating: nil, count: hydrogensToAtomsMap.count)
    
    // High-level specification of the algorithm structure.
    for hydrogenID in hydrogensToAtomsMap.indices {
      let atomList = hydrogensToAtomsMap[hydrogenID]
      guard atomList.count == 2 else {
        continue
      }
      
      let dimerGeometry = DimerGeometry(
        centerTypes: centerTypes,
        atomList: atomList)
      let initialLinkedList = createInitialLinkedList(
        hydrogenID: UInt32(hydrogenID),
        dimerGeometry: dimerGeometry)
      guard let initialLinkedList else {
        continue
      }
      
      // The IDs of the elements are interleaved. First a carbon, then the
      // connecting collision, then a carbon, then a collision, etc. The end
      // must always be a carbon.
      //
      // Even indices: carbon sites
      // Odd indices:  hydrogen sites / collision sites
      let linkedList = expandLinkedList(
        initialLinkedList: initialLinkedList)
      for i in linkedList.indices where i % 2 == 1 {
        let listElement = Int(linkedList[i])
        if updates[listElement] == nil {
          if i % 4 == 1 {
            updates[listElement] = .bond
          } else {
            updates[listElement] = .keep
          }
        }
      }
    }
    
    // Replace undefined updates with '.keep'.
    let definedUpdates = updates.map { $0 ?? .keep }
    updateCollisions(definedUpdates)
  }
  
  
  private func createInitialLinkedList(
    hydrogenID: UInt32,
    dimerGeometry: DimerGeometry
  ) -> SIMD3<UInt32>? {
    // This list is guaranteed to have 2 elements.
    let atomList = hydrogensToAtomsMap[Int(hydrogenID)]
    
    if dimerGeometry.bothBridgehead {
      let hydrogens = atomsToHydrogensMap[Int(atomList[0])]
      guard hydrogens.count == 1 else {
        fatalError("Bridgehead site did not have 1 hydrogen.")
      }
      
      return SIMD3(
        atomList[0],
        hydrogens[0],
        atomList[1])
    } else if dimerGeometry.bothSidewall {
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
        // 3: guaranteed to be impossible (?)
        guard oppositeAtomList.count == 1 else {
          continue
        }
        
        let oppositeAtomID = Self.opposite(
          original: atomID,
          unsortedList: atomList)
        
        // Making sense of all the opposites:
        // (oppositeHydrogenID, atomID, hydrogenID, oppositeAtomID)
        // (H, C, H, C, ...)
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
      return SIMD3(
        bridgeheadID,
        hydrogenID,
        sidewallID)
    } else {
      fatalError("This should never happen.")
    }
  }
  
  private func expandLinkedList(
    initialLinkedList: SIMD3<UInt32>
  ) -> [UInt32] {
    var linkedList = [
      initialLinkedList[0],
      initialLinkedList[1],
      initialLinkedList[2],
    ]
    
    // Iteratively search through the topology, seeing whether the chain
    // of linked center atoms finally ends.
    var iterationCount = 0
  outer:
    while true {
      defer {
        iterationCount += 1
        if iterationCount > 1000 {
          fatalError("Took too many iterations to find length of dimer chain.")
        }
      }
      let endOfList = linkedList[linkedList.count - 1]
      let existingHydrogen = linkedList[linkedList.count - 2]
      
      // Change this to a loop structure that calls a function, which returns
      // whether or not the chain terminated.
      var appendedListElements: [UInt32] = []
      var hydrogens = atomsToHydrogensMap[Int(endOfList)]
      switch hydrogens.count {
      case 1:
        // We found a bridgehead site.
        precondition(
          hydrogens[0] == existingHydrogen,
          "Unexpected hydrogen list.")
        
        break outer
        
      case 2:
        // We found a sidewall site.
        
        // TODO: This is very similar to a function for sorting.
        if hydrogens[0] == existingHydrogen {
          
        } else if hydrogens[1] == existingHydrogen {
          hydrogens = [hydrogens[1], hydrogens[0]]
        } else {
          fatalError("Unexpected hydrogen list.")
        }
        precondition(hydrogens.first! == existingHydrogen)
        precondition(hydrogens.last! != existingHydrogen)
        
        let nextHydrogen = hydrogens.last!
        var atomList = hydrogensToAtomsMap[Int(nextHydrogen)]
        if atomList.count == 1 {
          break outer
        }
        
        // Add to the linked list.
        appendedListElements.append(nextHydrogen)
        
        // TODO: What do you mean, "this may not always be true"?
        // Why is this conditional statement needed?
        precondition(atomList.count == 2) // this may not always be true
        
        precondition(atomList.contains(endOfList))
        atomList.removeAll(where: { $0 == endOfList })
        precondition(atomList.count == 1)
        appendedListElements.append(atomList[0])
        
        break
      
      default:
        fatalError("This should never happen.")
      }
      
      linkedList += appendedListElements
    }
    
    return linkedList
  }
  
  private mutating func updateCollisions(_ states: [CollisionState]) {
    var insertedBonds: [SIMD2<UInt32>] = []
    
  outer:
    for i in states.indices {
      precondition(
        i >= 0 && i < hydrogensToAtomsMap.count,
        "Hydrogen index out of bounds.")
      
      switch states[i] {
      case .keep:
        continue outer
      case .bond:
        break
      case .oneHydrogen(_):
        fatalError("This update is not reported yet.")
      }
      
      let atomList = hydrogensToAtomsMap[Int(i)]
      if atomList.count != 2 {
        // TODO: Important debug diagnostic? Is there a better way to handle
        // this than undesired console writing?
        //
        // Perhaps use a precondition.
        print("Not a two-way collision: \(atomList)")
      }
      
      precondition(atomList.count == 2, "Not a two-way collision: \(atomList)")
      hydrogensToAtomsMap[Int(i)] = []
      
      let bond = SIMD2(atomList[0], atomList[1])
      insertedBonds.append(bond)
      
      for j in atomList {
        precondition(
          j >= 0 && j < atomsToHydrogensMap.count,
          "Atom index is out of bounds.")
        var previous = atomsToHydrogensMap[Int(j)]
        precondition(previous.count > 0, "Hydrogen map already empty.")
        
        var matchIndex = -1
        for k in previous.indices {
          if previous[k] == UInt32(i) {
            matchIndex = k
            break
          }
        }
        precondition(matchIndex != -1, "Could not find a match.")
        previous.remove(at: matchIndex)
        atomsToHydrogensMap[Int(j)] = previous
      }
    }
    topology.insert(bonds: insertedBonds)
  }
}
