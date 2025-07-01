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
  
  private struct SiteGeometry {
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
  
  func validateCenterTypes() {
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
  
  mutating func resolveTwoWayCollisions() {
    var updates = [CollisionState?](
      repeating: nil, count: hydrogensToAtomsMap.count)
    
    // High-level specification of the algorithm structure.
    for i in hydrogensToAtomsMap.indices {
      let atomList = hydrogensToAtomsMap[i]
      guard atomList.count == 2 else {
        continue
      }
      
      let initialLinkedList = createInitialLinkedList(
        sourceAtomID: UInt32(i))
      guard let initialLinkedList else {
        continue
      }
      
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
    sourceAtomID: UInt32
  ) -> SIMD3<UInt32>? {
    let sourceAtomList = hydrogensToAtomsMap[Int(sourceAtomID)]
    let siteGeometry = SiteGeometry(
      centerTypes: centerTypes,
      atomList: sourceAtomList)
    
    // Extract this into an isolated function, initializeLinkedList().
    // Isolate it as much as possible from the local scope of the function
    // 'resolveTwoWayCollisions()'.
    if siteGeometry.bothBridgehead {
      let hydrogens = atomsToHydrogensMap[Int(sourceAtomList[0])]
      guard hydrogens.count == 1 else {
        fatalError("Bridgehead site did not have 1 hydrogen.")
      }
      
      return SIMD3(
        sourceAtomList[0],
        hydrogens[0],
        sourceAtomList[1])
    } else if siteGeometry.bothSidewall {
      // Sort the hydrogens, so the source atom appears first.
      func sort(hydrogens: [UInt32]) -> [UInt32] {
        if sourceAtomID == hydrogens[0] {
          return [hydrogens[0], hydrogens[1]]
        } else if sourceAtomID == hydrogens[1] {
          return [hydrogens[1], hydrogens[0]]
        } else {
          fatalError("Unexpected hydrogen list.")
        }
      }
      
      func createSidewallList() -> SIMD3<UInt32>? {
        for neighborAtomID in sourceAtomList {
          let unsortedHydrogens = atomsToHydrogensMap[Int(neighborAtomID)]
          guard unsortedHydrogens.count == 2 else {
            fatalError("Sidewall site did not have 2 hydrogens.")
          }
          
          let sortedHydrogens = sort(hydrogens: unsortedHydrogens)
          guard sourceAtomID == sortedHydrogens[0],
                sourceAtomID != sortedHydrogens[1] else {
            fatalError("Unexpected sorted hydrogen list.")
          }
          
          let nextHydrogen = Int(sortedHydrogens.last!)
          let nextHydrogenList = hydrogensToAtomsMap[nextHydrogen]
          if nextHydrogenList.count == 1 {
            // This is the end of the dimer chain.
            var listCopy = sourceAtomList
            listCopy.removeAll(where: { $0 == neighborAtomID })
            guard listCopy.count == 1 else {
              fatalError("This should never happen.")
            }
            
            return SIMD3(
              neighborAtomID,
              sortedHydrogens.first!,
              listCopy[0])
          }
        }
        
        return nil
      }
      
      let createdLinkedList = createSidewallList()
      guard let createdLinkedList else {
        // Edge case: middle of a bond chain. This is never handled. If there
        // is a self-referential ring, the entire ring is skipped.
        return nil
      }
      
      return createdLinkedList
    } else if let bridgeheadID = siteGeometry.bridgeheadID,
              let sidewallID = siteGeometry.sidewallID {
      // The IDs of the elements are interleaved. First a carbon, then the
      // connecting collision, then a carbon, then a collision, etc. The end
      // must always be a carbon.
      return SIMD3(
        bridgeheadID,
        sourceAtomID,
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
      let endOfList = linkedList.last!
      let existingHydrogen = linkedList[linkedList.count - 2]
      
      var hydrogens = atomsToHydrogensMap[Int(endOfList)]
      switch hydrogens.count {
      case 1:
        // If this happens, the end of the list is a bridgehead carbon.
        precondition(
          hydrogens[0] == UInt32(existingHydrogen),
          "Unexpected hydrogen list.")
        
        let centerType = centerTypes[Int(endOfList)]
        precondition(centerType == 3, "Must be a bridgehead carbon.")
        break outer
      case 2:
        if hydrogens[0] == UInt32(existingHydrogen) {
          
        } else if hydrogens[1] == UInt32(existingHydrogen) {
          hydrogens = [hydrogens[1], hydrogens[0]]
        } else {
          fatalError("Unexpected hydrogen list.")
        }
        precondition(hydrogens.first! == UInt32(existingHydrogen))
        precondition(hydrogens.last! != UInt32(existingHydrogen))
        
        let nextHydrogen = hydrogens.last!
        var atomList = hydrogensToAtomsMap[Int(nextHydrogen)]
        if atomList.count == 1 {
          // This is the end of the list.
          break outer
        }
        linkedList.append(nextHydrogen)
        precondition(atomList.count == 2) // this may not always be true
        
        precondition(atomList.contains(UInt32(endOfList)))
        atomList.removeAll(where: { $0 == UInt32(endOfList) })
        precondition(atomList.count == 1)
        linkedList.append(atomList[0])
        
        break
      case 3:
        fatalError("Chain terminated at 3-way collision.")
      default:
        fatalError("Unexpected hydrogen count: \(hydrogens.count)")
      }
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
