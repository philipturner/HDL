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
      initialTypeRawValues: [UInt8],
      atomList: [UInt32]
    ) {
      for atomID in atomList {
        switch initialTypeRawValues[Int(atomID)] {
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
  
  mutating func resolveTwoWayCollisions() {
    var updates = [CollisionState?](
      repeating: nil, count: hydrogensToAtomsMap.count)
    
    // Validate the integrity of 'initialTypeRawValues'.
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
      
      guard initialTypeRawValues[i] == expectedRawValue else {
        fatalError("Incorrect raw value.")
      }
    }
    
    // This needs to be a separate function because of strange control flow,
    // including a 'return' statement.
    func loopIteration(
      sourceAtomID: UInt32,
      sourceAtomList: [UInt32]
    ) {
      let siteGeometry = SiteGeometry(
        initialTypeRawValues: initialTypeRawValues,
        atomList: sourceAtomList)
      
      var linkedList: [UInt32] = []
      
      if siteGeometry.bothBridgehead {
        let hydrogens = atomsToHydrogensMap[Int(sourceAtomList[0])]
        guard hydrogens.count == 1 else {
          fatalError("Bridgehead site did not have 1 hydrogen.")
        }
        
        linkedList.append(sourceAtomList[0])
        linkedList.append(hydrogens[0])
        linkedList.append(sourceAtomList[1])
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
        
        // Wrap this in a nested function, to make the control flow more
        // legible.
        func createLinkedList() -> [UInt32]? {
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
            let atomList2 = hydrogensToAtomsMap[nextHydrogen]
            if atomList2.count == 1 {
              // This is the end of the list.
              var atomListCopy = sourceAtomList
              atomListCopy.removeAll(where: { $0 == neighborAtomID })
              guard atomListCopy.count == 1 else {
                fatalError("This should never happen.")
              }
              
              var output: [UInt32] = []
              output.append(neighborAtomID)
              output.append(sortedHydrogens.first!)
              output.append(atomListCopy[0])
              return output
            }
          }
          
          return nil
        }
        
//      outer:
//        for neighborAtomID in sourceAtomList {
//          let unsortedHydrogens = atomsToHydrogensMap[Int(neighborAtomID)]
//          guard unsortedHydrogens.count == 2 else {
//            fatalError("Sidewall site did not have 2 hydrogens.")
//          }
//          
//          let sortedHydrogens = sort(hydrogens: unsortedHydrogens)
//          guard sourceAtomID == sortedHydrogens[0],
//                sourceAtomID != sortedHydrogens[1] else {
//            fatalError("Unexpected sorted hydrogen list.")
//          }
//          
//          let nextHydrogen = Int(sortedHydrogens.last!)
//          let atomList2 = hydrogensToAtomsMap[nextHydrogen]
//          if atomList2.count == 1 {
//            // This is the end of the list.
//            linkedList.append(neighborAtomID)
//            linkedList.append(sortedHydrogens.first!)
//            
//            var atomListCopy = sourceAtomList
//            atomListCopy.removeAll(where: { $0 == neighborAtomID })
//            guard atomListCopy.count == 1 else {
//              fatalError("This should never happen.")
//            }
//            
//            linkedList.append(atomListCopy[0])
//            break outer
//          }
//        }
        
        let createdLinkedList = createLinkedList()
        guard let createdLinkedList else {
          // Edge case: middle of a bond chain. This is never handled. If there
          // is a self-referential ring, the entire ring is skipped.
          return
        }
        
        linkedList = createdLinkedList
        
//        if linkedList.count == 0 {
//          // Edge case: middle of a bond chain. This is never handled. If there
//          // is a self-referential ring, the entire ring is skipped.
//          //
//          // This is one of the most important parts of the algorithm for
//          // detecting dimer chains (?)
//          return
//        } else if linkedList.count == 3 {
//          // Proceed with the remainder of the loop iteration.
//        } else {
//          fatalError("This should never happen.")
//        }
      } else if let bridgeheadID = siteGeometry.bridgeheadID,
                let sidewallID = siteGeometry.sidewallID {
        // The IDs of the elements are interleaved. First a carbon, then the
        // connecting collision, then a carbon, then a collision, etc. The end
        // must always be a carbon.
        linkedList.append(bridgeheadID)
        linkedList.append(sourceAtomID)
        linkedList.append(sidewallID)
      } else {
        fatalError("This should never happen.")
      }
      
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
          
          let centerType = initialTypeRawValues[Int(endOfList)]
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
    
    // High-level specification of the algorithm structure.
    for i in hydrogensToAtomsMap.indices {
      let atomList = hydrogensToAtomsMap[i]
      guard atomList.count == 2 else {
        continue
      }
      
      loopIteration(
        sourceAtomID: UInt32(i),
        sourceAtomList: atomList)
    }
    updateCollisions(updates.map { $0 ?? .keep })
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
