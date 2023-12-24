//
//  Topology.swift
//  MolecularRenderer
//
//  Created by Philip Turner on 12/2/23.
//

public struct Topology {
  public internal(set) var atoms: [Entity] = []
  public internal(set) var bonds: [SIMD2<UInt32>] = []
  
  public init() {
    
  }
}

// MARK: - Insert and Remove

// TODO: Add test cases for these functions.
  
extension Topology {
  public mutating func insert(atoms: [Entity]) {
    for atom in atoms {
      guard atom.storage.w > 0 else {
        fatalError("Topology does not accept empty entities.")
      }
      guard atom.storage.w == Float(Int8(atom.storage.w)) else {
        fatalError("Entity must have an integer atomic number.")
      }
    }
    self.atoms += atoms
  }
  
  public mutating func insert(bonds: [SIMD2<UInt32>]) {
    let bondMax = UInt32(truncatingIfNeeded: atoms.count)
    for bond in bonds {
      guard all(bond .< bondMax) else {
        fatalError("Bond contained out-of-bounds atom index.")
      }
    }
    self.bonds += bonds
  }
  
  public mutating func remove(atoms indices: [UInt32]) {
    var marks = [Bool](repeating: true, count: atoms.count)
    guard indices.allSatisfy({ $0 < atoms.count }) else {
      fatalError("One or more atom indices were out of range.")
    }
    for index in indices {
      marks[Int(index)] = false
    }
    
    var newAtoms: [Entity] = []
    var reordering: [Int32] = []
    newAtoms.reserveCapacity(max(0, atoms.count - indices.count))
    reordering.reserveCapacity(atoms.count)
    for i in atoms.indices {
      var newIndex = Int32(newAtoms.count)
      if marks[i] {
        newAtoms.append(atoms[i])
      } else {
        newIndex = -1
      }
      reordering.append(newIndex)
    }
    self.atoms = newAtoms
    
    var removedBonds: [UInt32] = []
    for i in bonds.indices {
      let bond = SIMD2<Int>(truncatingIfNeeded: bonds[i])
      var newBond: SIMD2<Int32> = .zero
      newBond[0] = reordering[bond[0]]
      newBond[1] = reordering[bond[1]]
      bonds[i] = SIMD2(truncatingIfNeeded: newBond)
      
      if any(newBond .< 0) {
        removedBonds.append(UInt32(truncatingIfNeeded: i))
      }
    }
    self.remove(bonds: removedBonds)
  }
  
  public mutating func remove(bonds indices: [UInt32]) {
    var marks = [Bool](repeating: true, count: bonds.count)
    guard indices.allSatisfy({ $0 < bonds.count }) else {
      fatalError("One or more bond indices were out of range.")
    }
    for index in indices {
      marks[Int(index)] = false
    }
    
    var newBonds: [SIMD2<UInt32>] = []
    newBonds.reserveCapacity(max(0, bonds.count - indices.count))
    for i in bonds.indices {
      if marks[i] {
        newBonds.append(bonds[i])
      }
    }
    self.bonds = newBonds
  }
}

// MARK: - Sort

// TODO: Add test cases for Morton reordering and sorting.

extension Topology {
  @discardableResult
  public mutating func sort() -> [UInt32] {
    let grid = TopologyGrid(atoms: atoms)
    let reordering = grid.mortonReordering()
    
    for i in bonds.indices {
      let bond = SIMD2<Int>(truncatingIfNeeded: bonds[i])
      var newBond: SIMD2<UInt32> = .zero
      newBond[0] = reordering[bond[0]]
      newBond[1] = reordering[bond[1]]
      newBond = SIMD2(newBond.min(), newBond.max())
      bonds[i] = newBond
    }
    bonds.sort {
      if $0.x != $1.x {
        return $0.x < $1.x
      } else {
        return $0.y < $1.y
      }
    }
    return reordering
  }
}

// MARK: - Map

// TODO: Add test cases for mapping.

extension Topology {
  public enum MapType {
    case atoms
    case bonds
  }
  
  public func map(
    _ primaryType: MapType,
    to secondaryType: MapType
  ) -> [ArraySlice<UInt32>] {
    switch (primaryType, secondaryType) {
    case (.atoms, _):
      var connectionsMap = [SIMD8<Int32>](
        repeating: .init(repeating: -1), count: atoms.count)
      bonds.withUnsafeBufferPointer {
        let opaque = OpaquePointer($0.baseAddress.unsafelyUnwrapped)
        let casted = UnsafePointer<UInt32>(opaque)
        for i in 0..<2 * bonds.count {
          let atomID = Int(truncatingIfNeeded: casted[i])
          let bondID = Int32(truncatingIfNeeded: i / 2)
          var element = connectionsMap[atomID]
          for lane in 0..<8 {
            if element[lane] == -1 {
              element[lane] = bondID
              break
            }
          }
          connectionsMap[atomID] = element
        }
      }
      
      var outputArray: [UInt32] = []
      outputArray.reserveCapacity(atoms.count * 4)
      var outputRanges: [Range<UInt32>] = []
      outputRanges.reserveCapacity(atoms.count)
      var outputSlices: [ArraySlice<UInt32>] = []
      outputSlices.reserveCapacity(atoms.count)
      
      for atomID in atoms.indices {
        var element = connectionsMap[atomID]
        if secondaryType == .atoms {
          let atomID32 = UInt32(truncatingIfNeeded: atomID)
          for lane in 0..<8 {
            if element[lane] == -1 {
              break
            }
            
            let bondID = Int(truncatingIfNeeded: element[lane])
            let bond = bonds[bondID]
            guard any(bond .== atomID32) else {
              fatalError("This should never happen.")
            }
            let otherID32 = (bond[0] == atomID32) ? bond[1] : bond[0]
            element[lane] = Int32(truncatingIfNeeded: otherID32)
          }
        }
        
        let rangeStart = UInt32(truncatingIfNeeded: outputArray.count)
        for lane in 0..<8 {
          if element[lane] == -1 {
            break
          }
          outputArray.append(
            UInt32(truncatingIfNeeded: element[lane]))
        }
        let rangeEnd = UInt32(truncatingIfNeeded: outputArray.count)
        outputRanges.append(rangeStart..<rangeEnd)
      }
      
      
      for atomID in atoms.indices {
        let range = outputRanges[atomID]
        let rangeStart = Int(truncatingIfNeeded: range.lowerBound)
        let rangeEnd = Int(truncatingIfNeeded: range.upperBound)
        let slice = outputArray[rangeStart..<rangeEnd]
        outputSlices.append(slice)
      }
      return outputSlices
    case (.bonds, .atoms):
      var outputArray: [UInt32] = []
      outputArray.reserveCapacity(bonds.count * 2)
      for i in bonds.indices {
        let bond = bonds[i]
        outputArray.append(bond[0])
        outputArray.append(bond[1])
      }
      
      var outputSlices: [ArraySlice<UInt32>] = []
      outputSlices.reserveCapacity(bonds.count)
      for i in bonds.indices {
        let rangeStart = 2 * i
        let rangeEnd = 2 * (i + 1)
        let slice = outputArray[rangeStart..<rangeEnd]
        outputSlices.append(slice)
      }
      return outputSlices
    case (.bonds, .bonds):
      fatalError("Bonds to bonds map is not supported.")
    }
  }
}
