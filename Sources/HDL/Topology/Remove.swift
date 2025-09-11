//
//  Remove.swift
//  HDL
//
//  Created by Philip Turner on 7/5/25.
//

extension Topology {
  public mutating func remove(atoms indices: [UInt32]) {
    var marks = [Bool](repeating: true, count: atoms.count)
    guard indices.allSatisfy({ $0 < atoms.count }) else {
      fatalError("One or more atom indices were out of range.")
    }
    for index in indices {
      marks[Int(index)] = false
    }
    
    var newAtoms: [Atom] = []
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
