//
//  Sort.swift
//  HDL
//
//  Created by Philip Turner on 12/2/23.
//

import Dispatch
import QuartzCore

extension Topology {
  @discardableResult
  public mutating func sort() -> [UInt32] {
    let sorter = OctreeSorter(atoms: atoms)
    
//    let start = CACurrentMediaTime()
    let state = sorter.traverseHighLevels()
    let reordering = sorter.traverseLowLevels(state: state)
//    let end = CACurrentMediaTime()
//    debugProfile(start, end, "final")
    
//    let reordering = sorter.mortonReorderingDynamic()
    
//    let grid = sorter.oldCreateGrid()
//    let reordering = sorter.oldMortonReordering(grid: grid)
    
    let previousAtoms = atoms
    
    for i in reordering.indices {
      let mortonIndex = reordering[i]
      atoms[i] = previousAtoms[Int(mortonIndex)]
    }
    
    let inverted = OctreeSorter.invertOrder(reordering)
    for i in bonds.indices {
      let bond = SIMD2<Int>(truncatingIfNeeded: bonds[i])
      var newBond: SIMD2<UInt32> = .zero
      newBond[0] = inverted[bond[0]]
      newBond[1] = inverted[bond[1]]
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
    return inverted
  }
}
