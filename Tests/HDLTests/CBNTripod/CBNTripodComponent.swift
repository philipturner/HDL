//
//  CBNTripodComponent.swift
//  HDLTests
//
//  Created by Philip Turner on 12/30/23.
//

import HDL

protocol CBNTripodComponent {
  var topology: Topology { get }
  
  // Call the compiler passes here, instead of doing them post-initialization.
  init()
}

extension CBNTripodComponent {
  // Check that all bonds are correctly assigned.
  func createBondRecord() -> [SIMD2<UInt8>: Int] {
    var bondRecord: [SIMD2<UInt8>: Int] = [:]
    for bond in topology.bonds {
      var atomicNumbers: [UInt8] = [
        topology.atoms[Int(bond[0])].atomicNumber,
        topology.atoms[Int(bond[1])].atomicNumber,
      ]
      atomicNumbers.sort()
      
      let key = SIMD2(atomicNumbers[0], atomicNumbers[1])
      var value = bondRecord[key] ?? 0
      value += 1
      bondRecord[key] = value
    }
    return bondRecord
  }
}
