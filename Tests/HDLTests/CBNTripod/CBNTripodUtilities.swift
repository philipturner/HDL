//
//  CBNTripodUtilities.swift
//  HDL
//
//  Created by Philip Turner on 7/19/25.
//

import QuaternionModule

struct CBNTripodUtilities {
  // TODO: Add another utility to encapsulate all cross products
  
  static func quaternion(
    from start: SIMD3<Float>,
    to end: SIMD3<Float>
  ) -> Quaternion<Float> {
    func cross(_ _self: SIMD3<Float>, _ other: SIMD3<Float>) -> SIMD3<Float> {
      let yzx = SIMD3<Int>(1,2,0)
      let zxy = SIMD3<Int>(2,0,1)
      return (_self[yzx] * other[zxy]) - (_self[zxy] * other[yzx])
    }
    
    // Source: https://stackoverflow.com/a/1171995
    let a = cross(start, end)
    let xyz = a
    let v1LengthSq = (start * start).sum()
    let v2LengthSq = (end * end).sum()
    let w = (v1LengthSq + v2LengthSq).squareRoot() + (start * end).sum()
    let output = Quaternion(real: w, imaginary: xyz)
    
    guard let normalized = output.normalized else {
      fatalError("Could not normalize the quaternion.")
    }
    return normalized
  }
}
