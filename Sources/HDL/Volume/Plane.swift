//
//  Plane.swift
//  HDL
//
//  Created by Philip Turner on 10/29/23.
//

public struct Plane {
  @discardableResult
  public init(_ closure: () -> SIMD3<Float>) {
    guard GlobalScope.global == .lattice else {
      GlobalScope.throwUnrecognized(Self.self)
    }
    LatticeStack.touchGlobal()
    LatticeStack.global!.plane(normal: closure())
  }
}
