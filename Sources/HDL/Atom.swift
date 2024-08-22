//
//  Atom.swift
//  HDL
//
//  Created by Philip Turner on 10/22/23.
//

public typealias Atom = SIMD4<Float>

extension Atom {
  @inlinable @inline(__always)
  public var position: SIMD3<Float> {
    @_transparent
    get {
      unsafeBitCast(self, to: SIMD3<Float>.self)
    }
    @_transparent
    set {
      self = SIMD4(newValue, self.w)
    }
  }
  
  @inlinable @inline(__always)
  public var atomicNumber: UInt8 {
    @_transparent
    get {
      let atomicNumber = UInt8(self.w)
      return atomicNumber
    }
    @_transparent
    set {
      let atomicNumber = newValue
      self.w = Float(atomicNumber)
    }
  }
  
  @inlinable @inline(__always)
  public var element: Element {
    @_transparent
    get {
      return Element(atomicNumber: UInt8(self.w))
    }
    @_transparent
    set {
      self.w = Float(newValue.atomicNumber)
    }
  }
  
  @_transparent
  @inlinable @inline(__always)
  public init(position: SIMD3<Float>, atomicNumber: UInt8) {
    self = SIMD4(position, Float(atomicNumber))
  }
  
  @_transparent
  @inlinable @inline(__always)
  public init(position: SIMD3<Float>, element: Element) {
    self = SIMD4(position, Float(element.atomicNumber))
  }
}
