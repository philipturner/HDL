//
//  Bond.swift
//  MolecularRenderer
//
//  Created by Philip Turner on 10/22/23.
//

public enum Bond: RawRepresentable {
  /// A bond with order 1.
  case sigma
  
  @inlinable @inline(__always)
  public init(rawValue: Float) {
    switch rawValue {
    case 1: self = .sigma
    default:
      fatalError("Invalid raw value for bond: \(rawValue)")
    }
  }
  
  @inlinable @inline(__always)
  public var rawValue: Float {
    switch self {
    case .sigma: return 1
    }
  }
}
