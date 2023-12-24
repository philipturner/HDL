//
//  Element.swift
//  MolecularRenderer
//
//  Created by Philip Turner on 10/22/23.
//

/// The elements supported by the compiler.
public enum Element: UInt8, CustomStringConvertible {
  case hydrogen = 1
  case carbon = 6
  case nitrogen = 7
  case oxygen = 8
  case fluorine = 9
  case silicon = 14
  case phosphorus = 15
  case sulfur = 16
  case germanium = 32
  case gold = 79
  
  @inlinable @inline(__always)
  public init(_ atomicNumber: UInt8) {
    self.init(rawValue: atomicNumber)!
  }
  
  public var description: String {
    switch self {
    case .hydrogen: return ".hydrogen"
    case .carbon: return ".carbon"
    case .nitrogen: return ".nitrogen"
    case .oxygen: return ".oxygen"
    case .fluorine: return ".fluorine"
    case .silicon: return ".silicon"
    case .phosphorus: return ".phosphorus"
    case .sulfur: return ".sulfur"
    case .germanium: return ".germanium"
    case .gold: return ".gold"
    }
  }
  
  // Private API for approximating covalent bond length. This is meant to be
  // called in bulk, even when whole-module optimization is disabled.
  @usableFromInline
  static let covalentRadii: [Float] = {
    // Source: https://periodictable.com/Properties/A/CovalentRadius.v.log.html
    var output = [Float](repeating: -1, count: 127)
    output[0] = 0 / 1000
    output[1] = 31 / 1000
    output[6] = 76 / 1000
    output[7] = 71 / 1000
    output[8] = 66 / 1000
    output[9] = 57 / 1000
    output[14] = 111 / 1000
    output[15] = 107 / 1000
    output[16] = 105 / 1000
    output[32] = 120 / 1000
    output[79] = 136 / 1000
    return output
  }()
  
  @inlinable @inline(__always)
  public var covalentRadius: Float {
    Self.covalentRadii[Int(rawValue)]
  }
}
