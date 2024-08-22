//
//  Element.swift
//  HDL
//
//  Created by Philip Turner on 10/22/23.
//

/// The elements parameterized by the compiler.
public enum Element: UInt8, CustomStringConvertible {
  case hydrogen = 1
  
  case boron = 5
  case carbon = 6
  case nitrogen = 7
  case oxygen = 8
  case fluorine = 9
  
  case aluminum = 13
  case silicon = 14
  case phosphorus = 15
  case sulfur = 16
  case chlorine = 17
  
  case gallium = 31
  case germanium = 32
  case arsenic = 33
  case selenium = 34
  case bromine = 35
  
  case tin = 50
  case gold = 79
  case lead = 82
  
  @_transparent
  @inlinable @inline(__always)
  public init(atomicNumber: UInt8) {
    self.init(rawValue: atomicNumber)!
  }
  
  public var description: String {
    switch self {
    case .hydrogen: return ".hydrogen"
      
    case .boron: return ".boron"
    case .carbon: return ".carbon"
    case .nitrogen: return ".nitrogen"
    case .oxygen: return ".oxygen"
    case .fluorine: return ".fluorine"
      
    case .aluminum: return ".aluminum"
    case .silicon: return ".silicon"
    case .phosphorus: return ".phosphorus"
    case .sulfur: return ".sulfur"
    case .chlorine: return ".chlorine"
      
    case .gallium: return ".gallium"
    case .germanium: return ".germanium"
    case .arsenic: return ".arsenic"
    case .selenium: return ".selenium"
    case .bromine: return ".bromine"
      
    case .tin: return ".tin"
    case .gold: return ".gold"
    case .lead: return ".lead"
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
    
    output[5] = 85 / 1000
    output[6] = 76 / 1000
    output[7] = 71 / 1000
    output[8] = 66 / 1000
    output[9] = 57 / 1000
    
    output[13] = 121 / 1000
    output[14] = 111 / 1000
    output[15] = 107 / 1000
    output[16] = 105 / 1000
    output[17] = 102 / 1000
    
    output[31] = 122 / 1000
    output[32] = 120 / 1000
    output[33] = 119 / 1000
    output[34] = 120 / 1000
    output[35] = 120 / 1000
    
    output[50] = 139 / 1000
    output[79] = 136 / 1000
    output[82] = 146 / 1000
    return output
  }()
  
  @inlinable @inline(__always)
  public var covalentRadius: Float {
    Self.covalentRadii[Int(rawValue)]
  }
  
  @_transparent
  @inlinable @inline(__always)
  public var atomicNumber: UInt8 {
    rawValue
  }
}
