//
//  Material.swift
//  HDL
//
//  Created by Philip Turner on 10/26/23.
//

public enum MaterialType {
  case elemental(Element)
  case checkerboard(Element, Element)
}

@MainActor
public struct Material {
  @discardableResult
  public init(_ closure: () -> MaterialType) {
    guard GlobalScope.global == .lattice else {
      GlobalScope.throwUnrecognized(Self.self)
    }
    
    let materialType = closure()
    switch materialType {
    case .elemental(let element):
      switch element {
      case .carbon: break
      case .silicon: break
      case .germanium: break
      case .gold: break
      default: fatalError("Unrecognized material type: \(materialType)")
      }
    case .checkerboard(let a, let b):
      let minElement = (a.rawValue < b.rawValue) ? a : b
      let maxElement = (a.rawValue > b.rawValue) ? a : b
      switch (minElement, maxElement) {
      case (.boron, .nitrogen): break
      case (.boron, .phosphorus): break
      case (.boron, .arsenic): break
      case (.carbon, .silicon): break
      case (.carbon, .germanium): break
      case (.nitrogen, .aluminum): break
      case (.nitrogen, .gallium): break
      
      case (.aluminum, .phosphorus): break
      case (.aluminum, .arsenic): break
      case (.phosphorus, .gallium): break
      case (.gallium, .arsenic): break
      default: fatalError("Unrecognized material type: \(materialType)")
      }
    }
    
    guard LatticeStackDescriptor.global.materialType == nil else {
      fatalError("Already set material.")
    }
    LatticeStackDescriptor.global.materialType = materialType
  }
}

