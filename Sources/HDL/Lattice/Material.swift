//
//  Material.swift
//  HDL
//
//  Created by Philip Turner on 10/26/23.
//

// MARK: - Types

public enum MaterialType {
  case elemental(Element)
  case checkerboard(Element, Element)
}

public enum ConstantType {
  case hexagon
  case prism
  case square
  
  var latticeMultiplier: Float {
    switch self {
    case .hexagon:
      return Float(1.0 / 2).squareRoot()
    case .prism:
      return Float(4.0 / 3).squareRoot()
    case .square:
      return 1
    }
  }
}

// MARK: - Keywords

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

public typealias Constant = Float

extension Float {
  /// Do not use this initializer via the API `Float.init`. It is intended to be
  /// used through the API `Constant.init`.
  public init(
    _ constantType: ConstantType,
    _ closure: () -> MaterialType
  ) {
    let materialType = closure()
    var cubicSpacing: Float
    switch materialType {
    case .checkerboard(.boron, .nitrogen),
        .checkerboard(.nitrogen, .boron):
      cubicSpacing = 0.3615
    case .checkerboard(.boron, .phosphorus),
        .checkerboard(.phosphorus, .boron):
      cubicSpacing = 0.4538
    case .checkerboard(.boron, .arsenic),
        .checkerboard(.arsenic, .boron):
      cubicSpacing = 0.4778
      
    case .elemental(.carbon):
      cubicSpacing = 0.3567
    case .checkerboard(.carbon, .silicon),
        .checkerboard(.silicon, .carbon):
      cubicSpacing = 0.4360
    case .checkerboard(.carbon, .germanium),
        .checkerboard(.germanium, .carbon):
      cubicSpacing = 0.4523
      
    case .checkerboard(.nitrogen, .aluminum),
        .checkerboard(.aluminum, .nitrogen):
      cubicSpacing = 0.4358
    case .checkerboard(.nitrogen, .gallium),
        .checkerboard(.gallium, .nitrogen):
      cubicSpacing = 0.4498
      
    case .checkerboard(.aluminum, .phosphorus),
        .checkerboard(.phosphorus, .aluminum):
      cubicSpacing = 0.5464
    case .checkerboard(.aluminum, .arsenic),
        .checkerboard(.arsenic, .aluminum):
      cubicSpacing = 0.5660
      
    case .elemental(.silicon):
      cubicSpacing = 0.5431
      
    case .checkerboard(.phosphorus, .gallium),
        .checkerboard(.gallium, .phosphorus):
      cubicSpacing = 0.5450
      
    case .checkerboard(.gallium, .arsenic),
        .checkerboard(.arsenic, .gallium):
      cubicSpacing = 0.5653
      
    case .elemental(.germanium):
      cubicSpacing = 0.5658
      
    case .elemental(.gold):
      cubicSpacing = 0.4078
      guard constantType == .square else {
        fatalError("Hexagonal gold is unsupported.")
      }
      
    default:
      fatalError("Unrecognized material type: \(materialType)")
    }
    
    self = cubicSpacing * constantType.latticeMultiplier
  }
}
