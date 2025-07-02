//
//  Reconstruction.swift
//
//
//  Created by Philip Turner on 6/11/24.
//

// Revised user-facing API:
//
// struct Reconstruction {
//   var atoms: [SIMD4<Float>]?
//   var material: MaterialType?
//
//   func compile() -> Topology
// }
//
// var reconstruction = Reconstruction()
// reconstruction.atoms = ...
// reconstruction.material = ...
// let topology = reconstruction.compile()

public struct Reconstruction {
  public var atoms: [SIMD4<Float>]?
  public var material: MaterialType?
  
  public init() {
    
  }
  
  public func compile() -> Topology {
    guard let atoms else {
      fatalError("Did not specify atoms for surface reconstruction.")
    }
    guard let material else {
      fatalError("Did not specify material for surface reconstruction.")
    }
    
    var compilation = Compilation(
      atoms: atoms,
      material: material)
    compilation.compile()
    return compilation.topology
  }
}
