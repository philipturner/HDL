//
//  Reconstruction.swift
//
//
//  Created by Philip Turner on 6/11/24.
//

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
