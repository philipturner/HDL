//
//  Reconstruction.swift
//  HDL
//
//  Created by Philip Turner on 6/11/24.
//

public struct Reconstruction {
  /// Required.
  public var atoms: [SIMD4<Float>]?
  
  /// Required.
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
    let bonds = compilation.compile()
    
    var topology = Topology()
    topology.atoms = compilation.atoms
    topology.bonds = bonds
    return topology
  }
}
