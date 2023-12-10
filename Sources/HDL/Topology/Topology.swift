//
//  Topology.swift
//  MolecularRenderer
//
//  Created by Philip Turner on 12/2/23.
//

public struct Topology {
  var grid: TopologyGrid
  
  // Undocumented initializer to more ergonomically initialize a topology
  // without any filters.
  public init(_ entities: [Entity]) {
    self.init(entities) { }
  }
  
  public init(_ entities: [Entity], _ closure: () -> Void) {
    self.grid = TopologyGrid(entities: entities)
    closure()
  }
}
