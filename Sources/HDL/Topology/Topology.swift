//
//  Topology.swift
//  MolecularRenderer
//
//  Created by Philip Turner on 12/2/23.
//

public struct Topology {
  var grid: TopologyGrid
  
  public init(_ entities: [Entity]) {
    self.init(entities) { }
  }
  
  public init(_ entities: [Entity], _ closure: () -> Void) {
    var mappedLocations: [SIMD2<UInt32>] = []
    self.grid = TopologyGrid(
      entities: entities, mappedLocations: &mappedLocations)
    closure()
  }
}
