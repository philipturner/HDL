//
//  Topology+Match.swift
//
//
//  Created by Philip Turner on 12/10/23.
//

public typealias MatchType = (
  _ input: Entity,
  _ candidate: Entity
) -> Bool

extension Topology {
  public func match(
    _ input: [Entity]
  ) -> [UInt32?] {
    match(input) { _, _ in return true }
  }
  
  // TODO: Make a Swift package test for this function.
  public func match(
    _ input: [Entity],
    _ closure: MatchType
  ) -> [UInt32?] {
    fatalError("Not implemented.")
  }
}
