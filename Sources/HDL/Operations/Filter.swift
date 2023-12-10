//
//  Filter.swift
//  MolecularRenderer
//
//  Created by Philip Turner on 12/2/23.
//

public typealias FilterType = (
  _ atom: inout Entity,
  _ neighbors: inout [Entity]
) -> Void

public struct Filter {
  @discardableResult
  public init(_ closure: FilterType) {
    guard GlobalScope.global == .topology else {
      GlobalScope.throwUnrecognized(Self.self)
    }
    fatalError("Not implemented.")
  }
}
