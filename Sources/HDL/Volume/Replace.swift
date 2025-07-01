//
//  Replace.swift
//  HDL
//
//  Created by Philip Turner on 9/1/23.
//

public enum ReplaceType {
  case atom(Element)
  case empty
}

@MainActor
public struct Replace {
  @discardableResult
  public init(_ closure: () -> ReplaceType) {
    guard GlobalScope.global == .lattice else {
      GlobalScope.throwUnrecognized(Self.self)
    }
    LatticeStack.touchGlobal()
    
    switch closure() {
    case .atom(let element):
      let code = Int8(element.atomicNumber)
      LatticeStack.global!.replace(with: code)
    case .empty:
      LatticeStack.global!.replace(with: .zero)
    }
  }
}
