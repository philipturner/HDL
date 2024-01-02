//
//  GlobalScope.swift
//  HDL
//
//  Created by Philip Turner on 12/2/23.
//

/// Utility for checking that you're using each keyword in the correct scope.
enum GlobalScope {
  case lattice
  
  var description: String {
    switch self {
    case .lattice: return "Lattice"
    }
  }
  
  // You cannot instantiate an object while in the middle of defining
  // another one. This will cause an error because the global scope
  // doesn't have a FILO stack.
  static var global: GlobalScope? = nil
  
  static func push(_ newValue: GlobalScope) {
    guard global == nil else {
      fatalError("Pushed when already in a scope: \(global!)")
    }
    global = newValue
  }
  
  static func pop() -> GlobalScope {
    guard let global else {
      fatalError("Popped when not in a scope.")
    }
    Self.global = nil
    return global
  }
  
  static func throwUnrecognized(_ type: Any.Type) -> Never {
    if global == nil {
      fatalError("No global scope existed.")
    } else {
      let desc = global!.description
      fatalError("\(type) cannot be called in \(desc).")
    }
  }
}
