//
//  Lonsdaleite.swift
//  HDLTests
//
//  Created by Philip Turner on 12/25/23.
//

import HDL

struct Lonsdaleite {
  var atoms: [Entity]
  var bonds: [SIMD2<UInt32>]
  
  init() {
    self.atoms = lonsdaleiteAtoms.map(Entity.init(storage:))
    self.bonds = lonsdaleiteBonds
  }
}

private let lonsdaleiteAtoms: [SIMD4<Float>] = [
  SIMD4<Float>(0.126, 0.073, 0.000, 6),
  SIMD4<Float>(0.126, -0.030, 0.036, 1),
  SIMD4<Float>(0.126, 0.073, -0.109, 1),
  SIMD4<Float>(0.000, 0.146, 0.051, 6),
  SIMD4<Float>(-0.089, 0.094, 0.015, 1),
  SIMD4<Float>(0.000, 0.146, 0.206, 6),
  SIMD4<Float>(-0.089, 0.094, 0.242, 1),
  SIMD4<Float>(0.378, 0.073, 0.000, 6),
  SIMD4<Float>(0.378, -0.030, 0.036, 1),
  SIMD4<Float>(0.378, 0.073, -0.109, 1),
  SIMD4<Float>(0.252, 0.146, 0.051, 6),
  SIMD4<Float>(0.252, 0.146, 0.206, 6),
  SIMD4<Float>(0.126, 0.364, 0.051, 6),
  SIMD4<Float>(0.000, 0.291, 0.000, 6),
  SIMD4<Float>(0.000, 0.291, -0.109, 1),
  SIMD4<Float>(-0.089, 0.343, 0.036, 1),
  SIMD4<Float>(0.126, 0.364, 0.206, 6),
  SIMD4<Float>(0.252, 0.291, 0.000, 6),
  SIMD4<Float>(0.252, 0.291, -0.109, 1),
  SIMD4<Float>(0.378, 0.364, 0.051, 6),
  SIMD4<Float>(0.378, 0.364, 0.206, 6),
  SIMD4<Float>(0.126, 0.073, 0.257, 6),
  SIMD4<Float>(0.126, -0.030, 0.221, 1),
  SIMD4<Float>(0.126, 0.073, 0.412, 6),
  SIMD4<Float>(0.126, -0.030, 0.448, 1),
  SIMD4<Float>(0.000, 0.146, 0.463, 6),
  SIMD4<Float>(-0.089, 0.094, 0.427, 1),
  SIMD4<Float>(0.378, 0.073, 0.257, 6),
  SIMD4<Float>(0.378, -0.030, 0.221, 1),
  SIMD4<Float>(0.378, 0.073, 0.412, 6),
  SIMD4<Float>(0.378, -0.030, 0.448, 1),
  SIMD4<Float>(0.252, 0.146, 0.463, 6),
  SIMD4<Float>(0.000, 0.291, 0.257, 6),
  SIMD4<Float>(-0.089, 0.343, 0.221, 1),
  SIMD4<Float>(0.126, 0.364, 0.463, 6),
  SIMD4<Float>(0.000, 0.291, 0.412, 6),
  SIMD4<Float>(-0.089, 0.343, 0.448, 1),
  SIMD4<Float>(0.252, 0.291, 0.257, 6),
  SIMD4<Float>(0.252, 0.291, 0.412, 6),
  SIMD4<Float>(0.378, 0.364, 0.463, 6),
  SIMD4<Float>(0.631, 0.073, 0.000, 6),
  SIMD4<Float>(0.631, -0.030, 0.036, 1),
  SIMD4<Float>(0.631, 0.073, -0.109, 1),
  SIMD4<Float>(0.504, 0.146, 0.051, 6),
  SIMD4<Float>(0.504, 0.146, 0.206, 6),
  SIMD4<Float>(0.883, 0.073, 0.000, 6),
  SIMD4<Float>(0.883, -0.030, 0.036, 1),
  SIMD4<Float>(0.883, 0.073, -0.109, 1),
  SIMD4<Float>(0.757, 0.146, 0.051, 6),
  SIMD4<Float>(0.757, 0.146, 0.206, 6),
  SIMD4<Float>(0.504, 0.291, 0.000, 6),
  SIMD4<Float>(0.504, 0.291, -0.109, 1),
  SIMD4<Float>(0.631, 0.364, 0.051, 6),
  SIMD4<Float>(0.631, 0.364, 0.206, 6),
  SIMD4<Float>(0.883, 0.364, 0.051, 6),
  SIMD4<Float>(0.757, 0.291, 0.000, 6),
  SIMD4<Float>(0.757, 0.291, -0.109, 1),
  SIMD4<Float>(0.883, 0.364, 0.206, 6),
  SIMD4<Float>(0.631, 0.073, 0.257, 6),
  SIMD4<Float>(0.631, -0.030, 0.221, 1),
  SIMD4<Float>(0.631, 0.073, 0.412, 6),
  SIMD4<Float>(0.631, -0.030, 0.448, 1),
  SIMD4<Float>(0.504, 0.146, 0.463, 6),
  SIMD4<Float>(0.883, 0.073, 0.257, 6),
  SIMD4<Float>(0.883, -0.030, 0.221, 1),
  SIMD4<Float>(0.883, 0.073, 0.412, 6),
  SIMD4<Float>(0.883, -0.030, 0.448, 1),
  SIMD4<Float>(0.757, 0.146, 0.463, 6),
  SIMD4<Float>(0.504, 0.291, 0.257, 6),
  SIMD4<Float>(0.504, 0.291, 0.412, 6),
  SIMD4<Float>(0.631, 0.364, 0.463, 6),
  SIMD4<Float>(0.757, 0.291, 0.257, 6),
  SIMD4<Float>(0.883, 0.364, 0.463, 6),
  SIMD4<Float>(0.757, 0.291, 0.412, 6),
  SIMD4<Float>(0.126, 0.510, 0.000, 6),
  SIMD4<Float>(0.126, 0.510, -0.109, 1),
  SIMD4<Float>(0.000, 0.728, 0.000, 6),
  SIMD4<Float>(0.000, 0.728, -0.109, 1),
  SIMD4<Float>(-0.089, 0.779, 0.036, 1),
  SIMD4<Float>(0.000, 0.582, 0.051, 6),
  SIMD4<Float>(-0.089, 0.531, 0.015, 1),
  SIMD4<Float>(0.000, 0.582, 0.206, 6),
  SIMD4<Float>(-0.089, 0.531, 0.242, 1),
  SIMD4<Float>(0.378, 0.510, 0.000, 6),
  SIMD4<Float>(0.378, 0.510, -0.109, 1),
  SIMD4<Float>(0.252, 0.582, 0.051, 6),
  SIMD4<Float>(0.252, 0.728, 0.000, 6),
  SIMD4<Float>(0.252, 0.728, -0.109, 1),
  SIMD4<Float>(0.252, 0.582, 0.206, 6),
  SIMD4<Float>(0.126, 0.801, 0.051, 6),
  SIMD4<Float>(0.126, 0.801, 0.206, 6),
  SIMD4<Float>(0.126, 0.947, 0.000, 6),
  SIMD4<Float>(0.126, 0.947, -0.109, 1),
  SIMD4<Float>(0.378, 0.947, 0.000, 6),
  SIMD4<Float>(0.378, 0.947, -0.109, 1),
  SIMD4<Float>(0.378, 0.801, 0.051, 6),
  SIMD4<Float>(0.378, 0.801, 0.206, 6),
  SIMD4<Float>(0.126, 0.510, 0.257, 6),
  SIMD4<Float>(0.000, 0.728, 0.257, 6),
  SIMD4<Float>(-0.089, 0.779, 0.221, 1),
  SIMD4<Float>(0.126, 0.510, 0.412, 6),
  SIMD4<Float>(0.000, 0.728, 0.412, 6),
  SIMD4<Float>(-0.089, 0.779, 0.448, 1),
  SIMD4<Float>(0.000, 0.582, 0.463, 6),
  SIMD4<Float>(-0.089, 0.531, 0.427, 1),
  SIMD4<Float>(0.378, 0.510, 0.257, 6),
  SIMD4<Float>(0.252, 0.728, 0.257, 6),
  SIMD4<Float>(0.378, 0.510, 0.412, 6),
  SIMD4<Float>(0.252, 0.582, 0.463, 6),
  SIMD4<Float>(0.252, 0.728, 0.412, 6),
  SIMD4<Float>(0.126, 0.947, 0.257, 6),
  SIMD4<Float>(0.126, 0.801, 0.463, 6),
  SIMD4<Float>(0.126, 0.947, 0.412, 6),
  SIMD4<Float>(0.378, 0.947, 0.257, 6),
  SIMD4<Float>(0.378, 0.947, 0.412, 6),
  SIMD4<Float>(0.378, 0.801, 0.463, 6),
  SIMD4<Float>(0.631, 0.510, 0.000, 6),
  SIMD4<Float>(0.631, 0.510, -0.109, 1),
  SIMD4<Float>(0.504, 0.582, 0.051, 6),
  SIMD4<Float>(0.504, 0.582, 0.206, 6),
  SIMD4<Float>(0.504, 0.728, 0.000, 6),
  SIMD4<Float>(0.504, 0.728, -0.109, 1),
  SIMD4<Float>(0.883, 0.510, 0.000, 6),
  SIMD4<Float>(0.883, 0.510, -0.109, 1),
  SIMD4<Float>(0.757, 0.728, 0.000, 6),
  SIMD4<Float>(0.757, 0.728, -0.109, 1),
  SIMD4<Float>(0.757, 0.582, 0.051, 6),
  SIMD4<Float>(0.757, 0.582, 0.206, 6),
  SIMD4<Float>(0.631, 0.801, 0.051, 6),
  SIMD4<Float>(0.631, 0.947, 0.000, 6),
  SIMD4<Float>(0.631, 0.947, -0.109, 1),
  SIMD4<Float>(0.631, 0.801, 0.206, 6),
  SIMD4<Float>(0.883, 0.801, 0.051, 6),
  SIMD4<Float>(0.883, 0.801, 0.206, 6),
  SIMD4<Float>(0.883, 0.947, 0.000, 6),
  SIMD4<Float>(0.883, 0.947, -0.109, 1),
  SIMD4<Float>(0.631, 0.510, 0.257, 6),
  SIMD4<Float>(0.504, 0.728, 0.257, 6),
  SIMD4<Float>(0.631, 0.510, 0.412, 6),
  SIMD4<Float>(0.504, 0.582, 0.463, 6),
  SIMD4<Float>(0.504, 0.728, 0.412, 6),
  SIMD4<Float>(0.883, 0.510, 0.257, 6),
  SIMD4<Float>(0.757, 0.728, 0.257, 6),
  SIMD4<Float>(0.883, 0.510, 0.412, 6),
  SIMD4<Float>(0.757, 0.728, 0.412, 6),
  SIMD4<Float>(0.757, 0.582, 0.463, 6),
  SIMD4<Float>(0.631, 0.947, 0.257, 6),
  SIMD4<Float>(0.631, 0.801, 0.463, 6),
  SIMD4<Float>(0.631, 0.947, 0.412, 6),
  SIMD4<Float>(0.883, 0.947, 0.257, 6),
  SIMD4<Float>(0.883, 0.801, 0.463, 6),
  SIMD4<Float>(0.883, 0.947, 0.412, 6),
  SIMD4<Float>(0.000, 0.146, 0.618, 6),
  SIMD4<Float>(-0.089, 0.094, 0.654, 1),
  SIMD4<Float>(0.126, 0.073, 0.669, 6),
  SIMD4<Float>(0.126, -0.030, 0.633, 1),
  SIMD4<Float>(0.378, 0.073, 0.669, 6),
  SIMD4<Float>(0.378, -0.030, 0.633, 1),
  SIMD4<Float>(0.252, 0.146, 0.618, 6),
  SIMD4<Float>(0.126, 0.364, 0.618, 6),
  SIMD4<Float>(0.000, 0.291, 0.669, 6),
  SIMD4<Float>(-0.089, 0.343, 0.633, 1),
  SIMD4<Float>(0.252, 0.291, 0.669, 6),
  SIMD4<Float>(0.378, 0.364, 0.618, 6),
  SIMD4<Float>(0.126, 0.073, 0.824, 6),
  SIMD4<Float>(0.126, -0.030, 0.860, 1),
  SIMD4<Float>(0.000, 0.146, 0.875, 6),
  SIMD4<Float>(-0.089, 0.094, 0.839, 1),
  SIMD4<Float>(0.378, 0.073, 0.824, 6),
  SIMD4<Float>(0.378, -0.030, 0.860, 1),
  SIMD4<Float>(0.252, 0.146, 0.875, 6),
  SIMD4<Float>(0.126, 0.364, 0.875, 6),
  SIMD4<Float>(0.000, 0.291, 0.824, 6),
  SIMD4<Float>(-0.089, 0.343, 0.860, 1),
  SIMD4<Float>(0.252, 0.291, 0.824, 6),
  SIMD4<Float>(0.378, 0.364, 0.875, 6),
  SIMD4<Float>(0.631, 0.073, 0.669, 6),
  SIMD4<Float>(0.631, -0.030, 0.633, 1),
  SIMD4<Float>(0.504, 0.146, 0.618, 6),
  SIMD4<Float>(0.757, 0.146, 0.618, 6),
  SIMD4<Float>(0.883, 0.073, 0.669, 6),
  SIMD4<Float>(0.883, -0.030, 0.633, 1),
  SIMD4<Float>(0.504, 0.291, 0.669, 6),
  SIMD4<Float>(0.631, 0.364, 0.618, 6),
  SIMD4<Float>(0.883, 0.364, 0.618, 6),
  SIMD4<Float>(0.757, 0.291, 0.669, 6),
  SIMD4<Float>(0.631, 0.073, 0.824, 6),
  SIMD4<Float>(0.631, -0.030, 0.860, 1),
  SIMD4<Float>(0.504, 0.146, 0.875, 6),
  SIMD4<Float>(0.883, 0.073, 0.824, 6),
  SIMD4<Float>(0.883, -0.030, 0.860, 1),
  SIMD4<Float>(0.757, 0.146, 0.875, 6),
  SIMD4<Float>(0.504, 0.291, 0.824, 6),
  SIMD4<Float>(0.631, 0.364, 0.875, 6),
  SIMD4<Float>(0.883, 0.364, 0.875, 6),
  SIMD4<Float>(0.757, 0.291, 0.824, 6),
  SIMD4<Float>(0.000, 0.582, 0.618, 6),
  SIMD4<Float>(-0.089, 0.531, 0.654, 1),
  SIMD4<Float>(0.126, 0.510, 0.669, 6),
  SIMD4<Float>(0.000, 0.728, 0.669, 6),
  SIMD4<Float>(-0.089, 0.779, 0.633, 1),
  SIMD4<Float>(0.378, 0.510, 0.669, 6),
  SIMD4<Float>(0.252, 0.582, 0.618, 6),
  SIMD4<Float>(0.252, 0.728, 0.669, 6),
  SIMD4<Float>(0.126, 0.801, 0.618, 6),
  SIMD4<Float>(0.126, 0.947, 0.669, 6),
  SIMD4<Float>(0.378, 0.801, 0.618, 6),
  SIMD4<Float>(0.378, 0.947, 0.669, 6),
  SIMD4<Float>(0.126, 0.510, 0.824, 6),
  SIMD4<Float>(0.000, 0.728, 0.824, 6),
  SIMD4<Float>(-0.089, 0.779, 0.860, 1),
  SIMD4<Float>(0.000, 0.582, 0.875, 6),
  SIMD4<Float>(-0.089, 0.531, 0.839, 1),
  SIMD4<Float>(0.378, 0.510, 0.824, 6),
  SIMD4<Float>(0.252, 0.582, 0.875, 6),
  SIMD4<Float>(0.252, 0.728, 0.824, 6),
  SIMD4<Float>(0.126, 0.801, 0.875, 6),
  SIMD4<Float>(0.126, 0.947, 0.824, 6),
  SIMD4<Float>(0.378, 0.947, 0.824, 6),
  SIMD4<Float>(0.378, 0.801, 0.875, 6),
  SIMD4<Float>(0.631, 0.510, 0.669, 6),
  SIMD4<Float>(0.504, 0.582, 0.618, 6),
  SIMD4<Float>(0.504, 0.728, 0.669, 6),
  SIMD4<Float>(0.757, 0.582, 0.618, 6),
  SIMD4<Float>(0.883, 0.510, 0.669, 6),
  SIMD4<Float>(0.757, 0.728, 0.669, 6),
  SIMD4<Float>(0.631, 0.801, 0.618, 6),
  SIMD4<Float>(0.631, 0.947, 0.669, 6),
  SIMD4<Float>(0.883, 0.801, 0.618, 6),
  SIMD4<Float>(0.883, 0.947, 0.669, 6),
  SIMD4<Float>(0.631, 0.510, 0.824, 6),
  SIMD4<Float>(0.504, 0.582, 0.875, 6),
  SIMD4<Float>(0.504, 0.728, 0.824, 6),
  SIMD4<Float>(0.883, 0.510, 0.824, 6),
  SIMD4<Float>(0.757, 0.728, 0.824, 6),
  SIMD4<Float>(0.757, 0.582, 0.875, 6),
  SIMD4<Float>(0.631, 0.801, 0.875, 6),
  SIMD4<Float>(0.631, 0.947, 0.824, 6),
  SIMD4<Float>(0.883, 0.801, 0.875, 6),
  SIMD4<Float>(0.883, 0.947, 0.824, 6),
  SIMD4<Float>(1.009, 0.146, 0.051, 6),
  SIMD4<Float>(1.098, 0.094, 0.015, 1),
  SIMD4<Float>(1.009, 0.146, 0.206, 6),
  SIMD4<Float>(1.098, 0.094, 0.242, 1),
  SIMD4<Float>(1.009, 0.291, 0.000, 6),
  SIMD4<Float>(1.098, 0.343, 0.036, 1),
  SIMD4<Float>(1.009, 0.291, -0.109, 1),
  SIMD4<Float>(1.009, 0.146, 0.463, 6),
  SIMD4<Float>(1.098, 0.094, 0.427, 1),
  SIMD4<Float>(1.009, 0.291, 0.257, 6),
  SIMD4<Float>(1.098, 0.343, 0.221, 1),
  SIMD4<Float>(1.009, 0.291, 0.412, 6),
  SIMD4<Float>(1.098, 0.343, 0.448, 1),
  SIMD4<Float>(1.009, 0.582, 0.051, 6),
  SIMD4<Float>(1.098, 0.531, 0.015, 1),
  SIMD4<Float>(1.009, 0.728, 0.000, 6),
  SIMD4<Float>(1.098, 0.779, 0.036, 1),
  SIMD4<Float>(1.009, 0.728, -0.109, 1),
  SIMD4<Float>(1.009, 0.582, 0.206, 6),
  SIMD4<Float>(1.098, 0.531, 0.242, 1),
  SIMD4<Float>(1.009, 0.728, 0.257, 6),
  SIMD4<Float>(1.098, 0.779, 0.221, 1),
  SIMD4<Float>(1.009, 0.582, 0.463, 6),
  SIMD4<Float>(1.098, 0.531, 0.427, 1),
  SIMD4<Float>(1.009, 0.728, 0.412, 6),
  SIMD4<Float>(1.098, 0.779, 0.448, 1),
  SIMD4<Float>(1.009, 0.146, 0.618, 6),
  SIMD4<Float>(1.098, 0.094, 0.654, 1),
  SIMD4<Float>(1.009, 0.291, 0.669, 6),
  SIMD4<Float>(1.098, 0.343, 0.633, 1),
  SIMD4<Float>(1.009, 0.146, 0.875, 6),
  SIMD4<Float>(1.098, 0.094, 0.839, 1),
  SIMD4<Float>(1.009, 0.291, 0.824, 6),
  SIMD4<Float>(1.098, 0.343, 0.860, 1),
  SIMD4<Float>(1.009, 0.582, 0.618, 6),
  SIMD4<Float>(1.098, 0.531, 0.654, 1),
  SIMD4<Float>(1.009, 0.728, 0.669, 6),
  SIMD4<Float>(1.098, 0.779, 0.633, 1),
  SIMD4<Float>(1.009, 0.582, 0.875, 6),
  SIMD4<Float>(1.098, 0.531, 0.839, 1),
  SIMD4<Float>(1.009, 0.728, 0.824, 6),
  SIMD4<Float>(1.098, 0.779, 0.860, 1),
  SIMD4<Float>(0.126, 1.238, 0.051, 6),
  SIMD4<Float>(0.126, 1.341, 0.015, 1),
  SIMD4<Float>(0.000, 1.165, 0.000, 6),
  SIMD4<Float>(-0.089, 1.216, 0.036, 1),
  SIMD4<Float>(0.000, 1.165, -0.109, 1),
  SIMD4<Float>(0.000, 1.019, 0.051, 6),
  SIMD4<Float>(-0.089, 0.968, 0.015, 1),
  SIMD4<Float>(0.000, 1.019, 0.206, 6),
  SIMD4<Float>(-0.089, 0.968, 0.242, 1),
  SIMD4<Float>(0.126, 1.238, 0.206, 6),
  SIMD4<Float>(0.126, 1.341, 0.242, 1),
  SIMD4<Float>(0.252, 1.019, 0.051, 6),
  SIMD4<Float>(0.252, 1.165, 0.000, 6),
  SIMD4<Float>(0.252, 1.165, -0.109, 1),
  SIMD4<Float>(0.252, 1.019, 0.206, 6),
  SIMD4<Float>(0.378, 1.238, 0.051, 6),
  SIMD4<Float>(0.378, 1.341, 0.015, 1),
  SIMD4<Float>(0.378, 1.238, 0.206, 6),
  SIMD4<Float>(0.378, 1.341, 0.242, 1),
  SIMD4<Float>(0.000, 1.165, 0.257, 6),
  SIMD4<Float>(-0.089, 1.216, 0.221, 1),
  SIMD4<Float>(0.126, 1.238, 0.463, 6),
  SIMD4<Float>(0.126, 1.341, 0.427, 1),
  SIMD4<Float>(0.000, 1.165, 0.412, 6),
  SIMD4<Float>(-0.089, 1.216, 0.448, 1),
  SIMD4<Float>(0.000, 1.019, 0.463, 6),
  SIMD4<Float>(-0.089, 0.968, 0.427, 1),
  SIMD4<Float>(0.252, 1.165, 0.257, 6),
  SIMD4<Float>(0.252, 1.019, 0.463, 6),
  SIMD4<Float>(0.252, 1.165, 0.412, 6),
  SIMD4<Float>(0.378, 1.238, 0.463, 6),
  SIMD4<Float>(0.378, 1.341, 0.427, 1),
  SIMD4<Float>(0.504, 1.019, 0.051, 6),
  SIMD4<Float>(0.504, 1.019, 0.206, 6),
  SIMD4<Float>(0.504, 1.165, 0.000, 6),
  SIMD4<Float>(0.504, 1.165, -0.109, 1),
  SIMD4<Float>(0.631, 1.238, 0.051, 6),
  SIMD4<Float>(0.631, 1.341, 0.015, 1),
  SIMD4<Float>(0.631, 1.238, 0.206, 6),
  SIMD4<Float>(0.631, 1.341, 0.242, 1),
  SIMD4<Float>(0.883, 1.238, 0.051, 6),
  SIMD4<Float>(0.883, 1.341, 0.015, 1),
  SIMD4<Float>(0.757, 1.165, 0.000, 6),
  SIMD4<Float>(0.757, 1.165, -0.109, 1),
  SIMD4<Float>(0.757, 1.019, 0.051, 6),
  SIMD4<Float>(0.757, 1.019, 0.206, 6),
  SIMD4<Float>(0.883, 1.238, 0.206, 6),
  SIMD4<Float>(0.883, 1.341, 0.242, 1),
  SIMD4<Float>(0.504, 1.165, 0.257, 6),
  SIMD4<Float>(0.504, 1.019, 0.463, 6),
  SIMD4<Float>(0.504, 1.165, 0.412, 6),
  SIMD4<Float>(0.631, 1.238, 0.463, 6),
  SIMD4<Float>(0.631, 1.341, 0.427, 1),
  SIMD4<Float>(0.757, 1.165, 0.257, 6),
  SIMD4<Float>(0.883, 1.238, 0.463, 6),
  SIMD4<Float>(0.883, 1.341, 0.427, 1),
  SIMD4<Float>(0.757, 1.165, 0.412, 6),
  SIMD4<Float>(0.757, 1.019, 0.463, 6),
  SIMD4<Float>(0.000, 1.019, 0.618, 6),
  SIMD4<Float>(-0.089, 0.968, 0.654, 1),
  SIMD4<Float>(0.126, 1.238, 0.618, 6),
  SIMD4<Float>(0.126, 1.341, 0.654, 1),
  SIMD4<Float>(0.000, 1.165, 0.669, 6),
  SIMD4<Float>(-0.089, 1.216, 0.633, 1),
  SIMD4<Float>(0.252, 1.019, 0.618, 6),
  SIMD4<Float>(0.252, 1.165, 0.669, 6),
  SIMD4<Float>(0.378, 1.238, 0.618, 6),
  SIMD4<Float>(0.378, 1.341, 0.654, 1),
  SIMD4<Float>(0.126, 1.238, 0.875, 6),
  SIMD4<Float>(0.126, 1.341, 0.839, 1),
  SIMD4<Float>(0.000, 1.165, 0.824, 6),
  SIMD4<Float>(-0.089, 1.216, 0.860, 1),
  SIMD4<Float>(0.000, 1.019, 0.875, 6),
  SIMD4<Float>(-0.089, 0.968, 0.839, 1),
  SIMD4<Float>(0.252, 1.019, 0.875, 6),
  SIMD4<Float>(0.252, 1.165, 0.824, 6),
  SIMD4<Float>(0.378, 1.238, 0.875, 6),
  SIMD4<Float>(0.378, 1.341, 0.839, 1),
  SIMD4<Float>(0.504, 1.019, 0.618, 6),
  SIMD4<Float>(0.504, 1.165, 0.669, 6),
  SIMD4<Float>(0.631, 1.238, 0.618, 6),
  SIMD4<Float>(0.631, 1.341, 0.654, 1),
  SIMD4<Float>(0.757, 1.019, 0.618, 6),
  SIMD4<Float>(0.883, 1.238, 0.618, 6),
  SIMD4<Float>(0.883, 1.341, 0.654, 1),
  SIMD4<Float>(0.757, 1.165, 0.669, 6),
  SIMD4<Float>(0.504, 1.019, 0.875, 6),
  SIMD4<Float>(0.504, 1.165, 0.824, 6),
  SIMD4<Float>(0.631, 1.238, 0.875, 6),
  SIMD4<Float>(0.631, 1.341, 0.839, 1),
  SIMD4<Float>(0.883, 1.238, 0.875, 6),
  SIMD4<Float>(0.883, 1.341, 0.839, 1),
  SIMD4<Float>(0.757, 1.165, 0.824, 6),
  SIMD4<Float>(0.757, 1.019, 0.875, 6),
  SIMD4<Float>(1.009, 1.019, 0.051, 6),
  SIMD4<Float>(1.098, 0.968, 0.015, 1),
  SIMD4<Float>(1.009, 1.165, 0.000, 6),
  SIMD4<Float>(1.009, 1.165, -0.109, 1),
  SIMD4<Float>(1.098, 1.216, 0.036, 1),
  SIMD4<Float>(1.009, 1.019, 0.206, 6),
  SIMD4<Float>(1.098, 0.968, 0.242, 1),
  SIMD4<Float>(1.009, 1.165, 0.257, 6),
  SIMD4<Float>(1.098, 1.216, 0.221, 1),
  SIMD4<Float>(1.009, 1.019, 0.463, 6),
  SIMD4<Float>(1.098, 0.968, 0.427, 1),
  SIMD4<Float>(1.009, 1.165, 0.412, 6),
  SIMD4<Float>(1.098, 1.216, 0.448, 1),
  SIMD4<Float>(1.009, 1.019, 0.618, 6),
  SIMD4<Float>(1.098, 0.968, 0.654, 1),
  SIMD4<Float>(1.009, 1.165, 0.669, 6),
  SIMD4<Float>(1.098, 1.216, 0.633, 1),
  SIMD4<Float>(1.009, 1.019, 0.875, 6),
  SIMD4<Float>(1.098, 0.968, 0.839, 1),
  SIMD4<Float>(1.009, 1.165, 0.824, 6),
  SIMD4<Float>(1.098, 1.216, 0.860, 1),
  SIMD4<Float>(0.000, 0.146, 1.030, 6),
  SIMD4<Float>(-0.089, 0.094, 1.066, 1),
  SIMD4<Float>(0.126, 0.073, 1.081, 6),
  SIMD4<Float>(0.126, 0.073, 1.190, 1),
  SIMD4<Float>(0.126, -0.030, 1.045, 1),
  SIMD4<Float>(0.378, 0.073, 1.081, 6),
  SIMD4<Float>(0.378, 0.073, 1.190, 1),
  SIMD4<Float>(0.378, -0.030, 1.045, 1),
  SIMD4<Float>(0.252, 0.146, 1.030, 6),
  SIMD4<Float>(0.126, 0.364, 1.030, 6),
  SIMD4<Float>(0.000, 0.291, 1.081, 6),
  SIMD4<Float>(-0.089, 0.343, 1.045, 1),
  SIMD4<Float>(0.000, 0.291, 1.190, 1),
  SIMD4<Float>(0.252, 0.291, 1.081, 6),
  SIMD4<Float>(0.252, 0.291, 1.190, 1),
  SIMD4<Float>(0.378, 0.364, 1.030, 6),
  SIMD4<Float>(0.631, 0.073, 1.081, 6),
  SIMD4<Float>(0.631, 0.073, 1.190, 1),
  SIMD4<Float>(0.631, -0.030, 1.045, 1),
  SIMD4<Float>(0.504, 0.146, 1.030, 6),
  SIMD4<Float>(0.757, 0.146, 1.030, 6),
  SIMD4<Float>(0.883, 0.073, 1.081, 6),
  SIMD4<Float>(0.883, 0.073, 1.190, 1),
  SIMD4<Float>(0.883, -0.030, 1.045, 1),
  SIMD4<Float>(0.504, 0.291, 1.081, 6),
  SIMD4<Float>(0.504, 0.291, 1.190, 1),
  SIMD4<Float>(0.631, 0.364, 1.030, 6),
  SIMD4<Float>(0.883, 0.364, 1.030, 6),
  SIMD4<Float>(0.757, 0.291, 1.081, 6),
  SIMD4<Float>(0.757, 0.291, 1.190, 1),
  SIMD4<Float>(0.000, 0.582, 1.030, 6),
  SIMD4<Float>(-0.089, 0.531, 1.066, 1),
  SIMD4<Float>(0.126, 0.510, 1.081, 6),
  SIMD4<Float>(0.126, 0.510, 1.190, 1),
  SIMD4<Float>(0.000, 0.728, 1.081, 6),
  SIMD4<Float>(-0.089, 0.779, 1.045, 1),
  SIMD4<Float>(0.000, 0.728, 1.190, 1),
  SIMD4<Float>(0.378, 0.510, 1.081, 6),
  SIMD4<Float>(0.378, 0.510, 1.190, 1),
  SIMD4<Float>(0.252, 0.582, 1.030, 6),
  SIMD4<Float>(0.252, 0.728, 1.081, 6),
  SIMD4<Float>(0.252, 0.728, 1.190, 1),
  SIMD4<Float>(0.126, 0.801, 1.030, 6),
  SIMD4<Float>(0.126, 0.947, 1.081, 6),
  SIMD4<Float>(0.126, 0.947, 1.190, 1),
  SIMD4<Float>(0.378, 0.801, 1.030, 6),
  SIMD4<Float>(0.378, 0.947, 1.081, 6),
  SIMD4<Float>(0.378, 0.947, 1.190, 1),
  SIMD4<Float>(0.631, 0.510, 1.081, 6),
  SIMD4<Float>(0.631, 0.510, 1.190, 1),
  SIMD4<Float>(0.504, 0.582, 1.030, 6),
  SIMD4<Float>(0.504, 0.728, 1.081, 6),
  SIMD4<Float>(0.504, 0.728, 1.190, 1),
  SIMD4<Float>(0.757, 0.582, 1.030, 6),
  SIMD4<Float>(0.883, 0.510, 1.081, 6),
  SIMD4<Float>(0.883, 0.510, 1.190, 1),
  SIMD4<Float>(0.757, 0.728, 1.081, 6),
  SIMD4<Float>(0.757, 0.728, 1.190, 1),
  SIMD4<Float>(0.631, 0.801, 1.030, 6),
  SIMD4<Float>(0.631, 0.947, 1.081, 6),
  SIMD4<Float>(0.631, 0.947, 1.190, 1),
  SIMD4<Float>(0.883, 0.801, 1.030, 6),
  SIMD4<Float>(0.883, 0.947, 1.081, 6),
  SIMD4<Float>(0.883, 0.947, 1.190, 1),
  SIMD4<Float>(1.009, 0.146, 1.030, 6),
  SIMD4<Float>(1.098, 0.094, 1.066, 1),
  SIMD4<Float>(1.009, 0.291, 1.081, 6),
  SIMD4<Float>(1.009, 0.291, 1.190, 1),
  SIMD4<Float>(1.098, 0.343, 1.045, 1),
  SIMD4<Float>(1.009, 0.582, 1.030, 6),
  SIMD4<Float>(1.098, 0.531, 1.066, 1),
  SIMD4<Float>(1.009, 0.728, 1.081, 6),
  SIMD4<Float>(1.009, 0.728, 1.190, 1),
  SIMD4<Float>(1.098, 0.779, 1.045, 1),
  SIMD4<Float>(0.000, 1.019, 1.030, 6),
  SIMD4<Float>(-0.089, 0.968, 1.066, 1),
  SIMD4<Float>(0.126, 1.238, 1.030, 6),
  SIMD4<Float>(0.126, 1.341, 1.066, 1),
  SIMD4<Float>(0.000, 1.165, 1.081, 6),
  SIMD4<Float>(-0.089, 1.216, 1.045, 1),
  SIMD4<Float>(0.000, 1.165, 1.190, 1),
  SIMD4<Float>(0.252, 1.019, 1.030, 6),
  SIMD4<Float>(0.252, 1.165, 1.081, 6),
  SIMD4<Float>(0.252, 1.165, 1.190, 1),
  SIMD4<Float>(0.378, 1.238, 1.030, 6),
  SIMD4<Float>(0.378, 1.341, 1.066, 1),
  SIMD4<Float>(0.504, 1.019, 1.030, 6),
  SIMD4<Float>(0.504, 1.165, 1.081, 6),
  SIMD4<Float>(0.504, 1.165, 1.190, 1),
  SIMD4<Float>(0.631, 1.238, 1.030, 6),
  SIMD4<Float>(0.631, 1.341, 1.066, 1),
  SIMD4<Float>(0.757, 1.019, 1.030, 6),
  SIMD4<Float>(0.883, 1.238, 1.030, 6),
  SIMD4<Float>(0.883, 1.341, 1.066, 1),
  SIMD4<Float>(0.757, 1.165, 1.081, 6),
  SIMD4<Float>(0.757, 1.165, 1.190, 1),
  SIMD4<Float>(1.009, 1.019, 1.030, 6),
  SIMD4<Float>(1.098, 0.968, 1.066, 1),
  SIMD4<Float>(1.009, 1.165, 1.081, 6),
  SIMD4<Float>(1.098, 1.216, 1.045, 1),
  SIMD4<Float>(1.009, 1.165, 1.190, 1),
]

private let lonsdaleiteBonds: [SIMD2<UInt32>] = [
  SIMD2<UInt32>(0, 1),
  SIMD2<UInt32>(0, 2),
  SIMD2<UInt32>(0, 3),
  SIMD2<UInt32>(3, 4),
  SIMD2<UInt32>(3, 5),
  SIMD2<UInt32>(5, 6),
  SIMD2<UInt32>(7, 8),
  SIMD2<UInt32>(7, 9),
  SIMD2<UInt32>(0, 10),
  SIMD2<UInt32>(7, 10),
  SIMD2<UInt32>(10, 11),
  SIMD2<UInt32>(3, 13),
  SIMD2<UInt32>(12, 13),
  SIMD2<UInt32>(13, 14),
  SIMD2<UInt32>(13, 15),
  SIMD2<UInt32>(12, 16),
  SIMD2<UInt32>(10, 17),
  SIMD2<UInt32>(12, 17),
  SIMD2<UInt32>(17, 18),
  SIMD2<UInt32>(17, 19),
  SIMD2<UInt32>(19, 20),
  SIMD2<UInt32>(5, 21),
  SIMD2<UInt32>(11, 21),
  SIMD2<UInt32>(21, 22),
  SIMD2<UInt32>(21, 23),
  SIMD2<UInt32>(23, 24),
  SIMD2<UInt32>(23, 25),
  SIMD2<UInt32>(25, 26),
  SIMD2<UInt32>(11, 27),
  SIMD2<UInt32>(27, 28),
  SIMD2<UInt32>(27, 29),
  SIMD2<UInt32>(29, 30),
  SIMD2<UInt32>(23, 31),
  SIMD2<UInt32>(29, 31),
  SIMD2<UInt32>(5, 32),
  SIMD2<UInt32>(16, 32),
  SIMD2<UInt32>(32, 33),
  SIMD2<UInt32>(25, 35),
  SIMD2<UInt32>(32, 35),
  SIMD2<UInt32>(34, 35),
  SIMD2<UInt32>(35, 36),
  SIMD2<UInt32>(11, 37),
  SIMD2<UInt32>(16, 37),
  SIMD2<UInt32>(20, 37),
  SIMD2<UInt32>(31, 38),
  SIMD2<UInt32>(34, 38),
  SIMD2<UInt32>(37, 38),
  SIMD2<UInt32>(38, 39),
  SIMD2<UInt32>(40, 41),
  SIMD2<UInt32>(40, 42),
  SIMD2<UInt32>(7, 43),
  SIMD2<UInt32>(40, 43),
  SIMD2<UInt32>(43, 44),
  SIMD2<UInt32>(27, 44),
  SIMD2<UInt32>(45, 46),
  SIMD2<UInt32>(45, 47),
  SIMD2<UInt32>(40, 48),
  SIMD2<UInt32>(45, 48),
  SIMD2<UInt32>(48, 49),
  SIMD2<UInt32>(43, 50),
  SIMD2<UInt32>(19, 50),
  SIMD2<UInt32>(50, 51),
  SIMD2<UInt32>(50, 52),
  SIMD2<UInt32>(52, 53),
  SIMD2<UInt32>(48, 55),
  SIMD2<UInt32>(52, 55),
  SIMD2<UInt32>(54, 55),
  SIMD2<UInt32>(55, 56),
  SIMD2<UInt32>(54, 57),
  SIMD2<UInt32>(44, 58),
  SIMD2<UInt32>(49, 58),
  SIMD2<UInt32>(58, 59),
  SIMD2<UInt32>(58, 60),
  SIMD2<UInt32>(60, 61),
  SIMD2<UInt32>(29, 62),
  SIMD2<UInt32>(60, 62),
  SIMD2<UInt32>(49, 63),
  SIMD2<UInt32>(63, 64),
  SIMD2<UInt32>(63, 65),
  SIMD2<UInt32>(65, 66),
  SIMD2<UInt32>(60, 67),
  SIMD2<UInt32>(65, 67),
  SIMD2<UInt32>(44, 68),
  SIMD2<UInt32>(20, 68),
  SIMD2<UInt32>(53, 68),
  SIMD2<UInt32>(62, 69),
  SIMD2<UInt32>(39, 69),
  SIMD2<UInt32>(68, 69),
  SIMD2<UInt32>(69, 70),
  SIMD2<UInt32>(49, 71),
  SIMD2<UInt32>(53, 71),
  SIMD2<UInt32>(57, 71),
  SIMD2<UInt32>(67, 73),
  SIMD2<UInt32>(70, 73),
  SIMD2<UInt32>(71, 73),
  SIMD2<UInt32>(72, 73),
  SIMD2<UInt32>(12, 74),
  SIMD2<UInt32>(74, 75),
  SIMD2<UInt32>(76, 77),
  SIMD2<UInt32>(76, 78),
  SIMD2<UInt32>(74, 79),
  SIMD2<UInt32>(76, 79),
  SIMD2<UInt32>(79, 80),
  SIMD2<UInt32>(79, 81),
  SIMD2<UInt32>(81, 82),
  SIMD2<UInt32>(19, 83),
  SIMD2<UInt32>(83, 84),
  SIMD2<UInt32>(74, 85),
  SIMD2<UInt32>(83, 85),
  SIMD2<UInt32>(85, 86),
  SIMD2<UInt32>(86, 87),
  SIMD2<UInt32>(85, 88),
  SIMD2<UInt32>(76, 89),
  SIMD2<UInt32>(86, 89),
  SIMD2<UInt32>(89, 90),
  SIMD2<UInt32>(89, 91),
  SIMD2<UInt32>(91, 92),
  SIMD2<UInt32>(93, 94),
  SIMD2<UInt32>(86, 95),
  SIMD2<UInt32>(93, 95),
  SIMD2<UInt32>(95, 96),
  SIMD2<UInt32>(16, 97),
  SIMD2<UInt32>(81, 97),
  SIMD2<UInt32>(88, 97),
  SIMD2<UInt32>(81, 98),
  SIMD2<UInt32>(90, 98),
  SIMD2<UInt32>(98, 99),
  SIMD2<UInt32>(34, 100),
  SIMD2<UInt32>(97, 100),
  SIMD2<UInt32>(98, 101),
  SIMD2<UInt32>(101, 102),
  SIMD2<UInt32>(100, 103),
  SIMD2<UInt32>(101, 103),
  SIMD2<UInt32>(103, 104),
  SIMD2<UInt32>(20, 105),
  SIMD2<UInt32>(88, 105),
  SIMD2<UInt32>(88, 106),
  SIMD2<UInt32>(90, 106),
  SIMD2<UInt32>(96, 106),
  SIMD2<UInt32>(39, 107),
  SIMD2<UInt32>(105, 107),
  SIMD2<UInt32>(100, 108),
  SIMD2<UInt32>(107, 108),
  SIMD2<UInt32>(106, 109),
  SIMD2<UInt32>(108, 109),
  SIMD2<UInt32>(90, 110),
  SIMD2<UInt32>(101, 111),
  SIMD2<UInt32>(109, 111),
  SIMD2<UInt32>(110, 112),
  SIMD2<UInt32>(111, 112),
  SIMD2<UInt32>(96, 113),
  SIMD2<UInt32>(113, 114),
  SIMD2<UInt32>(109, 115),
  SIMD2<UInt32>(114, 115),
  SIMD2<UInt32>(52, 116),
  SIMD2<UInt32>(116, 117),
  SIMD2<UInt32>(83, 118),
  SIMD2<UInt32>(116, 118),
  SIMD2<UInt32>(118, 119),
  SIMD2<UInt32>(105, 119),
  SIMD2<UInt32>(118, 120),
  SIMD2<UInt32>(95, 120),
  SIMD2<UInt32>(120, 121),
  SIMD2<UInt32>(54, 122),
  SIMD2<UInt32>(122, 123),
  SIMD2<UInt32>(124, 125),
  SIMD2<UInt32>(116, 126),
  SIMD2<UInt32>(122, 126),
  SIMD2<UInt32>(124, 126),
  SIMD2<UInt32>(126, 127),
  SIMD2<UInt32>(120, 128),
  SIMD2<UInt32>(124, 128),
  SIMD2<UInt32>(128, 129),
  SIMD2<UInt32>(129, 130),
  SIMD2<UInt32>(128, 131),
  SIMD2<UInt32>(124, 132),
  SIMD2<UInt32>(132, 133),
  SIMD2<UInt32>(132, 134),
  SIMD2<UInt32>(134, 135),
  SIMD2<UInt32>(53, 136),
  SIMD2<UInt32>(119, 136),
  SIMD2<UInt32>(127, 136),
  SIMD2<UInt32>(119, 137),
  SIMD2<UInt32>(96, 137),
  SIMD2<UInt32>(131, 137),
  SIMD2<UInt32>(70, 138),
  SIMD2<UInt32>(136, 138),
  SIMD2<UInt32>(107, 139),
  SIMD2<UInt32>(138, 139),
  SIMD2<UInt32>(137, 140),
  SIMD2<UInt32>(139, 140),
  SIMD2<UInt32>(115, 140),
  SIMD2<UInt32>(57, 141),
  SIMD2<UInt32>(127, 141),
  SIMD2<UInt32>(127, 142),
  SIMD2<UInt32>(131, 142),
  SIMD2<UInt32>(133, 142),
  SIMD2<UInt32>(72, 143),
  SIMD2<UInt32>(141, 143),
  SIMD2<UInt32>(142, 144),
  SIMD2<UInt32>(138, 145),
  SIMD2<UInt32>(143, 145),
  SIMD2<UInt32>(144, 145),
  SIMD2<UInt32>(131, 146),
  SIMD2<UInt32>(140, 147),
  SIMD2<UInt32>(144, 147),
  SIMD2<UInt32>(146, 148),
  SIMD2<UInt32>(147, 148),
  SIMD2<UInt32>(133, 149),
  SIMD2<UInt32>(144, 150),
  SIMD2<UInt32>(149, 151),
  SIMD2<UInt32>(150, 151),
  SIMD2<UInt32>(25, 152),
  SIMD2<UInt32>(152, 153),
  SIMD2<UInt32>(152, 154),
  SIMD2<UInt32>(154, 155),
  SIMD2<UInt32>(156, 157),
  SIMD2<UInt32>(31, 158),
  SIMD2<UInt32>(154, 158),
  SIMD2<UInt32>(156, 158),
  SIMD2<UInt32>(34, 159),
  SIMD2<UInt32>(152, 160),
  SIMD2<UInt32>(159, 160),
  SIMD2<UInt32>(160, 161),
  SIMD2<UInt32>(158, 162),
  SIMD2<UInt32>(159, 162),
  SIMD2<UInt32>(39, 163),
  SIMD2<UInt32>(162, 163),
  SIMD2<UInt32>(154, 164),
  SIMD2<UInt32>(164, 165),
  SIMD2<UInt32>(164, 166),
  SIMD2<UInt32>(166, 167),
  SIMD2<UInt32>(156, 168),
  SIMD2<UInt32>(168, 169),
  SIMD2<UInt32>(164, 170),
  SIMD2<UInt32>(168, 170),
  SIMD2<UInt32>(160, 172),
  SIMD2<UInt32>(166, 172),
  SIMD2<UInt32>(171, 172),
  SIMD2<UInt32>(172, 173),
  SIMD2<UInt32>(162, 174),
  SIMD2<UInt32>(170, 174),
  SIMD2<UInt32>(171, 174),
  SIMD2<UInt32>(174, 175),
  SIMD2<UInt32>(176, 177),
  SIMD2<UInt32>(62, 178),
  SIMD2<UInt32>(156, 178),
  SIMD2<UInt32>(176, 178),
  SIMD2<UInt32>(67, 179),
  SIMD2<UInt32>(176, 179),
  SIMD2<UInt32>(179, 180),
  SIMD2<UInt32>(180, 181),
  SIMD2<UInt32>(178, 182),
  SIMD2<UInt32>(163, 182),
  SIMD2<UInt32>(70, 183),
  SIMD2<UInt32>(182, 183),
  SIMD2<UInt32>(72, 184),
  SIMD2<UInt32>(179, 185),
  SIMD2<UInt32>(183, 185),
  SIMD2<UInt32>(184, 185),
  SIMD2<UInt32>(176, 186),
  SIMD2<UInt32>(186, 187),
  SIMD2<UInt32>(168, 188),
  SIMD2<UInt32>(186, 188),
  SIMD2<UInt32>(180, 189),
  SIMD2<UInt32>(189, 190),
  SIMD2<UInt32>(186, 191),
  SIMD2<UInt32>(189, 191),
  SIMD2<UInt32>(182, 192),
  SIMD2<UInt32>(188, 192),
  SIMD2<UInt32>(175, 192),
  SIMD2<UInt32>(192, 193),
  SIMD2<UInt32>(185, 195),
  SIMD2<UInt32>(191, 195),
  SIMD2<UInt32>(193, 195),
  SIMD2<UInt32>(194, 195),
  SIMD2<UInt32>(103, 196),
  SIMD2<UInt32>(196, 197),
  SIMD2<UInt32>(159, 198),
  SIMD2<UInt32>(196, 198),
  SIMD2<UInt32>(196, 199),
  SIMD2<UInt32>(199, 200),
  SIMD2<UInt32>(163, 201),
  SIMD2<UInt32>(108, 202),
  SIMD2<UInt32>(198, 202),
  SIMD2<UInt32>(201, 202),
  SIMD2<UInt32>(202, 203),
  SIMD2<UInt32>(111, 204),
  SIMD2<UInt32>(199, 204),
  SIMD2<UInt32>(203, 204),
  SIMD2<UInt32>(204, 205),
  SIMD2<UInt32>(115, 206),
  SIMD2<UInt32>(203, 206),
  SIMD2<UInt32>(206, 207),
  SIMD2<UInt32>(198, 208),
  SIMD2<UInt32>(171, 208),
  SIMD2<UInt32>(199, 209),
  SIMD2<UInt32>(209, 210),
  SIMD2<UInt32>(208, 211),
  SIMD2<UInt32>(209, 211),
  SIMD2<UInt32>(211, 212),
  SIMD2<UInt32>(201, 213),
  SIMD2<UInt32>(175, 213),
  SIMD2<UInt32>(208, 214),
  SIMD2<UInt32>(213, 214),
  SIMD2<UInt32>(203, 215),
  SIMD2<UInt32>(214, 215),
  SIMD2<UInt32>(209, 216),
  SIMD2<UInt32>(215, 216),
  SIMD2<UInt32>(205, 217),
  SIMD2<UInt32>(216, 217),
  SIMD2<UInt32>(207, 218),
  SIMD2<UInt32>(215, 219),
  SIMD2<UInt32>(218, 219),
  SIMD2<UInt32>(183, 220),
  SIMD2<UInt32>(139, 221),
  SIMD2<UInt32>(201, 221),
  SIMD2<UInt32>(220, 221),
  SIMD2<UInt32>(221, 222),
  SIMD2<UInt32>(206, 222),
  SIMD2<UInt32>(145, 223),
  SIMD2<UInt32>(220, 223),
  SIMD2<UInt32>(184, 224),
  SIMD2<UInt32>(223, 224),
  SIMD2<UInt32>(223, 225),
  SIMD2<UInt32>(147, 226),
  SIMD2<UInt32>(222, 226),
  SIMD2<UInt32>(225, 226),
  SIMD2<UInt32>(226, 227),
  SIMD2<UInt32>(150, 228),
  SIMD2<UInt32>(225, 228),
  SIMD2<UInt32>(228, 229),
  SIMD2<UInt32>(220, 230),
  SIMD2<UInt32>(193, 230),
  SIMD2<UInt32>(213, 231),
  SIMD2<UInt32>(230, 231),
  SIMD2<UInt32>(222, 232),
  SIMD2<UInt32>(231, 232),
  SIMD2<UInt32>(219, 232),
  SIMD2<UInt32>(224, 233),
  SIMD2<UInt32>(194, 233),
  SIMD2<UInt32>(225, 234),
  SIMD2<UInt32>(230, 235),
  SIMD2<UInt32>(233, 235),
  SIMD2<UInt32>(234, 235),
  SIMD2<UInt32>(232, 236),
  SIMD2<UInt32>(234, 236),
  SIMD2<UInt32>(227, 237),
  SIMD2<UInt32>(236, 237),
  SIMD2<UInt32>(234, 238),
  SIMD2<UInt32>(229, 239),
  SIMD2<UInt32>(238, 239),
  SIMD2<UInt32>(45, 240),
  SIMD2<UInt32>(240, 241),
  SIMD2<UInt32>(240, 242),
  SIMD2<UInt32>(63, 242),
  SIMD2<UInt32>(242, 243),
  SIMD2<UInt32>(240, 244),
  SIMD2<UInt32>(54, 244),
  SIMD2<UInt32>(244, 245),
  SIMD2<UInt32>(244, 246),
  SIMD2<UInt32>(65, 247),
  SIMD2<UInt32>(247, 248),
  SIMD2<UInt32>(242, 249),
  SIMD2<UInt32>(57, 249),
  SIMD2<UInt32>(249, 250),
  SIMD2<UInt32>(247, 251),
  SIMD2<UInt32>(72, 251),
  SIMD2<UInt32>(249, 251),
  SIMD2<UInt32>(251, 252),
  SIMD2<UInt32>(122, 253),
  SIMD2<UInt32>(253, 254),
  SIMD2<UInt32>(253, 255),
  SIMD2<UInt32>(132, 255),
  SIMD2<UInt32>(255, 256),
  SIMD2<UInt32>(255, 257),
  SIMD2<UInt32>(253, 258),
  SIMD2<UInt32>(141, 258),
  SIMD2<UInt32>(258, 259),
  SIMD2<UInt32>(258, 260),
  SIMD2<UInt32>(133, 260),
  SIMD2<UInt32>(260, 261),
  SIMD2<UInt32>(143, 262),
  SIMD2<UInt32>(262, 263),
  SIMD2<UInt32>(260, 264),
  SIMD2<UInt32>(262, 264),
  SIMD2<UInt32>(150, 264),
  SIMD2<UInt32>(264, 265),
  SIMD2<UInt32>(247, 266),
  SIMD2<UInt32>(180, 266),
  SIMD2<UInt32>(266, 267),
  SIMD2<UInt32>(266, 268),
  SIMD2<UInt32>(184, 268),
  SIMD2<UInt32>(268, 269),
  SIMD2<UInt32>(189, 270),
  SIMD2<UInt32>(270, 271),
  SIMD2<UInt32>(268, 272),
  SIMD2<UInt32>(270, 272),
  SIMD2<UInt32>(194, 272),
  SIMD2<UInt32>(272, 273),
  SIMD2<UInt32>(262, 274),
  SIMD2<UInt32>(224, 274),
  SIMD2<UInt32>(274, 275),
  SIMD2<UInt32>(274, 276),
  SIMD2<UInt32>(228, 276),
  SIMD2<UInt32>(276, 277),
  SIMD2<UInt32>(233, 278),
  SIMD2<UInt32>(278, 279),
  SIMD2<UInt32>(276, 280),
  SIMD2<UInt32>(278, 280),
  SIMD2<UInt32>(238, 280),
  SIMD2<UInt32>(280, 281),
  SIMD2<UInt32>(282, 283),
  SIMD2<UInt32>(282, 284),
  SIMD2<UInt32>(284, 285),
  SIMD2<UInt32>(284, 286),
  SIMD2<UInt32>(91, 287),
  SIMD2<UInt32>(284, 287),
  SIMD2<UInt32>(287, 288),
  SIMD2<UInt32>(287, 289),
  SIMD2<UInt32>(110, 289),
  SIMD2<UInt32>(289, 290),
  SIMD2<UInt32>(282, 291),
  SIMD2<UInt32>(291, 292),
  SIMD2<UInt32>(91, 293),
  SIMD2<UInt32>(93, 293),
  SIMD2<UInt32>(282, 294),
  SIMD2<UInt32>(293, 294),
  SIMD2<UInt32>(294, 295),
  SIMD2<UInt32>(293, 296),
  SIMD2<UInt32>(110, 296),
  SIMD2<UInt32>(113, 296),
  SIMD2<UInt32>(294, 297),
  SIMD2<UInt32>(297, 298),
  SIMD2<UInt32>(297, 299),
  SIMD2<UInt32>(299, 300),
  SIMD2<UInt32>(289, 301),
  SIMD2<UInt32>(291, 301),
  SIMD2<UInt32>(301, 302),
  SIMD2<UInt32>(303, 304),
  SIMD2<UInt32>(301, 305),
  SIMD2<UInt32>(303, 305),
  SIMD2<UInt32>(305, 306),
  SIMD2<UInt32>(112, 307),
  SIMD2<UInt32>(305, 307),
  SIMD2<UInt32>(307, 308),
  SIMD2<UInt32>(291, 309),
  SIMD2<UInt32>(296, 309),
  SIMD2<UInt32>(299, 309),
  SIMD2<UInt32>(112, 310),
  SIMD2<UInt32>(114, 310),
  SIMD2<UInt32>(303, 311),
  SIMD2<UInt32>(309, 311),
  SIMD2<UInt32>(310, 311),
  SIMD2<UInt32>(311, 312),
  SIMD2<UInt32>(312, 313),
  SIMD2<UInt32>(93, 314),
  SIMD2<UInt32>(129, 314),
  SIMD2<UInt32>(314, 315),
  SIMD2<UInt32>(113, 315),
  SIMD2<UInt32>(146, 315),
  SIMD2<UInt32>(297, 316),
  SIMD2<UInt32>(314, 316),
  SIMD2<UInt32>(316, 317),
  SIMD2<UInt32>(316, 318),
  SIMD2<UInt32>(318, 319),
  SIMD2<UInt32>(318, 320),
  SIMD2<UInt32>(320, 321),
  SIMD2<UInt32>(322, 323),
  SIMD2<UInt32>(318, 324),
  SIMD2<UInt32>(322, 324),
  SIMD2<UInt32>(324, 325),
  SIMD2<UInt32>(129, 326),
  SIMD2<UInt32>(134, 326),
  SIMD2<UInt32>(324, 326),
  SIMD2<UInt32>(326, 327),
  SIMD2<UInt32>(146, 327),
  SIMD2<UInt32>(149, 327),
  SIMD2<UInt32>(322, 328),
  SIMD2<UInt32>(328, 329),
  SIMD2<UInt32>(299, 330),
  SIMD2<UInt32>(315, 330),
  SIMD2<UInt32>(320, 330),
  SIMD2<UInt32>(114, 331),
  SIMD2<UInt32>(148, 331),
  SIMD2<UInt32>(312, 332),
  SIMD2<UInt32>(330, 332),
  SIMD2<UInt32>(331, 332),
  SIMD2<UInt32>(332, 333),
  SIMD2<UInt32>(333, 334),
  SIMD2<UInt32>(320, 335),
  SIMD2<UInt32>(327, 335),
  SIMD2<UInt32>(328, 335),
  SIMD2<UInt32>(336, 337),
  SIMD2<UInt32>(333, 338),
  SIMD2<UInt32>(335, 338),
  SIMD2<UInt32>(336, 338),
  SIMD2<UInt32>(148, 339),
  SIMD2<UInt32>(151, 339),
  SIMD2<UInt32>(338, 339),
  SIMD2<UInt32>(307, 340),
  SIMD2<UInt32>(205, 340),
  SIMD2<UInt32>(340, 341),
  SIMD2<UInt32>(303, 342),
  SIMD2<UInt32>(342, 343),
  SIMD2<UInt32>(340, 344),
  SIMD2<UInt32>(342, 344),
  SIMD2<UInt32>(344, 345),
  SIMD2<UInt32>(310, 346),
  SIMD2<UInt32>(205, 346),
  SIMD2<UInt32>(207, 346),
  SIMD2<UInt32>(342, 347),
  SIMD2<UInt32>(346, 347),
  SIMD2<UInt32>(312, 348),
  SIMD2<UInt32>(347, 348),
  SIMD2<UInt32>(348, 349),
  SIMD2<UInt32>(350, 351),
  SIMD2<UInt32>(344, 352),
  SIMD2<UInt32>(350, 352),
  SIMD2<UInt32>(352, 353),
  SIMD2<UInt32>(217, 354),
  SIMD2<UInt32>(352, 354),
  SIMD2<UInt32>(354, 355),
  SIMD2<UInt32>(217, 356),
  SIMD2<UInt32>(218, 356),
  SIMD2<UInt32>(347, 357),
  SIMD2<UInt32>(350, 357),
  SIMD2<UInt32>(356, 357),
  SIMD2<UInt32>(357, 358),
  SIMD2<UInt32>(358, 359),
  SIMD2<UInt32>(331, 360),
  SIMD2<UInt32>(207, 360),
  SIMD2<UInt32>(227, 360),
  SIMD2<UInt32>(348, 361),
  SIMD2<UInt32>(360, 361),
  SIMD2<UInt32>(333, 362),
  SIMD2<UInt32>(361, 362),
  SIMD2<UInt32>(362, 363),
  SIMD2<UInt32>(339, 364),
  SIMD2<UInt32>(227, 364),
  SIMD2<UInt32>(229, 364),
  SIMD2<UInt32>(336, 365),
  SIMD2<UInt32>(365, 366),
  SIMD2<UInt32>(362, 367),
  SIMD2<UInt32>(364, 367),
  SIMD2<UInt32>(365, 367),
  SIMD2<UInt32>(218, 368),
  SIMD2<UInt32>(237, 368),
  SIMD2<UInt32>(361, 369),
  SIMD2<UInt32>(358, 369),
  SIMD2<UInt32>(368, 369),
  SIMD2<UInt32>(369, 370),
  SIMD2<UInt32>(370, 371),
  SIMD2<UInt32>(372, 373),
  SIMD2<UInt32>(367, 374),
  SIMD2<UInt32>(370, 374),
  SIMD2<UInt32>(372, 374),
  SIMD2<UInt32>(237, 375),
  SIMD2<UInt32>(239, 375),
  SIMD2<UInt32>(374, 375),
  SIMD2<UInt32>(134, 376),
  SIMD2<UInt32>(376, 377),
  SIMD2<UInt32>(322, 378),
  SIMD2<UInt32>(376, 378),
  SIMD2<UInt32>(378, 379),
  SIMD2<UInt32>(378, 380),
  SIMD2<UInt32>(376, 381),
  SIMD2<UInt32>(149, 381),
  SIMD2<UInt32>(381, 382),
  SIMD2<UInt32>(328, 383),
  SIMD2<UInt32>(381, 383),
  SIMD2<UInt32>(383, 384),
  SIMD2<UInt32>(151, 385),
  SIMD2<UInt32>(385, 386),
  SIMD2<UInt32>(336, 387),
  SIMD2<UInt32>(383, 387),
  SIMD2<UInt32>(385, 387),
  SIMD2<UInt32>(387, 388),
  SIMD2<UInt32>(385, 389),
  SIMD2<UInt32>(229, 389),
  SIMD2<UInt32>(389, 390),
  SIMD2<UInt32>(365, 391),
  SIMD2<UInt32>(389, 391),
  SIMD2<UInt32>(391, 392),
  SIMD2<UInt32>(239, 393),
  SIMD2<UInt32>(393, 394),
  SIMD2<UInt32>(391, 395),
  SIMD2<UInt32>(372, 395),
  SIMD2<UInt32>(393, 395),
  SIMD2<UInt32>(395, 396),
  SIMD2<UInt32>(166, 397),
  SIMD2<UInt32>(397, 398),
  SIMD2<UInt32>(397, 399),
  SIMD2<UInt32>(399, 400),
  SIMD2<UInt32>(399, 401),
  SIMD2<UInt32>(402, 403),
  SIMD2<UInt32>(402, 404),
  SIMD2<UInt32>(170, 405),
  SIMD2<UInt32>(399, 405),
  SIMD2<UInt32>(402, 405),
  SIMD2<UInt32>(171, 406),
  SIMD2<UInt32>(397, 407),
  SIMD2<UInt32>(406, 407),
  SIMD2<UInt32>(407, 408),
  SIMD2<UInt32>(407, 409),
  SIMD2<UInt32>(405, 410),
  SIMD2<UInt32>(406, 410),
  SIMD2<UInt32>(410, 411),
  SIMD2<UInt32>(175, 412),
  SIMD2<UInt32>(410, 412),
  SIMD2<UInt32>(413, 414),
  SIMD2<UInt32>(413, 415),
  SIMD2<UInt32>(188, 416),
  SIMD2<UInt32>(402, 416),
  SIMD2<UInt32>(413, 416),
  SIMD2<UInt32>(191, 417),
  SIMD2<UInt32>(413, 417),
  SIMD2<UInt32>(417, 418),
  SIMD2<UInt32>(418, 419),
  SIMD2<UInt32>(418, 420),
  SIMD2<UInt32>(416, 421),
  SIMD2<UInt32>(412, 421),
  SIMD2<UInt32>(421, 422),
  SIMD2<UInt32>(193, 423),
  SIMD2<UInt32>(421, 423),
  SIMD2<UInt32>(194, 424),
  SIMD2<UInt32>(417, 425),
  SIMD2<UInt32>(423, 425),
  SIMD2<UInt32>(424, 425),
  SIMD2<UInt32>(425, 426),
  SIMD2<UInt32>(211, 427),
  SIMD2<UInt32>(427, 428),
  SIMD2<UInt32>(406, 429),
  SIMD2<UInt32>(427, 429),
  SIMD2<UInt32>(429, 430),
  SIMD2<UInt32>(427, 431),
  SIMD2<UInt32>(431, 432),
  SIMD2<UInt32>(431, 433),
  SIMD2<UInt32>(412, 434),
  SIMD2<UInt32>(434, 435),
  SIMD2<UInt32>(214, 436),
  SIMD2<UInt32>(429, 436),
  SIMD2<UInt32>(434, 436),
  SIMD2<UInt32>(436, 437),
  SIMD2<UInt32>(437, 438),
  SIMD2<UInt32>(216, 439),
  SIMD2<UInt32>(431, 439),
  SIMD2<UInt32>(437, 439),
  SIMD2<UInt32>(439, 440),
  SIMD2<UInt32>(440, 441),
  SIMD2<UInt32>(219, 442),
  SIMD2<UInt32>(437, 442),
  SIMD2<UInt32>(442, 443),
  SIMD2<UInt32>(443, 444),
  SIMD2<UInt32>(423, 445),
  SIMD2<UInt32>(445, 446),
  SIMD2<UInt32>(231, 447),
  SIMD2<UInt32>(434, 447),
  SIMD2<UInt32>(445, 447),
  SIMD2<UInt32>(447, 448),
  SIMD2<UInt32>(442, 448),
  SIMD2<UInt32>(448, 449),
  SIMD2<UInt32>(235, 450),
  SIMD2<UInt32>(445, 450),
  SIMD2<UInt32>(424, 451),
  SIMD2<UInt32>(450, 451),
  SIMD2<UInt32>(451, 452),
  SIMD2<UInt32>(450, 453),
  SIMD2<UInt32>(453, 454),
  SIMD2<UInt32>(236, 455),
  SIMD2<UInt32>(448, 455),
  SIMD2<UInt32>(453, 455),
  SIMD2<UInt32>(455, 456),
  SIMD2<UInt32>(456, 457),
  SIMD2<UInt32>(238, 458),
  SIMD2<UInt32>(453, 458),
  SIMD2<UInt32>(458, 459),
  SIMD2<UInt32>(459, 460),
  SIMD2<UInt32>(270, 461),
  SIMD2<UInt32>(418, 461),
  SIMD2<UInt32>(461, 462),
  SIMD2<UInt32>(461, 463),
  SIMD2<UInt32>(424, 463),
  SIMD2<UInt32>(463, 464),
  SIMD2<UInt32>(463, 465),
  SIMD2<UInt32>(278, 466),
  SIMD2<UInt32>(451, 466),
  SIMD2<UInt32>(466, 467),
  SIMD2<UInt32>(466, 468),
  SIMD2<UInt32>(458, 468),
  SIMD2<UInt32>(468, 469),
  SIMD2<UInt32>(468, 470),
  SIMD2<UInt32>(354, 471),
  SIMD2<UInt32>(440, 471),
  SIMD2<UInt32>(471, 472),
  SIMD2<UInt32>(350, 473),
  SIMD2<UInt32>(473, 474),
  SIMD2<UInt32>(471, 475),
  SIMD2<UInt32>(473, 475),
  SIMD2<UInt32>(475, 476),
  SIMD2<UInt32>(475, 477),
  SIMD2<UInt32>(356, 478),
  SIMD2<UInt32>(440, 478),
  SIMD2<UInt32>(443, 478),
  SIMD2<UInt32>(473, 479),
  SIMD2<UInt32>(478, 479),
  SIMD2<UInt32>(479, 480),
  SIMD2<UInt32>(358, 481),
  SIMD2<UInt32>(479, 481),
  SIMD2<UInt32>(481, 482),
  SIMD2<UInt32>(368, 483),
  SIMD2<UInt32>(443, 483),
  SIMD2<UInt32>(456, 483),
  SIMD2<UInt32>(481, 484),
  SIMD2<UInt32>(483, 484),
  SIMD2<UInt32>(484, 485),
  SIMD2<UInt32>(370, 486),
  SIMD2<UInt32>(484, 486),
  SIMD2<UInt32>(486, 487),
  SIMD2<UInt32>(375, 488),
  SIMD2<UInt32>(456, 488),
  SIMD2<UInt32>(459, 488),
  SIMD2<UInt32>(372, 489),
  SIMD2<UInt32>(489, 490),
  SIMD2<UInt32>(486, 491),
  SIMD2<UInt32>(488, 491),
  SIMD2<UInt32>(489, 491),
  SIMD2<UInt32>(491, 492),
  SIMD2<UInt32>(393, 493),
  SIMD2<UInt32>(459, 493),
  SIMD2<UInt32>(493, 494),
  SIMD2<UInt32>(489, 495),
  SIMD2<UInt32>(493, 495),
  SIMD2<UInt32>(495, 496),
  SIMD2<UInt32>(495, 497),
]
