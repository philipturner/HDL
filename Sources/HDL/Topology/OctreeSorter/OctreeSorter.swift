//
//  OctreeSorter.swift
//  HDL
//
//  Created by Philip Turner on 7/14/25.
//

struct OctreeSorter {
  var atoms: [Atom] = []
  
  var origin: SIMD3<Float>
  var dimensions: SIMD3<Float>
  
  init(atoms: [Atom]) {
    self.atoms = atoms
    
    if atoms.count == 0 {
      origin = .zero
      dimensions = SIMD3(repeating: 1)
    } else {
      var minimum = SIMD3<Float>(repeating: .greatestFiniteMagnitude)
      var maximum = -minimum
      for atom in atoms {
        // @_transparent attribute is ineffective.
        let position = unsafeBitCast(atom, to: SIMD3<Float>.self)
        minimum.replace(with: position, where: position .< minimum)
        maximum.replace(with: position, where: position .> maximum)
      }
      minimum.round(.down)
      maximum.round(.up)
      
      origin = minimum
      dimensions = maximum - minimum
      dimensions.replace(
        with: SIMD3(repeating: 1),
        where: dimensions .< 1)
    }
  }
  
  static func invertOrder(_ input: [UInt32]) -> [UInt32] {
    return [UInt32](unsafeUninitializedCapacity: input.count) {
      $1 = input.count
      let baseAddress = $0.baseAddress.unsafelyUnwrapped
      
      for reorderedID in input.indices {
        let originalID = Int(input[reorderedID])
        baseAddress[originalID] = UInt32(reorderedID)
      }
    }
  }
  
  var highestLevelSize: Float {
    func truncatedDimensions() -> SIMD3<Float> {
      SIMD3<Float>(SIMD3<Int>(dimensions))
    }
    guard all(dimensions .== truncatedDimensions()) else {
      fatalError("Dimensions must be integers.")
    }
    
    var output: Float = 1
    for _ in 0..<100 {
      if output < dimensions.max() {
        output = 2 * output
      } else {
        return output
      }
    }
    
    fatalError("Took too many iterations to find highest level size.")
  }
}
