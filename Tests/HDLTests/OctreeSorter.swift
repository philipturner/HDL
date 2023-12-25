//
//  OctreeSorter.swift
//
//
//  Created by Philip Turner on 12/25/23.
//

import HDL

struct OctreeSorter {
  var atoms: [Entity] = []
  
  var origin: SIMD3<Float>
  var dimensions: SIMD3<Float>
  
  init(atoms: [Entity]) {
    self.atoms = atoms
    
    if atoms.count == 0 {
      origin = .zero
      dimensions = .init(repeating: 1)
    } else {
      var minimum = SIMD3<Float>(repeating: .greatestFiniteMagnitude)
      var maximum = -minimum
      for atom in atoms {
        let position = atom.position
        minimum.replace(with: position, where: position .< minimum)
        maximum.replace(with: position, where: position .> maximum)
      }
      for atom in atoms {
        let position = atom.position
        minimum.replace(with: position, where: position .< minimum)
      }
      
      origin = minimum
      dimensions = (maximum - minimum).rounded(.up)
      dimensions.replace(with: .init(repeating: 1), where: dimensions .<= 0)
      
      for lane in 0..<3 {
        let element = dimensions[lane]
        if element != element.binade {
          dimensions[lane] = 2 * element.binade
        }
      }
    }
  }
  
  func mortonReordering() -> [UInt32] {
    var output: [UInt32] = []
    
    func traverse(
      atomIDs: [UInt32],
      levelOrigin: SIMD3<Float>,
      levelSize: Float
    ) {
      // Optimization: perform this check before, not after, the function call.
      if levelSize <= 1 / 32 || atomIDs.count == 1 {
        output += atomIDs
        return
      }
      
      // Optimization: create an array of arrays, indexed by the key. The
      // underlying memory allocation is recycled each function call.
      var dictionary: [SIMD3<Int32>: [UInt32]] = [:]
      
      for atomID32 in atomIDs {
        var index: SIMD3<Int32> = .init(repeating: 1)
        let atomID = Int(truncatingIfNeeded: atomID32)
        let atomPosition = atoms[atomID].position - self.origin
        index.replace(
          with: .init(repeating: -1),
          where: atomPosition .< levelOrigin)
        
        // Optimization: modify the array in-place instead of copying it upon
        // every access.
        var list: [UInt32] = dictionary[index] ?? []
        list.append(UInt32(truncatingIfNeeded: atomID))
        dictionary[index] = list
      }
      
      // Optimization: remove the key sort operation, when the dictionary is
      // changed to a set of 8 arrays.
      let sortedKeys = dictionary.keys.sorted {
        if $0.z != $1.z {
          return $0.z < $1.z
        }
        if $0.y != $1.y {
          return $0.y < $1.y
        }
        if $0.x != $1.x {
          return $0.x < $1.x
        }
        return true
      }
      
      for key in sortedKeys {
        let newOrigin = levelOrigin + SIMD3<Float>(key) * levelSize / 2
        let values = dictionary[key]!
        traverse(
          atomIDs: values, levelOrigin: newOrigin, levelSize: levelSize / 2)
      }
    }
    
    let levelSize = self.dimensions.max() / 2
    let levelOrigin = SIMD3<Float>(repeating: levelSize)
    let initialArray = atoms.indices.map(UInt32.init(truncatingIfNeeded:))
    traverse(
      atomIDs: initialArray, levelOrigin: levelOrigin, levelSize: levelSize)
    precondition(output.count == atoms.count)
    
    var reordering = [UInt32](repeating: .max, count: atoms.count)
    for reorderedID in output.indices {
      let originalID32 = output[reorderedID]
      let originalID = Int(truncatingIfNeeded: originalID32)
      let reorderedID32 = UInt32(truncatingIfNeeded: reorderedID)
      reordering[originalID] = reorderedID32
    }
    precondition(!reordering.contains(.max))
    return reordering
  }
}
