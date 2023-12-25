//
//  OctreeSort.swift
//
//
//  Created by Philip Turner on 12/25/23.
//

/*
 func traverseOctree(atomIDs: [Int], origin: SIMD3<Float>, levelSize: Float, depth: Int) {
//    let prefix = String(repeating: "-", count: depth)
//    print(prefix, origin, atomIDs.count)
   if levelSize <= 1 / 32 || atomIDs.count <= 1 {
     octreeReordering += atomIDs
     return
   }
   
   // TODO: Speed up octree traversal after you know it makes correct outputs.
   // Profile the octree against the existing sorter for large crystals.
   var dictionary: [SIMD3<Int32>: [Int]] = [:]
   
   for atomID in atomIDs {
     var index: SIMD3<Int32> = .init(repeating: 1)
     let atomPosition = topology.atoms[atomID].position + SIMD3(1, 1, 1)
     index.replace(with: .init(repeating: -1), where: atomPosition .< origin)
     
     var list: [Int] = dictionary[index] ?? []
     list.append(atomID)
     dictionary[index] = list
   }
//    precondition(
//      dictionary.keys.count >= 1 && dictionary.keys.count <= 8)
   
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
     let newOrigin = origin + SIMD3<Float>(key) * levelSize / 2
     let values = dictionary[key]!
     traverseOctree(atomIDs: values, origin: newOrigin, levelSize: levelSize / 2, depth: depth + 1)
   }
 }
 traverseOctree(atomIDs: topology.atoms.indices.map { $0 }, origin: .zero, levelSize: 4, depth: 1)
 */
