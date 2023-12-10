//
//  Topology+Search.swift
//
//
//  Created by Philip Turner on 12/10/23.
//

// Efficient grid searching for match() and the bond generation algorithm.

// MARK: - Data Types

struct TopologyValue {
  var storage: SIMD3<Float>
  
  @inline(__always)
  init(cell: UInt32, offset: UInt32, distance: Float) {
    storage = SIMD3(
      Float(bitPattern: cell),
      Float(bitPattern: offset),
      distance)
  }
  
  var cell: UInt32 { storage.x.bitPattern }
  var offset: UInt32 { storage.y.bitPattern }
  var distance: Float { storage.z }
}

struct TopologyKey {
  var x: SIMD8<Float>
  var y: SIMD8<Float>
  var z: SIMD8<Float>
  var origin: SIMD3<Int32>
}

// A data structure to allocate once and recycle every outer loop iteration.
//
// Don't conditionally check whether each candidate is within range during the
// inner loop. Instead, sort first. Scan only a few chunks of the array that
// fall out of bounds, removing potentially large chunks all at once.
struct TopologyValues {
  var array0: [TopologyValue] = []
  var array1: [TopologyValue] = []
  var array2: [TopologyValue] = []
  var array3: [TopologyValue] = []
  var array4: [TopologyValue] = []
  var array5: [TopologyValue] = []
  var array6: [TopologyValue] = []
  var array7: [TopologyValue] = []
  
  init() {
    array0.reserveCapacity(1024)
    array1.reserveCapacity(1024)
    array2.reserveCapacity(1024)
    array3.reserveCapacity(1024)
    array4.reserveCapacity(1024)
    array5.reserveCapacity(1024)
    array6.reserveCapacity(1024)
    array7.reserveCapacity(1024)
  }
  
  mutating func append(cell: Int, offset: Int, distances: SIMD8<Float>) {
    let cell32 = UInt32(truncatingIfNeeded: cell)
    let offset32 = UInt32(truncatingIfNeeded: offset)
    
    @inline(__always)
    func makeValue(_ lane: Int) -> TopologyValue {
      TopologyValue(cell: cell32, offset: offset32, distance: distances[lane])
    }
    array0.append(makeValue(0))
    array1.append(makeValue(1))
    array2.append(makeValue(2))
    array3.append(makeValue(3))
    array4.append(makeValue(4))
    array5.append(makeValue(5))
    array6.append(makeValue(6))
    array7.append(makeValue(7))
  }
  
  mutating func sort(validKeys: Int) {
    func sortArray(_ keyID: Int, _ array: inout [TopologyValue]) {
      if keyID >= validKeys {
        return
      }
      array.sort(by: {
        $0.distance > $1.distance
      })
      
      let initialCount = array.count
      for i in (0..<16).reversed() {
        let startingPoint = i * initialCount / 16
        if startingPoint >= array.count {
          continue
        }
        
        if array[startingPoint].distance > 0.5 {
          array.removeLast(array.count - startingPoint)
        } else {
          while array.count > startingPoint,
                  array.last!.distance > 0.5 {
            array.removeLast()
          }
          return
        }
      }
    }
    sortArray(0, &array0)
    sortArray(1, &array1)
    sortArray(2, &array2)
    sortArray(3, &array3)
    sortArray(4, &array4)
    sortArray(5, &array5)
    sortArray(6, &array6)
    sortArray(7, &array7)
  }
  
  func withContents(
    validKeys: Int,
    _ closure: ([[TopologyValue]]) throws -> Void
  ) rethrows {
    var references = [
      array0, array1, array2, array3, array4, array5, array6, array7
    ]
    references.removeLast(8 - validKeys)
    try closure(references)
    references.removeAll()
  }
  
  mutating func clear() {
    array0.removeAll(keepingCapacity: true)
    array1.removeAll(keepingCapacity: true)
    array2.removeAll(keepingCapacity: true)
    array3.removeAll(keepingCapacity: true)
    array4.removeAll(keepingCapacity: true)
    array5.removeAll(keepingCapacity: true)
    array6.removeAll(keepingCapacity: true)
    array7.removeAll(keepingCapacity: true)
  }
}

// MARK: - Traversal

// TODO: Write a function that traverses a grid, given eight keys and an
// 'inout' values object.
