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
  
  @inline(__always)
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

extension TopologyGrid {
  func search(
    origin: SIMD3<Int32>,
    keys: ArraySlice<Entity>,
    values: inout TopologyValues
  ) {
    // Set up the values.
    values.clear()
    
    // Set up the keys.
    guard keys.count >= 0 else {
      fatalError("No method for handling the edge case where keys = 0.")
    }
    guard keys.count <= 8 else {
      fatalError("Too many keys.")
    }
    var x: SIMD8<Float> = .init(repeating: .nan)
    var y: SIMD8<Float> = .init(repeating: .nan)
    var z: SIMD8<Float> = .init(repeating: .nan)
    for keyID in keys.indices {
      let key = keys[keyID]
      x[keyID] = key.position.x
      y[keyID] = key.position.y
      z[keyID] = key.position.z
    }
    let minimum = SIMD3(x.min(), y.min(), z.min())
    let maximum = SIMD3(x.max(), y.max(), z.max())
    guard all(maximum - minimum .< 0.5 + 1e-3) else {
      fatalError("Points were not confined within 0.5 nm.")
    }
    let originDelta = (minimum / 0.5).rounded(.down) - SIMD3<Float>(origin)
    guard all(originDelta * originDelta .< 1e-3 * 1e-3) else {
      fatalError("Points did not match origin.")
    }
    
    // Iterate over all cells within the search radius.
    for xDelta in Int32(-1)...1 {
      for yDelta in Int32(-1)...1 {
        for zDelta in Int32(-1)...1 {
          // TODO: Profile the overhead of traversal against the setup code
          // surrounding it, such as array sorting. Determine whether the total
          // overhead of a reasonable size for 'match' is enough to parallelize
          // over multiple CPU cores.
          loop(SIMD3(xDelta, yDelta, zDelta))
        }
      }
    }
    @inline(__always)
    func loop(_ delta: SIMD3<Int32>) {
      // Return early if out of bounds.
      let originDelta = delta &+ (origin &- self.origin)
      guard all((originDelta .>= 0) .& (originDelta .< dimensions)) else {
        return
      }
      
      let cellID = self.createCellID(originDelta: originDelta)
      let cell = self.cells[cellID]
      guard let entities = cell.entities else {
        return
      }
      
      for atomID in entities.indices {
        let entity = entities[atomID]
        var deltaX = entity.storage.x - x
        var deltaY = entity.storage.y - y
        var deltaZ = entity.storage.z - z
        deltaX *= deltaX
        deltaY *= deltaY
        deltaZ *= deltaZ
        
        let distance = (deltaX + deltaY + deltaZ).squareRoot()
        values.append(cell: cellID, offset: atomID, distances: distance)
      }
    }
  }
}
