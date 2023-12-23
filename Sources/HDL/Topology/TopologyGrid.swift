//
//  TopologyGrid.swift
//  MolecularRenderer
//
//  Created by Philip Turner on 12/2/23.
//

#if arch(arm64)
typealias Half = Float16
#else
typealias Half = Float32
#endif

struct TopologyCell {
  var indices: [UInt32]?
  var covalentSearchRadius: Half = .zero
}

// Topology grids are transient and regenerated upon every function that
// requires them. This design choice decreases the complexity of state changes
// in Topology, although it may decrease performance. It may also increase
// algorithmic complexity for extremely tiny search lists. This performance
// issue can be fixed if the need arises, and there is a reasonable alternative.
struct TopologyGrid {
  let entities: [Entity]
  var cells: [TopologyCell] = []
  var entitiesToCellsMap: [SIMD2<UInt32>] = []
  
  // Origin and dimensions are in discrete multiples of cell width. To match
  // two grids to each other, they must have the exact same cell width.
  var cellWidth: Float
  var origin: SIMD3<Int32>
  var dimensions: SIMD3<Int32>
  
  init(entities: [Entity], cellWidth: Float) {
    self.entities = entities
    self.cellWidth = cellWidth
    precondition(cellWidth >= 0.1 && cellWidth <= 1, "Unreasonable cell width.")
    
    if entities.count == 0 {
      origin = .zero
      dimensions = .zero
    } else {
      var minimum: SIMD3<Float> = .init(repeating: .greatestFiniteMagnitude)
      var maximum: SIMD3<Float> = .init(repeating: -.greatestFiniteMagnitude)
      for entity in entities {
        let atom = entity.position
        minimum.replace(with: atom, where: atom .< minimum)
        maximum.replace(with: atom, where: atom .> maximum)
      }
      minimum /= cellWidth
      maximum /= cellWidth
      minimum.round(.down)
      maximum.round(.up)
      origin = SIMD3<Int32>(minimum)
      dimensions = SIMD3<Int32>(maximum - minimum)
    }
    
    // Use checking arithmetic to ensure the multiplication doesn't overflow.
    let cellCount = Int32(dimensions[0] * dimensions[1] * dimensions[2])
    cells = Array(repeating: TopologyCell(), count: Int(cellCount))
    entitiesToCellsMap = Array(repeating: .zero, count: entities.count)
    
    for atomID in entities.indices {
      let entity = entities[atomID]
      precondition(entity.storage.w != 0, "Empty entities are not permitted.")
      precondition(
        entity.storage.w == Float(Int8(entity.storage.w)),
        "Atomic number must be between 1 and 127.")
      
      var position = entity.position
      position /= cellWidth
      position.round(.down)
      let originDelta = SIMD3<Int32>(position) &- self.origin
      let cellID = self.createCellID(originDelta: originDelta)
      
      var mappedAtomID: Int = 0
      let atomID32 = UInt32(truncatingIfNeeded: atomID)
      if cells[cellID].indices == nil {
        cells[cellID].indices = Array(unsafeUninitializedCapacity: 8) {
          $0.baseAddress.unsafelyUnwrapped[0] = atomID32
          $1 = 1
        }
      } else {
        mappedAtomID = cells[cellID].indices!.count
        cells[cellID].indices!.append(atomID32)
      }
      let location = SIMD2<Int>(cellID, mappedAtomID)
      entitiesToCellsMap[atomID] = SIMD2(truncatingIfNeeded: location)
    }
    
    for cellID in cells.indices {
      guard let indices = cells[cellID].indices else {
        continue
      }
      var maxRadius: Float = .zero
      for atomID in indices {
        let entity = entities[Int(atomID)]
        let atomicNumber = Int(entity.storage.w)
        let covalentRadius = Element.covalentRadii[atomicNumber]
        maxRadius = max(maxRadius, covalentRadius)
      }
      cells[cellID].covalentSearchRadius = Half(maxRadius)
    }
  }
  
  // The input is the position relative to the grid's origin.
  @inline(__always)
  func createCellID(originDelta: SIMD3<Int32>) -> Int {
    var coords = originDelta
    coords.replace(with: SIMD3.zero, where: coords .< 0)
    coords.replace(with: dimensions &- 1, where: coords .>= dimensions)
    
    let x = coords.x
    let y = coords.y &* dimensions[0]
    let z = coords.z &* dimensions[0] &* dimensions[1]
    return Int(x &+ y &+ z)
  }
  
  static func createOptimalCellWidth(_ entities: [Entity]) -> Float {
    var sum: Float = 0
    var count: Int = 0
    for entity in entities where entity.storage.w > 1 {
      sum += Element.covalentRadii[Int(entity.storage.w)]
      count += 1
    }
    guard count > 0 else {
      print("WARNING: Zero entities when creating optimal cell width.")
      return 0.25
    }
    
    // Find the optimal heuristic for cell width when a bottleneck arises.
    var output = sum / Float(count)
    output *= 4
    if output < 0.25 { output = 0.25 }
    else if output > 0.50 { output = 0.50 }
    return output
  }
  
  func createMaxCovalentSearchRadius() -> Half {
    var searchRadius: Half = .zero
    for cell in cells {
      searchRadius = max(searchRadius, cell.covalentSearchRadius)
    }
    return searchRadius
  }
}
