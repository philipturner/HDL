//
//  TopologyGrid.swift
//  MolecularRenderer
//
//  Created by Philip Turner on 12/2/23.
//

// MARK: - Declaration

#if arch(arm64)
typealias Half = Float16
#else
typealias Half = Float32
#endif

struct TopologyCell {
  var indices: [UInt32] = []
}

// Topology grids are transient and regenerated upon every function that
// requires them. This design choice decreases the complexity of state changes
// in Topology, although it may decrease performance. It may also increase
// algorithmic complexity for extremely tiny search lists. This performance
// issue can be fixed if the need arises, and there is a reasonable alternative.
struct TopologyGrid {
  let atoms: [Entity]
  var cells: [TopologyCell] = []
  var atomsToCellsMap: [SIMD2<UInt32>] = []
  
  // Origin and dimensions are in discrete multiples of cell width.
  var cellWidth: Float
  var origin: SIMD3<Int32>
  var dimensions: SIMD3<Int32>
  
  init(atoms: [Entity], cellWidth: Float = 1.0) {
    self.atoms = atoms
    self.cellWidth = cellWidth
    
    if atoms.count == 0 {
      origin = .zero
      dimensions = .zero
    } else {
      var minimum: SIMD3<Float> = .init(repeating: .greatestFiniteMagnitude)
      var maximum: SIMD3<Float> = .init(repeating: -.greatestFiniteMagnitude)
      for atom in atoms {
        let position = atom.position
        minimum.replace(with: position, where: position .< minimum)
        maximum.replace(with: position, where: position .> maximum)
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
    atomsToCellsMap = Array(repeating: .zero, count: atoms.count)
    
    for atomID in atoms.indices {
      let atom = atoms[atomID]
      precondition(atom.storage.w != 0, "Empty entities are not permitted.")
      precondition(
        atom.storage.w == Float(Int8(atom.storage.w)),
        "Atomic number must be between 1 and 127.")
      
      var position = atom.position
      position /= cellWidth
      position.round(.down)
      let originDelta = SIMD3<Int32>(position) &- self.origin
      let cellID = self.createCellID(originDelta: originDelta)
      
      let mappedAtomID = cells[cellID].indices.count
      cells[cellID].indices.append(UInt32(truncatingIfNeeded: atomID))
      
      let location = SIMD2<Int>(cellID, mappedAtomID)
      atomsToCellsMap[atomID] = SIMD2(truncatingIfNeeded: location)
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
}

// MARK: - Morton Reordering

// Source: https://stackoverflow.com/a/18528775
@inline(__always)
private func morton_interleave(_ input: Int32) -> UInt64 {
  var x = UInt64(truncatingIfNeeded: input & 0x1fffff)
  x = (x | x &<< 32) & 0x1f00000000ffff
  x = (x | x &<< 16) & 0x1f0000ff0000ff
  x = (x | x &<< 8) & 0x100f00f00f00f00f
  x = (x | x &<< 4) & 0x10c30c30c30c30c3
  x = (x | x &<< 2) & 0x1249249249249249
  return x
}

// Source: https://stackoverflow.com/a/28358035
@inline(__always)
private func morton_deinterleave(_ input: UInt64) -> Int32 {
  var x = input & 0x9249249249249249
  x = (x | (x &>> 2))  & 0x30c30c30c30c30c3
  x = (x | (x &>> 4))  & 0xf00f00f00f00f00f
  x = (x | (x &>> 8))  & 0x00ff0000ff0000ff
  x = (x | (x &>> 16)) & 0xffff00000000ffff
  x = (x | (x &>> 32)) & 0x00000000ffffffff
  return Int32(truncatingIfNeeded: x)
}

// The inputs are 21-bit, unsigned, normalized integers.
@inline(__always)
private func morton_interleave(_ input: SIMD3<Int32>) -> UInt64 {
  let signed = SIMD3<Int32>(truncatingIfNeeded: input)
  let x = morton_interleave(signed.x) << 0
  let y = morton_interleave(signed.y) << 1
  let z = morton_interleave(signed.z) << 2
  return x | y | z
}

extension TopologyGrid {
  func mortonReordering() -> [UInt32] {
    var cellList: [SIMD2<UInt64>] = []
    for z in 0..<dimensions.z {
      for y in 0..<dimensions.y {
        for x in 0..<dimensions.x {
          let originDelta = SIMD3(x, y, z)
          let mortonCode = morton_interleave(originDelta)
          let cellID = createCellID(originDelta: originDelta)
          let cellID64 = UInt64(truncatingIfNeeded: cellID)
          cellList.append(SIMD2(cellID64, mortonCode))
        }
      }
    }
    cellList.sort { $0.y < $1.y }
    
    var output = [UInt32](repeating: .max, count: atoms.count)
    var mortonListSize: Int = 0
    
    var atomList: [SIMD2<UInt64>] = []
    for element in cellList {
      let cellID = Int(truncatingIfNeeded: element[0])
      let atomIndices = cells[cellID].indices
      for atomID in atomIndices {
        let atom = atoms[Int(atomID)]
        let scaledPosition = atom.position / cellWidth
        let floorPosition = scaledPosition.rounded(.down)
        let originDelta = SIMD3<Int32>(floorPosition) &- self.origin
        guard createCellID(originDelta: originDelta) == cellID else {
          fatalError("Atom not in the correct cell.")
        }
        
        let remainder = (scaledPosition - floorPosition) * Float(1 << 21)
        var remainderInt = SIMD3<Int32>(remainder)
        let maxValue = SIMD3<Int32>(repeating: (1 << 21 - 1))
        remainderInt.replace(with: 0, where: remainder .< 0)
        remainderInt.replace(with: maxValue, where: remainderInt .> maxValue)
        let mortonCode = morton_interleave(remainderInt)
        let atomID64 = UInt64(truncatingIfNeeded: atomID)
        atomList.append(SIMD2(atomID64, mortonCode))
      }
      atomList.sort { $0.y < $1.y }
      
      for element in atomList {
        let atomID = Int(truncatingIfNeeded: element[0])
        let mortonMapping = UInt32(truncatingIfNeeded: mortonListSize)
        output[atomID] = mortonMapping
        mortonListSize += 1
      }
      atomList.removeAll(keepingCapacity: true)
    }
    
    precondition(
      mortonListSize == atoms.count, "Morton reordered list was invalid.")
    if output.contains(.max) {
      fatalError("Morton reordered list was invalid.")
    }
    return output
  }
}
