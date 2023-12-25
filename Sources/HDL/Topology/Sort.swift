//
//  Sort.swift
//  MolecularRenderer
//
//  Created by Philip Turner on 12/2/23.
//

// Notes for when you get around to optimizing this:
//
// A test should compare the grid sorter to an alternative, recursive
// implementation. It first measures cold-start speed, then speed once the
// grid sorter can take advantage of the set already being sorted.
//
// The eventually chosen implementation might be a hybrid of these. It might
// temporarily generate Morton indices to check whether segments of the atoms
// are already sorted. It could also switch between different methods at
// different levels of the hierarchy.
//
// Note: random shuffling creates a different data distribution than the
// typical distribution from Lattice. It may unfairly increase execution time
// of the grid sorter.

extension Topology {
  @discardableResult
  public mutating func sort() -> [UInt32] {
    let grid = GridSorter(atoms: atoms)
    let reordering = grid.mortonReordering()
    let previousAtoms = atoms
    for originalID in reordering.indices {
      let reorderedID32 = reordering[originalID]
      let reorderedID = Int(truncatingIfNeeded: reorderedID32)
      atoms[reorderedID] = previousAtoms[originalID]
    }
    
    for i in bonds.indices {
      let bond = SIMD2<Int>(truncatingIfNeeded: bonds[i])
      var newBond: SIMD2<UInt32> = .zero
      newBond[0] = reordering[bond[0]]
      newBond[1] = reordering[bond[1]]
      newBond = SIMD2(newBond.min(), newBond.max())
      bonds[i] = newBond
    }
    bonds.sort {
      if $0.x != $1.x {
        return $0.x < $1.x
      } else {
        return $0.y < $1.y
      }
    }
    return reordering
  }
}

struct GridSorter {
  let atoms: [Entity]
  var cells: [[UInt32]] = []
  var atomsToCellsMap: [SIMD2<UInt32>] = []
  
  // Origin and dimensions are in discrete multiples of cell width.
  var cellWidth: Float
  var origin: SIMD3<Float>
  var dimensions: SIMD3<Int32>
  
  init(atoms: [Entity], cellWidth: Float = 1.0) {
    self.atoms = atoms
    self.cellWidth = cellWidth
    
    if atoms.count == 0 {
      origin = .zero
      dimensions = .zero
    } else {
      var minimum = SIMD3<Float>(repeating: .greatestFiniteMagnitude)
      var maximum = -minimum
      for atom in atoms {
        let position = atom.position
        minimum.replace(with: position, where: position .< minimum)
        maximum.replace(with: position, where: position .> maximum)
      }
      minimum /= cellWidth
      maximum /= cellWidth
      origin = minimum
      dimensions = SIMD3<Int32>((maximum - minimum).rounded(.up))
    }
    
    // Use checking arithmetic to ensure the multiplication doesn't overflow.
    let cellCount = Int32(dimensions[0] * dimensions[1] * dimensions[2])
    cells = Array(repeating: [], count: Int(cellCount))
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
      let originDelta = SIMD3<Int32>(position - origin)
      let cellID = self.createCellID(originDelta: originDelta)
      
      let mappedAtomID = cells[cellID].count
      cells[cellID].append(UInt32(truncatingIfNeeded: atomID))
      
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

extension GridSorter {
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
      for atomID in cells[cellID] {
        let atom = atoms[Int(atomID)]
        let scaledPosition = (atom.position - origin) / cellWidth
        let floorPosition = scaledPosition.rounded(.down)
        let originDelta = SIMD3<Int32>(floorPosition)
        guard createCellID(originDelta: originDelta) == cellID else {
          fatalError("Atom was not in the correct cell.")
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
