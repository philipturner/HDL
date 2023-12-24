//
//  Topology+Archive.swift
//
//
//  Created by Philip Turner on 12/23/23.
//

// An archive of old code to reference for operations like batched grid
// traversal and Morton reordering.

// MARK: - Dispatching Batched Neighbor Searches

#if false

// TODO: Remove this API, instead having the user rule out potential matches
// after fetching the candidates.
public typealias MatchType = (
  _ input: Entity,
  _ candidate: Entity
) -> Bool

extension Topology {
  // TODO: Match should take a radius for variable compute cost, with a default
  // of 0.5 nm. Therefore, grid cells are 0.25 nm by default for Match, but cell
  // size can vary. Cell size = max(0.25, radius / 2).
  public func match(
    _ input: [Entity]
  ) -> [UInt32?] {
    match(input) { _, _ in return true }
  }
  
  // TODO: Implement the linear prefix sum, use it instead of Morton order for
  // the time being. Make a Swift package test for this.
  //
  // TODO: Make a Swift package test for this function. Test a grid against
  // itself, another example where all atoms are shifted slightly, then another
  // example that rejects perfect matches.
  public func match(
    _ input: [Entity],
    _ closure: MatchType
  ) -> [UInt32?] {
    var matchResults: [SIMD2<UInt32>: UInt32] = [:]
    var mappedLocations: [SIMD2<UInt32>] = []
    let inputGrid = TopologyGrid(
      entities: input, mappedLocations: &mappedLocations)
    var values = TopologyValues()
    
    for z in 0..<inputGrid.dimensions.z {
      for y in 0..<inputGrid.dimensions.y {
        for x in 0..<inputGrid.dimensions.x {
          loop(SIMD3(x, y, z))
        }
      }
    }
    func loop(_ delta: SIMD3<Int32>) {
      let absoluteOrigin = inputGrid.origin &+ delta
      let cellID = inputGrid.createCellID(originDelta: delta)
      let cell = inputGrid.cells[cellID]
      guard let cellEntities = cell.entities else {
        return
      }
      
      var atomStart: Int = 0
      let atomEnd = cellEntities.count
      while atomStart < atomEnd {
        defer { atomStart += 8 }
        
        // TODO: Profile this (in release mode)...
        let validKeys = min(8, atomEnd - atomStart)
        let keys = cellEntities[atomStart..<atomStart + validKeys]
        self.grid.search(origin: absoluteOrigin, keys: keys, values: &values)
        
        // TODO: Against this (in release mode)...
        var candidateMatch: SIMD8<UInt32> = .init(repeating: .max)
        values.sort(validKeys: validKeys)
        values.withContents(validKeys: validKeys) { contents in
          for keyOffset in 0..<validKeys {
            let input = cellEntities[atomStart + keyOffset]
            let values = contents[keyOffset]
          inner:
            for valueID in values.indices {
              let value = values[valueID]
              let candidateCell = self.grid.cells[Int(value.cell)]
              guard let candidateCellEntities = candidateCell.entities else {
                fatalError("Candidate was against empty cell.")
              }
              let candidate = candidateCellEntities[Int(value.offset)]
              if closure(input, candidate) {
                let selfAddress = candidateCell.prefixSumMorton + value.offset
                candidateMatch[keyOffset] = selfAddress
                break inner
              }
            }
          }
        }
        for keyOffset in 0..<validKeys
        where candidateMatch[keyOffset] != .max {
          let location = SIMD2(cellID, atomStart + keyOffset)
          let location32 = SIMD2<UInt32>(truncatingIfNeeded: location)
          matchResults[location32] = candidateMatch[keyOffset]
        }
      }
    }
    
    var output = [UInt32?](repeating: nil, count: input.count)
    for atomID in input.indices {
      let location = mappedLocations[atomID]
      guard let match = matchResults[location] else {
        continue
      }
      output[atomID] = match
    }
    return output
  }
}

#endif

// MARK: - Batched Neighbor Searching

#if false

//struct TopologyValue {
//  var storage: SIMD3<Float>
//
//  @inline(__always)
//  init(cell: UInt32, offset: UInt32, distance: Float) {
//    storage = SIMD3(
//      Float(bitPattern: cell),
//      Float(bitPattern: offset),
//      distance)
//  }
//
//  var cell: UInt32 { storage.x.bitPattern }
//  var offset: UInt32 { storage.y.bitPattern }
//  var distance: Float { storage.z }
//}

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
  // The input is the absolute position in 3D space.
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
    for zDelta in Int32(-1)...1 {
      for yDelta in Int32(-1)...1 {
        for xDelta in Int32(-1)...1 {
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

#endif

// MARK: - Morton Reordering

#if false // This old code is a useful reference for Morton reordering.

/// Rounds an integer up to the nearest power of 2.
fileprivate func roundUpToPowerOf2(_ input: Int) -> Int {
  1 << (Int.bitWidth - max(0, input - 1).leadingZeroBitCount)
}

/// Rounds an integer down to the nearest power of 2.
fileprivate func roundDownToPowerOf2(_ input: Int) -> Int {
  1 << (Int.bitWidth - 1 - input.leadingZeroBitCount)
}

// Each atom stores the position in the first three lanes, and the atom type
// in the fourth. -1 corresponds to "sigma bond". When exporting, delete all
// atoms whose mask slot is zero.
struct Atom {
  var data: SIMD4<Float>
  var position: SIMD3<Float> { unsafeBitCast(data, to: SIMD3<Float>.self) }
  var atomicNumber: Float {
    get { data.w }
    set { data.w = newValue }
  }
}

struct Grid {
  // Tolerance for floating-point error or error in lattice constants when
  // comparing atom positions.
  // TODO: Override the nominal lonsdaleite lattice constants so they align
  // perfectly with diamond (111) surfaces, using 0.357 nm for diamond.
  static let epsilon: Float = 0.005
  
  // When appending an atom to a cell, add (1 << atoms.count) to the bitmask.
  // Measure the array count before, not after, adding the new atom.
  struct Cell {
    var mask: UInt64
    var atoms: [Atom]?
    
    func match(_ atom: Atom) -> Int32 {
      var match: Int32 = -1
      withUnsafeTemporaryAllocation(of: SIMD64<Int8>.self, capacity: 1) {
        let matchMemory = UnsafeMutablePointer<Int8>(
          OpaquePointer($0.baseAddress)).unsafelyUnwrapped
        for i in 0..<atoms.unsafelyUnwrapped.count {
          let delta = atoms.unsafelyUnwrapped[i].position - atom.position
          let distance = (delta * delta).sum()
          let index = Int8(truncatingIfNeeded: i)
          matchMemory[i] = (distance < epsilon) ? index : 0
        }
        
        let matchMask = $0.baseAddress.unsafelyUnwrapped
        match = Int32(truncatingIfNeeded: matchMask.pointee.max())
      }
      return match
    }
    
    mutating func append(_ atom: Atom) {
      if atoms == nil {
        atoms = []
      }
      mask |= 1 << atoms.unsafelyUnwrapped.count
      atoms?.append(atom)
    }
    
    mutating func merge(_ atom: Atom, match: Int32) {
      let index = Int(truncatingIfNeeded: match)
      let other = atoms.unsafelyUnwrapped[index]
      let selectMask: UInt64 = 1 << index
      
      if self.mask & selectMask == 0 ||
          other.atomicNumber < atom.atomicNumber {
        atoms?[index].atomicNumber = atom.atomicNumber
      }
      self.mask |= selectMask
    }
  }
  
  var cellWidth: Float
  var dimensions: SIMD3<Int32>
  var origin: SIMD3<Float>
  var cells: [Cell]
  
  func address(coords: SIMD3<Int32>) -> Int {
    var output = coords.z * dimensions.y * dimensions.x
    output &+= coords.y * dimensions.x
    output &+= coords.x
    return Int(truncatingIfNeeded: output)
  }
  
  // Incremental migration path: swap out the backend of Lattice and Solid with
  // a grid in lattice-space, with cell width 1. Then, change to use 0.357.
  init(
    cellWidth: Float,
    dimensions: SIMD3<Int32>,
    origin: SIMD3<Float> = .zero
  ) {
    self.cellWidth = cellWidth
    self.dimensions = dimensions
    self.origin = origin
    self.cells = []
    
    for _ in 0..<dimensions.z {
      for _ in 0..<dimensions.y {
        for _ in 0..<dimensions.x {
          cells.append(Cell(mask: 0, atoms: nil))
        }
      }
    }
  }
  
  // Merging appends are more expensive, as they check for and remove
  // duplicated atoms. When merging multiple grids, consider setting 'merging'
  // to 'false' when adding atoms from the first grid.
  mutating func append(contentsOf other: [Atom], merging: Bool = false) {
  outer:
    for atom in other {
      let offset = atom.position - origin
      var scaledOffset = offset / cellWidth
      if any(scaledOffset .< 0 - Self.epsilon) ||
          any(scaledOffset .> SIMD3<Float>(dimensions) + Self.epsilon) {
        fatalError("Atom was outside of grid bounds.")
      }
      
      var coords: SIMD3<Int32> = .init(scaledOffset.rounded(.down))
      coords.clamp(lowerBound: SIMD3.zero, upperBound: dimensions &- 1)
      if merging {
        let delta = scaledOffset - SIMD3<Float>(coords)
        var starts: SIMD3<Int32> = .init(repeating: -1)
        starts.replace(with: SIMD3.zero, where: delta .>= 0.5)
        
        for z in starts.z...starts.z + 1 {
          for y in starts.y...starts.y + 1 {
            for x in starts.x...starts.x + 1 {
              var adjusted = coords &+ SIMD3(x, y, z)
              if all(adjusted .>= 0 .& adjusted .<= dimensions &- 1) {
                let address = self.address(coords: coords)
                let match = cells[address].match(atom)
                if match > -1 {
                  cells[address].merge(atom, match: match)
                  continue outer
                }
              }
            }
          }
        }
      }
      
      let address = self.address(coords: coords)
      cells[address].append(atom)
    }
    
    // Always check that the number of atoms per cell never exceeds 64.
    for cell in cells where (cell.atoms?.count ?? 0) > 64 {
      fatalError("Cell had more than 64 atoms.")
    }
  }
}

extension Grid {
  private func cellsPrefixSum() -> [Int32] {
    var output: [Int32] = []
    var sum: Int = 0
    for cell in self.cells {
      output.append(Int32(truncatingIfNeeded: sum))
      sum += cell.atoms?.count ?? 0
    }
    output.append(Int32(truncatingIfNeeded: sum))
    return output
  }
  
  /// WARNING: This is a computed property. Access it sparingly.
  var atoms: [Atom] {
    var output: [Atom] = []
    for cell in self.cells where cell.atoms != nil {
      output.append(contentsOf: cell.atoms.unsafelyUnwrapped)
    }
    return output
  }
  
  /// WARNING: This is a computed property. Access it sparingly.
  ///
  /// Returns a map, of the atoms' new locations when arranged in Morton order.
  /// Usually, a 4/8 grid will be regenerated at 4/24 resolution to maximize
  /// spatial locality.
  var mortonReordering: [Int32] {
    // TODO: Check that results are correct by animating addition of atoms to
    // the scene, one at a time.
    
    // Interleave, sort, then deinterleave.
    //
    // Interleaving algorithm:
    // https://stackoverflow.com/a/18528775
    //
    // Deinterleaving algorithm:
    // https://stackoverflow.com/a/28358035
    
    @inline(__always)
    func morton_interleave(_ input: Int32) -> UInt64 {
      var x = UInt64(truncatingIfNeeded: input & 0x1fffff)
      x = (x | x &<< 32) & 0x1f00000000ffff
      x = (x | x &<< 16) & 0x1f0000ff0000ff
      x = (x | x &<< 8) & 0x100f00f00f00f00f
      x = (x | x &<< 4) & 0x10c30c30c30c30c3
      x = (x | x &<< 2) & 0x1249249249249249
      return x
    }
    
    @inline(__always)
    func morton_deinterleave(_ input: UInt64) -> Int32 {
      var x = input & 0x9249249249249249
      x = (x | (x &>> 2))  & 0x30c30c30c30c30c3
      x = (x | (x &>> 4))  & 0xf00f00f00f00f00f
      x = (x | (x &>> 8))  & 0x00ff0000ff0000ff
      x = (x | (x &>> 16)) & 0xffff00000000ffff
      x = (x | (x &>> 32)) & 0x00000000ffffffff
      return Int32(truncatingIfNeeded: x)
    }
    
    var mortonIndices: [UInt64] = []
    for z in 0..<dimensions.z {
      for y in 0..<dimensions.y {
        for x in 0..<dimensions.x {
          let address = self.address(coords: SIMD3(x, y, z))
          if cells[address].atoms?.count ?? 0 > 0 {
            let morton_x = morton_interleave(x)
            let morton_y = morton_interleave(y)
            let morton_z = morton_interleave(z)
            
            var mortonIndex = morton_x
            mortonIndex |= morton_y << 1
            mortonIndex |= morton_z << 2
            mortonIndices.append(mortonIndex)
          }
        }
      }
    }
    mortonIndices.sort()
    
    let prefixSum = cellsPrefixSum()
    var outputMap: [Int32] = .init(repeating: -1, count: Int(prefixSum.last!))
    var outputMapSize: Int = 0
    for mortonIndex in mortonIndices {
      let x = morton_deinterleave(mortonIndex)
      let y = morton_deinterleave(mortonIndex >> 1)
      let z = morton_deinterleave(mortonIndex >> 2)
      
      let address = self.address(coords: SIMD3(x, y, z))
      guard let atoms = cells[address].atoms else {
        fatalError("Morton indexing happened incorrectly.")
      }
      let inputPrefix = prefixSum[address]
      let outputPrefix = outputMapSize
      outputMapSize += atoms.count
      
      for i in 0..<atoms.count {
        outputMap[outputPrefix + i] = inputPrefix + 1
      }
    }
    return outputMap
  }
  
  /// WARNING: This is a computed property. Access it sparingly.
  ///
  /// If a carbon has more than four bonds, or less than two, this returns
  /// `[-2, -2, -2, -2]`. Therefore, anything checking for invalid bonds should
  /// simply check whether the bond index is less than zero.
  var bonds: [SIMD4<Int32>] {
    // 0.1545855 nm - carbon-carbon bond length, corrected to match diamond
    fatalError("Not implemented.")
  }
  
  /// Returns places where passivation must occur.
  ///
  /// The first 3 slots store the position. The last one stores the atom index
  /// as a floating point number. Passivations corresponding to the same atom
  /// are stored contiguously in memory. The atoms resulting from the
  /// passivations should be run through `.append` in a separate call, which
  /// prevents any atoms from being duplicated. Then, the bond topology should
  /// be regenerated, as the atoms are in a different order.
  ///
  /// If the atom is already passivated at a specific location, that location isn't
  /// returned. Instead, you can check whether a specific atom from the bond
  /// map is hydrogen or a single bond. Duplicated passivations will connect two
  /// carbons to the same hydrogen.
  ///
  /// NOTE: In the DSL, `Passivate` must override previous passivations.
  func passivations(bonds: [SIMD4<Int32>]) -> [SIMD4<Float>] {
    // 0.1545855 nm - carbon-carbon bond length, corrected to match diamond
    fatalError("Not implemented.")
  }
  
  // When adding support for 5-carbon rings (implementing the new MM4):
  //
  // Add a function that returns deltas to atom positions, to fix hydrogens when
  // reconstructing (100) surfaces. This function accepts the bond map as input.
}

#endif

// MARK: - Valence Orbitals

#if false
func createAtom(atoms: [MRAtom], atomID: Int) {
  precondition(atomID > -1)
  let thisAtom = atoms[atomID]
  
  let newAtomID = Int32(self.atoms.count)
  newIndicesMap[atomID] = newAtomID
  self.atoms.append(thisAtom)
  
  var neighborTypes: [Int] = []
  var neighborCenters: [SIMD3<Float>] = []
  for j in 0..<centerTypes[atomID] {
    let index = Int(centerNeighbors[atomID][j])
    neighborTypes.append(centerTypes[index])
    neighborCenters.append(atoms[index].origin)
    
    // Change this; store the bonds with indices being sorted inside the
    // bond, but only add a bond when the neighbor is already inside the
    // final list.
    let newNeighborID = newIndicesMap[index]
    guard newNeighborID > -1 else {
      continue
    }
    var newBond: SIMD2<Int32> = .zero
    newBond[0] = min(newAtomID, newNeighborID)
    newBond[1] = max(newAtomID, newNeighborID)
    bonds.append(newBond)
  }
  
  let valenceElectrons = Constants.valenceElectrons(
    element: thisAtom.element)
  if centerTypes[atomID] > valenceElectrons {
    fatalError("Too many bonds.")
  }
  
  var totalBonds = centerTypes[atomID]
  func addHydrogen(direction: SIMD3<Float>) {
    guard totalBonds < valenceElectrons else {
      return
    }
    totalBonds += 1
    
    let bondLength = Constants.bondLengths[
      [1, thisAtom.element]]!.average
    let hydrogenCenter = thisAtom.origin + bondLength * direction
    let hydrogenID = Int32(self.atoms.count)
    
    self.atoms.append(MRAtom(origin: hydrogenCenter, element: 1))
    self.bonds.append(SIMD2(Int32(newAtomID), hydrogenID))
  }
  
  switch centerTypes[atomID] {
  case 4:
    break
  case 3:
    let sideAB = neighborCenters[1] - neighborCenters[0]
    let sideAC = neighborCenters[2] - neighborCenters[0]
    var normal = _cross_platform_normalize(_cross_platform_cross(sideAB, sideAC))
    
    let deltaA = thisAtom.origin - neighborCenters[0]
    if _cross_platform_dot(normal, deltaA) < 0 {
      normal = -normal
    }
    
    addHydrogen(direction: normal)
  case 2:
    let midPoint = (neighborCenters[1] + neighborCenters[0]) / 2
    guard _cross_platform_distance(midPoint, thisAtom.origin) > 0.001 else {
      fatalError("sp3 carbons are too close to 180 degrees.")
    }
    
    let normal = _cross_platform_normalize(thisAtom.origin - midPoint)
    let axis = _cross_platform_normalize(neighborCenters[1] - midPoint)
    for angle in [-sp3BondAngle / 2, sp3BondAngle / 2] {
      let rotation = Quaternion<Float>(angle: angle, axis: axis)
      let direction = rotation.act(on: normal)
      addHydrogen(direction: direction)
    }
  case 1:
    guard neighborTypes[0] > 1 else {
      fatalError("Cannot determine structure of primary carbon.")
    }
    
    let j = Int(centerNeighbors[atomID][0])
    var referenceIndex: Int?
    for k in 0..<neighborTypes[0] {
      let index = Int(centerNeighbors[j][k])
      if atomID != index {
        referenceIndex = index
        break
      }
    }
    guard let referenceIndex else {
      fatalError("Could not find valid neighbor index.")
    }
    let referenceCenter = atoms[referenceIndex].origin
    let normal = _cross_platform_normalize(thisAtom.origin - atoms[j].origin)
    
    let referenceDelta = atoms[j].origin - referenceCenter
    var orthogonal = referenceDelta - normal * _cross_platform_dot(normal, referenceDelta)
    guard _cross_platform_length(orthogonal) > 0.001 else {
      fatalError("sp3 carbons are too close to 180 degrees.")
    }
    orthogonal = _cross_platform_normalize(orthogonal)
    let axis = _cross_platform_cross(normal, orthogonal)
    
    var directions: [SIMD3<Float>] = []
    let firstHydrogenRotation = Quaternion<Float>(
      angle: .pi - sp3BondAngle, axis: axis)
    directions.append(firstHydrogenRotation.act(on: normal))
    
    let secondHydrogenRotation = Quaternion<Float>(
      angle: 120 * .pi / 180, axis: normal)
    directions.append(secondHydrogenRotation.act(on: directions[0]))
    directions.append(secondHydrogenRotation.act(on: directions[1]))
    
    for direction in directions {
      addHydrogen(direction: direction)
    }
  default:
    fatalError("This should never happen.")
  }
}
#endif
