//
//  Topology+Match.swift
//
//
//  Created by Philip Turner on 12/10/23.
//

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
