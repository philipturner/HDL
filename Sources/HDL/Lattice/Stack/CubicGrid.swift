//
//  CubicGrid.swift
//  MolecularRenderer
//
//  Created by Philip Turner on 10/22/23.
//

struct CubicCell {
  // Multiply the plane's origin by [4, 4, 4].
  // Span: [0 -> h], [0 -> k], [0 -> l]
  static let x0 = SIMD8<Float>(0, 1, 0, 1, 2, 3, 2, 3)
  static let y0 = SIMD8<Float>(0, 1, 2, 3, 0, 1, 2, 3)
  static let z0 = SIMD8<Float>(0, 1, 2, 3, 2, 3, 0, 1)
  
  static let flags = SIMD8<UInt8>(
    1 << 0, 1 << 1, 1 << 2, 1 << 3,
    1 << 4, 1 << 5, 1 << 6, 1 << 7)
  
  // Binary mask corresponding to the plane's "one volume" and "zero volume".
  @inline(__always)
  static func intersect(
    origin: SIMD3<Float>,
    normal: SIMD3<Float>
  ) -> UInt8 {
    let scaledOrigin = origin * 4
    let scaledNormal = normal * 1
    
    let delta_x0 = x0 - scaledOrigin.x
    let delta_y0 = y0 - scaledOrigin.y
    let delta_z0 = z0 - scaledOrigin.z
    var dotProduct0 = delta_x0 * scaledNormal.x
    dotProduct0 += delta_y0 * scaledNormal.y
    dotProduct0 += delta_z0 * scaledNormal.z
    
    var mask0: SIMD8<Int32> = .zero
    mask0.replace(with: SIMD8(repeating: .max), where: dotProduct0 .> 0)
    let compressed = SIMD8<UInt8>(truncatingIfNeeded: mask0)
    return (compressed & CubicCell.flags).wrappedSum()
  }
}

struct CubicMask: LatticeMask {
  var mask: [UInt8]
  
  /// Create a mask using a plane.
  init(dimensions: SIMD3<Int32>, origin: SIMD3<Float>, normal: SIMD3<Float>) {
    // Initialize the mask with everything in the one volume, and filled. The
    // value should be overwritten somewhere in the inner loop.
    mask = Array(repeating: .max, count: Int(
      dimensions.x * dimensions.y * dimensions.z))
    if all(normal .== 0) {
      // This cannot be evaluated. It is a permissible escape hatch to create a
      // mask with no intersection.
      return
    }
    
    #if false
    // Slower version retained as a backup, in case the faster version has bugs.
    for z in 0..<dimensions.z {
      for y in 0..<dimensions.y {
        let baseAddress = (z &* dimensions.y &+ y) &* dimensions.x
        
        for x in 0..<dimensions.x {
          let lowerCorner = SIMD3<Float>(Float(x), Float(y), Float(z))
          
          let cellMask = CubicCell.intersect(
            origin: origin - lowerCorner, normal: normal)
          mask[Int(baseAddress &+ x)] = cellMask
        }
      }
    }
    #else
    mask.withUnsafeMutableBufferPointer { buffer in
      let opaque = OpaquePointer(buffer.baseAddress.unsafelyUnwrapped)
      let mask8 = UnsafeMutablePointer<UInt8>(opaque)
      let mask16 = UnsafeMutablePointer<UInt16>(opaque)
      let mask32 = UnsafeMutablePointer<UInt32>(opaque)
      let dims = dimensions
      
      @inline(__always)
      func intersect1(sector: SIMD3<Int32>) {
        let floatX = Float(sector.x)
        for z in sector.z..<min(dims.z, sector.z &+ 2) {
          for y in sector.y..<min(dims.y, sector.y &+ 2) {
            let baseAddress = (z &* dims.y &+ y) &* dims.x &+ sector.x
            do {
              let lowerCorner = SIMD3<Float>(floatX + 0, Float(y), Float(z))
              let cellMask = CubicCell.intersect(
                origin: origin - lowerCorner, normal: normal)
              mask8[Int(baseAddress &+ 0)] = cellMask
            }
            do {
              let lowerCorner = SIMD3<Float>(floatX + 1, Float(y), Float(z))
              let cellMask = CubicCell.intersect(
                origin: origin - lowerCorner, normal: normal)
              mask8[Int(baseAddress &+ 1)] = cellMask
            }
          }
        }
      }
      
      @inline(__always)
      func intersect2(sector: SIMD3<Int32>) {
        var loopBounds = dimensions &- sector
        loopBounds.replace(with: 2, where: loopBounds .> 2)
        for subSectorZ in 0..<Int32(loopBounds.z) {
          for subSectorY in 0..<Int32(loopBounds.y) {
            for subSectorX in 0..<Int32(loopBounds.x) {
              let sector2 = sector &+ SIMD3(
                subSectorX, subSectorY, subSectorZ) &* 2
              
              let permX: SIMD8<UInt8> = .init(0, 0, 0, 0, 2, 2, 2, 2)
              let permY: SIMD8<UInt8> = .init(0, 0, 2, 2, 0, 0, 2, 2)
              let permZ: SIMD8<UInt8> = .init(0, 2, 0, 2, 0, 2, 0, 2)
              var trialX = SIMD8(repeating: Float(sector2.x) - origin.x)
              var trialY = SIMD8(repeating: Float(sector2.y) - origin.y)
              var trialZ = SIMD8(repeating: Float(sector2.z) - origin.z)
              trialX += SIMD8<Float>(permX)
              trialY += SIMD8<Float>(permY)
              trialZ += SIMD8<Float>(permZ)
              
              var dotProduct = trialX * normal.x
              dotProduct += trialY * normal.y
              dotProduct += trialZ * normal.z
              let allNegative = all(dotProduct .< 0)
              let allPositive = all(dotProduct .> 0)
              
              if allPositive {
                // already initialized to 1111_1111
              } else if allNegative {
                for z in sector2.z..<min(dims.z, sector2.z &+ 2) {
                  for y in sector2.y..<min(dims.y, sector2.y &+ 2) {
                    var address = (z &* dims.y &+ y)
                    address = address &* (dimensions.x / 2) &+ (sector2.x / 2)
                    mask16[Int(address)] = .zero
                  }
                }
              } else {
                intersect1(sector: sector2)
              }
            }
          }
        }
      }
      
      let sectorsX = (dimensions.x &+ 3) / 4
      let sectorsY = (dimensions.y &+ 3) / 4
      let sectorsZ = (dimensions.z &+ 3) / 4
      for sectorZ in 0..<sectorsZ {
        for sectorY in 0..<sectorsY {
          for sectorX in 0..<sectorsX {
            let permX: SIMD8<UInt8> = .init(0, 0, 0, 0, 4, 4, 4, 4)
            let permY: SIMD8<UInt8> = .init(0, 0, 4, 4, 0, 0, 4, 4)
            let permZ: SIMD8<UInt8> = .init(0, 4, 0, 4, 0, 4, 0, 4)
            var trialX = SIMD8(repeating: Float(sectorX) * 4 - origin.x)
            var trialY = SIMD8(repeating: Float(sectorY) * 4 - origin.y)
            var trialZ = SIMD8(repeating: Float(sectorZ) * 4 - origin.z)
            trialX += SIMD8<Float>(permX)
            trialY += SIMD8<Float>(permY)
            trialZ += SIMD8<Float>(permZ)
            
            var dotProduct = trialX * normal.x
            dotProduct += trialY * normal.y
            dotProduct += trialZ * normal.z
            let allNegative = all(dotProduct .< 0)
            let allPositive = all(dotProduct .> 0)
            
            let sector = SIMD3(sectorX, sectorY, sectorZ) &* 4
            if allPositive {
              // already initialized to 1111_1111
            } else if allNegative {
              for z in sector.z..<min(dims.z, sector.z &+ 4) {
                for y in sector.y..<min(dims.y, sector.y &+ 4) {
                  var address = z &* dims.y &+ y
                  address = address &* (dimensions.x / 4) &+ sectorX
                  mask32[Int(address)] = .zero
                }
              }
            } else {
              intersect2(sector: sector)
            }
          }
        }
      }
    }
    #endif
  }
  
  static func &= (lhs: inout Self, rhs: Self) {
    guard lhs.mask.count == rhs.mask.count else {
      fatalError("Combined masks of different sizes.")
    }
    for elementID in lhs.mask.indices {
      lhs.mask[elementID] &= rhs.mask[elementID]
    }
  }
  
  static func |= (lhs: inout Self, rhs: Self) {
    guard lhs.mask.count == rhs.mask.count else {
      fatalError("Combined masks of different sizes.")
    }
    for elementID in lhs.mask.indices {
      lhs.mask[elementID] |= rhs.mask[elementID]
    }
  }
}

struct CubicGrid: LatticeGrid {
  var dimensions: SIMD3<Int32>
  var entityTypes: [SIMD8<Int8>]
  var squareSideLength: Float
  
  /// Create a mask using a plane.
  init(bounds: SIMD3<Float>, materialType: MaterialType) {
    guard all(bounds.rounded(.up) .== bounds) else {
      fatalError("Bounds were not integers.")
    }
    
    var repeatingUnit: SIMD8<Int8>
    switch materialType {
    case .elemental(let element):
      let scalar = Int8(clamping: element.rawValue)
      repeatingUnit = SIMD8(repeating: scalar)
      if element.rawValue == 79 {
        // Change gold to face-centered cubic by deleting some atoms in the
        // diamond cubic unit cell.
        repeatingUnit = SIMD8(79, 0, 79, 0, 79, 0, 79, 0)
      }
    case .checkerboard(let a, let b):
      let scalarA = Int8(clamping: a.rawValue)
      let scalarB = Int8(clamping: b.rawValue)
      let unit = unsafeBitCast(SIMD2(scalarA, scalarB), to: UInt16.self)
      let repeated = SIMD4<UInt16>(repeating: unit)
      repeatingUnit = unsafeBitCast(repeated, to: SIMD8<Int8>.self)
    }
    
    // Increase the bounds by a small amount, so atoms on the edge will be
    // present in the next cell.
    dimensions = SIMD3<Int32>((bounds + 0.001).rounded(.up))
    dimensions.replace(with: SIMD3.zero, where: dimensions .< 0)
    dimensions.x = (dimensions.x + 3) / 4 * 4
    entityTypes = Array(repeating: repeatingUnit, count: Int(
      dimensions.x * dimensions.y * dimensions.z))
    
    // Fetch the lattice constant using the 'Constant' API.
    squareSideLength = Constant(.square) { materialType }
    
    self.initializeBounds(bounds, normals: [
      SIMD3<Float>(-1, 0, 0),
      SIMD3<Float>(1, 0, 0),
      SIMD3<Float>(0, -1, 0),
      SIMD3<Float>(0, 1, 0),
      SIMD3<Float>(0, 0, -1),
      SIMD3<Float>(0, 0, 1),
    ])
  }
  
  mutating func replace(with other: Int8, where mask: CubicMask) {
    let newValue = SIMD8(repeating: other)
    
    for cellID in entityTypes.indices {
      let compressed = mask.mask[cellID]
      let flags = CubicCell.flags & compressed
      
      var codes = entityTypes[cellID]
      let select = codes .!= 0
      codes.replace(with: newValue, where: flags .> 0 .& select)
      entityTypes[cellID] = codes
    }
  }
  
  var entities: [Entity] {
    var output: [Entity] = []
    let outputTransform = (
      SIMD3<Float>(squareSideLength, 0, 0),
      SIMD3<Float>(0, squareSideLength, 0),
      SIMD3<Float>(0, 0, squareSideLength)
    )
    for z in 0..<dimensions.z {
      for y in 0..<dimensions.y {
        for x in 0..<dimensions.x {
          let lowerCorner = SIMD3<Float>(SIMD3(x, y, z))
          var cellID = z * dimensions.y + y
          cellID = cellID * dimensions.x + x
          
          let cell = entityTypes[Int(cellID)]
          if all(cell .== .zero) {
            continue
          }
          for lane in 0..<8 where cell[lane] != 0 {
            let x = CubicCell.x0[lane] / 4
            let y = CubicCell.y0[lane] / 4
            let z = CubicCell.z0[lane] / 4
            let type = EntityType(compactRepresentation: cell[lane])
            
            var position = SIMD3<Float>(x, y, z)
            position += lowerCorner
            position =
            outputTransform.0 * position.x +
            outputTransform.1 * position.y +
            outputTransform.2 * position.z
            
            let entity = Entity(
              position: position, type: type)
            output.append(entity)
          }
        }
      }
    }
    return output
  }
}
