//
//  Map.swift
//  HDLTests
//
//  Created by Philip Turner on 1/2/24.
//

import Atomics
import Dispatch

extension Topology {
  public enum MapNode: Sendable {
    case atoms
    case bonds
  }
  
  public struct MapStorage: Collection, Equatable, Sendable {
    @usableFromInline var storage: SIMD8<Int32>
    
    public typealias Index = Int
    
    public typealias Element = UInt32
    
    @_transparent
    public var startIndex: Int { 0 }
    
    @_transparent
    public var endIndex: Int {
      if storage[7] & 0x7FFF_FFF == storage[7] {
        return Int(storage[7])
      } else {
        return 8
      }
    }
    
    @_transparent
    public func index(after i: Int) -> Int {
      i &+ 1
    }
    
    public subscript(position: Int) -> UInt32 {
      @_transparent
      _read {
         yield UInt32(truncatingIfNeeded: storage[position] & 0x7FFF_FFFF)
      }
    }
  }
  
  public func map(
    _ sourceNode: MapNode,
    to targetNode: MapNode
  ) -> [MapStorage] {
    switch (sourceNode, targetNode) {
    case (.atoms, _):
      let connectionsMap = createConnectionsMap(targetNode: targetNode)
      return unsafeBitCast(connectionsMap, to: [_].self)
    case (.bonds, .atoms):
      var outputStorage: [MapStorage] = []
      outputStorage.reserveCapacity(bonds.count)
      for i in bonds.indices {
        let bond = bonds[i]
        var vector = SIMD8<Int32>.zero
        vector[0] = Int32(truncatingIfNeeded: bond[0])
        vector[1] = Int32(truncatingIfNeeded: bond[1])
        vector[7] = 2
        outputStorage.append(MapStorage(storage: vector))
      }
      return outputStorage
    case (.bonds, .bonds):
      fatalError("Bonds to bonds map is not supported.")
    }
  }
}

extension Topology {
  private func createConnectionsMap(
    targetNode: MapNode
  ) -> [SIMD8<Int32>] {
    var connectionsMap = [SIMD8<Int32>](
      repeating: .init(repeating: -1), count: atoms.count)
    guard atoms.count > 0 else {
      return []
    }
    
    // Optimization: use atomics with 16 bits instead of 32 bits
    nonisolated(unsafe)
    let atomicPointer: UnsafeMutablePointer<Int16.AtomicRepresentation> =
      .allocate(capacity: atoms.count)
    let atomicOpaque = OpaquePointer(atomicPointer)
    let atomicCasted = UnsafeMutablePointer<Int16>(atomicOpaque)
    atomicCasted.initialize(repeating: .zero, count: atoms.count)
    
    connectionsMap.withUnsafeMutableBufferPointer {
      let connectionsPointer = $0.baseAddress.unsafelyUnwrapped
      let connectionsOpaque = OpaquePointer(connectionsPointer)
      nonisolated(unsafe)
      let connectionsCasted = UnsafeMutablePointer<Int32>(connectionsOpaque)
      
      bonds.withUnsafeBufferPointer {
        let bondsPointer = $0.baseAddress.unsafelyUnwrapped
        let bondsOpaque = OpaquePointer(bondsPointer)
        nonisolated(unsafe)
        let bondsCasted = UnsafePointer<UInt32>(bondsOpaque)
        
        let taskSize: Int = 5_000
        let safeBondCount = bonds.count
        
        @Sendable
        func execute(taskID: Int) {
          let scalarStart = taskID &* taskSize
          let scalarEnd = min(scalarStart &+ taskSize, 2 * safeBondCount)
          
          for i in scalarStart..<scalarEnd {
            let atomID = Int(truncatingIfNeeded: bondsCasted[i])
            var idToWrite: Int32
            
            if targetNode == .atoms {
              var otherID: UInt32
              if i % 2 == 0 {
                otherID = bondsCasted[i &+ 1]
              } else {
                otherID = bondsCasted[i &- 1]
              }
              
              idToWrite = Int32(truncatingIfNeeded: otherID)
            } else {
              idToWrite = Int32(truncatingIfNeeded: i / 2)
            }
            
            let pointer = atomicPointer.advanced(by: atomID)
            let atomic = UnsafeAtomic<Int16>(at: pointer)
            let lane = atomic.loadThenWrappingIncrement(ordering: .relaxed)
            let lane64 = Int(truncatingIfNeeded: lane)
            
            if lane < 8 {
              connectionsCasted[atomID &* 8 &+ lane64] = Int32(
                truncatingIfNeeded: idToWrite)
            }
          }
        }
        
        // TODO: Unit test how the compiler behaves when it receives an empty
        // array, without adding any special checks/early returns for edge cases.
        let taskCount = (2 * bonds.count + taskSize - 1) / taskSize
        if taskCount <= 1 {
          for taskID in 0..<taskCount {
            execute(taskID: taskID)
          }
        } else {
          DispatchQueue.concurrentPerform(iterations: taskCount) { z in
            execute(taskID: z)
          }
        }
      }
    }
    
    for i in atoms.indices {
      let count = Int32(truncatingIfNeeded: atomicCasted[i])
      if count < 8 {
        connectionsMap[i][7] = count
      } else {
        connectionsMap[i][7] |= .init(truncatingIfNeeded: 0x8000_0000)
      }
    }
    
    atomicPointer.deallocate()
    return connectionsMap
  }
}
