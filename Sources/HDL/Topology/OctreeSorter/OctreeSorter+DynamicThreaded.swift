//
//  OctreeSorter+DynamicThreaded.swift
//  HDL
//
//  Created by Philip Turner on 7/27/25.
//

import Dispatch

// TODO: Implement this as just a copy of single threaded for now. Get tests
// working with the new API call.
//
// compute ideal task count
//   retrieve total atom count
//   retrieve levels remaining (7 @ 4 nm)
//   compute total latency from 7.5 ns/atom/level
//   ideal task count = round_to_nearest(total latency / 20 μs)
//   restrict ideal task count to 1 to 8
//
// early returns to disallow work splitting
//   level size is 1.0 nm or smaller
//   task count is 1
//
// child count = 1 to 8
// task count ≥ child count
//   every child
//   likely at highest level of tree
//   likely leaving breadcrumbs
// task count < child count
//   continue with algorithm

extension OctreeSorter {
  // Algorithm that adaptively uses multi-threading, when a subset of the
  // octree has enough atoms.
  func mortonReorderingDynamic() -> [UInt32] {
    return []
  }
}
