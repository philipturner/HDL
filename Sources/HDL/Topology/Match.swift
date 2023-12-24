//
//  Match.swift
//
//
//  Created by Philip Turner on 12/23/23.
//

// TODO: Start with an O(n^2) algorithm for now. Get all the tests finished,
// enabling a full-stack workflow. Later, implement the optimization that
// creates a grid for O(n) traversal. This approach allows an entire test suite
// to be created earlier on.
//
// Brilliant! The OpenMM approach without a grid, that just has a
// small-prefactor O(n^2) search.
// - Still use a grid for Morton reordering, there's no neighbor search there.
// - The block-wise approach is the most inherently parallel one.

// TODO: Make a general function that matches a group of SIMD4<Float> against
// another group of SIMD4<Float>, where the fourth vector index is the radius.
// It will be called at each level of recursion. Represent lower levels of atom
// blocks with a center point and radius for the bounding sphere.
