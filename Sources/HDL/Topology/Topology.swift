//
//  Topology.swift
//  MolecularRenderer
//
//  Created by Philip Turner on 12/2/23.
//

// TODO: Rewrite large sections of this codebase from scratch. Don't update the
// API documentation until the final result is completed. The grid size and
// search radius will be determined using some heuristics of atom covalent
// radius and population statistics of the target grid. The end product will be
// much more low-level than 'Diamondoid', but permits exotic use cases like
// (100) reconstruction, passivating strained shell structures, and Eric's
// bond formation algorithm. But without baking those algorithms into the
// compiler. Also, the '.bond(_)' entity needs to be removed. Topology may
// output a different IR than 'Entity/EntityType', which was simply a solution
// to match the API design language of the DSL.

// Always store a source of truth as a linear array in memory. The sorting into
// grid cells is transient, and frequently reconstructed (e.g. some added atoms
// may require the grid to expand its dimensions).

public struct Topology {
  public init(_ entities: [Entity]) {
    // Start by creating the stored properties of Topology, and the basic
    // computed properties. Add the functions for inserting and removing atoms
    // and bonds.
  }
}
