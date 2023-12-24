//
//  Orbitals.swift
//  
//
//  Created by Philip Turner on 12/24/23.
//

#if false
extension Topology {
  public enum MatchAlgorithm {
    // Search for neighbors within a fixed radius (in nanometers).
    case absoluteRadius(Float)
    
    // Search for neighbors within a multiple of covalent bond length.
    case covalentBondLength(Float)
  }
  
  public func match(
    _ input: [Entity],
    algorithm: MatchAlgorithm = .covalentBondLength(1.5)
  ) -> [ArraySlice<UInt32>] {
    
  }
}
#endif

extension Topology {
  public enum OrbitalHybridization {
    case sp1
    case sp2
    case sp3
  }
  
  public func nonbondingOrbitals(
    hybridization: OrbitalHybridization = .sp3
  ) -> [ArraySlice<SIMD3<Float>>] {
    fatalError("Not implemented.")
  }
}
