import XCTest
import HDL
import Numerics

// This file contains a performance test for topology.match(), which differs
// significantly from the other performance tests. It was moved into a separate
// file for better organization.

final class MatchTests: XCTestCase {
  static let printPerformanceSummary = false
  
  // We need to run performance tests of Topology.match, to ensure the
  // acceleration algorithm is working properly. One could imagine subtle bugs
  // that make it incorrect, resulting in O(n^2) scaling.
  
  // MARK: - Theory
  //
  // Proven lower bound: >11 scalar instructions per comparison
  // - (FP32) subtract X, Y, Z
  // - (FP16 or FP32) FMA X^2, Y^2, Z^2
  // - (FP16) add radii
  // - (FP16) radius^2
  // - (FP16) compare + select (2 instructions on CPU)
  // - (FP16) O(log(n)) reduction across vector lanes
  //
  // Latency on single M1 CPU core:
  // - 100 GFLOPS FP32 -> 50 GIPS FP32, 100 GIPS FP16
  // - effectively >14 instructions at FP16 bitwidth
  //
  // -     100 atoms ---        >1.4 µs
  // -   1,000 atoms ---      >140 µs
  // -  10,000 atoms ---    >14 ms
  // - 100,000 atoms --- >1.4 s
  //
  // What is the ratio of actual to theoretical minimum latency? Apply this rule
  // to different problem sizes and shapes, then report a percentage.
  // - The original implementation kept everything in FP32 before converting the
  //   delta squared to FP16. Compare to the number 1.7 instead of 1.4. Then,
  //   see whether switching to FP16 provides the expected speedup. This is a
  //   critical test of how well the theory matches reality.
  
  // MARK: - Experiment
  //
  // lattice size = 8, H-C, 856x4413
  //
  // Version        | Total Time | Ratio / 17n^2 | Compare % | Sort % |
  // -------------- | ---------- | ------------- | --------- | ------ |
  // Original       |       3212 |           5.0 |       59% |    40% |
  // Optimization 1 |       1316 |           2.0 |       94% |     4% |
  
  func testMatch() {
    // Accumulate statistics and sort by workload (size of a square representing
    // the number of comparisons). Also, report statistics for each problem
    // shape separately.
    var summary = MatchSummary()
    defer {
      print()
      print("Match Performance Report")
      print()
      print(summary.createReport(header: "C-C"))
      print()
      print(summary.createReport(header: "H-H"))
      print()
      print(summary.createReport(header: "H-C"))
      print()
    }
    
    // Compute cost scales with the sixth power of lattice width.
    #if RELEASE
    let latticeSizes: [Float] = [2, 3, 4, 5, 6, 7, 8]
    #else
    let latticeSizes: [Float] = [2, 3]
    #endif
    
    for latticeSize in latticeSizes {
      let lattice = Lattice<Cubic> { h, k, l in
        Bounds { latticeSize * (h + k + l) }
        Material { .elemental(.carbon) }
      }
      
      var topology = Topology()
      topology.atoms = lattice.atoms
      
      // Match carbons against carbons.
      do {
        let start = cross_platform_media_time()
        let matches = topology.match(topology.atoms)
        let end = cross_platform_media_time()
        summary.ccMetrics.append(
          SIMD3(
            topology.atoms.count,
            topology.atoms.count,
            Int((end - start) * 1e6)))
        
        var removedAtoms: [UInt32] = []
        var insertedBonds: [SIMD2<UInt32>] = []
        for i in topology.atoms.indices {
          let matchRange = matches[i]
          if matchRange.count <= 2 {
            removedAtoms.append(UInt32(i))
          }
          for j in matchRange where i < j {
            insertedBonds.append(SIMD2(UInt32(i), UInt32(j)))
          }
        }
        topology.insert(bonds: insertedBonds)
        topology.remove(atoms: removedAtoms)
      }
      
      // Define the hydrogen placement distance and the search radius for a
      // subsequent step.
      var ccBondLength = Constant(.square) { .elemental(.carbon) }
      ccBondLength *= Float(3).squareRoot() / 4
      
      // Generate hydrogens and de-duplicate them.
      var hydrogenTopology = Topology()
      do {
        let orbitals = topology.nonbondingOrbitals()
        
        var insertedAtoms: [Entity] = []
        for i in topology.atoms.indices {
          let atom = topology.atoms[i]
          for orbital in orbitals[i] {
            let position = atom.position + orbital * ccBondLength
            let hydrogen = Entity(position: position, type: .atom(.hydrogen))
            insertedAtoms.append(hydrogen)
          }
        }
        hydrogenTopology.insert(atoms: insertedAtoms)
        
        let start = cross_platform_media_time()
        let matches = hydrogenTopology.match(
          hydrogenTopology.atoms, algorithm: .absoluteRadius(0.020))
        let end = cross_platform_media_time()
        summary.hhMetrics.append(
          SIMD3(
            hydrogenTopology.atoms.count,
            hydrogenTopology.atoms.count,
            Int((end - start) * 1e6)))
        
        var removedAtoms: [UInt32] = []
        for i in hydrogenTopology.atoms.indices {
          let matchRange = matches[i]
          var survived = true
          for match in matchRange where match < i {
            survived = false
          }
          if !survived {
            removedAtoms.append(UInt32(i))
          }
        }
        hydrogenTopology.remove(atoms: removedAtoms)
      }
      
      // Locate the carbons attached to each hydrogen collision site.
      do {
        let start = cross_platform_media_time()
        _ = topology.match(
          hydrogenTopology.atoms, algorithm: .absoluteRadius(1.1 * ccBondLength))
        let end = cross_platform_media_time()
        summary.hcMetrics.append(
          SIMD3(
            hydrogenTopology.atoms.count,
            topology.atoms.count,
            Int((end - start) * 1e6)))
      }
    }
  }
}

private struct MatchSummary {
  var ccMetrics: [SIMD3<Int>] = []
  var hhMetrics: [SIMD3<Int>] = []
  var hcMetrics: [SIMD3<Int>] = []
  
  init() {
    
  }
  
  func createReport(header: String) -> String {
    var data: [SIMD3<Int>]
    switch header {
    case "C-C": data = ccMetrics
    case "H-H": data = hhMetrics
    case "H-C": data = hcMetrics
    default: fatalError("Unrecognized header: \(header)")
    }
    
    // MARK: - Gather Column Sizes
    
    var tableRows: [String] = []
    var maxColumnSizes: SIMD8<Int> = .zero
    
    for pass in 0..<2 {
      for measurement in data {
        let element1 = String(describing: measurement[0])
        let element2 = String(describing: measurement[1])
        let element3 = String(describing: measurement[2])
        
        var expected = Float(measurement[0]) * Float(measurement[1])
        expected *= 17 // cost in equivalent FP16 instructions
        expected /= 100e9 // CPU speed in FP16 instructions/second
        expected *= 1e6 // convert from s -> µs
        let element4 = String(describing: Int(expected))
        
        let ratio = Float(measurement[2]) / Float(expected)
        let element5 = String(format: "%.1f", ratio) + "x"
        
        var elements = [element1, element2, element3, element4, element5]
        if pass == 0 {
          var columnSizes: SIMD8<Int> = .zero
          for i in elements.indices {
            columnSizes[i] = elements[i].count
          }
          maxColumnSizes.replace(
            with: columnSizes, where: columnSizes .> maxColumnSizes)
        } else {
          for i in elements.indices {
            var element = elements[i]
            while element.count < maxColumnSizes[i] {
              element = " " + element
            }
            elements[i] = element
          }
          
          tableRows.append("""
            \(elements[0]) x \(elements[1]) | \
            \(elements[2]) | \(elements[3]) | \
            \(elements[4])
            """)
        }
      }
    }
    
    // MARK: - Create Table
    
    var dividerSections: [String] = []
    for lane in 0..<5 {
      let size = maxColumnSizes[lane]
      let divider = String(repeating: "-", count: size)
      dividerSections.append(divider)
    }
    
    var dividerRow = ""
    dividerRow += "\(dividerSections[0])---\(dividerSections[1])"
    dividerRow += " | "
    dividerRow += "\(dividerSections[2]) | \(dividerSections[3])"
    dividerRow += " | "
    dividerRow += "\(dividerSections[4])"
    
    var headerRow = ""
    let tableWidth = dividerRow.count
    do {
      let padding = UInt(tableWidth - header.count)
      let left = padding / 2
      headerRow += String(repeating: " ", count: Int(left))
      headerRow += header
      headerRow += String(repeating: " ", count: Int(padding - left))
    }
    
    var output = """
    \(headerRow)
    \(dividerRow)
    """
    for row in tableRows {
      output += "\n" + row
    }
    return output
  }
}
