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
  
  // MARK: - Experiment 1
  //
  // lattice size = 8, H-C, 856x4413
  //
  // Version        | Total Time | Ratio / 17n^2 | Compare % | Sort % |
  // -------------- | ---------- | ------------- | --------- | ------ |
  // Original       |       3212 |           5.0 |       59% |    40% |
  // Optimization 1 |       1316 |           2.0 |       94% |     4% |
  //
  // Some investigation revealed that Float16 does not improve performance at
  // all. Not even to reduce register pressure in the ISA. Therefore, we will
  // drop the metric of equivalent FP16 instructions. We have >67% ALU
  // utilization at FP32, which is impressive. We are currently at 15 FP32
  // instructions per comparison. Future benchmarks should consist of:
  // - FP32 instructions per O(n^2)
  // - FP32 instructions per O(n) at scalar granularity (#succeeding matches)
  //
  // Proven lower bound: >10 scalar instructions per comparison
  // - subtract X, Y, Z
  // - FMA X^2, Y^2, Z^2
  // - add radii
  // - radius^2
  // - compare + select (2 instructions on CPU)
  
  // MARK: - Experiment 2
  //
  // We will now experiment O(nlog(n)) algorithms. Performance will be
  // benchmarked in FP32 instructions per O(n^2). The optimizations are as
  // follows:
  // - Optimization 1: previous experiment
  // - Optimization 2: removing Float16
  // - Optimization 3: condition 1 from the OpenMM algorithm, 8x8 blocks
  //   - Condition 2 and sorting provided inconsistent results when applied at
  //     this stage; defer to later. Sorting still needs to be multithreaded.
  // Optimization 3 |   959 |   269 |   282 | 12.4 | 21.2 | 15.0 |
  //  cond 2 + sort |                         25.2 | 34.9 | 28.4 |
  //   sorting only |                         17.8 | 30.5 | 24.5 |
  //   cond. 2 only |                         16.0 | 27.4 | 18.2 |
  // Optimization 3 | 58062 |  5397 | 11226 |  2.5 |  7.6 |  5.0 |
  //  cond. 2 only  |                          3.0 |  8.3 |  5.5 |
  //  cond 2 + sort |                          3.0 |  5.1 |  3.8 |
  //  sorting only  |                          2.0 |  4.1 |  2.9 |
  //   - Changing block size from 8x8 to 8x32 improved performance, but only
  //     for the largest problem sizes, with sorting enabled. There is not
  //     sufficient evidence to change the block size.
  // - Optimization 4: change how the LHS is preprocessed. Although it decreases
  //   performance for the largest problem sizes, it is a necessary step before
  //   O(nlogn) and O(n) algorithms can be attempted.
  // - Optimization 5: use a 2nd hierarchy level and always enable sorting.
  //   Don't remove the code for OpenMM condition 2 yet. A future optimization
  //   could remove both that and 1/2 the overhead for generating block bounds.
  // C-C       |
  // --------- |
  // size2=8   | 398.9 | 101.4 |  50.8 |  28.8 |  18.8 |  13.0 |  10.2 |   6.9 |   4.2
  // size2=16  | 387.8 | 114.2 |  46.8 |  25.9 |  17.0 |  11.9 |   8.2 |   5.0 |   1.8
  // size2=32  | 393.4 | 110.3 |  46.8 |  30.3 |  16.6 |  11.0 |   8.1 |   4.6 |   1.4
  // size2=64  | 487.5 | 111.6 |  49.3 |  25.8 |  16.7 |  11.3 |   7.8 |   4.3 |   1.3
  // size2=128 | 404.4 | 111.6 |  47.5 |  26.6 |  16.4 |  10.9 |   7.7 |   4.5 |   1.3
  // size2=256 | 448.8 | 116.1 |  48.4 |  25.6 |  17.2 |  10.8 |   8.0 |   4.6 |   1.4
  //
  // H-H       |
  // --------- |
  // size2=8   | 458.8 | 137.3 |  83.9 |  45.1 |  28.9 |  21.3 |  16.7 |  11.7 |   6.1
  // size2=16  | 415.5 | 172.8 |  79.2 |  45.1 |  27.5 |  19.5 |  14.4 |   9.4 |   3.9
  // size2=32  | 493.4 | 165.4 |  75.3 |  41.9 |  27.5 |  19.0 |  14.4 |   8.9 |   3.4
  // size2=64  | 467.5 | 155.1 |  79.6 |  42.7 |  27.4 |  19.1 |  14.0 |   8.7 |   3.4
  // size2=128 | 441.5 | 175.7 |  73.1 |  39.4 |  29.0 |  18.6 |  14.0 |   8.8 |   3.4
  // size2=256 | 528.0 | 159.5 |  75.3 |  41.9 |  29.0 |  18.4 |  14.1 |   8.8 |   3.7
  //
  // H-C       |
  // --------- |
  // size2=8   | 479.2 | 136.4 |  73.2 |  39.8 |  25.3 |  17.5 |  13.8 |   8.9 |   5.3
  // size2=16  | 416.7 | 146.8 |  67.2 |  35.1 |  23.0 |  15.6 |  11.7 |   7.0 |   2.8
  // size2=32  | 791.7 | 134.9 |  66.5 |  36.6 |  23.0 |  15.1 |  12.0 |   7.0 |   2.5
  // size2=64  | 541.7 | 143.8 |  68.4 |  35.4 |  23.4 |  16.0 |  11.9 |   6.9 |   2.5
  // size2=128 | 458.3 | 139.4 |  68.7 |  36.1 |  23.0 |  15.4 |  11.4 |   6.7 |   2.6
  // size2=256 | 552.1 | 137.9 |  72.5 |  38.5 |  22.9 |  15.3 |  11.7 |   7.0 |   2.8
  //
  // lattice size = 6
  // - C-C (1963x1963)
  // - H-H (796x796)
  // - H-C (496x1895)
  //
  // Version        | Total Time            | Ratio / n^2        |
  // -------------- | --------------------- | ------------------ |
  //                | C-C   | H-H   | H-C   | C-C  | H-H  | H-C  |
  // -------------- | --------------------- | ------------------ |
  // Optimization 1 |       |       |       |      |      |      |
  // Optimization 2 |  1514 |   311 |   370 | 19.6 | 24.5 | 19.7 |
  // Optimization 3 |   959 |   269 |   282 | 12.4 | 21.2 | 15.0 |
  // Optimization 4 |   998 |   270 |   286 | 12.9 | 21.3 | 15.2 |
  // Optimization 5 |  1296 |   358 |   418 | 16.8 | 28.3 | 22.2 |
  //
  // lattice size = 16
  // - C-C (34353x34353)
  // - H-H (5956x5956)
  // - H-C (3256x34353)
  //
  // Version        | Total Time            | Ratio / n^2        |
  // -------------- | --------------------- | ------------------ |
  //                | C-C   | H-H   | H-C   | C-C  | H-H  | H-C  |
  // -------------- | --------------------- | ------------------ |
  // Optimization 1 |       |       | 34307 |      |      | 15.4 |
  // Optimization 2 |  3e+5 | 11358 | 33417 | 15.0 | 16.0 | 14.9 |
  // Optimization 3 | 58062 |  5397 | 11226 |  2.5 |  7.6 |  5.0 |
  // Optimization 4 | 71341 |  5612 | 12235 |  3.0 |  7.9 |  5.5 |
  // Optimization 5 | 34044 |  2513 |  5658 |  1.4 |  3.5 |  2.5 |
  
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
    let latticeSizes: [Float] = [2, 3, 4, 5, 6, 7, 8, 10, 16]
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
    var tableFooter: [String] = []
    var maxColumnSizes: SIMD8<Int> = .zero
    
    for pass in 0..<2 {
      for measurement in data {
        var instructionLatency = Float(measurement[0]) * Float(measurement[1])
        instructionLatency /= 50e9 // CPU speed in FP32 instructions/second
        instructionLatency *= 1e6 // convert from s -> µs
        
        var expected = instructionLatency
        expected *= 10 // theoretical minimum cost in FP32 instructions
        let instructions = Float(measurement[2]) / Float(instructionLatency)
        
        let element1 = String(describing: measurement[0])
        let element2 = String(describing: measurement[1])
        let element3 = String(describing: measurement[2])
        let element4 = String(describing: Int(expected))
        let element5 = String(format: "%.1f", instructions)
        
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
          tableFooter.append(elements[4])
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
    
    output += "\n"
    for i in tableFooter.indices {
      output += tableFooter[i]
      if i < tableFooter.count - 1 {
        output += " | "
      }
    }
    return output
  }
}
