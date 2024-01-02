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
  // - Optimization 6: added a 3rd hierarchy level, which required some major
  //   code changes. I have proven a small but measurable speedup for extremely
  //   large problem sizes. The code restructuring slightly slowed down smaller
  //   problems. For the largest problem sizes, we clearly have O(n) scaling.
  //   Going higher than 128x128 seems unnecessary.
  //
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
  // Optimization 6 |  1411 |   384 |   449 | 18.3 | 30.3 | 23.9 |
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
  // Optimization 6 | 32265 |  2476 |  6012 |  1.4 |  3.5 |  2.7 |
  
  // MARK: - Experiment 3
  //
  // This is the final round of optimization. Switch to an O(n) metric:
  // microseconds per RMS atom. This round will include multithreading and
  // kernel ensembles for different problem sizes.
  // - Optimization 7: change block bounds to a more efficient representation
  // - Optimization 8: use 32-bit integers instead of 64-bit integers
  // - Optimization 9: 3-way multithread the preparation stage
  // - Optimization 10: use 4 x SIMD8<Float> as the transformed format
  // - Optimization 11: fully optimize the preparation stage
  // - Optimization 12: multithread the searching stage
  // - Optimization 13: fuse Morton order mapping with preparation
  // - Optimization 14: remove a memory or object allocation bottleneck
  // - Optimization 15: fuse several post-processing stages into one pass
  // - Optimization 16: parallelize a tiny portion of pre-processing
  // - Optimization 17: disable sorting for the smallest problem sizes
  //
  // In the final state of optimization, asymptotically small problems only have
  // ~5 microseconds of latency. It is not clear exactly how such small latency
  // is being achieved with 'DispatchQueue.concurrentPerform', but it is a good
  // thing. The issue of O(1) overhead will not be investigated further.
  //
  // lattice size = 3
  //
  // Version         | Total Time               | Ratio / n^2           |
  // --------------- | ------------------------ | --------------------- |
  //                 | C-C    | H-H    | H-C    | C-C   | H-H   | H-C   |
  // --------------- | ------------------------ | --------------------- |
  // Optimization 6  |    178 |    115 |    115 | 0.64  | 0.62  | 0.63  |
  // Optimization 7  |    188 |    114 |     94 | 0.67  | 0.62  | 0.51  |
  // Optimization 8  |    175 |    111 |     92 | 0.62  | 0.60  | 0.50  |
  // Optimization 9  |    180 |    121 |     98 | 0.643 | 0.658 | 0.534 |
  // Optimization 10 |    172 |    114 |    100 | 0.614 | 0.620 | 0.545 |
  // Optimization 11 |    184 |    110 |     96 | 0.657 | 0.598 | 0.523 |
  // Optimization 12 |    142 |    112 |    100 | 0.507 | 0.609 | 0.545 |
  // Optimization 13 |    131 |     91 |     90 | 0.468 | 0.495 | 0.490 |
  // Optimization 14 |     97 |     70 |     74 | 0.346 | 0.380 | 0.403 |
  // Optimization 15 |     87 |     61 |     67 | 0.311 | 0.332 | 0.365 |
  // Optimization 16 |     86 |     62 |     64 | 0.307 | 0.337 | 0.348 |
  // Optimization 17 |     64 |     34 |     37 | 0.229 | 0.185 | 0.201 |
  //
  // lattice size = 6
  //
  // Version         | Total Time               | Ratio / n^2           |
  // --------------- | ------------------------ | --------------------- |
  //                 | C-C    | H-H    | H-C    | C-C   | H-H   | H-C   |
  // --------------- | ------------------------ | --------------------- |
  // Optimization 6  |   1505 |    362 |    474 | 0.77  | 0.45  | 0.49  |
  // Optimization 7  |   1397 |    381 |    467 | 0.71  | 0.48  | 0.48  |
  // Optimization 8  |   1394 |    363 |    466 | 0.71  | 0.46  | 0.48  |
  // Optimization 9  |   1287 |    359 |    435 | 0.656 | 0.451 | 0.449 |
  // Optimization 10 |   1277 |    354 |    449 | 0.651 | 0.445 | 0.463 |
  // Optimization 11 |   1252 |    349 |    464 | 0.638 | 0.438 | 0.479 |
  // Optimization 12 |    699 |    295 |    383 | 0.356 | 0.371 | 0.395 |
  // Optimization 13 |    712 |    298 |    366 | 0.363 | 0.374 | 0.378 |
  // Optimization 14 |    457 |    212 |    314 | 0.233 | 0.266 | 0.324 |
  // Optimization 15 |    383 |    190 |    293 | 0.195 | 0.239 | 0.302 |
  // Optimization 16 |    356 |    183 |    298 | 0.181 | 0.230 | 0.307 |
  // Optimization 17 |    248 |    131 |    185 | 0.126 | 0.165 | 0.191 |
  //
  // lattice size = 24
  //
  // Version         | Total Time               | Ratio / n^2           |
  // --------------- | ------------------------ | --------------------- |
  //                 | C-C    | H-H    | H-C    | C-C   | H-H   | H-C   |
  // --------------- | ------------------------ | --------------------- |
  // Optimization 6  | 128530 |   5680 |  18697 | 1.13  | 0.42  | 0.65  |
  // Optimization 7  | 113690 |   5379 |  17949 | 1.00  | 0.40  | 0.63  |
  // Optimization 8  | 109897 |   5121 |  17352 | 0.96  | 0.38  | 0.61  |
  // Optimization 9  | 105642 |   4637 |  17067 | 0.926 | 0.342 | 0.596 |
  // Optimization 10 | 102438 |   4503 |  16733 | 0.898 | 0.333 | 0.585 |
  // Optimization 11 | 101426 |   4566 |  16673 | 0.889 | 0.337 | 0.583 |
  // Optimization 12 |  36676 |   3037 |  10749 | 0.321 | 0.224 | 0.376 |
  // Optimization 13 |  38097 |   3333 |  10848 | 0.334 | 0.246 | 0.379 |
  // Optimization 14 |  24622 |   2171 |  10113 | 0.216 | 0.160 | 0.353 |
  // Optimization 15 |  18380 |   1709 |   9784 | 0.161 | 0.126 | 0.342 |
  // Optimization 16 |  17869 |   1621 |   9247 | 0.157 | 0.120 | 0.323 |
  // Optimization 17 |  17830 |   1715 |   9386 | 0.156 | 0.127 | 0.328 |
  
  func testMatch() {
    // Accumulate statistics and sort by workload (size of a square representing
    // the number of comparisons). Also, report statistics for each problem
    // shape separately.
    var summary = MatchSummary()
    defer {
      if Self.printPerformanceSummary {
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
    }
    
    // Compute cost scales with the sixth power of lattice width.
    #if RELEASE
    let latticeSizes: [Float] = [1, 1, 2, 3, 4, 5, 6, 7, 8, 10, 16, 24]
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
          hydrogenTopology.atoms,
          algorithm: .absoluteRadius(1.1 * ccBondLength))
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
        var element5 = String(format: "%.1f", instructions)
        
        // Revise so the output is µs/atom.
        var rmsAtoms = Float(measurement[0]) * Float(measurement[1])
        rmsAtoms.formSquareRoot()
        let µsPerAtom = Float(measurement[2]) / rmsAtoms
        element5 = String(format: "%.3f", µsPerAtom)
        
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
    
//    output += "\n"
//    for i in tableFooter.indices {
//      output += tableFooter[i]
//      if i < tableFooter.count - 1 {
//        output += " | "
//      }
//    }
    return output
  }
}

// MARK: - Raw Data

// With the match stage on single-core:
//
// atoms: 280 x 280
// total | sort | prepare | match | sort | map
// ----- | ---- | ------- | ----- | ---- | ---
//  100% |  31% |     16% |   35% |   9% |  9%
//   159 |   49 |      26 |    56 |   14 |  15
//
// atoms: 1963 x 1963
// total | sort | prepare | match | sort | map
// ----- | ---- | ------- | ----- | ---- | ---
//  100% |  12% |     10% |   62% |   8% |  8%
//  1237 |  153 |     120 |   765 |   96 | 101
//
// atoms: 114121 x 114121
// total | sort | prepare | match | sort |  map
// ----- | ---- | ------- | ----- | ---- | ----
//  100% |   7% |      5% |   76% |   6% |   7%
// 99795 | 6581 |    5180 | 75581 | 5881 | 6572

// With task=128 multithreading:
//
// atoms: 280 x 280
// total | sort | prepare | match | sort | map
// ----- | ---- | ------- | ----- | ---- | ---
//  100% |  38% |     16% |   22% |  12% | 11%
//   155 |   59 |      25 |    34 |   19 |  17
//
// atoms: 1963 x 1963
// total | sort | prepare | match | sort | map
// ----- | ---- | ------- | ----- | ---- | ---
//  100% |  23% |     18% |   26% |  16% | 16%
//   696 |  160 |     127 |   184 |  115 | 110
//
// atoms: 114121 x 114121
// total | sort | prepare | match | sort |  map
// ----- | ---- | ------- | ----- | ---- | ----
//  100% |  17% |     15% |   31% |  17% |  19%
// 36354 | 6325 |    5421 | 11406 | 6354 | 6849

// After several additional optimizations:
//
// atoms: 280 x 280
// total | sort | prepare | match | sort | map
// ----- | ---- | ------- | ----- | ---- | ---
//  100% |  68% |      3% |   29% |   0% |  0%
//   109 |   74 |       3 |    32 |    0 |   0
//
// atoms: 1963 x 1963
// total | sort | prepare | match | sort | map
// ----- | ---- | ------- | ----- | ---- | ---
//  100% |  59% |      1% |   40% |   0% |  0%
//   385 |  229 |       2 |   153 |    0 |   0
//
// atoms: 114121 x 114121
// total | sort | prepare | match | sort | map
// ----- | ---- | ------- | ----- | ---- | ---
//  100% |  51% |      0% |   48% |   1% |  0%
// 18399 | 9353 |       4 |  8889 |  152 |   0

// Final state before optimization is finished:
//
// atoms: 280 x 280
// total | sort | prepare | match | sort | map
// ----- | ---- | ------- | ----- | ---- | ---
//  100% |  43% |      2% |   54% |   1% |  0%
//    59 |   25 |       1 |    32 |    1 |   0
//
// atoms: 1963 x 1963
// total | sort | prepare | match | sort | map
// ----- | ---- | ------- | ----- | ---- | ---
//  100% |  30% |      2% |   68% |   0% |  0%
//   239 |   73 |       4 |   162 |    1 |   0
//
// atoms: 114121 x 114121
// total | sort | prepare | match | sort | map
// ----- | ---- | ------- | ----- | ---- | ---
//  100% |  49% |      0% |   50% |   1% |  0%
// 17747 | 8729 |       7 |  8861 |  149 |   1
//
//                  C-C
// --------------- | ----- | ------- | -----
//     18 x     18 |    19 |       0 | 1.056
//     18 x     18 |     7 |       0 | 0.389
//     95 x     95 |    20 |       1 | 0.211
//    280 x    280 |    66 |      15 | 0.236
//    621 x    621 |   110 |      77 | 0.177
//   1166 x   1166 |   173 |     271 | 0.148
//   1963 x   1963 |   241 |     770 | 0.123
//   3060 x   3060 |   369 |    1872 | 0.121
//   4505 x   4505 |   536 |    4059 | 0.119
//   8631 x   8631 |  1098 |   14898 | 0.127
//  34353 x  34353 |  5282 |  236025 | 0.154
// 114121 x 114121 | 17709 | 2604720 | 0.155
//
//                 H-H
// ------------- | ---- | ----- | -----
//    16 x    16 |    7 |     0 | 0.438
//    16 x    16 |    3 |     0 | 0.188
//    76 x    76 |   14 |     1 | 0.184
//   184 x   184 |   41 |     6 | 0.223
//   340 x   340 |   76 |    23 | 0.224
//   544 x   544 |  100 |    59 | 0.184
//   796 x   796 |  136 |   126 | 0.171
//  1096 x  1096 |  177 |   240 | 0.161
//  1444 x  1444 |  207 |   417 | 0.143
//  2284 x  2284 |  296 |  1043 | 0.130
//  5956 x  5956 |  972 |  7094 | 0.163
// 13540 x 13540 | 1630 | 36666 | 0.120
//
//               H-C
// ------------- | ---- | ------ | -----
//   16 x     10 |    4 |      0 | 0.316
//   16 x     10 |   14 |      0 | 1.107
//   64 x     75 |   14 |      0 | 0.202
//  136 x    248 |   46 |      6 | 0.250
//  232 x    577 |   67 |     26 | 0.183
//  352 x   1110 |  131 |     78 | 0.210
//  496 x   1895 |  200 |    187 | 0.206
//  664 x   2980 |  259 |    395 | 0.184
//  856 x   4413 |  327 |    755 | 0.168
// 1312 x   8515 |  519 |   2234 | 0.155
// 3256 x  34165 | 3070 |  22248 | 0.291
// 7192 x 113837 | 9510 | 163743 | 0.332

/*
 Here is code from the Swift Standard Library, explaining the internal layout of
 ArraySlice. MemoryLayout reports 32 bytes, which matches the data structure
 shown here.
 
 /// Buffer type for `ArraySlice<Element>`.
 @frozen
 @usableFromInline
 internal struct _SliceBuffer<Element>
   : _ArrayBufferProtocol,
     RandomAccessCollection
 {
   #if $Embedded
   @usableFromInline
   typealias AnyObject = Builtin.NativeObject
   #endif

   internal typealias NativeStorage = _ContiguousArrayStorage<Element>
   @usableFromInline
   internal typealias NativeBuffer = _ContiguousArrayBuffer<Element>

   /// An object that keeps the elements stored in this buffer alive.
   @usableFromInline
   internal var owner: AnyObject

   @usableFromInline
   internal let subscriptBaseAddress: UnsafeMutablePointer<Element>

   /// The position of the first element in a non-empty collection.
   ///
   /// In an empty collection, `startIndex == endIndex`.
   @usableFromInline
   internal var startIndex: Int

   /// [63:1: 63-bit index][0: has a native buffer]
   @usableFromInline
   internal var endIndexAndFlags: UInt
 */

// C-C       |
// --------- |
// task=128  | 0.716 | 0.543 | 0.440 | 0.417 | 0.361 | 0.318 | 0.321 | 0.314 | 0.315 | 0.318
// task=256  | 0.695 | 0.600 | 0.517 | 0.388 | 0.368 | 0.335 | 0.320 | 0.313 | 0.310 | 0.317
// task=512  | 0.726 | 0.586 | 0.572 | 0.474 | 0.388 | 0.346 | 0.332 | 0.311 | 0.322 | 0.324
// task=1024 | 0.705 | 0.586 | 0.625 | 0.569 | 0.471 | 0.391 | 0.365 | 0.337 | 0.314 | 0.318
// task=2048 | 0.695 | 0.582 | 0.605 | 0.607 | 0.648 | 0.534 | 0.452 | 0.361 | 0.317 | 0.318
// task=4096 | 0.747 | 0.604 | 0.614 | 0.615 | 0.643 | 0.671 | 0.645 | 0.477 | 0.329 | 0.331
//
// H-H       |
// --------- |
// task=128  | 0.671 | 0.603 | 0.512 | 0.434 | 0.394 | 0.330 | 0.333 | 0.279 | 0.236 | 0.224
// task=256  | 0.684 | 0.609 | 0.485 | 0.421 | 0.391 | 0.332 | 0.303 | 0.283 | 0.238 | 0.221
// task=512  | 0.658 | 0.636 | 0.515 | 0.458 | 0.388 | 0.356 | 0.319 | 0.285 | 0.260 | 0.223
// task=1024 | 0.592 | 0.598 | 0.500 | 0.441 | 0.425 | 0.378 | 0.342 | 0.313 | 0.245 | 0.222
// task=2048 | 0.618 | 0.592 | 0.521 | 0.436 | 0.412 | 0.389 | 0.386 | 0.361 | 0.259 | 0.229
// task=4096 | 0.592 | 0.500 | 0.485 | 0.445 | 0.428 | 0.387 | 0.371 | 0.373 | 0.311 | 0.243
//
// H-C       |
// --------- |
// task=128  | 0.664 | 0.446 | 0.528 | 0.429 | 0.392 | 0.331 | 0.323 | 0.308 | 0.318 | 0.372
// task=256  | 0.606 | 0.528 | 0.528 | 0.440 | 0.397 | 0.357 | 0.356 | 0.312 | 0.337 | 0.379
// task=512  | 0.606 | 0.528 | 0.547 | 0.474 | 0.456 | 0.429 | 0.398 | 0.351 | 0.336 | 0.389
// task=1024 | 0.577 | 0.555 | 0.549 | 0.469 | 0.461 | 0.459 | 0.456 | 0.437 | 0.376 | 0.390
// task=2048 | 0.635 | 0.528 | 0.528 | 0.467 | 0.458 | 0.457 | 0.464 | 0.486 | 0.451 | 0.428
// task=4096 | 0.664 | 0.550 | 0.533 | 0.480 | 0.485 | 0.465 | 0.479 | 0.476 | 0.530 | 0.483

// cutoff for whether to sort the atoms:
//
// C-C         |
// ----------- |
// cutoff=500  | 0.379 | 0.239 | 0.306 | 0.220 | 0.201 | 0.185 | 0.171 | 0.153 | 0.148 | 0.158
// cutoff=1000 | 0.274 | 0.289 | 0.180 | 0.230 | 0.201 | 0.180 | 0.161 | 0.158 | 0.147 | 0.160
// cutoff=2000 | 0.347 | 0.225 | 0.167 | 0.121 | 0.121 | 0.194 | 0.181 | 0.166 | 0.156 | 0.156
// cutoff=4000 | 0.368 | 0.264 | 0.174 | 0.160 | 0.125 | 0.125 | 0.168 | 0.168 | 0.146 | 0.157
// cutoff=8000 | 0.337 | 0.246 | 0.172 | 0.154 | 0.122 | 0.140 | 0.115 | 0.159 | 0.144 | 0.158
// cutoff=50k  | 0.295 | 0.243 | 0.167 | 0.161 | 0.120 | 0.119 | 0.112 | 0.133 | 0.185 | 0.157
// cutoff=150k | 0.326 | 0.239 | 0.180 | 0.149 | 0.126 | 0.224 | 0.114 | 0.126 | 0.183 | 0.318
//
// H-H         |
// ----------- |
// cutoff=500  | 0.197 | 0.174 | 0.171 | 0.298 | 0.259 | 0.210 | 0.212 | 0.173 | 0.134 | 0.136
// cutoff=1000 | 0.184 | 0.196 | 0.188 | 0.182 | 0.157 | 0.228 | 0.197 | 0.185 | 0.141 | 0.124
// cutoff=2000 | 0.184 | 0.168 | 0.191 | 0.178 | 0.157 | 0.131 | 0.131 | 0.170 | 0.144 | 0.125
// cutoff=4000 | 0.263 | 0.185 | 0.165 | 0.197 | 0.173 | 0.145 | 0.137 | 0.126 | 0.140 | 0.134
// cutoff=8000 | 0.237 | 0.196 | 0.171 | 0.193 | 0.177 | 0.156 | 0.163 | 0.128 | 0.149 | 0.128
// cutoff=50k  | 0.197 | 0.207 | 0.179 | 0.165 | 0.143 | 0.152 | 0.143 | 0.116 | 0.152 | 0.218
// cutoff=150k | 0.171 | 0.217 | 0.174 | 0.197 | 0.168 | 0.152 | 0.150 | 0.123 | 0.160 | 0.233
//
// H-C         |
// ----------- |
// cutoff=500  | 0.202 | 0.201 | 0.358 | 0.363 | 0.335 | 0.297 | 0.289 | 0.265 | 0.283 | 0.327
// cutoff=1000 | 0.159 | 0.191 | 0.153 | 0.389 | 0.343 | 0.338 | 0.314 | 0.285 | 0.280 | 0.332
// cutoff=2000 | 0.173 | 0.196 | 0.186 | 0.214 | 0.200 | 0.344 | 0.311 | 0.322 | 0.286 | 0.322
// cutoff=4000 | 0.202 | 0.229 | 0.238 | 0.211 | 0.196 | 0.194 | 0.312 | 0.334 | 0.392 | 0.327
// cutoff=8000 | 0.202 | 0.267 | 0.213 | 0.219 | 0.224 | 0.192 | 0.178 | 0.338 | 0.389 | 0.569
// cutoff=50k  | 0.188 | 0.250 | 0.200 | 0.216 | 0.175 | 0.171 | 0.167 | 0.187 | 0.214 | 0.584
// cutoff=150k | 0.159 | 0.201 | 0.197 | 0.222 | 0.192 | 0.175 | 0.178 | 0.156 | 0.225 | 0.360
