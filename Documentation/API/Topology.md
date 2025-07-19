# Topology

The following APIs are available for `Topology`.

```swift
extension Topology {
  enum MapNode {
    case atoms
    case bonds
  }
  
  // Underlying storage is a fixed-width vector, maximum capacity is 8.
  struct MapStorage: Collection {
    subscript(position: Int) -> UInt32
  }
}

func map(
  _ sourceNode: MapNode,
  to targetNode: MapNode
) -> [MapStorage]

// Example of usage.
var topology = Topology()
topology.atoms = atoms
topology.bonds = bonds
let atomsToAtomsMap = topology.map(.atoms, to: .atoms)
let atomsToBondsMap = topology.map(.atoms, to: .bonds)
```

Create a map that points from atoms/bonds to a list of connected atoms/bonds.

The results of this function follow a few particular rules:
- The source node must be `.atoms`.
- The number of targets for a given source node must not exceed 8.
- The indices are not guaranteed to be sorted in any particular order, for performance reasons. The underlying implementation uses multithreading + atomic synchronization when the bond count exceeds 2500.

```swift
extension Topology {
  enum MatchAlgorithm {
    // Search for neighbors within a fixed radius (in nanometers).
    case absoluteRadius(Float)
    
    // Search for neighbors within a multiple of covalent bond length.
    case covalentBondLength(Float)
  }
}

func match(
  _ source: [Atom], 
  algorithm: MatchAlgorithm = .covalentBondLength(1.5),
  maximumNeighborCount: Int = 8
) -> [ArraySlice<UInt32>]

// Example of usage.
var topology = Topology()
topology.atoms = atom1
let closeMatches = topology.match(atoms2)
let farMatches = topology.match(
  atoms2, algorithm: .covalentBondLength(2))
let angstromMatches = topology.match(
  atoms2, algorithm: .absoluteRadius(0.1))
```

Match the `source` atoms to target atom indices, where the targets reside in the `Topology`. The output has the same dimensions as the input array.

If the total number of atoms within the search radius exceeds the `maximumNeighborCount`, results are undefined. The internal matching algorithm relies on a given amount of memory being allocated. This memory temporarily holds references to other atoms. After all qualifying atoms are found, they are sorted in order of ascending distance. If the memory capacity is exceeded, there is nowhere to store the remaining references, even if one contains the closest atom.

If you specify a specific search radius, you should be able to calculate the expected neighbor count. Each material has a finite number of atoms per cubic nanometer. The density in atoms/nm<sup>3</sup> can be multiplied by the volume of a sphere with the given radius. If you are unsure exactly how many atoms to expect, you can provide a generous upper bound for memory capacity.

For `covalentBondScale`, bond length is determined by summing the covalent radii. The pairwise sum does not always equal the bond length predicted by the lattice constant. The default scale of 1.5 provides enough tolerance for 50% error in approximated bond length.

`maximumNeighborCount` is 8 by default. This means most `match()` invocations searching beyond a single covalent bond length will fail. You may increase the limit for returned neighbors, but doing so may increase compile time or memory consumption.

You are encouraged to sort the topology before calling `match()`. Otherwise, the search algorithm may degrade from $O(n)$ to $O(n^2)$. The overhead of sorting is significant and often takes more time than just running the match. Therefore, the internal implementation only performs sorting when the atom count is ~10,000. This is a performance sweet spot for highly ordered distributions (e.g. atoms directly fetched from a crystal lattice). However, it may not be a sweet spot for extremely disordered distributions.

```swift
extension Topology {
  enum OrbitalHybridization {
    case sp
    case sp2
    case sp3
  }
  
  // Underlying storage is a fixed-width vector, maximum capacity is 2.
  struct OrbitalStorage: Collection {
    subscript(position: Int) -> SIMD3<Float>
  }
}

func nonbondingOrbitals(
  hybridization: OrbitalHybridization = .sp3
) -> [OrbitalStorage]

// Example of usage.
let orbitalLists = topology.nonbondingOrbitals()
let orbitalList = orbitalLists[atomID]
let orbital = orbitalList[0]
```

Directions for nonbonding orbitals in the valence shell. The directions are represented by normalized vectors.

If the directions cannot be determined with absolute certainty from the immediate neighbors, no orbitals are reported. Examples are methane carbons, primary carbons in the sp3 hybridization, or sp2-bonded carbons with only a single bond filled. There are multiple reasonable heuristics for deciding where to place passivators in such edge cases. The choice of a heuristic is best deferred to the user. Furthermore, the appearance of such edge cases often signals malformed geometry or functional groups that ought to be removed.

Free radicals and lone pairs are treated the same way. This means a nitrogen with one missing passivator will return two options for N-H bond directions. Although the compiler cannot determine the passivator direction with certainty, it can narrow down a discrete set of choices.

For the carbon in an acetylene radical, only one orbital is reported. The reported orbital contains the free radical and is collinear with the two carbon atoms. It is also the only orbital known with absolute positional certainty. The two pi orbitals could be rotated into a infinite number of specific positions around the axis. It may be possible to exactly determine their orientation relative to other pi orbitals in a carbyne rod. However, that heuristic involves more than just immediate neighbors.

Another edge case is halogens. They have three sp3 lone pairs, which cannot be positionally constrained. Therefore, the compiler cannot generate a discrete set of orbital orientations for halogens.  Nonbonding orbitals of divalent oxygen and trivalent nitrogen can be computed analytically from the neighbor atom positions. The discrepancy between treatment of Group VII and Group V/VI has an analogue in purely hydrocarbon matter. Primary carbons report zero orbitals, while secondary and tertiary carbons report 1&ndash;2 orbitals.

```swift
mutating func remove(atoms indices: [UInt32])
mutating func remove(bonds indices: [UInt32])
```

Removes atoms or bonds at the specified indices. When removing an atom, the bonds connected to the atom are also removed.

An index may be specified multiple times in the input. The atom or bond will only be removed once.

The order of atoms and bonds is preserved after removal. The removed items are taken out of the list, and the remainder are compacted in place. This behavior is different from `sort()`, which scrambles the relative order of atoms.

```swift
// Sorts the atoms and returns the old atoms' indices in the new list.
@discardableResult
mutating func sort() -> [UInt32]
```

Sorts atoms in Morton order, then sorts bonds in ascending order based on atom indices.

The topology should be sorted before entering into a simulator. Sorting causes nearby atoms to appear in consecutive memory locations. OpenMM utilizes this spatial locality to compute nonbonded forces very fast. If you forget to sort, the algorithmic complexity may increase from $O(n)$ to $O(n^2)$.

Note that the spatial locality of the crystal lattice (8&ndash;12 wide atom blocks) typically achieves the same effect as sorting. But true Morton order results in the maximum possible degree of spatial locality.
