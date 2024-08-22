# Hardware Description Language

Domain-specific language for molecular nanotechnology. This repository includes a geometry and bond topology compiler.

Table of Contents
- [Overview](#overview)
- [Operations](#operations)
    - [Lattice](#lattice)
    - [Topology](#topology)
    - [Volume](#volume)
- [Testing](#testing)

## Overview

For an introduction, visit the [tutorial](./Documentation/GrapheneSiliceneBilayer.md).

```swift
typealias Atom = SIMD4<Float>

extension Atom {
  var position: SIMD3<Float>
  var atomicNumber: UInt8
  
  init(position: SIMD3<Float>, element: Element)
}
```

`Atom` is a typealias for a four-wide SIMD vector. Each vector lane contains an IEEE 754 single-precision floating point number. The first three lanes store the X, Y, and Z coordinates in real space. The fourth lane stores the atomic number.

```swift
enum Element: UInt8 {
  case hydrogen = 1
  
  case boron = 5
  case carbon = 6
  case nitrogen = 7
  case oxygen = 8
  case fluorine = 9
  
  case aluminum = 13
  case silicon = 14
  case phosphorus = 15
  case sulfur = 16
  case chlorine = 17
  
  case gallium = 31
  case germanium = 32
  case arsenic = 33
  case selenium = 34
  case bromine = 35
  
  case tin = 50
  case gold = 79
  case lead = 82
  
  var covalentRadius: Float { get }
}
```

Atomic numbers can be any element from Group III - VII, Period II - IV of the periodic table. In addition, a few heavy metals are parameterized.

```swift
// Specify lattice edits, if any, in the trailing closure.
Lattice<Basis> { h, k, l in
  Material { ... }
  Bounds { ... }
}

// Property to retrieve the geometry.
Lattice<Basis>.atoms
```

Object encapsulating crystal plane algebra.

Creates a lattice of crystal unit cells to edit. Coordinates are represented in numbers of crystal unit cells. The coordinate system may be mapped to a non-orthonormal coordinate system internally. Keep this in mind when processing `SIMD3<Float>` vectors. For example, avoid normalizing any vectors.

```swift
// Initialize an empty topology. Geometry will be added using member functions.
var topology = Topology()
```

Encapsulates low-level operations during bond topology formation. These include $O(n)$ neighbor searching, insertion/removal of atoms/bonds, and Morton reordering.

## Operations

### Lattice

The following keywords may be called inside a `Lattice`.

```swift
protocol Basis
Cubic: Basis
Hexagonal: Basis
```

Coordinate spaces for defining vectors in.

```swift
Bounds { SIMD3<Float> }
Bounds { 10 * h + 10 * k + 10 * l } // cubic
Bounds { 10 * h + 10 * (h + 2 * k) + 10 * l } // hexagonal
```

Sets the working set of crystal unit cells. The box spans from the world origin `[0, 0, 0]` the specified vector. This must be called in the top-level scope, before any `Volume` keywords.

For hexagonal crystals, the bounds are a cuboid defined by transforming the input vector. It is mapped from h/k/l space to h/h2k/l space. This allows the base lattice to be cartesian, similar to cubic. The quantity in each axis direction must be an integer.

```swift
Constant(ConstantType) { MaterialType }
ConstantType.hexagon // Hexagonal - hexagon side length
ConstantType.prism   // Hexagonal - prism height
ConstantType.square  // Cubic - square side length

// Query the lattice constant for diamond.
let latticeConstant = Constant(.square) { .elemental(.carbon) }
```

Values of the lattice constants, for use in custom geometry processing code.

The hexagonal lattice constants are not the exact energy minima. Rather, they are scaled to align with cubic (111) surfaces.

```swift
// Lattice vectors originate from the smallest repeatable unit of crystal. For
// cubic crystals, they are edges of a cube. For hexagonal crystals, they are
// sides of a hexagonal prism. The vectors aren't always orthogonal, so they are
// internally translated to nanometers before applying transforms.

// Hexagonal crystals are sometimes described with four unit vectors: h, k, i,
// and l. The 'i' vector is redundant and equals -h - k, creating a set of 3
// vectors symmetric around the perimeter of a hexagon. However, for the HDL, it
// can make representation more concise.
Lattice<Hexagonal> { h, k, l in
  // h + k + i == 0
  let i = -h - k
}

// Another helpful technique, which makes Hexagonal more similar to Cubic, is to
// replace 'k' with something orthogonal to 'h'. The coordinate basis has
// changed from h/k/l to h/h + 2k/l.
Lattice<Hexagonal> { h, k, l in
  // [3 * h, 2 * h2k, 2 * l] forms something close to a cube.
  let h2k = h + 2 * k
  ...
  
  Volume {
    Plane { -h }
    Plane { -h2k }
    Plane { -l }
    
    Origin { 3 * h + 2 * (h2k + l) }
    Plane { h }
    Plane { h2k }
    Plane { l }
    ...
  }
}
```

Unit vectors representing the crystal's basis.

```swift
Material { MaterialType }
```

Specifies the atom types to fill the lattice with, and the lattice constant. This must be called in the top-level scope.

### Topology

The following APIs are available for `Topology`.

```swift
var atoms: [Atom] { get set }
var bonds: [SIMD2<UInt32>] { get set }
```

The atoms and bonds backing the topology.

> WARNING: Although `atoms` and `bonds` are mutable, the setters do not enforce self-consistency. The intended use case is overwriting atoms/bonds in-place (e.g. adjusting an atom's position). Use `insert()` and `remove()` when the size of these arrays changes.

```swift
extension Topology {
  enum MapNode {
    case atoms
    case bonds
  }
}

func map(
  _ sourceNode: MapNode,
  to targetNode: MapNode
) -> [ArraySlice<UInt32>]

// Example of usage.
var topology = Topology()
topology.insert(atoms: atoms)
topology.insert(bonds: bonds)

let atomsToAtomsMap = topology.map(.atoms, to: .atoms)
let atomsToBondsMap = topology.map(.atoms, to: .bonds)
let bondsToAtomsMap = topology.map(.bonds, to: .atoms)
```

Create a map that points from atoms/bonds to a list of connected atoms/bonds.

The results of this function follow a few particular rules:
- The source and target nodes cannot both be `.bonds`.
- The number of targets for a given source node cannot exceed 8.
- If one of the nodes is `.bonds`, the target indices are sorted. Otherwise, the target indices (which are `.atoms`) correspond to bonds in ascending order.

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
topology.insert(atoms: atoms1)
let closeMatches = topology.match(atoms2)
let farMatches = topology.match(
  atoms2, algorithm: .covalentBondLength(2))
let angstromMatches = topology.match(
  atoms2, algorithm: .absoluteRadius(0.1))
  
// Another example of usage.
var atoms = [Atom(...), Atom(...)]
let atomsLocationsInStructure = topology.match(atoms)
someFunction(atoms[0], atomLocationsInStructure[0])
someFunction(atoms[1], atomLocationsInStructure[1])
```

Match the `source` atoms to target atom indices, where the targets reside in the `Topology`. The output has the same dimensions as the input array.

If the total number of atoms within the search radius exceeds the `maximumNeighborCount`, results are undefined. The internal matching algorithm relies on a given amount of memory being allocated. This memory temporarily holds references to other atoms. After all qualifying atoms are found, they are sorted in order of ascending distance. If the memory capacity is exceeded, there is nowhere to store the remaining references, even if one contains the closest atom.

If you specify a specific search radius, you should be able to calculate the expected neighbor count. Each material has a finite number of atoms per cubic nanometer. The density in atoms/nm<sup>3</sup> can be multiplied by the volume of a sphere with the given radius. If you are unsure exactly how many atoms to expect, you can provide a generous upper bound for memory capacity.

For `covalentBondScale`, bond length is determined by summing the covalent radii. The pairwise sum does not always equal the bond length from `MM4Parameters`; add some tolerance for such error. The default value of 1.5 provides enough tolerance for 50% error in approximated bond length.

`maximumNeighborCount` is 8 by default. This means most `match()` invocations searching beyond a single covalent bond length will fail. You may increase the limit for returned neighbors, but doing so may increase compile time or memory consumption.

You are encouraged to sort the topology before calling `match()`. Otherwise, the search algorithm may degrade from $O(n)$ to $O(n^2)$. The overhead of sorting is significant and often takes more time than just running the match. Therefore, the internal implementation only performs sorting when the atom count is ~10,000. This is a performance sweet spot for highly ordered distributions (e.g. atoms directly fetched from a crystal lattice). However, it may not be a sweet spot for extremely disordered distributions.

```swift
mutating func insert(atoms: [Atom])
mutating func insert(bonds: [SIMD2<UInt32>])
```

Adds new atoms/bonds to the topology.

The new atoms and bonds are added directly to the end of the list. The atoms and bonds are not checked for duplicates. For example, if you specify a bond multiple times, it will be added multiple times. This behavior is different from `remove()`, which checks for duplication before removing.

```swift
mutating func remove(atoms indices: [UInt32])
mutating func remove(bonds indices: [UInt32])
```

Removes atoms/bonds at the specified indices. For `removeAtoms`, bonds connected to the removed atoms are also removed.

An index may be specified multiple times in the input. The atom or bond will only be removed once.

The order of atoms and bonds is preserved after removal. The removed items are taken out of the list, and the remainder is compacted in place. This behavior is different from `sort()`, which scrambles the relative order of atoms.

```swift
extension Topology {
  enum OrbitalHybridization {
    case sp1
    case sp2
    case sp3
  }
}

func nonbondingOrbitals(
  hybridization: OrbitalHybridization = .sp3
) -> [ArraySlice<SIMD3<Float>>]
```

Directions for nonbonding orbitals in the valence shell. The directions are represented by normalized vectors.

If the directions cannot be determined with absolute certainty from the immediate neighbors, no orbitals are reported. Examples are methane carbons, primary carbons in the sp3 hybridization, or sp2-bonded carbons with only a single bond filled. There are multiple reasonable heuristics for deciding where to place passivators in such edge cases. The choice of a heuristic is best deferred to the user. Furthermore, the appearance of such edge cases often signals malformed geometry or functional groups that ought to be removed.

Free radicals and lone pairs are treated the same way. This means a nitrogen with one missing passivator will return two options for N-H bond directions. Although the compiler cannot determine the passivator direction with certainty, it can narrow down a discrete set of choices.

For the carbon in an acetylene radical, only one orbital is reported. The reported orbital contains the free radical and is collinear with the two carbon atoms. It is also the only orbital known with absolute positional certainty. The two pi orbitals could be rotated into a infinite number of specific positions around the axis. It may be possible to exactly determine their orientation relative to other pi orbitals in a carbyne rod. However, that heuristic involves more than just immediate neighbors.

Another edge case is halogens. They have three sp3 lone pairs, which cannot be positionally constrained. Therefore, the compiler cannot generate a discrete set of orbital orientations for halogens.  Nonbonding orbitals of divalent oxygen and trivalent nitrogen can be computed analytically from the neighbor atom positions. The discrepancy between treatment of Group VII and Group V/VI has an analogue in purely hydrocarbon matter. Primary carbons report zero orbitals, while secondary and tertiary carbons report 1&ndash;2 orbitals.

```swift
// Sorts the atoms and returns the old atoms' indices in the new list.
@discardableResult
mutating func sort() -> [UInt32]
```

Sorts atoms in Morton order, then sorts bonds in ascending order based on atom indices.

The topology should be sorted before entering into a simulator. Sorting causes nearby atoms to appear in consecutive memory locations. OpenMM utilizes this spatial locality to compute nonbonded forces very fast. If you forget to sort, the algorithmic complexity may increase from $O(n)$ to $O(n^2)$.

### Volume

The following keywords may be called inside a `Volume`.

```swift
Concave { }
```

Scope where every plane's "one" volume merges through AND in [disjunctive normal form](https://en.wikipedia.org/wiki/Disjunctive_normal_form). Upon exiting this scope, the added planes remain. This must be called inside a `Volume`.

```swift
Convex { }
```

Scope where every plane's "one" volume merges through OR in [disjunctive normal form](https://en.wikipedia.org/wiki/Disjunctive_normal_form). Upon exiting this scope, the added planes remain. This must be called inside a `Volume`.

```swift
Origin { SIMD3<Float> }
```

Translates the origin by a vector relative to the current origin. Modifications to the origin are undone when leaving the current scope. This may not be called in the top-level scope.

```swift
Plane { SIMD3<Float> }
```

Adds a plane to the stack. The plane will be combined with other planes, and used for selecting/deleting atoms.

A `Plane` divides the `Bounds` into two sections. The "one" volume is the side the normal vector points toward. The "zero" volume is the side the normal points away from. The "one" volume contains the atoms modified during a `Replace`. When planes combine into a `Concave`, only the crystal unit cells common to every plane's "one" volume are modifiable.

```swift
enum ReplaceType {
  case .atom(Element)
  case .empty
}

Replace { ReplaceType }
```

Transmute the atoms in the selected volume to a different atom type.

To delete atoms, use `Replace { .empty }`. Removed atoms cannot be restored by a subsequent `Replace`.

```swift
Volume { }
```

Encapsulates a set of planes, so that everything inside the scope is removed from the stack upon exiting. This must be called inside `Lattice` and may be called inside another `Volume`.

## Testing

Some unit tests are disabled by default. They take too much time to execute in debug mode. However, Swift release mode takes too long to compile the tests. The solution is to gate the tests under the macro `RELEASE`. Then, hack the debug-mode compiler to actually compile in release mode, but with incremental compilation.

Enter the following at command-line to run all of the tests:

```
swift test -Xswiftc -Ounchecked -Xswiftc -DRELEASE
```
