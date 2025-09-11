# Hardware Description Language

Domain-specific language for molecular nanotechnology. This repository includes a geometry and bond topology compiler.

Table of Contents
- [Overview](#overview)
- [API](#)
  - [Lattice](./Documentation/API/Lattice.md)
  - [Reconstruction](./Documentation/API/Reconstruction.md)
  - [Topology](./Documentation/API/Topology.md)
  - [Volume](./Documentation/API/Volume.md)
- [Testing](#testing)

## Overview

For an introduction, visit the [tutorial](./Documentation/Tutorial/GrapheneSiliceneBilayer.md).

```swift
typealias Atom = SIMD4<Float>

extension Atom {
  var position: SIMD3<Float>
  var element: Element
  
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
let lattice = Lattice<Basis> { h, k, l in
  Bounds { ... }
  Material { ... }
  
  Volume {
    ...
  }
}

// Property to retrieve the geometry.
let atoms = lattice.atoms
```

Creates a lattice of crystal unit cells to edit. Coordinates are represented in numbers of crystal unit cells. In the trailing closure of `Lattice`, a group of keywords can be written, which constitute a DSL.

```swift
struct Topology {
  var atoms: [Atom]
  var bonds: [SIMD2<UInt32>]
}
```

Encapsulates fundamental operations for bond topology formation:
- `map` - query the bonds that reference each atom
- `match` - $O(n)$ neighbor search
- `nonbondingOrbitals` - positions of unpaired electrons
- `remove` - update the bonds after atoms shift to fill a gap in the list
- `sort` - reorder the atoms to improve memory locality

## Testing

Some unit tests are disabled by default. They take too much time to execute in debug mode. However, Swift release mode takes too long to compile the tests. The solution is to gate the tests under the macro `RELEASE`. Then, hack the debug-mode compiler to actually compile in release mode, but with incremental compilation.

Enter the following at command-line to run all of the tests:

```
swift test -Xswiftc -Ounchecked -Xswiftc -DRELEASE
```
