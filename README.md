# Hardware Description Language

Domain-specific language for molecular nanotechnology.

Table of Contents
- [Overview](#overview)
    - [Objects](#objects)
    - [Scopes](#scopes)
- [Operations](#operations)
    - [Filter](#filter)
    - [Lattice](#lattice)
    - [Volume](#volume)
- [Tips](#tips)

## Overview

TODO: Provide a general explanation.

### Objects

```swift
enum Element { ... }

enum Bond {
  case sigma
}

enum EntityType {
  case atom(Element)
  case bond(Bond)
  case empty
}

struct Entity {
  var position: SIMD3<Float>
  var type: EntityType
}
```

`Entity` is a data structure that stores atoms and bond connectors. The position occupies 12 bytes and the entity type occupies 4 bytes. This format aligns the entity to a 16-byte vector word, improving compilation speed.

`EntityType` can extract information about the type, such as atomic number or bond order. Atomic numbers can be any element from the [MM4 force field](https://github.com/philipturner/MM4) (H, C, N, O, F, Si, P, S, and Ge). Zero for the type's raw value indicates `.empty`.

```swift
Lattice<Basis> { h, k, l in
  Material { ... }
  Bounds { ... }
}
Lattice<Basis>.entities
```

Object encapsulating crystal plane algebra.

Creates a lattice of crystal unit cells to edit. Coordinates are represented in numbers of crystal unit cells. The coordinate system may be mapped to a non-orthonormal coordinate system internally. Keep this in mind when processing `SIMD3<Float>` vectors. For example, avoid normalizing any vectors.

```swift
Topology { 
  Copy { [Entity] }
}
Topology.atomicNumbers: [UInt8]
Topology.bonds: [SIMD2<UInt32>]
Topology.positions: [SIMD3<Float>]
```

Object encapsulating bond topology filters.

Creates a list of atoms in Morton order, with sigma bonds connecting them. Free radicals are not yet passivated. Geometric overlap between potential passivators is detected.

The compiled topology is provided through a set of properties. These properties can be entered directly into `MM4ParametersDescriptor` or `MM4RigidBodyDescriptor`. If an entity exists where passivators collide, its atomic number is `0`.

### Scopes

```swift
Concave { }
```

Scope where every plane's "one" volume merges through AND in [disjunctive normal form](https://en.wikipedia.org/wiki/Disjunctive_normal_form). Upon exiting this scope, the added planes remain. This must be called inside a `Volume`.

```swift
Convex { }
```

Scope where every plane's "one" volume merges through OR in [disjunctive normal form](https://en.wikipedia.org/wiki/Disjunctive_normal_form). Upon exiting this scope, the added planes remain. This must be called inside a `Volume`.

```swift
Volume { }
```

Encapsulates a set of planes, so that everything inside the scope is removed from the stack upon exiting. This must be called inside `Lattice` and may be called inside another `Volume`.

## Operations

### Filter

The following documentation describes `Filter`, which may only be called inside a `Topology`.

```swift
Filter(FilterType)
Filter { atom, neighbors in ... }

// Shorthand for the function signature.
typealias FilterType = (
  atom: inout Entity, neighbors: inout [Entity]
) -> Void
```

Modify the topology using a custom filter. Three classes of filters are permitted:
- `add bonds` - change the entity type of empty neighbors to `.bond(.sigma)`.
- `remove atom` - change the entity type of `atom` to `.empty`.
- `passivate` - append atoms to the neighbors, keeping existing atoms intact.

After passivation, the neighbor list must equal the valence count. The valence shell includes only sp<sup>3</sup>-hybridized orbitals. In addition, added passivators must all have the same element. The following atom and passivator combinations are permitted:

| Atom      | Passivator | Valence |
| --------- | ---------- | ------- |
| carbon    | hydrogen   | 4       |
| carbon    | fluorine   | 4       |
| silicon   | hydrogen   | 4       |
| germanium | hydrogen   | 4       |

The atom's position can be adjusted during the filter. New atoms are copied into a separate list while the closure is called. The adjustment will not affect the value of existing neighbors during the function call. This functionality could be used to adjust carbon atom positions when reconstructing diamond (100) surfaces.

```swift
Filter.connectSharpCorners: FilterType
Filter.removePrimaryAtoms: FilterType
Filter.hydrogenPassivate: FilterType

Topology {
  Copy { ... }
  
  Filter(Filter.connectSharpCorners)
  Filter(Filter.removePrimaryAtoms)
  Filter(Filter.hydrogenPassivate)
}
```

A sequence of filters for cleaning up geometry.
1. Colliding passivators are replaced with sigma bonds, if both atoms have a single collision.
2. Primary carbons (methyl and trifluoromethyl groups) are removed.
3. All free radicals are passivated with hydrogen, except those with remaining passivator collisions.

> TODO: For (1), should connected atoms be displaced to shorten the bond? The passivator detection algorithm should be robust enough to handle slight displacements.

```swift
Filter.reconstructCubic100(SIMD3<Float>): FilterType

// The simplest way to call the filter.
let direction = SIMD3<Float>(1, 1, 1)
Filter(Filter.reconstructCubic100(x))

// You can exclude this filter from certain atoms.
// For example, you may want to reconstruct bonds in
// a different direction for different faces of a
// crystolecule.
Filter { atom, neighbors in
  let direction = SIMD3<Float>(1, 1, 1)
  guard condition(atom, neighbors) else {
    return
  }
  Filter.reconstructCubic100(direction)(atoms, neighbors)
}
```

A filter for cleaning up diamond (100) surfaces. Sigma bonds are generated approximately parallel to the specified direction. Surface atoms are displaced to shorten the bond. Passivating hydrogens are created at an uneven angle. This filter may fail to reconstruct bonds in certain cases.

 An atom located roughly at (0, 0, 0) will form a bond pointing in the specified direction. This fact can be used to change the parity of which atoms are connected. Flip the sign of the direction to alternate which atoms are connected.

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

> TODO: Change the hexagonal lattice constants to the true energy minima, rather than something scaled to align with cubic (111) surfaces.

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
    Plane { +h }
    Plane { +h2k }
    Plane { +l }
    ...
  }
}
```

Unit vectors representing the crystal's basis.

```swift
Material { MaterialType }
```

Specifies the atom types to fill the lattice with, and the lattice constant. This must be called in the top-level scope.

### Volume

The following keywords may be called inside a `Volume`.

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
Ridge(SIMD3<Float>) { SIMD3<Float> }
Ridge(normal) { reflector }
Valley(SIMD3<Float>) { SIMD3<Float> }
Valley(normal) { reflector }
```

Reflect `normal` across `reflector`. Generate a plane with the normal before reflection, then a second plane after the reflection. `Ridge` takes the union of the planes' "one" volumes, while `Valley` takes the intersection.

```swift
Replace { EntityType }
```

Replace all entities in the selected volume with a new entity.

To delete atoms, use `Replace { .empty }`. Removed atoms cannot be restored by a subsequent `Replace`.
