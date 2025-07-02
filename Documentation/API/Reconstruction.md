# Reconstruction

Surface reconstruction followed by hydrogen passivation. Useful for converting C(100)-(1×1) surfaces to C(100)-(2×1), fixing hydrogen collisions at concave corners, and removing methyl groups. Applies to any sp<sup>3</sup>-hybridized crystal. Relevant in both the `Cubic` and `Hexagonal` bases.

> TODO: Generalize this by making the hydrogen passivation component optional. Make a unit test for fluorine passivation, and another for no passivation at all. Implement this change once the internals of `Compilation` are more cleaned up and workable.

Revised user-facing API:

```swift
struct Reconstruction {
  var atoms: [SIMD4<Float>]?
  var material: MaterialType?
  var passivator: Element? = .hydrogen
  
  func compile() -> Topology
}

var reconstruction = Reconstruction()
reconstruction.atoms = ...
reconstruction.material = ...
let topology = reconstruction.compile()
```

The API doesn't provide the generality of mixed-element passivation across a structure, but it feels like the right design choice.
