# Reconstruction

Surface reconstruction followed by hydrogen passivation. Useful for converting C(100)-(1×1) surfaces to C(100)-(2×1), fixing hydrogen collisions at concave corners, and removing methyl groups. Applies to any sp<sup>3</sup>-hybridized crystal. Relevant in both the `Cubic` and `Hexagonal` bases.

> TODO: Generalize this by making the hydrogen passivation component optional. There should be a member of the `Reconstruction` descriptor-like structure that includes the passivator element (default value: `.hydrogen`). Make a unit test for fluorine passivation, and another for no passivation at all. This doesn't provide the generality of mixed-element passivation across a structure, but it feels like the right API design choice.

Revised user-facing API:

```swift
struct Reconstruction {
  var atoms: [SIMD4<Float>]?
  var material: MaterialType?
  
  func compile() -> Topology
}

var reconstruction = Reconstruction()
reconstruction.atoms = ...
reconstruction.material = ...
let topology = reconstruction.compile()
```
