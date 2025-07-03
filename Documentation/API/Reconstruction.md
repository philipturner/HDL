# Reconstruction

Surface reconstruction followed by passivation. Useful for converting C(100)-(1×1) surfaces to C(100)-(2×1), fixing hydrogen collisions at concave corners, and removing methyl groups. Applies to any sp<sup>3</sup>-hybridized crystal. Relevant in both the `Cubic` and `Hexagonal` bases.

> TODO: Generalize this by making the hydrogen passivation component optional. Make a unit test for fluorine passivation, and another for no passivation at all. Implement this change once the internals of `Compilation` are more cleaned up and workable.

Revised user-facing API:

```swift
struct Reconstruction {
  /// Required.
  var atoms: [SIMD4<Float>]?
  
  /// Required.
  var material: MaterialType?
  
  /// Optional. The default is `nil`. 
  var passivation: Element?
  
  func compile() -> Topology
}

var reconstruction = Reconstruction()
reconstruction.atoms = ...
reconstruction.material = ...
let topology = reconstruction.compile()
```

`passivation` places hydrogens according to lattice-aligned geometry before the reconstruction. Otherwise, you're stuck with manually generated hydrogens that have shifted.

The API doesn't provide the generality of mixed-element passivation across a structure, but it feels like the right design choice.
