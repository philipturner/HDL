# Reconstruction

> _Inspired by the PyTorch API choice, `torch.compile()`._

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

## Passivation?

Is it possible to completely isolate the H-passivation functionality from the surface reconstruction functionality?

`createHydrogenBonds()` does not depend on the material. So it could be factored out of Reconstruction and belong somewhere else. Maybe.

The two APIs could share some utility functions from the internal `Compilation` type. Perhaps re-computing the maps from atoms to hydrogen sites. These maps are where the main uncertainty lies. Need to clean up `Compilation` more, to get greater clarity.

Tentative API:

```swift
struct Passivation {
  /// Required.
  var atoms: [SIMD4<Float>]?
  
  /// Required.
  var bonds: [SIMD2<UInt32>]?
  
  /// Required.
  var element: Element?
  
  func compile() -> Topology
}

var passivation = Passivation()
passivation.atoms = topology.atoms
passivation.bonds = topology.bonds
passivation.element = .hydrogen
let passivatedTopology = passivation.compile()
```

Another realization: passivation would apply to 2D graphitic structures. And it would be less general-purpose, not satisfying many common use cases. Rather, it's a quick convenience function built into the library. Reduces the code size of typical client code. Reducing client code size alone is not a sufficient condition for inclusion; we care about minimalism and giving the user greater control.

Another possibility: passivation is an inherent side effect of surface reconstruction. An important piece of information that would be lost, if not exposed to the public API.
