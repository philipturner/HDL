# Reconstruction

> _Inspired by the PyTorch API choice, `torch.compile()`._

Surface reconstruction, _not_ followed by passivation. Useful for converting C(100)-(1×1) surfaces to C(100)-(2×1), fixing hydrogen collisions at concave corners, and removing methyl groups. Applies to any sp<sup>3</sup>-hybridized crystal. Relevant in both the `Cubic` and `Hexagonal` bases.

```swift
// Similar to 'Descriptor' APIs, the where the data required to initialize a
// complex object is specified, one property at a time. In line with these
// APIs, each property is 'Optional' and defaults to 'nil'. This contrasts with
// atoms and bonds in 'Topology', but it is a justifiable choice.
struct Reconstruction {
  /// Required.
  var atoms: [SIMD4<Float>]?
  
  /// Required.
  var material: MaterialType?
  
  func compile() -> Topology
}

// Example of usage.
var reconstruction = Reconstruction()
reconstruction.atoms = ...
reconstruction.material = ...
let topology = reconstruction.compile()
```
