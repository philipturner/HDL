# Volume

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

Translate the origin by a vector relative to the current origin. Modifications to the origin are undone when leaving the current scope. This may not be called in the top-level scope.

```swift
Plane { SIMD3<Float> }
```

Add a plane to the stack. The plane will be combined with other planes, and used for selecting/deleting atoms.

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

Encapsulate a set of planes, so that everything inside the scope is removed from the stack upon exiting. This may be called at the top level or inside another `Volume`. In most cases, it is best to only call once, at the top level.

## Difference Between Concave, Convex, and Volume

In the library code, a global singleton is notified when each scope starts or ends. For example, the singleton updates at the opening and closing brackets of `Volume`. This scope information allows a recursive stack to be held as the DSL is parsed. The three flavors of scope differ in how they affect destruction or merging of other scopes, for example at the closing bracket.

```swift
enum LatticeScopeType {
  case concave
  case convex
  case volume
  
  var modifiesPredecessor: Bool {
    self != .volume
  }
  
  var usesLogicalAnd: Bool {
    self == .concave
  }
}
```

The stack may only begin recording information (`Origin`, `Plane`, `Replace`) once the lattice constant, working grid size, and atom types are known. This is why `Bounds` and `Material` must be called before `Volume`. The global singleton keeps recording commands until it's notified that the `Lattice` reached its closing bracket. Then, the `.atoms` property materializes and a new `Lattice` can be started.

The embedded DSL design has two implications. First, it is not thread-safe. Multiple threads modifying the same global singleton will create a data race. Second, `Lattice` instances cannot be declared recursively inside of each other. Only one lattice can be recorded at a time.
