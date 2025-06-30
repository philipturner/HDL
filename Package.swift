// swift-tools-version: 5.8
// The swift-tools-version declares the minimum version of Swift required to build this package.

import PackageDescription

let package = Package(
  name: "HDL",
  products: [
    // Products define the executables and libraries a package produces, making them visible to other packages.
    .library(
      name: "HDL",
      targets: ["HDL"]),
  ],
  dependencies: [
    .package(url: "https://github.com/apple/swift-atomics.git", .upToNextMajor(from: "1.3.0")),
    .package(url: "https://github.com/philipturner/swift-numerics", branch: "Quaternions"),
  ],
  targets: [
    // Targets are the basic building blocks of a package, defining a module or a test suite.
    // Targets can depend on other targets in this package and products from dependencies.
    .target(
      name: "HDL",
      dependencies: [
        .product(name: "Atomics", package: "swift-atomics"),
      ]),
    .testTarget(
      name: "HDLTests",
      dependencies: [
        "HDL",
        .product(name: "Numerics", package: "swift-numerics"),
      ]),
  ]
)
