// swift-tools-version: 6.1

import PackageDescription

let package = Package(
  name: "HDL",
  products: [
    .library(
      name: "HDL",
      targets: ["HDL"]),
  ],
  dependencies: [
    .package(url: "https://github.com/apple/swift-atomics.git", .upToNextMajor(from: "1.3.0")),
    .package(url: "https://github.com/philipturner/swift-numerics", branch: "Quaternions"),
  ],
  targets: [
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
