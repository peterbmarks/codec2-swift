// swift-tools-version:5.9
import PackageDescription

let package = Package(
    name: "Codec2Swift",
    platforms: [
        .macOS(.v11)
    ],
    products: [
        .library(name: "Codec2", targets: ["Codec2"]),
        .executable(name: "c2enc", targets: ["c2enc"]),
        .executable(name: "c2dec", targets: ["c2dec"]),
    ],
    targets: [
        .target(
            name: "Codec2",
            path: "Sources/Codec2"
        ),
        .executableTarget(
            name: "c2enc",
            dependencies: ["Codec2"],
            path: "Sources/c2enc"
        ),
        .executableTarget(
            name: "c2dec",
            dependencies: ["Codec2"],
            path: "Sources/c2dec"
        ),
        .testTarget(
            name: "Codec2Tests",
            dependencies: ["Codec2"],
            path: "Tests/Codec2Tests",
            resources: [
                .copy("Reference")
            ]
        ),
    ]
)
