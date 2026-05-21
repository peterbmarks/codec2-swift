import XCTest
@testable import Codec2

final class PackTests: XCTestCase {

    func testRoundTripGray() {
        var bits = [UInt8](repeating: 0, count: 16)
        var idx: UInt32 = 0
        let values: [(Int32, UInt32)] = [
            (5, 4), (0, 1), (37, 6), (1, 1), (255, 8), (7, 3)
        ]
        for (v, w) in values { pack(&bits, bitIndex: &idx, field: v, fieldWidth: w) }

        var rIdx: UInt32 = 0
        for (v, w) in values {
            let got = unpack(bits, bitIndex: &rIdx, fieldWidth: w)
            XCTAssertEqual(got, v, "gray round-trip mismatch for value \(v) width \(w)")
        }
        XCTAssertEqual(idx, rIdx)
    }

    func testRoundTripNatural() {
        var bits = [UInt8](repeating: 0, count: 16)
        var idx: UInt32 = 0
        let values: [(Int32, UInt32)] = [
            (0xAB, 8), (0x3, 2), (0x1F, 5), (0x7FF, 11)
        ]
        for (v, w) in values {
            packNaturalOrGray(&bits, bitIndex: &idx, field: v, fieldWidth: w, gray: 0)
        }

        var rIdx: UInt32 = 0
        for (v, w) in values {
            let got = unpackNaturalOrGray(bits, bitIndex: &rIdx, fieldWidth: w, gray: 0)
            XCTAssertEqual(got, v)
        }
    }

    /// pack() with Gray = 1 must produce identical bytes to the C reference for
    /// the same input. The encoded bit pattern depends on the Gray mapping;
    /// these expected bytes were captured from the C implementation manually
    /// for a small fixed input so we can detect bit-level drift.
    func testKnownVectorGray() {
        // Pack the value 5 in a 4-bit field. Gray(5) = 5 ^ (5>>1) = 5 ^ 2 = 7 (0b0111).
        // With bitIndex starting at 0 we expect the high nibble of byte 0 to be 0x7.
        var bits = [UInt8](repeating: 0, count: 1)
        var idx: UInt32 = 0
        pack(&bits, bitIndex: &idx, field: 5, fieldWidth: 4)
        XCTAssertEqual(bits[0], 0x70)
        XCTAssertEqual(idx, 4)
    }
}
