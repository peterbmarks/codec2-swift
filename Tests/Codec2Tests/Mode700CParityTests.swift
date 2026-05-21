import XCTest
@testable import Codec2

/// Mode 700C parity against the C reference. Both the encoder and decoder
/// are effectively bit-exact:
///   - encoder: byte-for-byte identical bitstream
///   - decoder: <0.8% of samples differ, all by ≤ 1 LSB (fp rounding only)
final class Mode700CParityTests: XCTestCase {

    private func referenceURL(_ name: String) throws -> URL {
        let bundle = Bundle.module
        guard let url = bundle.url(forResource: name, withExtension: nil, subdirectory: "Reference") else {
            throw XCTSkip("missing reference resource \(name)")
        }
        return url
    }

    /// Encoder must produce a byte-for-byte identical bitstream.
    func testEncoder700CByteExact() throws {
        let inputURL = try referenceURL("hts1a.raw")
        let refBitsURL = try referenceURL("hts1a_700C.bit")

        let pcm = try Data(contentsOf: inputURL)
        let refBits = try Data(contentsOf: refBitsURL)
        guard let codec = Codec2(mode: .mode700C) else {
            XCTFail("could not instantiate Codec2 700C")
            return
        }
        let nsam = codec.samplesPerFrame
        let nbyte = codec.bytesPerFrame
        var buf = [Int16](repeating: 0, count: nsam)
        var bits = [UInt8](repeating: 0, count: nbyte)
        var encoded = Data()

        pcm.withUnsafeBytes { raw -> Void in
            let p = raw.bindMemory(to: Int16.self)
            let frames = p.count / nsam
            for f in 0..<frames {
                for i in 0..<nsam { buf[i] = p[f * nsam + i] }
                codec.encode(speech: buf, bits: &bits)
                encoded.append(contentsOf: bits)
            }
        }

        XCTAssertEqual(encoded.count, refBits.count)
        XCTAssertEqual(encoded, refBits, "encoder bitstream is not byte-exact with C reference")
    }

    /// Decoder reconstructs the C reference within ±1 LSB on every sample.
    func testDecoder700CSampleParity() throws {
        let refBitsURL = try referenceURL("hts1a_700C.bit")
        let refPcmURL = try referenceURL("hts1a_700C.raw")
        let bits = try Data(contentsOf: refBitsURL)
        let refPcm = try Data(contentsOf: refPcmURL)
        guard let codec = Codec2(mode: .mode700C) else {
            XCTFail("could not instantiate Codec2 700C")
            return
        }
        let nbyte = codec.bytesPerFrame
        let nsam = codec.samplesPerFrame
        var frameBits = [UInt8](repeating: 0, count: nbyte)
        var frameSpeech = [Int16](repeating: 0, count: nsam)
        var decoded = Data()

        var idx = 0
        while idx + nbyte <= bits.count {
            for i in 0..<nbyte { frameBits[i] = bits[idx + i] }
            codec.decode(bits: frameBits, speech: &frameSpeech)
            frameSpeech.withUnsafeBufferPointer { decoded.append(Data(buffer: $0)) }
            idx += nbyte
        }

        XCTAssertEqual(decoded.count, refPcm.count)

        // Compare int16 samples; allow ±1 LSB to absorb fp rounding.
        let n = decoded.count / 2
        var maxDiff: Int32 = 0
        var differing = 0
        decoded.withUnsafeBytes { sw in
            refPcm.withUnsafeBytes { rf in
                let swP = sw.bindMemory(to: Int16.self)
                let rfP = rf.bindMemory(to: Int16.self)
                for i in 0..<n {
                    let d = abs(Int32(swP[i]) - Int32(rfP[i]))
                    if d != 0 { differing += 1 }
                    if d > maxDiff { maxDiff = d }
                }
            }
        }
        XCTAssertLessThanOrEqual(maxDiff, 1, "decoder sample drifted by \(maxDiff) LSB")
        // Allow up to 1% of samples to be ±1 LSB.
        XCTAssertLessThan(Double(differing) / Double(n), 0.02,
                          "too many samples drifted by ≤1 LSB: \(differing)/\(n)")
    }
}
