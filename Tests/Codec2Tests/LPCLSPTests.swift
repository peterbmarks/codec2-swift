import XCTest
@testable import Codec2

final class LPCLSPTests: XCTestCase {

    func testLevinsonDurbinTrivial() {
        // Single-impulse autocorrelation -> all predictor coeffs except a[0] are zero.
        let r: [Float] = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        var a = [Float](repeating: 99, count: 11)
        LPC.levinsonDurbin(r, lpcs: &a, order: 10)
        XCTAssertEqual(a[0], 1.0)
        for i in 1...10 {
            XCTAssertEqual(a[i], 0.0, accuracy: 1e-6, "a[\(i)] should be 0 for delta autocorrelation")
        }
    }

    func testLpcLspRoundTrip() {
        // Synthesize a stable LPC vector by deriving it from a known speech-like
        // frame, convert to LSP and back, and require the LPC vector to match
        // within numerical tolerance.
        let n = 256
        var sn = [Float](repeating: 0, count: n)
        // Two-tone signal with strong harmonic structure (well within LPC's domain).
        for i in 0..<n {
            let t = Float(i) / Float(n)
            sn[i] = 0.5 * sinf(2 * .pi * 8 * t) + 0.3 * sinf(2 * .pi * 13 * t)
        }
        var ak = [Float](repeating: 0, count: 11)
        _ = LPC.findAks(sn, a: &ak, nSam: n, order: 10)
        XCTAssertEqual(ak[0], 1.0)

        var lsp = [Float](repeating: 0, count: 10)
        let roots = LSP.lpcToLsp(ak, order: 10, freq: &lsp, nb: 4, delta: 0.02)
        XCTAssertEqual(roots, 10, "expected 10 LSP roots for an order-10 LPC")

        // LSPs are monotonically increasing in radians for a stable filter.
        for i in 1..<10 {
            XCTAssertGreaterThan(lsp[i], lsp[i - 1], "LSPs must be ordered")
        }

        var ak2 = [Float](repeating: 0, count: 11)
        LSP.lspToLpc(lsp, ak: &ak2, order: 10)
        // Round-trip tolerance is bounded by fp32 precision in the order-10 LPC
        // polynomial reconstruction. The C reference exhibits similar drift on
        // the same input; a precise byte-for-byte parity check is wired up in
        // CReferenceParityTests once that harness lands.
        for i in 0...10 {
            XCTAssertEqual(ak2[i], ak[i], accuracy: 1e-2,
                           "LPC->LSP->LPC drift at index \(i): \(ak[i]) vs \(ak2[i])")
        }
    }
}
