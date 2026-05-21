import XCTest
@testable import Codec2

final class KissFFTRTests: XCTestCase {

    /// Real-FFT followed by real-IFFT should return the original signal
    /// scaled by N. This is the canonical KissFFTR round-trip identity.
    func testRealForwardInverseRoundTrip() {
        let n = 512
        let fwd = KissFFTRConfig(nfft: n, inverse: false)
        let inv = KissFFTRConfig(nfft: n, inverse: true)

        var time = [Float](repeating: 0, count: n)
        for i in 0..<n {
            let t = Float(i) / Float(n)
            time[i] = sinf(2 * .pi * 5 * t) + 0.4 * cosf(2 * .pi * 23 * t)
        }
        var freq = [COMP](repeating: COMP(), count: n / 2 + 1)
        KissFFTR.forward(fwd, timedata: time, freqdata: &freq)

        var back = [Float](repeating: 0, count: n)
        KissFFTR.inverse(inv, freqdata: freq, timedata: &back)

        let invN: Float = 1.0 / Float(n)
        for i in 0..<n {
            XCTAssertEqual(back[i] * invN, time[i], accuracy: 2e-4,
                           "index \(i) drift after FFTR round trip")
        }
    }

    /// Real FFT of a cosine at bin k peaks symmetrically: the standard test
    /// that catches sign/scaling/twiddle-table mistakes in the real-FFT
    /// rearrangement step.
    func testRealCosineBin() {
        let n = 64
        let k = 6
        let cfg = KissFFTRConfig(nfft: n, inverse: false)
        var time = [Float](repeating: 0, count: n)
        for i in 0..<n {
            time[i] = cosf(2 * .pi * Float(k) * Float(i) / Float(n))
        }
        var freq = [COMP](repeating: COMP(), count: n / 2 + 1)
        KissFFTR.forward(cfg, timedata: time, freqdata: &freq)

        for i in 0...(n / 2) {
            let expected: Float = (i == k) ? Float(n) / 2 : 0
            XCTAssertEqual(freq[i].real, expected, accuracy: 1e-3, "bin \(i) real")
            XCTAssertEqual(freq[i].imag, 0,         accuracy: 1e-3, "bin \(i) imag")
        }
    }
}
