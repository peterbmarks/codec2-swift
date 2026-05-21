import XCTest
@testable import Codec2

final class KissFFTTests: XCTestCase {

    /// KissFFT is unnormalised: forward followed by inverse multiplies by N.
    /// This verifies the Swift implementation's algebra against that property
    /// for the FFT sizes used by codec2 (FFT_ENC = FFT_DEC = 512).
    func testForwardInverseRoundTrip() {
        let n = 512
        let fwd = KissFFTConfig(nfft: n, inverse: false)
        let inv = KissFFTConfig(nfft: n, inverse: true)

        var input = [COMP](repeating: COMP(), count: n)
        for i in 0..<n {
            let t = Float(i) / Float(n)
            input[i] = COMP(sinf(2 * .pi * 5 * t) + 0.3 * sinf(2 * .pi * 21 * t),
                            0.0)
        }

        var spec = [COMP](repeating: COMP(), count: n)
        KissFFT.fft(fwd, input, &spec)

        var back = [COMP](repeating: COMP(), count: n)
        KissFFT.fft(inv, spec, &back)

        let inv_n = 1.0 / Float(n)
        for i in 0..<n {
            XCTAssertEqual(back[i].real * inv_n, input[i].real, accuracy: 1e-4)
            XCTAssertEqual(back[i].imag * inv_n, input[i].imag, accuracy: 1e-4)
        }
    }

    /// FFT of a unit impulse is a constant 1 across all bins. Independent of
    /// nfft factorisation it exercises the radix-2/4 mixed-radix path.
    func testImpulseSpectrum() {
        let n = 128
        let cfg = KissFFTConfig(nfft: n, inverse: false)
        var input = [COMP](repeating: COMP(), count: n)
        input[0] = COMP(1, 0)
        var spec = [COMP](repeating: COMP(), count: n)
        KissFFT.fft(cfg, input, &spec)
        for i in 0..<n {
            XCTAssertEqual(spec[i].real, 1.0, accuracy: 1e-5, "bin \(i) real")
            XCTAssertEqual(spec[i].imag, 0.0, accuracy: 1e-5, "bin \(i) imag")
        }
    }

    /// A pure complex sinusoid e^{j2πkn/N} has all energy at bin k.
    func testSingleBinSinusoid() {
        let n = 64
        let k = 7
        let cfg = KissFFTConfig(nfft: n, inverse: false)
        var input = [COMP](repeating: COMP(), count: n)
        for i in 0..<n {
            let phase = 2 * Float.pi * Float(k) * Float(i) / Float(n)
            input[i] = COMP(cosf(phase), sinf(phase))
        }
        var spec = [COMP](repeating: COMP(), count: n)
        KissFFT.fft(cfg, input, &spec)
        for i in 0..<n {
            let expected: Float = (i == k) ? Float(n) : 0
            XCTAssertEqual(spec[i].real, expected, accuracy: 1e-3, "bin \(i) real")
            XCTAssertEqual(spec[i].imag, 0.0, accuracy: 1e-3, "bin \(i) imag")
        }
    }
}
