import Foundation

// Swift port of codec2/src/kiss_fft.c (Mark Borgerding's KissFFT). Codec2
// builds KissFFT in float mode so we only port the float path. Macros are
// expanded inline; pointer arithmetic is replaced with explicit indices.

public final class KissFFTConfig {
    public let nfft: Int
    public let inverse: Bool
    /// nfft twiddle factors. The C code uses `twiddles[1]` as a flexible array
    /// member; we just hold them all in a Swift array.
    public var twiddles: [COMP]
    public var factors: [Int]    // pairs of (p, m), up to 2*MAXFACTORS entries

    public static let maxFactors = 32

    public init(nfft: Int, inverse: Bool) {
        self.nfft = nfft
        self.inverse = inverse
        self.twiddles = [COMP](repeating: COMP(), count: nfft)
        self.factors = [Int](repeating: 0, count: 2 * KissFFTConfig.maxFactors)

        // High-precision twiddles, matching the C double-precision pi.
        let pi = 3.141592653589793238462643383279502884197169399375105820974944
        for i in 0..<nfft {
            var phase = -2.0 * pi * Double(i) / Double(nfft)
            if inverse { phase = -phase }
            // cosf/sinf semantics: convert to Float at the same point as the C code.
            self.twiddles[i] = COMP(Float(cos(phase)), Float(sin(phase)))
        }
        KissFFTConfig.factor(nfft, into: &self.factors)
    }

    /// Mirrors `kf_factor` in the C code: pack as p1,m1,p2,m2,... where
    /// p_i * m_i = m_{i-1}.
    private static func factor(_ n: Int, into facbuf: inout [Int]) {
        var n = n
        var p = 4
        let floorSqrt = Int(floorf(sqrtf(Float(n))))
        var idx = 0
        repeat {
            while n % p != 0 {
                switch p {
                case 4: p = 2
                case 2: p = 3
                default: p += 2
                }
                if p > floorSqrt { p = n }   // no more factors, skip to end
            }
            n /= p
            facbuf[idx] = p; idx += 1
            facbuf[idx] = n; idx += 1
        } while n > 1
    }
}

public enum KissFFT {

    // MARK: - Public API

    public static func fft(_ cfg: KissFFTConfig, _ fin: [COMP], _ fout: inout [COMP]) {
        fftStride(cfg, fin: fin, fout: &fout, inStride: 1)
    }

    public static func fftStride(_ cfg: KissFFTConfig,
                                 fin: [COMP],
                                 fout: inout [COMP],
                                 inStride: Int) {
        // KissFFT is not truly in-place; the C code copies through a scratch
        // buffer when fin and fout alias. We don't expose aliasing in Swift,
        // so we always work straight from `fin` into `fout`.
        var tmp = fout
        kfWork(out: &tmp, outOffset: 0,
               f: fin, fOffset: 0,
               fstride: 1, inStride: inStride,
               factors: cfg.factors, factorsOffset: 0,
               st: cfg)
        fout = tmp
    }

    public static func nextFastSize(_ n: Int) -> Int {
        var n = n
        while true {
            var m = n
            while m % 2 == 0 { m /= 2 }
            while m % 3 == 0 { m /= 3 }
            while m % 5 == 0 { m /= 5 }
            if m <= 1 { break }
            n += 1
        }
        return n
    }

    // MARK: - Butterflies

    private static func kfBfly2(out: inout [COMP], outOffset: Int,
                                fstride: Int, st: KissFFTConfig, m mIn: Int) {
        var m = mIn
        var fOut = outOffset
        var fOut2 = outOffset + m
        var tw1Idx = 0
        repeat {
            // C_MUL(t, out2, tw1)
            let a = out[fOut2]
            let b = st.twiddles[tw1Idx]
            let t = COMP(a.real * b.real - a.imag * b.imag,
                         a.real * b.imag + a.imag * b.real)
            tw1Idx += fstride
            // C_SUB(*Fout2, *Fout, t); C_ADDTO(*Fout, t)
            out[fOut2] = COMP(out[fOut].real - t.real, out[fOut].imag - t.imag)
            out[fOut] = COMP(out[fOut].real + t.real, out[fOut].imag + t.imag)
            fOut2 += 1
            fOut += 1
            m -= 1
        } while m != 0
    }

    private static func kfBfly4(out: inout [COMP], outOffset: Int,
                                fstride: Int, st: KissFFTConfig, m: Int) {
        var k = m
        let m2 = 2 * m
        let m3 = 3 * m
        var tw1 = 0, tw2 = 0, tw3 = 0
        var fOut = outOffset
        var scratch = [COMP](repeating: COMP(), count: 6)

        repeat {
            // scratch[0] = Fout[m] * tw1
            let am = out[fOut + m]
            let tw1c = st.twiddles[tw1]
            scratch[0] = COMP(am.real * tw1c.real - am.imag * tw1c.imag,
                              am.real * tw1c.imag + am.imag * tw1c.real)
            // scratch[1] = Fout[m2] * tw2
            let a2 = out[fOut + m2]
            let tw2c = st.twiddles[tw2]
            scratch[1] = COMP(a2.real * tw2c.real - a2.imag * tw2c.imag,
                              a2.real * tw2c.imag + a2.imag * tw2c.real)
            // scratch[2] = Fout[m3] * tw3
            let a3 = out[fOut + m3]
            let tw3c = st.twiddles[tw3]
            scratch[2] = COMP(a3.real * tw3c.real - a3.imag * tw3c.imag,
                              a3.real * tw3c.imag + a3.imag * tw3c.real)

            // scratch[5] = Fout - scratch[1]
            scratch[5] = COMP(out[fOut].real - scratch[1].real,
                              out[fOut].imag - scratch[1].imag)
            // Fout += scratch[1]
            out[fOut] = COMP(out[fOut].real + scratch[1].real,
                             out[fOut].imag + scratch[1].imag)
            // scratch[3] = scratch[0] + scratch[2]; scratch[4] = scratch[0] - scratch[2]
            scratch[3] = COMP(scratch[0].real + scratch[2].real,
                              scratch[0].imag + scratch[2].imag)
            scratch[4] = COMP(scratch[0].real - scratch[2].real,
                              scratch[0].imag - scratch[2].imag)
            // Fout[m2] = Fout - scratch[3]
            out[fOut + m2] = COMP(out[fOut].real - scratch[3].real,
                                  out[fOut].imag - scratch[3].imag)
            tw1 += fstride
            tw2 += fstride * 2
            tw3 += fstride * 3
            // Fout += scratch[3]
            out[fOut] = COMP(out[fOut].real + scratch[3].real,
                             out[fOut].imag + scratch[3].imag)

            if st.inverse {
                out[fOut + m]  = COMP(scratch[5].real - scratch[4].imag,
                                      scratch[5].imag + scratch[4].real)
                out[fOut + m3] = COMP(scratch[5].real + scratch[4].imag,
                                      scratch[5].imag - scratch[4].real)
            } else {
                out[fOut + m]  = COMP(scratch[5].real + scratch[4].imag,
                                      scratch[5].imag - scratch[4].real)
                out[fOut + m3] = COMP(scratch[5].real - scratch[4].imag,
                                      scratch[5].imag + scratch[4].real)
            }
            fOut += 1
            k -= 1
        } while k != 0
    }

    private static func kfBfly3(out: inout [COMP], outOffset: Int,
                                fstride: Int, st: KissFFTConfig, m: Int) {
        var k = m
        let m2 = 2 * m
        let epi3 = st.twiddles[fstride * m]
        var tw1 = 0, tw2 = 0
        var fOut = outOffset
        var scratch = [COMP](repeating: COMP(), count: 5)

        repeat {
            // scratch[1] = Fout[m] * tw1
            let am = out[fOut + m]
            let tw1c = st.twiddles[tw1]
            scratch[1] = COMP(am.real * tw1c.real - am.imag * tw1c.imag,
                              am.real * tw1c.imag + am.imag * tw1c.real)
            // scratch[2] = Fout[m2] * tw2
            let a2 = out[fOut + m2]
            let tw2c = st.twiddles[tw2]
            scratch[2] = COMP(a2.real * tw2c.real - a2.imag * tw2c.imag,
                              a2.real * tw2c.imag + a2.imag * tw2c.real)
            // scratch[3] = scratch[1] + scratch[2], scratch[0] = scratch[1] - scratch[2]
            scratch[3] = COMP(scratch[1].real + scratch[2].real,
                              scratch[1].imag + scratch[2].imag)
            scratch[0] = COMP(scratch[1].real - scratch[2].real,
                              scratch[1].imag - scratch[2].imag)
            tw1 += fstride
            tw2 += fstride * 2

            // Fout[m].r = Fout.r - 0.5*scratch[3].r
            let foutR = out[fOut].real
            let foutI = out[fOut].imag
            out[fOut + m] = COMP(foutR - 0.5 * scratch[3].real,
                                 foutI - 0.5 * scratch[3].imag)

            // C_MULBYSCALAR(scratch[0], epi3.i)
            scratch[0] = COMP(scratch[0].real * epi3.imag,
                              scratch[0].imag * epi3.imag)

            // C_ADDTO(*Fout, scratch[3])
            out[fOut] = COMP(out[fOut].real + scratch[3].real,
                             out[fOut].imag + scratch[3].imag)

            // Fout[m2].r = Fout[m].r + scratch[0].i; .i = Fout[m].i - scratch[0].r
            out[fOut + m2] = COMP(out[fOut + m].real + scratch[0].imag,
                                  out[fOut + m].imag - scratch[0].real)

            // Fout[m].r -= scratch[0].i; .i += scratch[0].r
            out[fOut + m] = COMP(out[fOut + m].real - scratch[0].imag,
                                 out[fOut + m].imag + scratch[0].real)

            fOut += 1
            k -= 1
        } while k != 0
    }

    private static func kfBfly5(out: inout [COMP], outOffset: Int,
                                fstride: Int, st: KissFFTConfig, m: Int) {
        let ya = st.twiddles[fstride * m]
        let yb = st.twiddles[fstride * 2 * m]
        var scratch = [COMP](repeating: COMP(), count: 13)

        let f0Base = outOffset
        let f1Base = outOffset + m
        let f2Base = outOffset + 2 * m
        let f3Base = outOffset + 3 * m
        let f4Base = outOffset + 4 * m

        for u in 0..<m {
            scratch[0] = out[f0Base + u]

            // scratch[1] = *Fout1 * tw[u*fstride]
            let s1Tw = st.twiddles[u * fstride]
            let f1v = out[f1Base + u]
            scratch[1] = COMP(f1v.real * s1Tw.real - f1v.imag * s1Tw.imag,
                              f1v.real * s1Tw.imag + f1v.imag * s1Tw.real)
            let s2Tw = st.twiddles[2 * u * fstride]
            let f2v = out[f2Base + u]
            scratch[2] = COMP(f2v.real * s2Tw.real - f2v.imag * s2Tw.imag,
                              f2v.real * s2Tw.imag + f2v.imag * s2Tw.real)
            let s3Tw = st.twiddles[3 * u * fstride]
            let f3v = out[f3Base + u]
            scratch[3] = COMP(f3v.real * s3Tw.real - f3v.imag * s3Tw.imag,
                              f3v.real * s3Tw.imag + f3v.imag * s3Tw.real)
            let s4Tw = st.twiddles[4 * u * fstride]
            let f4v = out[f4Base + u]
            scratch[4] = COMP(f4v.real * s4Tw.real - f4v.imag * s4Tw.imag,
                              f4v.real * s4Tw.imag + f4v.imag * s4Tw.real)

            // 7 = 1+4, 10 = 1-4, 8 = 2+3, 9 = 2-3
            scratch[7]  = COMP(scratch[1].real + scratch[4].real,
                               scratch[1].imag + scratch[4].imag)
            scratch[10] = COMP(scratch[1].real - scratch[4].real,
                               scratch[1].imag - scratch[4].imag)
            scratch[8]  = COMP(scratch[2].real + scratch[3].real,
                               scratch[2].imag + scratch[3].imag)
            scratch[9]  = COMP(scratch[2].real - scratch[3].real,
                               scratch[2].imag - scratch[3].imag)

            // *Fout0 += scratch[7] + scratch[8]
            out[f0Base + u] = COMP(
                out[f0Base + u].real + scratch[7].real + scratch[8].real,
                out[f0Base + u].imag + scratch[7].imag + scratch[8].imag
            )

            scratch[5] = COMP(
                scratch[0].real + scratch[7].real * ya.real + scratch[8].real * yb.real,
                scratch[0].imag + scratch[7].imag * ya.real + scratch[8].imag * yb.real
            )
            scratch[6] = COMP(
                 scratch[10].imag * ya.imag + scratch[9].imag * yb.imag,
                -scratch[10].real * ya.imag - scratch[9].real * yb.imag
            )

            out[f1Base + u] = COMP(scratch[5].real - scratch[6].real,
                                   scratch[5].imag - scratch[6].imag)
            out[f4Base + u] = COMP(scratch[5].real + scratch[6].real,
                                   scratch[5].imag + scratch[6].imag)

            scratch[11] = COMP(
                scratch[0].real + scratch[7].real * yb.real + scratch[8].real * ya.real,
                scratch[0].imag + scratch[7].imag * yb.real + scratch[8].imag * ya.real
            )
            scratch[12] = COMP(
                -scratch[10].imag * yb.imag + scratch[9].imag * ya.imag,
                 scratch[10].real * yb.imag - scratch[9].real * ya.imag
            )

            out[f2Base + u] = COMP(scratch[11].real + scratch[12].real,
                                   scratch[11].imag + scratch[12].imag)
            out[f3Base + u] = COMP(scratch[11].real - scratch[12].real,
                                   scratch[11].imag - scratch[12].imag)
        }
    }

    private static func kfBflyGeneric(out: inout [COMP], outOffset: Int,
                                      fstride: Int, st: KissFFTConfig,
                                      m: Int, p: Int) {
        let nOrig = st.nfft
        var scratch = [COMP](repeating: COMP(), count: p)

        for u in 0..<m {
            var k = u
            for q1 in 0..<p {
                scratch[q1] = out[outOffset + k]
                k += m
            }
            k = u
            for q1 in 0..<p {
                var twidx = 0
                out[outOffset + k] = scratch[0]
                for q in 1..<p {
                    twidx += fstride * k
                    if twidx >= nOrig { twidx -= nOrig }
                    let tw = st.twiddles[twidx]
                    let t = COMP(scratch[q].real * tw.real - scratch[q].imag * tw.imag,
                                 scratch[q].real * tw.imag + scratch[q].imag * tw.real)
                    out[outOffset + k] = COMP(out[outOffset + k].real + t.real,
                                              out[outOffset + k].imag + t.imag)
                }
                k += m
                _ = q1
            }
        }
    }

    // MARK: - Recursive worker

    private static func kfWork(out: inout [COMP], outOffset: Int,
                               f: [COMP], fOffset: Int,
                               fstride: Int, inStride: Int,
                               factors: [Int], factorsOffset: Int,
                               st: KissFFTConfig) {
        let p = factors[factorsOffset]       // radix
        let m = factors[factorsOffset + 1]   // stage's fft length / p

        if m == 1 {
            // copy decimated input directly
            for i in 0..<p {
                out[outOffset + i] = f[fOffset + i * fstride * inStride]
            }
        } else {
            // recursive call: DFT of size m*p = p instances of size m
            for i in 0..<p {
                kfWork(out: &out, outOffset: outOffset + i * m,
                       f: f, fOffset: fOffset + i * fstride * inStride,
                       fstride: fstride * p, inStride: inStride,
                       factors: factors, factorsOffset: factorsOffset + 2,
                       st: st)
            }
        }

        // recombine the p smaller DFTs
        switch p {
        case 2: kfBfly2(out: &out, outOffset: outOffset, fstride: fstride, st: st, m: m)
        case 3: kfBfly3(out: &out, outOffset: outOffset, fstride: fstride, st: st, m: m)
        case 4: kfBfly4(out: &out, outOffset: outOffset, fstride: fstride, st: st, m: m)
        case 5: kfBfly5(out: &out, outOffset: outOffset, fstride: fstride, st: st, m: m)
        default: kfBflyGeneric(out: &out, outOffset: outOffset,
                               fstride: fstride, st: st, m: m, p: p)
        }
    }
}
