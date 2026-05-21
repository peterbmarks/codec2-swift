import Foundation

// Swift port of codec2/src/nlp.c.
//
// Non-Linear Pitch (NLP) estimation [Rowe PhD thesis, ch.4]. The decimating
// 16k→8k path is reachable from FreeDV but never from the Codec2-only API
// surface, so it is preserved here only as a precondition trap; we'll bring
// across the FIR resampler when/if FreeDV gets ported.

public final class NLPState {

    public static let pmaxM: Int = 320
    public static let coeff: Float = 0.95
    public static let peFFTSize: Int = 512
    public static let dec: Int = 5
    public static let sampleRate: Int = 8000
    public static let pi: Double = 3.141592654
    public static let cnlp: Float = 0.3
    public static let nlpNTap: Int = 48
    public static let f0Max: Int = 500

    /// 48-tap 600Hz low-pass FIR filter coefficients, copied verbatim from
    /// `nlp_fir[]` in nlp.c.
    public static let nlpFir: [Float] = [
        -1.0818124e-03, -1.1008344e-03, -9.2768838e-04, -4.2289438e-04,
         5.5034190e-04,  2.0029849e-03,  3.7058509e-03,  5.1449415e-03,
         5.5924666e-03,  4.3036754e-03,  8.0284511e-04, -4.8204610e-03,
        -1.1705810e-02, -1.8199275e-02, -2.2065282e-02, -2.0920610e-02,
        -1.2808831e-02,  3.2204775e-03,  2.6683811e-02,  5.5520624e-02,
         8.6305944e-02,  1.1480192e-01,  1.3674206e-01,  1.4867556e-01,
         1.4867556e-01,  1.3674206e-01,  1.1480192e-01,  8.6305944e-02,
         5.5520624e-02,  2.6683811e-02,  3.2204775e-03, -1.2808831e-02,
        -2.0920610e-02, -2.2065282e-02, -1.8199275e-02, -1.1705810e-02,
        -4.8204610e-03,  8.0284511e-04,  4.3036754e-03,  5.5924666e-03,
         5.1449415e-03,  3.7058509e-03,  2.0029849e-03,  5.5034190e-04,
        -4.2289438e-04, -9.2768838e-04, -1.1008344e-03, -1.0818124e-03,
    ]

    public var fs: Int
    public var m: Int
    public var w: [Float]
    public var sq: [Float]
    public var memX: Float = 0
    public var memY: Float = 0
    public var memFir: [Float]
    public var fftCfg: Codec2FFTConfig

    public init(c2const: C2Const) {
        precondition(c2const.fs == 8000, "NLP Fs=16000 path is FreeDV-only and not yet ported")
        self.fs = c2const.fs
        self.m = c2const.mPitch
        precondition(self.m <= NLPState.pmaxM)

        var w = [Float](repeating: 0, count: NLPState.pmaxM / NLPState.dec)
        let denom = Float((c2const.mPitch / NLPState.dec) - 1)
        for i in 0..<(c2const.mPitch / NLPState.dec) {
            // 0.5 - 0.5*cosf(2*PI*i/(m/DEC - 1)) with PI = 3.141592654
            w[i] = 0.5 - 0.5 * cosf(2.0 * Float(NLPState.pi) * Float(i) / denom)
        }
        self.w = w

        self.sq = [Float](repeating: 0, count: NLPState.pmaxM)
        self.memFir = [Float](repeating: 0, count: NLPState.nlpNTap)
        self.fftCfg = Codec2FFTConfig(nfft: NLPState.peFFTSize, inverse: false)
    }
}

public enum NLP {

    /// Pitch estimation. Returns F0 in Hz; writes estimated pitch period (in
    /// samples at the codec's current sample rate) to `pitch`. `prevF0` is
    /// read/written for pitch-tracking continuity.
    public static func run(state: NLPState,
                           sn: [Float],
                           n: Int,
                           pitch: inout Float,
                           prevF0: inout Float) -> Float {
        precondition(state.fs == 8000, "NLP Fs=16000 path not yet ported")
        let m = state.m
        let peN = NLPState.peFFTSize
        let dec = NLPState.dec

        // Square latest input samples.
        for i in (m - n)..<m {
            state.sq[i] = sn[i] * sn[i]
        }

        // Notch at DC. The "+1.0" trick on the last line matches the C source
        // and prevents a zero input from making kiss_fft pathologically slow.
        for i in (m - n)..<m {
            var notch: Float = state.sq[i] - state.memX
            notch += NLPState.coeff * state.memY
            state.memX = state.sq[i]
            state.memY = notch
            state.sq[i] = notch + 1.0
        }

        // 48-tap FIR low-pass.
        for i in (m - n)..<m {
            for j in 0..<(NLPState.nlpNTap - 1) {
                state.memFir[j] = state.memFir[j + 1]
            }
            state.memFir[NLPState.nlpNTap - 1] = state.sq[i]

            var acc: Float = 0
            for j in 0..<NLPState.nlpNTap {
                acc += state.memFir[j] * NLPState.nlpFir[j]
            }
            state.sq[i] = acc
        }

        // Decimate and DFT.
        var fw = [COMP](repeating: COMP(), count: peN)
        for i in 0..<(m / dec) {
            fw[i] = COMP(state.sq[i * dec] * state.w[i], 0)
        }
        Codec2FFT.fftInPlace(state.fftCfg, buffer: &fw)
        for i in 0..<peN {
            fw[i] = COMP(fw[i].real * fw[i].real + fw[i].imag * fw[i].imag, 0)
        }

        // Pitch search range.
        let pmin = Int(floor(Double(NLPState.sampleRate) * Codec2Constants.pMinSeconds))
        let pmax = Int(floor(Double(NLPState.sampleRate) * Codec2Constants.pMaxSeconds))

        var gmax: Float = 0
        var gmaxBin: Int = peN * dec / pmax
        let lo = peN * dec / pmax
        let hi = peN * dec / pmin
        for i in lo...hi {
            if fw[i].real > gmax {
                gmax = fw[i].real
                gmaxBin = i
            }
        }

        let bestF0 = postProcessSubMultiples(fw: fw, pmin: pmin, pmax: pmax,
                                             gmax: gmax, gmaxBin: gmaxBin,
                                             prevF0: prevF0)

        // Shift sq[] left by n samples to make room for the next frame.
        for i in 0..<(m - n) {
            state.sq[i] = state.sq[i + n]
        }

        pitch = Float(state.fs) / bestF0
        prevF0 = bestF0
        _ = m   // m is read but not mutated; silence the unused warning explicitly
        return bestF0
    }

    /// Submultiple post-processor — searches gmax_bin/2, /3, /4 ... for local
    /// maxima that exceed an empirical threshold and prefers lower F0 in noise.
    public static func postProcessSubMultiples(fw: [COMP],
                                               pmin: Int, pmax: Int,
                                               gmax: Float, gmaxBin: Int,
                                               prevF0: Float) -> Float {
        let peN = NLPState.peFFTSize
        let dec = NLPState.dec
        var mult = 2
        let minBin = peN * dec / pmax
        var cmaxBin = gmaxBin
        let prevF0Bin = Int(prevF0 * Float(peN * dec) / Float(NLPState.sampleRate))

        while gmaxBin / mult >= minBin {
            let bCenter = gmaxBin / mult
            var bmin = Int(0.8 * Float(bCenter))
            let bmax = Int(1.2 * Float(bCenter))
            if bmin < minBin { bmin = minBin }

            let thresh: Float
            if (prevF0Bin > bmin) && (prevF0Bin < bmax) {
                thresh = NLPState.cnlp * 0.5 * gmax
            } else {
                thresh = NLPState.cnlp * gmax
            }

            var lmax: Float = 0
            var lmaxBin = bmin
            for b in bmin...bmax {
                if fw[b].real > lmax {
                    lmax = fw[b].real
                    lmaxBin = b
                }
            }

            if lmax > thresh {
                if lmax > fw[lmaxBin - 1].real && lmax > fw[lmaxBin + 1].real {
                    cmaxBin = lmaxBin
                }
            }
            mult += 1
        }

        return Float(cmaxBin) * Float(NLPState.sampleRate) / Float(peN * dec)
    }
}
