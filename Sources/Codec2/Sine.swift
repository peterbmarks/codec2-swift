import Foundation

// Swift port of codec2/src/sine.c. Sinusoidal speech analysis and synthesis:
// time-domain windowing, DFT of speech, harmonic-sum pitch refinement,
// amplitude/phase estimation, MBE voicing decision, and trapezoidal-window
// overlap-add synthesis.

public enum Sine {

    public static let hpfBeta: Float = 0.125
    public static let fftEnc: Int = Codec2Constants.fftEnc      // 512
    public static let fftDec: Int = Codec2Constants.fftDec      // 512

    /// Builds the runtime constants table. Matches `c2const_create` exactly.
    public static func c2constCreate(fs: Int, framelengthS: Float) -> C2Const {
        precondition(fs == 8000 || fs == 16000)
        let nSamp = Int((Float(fs) * framelengthS).rounded())
        let pMin = Int(floor(Double(fs) * Codec2Constants.pMinSeconds))
        let pMax = Int(floor(Double(fs) * Codec2Constants.pMaxSeconds))
        let mPitch = Int(floor(Double(fs) * Codec2Constants.mPitchSeconds))
        let maxAmp = Int(floor(Double(fs) * Codec2Constants.pMaxSeconds / 2.0))
        let woMin = Float(Codec2Constants.twoPi / Double(pMax))
        let woMax = Float(Codec2Constants.twoPi / Double(pMin))
        let nw = (fs == 8000) ? 279 : 511
        let tw = Int(Double(fs) * Codec2Constants.twSeconds)

        return C2Const(fs: fs, nSamp: nSamp, maxAmp: maxAmp,
                       mPitch: mPitch, pMin: pMin, pMax: pMax,
                       woMin: woMin, woMax: woMax, nw: nw, tw: tw)
    }

    /// Generates the centred Hamming analysis window `w[]` and its frequency
    /// domain version `W[]`. Mirrors `make_analysis_window` with the same
    /// FFT-shift trick so W ends up real.
    public static func makeAnalysisWindow(c2const: C2Const,
                                          fftFwdCfg: Codec2FFTConfig,
                                          w: inout [Float],
                                          W: inout [Float]) {
        let mPitch = c2const.mPitch
        let nw = c2const.nw

        for i in 0..<(mPitch / 2 - nw / 2) { w[i] = 0 }
        var msum: Float = 0
        var j = 0
        for i in (mPitch / 2 - nw / 2)..<(mPitch / 2 + nw / 2) {
            w[i] = 0.5 - 0.5 * cosf(Float(Codec2Constants.twoPi) * Float(j) / Float(nw - 1))
            msum += w[i] * w[i]
            j += 1
        }
        for i in (mPitch / 2 + nw / 2)..<mPitch { w[i] = 0 }

        // Normalise so that frequency-domain amplitude estimation is direct.
        let scale = 1.0 / sqrtf(msum * Float(fftEnc))
        for i in 0..<mPitch { w[i] *= scale }

        // Build a circularly-shifted copy of w to make W[] real.
        var wshift = [COMP](repeating: COMP(), count: fftEnc)
        for i in 0..<(nw / 2) { wshift[i] = COMP(w[i + mPitch / 2], 0) }
        var src = mPitch / 2 - nw / 2
        for i in (fftEnc - nw / 2)..<fftEnc {
            wshift[i] = COMP(w[src], 0)
            src += 1
        }

        var temp = [COMP](repeating: COMP(), count: fftEnc)
        Codec2FFT.fft(fftFwdCfg, input: wshift, output: &temp)

        // Re-arrange to be symmetric about FFT_ENC/2.
        for i in 0..<(fftEnc / 2) {
            W[i] = temp[i + fftEnc / 2].real
            W[i + fftEnc / 2] = temp[i].real
        }
    }

    /// Single-pole high-pass filter, -3 dB at ~160 Hz.
    /// `states[0]` = y(n-1), `states[1]` = x(n-1).
    public static func hpf(_ x: Float, states: inout [Float]) -> Float {
        states[0] = -hpfBeta * states[0] + x - states[1]
        states[1] = x
        return states[0]
    }

    /// Forward DFT of the windowed centre of `Sn`. Output is the complex
    /// spectrum in `Sw`.
    public static func dftSpeech(c2const: C2Const,
                                 fftFwdCfg: Codec2FFTConfig,
                                 sw: inout [COMP],
                                 sn: [Float], w: [Float]) {
        let mPitch = c2const.mPitch
        let nw = c2const.nw

        for i in 0..<fftEnc { sw[i] = COMP() }
        // 2nd half of analysis window -> start of FFT buffer
        for i in 0..<(nw / 2) {
            let idx = i + mPitch / 2
            sw[i] = COMP(sn[idx] * w[idx], 0)
        }
        // 1st half -> end of FFT buffer (negative-time wrap)
        var srcIdx = mPitch / 2 - nw / 2
        for i in (fftEnc - nw / 2)..<fftEnc {
            sw[i] = COMP(sn[srcIdx] * w[srcIdx], 0)
            srcIdx += 1
        }
        Codec2FFT.fftInPlace(fftFwdCfg, buffer: &sw)
    }

    /// Two-stage harmonic-sum pitch refinement (coarse ±5 samples step 1.0,
    /// fine ±1 sample step 0.25). Sets `model.wo` and `model.l`.
    public static func twoStagePitchRefinement(c2const: C2Const,
                                               model: inout Codec2Model,
                                               sw: [COMP]) {
        let pi: Float = Float(Codec2Constants.pi)
        let twoPi: Float = Float(Codec2Constants.twoPi)

        var pmax = twoPi / model.wo + 5
        var pmin = twoPi / model.wo - 5
        var pstep: Float = 1.0
        hsPitchRefinement(model: &model, sw: sw, pmin: pmin, pmax: pmax, pstep: pstep)

        pmax = twoPi / model.wo + 1
        pmin = twoPi / model.wo - 1
        pstep = 0.25
        hsPitchRefinement(model: &model, sw: sw, pmin: pmin, pmax: pmax, pstep: pstep)

        if model.wo < twoPi / Float(c2const.pMax) { model.wo = twoPi / Float(c2const.pMax) }
        if model.wo > twoPi / Float(c2const.pMin) { model.wo = twoPi / Float(c2const.pMin) }

        model.l = Int(floorf(pi / model.wo))
        if model.wo * Float(model.l) >= 0.95 * pi {
            model.l -= 1
        }
        precondition(model.wo * Float(model.l) < pi)
    }

    /// Harmonic-sum pitch refinement step. Searches `Wo` in [pmin..pmax] at
    /// step `pstep` and keeps the value that maximises the harmonic energy
    /// sum.
    public static func hsPitchRefinement(model: inout Codec2Model,
                                         sw: [COMP],
                                         pmin: Float, pmax: Float, pstep: Float) {
        let pi: Float = Float(Codec2Constants.pi)
        let twoPi: Float = Float(Codec2Constants.twoPi)

        model.l = Int(pi / model.wo)
        var wom = model.wo
        var em: Float = 0
        let r = twoPi / Float(fftEnc)
        let oneOnR = 1.0 / r

        var p = pmin
        while p <= pmax {
            var e: Float = 0
            let wo = twoPi / p
            let bFloat = wo * oneOnR
            var current = bFloat
            for _ in 1...max(model.l, 1) {
                let b = Int(current + 0.5)
                if b < sw.count {
                    e += sw[b].real * sw[b].real + sw[b].imag * sw[b].imag
                }
                current += bFloat
            }
            if e > em { em = e; wom = wo }
            p += pstep
        }
        model.wo = wom
    }

    /// Estimates harmonic amplitudes (and phases if `estPhase != 0`) into
    /// `model.a[]` and `model.phi[]`.
    public static func estimateAmplitudes(model: inout Codec2Model,
                                          sw: [COMP], W: [Float],
                                          estPhase: Int) {
        let r = Float(Codec2Constants.twoPi) / Float(fftEnc)
        let oneOnR = 1.0 / r
        for m in 1...model.l {
            let am = Int((Float(m) - 0.5) * model.wo * oneOnR + 0.5)
            let bm = Int((Float(m) + 0.5) * model.wo * oneOnR + 0.5)
            var den: Float = 0
            for i in am..<bm {
                den += sw[i].real * sw[i].real + sw[i].imag * sw[i].imag
            }
            model.a[m] = sqrtf(den)
            if estPhase != 0 {
                let b = Int(Float(m) * model.wo / r + 0.5)
                model.phi[m] = atan2f(sw[b].imag, sw[b].real)
            }
            _ = W   // W used only by est_voicing_mbe; reference kept for parity
        }
    }

    /// MBE voicing decision. Returns SNR in dB; also sets `model.voiced`.
    public static func estVoicingMBE(c2const: C2Const,
                                     model: inout Codec2Model,
                                     sw: [COMP], W: [Float]) -> Float {
        let twoPi: Float = Float(Codec2Constants.twoPi)
        let l1000 = model.l * 1000 / (c2const.fs / 2)
        var sig: Float = 1e-4
        if l1000 >= 1 {
            for l in 1...l1000 { sig += model.a[l] * model.a[l] }
        }

        let wo = model.wo
        var error: Float = 1e-4
        if l1000 >= 1 {
            for l in 1...l1000 {
                var amR: Float = 0
                var amI: Float = 0
                var den: Float = 0
                let al = Int(ceilf((Float(l) - 0.5) * wo * Float(fftEnc) / twoPi))
                let bl = Int(ceilf((Float(l) + 0.5) * wo * Float(fftEnc) / twoPi))
                let offset = Int(Float(fftEnc) / 2 - Float(l) * wo * Float(fftEnc) / twoPi + 0.5)
                for m in al..<bl {
                    amR += sw[m].real * W[offset + m]
                    amI += sw[m].imag * W[offset + m]
                    den += W[offset + m] * W[offset + m]
                }
                amR /= den
                amI /= den
                for m in al..<bl {
                    let ewR = sw[m].real - amR * W[offset + m]
                    let ewI = sw[m].imag - amI * W[offset + m]
                    error += ewR * ewR + ewI * ewI
                }
            }
        }

        let snr = 10.0 * log10f(sig / error)
        model.voiced = snr > Codec2Constants.vThresh ? 1 : 0

        // Low/high band energy ratio post-correction.
        let l2000 = model.l * 2000 / (c2const.fs / 2)
        let l4000 = model.l * 4000 / (c2const.fs / 2)
        var elow: Float = 1e-4
        var ehigh: Float = 1e-4
        if l2000 >= 1 {
            for l in 1...l2000 { elow += model.a[l] * model.a[l] }
        }
        if l4000 >= l2000 {
            for l in l2000...l4000 { ehigh += model.a[l] * model.a[l] }
        }
        let eratio = 10.0 * log10f(elow / ehigh)

        if model.voiced == 0 {
            if eratio > 10.0 { model.voiced = 1 }
        }
        if model.voiced == 1 {
            if eratio < -10.0 { model.voiced = 0 }
            let sixty: Float = 60.0 * twoPi / Float(c2const.fs)
            if eratio < -4.0 && model.wo <= sixty { model.voiced = 0 }
        }
        return snr
    }

    /// Builds the trapezoidal (Parzen) synthesis window `Pn` of length 2*nSamp.
    public static func makeSynthesisWindow(c2const: C2Const, pn: inout [Float]) {
        let nSamp = c2const.nSamp
        let tw = c2const.tw

        for i in 0..<(nSamp / 2 - tw) { pn[i] = 0 }
        var win: Float = 0
        let step: Float = 1.0 / Float(2 * tw)
        for i in (nSamp / 2 - tw)..<(nSamp / 2 + tw) {
            pn[i] = win
            win += step
        }
        for i in (nSamp / 2 + tw)..<(3 * nSamp / 2 - tw) { pn[i] = 1.0 }
        win = 1.0
        for i in (3 * nSamp / 2 - tw)..<(3 * nSamp / 2 + tw) {
            pn[i] = win
            win -= step
        }
        for i in (3 * nSamp / 2 + tw)..<(2 * nSamp) { pn[i] = 0 }
    }

    /// Sinusoidal frequency-domain synthesis followed by trapezoidal overlap-add.
    /// `snOut` length must be 2*nSamp. `shift` indicates whether to roll the
    /// synthesis memory by nSamp (the standard transition between 10 ms
    /// frames).
    public static func synthesise(nSamp: Int,
                                  fftrInvCfg: Codec2FFTRConfig,
                                  snOut: inout [Float],
                                  model: Codec2Model,
                                  pn: [Float],
                                  shift: Int) {
        let nFft = fftDec
        let twoPi: Float = Float(Codec2Constants.twoPi)
        var swSpec = [COMP](repeating: COMP(), count: nFft / 2 + 1)
        var swTime = [Float](repeating: 0, count: nFft)

        if shift != 0 {
            for i in 0..<(nSamp - 1) {
                snOut[i] = snOut[i + nSamp]
            }
            snOut[nSamp - 1] = 0
        }

        // Place each harmonic at its DFT bin (clamped to the last positive bin).
        // The bin computation matches the C reference's mixed float/double:
        // `l * model.Wo * FFT_DEC` is float, divided by TWO_PI which is double,
        // so the (cast)int truncation happens on a double.
        let twoPiD = Codec2Constants.twoPi
        _ = twoPi
        for l in 1...model.l {
            let prod = Float(l) * model.wo * Float(nFft)
            var b = Int(Double(prod) / twoPiD + 0.5)
            if b > (nFft / 2) - 1 { b = nFft / 2 - 1 }
            swSpec[b] = COMP(model.a[l] * cosf(model.phi[l]),
                             model.a[l] * sinf(model.phi[l]))
        }

        Codec2FFT.fftri(fftrInvCfg, freq: swSpec, time: &swTime)

        // Overlap-add into the trapezoidal synthesis window.
        for i in 0..<(nSamp - 1) {
            snOut[i] += swTime[nFft - nSamp + 1 + i] * pn[i]
        }
        if shift != 0 {
            var j = 0
            for i in (nSamp - 1)..<(2 * nSamp) {
                snOut[i] = swTime[j] * pn[i]
                j += 1
            }
        } else {
            var j = 0
            for i in (nSamp - 1)..<(2 * nSamp) {
                snOut[i] += swTime[j] * pn[i]
                j += 1
            }
        }
    }

    // MARK: - codec2_rand

    public static let codec2RandMax: Int = 32767
    private static var randNext: UInt64 = 1
    /// Linear congruential PRNG, matches codec2_rand() exactly (same state and
    /// constants), used for unvoiced phase synthesis. The C reference holds
    /// this state in a file-scope static; in Swift it has to be process-global
    /// too, otherwise sample-by-sample parity with c2dec is lost. Tests that
    /// need to compare against the C reference should call
    /// `Sine.resetRandState()` before instantiating Codec2 (or simply pay the
    /// cost via `Codec2.init`, which does it for them).
    public static func codec2Rand() -> Int32 {
        randNext = randNext &* 1103515245 &+ 12345
        return Int32((randNext / 65536) % 32768)
    }

    /// Resets the PRNG to its at-process-start value (next = 1). Used by
    /// Codec2.init() so each fresh codec instance behaves like a fresh C
    /// process — necessary for sample-level parity.
    public static func resetRandState() {
        randNext = 1
    }
}
