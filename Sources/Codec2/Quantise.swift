import Foundation

// Swift port of codec2/src/quantise.c.
//
// This file currently covers the subset needed for Codec 2 mode 3200:
//   - speech_to_uq_lsps        (LPC analysis + LPC->LSP)
//   - encode_Wo / decode_Wo     (scalar pitch quantiser)
//   - encode_energy / decode_energy (scalar energy quantiser)
//   - encode_lspds_scalar / decode_lspds_scalar (LSP difference quantiser)
//   - lspd_bits
//   - aks_to_M2                (LPCs -> harmonic magnitudes via FFT)
//   - lpc_post_filter
//   - apply_lpc_correction
//   - quantise                 (inner Euclidean codebook search)
//
// The remaining encode/decode helpers (encode_lsps_scalar/vq, quantise_WoE,
// etc.) will follow when the larger modes (2400/1600/1400/1300/1200) get
// wired up.

public enum Quantise {

    public static let woBits: Int = 7
    public static let woLevels: Int = 1 << 7
    public static let eBits: Int = 5
    public static let eLevels: Int = 1 << 5
    public static let eMinDB: Float = -10.0
    public static let eMaxDB: Float =  40.0
    public static let lspScalarIndexes: Int = 10
    public static let lspdScalarIndexes: Int = 10
    public static let woEBits: Int = 8
    public static let lpcpfGamma: Float = 0.5
    public static let lpcpfBeta: Float = 0.2
    public static let lspDelta1: Float = 0.01

    /// Number of bits in the LSP-difference scalar codebook at index i.
    public static func lspdBits(_ i: Int) -> Int { LspdCodebooks.table[i].log2m }

    /// Pick the codebook entry minimising weighted squared error and add the
    /// best error to `se`. Returns the chosen index.
    @discardableResult
    public static func quantise(cb: [Float], vec: [Float], cbOffsetIntoVec: Int = 0,
                                w: [Float], k: Int, m: Int, se: inout Float) -> Int {
        var bestI = 0
        var bestE: Float = 1e32
        for j in 0..<m {
            var e: Float = 0
            for i in 0..<k {
                let diff = cb[j * k + i] - vec[cbOffsetIntoVec + i]
                e += diff * w[i] * diff * w[i]
            }
            if e < bestE { bestE = e; bestI = j }
        }
        se += bestE
        return bestI
    }

    // MARK: - Wo quantisation

    public static func encodeWo(c2const: C2Const, wo: Float, bits: Int) -> Int {
        let levels = 1 << bits
        let woMin = c2const.woMin
        let woMax = c2const.woMax
        let norm = (wo - woMin) / (woMax - woMin)
        var index = Int(floorf(Float(levels) * norm + 0.5))
        if index < 0 { index = 0 }
        if index > levels - 1 { index = levels - 1 }
        return index
    }

    public static func decodeWo(c2const: C2Const, index: Int, bits: Int) -> Float {
        let levels = 1 << bits
        let woMin = c2const.woMin
        let woMax = c2const.woMax
        let step = (woMax - woMin) / Float(levels)
        return woMin + step * Float(index)
    }

    /// Logarithmic Wo quantiser used by newamp1 (700C). Indexes are
    /// uniformly spaced in log10(Wo).
    public static func encodeLogWo(c2const: C2Const, wo: Float, bits: Int) -> Int {
        let levels = 1 << bits
        let woMin = c2const.woMin
        let woMax = c2const.woMax
        let norm = (log10f(wo) - log10f(woMin)) / (log10f(woMax) - log10f(woMin))
        var index = Int(floorf(Float(levels) * norm + 0.5))
        if index < 0 { index = 0 }
        if index > levels - 1 { index = levels - 1 }
        return index
    }

    public static func decodeLogWo(c2const: C2Const, index: Int, bits: Int) -> Float {
        let levels = 1 << bits
        let woMin = c2const.woMin
        let woMax = c2const.woMax
        let step = (log10f(woMax) - log10f(woMin)) / Float(levels)
        let wo = log10f(woMin) + step * Float(index)
        return pow10f(wo)
    }

    // MARK: - Energy quantisation

    public static func encodeEnergy(_ eIn: Float, bits: Int) -> Int {
        let levels = 1 << bits
        let e = 10.0 * log10f(eIn)
        let norm = (e - eMinDB) / (eMaxDB - eMinDB)
        var index = Int(floorf(Float(levels) * norm + 0.5))
        if index < 0 { index = 0 }
        if index > levels - 1 { index = levels - 1 }
        return index
    }

    public static func decodeEnergy(index: Int, bits: Int) -> Float {
        let levels = 1 << bits
        let step = (eMaxDB - eMinDB) / Float(levels)
        let eDb = eMinDB + step * Float(index)
        return pow10f(eDb / 10.0)
    }

    // MARK: - LSP-difference scalar quantiser (3200/1600/1400)

    public static func encodeLspdsScalar(indexes: inout [Int], lsp: [Float], order: Int) {
        let pi: Float = Float(Codec2Constants.pi)
        var lspHz = [Float](repeating: 0, count: order)
        var lspHzQ = [Float](repeating: 0, count: order)
        var dlsp = [Float](repeating: 0, count: order)
        var dlspQ = [Float](repeating: 0, count: order)
        let wt: [Float] = [1.0]

        for i in 0..<order { lspHz[i] = (4000.0 / pi) * lsp[i] }

        var se: Float = 0
        for i in 0..<order {
            if i > 0 {
                dlsp[i] = lspHz[i] - lspHzQ[i - 1]
            } else {
                dlsp[0] = lspHz[0]
            }
            let book = LspdCodebooks.table[i]
            let k = book.k
            let m = book.m
            indexes[i] = quantise(cb: book.cb, vec: dlsp, cbOffsetIntoVec: i,
                                  w: wt, k: k, m: m, se: &se)
            dlspQ[i] = book.cb[indexes[i] * k]
            if i > 0 {
                lspHzQ[i] = lspHzQ[i - 1] + dlspQ[i]
            } else {
                lspHzQ[0] = dlspQ[0]
            }
        }
    }

    public static func decodeLspdsScalar(lspOut: inout [Float], indexes: [Int], order: Int) {
        let pi: Float = Float(Codec2Constants.pi)
        var lspHzQ = [Float](repeating: 0, count: order)
        for i in 0..<order {
            let book = LspdCodebooks.table[i]
            let k = book.k
            let dlspQ = book.cb[indexes[i] * k]
            if i > 0 {
                lspHzQ[i] = lspHzQ[i - 1] + dlspQ
            } else {
                lspHzQ[0] = dlspQ
            }
            lspOut[i] = (pi / 4000.0) * lspHzQ[i]
        }
    }

    // MARK: - LPC -> harmonic magnitudes

    /// Compute the LPC synthesis-filter power spectrum, optionally apply a
    /// frequency-domain post filter, then integrate by harmonic to fill
    /// `model.a[]`. `Aw` is filled with the LPC spectrum (used later by phase
    /// synthesis).
    public static func aksToM2(fftrFwdCfg: Codec2FFTRConfig,
                               ak: [Float], order: Int,
                               model: inout Codec2Model,
                               E: Float,
                               snr: inout Float,
                               simPF: Int, pf: Int,
                               bassBoost: Int,
                               beta: Float, gamma: Float,
                               Aw: inout [COMP]) {
        let nFft = Sine.fftEnc
        let r = Float(Codec2Constants.twoPi) / Float(nFft)

        var aBuf = [Float](repeating: 0, count: nFft)
        for i in 0...order { aBuf[i] = ak[i] }
        Codec2FFT.fftr(fftrFwdCfg, time: aBuf, freq: &Aw)

        // |A(jw)|^2 -> reciprocal power spectrum P(w) = 1 / (|A|^2 + ε).
        var pw = [Float](repeating: 0, count: nFft / 2)
        for i in 0..<(nFft / 2) {
            pw[i] = 1.0 / (Aw[i].real * Aw[i].real + Aw[i].imag * Aw[i].imag + 1e-6)
        }

        if pf != 0 {
            lpcPostFilter(fftrFwdCfg: fftrFwdCfg, pw: &pw, ak: ak, order: order,
                          beta: beta, gamma: gamma, bassBoost: bassBoost, E: E)
        } else {
            for i in 0..<(nFft / 2) { pw[i] *= E }
        }

        var signal: Float = 1e-30
        var noise: Float = 1e-32
        for m in 1...model.l {
            let am = Int((Float(m) - 0.5) * model.wo / r + 0.5)
            var bm = Int((Float(m) + 0.5) * model.wo / r + 0.5)
            if bm > nFft / 2 { bm = nFft / 2 }
            var em: Float = 0
            for i in am..<bm { em += pw[i] }
            var Am = sqrtf(em)

            signal += model.a[m] * model.a[m]
            noise += (model.a[m] - Am) * (model.a[m] - Am)

            if simPF != 0 {
                if Am > model.a[m] { Am *= 0.7 }
                if Am < model.a[m] { Am *= 1.4 }
            }
            model.a[m] = Am
        }
        snr = 10.0 * log10f(signal / noise)
    }

    /// Frequency-domain LPC post filter as in `lpc_post_filter` in quantise.c.
    /// Suppresses inter-formant energy by multiplying `Pw` by R^β where R is
    /// the sqrt of (|W|² · Pw) and W is the γ-weighted denominator polynomial.
    public static func lpcPostFilter(fftrFwdCfg: Codec2FFTRConfig,
                                     pw: inout [Float], ak: [Float], order: Int,
                                     beta: Float, gamma: Float, bassBoost: Int, E: Float) {
        let nFft = Sine.fftEnc
        var x = [Float](repeating: 0, count: nFft)
        var ww = [COMP](repeating: COMP(), count: nFft / 2 + 1)

        x[0] = ak[0]
        var coeff = gamma
        for i in 1...order {
            x[i] = ak[i] * coeff
            coeff *= gamma
        }
        Codec2FFT.fftr(fftrFwdCfg, time: x, freq: &ww)

        // |W|^2 stored back in .real
        for i in 0..<(nFft / 2) {
            ww[i] = COMP(ww[i].real * ww[i].real + ww[i].imag * ww[i].imag, 0)
        }

        // R = sqrt(|W|^2 * Pw)
        var rw = [Float](repeating: 0, count: nFft / 2)
        for i in 0..<(nFft / 2) { rw[i] = sqrtf(ww[i].real * pw[i]) }

        var eBefore: Float = 1e-4
        for i in 0..<(nFft / 2) { eBefore += pw[i] }

        var eAfter: Float = 1e-4
        for i in 0..<(nFft / 2) {
            let pfw = powf(rw[i], beta)
            pw[i] *= pfw * pfw
            eAfter += pw[i]
        }
        var gain = eBefore / eAfter
        gain *= E
        for i in 0..<(nFft / 2) { pw[i] *= gain }

        if bassBoost != 0 {
            for i in 0..<(nFft / 8) {
                pw[i] *= 1.4 * 1.4
            }
        }
    }

    // MARK: - LPC analysis to LSP

    /// Windowed LPC analysis followed by LPC->LSP conversion. Returns the
    /// residual energy `E`. If the input is silent or LPC root finding fails,
    /// falls back to an evenly-spaced LSP vector to keep downstream code
    /// sane.
    public static func speechToUqLsps(lspOut: inout [Float], ak: inout [Float],
                                       sn: [Float], w: [Float],
                                       mPitch: Int, order: Int) -> Float {
        let pi: Float = Float(Codec2Constants.pi)
        var wn = [Float](repeating: 0, count: mPitch)
        var r = [Float](repeating: 0, count: order + 1)
        var e: Float = 0

        for i in 0..<mPitch {
            wn[i] = sn[i] * w[i]
            e += wn[i] * wn[i]
        }
        if e == 0 {
            for i in 0..<order { lspOut[i] = (pi / Float(order)) * Float(i) }
            return 0
        }

        LPC.autocorrelate(wn, rn: &r, nSam: mPitch, order: order)
        LPC.levinsonDurbin(r, lpcs: &ak, order: order)

        var energy: Float = 0
        for i in 0...order { energy += ak[i] * r[i] }

        // 15 Hz BW expansion to harden the LSP root finder.
        for i in 0...order { ak[i] *= powf(0.994, Float(i)) }

        let roots = LSP.lpcToLsp(ak, order: order, freq: &lspOut, nb: 5, delta: lspDelta1)
        if roots != order {
            for i in 0..<order { lspOut[i] = (pi / Float(order)) * Float(i) }
        }
        return energy
    }

    // MARK: - misc

    /// First-harmonic LPC correction for low-pitch males. Mirrors apply_lpc_correction.
    public static func applyLpcCorrection(_ model: inout Codec2Model) {
        if model.wo < Float(Codec2Constants.pi) * 150.0 / 4000.0 {
            model.a[1] *= 0.032
        }
    }
}
