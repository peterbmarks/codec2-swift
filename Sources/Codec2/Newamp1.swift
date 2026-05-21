import Foundation

// Swift port of codec2/src/newamp1.c.
//
// "newamp1" is the rate-K mel-spaced amplitude quantiser used by mode 700C.
// Encode side: variable rate L harmonic magnitudes -> rate K=20 mel-spaced
// dB vector, mean removed, optionally equalised, then 2-stage 9+9-bit VQ.
// Energy (mean) is scalar-quantised to 4 bits and Wo to 6 bits (with the
// smallest Wo index reserved as the UV flag).
// Decode side reconstructs the rate-K vector, linearly interpolates between
// 25 Hz (40 ms) frames back to 100 Hz (10 ms), and synthesises phase via
// the minimum-phase mag_to_phase cepstral trick.

public enum Newamp1 {

    public static let nIndexes: Int = 4         // vq1, vq2, energy, Wo
    public static let phaseNFFT: Int = 128
    public static let K: Int = 20
    public static let vqMBestDepth: Int = 5

    /// Pi as a Float, only used where the C reference also operates entirely
    /// in float. Inside newamp1 the C reference picks up `M_PI` from math.h
    /// (the standard double-precision value), so for those expressions we
    /// use `mPiD` below.
    @inlinable public static var mPi: Float { Float(Double.pi) }
    /// math.h M_PI as a Double. Matches the C macro newamp1.c picks up via
    /// math.h (NOT codec2's truncated PI = 3.141592654).
    public static let mPiD: Double = Double.pi

    // MARK: - Geometry helpers

    /// 2nd-order parabolic interpolation of (xp, yp) at points x. Used in
    /// place of cubic splines for the rate-L <-> rate-K resampling.
    public static func interpPara(y: inout [Float], yOffset: Int,
                                  xp: [Float], xpOffset: Int,
                                  yp: [Float], ypOffset: Int,
                                  np: Int, x: [Float], xOffset: Int, n: Int) {
        precondition(np >= 3)
        var k = 0
        for i in 0..<n {
            let xi = x[xOffset + i]
            while xp[xpOffset + k + 1] < xi && k < (np - 3) {
                k += 1
            }
            let x1 = xp[xpOffset + k]
            let y1 = yp[ypOffset + k]
            let x2 = xp[xpOffset + k + 1]
            let y2 = yp[ypOffset + k + 1]
            let x3 = xp[xpOffset + k + 2]
            let y3 = yp[ypOffset + k + 2]
            let a = ((y3 - y2) / (x3 - x2) - (y2 - y1) / (x2 - x1)) / (x3 - x1)
            let b = ((y3 - y2) / (x3 - x2) * (x2 - x1)
                   + (y2 - y1) / (x2 - x1) * (x3 - x2)) / (x3 - x1)
            y[yOffset + i] = a * (xi - x2) * (xi - x2) + b * (xi - x2) + y2
        }
    }

    /// Hz -> mel. Quantised to integer mels via the floor + 0.5 trick.
    public static func fToMel(_ fHz: Float) -> Float {
        return floorf(2595.0 * log10f(1.0 + fHz / 700.0) + 0.5)
    }

    /// Fill K sample frequencies (kHz) linearly spaced in mels between
    /// `melStart` and `melEnd`.
    public static func melSampleFreqsKHz(_ rateKSampleFreqsKHz: inout [Float],
                                         K: Int,
                                         melStart: Float, melEnd: Float) {
        let step = (melEnd - melStart) / Float(K - 1)
        var mel = melStart
        for k in 0..<K {
            rateKSampleFreqsKHz[k] = 0.7 * (pow10f(mel / 2595.0) - 1.0)
            mel += step
        }
    }

    // MARK: - Rate resampling

    /// Encoder side: model.A[1..L] in linear -> dB, resampled at K mel-spaced
    /// frequencies. Clips dynamic range to peak-50dB to keep VQ well-conditioned.
    public static func resampleConstRateF(c2const: C2Const,
                                          model: Codec2Model,
                                          rateKVec: inout [Float],
                                          rateKSampleFreqsKHz: [Float],
                                          K: Int) {
        var AmdB = [Float](repeating: 0, count: Codec2Constants.maxAmp + 1)
        var rateLSampleFreqsKHz = [Float](repeating: 0, count: Codec2Constants.maxAmp + 1)
        var AmdBPeak: Float = -100.0
        for m in 1...model.l {
            AmdB[m] = 20.0 * log10f(model.a[m] + 1e-16)
            if AmdB[m] > AmdBPeak { AmdBPeak = AmdB[m] }
            // C does this in mixed precision: `Fs/2000.0` is double, and
            // `/ M_PI` divides by math.h's double pi. The final assignment
            // truncates to float.
            let prod = Float(m) * model.wo            // float * float
            let scale = Double(c2const.fs) / 2000.0
            rateLSampleFreqsKHz[m] = Float(Double(prod) * scale / mPiD)
        }
        for m in 1...model.l {
            if AmdB[m] < AmdBPeak - 50.0 { AmdB[m] = AmdBPeak - 50.0 }
        }
        interpPara(y: &rateKVec, yOffset: 0,
                   xp: rateLSampleFreqsKHz, xpOffset: 1,
                   yp: AmdB, ypOffset: 1,
                   np: model.l,
                   x: rateKSampleFreqsKHz, xOffset: 0, n: K)
    }

    /// Decoder side: rate-K dB vector + 0 dB endpoints -> rate-L amplitudes in
    /// `model.a[1..L]`. Mirrors resample_rate_L() exactly.
    public static func resampleRateL(c2const: C2Const,
                                     model: inout Codec2Model,
                                     rateKVec: [Float],
                                     rateKSampleFreqsKHz: [Float], K: Int) {
        var rateKVecTerm = [Float](repeating: 0, count: K + 2)
        var rateKSampleFreqsKHzTerm = [Float](repeating: 0, count: K + 2)
        rateKVecTerm[0] = 0
        rateKVecTerm[K + 1] = 0
        rateKSampleFreqsKHzTerm[0] = 0
        rateKSampleFreqsKHzTerm[K + 1] = 4.0
        for k in 0..<K {
            rateKVecTerm[k + 1] = rateKVec[k]
            rateKSampleFreqsKHzTerm[k + 1] = rateKSampleFreqsKHz[k]
        }

        var AmdB = [Float](repeating: 0, count: Codec2Constants.maxAmp + 1)
        var rateLSampleFreqsKHz = [Float](repeating: 0, count: Codec2Constants.maxAmp + 1)
        for m in 1...model.l {
            // C does this in mixed precision: `Fs/2000.0` is double, and
            // `/ M_PI` divides by math.h's double pi. The final assignment
            // truncates to float.
            let prod = Float(m) * model.wo            // float * float
            let scale = Double(c2const.fs) / 2000.0
            rateLSampleFreqsKHz[m] = Float(Double(prod) * scale / mPiD)
        }

        interpPara(y: &AmdB, yOffset: 1,
                   xp: rateKSampleFreqsKHzTerm, xpOffset: 0,
                   yp: rateKVecTerm, ypOffset: 0,
                   np: K + 2,
                   x: rateLSampleFreqsKHz, xOffset: 1, n: model.l)
        for m in 1...model.l {
            model.a[m] = pow10f(AmdB[m] / 20.0)
        }
    }

    // MARK: - 2-stage VQ encoder

    /// 2-stage rate-K VQ search using MBest. Returns MSE; writes the two
    /// 9-bit indexes into `indexes[0..1]`.
    @discardableResult
    public static func rateKMBestEncode(indexes: inout [Int],
                                        x: [Float], xq: inout [Float],
                                        ndim: Int, mbestEntries: Int) -> Float {
        let cb1 = Newamp1VQCodebooks.table[0]
        let cb2 = Newamp1VQCodebooks.table[1]
        precondition(ndim == cb1.k)

        let stage1 = MBest(entries: mbestEntries)
        let stage2 = MBest(entries: mbestEntries)
        var index = [Int32](repeating: 0, count: MBEST_STAGES)

        MBest.search(cb: cb1.cb, cbOffset: 0, vec: x,
                     k: ndim, m: cb1.m, mbest: stage1, index: &index)

        var target = [Float](repeating: 0, count: ndim)
        for j in 0..<mbestEntries {
            let n1 = Int(stage1.list[j].index[0])
            index[1] = Int32(n1)
            for i in 0..<ndim { target[i] = x[i] - cb1.cb[ndim * n1 + i] }
            MBest.search(cb: cb2.cb, cbOffset: 0, vec: target,
                         k: ndim, m: cb2.m, mbest: stage2, index: &index)
        }

        let n1 = Int(stage2.list[0].index[1])
        let n2 = Int(stage2.list[0].index[0])
        var mse: Float = 0
        for i in 0..<ndim {
            let tmp = cb1.cb[ndim * n1 + i] + cb2.cb[ndim * n2 + i]
            mse += (x[i] - tmp) * (x[i] - tmp)
            xq[i] = tmp
        }
        indexes[0] = n1
        indexes[1] = n2
        return mse
    }

    // MARK: - Postfilter / equaliser / interpolation

    /// 20 dB/dec pre-emphasis postfilter that raises formants relative to
    /// anti-formants. Energy is normalised before and after the lift.
    public static func postFilterNewamp1(vec: inout [Float],
                                         sampleFreqKHz: [Float], K: Int,
                                         pfGain: Float) {
        var pre = [Float](repeating: 0, count: K)
        var eBefore: Float = 0
        var eAfter: Float = 0
        for k in 0..<K {
            pre[k] = 20.0 * log10f(sampleFreqKHz[k] / 0.3)
            vec[k] += pre[k]
            eBefore += pow10f(vec[k] / 10.0)
            vec[k] *= pfGain
            eAfter += pow10f(vec[k] / 10.0)
        }
        let gain = eAfter / eBefore
        let gainDb = 10.0 * log10f(gain)
        for k in 0..<K { vec[k] -= gainDb; vec[k] -= pre[k] }
    }

    /// Per-bin equaliser updated towards a fixed "ideal" mean-removed
    /// spectrum. Off by default in c2dec but lives here for tooling parity.
    public static func newamp1Eq(rateKVecNoMean: inout [Float],
                                 eq: inout [Float], K: Int, eqEn: Int) {
        let ideal: [Float] = [8, 10, 12, 14, 14, 14, 14, 14, 14, 14,
                              14, 14, 14, 14, 14, 14, 14, 14, 14, -20]
        let gain: Float = 0.02
        for k in 0..<K {
            let update = rateKVecNoMean[k] - ideal[k]
            eq[k] = (1.0 - gain) * eq[k] + gain * update
            if eq[k] < 0.0 { eq[k] = 0.0 }
            if eqEn != 0 { rateKVecNoMean[k] -= eq[k] }
        }
    }

    /// Decoder Wo + voicing interpolation between 25 Hz samples back to the
    /// internal 100 Hz frame rate (M=4 sub-frames per 40 ms).
    public static func interpWoV(Wo: inout [Float], L: inout [Int],
                                 voicing: inout [Int],
                                 Wo1: Float, Wo2: Float,
                                 voicing1: Int, voicing2: Int) {
        let M = 4
        for i in 0..<M { voicing[i] = 0 }
        // C uses M_PI (double); compute in double and cast to float.
        let uvWo: Float = Float(2.0 * mPiD / 100.0)

        if voicing1 == 0 && voicing2 == 0 {
            for i in 0..<M { Wo[i] = uvWo }
        }
        if voicing1 != 0 && voicing2 == 0 {
            Wo[0] = Wo1; Wo[1] = Wo1
            Wo[2] = uvWo; Wo[3] = uvWo
            voicing[0] = 1; voicing[1] = 1
        }
        if voicing1 == 0 && voicing2 != 0 {
            Wo[0] = uvWo; Wo[1] = uvWo
            Wo[2] = Wo2; Wo[3] = Wo2
            voicing[2] = 1; voicing[3] = 1
        }
        if voicing1 != 0 && voicing2 != 0 {
            var c: Float = 1.0
            for i in 0..<M {
                Wo[i] = Wo1 * c + Wo2 * (1.0 - c)
                voicing[i] = 1
                c -= 1.0 / Float(M)
            }
        }
        for i in 0..<M {
            // C: L_[i] = floorf(M_PI / Wo_[i]). M_PI is double, Wo_ is float;
            // division promotes to double, then floorf truncates to float.
            L[i] = Int(floorf(Float(mPiD / Double(Wo[i]))))
        }
    }

    /// Linear interpolation of two rate-K vectors over four 10ms sub-frames.
    public static func newamp1Interpolate(interpolatedSurface: inout [Float],
                                          leftVec: [Float], rightVec: [Float],
                                          K: Int) {
        let M = 4
        var c: Float = 1.0
        for i in 0..<M {
            for k in 0..<K {
                interpolatedSurface[i * K + k] = leftVec[k] * c + rightVec[k] * (1.0 - c)
            }
            c -= 1.0 / Float(M)
        }
    }

    // MARK: - Phase synthesis

    /// Determine the phase of each harmonic from the rate-L magnitude via the
    /// minimum-phase cepstral trick (`mag_to_phase`). Output goes into the
    /// excitation filter `H[1..L]`.
    public static func determinePhase(c2const: C2Const,
                                      H: inout [COMP], hOffset: Int,
                                      model: Codec2Model,
                                      Nfft: Int,
                                      fwdCfg: Codec2FFTConfig,
                                      invCfg: Codec2FFTConfig) {
        let ns = Nfft / 2 + 1
        var gdbfk = [Float](repeating: 0, count: ns)
        var sampleFreqsKHz = [Float](repeating: 0, count: ns)
        var phase = [Float](repeating: 0, count: ns)
        var AmdB = [Float](repeating: 0, count: Codec2Constants.maxAmp + 1)
        var rateLSampleFreqsKHz = [Float](repeating: 0, count: Codec2Constants.maxAmp + 1)

        for m in 1...model.l {
            precondition(model.a[m] != 0)
            AmdB[m] = 20.0 * log10f(model.a[m])
            // C does this in mixed precision: `Fs/2000.0` is double, and
            // `/ M_PI` divides by math.h's double pi. The final assignment
            // truncates to float.
            let prod = Float(m) * model.wo            // float * float
            let scale = Double(c2const.fs) / 2000.0
            rateLSampleFreqsKHz[m] = Float(Double(prod) * scale / mPiD)
        }
        for i in 0..<ns {
            sampleFreqsKHz[i] = (Float(c2const.fs) / 1000.0) * Float(i) / Float(Nfft)
        }

        interpPara(y: &gdbfk, yOffset: 0,
                   xp: rateLSampleFreqsKHz, xpOffset: 1,
                   yp: AmdB, ypOffset: 1,
                   np: model.l,
                   x: sampleFreqsKHz, xOffset: 0, n: ns)

        Phase.magToPhase(phase: &phase, gdbfk: gdbfk,
                         Nfft: Nfft, fftFwdCfg: fwdCfg, fftInvCfg: invCfg)

        for m in 1...model.l {
            // C does the division by (2.0 * M_PI) in double precision.
            let num = Float(m) * model.wo * Float(Nfft)
            let b = Int(floorf(0.5 + Float(Double(num) / (2.0 * mPiD))))
            H[hOffset + m] = COMP(cosf(phase[b]), sinf(phase[b]))
        }
    }

    // MARK: - Top-level encode/decode

    /// 700C encoder front end: model -> 4 indexes (vq1, vq2, energy, Wo).
    public static func modelToIndexes(c2const: C2Const,
                                      indexes: inout [Int],
                                      model: Codec2Model,
                                      rateKVec: inout [Float],
                                      rateKSampleFreqsKHz: [Float], K: Int,
                                      mean: inout Float,
                                      rateKVecNoMean: inout [Float],
                                      rateKVecNoMeanQ: inout [Float],
                                      se: inout Float,
                                      eq: inout [Float], eqEn: Int) {
        resampleConstRateF(c2const: c2const, model: model,
                           rateKVec: &rateKVec,
                           rateKSampleFreqsKHz: rateKSampleFreqsKHz, K: K)

        var sum: Float = 0
        for k in 0..<K { sum += rateKVec[k] }
        mean = sum / Float(K)
        for k in 0..<K { rateKVecNoMean[k] = rateKVec[k] - mean }

        newamp1Eq(rateKVecNoMean: &rateKVecNoMean, eq: &eq, K: K, eqEn: eqEn)
        rateKMBestEncode(indexes: &indexes,
                         x: rateKVecNoMean, xq: &rateKVecNoMeanQ,
                         ndim: K, mbestEntries: vqMBestDepth)

        for k in 0..<K {
            let d = rateKVecNoMean[k] - rateKVecNoMeanQ[k]
            se += d * d
        }

        let w: [Float] = [1.0]
        var seMean: Float = 0
        let energyCb = Newamp1EnergyCodebook.table[0]
        let meanArr = [mean]
        indexes[2] = Quantise.quantise(cb: energyCb.cb, vec: meanArr,
                                       cbOffsetIntoVec: 0,
                                       w: w, k: energyCb.k, m: energyCb.m, se: &seMean)

        if model.voiced != 0 {
            var index = Quantise.encodeLogWo(c2const: c2const, wo: model.wo, bits: 6)
            if index == 0 { index = 1 }
            indexes[3] = index
        } else {
            indexes[3] = 0
        }
    }

    /// Reconstruct the rate-K vector (mean restored, optionally postfiltered).
    public static func indexesToRateKVec(rateKVec: inout [Float],
                                         rateKVecNoMean: inout [Float],
                                         rateKSampleFreqsKHz: [Float], K: Int,
                                         meanOut: inout Float,
                                         indexes: [Int],
                                         userRateKVecNoMean: [Float]?,
                                         postFilterEn: Int) {
        let cb1 = Newamp1VQCodebooks.table[0]
        let cb2 = Newamp1VQCodebooks.table[1]
        let n1 = indexes[0]
        let n2 = indexes[1]

        if userRateKVecNoMean == nil {
            for k in 0..<K {
                rateKVecNoMean[k] = cb1.cb[K * n1 + k] + cb2.cb[K * n2 + k]
            }
        } else if let inj = userRateKVecNoMean {
            for k in 0..<K { rateKVecNoMean[k] = inj[k] }
        }

        if postFilterEn != 0 {
            postFilterNewamp1(vec: &rateKVecNoMean,
                              sampleFreqKHz: rateKSampleFreqsKHz, K: K,
                              pfGain: 1.5)
        }
        meanOut = Newamp1EnergyCodebook.table[0].cb[indexes[2]]
        for k in 0..<K { rateKVec[k] = rateKVecNoMean[k] + meanOut }
    }

    /// 700C decoder top level. Fills `models[0..3]` and `H[0..3]` (per-frame
    /// excitation filter spectrum) ready for `synthesise_one_frame`.
    public static func indexesToModel(c2const: C2Const,
                                      models: inout [Codec2Model],
                                      H: inout [COMP],     // length 4 * (MAX_AMP+1)
                                      interpolatedSurface: inout [Float], // 4*K
                                      prevRateKVec: inout [Float],
                                      WoLeft: inout Float,
                                      voicingLeft: inout Int,
                                      rateKSampleFreqsKHz: [Float], K: Int,
                                      fwdCfg: Codec2FFTConfig,
                                      invCfg: Codec2FFTConfig,
                                      indexes: [Int],
                                      userRateKVecNoMean: [Float]?,
                                      postFilterEn: Int) {
        let M = 4
        var rateKVec = [Float](repeating: 0, count: K)
        var rateKVecNoMean = [Float](repeating: 0, count: K)
        var mean: Float = 0

        indexesToRateKVec(rateKVec: &rateKVec,
                          rateKVecNoMean: &rateKVecNoMean,
                          rateKSampleFreqsKHz: rateKSampleFreqsKHz, K: K,
                          meanOut: &mean,
                          indexes: indexes,
                          userRateKVecNoMean: userRateKVecNoMean,
                          postFilterEn: postFilterEn)

        let WoRight: Float
        let voicingRight: Int
        if indexes[3] != 0 {
            WoRight = Quantise.decodeLogWo(c2const: c2const, index: indexes[3], bits: 6)
            voicingRight = 1
        } else {
            WoRight = Float(2.0 * mPiD / 100.0)
            voicingRight = 0
        }

        newamp1Interpolate(interpolatedSurface: &interpolatedSurface,
                           leftVec: prevRateKVec, rightVec: rateKVec, K: K)

        var aWo = [Float](repeating: 0, count: M)
        var aL = [Int](repeating: 0, count: M)
        var aVoiced = [Int](repeating: 0, count: M)
        interpWoV(Wo: &aWo, L: &aL, voicing: &aVoiced,
                  Wo1: WoLeft, Wo2: WoRight,
                  voicing1: voicingLeft, voicing2: voicingRight)

        for i in 0..<M {
            models[i].wo = aWo[i]
            models[i].l = aL[i]
            models[i].voiced = Int32(aVoiced[i])
            resampleRateL(c2const: c2const, model: &models[i],
                          rateKVec: Array(interpolatedSurface[(K * i)..<(K * i + K)]),
                          rateKSampleFreqsKHz: rateKSampleFreqsKHz, K: K)
            determinePhase(c2const: c2const, H: &H, hOffset: (Codec2Constants.maxAmp + 1) * i,
                           model: models[i], Nfft: phaseNFFT,
                           fwdCfg: fwdCfg, invCfg: invCfg)
        }

        for k in 0..<K { prevRateKVec[k] = rateKVec[k] }
        WoLeft = WoRight
        voicingLeft = voicingRight
    }
}
