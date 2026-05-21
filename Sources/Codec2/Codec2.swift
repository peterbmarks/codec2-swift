import Foundation

// Swift port of codec2/src/codec2.c.
//
// Currently implements the public API plus the codepath for mode 3200.
// Other modes (2400, 1600, 1400, 1300, 1200, 700C) share the same state
// machine; their per-mode encoders/decoders will land alongside the
// codebooks they need.

public enum Codec2Mode: Int {
    case mode3200 = 0
    case mode2400 = 1
    case mode1600 = 2
    case mode1400 = 3
    case mode1300 = 4
    case mode1200 = 5
    case mode700C = 8

    public static func fromString(_ s: String) -> Codec2Mode? {
        switch s.uppercased() {
        case "3200": return .mode3200
        case "2400": return .mode2400
        case "1600": return .mode1600
        case "1400": return .mode1400
        case "1300": return .mode1300
        case "1200": return .mode1200
        case "700C": return .mode700C
        default: return nil
        }
    }
}

/// Full state of a single codec2 encoder/decoder. One instance is sufficient
/// for full duplex use, matching the C `struct CODEC2`.
public final class Codec2 {

    public let mode: Codec2Mode
    public var c2const: C2Const
    public var fs: Int { c2const.fs }
    public var nSamp: Int { c2const.nSamp }
    public var mPitch: Int { c2const.mPitch }

    // Analysis buffers and synthesis memory.
    public var sn: [Float]                   // m_pitch samples of speech history
    public var w: [Float]                    // m_pitch-sample analysis window
    public var W: [Float]                    // FFT_ENC samples, freq-domain window
    public var snSynth: [Float]              // 2*n_samp synthesis memory (Sn_)
    public var pn: [Float]                   // 2*n_samp Parzen synthesis window

    public var fftFwdCfg: Codec2FFTConfig
    public var fftrFwdCfg: Codec2FFTRConfig
    public var fftrInvCfg: Codec2FFTRConfig
    public var nlp: NLPState

    public var prevF0Enc: Float
    public var bgEst: Float
    public var exPhase: Float

    public var prevModelDec: Codec2Model
    public var prevLspsDec: [Float]
    public var prevEDec: Float

    public var lpcPF: Int = 1
    public var bassBoost: Int = 1
    public var beta: Float = Quantise.lpcpfBeta
    public var gamma: Float = Quantise.lpcpfGamma

    public var hpfStates: [Float] = [0, 0]
    public var gray: Int = 1

    // -- newamp1 (mode 700C) state --
    public var rateKSampleFreqsKHz: [Float] = []
    public var prevRateKVec: [Float] = []
    public var eq: [Float] = []
    public var eqEn: Int = 0
    public var WoLeft: Float = 0
    public var voicingLeft: Int = 0
    public var phaseFftFwdCfg: Codec2FFTConfig?
    public var phaseFftInvCfg: Codec2FFTConfig?
    public var se: Float = 0
    public var nse: Int = 0
    public var postFilterEn: Int = 1

    public init?(mode: Codec2Mode) {
        self.mode = mode
        // Match C process-start state for the unvoiced-phase PRNG. The C
        // reference uses a file-scope static initialised to 1 — equivalent to
        // resetting once per c2enc/c2dec process launch.
        Sine.resetRandState()
        let c = Sine.c2constCreate(fs: 8000, framelengthS: Float(Codec2Constants.nSeconds))
        self.c2const = c

        self.sn = [Float](repeating: 1.0, count: c.mPitch)
        self.w = [Float](repeating: 0, count: c.mPitch)
        self.W = [Float](repeating: 0, count: Codec2Constants.fftEnc)
        self.snSynth = [Float](repeating: 0, count: 2 * c.nSamp)
        self.pn = [Float](repeating: 0, count: 2 * c.nSamp)

        self.fftFwdCfg = Codec2FFTConfig(nfft: Codec2Constants.fftEnc, inverse: false)
        self.fftrFwdCfg = Codec2FFTRConfig(nfft: Codec2Constants.fftEnc, inverse: false)
        self.fftrInvCfg = Codec2FFTRConfig(nfft: Codec2Constants.fftDec, inverse: true)

        self.nlp = NLPState(c2const: c)

        self.prevF0Enc = 1.0 / Float(Codec2Constants.pMaxSeconds)
        self.bgEst = 0
        self.exPhase = 0

        var prev = Codec2Model()
        prev.wo = Float(Codec2Constants.twoPi) / Float(c.pMax)
        prev.l = Int(Float(Codec2Constants.pi) / prev.wo)
        prev.voiced = 0
        self.prevModelDec = prev

        var lsps = [Float](repeating: 0, count: Codec2Constants.lpcOrd)
        for i in 0..<Codec2Constants.lpcOrd {
            lsps[i] = Float(i) * Float(Codec2Constants.pi) / Float(Codec2Constants.lpcOrd + 1)
        }
        self.prevLspsDec = lsps
        self.prevEDec = 1

        Sine.makeAnalysisWindow(c2const: c, fftFwdCfg: fftFwdCfg, w: &self.w, W: &self.W)
        Sine.makeSynthesisWindow(c2const: c, pn: &self.pn)

        // newamp1 (700C) setup.
        if mode == .mode700C {
            self.rateKSampleFreqsKHz = [Float](repeating: 0, count: Newamp1.K)
            Newamp1.melSampleFreqsKHz(&self.rateKSampleFreqsKHz, K: Newamp1.K,
                                      melStart: Newamp1.fToMel(200.0),
                                      melEnd: Newamp1.fToMel(3700.0))
            self.prevRateKVec = [Float](repeating: 0, count: Newamp1.K)
            self.eq = [Float](repeating: 0, count: Newamp1.K)
            self.phaseFftFwdCfg = Codec2FFTConfig(nfft: Newamp1.phaseNFFT, inverse: false)
            self.phaseFftInvCfg = Codec2FFTConfig(nfft: Newamp1.phaseNFFT, inverse: true)
        }
    }

    /// Bits per encoded frame for this mode.
    public var bitsPerFrame: Int {
        switch mode {
        case .mode3200: return 64
        case .mode2400: return 48
        case .mode1600: return 64
        case .mode1400: return 56
        case .mode1300: return 52
        case .mode1200: return 48
        case .mode700C: return 28
        }
    }

    /// Bytes needed to hold one encoded frame.
    public var bytesPerFrame: Int { (bitsPerFrame + 7) / 8 }

    /// Speech samples (8 kHz, int16) per frame.
    public var samplesPerFrame: Int {
        switch mode {
        case .mode3200, .mode2400: return 160
        default: return 320
        }
    }

    // MARK: - Encode/decode entry points

    public func encode(speech: [Int16], bits: inout [UInt8]) {
        precondition(speech.count >= samplesPerFrame)
        precondition(bits.count >= bytesPerFrame)
        switch mode {
        case .mode3200: encode3200(speech: speech, bits: &bits)
        case .mode700C: encode700C(speech: speech, bits: &bits)
        default:
            fatalError("Codec2.encode: mode \(mode) not yet ported")
        }
    }

    public func decode(bits: [UInt8], speech: inout [Int16]) {
        precondition(speech.count >= samplesPerFrame)
        switch mode {
        case .mode3200: decode3200(bits: bits, speech: &speech)
        case .mode700C: decode700C(bits: bits, speech: &speech)
        default:
            fatalError("Codec2.decode: mode \(mode) not yet ported")
        }
    }

    // MARK: - Mode 3200

    private func encode3200(speech: [Int16], bits: inout [UInt8]) {
        for i in 0..<bytesPerFrame { bits[i] = 0 }

        var model = Codec2Model()
        var ak = [Float](repeating: 0, count: Codec2Constants.lpcOrd + 1)
        var lsps = [Float](repeating: 0, count: Codec2Constants.lpcOrd)
        var nbit: UInt32 = 0
        var lspdIndexes = [Int](repeating: 0, count: Codec2Constants.lpcOrd)

        // First 10ms frame: emit voicing bit only.
        analyseOneFrame(model: &model, speechSlice: Array(speech[0..<nSamp]))
        pack(&bits, bitIndex: &nbit, field: model.voiced, fieldWidth: 1)

        // Second 10ms frame: full parameter set.
        analyseOneFrame(model: &model, speechSlice: Array(speech[nSamp..<2 * nSamp]))
        pack(&bits, bitIndex: &nbit, field: model.voiced, fieldWidth: 1)

        let woIndex = Quantise.encodeWo(c2const: c2const, wo: model.wo, bits: Quantise.woBits)
        pack(&bits, bitIndex: &nbit, field: Int32(woIndex),
             fieldWidth: UInt32(Quantise.woBits))

        let e = Quantise.speechToUqLsps(lspOut: &lsps, ak: &ak,
                                        sn: sn, w: w,
                                        mPitch: mPitch, order: Codec2Constants.lpcOrd)
        let eIndex = Quantise.encodeEnergy(e, bits: Quantise.eBits)
        pack(&bits, bitIndex: &nbit, field: Int32(eIndex),
             fieldWidth: UInt32(Quantise.eBits))

        Quantise.encodeLspdsScalar(indexes: &lspdIndexes, lsp: lsps, order: Codec2Constants.lpcOrd)
        for i in 0..<Quantise.lspdScalarIndexes {
            pack(&bits, bitIndex: &nbit,
                 field: Int32(lspdIndexes[i]),
                 fieldWidth: UInt32(Quantise.lspdBits(i)))
        }
        precondition(nbit == UInt32(bitsPerFrame))
    }

    private func decode3200(bits: [UInt8], speech: inout [Int16]) {
        var nbit: UInt32 = 0
        var model = [Codec2Model(), Codec2Model()]
        var lsps: [[Float]] = [
            [Float](repeating: 0, count: Codec2Constants.lpcOrd),
            [Float](repeating: 0, count: Codec2Constants.lpcOrd)
        ]
        var ak: [[Float]] = [
            [Float](repeating: 0, count: Codec2Constants.lpcOrd + 1),
            [Float](repeating: 0, count: Codec2Constants.lpcOrd + 1)
        ]
        var e = [Float](repeating: 0, count: 2)
        var snr: Float = 0
        var Aw = [COMP](repeating: COMP(), count: Codec2Constants.fftEnc / 2 + 1)
        var lspdIndexes = [Int](repeating: 0, count: Codec2Constants.lpcOrd)

        model[0].voiced = unpack(bits, bitIndex: &nbit, fieldWidth: 1)
        model[1].voiced = unpack(bits, bitIndex: &nbit, fieldWidth: 1)

        let woIndex = Int(unpack(bits, bitIndex: &nbit, fieldWidth: UInt32(Quantise.woBits)))
        model[1].wo = Quantise.decodeWo(c2const: c2const, index: woIndex, bits: Quantise.woBits)
        model[1].l = Int(Float(Codec2Constants.pi) / model[1].wo)

        let eIndex = Int(unpack(bits, bitIndex: &nbit, fieldWidth: UInt32(Quantise.eBits)))
        e[1] = Quantise.decodeEnergy(index: eIndex, bits: Quantise.eBits)

        for i in 0..<Quantise.lspdScalarIndexes {
            lspdIndexes[i] = Int(unpack(bits, bitIndex: &nbit,
                                        fieldWidth: UInt32(Quantise.lspdBits(i))))
        }
        Quantise.decodeLspdsScalar(lspOut: &lsps[1], indexes: lspdIndexes,
                                    order: Codec2Constants.lpcOrd)

        Interp.interpWo(&model[0], prev: prevModelDec, next: model[1], woMin: c2const.woMin)
        e[0] = Interp.interpEnergy(prevEDec, e[1])
        Interp.interpolateLspVer2(&lsps[0], prev: prevLspsDec, next: lsps[1],
                                  weight: 0.5, order: Codec2Constants.lpcOrd)

        for i in 0..<2 {
            LSP.lspToLpc(lsps[i], ak: &ak[i], order: Codec2Constants.lpcOrd)
            Quantise.aksToM2(fftrFwdCfg: fftrFwdCfg,
                             ak: ak[i], order: Codec2Constants.lpcOrd,
                             model: &model[i], E: e[i], snr: &snr,
                             simPF: 0, pf: lpcPF,
                             bassBoost: bassBoost,
                             beta: beta, gamma: gamma,
                             Aw: &Aw)
            Quantise.applyLpcCorrection(&model[i])
            synthesiseOneFrame(speech: &speech, speechOffset: nSamp * i,
                               model: &model[i], Aw: Aw, gain: 1.0)
        }

        prevModelDec = model[1]
        prevEDec = e[1]
        for i in 0..<Codec2Constants.lpcOrd { prevLspsDec[i] = lsps[1][i] }
    }

    // MARK: - Mode 700C (newamp1)

    private func encode700C(speech: [Int16], bits: inout [UInt8]) {
        for i in 0..<bytesPerFrame { bits[i] = 0 }

        var model = Codec2Model()
        // Analyse 4 internal 10ms frames; only the last frame's model
        // parameters get quantised — the C reference does the same.
        for i in 0..<4 {
            let slice = Array(speech[(i * nSamp)..<((i + 1) * nSamp)])
            analyseOneFrame(model: &model, speechSlice: slice)
        }

        var indexes = [Int](repeating: 0, count: Newamp1.nIndexes)
        var rateKVec = [Float](repeating: 0, count: Newamp1.K)
        var rateKVecNoMean = [Float](repeating: 0, count: Newamp1.K)
        var rateKVecNoMeanQ = [Float](repeating: 0, count: Newamp1.K)
        var meanVal: Float = 0

        Newamp1.modelToIndexes(c2const: c2const,
                               indexes: &indexes,
                               model: model,
                               rateKVec: &rateKVec,
                               rateKSampleFreqsKHz: rateKSampleFreqsKHz,
                               K: Newamp1.K,
                               mean: &meanVal,
                               rateKVecNoMean: &rateKVecNoMean,
                               rateKVecNoMeanQ: &rateKVecNoMeanQ,
                               se: &se,
                               eq: &eq, eqEn: eqEn)
        nse += Newamp1.K

        var nbit: UInt32 = 0
        packNaturalOrGray(&bits, bitIndex: &nbit, field: Int32(indexes[0]),
                          fieldWidth: 9, gray: 0)
        packNaturalOrGray(&bits, bitIndex: &nbit, field: Int32(indexes[1]),
                          fieldWidth: 9, gray: 0)
        packNaturalOrGray(&bits, bitIndex: &nbit, field: Int32(indexes[2]),
                          fieldWidth: 4, gray: 0)
        packNaturalOrGray(&bits, bitIndex: &nbit, field: Int32(indexes[3]),
                          fieldWidth: 6, gray: 0)
        precondition(nbit == UInt32(bitsPerFrame))
    }

    private func decode700C(bits: [UInt8], speech: inout [Int16]) {
        var nbit: UInt32 = 0
        var indexes = [Int](repeating: 0, count: Newamp1.nIndexes)
        indexes[0] = Int(unpackNaturalOrGray(bits, bitIndex: &nbit, fieldWidth: 9, gray: 0))
        indexes[1] = Int(unpackNaturalOrGray(bits, bitIndex: &nbit, fieldWidth: 9, gray: 0))
        indexes[2] = Int(unpackNaturalOrGray(bits, bitIndex: &nbit, fieldWidth: 4, gray: 0))
        indexes[3] = Int(unpackNaturalOrGray(bits, bitIndex: &nbit, fieldWidth: 6, gray: 0))

        let M = 4
        var models = [Codec2Model](repeating: Codec2Model(), count: M)
        // HH is M × (MAX_AMP+1) excitation-filter spectra.
        var HH = [COMP](repeating: COMP(), count: M * (Codec2Constants.maxAmp + 1))
        var interpolatedSurface = [Float](repeating: 0, count: M * Newamp1.K)

        guard let fwdCfg = phaseFftFwdCfg, let invCfg = phaseFftInvCfg else {
            fatalError("700C phase FFT configs missing — Codec2 must be created with mode .mode700C")
        }

        Newamp1.indexesToModel(c2const: c2const,
                               models: &models,
                               H: &HH,
                               interpolatedSurface: &interpolatedSurface,
                               prevRateKVec: &prevRateKVec,
                               WoLeft: &WoLeft,
                               voicingLeft: &voicingLeft,
                               rateKSampleFreqsKHz: rateKSampleFreqsKHz,
                               K: Newamp1.K,
                               fwdCfg: fwdCfg, invCfg: invCfg,
                               indexes: indexes,
                               userRateKVecNoMean: nil,
                               postFilterEn: postFilterEn)

        // 700C is a little quieter; the C reference applies a 1.5x gain.
        for i in 0..<M {
            let base = i * (Codec2Constants.maxAmp + 1)
            let Hi = Array(HH[base..<(base + Codec2Constants.maxAmp + 1)])
            synthesiseOneFrame(speech: &speech, speechOffset: nSamp * i,
                               model: &models[i], Aw: Hi, gain: 1.5)
        }
    }

    // MARK: - Frame helpers

    /// Sinusoidal analysis of one 10 ms (nSamp) speech frame. Shifts the
    /// internal `sn[]` history, appends the new samples, runs the pitch +
    /// DFT + amplitude + voicing pipeline.
    public func analyseOneFrame(model: inout Codec2Model, speechSlice: [Int16]) {
        var sw = [COMP](repeating: COMP(), count: Codec2Constants.fftEnc)
        var pitch: Float = 0

        for i in 0..<(mPitch - nSamp) { sn[i] = sn[i + nSamp] }
        for i in 0..<nSamp { sn[i + mPitch - nSamp] = Float(speechSlice[i]) }

        Sine.dftSpeech(c2const: c2const, fftFwdCfg: fftFwdCfg, sw: &sw, sn: sn, w: w)
        _ = NLP.run(state: nlp, sn: sn, n: nSamp, pitch: &pitch, prevF0: &prevF0Enc)
        model.wo = Float(Codec2Constants.twoPi) / pitch
        model.l = Int(Float(Codec2Constants.pi) / model.wo)

        Sine.twoStagePitchRefinement(c2const: c2const, model: &model, sw: sw)
        Sine.estimateAmplitudes(model: &model, sw: sw, W: W, estPhase: 0)
        _ = Sine.estVoicingMBE(c2const: c2const, model: &model, sw: sw, W: W)
    }

    /// Synthesise a single 10 ms frame. Applies LPC phase synthesis,
    /// background-noise postfilter, overlap-add, gain, ear protection and
    /// clipping. Writes nSamp Int16 samples starting at `speechOffset`.
    public func synthesiseOneFrame(speech: inout [Int16], speechOffset: Int,
                                   model: inout Codec2Model, Aw: [COMP], gain: Float) {
        var H = [COMP](repeating: COMP(), count: Codec2Constants.maxAmp + 1)
        if mode == .mode700C {
            // newamp1 — rate-L phase already computed by determine_phase.
            for i in 0..<min(H.count, Aw.count) { H[i] = Aw[i] }
        } else {
            Phase.samplePhase(model: model, H: &H, A: Aw)
        }
        Phase.phaseSynthZeroOrder(nSamp: nSamp, model: &model, exPhase: &exPhase, H: H)

        Postfilter.postfilter(&model, bgEst: &bgEst)
        Sine.synthesise(nSamp: nSamp, fftrInvCfg: fftrInvCfg,
                        snOut: &snSynth, model: model, pn: pn, shift: 1)

        for i in 0..<nSamp { snSynth[i] *= gain }
        Codec2.earProtection(buffer: &snSynth, n: nSamp)

        for i in 0..<nSamp {
            let v = snSynth[i]
            let clipped: Float = v > 32767 ? 32767 : (v < -32767 ? -32767 : v)
            speech[speechOffset + i] = Int16(clipped)
        }
    }

    /// Limit the peak sample to ~30000 to mask bit-error excursions.
    /// Mirrors `ear_protection` in codec2.c.
    public static func earProtection(buffer: inout [Float], n: Int) {
        var maxSample: Float = 0
        for i in 0..<n { if buffer[i] > maxSample { maxSample = buffer[i] } }
        let over = maxSample / 30000.0
        if over > 1.0 {
            let g = 1.0 / (over * over)
            for i in 0..<n { buffer[i] *= g }
        }
    }
}
