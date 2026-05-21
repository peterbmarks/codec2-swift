import Foundation

// Swift port of codec2/src/postfilter.c.
//
// Background-noise adaptive postfilter. Tracks the low-level noise floor in
// unvoiced frames; in voiced frames any harmonic whose amplitude falls below
// (bg_est + BG_MARGIN) dB has its phase randomised to mask "clicky" noise.

public enum Postfilter {

    public static let bgThresh: Float = 40.0
    public static let bgBeta: Float = 0.1
    public static let bgMargin: Float = 6.0

    /// Postfilter a single 10 ms model frame, updating the running background
    /// estimate `bgEst`. Mirrors `postfilter()` exactly, including the side
    /// effect of randomising sub-threshold harmonic phases.
    public static func postfilter(_ model: inout Codec2Model, bgEst: inout Float) {
        var e: Float = 1e-12
        for m in 1...model.l {
            e += model.a[m] * model.a[m]
        }
        precondition(e > 0)
        e = 10.0 * log10f(e / Float(model.l))

        if e < bgThresh && model.voiced == 0 {
            bgEst = bgEst * (1.0 - bgBeta) + e * bgBeta
        }

        let thresh = pow10f((bgEst + bgMargin) / 20.0)
        if model.voiced != 0 {
            for m in 1...model.l {
                if model.a[m] < thresh {
                    model.phi[m] = (Float(Codec2Constants.twoPi) /
                                    Float(Sine.codec2RandMax)) *
                                   Float(Sine.codec2Rand())
                }
            }
        }
    }
}
