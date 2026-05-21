import Foundation

// Swift port of codec2/src/interp.c.
//
// Interpolates the centre 10ms of two model frames sampled 20ms apart.
// Used by every codec2 mode at decode time: Wo, energy and LSPs are all
// transmitted at the 20ms cadence and the in-between frame is recovered by
// interpolation.

public enum Interp {

    /// Samples the log-domain amplitude envelope at arbitrary `w` (radians).
    /// Implements the linear interpolation between adjacent harmonic
    /// amplitudes from `sample_log_amp` in interp.c.
    public static func sampleLogAmp(model: Codec2Model, w: Float) -> Float {
        let pi: Float = Float(Codec2Constants.pi)
        precondition(w > 0.0)
        precondition(w <= pi)

        let m = Int(floorf(w / model.wo + 0.5))
        let f = (w - Float(m) * model.wo) / w
        precondition(f <= 1.0)

        if m < 1 {
            return f * log10f(model.a[1] + 1e-6)
        } else if (m + 1) > model.l {
            return (1.0 - f) * log10f(model.a[model.l] + 1e-6)
        } else {
            return (1.0 - f) * log10f(model.a[m] + 1e-6)
                 + f * log10f(model.a[m + 1] + 1e-6)
        }
    }

    /// Equal-weight Wo interpolation.
    public static func interpWo(_ interp: inout Codec2Model,
                                prev: Codec2Model, next: Codec2Model,
                                woMin: Float) {
        interpWo2(&interp, prev: prev, next: next, weight: 0.5, woMin: woMin)
    }

    /// Weighted Wo interpolation. `weight` ∈ [0,1] is the next-frame mix.
    public static func interpWo2(_ interp: inout Codec2Model,
                                 prev: Codec2Model, next: Codec2Model,
                                 weight: Float, woMin: Float) {
        // Trap noisy V/UV pattern.
        if interp.voiced != 0 && prev.voiced == 0 && next.voiced == 0 {
            interp.voiced = 0
        }
        if interp.voiced != 0 {
            if prev.voiced != 0 && next.voiced != 0 {
                interp.wo = (1.0 - weight) * prev.wo + weight * next.wo
            }
            if prev.voiced == 0 && next.voiced != 0 { interp.wo = next.wo }
            if prev.voiced != 0 && next.voiced == 0 { interp.wo = prev.wo }
        } else {
            interp.wo = woMin
        }
        interp.l = Int(Float(Codec2Constants.pi) / interp.wo)
    }

    /// Geometric mean energy interpolation — matches `interp_energy` (sqrt
    /// product form, mathematically equivalent to the log-domain average).
    public static func interpEnergy(_ prev: Float, _ next: Float) -> Float {
        return sqrtf(prev * next)
    }

    /// Weighted log-domain energy interpolation. Matches `interp_energy2`.
    public static func interpEnergy2(_ prev: Float, _ next: Float, weight: Float) -> Float {
        return pow10f((1.0 - weight) * log10f(prev) + weight * log10f(next))
    }

    /// Weighted LSP interpolation between `prev` and `next`. Matches
    /// `interpolate_lsp_ver2`.
    public static func interpolateLspVer2(_ interp: inout [Float],
                                          prev: [Float], next: [Float],
                                          weight: Float, order: Int) {
        for i in 0..<order {
            interp[i] = (1.0 - weight) * prev[i] + weight * next[i]
        }
    }
}
