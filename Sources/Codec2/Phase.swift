import Foundation

// Swift port of codec2/src/phase.c. Phase modelling: sample the LPC filter
// phase response at each harmonic, the zero-order phase synthesis rule, and
// minimum-phase reconstruction from a magnitude spectrum.

public enum Phase {

    /// Sample the LPC analysis-filter spectrum at each harmonic centre and
    /// store the conjugate (synthesis filter has opposite phase to analysis).
    public static func samplePhase(model: Codec2Model, H: inout [COMP], A: [COMP]) {
        let r = Float(Codec2Constants.twoPi) / Float(Codec2Constants.fftEnc)
        for m in 1...model.l {
            let b = Int(Float(m) * model.wo / r + 0.5)
            H[m] = cconj(A[b])
        }
    }

    /// Zero-order phase synthesis. For voiced frames, harmonic m's excitation
    /// phase is m * ex_phase[0]; for unvoiced frames a uniform random phase
    /// is used. The result is passed through the LPC filter and atan2'd to
    /// recover the harmonic phase.
    public static func phaseSynthZeroOrder(nSamp: Int,
                                           model: inout Codec2Model,
                                           exPhase: inout Float,
                                           H: [COMP]) {
        var ex = [COMP](repeating: COMP(), count: Codec2Constants.maxAmp + 1)
        var aBar = [COMP](repeating: COMP(), count: Codec2Constants.maxAmp + 1)

        // The C reference uses the double-precision TWO_PI macro inside this
        // wrap, so the intermediate division and subtraction promote to
        // double. Doing the same in Swift matches the C bit-for-bit.
        let twoPiD: Double = Codec2Constants.twoPi
        exPhase += model.wo * Float(nSamp)
        let arg = Float(Double(exPhase) / twoPiD + 0.5)
        let floored = floorf(arg)
        exPhase = Float(Double(exPhase) - twoPiD * Double(floored))

        for m in 1...model.l {
            if model.voiced != 0 {
                ex[m] = COMP(cosf(exPhase * Float(m)), sinf(exPhase * Float(m)))
            } else {
                let phi = Float(Codec2Constants.twoPi) *
                          Float(Sine.codec2Rand()) / Float(Sine.codec2RandMax)
                ex[m] = COMP(cosf(phi), sinf(phi))
            }
            aBar[m] = COMP(H[m].real * ex[m].real - H[m].imag * ex[m].imag,
                           H[m].imag * ex[m].real + H[m].real * ex[m].imag)
            model.phi[m] = atan2f(aBar[m].imag, aBar[m].real + 1e-12)
        }
    }

    /// Given a dB magnitude spectrum, compute the corresponding minimum-phase
    /// phase spectrum via the real cepstrum trick from
    /// http://www.dsprelated.com/showcode/20.php (and octave/mag_to_phase.m).
    /// `phase` and `gdbfk` are both length Nfft/2+1.
    public static func magToPhase(phase: inout [Float],
                                  gdbfk: [Float],
                                  Nfft: Int,
                                  fftFwdCfg: Codec2FFTConfig,
                                  fftInvCfg: Codec2FFTConfig) {
        let ns = Nfft / 2 + 1
        var sdb = [COMP](repeating: COMP(), count: Nfft)
        var c = [COMP](repeating: COMP(), count: Nfft)
        var cf = [COMP](repeating: COMP(), count: Nfft)
        var cF = [COMP](repeating: COMP(), count: Nfft)

        // Install negative-frequency mirror; 1/Nfft is the missing ifft scaling.
        sdb[0] = COMP(gdbfk[0], 0)
        for i in 1..<ns {
            sdb[i] = COMP(gdbfk[i], 0)
            sdb[Nfft - i] = COMP(gdbfk[i], 0)
        }

        Codec2FFT.fft(fftInvCfg, input: sdb, output: &c)
        let invN: Float = 1.0 / Float(Nfft)
        for i in 0..<Nfft {
            c[i] = COMP(c[i].real * invN, c[i].imag * invN)
        }

        // Fold cepstrum so non-min-phase zeros end up inside the unit circle.
        cf[0] = c[0]
        for i in 1..<(ns - 1) {
            cf[i] = cadd(c[i], c[Nfft - i])
        }
        cf[ns - 1] = c[ns - 1]
        for i in ns..<Nfft { cf[i] = COMP() }

        Codec2FFT.fft(fftFwdCfg, input: cf, output: &cF)

        // log(x) vs 20*log10(x) scaling.
        let scale: Float = 20.0 / logf(10.0)
        for i in 0..<ns { phase[i] = cF[i].imag / scale }
    }
}
