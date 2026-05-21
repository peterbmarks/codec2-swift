import Foundation

// Swift port of codec2/src/kiss_fftr.c. Real-input/output FFT layered on top
// of a half-size complex FFT, using the standard Borgerding super-twiddle
// trick. Float mode only — codec2 doesn't build the fixed-point path.

public final class KissFFTRConfig {
    public let nfft: Int          // real length (must be even)
    public let ncfft: Int         // half-size complex FFT length
    public let inverse: Bool
    public let substate: KissFFTConfig
    public var superTwiddles: [COMP]  // ncfft/2 entries

    public init(nfft: Int, inverse: Bool) {
        precondition((nfft & 1) == 0, "Real FFT optimization must be even.")
        self.nfft = nfft
        self.inverse = inverse
        self.ncfft = nfft >> 1
        self.substate = KissFFTConfig(nfft: ncfft, inverse: inverse)

        var twids = [COMP](repeating: COMP(), count: ncfft / 2)
        let pi = 3.14159265358979323846264338327
        for i in 0..<(ncfft / 2) {
            var phase = -pi * (Double(i + 1) / Double(ncfft) + 0.5)
            if inverse { phase = -phase }
            // The C code casts the double phase to float before calling cosf/sinf
            // via the kf_cexp macro. Match the precision step exactly.
            let f = Float(phase)
            twids[i] = COMP(cosf(f), sinf(f))
        }
        self.superTwiddles = twids
    }
}

public enum KissFFTR {

    /// Forward real-to-complex FFT. `timedata` length == nfft, `freqdata` length == nfft/2 + 1.
    public static func forward(_ st: KissFFTRConfig, timedata: [Float], freqdata: inout [COMP]) {
        precondition(!st.inverse, "kiss_fftr called with inverse cfg")
        let ncfft = st.ncfft

        // Pack real time series as a complex sequence of length ncfft.
        var packed = [COMP](repeating: COMP(), count: ncfft)
        for k in 0..<ncfft {
            packed[k] = COMP(timedata[2 * k], timedata[2 * k + 1])
        }
        var tmpbuf = [COMP](repeating: COMP(), count: ncfft)
        KissFFT.fft(st.substate, packed, &tmpbuf)

        // DC + Nyquist split.
        let tdc = tmpbuf[0]
        freqdata[0]     = COMP(tdc.real + tdc.imag, 0)
        freqdata[ncfft] = COMP(tdc.real - tdc.imag, 0)

        if ncfft >= 2 {
            for k in 1...(ncfft / 2) {
                let fpk = tmpbuf[k]
                let fpnk = COMP(tmpbuf[ncfft - k].real, -tmpbuf[ncfft - k].imag)
                let f1k = COMP(fpk.real + fpnk.real, fpk.imag + fpnk.imag)
                let f2k = COMP(fpk.real - fpnk.real, fpk.imag - fpnk.imag)
                let tw0 = st.superTwiddles[k - 1]
                let tw = COMP(f2k.real * tw0.real - f2k.imag * tw0.imag,
                              f2k.real * tw0.imag + f2k.imag * tw0.real)

                freqdata[k]         = COMP((f1k.real + tw.real) * 0.5,
                                           (f1k.imag + tw.imag) * 0.5)
                freqdata[ncfft - k] = COMP((f1k.real - tw.real) * 0.5,
                                           (tw.imag - f1k.imag) * 0.5)
            }
        }
    }

    /// Inverse complex-to-real FFT. `freqdata` length == nfft/2 + 1, `timedata` length == nfft.
    public static func inverse(_ st: KissFFTRConfig, freqdata: [COMP], timedata: inout [Float]) {
        precondition(st.inverse, "kiss_fftri called with non-inverse cfg")
        let ncfft = st.ncfft

        var tmpbuf = [COMP](repeating: COMP(), count: ncfft)
        tmpbuf[0] = COMP(freqdata[0].real + freqdata[ncfft].real,
                         freqdata[0].real - freqdata[ncfft].real)

        if ncfft >= 2 {
            for k in 1...(ncfft / 2) {
                let fk = freqdata[k]
                let fnkc = COMP(freqdata[ncfft - k].real, -freqdata[ncfft - k].imag)
                let fek = COMP(fk.real + fnkc.real, fk.imag + fnkc.imag)
                let tmp = COMP(fk.real - fnkc.real, fk.imag - fnkc.imag)
                let tw0 = st.superTwiddles[k - 1]
                let fok = COMP(tmp.real * tw0.real - tmp.imag * tw0.imag,
                               tmp.real * tw0.imag + tmp.imag * tw0.real)
                tmpbuf[k] = COMP(fek.real + fok.real, fek.imag + fok.imag)
                tmpbuf[ncfft - k] = COMP(fek.real - fok.real, -(fek.imag - fok.imag))
            }
        }

        var packed = [COMP](repeating: COMP(), count: ncfft)
        KissFFT.fft(st.substate, tmpbuf, &packed)
        for k in 0..<ncfft {
            timedata[2 * k]     = packed[k].real
            timedata[2 * k + 1] = packed[k].imag
        }
    }
}
