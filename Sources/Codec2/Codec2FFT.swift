import Foundation

// Swift port of codec2/src/codec2_fft.c. In the C reference this file is a thin
// dispatch layer between codec2 callers and Kiss FFT (USE_KISS_FFT) or CMSIS
// DSP (FDV_ARM_MATH). Codec2's desktop build always uses the Kiss path; we
// match that here so callers can write `Codec2FFT.fft(...)` instead of
// reaching for KissFFT directly.

public typealias Codec2FFTConfig = KissFFTConfig
public typealias Codec2FFTRConfig = KissFFTRConfig

public enum Codec2FFT {
    @inlinable
    public static func fft(_ cfg: Codec2FFTConfig, input: [COMP], output: inout [COMP]) {
        KissFFT.fft(cfg, input, &output)
    }

    @inlinable
    public static func fftInPlace(_ cfg: Codec2FFTConfig, buffer: inout [COMP]) {
        var tmp = buffer
        KissFFT.fft(cfg, buffer, &tmp)
        buffer = tmp
    }

    @inlinable
    public static func fftr(_ cfg: Codec2FFTRConfig, time: [Float], freq: inout [COMP]) {
        KissFFTR.forward(cfg, timedata: time, freqdata: &freq)
    }

    @inlinable
    public static func fftri(_ cfg: Codec2FFTRConfig, freq: [COMP], time: inout [Float]) {
        KissFFTR.inverse(cfg, freqdata: freq, timedata: &time)
    }
}
