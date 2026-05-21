import Foundation

// Direct Swift port of codec2/src/defines.h.
// Constant values are kept identical to the C reference so floating-point
// arithmetic matches bit-for-bit where possible.

public enum Codec2Constants {
    public static let nSeconds: Double = 0.01      // internal proc frame length in secs
    public static let twSeconds: Double = 0.005    // trapezoidal synth window overlap
    public static let maxAmp: Int = 160            // maximum number of harmonics
    public static let pi: Double = 3.141592654
    public static let twoPi: Double = 6.283185307
    public static let mPiFloat: Float = 3.14159265358979323846
    public static let maxStr: Int = 2048

    public static let fftEnc: Int = 512            // FFT size used for encoder
    public static let fftDec: Int = 512            // FFT size used in decoder
    public static let vThresh: Float = 6.0         // voicing threshold in dB
    public static let lpcOrd: Int = 10             // LPC order
    public static let lpcOrdLow: Int = 6           // LPC order for lower rates

    public static let mPitchSeconds: Double = 0.0400
    public static let pMinSeconds: Double = 0.0025
    public static let pMaxSeconds: Double = 0.0200
}

/// Constants derived at runtime from the sample rate, mirroring `C2CONST`.
public struct C2Const {
    public var fs: Int          // sample rate
    public var nSamp: Int       // samples per 10ms frame at fs
    public var maxAmp: Int      // maximum number of harmonics
    public var mPitch: Int      // pitch estimation window size in samples
    public var pMin: Int        // minimum pitch period in samples
    public var pMax: Int        // maximum pitch period in samples
    public var woMin: Float
    public var woMax: Float
    public var nw: Int          // analysis window size in samples
    public var tw: Int          // trapezoidal synthesis window overlap
}

/// Sinusoidal model parameters for a single 10ms frame, mirroring `MODEL`.
public struct Codec2Model {
    public var wo: Float = 0                              // fundamental frequency, radians
    public var l: Int = 0                                 // number of harmonics
    public var a: [Float]                                 // amplitude of each harmonic, size maxAmp+1
    public var phi: [Float]                               // phase of each harmonic, size maxAmp+1
    public var voiced: Int32 = 0                          // non-zero if voiced

    public init(maxAmp: Int = Codec2Constants.maxAmp) {
        self.a = Array(repeating: 0, count: maxAmp + 1)
        self.phi = Array(repeating: 0, count: maxAmp + 1)
    }
}

/// Vector quantisation codebook entry, mirroring `struct lsp_codebook`.
public struct LSPCodebook {
    public let k: Int        // dimension of vector
    public let log2m: Int    // number of bits in m
    public let m: Int        // elements in codebook
    public let cb: [Float]   // flattened m*k entries

    public init(k: Int, log2m: Int, m: Int, cb: [Float]) {
        self.k = k
        self.log2m = log2m
        self.m = m
        self.cb = cb
    }
}

@inlinable
public func pow10f(_ x: Float) -> Float {
    // Matches the codec's POW10F definition exactly.
    return Foundation.expf(2.302585092994046 as Float * x)
}
