import Foundation

// Swift port of codec2/src/lpc.c. Constants and algorithms match the C
// reference (1990-vintage Makhoul Levinson-Durbin, ALPHA/BETA pre/de-emphasis).

public enum LPC {
    public static let maxOrder: Int = 20
    private static let alpha: Float = 1.0
    private static let beta: Float = 0.94

    /// Pre-emphasise (high-pass) a frame of speech samples.
    /// `mem` is a single-sample state, the previous Sn[-1].
    public static func preEmp(_ snPre: inout [Float], sn: [Float], mem: inout Float, nSam: Int) {
        for i in 0..<nSam {
            snPre[i] = sn[i] - alpha * mem
            mem = sn[i]
        }
    }

    /// De-emphasis (low-pass) inverse of `preEmp`.
    public static func deEmp(_ snDe: inout [Float], sn: [Float], mem: inout Float, nSam: Int) {
        for i in 0..<nSam {
            snDe[i] = sn[i] + beta * mem
            mem = snDe[i]
        }
    }

    /// Hanning window of length `nSam` applied sample-by-sample.
    public static func hanningWindow(_ sn: [Float], wn: inout [Float], nSam: Int) {
        let denom = Float(nSam - 1)
        for i in 0..<nSam {
            // matches: 0.5 - 0.5*cosf(2*PI*i/(Nsam-1)) where PI is 3.141592654
            let w = 0.5 - 0.5 * cosf(2.0 * Float(Codec2Constants.pi) * Float(i) / denom)
            wn[i] = sn[i] * w
        }
    }

    /// First `order+1` autocorrelation coefficients of `Sn`.
    public static func autocorrelate(_ sn: [Float], rn: inout [Float], nSam: Int, order: Int) {
        for j in 0...order {
            var acc: Float = 0
            let end = nSam - j
            var i = 0
            while i < end {
                acc += sn[i] * sn[i + j]
                i += 1
            }
            rn[j] = acc
        }
    }

    /// Levinson-Durbin recursion. Produces `lpcs[0..order]` with `lpcs[0] = 1.0`.
    public static func levinsonDurbin(_ r: [Float], lpcs: inout [Float], order: Int) {
        // 2-D scratch laid out [order+1][order+1] to match the C VLA.
        let stride = order + 1
        var a = [Float](repeating: 0, count: stride * stride)
        var e: Float = r[0]

        for i in 1...order {
            var sum: Float = 0
            if i >= 2 {
                for j in 1...(i - 1) {
                    sum += a[(i - 1) * stride + j] * r[i - j]
                }
            }
            var k: Float = -1.0 * (r[i] + sum) / e
            if fabsf(k) > 1.0 { k = 0.0 }

            a[i * stride + i] = k
            if i >= 2 {
                for j in 1...(i - 1) {
                    a[i * stride + j] = a[(i - 1) * stride + j] + k * a[(i - 1) * stride + (i - j)]
                }
            }
            e *= (1 - k * k)
        }

        for i in 1...order {
            lpcs[i] = a[order * stride + i]
        }
        lpcs[0] = 1.0
    }

    /// Inverse filter A(z). `sn` must contain `order` samples of memory
    /// preceding the frame at offset `snOffset`, mirroring the C convention
    /// where the caller passes a pointer into a larger buffer.
    public static func inverseFilter(_ sn: [Float], snOffset: Int,
                                     a: [Float], nSam: Int,
                                     res: inout [Float], resOffset: Int,
                                     order: Int) {
        for i in 0..<nSam {
            var acc: Float = 0
            for j in 0...order {
                acc += sn[snOffset + i - j] * a[j]
            }
            res[resOffset + i] = acc
        }
    }

    /// IIR synthesis filter 1/A(z). `snOut` must carry `order` samples of
    /// memory before the frame at `snOffset` (caller updates the memory).
    public static func synthesisFilter(_ res: [Float], resOffset: Int,
                                       a: [Float], nSam: Int, order: Int,
                                       snOut: inout [Float], snOffset: Int) {
        for i in 0..<nSam {
            var v: Float = res[resOffset + i] * a[0]
            for j in 1...order {
                v -= snOut[snOffset + i - j] * a[j]
            }
            snOut[snOffset + i] = v
        }
    }

    /// Combined Hanning -> autocorrelate -> Levinson-Durbin -> residual energy.
    public static func findAks(_ sn: [Float], a: inout [Float], nSam: Int, order: Int) -> Float {
        precondition(nSam < 512, "Nsam must be < LPC_MAX_N (512)")
        var wn = [Float](repeating: 0, count: nSam)
        var r = [Float](repeating: 0, count: order + 1)

        hanningWindow(sn, wn: &wn, nSam: nSam)
        autocorrelate(wn, rn: &r, nSam: nSam, order: order)
        levinsonDurbin(r, lpcs: &a, order: order)

        var energy: Float = 0
        for i in 0...order { energy += a[i] * r[i] }
        if energy < 0.0 { energy = 1e-12 }
        return energy
    }

    /// Bandwidth-expansion weighting of an LPC vector.
    public static func weight(_ ak: [Float], gamma: Float, order: Int, akw: inout [Float]) {
        for i in 1...order {
            akw[i] = ak[i] * powf(gamma, Float(i))
        }
    }
}
