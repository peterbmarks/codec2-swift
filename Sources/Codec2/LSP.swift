import Foundation

// Swift port of codec2/src/lsp.c. The LPC <-> LSP conversion is at the heart
// of every codec2 mode that quantises spectrum, so the port reproduces the
// C arithmetic step-for-step. Pointer arithmetic is replaced with explicit
// index variables to keep the algorithm identical.

public enum LSP {

    /// Evaluates a series of Chebyshev polynomials. Matches `cheb_poly_eva`.
    private static func chebPolyEva(_ coef: [Float], coefOffset: Int, x: Float, order: Int) -> Float {
        let halfOrder = order / 2
        var t = [Float](repeating: 0, count: halfOrder + 1)
        t[0] = 1.0
        t[1] = x
        if halfOrder >= 2 {
            for i in 2...halfOrder {
                t[i] = (2 * x) * t[i - 1] - t[i - 2]
            }
        }
        var sum: Float = 0
        for i in 0...halfOrder {
            sum += coef[coefOffset + (halfOrder - i)] * t[i]
        }
        return sum
    }

    /// LPC -> LSP conversion. Returns the number of LSP roots found.
    /// - Parameters:
    ///   - a: LPC coefficients, length `order + 1`.
    ///   - order: LPC order (typically 10).
    ///   - freq: output LSP frequencies in radians, length `order`.
    ///   - nb: number of bisection sub-intervals (codec uses 4).
    ///   - delta: grid spacing (codec uses 0.02).
    @discardableResult
    public static func lpcToLsp(_ a: [Float], order: Int, freq: inout [Float], nb: Int, delta: Float) -> Int {
        let m = order / 2
        var p = [Float](repeating: 0, count: order + 1)
        var q = [Float](repeating: 0, count: order + 1)

        // Build P' and Q' from A: P'(z) = P(z)/(1 + z^-1), Q'(z) = Q(z)/(1 - z^-1).
        p[0] = 1.0
        q[0] = 1.0
        for i in 1...m {
            p[i] = a[i] + a[order + 1 - i] - p[i - 1]
            q[i] = a[i] - a[order + 1 - i] + q[i - 1]
        }
        for i in 0..<m {
            p[i] = 2 * p[i]
            q[i] = 2 * q[i]
        }

        var roots = 0
        var xr: Float = 0
        var xl: Float = 1.0
        var xm: Float = 0

        for j in 0..<order {
            let usingQ = (j % 2) != 0
            // Helper to evaluate the active polynomial.
            func eva(_ x: Float) -> Float {
                return usingQ
                    ? chebPolyEva(q, coefOffset: 0, x: x, order: order)
                    : chebPolyEva(p, coefOffset: 0, x: x, order: order)
            }

            var psuml = eva(xl)
            var flag = true
            while flag && (xr >= -1.0) {
                xr = xl - delta
                var psumr = eva(xr)
                let tempPsumr = psumr
                let tempXr = xr

                if (psumr * psuml) < 0.0 || psumr == 0.0 {
                    roots += 1
                    var psumm = psuml
                    for _ in 0...nb {
                        xm = (xl + xr) / 2
                        psumm = eva(xm)
                        if psumm * psuml > 0.0 {
                            psuml = psumm
                            xl = xm
                        } else {
                            psumr = psumm
                            xr = xm
                        }
                    }
                    freq[j] = xm
                    xl = xm
                    flag = false
                } else {
                    psuml = tempPsumr
                    xl = tempXr
                }
            }
        }

        // x-domain -> radians
        for i in 0..<order {
            freq[i] = acosf(freq[i])
        }
        return roots
    }

    /// LSP -> LPC conversion. Reconstructs ak[0..order] from LSP frequencies.
    public static func lspToLpc(_ lsp: [Float], ak: inout [Float], order: Int) {
        var freq = [Float](repeating: 0, count: order)
        for i in 0..<order { freq[i] = cosf(lsp[i]) }

        let wpLen = (order * 4) + 2
        var wp = [Float](repeating: 0, count: wpLen)

        var xin1: Float = 1.0
        var xin2: Float = 1.0
        var lastN4: Int = 0   // index into wp; needed after the inner loop, matches C's `n4`

        for j in 0...order {
            for i in 0..<(order / 2) {
                let n1 = i * 4
                let n2 = n1 + 1
                let n3 = n2 + 1
                let n4 = n3 + 1
                let xout1 = xin1 - 2 * freq[2 * i] * wp[n1] + wp[n2]
                let xout2 = xin2 - 2 * freq[2 * i + 1] * wp[n3] + wp[n4]
                wp[n2] = wp[n1]
                wp[n4] = wp[n3]
                wp[n1] = xin1
                wp[n3] = xin2
                xin1 = xout1
                xin2 = xout2
                lastN4 = n4
            }
            let xout1 = xin1 + wp[lastN4 + 1]
            let xout2 = xin2 - wp[lastN4 + 2]
            ak[j] = (xout1 + xout2) * 0.5
            wp[lastN4 + 1] = xin1
            wp[lastN4 + 2] = xin2

            xin1 = 0.0
            xin2 = 0.0
        }
    }
}
