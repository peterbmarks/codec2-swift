import Foundation

// Swift port of codec2/src/comp.h and codec2/src/comp_prim.h.
//
// codec2 uses its own COMP type rather than C99 complex.h so the same code can
// compile under non-C99 toolchains. We mirror the field layout exactly:
// real part first, imaginary part second.

public struct COMP: Equatable {
    public var real: Float
    public var imag: Float

    @inlinable
    public init(real: Float = 0, imag: Float = 0) {
        self.real = real
        self.imag = imag
    }

    @inlinable
    public init(_ real: Float, _ imag: Float) {
        self.real = real
        self.imag = imag
    }
}

@inlinable
public func cneg(_ a: COMP) -> COMP {
    COMP(-a.real, -a.imag)
}

@inlinable
public func cconj(_ a: COMP) -> COMP {
    COMP(a.real, -a.imag)
}

@inlinable
public func fcmult(_ k: Float, _ a: COMP) -> COMP {
    COMP(k * a.real, k * a.imag)
}

@inlinable
public func cmult(_ a: COMP, _ b: COMP) -> COMP {
    COMP(a.real * b.real - a.imag * b.imag,
         a.real * b.imag + a.imag * b.real)
}

@inlinable
public func cadd(_ a: COMP, _ b: COMP) -> COMP {
    COMP(a.real + b.real, a.imag + b.imag)
}

@inlinable
public func csub(_ a: COMP, _ b: COMP) -> COMP {
    COMP(a.real - b.real, a.imag - b.imag)
}

@inlinable
public func cabsolute(_ a: COMP) -> Float {
    return sqrtf(a.real * a.real + a.imag * a.imag)
}

/// Complex exponential e^{jθ} — matches the codec's `comp_exp_j`.
@inlinable
public func compExpJ(_ theta: Float) -> COMP {
    COMP(cosf(theta), sinf(theta))
}
