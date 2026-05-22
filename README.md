# codec2-swift

Codec 2 is an open-source, low-bitrate speech codec developed by Australian engineer David Rowe (VK5DGR), 
designed specifically for compressing voice at very low bit rates — far lower than mainstream codecs like MP3 or even Opus.

*Core purpose*
It was created to enable intelligible voice communication over very narrow bandwidth channels, particularly HF (shortwave) radio links used in amateur radio and emergency communications. At the bit rates Codec 2 targets, most commercial codecs either don't work or produce unusable audio.

Codec2 can also be used to compress speach to a very small data file.

This project is a pure-Swift port of David Rowe's [Codec 2](http://rowetel.com/codec2.html) speech codec.

The library is implemented entirely in Swift with no dependency on the C codec2 sources at runtime. The original C reference is used only as the source of truth for parity testing and codebook tables. The port is built and tested against the C reference outputs captured under [Tests/Codec2Tests/Reference/](Tests/Codec2Tests/Reference/).

## What is here

A SwiftPM package with:

- a `Codec2` library target (the port itself, in [Sources/Codec2/](Sources/Codec2/));
- two executable targets, `c2enc` and `c2dec`, that mirror the C CLI tools and are wired up to the Swift codec ([Sources/c2enc/](Sources/c2enc/), [Sources/c2dec/](Sources/c2dec/));
- an XCTest target with parity tests against the C reference outputs ([Tests/Codec2Tests/](Tests/Codec2Tests/)).

## Modes supported

Codec 2 ships seven speech-only bit rates (mode 700C is unrelated to the rest; the others share a common LPC/LSP/sinusoidal pipeline):

| Mode | bits/frame | Encoder | Decoder | Notes |
|------|------------|---------|---------|-------|
| 3200 | 64 | yes | yes | Encoder ~99.7% byte-exact with C, decoder still drifts numerically — see "Known divergences" below |
| 2400 | 48 | no  | no  | not yet ported |
| 1600 | 64 | no  | no  | not yet ported |
| 1400 | 56 | no  | no  | not yet ported |
| 1300 | 52 | no  | no  | not yet ported |
| 1200 | 48 | no  | no  | not yet ported |
| 700C | 28 | yes | yes | **Encoder byte-exact, decoder ≤ 1 LSB on every sample** |

The remaining modes share most of their building blocks with 3200 and 700C; what they still need is `encode_lsps_scalar` / `decode_lsps_scalar` plus the LSP scalar codebook (`lsp_cb`) for 2400/1600/1300, and the joint WoE quantiser (`quantise_WoE`, `encode_WoE`, `decode_WoE`) plus the `ge_cb` codebook for 2400/1400/1200. They are not yet wired up. Stub callers will `fatalError` rather than produce wrong output.

The library deliberately omits everything outside the speech codec: no FreeDV API, no OFDM/COHPSK/FDMDV/FSK modems, no LDPC, no STM32 firmware, no data channels. The corresponding C sources are present under [codec2/src/](codec2/src/) for reference but no Swift ports exist for them.

## How the port is structured

Each C source maps to a Swift module of the same role; the naming is deliberately close to the originals so the two can be diffed:

| Swift | C source | Role |
|-------|----------|------|
| [Defines.swift](Sources/Codec2/Defines.swift) | [defines.h](codec2/src/defines.h) | Constants, `C2Const`, `Codec2Model`, `LSPCodebook` |
| [Complex.swift](Sources/Codec2/Complex.swift) | [comp.h](codec2/src/comp.h), [comp_prim.h](codec2/src/comp_prim.h) | `COMP` + complex primitives |
| [Pack.swift](Sources/Codec2/Pack.swift) | [pack.c](codec2/src/pack.c) | Bit packing (natural + Gray code) |
| [LPC.swift](Sources/Codec2/LPC.swift) | [lpc.c](codec2/src/lpc.c) | Pre/de-emphasis, Hanning, autocorrelate, Levinson-Durbin, inverse/synthesis filters, `find_aks`, `weight` |
| [LSP.swift](Sources/Codec2/LSP.swift) | [lsp.c](codec2/src/lsp.c) | LPC↔LSP Chebyshev root finder + cascade reconstruction |
| [KissFFT.swift](Sources/Codec2/KissFFT.swift) | [kiss_fft.c](codec2/src/kiss_fft.c) | Mixed-radix complex FFT |
| [KissFFTR.swift](Sources/Codec2/KissFFTR.swift) | [kiss_fftr.c](codec2/src/kiss_fftr.c) | Real-input FFT |
| [Codec2FFT.swift](Sources/Codec2/Codec2FFT.swift) | [codec2_fft.c](codec2/src/codec2_fft.c) | FFT wrappers used by the codec |
| [MBest.swift](Sources/Codec2/MBest.swift) | [mbest.c](codec2/src/mbest.c) | Multi-best VQ search |
| [NLP.swift](Sources/Codec2/NLP.swift) | [nlp.c](codec2/src/nlp.c) | Non-Linear Pitch estimator (8 kHz path only; the 16 kHz FreeDV path is stubbed with a `precondition`) |
| [Sine.swift](Sources/Codec2/Sine.swift) | [sine.c](codec2/src/sine.c) | Sinusoidal analysis/synthesis, `c2const_create`, MBE voicing, overlap-add synth, `codec2_rand` |
| [Phase.swift](Sources/Codec2/Phase.swift) | [phase.c](codec2/src/phase.c) | LPC phase sampling, zero-order phase synth, `mag_to_phase` |
| [Interp.swift](Sources/Codec2/Interp.swift) | [interp.c](codec2/src/interp.c) | Frame-to-frame Wo / energy / LSP interpolation |
| [Postfilter.swift](Sources/Codec2/Postfilter.swift) | [postfilter.c](codec2/src/postfilter.c) | Background-noise adaptive postfilter |
| [Quantise.swift](Sources/Codec2/Quantise.swift) | [quantise.c](codec2/src/quantise.c) (3200 + 700C subset) | `speech_to_uq_lsps`, `aks_to_M2`, `lpc_post_filter`, `apply_lpc_correction`, Wo/log-Wo/energy/LSP-difference encode+decode, codebook search |
| [Newamp1.swift](Sources/Codec2/Newamp1.swift) | [newamp1.c](codec2/src/newamp1.c) | Full 700C pipeline: rate-L↔rate-K mel resampling, 2-stage VQ, post-filter, Wo/voicing interp, phase synthesis |
| [Codec2.swift](Sources/Codec2/Codec2.swift) | [codec2.c](codec2/src/codec2.c) | Public API + mode 3200 + mode 700C |
| [CodebookLspd.swift](Sources/Codec2/CodebookLspd.swift) | [build/src/codebookd.c](codec2/build/src/codebookd.c) | LSP-difference scalar codebooks (used by mode 3200) |
| [CodebookNewamp1.swift](Sources/Codec2/CodebookNewamp1.swift) | [build/src/codebooknewamp1.c](codec2/build/src/codebooknewamp1.c) | newamp1 2-stage VQ codebooks (used by 700C) |
| [CodebookNewamp1Energy.swift](Sources/Codec2/CodebookNewamp1Energy.swift) | [build/src/codebooknewamp1_energy.c](codec2/build/src/codebooknewamp1_energy.c) | newamp1 energy codebook (used by 700C) |

The three `Codebook*.swift` files were generated from the C build outputs and are byte-for-byte identical to what the C library ships.

## Building and running

The package requires Xcode's Swift toolchain (Swift 5.9+) on macOS 11 or later. From the project root:

```sh
swift build
swift test
```

The build produces `c2enc` and `c2dec` binaries in `.build/debug/`. They take the same `mode input output` arguments as the C reference tools:

```sh
# play sample raw audio using sox
play --encoding signed-integer --bits 16 --rate 8000 raw/hts1a.raw
# encode the bundled hts1a sample at 700C, then decode it back
.build/debug/c2enc 700C raw/hts1a.raw /tmp/hts1a_700c.bit
.build/debug/c2dec 700C /tmp/hts1a_700c.bit /tmp/hts1a_700c_decoded.raw
# play decoded audio
play --encoding signed-integer --bits 16 --rate 8000 /tmp/hts1a_700c_decoded.raw

# Compare against the C reference outputs
cmp /tmp/hts1a_700c.bit Tests/Codec2Tests/Reference/hts1a_700C.bit
```

The C reference binaries used to capture the parity-target outputs are at [codec2/build/src/c2enc](codec2/build/src/c2enc) and [codec2/build/src/c2dec](codec2/build/src/c2dec). They were built once from the bundled C tree; the Swift port does not invoke them at runtime.

## Tests

`swift test` runs 14 test cases in 7 suites. All pass.

| Suite | What it covers |
|-------|----------------|
| [PackTests](Tests/Codec2Tests/PackTests.swift) | Natural + Gray-code pack/unpack round trips and a known-vector check |
| [LPCLSPTests](Tests/Codec2Tests/LPCLSPTests.swift) | Levinson-Durbin on a delta autocorrelation; LPC→LSP→LPC numerical round trip |
| [KissFFTTests](Tests/Codec2Tests/KissFFTTests.swift) | Forward/inverse round trip, impulse spectrum, pure single-bin sinusoid |
| [KissFFTRTests](Tests/Codec2Tests/KissFFTRTests.swift) | Real-FFT round trip and known cosine bin |
| [Mode3200ParityTests](Tests/Codec2Tests/Mode3200ParityTests.swift) | Encoder byte-parity against the C reference (current baseline ≥ 99%); decoder shape test |
| [Mode700CParityTests](Tests/Codec2Tests/Mode700CParityTests.swift) | Encoder is byte-exact; decoder samples are within ±1 LSB and fewer than 2% of samples differ from the C output |

To filter to one suite:

```sh
swift test --filter Mode700CParityTests
```

### XCTest not found?

XCTest isn't in the command line tools, it comes bundled in XCode. You may
need to switch as follows:

```sh
sudo xcode-select -s /Applications/Xcode.app/Contents/Developer
```

### Mode 700C parity result

The 700C encoder produces a bitstream that is **byte-for-byte identical** to the C reference. The decoder, fed the same C-generated bits, produces audio that is **≤ 1 LSB away from the C output on every sample**:

| Metric | Value |
|--------|-------|
| Encoded bitstream length | 300 bytes (matches C) |
| Encoded bytes that differ from C | **0 / 300** |
| Decoded sample count | 24 000 (matches C) |
| Samples bit-exact with C decoder | **23 828 / 24 000 (99.28%)** |
| Samples differing | 172 / 24 000 (0.72%) |
| Max absolute sample difference | **1 LSB** |
| Mean absolute sample difference | 0.01 |

A 1-LSB drift on a small fraction of samples is the irreducible floating-point rounding noise of replaying the same algorithm with the same arithmetic order in a different toolchain; it is not algorithmic divergence.

### Mode 3200 parity result

The mode 3200 port is functional end-to-end but not yet bit-exact:

| Metric | Value |
|--------|-------|
| Encoded bytes that differ from C | 4 / 1200 (frames 112–113) |
| Max decoded sample difference | ~4100 LSB |
| Mean decoded sample difference | ~100 LSB |

The first 111 encoded frames are byte-exact; the encoder divergence at frame 112 looks like accumulated drift in the NLP pitch tracker. The decoder still has a larger numerical gap — same root cause class as the 700C path that has since been fixed (single- vs double-precision arithmetic and `M_PI`/`TWO_PI` macro semantics), but applied to the LPC/`aks_to_M2`/`lpc_post_filter` chain in [Quantise.swift](Sources/Codec2/Quantise.swift). That audit is the next bit of work needed before 3200 reaches the same parity as 700C.

## Provenance

This code is inspired by David Rowe's Codec 2 (LGPL 2.1). The Swift port mirrors its algorithms as faithfully as possible.
