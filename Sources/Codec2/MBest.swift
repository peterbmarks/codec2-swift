import Foundation

// Swift port of codec2/src/mbest.c.
//
// Multistage VQ search that keeps the m best candidates from each stage.
// Used by codec2 wherever a tree-structured / multistage VQ codebook is
// searched (notably newamp1 for 700C, but generally useful).
//
// The C reference uses memmove/memcpy on a packed struct array. We use a
// Swift array of value-typed entries and an explicit shift-down loop that
// matches the same logical behaviour.

public let MBEST_STAGES = 4

public struct MBestEntry {
    public var index: [Int32]            // length MBEST_STAGES
    public var error: Float

    public init() {
        self.index = [Int32](repeating: 0, count: MBEST_STAGES)
        self.error = 1e32
    }
}

public final class MBest {
    public var entries: Int
    public var list: [MBestEntry]

    public init(entries: Int) {
        precondition(entries > 0)
        self.entries = entries
        self.list = Array(repeating: MBestEntry(), count: entries)
    }

    /// In-place scale of each codebook entry by the per-dimension weight vector
    /// `w[]`. After this, an unweighted Euclidean search produces a weighted
    /// match (see comment in mbest.c for the algebraic derivation).
    public static func precomputeWeight(cb: inout [Float], w: [Float], k: Int, m: Int) {
        for j in 0..<m {
            for i in 0..<k {
                cb[k * j + i] *= w[i]
            }
        }
    }

    /// Inserts an (index, error) pair into the sorted m-best list. The list is
    /// ordered by ascending error.
    public func insert(index: [Int32], error: Float) {
        var i = 0
        while i < entries {
            if error < list[i].error {
                // Shift entries [i, entries-2] right by one, dropping the last.
                if i + 1 < entries {
                    for j in stride(from: entries - 1, to: i, by: -1) {
                        list[j] = list[j - 1]
                    }
                }
                var newEntry = MBestEntry()
                for s in 0..<MBEST_STAGES { newEntry.index[s] = index[s] }
                newEntry.error = error
                list[i] = newEntry
                return
            }
            i += 1
        }
    }

    /// Compares `vec[]` against `m` rows of length `k` in `cb`, retaining the
    /// m-best closest matches inside `mbest`. `index[]` holds the stage path
    /// that led us here; the leaf stage's chosen row is written to `index[0]`.
    public static func search(cb: [Float], cbOffset: Int,
                              vec: [Float],
                              k: Int, m: Int,
                              mbest: MBest,
                              index: inout [Int32]) {
        var ptr = cbOffset
        for j in 0..<m {
            var e: Float = 0
            for i in 0..<k {
                let diff = cb[ptr] - vec[i]
                ptr += 1
                e += diff * diff
            }
            index[0] = Int32(j)
            if e < mbest.list[mbest.entries - 1].error {
                mbest.insert(index: index, error: e)
            }
        }
    }

    /// 450 mode variant: search only the first `shorterK` dimensions, weighting
    /// each contribution by `w[i]`. Behaviour matches `mbest_search450` from
    /// mbest.c byte-for-byte.
    public static func search450(cb: [Float], vec: [Float], w: [Float],
                                 k: Int, shorterK: Int, m: Int,
                                 mbest: MBest, index: inout [Int32]) {
        for j in 0..<m {
            var e: Float = 0
            for i in 0..<k {
                if i < shorterK {
                    let diff = cb[j * k + i] - vec[i]
                    e += diff * w[i] * diff * w[i]
                }
            }
            index[0] = Int32(j)
            mbest.insert(index: index, error: e)
        }
    }
}
