import Foundation

// Swift port of codec2/src/pack.c.
//
// Bit-packs and unpacks integer fields into a byte array. Used by the encoder
// and decoder to serialise quantised parameters. The Gray-code variants
// preserve the C behaviour: a single-bit channel error only changes the
// decoded value by one step.
//
// Conventions inherited from the C code:
//   * bitIndex counts BITS within bitArray, not bytes.
//   * bitArray bytes must be zero-initialised before the first pack() call.
//   * Field widths are <= 8 bits when Gray coding is used.

private let wordSize: UInt32 = 8
private let indexMask: UInt32 = 0x7
private let shiftRight: UInt32 = 3

@inlinable
public func pack(_ bitArray: inout [UInt8], bitIndex: inout UInt32, field: Int32, fieldWidth: UInt32) {
    packNaturalOrGray(&bitArray, bitIndex: &bitIndex, field: field, fieldWidth: fieldWidth, gray: 1)
}

public func packNaturalOrGray(_ bitArray: inout [UInt8],
                              bitIndex: inout UInt32,
                              field rawField: Int32,
                              fieldWidth: UInt32,
                              gray: UInt32) {
    var fieldWidth = fieldWidth
    // The C code treats `field` as signed int but operates on it bitwise.
    // We do the same conversion through UInt32 for two's-complement masking.
    var field = UInt32(bitPattern: rawField)
    if gray != 0 {
        field = (field >> 1) ^ field
    }
    repeat {
        let bI = bitIndex
        let bitsLeft = wordSize - (bI & indexMask)
        let sliceWidth = bitsLeft < fieldWidth ? bitsLeft : fieldWidth
        let wordIndex = Int(bI >> shiftRight)

        let shifted = (field >> (fieldWidth - sliceWidth)) << (bitsLeft - sliceWidth)
        bitArray[wordIndex] |= UInt8(truncatingIfNeeded: shifted)

        bitIndex = bI + sliceWidth
        fieldWidth -= sliceWidth
    } while fieldWidth != 0
}

@inlinable
public func unpack(_ bitArray: [UInt8], bitIndex: inout UInt32, fieldWidth: UInt32) -> Int32 {
    unpackNaturalOrGray(bitArray, bitIndex: &bitIndex, fieldWidth: fieldWidth, gray: 1)
}

public func unpackNaturalOrGray(_ bitArray: [UInt8],
                                bitIndex: inout UInt32,
                                fieldWidth: UInt32,
                                gray: UInt32) -> Int32 {
    var fieldWidth = fieldWidth
    var field: UInt32 = 0
    repeat {
        let bI = bitIndex
        let bitsLeft = wordSize - (bI & indexMask)
        let sliceWidth = bitsLeft < fieldWidth ? bitsLeft : fieldWidth

        let byte = UInt32(bitArray[Int(bI >> shiftRight)])
        let mask = (UInt32(1) << sliceWidth) - 1
        field |= ((byte >> (bitsLeft - sliceWidth)) & mask) << (fieldWidth - sliceWidth)

        bitIndex = bI + sliceWidth
        fieldWidth -= sliceWidth
    } while fieldWidth != 0

    let t: UInt32
    if gray != 0 {
        // Gray -> binary, valid for <= 8-bit fields (matches C comment).
        var u = field ^ (field >> 8)
        u ^= (u >> 4)
        u ^= (u >> 2)
        u ^= (u >> 1)
        t = u
    } else {
        t = field
    }
    return Int32(bitPattern: t)
}
