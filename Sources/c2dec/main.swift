import Foundation
import Codec2

// Swift port of codec2/src/c2dec.c. Reads packed bit frames and writes
// reconstructed int16 PCM samples. Only the headerless / no-bit-error
// path is currently implemented; --softdec, --ber, --error_pattern etc.
// from the C tool will follow when the modes that exercise them land.

let argv = CommandLine.arguments

guard argv.count >= 4 else {
    FileHandle.standardError.write(Data("""
        usage: c2dec 3200|2400|1600|1400|1300|1200|700C InputBitFile OutputRawSpeechFile
        e.g. c2dec 3200 hts1a.bit hts1a_3200.raw

        """.utf8))
    exit(1)
}

guard let mode = Codec2Mode.fromString(argv[1]) else {
    FileHandle.standardError.write(Data("Error: mode must be 3200, 2400, 1600, 1400, 1300, 1200, 700C\n".utf8))
    exit(1)
}

guard let codec = Codec2(mode: mode) else {
    FileHandle.standardError.write(Data("Error: could not initialise Codec 2\n".utf8))
    exit(1)
}

let inHandle: FileHandle
if argv[2] == "-" {
    inHandle = FileHandle.standardInput
} else {
    guard let h = FileHandle(forReadingAtPath: argv[2]) else {
        FileHandle.standardError.write(Data("Error opening input bit file: \(argv[2])\n".utf8))
        exit(1)
    }
    inHandle = h
}

let outURL = URL(fileURLWithPath: argv[3])
if argv[3] != "-" {
    FileManager.default.createFile(atPath: argv[3], contents: nil, attributes: nil)
}
let outHandle: FileHandle
if argv[3] == "-" {
    outHandle = FileHandle.standardOutput
} else {
    guard let h = FileHandle(forWritingAtPath: outURL.path) else {
        FileHandle.standardError.write(Data("Error opening output speech file: \(argv[3])\n".utf8))
        exit(1)
    }
    outHandle = h
}

let nsam = codec.samplesPerFrame
let nbyte = codec.bytesPerFrame
var bits = [UInt8](repeating: 0, count: nbyte)
var speech = [Int16](repeating: 0, count: nsam)

while true {
    let data = inHandle.readData(ofLength: nbyte)
    if data.count < nbyte { break }
    for i in 0..<nbyte { bits[i] = data[i] }
    codec.decode(bits: bits, speech: &speech)

    speech.withUnsafeBufferPointer { ptr in
        outHandle.write(Data(buffer: ptr))
    }
}

if argv[2] != "-" { try? inHandle.close() }
if argv[3] != "-" { try? outHandle.close() }
