import Foundation
import Codec2

// Swift port of codec2/src/c2enc.c. Streams int16 mono speech samples and
// writes packed codec2 bits. Supports only the subset of options that the
// Codec2 Swift library currently implements (no header, no softdec, no
// codebook loading); flags are surfaced in the usage text but not honoured.

let argv = CommandLine.arguments

guard argv.count >= 4 else {
    FileHandle.standardError.write(Data("""
        usage: c2enc 3200|2400|1600|1400|1300|1200|700C InputRawspeechFile OutputBitFile
        e.g. c2enc 3200 ../raw/hts1a.raw hts1a.bit

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
        FileHandle.standardError.write(Data("Error opening input speech file: \(argv[2])\n".utf8))
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
        FileHandle.standardError.write(Data("Error opening output bit file: \(argv[3])\n".utf8))
        exit(1)
    }
    outHandle = h
}

let nsam = codec.samplesPerFrame
let nbyte = codec.bytesPerFrame
var buf = [Int16](repeating: 0, count: nsam)
var bits = [UInt8](repeating: 0, count: nbyte)

while true {
    let bytesPerSample = MemoryLayout<Int16>.size
    let needed = nsam * bytesPerSample
    let data = inHandle.readData(ofLength: needed)
    if data.count < needed { break }
    data.withUnsafeBytes { raw in
        let ptr = raw.bindMemory(to: Int16.self)
        for i in 0..<nsam { buf[i] = ptr[i] }
    }
    codec.encode(speech: buf, bits: &bits)
    outHandle.write(Data(bits))
}

if argv[2] != "-" { try? inHandle.close() }
if argv[3] != "-" { try? outHandle.close() }
