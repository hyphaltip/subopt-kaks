//! Genetic codes and codon encoding, matching the encoding used by PAML/yn00.
//!
//! Nucleotides are encoded T=0, C=1, A=2, G=3 (the `BASEs` order in tools.c).
//! A codon is encoded as `base0*16 + base1*4 + base2`, giving an index 0..=63.
//! `GeneticCode[icode][codon]` is the amino-acid index (0..=19) or -1 for a stop.

/// Amino-acid single-letter codes indexed by the values stored in [`GENETIC_CODES`].
pub const AAS: &[u8; 20] = b"ARNDCQEGHILKMFPSTWYV";

/// Position weights so that `codon = b0*BY[0] + b1*BY[1] + b2*BY[2]`.
pub const BY: [usize; 3] = [16, 4, 1];

/// The subset of NCBI translation tables carried over from PAML (`GeneticCode` in tools.c).
/// Index by codon 0..=63; value is the amino-acid index or -1 for a stop codon.
pub const GENETIC_CODES: [[i8; 64]; 6] = [
    // 0: universal / standard
    [
        13, 13, 10, 10, 15, 15, 15, 15, 18, 18, -1, -1, 4, 4, -1, 17, 10, 10, 10, 10, 14, 14, 14,
        14, 8, 8, 5, 5, 1, 1, 1, 1, 9, 9, 9, 12, 16, 16, 16, 16, 2, 2, 11, 11, 15, 15, 1, 1, 19,
        19, 19, 19, 0, 0, 0, 0, 3, 3, 6, 6, 7, 7, 7, 7,
    ],
    // 1: vertebrate mitochondrial
    [
        13, 13, 10, 10, 15, 15, 15, 15, 18, 18, -1, -1, 4, 4, 17, 17, 10, 10, 10, 10, 14, 14, 14,
        14, 8, 8, 5, 5, 1, 1, 1, 1, 9, 9, 12, 12, 16, 16, 16, 16, 2, 2, 11, 11, 15, 15, -1, -1, 19,
        19, 19, 19, 0, 0, 0, 0, 3, 3, 6, 6, 7, 7, 7, 7,
    ],
    // 2: yeast mitochondrial
    [
        13, 13, 10, 10, 15, 15, 15, 15, 18, 18, -1, -1, 4, 4, -1, 17, 16, 16, 16, 16, 14, 14, 14,
        14, 8, 8, 5, 5, 1, 1, 1, 1, 9, 9, 12, 12, 16, 16, 16, 16, 2, 2, 11, 11, 15, 15, 1, 1, 19,
        19, 19, 19, 0, 0, 0, 0, 3, 3, 6, 6, 7, 7, 7, 7,
    ],
    // 3: mold, protozoan, and coelenterate mitochondrial
    [
        13, 13, 10, 10, 15, 15, 15, 15, 18, 18, -1, -1, 4, 4, 17, 17, 10, 10, 10, 10, 14, 14, 14,
        14, 8, 8, 5, 5, 1, 1, 1, 1, 9, 9, 9, 12, 16, 16, 16, 16, 2, 2, 11, 11, 15, 15, 1, 1, 19,
        19, 19, 19, 0, 0, 0, 0, 3, 3, 6, 6, 7, 7, 7, 7,
    ],
    // 4: invertebrate mitochondrial
    [
        13, 13, 10, 10, 15, 15, 15, 15, 18, 18, -1, -1, 4, 4, 17, 17, 10, 10, 10, 10, 14, 14, 14,
        14, 8, 8, 5, 5, 1, 1, 1, 1, 9, 9, 12, 12, 16, 16, 16, 16, 2, 2, 11, 11, 15, 15, 15, 15, 19,
        19, 19, 19, 0, 0, 0, 0, 3, 3, 6, 6, 7, 7, 7, 7,
    ],
    // 5: ciliate nuclear
    [
        13, 13, 10, 10, 15, 15, 15, 15, 18, 18, 5, 5, 4, 4, -1, 17, 10, 10, 10, 10, 14, 14, 14, 14,
        8, 8, 5, 5, 1, 1, 1, 1, 9, 9, 9, 12, 16, 16, 16, 16, 2, 2, 11, 11, 15, 15, 1, 1, 19, 19, 19,
        19, 0, 0, 0, 0, 3, 3, 6, 6, 7, 7, 7, 7,
    ],
];

/// Number of supported genetic codes.
pub const NUM_CODES: usize = GENETIC_CODES.len();

/// Map a nucleotide character to its base index (T/U=0, C=1, A=2, G=3).
/// Returns `None` for gaps or ambiguous/unknown characters.
#[inline]
pub fn base_index(c: u8) -> Option<u8> {
    match c.to_ascii_uppercase() {
        b'T' | b'U' => Some(0),
        b'C' => Some(1),
        b'A' => Some(2),
        b'G' => Some(3),
        _ => None,
    }
}

/// Amino acid for a codon index under the given code, or `'*'` for a stop, `'X'` if invalid.
pub fn amino_acid(code: usize, codon: usize) -> u8 {
    if codon >= 64 {
        return b'X';
    }
    match GENETIC_CODES[code][codon] {
        -1 => b'*',
        aa => AAS[aa as usize],
    }
}

/// Translate a nucleotide slice (already base-indexed, no gaps) into amino-acid codes.
/// Trailing partial codons are ignored.
pub fn translate(bases: &[u8], code: usize) -> Vec<u8> {
    bases
        .chunks_exact(3)
        .map(|c| {
            let codon = c[0] as usize * 16 + c[1] as usize * 4 + c[2] as usize;
            amino_acid(code, codon)
        })
        .collect()
}
