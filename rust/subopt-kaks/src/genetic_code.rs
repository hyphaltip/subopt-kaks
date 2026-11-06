//! Genetic codes and codon encoding, matching the encoding used by PAML/yn00.
//!
//! Nucleotides are encoded T=0, C=1, A=2, G=3 (the `BASEs` order in tools.c).
//! A codon is encoded as `base0*16 + base1*4 + base2`, giving an index 0..=63.
//! This is the same TCAG order NCBI uses for its 64-character `ncbieaa`
//! strings, so the tables below are the NCBI strings verbatim.
//!
//! Source: NCBI `gc.prt` version 4.6
//! (https://ftp.ncbi.nlm.nih.gov/entrez/misc/data/gc.prt), all 27 tables.
//! Codes are selected by their NCBI translation-table id (1 = standard).

/// Amino-acid single-letter codes indexed by the values stored in [`GENETIC_CODES`].
pub const AAS: &[u8; 20] = b"ARNDCQEGHILKMFPSTWYV";

/// Position weights so that `codon = b0*BY[0] + b1*BY[1] + b2*BY[2]`.
pub const BY: [usize; 3] = [16, 4, 1];

/// One NCBI translation table.
pub struct GeneticCode {
    /// NCBI translation-table id (`transl_table` in GenBank records).
    pub id: u8,
    /// NCBI name.
    pub name: &'static str,
    /// Amino acid for each codon 0..=63 in TCAG order; `*` is a stop.
    pub aas: &'static [u8; 64],
}

/// All NCBI translation tables, sorted by id.
pub const CODES: &[GeneticCode] = &[
    GeneticCode { id: 1, name: "Standard", aas: b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG" },
    GeneticCode { id: 2, name: "Vertebrate Mitochondrial", aas: b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG" },
    GeneticCode { id: 3, name: "Yeast Mitochondrial", aas: b"FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG" },
    GeneticCode { id: 4, name: "Mold Mitochondrial; Protozoan Mitochondrial; Coelenterate Mitochondrial; Mycoplasma; Spiroplasma", aas: b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG" },
    GeneticCode { id: 5, name: "Invertebrate Mitochondrial", aas: b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG" },
    GeneticCode { id: 6, name: "Ciliate Nuclear; Dasycladacean Nuclear; Hexamita Nuclear", aas: b"FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG" },
    GeneticCode { id: 9, name: "Echinoderm Mitochondrial; Flatworm Mitochondrial", aas: b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG" },
    GeneticCode { id: 10, name: "Euplotid Nuclear", aas: b"FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG" },
    GeneticCode { id: 11, name: "Bacterial, Archaeal and Plant Plastid", aas: b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG" },
    GeneticCode { id: 12, name: "Alternative Yeast Nuclear", aas: b"FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG" },
    GeneticCode { id: 13, name: "Ascidian Mitochondrial", aas: b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG" },
    GeneticCode { id: 14, name: "Alternative Flatworm Mitochondrial", aas: b"FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG" },
    GeneticCode { id: 15, name: "Blepharisma Macronuclear", aas: b"FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG" },
    GeneticCode { id: 16, name: "Chlorophycean Mitochondrial", aas: b"FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG" },
    GeneticCode { id: 21, name: "Trematode Mitochondrial", aas: b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG" },
    GeneticCode { id: 22, name: "Scenedesmus obliquus Mitochondrial", aas: b"FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG" },
    GeneticCode { id: 23, name: "Thraustochytrium Mitochondrial", aas: b"FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG" },
    GeneticCode { id: 24, name: "Rhabdopleuridae Mitochondrial", aas: b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG" },
    GeneticCode { id: 25, name: "Candidate Division SR1 and Gracilibacteria", aas: b"FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG" },
    GeneticCode { id: 26, name: "Pachysolen tannophilus Nuclear", aas: b"FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG" },
    GeneticCode { id: 27, name: "Karyorelict Nuclear", aas: b"FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG" },
    GeneticCode { id: 28, name: "Condylostoma Nuclear", aas: b"FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG" },
    GeneticCode { id: 29, name: "Mesodinium Nuclear", aas: b"FFLLSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG" },
    GeneticCode { id: 30, name: "Peritrich Nuclear", aas: b"FFLLSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG" },
    GeneticCode { id: 31, name: "Blastocrithidia Nuclear", aas: b"FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG" },
    GeneticCode { id: 32, name: "Balanophoraceae Plastid", aas: b"FFLLSSSSYY*WCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG" },
    GeneticCode { id: 33, name: "Cephalodiscidae Mitochondrial", aas: b"FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG" },
];

/// Number of supported genetic codes.
pub const NUM_CODES: usize = CODES.len();

const fn aa_index(c: u8) -> i8 {
    if c == b'*' {
        return -1;
    }
    let mut k = 0;
    while k < AAS.len() {
        if AAS[k] == c {
            return k as i8;
        }
        k += 1;
    }
    -1
}

const fn build_table(aas: &[u8; 64]) -> [i8; 64] {
    let mut t = [0i8; 64];
    let mut i = 0;
    while i < 64 {
        t[i] = aa_index(aas[i]);
        i += 1;
    }
    t
}

const fn build_all() -> [[i8; 64]; NUM_CODES] {
    let mut all = [[0i8; 64]; NUM_CODES];
    let mut i = 0;
    while i < NUM_CODES {
        all[i] = build_table(CODES[i].aas);
        i += 1;
    }
    all
}

/// `GENETIC_CODES[index][codon]` is the amino-acid index (0..=19) or -1 for a
/// stop, where `index` is a position in [`CODES`] (see [`index_of`]).
pub const GENETIC_CODES: [[i8; 64]; NUM_CODES] = build_all();

/// Position in [`CODES`] / [`GENETIC_CODES`] of the table with NCBI id `id`.
pub fn index_of(id: u8) -> Option<usize> {
    CODES.iter().position(|c| c.id == id)
}

/// NCBI id of the table at `index`.
pub fn ncbi_id(index: usize) -> u8 {
    CODES[index].id
}

/// One line per table: `id<TAB>name`.
pub fn list() -> String {
    CODES
        .iter()
        .map(|c| format!("{:>2}\t{}\n", c.id, c.name))
        .collect()
}

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

/// Amino acid for a codon index under the given code index, or `'*'` for a stop, `'X'` if invalid.
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

#[cfg(test)]
mod tests {
    use super::*;

    /// The six tables PAML ships (`GeneticCode` in tools.c), in PAML order
    /// 0..=5, correspond to NCBI tables 1,2,3,4,5,6. They are identical except
    /// for one codon: PAML's yeast-mitochondrial table (index 2) keeps TGA as
    /// a stop, while NCBI table 3 translates TGA as Trp. NCBI is used here.
    #[test]
    fn matches_paml_tables() {
        const PAML: [[i8; 64]; 6] = [
            [13, 13, 10, 10, 15, 15, 15, 15, 18, 18, -1, -1, 4, 4, -1, 17, 10, 10, 10, 10, 14, 14, 14,
             14, 8, 8, 5, 5, 1, 1, 1, 1, 9, 9, 9, 12, 16, 16, 16, 16, 2, 2, 11, 11, 15, 15, 1, 1, 19,
             19, 19, 19, 0, 0, 0, 0, 3, 3, 6, 6, 7, 7, 7, 7],
            [13, 13, 10, 10, 15, 15, 15, 15, 18, 18, -1, -1, 4, 4, 17, 17, 10, 10, 10, 10, 14, 14, 14,
             14, 8, 8, 5, 5, 1, 1, 1, 1, 9, 9, 12, 12, 16, 16, 16, 16, 2, 2, 11, 11, 15, 15, -1, -1, 19,
             19, 19, 19, 0, 0, 0, 0, 3, 3, 6, 6, 7, 7, 7, 7],
            [13, 13, 10, 10, 15, 15, 15, 15, 18, 18, -1, -1, 4, 4, -1, 17, 16, 16, 16, 16, 14, 14, 14,
             14, 8, 8, 5, 5, 1, 1, 1, 1, 9, 9, 12, 12, 16, 16, 16, 16, 2, 2, 11, 11, 15, 15, 1, 1, 19,
             19, 19, 19, 0, 0, 0, 0, 3, 3, 6, 6, 7, 7, 7, 7],
            [13, 13, 10, 10, 15, 15, 15, 15, 18, 18, -1, -1, 4, 4, 17, 17, 10, 10, 10, 10, 14, 14, 14,
             14, 8, 8, 5, 5, 1, 1, 1, 1, 9, 9, 9, 12, 16, 16, 16, 16, 2, 2, 11, 11, 15, 15, 1, 1, 19,
             19, 19, 19, 0, 0, 0, 0, 3, 3, 6, 6, 7, 7, 7, 7],
            [13, 13, 10, 10, 15, 15, 15, 15, 18, 18, -1, -1, 4, 4, 17, 17, 10, 10, 10, 10, 14, 14, 14,
             14, 8, 8, 5, 5, 1, 1, 1, 1, 9, 9, 12, 12, 16, 16, 16, 16, 2, 2, 11, 11, 15, 15, 15, 15, 19,
             19, 19, 19, 0, 0, 0, 0, 3, 3, 6, 6, 7, 7, 7, 7],
            [13, 13, 10, 10, 15, 15, 15, 15, 18, 18, 5, 5, 4, 4, -1, 17, 10, 10, 10, 10, 14, 14, 14, 14,
             8, 8, 5, 5, 1, 1, 1, 1, 9, 9, 9, 12, 16, 16, 16, 16, 2, 2, 11, 11, 15, 15, 1, 1, 19, 19, 19,
             19, 0, 0, 0, 0, 3, 3, 6, 6, 7, 7, 7, 7],
        ];
        const TGA: usize = 3 * 4 + 2; // T=0, G=3, A=2
        for (paml_idx, ncbi_id) in [(0u8, 1u8), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6)] {
            let idx = index_of(ncbi_id).unwrap();
            let mut expected = PAML[paml_idx as usize];
            if ncbi_id == 3 {
                assert_eq!(expected[TGA], -1);
                expected[TGA] = aa_index(b'W');
            }
            assert_eq!(GENETIC_CODES[idx], expected, "NCBI table {ncbi_id}");
        }
    }

    #[test]
    fn ids_sorted_unique_and_yeast_mito_differs() {
        for w in CODES.windows(2) {
            assert!(w[0].id < w[1].id);
        }
        assert!(index_of(0).is_none());
        assert!(index_of(7).is_none());
        // CTN is Thr in yeast mitochondria (table 3), Leu in the standard code
        let ctt = 1 * 16 + 0 * 4 + 0;
        assert_eq!(amino_acid(index_of(1).unwrap(), ctt), b'L');
        assert_eq!(amino_acid(index_of(3).unwrap(), ctt), b'T');
    }
}
