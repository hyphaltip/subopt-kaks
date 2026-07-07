//! Protein-guided codon alignment for a pair of CDS sequences.
//!
//! CDS are translated, the proteins are aligned with an affine-gap Gotoh algorithm
//! (global or local), and the codon alignment is recovered by back-translation
//! (the role of `mrtrans` in the C code). The aligner is a standard implementation
//! and is not expected to reproduce the C Zuker aligner byte-for-byte.

use crate::genetic_code::{base_index, translate};
use crate::matrix::ScoreMatrix;

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mode {
    Global,
    Local,
}

/// Align two amino-acid strings with affine gaps.
///
/// Returns the aligned strings (with `-` gaps) plus the residue offset in `a` and `b`
/// of the first aligned column. For global alignment both offsets are 0; for local
/// alignment they mark where the aligned block begins, which the caller needs to keep
/// the back-translation in frame.
fn align_proteins(
    a: &[u8],
    b: &[u8],
    mat: &ScoreMatrix,
    mode: Mode,
    gap_open: i32,
    gap_ext: i32,
) -> (Vec<u8>, Vec<u8>, usize, usize) {
    let n = a.len();
    let m = b.len();
    let neg = i32::MIN / 4;
    // M/Ix/Iy matrices
    let idx = |i: usize, j: usize| i * (m + 1) + j;
    let mut mm = vec![neg; (n + 1) * (m + 1)];
    let mut ix = vec![neg; (n + 1) * (m + 1)]; // gap in b (consume a)
    let mut iy = vec![neg; (n + 1) * (m + 1)]; // gap in a (consume b)
    // traceback: 0=M,1=Ix,2=Iy
    let mut tb = vec![0u8; (n + 1) * (m + 1) * 3];

    mm[idx(0, 0)] = 0;
    if mode == Mode::Global {
        for i in 1..=n {
            ix[idx(i, 0)] = -(gap_open + (i as i32 - 1) * gap_ext);
        }
        for j in 1..=m {
            iy[idx(0, j)] = -(gap_open + (j as i32 - 1) * gap_ext);
        }
    } else {
        for i in 0..=n {
            mm[idx(i, 0)] = 0;
        }
        for j in 0..=m {
            mm[idx(0, j)] = 0;
        }
    }

    for i in 1..=n {
        for j in 1..=m {
            let s = mat.score(a[i - 1], b[j - 1]);
            // M
            let diag = *[mm[idx(i - 1, j - 1)], ix[idx(i - 1, j - 1)], iy[idx(i - 1, j - 1)]]
                .iter()
                .max()
                .unwrap();
            let mut best_m = diag + s;
            let mut src_m = if diag == mm[idx(i - 1, j - 1)] {
                0
            } else if diag == ix[idx(i - 1, j - 1)] {
                1
            } else {
                2
            };
            if mode == Mode::Local && best_m < 0 {
                best_m = 0;
                src_m = 0;
            }
            mm[idx(i, j)] = best_m;
            tb[idx(i, j) * 3] = src_m;
            // Ix (gap in b)
            let open_x = mm[idx(i - 1, j)] - gap_open;
            let ext_x = ix[idx(i - 1, j)] - gap_ext;
            if open_x >= ext_x {
                ix[idx(i, j)] = open_x;
                tb[idx(i, j) * 3 + 1] = 0;
            } else {
                ix[idx(i, j)] = ext_x;
                tb[idx(i, j) * 3 + 1] = 1;
            }
            // Iy (gap in a)
            let open_y = mm[idx(i, j - 1)] - gap_open;
            let ext_y = iy[idx(i, j - 1)] - gap_ext;
            if open_y >= ext_y {
                iy[idx(i, j)] = open_y;
                tb[idx(i, j) * 3 + 2] = 0;
            } else {
                iy[idx(i, j)] = ext_y;
                tb[idx(i, j) * 3 + 2] = 2;
            }
        }
    }

    // find start of traceback
    let (mut i, mut j, mut state);
    if mode == Mode::Global {
        i = n;
        j = m;
        let ends = [mm[idx(n, m)], ix[idx(n, m)], iy[idx(n, m)]];
        state = ends.iter().enumerate().max_by_key(|(_, v)| **v).unwrap().0 as u8;
    } else {
        // best M cell
        let (mut bi, mut bj, mut bv) = (0usize, 0usize, 0i32);
        for x in 0..=n {
            for y in 0..=m {
                if mm[idx(x, y)] >= bv {
                    bv = mm[idx(x, y)];
                    bi = x;
                    bj = y;
                }
            }
        }
        i = bi;
        j = bj;
        state = 0;
    }

    let mut ra = Vec::new();
    let mut rb = Vec::new();
    while i > 0 || j > 0 {
        if mode == Mode::Local && state == 0 && mm[idx(i, j)] == 0 {
            break;
        }
        match state {
            0 if i > 0 && j > 0 => {
                ra.push(a[i - 1]);
                rb.push(b[j - 1]);
                state = tb[idx(i, j) * 3];
                i -= 1;
                j -= 1;
            }
            // At an edge in state M (global): flush the remaining prefix as terminal gaps
            // so the traceback reaches (0,0) and every residue is accounted for.
            0 if i > 0 => {
                ra.push(a[i - 1]);
                rb.push(b'-');
                i -= 1;
            }
            0 => {
                ra.push(b'-');
                rb.push(b[j - 1]);
                j -= 1;
            }
            1 => {
                ra.push(a[i - 1]);
                rb.push(b'-');
                let next = tb[idx(i, j) * 3 + 1];
                i -= 1;
                state = if next == 0 { 0 } else { 1 };
            }
            _ => {
                ra.push(b'-');
                rb.push(b[j - 1]);
                let next = tb[idx(i, j) * 3 + 2];
                j -= 1;
                state = if next == 0 { 0 } else { 2 };
            }
        }
    }
    ra.reverse();
    rb.reverse();
    // `i`/`j` now hold the number of leading residues in `a`/`b` left out of the alignment.
    (ra, rb, i, j)
}

/// Result of aligning a pair of CDS: aligned codon (nucleotide) strings.
pub struct CodonAlignment {
    pub names: [String; 2],
    pub aln: [String; 2],
}

/// Translate, align proteins, and back-translate to a codon alignment.
pub fn align_cds_pair(
    name1: &str,
    cds1: &str,
    name2: &str,
    cds2: &str,
    code: usize,
    mat: &ScoreMatrix,
    mode: Mode,
    gap_open: i32,
    gap_ext: i32,
) -> Result<CodonAlignment, String> {
    let b1: Vec<u8> = cds1.bytes().filter_map(base_index).collect();
    let b2: Vec<u8> = cds2.bytes().filter_map(base_index).collect();
    let mut p1 = translate(&b1, code);
    let mut p2 = translate(&b2, code);
    // drop a trailing stop, matching the tool's expectation of stop-stripped CDS
    if p1.last() == Some(&b'*') {
        p1.pop();
    }
    if p2.last() == Some(&b'*') {
        p2.pop();
    }

    let (a1, a2, off1, off2) = align_proteins(&p1, &p2, mat, mode, gap_open, gap_ext);

    let mut out1 = String::new();
    let mut out2 = String::new();
    // Back-translate using ungapped codon indices into the original CDS. Local
    // alignments start partway in, so the codon cursors begin at the block offsets.
    let mut codon_index1 = off1;
    let mut codon_index2 = off2;
    let ungap1: Vec<u8> = cds1.bytes().filter(|c| base_index(*c).is_some()).collect();
    let ungap2: Vec<u8> = cds2.bytes().filter(|c| base_index(*c).is_some()).collect();
    for k in 0..a1.len() {
        if a1[k] == b'-' {
            out1.push_str("---");
        } else {
            let s = codon_index1 * 3;
            out1.push_str(std::str::from_utf8(&ungap1[s..s + 3]).unwrap_or("NNN"));
            codon_index1 += 1;
        }
        if a2[k] == b'-' {
            out2.push_str("---");
        } else {
            let s = codon_index2 * 3;
            out2.push_str(std::str::from_utf8(&ungap2[s..s + 3]).unwrap_or("NNN"));
            codon_index2 += 1;
        }
    }

    Ok(CodonAlignment {
        names: [name1.to_string(), name2.to_string()],
        aln: [out1, out2],
    })
}
