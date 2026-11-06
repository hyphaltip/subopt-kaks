//! dN/dS estimation ported from the subopt-kaks C code (Yang & Nielsen 2000
//! single-pass estimator plus Nei & Gojobori 1986).
//!
//! This reproduces the behaviour of `yn00_cds_prealigned` / `yn00_cds_optimal`:
//! the shipped C tool runs YN00 with `com.iteration == 0`, i.e. a single pass with
//! equal weighting of substitution pathways (the codon transition matrix P(t) is
//! taken as all-ones), so no codon-model eigendecomposition is required here.

use crate::genetic_code::{AAS, BY, GENETIC_CODES};

const NCODE: usize = 64;
const ACCURACY: f64 = 1e-9;

/// Cleaned, codon-encoded alignment.
pub struct CodonData {
    pub names: Vec<String>,
    pub ns: usize,
    pub ls: usize,
    pub alen: usize,
    /// `z[i][h]` is the codon index (0..=63) of codon `h` in sequence `i`.
    pub z: Vec<Vec<u8>>,
    /// Index into `genetic_code::CODES` (not the NCBI id).
    pub code: usize,
}

/// Result of one pairwise comparison.
#[derive(Clone, Debug)]
pub struct PairResult {
    pub seq1: String,
    pub seq2: String,
    // Yang & Nielsen (single pass)
    pub dn: f64,
    pub ds: f64,
    pub omega: f64,
    pub n_sites: f64,
    pub s_sites: f64,
    pub kappa: f64,
    pub t: f64,
    pub se_dn: f64,
    pub se_ds: f64,
    // Nei & Gojobori 1986
    pub ng_dn: f64,
    pub ng_ds: f64,
    pub ng_omega: f64,
    pub length: usize,
}

/// Build codon-encoded data from an aligned set of nucleotide records.
///
/// Whole codons (triplets) containing a gap or ambiguous base in any sequence are
/// dropped — a frame-safe variant of the column dropping done by the C code (identical
/// for gapless alignments, correct for triplet-gapped ones).
pub fn clean_and_encode(
    names: Vec<String>,
    seqs: &[Vec<u8>],
    code: usize,
) -> Result<CodonData, String> {
    let ns = seqs.len();
    if ns < 2 {
        return Err("need at least 2 sequences".into());
    }
    let alen = seqs[0].len();
    if seqs.iter().any(|s| s.len() != alen) {
        return Err("aligned sequences differ in length".into());
    }
    let ncodons = alen / 3;
    let table = &GENETIC_CODES[code];

    // Base-index every character once (None = gap/ambiguous).
    let coded: Vec<Vec<Option<u8>>> = seqs
        .iter()
        .map(|s| s.iter().map(|&c| crate::genetic_code::base_index(c)).collect())
        .collect();

    let mut z: Vec<Vec<u8>> = vec![Vec::with_capacity(ncodons); ns];
    for cx in 0..ncodons {
        let base = cx * 3;
        // keep the codon only if every sequence has three clean bases here
        let keep = coded
            .iter()
            .all(|s| s[base].is_some() && s[base + 1].is_some() && s[base + 2].is_some());
        if !keep {
            continue;
        }
        for i in 0..ns {
            let b0 = coded[i][base].unwrap() as usize;
            let b1 = coded[i][base + 1].unwrap() as usize;
            let b2 = coded[i][base + 2].unwrap() as usize;
            let codon = b0 * 16 + b1 * 4 + b2;
            if table[codon] == -1 {
                return Err(format!(
                    "stop codon in sequence '{}' at codon {}",
                    names[i],
                    cx + 1
                ));
            }
            z[i].push(codon as u8);
        }
    }
    let ls = z[0].len();
    if ls == 0 {
        return Err("no usable codons after cleaning".into());
    }
    Ok(CodonData { names, ns, ls, alen, z, code })
}

/// Per-sequence 3x4 base composition table (position x base), matching `Statistics()`.
fn f3x4s(data: &CodonData) -> Vec<[f64; 12]> {
    let ls = data.ls as f64;
    data.z
        .iter()
        .map(|zi| {
            let mut f = [0.0f64; 12];
            for &codon in zi {
                let c = codon as usize;
                let b = [c / 16, (c % 16) / 4, c % 4];
                for pos in 0..3 {
                    f[pos * 4 + b[pos]] += 1.0 / ls;
                }
            }
            f
        })
        .collect()
}

/// Codon frequencies pi[64] from the averaged F3x4 table (`GetFreqs`, fcommon=1).
fn get_freqs(data: &CodonData, f3: &[[f64; 12]]) -> Result<[f64; NCODE], String> {
    let table = &GENETIC_CODES[data.code];
    let mut f3x4 = [0.0f64; 12];
    for fi in f3 {
        for k in 0..12 {
            f3x4[k] += fi[k] / data.ns as f64;
        }
    }
    let mut pi = [0.0f64; NCODE];
    let mut fstop = 0.0;
    for j in 0..NCODE {
        let v = f3x4[j / 16] * f3x4[4 + (j % 16) / 4] * f3x4[8 + j % 4];
        if table[j] == -1 {
            fstop += v;
        } else {
            pi[j] = v;
        }
    }
    for v in pi.iter_mut() {
        *v /= 1.0 - fstop;
    }
    if (1.0 - pi.iter().sum::<f64>()).abs() > 1e-6 {
        return Err("codon frequencies do not sum to 1".into());
    }
    Ok(pi)
}

struct F84 {
    kappa: f64,
    t: f64,
    se: f64,
}

/// F84 distance from transition/transversion proportions (`DistanceF84`).
/// Falls back F84 -> K80 -> JC69 exactly as the C code does.
fn distance_f84(n: f64, mut p: f64, q: f64, pi: &[f64; 4], want_se: bool) -> F84 {
    let (mut fail_f84, mut fail_k80, mut fail_jc69) = (false, false, false);
    let (mut a, mut b) = (0.0f64, 0.0f64);
    let qsmall = (1e-10f64).min(0.1 / n);
    let maxkappa = 999.0;
    let maxt = 99.0;
    let mut khky = -1.0f64;
    let mut t = 0.0f64;
    let mut se = 0.0f64;

    let y = pi[0] + pi[1];
    let r = pi[2] + pi[3];
    let tc = pi[0] * pi[1];
    let ag = pi[2] * pi[3];

    if p + q > 1.0 {
        return F84 { kappa: 1.0, t: maxt, se: 0.0 };
    }
    let (aa, bb, cc);
    if q < qsmall {
        fail_f84 = true;
        fail_k80 = true;
        aa = 0.0;
        bb = 0.0;
        cc = 0.0;
    } else if y <= 0.0 || r <= 0.0 || (tc <= 0.0 && ag <= 0.0) {
        fail_f84 = true;
        aa = 0.0;
        bb = 0.0;
        cc = 0.0;
    } else {
        aa = tc / y + ag / r;
        bb = tc + ag;
        cc = y * r;
        a = (2.0 * bb + 2.0 * (tc * r / y + ag * y / r) * (1.0 - q / (2.0 * cc)) - p) / (2.0 * aa);
        b = 1.0 - q / (2.0 * cc);
        if a <= 0.0 || b <= 0.0 {
            fail_f84 = true;
        }
    }

    if !fail_f84 {
        a = -0.5 * a.ln();
        b = -0.5 * b.ln();
        if b <= 0.0 {
            fail_f84 = true;
        } else {
            let k_f84 = a / b - 1.0;
            t = 4.0 * b * (tc * (1.0 + k_f84 / y) + ag * (1.0 + k_f84 / r) + cc);
            khky = (bb + (tc / y + ag / r) * k_f84) / bb;
            if want_se {
                let sa = aa * cc / (aa * cc - cc * p / 2.0 - (aa - bb) * q / 2.0);
                let sb = aa * (aa - bb) / (aa * cc - cc * p / 2.0 - (aa - bb) * q / 2.0)
                    - (aa - bb - cc) / (cc - q / 2.0);
                se = ((sa * sa * p + sb * sb * q - (sa * p + sb * q).powi(2)) / n).sqrt();
            }
        }
    }
    if fail_f84 && !fail_k80 {
        let mut ka = 1.0 - 2.0 * p - q;
        let mut kb = 1.0 - 2.0 * q;
        if ka <= 0.0 || kb <= 0.0 {
            fail_k80 = true;
        } else {
            ka = -ka.ln();
            kb = -kb.ln();
            if kb <= 0.0 {
                fail_k80 = true;
            } else {
                khky = (0.5 * ka - 0.25 * kb) / (0.25 * kb);
                t = 0.5 * ka + 0.25 * kb;
            }
            if want_se {
                let sa = 1.0 / (1.0 - 2.0 * p - q);
                let sb = (sa + 1.0 / (1.0 - 2.0 * q)) / 2.0;
                se = ((sa * sa * p + sb * sb * q - (sa * p + sb * q).powi(2)) / n).sqrt();
            }
        }
    }
    if fail_k80 {
        p += q;
        if p >= 0.75 {
            fail_jc69 = true;
            p = 0.75 * (n - 1.0) / n;
        }
        t = -0.75 * (1.0 - p * 4.0 / 3.0).ln();
        if t > maxt {
            t = maxt;
        }
        if want_se {
            se = (9.0 * p * (1.0 - p) / n).sqrt() / (3.0 - 4.0 * p);
        }
    }
    if khky > maxkappa {
        khky = maxkappa;
    }
    let _ = (fail_jc69,);
    F84 { kappa: khky, t, se }
}

/// Transition/transversion rate ratio kappa from degenerate sites (`GetKappa`).
fn get_kappa(data: &CodonData, kappa_init: f64) -> f64 {
    let table = &GENETIC_CODES[data.code];
    let kdefault = if kappa_init > 0.0 {
        kappa_init
    } else if crate::genetic_code::ncbi_id(data.code) == 2 {
        10.0
    } else {
        2.0
    };
    let mut kappa = 0.0f64;
    let npairs = data.ns as f64 * (data.ns as f64 - 1.0) / 2.0;

    for is in 0..data.ns {
        for js in 0..is {
            let mut f = [[0.0f64; 16]; 2]; // [0]=non-degenerate, [1]=4-fold
            for h in 0..data.ls {
                let c = [data.z[is][h] as usize, data.z[js][h] as usize];
                let b = [
                    [c[0] / 16, (c[0] % 16) / 4, c[0] % 4],
                    [c[1] / 16, (c[1] % 16) / 4, c[1] % 4],
                ];
                let aa = [table[c[0]], table[c[1]]];

                // non-degenerate sites
                for pos in 0..3 {
                    let mut ndeg = 0;
                    for k in 0..2 {
                        let mut i1 = 0;
                        while i1 < 4 {
                            if i1 != b[k][pos] {
                                let idx = c[k] as isize
                                    + (i1 as isize - b[k][pos] as isize) * BY[pos] as isize;
                                if table[idx as usize] == aa[k] {
                                    break;
                                }
                            }
                            i1 += 1;
                        }
                        if i1 == 4 {
                            ndeg += 1;
                        }
                    }
                    if ndeg == 2 {
                        f[0][b[0][pos] * 4 + b[1][pos]] += 0.5;
                        f[0][b[1][pos] * 4 + b[0][pos]] += 0.5;
                    }
                }

                // 4-fold degenerate third positions
                let mut ndeg = 0;
                for k in 0..2 {
                    let i1 = c[k] - b[k][2]; // codon with 3rd base 0
                    let mut j = 0;
                    while j < 4 {
                        if j != b[k][2] && table[i1 + j] != aa[k] {
                            break;
                        }
                        j += 1;
                    }
                    if aa[0] == aa[1] && j == 4 {
                        ndeg += 1;
                    }
                }
                if ndeg < 2 {
                    continue;
                }
                f[1][b[0][2] * 4 + b[1][2]] += 0.5;
                f[1][b[1][2] * 4 + b[0][2]] += 0.5;
            }

            let mut ka = [0.0f64; 2];
            let mut wk = [0.0f64; 2];
            for k in 0..2 {
                let s: f64 = f[k].iter().sum();
                if s <= 0.0 {
                    wk[k] = 0.0;
                    continue;
                }
                for v in f[k].iter_mut() {
                    *v /= s;
                }
                // transition pairs: (T,C) at row 0 col 1 and (A,G) at row 2 col 3
                let p = (f[k][1] + f[k][2 * 4 + 3]) * 2.0;
                let qq = 1.0 - (f[k][0] + f[k][5] + f[k][10] + f[k][15]) - p;
                let pi4 = [
                    f[k][0] + f[k][1] + f[k][2] + f[k][3],
                    f[k][4] + f[k][5] + f[k][6] + f[k][7],
                    f[k][8] + f[k][9] + f[k][10] + f[k][11],
                    f[k][12] + f[k][13] + f[k][14] + f[k][15],
                ];
                let r = distance_f84(s, p, qq, &pi4, false);
                ka[k] = r.kappa;
                wk[k] = if ka[k] > 0.0 { s } else { 0.0 };
            }
            let pair_kappa = if wk[0] + wk[1] == 0.0 {
                kdefault
            } else {
                (ka[0] * wk[0] + ka[1] * wk[1]) / (wk[0] + wk[1])
            };
            kappa += pair_kappa / npairs;
        }
    }
    kappa
}

/// Synonymous / nonsynonymous site counts for one sequence (`CountSites`).
fn count_sites(z: &[u8], pi: &[f64; NCODE], kappa: f64, ls: usize, code: usize) -> (f64, f64, [f64; 4], [f64; 4]) {
    let table = &GENETIC_CODES[code];
    let (mut stot, mut ntot) = (0.0, 0.0);
    let mut fbs = [0.0f64; 4];
    let mut fbn = [0.0f64; 4];
    for &codon in z {
        let c0 = codon as usize;
        let b = [c0 / 16, (c0 % 16) / 4, c0 % 4];
        let aa0 = table[c0];
        // stop codons were rejected during encoding; guard anyway
        if aa0 == -1 {
            continue;
        }
        let (mut s, mut n) = (0.0f64, 0.0f64);
        for j in 0..3 {
            for k in 0..4u32 {
                if k as usize == b[j] {
                    continue;
                }
                let c1 = (c0 as isize + (k as isize - b[j] as isize) * BY[j] as isize) as usize;
                let aa1 = table[c1];
                if aa1 == -1 {
                    continue;
                }
                let mut r = pi[c1];
                if k + b[j] as u32 == 1 || k + b[j] as u32 == 5 {
                    r *= kappa; // transition
                }
                if aa0 == aa1 {
                    s += r;
                    fbs[b[j]] += r;
                } else {
                    n += r;
                    fbn[b[j]] += r;
                }
            }
        }
        stot += s;
        ntot += n;
    }
    let scale = 3.0 * ls as f64 / (stot + ntot);
    stot *= scale;
    ntot *= scale;
    let rs: f64 = fbs.iter().sum();
    let rn: f64 = fbn.iter().sum();
    for v in fbs.iter_mut() {
        *v /= rs;
    }
    for v in fbn.iter_mut() {
        *v /= rn;
    }
    (stot, ntot, fbs, fbn)
}

/// Counts of synonymous/nonsynonymous transitions & transversions between two
/// sequences with equal pathway weighting (`CountDiffs`, PMat = all ones).
fn count_diffs(z1: &[u8], z2: &[u8], code: usize) -> (f64, f64, f64, f64) {
    let table = &GENETIC_CODES[code];
    let (mut sdts, mut sdtv, mut ndts, mut ndtv) = (0.0, 0.0, 0.0, 0.0);
    for h in 0..z1.len() {
        let c = [z1[h] as usize, z2[h] as usize];
        if c[0] == c[1] {
            continue;
        }
        let b = [
            [c[0] / 16, (c[0] % 16) / 4, c[0] % 4],
            [c[1] / 16, (c[1] % 16) / 4, c[1] % 4],
        ];
        let aa = [table[c[0]], table[c[1]]];
        let (mut sts, mut stv, mut nts, mut ntv) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
        let mut dmark = [-1i32; 3];
        let mut ndiff = 0usize;
        for (k, (x, y)) in b[0].iter().zip(b[1].iter()).enumerate() {
            if x != y {
                dmark[ndiff] = k as i32;
                ndiff += 1;
            }
        }
        if ndiff == 1 {
            let d = dmark[0] as usize;
            let transi = b[0][d] + b[1][d] == 1 || b[0][d] + b[1][d] == 5;
            if aa[0] == aa[1] {
                if transi {
                    sts += 1.0;
                } else {
                    stv += 1.0;
                }
            } else if transi {
                nts += 1.0;
            } else {
                ntv += 1.0;
            }
        } else {
            // ndiff == 2 or 3: enumerate pathways, equal weight (PMat = 1)
            let npath = if ndiff == 2 { 2 } else { 6 };
            let mut nstop = 0usize;
            let mut stspath = [0i32; 6];
            let mut stvpath = [0i32; 6];
            let mut ntspath = [0i32; 6];
            let mut ntvpath = [0i32; 6];
            let mut ok = [false; 6];
            for k in 0..npath {
                let mut step = [0usize; 3];
                if ndiff == 2 {
                    step[0] = dmark[k] as usize;
                    step[1] = dmark[1 - k] as usize;
                } else {
                    step[0] = k / 2;
                    step[1] = k % 2;
                    if step[0] <= step[1] {
                        step[1] += 1;
                    }
                    step[2] = 3 - step[0] - step[1];
                }
                let mut bt1 = b[0];
                let mut bt2 = b[0];
                let mut path_ok = true;
                for i1 in 0..ndiff {
                    bt2[step[i1]] = b[1][step[i1]];
                    let ct1 = bt1[0] * 16 + bt1[1] * 4 + bt1[2];
                    let ct2 = bt2[0] * 16 + bt2[1] * 4 + bt2[2];
                    let a0 = table[ct1];
                    let a1 = table[ct2];
                    if a1 == -1 {
                        nstop += 1;
                        path_ok = false;
                        break;
                    }
                    let sp = step[i1];
                    let transi = b[0][sp] + b[1][sp] == 1 || b[0][sp] + b[1][sp] == 5;
                    if a0 == a1 {
                        if transi {
                            stspath[k] += 1;
                        } else {
                            stvpath[k] += 1;
                        }
                    } else if transi {
                        ntspath[k] += 1;
                    } else {
                        ntvpath[k] += 1;
                    }
                    bt1 = bt2;
                }
                ok[k] = path_ok;
            }
            if npath == nstop {
                if ndiff == 2 {
                    nts = 0.5;
                    ntv = 1.5;
                } else {
                    nts = 0.5;
                    ntv = 2.5;
                }
            } else {
                let valid = (npath - nstop) as f64;
                for k in 0..npath {
                    if !ok[k] {
                        continue;
                    }
                    let p = 1.0 / valid; // equal weighting
                    sts += stspath[k] as f64 * p;
                    stv += stvpath[k] as f64 * p;
                    nts += ntspath[k] as f64 * p;
                    ntv += ntvpath[k] as f64 * p;
                }
            }
        }
        sdts += sts;
        sdtv += stv;
        ndts += nts;
        ndtv += ntv;
    }
    (sdts, sdtv, ndts, ndtv)
}

/// One pairwise YN00 (single-pass) estimate (`DistanceYN00`).
#[allow(clippy::too_many_arguments)]
fn distance_yn00(
    is: usize,
    js: usize,
    data: &CodonData,
    pi: &[f64; NCODE],
    kappa: f64,
) -> (f64, f64, f64, f64, f64, f64, f64, f64) {
    // average site counts over the two sequences
    let (mut s, mut n) = (0.0f64, 0.0f64);
    let mut fbs = [0.0f64; 4];
    let mut fbn = [0.0f64; 4];
    for &seq_idx in &[is, js] {
        let (st, nt, fbst, fbnt) = count_sites(&data.z[seq_idx], pi, kappa, data.ls, data.code);
        s += st / 2.0;
        n += nt / 2.0;
        for k in 0..4 {
            fbs[k] += fbst[k] / 2.0;
            fbn[k] += fbnt[k] / 2.0;
        }
    }

    let (sdts, sdtv, ndts, ndtv) = count_diffs(&data.z[is], &data.z[js], data.code);

    let syn = distance_f84(s, sdts / s, sdtv / s, &fbs, true);
    let non = distance_f84(n, ndts / n, ndtv / n, &fbn, true);
    let mut ds = syn.t;
    let mut dn = non.t;
    if dn < ACCURACY {
        dn = 0.0;
    }
    let omega;
    if ds < ACCURACY {
        omega = 99.0;
        ds = 0.0;
    } else {
        omega = (1e-5f64).max(dn / ds);
    }
    let t = ds * 3.0 * s / (s + n) + dn * 3.0 * n / (s + n);
    (s, n, ds, dn, syn.se, non.se, t, omega)
}

/// Nei & Gojobori (1986) distance for one pair (`DistanceMatNG86` + `difcodonNG`).
fn ng86_pair(z1: &[u8], z2: &[u8], ls: usize, code: usize) -> (f64, f64, f64) {
    let (mut nst, mut nat, mut st_tot, mut nt_tot) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
    for h in 0..z1.len() {
        let (s, n, ns, na) = difcodon_ng(z1[h] as usize, z2[h] as usize, code);
        st_tot += s;
        nt_tot += n;
        nst += ns;
        nat += na;
    }
    let y = 3.0 * ls as f64 / (st_tot + nt_tot);
    st_tot *= y;
    nt_tot *= y;
    let mut ds = 1.0 - 4.0 / 3.0 * nst / st_tot;
    let mut dn = 1.0 - 4.0 / 3.0 * nat / nt_tot;
    ds = if ds == 1.0 {
        0.0
    } else if ds < 0.0 {
        -1.0
    } else {
        3.0 / 4.0 * (-ds.ln())
    };
    dn = if dn == 1.0 {
        0.0
    } else if dn < 0.0 {
        -1.0
    } else {
        -3.0 / 4.0 * dn.ln()
    };
    let omega = if ds > 0.0 { dn / ds } else { -1.0 };
    (dn, ds, omega)
}

/// Sites and differences for a codon pair under NG86 (`difcodonNG`, transfed=1).
fn difcodon_ng(iy0: usize, iy1: usize, code: usize) -> (f64, f64, f64, f64) {
    let table = &GENETIC_CODES[code];
    let iy = [iy0, iy1];
    let b = [
        [iy0 / 16, (iy0 % 16) / 4, iy0 % 4],
        [iy1 / 16, (iy1 % 16) / 4, iy1 % 4],
    ];
    let iaa = [table[iy0], table[iy1]];

    let mut syn_site = 0.0f64;
    let mut nstop = 0usize;
    for i in 0..2 {
        for j in 0..3 {
            for k in 0..4u32 {
                if k as usize == b[i][j] {
                    continue;
                }
                let idx = iy[i] as isize + (k as isize - b[i][j] as isize) * BY[j] as isize;
                let i1 = table[idx as usize];
                if i1 == -1 {
                    nstop += 1;
                } else if i1 == iaa[i] {
                    syn_site += 1.0;
                }
            }
        }
    }
    syn_site *= 3.0 / 18.0;
    let asyn_site = 3.0 * (1.0 - nstop as f64 / 18.0) - syn_site;

    let mut dmark = [-1i32; 3];
    let mut ndiff = 0usize;
    for (k, (x, y)) in b[0].iter().zip(b[1].iter()).enumerate() {
        if x != y {
            dmark[ndiff] = k as i32;
            ndiff += 1;
        }
    }
    let (mut syn_dif, mut asyn_dif) = (0.0f64, 0.0f64);
    if ndiff == 0 {
        return (syn_site, asyn_site, 0.0, 0.0);
    }
    let npath = if ndiff == 1 { 1 } else if ndiff == 2 { 2 } else { 6 };
    let mut nstop2 = 0usize;
    if ndiff == 1 {
        if iaa[0] == iaa[1] {
            syn_dif += 1.0;
        } else {
            asyn_dif += 1.0;
        }
    } else {
        for k in 0..npath {
            let mut step = [0usize; 3];
            if ndiff == 2 {
                step[0] = dmark[k] as usize;
                step[1] = dmark[1 - k] as usize;
            } else {
                step[0] = k / 2;
                step[1] = k % 2;
                if step[0] <= step[1] {
                    step[1] += 1;
                }
                step[2] = 3 - step[0] - step[1];
            }
            let mut bt1 = b[0];
            let mut bt2 = b[0];
            let (mut sdpath, mut ndpath) = (0i32, 0i32);
            let mut broke = false;
            for i1 in 0..ndiff {
                bt2[step[i1]] = b[1][step[i1]];
                let ic0 = bt1[0] * 16 + bt1[1] * 4 + bt1[2];
                let ic1 = bt2[0] * 16 + bt2[1] * 4 + bt2[2];
                let a0 = table[ic0];
                let a1 = table[ic1];
                if a1 == -1 {
                    nstop2 += 1;
                    sdpath = 0;
                    ndpath = 0;
                    broke = true;
                    break;
                }
                if a0 == a1 {
                    sdpath += 1;
                } else {
                    ndpath += 1;
                }
                bt1 = bt2;
            }
            let _ = broke;
            syn_dif += sdpath as f64;
            asyn_dif += ndpath as f64;
        }
    }
    if npath == nstop2 {
        if ndiff == 2 {
            syn_dif = 0.0;
            asyn_dif = 2.0;
        } else {
            syn_dif = 1.0;
            asyn_dif = 2.0;
        }
    } else {
        syn_dif /= (npath - nstop2) as f64;
        asyn_dif /= (npath - nstop2) as f64;
    }
    (syn_site, asyn_site, syn_dif, asyn_dif)
}

/// Run all pairwise comparisons.
pub fn run(data: &CodonData) -> Vec<PairResult> {
    let f3 = f3x4s(data);
    let pi = get_freqs(data, &f3).unwrap_or([0.0; NCODE]);
    let kappa = get_kappa(data, 4.6);

    let mut out = Vec::new();
    for i in 0..data.ns {
        for j in (i + 1)..data.ns {
            let (s, n, ds, dn, se_ds, se_dn, t, omega) = distance_yn00(i, j, data, &pi, kappa);
            let (ng_dn, ng_ds, ng_omega) = ng86_pair(&data.z[i], &data.z[j], data.ls, data.code);
            out.push(PairResult {
                seq1: data.names[i].clone(),
                seq2: data.names[j].clone(),
                dn,
                ds,
                omega,
                n_sites: n,
                s_sites: s,
                kappa,
                t,
                se_dn,
                se_ds,
                ng_dn,
                ng_ds,
                ng_omega,
                length: data.alen,
            });
        }
    }
    out
}

// silence unused import warning if AAS is only referenced conditionally
#[allow(dead_code)]
const _AAS_KEEPALIVE: &[u8; 20] = AAS;
