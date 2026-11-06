//! Pairwise maximum-likelihood dN/dS under the Goldman & Yang (1994) codon
//! model, as computed by PAML `codeml` with `runmode = -2` (`PairwiseCodon()`
//! in codeml.c). Three parameters are estimated per pair: the distance `t`
//! (expected substitutions per codon), the transition/transversion ratio
//! `kappa`, and `omega` = dN/dS. S, N, dS and dN are then derived from the
//! rate matrix at the MLE exactly as `eigenQcodon(mode = 2)` does.
//!
//! Numerics: the reversible rate matrix is symmetrized with sqrt(pi) and
//! diagonalized with cyclic Jacobi rotations (PAML uses Householder + QL; both
//! are exact to rounding). The optimizer is the `ming2` port in `optim.rs`.

use crate::genetic_code::{ncbi_id, GENETIC_CODES};
use crate::optim::{hessian, invert, ming2, SMALL_DIFF};

/// Equilibrium codon-frequency model (PAML `CodonFreq`).
#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub enum CodonFreq {
    /// 1/n for every sense codon (CodonFreq = 0)
    Fequal,
    /// product of the overall base frequencies (CodonFreq = 1)
    F1x4,
    /// product of the position-specific base frequencies (CodonFreq = 2)
    F3x4,
    /// observed codon frequencies (CodonFreq = 3)
    Fcodon,
}

/// Result of one pairwise ML estimate.
#[derive(Clone, Debug)]
pub struct MlResult {
    pub t: f64,
    pub kappa: f64,
    pub omega: f64,
    /// log-likelihood at the MLE
    pub lnl: f64,
    pub s_sites: f64,
    pub n_sites: f64,
    pub ds: f64,
    pub dn: f64,
    pub se_t: f64,
    pub se_kappa: f64,
    pub se_omega: f64,
    pub se_dn: f64,
    pub se_ds: f64,
    pub converged: bool,
    /// likelihood evaluations used by the optimizer
    pub nfun: usize,
}

/// One single-nucleotide neighbour pair in the sense-codon index space.
struct Edge {
    i: usize,
    j: usize,
    transition: bool,
    synonymous: bool,
}

/// The codon model for one genetic code: sense codons, their amino acids and
/// the neighbour graph. Independent of the data.
struct CodonModel {
    n: usize,
    /// sense index -> codon 0..=63
    codon64: Vec<usize>,
    /// codon 0..=63 -> sense index (or usize::MAX for a stop)
    sense: [usize; 64],
    edges: Vec<Edge>,
    code: usize,
}

impl CodonModel {
    fn new(code: usize) -> Self {
        let table = &GENETIC_CODES[code];
        let mut codon64 = Vec::new();
        let mut sense = [usize::MAX; 64];
        for c in 0..64 {
            if table[c] != -1 {
                sense[c] = codon64.len();
                codon64.push(c);
            }
        }
        let n = codon64.len();
        let mut edges = Vec::new();
        for (i, &c1) in codon64.iter().enumerate().skip(1) {
            let from = [c1 / 16, (c1 / 4) % 4, c1 % 4];
            for (j, &c2) in codon64.iter().enumerate().take(i) {
                let to = [c2 / 16, (c2 / 4) % 4, c2 % 4];
                let mut ndiff = 0;
                let mut pos = 0;
                for k in 0..3 {
                    if from[k] != to[k] {
                        ndiff += 1;
                        pos = k;
                    }
                }
                if ndiff != 1 {
                    continue;
                }
                // T=0,C=1,A=2,G=3: TC (sum 1) and AG (sum 5) are transitions
                let sum = from[pos] + to[pos];
                edges.push(Edge {
                    i,
                    j,
                    transition: sum == 1 || sum == 5,
                    synonymous: table[c1] == table[c2],
                });
            }
        }
        CodonModel { n, codon64, sense, edges, code }
    }

    /// Expected codon frequencies for one pair from the observed frequencies
    /// `obs` (PAML `GetCodonFreqs2`).
    fn expected_pi(&self, obs: &[f64], cf: CodonFreq) -> Vec<f64> {
        let n = self.n;
        match cf {
            CodonFreq::Fequal => vec![1.0 / n as f64; n],
            CodonFreq::Fcodon => obs.to_vec(),
            CodonFreq::F1x4 | CodonFreq::F3x4 => {
                let mut fb3x4 = [0.0f64; 12];
                let mut fb4 = [0.0f64; 4];
                for (&c, &o) in self.codon64.iter().zip(obs) {
                    let b = [c / 16, (c / 4) % 4, c % 4];
                    for (j, &bj) in b.iter().enumerate() {
                        fb3x4[j * 4 + bj] += o;
                        fb4[bj] += o / 3.0;
                    }
                }
                let mut pi = vec![0.0; n];
                for (p, &c) in pi.iter_mut().zip(&self.codon64) {
                    let b = [c / 16, (c / 4) % 4, c % 4];
                    *p = if cf == CodonFreq::F3x4 {
                        fb3x4[b[0]] * fb3x4[4 + b[1]] * fb3x4[8 + b[2]]
                    } else {
                        fb4[b[0]] * fb4[b[1]] * fb4[b[2]]
                    };
                }
                let s: f64 = pi.iter().sum();
                for v in pi.iter_mut() {
                    *v /= s;
                }
                pi
            }
        }
    }

    /// Build the unscaled rate matrix Q (row-major n*n) and return it with
    /// its mean rate `mr = -sum_i pi_i Q_ii`.
    fn rate_matrix(&self, pi: &[f64], kappa: f64, omega: f64) -> (Vec<f64>, f64) {
        let n = self.n;
        let mut q = vec![0.0; n * n];
        for e in &self.edges {
            let mut r = 1.0;
            if e.transition {
                r *= kappa;
            }
            if !e.synonymous {
                r *= omega;
            }
            q[e.i * n + e.j] = r * pi[e.j];
            q[e.j * n + e.i] = r * pi[e.i];
        }
        let mut mr = 0.0;
        for i in 0..n {
            let row: f64 = (0..n).filter(|&j| j != i).map(|j| q[i * n + j]).sum();
            q[i * n + i] = -row;
            mr += pi[i] * row;
        }
        (q, mr)
    }

    /// S, dS and dN from the rate matrix at distance `t` (PAML `eigenQcodon`,
    /// mode 2). Returns (S, dS, dN); N = 3*ls - S.
    fn statistics(&self, pi: &[f64], kappa: f64, omega: f64, t: f64, ls: usize) -> (f64, f64, f64) {
        let n = self.n;
        let (q, mr) = self.rate_matrix(pi, kappa, omega);
        let (mut rs, mut ra, mut ra0) = (0.0, 0.0, 0.0);
        for e in &self.edges {
            let (i, j) = (e.i, e.j);
            let flux = pi[i] * q[i * n + j] + pi[j] * q[j * n + i];
            if e.synonymous {
                rs += flux;
            } else {
                ra += flux;
                ra0 += flux / omega;
            }
        }
        let rs0 = rs;
        let w = rs0 + ra0;
        let rs0n = rs0 / w;
        let ra0n = ra0 / w;
        let s = rs0n * 3.0 * ls as f64;
        if t == 0.0 {
            return (s, 0.0, 0.0);
        }
        let ds = t * (rs / mr) / (3.0 * rs0n);
        let dn = t * (ra / mr) / (3.0 * ra0n);
        (s, ds, dn)
    }

    /// Eigen-decomposition of the scaled reversible rate matrix (PAML
    /// `eigenQREV` + scaling by the mean rate): returns (root, u, v) with
    /// P(t) = U diag(exp(root t)) V, root in units of substitutions per codon.
    fn eigen(&self, pi: &[f64], kappa: f64, omega: f64) -> Eigen {
        let n = self.n;
        let (q, mr) = self.rate_matrix(pi, kappa, omega);
        let smallv = 1e-100;
        let live: Vec<usize> = (0..n).filter(|&i| pi[i] > smallv).collect();
        let m = live.len();
        let sq: Vec<f64> = live.iter().map(|&i| pi[i].sqrt()).collect();
        // symmetric S = D^1/2 Q D^-1/2 on the live codons
        let mut a = vec![0.0; m * m];
        for (ii, &i) in live.iter().enumerate() {
            a[ii * m + ii] = q[i * n + i];
            for (jj, &j) in live.iter().enumerate().take(ii) {
                let val = q[i * n + j] * sq[ii] / sq[jj];
                a[ii * m + jj] = val;
                a[jj * m + ii] = val;
            }
        }
        let (d, vec) = jacobi_eigen(&mut a, m);
        let mut root = vec![0.0; n];
        let mut u = vec![0.0; n * n];
        let mut v = vec![0.0; n * n];
        for i in 0..n {
            u[i * n + i] = 1.0;
            v[i * n + i] = 1.0;
        }
        for (kk, &k) in live.iter().enumerate() {
            root[k] = d[kk] / mr;
            for (ii, &i) in live.iter().enumerate() {
                // vec[ii*m + kk] = component ii of eigenvector kk
                u[i * n + k] = vec[ii * m + kk] / sq[ii];
                v[k * n + i] = vec[ii * m + kk] * sq[ii];
            }
        }
        Eigen { root, u, v }
    }
}

struct Eigen {
    root: Vec<f64>,
    u: Vec<f64>,
    v: Vec<f64>,
}

/// Cyclic Jacobi eigen-decomposition of a symmetric matrix (row-major `m*m`,
/// destroyed). Returns eigenvalues and eigenvectors stored by column:
/// `vec[i*m + k]` is component `i` of eigenvector `k`.
fn jacobi_eigen(a: &mut [f64], m: usize) -> (Vec<f64>, Vec<f64>) {
    let mut v = vec![0.0; m * m];
    for i in 0..m {
        v[i * m + i] = 1.0;
    }
    let mut d: Vec<f64> = (0..m).map(|i| a[i * m + i]).collect();
    let mut b = d.clone();
    let mut z = vec![0.0; m];
    for sweep in 0..60 {
        let mut sm = 0.0;
        for ip in 0..m {
            for iq in ip + 1..m {
                sm += a[ip * m + iq].abs();
            }
        }
        if sm == 0.0 {
            break;
        }
        let tresh = if sweep < 3 { 0.2 * sm / (m * m) as f64 } else { 0.0 };
        for ip in 0..m {
            for iq in ip + 1..m {
                let apq = a[ip * m + iq];
                let g = 100.0 * apq.abs();
                if sweep > 3 && d[ip].abs() + g == d[ip].abs() && d[iq].abs() + g == d[iq].abs() {
                    a[ip * m + iq] = 0.0;
                } else if apq.abs() > tresh {
                    let mut h = d[iq] - d[ip];
                    let t = if h.abs() + g == h.abs() {
                        apq / h
                    } else {
                        let theta = 0.5 * h / apq;
                        let mut t = 1.0 / (theta.abs() + (1.0 + theta * theta).sqrt());
                        if theta < 0.0 {
                            t = -t;
                        }
                        t
                    };
                    let c = 1.0 / (1.0 + t * t).sqrt();
                    let s = t * c;
                    let tau = s / (1.0 + c);
                    h = t * apq;
                    z[ip] -= h;
                    z[iq] += h;
                    d[ip] -= h;
                    d[iq] += h;
                    a[ip * m + iq] = 0.0;
                    let rot = |a: &mut [f64], i: usize, j: usize, k: usize, l: usize| {
                        let g = a[i * m + j];
                        let h = a[k * m + l];
                        a[i * m + j] = g - s * (h + g * tau);
                        a[k * m + l] = h + s * (g - h * tau);
                    };
                    for j in 0..ip {
                        rot(a, j, ip, j, iq);
                    }
                    for j in ip + 1..iq {
                        rot(a, ip, j, j, iq);
                    }
                    for j in iq + 1..m {
                        rot(a, ip, j, iq, j);
                    }
                    for j in 0..m {
                        rot(&mut v, j, ip, j, iq);
                    }
                }
            }
        }
        for ip in 0..m {
            b[ip] += z[ip];
            d[ip] = b[ip];
            z[ip] = 0.0;
        }
    }
    (d, v)
}

/// Site-pattern table for one pair: (codon index i, codon index j, count),
/// symmetrized as PAML does (i >= j).
struct PairData {
    patt: Vec<(usize, usize, f64)>,
    ls: usize,
    /// observed codon frequencies over both sequences
    obs: Vec<f64>,
}

impl PairData {
    fn new(model: &CodonModel, z1: &[u8], z2: &[u8]) -> PairData {
        let n = model.n;
        let mut fp = vec![0.0f64; n * n];
        for h in 0..z1.len() {
            let a = model.sense[z1[h] as usize];
            let b = model.sense[z2[h] as usize];
            let (j, k) = if a >= b { (a, b) } else { (b, a) };
            fp[j * n + k] += 1.0;
        }
        let ls = z1.len();
        let mut patt = Vec::new();
        let mut obs = vec![0.0; n];
        for j in 0..n {
            for k in 0..=j {
                let c = fp[j * n + k];
                if c > 0.0 {
                    patt.push((j, k, c));
                    obs[j] += c / (2.0 * ls as f64);
                    obs[k] += c / (2.0 * ls as f64);
                }
            }
        }
        PairData { patt, ls, obs }
    }
}

/// Negative log-likelihood of x = [t, kappa, omega], with the eigen-solution
/// cached across calls that only change `t`.
struct LnL<'a> {
    model: &'a CodonModel,
    data: &'a PairData,
    pi: Vec<f64>,
    cache: Option<(f64, f64, Eigen)>,
}

impl LnL<'_> {
    fn eval(&mut self, x: &[f64]) -> f64 {
        let (t, kappa, omega) = (x[0], x[1], x[2]);
        let needs = match &self.cache {
            Some((k, w, _)) => *k != kappa || *w != omega,
            None => true,
        };
        if needs {
            let e = self.model.eigen(&self.pi, kappa, omega);
            self.cache = Some((kappa, omega, e));
        }
        let e = &self.cache.as_ref().unwrap().2;
        let n = self.model.n;
        let expt: Vec<f64> = e.root.iter().map(|r| (t * r).exp()).collect();
        let mut lnl = 0.0;
        for &(z0, z1, cnt) in &self.data.patt {
            let mut fh = 0.0;
            for (k, &ek) in expt.iter().enumerate() {
                fh += e.u[z0 * n + k] * ek * e.v[k * n + z1];
            }
            fh *= self.pi[z0];
            if fh <= 0.0 {
                fh = 1e-70;
            }
            lnl -= fh.ln() * cnt;
        }
        lnl
    }
}

/// Pairwise ML estimate for two codon-encoded sequences (codon indices
/// 0..=63, equal length, no gaps or stops). `t_init` seeds the distance (use
/// the NG86 distance when available, else a negative value).
pub fn pairwise_ml(
    z1: &[u8],
    z2: &[u8],
    code: usize,
    cf: CodonFreq,
    t_init: f64,
) -> Result<MlResult, String> {
    if z1.len() != z2.len() || z1.is_empty() {
        return Err("pairwise_ml: sequences must be non-empty and equal in length".into());
    }
    let model = CodonModel::new(code);
    let data = PairData::new(&model, z1, z2);
    let pi = model.expected_pi(&data.obs, cf);
    let ls = data.ls;

    // starting values and bounds as in PairwiseCodon(), without the random jitter
    let mut t0 = if t_init.is_finite() && t_init > 0.0 { t_init } else { 0.5 };
    if t0 > 3.0 {
        t0 = 2.0;
    }
    if t0 < 1e-6 {
        t0 = 0.25;
    }
    let kappa0 = if ncbi_id(model.code) == 2 { 4.5 } else { 2.0 };
    let mut x = vec![t0, kappa0, 0.3];
    let xb = [[1e-5, 50.0], [0.4, 999.0], [0.001, 99.0]];

    let mut lnl = LnL { model: &model, data: &data, pi: pi.clone(), cache: None };
    let mut fun = |x: &[f64]| lnl.eval(x);
    let res = ming2(&mut fun, &mut x, &xb, 1e-7);
    let (t, kappa, omega) = (x[0], x[1], x[2]);
    let (s, ds, dn) = model.statistics(&pi, kappa, omega, t, ls);

    // standard errors: inverse Hessian of -lnL, then the delta method for dN, dS
    let mut se = [0.0f64; 3];
    let (mut se_dn, mut se_ds) = (0.0, 0.0);
    let hm = hessian(&mut fun, &x, res.f);
    if let Some(var) = invert(&hm, 3) {
        for k in 0..3 {
            let v = var[k * 3 + k];
            se[k] = if v > 0.0 { v.sqrt() } else { 0.0 };
        }
        // Jacobian rows: dS and dN with respect to each parameter
        let mut jac = [[0.0f64; 3]; 2];
        for k in 0..3 {
            let eh = (x[k] + 1.0) * SMALL_DIFF;
            let mut y = x.clone();
            y[k] += eh;
            let (_, ds1, dn1) = model.statistics(&pi, y[1], y[2], y[0], ls);
            y[k] -= 2.0 * eh;
            let (_, ds2, dn2) = model.statistics(&pi, y[1], y[2], y[0], ls);
            jac[0][k] = (ds1 - ds2) / (2.0 * eh);
            jac[1][k] = (dn1 - dn2) / (2.0 * eh);
        }
        let mut vsn = [[0.0f64; 2]; 2];
        for a in 0..2 {
            for b in 0..2 {
                let mut acc = 0.0;
                for k in 0..3 {
                    for l in 0..3 {
                        acc += jac[a][k] * var[k * 3 + l] * jac[b][l];
                    }
                }
                vsn[a][b] = acc;
            }
        }
        se_ds = if vsn[0][0] > 0.0 { vsn[0][0].sqrt() } else { 0.0 };
        se_dn = if vsn[1][1] > 0.0 { vsn[1][1].sqrt() } else { 0.0 };
    }

    Ok(MlResult {
        t,
        kappa,
        omega,
        lnl: -res.f,
        s_sites: s,
        n_sites: 3.0 * ls as f64 - s,
        ds,
        dn,
        se_t: se[0],
        se_kappa: se[1],
        se_omega: se[2],
        se_dn,
        se_ds,
        converged: res.converged,
        nfun: res.nfun,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn jacobi_recovers_known_eigenvalues() {
        // symmetric 3x3 with eigenvalues 1, 2, 3
        let mut a = vec![2.0, -1.0, 0.0, -1.0, 2.0, -1.0, 0.0, -1.0, 2.0];
        let (mut d, v) = jacobi_eigen(&mut a, 3);
        d.sort_by(|x, y| x.partial_cmp(y).unwrap());
        let want = [2.0 - 2f64.sqrt(), 2.0, 2.0 + 2f64.sqrt()];
        for (x, y) in d.iter().zip(want.iter()) {
            assert!((x - y).abs() < 1e-12, "{x} vs {y}");
        }
        // orthonormal eigenvectors
        for k in 0..3 {
            let nrm: f64 = (0..3).map(|i| v[i * 3 + k] * v[i * 3 + k]).sum();
            assert!((nrm - 1.0).abs() < 1e-12);
        }
    }

    #[test]
    fn standard_code_has_61_sense_codons_and_transition_edges() {
        let m = CodonModel::new(crate::genetic_code::index_of(1).unwrap());
        assert_eq!(m.n, 61);
        // every codon has at most 9 single-step neighbours; 61*9/2 minus
        // edges that would go to stops
        assert!(m.edges.len() < 61 * 9 / 2);
        assert!(m.edges.iter().any(|e| e.transition && e.synonymous));
    }
}
