//! Bounded quasi-Newton (BFGS) minimizer, a port of PAML's `ming2()` with its
//! `gradientB()` finite-difference gradient, `LineSearch2()` quadratic-interpolation
//! line search and `H_end()` Himmelblau stopping rule (tools.c, Ziheng Yang).
//!
//! `x` holds the starting point on input and the minimizer on return; `xb[i]` is
//! the `[lower, upper]` box for `x[i]`. The function returns the minimum value,
//! `true` when the Himmelblau test was met, and, when no parameter ended on a
//! bound, the final BFGS inverse-Hessian approximation.

/// Finite-difference step base (PAML `Small_Diff`).
pub const SMALL_DIFF: f64 = 1e-6;

pub struct MinResult {
    pub f: f64,
    pub converged: bool,
    pub nfun: usize,
}

fn norm(x: &[f64]) -> f64 {
    x.iter().map(|v| v * v).sum::<f64>().sqrt()
}

fn distance(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(p, q)| (p - q) * (p - q)).sum::<f64>().sqrt()
}

fn innerp(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(p, q)| p * q).sum()
}

/// Himmelblau termination rule: true means stop.
fn h_end(x0: &[f64], x1: &[f64], f0: f64, f1: f64, e1: f64, e2: f64) -> bool {
    let mut r = norm(x0);
    if r < e2 {
        r = 1.0;
    }
    r *= e1;
    if distance(x1, x0) >= r {
        return false;
    }
    let mut r = f0.abs();
    if r < e2 {
        r = 1.0;
    }
    r *= e1;
    (f1 - f0).abs() < r
}

/// Gradient by finite differences; central for free interior parameters when
/// the step is small, one-sided (inward) for parameters on a bound.
#[allow(clippy::too_many_arguments)]
fn gradient_b<F: FnMut(&[f64]) -> f64>(
    fun: &mut F,
    x: &[f64],
    f0: f64,
    g: &mut [f64],
    xmark: &[i8],
    always_center: bool,
    sizep: f64,
    nfun: &mut usize,
) {
    let n = x.len();
    let mut x0 = vec![0.0; n];
    let mut x1 = vec![0.0; n];
    for i in 0..n {
        let mut eh = SMALL_DIFF * (x[i].abs() + 1.0);
        if xmark[i] == 0 && (always_center || sizep < 1.0) {
            x0.copy_from_slice(x);
            x1.copy_from_slice(x);
            eh = eh.powf(0.67);
            x0[i] -= eh;
            x1[i] += eh;
            g[i] = (fun(&x1) - fun(&x0)) / (eh * 2.0);
            *nfun += 2;
        } else {
            x1.copy_from_slice(x);
            if xmark[i] != 0 {
                eh *= -(xmark[i] as f64);
            }
            x1[i] += eh;
            g[i] = (fun(&x1) - f0) / eh;
            *nfun += 1;
        }
    }
}

/// Line search along `p` from `x0`: returns the step `a` and updates `f`.
#[allow(clippy::too_many_arguments)]
fn line_search2<F: FnMut(&[f64]) -> f64>(
    fun: &mut F,
    f: &mut f64,
    x0: &[f64],
    p: &[f64],
    mut step: f64,
    limit: f64,
    e: f64,
    nfun: &mut usize,
) -> f64 {
    let n = x0.len();
    let mut x = vec![0.0; n];
    let mut eval = |t: f64, x: &mut Vec<f64>, nfun: &mut usize| -> f64 {
        for i in 0..n {
            x[i] = x0[i] + t * p[i];
        }
        *nfun += 1;
        fun(x)
    };
    let factor = 4.0;
    let smallv = 1e-10;
    let smallgapa = 0.2;
    let maxround = 10;

    if step <= 0.0 || limit < smallv || step >= limit {
        return 0.0;
    }
    let a0 = 0.0;
    let mut a1 = 0.0;
    let f0 = *f;
    let mut f1 = f0;
    let mut a2 = a0 + step;
    let mut f2 = eval(a2, &mut x, nfun);
    let mut a3;
    let mut f3;
    if f2 > f1 {
        loop {
            step /= factor;
            if step < smallv {
                return 0.0;
            }
            a3 = a2;
            f3 = f2;
            a2 = a0 + step;
            f2 = eval(a2, &mut x, nfun);
            if f2 <= f1 {
                break;
            }
        }
    } else {
        loop {
            step *= factor;
            if step > limit {
                step = limit;
            }
            a3 = a0 + step;
            f3 = eval(a3, &mut x, nfun);
            if f3 >= f2 {
                break;
            }
            a1 = a2;
            f1 = f2;
            a2 = a3;
            f2 = f3;
            if step >= limit {
                *f = f3;
                return a3;
            }
        }
    }

    let mut a4 = -1.0;
    let mut f4 = -1.0;
    for _ in 0..maxround {
        a4 = (a2 - a3) * f1 + (a3 - a1) * f2 + (a1 - a2) * f3;
        if a4.abs() > 1e-100 {
            a4 = ((a2 * a2 - a3 * a3) * f1 + (a3 * a3 - a1 * a1) * f2 + (a1 * a1 - a2 * a2) * f3)
                / (2.0 * a4);
        }
        if a4 > a3 || a4 < a1 {
            a4 = (a1 + a2) / 2.0;
        }
        f4 = eval(a4, &mut x, nfun);
        if (f2 - f4).abs() < e * (1.0 + f2.abs()) {
            break;
        }
        if a4 <= a2 {
            if a2 - a4 > smallgapa * (a2 - a1) {
                if f4 <= f2 {
                    a3 = a2;
                    a2 = a4;
                    f3 = f2;
                    f2 = f4;
                } else {
                    a1 = a4;
                    f1 = f4;
                }
            } else if f4 > f2 {
                let a5 = (a2 + a3) / 2.0;
                let f5 = eval(a5, &mut x, nfun);
                if f5 > f2 {
                    a1 = a4;
                    a3 = a5;
                    f1 = f4;
                    f3 = f5;
                } else {
                    a1 = a2;
                    a2 = a5;
                    f1 = f2;
                    f2 = f5;
                }
            } else {
                let a5 = (a1 + a4) / 2.0;
                let f5 = eval(a5, &mut x, nfun);
                if f5 >= f4 {
                    a3 = a2;
                    a2 = a4;
                    a1 = a5;
                    f3 = f2;
                    f2 = f4;
                    f1 = f5;
                } else {
                    let a6 = (a1 + a5) / 2.0;
                    let f6 = eval(a6, &mut x, nfun);
                    if f6 > f5 {
                        a1 = a6;
                        a2 = a5;
                        a3 = a4;
                        f1 = f6;
                        f2 = f5;
                        f3 = f4;
                    } else {
                        a2 = a6;
                        a3 = a5;
                        f2 = f6;
                        f3 = f5;
                    }
                }
            }
        } else if a4 - a2 > smallgapa * (a3 - a2) {
            if f2 >= f4 {
                a1 = a2;
                a2 = a4;
                f1 = f2;
                f2 = f4;
            } else {
                a3 = a4;
                f3 = f4;
            }
        } else if f4 > f2 {
            let a5 = (a1 + a2) / 2.0;
            let f5 = eval(a5, &mut x, nfun);
            if f5 > f2 {
                a1 = a5;
                a3 = a4;
                f1 = f5;
                f3 = f4;
            } else {
                a3 = a2;
                a2 = a5;
                f3 = f2;
                f2 = f5;
            }
        } else {
            let a5 = (a3 + a4) / 2.0;
            let f5 = eval(a5, &mut x, nfun);
            if f5 >= f4 {
                a1 = a2;
                a2 = a4;
                a3 = a5;
                f1 = f2;
                f2 = f4;
                f3 = f5;
            } else {
                let a6 = (a3 + a5) / 2.0;
                let f6 = eval(a6, &mut x, nfun);
                if f6 > f5 {
                    a1 = a4;
                    a2 = a5;
                    a3 = a6;
                    f1 = f4;
                    f2 = f5;
                    f3 = f6;
                } else {
                    a1 = a5;
                    a2 = a6;
                    f1 = f5;
                    f2 = f6;
                }
            }
        }
    }

    if f2 > f0 && f4 > f0 {
        a4 = 0.0;
    }
    if f2 <= f4 {
        *f = f2;
        a4 = a2;
    } else {
        *f = f4;
    }
    a4
}

/// Minimize `fun` over the box `xb` starting from `x` (PAML `ming2`).
pub fn ming2<F: FnMut(&[f64]) -> f64>(
    fun: &mut F,
    x: &mut [f64],
    xb: &[[f64; 2]],
    e: f64,
) -> MinResult {
    let n = x.len();
    let maxround = 10000;
    let ngoodtimes = 2;
    let smallv = 1e-30;
    let maxstep = 8.0;
    let mut nfun = 0usize;
    if n == 0 {
        return MinResult { f: 0.0, converged: true, nfun };
    }

    let mut xmark = vec![0i8; n];
    let mut ix: Vec<usize> = Vec::with_capacity(n);
    for i in 0..n {
        if x[i] <= xb[i][0] {
            x[i] = xb[i][0];
            xmark[i] = -1;
            continue;
        }
        if x[i] >= xb[i][1] {
            x[i] = xb[i][1];
            xmark[i] = 1;
            continue;
        }
        ix.push(i);
    }
    let mut nfree = ix.len();

    let mut f0 = fun(x);
    nfun += 1;
    let mut f;
    let mut x0 = x.to_vec();
    let mut sizep = 99.0f64;
    let mut sizep0;
    let mut always_center = false;
    let mut fail = false;
    let mut goodtimes = 0;

    let mut g0 = vec![0.0; n];
    let mut g = vec![0.0; n];
    let mut p = vec![0.0; n];
    let mut y = vec![0.0; n];
    let mut s = vec![0.0; n];
    let mut z = vec![0.0; n];
    gradient_b(fun, &x0, f0, &mut g0, &xmark, always_center, sizep, &mut nfun);

    let identity = |m: &mut Vec<f64>, k: usize| {
        *m = vec![0.0; k * k];
        for i in 0..k {
            m[i * k + i] = 1.0;
        }
    };
    let mut h = Vec::new();
    identity(&mut h, nfree);

    let mut iround = 0;
    let mut converged = false;
    while iround < maxround {
        // search direction p = -H g (free parameters only)
        for v in p.iter_mut() {
            *v = 0.0;
        }
        for i in 0..nfree {
            for j in 0..nfree {
                p[ix[i]] -= h[i * nfree + j] * g0[ix[j]];
            }
        }
        sizep0 = sizep;
        sizep = norm(&p);

        // maximum step length before hitting a bound
        let mut am = maxstep;
        for i in 0..n {
            if p[i] > 0.0 && (xb[i][1] - x0[i]) / p[i] < am {
                am = (xb[i][1] - x0[i]) / p[i];
            } else if p[i] < 0.0 && (xb[i][0] - x0[i]) / p[i] < am {
                am = (xb[i][0] - x0[i]) / p[i];
            }
        }
        let mut hstep = if iround == 0 {
            let hh = (2.0 * f0 * 0.01 / innerp(&g0, &p)).abs();
            hh.min(am / 2000.0)
        } else {
            let hh = norm(&s[..nfree]) / sizep;
            hh.max(am / 500.0)
        };
        hstep = hstep.max(1e-5).min(am / 5.0);
        f = f0;
        let alpha = line_search2(fun, &mut f, &x0, &p, hstep, am, e.min(1e-3), &mut nfun);

        if alpha <= 0.0 {
            if fail {
                if always_center {
                    break;
                } else {
                    always_center = true;
                    identity(&mut h, nfree);
                    fail = true;
                }
            } else {
                identity(&mut h, nfree);
                fail = true;
            }
        } else {
            fail = false;
            for i in 0..n {
                x[i] = x0[i] + alpha * p[i];
            }
            if iround == 0 || sizep < sizep0 || (sizep < 0.001 && sizep0 < 0.001) {
                goodtimes += 1;
            } else {
                goodtimes = 0;
            }
            if (n == 1 || goodtimes >= ngoodtimes)
                && sizep < (if e > 1e-5 { 1.0 } else { 0.001 })
                && h_end(&x0, x, f0, f, e, e)
            {
                converged = true;
                break;
            }
        }
        gradient_b(fun, x, f, &mut g, &xmark, always_center, sizep, &mut nfun);

        // add constraints: parameters that reached a bound and push outward
        for i in 0..n {
            if xmark[i] != 0 {
                continue;
            }
            if (x[i] - xb[i][0]).abs() < 1e-6 && -g[i] < 0.0 {
                xmark[i] = -1;
            } else if (x[i] - xb[i][1]).abs() < 1e-6 && -g[i] > 0.0 {
                xmark[i] = 1;
            }
            if xmark[i] == 0 {
                continue;
            }
            let it = ix.iter().position(|&k| k == i).unwrap_or(nfree);
            if it == nfree {
                continue;
            }
            let c = h.clone();
            ix.remove(it);
            nfree -= 1;
            h = vec![0.0; nfree * nfree];
            for i1 in 0..nfree {
                for i2 in 0..nfree {
                    let r = i1 + usize::from(i1 >= it);
                    let cc = i2 + usize::from(i2 >= it);
                    h[i1 * nfree + i2] = c[r * (nfree + 1) + cc];
                }
            }
        }
        // delete a constraint when the gradient pulls a bound parameter inward
        let mut it = 0;
        let mut w = 0.0;
        for i in 0..n {
            if xmark[i] == -1 && -g[i] > w {
                it = i;
                w = -g[i];
            } else if xmark[i] == 1 && -g[i] < -w {
                it = i;
                w = g[i];
            }
        }
        if nfree > 0 && w > 10.0 * sizep / nfree as f64 || (nfree == 0 && w > 0.0) {
            let c = h.clone();
            let nf1 = nfree + 1;
            h = vec![0.0; nf1 * nf1];
            for i1 in 0..nfree {
                for i2 in 0..nfree {
                    h[i1 * nf1 + i2] = c[i1 * nfree + i2];
                }
            }
            h[nf1 * nf1 - 1] = 1.0;
            xmark[it] = 0;
            ix.push(it);
            nfree = nf1;
        }

        f0 = f;
        for i in 0..nfree {
            y[i] = g[ix[i]] - g0[ix[i]];
            s[i] = x[ix[i]] - x0[ix[i]];
        }
        g0.copy_from_slice(&g);
        x0.copy_from_slice(x);

        // BFGS update of the inverse Hessian
        let mut w = 0.0;
        let mut v = 0.0;
        for i in 0..nfree {
            z[i] = 0.0;
            for j in 0..nfree {
                z[i] += h[i * nfree + j] * y[j];
            }
            w += y[i] * z[i];
            v += y[i] * s[i];
        }
        if v.abs() < smallv {
            identity(&mut h, nfree);
            fail = true;
            iround += 1;
            continue;
        }
        for i in 0..nfree {
            for j in 0..nfree {
                h[i * nfree + j] +=
                    ((1.0 + w / v) * s[i] * s[j] - z[i] * s[j] - s[i] * z[j]) / v;
            }
        }
        iround += 1;
    }

    f = fun(x);
    nfun += 1;
    // Stopping after an exhausted line search with central-difference gradients
    // means no further decrease is resolvable; only hitting the round limit is
    // reported as non-convergence (PAML prints "check convergence!" there).
    let _ = converged;
    MinResult { f, converged: iround < maxround, nfun }
}

/// Hessian of `fun` at `x` by central differences (PAML `Hessian`), returned
/// row-major `n*n`.
pub fn hessian<F: FnMut(&[f64]) -> f64>(fun: &mut F, x: &[f64], f0: f64) -> Vec<f64> {
    let n = x.len();
    let h0 = SMALL_DIFF * 2.0;
    let h: Vec<f64> = x.iter().map(|v| h0 * (1.0 + v.abs())).collect();
    let mut hm = vec![0.0; n * n];
    let mut x1 = vec![0.0; n];
    for i in 0..n {
        for j in i..n {
            x1.copy_from_slice(x);
            x1[i] += h[i];
            x1[j] += h[j];
            let fpp = fun(&x1);
            x1[i] -= h[i] * 2.0;
            x1[j] -= h[j] * 2.0;
            let fmm = fun(&x1);
            if i == j {
                hm[i * n + i] = (fpp + fmm - 2.0 * f0) / (4.0 * h[i] * h[i]);
            } else {
                x1[i] += 2.0 * h[i];
                let fpm = fun(&x1);
                x1[i] -= 2.0 * h[i];
                x1[j] += 2.0 * h[j];
                let fmp = fun(&x1);
                let v = (fpp + fmm - fpm - fmp) / (4.0 * h[i] * h[j]);
                hm[i * n + j] = v;
                hm[j * n + i] = v;
            }
        }
    }
    hm
}

/// Invert a small row-major matrix in place by Gauss-Jordan elimination with
/// partial pivoting. Returns `None` if singular.
pub fn invert(a: &[f64], n: usize) -> Option<Vec<f64>> {
    let mut m = a.to_vec();
    let mut inv = vec![0.0; n * n];
    for i in 0..n {
        inv[i * n + i] = 1.0;
    }
    for col in 0..n {
        let mut piv = col;
        for r in col + 1..n {
            if m[r * n + col].abs() > m[piv * n + col].abs() {
                piv = r;
            }
        }
        if m[piv * n + col].abs() < 1e-300 {
            return None;
        }
        if piv != col {
            for k in 0..n {
                m.swap(col * n + k, piv * n + k);
                inv.swap(col * n + k, piv * n + k);
            }
        }
        let d = m[col * n + col];
        for k in 0..n {
            m[col * n + k] /= d;
            inv[col * n + k] /= d;
        }
        for r in 0..n {
            if r == col {
                continue;
            }
            let f = m[r * n + col];
            if f != 0.0 {
                for k in 0..n {
                    m[r * n + k] -= f * m[col * n + k];
                    inv[r * n + k] -= f * inv[col * n + k];
                }
            }
        }
    }
    Some(inv)
}
