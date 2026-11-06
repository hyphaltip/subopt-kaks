//! Rendering of [`PairResult`]s as TSV, an aligned table, or JSON.

use crate::estimate::PairResult;
use std::io::Write;

/// Which estimator(s) to report.
#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub enum Method {
    /// Yang & Nielsen 2000 counting method (PAML yn00, weighting = 0)
    Yn00,
    /// Nei & Gojobori 1986
    Ng86,
    /// yn00 + ng86
    Both,
    /// Pairwise maximum likelihood, Goldman & Yang 1994 (codeml runmode = -2)
    Ml,
    /// yn00 + ng86 + ml
    All,
}

impl Method {
    pub fn parse(s: &str) -> Option<Method> {
        match s.to_ascii_lowercase().as_str() {
            "yn00" => Some(Method::Yn00),
            "ng86" => Some(Method::Ng86),
            "both" => Some(Method::Both),
            "ml" => Some(Method::Ml),
            "all" => Some(Method::All),
            _ => None,
        }
    }
}

/// Output layout.
#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub enum OutFormat {
    Tsv,
    Table,
    Json,
}

impl OutFormat {
    pub fn parse(s: &str) -> Option<OutFormat> {
        match s.to_ascii_lowercase().as_str() {
            "tsv" => Some(OutFormat::Tsv),
            "table" => Some(OutFormat::Table),
            "json" => Some(OutFormat::Json),
            _ => None,
        }
    }
}

const YN_COLS: [&str; 8] = ["dN", "dS", "OMEGA", "N", "S", "kappa", "t", "LENGTH"];
const NG_COLS: [&str; 3] = ["NG_dN", "NG_dS", "NG_OMEGA"];
const ML_COLS: [&str; 11] = [
    "ML_dN", "ML_dS", "ML_OMEGA", "ML_N", "ML_S", "ML_kappa", "ML_t", "lnL", "SE_dN", "SE_dS",
    "LENGTH",
];

/// Column names for `method`.
pub fn columns(method: Method) -> Vec<&'static str> {
    let mut c = vec!["SEQ1", "SEQ2"];
    match method {
        Method::Yn00 => c.extend(YN_COLS),
        Method::Ng86 => {
            c.extend(NG_COLS);
            c.push("LENGTH");
        }
        Method::Both => {
            c.extend(&YN_COLS[..7]);
            c.extend(NG_COLS);
            c.push("LENGTH");
        }
        // alone, the ML estimate uses the plain column names
        Method::Ml => c.extend(ML_COLS.iter().map(|s| s.strip_prefix("ML_").unwrap_or(s))),
        Method::All => {
            c.extend(&YN_COLS[..7]);
            c.extend(NG_COLS);
            c.extend(ML_COLS);
        }
    }
    c
}

/// Formatted cell values for one result, in the order of [`columns`].
pub fn row_values(r: &PairResult, method: Method) -> Vec<String> {
    let mut v = vec![r.seq1.clone(), r.seq2.clone()];
    let f = |x: f64| format!("{x:.6}");
    let yn = |v: &mut Vec<String>| {
        v.extend([
            f(r.dn),
            f(r.ds),
            f(r.omega),
            f(r.n_sites),
            f(r.s_sites),
            f(r.kappa),
            f(r.t),
        ])
    };
    let ng = |v: &mut Vec<String>| v.extend([f(r.ng_dn), f(r.ng_ds), f(r.ng_omega)]);
    let mlv = |v: &mut Vec<String>| match &r.ml {
        Some(m) => v.extend([
            f(m.dn),
            f(m.ds),
            f(m.omega),
            f(m.n_sites),
            f(m.s_sites),
            f(m.kappa),
            f(m.t),
            f(m.lnl),
            f(m.se_dn),
            f(m.se_ds),
        ]),
        None => v.extend(std::iter::repeat_n("NA".to_string(), 10)),
    };
    match method {
        Method::Yn00 => yn(&mut v),
        Method::Ng86 => ng(&mut v),
        Method::Both => {
            yn(&mut v);
            ng(&mut v);
        }
        Method::Ml => mlv(&mut v),
        Method::All => {
            yn(&mut v);
            ng(&mut v);
            mlv(&mut v);
        }
    }
    v.push(r.length.to_string());
    v
}

/// Write `results` to `out` in `format`.
pub fn write(
    out: &mut dyn Write,
    results: &[PairResult],
    method: Method,
    format: OutFormat,
    header: bool,
) -> std::io::Result<()> {
    match format {
        OutFormat::Json => write_json(out, results, method),
        OutFormat::Tsv => write_delimited(out, results, method, '\t', header),
        OutFormat::Table => write_table(out, results, method),
    }
}

/// Render `results` to a `String`.
pub fn to_string(results: &[PairResult], method: Method, format: OutFormat, header: bool) -> String {
    let mut buf = Vec::new();
    write(&mut buf, results, method, format, header).expect("writing to a Vec cannot fail");
    String::from_utf8(buf).expect("report is ASCII")
}

fn write_delimited(
    out: &mut dyn Write,
    results: &[PairResult],
    method: Method,
    sep: char,
    header: bool,
) -> std::io::Result<()> {
    if header {
        writeln!(out, "{}", columns(method).join(&sep.to_string()))?;
    }
    for r in results {
        writeln!(out, "{}", row_values(r, method).join(&sep.to_string()))?;
    }
    Ok(())
}

fn write_table(out: &mut dyn Write, results: &[PairResult], method: Method) -> std::io::Result<()> {
    let cols = columns(method);
    let mut widths: Vec<usize> = cols.iter().map(|c| c.len()).collect();
    let rows: Vec<Vec<String>> = results.iter().map(|r| row_values(r, method)).collect();
    for row in &rows {
        for (k, cell) in row.iter().enumerate() {
            widths[k] = widths[k].max(cell.len());
        }
    }
    let fmt_row = |cells: &[String]| {
        cells
            .iter()
            .enumerate()
            .map(|(k, c)| format!("{:>w$}", c, w = widths[k]))
            .collect::<Vec<_>>()
            .join("  ")
    };
    let header: Vec<String> = cols.iter().map(|s| s.to_string()).collect();
    writeln!(out, "{}", fmt_row(&header))?;
    for row in &rows {
        writeln!(out, "{}", fmt_row(row))?;
    }
    Ok(())
}

fn write_json(out: &mut dyn Write, results: &[PairResult], method: Method) -> std::io::Result<()> {
    writeln!(out, "[")?;
    for (i, r) in results.iter().enumerate() {
        let comma = if i + 1 < results.len() { "," } else { "" };
        let mut fields = format!(
            "\"seq1\":{:?},\"seq2\":{:?},\"length\":{}",
            r.seq1, r.seq2, r.length
        );
        if matches!(method, Method::Yn00 | Method::Both | Method::All) {
            fields += &format!(
                ",\"dN\":{:.6},\"dS\":{:.6},\"omega\":{:.6},\"N\":{:.6},\"S\":{:.6},\"kappa\":{:.6},\"t\":{:.6},\"SE_dN\":{:.6},\"SE_dS\":{:.6}",
                r.dn, r.ds, r.omega, r.n_sites, r.s_sites, r.kappa, r.t, r.se_dn, r.se_ds
            );
        }
        if matches!(method, Method::Ng86 | Method::Both | Method::All) {
            fields += &format!(
                ",\"NG_dN\":{:.6},\"NG_dS\":{:.6},\"NG_omega\":{:.6}",
                r.ng_dn, r.ng_ds, r.ng_omega
            );
        }
        if matches!(method, Method::Ml | Method::All) {
            if let Some(m) = &r.ml {
                fields += &format!(
                    ",\"ML_dN\":{:.6},\"ML_dS\":{:.6},\"ML_omega\":{:.6},\"ML_N\":{:.6},\"ML_S\":{:.6},\"ML_kappa\":{:.6},\"ML_t\":{:.6},\"lnL\":{:.6},\"ML_SE_dN\":{:.6},\"ML_SE_dS\":{:.6},\"ML_SE_t\":{:.6},\"ML_SE_kappa\":{:.6},\"ML_SE_omega\":{:.6},\"ML_converged\":{},\"ML_nfun\":{}",
                    m.dn, m.ds, m.omega, m.n_sites, m.s_sites, m.kappa, m.t, m.lnl, m.se_dn, m.se_ds, m.se_t, m.se_kappa, m.se_omega, m.converged, m.nfun
                );
            }
        }
        writeln!(out, "  {{{fields}}}{comma}")?;
    }
    writeln!(out, "]")?;
    Ok(())
}
