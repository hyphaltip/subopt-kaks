//! `kaks` — estimate dN/dS (Ka/Ks) from coding-sequence alignments.
//!
//! A Rust port of the subopt-kaks C tools with a single command-line-driven binary
//! and choice of output format. Two input modes:
//! * `--prealigned` : input is an in-frame codon alignment (like `yn00_cds_prealigned`)
//! * default        : input is unaligned CDS; sequences are translated, aligned on the
//!   protein, and back-translated (like `yn00_cds_optimal`).

mod align;
mod estimate;
mod genetic_code;
mod matrix;
mod seq;

use clap::{Parser, ValueEnum};
use estimate::PairResult;
use std::io::Write;

#[derive(Copy, Clone, PartialEq, Eq, ValueEnum, Debug)]
enum InFormat {
    Fasta,
    Phylip,
}

#[derive(Copy, Clone, PartialEq, Eq, ValueEnum, Debug)]
enum OutFormat {
    /// Tab-separated table (default)
    Tsv,
    /// Aligned, human-readable columns
    Table,
    /// One JSON object per pair (JSON array)
    Json,
}

#[derive(Copy, Clone, PartialEq, Eq, ValueEnum, Debug)]
enum Method {
    /// Yang & Nielsen 2000 (single pass) — the tool's primary estimate
    Yn00,
    /// Nei & Gojobori 1986
    Ng86,
    /// Both methods
    Both,
}

#[derive(Copy, Clone, PartialEq, Eq, ValueEnum, Debug)]
enum AlnMode {
    Global,
    Local,
}

#[derive(Parser, Debug)]
#[command(
    name = "kaks",
    version,
    about = "Estimate dN/dS (Ka/Ks) from coding-sequence alignments",
    long_about = None
)]
struct Cli {
    /// Input sequence file (`-` for stdin)
    #[arg(required_unless_present = "list_codes")]
    input: Option<String>,

    /// Treat input as an in-frame codon alignment (skip the alignment step)
    #[arg(long)]
    prealigned: bool,

    /// Input file format (auto-detected if omitted)
    #[arg(long, value_enum)]
    informat: Option<InFormat>,

    /// Output format
    #[arg(long, value_enum, default_value_t = OutFormat::Tsv)]
    format: OutFormat,

    /// Which estimator(s) to report
    #[arg(long, value_enum, default_value_t = Method::Yn00)]
    method: Method,

    /// NCBI translation-table id (1=standard, 2=vertebrate mito, 3=yeast mito,
    /// 4=mold/protozoan mito, 5=invertebrate mito, 6=ciliate nuclear, 9..33);
    /// see --list-codes
    #[arg(long, default_value_t = 1, value_name = "ID")]
    genetic_code: u8,

    /// Print all supported genetic codes (id and name) and exit
    #[arg(long)]
    list_codes: bool,

    /// Alignment mode when starting from unaligned CDS
    #[arg(long, value_enum, default_value_t = AlnMode::Global)]
    aln_mode: AlnMode,

    /// Substitution matrix for the from-CDS aligner: a file path, or a bare name
    /// (e.g. BLOSUM62, PAM250) resolved within a matrix directory
    #[arg(long)]
    matrix: Option<String>,

    /// Directory to search for named matrices (else $SUBOPT_MATRIX_DIR, $SUBOPTDIR,
    /// ./matricies/aa, ...). BLOSUM62 falls back to a built-in copy.
    #[arg(long)]
    matrix_dir: Option<String>,

    /// Gap-open penalty (protein alignment)
    #[arg(long, default_value_t = 11)]
    gap_open: i32,

    /// Gap-extension penalty (protein alignment)
    #[arg(long, default_value_t = 1)]
    gap_ext: i32,

    /// Suppress the header line
    #[arg(long)]
    no_header: bool,

    /// Write output to FILE instead of stdout
    #[arg(long, short)]
    output: Option<String>,
}

fn main() {
    let cli = Cli::parse();
    if cli.list_codes {
        print!("{}", genetic_code::list());
        return;
    }
    let code = match genetic_code::index_of(cli.genetic_code) {
        Some(i) => i,
        None => {
            eprintln!(
                "error: unknown genetic code {}; NCBI ids are: {}",
                cli.genetic_code,
                genetic_code::CODES
                    .iter()
                    .map(|c| c.id.to_string())
                    .collect::<Vec<_>>()
                    .join(",")
            );
            std::process::exit(2);
        }
    };
    if let Err(e) = run(&cli, code) {
        eprintln!("error: {e}");
        std::process::exit(1);
    }
}

fn run(cli: &Cli, code: usize) -> Result<(), Box<dyn std::error::Error>> {
    let fmt = cli.informat.map(|f| match f {
        InFormat::Fasta => seq::Format::Fasta,
        InFormat::Phylip => seq::Format::Phylip,
    });
    let input = cli.input.as_deref().ok_or("missing INPUT")?;
    let records = seq::read_path(input, fmt)?;

    let results: Vec<PairResult> = if cli.prealigned {
        let names: Vec<String> = records.iter().map(|r| r.name.clone()).collect();
        let seqs: Vec<Vec<u8>> = records.iter().map(|r| r.seq.clone().into_bytes()).collect();
        let data = estimate::clean_and_encode(names, &seqs, code)?;
        estimate::run(&data)
    } else {
        align_and_estimate(cli, &records, code)?
    };

    let out: Box<dyn Write> = match &cli.output {
        Some(p) if p != "-" => Box::new(std::fs::File::create(p)?),
        _ => Box::new(std::io::stdout()),
    };
    write_results(out, &results, cli)?;
    Ok(())
}

/// Align every pair of CDS records and estimate.
fn align_and_estimate(
    cli: &Cli,
    records: &[seq::Record],
    code: usize,
) -> Result<Vec<PairResult>, Box<dyn std::error::Error>> {
    let mode = match cli.aln_mode {
        AlnMode::Global => align::Mode::Global,
        AlnMode::Local => align::Mode::Local,
    };
    let mat = matrix::load(cli.matrix.as_deref(), cli.matrix_dir.as_deref())?;
    let mut results = Vec::new();
    for i in 0..records.len() {
        for j in (i + 1)..records.len() {
            let ca = align::align_cds_pair(
                &records[i].name,
                &records[i].seq,
                &records[j].name,
                &records[j].seq,
                code,
                &mat,
                mode,
                cli.gap_open,
                cli.gap_ext,
            )?;
            let names = vec![ca.names[0].clone(), ca.names[1].clone()];
            let seqs = vec![ca.aln[0].clone().into_bytes(), ca.aln[1].clone().into_bytes()];
            let data = estimate::clean_and_encode(names, &seqs, code)?;
            results.extend(estimate::run(&data));
        }
    }
    Ok(results)
}

fn write_results(
    mut out: Box<dyn Write>,
    results: &[PairResult],
    cli: &Cli,
) -> std::io::Result<()> {
    match cli.format {
        OutFormat::Json => write_json(&mut out, results, cli.method),
        OutFormat::Tsv => write_delimited(&mut out, results, cli, '\t'),
        OutFormat::Table => write_table(&mut out, results, cli.method),
    }
}

fn columns(method: Method) -> Vec<&'static str> {
    let mut c = vec!["SEQ1", "SEQ2"];
    match method {
        Method::Yn00 => c.extend(["dN", "dS", "OMEGA", "N", "S", "kappa", "t", "LENGTH"]),
        Method::Ng86 => c.extend(["NG_dN", "NG_dS", "NG_OMEGA", "LENGTH"]),
        Method::Both => c.extend([
            "dN", "dS", "OMEGA", "N", "S", "kappa", "t", "NG_dN", "NG_dS", "NG_OMEGA", "LENGTH",
        ]),
    }
    c
}

fn row_values(r: &PairResult, method: Method) -> Vec<String> {
    let mut v = vec![r.seq1.clone(), r.seq2.clone()];
    let f = |x: f64| format!("{x:.6}");
    match method {
        Method::Yn00 => v.extend([
            f(r.dn),
            f(r.ds),
            f(r.omega),
            f(r.n_sites),
            f(r.s_sites),
            f(r.kappa),
            f(r.t),
            r.length.to_string(),
        ]),
        Method::Ng86 => v.extend([
            f(r.ng_dn),
            f(r.ng_ds),
            f(r.ng_omega),
            r.length.to_string(),
        ]),
        Method::Both => v.extend([
            f(r.dn),
            f(r.ds),
            f(r.omega),
            f(r.n_sites),
            f(r.s_sites),
            f(r.kappa),
            f(r.t),
            f(r.ng_dn),
            f(r.ng_ds),
            f(r.ng_omega),
            r.length.to_string(),
        ]),
    }
    v
}

fn write_delimited(
    out: &mut Box<dyn Write>,
    results: &[PairResult],
    cli: &Cli,
    sep: char,
) -> std::io::Result<()> {
    if !cli.no_header {
        writeln!(out, "{}", columns(cli.method).join(&sep.to_string()))?;
    }
    for r in results {
        writeln!(out, "{}", row_values(r, cli.method).join(&sep.to_string()))?;
    }
    Ok(())
}

fn write_table(out: &mut Box<dyn Write>, results: &[PairResult], method: Method) -> std::io::Result<()> {
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

fn write_json(out: &mut Box<dyn Write>, results: &[PairResult], method: Method) -> std::io::Result<()> {
    writeln!(out, "[")?;
    for (i, r) in results.iter().enumerate() {
        let comma = if i + 1 < results.len() { "," } else { "" };
        let mut fields = format!(
            "\"seq1\":{:?},\"seq2\":{:?},\"length\":{}",
            r.seq1, r.seq2, r.length
        );
        if matches!(method, Method::Yn00 | Method::Both) {
            fields += &format!(
                ",\"dN\":{:.6},\"dS\":{:.6},\"omega\":{:.6},\"N\":{:.6},\"S\":{:.6},\"kappa\":{:.6},\"t\":{:.6},\"SE_dN\":{:.6},\"SE_dS\":{:.6}",
                r.dn, r.ds, r.omega, r.n_sites, r.s_sites, r.kappa, r.t, r.se_dn, r.se_ds
            );
        }
        if matches!(method, Method::Ng86 | Method::Both) {
            fields += &format!(
                ",\"NG_dN\":{:.6},\"NG_dS\":{:.6},\"NG_omega\":{:.6}",
                r.ng_dn, r.ng_ds, r.ng_omega
            );
        }
        writeln!(out, "  {{{fields}}}{comma}")?;
    }
    writeln!(out, "]")?;
    Ok(())
}
