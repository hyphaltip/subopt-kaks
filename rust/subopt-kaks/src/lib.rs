//! `subopt_kaks`: dN/dS (Ka/Ks) estimation from coding sequences.
//!
//! The library behind the `kaks` command-line tool and the browser build.
//! [`analyze`] takes sequence text plus [`Options`] and returns one
//! [`PairResult`] per sequence pair; [`report::write`] renders them as TSV,
//! an aligned table, or JSON.

pub mod align;
pub mod estimate;
pub mod genetic_code;
pub mod matrix;
pub mod ml;
pub mod optim;
pub mod report;
pub mod seq;
#[cfg(feature = "wasm")]
pub mod wasm;

pub use estimate::PairResult;
pub use ml::CodonFreq;
pub use report::{Method, OutFormat};

/// Analysis settings shared by every front end.
pub struct Options {
    /// Input is already an in-frame codon alignment.
    pub prealigned: bool,
    /// Input format; auto-detected when `None`.
    pub informat: Option<seq::Format>,
    pub method: Method,
    /// yn00: pool F3x4 frequencies over all sequences (PAML `commonf3x4 = 1`).
    pub common_f3x4: bool,
    /// Codon-frequency model for the ML estimate.
    pub codon_freq: CodonFreq,
    /// NCBI translation-table id.
    pub genetic_code: u8,
    pub aln_mode: align::Mode,
    /// Scoring matrix for the protein aligner: `None` uses the embedded BLOSUM62.
    pub matrix: Option<matrix::ScoreMatrix>,
    pub gap_open: i32,
    pub gap_ext: i32,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            prealigned: false,
            informat: None,
            method: Method::Yn00,
            common_f3x4: false,
            codon_freq: CodonFreq::F3x4,
            genetic_code: 1,
            aln_mode: align::Mode::Global,
            matrix: None,
            gap_open: 11,
            gap_ext: 1,
        }
    }
}

/// Parse `text` (FASTA or PHYLIP) and estimate every pair.
pub fn analyze(text: &str, opts: &Options) -> Result<Vec<PairResult>, String> {
    let fmt = opts.informat.unwrap_or_else(|| seq::detect(text));
    let records = match fmt {
        seq::Format::Fasta => seq::parse_fasta(text),
        seq::Format::Phylip => seq::parse_phylip(text),
    }
    .map_err(|e| e.to_string())?;
    analyze_records(&records, opts)
}

/// Estimate every pair of already-parsed records.
pub fn analyze_records(records: &[seq::Record], opts: &Options) -> Result<Vec<PairResult>, String> {
    let code = genetic_code::index_of(opts.genetic_code).ok_or_else(|| {
        format!(
            "unknown genetic code {}; NCBI ids are: {}",
            opts.genetic_code,
            genetic_code::CODES
                .iter()
                .map(|c| c.id.to_string())
                .collect::<Vec<_>>()
                .join(",")
        )
    })?;
    let ml = match opts.method {
        Method::Ml | Method::All => Some(opts.codon_freq),
        _ => None,
    };
    if opts.prealigned {
        let names: Vec<String> = records.iter().map(|r| r.name.clone()).collect();
        let seqs: Vec<Vec<u8>> = records.iter().map(|r| r.seq.clone().into_bytes()).collect();
        let data = estimate::clean_and_encode(names, &seqs, code)?;
        return estimate::run(&data, ml, opts.common_f3x4);
    }
    let embedded;
    let mat = match &opts.matrix {
        Some(m) => m,
        None => {
            embedded = matrix::ScoreMatrix::embedded_blosum62();
            &embedded
        }
    };
    let mut results = Vec::new();
    for i in 0..records.len() {
        for j in (i + 1)..records.len() {
            let ca = align::align_cds_pair(
                &records[i].name,
                &records[i].seq,
                &records[j].name,
                &records[j].seq,
                code,
                mat,
                opts.aln_mode,
                opts.gap_open,
                opts.gap_ext,
            )?;
            let names = vec![ca.names[0].clone(), ca.names[1].clone()];
            let seqs = vec![ca.aln[0].clone().into_bytes(), ca.aln[1].clone().into_bytes()];
            let data = estimate::clean_and_encode(names, &seqs, code)?;
            results.extend(estimate::run(&data, ml, opts.common_f3x4)?);
        }
    }
    Ok(results)
}
