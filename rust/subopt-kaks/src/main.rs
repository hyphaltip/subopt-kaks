//! `kaks` — estimate dN/dS (Ka/Ks) from coding-sequence alignments.
//!
//! Command-line front end over the `subopt_kaks` library. Two input modes:
//! * `--prealigned` : input is an in-frame codon alignment (like `yn00_cds_prealigned`)
//! * default        : input is unaligned CDS; sequences are translated, aligned on the
//!   protein, and back-translated (like `yn00_cds_optimal`).

use clap::{Parser, ValueEnum};
use std::io::Write;
use subopt_kaks::{align, genetic_code, matrix, report, seq, CodonFreq, Method, Options, OutFormat};

#[derive(Copy, Clone, PartialEq, Eq, ValueEnum, Debug)]
enum InFormat {
    Fasta,
    Phylip,
}

#[derive(Copy, Clone, PartialEq, Eq, ValueEnum, Debug)]
enum OutFormatArg {
    /// Tab-separated table (default)
    Tsv,
    /// Aligned, human-readable columns
    Table,
    /// One JSON object per pair (JSON array)
    Json,
}

#[derive(Copy, Clone, PartialEq, Eq, ValueEnum, Debug)]
enum MethodArg {
    /// Yang & Nielsen 2000 (single pass) — the tool's primary estimate
    Yn00,
    /// Nei & Gojobori 1986
    Ng86,
    /// Both counting methods (yn00 + ng86)
    Both,
    /// Pairwise maximum likelihood under the Goldman & Yang 1994 codon model
    /// (PAML codeml runmode = -2)
    Ml,
    /// yn00, ng86 and ml
    All,
}

#[derive(Copy, Clone, PartialEq, Eq, ValueEnum, Debug)]
enum CodonFreqArg {
    /// Equal frequency for every sense codon (codeml CodonFreq = 0)
    Fequal,
    /// Product of overall base frequencies (CodonFreq = 1)
    F1x4,
    /// Product of position-specific base frequencies (CodonFreq = 2)
    F3x4,
    /// Observed codon frequencies (CodonFreq = 3)
    Fcodon,
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
    #[arg(long, value_enum, default_value_t = OutFormatArg::Tsv)]
    format: OutFormatArg,

    /// Which estimator(s) to report
    #[arg(long, value_enum, default_value_t = MethodArg::Yn00)]
    method: MethodArg,

    /// yn00: pool F3x4 codon frequencies over all sequences instead of each
    /// pair (PAML yn00 commonf3x4 = 1; the C yn00_cds_* tools' behaviour).
    /// No effect with two sequences.
    #[arg(long)]
    common_f3x4: bool,

    /// Codon-frequency model for --method ml (codeml CodonFreq)
    #[arg(long, value_enum, default_value_t = CodonFreqArg::F3x4)]
    codon_freq: CodonFreqArg,

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
    if let Err(e) = run(&cli) {
        eprintln!("error: {e}");
        std::process::exit(1);
    }
}

fn run(cli: &Cli) -> Result<(), Box<dyn std::error::Error>> {
    let method = match cli.method {
        MethodArg::Yn00 => Method::Yn00,
        MethodArg::Ng86 => Method::Ng86,
        MethodArg::Both => Method::Both,
        MethodArg::Ml => Method::Ml,
        MethodArg::All => Method::All,
    };
    let format = match cli.format {
        OutFormatArg::Tsv => OutFormat::Tsv,
        OutFormatArg::Table => OutFormat::Table,
        OutFormatArg::Json => OutFormat::Json,
    };
    // Only resolve a matrix when it will be used; a missing named matrix is
    // then reported before any sequence work.
    let matrix = if cli.prealigned {
        None
    } else {
        Some(matrix::load(cli.matrix.as_deref(), cli.matrix_dir.as_deref())?)
    };
    let opts = Options {
        prealigned: cli.prealigned,
        informat: cli.informat.map(|f| match f {
            InFormat::Fasta => seq::Format::Fasta,
            InFormat::Phylip => seq::Format::Phylip,
        }),
        method,
        common_f3x4: cli.common_f3x4,
        codon_freq: match cli.codon_freq {
            CodonFreqArg::Fequal => CodonFreq::Fequal,
            CodonFreqArg::F1x4 => CodonFreq::F1x4,
            CodonFreqArg::F3x4 => CodonFreq::F3x4,
            CodonFreqArg::Fcodon => CodonFreq::Fcodon,
        },
        genetic_code: cli.genetic_code,
        aln_mode: match cli.aln_mode {
            AlnMode::Global => align::Mode::Global,
            AlnMode::Local => align::Mode::Local,
        },
        matrix,
        gap_open: cli.gap_open,
        gap_ext: cli.gap_ext,
    };
    let input = cli.input.as_deref().ok_or("missing INPUT")?;
    let records = seq::read_path(input, opts.informat)?;
    let results = subopt_kaks::analyze_records(&records, &opts)?;

    let mut out: Box<dyn Write> = match &cli.output {
        Some(p) if p != "-" => Box::new(std::fs::File::create(p)?),
        _ => Box::new(std::io::stdout()),
    };
    report::write(&mut *out, &results, method, format, !cli.no_header)?;
    Ok(())
}
