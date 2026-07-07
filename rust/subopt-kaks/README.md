# kaks — Rust port of subopt-kaks

A single self-contained binary that estimates dN/dS (Ka/Ks) from coding-sequence
alignments, reimplementing the numerical core of the `yn00_cds_*` C tools with a
modern command-line interface: **no config files, command-line flags only, and a
choice of output formats.**

## Why a port

The original C tools depend on PAML-style control files and, for the
alignment-from-CDS path, external substitution-matrix files loaded from a hardcoded
relative path (`../matricies/aa/BLOSUM50`, which does not even ship with the repo).
This port:

* takes **all** configuration from command-line flags,
* embeds its scoring matrix (BLOSUM62) and genetic-code tables — nothing to install,
* reads FASTA or PHYLIP and writes TSV, an aligned table, or JSON,
* has **zero** required third-party runtime dependencies (only `clap` for CLI parsing).

## Build

```sh
cargo build --release
# binary: target/release/kaks
```

## Usage

```sh
# Unaligned CDS: translate → align proteins → back-translate → estimate (default)
kaks input.cds.fasta

# Already an in-frame codon alignment: skip the alignment step
kaks --prealigned input.cds.aln

# Nei-Gojobori 1986 as well as Yang-Nielsen 2000
kaks --method both input.cds.fasta

# JSON output, read from stdin
cat input.cds.fasta | kaks --format json -
```

### Key options

| Flag | Default | Meaning |
|------|---------|---------|
| `--prealigned` | off | Treat input as an in-frame codon alignment |
| `--informat {fasta,phylip}` | auto | Input format (auto-detected otherwise) |
| `--format {tsv,table,json}` | `tsv` | Output format |
| `--method {yn00,ng86,both}` | `yn00` | Estimator(s) to report |
| `--genetic-code N` | `0` | NCBI code: 0=universal, 1=vert-mito, 2=yeast-mito, 3=mold-mito, 4=invert-mito, 5=ciliate |
| `--aln-mode {global,local}` | `global` | Aligner mode for the from-CDS path |
| `--gap-open` / `--gap-ext` | `11` / `1` | Affine gap penalties (protein alignment) |
| `--no-header` | off | Suppress the header row |
| `-o, --output FILE` | stdout | Write results to a file |

Every unordered pair of input sequences is estimated.

## Validation

The **prealigned** path reproduces the C `yn00_cds_prealigned` output byte-for-byte
(dN, dS, ω, N, S, κ, t and alignment length all match to six decimals) on the bundled
`tst/F35G8.1_CBG16160.cds.aln` test pair. This is the scientifically critical path —
the estimator is a faithful reimplementation of the YN00 (Yang & Nielsen 2000) and
NG86 (Nei & Gojobori 1986) methods as used by the C tool (single-pass, equal path
weighting; `com.iteration = 0`).

The **from-CDS** path uses a standard Gotoh affine-gap protein aligner (BLOSUM62)
rather than the C Zuker suboptimal aligner, so the alignment itself may differ; for
well-conserved pairs it reproduces the reference alignment and estimate exactly.
