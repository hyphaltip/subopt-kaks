# kaks — Rust port of subopt-kaks

A single self-contained binary that estimates dN/dS (Ka/Ks) from coding
sequences. It reimplements the numerical core of the `yn00_cds_optimal` and
`yn00_cds_prealigned` C tools with command-line flags only, no control files,
and a choice of output formats.

* takes **all** configuration from command-line flags
* embeds BLOSUM62 and all 27 NCBI genetic codes (`gc.prt` 4.6), so nothing
  else needs installing (other matrices can still be loaded from a file or a
  matrix directory)
* reads FASTA or PHYLIP, from a file or stdin
* writes TSV, an aligned table, or JSON
* no runtime dependencies (only `clap` at build time)

## Build and install

Requires a Rust toolchain (`cargo`). From this directory:

```sh
make                       # cargo build --release -> target/release/kaks
make test                  # unit tests + regression tests against ../../tst
make install               # copies kaks into ~/.local/bin
make install PREFIX=/usr/local       # -> /usr/local/bin/kaks
make install BINDIR=/path/to/bin     # any bin directory
make uninstall
```

`DESTDIR` is honoured for staged installs. The same targets are reachable
from the repository root as `make rust`, `make test-rust`,
`make install-rust PREFIX=...` and `make uninstall-rust`.

Without `make`:

```sh
cargo build --release
install -m 755 target/release/kaks /path/to/bin/
```

## Usage

```
kaks [OPTIONS] <INPUT>          INPUT is a file, or - for stdin
```

Every unordered pair of sequences in the input is estimated.

```sh
# Unaligned CDS: translate -> align proteins -> back-translate -> estimate
kaks input.cds.fasta

# Already an in-frame codon alignment: skip the alignment step
kaks --prealigned input.cds.aln

# Nei-Gojobori 1986 as well as Yang-Nielsen 2000, as a readable table
kaks --method both --format table input.cds.fasta

# JSON, reading from stdin, written to a file
cat input.cds.fasta | kaks --format json -o result.json -

# Local protein alignment with a different matrix and gap penalties
kaks --aln-mode local --matrix BLOSUM45 --gap-open 10 --gap-ext 2 input.cds.fasta

# Vertebrate mitochondrial code (NCBI transl_table 2)
kaks --genetic-code 2 input.cds.fasta

# All supported genetic codes
kaks --list-codes
```

### Options

| Flag | Default | Meaning |
|------|---------|---------|
| `--prealigned` | off | Treat input as an in-frame codon alignment |
| `--informat {fasta,phylip}` | auto | Input format (auto-detected otherwise) |
| `--format {tsv,table,json}` | `tsv` | Output format |
| `--method {yn00,ng86,both}` | `yn00` | Estimator(s) to report |
| `--genetic-code ID` | `1` | NCBI translation table id (`transl_table`): any of the 27 NCBI tables, e.g. 1 standard, 2 vertebrate mito, 3 yeast mito, 4 mold/protozoan mito, 5 invertebrate mito, 6 ciliate, 11 bacterial/plastid |
| `--list-codes` | | Print all supported genetic codes (id and name) and exit |
| `--aln-mode {global,local}` | `global` | Protein aligner mode for the from-CDS path |
| `--matrix NAME\|FILE` | built-in BLOSUM62 | Scoring matrix: a path, or a bare name (e.g. `PAM250`) |
| `--matrix-dir DIR` | | Directory searched for bare matrix names (else `$SUBOPT_MATRIX_DIR`, `$SUBOPTDIR/matricies/aa`, `./matricies/aa`) |
| `--gap-open` / `--gap-ext` | `11` / `1` | Affine gap penalties (protein alignment) |
| `--no-header` | off | Suppress the header row |
| `-o, --output FILE` | stdout | Write results to a file |

### Input requirements

* CDS sequences with introns removed, in frame, starting at the first codon.
* A trailing stop codon is removed automatically on the from-CDS path.
* With `--prealigned`, gaps must be in codon triplets. Codons containing a
  non-ACGT character are skipped.

### Output columns

`SEQ1 SEQ2 dN dS OMEGA N S kappa t LENGTH`, plus `SE_dN SE_dS` and the NG86
columns `NG_dN NG_dS NG_omega` when `--method both`. `LENGTH` is the codon
alignment length in nucleotides.

Example (bundled test pair):

```
SEQ1     SEQ2      dN        dS        OMEGA     N           S           kappa     t         LENGTH
F35G8.1  CBG16160  0.066778  1.472123  0.045362  832.180782  289.819218  1.564887  1.289361  1227
```

### Genetic codes

`--genetic-code` takes the NCBI translation-table id used in GenBank
`transl_table` qualifiers. All 27 tables in NCBI `gc.prt` version 4.6 are
embedded (1-6, 9-16, 21-33); `kaks --list-codes` prints them. The default is 1
(standard). Note for users of the C tools or PAML: PAML numbers its six codes
0-5, and its yeast-mitochondrial table keeps TGA as a stop, whereas NCBI
table 3 translates TGA as Trp. `kaks` follows NCBI.

## Validation

The **prealigned** path reproduces the C `yn00_cds_prealigned` output
byte-for-byte (dN, dS, ω, N, S, κ, t and alignment length) on the bundled
`tst/F35G8.1_CBG16160.cds.aln` pair. The estimator is a faithful
reimplementation of YN00 (Yang & Nielsen 2000) and NG86 (Nei & Gojobori 1986)
as used by the C tool (single pass, equal path weighting; `com.iteration = 0`).

The **from-CDS** path uses a Gotoh affine-gap protein aligner, not the Zuker
suboptimal aligner of the C `*zuker*` tools. On the bundled pair it gives the
same alignment and the same estimates as `yn00_cds_optimal`. The Zuker
suboptimal-alignment tools exist only in C.
