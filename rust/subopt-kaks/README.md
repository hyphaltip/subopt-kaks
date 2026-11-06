# kaks — Rust port of subopt-kaks

A single self-contained binary that estimates dN/dS (Ka/Ks) from coding
sequences. It reimplements the numerical core of the `yn00_cds_optimal` and
`yn00_cds_prealigned` C tools with command-line flags only, no control files,
and a choice of output formats, and adds the pairwise maximum-likelihood
estimate of PAML `codeml` (`runmode = -2`).

Estimators:

| `--method` | What it is | PAML equivalent |
|---|---|---|
| `yn00` (default) | Yang & Nielsen 2000 counting method, single pass, equal path weighting | `yn00` with `weighting = 0` |
| `ng86` | Nei & Gojobori 1986 | `yn00` output section A |
| `ml` | Goldman & Yang 1994 codon model, t / kappa / omega by maximum likelihood, S.E. by the curvature of the likelihood | `codeml` with `runmode = -2`, `model = 0`, `NSsites = 0` |
| `both` | yn00 + ng86 | |
| `all` | yn00 + ng86 + ml | |

Only pairwise estimation is implemented. There is no tree-based likelihood,
no site or branch models, and no NSsites tests; for those use `codeml`.

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

# Maximum likelihood (codeml pairwise), F3x4 codon frequencies
kaks --method ml --prealigned input.cds.aln
kaks --method all --codon-freq fcodon --prealigned input.cds.aln

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
| `--method {yn00,ng86,ml,both,all}` | `yn00` | Estimator(s) to report (see table above) |
| `--codon-freq {fequal,f1x4,f3x4,fcodon}` | `f3x4` | Equilibrium codon frequencies for `--method ml` (codeml `CodonFreq` 0-3) |
| `--common-f3x4` | off | yn00: pool F3x4 frequencies over all sequences (PAML `commonf3x4 = 1`, the C tools' behaviour) instead of per pair |
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
alignment length in nucleotides. With `--method ml` the same columns hold the
ML values plus `lnL SE_dN SE_dS`; with `--method all` they are prefixed
`ML_`. JSON output also carries `ML_SE_t`, `ML_SE_kappa`, `ML_SE_omega`,
`ML_converged` and `ML_nfun` (likelihood evaluations).

Example (bundled test pair):

```
SEQ1     SEQ2      dN        dS        OMEGA     N           S           kappa     t         LENGTH
F35G8.1  CBG16160  0.066778  1.472123  0.045362  832.180782  289.819218  1.564887  1.289361  1227
```

### Genetic codes

`--genetic-code` takes the NCBI translation-table id used in GenBank
`transl_table` qualifiers. All 27 tables in NCBI `gc.prt` version 4.6 are
embedded (1-6, 9-16, 21-33); `kaks --list-codes` prints them. The default is 1
(standard). Note for users of the C tools: PAML numbers its codes 0-10
(`icode`), which map to NCBI ids 1, 2, 3, 4, 5, 6, 9, 10, 12, 13, 15; current
PAML 4.10 tables are identical to NCBI's. The 2004 PAML copy bundled in `src/`
keeps TGA as a stop in its yeast-mitochondrial table, whereas NCBI table 3
translates TGA as Trp. `kaks` follows NCBI.

## Validation against PAML

`../../tst/paml_compare/compare_paml.py` ran PAML 4.10.10 `codeml` (pairwise,
`runmode = -2`, all four `CodonFreq` models) and `yn00` next to `kaks` on
PAML's own example alignments: 2 to 25 sequences, 273 to 10,707 nt,
interleaved and sequential files, identity dots and ambiguity codes, under 11
genetic codes including vertebrate mitochondrial on mitochondrial data.

| | |
|---|---|
| dataset x code combinations | 17 |
| yn00 pairs compared | 794, 0 differences beyond printed precision |
| codeml ML pair x CodonFreq comparisons | 3176, 3 differences |

* lnL agrees to all six decimals codeml prints; t, kappa, omega, dN, dS, S
  and N agree to codeml's five decimals on every pair except three.
* Of those three, one is a codeml local optimum (its lnL is 1.04 units worse
  than the `kaks` value; codeml starts from random values, `kaks` from the
  NG86 distance). Two are saturated pairs with dS above 25 where the
  likelihood is flat in t and the two programs differ by under 0.5 %.
* Standard errors agree within 15 % on 95 % of pairs. Both programs use a
  finite-difference Hessian; the `kaks` values are stable across step sizes,
  codeml's change between runs. Treat S.E. values from either program as
  approximate.
* All 11 PAML genetic-code tables are identical to the NCBI tables `kaks`
  embeds.

Conclusion: for pairwise estimation `kaks --method ml` gives the same answer
as `codeml`, and where they differ `kaks` has the higher likelihood. For more
than two sequences both programs produce independent pairwise estimates; the
tree-based joint estimate of `codeml` is not implemented. Full report and
method: `../../tst/paml_compare/README.md`.

The **prealigned** path reproduces the C `yn00_cds_prealigned` output
byte-for-byte (dN, dS, ω, N, S, κ, t and alignment length) on the bundled
`tst/F35G8.1_CBG16160.cds.aln` pair. The estimator is a faithful
reimplementation of YN00 (Yang & Nielsen 2000) and NG86 (Nei & Gojobori 1986)
as used by the C tool (single pass, equal path weighting; `com.iteration = 0`).

The **from-CDS** path uses a Gotoh affine-gap protein aligner, not the Zuker
suboptimal aligner of the C `*zuker*` tools. On the bundled pair it gives the
same alignment and the same estimates as `yn00_cds_optimal`. The Zuker
suboptimal-alignment tools exist only in C.
