# subopt-kaks

Suboptimal alignments and Ka/Ks (dN/dS) calculations for pairs of coding
sequences.

This source is freely available, but should be considered a personal
communication: it has not been properly packaged and does not yet give
adequate credit to the code it reuses. Please contact me if you have
questions or intend to use this for large-scale analyses so I can be sure
it is performing properly for you.

* **SQUID** is Sean Eddy's sequence library (bundled in `squid/`).
* The alignment code is based on design ideas from Michael Smoot.
* dN/dS estimation reuses Z. Yang's `yn00` implementation from PAML, with
  modifications to detach it from the PAML data structures.

## Build

Requirements: a C compiler (GCC 8 through 16 are tested), `make`, and a
POSIX shell. Optional: a Rust toolchain for the `kaks` port.

```sh
make            # configures and builds squid, then the C tools in src/
make test       # runs the regression tests in tst/
make install    # copies the tools into ./bin
```

`make install PREFIX=/usr/local` installs into `$PREFIX/bin`, puts the
substitution matrices in `$PREFIX/share/subopt-kaks/matricies`, and compiles
that data directory into the binaries, so no environment variable is needed
at run time. `DESTDIR` is honoured for staged installs. `make uninstall`
removes what `make install` added.

## Programs

| Program | Purpose |
|---------|---------|
| `yn00_cds_optimal` | Translate two or more CDS, align the proteins (local by default, `--global` for global), back-translate, and report dN/dS for every pair. Use `--showtable` for the tab-delimited summary. |
| `yn00_cds_prealigned` | Compute dN/dS for every pair in an existing in-frame codon alignment (one or more files). |
| `zuker_align` | Zuker suboptimal alignments of a pair of sequences. |
| `cds_zuker_align` | Zuker suboptimal alignments of a pair of CDS, back-translated to codons. |
| `yn00_zuker_align` | dN/dS for each Zuker suboptimal alignment of a CDS pair (writes `summary.dat` unless `--output` is given; `--output -` writes to stdout). |

Every program prints its options with `-h`.

```sh
# Ka/Ks table from unaligned CDS (stop codon removed automatically)
./bin/yn00_cds_optimal --showtable tst/F35G8.1_CBG16160.cds.fas

# Ka/Ks from a pre-aligned codon alignment
./bin/yn00_cds_prealigned tst/F35G8.1_CBG16160.cds.aln
```

Input files must contain CDS sequences with introns removed. A trailing
stop codon is stripped by `yn00_cds_optimal`; the other tools expect it to be
absent already (this is a PAML requirement).

### Substitution matrices

The default matrix is `BLOSUM62`. `--matrix` accepts either a path or a bare
name such as `BLOSUM45` or `PAM250`. A bare name is searched for, in order, in
`--matrix-dir`, `$SUBOPT_MATRIX_DIR`, `$SUBOPTDIR/matricies/aa`,
`./matricies/aa`, `../matricies/aa`, and the directory compiled in at build
time. The bundled matrices are in `matricies/aa` and `matricies/nt`.

### Gap penalties

`--gapopen` and `--gapext` take a penalty magnitude. `--gapopen 10` and
`--gapopen -10` mean the same thing (a penalty of 10). Defaults are 10 and 2.

## Rust port: `kaks` (recommended for dN/dS)

`rust/subopt-kaks` is a single-binary reimplementation of `yn00_cds_optimal`
and `yn00_cds_prealigned`: flags only, embedded BLOSUM62 and genetic-code
tables, FASTA/PHYLIP input from a file or stdin, TSV/table/JSON output, all 27
NCBI genetic codes. Estimators: YN00 (the counting method of PAML `yn00`,
default), NG86, and pairwise maximum likelihood under the Goldman-Yang codon
model (`--method ml`, equivalent to PAML `codeml` with `runmode = -2`). It is
validated against PAML's own `codeml` and `yn00` on PAML's example data
(`tst/paml_compare/`). Its `--prealigned` path reproduces
`yn00_cds_prealigned` byte-for-byte on the bundled test pair. The Zuker
suboptimal-alignment tools exist only in C.

### Build and install

Requires a Rust toolchain (`cargo`).

```sh
make rust                              # build -> rust/subopt-kaks/target/release/kaks
make test-rust                         # cargo test + regression tests
make install-rust                      # copy kaks into ./bin
make install-rust PREFIX=/usr/local    # -> /usr/local/bin/kaks
make install-rust BINDIR=$HOME/bin     # any bin directory
make uninstall-rust
```

Or from `rust/subopt-kaks/`: `make`, `make test`, `make install [PREFIX=|BINDIR=]`.
Only the one binary is installed; it needs no data files or environment variables.

WebAssembly: `make -C rust/subopt-kaks wasi` builds the same CLI for any WASI
runtime, and `make web` + `make serve-web` build and serve a browser page
(`web/index.html`) that runs every estimator locally on one or many uploaded
files. It is published at https://hyphaltip.github.io/subopt-kaks/ by the
Pages workflow on every push. See `rust/subopt-kaks/README.md`.

### Usage

```sh
kaks input.cds.fasta                        # unaligned CDS -> align -> dN/dS
kaks --prealigned input.cds.aln             # in-frame codon alignment
kaks --method both --format table in.fasta  # YN00 + NG86, readable columns
kaks --method ml --prealigned in.aln        # codeml-style pairwise ML (F3x4)
cat in.fasta | kaks --format json -         # JSON from stdin
kaks --aln-mode local --matrix BLOSUM45 --gap-open 10 --gap-ext 2 in.fasta
kaks --genetic-code 2 in.fasta              # NCBI transl_table 2, vertebrate mito
kaks --list-codes                           # all 27 supported genetic codes
```

Every pair of sequences in the input is reported. Output columns are
`SEQ1 SEQ2 dN dS OMEGA N S kappa t LENGTH` (plus `SE_dN SE_dS NG_dN NG_dS
NG_omega` with `--method both`). `kaks --help` lists all options; see
`rust/subopt-kaks/README.md` for details.

### Validation against PAML (2026-11-05)

`tst/paml_compare/compare_paml.py` ran PAML 4.10.10 `codeml` (pairwise,
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
method: `tst/paml_compare/README.md`.

## Tests

`tst/run_tests.sh` runs each tool on the bundled *C. elegans* / *C. briggsae*
pair and HIV pair and compares the output with `tst/expected/`. Reference
values for the bundled pair (universal code, YN00):

```
dN=0.066778  dS=1.472123  omega=0.045362  N=832.180782  S=289.819218  kappa=1.564887  t=1.289361
```

## Changes

See `Changelog`.

Jason Stajich
