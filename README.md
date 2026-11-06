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

## Rust port

`rust/subopt-kaks` is a single-binary reimplementation (`kaks`) with flags
only, embedded BLOSUM62 and genetic-code tables, FASTA/PHYLIP input, and
TSV/table/JSON output. Its `--prealigned` path reproduces
`yn00_cds_prealigned` byte-for-byte on the bundled test pair. See
`rust/subopt-kaks/README.md`.

```sh
make rust                      # cargo build --release
make test-rust                 # cargo test + regression tests
make install-rust PREFIX=~/.local
```

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
