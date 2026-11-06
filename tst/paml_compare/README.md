# kaks versus PAML: pairwise comparison sweep

`compare_paml.py` runs PAML `codeml` (pairwise ML, `runmode = -2`, all four
`CodonFreq` models, `getSE = 1`, `Small_Diff = 1e-6`, `cleandata = 1`) and
`yn00` (`weighting = 0`, `commonf3x4 = 0`) next to `kaks --method all` on
PAML's own example alignments, and compares every sequence pair. It also
checks PAML's `GeneticCode[][64]` table against the NCBI `gc.prt` tables that
`kaks` embeds.

## Running it

```sh
git clone https://github.com/abacus-gene/paml /some/dir/paml
make -C /some/dir/paml/src codeml yn00
curl -o gc.prt https://ftp.ncbi.nlm.nih.gov/entrez/misc/data/gc.prt

tst/paml_compare/compare_paml.py \
    --paml-src /some/dir/paml \
    --codeml /some/dir/paml/src/codeml --yn00 /some/dir/paml/src/yn00 \
    --kaks rust/subopt-kaks/target/release/kaks \
    --workdir /scratch/kaks-sweep --gcprt gc.prt --out report.tsv
```

Runtime is about 25 minutes, almost all of it `codeml`. `--datasets` takes a
comma-separated list of labels to run a subset.

## Datasets

| Label | File (paml/examples) | Sequences | PAML icode | NCBI id | Notes |
|---|---|---|---|---|---|
| abglobin_5seq_interleaved_dots | abglobin.nuc | 5 | 0 | 1 | interleaved, `.` identity dots |
| lysozymeSmall_7seq_ambiguity | lysozyme/lysozymeSmall.nuc | 7 | 0 | 1 | ambiguity codes present |
| HIVenv_13seq | HIVNSsites/HIVenvSweden.txt | 13 | 0 | 1 | |
| myxovirus_12seq | myxo/myxovirus.aln | 12 | 0 | 1 | several saturated pairs |
| ECP_EDN_15seq | CladeModelCD/ECP_EDN_15.nuc | 15 | 0 | 1 | |
| lysin_25seq | lysin/lysin.nuc | 25 | 0 | 1 | 300 pairs |
| mtCDNApri_7seq_vertmito | mtCDNA/mtCDNApri.nuc | 7 | 1 | 2 | 9993 nt, vertebrate mitochondrial |
| mtCDNAape_6seq_vertmito | mtCDNAape/mtCDNAape.txt | 6 | 1 | 2 | 10707 nt |
| lysozymeSmall_{moldmito,invertmito,ciliate,echinoderm,euplotid,altyeast,ascidian,blepharisma} | lysozyme/lysozymeSmall.nuc | 7 | 3,4,5,6,7,8,9,10 | 4,5,6,9,10,12,13,15 | same alignment under codes that introduce no stops |
| abglobin_yeastmito | abglobin.nuc | 5 | 2 | 3 | |

## Tolerances

codeml prints 5 decimals for t, kappa, omega, dN and dS, 6 for lnL and 1 for
S and N. A value fails when it differs by more than the printed precision
(2e-4; lnL 1e-3; S and N 0.06) **and** by more than 0.2 % relative. yn00
prints 4 decimals; its tolerance is 1.5e-3. Standard errors are compared at
15 % relative and reported as `SE_WARN`, never as failures (see below).

## Result, 2026-11-05 (PAML 4.10.10, commit 4c7902f; kaks at this commit)

Full report: `report-2026-11-05.tsv`.

| | |
|---|---|
| dataset x code combinations | 17 |
| pair comparisons (yn00 + 4 x codeml) | 3970 |
| failures | 3 |

All 11 PAML genetic-code tables (icode 0-10) are identical to NCBI tables
1, 2, 3, 4, 5, 6, 9, 10, 12, 13 and 15. PAML's icode 11 ("Ziheng's regular
code") is an artificial table and has no NCBI counterpart.

yn00: 0 failures on 794 pairs across all codes once `kaks` adopted PAML's
default of per-pair F3x4 frequencies (`commonf3x4 = 0`). The C tools in
`src/` pool the frequencies over all sequences; `kaks --common-f3x4`
reproduces that.

codeml ML: 3 failures on 3176 pair x CodonFreq comparisons.

1. `abglobin`, marsupial-rat, F3x4: codeml reported lnL = -1843.676403
   (t = 2.575, kappa = 0.824, omega = 0.063); kaks reports lnL = -1842.638722
   (t = 1.663, kappa = 1.021, omega = 0.123). kaks's likelihood is higher by
   1.04 log units, so codeml stopped at a local optimum. codeml starts each
   pair from randomly jittered values; kaks starts deterministically from the
   NG86 distance.
2. `myxovirus`, Chicken_Mx-Pig_Mx, F3x4: t 26.908 vs 26.972, dS 32.53 vs
   32.61 (0.24 % relative).
3. `myxovirus`, Duck_Mx-Pig_Mx, Fcodon: t 22.682 vs 22.794, dS 26.85 vs
   26.99 (0.5 % relative).

Cases 2 and 3 are saturated pairs (dS above 25) where the likelihood is flat
in t; both programs stop at the same lnL to within their 1e-7 convergence
tolerance.

Standard errors: 165 of the 3176 ML comparisons had at least one S.E. more
than 15 % from codeml's, mostly `se_omega`, `se_kappa` and `se_dS` on short or
saturated pairs. Both programs compute the S.E. from a finite-difference
Hessian with step 2e-6 x (1 + |x|). kaks's values are stable when that step
is changed between 2e-6 and 1e-4; codeml's values for the same pair change
between runs (0.00807 vs 0.00812 for the test pair's SE(dN) with
`Small_Diff` 5e-7 vs 1e-6). Point estimates are unaffected.

## What this does not cover

Only pairwise estimates. Tree-based `codeml` (branch lengths, site models,
branch models, LRTs) is not implemented in `kaks`.
