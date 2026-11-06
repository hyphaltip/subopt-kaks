#!/usr/bin/env python3
"""Compare kaks against PAML codeml (runmode = -2) and yn00 on a set of alignments.

    compare_paml.py --paml-src DIR --codeml BIN --yn00 BIN --kaks BIN \
                    --workdir DIR [--out report.tsv] [--datasets NAME,...]

For every dataset x genetic code x CodonFreq the script
  1. converts the PAML sequence file to FASTA (expanding '.' identity dots),
  2. runs codeml pairwise ML (CodonFreq 0-3, getSE = 1, Small_Diff = 1e-6),
     yn00 (weighting = 0), and `kaks --method all` with the matching
     --genetic-code and --codon-freq,
  3. joins the pairwise results by sequence name and compares lnL, t, kappa,
     omega, dN, dS, S, N and the standard errors against tolerances.

It also checks PAML's GeneticCode[][64] table (tools.c) against the NCBI
ncbieaa strings in gc.prt to document which PAML icode maps to which NCBI id
and where the two disagree.

Exit status is 1 if any value outside tolerance was found.
"""
import argparse
import json
import os
import re
import subprocess
import sys
from pathlib import Path

# (relative path in paml/examples, PAML icode, NCBI transl_table id, label)
DATASETS = [
    ("abglobin.nuc", 0, 1, "abglobin_5seq_interleaved_dots"),
    ("lysozyme/lysozymeSmall.nuc", 0, 1, "lysozymeSmall_7seq_ambiguity"),
    ("HIVNSsites/HIVenvSweden.txt", 0, 1, "HIVenv_13seq"),
    ("myxo/myxovirus.aln", 0, 1, "myxovirus_12seq"),
    ("CladeModelCD/ECP_EDN_15.nuc", 0, 1, "ECP_EDN_15seq"),
    ("lysin/lysin.nuc", 0, 1, "lysin_25seq"),
    ("mtCDNA/mtCDNApri.nuc", 1, 2, "mtCDNApri_7seq_vertmito"),
    ("mtCDNAape/mtCDNAape.txt", 1, 2, "mtCDNAape_6seq_vertmito"),
    # the same nuclear alignment re-analysed under other codes that add no stops
    ("lysozyme/lysozymeSmall.nuc", 3, 4, "lysozymeSmall_moldmito"),
    ("lysozyme/lysozymeSmall.nuc", 4, 5, "lysozymeSmall_invertmito"),
    ("lysozyme/lysozymeSmall.nuc", 5, 6, "lysozymeSmall_ciliate"),
    ("lysozyme/lysozymeSmall.nuc", 6, 9, "lysozymeSmall_echinoderm"),
    ("lysozyme/lysozymeSmall.nuc", 7, 10, "lysozymeSmall_euplotid"),
    ("lysozyme/lysozymeSmall.nuc", 8, 12, "lysozymeSmall_altyeast"),
    ("lysozyme/lysozymeSmall.nuc", 9, 13, "lysozymeSmall_ascidian"),
    ("lysozyme/lysozymeSmall.nuc", 10, 15, "lysozymeSmall_blepharisma"),
    ("abglobin.nuc", 2, 3, "abglobin_yeastmito"),
]

CODONFREQ = {0: "fequal", 1: "f1x4", 2: "f3x4", 3: "fcodon"}

# absolute tolerances on codeml's printed precision
TOL = {
    "lnL": 1e-3,
    "t": 2e-4,
    "kappa": 2e-4,
    "omega": 2e-4,
    "dN": 2e-4,
    "dS": 2e-4,
    "S": 0.06,
    "N": 0.06,
}
SE_REL_TOL = 0.15  # standard errors: finite-difference Hessians, reported only
YN_TOL = 1.5e-3    # yn00 prints 4 decimals


# --------------------------------------------------------------------------- IO
def parse_paml_seqs(path):
    """Return [(name, seq)] from a PAML/PHYLIP sequence file (sequential or
    interleaved, optional '.' identity dots, spaces between codons)."""
    lines = Path(path).read_text(errors="replace").splitlines()
    while lines and not lines[0].strip():
        lines.pop(0)
    head = lines[0].split()
    ns, ls = int(head[0]), int(head[1])
    opts = "".join(head[2:]).upper()
    body = [l for l in lines[1:] if l.strip()]
    names, seqs = [], []

    def clean(chunk):
        return re.sub(r"[^A-Za-z.\-?]", "", chunk)

    if "I" in opts:
        # first block: names (sequence may follow on the same line)
        for l in body[:ns]:
            parts = l.split(None, 1)
            names.append(parts[0])
            seqs.append(clean(parts[1]) if len(parts) > 1 else "")
        k = 0
        for l in body[ns:]:
            if re.fullmatch(r"[\d\s]+", l):
                continue  # site-number lines
            seqs[k % ns] += clean(l)
            k += 1
            if all(len(s) >= ls for s in seqs):
                break
    else:
        tokens = []
        for l in body:
            tokens.extend(l.split())
        i = 0
        for _ in range(ns):
            names.append(tokens[i])
            i += 1
            s = ""
            while len(s) < ls and i < len(tokens):
                if re.fullmatch(r"\d+", tokens[i]):
                    i += 1
                    continue
                s += clean(tokens[i])
                i += 1
            seqs.append(s)
    out = []
    for n, s in zip(names, seqs):
        s = s[:ls]
        if len(s) != ls:
            raise ValueError(f"{path}: {n} has {len(s)} chars, expected {ls}")
        out.append((n, s))
    ref = out[0][1]
    expanded = []
    for n, s in out:
        s = "".join(ref[i] if c == "." else c for i, c in enumerate(s))
        expanded.append((n, s.upper()))
    return expanded


def write_fasta(recs, path):
    with open(path, "w") as fh:
        for n, s in recs:
            fh.write(f">{n}\n")
            for i in range(0, len(s), 60):
                fh.write(s[i:i + 60] + "\n")


def run(cmd, cwd=None, timeout=3600):
    r = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, timeout=timeout)
    if r.returncode != 0:
        raise RuntimeError(f"{' '.join(map(str, cmd))} failed ({r.returncode}):\n{r.stderr[-2000:]}")
    return r.stdout


# ------------------------------------------------------------------ PAML runs
def run_codeml(codeml, seqfile, icode, codonfreq, wd):
    wd.mkdir(parents=True, exist_ok=True)
    (wd / "codeml.ctl").write_text(
        f"""seqfile = {seqfile}
outfile = mlc
noisy = 0
verbose = 0
runmode = -2
seqtype = 1
CodonFreq = {codonfreq}
model = 0
NSsites = 0
icode = {icode}
fix_kappa = 0
kappa = 2
fix_omega = 0
omega = 0.4
fix_alpha = 1
alpha = 0
getSE = 1
RateAncestor = 0
Small_Diff = 1e-6
cleandata = 1
"""
    )
    run([codeml, "codeml.ctl"], cwd=wd)
    return parse_mlc(wd / "mlc")


def parse_mlc(path):
    """{(name_i, name_j): {...}} from the pairwise section of codeml's mlc."""
    txt = Path(path).read_text()
    res = {}
    blocks = re.split(r"\n(?=\d+ \([^)]*\) \.\.\. \d+ \([^)]*\))", txt)
    for b in blocks[1:]:
        m = re.match(r"\d+ \(([^)]*)\) \.\.\. \d+ \(([^)]*)\)", b)
        n1, n2 = m.group(1).strip(), m.group(2).strip()
        lnl = float(re.search(r"lnL =\s*([-\d.]+)", b).group(1))
        pars = re.search(r"lnL =.*\n\s*([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)", b)
        t, kappa, omega = map(float, pars.groups())
        se = re.search(r"SEs for parameters:\n\s*([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)", b)
        se_t, se_k, se_w = map(float, se.groups()) if se else (None,) * 3
        st = re.search(r"t=\s*([-\d.]+)\s+S=\s*([-\d.]+)\s+N=\s*([-\d.]+)\s+dN/dS=\s*([-\d.]+)\s+dN =\s*([-\d.]+)\s+dS =\s*([-\d.]+)", b)
        S, N = float(st.group(2)), float(st.group(3))
        m1 = re.search(r"dN =\s*([-\d.]+) \+- ([-\d.]+)\s+dS =\s*([-\d.]+) \+- ([-\d.]+) \(by method 1\)", b)
        dN, se_dN, dS, se_dS = map(float, m1.groups())
        res[(n1, n2)] = dict(lnL=lnl, t=t, kappa=kappa, omega=omega, S=S, N=N, dN=dN, dS=dS,
                             se_t=se_t, se_kappa=se_k, se_omega=se_w, se_dN=se_dN, se_dS=se_dS)
    return res


def run_yn00(yn00, seqfile, icode, wd):
    wd.mkdir(parents=True, exist_ok=True)
    (wd / "yn00.ctl").write_text(
        f"seqfile = {seqfile}\noutfile = yn\nverbose = 0\nicode = {icode}\nweighting = 0\ncommonf3x4 = 0\n"
    )
    run([yn00, "yn00.ctl"], cwd=wd)
    txt = (wd / "yn").read_text()
    names = re.findall(r"^\s*(\S+)\s+\d+\s+[\d.]+", txt, re.M)  # not reliable; use rst-style table below
    # table: "seq. seq.     S       N        t   kappa   omega     dN +- SE    dS +- SE"
    m = re.search(r"seq\. seq\.\s+S\s+N\s+t\s+kappa\s+omega\s+dN \+- SE\s+dS \+- SE\n(.*?)(?:\n\n|\Z)", txt, re.S)
    rows = {}
    if m:
        for line in m.group(1).strip().splitlines():
            f = line.split()
            if len(f) < 11:
                continue
            i, j = int(f[0]), int(f[1])
            rows[(i, j)] = dict(S=float(f[2]), N=float(f[3]), t=float(f[4]), kappa=float(f[5]),
                                omega=float(f[6]), dN=float(f[7]), se_dN=float(f[9]),
                                dS=float(f[10]), se_dS=float(f[12]))
    return rows


# ------------------------------------------------------------------ kaks run
def run_kaks(kaks, fasta, ncbi, codonfreq):
    out = run([kaks, "--prealigned", "--method", "all", "--format", "json",
               "--genetic-code", str(ncbi), "--codon-freq", CODONFREQ[codonfreq], str(fasta)])
    rows = {}
    for r in json.loads(out):
        rows[(r["seq1"], r["seq2"])] = r
    return rows


# ------------------------------------------------------------- comparisons
def cmp_val(label, a, b, tol, fails, notes, rel=2e-3):
    """Fail when |a-b| exceeds both the absolute tolerance and `rel`*|a|
    (the relative part covers saturated pairs with t or dS in the tens)."""
    if a is None or b is None:
        return
    d = abs(a - b)
    if d > tol and d > rel * abs(a):
        fails.append(f"{label}: codeml={a} kaks={b} diff={d:.3g} > {tol}")


def compare_dataset(args, ds, report):
    rel, icode, ncbi, label = ds
    seqfile = Path(args.paml_src) / "examples" / rel
    wd = Path(args.workdir) / label
    wd.mkdir(parents=True, exist_ok=True)
    recs = parse_paml_seqs(seqfile)
    fasta = wd / "aln.fasta"
    write_fasta(recs, fasta)
    order = [n for n, _ in recs]
    nfail = 0
    npairs = 0

    # yn00
    try:
        yn = run_yn00(args.yn00, seqfile.resolve(), icode, wd / "yn00")
    except RuntimeError as e:
        report.append([label, "yn00", "-", "ERROR", str(e).splitlines()[-1][:120]])
        yn = {}
    kk = run_kaks(args.kaks, fasta, ncbi, 2)
    for (i, j), y in yn.items():
        n1, n2 = order[i - 1], order[j - 1]
        k = kk.get((n2, n1)) or kk.get((n1, n2))
        if k is None:
            report.append([label, "yn00", f"{n1}-{n2}", "MISSING", ""])
            nfail += 1
            continue
        fails = []
        for key, kkey in [("dN", "dN"), ("dS", "dS"), ("omega", "omega"), ("S", "S"), ("N", "N"),
                          ("t", "t"), ("kappa", "kappa"), ("se_dN", "SE_dN"), ("se_dS", "SE_dS")]:
            tol = YN_TOL if key in ("dN", "dS", "omega", "t", "kappa", "se_dN", "se_dS") else 0.06
            cmp_val(key, y[key], k[kkey], tol, fails, None)
        npairs += 1
        if fails:
            nfail += 1
            report.append([label, "yn00", f"{n1}-{n2}", "FAIL", "; ".join(fails)])

    # codeml for each CodonFreq
    for cf in range(4):
        try:
            ml = run_codeml(args.codeml, seqfile.resolve(), icode, cf, wd / f"codeml_cf{cf}")
        except RuntimeError as e:
            report.append([label, f"ml_cf{cf}", "-", "ERROR", str(e).splitlines()[-1][:120]])
            nfail += 1
            continue
        try:
            kk = run_kaks(args.kaks, fasta, ncbi, cf)
        except RuntimeError as e:
            report.append([label, f"ml_cf{cf}", "-", "KAKS_ERROR", str(e).splitlines()[-1][:160]])
            nfail += 1
            continue
        for (n1, n2), c in ml.items():
            k = kk.get((n2, n1)) or kk.get((n1, n2))
            if k is None:
                report.append([label, f"ml_cf{cf}", f"{n1}-{n2}", "MISSING", ""])
                nfail += 1
                continue
            fails, warn = [], []
            cmp_val("lnL", c["lnL"], k["lnL"], TOL["lnL"], fails, None)
            cmp_val("t", c["t"], k["ML_t"], TOL["t"], fails, None)
            cmp_val("kappa", c["kappa"], k["ML_kappa"], TOL["kappa"], fails, None)
            cmp_val("omega", c["omega"], k["ML_omega"], TOL["omega"], fails, None)
            cmp_val("dN", c["dN"], k["ML_dN"], TOL["dN"], fails, None)
            cmp_val("dS", c["dS"], k["ML_dS"], TOL["dS"], fails, None)
            cmp_val("S", c["S"], k["ML_S"], TOL["S"], fails, None)
            cmp_val("N", c["N"], k["ML_N"], TOL["N"], fails, None)
            for key, kkey in [("se_dN", "ML_SE_dN"), ("se_dS", "ML_SE_dS"), ("se_t", "ML_SE_t"),
                              ("se_kappa", "ML_SE_kappa"), ("se_omega", "ML_SE_omega")]:
                a, b = c[key], k[kkey]
                if a is not None and a > 0 and abs(a - b) / a > SE_REL_TOL:
                    warn.append(f"{key}: codeml={a} kaks={b}")
            npairs += 1
            if fails:
                nfail += 1
                report.append([label, f"ml_cf{cf}", f"{n1}-{n2}", "FAIL", "; ".join(fails)])
            elif warn:
                report.append([label, f"ml_cf{cf}", f"{n1}-{n2}", "SE_WARN", "; ".join(warn)])
    report.append([label, "summary", f"{len(order)} seqs", f"{npairs} comparisons", f"{nfail} failures"])
    return nfail


def check_genetic_codes(paml_src, gcprt, report):
    """Compare PAML GeneticCode[][64] with NCBI gc.prt tables."""
    aas = "ARNDCQEGHILKMFPSTWYV"
    tools = (Path(paml_src) / "src" / "tools.c").read_text(errors="replace")
    m = re.search(r"int GeneticCode\[\]\[64\] =\s*\{(.*?)\n\};", tools, re.S)
    body = m.group(1)
    tables = []
    for tm in re.finditer(r"\{([^{}]*)\}\s*,?\s*/\*\s*(\d+):([^*]*)\*/", body):
        nums = [int(x) for x in tm.group(1).replace("\n", " ").split(",") if x.strip()]
        s = "".join("*" if v == -1 else aas[v] for v in nums)
        tables.append((int(tm.group(2)), tm.group(3).strip(), s))
    ncbi = {}
    txt = Path(gcprt).read_text()
    for b in re.findall(r"\{\s*(name\s+\"[^}]*?)\}", txt, re.S):
        ident = int(re.search(r"\bid\s+(\d+)", b).group(1))
        ncbi[ident] = re.search(r"ncbieaa\s+\"([^\"]+)\"", b).group(1)
    nfail = 0
    for icode, name, s in tables:
        best = [(sum(a != b for a, b in zip(s, t)), i) for i, t in ncbi.items()]
        best.sort()
        ndiff, nid = best[0]
        diffs = []
        if ndiff:
            t = ncbi[nid]
            for k in range(64):
                if s[k] != t[k]:
                    codon = "TCAG"[k // 16] + "TCAG"[(k // 4) % 4] + "TCAG"[k % 4]
                    diffs.append(f"{codon}:PAML={s[k]},NCBI={t[k]}")
        status = "OK" if ndiff == 0 else "DIFF"
        report.append(["genetic_code", f"paml_icode_{icode}", name, f"NCBI {nid} {status}", " ".join(diffs)])
    return nfail


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--paml-src", required=True, help="PAML checkout (with src/ and examples/)")
    ap.add_argument("--codeml", required=True)
    ap.add_argument("--yn00", required=True)
    ap.add_argument("--kaks", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--gcprt", help="NCBI gc.prt to check PAML's code tables against")
    ap.add_argument("--out", default="-")
    ap.add_argument("--datasets", help="comma-separated labels to run (default: all)")
    args = ap.parse_args()

    report = [["dataset", "test", "pair", "status", "detail"]]
    total_fail = 0
    if args.gcprt:
        total_fail += check_genetic_codes(args.paml_src, args.gcprt, report)
    wanted = set(args.datasets.split(",")) if args.datasets else None
    for ds in DATASETS:
        if wanted and ds[3] not in wanted:
            continue
        print(f"== {ds[3]}", file=sys.stderr, flush=True)
        try:
            total_fail += compare_dataset(args, ds, report)
        except Exception as e:  # keep going; record
            report.append([ds[3], "dataset", "-", "ERROR", repr(e)[:200]])
            total_fail += 1

    out = sys.stdout if args.out == "-" else open(args.out, "w")
    for row in report:
        out.write("\t".join(row) + "\n")
    if out is not sys.stdout:
        out.close()
    print(f"total failures: {total_fail}", file=sys.stderr)
    sys.exit(1 if total_fail else 0)


if __name__ == "__main__":
    main()
