#!/usr/bin/env bash
# Regression tests for subopt-kaks.
#
#   tst/run_tests.sh [BINDIR]          test the C tools (default: ../src)
#   tst/run_tests.sh --rust [BINDIR]   test the Rust `kaks` binary
#
# Every test runs a tool on the bundled F35G8.1/CBG16160 pair and compares
# the output byte-for-byte with tst/expected/<name>.tsv.
set -u

TST=$(cd "$(dirname "$0")" && pwd)
mode=c
if [ "${1:-}" = "--rust" ]; then mode=rust; shift; fi
if [ "$mode" = rust ]; then
    BIN=${1:-$TST/../rust/subopt-kaks/target/release}
else
    BIN=${1:-$TST/../src}
fi
BIN=$(cd "$BIN" && pwd) || { echo "no such directory: $BIN" >&2; exit 2; }

tmp=$(mktemp -d "${SCRATCH:-${TMPDIR:-/tmp}}/subopt-test.XXXXXX")
trap 'rm -rf "$tmp"' EXIT

pass=0; fail=0
check() {   # check NAME COMMAND...
    local name=$1; shift
    local exp="$TST/expected/$name.tsv"
    if "$@" > "$tmp/$name.out" 2> "$tmp/$name.err"; then
        if cmp -s "$tmp/$name.out" "$exp"; then
            echo "PASS $name"; pass=$((pass+1)); return
        fi
        echo "FAIL $name: output differs from $exp"
        diff "$exp" "$tmp/$name.out" | head -20
    else
        echo "FAIL $name: exit status $? ($*)"
        head -20 "$tmp/$name.err"
    fi
    fail=$((fail+1))
}

cd "$TST" || exit 2
# never let an inherited environment pick a different matrix directory
unset SUBOPTDIR SUBOPT_DIR SUBOPT_MATRIX_DIR

if [ "$mode" = c ]; then
    check prealigned      "$BIN/yn00_cds_prealigned" F35G8.1_CBG16160.cds.aln
    check prealigned_x2   "$BIN/yn00_cds_prealigned" F35G8.1_CBG16160.cds.aln F35G8.1_CBG16160.cds.aln
    check optimal_local   "$BIN/yn00_cds_optimal" --showtable F35G8.1_CBG16160.cds.fas
    check optimal_global  "$BIN/yn00_cds_optimal" --showtable --global F35G8.1_CBG16160.cds.fas
    check optimal_gapsign "$BIN/yn00_cds_optimal" --showtable --noheader --gapopen 10 --gapext 2 F35G8.1_CBG16160.cds.fas
    check optimal_hiv     "$BIN/yn00_cds_optimal" --showtable HIV_test.fa
    check optimal_aln     "$BIN/yn00_cds_optimal" F35G8.1_CBG16160.cds.fas
    check zuker_pep       "$BIN/zuker_align" --max_to_show 2 F35G8.1_CBG16160.pep
    check cds_zuker       "$BIN/cds_zuker_align" --max_to_show 2 --output - F35G8.1_CBG16160.cds.fas
    check yn00_zuker      "$BIN/yn00_zuker_align" --sorted --max_to_show 5 --output - F35G8.1_CBG16160.cds.fas
else
    check kaks_prealigned "$BIN/kaks" --prealigned F35G8.1_CBG16160.cds.aln
    check kaks_cds        "$BIN/kaks" F35G8.1_CBG16160.cds.fas
fi

echo "$pass passed, $fail failed"
[ "$fail" -eq 0 ]
