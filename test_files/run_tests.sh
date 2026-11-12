#!/bin/bash
# Run the mstar regression tests and compare with the reference outputs.
#
# Usage (from any directory):
#   test_files/run_tests.sh [INPUT_DIR]
#
# INPUT_DIR holds the DFT input files in the subfolders WIEN2k and VASP.
# The CoSi inputs are shared with BerryCPT, so the default is the BerryCPT
# test folder next to the mstar repository:
#   ../berrycpt/test_files   (relative to the mstar repository root)
#
# The mstar executable is taken from the repository root. Each test runs in
# test_files/<test>/run, which is removed and recreated.

set -u
here=$(cd "$(dirname "$0")" && pwd)
root=$(cd "$here/.." && pwd)
mstar="$root/mstar"
inputs=$(cd "${1:-$root/../berrycpt/test_files}" 2>/dev/null && pwd) || {
    echo "Input folder not found: ${1:-$root/../berrycpt/test_files}"
    exit 1
}

if [ ! -x "$mstar" ]; then
    echo "mstar executable not found: $mstar (run make first)"
    exit 1
fi

status=0

# run_test NAME INPUT_SUBDIR REFERENCE_SUBDIR "INPUT FILES" MSTAR_ARGUMENTS...
run_test() {
    local name=$1 ; shift
    local input_subdir=$1 ; shift
    local reference_subdir=$1 ; shift
    local input_files=$1 ; shift
    local rundir="$here/$name/run"
    echo "=== $name"
    rm -rf "$rundir"
    mkdir -p "$rundir"
    for f in $input_files; do
        if [ ! -f "$inputs/$input_subdir/$f" ]; then
            echo "Missing input file: $inputs/$input_subdir/$f"
            status=1
            return
        fi
        ln -s "$inputs/$input_subdir/$f" "$rundir/$f" 2>/dev/null || \
            cp "$inputs/$input_subdir/$f" "$rundir/$f"
    done
    ( cd "$rundir" && "$mstar" "$@" > mstar.log 2>&1 )
    local rc=$?
    if [ $rc -ne 0 ]; then
        echo "FAIL mstar exited with status $rc (see $rundir/mstar.log)"
        status=1
        return
    fi
    python3 "$here/compare_outputs.py" "$here/$reference_subdir/expected_output" "$rundir" || status=1
}

run_test WIEN2k WIEN2k WIEN2k "CoSi-berrycpt.mommat2" CoSi-berrycpt.mommat2 1e-5
run_test VASP VASP VASP "WAVEDER EIGENVAL" WAVEDER
# --nbmax 3 cuts through the fourfold block of bands 1-4 at k-point 1
run_test WIEN2k_nbmax3 WIEN2k WIEN2k_nbmax3 "CoSi-berrycpt.mommat2" \
    CoSi-berrycpt.mommat2 1e-5 --nbmax 3
# --nbmax above NBANDS_CDER (156) must not change the VASP results; the input
# is given with its path (EIGENVAL is taken from the same directory)
run_test VASP_nbmax1000 VASP VASP "" --nbmax 1000 "$inputs/VASP/WAVEDER"

if [ $status -eq 0 ]; then
    echo "All tests passed."
else
    echo "Some tests failed."
fi
exit $status
