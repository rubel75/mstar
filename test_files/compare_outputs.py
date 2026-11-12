#!/usr/bin/env python3
"""Compare mstar output files with reference outputs.

Usage:
    compare_outputs.py REFERENCE_DIR RESULT_DIR [--rtol R] [--atol A]

All minv_*.dat files found in REFERENCE_DIR are compared with the files of
the same name in RESULT_DIR. Records are matched by k-point and band index.
Numerical columns agree when |a - b| <= atol + rtol * max(|a|, |b|).
Comment lines other than the "# KP:" headers are ignored, so changes in the
descriptive file header do not cause a failure. Text annotations on a data
line (for example the hyperbolic-dispersion warning in minv_d.dat) must match.

The defaults (rtol = 1e-5, atol = 1e-6) allow for rounding in the last
printed digit (7 significant digits) and for differences between compilers
and LAPACK libraries. Use larger values when comparing with outputs of
mstar versions before 2026-11, which were written with 4 significant digits
and accumulated in single precision (e.g., --rtol 2e-3 --atol 1e-5).
Exit status is 0 if all files agree and 1 otherwise.
"""

import argparse
import glob
import os
import re
import sys

KP_RE = re.compile(r"#\s*KP:\s*(\d+)")


def parse(fname):
    """Return {(kpoint, band): (numbers, text)} and the list of KP headers."""
    records = {}
    headers = []
    kpt = None
    with open(fname) as f:
        for line in f:
            m = KP_RE.match(line)
            if m:
                kpt = int(m.group(1))
                headers.append(" ".join(line.split()))
                continue
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split()
            band = int(fields[0])
            numbers, text = [], []
            for x in fields[1:]:
                try:
                    numbers.append(float(x))
                except ValueError:
                    text.append(x)
            records[(kpt, band)] = (numbers, " ".join(text))
    return records, headers


def compare(ref_file, new_file, rtol, atol):
    ref, ref_headers = parse(ref_file)
    new, new_headers = parse(new_file)
    problems = []
    if ref_headers != new_headers:
        problems.append("k-point headers differ")
    missing = sorted(set(ref) - set(new))
    extra = sorted(set(new) - set(ref))
    if missing:
        problems.append(f"{len(missing)} records missing, first {missing[0]}")
    if extra:
        problems.append(f"{len(extra)} unexpected records, first {extra[0]}")
    worst = (0.0, None)
    nbad = 0
    for key in sorted(set(ref) & set(new)):
        a, ta = ref[key]
        b, tb = new[key]
        if len(a) != len(b) or ta != tb:
            nbad += 1
            problems.append(f"record {key} has a different structure")
            continue
        for i, (x, y) in enumerate(zip(a, b)):
            diff = abs(x - y)
            tol = atol + rtol * max(abs(x), abs(y))
            if diff > tol:
                nbad += 1
                excess = diff / tol
                if excess > worst[0]:
                    worst = (excess, (key, i + 1, x, y))
    if nbad:
        (k, b), col, x, y = worst[1] if worst[1] else ((None, None), 0, 0, 0)
        problems.append(
            f"{nbad} values outside tolerance; worst at k-point {k}, band {b}, "
            f"column {col}: reference {x:.6e}, result {y:.6e}")
    return problems


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("reference_dir")
    p.add_argument("result_dir")
    p.add_argument("--rtol", type=float, default=1e-5)
    p.add_argument("--atol", type=float, default=1e-6)
    args = p.parse_args()

    ref_files = sorted(glob.glob(os.path.join(args.reference_dir, "minv_*.dat")))
    if not ref_files:
        print(f"No minv_*.dat files found in {args.reference_dir}")
        return 1
    failed = False
    for ref_file in ref_files:
        name = os.path.basename(ref_file)
        new_file = os.path.join(args.result_dir, name)
        if not os.path.isfile(new_file):
            print(f"FAIL {name}: not produced")
            failed = True
            continue
        problems = compare(ref_file, new_file, args.rtol, args.atol)
        if problems:
            failed = True
            print(f"FAIL {name}")
            for msg in problems:
                print(f"     {msg}")
        else:
            print(f"OK   {name}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
