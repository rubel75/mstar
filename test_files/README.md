# mstar regression tests

Two tests use CoSi, with the same DFT inputs as the BerryCPT tests
(https://github.com/rubel75/berrycpt, folder `test_files`). The inputs are
not duplicated here. By default the scripts expect the BerryCPT repository
next to the mstar repository:

```
<parent>/mstar
<parent>/berrycpt/test_files/WIEN2k
<parent>/berrycpt/test_files/VASP
```

| Test | Input files | Command |
|---|---|---|
| `WIEN2k` | `CoSi-berrycpt.mommat2` (SOC, 12 k-points) | `mstar CoSi-berrycpt.mommat2 1e-5` |
| `VASP` | `WAVEDER`, `EIGENVAL` (17 k-points, 422 bands, NBANDS_CDER = 156) | `mstar WAVEDER` |

The DFT input and output files used to generate these inputs are in the
`etc` folders of the BerryCPT tests.

## Running

Build mstar in the repository root, then run

```bash
test_files/run_tests.sh [INPUT_DIR]
```

`INPUT_DIR` overrides the location of the input folders. Each test runs in
`test_files/<test>/run`, and its `minv_*.dat` files are compared with
`test_files/<test>/expected_output` by `compare_outputs.py`. Values agree when
`|a - b| <= 1e-6 + 1e-5 * max(|a|, |b|)`, which allows for rounding in the
last printed digit and for differences between compilers and LAPACK
libraries. The tolerances can be changed with `--rtol` and `--atol`.

## History of the reference outputs

- 2026-11-12 (pass 3): regenerated after adding the degenerate-block ID
  column and changing to chain grouping of degenerate bands. All values are
  identical to the pass-2 references; the groups are the same for these
  inputs.
- 2026-11-12 (pass 2): regenerated with gfortran 13.3 and reference LAPACK
  after the change to double precision with compensated summation, the
  free-format reading of EIGENVAL, and the 7-digit output format. WIEN2k:
  agrees with the previous reference within 2e-3; the only larger changes are
  off-diagonal elements below 1e-2 at k-point 1, where the change is below
  1e-7 of the largest diagonal element of the same band (single-precision
  rounding in the previous version). VASP: the previous reference is invalid
  because EIGENVAL was truncated to 0.01 eV. With that reading fixed in the
  previous code, the new results agree within 2e-3 except two off-diagonal
  elements (k-point 4, bands 89 and 90), where the change is below 1e-7 of
  the largest diagonal element.
- 2026-11-12 (pass 1, commit c737808): generated with gfortran 13.3 and
  reference LAPACK, after the fix of the VASP reader that omitted the
  couplings of band NBANDS_CDER to all higher bands.
