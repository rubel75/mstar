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
`test_files/<test>/run` with `OMP_NUM_THREADS=1`, and its `minv_*.dat` files
are compared with `test_files/<test>/expected_output` by
`compare_outputs.py`. Values agree when
`|a - b| <= 1e-5 + 2e-3 * max(|a|, |b|)`, which allows for rounding in the last
printed digit and for differences between compilers and LAPACK libraries.

With more than one OpenMP thread, the last printed digit can change because
the band sum is accumulated in single precision in a thread-dependent order.

## History of the reference outputs

- 2026-11-12: generated with gfortran 13.3 and reference LAPACK, after the
  fix of the VASP reader that omitted the couplings of band NBANDS_CDER to all
  higher bands. Relative to commit 1a5c1c0, only band 156 and the bands
  degenerate with it change in the VASP test. The WIEN2k test is unchanged.
