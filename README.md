# mstar

mstar calculates band-resolved effective mass tensors from density-functional-theory momentum (velocity) matrix elements using degenerate k·p perturbation theory.

Supported electronic-structure codes:

- [WIEN2k](https://www.wien2k.at)
- [VASP](https://www.vasp.at)

mstar is written in Fortran and is intended for Linux-based systems. The method and its validation for WIEN2k are described in [Comp. Phys. Commun. 261, 107648 (2021)](https://doi.org/10.1016/j.cpc.2020.107648).

## Status

The `master` branch is the reference version. Every output file records the Git revision of the build in its first header line, and the same revision is printed at the start and end of each run. Please report results together with that revision.

**VASP results obtained with versions before November 2026 should be recalculated.** Earlier versions read `EIGENVAL` with a fixed VASP 5 format, which truncated the eigenvalues to 0.01 eV for files written by VASP 6.4.3 (and possibly other VASP 6 versions), and they omitted the couplings of band `NBANDS_CDER` to all higher bands. WIEN2k results of earlier versions are not affected, apart from rounding in the last digits. See [Changes](#changes) for details.

## Installation

Clone the repository:

```bash
git clone https://github.com/rubel75/mstar
cd mstar
```

The `makefile` is set up for the Intel Fortran compiler (`ifort`) and Intel MKL, and requires the environment variable `MKLROOT`. Settings for GNU Fortran with MKL or OpenBLAS, and a debugging build, are included as commented lines. With oneAPI 2025 or later, use `ifx` (`make FC=ifx`). Compile with:

```bash
make
```

Useful build targets are `make clean` (object and module files) and `make veryclean` (also the executable and the generated `version_mod.f90`). Only a Fortran compiler and a LAPACK library are required.

## Input preparation

mstar needs momentum (WIEN2k) or velocity (VASP) matrix elements between all pairs of bands at each k-point, from a calculation with _plenty_ of empty bands. Step-by-step instructions are available on the Wiki:

- [Generating WIEN2k `case.mommat2` files](https://github.com/rubel75/mstar/wiki/Generate-case.mommat2-file-in-WIEN2k)
- [Generating VASP `WAVEDER` files](https://github.com/rubel75/mstar/wiki/Generate-WAVEDER-file-in-VASP)

### WIEN2k

- Writing of matrix elements is not the default. In `case.inop`, change `OFF` to `ON` in the line `ON           ON/OFF   writes MME to unit 4`.
- Include empty states up to about 5 Ry above the Fermi energy: increase `de` in `case.in1(c)` and, with spin-orbit coupling (SOC), `Emax` in `case.inso`. Set `Emax` in `case.inop` to at least the same value.
- The band range of `case.mommat2` must start at band 1 (`NEMIN = 1`). If it does not, there are valence states below `Emin` of `case.inop`. Lower `Emin` in that case (e.g., to -8 Ry).
- Only the free format of `case.mommat2` written by WIEN2k v20.1 and later is supported. For older WIEN2k versions, use an older release of mstar.
- For states with strong d character, high-energy local orbitals (HELOs) may be needed to converge the masses (see the CuI example in the reference).

### VASP

- Use `LOPTICS = .TRUE.` and `LPEAD = .FALSE.`. With `LPEAD = .TRUE.`, the derivatives are obtained by finite differences on the k-point mesh, which is usually too coarse for effective masses.
- Increase `NBANDS` (at least three times the number of occupied bands as a starting point) and test the convergence.
- GW PAW potentials are recommended, since they describe the high-energy empty states entering the perturbation sums better.
- `WAVEDER` and `EIGENVAL` must come from the same calculation and must be in the same directory.

### Convergence

The masses converge slowly with the number of empty bands, in particular for weakly dispersive bands (e.g., the longitudinal conduction-band mass of Si) and for states with d character. The number of bands, the energy range of the matrix elements, and, in WIEN2k, the local orbitals must be convergence-tested for each material.

## Command-line syntax

```bash
mstar MATRIX_FILE [TOLERANCE] [--nbmax N]
```

- `MATRIX_FILE` is `case.mommat2`, `case.mommat2up`, or `case.mommat2dn` for WIEN2k, or `WAVEDER` (or `WAVEDERF`) for VASP. The file may be given with a path. VASP is recognized by the file name. `EIGENVAL` is read from the same directory as `WAVEDER`, and the binary `WAVEDER` is read when `WAVEDERF` is given.
- `TOLERANCE` (optional) is the degeneracy energy tolerance in Ha (default `1e-6`). Adjacent bands whose energies differ by no more than this value are treated as degenerate (see [Degenerate bands](#degenerate-bands)).
- `--nbmax N` (optional, anywhere on the command line) restricts the output to bands 1 to N. In VASP, the output is limited to the `NBANDS_CDER` bands of `WAVEDER` in any case. All bands still enter the perturbation sums, so the masses of the bands written do not change. A degenerate block cut by N is excluded. This avoids writing masses of high-lying bands, which are not converged with respect to the number of empty states.

The built-in help is displayed with `mstar --help`.

Examples:

```bash
# WIEN2k
/path/to/mstar/mstar case.mommat2
/path/to/mstar/mstar case.mommat2up 1e-5 --nbmax 60

# VASP
/path/to/mstar/mstar WAVEDER
/path/to/mstar/mstar /path/to/calculation/WAVEDER 1e-5
```

WIEN2k (v20.1 and later) also includes a version of mstar that is run as `x mstar [-up/-dn] [-settol 1.0e-5]`. That version does not contain the changes listed under [Changes](#changes).

In spin-polarized WIEN2k calculations, and in WIEN2k calculations with SOC (which are run as spin-polarized for the OPTIC step, see the Wiki), each `case.mommat2up/dn` file is processed in a separate run. A spin-polarized VASP `WAVEDER` is processed in one run.

## Method

### Effective mass of a non-degenerate band

The inverse effective mass tensor of band n follows from second-order k·p perturbation theory [Eq. (5) of the reference],

```math
\frac{m_0}{m^*_{\alpha\beta,n}} = \delta_{\alpha\beta}
+ \frac{1}{m_0} \sum_{l \neq n}
\frac{p^{(\alpha)}_{nl}\, p^{(\beta)}_{ln} + p^{(\beta)}_{nl}\, p^{(\alpha)}_{ln}}{E_n - E_l} ,
```

where α, β are Cartesian directions, $p^{(\alpha)}_{nl}$ are momentum matrix elements between cell-periodic Bloch states at the same k, and the sum runs over all bands available in the input. Negative values correspond to a downward curvature (e.g., holes at a valence band maximum).

### Matrix elements and limitations

The general form of the first term is $(m_0/\hbar^2)\langle u_n | \partial^2 \hat{H}_\mathbf{k} / \partial k_\alpha \partial k_\beta | u_n \rangle$, and the momentum operator in the sum is, more generally, the velocity operator $m_0\hat{\mathbf{v}} = (m_0/\hbar)\,\partial \hat{H}_\mathbf{k}/\partial \mathbf{k}$. This leads to the following differences between the two codes.

- **WIEN2k** (`OPTIC`) provides momentum matrix elements. For a local potential (LDA, GGA, TB-mBJ), momentum and velocity matrix elements coincide, and the first term is exactly $\delta_{\alpha\beta}$. A nonlocal potential adds $i[\hat{V}_\mathrm{NL}, \hat{\mathbf{r}}]$ to the velocity operator, which WIEN2k does not provide. Hybrid functionals are therefore not supported (with momentum matrix elements, the masses would be systematically too heavy, by about 25% for YS-PBE0). No problems were observed for DFT+U. With SOC, the matrix elements are taken between spinor states, but the spin-orbit correction to the velocity operator itself is not included.
- **VASP** (`WAVEDER`) provides derivative overlaps that include the nonlocal PAW terms, so the matrix elements in the sum are velocity matrix elements. The conventions are documented in `read_mommat_pij_vasp_mod.f90`. The phase of the matrix elements does not affect the effective masses.
- **Free-electron term.** The VASP PAW Hamiltonian contains nonlocal projectors that depend on k. The velocity matrix elements in `WAVEDER` include their first derivative. The second derivative $\langle u_n | \partial^2 \hat{V}_\mathrm{NL} / \partial k_\alpha \partial k_\beta | u_n \rangle$ and the related terms of the PAW overlap operator are not available and are approximated by $\delta_{\alpha\beta}$ [see Phys. Rev. B 93, 205147 (2016)]. Effective masses obtained with VASP have agreed well with WIEN2k results in our comparisons, which suggests that this approximation is minor for the systems studied so far.
- **Truncated band sum.** Only the bands available in the input enter the sums (see [Convergence](#convergence)). Core states are not included.

## Degenerate bands

### Degenerate perturbation theory

For a block D of degenerate bands, mstar constructs for every tensor component αβ the Hermitian matrix [Eq. (6) of the reference]

```math
M^{(\alpha\beta)}_{nm} = \frac{1}{2m_0} \sum_{l \notin D}
\left( p^{(\alpha)}_{nl}\, p^{(\beta)}_{lm} + p^{(\beta)}_{nl}\, p^{(\alpha)}_{lm} \right)
\left[ \frac{1}{E_n - E_l} + \frac{1}{E_m - E_l} \right] ,
\qquad n, m \in D ,
```

and the inverse masses of the block are $\delta_{\alpha\beta} + \mathrm{eig}\, M^{(\alpha\beta)}$ [Eq. (7)], listed in ascending order. For a single band this is the non-degenerate expression. For exactly degenerate states, the averaged inverse energy difference equals the expression of the reference. The sums are accumulated in double precision with compensated (Kahan) summation.

The theory assumes that the degeneracy is not lifted at first order in k, i.e., that the bands of the block have zero velocity at k within the block. At points where bands cross linearly (e.g., Weyl or multifold fermions, such as CoSi at Γ), the dispersion is not quadratic and the masses of these blocks do not describe it.

### Grouping and tolerance

Bands are grouped into blocks at every k-point by chaining adjacent bands whose energies differ by no more than the tolerance (as in BerryCPT). Eigenvalues of states that are degenerate by symmetry differ slightly for numerical reasons, and treating them as non-degenerate leads to large errors through the energy denominators. Inspect the eigenvalues of the states that should be degenerate (e.g., the top valence bands of GaAs at Γ) and set the tolerance slightly above their spread. The approach is stable as long as the block is well separated from other bands compared with the tolerance.

### Band warping and the `# non-commuting block` flag

The six matrices $M^{(\alpha\beta)}$ of a block generally do not commute, so their eigenvalues belong to different combinations of the degenerate states. A record of `minv_ij.dat` then combines values that belong to different states, and a mass tensor cannot be assigned to an individual band (band warping, e.g., heavy and light holes) [Phys. Rev. B 93, 205147 (2016), Sec. III].

For every block, mstar evaluates

```math
r = \max_p \frac{\left\| \mathrm{offdiag}\left( U^\dagger M_p U \right) \right\|_F}{\max_p \| M_p \|_F} ,
```

where $M_p$ are the six matrices of the block and $U$ diagonalizes a generic linear combination of them. If the matrices commute (e.g., for a Kramers pair that remains degenerate around k), r vanishes up to rounding. Otherwise r estimates the relative error of a per-band tensor. Blocks with r > 1e-3 are marked `# non-commuting block`, and `minv_ij.dat` also gives r. For records of marked blocks:

- the diagonal components xx, yy, zz in `minv_ij.dat` are the inverse masses along the Cartesian axes (direction-resolved masses, labelled by ascending order),
- the off-diagonal components yz, xz, xy and the values in `minv_c.dat`, `minv_pr.dat`, and `minv_d.dat` are not defined for individual bands.

## Output

Output files (with the suffix `-up`/`-dn` for spin-polarized calculations):

| File | Content |
|---|---|
| `minv_ij.dat` | components of the inverse effective mass tensor m0/m*_ij (xx, yy, zz, yz, xz, xy) |
| `minv_c.dat` | conductivity inverse effective mass m0/m*_c = (1/3) Tr(m0/m*_ij) |
| `minv_pr.dat` | principal components (eigenvalues) of m0/m*_ij |
| `minv_d.dat` | density-of-states inverse effective mass m0/m*_d = m0/(m_1 m_2 m_3)^(1/3), signed |

All values are dimensionless inverse masses in units of 1/m0. Each k-point block starts with

```text
# KP: k NBCDER: Nout NEMAX: Nbands
```

where k is the k-point index in the order of the input file, Nout is the last band written, and Nbands is the number of bands at that k-point. Each record starts with the band index and the ID of the degenerate block (numbered from 1 at each k-point), followed by the values. Annotations follow a `#` at the end of a record (`# non-commuting block`, and `# hyperbolic dispersion` in `minv_d.dat` when the principal components have different signs, where the geometric average is of limited use). The files can therefore be read directly by gnuplot or `numpy.loadtxt`. Every file has a header that describes its content.

Existing output files are replaced. The input is checked before the output files are opened, so a run that stops on unusable input does not overwrite the results of a previous run. All errors terminate with a nonzero exit status.

The script `utils/dat2csv.sh` combines band energies along a WIEN2k band-structure path with selected mass components into a CSV file.

## Tests

The folder `test_files` contains regression tests for WIEN2k and VASP (CoSi with spin-orbit coupling), including the DFT input files and the calculation details in the `etc` folders. After compiling, run

```bash
test_files/run_tests.sh
```

See `test_files/README.md` for details. The results agree with the reference outputs up to rounding in the last digits, which depends on the compiler and the LAPACK library (tested with gfortran and reference LAPACK, and with Intel Fortran and MKL).

## Reference

If you find the results useful and publishable, please cite:

- O. Rubel, F. Tran, X. Rocquefelte, and P. Blaha, "Perturbation approach to _ab initio_ effective mass calculations", [Comp. Phys. Commun. 261, 107648 (2021)](https://doi.org/10.1016/j.cpc.2020.107648)

together with the Git revision printed by mstar. Examples of applications can be found in the reference ([arXiv:2007.03816](https://arxiv.org/abs/2007.03816)). Please communicate feedback, support, or feature requests via the WIEN2k [mailing list](http://www.wien2k.at/reg_user/mailing_list).

## Changes

### November 2026

Corrections:

- VASP: `EIGENVAL` was read with a fixed VASP 5 format. For files written by VASP 6.4.3 (and possibly other VASP 6 versions), with a wider band-index field, the eigenvalues were truncated to 0.01 eV. This corrupted the energy denominators and could merge non-degenerate states into false degenerate groups. `EIGENVAL` is now read in free format.
- VASP: the couplings of band `NBANDS_CDER` (the last band with derivative matrix elements in `WAVEDER`) to all higher bands were omitted, which affected that band and the bands degenerate with it.
- Errors now terminate with a nonzero exit status.

Numerics:

- The band sums are accumulated in double precision with compensated summation, and all matrices are diagonalized in double precision. Values are written with 7 significant digits. The code no longer uses OpenMP, and the results do not depend on the summation order of a parallel run.
- For numerically near-degenerate states inside a block, the effective matrix uses the average of the inverse energy differences, which keeps it Hermitian. For exactly degenerate states the result is unchanged.

Degenerate bands and output:

- Each record gives the ID of the degenerate block after the band index. Scripts that read the output by column need to account for the extra column (`utils/dat2csv.sh` is updated).
- Blocks are formed by chaining adjacent bands within the tolerance (previously, all bands of a block had to lie within the tolerance of its first band).
- Blocks whose tensor components do not commute are marked `# non-commuting block` (see [Degenerate bands](#degenerate-bands)).
- Annotations are written as trailing `#` comments (`warning: hyperbolic dispersion` became `# hyperbolic dispersion`).
- The Git revision of the build is written to all output files and to the screen.

Input and interface:

- New option `--nbmax N`.
- VASP and WIEN2k input is read with the routines of [BerryCPT](https://github.com/rubel75/berrycpt). Only the rows of the target bands are stored. WIEN2k files are read in a single pass. The VASP eigenvalues are kept in double precision.
- The input file may be given with a path. For VASP, `EIGENVAL` is read from the directory of `WAVEDER`.
- The input is checked before the output files are opened.
- Regression tests in `test_files`.

### Earlier versions

- October 13, 2021
- March 24, 2021 (gfortran compatibility)
