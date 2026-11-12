# mstar
Effective mass calculation with DFT using a perturbation theory. Currently supported codes:
* [WIEN2k](http://www.wien2k.at)
* [VASP](https://www.vasp.at)

It is written in Fortran and intended for Linux OS

### Current Version

Development version, November 2026 (changes since Oct 13, 2021):

* **VASP results obtained with earlier versions should be recalculated.** (1) The `EIGENVAL` file was read with a fixed format of VASP 5. With the wider band-index field written by VASP 6.4.3 (and possibly other VASP 6 versions), the eigenvalues were truncated to 0.01 eV. This corrupts the energy denominators and can merge non-degenerate states into false degenerate groups. `EIGENVAL` is now read in free format. (2) The couplings of band `NBANDS_CDER` (the last band with derivative matrix elements in `WAVEDER`) to all higher bands were omitted, which affected that band and the bands degenerate with it.
* The band sums are accumulated in double precision with compensated (Kahan) summation, and the effective matrices are diagonalized in double precision. The code no longer uses OpenMP, so it compiles without OpenMP support, and the results do not depend on the summation order of a parallel run. Values are written with 7 significant digits.
* For numerically near-degenerate states inside a degenerate group, the effective matrix uses the average of the inverse energy differences, which keeps it Hermitian. For exactly degenerate states the result is unchanged.
* VASP input is read with the routines of [BerryCPT](https://github.com/rubel75/berrycpt). Only the derivative-band block is stored, which reduces the memory needed for `WAVEDER` by more than half.
* Each record of the output files now gives the band index followed by the ID of the degenerate block to which the band belongs (numbered from 1 at each k-point). Scripts that read the output by column need to account for the extra column (`utils/dat2csv.sh` is updated).
* Degenerate blocks are formed by chaining adjacent bands whose energies differ by no more than the tolerance, as in BerryCPT. Previously, all bands of a block had to lie within the tolerance of the first band of the block.
* Records of degenerate blocks whose tensor components do not share eigenvectors (band warping, e.g., heavy and light holes) are marked `# non-commuting block`. For such bands, a mass tensor cannot be assigned to an individual band [Laflamme Janssen et al., Phys. Rev. B 93, 205147 (2016), Sec. III]: only the diagonal components xx, yy, zz in `minv_ij` (inverse masses along the Cartesian axes) are meaningful, while the off-diagonal components and the values in `minv_c`, `minv_pr`, and `minv_d` are not. In `minv_ij`, the flag also gives `r`, the estimated relative error of a per-band tensor; blocks are marked when r > 1e-3. Blocks with commuting components (e.g., Kramers pairs) are not marked.
* All annotations follow a `#` at the end of a record (e.g., `# hyperbolic dispersion` in `minv_d`, formerly `warning: hyperbolic dispersion`), so that the files can be read directly by gnuplot or `numpy.loadtxt`.
* New option `--nbmax N` limits the output to bands 1..N (in VASP, at most `NBANDS_CDER`). All bands still enter the perturbation sums, so the masses of the bands written are unchanged. A degenerate block cut by N is excluded. This avoids writing masses of high-lying bands, which are not converged with respect to the number of empty states.
* WIEN2k input is read with the routines of BerryCPT in a single pass, one k-point at a time. Only the rows of the target bands are stored, so `--nbmax` also reduces the memory. The input is checked before the output files are opened, so that the output of a previous run is not overwritten when the input is not usable.
* The input file may be given with a path. For VASP, `EIGENVAL` is read from the directory of `WAVEDER`.
* The Git revision of the build is written in the first line of every output file and to the screen (see `makefile`). Please report it together with results.
* Errors terminate with a nonzero exit status.
* Regression tests are in `test_files` (see `test_files/README.md`).

WIEN2k results are unchanged apart from rounding in the last digits.

Version Oct 13, 2021

Version Mar 24, 2021 (+gfortran compatibility)


### Installation:
First clone the GitHub repository

`$ git clone https://github.com/rubel75/mstar`

**WIEN2k compatibility note**:
The format of the `case.mommat2` file changed in WIEN2k v20.1 (free format, to allow more than 9999 bands). This version of `mstar` reads only the new format. Files written by WIEN2k versions _prior_ to v20.1 are not supported; use an older release of `mstar` for them.

The `makefile` is set up for Intel Fortran compiler `ifort` and Intel MKL. (To compile with `gfortran` you need to uncoment corresponding `FC`, `FCFLAGS`, and `FLFLAGS` variables in the `makefile`). To compile, simply execute

`$ cd mstar; make`


### Execution
First, you need to perform a standard SCF (self-consistent field) calculation and generate a file that contains optical matrix elements (`case.mommat2[up/dn]` in WIEN2k or `WAVEDER` in VASP). Tips on how to do this can be found at this [Wiki page](https://github.com/rubel75/mstar/wiki). Once the file is ready, execute

`x mstar [-up/-dn] [-settol 1.0e-5] # if you use the version built into WIEN2k starting with v20.1`

`/path/to/mstar/mstar case.mommat2[up/dn] [1e-5] # if you use this GitHub version and WIEN2k (see the compatibility note above)`

`/path/to/mstar/mstar WAVEDER [1e-5] # VASP`

The option `--nbmax N` (anywhere on the command line) limits the output to bands 1..N.

Options:

  * `[-up/-dn]` tells `mstar` to read `case.mommat2[up/dn]` files (needed for spin-polarized calculations or calculations with spin-orbit coupling).

  * `[-settol 1.0e-5]` or `[1e-5]` is (optional) degeneracy energy tolerance [Ha], which is the max dE between adjacent states for them to be grouped as degenerate (default value is 1.0e-6 Ha). Eigenvalues of states that are physically expected to be degenerate (for example light and heavy holes in GaAs at Gamma point) are not exactly the same for numerical reasons. If we treat them as non-degenerate states, a large error will occur due to the term p^2/dE when dE -> 0. These states need to be "bundled" together and treated as degenerate. It is recommended to inspect eigenvalues (in the case of GaAs it will be the top valence bands) and evaluate the energy difference dE between states that should be treated as degenerate. Then set the tolerance slightly greater than this value.


### Output

`minv_ij.dat` - contains elements of the (m0/m*_ij) tensor for each k point and band.

`minv_c.dat` - conductivity effective mass (m0/m*_c).

`minv_pr.dat` - contains principal components of the inverse eff. mass tensor eig(m0/m*_ij).

`minv_d.dat` - density of states inverse effective mass m0/m*_d = m0/(m_1 *m_2 *m_3)**(1/3).

The file headers explain the content.


### Documentation

For tutorials and tips on how to generate a file that contains optical matrix elements (`case.mommat2[up/dn]` in WIEN2k or `WAVEDER` in VASP), please refer to the [Wiki page](https://github.com/rubel75/mstar/wiki).

Examples of "real life" applications can be found in [arXiv:2007.03816 [cond-mat.mtrl-sci]](https://arxiv.org/abs/2007.03816)

Please communicate your feedback, support or feature requests via WIEN2k [mailing list](http://www.wien2k.at/reg_user/mailing_list)

### Reference

If you find the results useful and publishable, we will appreciate citing the following paper:

* O. Rubel, F. Tran, X. Rocquefelte, and P. Blaha "Perturbation approach to _ab initio_ effective mass calculations" [Comp. Phys. Commun. 261, 107648 (2021)](https://doi.org/10.1016/j.cpc.2020.107648).
