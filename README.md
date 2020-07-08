# mstar
Effective mass calculation with DFT using a perturbation theory. Currently supported codes:
* [WIEN2k](http://www.wien2k.at)
* [VASP](https://www.vasp.at)

It is written in Fortran and intended for Linux OS


### Installation:
First clone the GitHub repository

`$ git clone https://github.com/rubel75/mstar`

The `makefile` is set up for Intel Fortran compiler `ifort`. To compile, simply execute

`$ cd mstar; make`


### Execution
First, you need to perform a standard SCF calculation and generate a file that contains optical matrix elements (`case.mommat2[up/dn]` in WIEN2k or `WAVEDER` in VASP). Tips for this can be found at... Once the file is ready, execute

`x mstar [-up] [-settol 1.0e-5] # if you use the version built into WIEN2k starting with v20.1`

`/path/to/mstar case.mommat2[up/dn] [1e-5] # if you use this GitHub version and WIEN2k prior to v20.1 (see the compatibility note above)`

`/path/to/mstar WAVEDER [1e-5] # VASP`

Options:

  `[1e-5]` is (optional) degeneracy energy tolerance [Ha], which is max dE for 2 states to be considered as degenerate (default value is 1.0e-6 Ha).


### Output

`minv_ij.dat` - contains elements of the (m0/m*_ij) tensor for each k point and band.

`min_c.dat` - conductivity effective mass (m0/m*_c).

`minv_pr.dat` - contains principal components of the inverse eff. mass tensor eig(m0/m*_ij).

`minv_d.dat` - density of states inverse effective mass m0/m*_d = m0/(m_1 *m_2 *m_3)**(1/3).

The file headers explain the content.


For tutorials, please refer to ...

Examples of "real life" applications can be found in ...

Please communicate your feedback, support or feature requests via WIEN2k [mailing list](http://www.wien2k.at/reg_user/mailing_list)

### Reference

If you find the results useful and publishable, we will appreciate citing the following paper:

* O. Rubel, F. Tran, X. Rocquefelte, and P. Blaha "Perturbation approach to _ab initio_ effective mass calculations" [to appear on arXiv.org](https://arxiv.org).
