PROGRAM mstar
! includes modules and subprograms:
!     degenmass_mod.f90 -- Effective-matrix inverse mass for a block of
!       degenerate bands (double precision, Kahan summation)
!     dgenen_mod.f90 -- Group degenerate states
!     eigvs_mod.f90 -- Eigenvalues of a real symmetric matrix
!     kdelta.f90 -- Kronecker delta function
!     read_mommat_nb.f90 -- Determine number of bands in a k-point block of
!       mommat file
!     read_mommat_pij.f90 -- Read momentum matrix elements <i|p_a|j> (a=x,y,z)
!       and energy differences E_i - E_j in a k-point block of mommat file
!     read_numlines.f90 -- Read number of lines in a file
!   shared with BerryCPT (https://github.com/rubel75/berrycpt):
!     precision_mod.f90, eigvz_mod.f90, is_hermitian_mod.f90,
!     read_waveder_header_mod.f90, read_eigenvalues_vasp_mod.f90,
!     read_mommat_pij_vasp_mod.f90
!
! (c) Oleg Rubel, Mar 2021
!
! Execution:
!   $ ./mstar arg1 [arg2]
! arg1 - input mommat file name
! arg2 - (optional) degeneracy energy tolerance [Ha]
!        max dE between adjacent states to be grouped as degenerate
!        default value is 1.0e-6 Ha
!
! Output:
! minv_ij.dat - contains elements of the (m0/m*_ij) tensor.               
! minv_c.dat - conductivity effective mass (m0/m*_c).
! minv_pr.dat - contains principal components of the inverse eff. mass 
!               tensor eig(m0/m*_ij)
! minv_d.dat - density of states inverse effective mass 
!              m0/m*_d = m0/(m_1*m_2*m_3)**(1/3)
! If the files exist from a previous run, they will be removed
!
! Tips:
! (1): Writing of the mommat file is _not_ default in WIEN2k.
!      To enable writing, edit the case.inop file and change
!      OFF to ON in the following line:
!      ON           ON/OFF   writes MME to unit 4
!      -^
! (2): Make sure to get _plenty_ of empty bands during SCF.
!      This requires modification in several input files:
!      (a) extend "de" in case.in1(c) above 5 Ry
!          K-VECTORS FROM UNIT:4   -9.0      10.0    10   emin / de (emax=Ef+de) / nband
!          -----------------------------------^
!      (b) if you do SOC calculation, extend "Emax" in
!          case.inso up to 5 Ry
!          -10 5.0                Emin, Emax
!          -----^
!      (c) default "Emax" in case.inop is 3 Ry, which should
!          be OK, but it can be good to test the convergence
!          and push this parameter up to 5 Ry
!          -5.0 3.5 9999 Emin, Emax for matrix elements, NBvalMAX
!          ------^
! (3): In VASP calculations use LOPTICS = .TRUE., increase (at least x3) 
!      the number NBANDS = XXXX, and disable a finite differences derivative 
!      of the cell-periodic part of the orbitals LPEAD =.FALSE. Tests show 
!      that m* values calculated with LPEAD =.TRUE. make no sense.


!! Variables

USE precision_mod, ONLY: sp, dp
USE dgenen_mod, ONLY: dgenen
USE degenmass_mod, ONLY: degenmass, block_inconsistency
USE eigvs_mod, ONLY: eigvs
USE read_waveder_header_mod, ONLY: read_waveder_header
USE read_eigenvalues_vasp_mod, ONLY: read_eigenvalues_vasp
USE read_mommat_pij_vasp_mod, ONLY: read_mommat_pij_vasp

implicit none
CHARACTER(len=256) :: &
    arg1, arg2, & ! command line input arguments
    fnameinp, fnameout2, &
    fnameout3, fnameout4, fnameout5, & ! input/output file names
    wformat2, wformat3, &
    wformat4, wformat5, & ! format for writing/reading data
    charspin ! spin component for spin-polarized calculations
INTEGER :: &
    nltot, & ! total number of lines in mommat file
    nktot, & ! total number of k-points in WAVEDER file
    nstot, & ! total number of spins in WAVEDER file
    iline, & ! current line number during reading of mommat file
    ispin, & ! current spin (1/2)
    nb, & ! number of bands for a current k-point
    nbb, & ! number of band-to-band transition
    nbcder, & ! the max number of bands for which m* is calculated
              ! In WIEN2k we set it = to the total number of bands
              ! In VASP it is set in WAVEDER file as NBANDS_CDER
    i, ikpt, & ! counters
    ierr ! error code
INTEGER, ALLOCATABLE :: &
    nbk(:,:) ! number of bands at each k-point and spin (VASP)
INTEGER, EXTERNAL :: &
    kdelta ! Kronecker delta function
REAL(kind=sp) :: &
    dEtol ! max dE between adjacent states to be grouped as degenerate [Ha]
REAL(kind=sp), ALLOCATABLE :: &
    dEij(:,:), & ! energy differences E_j-E_i [Ha] (WIEN2k, one k-point)
    dEijks(:,:,:,:), & ! k- and spin-dependent E_j-E_i [Ha] (VASP)
    occup(:,:,:) ! occupations (VASP, not used)
REAL(kind=dp), ALLOCATABLE :: &
    ene(:,:,:) ! eigenvalues [Ha] (VASP)
COMPLEX(kind=sp), ALLOCATABLE :: &
    pij(:,:,:), & ! momentum matrix elements [at.u.] (WIEN2k, one k-point)
    pijks(:,:,:,:,:) ! k- and spin-dependent velocity matrix elements
                     ! [at.u.] (VASP), rectangular (3,nbcder,nb,nktot,nstot)
INTEGER :: &
    nflagged = 0 ! number of non-commuting degenerate blocks
REAL(kind=dp), PARAMETER :: &
    block_tol = 1.0e-3_dp ! threshold on block_inconsistency for flagging
LOGICAL :: &
    fmommatend, & ! end of mommat file
    file_exists, &
    wien2k, & ! true if this is a WIEN2k calculation
    warnHyperb ! warning on hyperbolic band dispersion

!! Get command line input arguments

DO i = 1, iargc() ! loop through all command line arguments to search for help
    CALL GETARG(i,arg1)
    IF ( TRIM(arg1)=='-h' .or.  TRIM(arg1)=='--h' .or. & !...
            TRIM(arg1)=='-help' .or. TRIM(arg1)=='--help') THEN
        GOTO 911 ! print help and STOP
    END IF
END DO
write(*,'(A,I0)') ' Detected input arguments = ', iargc()
IF (iargc() == 1) THEN ! check number of input arguments
    CALL GETARG(1,arg1) ! mommat file name
    write (*,*) 'Input mommat file = ', TRIM(arg1)
    fnameinp = TRIM(arg1)
    dEtol = 1.0e-6 ! default degeneracy tolerance [Ha]
    write (*,*) 'Degeneracy tolerance is not specified,'
    write (*,'(A,es12.5,A)') '  default will be used Etol = ', dEtol, ' [Ha]'
ELSEIF (iargc() == 2) THEN ! for 2 arguments
    CALL GETARG(1,arg1) ! mommat file name
    write (*,*) 'Input mommat file = ', TRIM(arg1)
    fnameinp = TRIM(arg1)
    CALL GETARG(2,arg2) ! degeneracy tolerance [Ha]
    write (*,*) 'Degeneracy tolerance dEtol = ', TRIM(arg2), ' [Ha]'
    read(arg2,*,IOSTAT=ierr) dEtol
    IF (ierr /= 0) THEN ! input error for 2nd argument
        write(*,*) 'Error detected for the 2nd input argument (must be a number)'
        write(*,*) 'Suggested execution:'
        write(*,*) '$ ./mstar mass.mommat2up 1e-6'
        write(*,*) 'or'
        write(*,*) '$ ./mstar mass.mommat2up'
        ERROR STOP 1
    END IF
    write (*,'(A,es12.5,A)') & !...
        ' Confirming text-to-number conversion dEtol = ', dEtol, ' [Ha]'
ELSE ! impossible
    GOTO 912 ! print error and STOP
END IF

!! WIEN2k or VASP?

IF ( TRIM(fnameinp) == 'WAVEDER' .OR. TRIM(fnameinp) == 'WAVEDERF' ) THEN ! VASP
    wien2k = .false.
    write (*,*) 'Assume VASP calculation'
    ! change to a binary file in case it is pointed at the formatted
    ! file WAVEDER
    fnameinp = 'WAVEDER'
    write(*,*) 'Assumed VASP calculation based on the input file name.'
ELSE ! WIEN2k (default)
    wien2k = .true.
    write(*,*) 'Assumed WIEN2k calculation based on the input file name.'
    write(*,*) '(If you would like to read VASP file, the input file should'
    write(*,*) 'have the exact name WAVEDER.)'
END IF

!! Determine number of lines in mommat file & read VASP input

! check if the file exists
INQUIRE(FILE=fnameinp, EXIST=file_exists)
IF ( file_exists ) THEN
    write(*,*) 'The input file ', TRIM(fnameinp), ' was found.'
ELSE IF ( .not.(file_exists) ) THEN
    write(*,*) 'The input file ', TRIM(fnameinp), ' does not exist. Exiting'
    ERROR STOP 1
END IF
IF ( file_exists .AND. (.not.(wien2k)) ) THEN ! check EIGENVAL for VASP
    INQUIRE(FILE='EIGENVAL', EXIST=file_exists)
    IF ( .not.(file_exists) ) THEN
        write(*,*) 'The file EIGENVAL is also required, but it ', &
            'does not exist. Exiting'
        ERROR STOP 1
    ELSE
        write(*,*) 'The file EIGENVAL is also required. ', &
            'It is found.'
    END IF
END IF
IF (wien2k) THEN
    CALL read_numlines(fnameinp, 1, & ! <- args in
                nltot) ! -> args out
    write (*,'(A,I0)') '  number of lines in mommat file = ', nltot
    IF (nltot < 10) THEN
        write(*,*) 'The file ', TRIM(fnameinp), ' is too short ', & !...
            'and most likely useless. Stopping'
        ERROR STOP 1
    END IF
    ! assume 1 spin since case.momat2up and case.momat2dn files should be
    ! read one after the other any way
    nstot = 1
    ! open input mommat file for reading (file ID=1)
    OPEN (1, file = TRIM(fnameinp), status = 'old')
ELSE ! VASP
    ! WAVEDER header: NB_TOT, NBANDS_CDER, NKPTS, ISPIN
    CALL read_waveder_header( &
        fnameinp, &               ! in
        nstot, nktot, nbcder, nb) ! out
    ! eigenvalues [Ha] from EIGENVAL (occupations are not used)
    CALL read_eigenvalues_vasp( &
        'EIGENVAL', &             ! in
        nstot, nktot, nb, &       ! in
        nbk, ene, occup)          ! out
    ! velocity matrix elements <i|v|j> [a.u.] and E_j-E_i [Ha] in
    ! rectangular storage (3,nbcder,nb,nktot,nstot) and
    ! (nbcder,nb,nktot,nstot). The VASP phase convention is documented in
    ! read_mommat_pij_vasp_mod.f90. Only products of the form
    ! p_a(i,l) p_b(l,j) enter the effective mass, and for bands i, j of one
    ! contiguous degenerate block (l outside the block) the phase factors of
    ! the two elements cancel. The effective masses therefore do not depend
    ! on that convention.
    CALL read_mommat_pij_vasp( &
        fnameinp, nstot, nktot, nbcder, nb, ene, & ! in
        pijks, dEijks)                             ! out
    DEALLOCATE( ene, occup, nbk )
END IF

! Loop over spins (for VASP only)
! In WIEN2k nstot=1 since case.momat2up and case.momat2dn files should be
! read one after the other any way
warnHyperb = .false.
write(*,*) 'Entering the main loop...'
DO ispin = 1, nstot

    !! Spin suffix

    IF (wien2k) THEN
        ! take last two letters of the input file name
        charspin = fnameinp( LEN(TRIM(fnameinp))-1 : LEN(TRIM(fnameinp)) )
        IF ( TRIM(charspin) == 'up' ) THEN
            charspin = '-up'
        ELSE IF ( TRIM(charspin) == 'dn' ) THEN
            charspin = '-dn'
        ELSE
            charspin = '' ! no spin identity
        END IF
    ELSE ! VASP
        IF (ispin==1 .AND. nstot==2) THEN
            charspin = '-up'
        ELSE IF (ispin==2 .AND. nstot==2) THEN
            charspin = '-dn'
        ELSE
            charspin = '' ! no spin identity
        END IF
    END IF

    !! Prepare output files

    fnameout2 = 'minv_ij'//TRIM(charspin)//'.dat'
    OPEN (2, file = TRIM(fnameout2), status = 'UNKNOWN') ! output file 1
    write(2,'(A33)') '# This file is generated by mstar'
    write(2,'(A62)') '# the output contains inverse effective masses (m0/m_ij*) that'
    write(2,'(A57)') '# are grouped by k-point index and then by the band index'
    write(2,'(A)') '# Each band record also gives the ID of its degenerate block (numbered'
    write(2,'(A)') '# from 1 at each k-point). Bands with the same block ID were treated'
    write(2,'(A)') '# together by degenerate perturbation theory (Eqs. 6 and 7 of the'
    write(2,'(A)') '# reference); adjacent bands belong to one block when their energies'
    write(2,'(A)') '# differ by no more than the degeneracy tolerance.'
    write(2,'(A)') '# A block is marked "# non-commuting block" when the effective matrices'
    write(2,'(A)') '# of its tensor components (Eq. 6) do not share eigenvectors (band'
    write(2,'(A)') '# warping). Each component is diagonalized separately and its eigenvalues'
    write(2,'(A)') '# are listed in ascending order, so a record then combines values that'
    write(2,'(A)') '# belong to different states [Phys. Rev. B 93, 205147 (2016), Sec. III].'
    write(2,'(A)') '# In such blocks, xx, yy, zz are the inverse masses along the Cartesian'
    write(2,'(A)') '# axes, while yz, xz, xy are not defined for individual bands. r is the'
    write(2,'(A)') '# estimated relative error of a per-band tensor (blocks are marked if'
    write(2,'(A)') '# r > 1e-3).'
    write(2,'(A53)') '# columns correspond to Cartesian directions for m_ij'
    write(2,'(A)') '# band block 1=xx;     2=yy;      3=zz;     4=yz;'//&!...
        '       5=xz;      6=xy'

    fnameout3 = 'minv_c'//TRIM(charspin)//'.dat'
    OPEN (3, file = TRIM(fnameout3), status = 'UNKNOWN') ! output file 2
    write(3,'(A)') '# This file is generated by mstar'
    write(3,'(A)') '# the output contains inverse conductivity eff. masses'
    write(3,'(A)') '# m0/m_c = 1/3*Tr(m0/m*_ij) that'
    write(3,'(A)') '# are grouped by k-point index and then by the band index'
    write(3,'(A)') '# Each band record also gives the ID of its degenerate block (numbered'
    write(3,'(A)') '# from 1 at each k-point). Bands with the same block ID were treated'
    write(3,'(A)') '# together by degenerate perturbation theory (Eqs. 6 and 7 of the'
    write(3,'(A)') '# reference); adjacent bands belong to one block when their energies'
    write(3,'(A)') '# differ by no more than the degeneracy tolerance.'
    write(3,'(A)') '# A block is marked "# non-commuting block" when the effective matrices'
    write(3,'(A)') '# of its tensor components (Eq. 6) do not share eigenvectors (band'
    write(3,'(A)') '# warping). Each component is diagonalized separately and its eigenvalues'
    write(3,'(A)') '# are listed in ascending order, so a record then combines values that'
    write(3,'(A)') '# belong to different states [Phys. Rev. B 93, 205147 (2016), Sec. III].'
    write(3,'(A)') '# The values in such blocks are not defined for individual bands.'
    write(3,'(A)') '# columns correspond to'
    write(3,'(A)') '# band block m0/m_c'

    fnameout4 = 'minv_pr'//TRIM(charspin)//'.dat'
    OPEN (4, file = TRIM(fnameout4), status = 'UNKNOWN') ! output file 3
    write(4,'(A)') '# This file is generated by mstar'
    write(4,'(A)') '# the output contains principal components of the inverse'
    write(4,'(A)') '# eff. mass tensor eig(m0/m*_ij) that'
    write(4,'(A)') '# are grouped by k-point index and then by the band index'
    write(4,'(A)') '# Each band record also gives the ID of its degenerate block (numbered'
    write(4,'(A)') '# from 1 at each k-point). Bands with the same block ID were treated'
    write(4,'(A)') '# together by degenerate perturbation theory (Eqs. 6 and 7 of the'
    write(4,'(A)') '# reference); adjacent bands belong to one block when their energies'
    write(4,'(A)') '# differ by no more than the degeneracy tolerance.'
    write(4,'(A)') '# A block is marked "# non-commuting block" when the effective matrices'
    write(4,'(A)') '# of its tensor components (Eq. 6) do not share eigenvectors (band'
    write(4,'(A)') '# warping). Each component is diagonalized separately and its eigenvalues'
    write(4,'(A)') '# are listed in ascending order, so a record then combines values that'
    write(4,'(A)') '# belong to different states [Phys. Rev. B 93, 205147 (2016), Sec. III].'
    write(4,'(A)') '# The values in such blocks are not defined for individual bands.'
    write(4,'(A)') '# columns correspond to'
    write(4,'(A)') '# band block m0/m_1  m0/m_2  m0/m_3'

    fnameout5 = 'minv_d'//TRIM(charspin)//'.dat'
    OPEN (5, file = TRIM(fnameout5), status = 'UNKNOWN') ! output file 4
    write(5,'(A)') '# This file is generated by mstar'
    write(5,'(A)') '# the output contains density of states inverse effective mass'
    write(5,'(A)') '# m0/m_d = m0/(m_1*m_2*m_3)**(1/3) that'
    write(5,'(A)') '# are grouped by k-point index and then by the band index. Here'
    write(5,'(A)') '# m0/m_1  m0/m_2  m0/m_3 are principal components of the inverse'
    write(5,'(A)') '# eff. mass tensor eig(m0/m*_ij).'
    write(5,'(A)') '# Each band record also gives the ID of its degenerate block (numbered'
    write(5,'(A)') '# from 1 at each k-point). Bands with the same block ID were treated'
    write(5,'(A)') '# together by degenerate perturbation theory (Eqs. 6 and 7 of the'
    write(5,'(A)') '# reference); adjacent bands belong to one block when their energies'
    write(5,'(A)') '# differ by no more than the degeneracy tolerance.'
    write(5,'(A)') '# A block is marked "# non-commuting block" when the effective matrices'
    write(5,'(A)') '# of its tensor components (Eq. 6) do not share eigenvectors (band'
    write(5,'(A)') '# warping). Each component is diagonalized separately and its eigenvalues'
    write(5,'(A)') '# are listed in ascending order, so a record then combines values that'
    write(5,'(A)') '# belong to different states [Phys. Rev. B 93, 205147 (2016), Sec. III].'
    write(5,'(A)') '# The values in such blocks are not defined for individual bands.'
    write(5,'(A)') '# Records with principal components of different signs are marked'
    write(5,'(A)') '# "# hyperbolic dispersion", where the geometric average is of limited'
    write(5,'(A)') '# validity.'
    write(5,'(A)') '# Columns correspond to'
    write(5,'(A)') '# band block m0/m_d'

    !! Main part nested loops

    fmommatend = .false. ! end of mommat file is not reached
    ikpt = 0 ! initialize the counter for k-points
    nbb = 0 ! initialize the number of band-to-band transitions
    iline = 0 ! initialize the number of lines to skip

    DO WHILE (.not.(fmommatend)) ! loop over k-points until the file ends
        ikpt = ikpt + 1 ! count number of k-points

        IF (wien2k) THEN

            !! Determine number of bands in mommat file

            CALL read_mommat_nb(1, & ! <- args in
                    iline, & ! <-> args in-out
                    nb) ! -> args out
            ! the max number of bands for which m* is calculated
            ! here (WIEN2k) we set it = to the total number of bands
            nbcder = nb

            !! Read <b_i|p_a|b_j> from mommat file, a=1,2,3 (x,y,z)

            ALLOCATE( pij(3,nb,nb), dEij(nb,nb) )
            nbb = (nb+nb**2)/2 ! number of band-to-band transitions
            CALL read_mommat_pij (1, nb, nbb, & ! <- args in
                    iline, & ! <-> args in-out
                    pij, dEij) ! -> args out
            CALL process_kpoint(pij, dEij)
            DEALLOCATE( pij, dEij )
        ELSE ! VASP
            CALL process_kpoint(pijks(:,:,:,ikpt,ispin), &
                dEijks(:,:,ikpt,ispin))
        END IF

        !! Output progress to screen

        IF (wien2k) THEN
            write(*,'(A,I0,A,I0,A,I3,A)') & !...
                ' KP: ', ikpt, ' bands: ', nb, & !...
                ' progress: ', INT(100*iline/nltot), '%'
        ELSE ! VASP
            IF (ispin==1 .AND. nstot==2) THEN
                write(*,'(A,I0,A,I0,A,I0,A,I3,A)') & !...
                    ' KP: ', ikpt, ' spin: ', ispin, ' bands: ', nb, & !...
                    ' progress: ', INT(100*ikpt/nktot/2), '%'
            ELSE IF (ispin==2 .AND. nstot==2) THEN
                write(*,'(A,I0,A,I0,A,I0,A,I3,A)') & !...
                    ' KP: ', ikpt, ' spin: ', ispin, ' bands: ', nb, & !...
                    ' progress: ', INT(100*ikpt/nktot/2 + 50), '%'
            ELSE
                write(*,'(A,I0,A,I0,A,I3,A)') & !...
                    ' KP: ', ikpt, ' bands: ', nb, & !...
                    ' progress: ', INT(100*ikpt/nktot), '%'
            END IF
        END IF

        !! Check if the end of file mommat is reached

        IF (wien2k) THEN
            IF (iline == nltot) THEN ! end of mommat file, exit WHILE loop
                fmommatend = .true.
            END IF
        ELSE ! VASP
            IF (ikpt == nktot) THEN ! end of WAVEDER file, exit WHILE loop
                fmommatend = .true.
            END IF
        END IF

    END DO ! loop over k-points

    !! Close output files

    CLOSE (2) ! output file 1
    CLOSE (3) ! output file 2
    CLOSE (4) ! output file 3
    CLOSE (5) ! output file 4

END DO ! loop over spins

IF (wien2k) CLOSE (1) ! input  file 1
IF (ALLOCATED(pijks)) DEALLOCATE( pijks, dEijks )

write(*,*) 'Summary of the output:'
write(*,*) '(1) Components of the (m0/m*_ij) tensor are stored in file ', &!...
    TRIM(fnameout2)
write(*,*) '(2) The conductivity inverse eff. masses m0/m_c = 1/3*Tr(m0/m*_ij)'
write(*,*) '    are stored in file ', TRIM(fnameout3)
write(*,*) '(3) Principal components of the (m0/m*_ij) tensor'
write(*,*) '    are stored in file ', TRIM(fnameout4)
write(*,*) '(4) Density of states inverse effective mass'
write(*,*) '    m0/m*_d = m0/(m_1*m_2*m_3)**(1/3) are stored in file ', &!...
    TRIM(fnameout5)
IF (warnHyperb) THEN
    write(*,*) '    Some records are marked "# hyperbolic dispersion" in this file.'
    write(*,*) '    This indicates a limited validity of the geometric average'
    write(*,*) '    (m_1*m_2*m_3)**(1/3) for the DOS effective mass in the situations when'
    write(*,*) '    the band dispersion is _not_ parabolic.'
END IF
IF (nflagged > 0) THEN
    write(*,'(A,I0,A)') '  Records of ', nflagged, ' degenerate blocks are marked'
    write(*,*) '    "# non-commuting block". In these blocks, only the diagonal'
    write(*,*) '    components xx, yy, zz in minv_ij are defined for individual bands'
    write(*,*) '    (directional masses along the Cartesian axes). The off-diagonal'
    write(*,*) '    components and the values in minv_c, minv_pr, and minv_d are not.'
END IF
write(*,*) 'See the file header for the description'
write(*,*) 'Suggested reference:'
write(*,*) '[1] O. Rubel, F. Tran, X. Rocquefelte, and P. Blaha "Perturbation'
write(*,*) '    approach to ab initio effective mass calculations"'
write(*,*) '    Comp. Phys. Commun. 261, 107648 (2021).'
write(*,*) '    https://doi.org/10.1016/j.cpc.2020.107648'
STOP ! end of the main code

!! Help section

911 & ! label for GOTO statement
write(*,*) 'Execution:'
write(*,*) '  $ ./mstar arg1 [arg2]'
write(*,*) 'arg1 - input case.mommat2 (WIEN2k) or WAVEDER (VASP) file name'
write(*,*) 'arg2 - (optional) degeneracy energy tolerance [Ha]'
write(*,*) '       max dE between adjacent states to be grouped as degenerate'
write(*,*) '       default value is 1.0e-6 Ha'
write(*,*) ''
write(*,*) 'Output:'
write(*,*) 'minv_ij.dat - contains elements of the (m0/m*_ij) tensor.'
write(*,*) 'minv_c.dat - conductivity effective mass (m0/m*_c).'
write(*,*) 'minv_pr.dat - contains principal components of the inverse'
write(*,*) '              eff. mass tensor eig(m0/m*_ij)'
write(*,*) 'minv_d.dat - density of states effective mass (m0/m*_d).'
write(*,*) 'Each data record starts with the band index and the ID of the'
write(*,*) 'degenerate block (numbered from 1 at each k-point) to which the band belongs.'
write(*,*) 'If the files exist from a previous run, they will be removed'
write(*,*) ''
write(*,*) 'Tips:'
write(*,*) '(1): Writing of the mommat file is _not_ default in WIEN2k.'
write(*,*) '     To enable writing, edit the case.inop file and change'
write(*,*) '     OFF to ON in the following line:'
write(*,*) '     ON           ON/OFF   writes MME to unit 4'
write(*,*) '     -^'
write(*,*) '(2): Make sure to get _plenty_ of empty bands during SCF.'
write(*,*) '     This requires modification in several input files:'
write(*,*) '     (a) extend "de" in case.in1(c) above 5 Ry'
write(*,'(A)') '          K-VECTORS FROM UNIT:4   -9.0      10.0'//&!...
    '    10   emin / de (emax=Ef+de) / nband'
write(*,*) '         -----------------------------------^'
write(*,*) '     (b) if you do SOC calculation, extend "Emax" in'
write(*,*) '         case.inso up to 5 Ry'
write(*,*) '         -10 5.0                Emin, Emax'
write(*,*) '         -----^'
write(*,*) '     (c) default "Emax" in case.inop is 3 Ry, which should'
write(*,*) '         be OK, but it can be good to test the convergence'
write(*,*) '         and push this parameter up to 5 Ry'
write(*,*) '         -5.0 3.5 9999 Emin, Emax for matrix elements, NBvalMAX'
write(*,*) '         ------^'
write(*,*) '(3): In VASP calculations use LOPTICS = .TRUE., increase'
write(*,*) '     (at least x3) the number NBANDS = XXXX, and disable a finite'
write(*,*) '     differences derivative of the cell-periodic part of'
write(*,*) '     the orbitals LPEAD =.FALSE. Tests show that m* values'
write(*,*) '     calculated with LPEAD =.TRUE. make no sense.'
STOP

!! Error section

912 & ! label for GOTO statement
write(*,*) 'Error detected for the number of input arguments.'
write(*,*) 'There should be 1 or 2 arguments'
write(*,*) 'Suggested execution (WIEN2k):'
write(*,*) '$ ./mstar mass.mommat2 1e-6'
write(*,*) 'or'
write(*,*) '$ ./mstar mass.mommat2'
write(*,*) ' '
write(*,*) 'Suggested execution (VASP):'
write(*,*) '$ ./mstar WAVEDER 1e-6'
write(*,*) 'or'
write(*,*) '$ ./mstar WAVEDER'
ERROR STOP 1

CONTAINS

SUBROUTINE process_kpoint(pij_k, dEij_k)
! Calculate and write inverse effective masses for one k-point and spin.
! Host variables used: nb, ikpt, dEtol, warnHyperb, and output units 2-5.
!
! pij_k(a,i,l) = <i|p_a|l> [a.u.] and dEij_k(i,l) = E_l - E_i [Ha] use
! rectangular storage (component, target band, all bands). Masses are
! calculated for the target bands 1..SIZE(pij_k,2) (= nbcder), truncated at
! the last complete degenerate group.

COMPLEX(kind=sp), INTENT(IN) :: pij_k(:,:,:)
REAL(kind=sp), INTENT(IN) :: dEij_k(:,:)

INTEGER :: &
    nbt, & ! number of target bands
    ndblen, & ! number of complete degenerate groups among target bands
    nbout, & ! last band for which masses are written
    id, nd, j, n, & ! counters, size of a degenerate block, band index
    ivoigt, & ! Voigt index (1..6)
    alpha, beta ! Cartesian directions 1,2,3 = x,y,z
INTEGER, ALLOCATABLE :: &
    nmlist(:,:) ! list of degeneracy indices, e.g. [1 1; 2 4; ...]
REAL(kind=dp) :: &
    minv_det, & ! determinant of the inverse eff. mass tensor
    minv_d, & ! density of states inverse eff. mass
    minv_U(3,3), & ! the upper triangular part of the m0/m*_ij tensor
    eigv_pr(3) ! principal components of m0/m*_ij
REAL(kind=dp), ALLOCATABLE :: &
    eigv(:), & ! eigenvalues of the effective matrix M
    minv(:,:) ! inverse masses of a degenerate block (band, Voigt index)
COMPLEX(kind=dp), ALLOCATABLE :: &
    Mall(:,:,:) ! effective matrices M^(ab) of a block (band, band, Voigt)
REAL(kind=dp) :: &
    rblock ! inconsistency of the block (see block_inconsistency)
CHARACTER(len=64) :: &
    flag_ij, flag ! annotations of records in non-commuting blocks
CHARACTER(len=16) :: cw ! width of the band-index field
! Cartesian pairs for Voigt indices 1=xx; 2=yy; 3=zz; 4=yz; 5=xz; 6=xy
INTEGER, PARAMETER :: &
    alpha_voigt(6) = [1, 2, 3, 2, 1, 1], &
    beta_voigt(6)  = [1, 2, 3, 3, 3, 2]

nbt = SIZE(pij_k,2)
IF (SIZE(pij_k,3) /= nb .OR. SIZE(dEij_k,1) /= nbt .OR. &
        SIZE(dEij_k,2) /= nb) THEN
    write(*,'(A,I0)') ' ERROR: inconsistent matrix-element arrays at k-point ', &
        ikpt
    ERROR STOP 1
END IF

!! Find degenerate states within the energy tolerance

ALLOCATE( nmlist(nbt,2) )
CALL dgenen(nb, nbt, dEij_k, dEtol, & ! <- args in
        nmlist, ndblen) ! -> args out
! NOTE: nbt (= nbcder) can cut through a degenerate group of states. Then
!       the entire group is excluded from the output, and the output stops
!       at nmlist(ndblen,2).
IF (ndblen > 0) THEN
    nbout = nmlist(ndblen,2)
ELSE
    nbout = 0
END IF

!! Write information about the current k-point

write(2,'(A,I0,1X,A,I0,1X,A,I0)') & !...
    '# KP: ', ikpt, 'NBCDER: ', nbout, 'NEMAX: ', nb
write(3,'(A,I0,1X,A,I0,1X,A,I0)') & !...
    '# KP: ', ikpt, 'NBCDER: ', nbout, 'NEMAX: ', nb
write(4,'(A,I0,1X,A,I0,1X,A,I0)') & !...
    '# KP: ', ikpt, 'NBCDER: ', nbout, 'NEMAX: ', nb
write(5,'(A,I0,1X,A,I0,1X,A,I0)') & !...
    '# KP: ', ikpt, 'NBCDER: ', nbout, 'NEMAX: ', nb

!! prepare output formats for effective masses

write(cw,'(I0)') nbt ! band index field as wide as the largest target band
write(cw,'(I0)') LEN_TRIM(cw)
! band index and block ID use the same field width
! a trailing annotation (comment starting with #) may follow the values
wformat2 = '(I' // TRIM(cw) // ',1X,I' // TRIM(cw) // ',6(1X,ES14.6E3),A)'
wformat3 = '(I' // TRIM(cw) // ',1X,I' // TRIM(cw) // ',1X,ES14.6E3,A)'
wformat4 = '(I' // TRIM(cw) // ',1X,I' // TRIM(cw) // ',3(1X,ES14.6E3),A)'
wformat5 = '(I' // TRIM(cw) // ',1X,I' // TRIM(cw) // ',1X,ES14.6E3,A)'

!! Loop through blocks of degenerate states

DO id = 1, ndblen
    ! number of degenerate states in the block
    nd = nmlist(id,2) - nmlist(id,1) + 1
    ! array to store inverse masses; size 6 is because of 6 Voigt indices
    ALLOCATE( minv(nd,6), Mall(nd,nd,6) )
    ! loop over Voigt indices: 1=xx; 2=yy; 3=zz; 4=yz; 5=xz; 6=xy
    DO ivoigt = 1, 6
        alpha = alpha_voigt(ivoigt)
        beta = beta_voigt(ivoigt)
        ! eigenvalues of the effective matrix M, Eqs. (6) and (7)
        CALL degenmass(nmlist(id,1), nmlist(id,2), alpha, beta, & ! in
            pij_k, dEij_k, & ! in
            eigv, Mall(:,:,ivoigt)) ! out
        ! store the inverse effective masses (+ 1 if alpha=beta)
        minv(:,ivoigt) = eigv + REAL(kdelta(alpha,beta), dp)
        DEALLOCATE( eigv )
    END DO ! loop over 'ivoigt'
    ! non-commuting (warped) degenerate block?
    rblock = block_inconsistency(Mall)
    IF (rblock > block_tol) THEN
        nflagged = nflagged + 1
        write(flag_ij,'(A,ES9.2)') '  # non-commuting block, r =', rblock
        flag = '  # non-commuting block'
    ELSE
        flag_ij = ''
        flag = ''
    END IF
    ! store the inverse effective masses
    DO j = 1, nd
        n = nmlist(id,1) + j - 1
        write(2,TRIM(wformat2)) n, id, minv(j,1:6), TRIM(flag_ij)
        ! the trace of the matrix does not change after diagonalization
        ! so the sum of eigenvalues will be just a trace
        write(3,TRIM(wformat3)) n, id, SUM(minv(j,1:3))/3.0_dp, TRIM(flag)
        ! determine principal components of the eff. mass tensor
        minv_U = 0.0_dp
        minv_U(1,1) = minv(j,1) ! the upper triangular part of minv tensor
        minv_U(2,2) = minv(j,2) ! 1=xx; 2=yy; 3=zz; 4=yz; 5=xz; 6=xy
        minv_U(3,3) = minv(j,3)
        minv_U(2,3) = minv(j,4)
        minv_U(1,3) = minv(j,5)
        minv_U(1,2) = minv(j,6)
        CALL eigvs(3, minv_U, & ! <- args in
                eigv_pr) ! -> args out
        write(4,TRIM(wformat4)) n, id, eigv_pr(1:3), TRIM(flag)
        ! DOS effective mass
        ! determinant of inv. eff. mass tensor
        minv_det = eigv_pr(1)*eigv_pr(2)*eigv_pr(3)
        minv_d = SIGN(ABS(minv_det)**(1.0_dp/3.0_dp), minv_det)
        IF (ALL(eigv_pr < 0.0_dp) .OR. ALL(eigv_pr > 0.0_dp)) THEN
            ! parabolic band dispersion
            write(5,TRIM(wformat5)) n, id, minv_d, TRIM(flag)
        ELSE
            ! hyperbolic band dispersion (signs of masses are different)
            IF (LEN_TRIM(flag) > 0) THEN
                write(5,TRIM(wformat5)) n, id, minv_d, &
                    TRIM(flag)//', hyperbolic dispersion'
            ELSE
                write(5,TRIM(wformat5)) n, id, minv_d, &
                    '  # hyperbolic dispersion'
            END IF
            warnHyperb = .true.
        END IF
    END DO
    DEALLOCATE( minv, Mall )
END DO ! loop over 'id'
DEALLOCATE( nmlist )

END SUBROUTINE process_kpoint

END PROGRAM mstar
