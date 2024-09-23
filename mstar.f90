PROGRAM mstar
! includes subprograms:
!     dgenen.f90 -- Group degenerate states
!     eigvz.f90 -- Solve a complex eigenvalue problem for a Hermitian matrix
!     kdelta.f90 -- Kronecker delta function
!     read_mommat_nb.f90 -- Determine number of bands in a k-point block of 
!       mommat file
!     read_mommat_pij.f90 -- Read momentum matrix elements <i|p_a|j> (a=x,y,z) 
!       and energy differences E_i - E_j in a k-point block of mommat file
!     read_numlines.f90 -- Read number of lines in a file
!
! (c) Oleg Rubel, Mar 2021
!
! Execution:
!   $ ./mstar arg1 [arg2]
! arg1 - input mommat file name
! arg2 - (optional) degeneracy energy tolerance [Ha]
!        max dE for 2 states to be considered as degenerate
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

USE OMP_LIB

implicit none
CHARACTER(len=256) :: &
    arg1, arg2, & ! command line input arguments
    fnameinp, fnameout2, &
    fnameout3, fnameout4, fnameout5, & ! input/output file names
    wformat2, wformat3, &
    wformat4, wformat5, wformat5w, & ! format for writing/reading data
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
    id, in, im, i, j, ikpt, & ! counters
    n, m, k, & ! band indices
    ndblen, & ! length of the useful part in 'nmlist' array
    nd, & ! number of degenerate states in the block
    ivoigt, & ! Voigt index (1..6)
    alpha, beta, & ! Cartesian directions 1,2,3 = x,y,z
    kdelta, & ! Kronecker delta function
    ierr ! error code
INTEGER, ALLOCATABLE :: &
    nmlist(:,:) ! list of degeneracy indices, e.g. [1 1; 2 4; ...]
REAL(kind=4) :: &
    dE, & ! energy difference [Ha]
    dEtol, & ! max dE for 2 states to be considered as degenerate [Ha]
    minv_det, & ! determinant of the inverse eff. mass tensor
    minv_d, & ! density of states eff. mass
    minv_U(3,3) ! the upper triangular part of the m0/m*_ij tensor
REAL(kind=4), ALLOCATABLE :: &
    dEij(:,:), & ! energy differences E_i-E_j [Ha]
    dEijks(:,:,:,:), & ! k- and spin-dependent energy differences E_i-E_j [Ha]
    EIGV(:), & ! eigenvalues of Mnm or 1/m*
    EIGFR(:,:), & ! real eigenvectors of 1/m*
    minv(:,:) ! array to store inverse masses
COMPLEX(kind=4) :: &
    dM, & ! contribution to the Mnm matrix
    p2 ! product of momentum matrix elements
COMPLEX(kind=4), ALLOCATABLE :: &
    pij(:,:,:), & ! momentum matrix elements [at.u.]
    pijks(:,:,:,:,:), & ! k- and spin-dependent momentum matrix elements [at.u.]
    Mnm(:,:), & ! momentum matrix elements
    EIGFC(:,:) ! complex eigenvectors of Mnm
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
        STOP
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

!! Determine number of lines in mommat file & read VASP WAVEDER header

! check if the file exists
INQUIRE(FILE=fnameinp, EXIST=file_exists)
IF ( file_exists ) THEN
    write(*,*) 'The input file ', TRIM(fnameinp), ' was found.'
ELSE IF ( .not.(file_exists) ) THEN
    write(*,*) 'The input file ', TRIM(fnameinp), ' does not exist. Exiting'
    STOP
END IF
IF ( file_exists .AND. (.not.(wien2k)) ) THEN ! check EIGENVAL for VASP
    INQUIRE(FILE='EIGENVAL', EXIST=file_exists)
    IF ( .not.(file_exists) ) THEN
        write(*,*) 'The file EIVENVAL is also required, but it ', &
            'does not exist. Exiting'
        STOP
    ELSE
        write(*,*) 'The file EIVENVAL is also required. ', &
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
        STOP
    END IF
    ! assume 1 spin since case.momat2up and case.momat2dn files should be
    ! read one after the other any way
    nstot = 1 
ELSE ! VASP
    CALL read_numlines_vasp(fnameinp, 1, & ! <- args in 
        nstot, nktot, nbcder, nb) ! -> args out
    write (*,'(A,I0)') '  number of spins in WAVEDER file = ', nstot
    write (*,'(A,I0)') '  number of k-points in WAVEDER file = ', nktot
    write (*,'(A,I0)') '  smaller number of bands in WAVEDER file = ', nbcder
    write (*,'(A,I0)') '  number of bands per k-point in EIGENVAL file = ', &
        nb
    ! Memory estimate: 
    ! 4 bytes * 2 (complex) * 3 (3D x,y,z) * 10e9 (Bytes -> GB)
    write (*,'(A,1X,F5.1,1X,A)') ' Memory required to store matrix elements'//&
        ' from WAVEDER file', nstot*nktot*(nb**2.)*4*2*3/10.0**9., 'GB'
    write (*,'(A,1X,F5.1,1X,A)') ' Memory required to store dE_ij'//&
        ' from EIGENVAL file', nstot*nktot*(nb**2.)*4*1/10.0**9., 'GB'
    write (*,'(A,1X,F5.1,1X,A)') ' It will take additional', &
        1*1*(nb**2.)*4*(2*3+1)/10.0**9., 'GB overhead'//&
        ' to run the calculation'
    write (*,'(A,1X,F5.1,1X,A)') ' Overall, you will need at least', &
        nstot*nktot*(nb**2.)*4*2*3/10.0**9. + &
        nstot*nktot*(nb**2.)*4*1/10.0**9. + &
        1*1*(nb**2.)*4*(2*3+1)/10.0**9., 'GB (+ 20% incidental) of RAM'//&
        ' to run the calculation'
END IF
! open input mommat/WAVEDER file for reading (file ID=1)
IF (wien2k) THEN
    OPEN (1, file = TRIM(fnameinp), status = 'old')
ELSE ! VASP
    OPEN (1, file = TRIM(fnameinp), form = 'unformatted', status = 'old')
END IF

!! Read VASP WAVEDER and EIGENVAL files to determine matrix elements

IF (.not.(wien2k)) THEN
    ALLOCATE( pijks(3,nb,nb,nktot,nstot), dEijks(nb,nb,nktot,nstot) )
    CALL read_mommat_pij_vasp (1, nstot, nktot, nbcder, nb, & ! <- args in
            pijks, dEijks) ! -> args out
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
    write(2,'(A53)') '# columns correspond to Cartesian directions for m_ij'
    write(2,'(A)') '# band 1=xx;     2=yy;      3=zz;     4=yz;'//&!...
        '       5=xz;      6=xy'

    fnameout3 = 'minv_c'//TRIM(charspin)//'.dat'
    OPEN (3, file = TRIM(fnameout3), status = 'UNKNOWN') ! output file 2
    write(3,'(A)') '# This file is generated by mstar'
    write(3,'(A)') '# the output contains inverse conductivity eff. masses'
    write(3,'(A)') '# m0/m_c = 1/3*Tr(m0/m*_ij) that'
    write(3,'(A)') '# are grouped by k-point index and then by the band index'
    write(3,'(A)') '# columns correspond to'
    write(3,'(A)') '# band m0/m_c'

    fnameout4 = 'minv_pr'//TRIM(charspin)//'.dat'
    OPEN (4, file = TRIM(fnameout4), status = 'UNKNOWN') ! output file 3
    write(4,'(A)') '# This file is generated by mstar'
    write(4,'(A)') '# the output contains principal components of the inverse'
    write(4,'(A)') '# eff. mass tensor eig(m0/m*_ij) that'
    write(4,'(A)') '# are grouped by k-point index and then by the band index'
    write(4,'(A)') '# columns correspond to'
    write(4,'(A)') '# band m0/m_1  m0/m_2  m0/m_3'

    fnameout5 = 'minv_d'//TRIM(charspin)//'.dat'
    OPEN (5, file = TRIM(fnameout5), status = 'UNKNOWN') ! output file 4
    write(5,'(A)') '# This file is generated by mstar'
    write(5,'(A)') '# the output contains density of states inverse effective mass' 
    write(5,'(A)') '# m0/m_d = m0/(m_1*m_2*m_3)**(1/3) that'
    write(5,'(A)') '# are grouped by k-point index and then by the band index. Here'
    write(5,'(A)') '# m0/m_1  m0/m_2  m0/m_3 are principal components of the inverse'
    write(5,'(A)') '# eff. mass tensor eig(m0/m*_ij).'
    write(5,'(A)') '# Columns correspond to'
    write(5,'(A)') '# band m0/m_d'

    !! Main part nested loops
    
    fmommatend = .false. ! end of mommat file is not reached
    ikpt = 0 ! initialize the counter for k-points
    nbb = 0 ! initialize the number of band-to-band transitions
    iline = 0 ! initialize the number of lines to skip

    DO WHILE (.not.(fmommatend)) ! loop over k-points until the file ends
        ikpt = ikpt + 1 ! count number of k-points

        !! Determine number of bands in mommat file

        IF (wien2k) THEN ! do for WIEN2k only (nb was read above in case of VASP)
            CALL read_mommat_nb(1, & ! <- args in
                    iline, & ! <-> args in-out
                    nb) ! -> args out
            ! the max number of bands for which m* is calculated
            ! here (WIEN2k) we set it = to the total number of bands
            nbcder = nb
        END IF

        !! Read <b_i|p_a|b_j> from mommat file, a=1,2,3 (x,y,z)

        ALLOCATE( pij(3,nb,nb), dEij(nb,nb) )
        IF (wien2k) THEN
            nbb = (nb+nb**2)/2 ! number of band-to-band transitions
            CALL read_mommat_pij (1, nb, nbb, & ! <- args in
                    iline, & ! <-> args in-out
                    pij, dEij) ! -> args out
        ELSE ! VASP
            pij = pijks(:,:,:,ikpt,ispin)
            dEij = dEijks(:,:,ikpt,ispin)
        END IF
        
        !! Find degenerate states within the energy tolerance

        ALLOCATE( nmlist(nb,2) )
        CALL dgenen(nb, nbcder, dEij, dEtol, & ! <- args in 
                nmlist, ndblen) ! -> args out
        ! NOTE: nbcder is not used from this point down the code
        !       it is substituted by nmlist(ndblen,2). The reason is that
        !       nbcder can cut through a degenerate group of states. Then we
        !       have to exclude the entire group of bands from the output of m*

        !! Write information about the current k-point

        write(2,'(A,I0,1X,A,I0,1X,A,I0)') & !...
            '# KP: ', ikpt, 'NBCDER: ', nmlist(ndblen,2), 'NEMAX: ', nb
        write(3,'(A,I0,1X,A,I0,1X,A,I0)') & !...
            '# KP: ', ikpt, 'NBCDER: ', nmlist(ndblen,2), 'NEMAX: ', nb
        write(4,'(A,I0,1X,A,I0,1X,A,I0)') & !...
            '# KP: ', ikpt, 'NBCDER: ', nmlist(ndblen,2), 'NEMAX: ', nb
        write(5,'(A,I0,1X,A,I0,1X,A,I0)') & !...
            '# KP: ', ikpt, 'NBCDER: ', nmlist(ndblen,2), 'NEMAX: ', nb
        
        !! prepare output formats for effective masses
        
        write(wformat2,'(I0)') nbcder ! make a character of the length 'nbcder'
        ! format line to write inverse effective masses
        write(wformat2,'(I0)') LEN(TRIM(wformat2))
        wformat2 = '(I' // TRIM(wformat2) // ',1X,5(es10.3,1X),es10.3)'
        write(wformat3,'(I0)') nbcder
        write(wformat3,'(I0)') LEN(TRIM(wformat3))
        wformat3 = '(I' // TRIM(wformat3) // ',1X,es10.3)'
        write(wformat4,'(I0)') nbcder
        write(wformat4,'(I0)') LEN(TRIM(wformat4))
        wformat4 = '(I' // TRIM(wformat4) // ',1X,2(es10.3,1X),es10.3)'
        write(wformat5,'(I0)') nbcder
        write(wformat5,'(I0)') LEN(TRIM(wformat5))
        wformat5 = '(I' // TRIM(wformat5) // ',1X,es10.3)'
        write(wformat5w,'(I0)') nbcder
        write(wformat5w,'(I0)') LEN(TRIM(wformat5w))
        wformat5w = '(I' // TRIM(wformat5w) // ',1X,es10.3,1X,A)'
            
        !! Loop through blocks of degenerate states
                
        DO id = 1,ndblen
            ! number of degenerate states in the block
            nd = nmlist(id,2) - nmlist(id,1) + 1
            ! array to store inverse masses; size 6 is because of 6 Voigt indices
            ALLOCATE( minv(nd,6) ) 
            
            ! loop over Voigt indices: 1=xx; 2=yy; 3=zz; 4=yz; 5=xz; 6=xy
            DO ivoigt = 1,6
                ! handle Voigt notations
                SELECT CASE (ivoigt) ; ; 
                    CASE (1) ! 1=xx
                        alpha = 1
                        beta = 1
                    CASE (2) ! 2=yy
                        alpha = 2
                        beta = 2
                    CASE (3) ! 3=zz
                        alpha = 3
                        beta = 3
                    CASE (4) ! 4=yz
                        alpha = 2
                        beta = 3
                    CASE (5) ! 5=xz
                        alpha = 1
                        beta = 3
                    CASE (6) ! 6=xy
                        alpha = 1
                        beta = 2
                    CASE DEFAULT
                END SELECT
                ALLOCATE( Mnm(nd,nd) )
                Mnm = (0.0, 0.0) ! initialize with complex 0's
                DO in = 1,nd
                    n = nmlist(id,1) + in - 1
                    DO im = 1,nd
                        m = nmlist(id,1) + im - 1
                        !$OMP PARALLEL DO REDUCTION(+:Mnm) &
                        !$OMP PRIVATE(dM,p2,dE,k)
                        DO k = 1,nb ! all bands (also above nbcder) are included
                            ! 'k' is outside of the degeneracy list
                            IF (k<nmlist(id,1) .or. k > nmlist(id,2)) THEN
                                p2 = pij(alpha,n,k)*pij(beta,k,m) + & !...
                                    pij(beta,n,k)*pij(alpha,k,m)
                                ! account for slight violation of degeneracy by
                                ! taking an average energy difference
                                dE = -dEij(n,k)/2 -dEij(m,k)/2
                                dM = p2/dE
                                ! make sure dM is finite (not NaN and not Inf)
                                IF (dM /= dM .or. abs(dM) > HUGE(abs(dM))) THEN
                                    write(*,*) 'ikpt =', ikpt
                                    write(*,*) 'n =', n
                                    write(*,*) 'k =', k
                                    write(*,*) 'm =', m
                                    write(*,*) 'alpha =', alpha
                                    write(*,*) 'beta =', beta
                                    write(*,*) 'pij(alpha,n,k) =', &
                                        pij(alpha,n,k)
                                    write(*,*) 'pij(beta,k,m) =', &
                                        pij(beta,k,m)
                                    write(*,*) 'pij(beta,n,k) =', &
                                        pij(beta,n,k)
                                    write(*,*) 'pij(alpha,k,m) =', &
                                        pij(alpha,k,m)
                                    write(*,*) 'dEij(n,k) =', dEij(n,k)
                                    write(*,*) 'dEij(m,k) =', dEij(m,k)
                                    write(*,*) 'dM = ', dM
                                    write(*,*) 'dE = ', dE
                                    write(*,*) 'p2 = ', p2
                                    STOP 'Error: dM is not finite'
                                END IF
                                ! update Mnm array
                                Mnm(in,im) = Mnm(in,im) + dM
                            END IF
                        END DO ! loop over 'ik'
                        !$OMP END PARALLEL DO
                    END DO ! loop over 'im'
                END DO ! loop over 'in'
                ! find eigenvalues of Mnm
                ALLOCATE( EIGV(nd), EIGFC(nd,nd) )
                CALL eigvz(nd, Mnm, & ! <- args in 
                        EIGV, EIGFC) ! -> args out
                ! store the inverse effective masses (+ 1 if alpha=beta)
                minv(:,ivoigt) = EIGV + kdelta(alpha,beta)
                DEALLOCATE( Mnm, EIGV, EIGFC ) ! EIGFC is not used
            END DO ! loop over 'ivoigt'
            ! store the inverse effective masses
            DO i = 1,nd
                write(2,TRIM(wformat2)) nmlist(id,1)+i-1, (minv(i,j), j=1,6)
                ! the trace of the matrix does not change after diagonalization
                ! so the sum of eigenvalues will be just a trace
                write(3,TRIM(wformat3)) nmlist(id,1)+i-1, & !...
                    (1.0/3.0)*SUM(minv(i,1:3))
                ! determine principal components of the eff. mass tensor
                minv_U(1,1) = minv(i,1) ! the upper triangular part of minv tensor
                minv_U(2,2) = minv(i,2) ! 1=xx; 2=yy; 3=zz; 4=yz; 5=xz; 6=xy
                minv_U(3,3) = minv(i,3)
                minv_U(2,3) = minv(i,4)
                minv_U(1,3) = minv(i,5)
                minv_U(1,2) = minv(i,6)
                ALLOCATE( EIGV(3), EIGFR(3,3) )
                CALL eigvs(3, minv_U, & ! <- args in 
                        EIGV, EIGFR) ! -> args out
                write(4,TRIM(wformat4)) nmlist(id,1)+i-1, EIGV(1:3)
                ! DOS effective mass
                ! determinant of inv. eff. mass tensor
                minv_det = EIGV(1)*EIGV(2)*EIGV(3)
                minv_d = SIGN(abs(minv_det)**(1.0/3.0) , minv_det)
                IF (((EIGV(1) .LT. 0.0) .AND. (EIGV(2) .LT. 0.0) &!...
                        .AND. (EIGV(3) .LT. 0.0)) .OR. &!...
                        ((EIGV(1) .GT. 0.0) .AND. (EIGV(2) .GT. 0.0) &!...
                        .AND. (EIGV(3) .GT. 0.0))) THEN
                    ! parabolic band dispersion
                    write(5,TRIM(wformat5)) nmlist(id,1)+i-1, &!...
                        minv_d
                ELSE
                    ! hyperbolic band dispersion (signs of masses are different)
                    write(5,TRIM(wformat5w)) nmlist(id,1)+i-1, &!...
                        minv_d, 'warning: hyperbolic dispersion'
                    warnHyperb = .true.
                END IF
                DEALLOCATE( EIGV, EIGFR ) ! EIGFR is not used
            END DO
            DEALLOCATE( minv )
        END DO ! loop over 'id'
        DEALLOCATE( nmlist, pij, dEij ) ! all k-point specific variables
        
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

CLOSE (1) ! input  file 1

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
    write(*,*) '    There was a warning about hyperbolic band dispersion in the output file.'
    write(*,*) '    This indicates a limited validity of the geometric average'
    write(*,*) '    (m_1*m_2*m_3)**(1/3) for the DOS effective mass in the situations when'
    write(*,*) '    the band dispersion is _not_ parabolic.'
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
write(*,*) '       max dE for 2 states to be considered as degenerate'
write(*,*) '       default value is 1.0e-6 Ha'
write(*,*) ''
write(*,*) 'Output:'
write(*,*) 'minv_ij.dat - contains elements of the (m0/m*_ij) tensor.'
write(*,*) 'minv_c.dat - conductivity effective mass (m0/m*_c).'
write(*,*) 'minv_pr.dat - contains principal components of the inverse'
write(*,*) '              eff. mass tensor eig(m0/m*_ij)'
write(*,*) 'minv_d.dat - density of states effective mass (m0/m*_d).'
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
STOP

END PROGRAM mstar
