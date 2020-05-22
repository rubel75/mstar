SUBROUTINE read_mommat_pij_vasp (fid, nstot, nktot, nbcder, nb, & ! <- args in
                pijks, dEijks) ! -> args out

! Read dipole matrix elements <i|r_a|j> (a=x,y,z) from WAVEDER and 
! energy differences E_i - E_j from EIGENVAL files. Then convert
! dipole matrix elements <i|r_a|j> into the momentum matrix elements <i|p_a|j>

!! Variables in-out

implicit none
INTEGER, intent(in) :: &
    fid, & ! file ID
    nstot, & ! total number of spins in WAVEDER file
    nktot, & ! total number of k-points in WAVEDER file
    nbcder, & ! max number of bands for which m* is calculated
    nb ! max number of bands
REAL(kind=4), intent(out) :: &
    dEijks(nb,nb,nktot,nstot) ! k- and spin-dependent energy differences 
                              ! E_i-E_j (at.u.)
COMPLEX(kind=4), intent(out) :: &
    pijks(3,nb,nb,nktot,nstot) ! k- and spin-dependent momentum matrix
                               ! elements (at.u.)

!! Variables internal

CHARACTER(len=256) :: cline ! text in the current line
INTEGER :: ispin, ikpt, iband, idim, im, in, i ! counters
INTEGER(kind=4) :: dum1, dum2, dum3, dum4 ! not used
REAL(kind=8) :: NODES_IN_DIELECTRIC_FUNCTION, WPLASMON(3,3) ! not used
REAL(kind=4) :: occup(nstot) ! not used
REAL(kind=8) :: ene(nb,nktot,nstot) ! energy eigenvalues (eV)
COMPLEX(kind=4) ::  rijks(nb,nbcder,nktot,nstot,3) ! dipole matrix elements (Ang)


!! Read dipole matrix element from WAVEDER file

! skip through header lines
READ(fid) dum1, dum2, dum3, dum4 ! WDES%NB_TOT, NBANDS_CDER, WDES%NKPTS, WDES%ISPIN
READ(fid) NODES_IN_DIELECTRIC_FUNCTION ! not used
READ(fid) WPLASMON ! not used
! read matrix elements
READ(fid) rijks ! VASP calls it CDER_BETWEEN_STATES or CHAM

!! Read eigenvalues from EIGENVAL file

OPEN (99, file = 'EIGENVAL', status = 'old')
DO i = 1,6 ! skip first 6 lines
    READ (99,*, END=10)
END DO
! initialize array
ene=0.0
! read eigenvalues into ene(iband,ikpt)
DO ikpt = 1, nktot ! loop over k-points
    READ (99,*, END=10) ! skip 2 lines
    READ (99,*, END=10)
    DO iband = 1, nb ! loop over bands
        READ(99,'(A)',END=10) cline
        IF (nstot == 1) THEN
           ! Standard VASP format for EIGENVAL is (1X,I4,4X,F12.6,2X,F9.6)
           ! It includes band index, ene(iband,ikpt,1), occup
           ! When the number of bands > 9999 it dies not fit in (I4) and
           ! writes '****' in the EIGENVAL for the band index.
           ! Since we do not use the number of bands from EIGENVAL, we can
           ! just skip the first 9 characters, i.e., merge (1X,I4,4X) => (9X)
            READ(cline,'(9X,F12.6,2X,F9.6)',ERR=20) &
                ene(iband,ikpt,1), occup
        ELSE
            ! (1) spin up, (2) spin down (see the above note on the format)
            READ(cline,'(9X,F12.6,2X,F12.6,2X,F9.6,2X,F9.6)',ERR=30) &
                ene(iband,ikpt,1), ene(iband,ikpt,2), occup(1), occup(2)
        END IF
    END DO ! bands
END DO ! k-points
10 CLOSE(99) ! EIGENVAL
! cross check
IF (iband-1==nb .AND. ikpt-1==nktot) THEN
    write (*,*) 'WAVEDER and EIGENVAL files were read successfully'
ELSE
    write (*,*) 'Error in reading EIGENVAL file'
    write(*,'(I0,A,I0)') ikpt-1, &
        ' k-points found in EIGENVAL file, while expecting ', nktot
    write(*,'(I0,A,I0)') iband-1, &
        ' bands found in EIGENVAL file, while expecting ', nb
    STOP
END IF

!! Determine energy difference between pairs of states in [Ha]

dEijks=0.0
DO ispin = 1, nstot ! loop over spins
    DO ikpt = 1, nktot ! loop over k-points
        DO in = 1,nb-1
            DO im = in+1,nb
                ! 27.211386245988 [eV -> Ha]
                dEijks(in,im,ikpt,ispin) = &
                    ( ene(im,ikpt,ispin) - ene(in,ikpt,ispin) )/27.211386245988
                dEijks(im,in,ikpt,ispin) = -dEijks(in,im,ikpt,ispin)
            END DO
        END DO
    END DO
END DO

!! Convert dipole matrix elements in [Ang] into momentum matrix elements
!! in [at.u.]

pijks = (0.0, 0.0)
DO ispin = 1, nstot ! loop over spins
    DO ikpt = 1, nktot ! loop over k-points
        DO im = 1,nbcder-1
            DO in = im+1,nb
                DO idim = 1, 3
                    ! 1.88972623847507 [Ang -> bohr]
                    pijks(idim,in,im,ikpt,ispin) = & !...
                        rijks(in,im,ikpt,ispin,idim)*1.88972623847507* & !...
                        dEijks(in,im,ikpt,ispin)
                    ! ensure Hermitian properties of pij
                    pijks(idim,im,in,ikpt,ispin) = & !...
                        CONJG( pijks(idim,in,im,ikpt,ispin) )
                    ! debug
                    ! write(*,'(I0,1X,I0,1X,I0,1X,I0,1X,I0,1X,E15.6,1X,E15.6)') &
                    !     idim,in,im,ikpt,ispin,pijks(idim,in,im,ikpt,ispin)
                END DO
            END DO
        END DO
    END DO
END DO
RETURN

!! Read error handling

20 & ! read error in EIGENVAL (single spin)
write(*,*) 'ERROR while reading EIGENVAL (single spin mode)'
write(*,*) 'Here is the line that caused the error'
write(*,'(A)') TRIM(cline)
write(*,*) 'while expected'
write(*,'(A)') ' 9999       20.826607   0.000000'
write(*,*) 'Execution terminated'
STOP

30 & ! read error in EIGENVAL (spin up & down)
write(*,*) 'ERROR while reading EIGENVAL (2 spins mode)'
write(*,*) 'Here is the line that caused the error'
write(*,'(A)') TRIM(cline)
write(*,*) 'while expected'
write(*,'(A)') ' 9999       20.826607        20.826607   0.000000     0.000000'
write(*,*) 'Execution terminated'
STOP

END SUBROUTINE read_mommat_pij_vasp
