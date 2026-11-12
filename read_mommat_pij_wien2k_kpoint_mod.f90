MODULE read_mommat_pij_wien2k_kpoint_mod

USE precision_mod, ONLY: sp

IMPLICIT NONE

PRIVATE
PUBLIC :: read_mommat_pij_wien2k_kpoint

CONTAINS

SUBROUTINE read_mommat_pij_wien2k_kpoint( &
        unitinp, unitinpUP, unitinpDN, soc_sp_resolv_vij, &
        ikpt, nbki, nbcderki, &
        vij, dEij, vijUP, vijDN)

! Convention of case.mommat2[up/dn]
! ---------------------------------
! Taken from the WIEN2k 24.1 source of OPTIC (SRC_optic):
!
!   sph-UP.F, lines 690-709: the spin blocks OPMX(:,s), s = 1, 2, hold
!   <u_NB1|grad|u_NB2> restricted to the up-up and down-down spinor
!   components (zero-initialized, so s = 2 stays zero without SOC), and
!   the total is OPMATX = -i (OPMX(:,1) + OPMX(:,2)) = <u_NB1|-i grad|u_NB2>.
!   Line 715 documents the same quantity as <j,k|-i*grad|j',k>.
!   The interstitial part is added in planew.F, line 282.
!
!   outmat.f, lines 146-154: each record (NB1,NB2) with NB1 <= NB2 holds
!   the bra band NB1, the ket band NB2, the three Cartesian components
!   of <u_NB1|-i grad|u_NB2>, and E(NB2) - E(NB1) in Ry. The total is
!   written to unit 4 (case.mommat2). With SOC, the spin blocks are
!   written as -i OPMX(:,1) and -i OPMX(:,2) to units 32 and 33
!   (case.mommat2up and case.mommat2dn of an OPTIC run with SOC).
!   Without SOC, a spin-polarized calculation writes each spin channel
!   to its own unit-4 file (case.mommat2up or case.mommat2dn).
!
! Hence each file holds momentum matrix elements p = -i hbar grad in units
! of hbar/bohr with the same phase convention, and the total equals the
! sum of the two spin-diagonal blocks, p = p^upup + p^dndn. BerryCPT works
! in Hartree atomic units (hbar = m0 = 1), where the velocity matrix element
! v = p/m0 has the same numerical value as p. The conversion v = p/m0 is
! valid only for a local potential (see README.md). The reverse element
! is <u_NB2|p|u_NB1> = CONJG(<u_NB1|p|u_NB2>), since p is Hermitian.

! File units for the total and optional spin-resolved matrix elements
INTEGER, INTENT(IN) :: &
    unitinp, &
    unitinpUP, &
    unitinpDN

! Spin-resolved matrix-element input
LOGICAL, INTENT(IN) :: soc_sp_resolv_vij

! Current k-point metadata
INTEGER, INTENT(IN) :: &
    ikpt, &
    nbki, &
    nbcderki

! Total velocity matrix elements v = p/m0 [a.u.] and energy differences [Ha]
COMPLEX(KIND=sp), ALLOCATABLE, INTENT(OUT) :: vij(:,:,:)
REAL(KIND=sp), ALLOCATABLE, INTENT(OUT) :: dEij(:,:)

! Optional spin-diagonal velocity matrix elements [a.u.]
COMPLEX(KIND=sp), ALLOCATABLE, INTENT(OUT) :: &
    vijUP(:,:,:), &
    vijDN(:,:,:)

CALL read_single_mommat_kpoint( &
    unitinp, 'total', ikpt, nbki, nbcderki, & ! in
    vij, dEij)                                ! out

IF (soc_sp_resolv_vij) THEN
    CALL read_single_mommat_kpoint( &
        unitinpUP, 'spin-up', ikpt, nbki, nbcderki, & ! in
        vijUP)                                      ! out

    CALL read_single_mommat_kpoint( &
        unitinpDN, 'spin-down', ikpt, nbki, nbcderki, & ! in
        vijDN)                                        ! out
END IF

END SUBROUTINE read_mommat_pij_wien2k_kpoint


SUBROUTINE read_single_mommat_kpoint( &
        unitinp, file_label, ikpt, nbki, nbcderki, &
        vij, dEij)

! Matrix-element file and its role
INTEGER, INTENT(IN) :: unitinp
CHARACTER(LEN=*), INTENT(IN) :: file_label

! Expected k-point metadata
INTEGER, INTENT(IN) :: &
    ikpt, &
    nbki, &
    nbcderki

! Rectangular momentum matrix elements
COMPLEX(KIND=sp), ALLOCATABLE, INTENT(OUT) :: vij(:,:,:)

! Rectangular energy differences, required only for the total file
REAL(KIND=sp), ALLOCATABLE, OPTIONAL, INTENT(OUT) :: dEij(:,:)

! Current input record and I/O diagnostics
CHARACTER(LEN=256) :: cline, iomsg
INTEGER :: ierr

! Matrix-record indices and counters
INTEGER :: &
    bii, &
    bjj, &
    irow, &
    nbb

! Matrix elements and energy difference read from one record
REAL(KIND=sp) :: &
    p1_re, p1_im, &
    p2_re, p2_im, &
    p3_re, p3_im, &
    energy_difference

! Momentum matrix element <u_bii|p|u_bjj> read from one record
COMPLEX(KIND=sp) :: pvalue(3)

IF (ikpt < 1) THEN
    CALL mommat_error( &
        file_label, ikpt, & ! in
        'Invalid k-point index supplied to the WIEN2k reader') ! in
END IF

IF (nbki < 1) THEN
    CALL mommat_error( &
        file_label, ikpt, & ! in
        'Invalid number of bands supplied to the WIEN2k reader') ! in
END IF

IF (nbcderki < 1 .OR. nbcderki > nbki) THEN
    CALL mommat_error( &
        file_label, ikpt, & ! in
        'Invalid derivative-band limit supplied to the WIEN2k reader') ! in
END IF

CALL read_mommat_kpoint_header( &
    unitinp, file_label, ikpt, nbki) ! in

CALL allocate_vij( &
    file_label, ikpt, nbki, nbcderki, & ! in
    vij)                               ! out

IF (PRESENT(dEij)) THEN
    CALL allocate_dEij( &
        file_label, ikpt, nbki, nbcderki, & ! in
        dEij)                              ! out
END IF

vij = CMPLX(0.0_sp,0.0_sp,KIND=sp)
IF (PRESENT(dEij)) dEij = 0.0_sp

nbb = nbki * (nbki + 1) / 2

irow = 0
DO WHILE (irow < nbb)
    READ(unitinp,'(A)',IOSTAT=ierr,IOMSG=iomsg) cline

    IF (ierr < 0) THEN
        CALL mommat_error( &
            file_label, ikpt, & ! in
            'Reached the end of the file before all matrix elements were read') ! in
    ELSE IF (ierr > 0) THEN
        CALL mommat_error( &
            file_label, ikpt, & ! in
            'Unable to read a matrix-element record: '//TRIM(iomsg)) ! in
    END IF

    IF (LEN_TRIM(cline) == 0) CYCLE

    IF (INDEX(cline,'KP:') /= 0) THEN
        CALL mommat_error( &
            file_label, ikpt, & ! in
            'Encountered the next k-point header before all matrix elements were read') ! in
    END IF

    irow = irow + 1

    READ(cline,*,IOSTAT=ierr,IOMSG=iomsg) &
        bii, bjj, &
        p1_re, p1_im, p2_re, p2_im, p3_re, p3_im, &
        energy_difference

    IF (ierr /= 0) THEN
        WRITE(*,'(/,A)') ' ERROR: Unable to parse a WIEN2k matrix-element record'
        WRITE(*,'(A)') ' File role = '//TRIM(file_label)
        WRITE(*,'(A,I0)') ' K-point = ', ikpt
        WRITE(*,'(A,I0)') ' Matrix record = ', irow
        WRITE(*,'(A)') ' Record = '//TRIM(cline)
        WRITE(*,'(A)') ' I/O message = '//TRIM(iomsg)
        ERROR STOP 1
    END IF

    IF (bii < 1 .OR. bii > nbki .OR. bjj < 1 .OR. bjj > nbki) THEN
        WRITE(*,'(/,A)') ' ERROR: Invalid band index in a WIEN2k matrix-element record'
        WRITE(*,'(A)') ' File role = '//TRIM(file_label)
        WRITE(*,'(A,I0)') ' K-point = ', ikpt
        WRITE(*,'(A,I0)') ' First band index = ', bii
        WRITE(*,'(A,I0)') ' Second band index = ', bjj
        WRITE(*,'(A,I0)') ' Number of bands = ', nbki
        ERROR STOP 1
    END IF

    IF (bii > bjj) THEN
        WRITE(*,'(/,A)') ' ERROR: WIEN2k matrix records are not upper triangular'
        WRITE(*,'(A)') ' File role = '//TRIM(file_label)
        WRITE(*,'(A,I0)') ' K-point = ', ikpt
        WRITE(*,'(A,I0)') ' First band index = ', bii
        WRITE(*,'(A,I0)') ' Second band index = ', bjj
        ERROR STOP 1
    END IF

    IF (irow == 1) THEN
        IF (bii /= 1 .OR. bjj /= 1) THEN
            CALL mommat_error( &
                file_label, ikpt, & ! in
                'The first matrix-element record is not for bands 1 and 1') ! in
        END IF
    ELSE IF (irow == nbb) THEN
        IF (bii /= nbki .OR. bjj /= nbki) THEN
            CALL mommat_error( &
                file_label, ikpt, & ! in
                'The last matrix-element record does not contain the final band') ! in
        END IF
    END IF

    ! <u_bii|p|u_bjj> [hbar/bohr], numerically equal to the velocity
    ! matrix element v = p/m0 in Hartree atomic units
    pvalue(1) = CMPLX(p1_re,p1_im,KIND=sp)
    pvalue(2) = CMPLX(p2_re,p2_im,KIND=sp)
    pvalue(3) = CMPLX(p3_re,p3_im,KIND=sp)

    ! WIEN2k prints E_j-E_i in Ry for the upper-triangular record (i,j).
    ! BerryCPT stores dEij(i,j)=E_j-E_i in Ha.
    energy_difference = energy_difference / 2.0_sp

    IF (bii <= nbcderki) THEN
        vij(:,bii,bjj) = pvalue
        IF (PRESENT(dEij)) dEij(bii,bjj) = energy_difference
    END IF

    IF (bjj <= nbcderki) THEN
        vij(:,bjj,bii) = CONJG(pvalue)
        IF (PRESENT(dEij)) dEij(bjj,bii) = -energy_difference
    END IF
END DO

END SUBROUTINE read_single_mommat_kpoint


SUBROUTINE read_mommat_kpoint_header( &
        unitinp, file_label, ikpt, nbki)

! Matrix-element file and its role
INTEGER, INTENT(IN) :: unitinp
CHARACTER(LEN=*), INTENT(IN) :: file_label

! Expected k-point metadata
INTEGER, INTENT(IN) :: ikpt, nbki

! Current input record and I/O diagnostics
CHARACTER(LEN=256) :: cline, iomsg
INTEGER :: ierr

! Metadata read from the k-point header
INTEGER :: ikpt_file, nbmin_file, nbmax_file

DO
    READ(unitinp,'(A)',IOSTAT=ierr,IOMSG=iomsg) cline

    IF (ierr < 0) THEN
        CALL mommat_error( &
            file_label, ikpt, & ! in
            'Reached the end of the file while searching for the k-point header') ! in
    ELSE IF (ierr > 0) THEN
        CALL mommat_error( &
            file_label, ikpt, & ! in
            'Unable to read while searching for the k-point header: '//TRIM(iomsg)) ! in
    END IF

    IF (INDEX(cline,'KP:') /= 0) EXIT
END DO

READ(cline,'(7X,I5)',IOSTAT=ierr,IOMSG=iomsg) ikpt_file
IF (ierr /= 0) THEN
    CALL mommat_error( &
        file_label, ikpt, & ! in
        'Unable to read the k-point index from the header: '//TRIM(iomsg)) ! in
END IF

READ(cline,'(27X,2I5)',IOSTAT=ierr,IOMSG=iomsg) &
    nbmin_file, nbmax_file
IF (ierr /= 0) THEN
    CALL mommat_error( &
        file_label, ikpt, & ! in
        'Unable to read the band range from the k-point header: '//TRIM(iomsg)) ! in
END IF

IF (ikpt_file /= ikpt) THEN
    WRITE(*,'(/,A)') ' ERROR: Unexpected k-point index in a WIEN2k matrix-element file'
    WRITE(*,'(A)') ' File role = '//TRIM(file_label)
    WRITE(*,'(A,I0)') ' Expected k-point = ', ikpt
    WRITE(*,'(A,I0)') ' K-point in file = ', ikpt_file
    ERROR STOP 1
END IF

IF (nbmin_file /= 1 .OR. nbmax_file /= nbki) THEN
    WRITE(*,'(/,A)') ' ERROR: Unexpected band range in a WIEN2k matrix-element file'
    WRITE(*,'(A)') ' File role = '//TRIM(file_label)
    WRITE(*,'(A,I0)') ' K-point = ', ikpt
    WRITE(*,'(A,I0)') ' First band in file = ', nbmin_file
    WRITE(*,'(A,I0)') ' Last band in file = ', nbmax_file
    WRITE(*,'(A,I0)') ' Expected number of bands = ', nbki
    ERROR STOP 1
END IF

END SUBROUTINE read_mommat_kpoint_header


SUBROUTINE allocate_vij( &
        file_label, ikpt, nbki, nbcderki, &
        vij)

! Matrix-element file role and current dimensions
CHARACTER(LEN=*), INTENT(IN) :: file_label
INTEGER, INTENT(IN) :: ikpt, nbki, nbcderki

! Momentum matrix allocation
COMPLEX(KIND=sp), ALLOCATABLE, INTENT(OUT) :: vij(:,:,:)

! Allocation status and message
INTEGER :: ierr
CHARACTER(LEN=256) :: iomsg

ALLOCATE( &
    vij(3,nbcderki,nbki), &
    STAT=ierr, ERRMSG=iomsg)

IF (ierr /= 0) THEN
    WRITE(*,'(/,A)') ' ERROR: Unable to allocate a WIEN2k momentum matrix'
    WRITE(*,'(A)') ' File role = '//TRIM(file_label)
    WRITE(*,'(A,I0)') ' K-point = ', ikpt
    WRITE(*,'(A,3(I0,1X))') ' Requested shape = ', 3, nbcderki, nbki
    WRITE(*,'(A)') ' Allocation message = '//TRIM(iomsg)
    ERROR STOP 1
END IF

END SUBROUTINE allocate_vij


SUBROUTINE allocate_dEij( &
        file_label, ikpt, nbki, nbcderki, &
        dEij)

! Matrix-element file role and current dimensions
CHARACTER(LEN=*), INTENT(IN) :: file_label
INTEGER, INTENT(IN) :: ikpt, nbki, nbcderki

! Energy-difference matrix allocation
REAL(KIND=sp), ALLOCATABLE, INTENT(OUT) :: dEij(:,:)

! Allocation status and message
INTEGER :: ierr
CHARACTER(LEN=256) :: iomsg

ALLOCATE( &
    dEij(nbcderki,nbki), &
    STAT=ierr, ERRMSG=iomsg)

IF (ierr /= 0) THEN
    WRITE(*,'(/,A)') ' ERROR: Unable to allocate a WIEN2k energy-difference matrix'
    WRITE(*,'(A)') ' File role = '//TRIM(file_label)
    WRITE(*,'(A,I0)') ' K-point = ', ikpt
    WRITE(*,'(A,2(I0,1X))') ' Requested shape = ', nbcderki, nbki
    WRITE(*,'(A)') ' Allocation message = '//TRIM(iomsg)
    ERROR STOP 1
END IF

END SUBROUTINE allocate_dEij


SUBROUTINE mommat_error( &
        file_label, ikpt, message)

! Matrix-element file role and current k-point
CHARACTER(LEN=*), INTENT(IN) :: file_label
INTEGER, INTENT(IN) :: ikpt

! Error description
CHARACTER(LEN=*), INTENT(IN) :: message

WRITE(*,'(/,A)') ' ERROR: '//TRIM(message)
WRITE(*,'(A)') ' File role = '//TRIM(file_label)
WRITE(*,'(A,I0)') ' K-point = ', ikpt
ERROR STOP 1

END SUBROUTINE mommat_error

END MODULE read_mommat_pij_wien2k_kpoint_mod
