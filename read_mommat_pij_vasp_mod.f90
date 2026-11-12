MODULE read_mommat_pij_vasp_mod

! Copied from BerryCPT (https://github.com/rubel75/berrycpt, commit 993047e)
! with one modification for mstar: the eigenvalues ene are double precision,
! and energy differences are formed in double precision before they are
! stored in single precision (see read_eigenvalues_vasp_mod.f90).

USE, INTRINSIC :: iso_fortran_env, ONLY: int32, int64, output_unit
USE precision_mod, ONLY: sp, dp

IMPLICIT NONE

PRIVATE
PUBLIC :: read_mommat_pij_vasp

CONTAINS

SUBROUTINE read_mommat_pij_vasp( &
        fnameinp, nstot, nktot, nbcder, nb, ene, &
        vijks, dEijks)

! Name of the VASP WAVEDER file
CHARACTER(LEN=*), INTENT(IN) :: fnameinp

! WAVEDER metadata
INTEGER, INTENT(IN) :: &
    nstot, &  ! total number of spin channels
    nktot, &  ! total number of k-points
    nbcder, & ! number of bands included in derivative overlaps
    nb         ! total number of bands

! Band-, k-point-, and spin-dependent eigenvalues [Ha]
REAL(KIND=dp), INTENT(IN) :: ene(:,:,:)

! Rectangular velocity matrix elements [a.u.]
! vijks(:,i,j,k,s) stores <i|v|j>, where i <= nbcder and j <= nb.
COMPLEX(KIND=sp), ALLOCATABLE, INTENT(OUT) :: vijks(:,:,:,:,:)

! Rectangular energy differences [Ha]
! dEijks(i,j,k,s) stores E_j-E_i, where i <= nbcder and j <= nb.
REAL(KIND=sp), ALLOCATABLE, INTENT(OUT) :: dEijks(:,:,:,:)

! WAVEDER metadata read from the first record
INTEGER(KIND=int32) :: &
    nb_file, &
    nbcder_file, &
    nktot_file, &
    nstot_file

! WAVEDER records not used by BerryCPT
REAL(KIND=dp) :: &
    nodes_in_dielectric_function, &
    wplasmon(3,3)

! Input unit and I/O status
INTEGER :: unitinp, ierr
CHARACTER(LEN=256) :: iomsg

! Band, k-point, spin, and Cartesian counters
INTEGER :: &
    iband, iband_cder, ikpt, ispin, idim

! Energy difference for the current band pair [Ha]
REAL(KIND=sp) :: energy_difference

! Memory required by the returned arrays [bytes]
INTEGER(KIND=int64) :: &
    vijks_bytes, &
    dEijks_bytes

! Unit conversion from Angstrom to bohr
REAL(KIND=sp), PARAMETER :: angstrom_to_bohr = 1.8897261_sp

CALL validate_dimensions( &
    fnameinp, nstot, nktot, nbcder, nb, ene) ! in

CALL validate_waveder_kinds()

OPEN(NEWUNIT=unitinp,FILE=TRIM(fnameinp),FORM='UNFORMATTED', &
    ACCESS='SEQUENTIAL',STATUS='OLD',ACTION='READ', &
    IOSTAT=ierr,IOMSG=iomsg)
IF (ierr /= 0) THEN
    CALL waveder_file_error( &
        fnameinp, & ! in
        'Unable to open the WAVEDER file: '//TRIM(iomsg)) ! in
END IF

READ(unitinp,IOSTAT=ierr,IOMSG=iomsg) &
    nb_file, nbcder_file, nktot_file, nstot_file
IF (ierr /= 0) THEN
    CLOSE(unitinp)
    CALL waveder_file_error( &
        fnameinp, & ! in
        'Unable to read the WAVEDER metadata: '//TRIM(iomsg)) ! in
END IF

IF (nb_file /= nb .OR. nbcder_file /= nbcder .OR. &
        nktot_file /= nktot .OR. nstot_file /= nstot) THEN
    CLOSE(unitinp)
    WRITE(*,'(/,A)') ' ERROR: Inconsistent WAVEDER metadata'
    WRITE(*,'(A,I0,A,I0)') &
        ' Number of bands: expected ', nb, ', found ', nb_file
    WRITE(*,'(A,I0,A,I0)') &
        ' Number of derivative bands: expected ', nbcder, &
        ', found ', nbcder_file
    WRITE(*,'(A,I0,A,I0)') &
        ' Number of k-points: expected ', nktot, ', found ', nktot_file
    WRITE(*,'(A,I0,A,I0)') &
        ' Number of spin channels: expected ', nstot, &
        ', found ', nstot_file
    ERROR STOP 1
END IF

READ(unitinp,IOSTAT=ierr,IOMSG=iomsg) nodes_in_dielectric_function
IF (ierr /= 0) THEN
    CLOSE(unitinp)
    CALL waveder_file_error( &
        fnameinp, & ! in
        'Unable to read the second WAVEDER record: '//TRIM(iomsg)) ! in
END IF

READ(unitinp,IOSTAT=ierr,IOMSG=iomsg) wplasmon
IF (ierr /= 0) THEN
    CLOSE(unitinp)
    CALL waveder_file_error( &
        fnameinp, & ! in
        'Unable to read the third WAVEDER record: '//TRIM(iomsg)) ! in
END IF

CALL estimate_output_memory( &
    nstot, nktot, nbcder, nb, & ! in
    vijks_bytes, dEijks_bytes)   ! out

CALL report_output_memory( &
    vijks_bytes, dEijks_bytes) ! in

ALLOCATE( &
    vijks(3,nbcder,nb,nktot,nstot), &
    STAT=ierr, ERRMSG=iomsg)
IF (ierr /= 0) THEN
    CLOSE(unitinp)
    CALL allocation_error( &
        fnameinp, 'vijks', vijks_bytes, iomsg) ! in
END IF

ALLOCATE( &
    dEijks(nbcder,nb,nktot,nstot), &
    STAT=ierr, ERRMSG=iomsg)
IF (ierr /= 0) THEN
    DEALLOCATE(vijks)
    CLOSE(unitinp)
    CALL allocation_error( &
        fnameinp, 'dEijks', dEijks_bytes, iomsg) ! in
END IF

! WAVEDER stores D^a_ji = <u_j| -i d/dk_a u_i> with the all-band (bra)
! index j first and the derivative-band (ket) index i second (see the
! convention notes below). Read the complete record directly into
! vijks(:,i,j,k,s), avoiding a second full-size derivative-overlap array.
READ(unitinp,IOSTAT=ierr,IOMSG=iomsg) &
    (((((vijks(idim,iband_cder,iband,ikpt,ispin), &
        iband = 1, nb), &
        iband_cder = 1, nbcder), &
        ikpt = 1, nktot), &
        ispin = 1, nstot), &
        idim = 1, 3)
IF (ierr /= 0) THEN
    CLOSE(unitinp)
    CALL waveder_file_error( &
        fnameinp, & ! in
        'Unable to read the WAVEDER derivative overlaps: '//TRIM(iomsg)) ! in
END IF

CLOSE(unitinp,IOSTAT=ierr,IOMSG=iomsg)
IF (ierr /= 0) THEN
    CALL waveder_file_error( &
        fnameinp, & ! in
        'Unable to close the WAVEDER file: '//TRIM(iomsg)) ! in
END IF

! Convert D^a_ji [Angstrom] to v^a_ij [a.u.] using Hermiticity.
!
! Convention of WAVEDER
! ---------------------
! The two available descriptions of WAVEDER differ by a factor -i.
!
! (1) VASP wiki, page "WAVEDER" (https://www.vasp.at/wiki/index.php/WAVEDER),
!     revision oldid=35849, last edited 15 April 2026, read 11 November 2026.
!     It defines the stored element as <u_n'|d u_n/d k_a> in Angstrom,
!     without a factor -i, with the array order (n', n, k, a).
!
! (2) VASP 6.6.1 source, linear_optics.F, module mlr_optic.
!     The module header (lines 13-14 and 87-88) defines
!         CDER_BETWEEN_STATES(m,n,k,s,a) = <u_m| -i d/dk_a u_n>
!                                        = -<u_m|r_a|u_n>,
!     LRF_EPSILON_MELS (lines 737 and 751) builds the first-order change
!     -i d/dk u_n = -r u_n and its overlaps with all bands m, and INPROD_W
!     (line 1115) stores CHAM(j,i) = <W0(j)|W1(i)>, so the bra band m
!     (1:NB_TOT) is the first index and the ket band n (1:NBANDS_CDER) the
!     second. WRT_CDER_BETWEEN_STATES (line 3118) writes the whole array
!     CDER_BETWEEN_STATES(NB_TOT, NBANDS_CDER, NKPTS, ISPIN, 3) as one
!     record. The routine LRF_RPHI0, which solves for the first-order
!     change, is in another source file and was not inspected.
!
! BerryCPT follows the source code (2). The record is read here as
!
!   D^a_ji = <u_j| -i d/dk_a u_i> = -r^a_ji,
!
! with j = iband (all bands) and i = iband_cder (i <= NBANDS_CDER), in
! Angstrom.
!
! Velocity matrix elements
! ------------------------
! With the dipole matrix element r^a_ji = i <u_j|d u_i/d k_a> = -D^a_ji,
! the expression used in README.md and in the paper,
!
!   v^a_ji = (i/hbar) (E_j - E_i) r^a_ji,
!
! becomes
!
!   v^a_ji = (i/hbar) (E_i - E_j) D^a_ji,
!   v^a_ij = CONJG(v^a_ji) = (i/hbar) (E_j - E_i) CONJG(D^a_ji).
!
! The result is in Hartree atomic units (hbar = m0 = 1), with energies in
! Ha and D converted from Angstrom to bohr. Because the derivative overlap
! includes the nonlocal PAW terms, these are velocity matrix elements, not
! momentum matrix elements divided by m0.
!
! Every element carries its proper phase. BerryCPT releases up to commit
! a564a1e omitted the factor i. That omission cancelled in the
! band-diagonal products v_a(n,l) v_b(l,n) used for the Berry curvature
! and OAM and in the effective matrices of contiguous degenerate blocks,
! but not in the off-diagonal elements of the generalized OAM matrix.
!
! Only the rectangular block i <= nbcder is retained. For pairs entirely
! within that block, the reverse element and energy difference are set
! explicitly rather than taken independently from WAVEDER.
dEijks = 0.0_sp

DO ispin = 1, nstot
    DO ikpt = 1, nktot
        DO iband_cder = 1, nbcder
            vijks(:,iband_cder,iband_cder,ikpt,ispin) = &
                CMPLX(0.0_sp,0.0_sp,KIND=sp)

            DO iband = iband_cder + 1, nb
                energy_difference = REAL(ene(iband,ikpt,ispin) - &
                    ene(iband_cder,ikpt,ispin), KIND=sp)

                dEijks(iband_cder,iband,ikpt,ispin) = &
                    energy_difference

                DO idim = 1, 3
                    vijks(idim,iband_cder,iband,ikpt,ispin) = &
                        CMPLX(0.0_sp,1.0_sp,KIND=sp) * &
                        angstrom_to_bohr * energy_difference * &
                        CONJG(vijks(idim,iband_cder,iband,ikpt,ispin))
                END DO

                IF (iband <= nbcder) THEN
                    dEijks(iband,iband_cder,ikpt,ispin) = &
                        -energy_difference
                    vijks(:,iband,iband_cder,ikpt,ispin) = &
                        CONJG(vijks(:,iband_cder,iband,ikpt,ispin))
                END IF
            END DO
        END DO
    END DO
END DO

WRITE(*,'(A)') ' VASP velocity matrix elements were read successfully.'
WRITE(*,'(A,5(I0,1X),A)') &
    ' Shape of vijks = (', SHAPE(vijks), &
    ') (dim_3D, nbands_smaller, nbands_larger, nkpts, nspins)'
WRITE(*,'(A,4(I0,1X),A)') &
    ' Shape of dEijks = (', SHAPE(dEijks), &
    ') (nbands_smaller, nbands_larger, nkpts, nspins)'

END SUBROUTINE read_mommat_pij_vasp


SUBROUTINE validate_waveder_kinds()

! WAVEDER uses 32-bit integers, 64-bit real header values, and
! 32-bit real components in the complex derivative-overlap record.
IF (STORAGE_SIZE(0_int32) /= 32) THEN
    WRITE(*,'(/,A)') ' ERROR: Unsupported integer storage size for WAVEDER'
    WRITE(*,'(A,I0)') ' INTEGER(KIND=int32) storage size [bits] = ', &
        STORAGE_SIZE(0_int32)
    ERROR STOP 1
END IF

IF (STORAGE_SIZE(0.0_dp) /= 64) THEN
    WRITE(*,'(/,A)') ' ERROR: Unsupported double-precision storage size for WAVEDER'
    WRITE(*,'(A,I0)') ' REAL(KIND=dp) storage size [bits] = ', &
        STORAGE_SIZE(0.0_dp)
    ERROR STOP 1
END IF

IF (STORAGE_SIZE(0.0_sp) /= 32) THEN
    WRITE(*,'(/,A)') ' ERROR: Unsupported single-precision storage size for WAVEDER'
    WRITE(*,'(A,I0)') ' REAL(KIND=sp) storage size [bits] = ', &
        STORAGE_SIZE(0.0_sp)
    ERROR STOP 1
END IF

IF (STORAGE_SIZE(CMPLX(0.0_sp,0.0_sp,KIND=sp)) /= 64) THEN
    WRITE(*,'(/,A)') ' ERROR: Unsupported complex storage size for WAVEDER'
    WRITE(*,'(A,I0)') ' COMPLEX(KIND=sp) storage size [bits] = ', &
        STORAGE_SIZE(CMPLX(0.0_sp,0.0_sp,KIND=sp))
    ERROR STOP 1
END IF

END SUBROUTINE validate_waveder_kinds


SUBROUTINE estimate_output_memory( &
        nstot, nktot, nbcder, nb, &
        vijks_bytes, dEijks_bytes)

! WAVEDER dimensions
INTEGER, INTENT(IN) :: nstot, nktot, nbcder, nb

! Memory required by the returned arrays [bytes]
INTEGER(KIND=int64), INTENT(OUT) :: &
    vijks_bytes, &
    dEijks_bytes

! Number of band-pair, k-point, and spin entries
INTEGER(KIND=int64) :: matrix_entries

! Storage required by one array element [bytes]
INTEGER(KIND=int64) :: &
    complex_sp_bytes, &
    real_sp_bytes

matrix_entries = &
    INT(nbcder,KIND=int64) * INT(nb,KIND=int64) * &
    INT(nktot,KIND=int64) * INT(nstot,KIND=int64)

complex_sp_bytes = INT( &
    STORAGE_SIZE(CMPLX(0.0_sp,0.0_sp,KIND=sp)) / 8, &
    KIND=int64)
real_sp_bytes = INT(STORAGE_SIZE(0.0_sp) / 8,KIND=int64)

vijks_bytes = 3_int64 * matrix_entries * complex_sp_bytes
dEijks_bytes = matrix_entries * real_sp_bytes

END SUBROUTINE estimate_output_memory


SUBROUTINE report_output_memory( &
        vijks_bytes, dEijks_bytes)

! Memory required by the returned arrays [bytes]
INTEGER(KIND=int64), INTENT(IN) :: &
    vijks_bytes, &
    dEijks_bytes

! Combined memory required by the returned arrays [bytes]
INTEGER(KIND=int64) :: total_bytes

! Conversion from bytes to gibibytes
REAL(KIND=dp), PARAMETER :: bytes_per_gib = 1024.0_dp**3

total_bytes = vijks_bytes + dEijks_bytes

WRITE(*,'(/,A)') ' VASP velocity-matrix memory estimate'
WRITE(*,'(A,F12.3,A)') &
    ' vijks   = ', REAL(vijks_bytes,KIND=dp) / bytes_per_gib, ' GiB'
WRITE(*,'(A,F12.3,A)') &
    ' dEijks  = ', REAL(dEijks_bytes,KIND=dp) / bytes_per_gib, ' GiB'
WRITE(*,'(A,F12.3,A)') &
    ' combined = ', REAL(total_bytes,KIND=dp) / bytes_per_gib, ' GiB'
WRITE(*,'(A)') &
    ' This estimate excludes eigenvalues, occupations, work arrays, and runtime overhead.'
FLUSH(OUTPUT_UNIT)

END SUBROUTINE report_output_memory


SUBROUTINE allocation_error( &
        fnameinp, array_name, requested_bytes, message)

! File name, array name, and allocation error description
CHARACTER(LEN=*), INTENT(IN) :: fnameinp, array_name, message

! Requested allocation size [bytes]
INTEGER(KIND=int64), INTENT(IN) :: requested_bytes

! Conversion from bytes to gibibytes
REAL(KIND=dp), PARAMETER :: bytes_per_gib = 1024.0_dp**3

WRITE(*,'(/,A)') ' ERROR: Unable to allocate a VASP velocity-matrix array'
WRITE(*,'(A)') ' File = '//TRIM(fnameinp)
WRITE(*,'(A)') ' Array = '//TRIM(array_name)
WRITE(*,'(A,F12.3,A)') &
    ' Requested allocation = ', &
    REAL(requested_bytes,KIND=dp) / bytes_per_gib, ' GiB'
WRITE(*,'(A)') ' '//TRIM(message)
ERROR STOP 1

END SUBROUTINE allocation_error


SUBROUTINE validate_dimensions( &
        fnameinp, nstot, nktot, nbcder, nb, ene)

! File name and WAVEDER metadata
CHARACTER(LEN=*), INTENT(IN) :: fnameinp
INTEGER, INTENT(IN) :: nstot, nktot, nbcder, nb

! Band-, k-point-, and spin-dependent eigenvalues [Ha]
REAL(KIND=dp), INTENT(IN) :: ene(:,:,:)

IF (nstot < 1 .OR. nstot > 2 .OR. nktot <= 0 .OR. nb <= 0 .OR. &
        nbcder <= 0 .OR. nbcder > nb) THEN
    WRITE(*,'(/,A)') ' ERROR: Invalid WAVEDER dimensions were supplied'
    WRITE(*,'(A)') ' File = '//TRIM(fnameinp)
    WRITE(*,'(A,I0)') ' Number of spin channels = ', nstot
    WRITE(*,'(A,I0)') ' Number of k-points = ', nktot
    WRITE(*,'(A,I0)') ' Number of derivative bands = ', nbcder
    WRITE(*,'(A,I0)') ' Number of bands = ', nb
    ERROR STOP 1
END IF

IF (ANY(SHAPE(ene) /= [nb,nktot,nstot])) THEN
    WRITE(*,'(/,A)') ' ERROR: Inconsistent dimensions of ene'
    WRITE(*,'(A,3(I0,1X))') ' Actual shape = ', SHAPE(ene)
    WRITE(*,'(A,3(I0,1X))') &
        ' Expected shape = ', nb, nktot, nstot
    ERROR STOP 1
END IF

END SUBROUTINE validate_dimensions


SUBROUTINE waveder_file_error( &
        fnameinp, message)

! WAVEDER file name and error description
CHARACTER(LEN=*), INTENT(IN) :: fnameinp, message

WRITE(*,'(/,A)') ' ERROR: Unable to read VASP derivative matrix elements'
WRITE(*,'(A)') ' File = '//TRIM(fnameinp)
WRITE(*,'(A)') ' '//TRIM(message)
ERROR STOP 1

END SUBROUTINE waveder_file_error

END MODULE read_mommat_pij_vasp_mod
