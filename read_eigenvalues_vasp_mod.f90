MODULE read_eigenvalues_vasp_mod

! Copied from BerryCPT (https://github.com/rubel75/berrycpt, commit 993047e)
! with one modification for mstar: eigenvalues are returned in double
! precision. EIGENVAL gives eigenvalues in eV with 6 decimals. In single
! precision, the spacing of representable numbers exceeds 1e-6 eV for
! |E| > 8 eV (about 4e-6 eV at 40 eV), so the last printed digit is lost.
! In the CoSi test, single precision changes m0/m_c by a median relative
! amount of 4e-6 (99th percentile 6e-4), about 2-3 times the change caused
! by perturbing the eigenvalues within their printed resolution of 1e-6 eV
! (median 2e-6, 99th percentile 2e-4). The largest changes (several percent)
! occur where the trace results from cancelling components of order 100 and
! are limited by the resolution of EIGENVAL itself.

USE precision_mod, ONLY: sp, dp

IMPLICIT NONE

PRIVATE
PUBLIC :: read_eigenvalues_vasp


CONTAINS

SUBROUTINE read_eigenvalues_vasp( &
        fnameene, &
        nstot, nktot, nb, &
        nbk, ene, occup)

! Name of the VASP EIGENVAL file
CHARACTER(LEN=*), INTENT(IN) :: fnameene

! Metadata obtained from the WAVEDER header
INTEGER, INTENT(IN) :: nstot, nktot, nb

! Number of eigenvalues at each k-point and spin
INTEGER, ALLOCATABLE, INTENT(OUT) :: nbk(:,:)

! Band-, k-point-, and spin-dependent eigenvalues [Ha]
REAL(KIND=dp), ALLOCATABLE, INTENT(OUT) :: ene(:,:,:)

! Band-, k-point-, and spin-dependent occupations
REAL(KIND=sp), ALLOCATABLE, INTENT(OUT) :: occup(:,:,:)

! Text from the current record and the band-index field
CHARACTER(LEN=1024) :: cline
CHARACTER(LEN=256) :: iband_token

! Input unit, metadata, counters, and I/O status
INTEGER :: &
    unitene, &
    nions1, nions2, nions3, nstot_file, &
    nktot_file, nb_file, &
    ikpt, iband, iband_file, i, &
    ierr, ierr_token

! Electron count and k-point metadata
REAL(KIND=sp) :: &
    nelect, &
    kvec(3), kweight

! I/O error message
CHARACTER(LEN=256) :: iomsg

IF (nstot < 1 .OR. nstot > 2 .OR. nktot <= 0 .OR. nb <= 0) THEN
    CALL eigenvalue_file_error( &
        fnameene, & ! in
        'Invalid VASP metadata were supplied') ! in
END IF

CALL open_energy_file( &
    fnameene, & ! in
    unitene)    ! out

! The first EIGENVAL record contains ISPIN
READ(unitene,'(A)',IOSTAT=ierr,IOMSG=iomsg) cline
IF (ierr /= 0) THEN
    CLOSE(unitene)
    CALL eigenvalue_file_error( &
        fnameene, & ! in
        'Unable to read the first EIGENVAL record: '//TRIM(iomsg)) ! in
END IF

READ(cline,*,IOSTAT=ierr) nions1, nions2, nions3, nstot_file
IF (ierr /= 0) THEN
    CLOSE(unitene)
    CALL eigenvalue_file_error( &
        fnameene, & ! in
        'Unable to interpret the first EIGENVAL record') ! in
END IF

IF (nstot_file /= nstot) THEN
    CLOSE(unitene)
    WRITE(*,'(/,A)') ' ERROR: Inconsistent WAVEDER and EIGENVAL metadata'
    WRITE(*,'(A,I0,A,I0)') &
        ' Number of spin channels: WAVEDER = ', nstot, &
        ', EIGENVAL = ', nstot_file
    ERROR STOP 1
END IF

! Skip EIGENVAL records two through five
DO i = 2, 5
    READ(unitene,'(A)',IOSTAT=ierr,IOMSG=iomsg) cline
    IF (ierr /= 0) THEN
        CLOSE(unitene)
        CALL eigenvalue_file_error( &
            fnameene, & ! in
            'Unable to read the EIGENVAL header: '//TRIM(iomsg)) ! in
    END IF
END DO

! The sixth EIGENVAL record contains NELECT, NKPTS, and NBANDS
READ(unitene,'(A)',IOSTAT=ierr,IOMSG=iomsg) cline
IF (ierr /= 0) THEN
    CLOSE(unitene)
    CALL eigenvalue_file_error( &
        fnameene, & ! in
        'Unable to read the EIGENVAL metadata: '//TRIM(iomsg)) ! in
END IF

READ(cline,*,IOSTAT=ierr) nelect, nktot_file, nb_file
IF (ierr /= 0) THEN
    CLOSE(unitene)
    CALL eigenvalue_file_error( &
        fnameene, & ! in
        'Unable to interpret the sixth EIGENVAL record') ! in
END IF

IF (nktot_file /= nktot .OR. nb_file /= nb) THEN
    CLOSE(unitene)
    WRITE(*,'(/,A)') ' ERROR: Inconsistent WAVEDER and EIGENVAL metadata'
    WRITE(*,'(A,I0,A,I0)') &
        ' Number of k-points: WAVEDER = ', nktot, &
        ', EIGENVAL = ', nktot_file
    WRITE(*,'(A,I0,A,I0)') &
        ' Number of bands: WAVEDER = ', nb, &
        ', EIGENVAL = ', nb_file
    ERROR STOP 1
END IF

ALLOCATE(nbk(nktot,nstot))
ALLOCATE(ene(nb,nktot,nstot))
ALLOCATE(occup(nb,nktot,nstot))

nbk = nb
ene = 0.0_dp
occup = 0.0_sp

DO ikpt = 1, nktot
    CALL read_next_nonempty_line( &
        unitene, fnameene, & ! in
        cline)               ! out

    READ(cline,*,IOSTAT=ierr) kvec, kweight
    IF (ierr /= 0) THEN
        CLOSE(unitene)
        CALL eigenvalue_file_error( &
            fnameene, & ! in
            'Unable to interpret a VASP k-point record') ! in
    END IF

    DO iband = 1, nb
        READ(unitene,'(A)',IOSTAT=ierr,IOMSG=iomsg) cline
        IF (ierr /= 0) THEN
            CLOSE(unitene)
            CALL eigenvalue_file_error( &
                fnameene, & ! in
                'Unable to read a VASP eigenvalue record: '//TRIM(iomsg)) ! in
        END IF

        IF (nstot == 1) THEN
            READ(cline,*,IOSTAT=ierr) &
                iband_token, ene(iband,ikpt,1), occup(iband,ikpt,1)
        ELSE
            READ(cline,*,IOSTAT=ierr) &
                iband_token, ene(iband,ikpt,1), ene(iband,ikpt,2), &
                occup(iband,ikpt,1), occup(iband,ikpt,2)
        END IF

        IF (ierr /= 0) THEN
            CLOSE(unitene)
            WRITE(*,'(/,A)') &
                ' ERROR: Unable to interpret a VASP eigenvalue record'
            WRITE(*,'(A)') ' '//TRIM(cline)
            ERROR STOP 1
        END IF

        ! A numeric band index must agree with the record position.
        ! A field containing only asterisks is accepted for large band indices.
        READ(iband_token,*,IOSTAT=ierr_token) iband_file

        IF (ierr_token == 0) THEN
            IF (iband_file /= iband) THEN
                CLOSE(unitene)
                WRITE(*,'(/,A)') ' ERROR: Inconsistent VASP band index'
                WRITE(*,'(A,I0,A,A)') &
                    ' Expected band ', iband, ', found ', TRIM(iband_token)
                ERROR STOP 1
            END IF
        ELSE IF (LEN_TRIM(iband_token) == 0 .OR. &
                    VERIFY(TRIM(iband_token),'*') /= 0) THEN
            CLOSE(unitene)
            WRITE(*,'(/,A)') ' ERROR: Invalid VASP band-index field'
            WRITE(*,'(A)') ' '//TRIM(cline)
            ERROR STOP 1
        END IF
    END DO
END DO

! Convert VASP eigenvalues from eV to Hartree
ene = ene / 27.211386245988_dp ! Hartree energy in eV (CODATA 2018)

CLOSE(unitene,IOSTAT=ierr,IOMSG=iomsg)

IF (ierr /= 0) THEN
    CALL eigenvalue_file_error( &
        fnameene, & ! in
        'Unable to close the energy file: '//TRIM(iomsg)) ! in
END IF

WRITE(*,'(A)') ' VASP eigenvalues and occupations were read successfully.'
WRITE(*,'(A,I0)') ' Number of spin channels = ', nstot
WRITE(*,'(A,I0)') ' Number of k-points = ', nktot
WRITE(*,'(A,I0)') ' Number of bands = ', nb

END SUBROUTINE read_eigenvalues_vasp


SUBROUTINE read_next_nonempty_line( &
        unitene, fnameene, &
        cline)

! Input unit and file name
INTEGER, INTENT(IN) :: unitene
CHARACTER(LEN=*), INTENT(IN) :: fnameene

! First nonempty record
CHARACTER(LEN=*), INTENT(OUT) :: cline

! I/O status and message
INTEGER :: ierr
CHARACTER(LEN=256) :: iomsg

DO
    READ(unitene,'(A)',IOSTAT=ierr,IOMSG=iomsg) cline

    IF (ierr /= 0) THEN
        CALL eigenvalue_file_error( &
            fnameene, & ! in
            'Unable to read the next record: '//TRIM(iomsg)) ! in
    END IF

    IF (LEN_TRIM(cline) > 0) EXIT
END DO

END SUBROUTINE read_next_nonempty_line


SUBROUTINE open_energy_file( &
        fnameene, &
        unitene)

! Name of the energy file
CHARACTER(LEN=*), INTENT(IN) :: fnameene

! Input unit
INTEGER, INTENT(OUT) :: unitene

! I/O status and message
INTEGER :: ierr
CHARACTER(LEN=256) :: iomsg

OPEN(NEWUNIT=unitene,FILE=TRIM(fnameene),FORM='FORMATTED', &
    STATUS='OLD',ACTION='READ',IOSTAT=ierr,IOMSG=iomsg)

IF (ierr /= 0) THEN
    CALL eigenvalue_file_error( &
        fnameene, & ! in
        'Unable to open the energy file: '//TRIM(iomsg)) ! in
END IF

END SUBROUTINE open_energy_file


SUBROUTINE eigenvalue_file_error( &
        fnameene, message)

! Name of the energy file and description of the error
CHARACTER(LEN=*), INTENT(IN) :: fnameene, message

WRITE(*,'(/,A)') ' ERROR: '//TRIM(message)
WRITE(*,'(A)') ' Energy file: '//TRIM(fnameene)
ERROR STOP 1

END SUBROUTINE eigenvalue_file_error

END MODULE read_eigenvalues_vasp_mod