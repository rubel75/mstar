MODULE read_waveder_header_mod

USE precision_mod, ONLY: sp

IMPLICIT NONE

PRIVATE
PUBLIC :: read_waveder_header

CONTAINS

SUBROUTINE read_waveder_header( &
        fnameinp, &
        nstot, nktot, nbcder, nb)

! Name of the VASP WAVEDER file
CHARACTER(LEN=*), INTENT(IN) :: fnameinp

! Metadata stored in the first WAVEDER record
INTEGER, INTENT(OUT) :: &
    nstot, &  ! total number of spins
    nktot, &  ! total number of k-points
    nbcder, & ! number of bands included in derivative overlaps
    nb         ! total number of bands

! The first WAVEDER record contains four 4-byte integers
INTEGER(KIND=sp) :: &
    xnb, &      ! WDES%NB_TOT
    xnbcder, &  ! NBANDS_CDER
    xnktot, &   ! WDES%NKPTS
    xnstot      ! WDES%ISPIN

! Input unit and I/O status
INTEGER :: unitinp, ierr
CHARACTER(LEN=256) :: iomsg

! INTEGER(KIND=sp) must occupy four bytes to match the WAVEDER header
IF (STORAGE_SIZE(xnb) /= 32) THEN
    WRITE(*,'(/,A)') ' ERROR: INTEGER(KIND=sp) does not occupy four bytes'
    ERROR STOP 1
END IF

OPEN(NEWUNIT=unitinp, FILE=TRIM(fnameinp), FORM='UNFORMATTED', &
    ACCESS='SEQUENTIAL', STATUS='OLD', ACTION='READ', &
    IOSTAT=ierr, IOMSG=iomsg)
IF (ierr /= 0) THEN
    WRITE(*,'(/,A)') ' ERROR: Unable to open the WAVEDER file '//TRIM(fnameinp)
    WRITE(*,'(A)') ' '//TRIM(iomsg)
    ERROR STOP 1
END IF

READ(unitinp, IOSTAT=ierr, IOMSG=iomsg) xnb, xnbcder, xnktot, xnstot
IF (ierr /= 0) THEN
    CLOSE(unitinp)
    WRITE(*,'(/,A)') ' ERROR: Unable to read the WAVEDER header from '// &
        TRIM(fnameinp)
    WRITE(*,'(A)') ' '//TRIM(iomsg)
    ERROR STOP 1
END IF

CLOSE(unitinp, IOSTAT=ierr, IOMSG=iomsg)
IF (ierr /= 0) THEN
    WRITE(*,'(/,A)') ' ERROR: Unable to close the WAVEDER file '//TRIM(fnameinp)
    WRITE(*,'(A)') ' '//TRIM(iomsg)
    ERROR STOP 1
END IF

! Convert the four-byte header integers to the default integer kind
nb = xnb
nbcder = xnbcder
nktot = xnktot
nstot = xnstot

! Check that the header values are physically meaningful
IF (nb <= 0 .OR. nbcder <= 0 .OR. nktot <= 0 .OR. &
        nstot < 1 .OR. nstot > 2 .OR. nbcder > nb) THEN
    WRITE(*,'(/,A)') ' ERROR: Invalid metadata were read from the WAVEDER header'
    WRITE(*,'(A,I0)') ' NB_TOT = ', nb
    WRITE(*,'(A,I0)') ' NBANDS_CDER = ', nbcder
    WRITE(*,'(A,I0)') ' NKPTS = ', nktot
    WRITE(*,'(A,I0)') ' ISPIN = ', nstot
    ERROR STOP 1
ELSE
    WRITE(*,'(/,A)') ' Valid metadata were read from the WAVEDER header'
    WRITE(*,'(A,I0)') ' NB_TOT = ', nb
    WRITE(*,'(A,I0)') ' NBANDS_CDER = ', nbcder
    WRITE(*,'(A,I0)') ' NKPTS = ', nktot
    WRITE(*,'(A,I0)') ' ISPIN = ', nstot
END IF

END SUBROUTINE read_waveder_header

END MODULE read_waveder_header_mod
