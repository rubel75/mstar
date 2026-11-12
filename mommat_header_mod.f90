MODULE mommat_header_mod
! Locate the next k-point record of a WIEN2k case.mommat2[up/dn] file
IMPLICIT NONE
PRIVATE
PUBLIC :: peek_mommat_kpoint

CONTAINS

SUBROUTINE peek_mommat_kpoint(unitinp, fname, & ! in
        found, ikpt_file, nb)                     ! out
! Read forward to the next line containing "KP:", e.g.
!    KP:     1 NEMIN NEMAX :     1  116 dE:  -5.0   5.0 K:
! and return the k-point index and the number of bands NEMAX. The file is
! then positioned back at the start of that line (BACKSPACE), so that the
! complete k-point record can be read by read_mommat_pij_wien2k_kpoint.
! found = .FALSE. when the end of the file is reached without a header.
!
! The band range must start at NEMIN = 1, since band indices are used to
! identify degenerate groups and the mass is computed with all bands.

INTEGER, INTENT(IN) :: unitinp ! connected unit of the mommat file
CHARACTER(LEN=*), INTENT(IN) :: fname ! file name (for messages)
LOGICAL, INTENT(OUT) :: found ! a k-point header was found
INTEGER, INTENT(OUT) :: &
    ikpt_file, & ! k-point index in the file
    nb ! number of bands (NEMAX)

CHARACTER(LEN=256) :: cline, iomsg
INTEGER :: ierr, nbmin

found = .FALSE.
ikpt_file = 0
nb = 0
DO
    READ(unitinp,'(A)',IOSTAT=ierr,IOMSG=iomsg) cline
    IF (ierr < 0) RETURN ! end of file, no further k-point
    IF (ierr > 0) THEN
        WRITE(*,'(A)') ' ERROR: unable to read the file '//TRIM(fname)
        WRITE(*,'(A)') ' '//TRIM(iomsg)
        ERROR STOP 1
    END IF
    IF (INDEX(cline,'KP:') /= 0) EXIT
END DO

READ(cline,'(7X,I5)',IOSTAT=ierr) ikpt_file
IF (ierr == 0) READ(cline,'(27X,2I5)',IOSTAT=ierr) nbmin, nb
IF (ierr /= 0) THEN
    WRITE(*,'(A)') ' ERROR: unable to read the k-point index and band range'// &
        ' from the following line of '//TRIM(fname)
    WRITE(*,'(A)') TRIM(cline)
    WRITE(*,'(A)') ' The line is expected to look like this'
    WRITE(*,'(A)') '   KP:     1 NEMIN NEMAX :     1  116 dE:  -5.0   5.0 K:'
    ERROR STOP 1
END IF
IF (nbmin /= 1 .OR. nb < 2) THEN
    WRITE(*,'(A,I0,A,I0)') ' ERROR: unexpected band range NEMIN = ', nbmin, &
        ', NEMAX = ', nb
    WRITE(*,'(A)') ' in the following line of '//TRIM(fname)
    WRITE(*,'(A)') TRIM(cline)
    WRITE(*,'(A)') ' NEMIN must be 1. NEMIN > 1 usually means that there are'// &
        ' valence states below Emin of the optic window (check case.scf).'
    WRITE(*,'(A)') ' In that case, set a lower Emin in case.inop (e.g., -8 Ry).'
    ERROR STOP 1
END IF
found = .TRUE.
BACKSPACE(unitinp, IOSTAT=ierr)
IF (ierr /= 0) THEN
    WRITE(*,'(A)') ' ERROR: unable to reposition the file '//TRIM(fname)
    ERROR STOP 1
END IF

END SUBROUTINE peek_mommat_kpoint

END MODULE mommat_header_mod
