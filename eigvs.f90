SUBROUTINE eigvs(n, H, & ! <- args in 
            EIGV, EIGF) ! -> args out

! Solve a _real_ matrix eigenvalue problem

!! Variables in-out

implicit none
INTEGER, intent(in) :: &
    n ! size of H matrixs
REAL(kind=4), intent(in) :: &
    H(n,n) ! square matrix
REAL(kind=4), intent(out) :: &
    EIGV(n), & ! eigenvalues
    EIGF(n,n) ! eigenvectors

!! Internal variables

INTEGER :: &
    lwork, liwork, & ! The size of the work arrays
    info ! status argument
INTEGER, ALLOCATABLE :: iwork(:) ! The size of the work array (lwork>=n)
REAL(kind=4) :: &
    A(n,n) ! used for H matrix and then overwritten by eigenvectors
REAL(kind=4), ALLOCATABLE :: &
    work(:) ! workspace array

!! Parameters

CHARACTER(len=1), PARAMETER :: &
    uplo='U', & ! stores the upper triangular part of A
    jobz = 'V' ! compute eigenvectors too

!! Reassign H to another variable to avoid it being overwritten

A = H

!! Determine size of work arrays

lwork = -1
liwork = -1
allocate( work(1), iwork(1) )
!call dsyevd(jobz, uplo, n, A, n, EIGV, work, lwork, iwork, liwork, info) ! double precision (kind=8)
call ssyevd(jobz, uplo, n, A, n, EIGV, work, lwork, iwork, liwork, info) ! single precision (kind=4)
IF (info .ne. 0) THEN
    write (*,*) "ERROR in ssyevd (1st call): info = ", info
    write (*,*) "If info = -i, the i-th parameter had an illegal value."
    write (*,*) "If info = i, and jobz = 'V', then the algorithm failed to "//&
        "compute an eigenvalue while working on the submatrix lying in "//&
        "rows and columns info/(n+1) through mod(info,n+1)."
    stop
END IF
lwork = INT( work(1) )
liwork = iwork(1)
deallocate ( work, iwork )

!! Solve eigenvalues problem

allocate( work(lwork), iwork(liwork) )
!call dsyevd(jobz, uplo, n, A, n, EIGV, work, lwork, iwork, liwork, info) ! double precision (kind=8)
call ssyevd(jobz, uplo, n, A, n, EIGV, work, lwork, iwork, liwork, info) ! single precision (kind=4)
deallocate ( work, iwork )
IF (info .ne. 0) THEN
    write (*,*) "ERROR in ssyevd (2nd call): info = ", info
    write (*,*) "If info = -i, the i-th parameter had an illegal value."
    write (*,*) "If info = i, and jobz = 'V', then the algorithm failed to "//&
        "compute an eigenvalue while working on the submatrix lying in rows "//&
        "and columns info/(n+1) through mod(info,n+1)."
    stop
END IF
EIGF = A ! eigenvectors are returned through A-matrix

RETURN
END SUBROUTINE eigvs
