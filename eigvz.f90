SUBROUTINE eigvz(n, H, & ! <- args in 
            EIGV, EIGF) ! -> args out

! Solve a _complex_ eigenvalue problem for a Hermitian matrix

!! Variables in-out

implicit none
INTEGER, intent(in) :: &
    n ! size of H matrices
COMPLEX(kind=4), intent(in) :: &
    H(n,n) ! square Hermitian matrix
REAL(kind=4), intent(out) :: &
    EIGV(n), & ! eigenvalues
    EIGF(n,n) ! eigenvectors

!! Internal variables

INTEGER :: &
    info, & ! status argument
    lwork, lrwork, liwork, & ! The size of the work arrays
    I, J ! counters
INTEGER, ALLOCATABLE :: &
    iwork(:) ! The size of the work array (lwork>=n)
REAL(kind=4), ALLOCATABLE :: &
    rwork(:) ! workspace array
COMPLEX(kind=4) :: &
    A(n,n) ! used for H matrix and then overwritten by eigenvectors
COMPLEX(kind=4), ALLOCATABLE :: &
    work(:) ! workspace array

!! Parameters

CHARACTER(len=1), PARAMETER :: &
    uplo='U', & ! stores the upper triangular part of A
    jobz = 'V' ! compute eigenvectors too

!! Reassign H to another variable to avoid it being overwritten

A = H

!! Determine size of work arrays

lwork = -1
lrwork = -1
liwork = -1
allocate( work(1), rwork(1), iwork(1) )
call CHEEVD(jobz, uplo, n, A, n, EIGV, work, lwork, rwork, lrwork, & !...
    iwork, liwork, info) ! single precision (kind=4)
IF (info .ne. 0) THEN
    write (*,*) "ERROR in dsyevd: info = ", info
    write (*,*) "If info = -i, the i-th parameter had an illegal value."
    write (*,*) "If info = i, and jobz = 'V', then the algorithm failed to"
    write (*,*) "compute an eigenvalue while working on the submatrix lying"
    write (*,*) "in rows and columns info/(n+1) through mod(info,n+1)."
    stop
END IF
lwork = INT( work(1) )
lrwork = INT( rwork(1) )
liwork = iwork(1)
deallocate ( work, rwork, iwork )

!! Solve eigenvalues problem

allocate( work(lwork), rwork(lrwork), iwork(liwork) )
call CHEEVD(jobz, uplo, n, A, n, EIGV, work, lwork, rwork, lrwork, & !...
    iwork, liwork, info) ! single precision (kind=4)
deallocate ( work, rwork, iwork )
IF (info .ne. 0) THEN
    write (*,*) "ERROR in dsyevd: info = ", info
    write (*,*) "If info = -i, the i-th parameter had an illegal value."
    write (*,*) "If info = i, and jobz = 'V', then the algorithm failed to"
    write (*,*) "compute an eigenvalue while working on the submatrix lying"
    write (*,*) "in rows and columns info/(n+1) through mod(info,n+1)."
    stop
END IF
EIGF = A ! eigenvectors are returned through A-matrix

RETURN
END SUBROUTINE eigvz
