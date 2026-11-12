MODULE eigvz_mod
CONTAINS

SUBROUTINE eigvz(n, H, & ! <- args in 
            EIGV) ! -> args out

! Solve a _real_ matrix eigenvalue problem double precision

USE precision_mod, ONLY: dp
implicit none

!! Variables in-out

INTEGER, intent(in) :: &
    n ! size of H matrixs
COMPLEX(kind=dp), intent(in) :: &
    H(:,:) ! square matrix (size n x n)
REAL(kind=dp), ALLOCATABLE, intent(out) :: &
    EIGV(:) ! eigenvalues

!! Internal variables

INTEGER :: &
    lwork, lrwork, liwork, & ! The size of the work arrays
    info ! status argument
INTEGER, ALLOCATABLE :: &
    iwork(:) ! The size of the work array (lwork>=n)
COMPLEX(kind=dp), ALLOCATABLE :: &
    A(:,:), & ! used for H matrix and then overwritten by eigenvectors
    work(:) ! workspace array
REAL(kind=dp), ALLOCATABLE :: &
    rwork(:) ! workspace array

!! Parameters

CHARACTER(len=1), PARAMETER :: &
    uplo='U', & ! stores the upper triangular part of A
    jobz = 'N' ! do not compute eigenvectors

!! External subroutines

EXTERNAL :: &
    zheevd ! Intel MKL

ALLOCATE( EIGV(n) ) ! allocate intent(out) array, but not deallocate here
ALLOCATE( A(n,n) )

!! Check dimension and reassign H to another variable to avoid it being 
!! overwritten

IF (SIZE(H,1) /= n .OR. SIZE(H,2) /= n) THEN
    WRITE(*,'(A,I0,A,I0,A,I0,A,I0,A)') &
        'Error: input matrix H has shape (', SIZE(H,1), ',', SIZE(H,2), &
        ') but expected (', n, ',', n, ')'
    ERROR STOP 'Invalid matrix dimensions in eigvd'
END IF
A = H

!! Determine size of work arrays

lwork = -1
lrwork = -1
liwork = -1
allocate( work(1), rwork(1), iwork(1) )
call zheevd(jobz, uplo, n, A, n, EIGV, work, lwork, rwork, lrwork, & !...
    iwork, liwork, info) ! double precision (kind=8)
IF (info .ne. 0) THEN
    write (*,*) "ERROR in zheevd (1st call): info = ", info
    write (*,*) "If info = -i, the i-th parameter had an illegal value."
    write (*,*) "If info = i, and jobz = 'V', then the algorithm failed to "//&
        "compute an eigenvalue while working on the submatrix lying in "//&
        "rows and columns info/(n+1) through mod(info,n+1)."
    error stop
END IF
lwork = INT( work(1) )
lrwork = INT( rwork(1) )
liwork = iwork(1)
deallocate ( work, rwork, iwork )

!! Solve eigenvalues problem

allocate( work(lwork), rwork(lrwork), iwork(liwork) )
call zheevd(jobz, uplo, n, A, n, EIGV, work, lwork, rwork, lrwork, & !...
    iwork, liwork, info) ! double precision (kind=8)
deallocate ( work, rwork, iwork, A )
IF (info .ne. 0) THEN
    write (*,*) "ERROR in zheevd (2nd call): info = ", info
    write (*,*) "If info = -i, the i-th parameter had an illegal value."
    write (*,*) "If info = i, and jobz = 'V', then the algorithm failed to "//&
        "compute an eigenvalue while working on the submatrix lying in rows "//&
        "and columns info/(n+1) through mod(info,n+1)."
    error stop
END IF
! eigenvectors are returned through A-matrix are not used

RETURN
END SUBROUTINE eigvz

END MODULE eigvz_mod