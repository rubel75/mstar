MODULE eigvs_mod
! Eigenvalues of a real symmetric matrix (double precision)
USE precision_mod, ONLY: dp
IMPLICIT NONE
PRIVATE
PUBLIC :: eigvs

CONTAINS

SUBROUTINE eigvs(n, H, & ! <- args in
            EIGV) ! -> args out
! Eigenvalues (ascending) of a real symmetric matrix. Only the upper
! triangle of H is referenced. Eigenvectors are not computed.

!! Variables in-out
INTEGER, INTENT(IN) :: &
    n ! size of H matrix
REAL(KIND=dp), INTENT(IN) :: &
    H(:,:) ! square symmetric matrix (n x n), upper triangle used
REAL(KIND=dp), INTENT(OUT) :: &
    EIGV(:) ! eigenvalues (size n)

!! Internal variables
INTEGER :: &
    lwork, liwork, & ! The size of the work arrays
    info ! status argument
INTEGER, ALLOCATABLE :: iwork(:) ! integer workspace
REAL(KIND=dp), ALLOCATABLE :: &
    A(:,:), & ! copy of H, overwritten by DSYEVD
    work(:) ! workspace array

!! Parameters
CHARACTER(len=1), PARAMETER :: &
    uplo = 'U', & ! use the upper triangular part of A
    jobz = 'N' ! eigenvalues only

!! External subroutines
EXTERNAL :: dsyevd ! LAPACK

IF (SIZE(H,1) /= n .OR. SIZE(H,2) /= n .OR. SIZE(EIGV) /= n) THEN
    WRITE(*,'(A,I0)') ' ERROR in eigvs: inconsistent dimensions for n = ', n
    ERROR STOP 1
END IF
ALLOCATE( A(n,n) )
A = H

!! Determine size of work arrays
lwork = -1
liwork = -1
ALLOCATE( work(1), iwork(1) )
CALL dsyevd(jobz, uplo, n, A, n, EIGV, work, lwork, iwork, liwork, info)
IF (info /= 0) THEN
    WRITE(*,*) "ERROR in dsyevd (workspace query): info = ", info
    ERROR STOP 1
END IF
lwork = INT( work(1) )
liwork = iwork(1)
DEALLOCATE( work, iwork )

!! Solve eigenvalue problem
ALLOCATE( work(lwork), iwork(liwork) )
CALL dsyevd(jobz, uplo, n, A, n, EIGV, work, lwork, iwork, liwork, info)
DEALLOCATE( work, iwork, A )
IF (info /= 0) THEN
    WRITE(*,*) "ERROR in dsyevd: info = ", info
    WRITE(*,*) "If info = -i, the i-th parameter had an illegal value."
    WRITE(*,*) "If info = i, the algorithm failed to converge."
    ERROR STOP 1
END IF

END SUBROUTINE eigvs

END MODULE eigvs_mod
