MODULE degenmass_mod
! Inverse effective mass for a block of degenerate bands (degenerate k.p
! perturbation theory, Eqs. (6) and (7) of Comp. Phys. Commun. 261, 107648
! (2021)). The structure follows degenbc_mod and degenoam_mod of BerryCPT.
USE precision_mod, ONLY: sp, dp
USE eigvz_mod, ONLY: eigvz
USE is_hermitian_mod, ONLY: is_hermitian
USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
IMPLICIT NONE
PRIVATE
PUBLIC :: degenmass

CONTAINS

SUBROUTINE degenmass( &
    idg1, idg2, alpha, beta, pij, dEij, & ! in
    minv)                                 ! out
! Effective-matrix contribution to the inverse effective mass m0/m*_ab for
! the degenerate block D = {idg1, ..., idg2}. The free-electron term
! delta_ab is NOT included; it is added by the caller.
!
! Array conventions (Hartree atomic units, hbar = m0 = 1)
! -------------------------------------------------------
! pij(a,i,l) = <i|p_a|l>, rectangular storage (component, target band,
!              all bands). Target bands run 1..SIZE(pij,2) >= idg2,
!              transition bands run 1..SIZE(pij,3) = nb.
! dEij(i,l)  = E_l - E_i [Ha], same rectangular storage.
! The reverse element <l|p_b|j> is reconstructed as CONJG(<j|p_b|l>)
! (p is Hermitian).
!
! Effective matrix
! ----------------
! For i, j in D,
!
!   M_ij = SUM_{l not in D} (1/2) [ p_a(i,l) p_b(l,j) + p_b(i,l) p_a(l,j) ]
!          x [ 1/(E_i - E_l) + 1/(E_j - E_l) ] .
!
! For i = j this is exactly the non-degenerate expression, Eq. (5),
! SUM_l [p_a(i,l) p_b(l,i) + p_b(i,l) p_a(l,i)] / (E_i - E_l).
! For E_i = E_j it is Eq. (6) of the paper. Versions of mstar before
! 2026-11 used 1/[(E_i - E_l + E_j - E_l)/2] (inverse of the averaged gap)
! instead of the averaged inverse gap. The two coincide for exact
! degeneracy and differ only for numerically near-degenerate states. The
! averaged inverse gap keeps M Hermitian and matches BerryCPT (degenoam_mod).
! The inverse masses of the block are the eigenvalues of M, Eq. (7).
!
! Precision
! ---------
! Matrix elements and energy differences are stored in sp. Each factor is
! promoted to dp before the products are formed. Every term is accumulated
! with complex Kahan compensated summation in dp. M and its eigenvalues are
! dp. The result does not depend on the summation order of a parallel run.

INTEGER, INTENT(IN) :: &
    idg1, idg2, & ! first and last band of the degenerate block
    alpha, beta   ! Cartesian directions (1, 2, 3 = x, y, z)
COMPLEX(KIND=sp), INTENT(IN) :: &
    pij(:,:,:)    ! <i|p_a|l> [a.u.] (component, target band, all bands)
REAL(KIND=sp), INTENT(IN) :: &
    dEij(:,:)     ! E_l - E_i [Ha] (target band, all bands)
REAL(KIND=dp), ALLOCATABLE, INTENT(OUT) :: &
    minv(:)       ! eigenvalues of M (m0/m* without delta_ab)

COMPLEX(KIND=dp), ALLOCATABLE :: &
    M(:,:), &     ! effective matrix
    Mcorr(:,:)    ! Kahan compensation
COMPLEX(KIND=dp) :: &
    pa_il, pb_il, pa_jl, pb_jl, & ! matrix elements promoted to dp
    p2, &         ! p_a(i,l) p_b(l,j) + p_b(i,l) p_a(l,j)
    term, &       ! contribution of one intermediate state l
    corrected_term, temp ! Kahan intermediates
REAL(KIND=dp) :: &
    dEi, dEj, &   ! E_i - E_l and E_j - E_l [Ha]
    inverse_gap   ! averaged inverse gap [1/Ha]
INTEGER :: &
    nb, ndg, i, j, ib, jb, l
CHARACTER(LEN=*), PARAMETER :: &
    wformat_m = '(*( :, "(", ES12.4E3, ",", ES12.4E3, ")", : , 1X ))'

nb = SIZE(pij,3)
ndg = idg2 - idg1 + 1

!! Check the interface
IF (idg1 < 1 .OR. idg2 < idg1 .OR. idg2 > SIZE(pij,2)) THEN
    WRITE(*,'(A,I0,A,I0,A,I0)') ' ERROR in degenmass: block ', idg1, ':', &
        idg2, ' is outside the target range 1:', SIZE(pij,2)
    ERROR STOP 1
END IF
IF (SIZE(pij,1) /= 3 .OR. SIZE(dEij,1) /= SIZE(pij,2) .OR. &
        SIZE(dEij,2) /= nb) THEN
    WRITE(*,'(A)') ' ERROR in degenmass: inconsistent shapes of pij and dEij'
    WRITE(*,'(A,3(I0,1X))') ' SHAPE(pij)  = ', SHAPE(pij)
    WRITE(*,'(A,2(I0,1X))') ' SHAPE(dEij) = ', SHAPE(dEij)
    ERROR STOP 1
END IF
IF (alpha < 1 .OR. alpha > 3 .OR. beta < 1 .OR. beta > 3) THEN
    WRITE(*,'(A,I0,1X,I0)') ' ERROR in degenmass: invalid alpha, beta = ', &
        alpha, beta
    ERROR STOP 1
END IF

ALLOCATE( M(ndg,ndg), Mcorr(ndg,ndg) )
M = (0.0_dp, 0.0_dp)
Mcorr = (0.0_dp, 0.0_dp)

!! Construct the upper triangle of M
DO i = 1, ndg
    ib = idg1 + i - 1
    DO j = i, ndg
        jb = idg1 + j - 1
        DO l = 1, nb
            IF (l >= idg1 .AND. l <= idg2) CYCLE ! exclude the block
            pa_il = CMPLX(pij(alpha,ib,l), KIND=dp)
            pb_il = CMPLX(pij(beta,ib,l), KIND=dp)
            pa_jl = CMPLX(pij(alpha,jb,l), KIND=dp)
            pb_jl = CMPLX(pij(beta,jb,l), KIND=dp)
            ! <l|p_b|j> = CONJG(<j|p_b|l>), <l|p_a|j> = CONJG(<j|p_a|l>)
            p2 = pa_il*CONJG(pb_jl) + pb_il*CONJG(pa_jl)
            dEi = -REAL(dEij(ib,l), KIND=dp) ! E_i - E_l
            dEj = -REAL(dEij(jb,l), KIND=dp) ! E_j - E_l
            IF (i == j) THEN
                inverse_gap = 1.0_dp/dEi
            ELSE
                inverse_gap = 0.5_dp*(1.0_dp/dEi + 1.0_dp/dEj)
            END IF
            term = p2*inverse_gap
            IF (.NOT. IEEE_IS_FINITE(REAL(term,dp)) .OR. &
                    .NOT. IEEE_IS_FINITE(AIMAG(term))) THEN
                WRITE(*,'(A)') ' ERROR in degenmass: non-finite contribution'
                WRITE(*,'(A,I0,A,I0,A,I0)') ' bands i = ', ib, ', j = ', jb, &
                    ', l = ', l
                WRITE(*,'(A,I0,1X,I0)') ' alpha, beta = ', alpha, beta
                WRITE(*,*) 'p_a(i,l) = ', pa_il, ' p_b(i,l) = ', pb_il
                WRITE(*,*) 'p_a(j,l) = ', pa_jl, ' p_b(j,l) = ', pb_jl
                WRITE(*,*) 'E_i - E_l = ', dEi, ' E_j - E_l = ', dEj
                WRITE(*,'(A)') ' A zero energy difference outside a '// &
                    'degenerate block usually means that the degeneracy '// &
                    'tolerance is too small.'
                ERROR STOP 1
            END IF
            ! Kahan compensated summation into M(i,j)
            corrected_term = term - Mcorr(i,j)
            temp = M(i,j) + corrected_term
            Mcorr(i,j) = (temp - M(i,j)) - corrected_term
            M(i,j) = temp
        END DO ! l
    END DO ! j
    ! the diagonal is real in exact arithmetic
    M(i,i) = CMPLX(REAL(M(i,i),dp), 0.0_dp, KIND=dp)
END DO ! i
DO i = 1, ndg
    DO j = i+1, ndg
        M(j,i) = CONJG(M(i,j))
    END DO
END DO

!! Inverse masses are the eigenvalues of M
IF (ndg > 1) THEN
    IF (.NOT. is_hermitian(M)) THEN
        WRITE(*,'(A)') ' M ='
        DO i = 1, ndg
            WRITE(*,wformat_m) (M(i,j), j = 1, ndg)
        END DO
        ERROR STOP 'in degenmass: matrix M is not Hermitian'
    END IF
    CALL eigvz(ndg, M, & ! in
        minv)            ! out (allocated inside)
ELSE
    ALLOCATE( minv(1) )
    minv(1) = REAL(M(1,1), dp)
END IF
DEALLOCATE( M, Mcorr )

END SUBROUTINE degenmass

END MODULE degenmass_mod
