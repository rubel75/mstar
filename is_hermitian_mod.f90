module is_hermitian_mod
contains
  
pure logical function is_hermitian(A) result(ok)
    USE precision_mod, ONLY: dp
    IMPLICIT NONE

    complex(dp), intent(in) :: A(:,:)

    integer  :: n, i, j
    real(dp) :: eps, scale, tol, err

    n = size(A,1)
    if (n /= size(A,2)) then
        ok = .false.; return
    end if

    ! Scale by the largest magnitude entry to make the test relative.
    scale = max(1.0_dp, maxval(abs(A)))
    eps   = epsilon(1.0_dp)

    ! Heuristic tolerance: O(eps * scale), inflated a bit for safety.
    tol = 20.0_dp * eps * scale

    err = 0.0_dp

    ! Check diagonals have negligible imaginary part
    do i = 1, n
        err = max(err, abs(aimag(A(i,i))))
        if (err > tol) then
            ok = .false.; return
        end if
    end do

    ! Check upper triangle against conjugate of lower
    do j = 1, n
        do i = j+1, n
            err = max(err, abs(A(i,j) - conjg(A(j,i))))
            if (err > tol) then
                ok = .false.; return
            end if
        end do
    end do

    ok = .true.
end function is_hermitian
  
end module is_hermitian_mod
  