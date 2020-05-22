FUNCTION kdelta(i,j)

! Kronecker delta function

implicit none
INTEGER :: i,j,kdelta
kdelta=0
IF (i == j) kdelta = 1

RETURN
END FUNCTION kdelta
