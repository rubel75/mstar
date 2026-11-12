MODULE dgenen_mod
! Group degenerate states
USE precision_mod, ONLY: sp
IMPLICIT NONE
PRIVATE
PUBLIC :: dgenen

CONTAINS

SUBROUTINE dgenen(nb, nbcder, dEij, dEtol, & ! <- args in
            nmlist, ndblen) ! -> args out
! Group degenerate states among the target bands 1..nbcder.
!
! A group starts at band i and extends over the following bands j for which
! |E_j - E_i| <= dEtol, i.e., all bands of a group lie within dEtol of its
! first band. Only complete groups are returned: if the last group that
! starts at or below nbcder extends beyond nbcder, it is excluded, since
! derivative matrix elements are not available for all of its members.
!
! dEij may use rectangular storage (target band, all bands) with at least
! nbcder rows; only rows 1..nbcder are accessed. Before 2026-11 this routine
! required the full square array.

!! Variables in-out
INTEGER, INTENT(IN) :: &
    nb, & ! number of bands
    nbcder ! number of low energy bands for which the mass will be calculated
REAL(KIND=sp), INTENT(IN) :: &
    dEij(:,:), & ! energy differences E_j-E_i (target band i, all bands j)
    dEtol ! max dE for 2 states to be considered as degenerate
INTEGER, INTENT(OUT) :: &
    nmlist(:,:), & ! list of degeneracy indices, e.g. [1 1; 2 4; ...]
    ndblen ! length of the useful part in 'nmlist' array

!! Variables internal
INTEGER :: &
    i, j ! first and last band of the current group

IF (nbcder < 1 .OR. nbcder > nb .OR. SIZE(dEij,1) < nbcder .OR. &
        SIZE(dEij,2) < nb .OR. SIZE(nmlist,1) < nbcder .OR. &
        SIZE(nmlist,2) /= 2) THEN
    WRITE(*,'(A)') ' ERROR in dgenen: inconsistent array dimensions'
    WRITE(*,'(A,I0,A,I0)') ' nb = ', nb, ', nbcder = ', nbcder
    WRITE(*,'(A,2(I0,1X))') ' SHAPE(dEij) = ', SHAPE(dEij)
    WRITE(*,'(A,2(I0,1X))') ' SHAPE(nmlist) = ', SHAPE(nmlist)
    ERROR STOP 1
END IF

nmlist = 0
ndblen = 0
i = 1
DO WHILE (i <= nbcder)
    ! extend the group as long as the next band is within dEtol of band i
    j = i
    DO WHILE (j < nb)
        IF (ABS(dEij(i,j+1)) > dEtol) EXIT
        j = j + 1
    END DO
    ! exclude a group that is cut by nbcder (and stop)
    IF (j > nbcder) EXIT
    ndblen = ndblen + 1
    nmlist(ndblen,1) = i
    nmlist(ndblen,2) = j
    i = j + 1
END DO

END SUBROUTINE dgenen

END MODULE dgenen_mod
