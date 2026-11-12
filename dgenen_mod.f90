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
! Grouping rule (the same as in BerryCPT, find_degenerate_groups_mod.f90):
! bands are ordered by energy, and adjacent bands j and j+1 belong to the
! same group when E_{j+1} - E_j <= dEtol. Groups are therefore chains of
! adjacent bands, and the total width of a group can exceed dEtol.
! (Versions before 2026-11 required every band of a group to lie within
! dEtol of the first band of the group.)
!
! Only complete groups are returned: if the last group that starts at or
! below nbcder extends beyond nbcder, it is excluded, since derivative
! matrix elements are not available for all of its members.
!
! dEij may use rectangular storage (target band, all bands) with at least
! nbcder rows. Only the adjacent-band differences dEij(j,j+1), j <= nbcder,
! are accessed.

!! Variables in-out
INTEGER, INTENT(IN) :: &
    nb, & ! number of bands
    nbcder ! number of low energy bands for which the mass will be calculated
REAL(KIND=sp), INTENT(IN) :: &
    dEij(:,:), & ! energy differences E_j-E_i (target band i, all bands j)
    dEtol ! max dE for 2 adjacent states to be considered as degenerate
INTEGER, INTENT(OUT) :: &
    nmlist(:,:), & ! list of degeneracy indices, e.g. [1 1; 2 4; ...]
    ndblen ! length of the useful part in 'nmlist' array

!! Variables internal
INTEGER :: &
    i, j ! first and last band of the current group
REAL(KIND=sp) :: &
    gap ! E_{j+1} - E_j

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
    ! extend the chain while the next band is within dEtol of the last one
    j = i
    DO WHILE (j < nb)
        gap = dEij(j,j+1)
        IF (gap < -dEtol) THEN
            WRITE(*,'(A)') ' ERROR in dgenen: eigenvalues are not ordered '// &
                'by band index'
            WRITE(*,'(A,I0,A,I0,A,ES14.6E3,A)') ' E(', j+1, ') - E(', j, &
                ') = ', gap, ' Ha'
            ERROR STOP 1
        END IF
        IF (gap > dEtol) EXIT
        j = j + 1
        ! rows above nbcder are not available (and not needed): a chain
        ! that reaches beyond nbcder is excluded below
        IF (j > nbcder) EXIT
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
