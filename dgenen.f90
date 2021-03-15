SUBROUTINE dgenen(nb, nbcder, dEij, dEtol, & ! <- args in 
            nmlist, ndblen) ! -> args out 

! Group degenerate states

!! Variables in-out

implicit none
INTEGER, intent(in) :: &
    nb, & ! number of bands
    nbcder ! number of low energy bands for which the mass will be calculated
REAL(kind=4), intent(in) :: &
    dEij(nb,nb), & ! energy differences E_i-E_j
    dEtol ! max dE for 2 states to be considered as degenerate
INTEGER, intent(out) :: &
    nmlist(nb,2), & ! list of degeneracy indices, e.g. [1 1; 2 4; ...]
    ndblen ! length of the useful part in 'nmlist' array

!! Variables internal

INTEGER :: &
    i,j, & ! counters
    uboundnmlist, & ! length of array nmlist
    nmlocal(2*nb) ! 1D array of list of degeneracy indices
REAL(kind=4) :: &
    dE ! energy difference

!! create list of degenerate bands (1,3,3,7,...,nb)

i = 1
uboundnmlist = 0
loop1: DO WHILE (.TRUE.)
    loop2: DO j = i+1,nb
        dE = abs(dEij(i,j))
        IF (dE > dEtol) THEN
            nmlocal(uboundnmlist+1) = i
            nmlocal(uboundnmlist+2) = j-1
            uboundnmlist = uboundnmlist + 2
            i = j
            exit loop2
        ELSEIF (j==nb) THEN
            nmlocal(uboundnmlist+1) = i
            nmlocal(uboundnmlist+2) = j
            uboundnmlist = uboundnmlist + 2
            i = j
            exit loop2
        END IF
    END DO loop2
    IF (j == nb) THEN
        IF (nmlocal(uboundnmlist)==nb) THEN
            exit loop1
        ELSE
            nmlocal(uboundnmlist+1) = i
            nmlocal(uboundnmlist+2) = j
            uboundnmlist = uboundnmlist + 2
            exit loop1
        END IF
    END IF
END DO loop1

!! reshape the list
!! (1,3)
!! (3,7) ...

nmlist = 0 ! initialize with 0's
DO i = 2, uboundnmlist, 2
    nmlist(i/2,1) = nmlocal(i-1)
    nmlist(i/2,2) = nmlocal(i)
END DO
ndblen = uboundnmlist/2

!! truncate the list at nbcder

IF (nbcder < nb) THEN ! only if nb > nbcder
    loop3: DO i = 1, ndblen
        ! locate the degeneracy range that includes nbcder
        IF ( nmlist(i,1)<=nbcder .AND. nmlist(i,2)>=nbcder ) THEN
            IF ( nmlist(i,2)==nbcder ) THEN ! the range ends with nbcder
                ndblen = i ! stop degeneracy list at this point
            ELSE ! nbcder is within the range, but not at the end
                ndblen = i-1 ! exclude the current range
            END IF
            exit loop3
        END IF
    END DO loop3
END IF

RETURN
END SUBROUTINE dgenen
