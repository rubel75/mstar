SUBROUTINE read_numlines_vasp(fname, fid, & ! <- args in 
                nstot, nktot, nbcder, nb) ! -> args out

! Read header of the WAVEDER file

!! Variables in-out

implicit none
CHARACTER(len=256), intent(in) :: &
    fname ! command line input arguments
INTEGER, intent(in) :: &
    fid ! file ID
INTEGER, intent(out) :: &
    nktot, & ! total number of k-points in WAVEDER file
    nstot, & ! total number of spins in WAVEDER file
    nbcder, & ! smaller number of bands in WAVEDER file
    nb ! total number of bands in VASP calculation

!! Variables internal

INTEGER(kind=4) :: &
    xnb, xnbcder, xnktot, xnstot ! temp variable

!! Read the head record of WAVEDER file to determine 
!! number of bands, k-points, spins

OPEN (fid, file = TRIM(fname), form = 'unformatted', status = 'old')
READ(fid) xnb, xnbcder, xnktot, xnstot
CLOSE(fid)

!! convert INTEGER(kind=4) into INTEGER

nb=xnb
nbcder=xnbcder
nktot=xnktot
nstot=xnstot

RETURN
END SUBROUTINE read_numlines_vasp
