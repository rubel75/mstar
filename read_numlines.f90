SUBROUTINE read_numlines(fname, fid, & ! <- args in 
            nltot) ! -> args out 

! Read number of lines in a file

!! Variables

implicit none
CHARACTER(len=256), intent(in) :: &
    fname ! command line imput arguments
INTEGER, intent(in) :: &
    fid ! file ID
INTEGER, intent(out) :: &
    nltot ! number of lines in file

!! Determine the number of rows

OPEN (fid, file = TRIM(fname), status = 'old')
nltot = 0
DO
    READ (1,*, END=10)
    nltot = nltot + 1
END DO
10 CLOSE (fid)

RETURN
END SUBROUTINE read_numlines
