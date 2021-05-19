SUBROUTINE read_mommat_nb(fid, & ! <- args in
            iline, & ! <-> args in-out
            nb) ! -> args out 

! Determine number of bands in a k-point block of mommat file

!! Variables in-out

implicit none
INTEGER, intent(in) :: &
    fid ! file ID
INTEGER, intent(inout) :: &
    iline ! current line number during reading of mommat file
INTEGER, intent(out) :: &
    nb ! number of bands

!! Variables internal

CHARACTER(len=256) :: &
    cline, & ! text in the current line
    dum1, dum2 ! dummy variables required for reading
INTEGER :: &
    nbmin, nbmax ! number of bands range
LOGICAL :: &
    ldone ! done reading?

!! Read line-by-line trying to find "KP:" keyword

ldone = .false.
DO WHILE ( .not.(ldone) )
    READ(fid,'(A)',END=10) cline
    iline = iline + 1 ! count total lines read since opening the file
    IF (INDEX(cline,'KP:') .ne. 0) THEN ! check for "KP:" in the line
        ! extract the range of bands
        READ(cline,'(A27,2I5,A)',ERR=20) dum1, nbmin, nbmax, dum2
        ! check for the range to start with 1st band
        IF (nbmin == 1 .and. nbmax > nbmin) THEN
            nb = nbmax
            RETURN ! done
        ELSE ! read or data error, abort
            write(*,'(A,I0,A)') 'nbmin = ', nbmin, ' while expected = 1'
            write(*,'(A,I0,A)') 'nbmax = ', nbmax, ' expected nbmax > nbmin'
            GOTO 20 
        END IF
        write(*,*) TRIM(dum1), nbmin, nbmax, TRIM(dum2)
        STOP
    END IF
END DO

! unexpected end of file
10 write(*,*) 'reached end of mommat file, but unable to find a line with "KP:"'
write(*,*) 'Sometimes it happens when VASP file is read assuming WIEN2k', &
    ' format. If you read VASP file, the input file should have exact name', &
    ' WAVEDER.'
CLOSE(fid)
STOP

! read or data error, abort
20 write(*,*) 'unable to read the number of bands from the following line'
write(*,*) TRIM(cline)
write(*,*) 'the line is expected to look like this'
write(*,*) '   KP:     1 NEMIN NEMAX :     1  116 dE:  -5.0   5.0 K:'
write(*,*) 'If NEMIN > 1, it most likely implies that you have valence states', &
    ' below -5 Ry (check case.scf). If this is the case, set a lower value of', &
    ' Emin in optic.inop (e.g., -8 Ry)'
write(*,*) 'If the line printed above is not readable, most likely you read', &
    ' the VASP binary file assuming WIEN2k format. If you read VASP file,', &
    ' the input file should have the exact name WAVEDER.'
CLOSE(fid)
STOP

END SUBROUTINE read_mommat_nb
