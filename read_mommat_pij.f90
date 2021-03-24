SUBROUTINE read_mommat_pij (fid, nb, nbb, & ! <- args in
            iline, & ! <-> args in-out
            pij, dEij) ! -> args out 

! Read momentum matrix elements <i|p_a|j> (a=x,y,z) and energy differences 
! E_i - E_j in a k-point block of mommat file

!! Variables in-out

implicit none
INTEGER, intent(in) :: &
    fid, & ! file ID
    nb, & ! number of bands
    nbb ! number of band-to-band transition
INTEGER, intent(inout) :: &
    iline ! current line number during reading of mommat file
REAL(kind=4), intent(out) :: &
    dEij(nb,nb) ! energy differences E_i-E_j [Ha]
COMPLEX(kind=4), intent(out) :: &
    pij(3,nb,nb) ! momentum matrix elements [at. units]

!! Variables internal

CHARACTER(len=256) :: cline ! text in the current line
INTEGER :: row  ! counter for rows
INTEGER :: bii, bjj  ! band indices used to determine the total number of bands
REAL(kind=4) :: p1_Re, p1_Im, p2_Re, p2_Im, p3_Re, p3_Im

!! Read momentum matrix element

row = 0
DO WHILE ( .true. )
    READ(fid,'(A)',END=10) cline ! read current line as text
    iline = iline + 1 ! count total lines read since opening the file
    IF (LEN(TRIM(cline)) == 0) CYCLE ! empty line, skip
    ! check for "KP:" in the line; it should NOT be there
    IF (INDEX(cline,'KP:') .ne. 0) GOTO 20 ! corrupted record (too short)
    row = row + 1 ! useful line
    ! WIEN2k prior May 2020 (the case.mommat2 is formated, there will be 
    ! an error for nb > 9999!)
    ! uncomment next 2 lines and coment the the following READ statement
    ! if you use an older version
    !READ(cline,'(3X,2I4,6E13.6,F13.8)',ERR=20) bii ,bjj, & !...
    !    p1_Re, p1_Im, p2_Re, p2_Im, p3_Re, p3_Im, dEij(bii,bjj)
    
    ! WIEN2k after May 2020 (the case.mommat2 file is not a fixed format)
    READ(cline,*,ERR=30) bii ,bjj, & !...
        p1_Re, p1_Im, p2_Re, p2_Im, p3_Re, p3_Im, dEij(bii,bjj)
    pij(1,bii,bjj) = CMPLX(p1_Re,p1_Im) ! convert to complex
    pij(2,bii,bjj) = CMPLX(p2_Re,p2_Im)
    pij(3,bii,bjj) = CMPLX(p3_Re,p3_Im)
    ! ensure Hermitian properties of pij
    pij(1,bjj,bii) = CONJG( pij(1,bii,bjj) ) ! conjugate
    pij(2,bjj,bii) = CONJG( pij(2,bii,bjj) )
    pij(3,bjj,bii) = CONJG( pij(3,bii,bjj) )
    ! convert energies [Ry] -> [Ha]
    dEij(bii,bjj) = dEij(bii,bjj)/2
    dEij(bjj,bii) = -dEij(bii,bjj)
    ! checks
    IF (row==1) THEN ! beginning of the record bii=bjj=1
        IF (bii/=1 .or. bjj/=1) THEN
            write(*,*) 'error reading mommat file'
            write(*,*) 'expected bii=bjj=1, but got', bii, ' and', bjj
            write(*,'(A,I0)') 'line number ', iline
            STOP
        END IF
    ELSEIF (row==nbb) THEN ! end of the record bii=bjj=nb
        IF (bii/=bjj .or. bii/=nb .or. bjj /=nb) THEN
            write(*,*) 'error reading mommat file'
            write(*,*) 'expected bii=bjj, but got', bii, ' /=', bjj
            write(*,'(A,I0)') 'line number ', iline
            STOP
        END IF
        RETURN ! done
    END IF
END DO

! illegal (unexpected) end of mommat file
10 write(*,*) 'ERROR: reached end of mommat file too early'
write(*,'(A,I0,A,I0)') 'expected bii=bjj, but got ', bii, ' /= ', bjj
write(*,'(A,I0)') 'line number ', iline
CLOSE(fid)
STOP

! error reading the line (old format)
20 write(*,*) 'ERROR reading the line from mommat file containing'
write(*,'(A)') TRIM(cline)
write(*,'(A,I0)') 'line number ', iline
write(*,*) 'while expecting a _fixed_ format line like this'
write(*,'(A)')  '      1  32 0.270604E-02 0.659085E-03 0.555417E-03'//&
    '-0.308709E-02-0.132694E-11-0.600436E-11   2.35217055'
CLOSE(fid)
STOP

! error reading the line (new format)
30 write(*,*) 'ERROR reading the line from mommat file containing'
write(*,'(A)') TRIM(cline)
write(*,'(A,I0)') 'line number ', iline
write(*,*) 'while expecting a _free_ format line like this'
write(*,'(A)')  '  1  32 0.270604E-02 0.659085E-03 0.555417E-03'//&
    ' -0.308709E-02 -0.132694E-11 -0.600436E-11 2.35217055'
CLOSE(fid)
STOP

END SUBROUTINE read_mommat_pij
