SUBROUTINE surfstat
     ! surfstat calculates surface state using             !
     ! green's function method  ---  J.Phys.F.Met.Phys.15(1985)851-858      !
     !
     ! History:
     !
     !         by Quan Sheng Wu on 4/20/2010                                !
     !
     !            mpi version      4/21/2010
     !
     !            Quansheng Wu on Jan 30 2015 at ETH Zurich
     !
     ! Copyright (c) 2010 QuanSheng Wu. All rights reserved.

    USE wmpi
    USE para
    IMPLICIT NONE

    INTEGER :: ierr

    INTEGER :: doslfile, dosrfile, dosbulkfile

    ! general loop index
    INTEGER :: i, j, io

    ! kpoint loop index
    INTEGER :: ikp

    REAL(dp) :: emin
    REAL(dp) :: emax
    REAL(dp) :: k(2), w

    REAL(dp) :: time_start, time_end, time_togo

    REAL(dp), ALLOCATABLE :: omega(:)

    REAL(dp), ALLOCATABLE :: dos_l(:,:)
    REAL(dp), ALLOCATABLE :: dos_r(:,:)
    REAL(dp), ALLOCATABLE :: dos_l_only(:,:)
    REAL(dp), ALLOCATABLE :: dos_r_only(:,:)
    REAL(dp), ALLOCATABLE :: dos_l_mpi(:,:)
    REAL(dp), ALLOCATABLE :: dos_r_mpi(:,:)
    REAL(dp), ALLOCATABLE :: dos_bulk(:,:)
    REAL(dp), ALLOCATABLE :: dos_bulk_mpi(:,:)

    COMPLEX(dp), ALLOCATABLE :: GLL(:,:)
    COMPLEX(dp), ALLOCATABLE :: GRR(:,:)
    COMPLEX(dp), ALLOCATABLE :: GB (:,:)
    COMPLEX(dp), ALLOCATABLE :: H00(:,:)
    COMPLEX(dp), ALLOCATABLE :: H01(:,:)
    COMPLEX(dp), ALLOCATABLE :: ones(:,:)

    ! string for integer
    CHARACTER(LEN=140) :: ichar, jchar, kchar, fmt

    ALLOCATE( omega(omeganum))
    ALLOCATE( dos_l(knv2, omeganum))
    ALLOCATE( dos_r(knv2, omeganum))
    ALLOCATE( dos_l_only(knv2, omeganum))
    ALLOCATE( dos_r_only(knv2, omeganum))
    ALLOCATE( dos_l_mpi(knv2, omeganum))
    ALLOCATE( dos_r_mpi(knv2, omeganum))
    ALLOCATE( dos_bulk(knv2, omeganum))
    ALLOCATE( dos_bulk_mpi(knv2, omeganum))

    omega=0d0
    dos_l=0d0
    dos_r=0d0
    dos_l_only=0d0
    dos_r_only=0d0
    dos_l_mpi=0d0
    dos_r_mpi=0d0
    dos_bulk=0d0
    dos_bulk_mpi=0d0

    eta=(omegamax- omegamin)/dble(omeganum)*1.5d0

    DO i= 1, omeganum
        omega(i)=omegamin+(i-1)*(omegamax-omegamin)/dble(omeganum)
    ENDDO

    !> deal with phonon system
    !> for phonon system, omega should be changed to omega^2
    IF (index(Particle,'phonon')/=0) then
        DO i= 1, omeganum
            omega(i)= omega(i)*omega(i)
        ENDDO
    ENDIF

    ALLOCATE(GLL(Ndim, Ndim))
    ALLOCATE(GRR(Ndim, Ndim))
    ALLOCATE(GB (Ndim, Ndim))
    ALLOCATE(H00(Ndim, Ndim))
    ALLOCATE(H01(Ndim, Ndim))
    ALLOCATE(ones(Ndim, Ndim))
    GLL= 0d0
    GRR= 0d0
    GB = 0d0
    H00= 0d0
    H01= 0d0
    ones= 0d0

    DO i=1,Ndim
        ones(i,i)=1.0d0
    ENDDO

    time_start= 0d0
    time_end= 0d0
    DO ikp= 1+cpuid, knv2, num_cpu

        k= k2_path(ikp,:)

        call now(time_start)
        !> deal with phonon system
        !> get the hopping matrix between two principle layers
        if (index(Particle,'phonon')/=0.and.LOTO_correction) then
            call ham_qlayer2qlayer_LOTO(k,H00,H01)
        else
            call ham_qlayer2qlayer(k,H00,H01)
        endif

        DO j = 1, omeganum

            w=omega(j)

            !> calculate surface green function
            ! there are two method to calculate surface green's function
            ! the method in 1985 is better, you can find the ref in the
            ! subroutine
            CALL surfgreen_1985(w,GLL,GRR,GB,H00,H01,ones)
            ! CALL surfgreen_1984(w,GLL,GRR,H00,H01,ones)

            ! calculate spectral function
            DO i= 1, NtopOrbitals
                io= TopOrbitals(i)
                dos_l(ikp, j)=dos_l(ikp,j)- aimag(GLL(io,io))
            ENDDO ! i
            DO i= 1, NBottomOrbitals
                io= Ndim- Num_wann+ BottomOrbitals(i)
                dos_r(ikp, j)=dos_r(ikp,j)- aimag(GRR(io,io))
            ENDDO ! i
            DO i= 1, Ndim
                dos_bulk(ikp, j)=dos_bulk(ikp,j)- aimag(GB(i,i))
            ENDDO ! i

        ENDDO ! j
        CALL now(time_end)

        time_togo = (knv2-ikp)*(time_end-time_start)/num_cpu
        IF (ikp == 1) THEN
            WRITE(ichar, *) knv2
            WRITE(jchar, '("I", I0)') LEN(TRIM(ADJUSTL(ichar)))
            WRITE(ichar, *) int(time_togo)
            WRITE(kchar, '("F", I0, ".2")') LEN(TRIM(ADJUSTL(ichar)))+4
            fmt = "(' SSs ik/nk : ',"//TRIM(jchar)//",'/',"//TRIM(jchar)//",' , &
            & Estimated time remaining : ', "//TRIM(kchar)//", ' s')"
        ENDIF

        IF ( cpuid==0 .AND. MOD(ikp/num_cpu, 10)==0 .AND. ikp/=1 ) &
        WRITE( stdout, fmt ) ikp, knv2, time_togo

    ENDDO ! ikp

!> we don't have to do allreduce operation
#if defined (MPI)
    CALL mpi_reduce(dos_l, dos_l_mpi, size(dos_l), mpi_double_precision,&
         mpi_sum, 0, mpi_comm_world, ierr)
    CALL mpi_reduce(dos_r, dos_r_mpi, size(dos_r), mpi_double_precision,&
         mpi_sum, 0, mpi_comm_world, ierr)
    CALL mpi_reduce(dos_bulk, dos_bulk_mpi, size(dos_bulk), mpi_double_precision,&
         mpi_sum, 0, mpi_comm_world, ierr)
#else
    dos_l_mpi= dos_l
    dos_r_mpi= dos_r
    dos_bulk_mpi= dos_bulk
#endif

    dos_l=log(abs(dos_l_mpi))
    dos_r=log(abs(dos_r_mpi))
    dos_bulk=log(abs(dos_bulk_mpi)+eps9)
    DO ikp=1, knv2
        DO j=1, omeganum
            dos_l_only(ikp, j)= dos_l_mpi(ikp, j)- dos_bulk_mpi(ikp, j)
            IF (dos_l_only(ikp, j)<0) dos_l_only(ikp, j)=eps9
            dos_r_only(ikp, j)= dos_r_mpi(ikp, j)- dos_bulk_mpi(ikp, j)
            IF (dos_r_only(ikp, j)<0) dos_r_only(ikp, j)=eps9
        ENDDO
    ENDDO

    outfileindex= outfileindex+ 1
    doslfile= outfileindex
    outfileindex= outfileindex+ 1
    dosrfile= outfileindex
    outfileindex= outfileindex+ 1
    dosbulkfile= outfileindex

!> deal with phonon system
!> for phonon system, omega should be changed to omega^2
    IF (index(Particle,'phonon')/=0) then
        DO i= 1, omeganum
            omega(i)= sqrt(omega(i))
        ENDDO
    ENDIF

    IF (cpuid.eq.0)then
        OPEN (unit=doslfile, file='dos.dat_l')
        OPEN (unit=dosrfile, file='dos.dat_r')
        OPEN (unit=dosbulkfile, file='dos.dat_bulk')
        DO ikp=1, knv2
            DO j=1, omeganum
                WRITE(doslfile, '(30f16.8)')k2len(ikp), omega(j), dos_l(ikp, j), log(dos_l_only(ikp, j))
                WRITE(dosrfile, '(30f16.8)')k2len(ikp), omega(j), dos_r(ikp, j), log(dos_r_only(ikp, j))
                WRITE(dosbulkfile, '(30f16.8)')k2len(ikp), omega(j), dos_bulk(ikp, j)
            ENDDO
            WRITE(doslfile, *) ' '
            WRITE(dosrfile, *) ' '
            WRITE(dosbulkfile, *) ' '
        ENDDO
        CLOSE(doslfile)
        CLOSE(dosrfile)
        CLOSE(dosbulkfile)

        WRITE(stdout,*)'ndim',ndim
        WRITE(stdout,*) 'knv2,omeganum,eta',knv2, omeganum, eta
        WRITE(stdout,*)'calculate density of state successfully'
    ENDIF


    emin= minval(omega)
    emax= maxval(omega)
    !> write script for gnuplot
    outfileindex= outfileindex+ 1
    IF (cpuid==0) THEN
        OPEN(unit=outfileindex, file='surfdos_l.gnu')
        WRITE(outfileindex, '(a)')"set encoding iso_8859_1"
        WRITE(outfileindex, '(a)')'#set terminal  postscript enhanced color'
        WRITE(outfileindex, '(a)')"#set output 'surfdos_l.eps'"
        WRITE(outfileindex, '(3a)')'#set terminal  pngcairo truecolor enhanced', &
        '  font ", 60" size 1920, 1680'
        WRITE(outfileindex, '(3a)')'set terminal  png truecolor enhanced', &
        ' font ", 60" size 1920, 1680'
        WRITE(outfileindex, '(a)')"set output 'surfdos_l.png'"
        WRITE(outfileindex,'(2a)') 'set palette defined (-10 "#194eff", ', &
        '0 "white", 10 "red" )'
        WRITE(outfileindex, '(a)')'#set palette rgbformulae 33,13,10'
        WRITE(outfileindex, '(a)')'set style data linespoints'
        WRITE(outfileindex, '(a)')'set size 0.8, 1'
        WRITE(outfileindex, '(a)')'set origin 0.1, 0'
        WRITE(outfileindex, '(a)')'unset ztics'
        WRITE(outfileindex, '(a)')'unset key'
        WRITE(outfileindex, '(a)')'set pointsize 0.8'
        WRITE(outfileindex, '(a)')'set pm3d'
        WRITE(outfileindex, '(a)')'#set view equal xyz'
        WRITE(outfileindex, '(a)')'set view map'
        WRITE(outfileindex, '(a)')'set border lw 3'
        WRITE(outfileindex, '(a)')'#set cbtics font ",48"'
        WRITE(outfileindex, '(a)')'#set xtics font ",48"'
        WRITE(outfileindex, '(a)')'#set ytics font ",48"'
        WRITE(outfileindex, '(a)')'#set ylabel font ",48"'
        WRITE(outfileindex, '(a)')'set ylabel "Energy (eV)"'
        WRITE(outfileindex, '(a)')'#set xtics offset 0, -1'
        WRITE(outfileindex, '(a)')'#set ylabel offset -6, 0 '
        WRITE(outfileindex, '(a, f18.5, a)')'set xrange [0: ', maxval(k2len), ']'
        WRITE(outfileindex, '(a, f18.5, a, f18.5, a)')'set yrange [', emin, ':', emax, ']'
        WRITE(outfileindex, 202, advance="no") (k2line_name(i), k2line_stop(i), i=1, nk2lines)
        WRITE(outfileindex, 203)k2line_name(nk2lines+1), k2line_stop(nk2lines+1)
        DO i=1, nk2lines-1
            WRITE(outfileindex, 204)k2line_stop(i+1), emin, k2line_stop(i+1), emax
        ENDDO
        WRITE(outfileindex, '(a)')'set pm3d interpolate 2,2'
        WRITE(outfileindex, '(2a)')"splot 'dos.dat_l' u 1:2:3 w pm3d"
        CLOSE(outfileindex)
    ENDIF

    outfileindex= outfileindex+ 1
    IF (cpuid==0) THEN
        OPEN(unit=outfileindex, file='surfdos_l_only.gnu')
        WRITE(outfileindex, '(a)')"set encoding iso_8859_1"
        WRITE(outfileindex, '(a)')'#set terminal  postscript enhanced color'
        WRITE(outfileindex, '(a)')"#set output 'surfdos_l.eps'"
        WRITE(outfileindex, '(3a)')'#set terminal  pngcairo truecolor enhanced', &
        '  font ", 60" size 1920, 1680'
        WRITE(outfileindex, '(3a)')'set terminal  png truecolor enhanced', &
        ' font ", 60" size 1920, 1680'
        WRITE(outfileindex, '(a)')"set output 'surfdos_l_only.png'"
        WRITE(outfileindex,'(2a)') 'set palette defined (0  "white", ', &
        '6 "red", 20 "black" )'
        WRITE(outfileindex, '(a)')'#set palette rgbformulae 33,13,10'
        WRITE(outfileindex, '(a)')'set style data linespoints'
        WRITE(outfileindex, '(a)')'set size 0.8, 1'
        WRITE(outfileindex, '(a)')'set origin 0.1, 0'
        WRITE(outfileindex, '(a)')'unset ztics'
        WRITE(outfileindex, '(a)')'unset key'
        WRITE(outfileindex, '(a)')'set pointsize 0.8'
        WRITE(outfileindex, '(a)')'set pm3d'
        WRITE(outfileindex, '(a)')'#set view equal xyz'
        WRITE(outfileindex, '(a)')'set view map'
        WRITE(outfileindex, '(a)')'set border lw 3'
        WRITE(outfileindex, '(a)')'#set cbtics font ",48"'
        WRITE(outfileindex, '(a)')'#set xtics font ",48"'
        WRITE(outfileindex, '(a)')'#set ytics font ",48"'
        WRITE(outfileindex, '(a)')'unset cbtics'
        WRITE(outfileindex, '(a)')'#set ylabel font ",48"'
        WRITE(outfileindex, '(a)')'set ylabel "Energy (eV)"'
        WRITE(outfileindex, '(a)')'#set xtics offset 0, -1'
        WRITE(outfileindex, '(a)')'#set ylabel offset -6, 0 '
        WRITE(outfileindex, '(a, f18.5, a)')'set xrange [0: ', maxval(k2len), ']'
        WRITE(outfileindex, '(a, f18.5, a, f18.5, a)')'set yrange [', emin, ':', emax, ']'
        WRITE(outfileindex, 202, advance="no") (k2line_name(i), k2line_stop(i), i=1, nk2lines)
        WRITE(outfileindex, 203)k2line_name(nk2lines+1), k2line_stop(nk2lines+1)

        DO i=1, nk2lines-1
            WRITE(outfileindex, 204)k2line_stop(i+1), emin, k2line_stop(i+1), emax
        ENDDO
        WRITE(outfileindex, '(a)')'set pm3d interpolate 2,2'
        WRITE(outfileindex, '(2a)')"splot 'dos.dat_l' u 1:2:(exp($4)) w pm3d"
        CLOSE(outfileindex)
    ENDIF

    !> WRITE script for gnuplot
    outfileindex= outfileindex+ 1
    IF (cpuid==0) THEN
        OPEN(unit=outfileindex, file='surfdos_r_only.gnu')
        WRITE(outfileindex, '(a)')"set encoding iso_8859_1"
        WRITE(outfileindex, '(a)')'#set terminal  postscript enhanced color'
        WRITE(outfileindex, '(a)')"#set output 'surfdos_r.eps'"
        WRITE(outfileindex, '(3a)')'#set terminal pngcairo truecolor enhanced', &
        '  font ", 60" size 1920, 1680'
        WRITE(outfileindex, '(3a)')'set terminal png truecolor enhanced', &
        ' font ", 60" size 1920, 1680'
        WRITE(outfileindex, '(a)')"set output 'surfdos_r_only.png'"
        WRITE(outfileindex,'(2a)') 'set palette defined (0  "white", ', &
        '6 "red", 20 "black" )'
        WRITE(outfileindex, '(a)')'#set palette rgbformulae 33,13,10'
        WRITE(outfileindex, '(a)')'set style data linespoints'
        WRITE(outfileindex, '(a)')'unset ztics'
        WRITE(outfileindex, '(a)')'unset key'
        WRITE(outfileindex, '(a)')'set pointsize 0.8'
        WRITE(outfileindex, '(a)')'set pm3d'
        WRITE(outfileindex, '(a)')'set border lw 3'
        WRITE(outfileindex, '(a)')'set size 0.8, 1'
        WRITE(outfileindex, '(a)')'set origin 0.1, 0'
        WRITE(outfileindex, '(a)')'#set size ratio -1'
        WRITE(outfileindex, '(a)')'#set view equal xyz'
        WRITE(outfileindex, '(a)')'set view map'
        WRITE(outfileindex, '(a)')'#set cbtics font ",48"'
        WRITE(outfileindex, '(a)')'#set xtics font ",48"'
        WRITE(outfileindex, '(a)')'#set ytics font ",48"'
        WRITE(outfileindex, '(a)')'#set ylabel font ",48"'
        WRITE(outfileindex, '(a)')'set ylabel "Energy (eV)"'
        WRITE(outfileindex, '(a)')'unset cbtics'
        WRITE(outfileindex, '(a)')'#set xtics offset 0, -1'
        WRITE(outfileindex, '(a)')'#set ylabel offset -6, 0 '
        WRITE(outfileindex, '(a, f18.5, a)')'set xrange [0: ', maxval(k2len), ']'
        WRITE(outfileindex, '(a, f18.5, a, f18.5, a)')'set yrange [', emin, ':', emax, ']'
        WRITE(outfileindex, 202, advance="no") (k2line_name(i), k2line_stop(i), i=1, nk2lines)
        WRITE(outfileindex, 203)k2line_name(nk2lines+1), k2line_stop(nk2lines+1)

        DO i=1, nk2lines-1
            WRITE(outfileindex, 204)k2line_stop(i+1), emin, k2line_stop(i+1), emax
        ENDDO
        WRITE(outfileindex, '(a)')'set pm3d interpolate 2,2'
        WRITE(outfileindex, '(2a)')"splot 'dos.dat_r' u 1:2:(exp($4)) w pm3d"
        CLOSE(outfileindex)
    ENDIF


    !> WRITE script for gnuplot
    outfileindex= outfileindex+ 1
    IF (cpuid==0) THEN
        OPEN(unit=outfileindex, file='surfdos_r.gnu')
        WRITE(outfileindex, '(a)')"set encoding iso_8859_1"
        WRITE(outfileindex, '(a)')'#set terminal  postscript enhanced color'
        WRITE(outfileindex, '(a)')"#set output 'surfdos_r.eps'"
        WRITE(outfileindex, '(3a)')'#set terminal pngcairo truecolor enhanced', &
        '  font ", 60" size 1920, 1680'
        WRITE(outfileindex, '(3a)')'set terminal png truecolor enhanced', &
        ' font ", 60" size 1920, 1680'
        WRITE(outfileindex, '(a)')"set output 'surfdos_r.png'"
        WRITE(outfileindex,'(2a)') 'set palette defined (-10 "#194eff", ', &
        '0 "white", 10 "red" )'
        WRITE(outfileindex, '(a)')'#set palette rgbformulae 33,13,10'
        WRITE(outfileindex, '(a)')'set style data linespoints'
        WRITE(outfileindex, '(a)')'unset ztics'
        WRITE(outfileindex, '(a)')'unset key'
        WRITE(outfileindex, '(a)')'set pointsize 0.8'
        WRITE(outfileindex, '(a)')'set pm3d'
        WRITE(outfileindex, '(a)')'set border lw 3'
        WRITE(outfileindex, '(a)')'set size 0.8, 1'
        WRITE(outfileindex, '(a)')'set origin 0.1, 0'
        WRITE(outfileindex, '(a)')'#set size ratio -1'
        WRITE(outfileindex, '(a)')'#set view equal xyz'
        WRITE(outfileindex, '(a)')'set view map'
        WRITE(outfileindex, '(a)')'#set cbtics font ",48"'
        WRITE(outfileindex, '(a)')'#set xtics font ",48"'
        WRITE(outfileindex, '(a)')'#set ytics font ",48"'
        WRITE(outfileindex, '(a)')'#set ylabel font ",48"'
        WRITE(outfileindex, '(a)')'set ylabel "Energy (eV)"'
        WRITE(outfileindex, '(a)')'#set xtics offset 0, -1'
        WRITE(outfileindex, '(a)')'#set ylabel offset -6, 0 '
        WRITE(outfileindex, '(a, f18.5, a)')'set xrange [0: ', maxval(k2len), ']'
        WRITE(outfileindex, '(a, f18.5, a, f18.5, a)')'set yrange [', emin, ':', emax, ']'
        WRITE(outfileindex, 202, advance="no") (k2line_name(i), k2line_stop(i), i=1, nk2lines)
        WRITE(outfileindex, 203)k2line_name(nk2lines+1), k2line_stop(nk2lines+1)
        DO i=1, nk2lines-1
            WRITE(outfileindex, 204)k2line_stop(i+1), emin, k2line_stop(i+1), emax
        ENDDO
        WRITE(outfileindex, '(a)')'set pm3d interpolate 2,2'
        WRITE(outfileindex, '(2a)')"splot 'dos.dat_r' u 1:2:3 w pm3d"
        CLOSE(outfileindex)
    ENDIF

    !> WRITE script for gnuplot
    outfileindex= outfileindex+ 1
    IF (cpuid==0) THEN
        OPEN(unit=outfileindex, file='surfdos_bulk.gnu')
        WRITE(outfileindex, '(a)')"set encoding iso_8859_1"
        WRITE(outfileindex, '(a)')'#set terminal  postscript enhanced color'
        WRITE(outfileindex, '(a)')"#set output 'surfdos_bulk.eps'"
        WRITE(outfileindex, '(3a)')'#set terminal  pngcairo truecolor enhanced', &
        '  font ", 60" size 1920, 1680'
        WRITE(outfileindex, '(3a)')'set terminal  png truecolor enhanced', &
        ' font ", 60" size 1920, 1680'
        WRITE(outfileindex, '(a)')"set output 'surfdos_bulk.png'"
        WRITE(outfileindex,'(2a)') 'set palette defined (-10 "#194eff", ', &
        '0 "white", 10 "red" )'
        WRITE(outfileindex, '(a)')'#set palette rgbformulae 33,13,10'
        WRITE(outfileindex, '(a)')'set style data linespoints'
        WRITE(outfileindex, '(a)')'set size 0.8, 1'
        WRITE(outfileindex, '(a)')'set origin 0.1, 0'
        WRITE(outfileindex, '(a)')'unset ztics'
        WRITE(outfileindex, '(a)')'unset key'
        WRITE(outfileindex, '(a)')'set pointsize 0.8'
        WRITE(outfileindex, '(a)')'set pm3d'
        WRITE(outfileindex, '(a)')'#set view equal xyz'
        WRITE(outfileindex, '(a)')'set view map'
        WRITE(outfileindex, '(a)')'set border lw 3'
        WRITE(outfileindex, '(a)')'#set cbtics font ",48"'
        WRITE(outfileindex, '(a)')'#set xtics font ",48"'
        WRITE(outfileindex, '(a)')'#set ytics font ",48"'
        WRITE(outfileindex, '(a)')'#set ylabel font ",48"'
        WRITE(outfileindex, '(a)')'unset cbtics'
        WRITE(outfileindex, '(a)')'set ylabel "Energy (eV)"'
        WRITE(outfileindex, '(a)')'#set xtics offset 0, -1'
        WRITE(outfileindex, '(a)')'#set ylabel offset -6, 0 '
        WRITE(outfileindex, '(a, f18.5, a)')'set xrange [0: ', maxval(k2len), ']'
        WRITE(outfileindex, '(a, f18.5, a, f18.5, a)')'set yrange [', emin, ':', emax, ']'
        WRITE(outfileindex, 202, advance="no") (k2line_name(i), k2line_stop(i), i=1, nk2lines)
        WRITE(outfileindex, 203)k2line_name(nk2lines+1), k2line_stop(nk2lines+1)

        DO i=1, nk2lines-1
            WRITE(outfileindex, 204)k2line_stop(i+1), emin, k2line_stop(i+1), emax
        ENDDO
        WRITE(outfileindex, '(a)')'set pm3d interpolate 2,2'
        WRITE(outfileindex, '(2a)')"splot 'dos.dat_bulk' u 1:2:(exp($3)) w pm3d"
        CLOSE(outfileindex)
    ENDIF

202 FORMAT('set xtics (',:20('"',A1,'" ',F8.5,','))
203 FORMAT(A1,'" ',F8.5,')')
204 FORMAT('set arrow from ',F8.5,',',F10.5,' to ',F8.5,',',F10.5, ' nohead front lw 3')

#if defined (MPI)
    CALL mpi_barrier(mpi_cmw, ierr)
#endif

    DEALLOCATE( omega)
    DEALLOCATE( dos_l)
    DEALLOCATE( dos_r)
    DEALLOCATE( dos_l_only)
    DEALLOCATE( dos_r_only)
    DEALLOCATE( dos_l_mpi)
    DEALLOCATE( dos_r_mpi)
    DEALLOCATE( dos_bulk)
    DEALLOCATE( dos_bulk_mpi)
    DEALLOCATE(GLL)
    DEALLOCATE(GRR)
    DEALLOCATE(GB )
    DEALLOCATE(H00)
    DEALLOCATE(H01)
    DEALLOCATE(ones)

    RETURN

END SUBROUTINE surfstat
