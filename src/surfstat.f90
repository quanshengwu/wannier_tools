  subroutine surfstat
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

     use wmpi
     use para
     implicit none

     integer :: ierr, doslfile, dosrfile, dosbulkfile

     ! general loop index
     integer :: i, j, io, ikp, nw_half, spindoslfile, spindosrfile

     real(dp) :: emin, emax, w, k(2), time_start, time_end, s0(3), s1(3), eta_broadening

     real(dp), allocatable :: omega(:)

     real(dp), allocatable :: dos_l(:,:), dos_r(:,:), dos_l_only(:, :), dos_r_only(:,:)
     real(dp), allocatable :: dos_l_mpi(:,:), dos_r_mpi(:,:), dos_bulk(:,:), dos_bulk_mpi(:,:)

     complex(dp), allocatable :: GLL(:,:), GRR(:,:), GB (:,:), H00(:,:), H01(:,:), ones(:,:)
 
     ! Spin resolved component
     REAL(DP),    ALLOCATABLE  :: sx_l(:, :), sy_l(:, :), sz_l(:, :)
     REAL(DP),    ALLOCATABLE  :: sx_r(:, :), sy_r(:, :), sz_r(:, :)
     REAL(DP),    ALLOCATABLE  :: sx_l_mpi(:, :), sy_l_mpi(:, :), sz_l_mpi(:, :)
     REAL(DP),    ALLOCATABLE  :: sx_r_mpi(:, :), sy_r_mpi(:, :), sz_r_mpi(:, :)
     COMPLEX(DP), ALLOCATABLE  :: sigma_x(:,:), sigma_y(:,:), sigma_z(:,:)
     COMPLEX(DP), ALLOCATABLE  :: ctemp(:,:)


     !> for special line
    !ke(1,:)=(/0.0d0, 0.00d0/)
    !kp(1, 2)= surf_onsite
    !ke(1, 2)= surf_onsite
    !kp(2,:)=(/0.0d0, 0.00d0/)  ; k2line_name(2)= 'G'
    !k1= -0.5d0*ka2+0.5d0*kb2
    !t1= 0d0
    !do i=1, NN
    !   s= (i-1d0)/(NN-1d0)-0.5d0
    !   k2= k1
    !   kpoint(i, 1)= s
    !   kpoint(i, 2)= -2.818d0*s**3-0.2955*s
    !   k1= kpoint(i, 1)*ka2+ kpoint(i, 2)*kb2
    !   temp= dsqrt((k2(1)- k1(1))**2 &
    !         +(k2(2)- k1(2))**2)

    !   if (i.gt.1) THEN
    !      t1=t1+temp
    !   endif
    !   k2len(i)= t1
    !   print *, i, t1
    !ENDDO

     allocate( omega(omeganum), dos_l(knv2, omeganum), dos_r(knv2, omeganum))
     allocate( dos_l_only(knv2, omeganum), dos_r_only(knv2, omeganum))
     allocate( dos_l_mpi(knv2, omeganum), dos_r_mpi(knv2, omeganum))
     allocate( dos_bulk(knv2, omeganum), dos_bulk_mpi(knv2, omeganum))
     omega=0d0;  dos_l=0d0;  dos_r=0d0;  dos_l_only=0d0;  dos_r_only=0d0
     dos_l_mpi=0d0;  dos_r_mpi=0d0;  dos_bulk=0d0;  dos_bulk_mpi=0d0

     ALLOCATE( sx_l(knv2, omeganum), sy_l(knv2, omeganum), sz_l(knv2, omeganum))
     ALLOCATE( sx_l_mpi(knv2, omeganum), sy_l_mpi(knv2, omeganum), sz_l_mpi(knv2, omeganum))
     ALLOCATE( sx_r(knv2, omeganum), sy_r(knv2, omeganum), sz_r(knv2, omeganum))
     ALLOCATE( sx_r_mpi(knv2, omeganum), sy_r_mpi(knv2, omeganum), sz_r_mpi(knv2, omeganum))
     ALLOCATE( sigma_x(ndim,ndim), sigma_y(ndim,ndim), sigma_z(ndim,ndim), ctemp(ndim,ndim))
     sigma_x      = 0d0;      sigma_y      = 0d0;      sigma_z      = 0d0
     sx_l         = 0d0;      sy_l         = 0d0;      sz_l         = 0d0
     sx_r         = 0d0;      sy_r         = 0d0;      sz_r         = 0d0
     sx_l_mpi     = 0d0;      sy_l_mpi     = 0d0;      sz_l_mpi     = 0d0
     sx_r_mpi     = 0d0;      sy_r_mpi     = 0d0;      sz_r_mpi     = 0d0


     eta_broadening=(omegamax- omegamin)/dble(omeganum)*3.0d0

     do i= 1, omeganum
        omega(i)=omegamin+(i-1)*(omegamax-omegamin)/dble(omeganum)
     enddo

     !> deal with phonon system
     !> for phonon system, omega should be changed to omega^2
     if (index(Particle,'phonon')/=0) then
        do i= 1, omeganum
           omega(i)= omega(i)*omega(i)
        enddo
     endif

     allocate(GLL(Ndim, Ndim), GRR(Ndim, Ndim), GB (Ndim, Ndim))
     allocate(H00(Ndim, Ndim), H01(Ndim, Ndim), ones(Ndim, Ndim))
     GLL= 0d0; GRR= 0d0; GB = 0d0; H00= 0d0; H01= 0d0; ones= 0d0

     !> spin operator matrix
     !> Note, the basis here should be |↑↑↓↓>
     nw_half = Num_wann/2
     do i=1, Np
        do j=1, nw_half
           sigma_x( Num_wann*(i-1)+j        , Num_wann*(i-1)+j+nw_half ) =  1.0d0
           sigma_x( Num_wann*(i-1)+j+nw_half, Num_wann*(i-1)+j         ) =  1.0d0
           sigma_y( Num_wann*(i-1)+j        , Num_wann*(i-1)+j+nw_half ) = -zi
           sigma_y( Num_wann*(i-1)+j+nw_half, Num_wann*(i-1)+j         ) =  zi
           sigma_z( Num_wann*(i-1)+j        , Num_wann*(i-1)+j         ) =  1.0d0
           sigma_z( Num_wann*(i-1)+j+nw_half, Num_wann*(i-1)+j+nw_half ) = -1.0d0
        enddo 
     enddo


     do i=1,Ndim
        ones(i,i)=1.0d0
     enddo

     time_start= 0d0
     time_end= 0d0
     do ikp= 1+cpuid, knv2, num_cpu
        if (cpuid==0.and. mod(ikp/num_cpu, 10)==0) &
           WRITE(stdout, *) 'SurfaceSS, ik', ikp, 'Nk', knv2, 'time left', &
           (knv2-ikp)*(time_end- time_start)/num_cpu, ' s'

        k= k2_path(ikp,:)

        call now(time_start)
        !> deal with phonon system
        !> get the hopping matrix between two principle layers
        if (index(Particle,'phonon')/=0.and.LOTO_correction) then
           call ham_qlayer2qlayer_LOTO(k,H00,H01)
        else
           call ham_qlayer2qlayer(k,H00,H01)
        endif

        do j = 1, omeganum
           w=omega(j)

           !> calculate surface green function
           ! there are two method to calculate surface green's function
           ! the method in 1985 is better, you can find the ref in the
           ! subroutine
            call surfgreen_1985(w,GLL,GRR,GB,H00,H01,ones, eta_broadening)
           ! call surfgreen_1984(w,GLL,GRR,H00,H01,ones, eta_broadening)

           ! calculate spectral function
           do i= 1, NtopOrbitals
              io= TopOrbitals(i)
              dos_l(ikp, j)=dos_l(ikp,j)- aimag(GLL(io,io))
           enddo ! i
           do i= 1, NBottomOrbitals
              io= Ndim- Num_wann+ BottomOrbitals(i)
              dos_r(ikp, j)=dos_r(ikp,j)- aimag(GRR(io,io))
           enddo ! i
           do i= 1, Ndim
              dos_bulk(ikp, j)=dos_bulk(ikp,j)- aimag(GB(i,i))
           enddo ! i

            ! Spin resolved sprectrum
            call mat_mul(ndim,gll,sigma_x,ctemp)
            do i = 1, NtopOrbitals
                io = TopOrbitals(i)
                sx_l_mpi(ikp, j) = sx_l_mpi(ikp, j)- AIMAG(ctemp(io,io))
            enddo ! i
            call mat_mul(ndim,gll,sigma_y,ctemp)
            do i = 1, NtopOrbitals
                io = TopOrbitals(i)
                sy_l_mpi(ikp, j) = sy_l_mpi(ikp, j)- AIMAG(ctemp(io,io))
            enddo !
            call mat_mul(ndim,gll,sigma_z,ctemp)
            do i = 1, NtopOrbitals
                io = TopOrbitals(i)
                sz_l_mpi(ikp, j) = sz_l_mpi(ikp, j)- AIMAG(ctemp(io,io))
            enddo ! i
            call mat_mul(ndim,grr,sigma_x,ctemp)
            do i = 1, NtopOrbitals
                io = TopOrbitals(i)
                sx_r_mpi(ikp, j) = sx_r_mpi(ikp, j)- AIMAG(ctemp(io,io))
            enddo ! i
            call mat_mul(ndim,grr,sigma_y,ctemp)
            do i = 1, NtopOrbitals
                io = TopOrbitals(i)
                sy_r_mpi(ikp, j) = sy_r_mpi(ikp, j)- AIMAG(ctemp(io,io))
            enddo !
            call mat_mul(ndim,grr,sigma_z,ctemp)
            do i = 1, NtopOrbitals
                io = TopOrbitals(i)
                sz_r_mpi(ikp, j) = sz_r_mpi(ikp, j)- AIMAG(ctemp(io,io))
            enddo ! i


        enddo ! j
        call now(time_end)
     enddo ! ikp

!> we do have to do allreduce operation
#ifdef MPI
    call mpi_allreduce(sx_l_mpi, sx_l, SIZE(sx_l), mpi_double_precision,&
    mpi_sum, mpi_comm_world, ierr)
    call mpi_allreduce(sy_l_mpi, sy_l, SIZE(sy_l), mpi_double_precision,&
    mpi_sum, mpi_comm_world, ierr)
    call mpi_allreduce(sz_l_mpi, sz_l, SIZE(sz_l), mpi_double_precision,&
    mpi_sum, mpi_comm_world, ierr)
    call mpi_allreduce(sx_r_mpi, sx_r, SIZE(sx_r), mpi_double_precision,&
    mpi_sum, mpi_comm_world, ierr)
    call mpi_allreduce(sy_r_mpi, sy_r, SIZE(sy_r), mpi_double_precision,&
    mpi_sum, mpi_comm_world, ierr)
    call mpi_allreduce(sz_r_mpi, sz_r, SIZE(sz_r), mpi_double_precision,&
    mpi_sum, mpi_comm_world, ierr)
#else
    sx_l     = sx_l_mpi;      sy_l     = sy_l_mpi;      sz_l     = sz_l_mpi
    sx_r     = sx_r_mpi;      sy_r     = sy_r_mpi;      sz_r     = sz_r_mpi
#endif

    deallocate( sx_l_mpi, sy_l_mpi, sz_l_mpi )
    deallocate( sx_r_mpi, sy_r_mpi, sz_r_mpi )
    deallocate( sigma_x, sigma_y, sigma_z, ctemp )


     !> we don't have to do allreduce operation
#if defined (MPI)
     call mpi_reduce(dos_l, dos_l_mpi, size(dos_l), mpi_double_precision,&
                     mpi_sum, 0, mpi_comm_world, ierr)
     call mpi_reduce(dos_r, dos_r_mpi, size(dos_r), mpi_double_precision,&
                     mpi_sum, 0, mpi_comm_world, ierr)
     call mpi_reduce(dos_bulk, dos_bulk_mpi, size(dos_bulk), mpi_double_precision,&
                     mpi_sum, 0, mpi_comm_world, ierr)
#else
     dos_l_mpi= dos_l
     dos_r_mpi= dos_r
     dos_bulk_mpi= dos_bulk
#endif

     dos_l=log(abs(dos_l_mpi))
     dos_r=log(abs(dos_r_mpi))
     dos_bulk=log(abs(dos_bulk_mpi)+eps9)
     do ikp=1, knv2
        do j=1, omeganum
           dos_l_only(ikp, j)= dos_l_mpi(ikp, j)- dos_bulk_mpi(ikp, j)
           if (dos_l_only(ikp, j)<0) dos_l_only(ikp, j)=eps9
           dos_r_only(ikp, j)= dos_r_mpi(ikp, j)- dos_bulk_mpi(ikp, j)
           if (dos_r_only(ikp, j)<0) dos_r_only(ikp, j)=eps9
        enddo
     enddo

     outfileindex= outfileindex+ 1
     doslfile= outfileindex
     outfileindex= outfileindex+ 1
     dosrfile= outfileindex
     outfileindex= outfileindex+ 1
     dosbulkfile= outfileindex
     outfileindex = outfileindex+ 1
     spindoslfile = outfileindex
     outfileindex = outfileindex+ 1
     spindosrfile = outfileindex

     !> deal with phonon system
     !> for phonon system, omega should be changed to omega^2
     if (index(Particle,'phonon')/=0) then
        do i= 1, omeganum
           omega(i)= sqrt(omega(i))
        enddo
     endif

 
     ! Write surface state to files
     IF (cpuid .eq. 0) THEN
         open(unit=doslfile    , file='dos.dat_l')
         open(unit=dosrfile    , file='dos.dat_r')
         open(unit=dosbulkfile , file='dos.dat_bulk')
         open(unit=spindoslfile, file='spindos.dat_l')
         open(unit=spindosrfile, file='spindos.dat_r')
         write(doslfile, '("#", a12, 6a17)') ' k(1/Ang)', ' E(eV)', 'dos_l', 'dos_l_only'
         write(dosrfile, '("#", a12, 6a17)') ' k(1/Ang)', ' E(eV)', 'dos_r', 'dos_r_only'
         write(dosbulkfile, '("#", a12, 6a17)') ' k(1/Ang)', ' E(eV)', 'dos_bulk'
         write(spindoslfile, '("#", a)') ' spin dos_l, the axis is rotated as '
         write(spindoslfile, '("#", a)') " x is along R1', z is along R1'xR2', y is along z x y"
         write(spindoslfile, '("#", a12, 6a17)') ' k(1/Ang)', ' E(eV)', 'sx', 'sy', 'sz'
         write(spindosrfile, '("#", a)') ' spin dos_r, the axis is rotated as '
         write(spindosrfile, '("#", a)') " x is along R1', z is along R1'xR2', y is along z x y"
         write(spindosrfile, '("#", a12, 6a17)') ' k(1/Ang)', ' E(eV)', 'sx', 'sy', 'sz'
         do ikp = 1, knv2
            do j = 1, omeganum
                write(doslfile,    2002) k2len(ikp)*Angstrom2atomic, omega(j)/eV2Hartree, dos_l(ikp, j), dos_l_only(ikp, j)
                write(dosrfile,    2002) k2len(ikp)*Angstrom2atomic, omega(j)/eV2Hartree, dos_r(ikp, j), dos_r_only(ikp, j)
                write(dosbulkfile, 2003) k2len(ikp)*Angstrom2atomic, omega(j)/eV2Hartree, dos_bulk(ikp, j)
                s0(1)=sx_l(ikp, j); s0(2)=sy_l(ikp, j); s0(3)=sz_l(ikp, j); 
                call rotate(s0, s1)
                write(spindoslfile,2001) k2len(ikp)*Angstrom2atomic, omega(j)/eV2Hartree, s1
                s0(1)=sx_r(ikp, j); s0(2)=sy_r(ikp, j); s0(3)=sz_r(ikp, j); 
                call rotate(s0, s1)
                write(spindosrfile,2001) k2len(ikp)*Angstrom2atomic, omega(j)/eV2Hartree, s1
             ENDDO
             write(doslfile    , *)
             write(dosrfile    , *)
             write(dosbulkfile , *)
             write(spindoslfile, *)
             write(spindosrfile, *)
         ENDDO
         close(doslfile)
         close(dosrfile)
         close(dosbulkfile)
         close(spindoslfile)
         close(spindosrfile)
         write(stdout,*)'ndim',ndim
         write(stdout,*) 'knv2,omeganum,eta_broadening',knv2, omeganum, eta_broadening/eV2Hartree
         write(stdout,*)'calculate density of state successfully'
     ENDIF



2001 FORMAT(5(1X,E16.8))
2002 FORMAT(4(1X,E16.8))
2003 FORMAT(3(1X,E16.8))

     emin= minval(omega)/eV2Hartree
     emax= maxval(omega)/eV2Hartree
     !> write script for gnuplot
     outfileindex= outfileindex+ 1
     if (cpuid==0) then
        open(unit=outfileindex, file='surfdos_l.gnu')
        write(outfileindex, '(a)')"set encoding iso_8859_1"
        write(outfileindex, '(a)')'#set terminal  postscript enhanced color'
        write(outfileindex, '(a)')"#set output 'surfdos_l.eps'"
        write(outfileindex, '(3a)')'#set terminal  pngcairo truecolor enhanced', &
           '  font ", 60" size 1920, 1680'
        write(outfileindex, '(3a)')'set terminal  png truecolor enhanced', &
           ' font ", 60" size 1920, 1680'
        write(outfileindex, '(a)')"set output 'surfdos_l.png'"
        write(outfileindex,'(2a)') 'set palette defined (-10 "#194eff", ', &
           '0 "white", 10 "red" )'
        write(outfileindex, '(a)')'#set palette rgbformulae 33,13,10'
        write(outfileindex, '(a)')'set style data linespoints'
        write(outfileindex, '(a)')'set size 0.8, 1'
        write(outfileindex, '(a)')'set origin 0.1, 0'
        write(outfileindex, '(a)')'unset ztics'
        write(outfileindex, '(a)')'unset key'
        write(outfileindex, '(a)')'set pointsize 0.8'
        write(outfileindex, '(a)')'set pm3d'
        write(outfileindex, '(a)')'#set view equal xyz'
        write(outfileindex, '(a)')'set view map'
        write(outfileindex, '(a)')'set border lw 3'
        write(outfileindex, '(a)')'#set cbtics font ",48"'
        write(outfileindex, '(a)')'#set xtics font ",48"'
        write(outfileindex, '(a)')'#set ytics font ",48"'
        write(outfileindex, '(a)')'#set ylabel font ",48"'
        write(outfileindex, '(a)')'set ylabel "Energy (eV)"'
        write(outfileindex, '(a)')'#set xtics offset 0, -1'
        write(outfileindex, '(a)')'#set ylabel offset -6, 0 '
        write(outfileindex, '(a, f18.5, a)')'set xrange [0: ', maxval(k2len)*Angstrom2atomic, ']'
        write(outfileindex, '(a, f18.5, a, f18.5, a)')'set yrange [', emin, ':', emax, ']'
        write(outfileindex, 202, advance="no") (k2line_name(i), k2line_stop(i)*Angstrom2atomic, i=1, nk2lines)
        write(outfileindex, 203)k2line_name(nk2lines+1), k2line_stop(nk2lines+1)*Angstrom2atomic

        do i=1, nk2lines-1
           write(outfileindex, 204)k2line_stop(i+1)*Angstrom2atomic, emin, k2line_stop(i+1)*Angstrom2atomic, emax
        enddo
        write(outfileindex, '(a)')'set pm3d interpolate 2,2'
        write(outfileindex, '(2a)')"splot 'dos.dat_l' u 1:2:3 w pm3d"

        write(outfileindex, '(3a)')'set terminal png truecolor enhanced', &
           ' font ", 30" size 1920, 1680'
        write(outfileindex, '(a)')"set output 'spindos_l.png'"
        write(outfileindex, '(a)')"set multiplot layout 3, 1"
        write(outfileindex, '(a)')"set title 'sx'"
        write(outfileindex, '(2a)')"splot 'spindos.dat_l' u 1:2:3 w pm3d "
        write(outfileindex, '(a)')"set title 'sy'"
        write(outfileindex, '(2a)')"splot 'spindos.dat_l' u 1:2:4 w pm3d"
        write(outfileindex, '(a)')"set title 'sz'"
        write(outfileindex, '(2a)')"splot 'spindos.dat_l' u 1:2:5 w pm3d"
 

        close(outfileindex)

     endif

     outfileindex= outfileindex+ 1
     if (cpuid==0) then
        open(unit=outfileindex, file='surfdos_l_only.gnu')
        write(outfileindex, '(a)')"set encoding iso_8859_1"
        write(outfileindex, '(a)')'#set terminal  postscript enhanced color'
        write(outfileindex, '(a)')"#set output 'surfdos_l.eps'"
        write(outfileindex, '(3a)')'#set terminal  pngcairo truecolor enhanced', &
           '  font ", 60" size 1920, 1680'
        write(outfileindex, '(3a)')'set terminal  png truecolor enhanced', &
           ' font ", 60" size 1920, 1680'
        write(outfileindex, '(a)')"set output 'surfdos_l_only.png'"
        write(outfileindex,'(2a)') 'set palette defined (0  "white", ', &
           '6 "red", 20 "black" )'
        write(outfileindex, '(a)')'#set palette rgbformulae 33,13,10'
        write(outfileindex, '(a)')'set style data linespoints'
        write(outfileindex, '(a)')'set size 0.8, 1'
        write(outfileindex, '(a)')'set origin 0.1, 0'
        write(outfileindex, '(a)')'unset ztics'
        write(outfileindex, '(a)')'unset key'
        write(outfileindex, '(a)')'set pointsize 0.8'
        write(outfileindex, '(a)')'set pm3d'
        write(outfileindex, '(a)')'#set view equal xyz'
        write(outfileindex, '(a)')'set view map'
        write(outfileindex, '(a)')'set border lw 3'
        write(outfileindex, '(a)')'#set cbtics font ",48"'
        write(outfileindex, '(a)')'#set xtics font ",48"'
        write(outfileindex, '(a)')'#set ytics font ",48"'
        write(outfileindex, '(a)')'unset cbtics'
        write(outfileindex, '(a)')'#set ylabel font ",48"'
        write(outfileindex, '(a)')'set ylabel "Energy (eV)"'
        write(outfileindex, '(a)')'#set xtics offset 0, -1'
        write(outfileindex, '(a)')'#set ylabel offset -6, 0 '
        write(outfileindex, '(a, f18.5, a)')'set xrange [0: ', maxval(k2len)*Angstrom2atomic, ']'
        write(outfileindex, '(a, f18.5, a, f18.5, a)')'set yrange [', emin, ':', emax, ']'
        write(outfileindex, 202, advance="no") (k2line_name(i), k2line_stop(i)*Angstrom2atomic, i=1, nk2lines)
        write(outfileindex, 203)k2line_name(nk2lines+1), k2line_stop(nk2lines+1)*Angstrom2atomic

        do i=1, nk2lines-1
           write(outfileindex, 204)k2line_stop(i+1)*Angstrom2atomic, emin, k2line_stop(i+1)*Angstrom2atomic, emax
        enddo
        write(outfileindex, '(a)')'set pm3d interpolate 2,2'
        write(outfileindex, '(2a)')"splot 'dos.dat_l' u 1:2:(exp($4)) w pm3d"
        CLOSE(outfileindex)

     endif

     !> write script for gnuplot
     outfileindex= outfileindex+ 1
     if (cpuid==0) then
        open(unit=outfileindex, file='surfdos_r_only.gnu')
        write(outfileindex, '(a)')"set encoding iso_8859_1"
        write(outfileindex, '(a)')'#set terminal  postscript enhanced color'
        write(outfileindex, '(a)')"#set output 'surfdos_r.eps'"
        write(outfileindex, '(3a)')'#set terminal pngcairo truecolor enhanced', &
           '  font ", 60" size 1920, 1680'
        write(outfileindex, '(3a)')'set terminal png truecolor enhanced', &
           ' font ", 60" size 1920, 1680'
        write(outfileindex, '(a)')"set output 'surfdos_r_only.png'"
        write(outfileindex,'(2a)') 'set palette defined (0  "white", ', &
           '6 "red", 20 "black" )'
        write(outfileindex, '(a)')'#set palette rgbformulae 33,13,10'
        write(outfileindex, '(a)')'set style data linespoints'
        write(outfileindex, '(a)')'unset ztics'
        write(outfileindex, '(a)')'unset key'
        write(outfileindex, '(a)')'set pointsize 0.8'
        write(outfileindex, '(a)')'set pm3d'
        write(outfileindex, '(a)')'set border lw 3'
        write(outfileindex, '(a)')'set size 0.8, 1'
        write(outfileindex, '(a)')'set origin 0.1, 0'
        write(outfileindex, '(a)')'#set size ratio -1'
        write(outfileindex, '(a)')'#set view equal xyz'
        write(outfileindex, '(a)')'set view map'
        write(outfileindex, '(a)')'#set cbtics font ",48"'
        write(outfileindex, '(a)')'#set xtics font ",48"'
        write(outfileindex, '(a)')'#set ytics font ",48"'
        write(outfileindex, '(a)')'#set ylabel font ",48"'
        write(outfileindex, '(a)')'set ylabel "Energy (eV)"'
        write(outfileindex, '(a)')'unset cbtics'
        write(outfileindex, '(a)')'#set xtics offset 0, -1'
        write(outfileindex, '(a)')'#set ylabel offset -6, 0 '
        write(outfileindex, '(a, f18.5, a)')'set xrange [0: ', maxval(k2len)*Angstrom2atomic, ']'
        write(outfileindex, '(a, f18.5, a, f18.5, a)')'set yrange [', emin, ':', emax, ']'
        write(outfileindex, 202, advance="no") (k2line_name(i), k2line_stop(i)*Angstrom2atomic, i=1, nk2lines)
        write(outfileindex, 203)k2line_name(nk2lines+1), k2line_stop(nk2lines+1)*Angstrom2atomic

        do i=1, nk2lines-1
           write(outfileindex, 204)k2line_stop(i+1)*Angstrom2atomic, emin, k2line_stop(i+1)*Angstrom2atomic, emax
        enddo
        write(outfileindex, '(a)')'set pm3d interpolate 2,2'
        write(outfileindex, '(2a)')"splot 'dos.dat_r' u 1:2:(exp($4)) w pm3d"
        CLOSE(outfileindex)

     endif


     !> write script for gnuplot
     outfileindex= outfileindex+ 1
     if (cpuid==0) then
        open(unit=outfileindex, file='surfdos_r.gnu')
        write(outfileindex, '(a)')"set encoding iso_8859_1"
        write(outfileindex, '(a)')'#set terminal  postscript enhanced color'
        write(outfileindex, '(a)')"#set output 'surfdos_r.eps'"
        write(outfileindex, '(3a)')'#set terminal pngcairo truecolor enhanced', &
           '  font ", 60" size 1920, 1680'
        write(outfileindex, '(3a)')'set terminal png truecolor enhanced', &
           ' font ", 60" size 1920, 1680'
        write(outfileindex, '(a)')"set output 'surfdos_r.png'"
        write(outfileindex,'(2a)') 'set palette defined (-10 "#194eff", ', &
           '0 "white", 10 "red" )'
        write(outfileindex, '(a)')'#set palette rgbformulae 33,13,10'
        write(outfileindex, '(a)')'set style data linespoints'
        write(outfileindex, '(a)')'unset ztics'
        write(outfileindex, '(a)')'unset key'
        write(outfileindex, '(a)')'set pointsize 0.8'
        write(outfileindex, '(a)')'set pm3d'
        write(outfileindex, '(a)')'set border lw 3'
        write(outfileindex, '(a)')'set size 0.8, 1'
        write(outfileindex, '(a)')'set origin 0.1, 0'
        write(outfileindex, '(a)')'#set size ratio -1'
        write(outfileindex, '(a)')'#set view equal xyz'
        write(outfileindex, '(a)')'set view map'
        write(outfileindex, '(a)')'#set cbtics font ",48"'
        write(outfileindex, '(a)')'#set xtics font ",48"'
        write(outfileindex, '(a)')'#set ytics font ",48"'
        write(outfileindex, '(a)')'#set ylabel font ",48"'
        write(outfileindex, '(a)')'set ylabel "Energy (eV)"'
        write(outfileindex, '(a)')'#set xtics offset 0, -1'
        write(outfileindex, '(a)')'#set ylabel offset -6, 0 '
        write(outfileindex, '(a, f18.5, a)')'set xrange [0: ', maxval(k2len)*Angstrom2atomic, ']'
        write(outfileindex, '(a, f18.5, a, f18.5, a)')'set yrange [', emin, ':', emax, ']'
        write(outfileindex, 202, advance="no") (k2line_name(i), k2line_stop(i)*Angstrom2atomic, i=1, nk2lines)
        write(outfileindex, 203)k2line_name(nk2lines+1), k2line_stop(nk2lines+1)*Angstrom2atomic

        do i=1, nk2lines-1
           write(outfileindex, 204)k2line_stop(i+1)*Angstrom2atomic, emin, k2line_stop(i+1)*Angstrom2atomic, emax
        enddo
        write(outfileindex, '(a)')'set pm3d interpolate 2,2'
        write(outfileindex, '(2a)')"splot 'dos.dat_r' u 1:2:3 w pm3d"

        write(outfileindex, '(3a)')'set terminal png truecolor enhanced', &
           ' font ", 30" size 1920, 1680'
        write(outfileindex, '(a)')"set output 'spindos_r.png'"
        write(outfileindex, '(a)')"set multiplot layout 3, 1"
        write(outfileindex, '(a)')"set title 'sx'"
        write(outfileindex, '(2a)')"splot 'spindos.dat_r' u 1:2:3 w pm3d "
        write(outfileindex, '(a)')"set title 'sy'"
        write(outfileindex, '(2a)')"splot 'spindos.dat_r' u 1:2:4 w pm3d"
        write(outfileindex, '(a)')"set title 'sz'"
        write(outfileindex, '(2a)')"splot 'spindos.dat_r' u 1:2:5 w pm3d"
        close(outfileindex)

     endif

     !> write script for gnuplot
     outfileindex= outfileindex+ 1
     if (cpuid==0) then
        open(unit=outfileindex, file='surfdos_bulk.gnu')
        write(outfileindex, '(a)')"set encoding iso_8859_1"
        write(outfileindex, '(a)')'#set terminal  postscript enhanced color'
        write(outfileindex, '(a)')"#set output 'surfdos_bulk.eps'"
        write(outfileindex, '(3a)')'#set terminal  pngcairo truecolor enhanced', &
           '  font ", 60" size 1920, 1680'
        write(outfileindex, '(3a)')'set terminal  png truecolor enhanced', &
           ' font ", 60" size 1920, 1680'
        write(outfileindex, '(a)')"set output 'surfdos_bulk.png'"
        write(outfileindex,'(2a)') 'set palette defined (-10 "#194eff", ', &
           '0 "white", 10 "red" )'
        write(outfileindex, '(a)')'#set palette rgbformulae 33,13,10'
        write(outfileindex, '(a)')'set style data linespoints'
        write(outfileindex, '(a)')'set size 0.8, 1'
        write(outfileindex, '(a)')'set origin 0.1, 0'
        write(outfileindex, '(a)')'unset ztics'
        write(outfileindex, '(a)')'unset key'
        write(outfileindex, '(a)')'set pointsize 0.8'
        write(outfileindex, '(a)')'set pm3d'
        write(outfileindex, '(a)')'#set view equal xyz'
        write(outfileindex, '(a)')'set view map'
        write(outfileindex, '(a)')'set border lw 3'
        write(outfileindex, '(a)')'#set cbtics font ",48"'
        write(outfileindex, '(a)')'#set xtics font ",48"'
        write(outfileindex, '(a)')'#set ytics font ",48"'
        write(outfileindex, '(a)')'#set ylabel font ",48"'
        write(outfileindex, '(a)')'unset cbtics'
        write(outfileindex, '(a)')'set ylabel "Energy (eV)"'
        write(outfileindex, '(a)')'#set xtics offset 0, -1'
        write(outfileindex, '(a)')'#set ylabel offset -6, 0 '
        write(outfileindex, '(a, f18.5, a)')'set xrange [0: ', maxval(k2len)*Angstrom2atomic, ']'
        write(outfileindex, '(a, f18.5, a, f18.5, a)')'set yrange [', emin, ':', emax, ']'
        write(outfileindex, 202, advance="no") (k2line_name(i), k2line_stop(i)*Angstrom2atomic, i=1, nk2lines)
        write(outfileindex, 203)k2line_name(nk2lines+1), k2line_stop(nk2lines+1)*Angstrom2atomic

        do i=1, nk2lines-1
           write(outfileindex, 204)k2line_stop(i+1)*Angstrom2atomic, emin, k2line_stop(i+1)*Angstrom2atomic, emax
        enddo
        write(outfileindex, '(a)')'set pm3d interpolate 2,2'
        write(outfileindex, '(2a)')"splot 'dos.dat_bulk' u 1:2:(exp($3)) w pm3d"
        CLOSE(outfileindex)

     endif

     202 format('set xtics (',:20('"',A1,'" ',F8.5,','))
     203 format(A1,'" ',F8.5,')')
     204 format('set arrow from ',F8.5,',',F10.5,' to ',F8.5,',',F10.5, ' nohead front lw 3')

#if defined (MPI)
     call mpi_barrier(mpi_cmw, ierr)
#endif

     deallocate( omega)
     deallocate( dos_l)
     deallocate( dos_r)
     deallocate( dos_l_only)
     deallocate( dos_r_only)
     deallocate( dos_l_mpi)
     deallocate( dos_r_mpi)
     deallocate( dos_bulk)
     deallocate( dos_bulk_mpi)
     deallocate(GLL)
     deallocate(GRR)
     deallocate(GB )
     deallocate(H00)
     deallocate(H01)
     deallocate(ones)

  return
  end subroutine surfstat


SUBROUTINE surfstat_jdos

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
!> The jdos part in this subroutine was added by Jianzhou Zhao
!
! Licence : GPL v3

    USE wmpi
    USE para, ONLY: omeganum, omegamin, omegamax, ndim, knv2, k2_path, outfileindex, &
                    BottomOrbitals, TopOrbitals, NBottomOrbitals, NtopOrbitals, stdout, &
                    k2len, Num_wann, eps9, zi, Np, eV2Hartree, Angstrom2atomic
    IMPLICIT NONE

    ! MPI error code
    INTEGER  :: ierr

    ! file id
    INTEGER  :: jdoslfile, jdosrfile, dosbulkfile
    INTEGER  :: doslfile, dosrfile, spindoslfile, spindosrfile

    ! general loop index
    INTEGER  :: i, j, io, iq, iq1, ik1, ikp
    INTEGER  :: Nk_half, imin1, imax1, nw_half
    REAL(DP) :: ktmp(2), eta_broadening, s0(3), s1(3)
    ! string for integer
    CHARACTER(LEN=140) :: ichar, jchar, kchar, fmt

    ! start & stop time
    REAL(DP) :: time_start, time_end, time_togo

    ! Omega
    REAL(DP),    ALLOCATABLE  :: omega(:)
    ! DOS for surface and bulk
    REAL(DP),    ALLOCATABLE  :: dos_l(:,:), dos_r(:,:), dos_bulk(:,:)
    ! DOS for surface ONLY
    REAL(DP),    ALLOCATABLE  :: dos_l_only(:,:), dos_r_only(:,:)
    ! JDOS for both side
    REAL(DP),    ALLOCATABLE  :: jdos_l(:,:), jdos_r(:,:), jdos_l_only(:,:), jdos_r_only(:,:)
    ! Array for MPI
    REAL(DP),    ALLOCATABLE  :: dos_l_mpi(:,:), dos_r_mpi(:,:), dos_bulk_mpi(:,:)
    REAL(DP),    ALLOCATABLE  :: jdos_l_mpi(:,:), jdos_r_mpi(:,:), jdos_l_only_mpi(:,:), jdos_r_only_mpi(:,:)
    ! Green's function
    COMPLEX(DP), ALLOCATABLE  :: GLL(:,:), GRR(:,:), GB (:,:)
    COMPLEX(DP), ALLOCATABLE  :: H00(:,:), H01(:,:)
    ! Unit array
    COMPLEX(DP), ALLOCATABLE  :: ones(:,:)
    ! Spin resolved component
    REAL(DP),    ALLOCATABLE  :: sx_l(:, :), sy_l(:, :), sz_l(:, :)
    REAL(DP),    ALLOCATABLE  :: sx_r(:, :), sy_r(:, :), sz_r(:, :)
    REAL(DP),    ALLOCATABLE  :: sx_l_mpi(:, :), sy_l_mpi(:, :), sz_l_mpi(:, :)
    REAL(DP),    ALLOCATABLE  :: sx_r_mpi(:, :), sy_r_mpi(:, :), sz_r_mpi(:, :)
    COMPLEX(DP), ALLOCATABLE  :: sigma_x(:,:), sigma_y(:,:), sigma_z(:,:)
    COMPLEX(DP), ALLOCATABLE  :: ctemp(:,:)

    ALLOCATE( omega(omeganum) )
    ALLOCATE( dos_l(knv2, omeganum), dos_r(knv2, omeganum), dos_bulk(knv2, omeganum) )
    ALLOCATE( dos_l_only(knv2, omeganum), dos_r_only(knv2, omeganum) )
    ALLOCATE( dos_l_mpi(knv2, omeganum), dos_r_mpi(knv2, omeganum) )
    ALLOCATE( jdos_l(knv2, omeganum), jdos_r(knv2, omeganum) )
    ALLOCATE( jdos_l_only(knv2, omeganum), jdos_r_only(knv2, omeganum) )
    ALLOCATE( dos_bulk_mpi(knv2, omeganum))
    ALLOCATE( jdos_l_mpi(knv2, omeganum), jdos_r_mpi(knv2, omeganum) )
    ALLOCATE( jdos_l_only_mpi(knv2, omeganum), jdos_r_only_mpi(knv2, omeganum) )
    ALLOCATE( GLL(Ndim, Ndim),GRR(Ndim, Ndim),GB (Ndim, Ndim) )
    ALLOCATE( H00(Ndim, Ndim),H01(Ndim, Ndim) )
    ALLOCATE( ones(Ndim, Ndim) )
    ALLOCATE( sx_l(knv2, omeganum), sy_l(knv2, omeganum), sz_l(knv2, omeganum))
    ALLOCATE( sx_l_mpi(knv2, omeganum), sy_l_mpi(knv2, omeganum), sz_l_mpi(knv2, omeganum))
    ALLOCATE( sx_r(knv2, omeganum), sy_r(knv2, omeganum), sz_r(knv2, omeganum))
    ALLOCATE( sx_r_mpi(knv2, omeganum), sy_r_mpi(knv2, omeganum), sz_r_mpi(knv2, omeganum))
    ALLOCATE( sigma_x(ndim,ndim), sigma_y(ndim,ndim), sigma_z(ndim,ndim), ctemp(ndim,ndim))

    omega        = 0d0;      ones         = 0d0
    dos_l        = 0d0;      dos_r        = 0d0
    dos_l_only   = 0d0;      dos_r_only   = 0d0
    dos_l_mpi    = 0d0;      dos_r_mpi    = 0d0
    dos_bulk     = 0d0;      dos_bulk_mpi = 0d0
    GLL          = 0d0;      GRR          = 0d0;      GB           = 0d0
    H00          = 0d0;      H01          = 0d0
    sigma_x      = 0d0;      sigma_y      = 0d0;      sigma_z      = 0d0
    sx_l         = 0d0;      sy_l         = 0d0;      sz_l         = 0d0
    sx_r         = 0d0;      sy_r         = 0d0;      sz_r         = 0d0
    sx_l_mpi     = 0d0;      sy_l_mpi     = 0d0;      sz_l_mpi     = 0d0
    sx_r_mpi     = 0d0;      sy_r_mpi     = 0d0;      sz_r_mpi     = 0d0

    ! Broaden coeffient
    eta_broadening=(omegamax- omegamin)/DBLE(omeganum)
    ! omega list
    DO i = 1, omeganum
        omega(i) = omegamin+(i-1)*eta_broadening
    ENDDO
    eta_broadening = eta_broadening * 3.0d0

    DO i=1,Ndim
        ones(i,i) = 1.0d0
    ENDDO

    !> spin operator matrix
    !> Note, the basis here should be |↑↑↓↓>
    nw_half = Num_wann/2
    DO i=1, Np
        DO j=1, nw_half
            sigma_x( Num_wann*(i-1)+j        , Num_wann*(i-1)+j+nw_half ) =  1.0d0
            sigma_x( Num_wann*(i-1)+j+nw_half, Num_wann*(i-1)+j         ) =  1.0d0
            sigma_y( Num_wann*(i-1)+j        , Num_wann*(i-1)+j+nw_half ) = -zi
            sigma_y( Num_wann*(i-1)+j+nw_half, Num_wann*(i-1)+j         ) =  zi
            sigma_z( Num_wann*(i-1)+j        , Num_wann*(i-1)+j         ) =  1.0d0
            sigma_z( Num_wann*(i-1)+j+nw_half, Num_wann*(i-1)+j+nw_half ) = -1.0d0
        ENDDO
    ENDDO

    time_start = 0d0
    time_end   = 0d0
    DO ikp = 1+cpuid, knv2, num_cpu

        CALL now(time_start)

        ktmp = k2_path(ikp,:)

        CALL ham_qlayer2qlayer(ktmp,h00,h01)

        DO j = 1, omeganum
            !> calculate surface green function
            ! there are two method to calculate surface green's function
            ! the method in 1985 is better, you can find the ref in the
            ! subroutine
            CALL surfgreen_1985(omega(j),GLL,GRR,GB,H00,H01,ones, eta_broadening)
            ! call surfgreen_1984(w,GLL,GRR,H00,H01,ones, eta_broadening)
            ! calculate spectral function
            DO i = 1, NtopOrbitals
                io = TopOrbitals(i)
                dos_l_mpi(ikp, j) = dos_l_mpi(ikp,j) - AIMAG(GLL(io,io))
            ENDDO ! i
            DO i = 1, NBottomOrbitals
                io = Ndim- Num_wann+ BottomOrbitals(i)
                dos_r_mpi(ikp, j) = dos_r_mpi(ikp,j) - AIMAG(GRR(io,io))
            ENDDO ! i
            DO i = 1, Ndim
                dos_bulk_mpi(ikp, j) = dos_bulk_mpi(ikp,j) - AIMAG(GB(i,i))
            ENDDO ! i

            ! Spin resolved sprectrum
            CALL mat_mul(ndim,gll,sigma_x,ctemp)
            DO i = 1, NtopOrbitals
                io = TopOrbitals(i)
                sx_l_mpi(ikp, j) = sx_l_mpi(ikp, j)- AIMAG(ctemp(io,io))
            ENDDO ! i
            CALL mat_mul(ndim,gll,sigma_y,ctemp)
            DO i = 1, NtopOrbitals
                io = TopOrbitals(i)
                sy_l_mpi(ikp, j) = sy_l_mpi(ikp, j)- AIMAG(ctemp(io,io))
            ENDDO !
            CALL mat_mul(ndim,gll,sigma_z,ctemp)
            DO i = 1, NtopOrbitals
                io = TopOrbitals(i)
                sz_l_mpi(ikp, j) = sz_l_mpi(ikp, j)- AIMAG(ctemp(io,io))
            ENDDO ! i
            CALL mat_mul(ndim,grr,sigma_x,ctemp)
            DO i = 1, NtopOrbitals
                io = TopOrbitals(i)
                sx_r_mpi(ikp, j) = sx_r_mpi(ikp, j)- AIMAG(ctemp(io,io))
            ENDDO ! i
            CALL mat_mul(ndim,grr,sigma_y,ctemp)
            DO i = 1, NtopOrbitals
                io = TopOrbitals(i)
                sy_r_mpi(ikp, j) = sy_r_mpi(ikp, j)- AIMAG(ctemp(io,io))
            ENDDO !
            CALL mat_mul(ndim,grr,sigma_z,ctemp)
            DO i = 1, NtopOrbitals
                io = TopOrbitals(i)
                sz_r_mpi(ikp, j) = sz_r_mpi(ikp, j)- AIMAG(ctemp(io,io))
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

!> we do have to do allreduce operation
#ifdef MPI
    CALL mpi_allreduce(dos_l_mpi, dos_l, SIZE(dos_l), mpi_double_precision,&
    mpi_sum, mpi_comm_world, ierr)
    CALL mpi_allreduce(dos_r_mpi, dos_r, SIZE(dos_r), mpi_double_precision,&
    mpi_sum, mpi_comm_world, ierr)
    CALL mpi_allreduce(dos_bulk_mpi, dos_bulk, SIZE(dos_bulk), mpi_double_precision,&
    mpi_sum, mpi_comm_world, ierr)
    CALL mpi_allreduce(sx_l_mpi, sx_l, SIZE(sx_l), mpi_double_precision,&
    mpi_sum, mpi_comm_world, ierr)
    CALL mpi_allreduce(sy_l_mpi, sy_l, SIZE(sy_l), mpi_double_precision,&
    mpi_sum, mpi_comm_world, ierr)
    CALL mpi_allreduce(sz_l_mpi, sz_l, SIZE(sz_l), mpi_double_precision,&
    mpi_sum, mpi_comm_world, ierr)
    CALL mpi_allreduce(sx_r_mpi, sx_r, SIZE(sx_r), mpi_double_precision,&
    mpi_sum, mpi_comm_world, ierr)
    CALL mpi_allreduce(sy_r_mpi, sy_r, SIZE(sy_r), mpi_double_precision,&
    mpi_sum, mpi_comm_world, ierr)
    CALL mpi_allreduce(sz_r_mpi, sz_r, SIZE(sz_r), mpi_double_precision,&
    mpi_sum, mpi_comm_world, ierr)
#else
    dos_l    = dos_l_mpi;     dos_r    = dos_r_mpi;     dos_bulk = dos_bulk_mpi
    sx_l     = sx_l_mpi;      sy_l     = sy_l_mpi;      sz_l     = sz_l_mpi
    sx_r     = sx_r_mpi;      sy_r     = sy_r_mpi;      sz_r     = sz_r_mpi
#endif

    DEALLOCATE( dos_l_mpi, dos_r_mpi, dos_bulk_mpi )
    DEALLOCATE( sx_l_mpi, sy_l_mpi, sz_l_mpi )
    DEALLOCATE( sx_r_mpi, sy_r_mpi, sz_r_mpi )
    DEALLOCATE( sigma_x, sigma_y, sigma_z, ctemp )

    DO ikp=1, knv2
        DO j=1, omeganum
            dos_l_only(ikp, j) = dos_l(ikp, j)- dos_bulk(ikp, j)
            IF (dos_l_only(ikp, j)<0) dos_l_only(ikp, j) = eps9
            dos_r_only(ikp, j) = dos_r(ikp, j)- dos_bulk(ikp, j)
            IF (dos_r_only(ikp, j)<0) dos_r_only(ikp, j) = eps9
        ENDDO
    ENDDO

    outfileindex = outfileindex+ 1
    doslfile     = outfileindex
    outfileindex = outfileindex+ 1
    dosrfile     = outfileindex
    outfileindex = outfileindex+ 1
    dosbulkfile  = outfileindex
    outfileindex = outfileindex+ 1
    spindoslfile = outfileindex
    outfileindex = outfileindex+ 1
    spindosrfile = outfileindex

    ! Write surface state to files
    IF (cpuid .eq. 0) THEN
        OPEN(unit=doslfile    , file='dos.dat_l')
        OPEN(unit=dosrfile    , file='dos.dat_r')
        OPEN(unit=dosbulkfile , file='dos.dat_bulk')
        OPEN(unit=spindoslfile, file='spindos.dat_l')
        OPEN(unit=spindosrfile, file='spindos.dat_r')
        DO ikp = 1, knv2
            DO j = 1, omeganum
                WRITE(doslfile,    2002) k2len(ikp)*Angstrom2atomic, omega(j)/eV2Hartree, dos_l(ikp, j), dos_l_only(ikp, j)
                WRITE(dosrfile,    2002) k2len(ikp)*Angstrom2atomic, omega(j)/eV2Hartree, dos_r(ikp, j), dos_r_only(ikp, j)
                WRITE(dosbulkfile, 2003) k2len(ikp)*Angstrom2atomic, omega(j)/eV2Hartree, dos_bulk(ikp, j)
                s0(1)=sx_l(ikp, j); s0(2)=sy_l(ikp, j); s0(3)=sz_l(ikp, j); 
                call rotate(s0, s1)
                WRITE(spindoslfile,2001) k2len(ikp)*Angstrom2atomic, omega(j)/eV2Hartree, s1
                s0(1)=sx_r(ikp, j); s0(2)=sy_r(ikp, j); s0(3)=sz_r(ikp, j); 
                call rotate(s0, s1)
                WRITE(spindosrfile,2001) k2len(ikp)*Angstrom2atomic, omega(j)/eV2Hartree, s1
            ENDDO
            WRITE(doslfile    , *)
            WRITE(dosrfile    , *)
            WRITE(dosbulkfile , *)
            WRITE(spindoslfile, *)
            WRITE(spindosrfile, *)
        ENDDO
        CLOSE(doslfile)
        CLOSE(dosrfile)
        CLOSE(dosbulkfile)
        CLOSE(spindoslfile)
        CLOSE(spindosrfile)
    ENDIF

    ! Begin to calculate JDOS from DOS before
    Nk_half         = (knv2-1)/2
    jdos_l_mpi      = 0.0d0
    jdos_r_mpi      = 0.0d0
    jdos_l_only_mpi = 0.0d0
    jdos_r_only_mpi = 0.0d0
    DO j = 1+cpuid, omeganum,num_cpu
        DO iq= 1, knv2
            iq1 = iq - Nk_half
            imin1 = MAX(-Nk_half-iq1, -Nk_half)+ Nk_half+ 1
            imax1 = MIN( Nk_half-iq1,  Nk_half)+ Nk_half+ 1
            DO ik1 = imin1, imax1
                jdos_l_mpi(iq,j) = jdos_l_mpi(iq,j)+ dos_l(ik1, j)* dos_l(ik1+iq1, j)
                jdos_r_mpi(iq,j) = jdos_r_mpi(iq,j)+ dos_r(ik1, j)* dos_r(ik1+iq1, j)
                jdos_l_only_mpi(iq,j) = jdos_l_only_mpi(iq,j)+ dos_l_only(ik1, j)* dos_l_only(ik1+iq1, j)
                jdos_r_only_mpi(iq,j) = jdos_r_only_mpi(iq,j)+ dos_r_only(ik1, j)* dos_r_only(ik1+iq1, j)
            ENDDO !ik1
        ENDDO !iq
    ENDDO

#ifdef MPI
    CALL mpi_reduce(jdos_l_mpi, jdos_l, SIZE(jdos_l), mpi_double_precision,&
    mpi_sum, 0, mpi_comm_world, ierr)
    CALL mpi_reduce(jdos_r_mpi, jdos_r, SIZE(jdos_r), mpi_double_precision,&
    mpi_sum, 0, mpi_comm_world, ierr)
    CALL mpi_reduce(jdos_l_only_mpi, jdos_l_only, SIZE(jdos_l_only), mpi_double_precision,&
    mpi_sum, 0, mpi_comm_world, ierr)
    CALL mpi_reduce(jdos_r_only_mpi, jdos_r_only, SIZE(jdos_r_only), mpi_double_precision,&
    mpi_sum, 0, mpi_comm_world, ierr)
#else
    jdos_l      = jdos_l_mpi;         jdos_r      = jdos_r_mpi
    jdos_l_only = jdos_l_only_mpi;    jdos_r_only = jdos_r_only_mpi
#endif

    outfileindex = outfileindex+ 1
    jdoslfile    = outfileindex
    outfileindex = outfileindex+ 1
    jdosrfile    = outfileindex

    IF(cpuid.eq.0) THEN
        OPEN(unit=jdoslfile, file='dos.jdat_l')
        OPEN(unit=jdosrfile, file='dos.jdat_r')
        write(jdoslfile, '("#", a12, 3a17)') ' k(1/Ang)', ' E(eV)', 'jdos_l', 'jdos_l_only'
        write(jdosrfile, '("#", a12, 3a17)') ' k(1/Ang)', ' E(eV)', 'jdos_l', 'jdos_l_only'
        DO ikp = 1, knv2
            DO j = 1, omeganum
                WRITE(jdoslfile, 2002) k2len(ikp)*Angstrom2atomic, omega(j)/eV2Hartree, jdos_l(ikp, j), jdos_l_only(ikp, j)
                WRITE(jdosrfile, 2002) k2len(ikp)*Angstrom2atomic, omega(j)/eV2Hartree, jdos_r(ikp, j), jdos_r_only(ikp, j)
            ENDDO
            WRITE(jdoslfile, *)
            WRITE(jdosrfile, *)
        ENDDO
        CLOSE(jdoslfile)
        CLOSE(jdosrfile)
        WRITE(stdout,*)
        WRITE(stdout,*)'calculate joint density of state successfully'
    ENDIF

    DEALLOCATE( ones, omega )
    DEALLOCATE( dos_l, dos_r, dos_l_only, dos_r_only, dos_bulk )
    DEALLOCATE( GLL, GRR, GB )
    DEALLOCATE( H00, H01 )
    DEALLOCATE( jdos_l, jdos_r, jdos_l_only, jdos_r_only )
    DEALLOCATE( jdos_l_mpi, jdos_r_mpi, jdos_l_only_mpi, jdos_r_only_mpi )
    DEALLOCATE( sx_l, sy_l, sz_l )
    DEALLOCATE( sx_r, sy_r, sz_r )

2001 FORMAT(5(1X,E16.8))
2002 FORMAT(4(1X,E16.8))
2003 FORMAT(3(1X,E16.8))

END SUBROUTINE surfstat_jdos
