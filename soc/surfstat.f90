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
     
     integer :: ierr

     integer :: doslfile, dosrfile, dosbulkfile

     ! general loop index
     integer :: i, j, io

     ! kpoint loop index
     integer :: ikp

     real(dp) :: emin
     real(dp) :: emax
     real(dp) :: k(2), w

     real(dp) :: time_start, time_end

     real(dp), allocatable :: omega(:)
     
     real(dp), allocatable :: dos_l(:,:)
     real(dp), allocatable :: dos_r(:,:)
     real(dp), allocatable :: dos_l_only(:,:)
     real(dp), allocatable :: dos_r_only(:,:)
     real(dp), allocatable :: dos_l_mpi(:,:)
     real(dp), allocatable :: dos_r_mpi(:,:)
     real(dp), allocatable :: dos_bulk(:,:)
     real(dp), allocatable :: dos_bulk_mpi(:,:)

     complex(dp), allocatable :: GLL(:,:)
     complex(dp), allocatable :: GRR(:,:)
     complex(dp), allocatable :: GB (:,:)
     complex(dp), allocatable :: H00(:,:)
     complex(dp), allocatable :: H01(:,:)
     complex(dp), allocatable :: ones(:,:)

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

    !   if (i.gt.1) then
    !      t1=t1+temp
    !   endif
    !   k2len(i)= t1
    !   print *, i, t1
    !enddo



     allocate( omega(omeganum))
     allocate( dos_l(knv2, omeganum))
     allocate( dos_r(knv2, omeganum))
     allocate( dos_l_only(knv2, omeganum))
     allocate( dos_r_only(knv2, omeganum))
     allocate( dos_l_mpi(knv2, omeganum))
     allocate( dos_r_mpi(knv2, omeganum))
     allocate( dos_bulk(knv2, omeganum))
     allocate( dos_bulk_mpi(knv2, omeganum))
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


     allocate(GLL(Ndim, Ndim))
     allocate(GRR(Ndim, Ndim))
     allocate(GB (Ndim, Ndim))
     allocate(H00(Ndim, Ndim))
     allocate(H01(Ndim, Ndim))
     allocate(ones(Ndim, Ndim))
     GLL= 0d0
     GRR= 0d0
     GB = 0d0
     H00= 0d0
     H01= 0d0
     ones= 0d0

     do i=1,Ndim
        ones(i,i)=1.0d0
     enddo

     time_start= 0d0
     time_end= 0d0
     do ikp= 1+cpuid, knv2, num_cpu
        if (cpuid==0.and. mod(ikp/num_cpu, 10)==0) &
           write(stdout, *) 'SurfaceSS, ik', ikp, 'Nk', knv2, 'time left', &
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
             call surfgreen_1985(w,GLL,GRR,GB,H00,H01,ones)
           ! call surfgreen_1984(w,GLL,GRR,H00,H01,ones)

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
        enddo ! j
        call now(time_end)
     enddo ! ikp

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

     !> deal with phonon system
     !> for phonon system, omega should be changed to omega^2
     if (index(Particle,'phonon')/=0) then
        do i= 1, omeganum
           omega(i)= sqrt(omega(i))
        enddo
     endif

     if (cpuid.eq.0)then
        open (unit=doslfile, file='dos.dat_l')
        open (unit=dosrfile, file='dos.dat_r')
        open (unit=dosbulkfile, file='dos.dat_bulk')
        do ikp=1, knv2
           do j=1, omeganum 
              write(doslfile, '(30f16.8)')k2len(ikp), omega(j), dos_l(ikp, j), log(dos_l_only(ikp, j))
              write(dosrfile, '(30f16.8)')k2len(ikp), omega(j), dos_r(ikp, j), log(dos_r_only(ikp, j))
              write(dosbulkfile, '(30f16.8)')k2len(ikp), omega(j), log(abs(dos_bulk(ikp, j))+eps9)
           enddo
           write(doslfile, *) ' '
           write(dosrfile, *) ' '
           write(dosbulkfile, *) ' '
        enddo
        close(doslfile)
        close(dosrfile)
        close(dosbulkfile)

        write(stdout,*)'ndim',ndim
        write(stdout,*) 'knv2,omeganum,eta',knv2, omeganum, eta
        write(stdout,*)'calculate density of state successfully'    
     endif


     emin= minval(omega)
     emax= maxval(omega)
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
        write(outfileindex, '(a, f18.5, a)')'set xrange [0: ', maxval(k2len), ']'
        write(outfileindex, '(a, f18.5, a, f18.5, a)')'set yrange [', emin, ':', emax, ']'
        write(outfileindex, 202, advance="no") (k2line_name(i), k2line_stop(i), i=1, nk2lines)
        write(outfileindex, 203)k2line_name(nk2lines+1), k2line_stop(nk2lines+1)

        do i=1, nk2lines-1
           write(outfileindex, 204)k2line_stop(i+1), emin, k2line_stop(i+1), emax
        enddo
        write(outfileindex, '(a)')'set pm3d interpolate 2,2'
        write(outfileindex, '(2a)')"splot 'dos.dat_l' u 1:2:3 w pm3d"
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
        write(outfileindex, '(a, f18.5, a)')'set xrange [0: ', maxval(k2len), ']'
        write(outfileindex, '(a, f18.5, a, f18.5, a)')'set yrange [', emin, ':', emax, ']'
        write(outfileindex, 202, advance="no") (k2line_name(i), k2line_stop(i), i=1, nk2lines)
        write(outfileindex, 203)k2line_name(nk2lines+1), k2line_stop(nk2lines+1)

        do i=1, nk2lines-1
           write(outfileindex, 204)k2line_stop(i+1), emin, k2line_stop(i+1), emax
        enddo
        write(outfileindex, '(a)')'set pm3d interpolate 2,2'
        write(outfileindex, '(2a)')"splot 'dos.dat_l' u 1:2:(exp($4)) w pm3d"
        close(outfileindex)

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
        write(outfileindex, '(a, f18.5, a)')'set xrange [0: ', maxval(k2len), ']'
        write(outfileindex, '(a, f18.5, a, f18.5, a)')'set yrange [', emin, ':', emax, ']'
        write(outfileindex, 202, advance="no") (k2line_name(i), k2line_stop(i), i=1, nk2lines)
        write(outfileindex, 203)k2line_name(nk2lines+1), k2line_stop(nk2lines+1)

        do i=1, nk2lines-1
           write(outfileindex, 204)k2line_stop(i+1), emin, k2line_stop(i+1), emax
        enddo
        write(outfileindex, '(a)')'set pm3d interpolate 2,2'
        write(outfileindex, '(2a)')"splot 'dos.dat_r' u 1:2:(exp($4)) w pm3d"
        close(outfileindex)

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
        write(outfileindex, '(a, f18.5, a)')'set xrange [0: ', maxval(k2len), ']'
        write(outfileindex, '(a, f18.5, a, f18.5, a)')'set yrange [', emin, ':', emax, ']'
        write(outfileindex, 202, advance="no") (k2line_name(i), k2line_stop(i), i=1, nk2lines)
        write(outfileindex, 203)k2line_name(nk2lines+1), k2line_stop(nk2lines+1)

        do i=1, nk2lines-1
           write(outfileindex, 204)k2line_stop(i+1), emin, k2line_stop(i+1), emax
        enddo
        write(outfileindex, '(a)')'set pm3d interpolate 2,2'
        write(outfileindex, '(2a)')"splot 'dos.dat_r' u 1:2:3 w pm3d"
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
        write(outfileindex, '(a, f18.5, a)')'set xrange [0: ', maxval(k2len), ']'
        write(outfileindex, '(a, f18.5, a, f18.5, a)')'set yrange [', emin, ':', emax, ']'
        write(outfileindex, 202, advance="no") (k2line_name(i), k2line_stop(i), i=1, nk2lines)
        write(outfileindex, 203)k2line_name(nk2lines+1), k2line_stop(nk2lines+1)

        do i=1, nk2lines-1
           write(outfileindex, 204)k2line_stop(i+1), emin, k2line_stop(i+1), emax
        enddo
        write(outfileindex, '(a)')'set pm3d interpolate 2,2'
        write(outfileindex, '(2a)')"splot 'dos.dat_bulk' u 1:2:(exp($3)) w pm3d"
        close(outfileindex)

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


