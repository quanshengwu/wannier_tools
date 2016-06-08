!+---------+---------+---------+---------+---------+---------+--------+!
! this subroutine is used to calculate surface state using             !
! green's function method  ---  J.Phys.F.Met.Phys.15(1985)851-858      !
! 
! History:
!         by Quan Sheng Wu on 4/20/2010                                !
!            mpi version      4/21/2010
!            change Kb to K=(Ka+Kb)/3 direction 4/22/2010
!            Quansheng Wu on Jan 30 2015 at ETH Zurich
! Copyright (c) 2010 QuanSheng Wu. All rights reserved.
!+---------+---------+---------+---------+---------+---------+--------+!
  subroutine surfstat

     use wmpi
     use para
     implicit none
     
     integer :: ierr

     ! general loop index
     integer :: i,j 

     ! kpoint loop index
     integer :: ikp

     real(dp) :: emin
     real(dp) :: emax
     real(dp) :: k(2), w

     real(dp), allocatable :: omega(:)
     
     real(dp), allocatable :: dos_l(:,:)
     real(dp), allocatable :: dos_r(:,:)
     real(dp), allocatable :: dos_l_mpi(:,:)
     real(dp), allocatable :: dos_r_mpi(:,:)

     complex(dp), allocatable :: GLL(:,:)
     complex(dp), allocatable :: GRR(:,:)
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
     allocate( dos_l_mpi(knv2, omeganum))
     allocate( dos_r_mpi(knv2, omeganum))
     omega=0d0
     dos_l=0d0
     dos_r=0d0
     dos_l_mpi=0d0
     dos_r_mpi=0d0

     eta=(omegamax- omegamin)/dble(omeganum)*1.5d0

     do i= 1, omeganum
        omega(i)=omegamin+(i-1)*(omegamax-omegamin)/dble(omeganum)
     enddo

     allocate(GLL(Ndim, Ndim))
     allocate(GRR(Ndim, Ndim))
     allocate(H00(Ndim, Ndim))
     allocate(H01(Ndim, Ndim))
     allocate(ones(Ndim, Ndim))
     GLL= 0d0
     GRR= 0d0
     H00= 0d0
     H01= 0d0
     ones= 0d0

     do i=1,Ndim
        ones(i,i)=1.0d0
     enddo

     do ikp= 1+cpuid, knv2, num_cpu
        if (cpuid==0) write(stdout, *) 'SurfaceSS, ik', ikp, 'Nk', knv2
        k= k2_path(ikp,:)

        !> get the hopping matrix between two principle layers
        call ham_qlayer2qlayer(k,H00,H01) 

        do j = 1, omeganum
           w=omega(j)

           !> calculate surface green function
           ! there are two method to calculate surface green's function 
           ! the method in 1985 is better, you can find the ref in the
           ! subroutine
           call surfgreen_1985(w,GLL,GRR,H00,H01,ones)
           !call surfgreen_1984(w,GLL,GRR,H00,H01,ones)

           ! calculate spectral function
           do i= 1, ndim
              dos_l(ikp, j)=dos_l(ikp,j)- aimag(GLL(i,i))
              dos_r(ikp, j)=dos_r(ikp,j)- aimag(GRR(i,i))
           enddo ! i
        enddo ! j
     enddo ! ikp

     !> we don't have to do allreduce operation
     call mpi_reduce(dos_l, dos_l_mpi, size(dos_l), mpi_double_precision,&
                     mpi_sum, 0, mpi_comm_world, ierr)
     call mpi_reduce(dos_r, dos_r_mpi, size(dos_r), mpi_double_precision,&
                     mpi_sum, 0, mpi_comm_world, ierr)
     dos_l=log(dos_l_mpi)
     dos_r=log(dos_r_mpi)

     if (cpuid.eq.0)then
        open (unit=12, file='dos.dat_l')
        open (unit=13, file='dos.dat_r')
        do ikp=1, knv2
           do j=1, omeganum 
              write(12, '(3f16.8)')k2len(ikp), omega(j), dos_l(ikp, j)
              write(13, '(3f16.8)')k2len(ikp), omega(j), dos_r(ikp, j)
           enddo
           write(12, *) ' '
           write(13, *) ' '
        enddo
        close(12)
        close(13)
        write(stdout,*)'ndim',ndim
        write(stdout,*) 'knv2,omeganum,eta',knv2, omeganum, eta
        write(stdout,*)'calculate density of state successfully'    
     endif

     emin= minval(omega)
     emax= maxval(omega)
     !> write script for gnuplot
     if (cpuid==0) then
        open(unit=116, file='surfdos_l.gnu')
        write(116, '(a)')"set encoding iso_8859_1"
        write(116, '(a)')'#set terminal  postscript enhanced color'
        write(116, '(a)')"#set output 'surfdos_l.eps'"
        write(116, '(3a)')'#set terminal  pngcairo truecolor enhanced', &
           '# font ", 60" size 1920, 1680'
        write(116, '(3a)')'set terminal  png truecolor enhanced', &
           ' font ", 60" size 1920, 1680'
        write(116, '(a)')"set output 'surfdos_l.png'"
        write(116,'(2a)') 'set palette defined (-10 "#194eff", ', &
           '0 "white", 10 "red" )'
        write(116, '(a)')'#set palette rgbformulae 33,13,10'
        write(116, '(a)')'set style data linespoints'
        write(116, '(a)')'set size 0.8, 1'
        write(116, '(a)')'set origin 0.1, 0'
        write(116, '(a)')'unset ztics'
        write(116, '(a)')'unset key'
        write(116, '(a)')'set pointsize 0.8'
        write(116, '(a)')'set pm3d'
        write(116, '(a)')'#set view equal xyz'
        write(116, '(a)')'set view map'
        write(116, '(a)')'set border lw 3'
        write(116, '(a)')'#set cbtics font ",48"'
        write(116, '(a)')'#set xtics font ",48"'
        write(116, '(a)')'#set ytics font ",48"'
        write(116, '(a)')'#set ylabel font ",48"'
        write(116, '(a)')'set ylabel "Energy (eV)"'
        write(116, '(a)')'#set xtics offset 0, -1'
        write(116, '(a)')'#set ylabel offset -6, 0 '
        write(116, '(a, f8.5, a)')'set xrange [0: ', maxval(k2len), ']'
        write(116, '(a, f8.5, a, f8.5, a)')'set yrange [', emin, ':', emax, ']'
        write(116, 202, advance="no") (k2line_name(i), k2line_stop(i), i=1, nk2lines)
        write(116, 203)k2line_name(nk2lines+1), k2line_stop(nk2lines+1)

        do i=1, nk2lines-1
           write(116, 204)k2line_stop(i+1), emin, k2line_stop(i+1), emax
        enddo
        write(116, '(a)')'set pm3d interpolate 2,2'
        write(116, '(2a)')"splot 'dos.dat_l' u 1:2:3 w pm3d"
        close(116)

     endif

     !> write script for gnuplot
     if (cpuid==0) then
        open(unit=117, file='surfdos_r.gnu')
        write(117, '(a)')"set encoding iso_8859_1"
        write(117, '(a)')'#set terminal  postscript enhanced color'
        write(117, '(a)')"#set output 'surfdos_r.eps'"
        write(117, '(3a)')'#set terminal pngcairo truecolor enhanced', &
           '# font ", 60" size 1920, 1680'
        write(117, '(3a)')'set terminal png truecolor enhanced', &
           ' font ", 60" size 1920, 1680'
        write(117, '(a)')"set output 'surfdos_r.png'"
        write(117,'(2a)') 'set palette defined (-10 "#194eff", ', &
           '0 "white", 10 "red" )'
        write(117, '(a)')'#set palette rgbformulae 33,13,10'
        write(117, '(a)')'set style data linespoints'
        write(117, '(a)')'unset ztics'
        write(117, '(a)')'unset key'
        write(117, '(a)')'set pointsize 0.8'
        write(117, '(a)')'set pm3d'
        write(117, '(a)')'set border lw 3'
        write(117, '(a)')'set size 0.8, 1'
        write(117, '(a)')'set origin 0.1, 0'
        write(117, '(a)')'#set size ratio -1'
        write(117, '(a)')'#set view equal xyz'
        write(117, '(a)')'set view map'
        write(117, '(a)')'#set cbtics font ",48"'
        write(117, '(a)')'#set xtics font ",48"'
        write(117, '(a)')'#set ytics font ",48"'
        write(117, '(a)')'#set ylabel font ",48"'
        write(117, '(a)')'set ylabel "Energy (eV)"'
        write(117, '(a)')'#set xtics offset 0, -1'
        write(117, '(a)')'#set ylabel offset -6, 0 '
        write(117, '(a, f8.5, a)')'set xrange [0: ', maxval(k2len), ']'
        write(117, '(a, f8.5, a, f8.5, a)')'set yrange [', emin, ':', emax, ']'
        write(117, 202, advance="no") (k2line_name(i), k2line_stop(i), i=1, nk2lines)
        write(117, 203)k2line_name(nk2lines+1), k2line_stop(nk2lines+1)

        do i=1, nk2lines-1
           write(117, 204)k2line_stop(i+1), emin, k2line_stop(i+1), emax
        enddo
        write(117, '(a)')'set pm3d interpolate 2,2'
        write(117, '(2a)')"splot 'dos.dat_r' u 1:2:3 w pm3d"
        close(117)

     endif


     202 format('set xtics (',:20('"',A1,'" ',F8.5,','))
     203 format(A1,'" ',F8.5,')')
     204 format('set arrow from ',F8.5,',',F10.5,' to ',F8.5,',',F10.5, ' nohead front lw 3')
 

  return   
  end subroutine
