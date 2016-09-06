!+---------+---------+---------+---------+---------+---------+--------+!
! this subroutine is used to calculate surface state using             !
! green's function method  ---  J.Phys.F.Met.Phys.15(1985)851-858      !
! 
! History:
!         by Quan Sheng Wu on 4/20/2010                                !
!            mpi version      4/21/2010
!         Ca3PbO version     11/12/2011
!+---------+---------+---------+---------+---------+---------+--------+!
  subroutine fermiarc


     use wmpi
     use para
     implicit none
     


     integer :: ierr

     ! general loop index
     integer :: i,j 

     ! kpoint loop index
     integer :: ikp

     real(dp) :: k(2)

     real(dp) :: omega
     real(dp) :: k1min_shape, k1max_shape, k2min_shape, k2max_shape

     real(dp), allocatable :: k12(:,:)
     real(dp), allocatable :: k12_shape(:,:)
     
     real(dp), allocatable :: dos_l(:)
     real(dp), allocatable :: dos_l_mpi(:)
     real(dp), allocatable :: dos_r(:)
     real(dp), allocatable :: dos_r_mpi(:)

     complex(dp), allocatable :: ones(:,:)
     complex(dp), allocatable :: GLL(:,:), GRR(:,:)
     complex(dp), allocatable :: H00(:,:), H01(:,:)


     allocate( k12(2, nk1*nk2))
     allocate( k12_shape(2, nk1*nk2))
     allocate( dos_l(nk1*nk2))
     allocate( dos_l_mpi(nk1*nk2))
     allocate( dos_r(nk1*nk2))
     allocate( dos_r_mpi(nk1*nk2))
     allocate( GLL(ndim,ndim), GRR(ndim,ndim))
     k12=0d0
     k12_shape=0d0
     dos_l=0d0
     dos_l_mpi=0d0
     dos_r=0d0
     dos_r_mpi=0d0

     ikp=0
     do i= 1, nk1
        do j= 1, nk2
           ikp=ikp+1
           k12(:, ikp)=K2D_start+ (i-1)*K2D_vec1/dble(nk1-1) &
                      + (j-1)*K2D_vec2/dble(nk2-1)
           k12_shape(:, ikp)= k12(1, ikp)* Ka2+ k12(2, ikp)* Kb2
        enddo
     enddo

     k1min_shape= minval(k12_shape(1, :))
     k2min_shape= minval(k12_shape(2, :))
     k1max_shape= maxval(k12_shape(1, :))
     k2max_shape= maxval(k12_shape(2, :))

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


     omega = E_arc
     eta= eta_arc

     do ikp= 1+cpuid, nk1*nk2, num_cpu
        if (cpuid==0) write(stdout, *) 'Arc, ik, knv2', ikp, nk1*nk2
        k(1)= k12(1, ikp)
        k(2)= k12(2, ikp)

        call ham_qlayer2qlayer(k,H00,H01) 

        !> calculate surface green function
        ! there are two method to calculate surface green's function 
        ! the method in 1985 is better, you can find the ref in the
        ! subroutine
        call surfgreen_1985(omega,GLL,GRR,H00,H01,ones)
        ! call surfgreen_1984(omega,GLL,GRR,H00,H01,ones)


        ! calculate spectral function
        do i= 1, ndim
           dos_l(ikp)=dos_l(ikp)- aimag(GLL(i,i))
           dos_r(ikp)=dos_r(ikp)- aimag(GRR(i,i))
        enddo
     enddo

     !> we don't have to do allreduce operation
     call mpi_reduce(dos_l, dos_l_mpi, size(dos_l),mpi_double_precision,&
                     mpi_sum, 0, mpi_comm_world, ierr)
     call mpi_reduce(dos_r, dos_r_mpi, size(dos_r),mpi_double_precision,&
                     mpi_sum, 0, mpi_comm_world, ierr)

     if (cpuid.eq.0)then
        open (unit=12, file='arc.dat_l')
        open (unit=13, file='arc.dat_r')
        do ikp=1, nk1*nk2
           write(12, '(3f16.8)')k12_shape(:, ikp), log(dos_l_mpi(ikp))
           if (mod(ikp, nk2)==0) write(12, *)' '
           write(13, '(3f16.8)')k12_shape(:, ikp), log(dos_r_mpi(ikp))
           if (mod(ikp, nk2)==0) write(13, *)' '
        enddo
        close(12)
        close(13)
        write(stdout,*)'ndim',ndim
        write(stdout,*) 'Nk1,Nk2,eta',Nk1, Nk2, eta
        write(stdout,*)'calculate density of state successfully'    
     endif

     !> write script for gnuplot
     if (cpuid==0) then
        open(unit=114, file='arc_l.gnu')
        write(114, '(a)')"set encoding iso_8859_1"
        write(114, '(a)')'#set terminal  postscript enhanced color'
        write(114, '(a)')"#set output 'arc_l.eps'"
        write(114, '(3a)')'#set terminal pngcairo truecolor enhanced', &
           '  font ",60" size 1920, 1680'
        write(114, '(3a)')'set terminal png truecolor enhanced', &
           ' font ",60" size 1920, 1680'
        write(114, '(a)')"set output 'arc_l.png'"
        write(114,'(2a)') 'set palette defined ( -10 "#194eff", ', &
           '0 "white", 10 "red" )'
        write(114, '(a)')'#set palette rgbformulae 33,13,10'
        write(114, '(a)')'unset ztics'
        write(114, '(a)')'unset key'
        write(114, '(a)')'set pm3d'
        write(114, '(a)')'set border lw 6'
        write(114, '(a)')'set size ratio -1'
        write(114, '(a)')'set view map'
        write(114, '(a)')'set xtics'
        write(114, '(a)')'set ytics'
        write(114, '(a)')'set xlabel "K_1 (1/{\305})"'
        write(114, '(a)')'set ylabel "K_2 (1/{\305})"'
        write(114, '(a)')'set colorbox'
        write(114, '(a, f8.5, a, f8.5, a)')'set xrange [', k1min_shape, ':', k1max_shape, ']'
        write(114, '(a, f8.5, a, f8.5, a)')'set yrange [', k2min_shape, ':', k2max_shape, ']'
        write(114, '(a)')'set pm3d interpolate 2,2'
        write(114, '(2a)')"splot 'arc.dat_l' u 1:2:3 w pm3d"

        close(114)
     endif

     !> write script for gnuplot
     if (cpuid==0) then
        open(unit=115, file='arc_r.gnu')
        write(115, '(a)')"set encoding iso_8859_1"
        write(115, '(a)')'#set terminal  postscript enhanced color'
        write(115, '(a)')"#set output 'arc_r.eps'"
        write(115, '(3a)')'#set terminal pngcairo truecolor enhanced', &
           '  font ",60" size 1920, 1680'
        write(115, '(3a)')'set terminal png truecolor enhanced', &
           ' font ",60" size 1920, 1680'
        write(115, '(a)')"set output 'arc_r.png'"
        write(115,'(2a)') 'set palette defined ( -10 "#194eff", ', &
           '0 "white", 10 "red" )'
        write(115, '(a)')'#set palette rgbformulae 33,13,10'
        write(115, '(a)')'unset ztics'
        write(115, '(a)')'unset key'
        write(115, '(a)')'set pm3d'
        write(115, '(a)')'set border lw 6'
        write(115, '(a)')'set size ratio -1'
        write(115, '(a)')'set view map'
        write(115, '(a)')'set xtics'
        write(115, '(a)')'set ytics'
        write(115, '(a)')'set xlabel "K_1 (1/{\305})"'
        write(115, '(a)')'set ylabel "K_2 (1/{\305})"'
        write(115, '(a)')'set colorbox'
        write(115, '(a, f8.5, a, f8.5, a)')'set xrange [', k1min_shape, ':', k1max_shape, ']'
        write(115, '(a, f8.5, a, f8.5, a)')'set yrange [', k2min_shape, ':', k2max_shape, ']'
        write(115, '(a)')'set pm3d interpolate 2,2'
        write(115, '(2a)')"splot 'arc.dat_r' u 1:2:3 w pm3d"

        close(115)
     endif


  return   
  end subroutine fermiarc


  !> calculation joint density of state 
  !> refs :http://science.sciencemag.org/content/sci/suppl/2016/03/09/351.6278.1184.DC1/Inoue.SM.pdf

  subroutine fermiarc_jdos


     use wmpi
     use para
     implicit none
     


     integer :: ierr

     ! general loop index
     integer :: i,j 

     ! kpoint loop index
     integer :: ikp, ik1, ik2, iq, Nk1_half, Nk2_half
     integer :: imin1, imax1, imin2, imax2, iq1, iq2

     real(dp) :: k(2), dis

     real(dp) :: omega
     real(dp) :: k1min_shape, k1max_shape, k2min_shape, k2max_shape

     integer , allocatable :: ik12(:,:)
     real(dp), allocatable :: k12(:,:)
     real(dp), allocatable :: k12_shape(:,:)
     
     real(dp), allocatable :: dos_l(:,:)
     real(dp), allocatable :: dos_l_mpi(:,:)
     real(dp), allocatable :: dos_r(:,:)
     real(dp), allocatable :: dos_r_mpi(:,:)
     real(dp), allocatable :: jdos_l(:)
     real(dp), allocatable :: jdos_l_mpi(:)
     real(dp), allocatable :: jdos_r(:)
     real(dp), allocatable :: jdos_r_mpi(:)

     complex(dp), allocatable :: ones(:,:)
     complex(dp), allocatable :: GLL(:,:), GRR(:,:)
     complex(dp), allocatable :: H00(:,:), H01(:,:)

     Nk1_half= (Nk1-1)/2
     Nk2_half= (Nk2-1)/2

     allocate( ik12(2, nk1*nk2))
     allocate( k12(2, nk1*nk2))
     allocate( k12_shape(2, nk1*nk2))
     allocate( dos_l(nk1,nk2))
     allocate( dos_l_mpi(nk1,nk2))
     allocate( dos_r(nk1,nk2))
     allocate( dos_r_mpi(nk1,nk2))
     allocate( jdos_l(nk1*nk2))
     allocate( jdos_l_mpi(nk1*nk2))
     allocate( jdos_r(nk1*nk2))
     allocate( jdos_r_mpi(nk1*nk2))
     allocate( GLL(ndim,ndim), GRR(ndim,ndim))
     ik12=0
     k12=0d0
     k12_shape=0d0
     dos_l=0d0
     dos_l_mpi=0d0
     dos_r=0d0
     dos_r_mpi=0d0
     jdos_l=0d0
     jdos_l_mpi=eps9
     jdos_r=0d0
     jdos_r_mpi=eps9

     ikp=0
     do i= 1, nk1
        do j= 1, nk2
           ikp=ikp+1
           ik12(1, ikp)= i
           ik12(2, ikp)= j
           k12(:, ikp)=K2D_start+ (i-1)*K2D_vec1/dble(nk1-1) &
                      + (j-1)*K2D_vec2/dble(nk2-1)
           k12_shape(:, ikp)= k12(1, ikp)* Ka2+ k12(2, ikp)* Kb2
        enddo
     enddo

     k1min_shape= minval(k12_shape(1, :))
     k2min_shape= minval(k12_shape(2, :))
     k1max_shape= maxval(k12_shape(1, :))
     k2max_shape= maxval(k12_shape(2, :))

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


     omega = E_arc
     eta= eta_arc

     do ikp= 1+cpuid, nk1*nk2, num_cpu
        if (cpuid==0) write(stdout, *) 'Arc, ik, knv2', ikp, nk1*nk2
        k(1)= k12(1, ikp)
        k(2)= k12(2, ikp)

        call ham_qlayer2qlayer(k,H00,H01) 

        !> calculate surface green function
        ! there are two method to calculate surface green's function 
        ! the method in 1985 is better, you can find the ref in the
        ! subroutine
        call surfgreen_1985(omega,GLL,GRR,H00,H01,ones)
        ! call surfgreen_1984(omega,GLL,GRR,H00,H01,ones)

        ik1= ik12(1, ikp)
        ik2= ik12(2, ikp)

        ! calculate spectral function
        do i= 1, ndim
           dos_l(ik1, ik2)=dos_l(ik1, ik2)- aimag(GLL(i,i))
           dos_r(ik1, ik2)=dos_r(ik1, ik2)- aimag(GRR(i,i))
        enddo

     enddo

     !> we don't have to do allreduce operation
     call mpi_allreduce(dos_l, dos_l_mpi, size(dos_l),mpi_double_precision,&
                     mpi_sum, mpi_comm_world, ierr)
     call mpi_allreduce(dos_r, dos_r_mpi, size(dos_r),mpi_double_precision,&
                     mpi_sum, mpi_comm_world, ierr)


     !> calculate jdos
     do iq= 1+ cpuid, Nk1*Nk2, num_cpu
        iq1= ik12(1, iq)- Nk1_half
        iq2= ik12(2, iq)- Nk2_half
        if (cpuid==0) write(stdout, *) 'Jdos, iq, knv2', iq, nk1*nk2
        imin1= max(-Nk1_half-iq1, -Nk1_half)+ Nk1_half+ 1
        imax1= min(Nk1_half-iq1, Nk1_half)+ Nk1_half+ 1
        imin2= max(-Nk2_half-iq2, -Nk2_half)+ Nk2_half+ 1
        imax2= min(Nk2_half-iq2, Nk2_half)+ Nk2_half+ 1
        do ik1= imin1, imax1
           do ik2= imin2, imax2
              jdos_l(iq)= jdos_l(iq)+ dos_l_mpi(ik1, ik2)* dos_l_mpi(ik1+iq1, ik2+iq2)
              jdos_r(iq)= jdos_r(iq)+ dos_r_mpi(ik1, ik2)* dos_r_mpi(ik1+iq1, ik2+iq2)
           enddo !ik
        enddo !ik
     enddo !ik

     jdos_l_mpi=eps9
     jdos_r_mpi=eps9
     call mpi_reduce(jdos_l, jdos_l_mpi, size(jdos_l),mpi_double_precision,&
                     mpi_sum, 0, mpi_comm_world, ierr)
     call mpi_reduce(jdos_r, jdos_r_mpi, size(jdos_r),mpi_double_precision,&
                     mpi_sum, 0, mpi_comm_world, ierr)


     if (cpuid.eq.0)then
        open (unit=12, file='arc.dat_l')
        open (unit=13, file='arc.dat_r')
        do ikp=1, nk1*nk2
           write(12, '(3f16.8)')k12_shape(:, ikp), log(dos_l_mpi(ik12(1, ikp), ik12(2, ikp)))
           if (mod(ikp, nk2)==0) write(12, *)' '
           write(13, '(3f16.8)')k12_shape(:, ikp), log(dos_r_mpi(ik12(1, ikp), ik12(2, ikp)))
           if (mod(ikp, nk2)==0) write(13, *)' '
        enddo
        close(12)
        close(13)
        write(stdout,*)'ndim',ndim
        write(stdout,*) 'Nk1,Nk2,eta',Nk1, Nk2, eta
        write(stdout,*)'calculate density of state successfully'    
     endif

     if (cpuid.eq.0)then
        open (unit=1212, file='arc.jdat_l')
        open (unit=1313, file='arc.jdat_r')
        do ikp=1, nk1*nk2
           write(1212, '(30f16.8)')k12_shape(:, ikp), log(jdos_l_mpi(ikp))
           if (mod(ikp, nk2)==0) write(1212, *)' '
           write(1313, '(30f16.8)')k12_shape(:, ikp), log(jdos_r_mpi(ikp))
           if (mod(ikp, nk2)==0) write(1313, *)' '
        enddo
        close(1212)
        close(1313)
        write(stdout,*)'calculate joint density of state successfully'    
     endif

     !> write script for gnuplot
     if (cpuid==0) then
        open(unit=114, file='arc_l.gnu')
        write(114, '(a)')"set encoding iso_8859_1"
        write(114, '(a)')'#set terminal  postscript enhanced color'
        write(114, '(a)')"#set output 'arc_l.eps'"
        write(114, '(3a)')'#set terminal pngcairo truecolor enhanced', &
           '  font ",60" size 1920, 1680'
        write(114, '(3a)')'set terminal png truecolor enhanced', &
           ' font ",60" size 1920, 1680'
        write(114, '(a)')"set output 'arc_l.png'"
        write(114,'(2a)') 'set palette defined ( -10 "#194eff", ', &
           '0 "white", 10 "red" )'
        write(114, '(a)')'#set palette rgbformulae 33,13,10'
        write(114, '(a)')'unset ztics'
        write(114, '(a)')'unset key'
        write(114, '(a)')'set pm3d'
        write(114, '(a)')'set border lw 6'
        write(114, '(a)')'set size ratio -1'
        write(114, '(a)')'set view map'
        write(114, '(a)')'set xtics'
        write(114, '(a)')'set ytics'
        write(114, '(a)')'set xlabel "K_1 (1/{\305})"'
        write(114, '(a)')'set ylabel "K_2 (1/{\305})"'
        write(114, '(a)')'set colorbox'
        write(114, '(a, f8.5, a, f8.5, a)')'set xrange [', k1min_shape, ':', k1max_shape, ']'
        write(114, '(a, f8.5, a, f8.5, a)')'set yrange [', k2min_shape, ':', k2max_shape, ']'
        write(114, '(a)')'set pm3d interpolate 2,2'
        write(114, '(2a)')"splot 'arc.dat_l' u 1:2:3 w pm3d"

        close(114)
     endif

     !> write script for gnuplot
     if (cpuid==0) then
        open(unit=115, file='arc_r.gnu')
        write(115, '(a)')"set encoding iso_8859_1"
        write(115, '(a)')'#set terminal  postscript enhanced color'
        write(115, '(a)')"#set output 'arc_r.eps'"
        write(115, '(3a)')'#set terminal pngcairo truecolor enhanced', &
           '  font ",60" size 1920, 1680'
        write(115, '(3a)')'set terminal png truecolor enhanced', &
           ' font ",60" size 1920, 1680'
        write(115, '(a)')"set output 'arc_r.png'"
        write(115,'(2a)') 'set palette defined ( -10 "#194eff", ', &
           '0 "white", 10 "red" )'
        write(115, '(a)')'#set palette rgbformulae 33,13,10'
        write(115, '(a)')'unset ztics'
        write(115, '(a)')'unset key'
        write(115, '(a)')'set pm3d'
        write(115, '(a)')'set border lw 6'
        write(115, '(a)')'set size ratio -1'
        write(115, '(a)')'set view map'
        write(115, '(a)')'set xtics'
        write(115, '(a)')'set ytics'
        write(115, '(a)')'set xlabel "K_1 (1/{\305})"'
        write(115, '(a)')'set ylabel "K_2 (1/{\305})"'
        write(115, '(a)')'set colorbox'
        write(115, '(a, f8.5, a, f8.5, a)')'set xrange [', k1min_shape, ':', k1max_shape, ']'
        write(115, '(a, f8.5, a, f8.5, a)')'set yrange [', k2min_shape, ':', k2max_shape, ']'
        write(115, '(a)')'set pm3d interpolate 2,2'
        write(115, '(2a)')"splot 'arc.dat_r' u 1:2:3 w pm3d"

        close(115)
     endif


  return   
  end subroutine fermiarc_jdos
