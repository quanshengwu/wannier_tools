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


     use mpi
     use para
     implicit none
     


     integer :: ierr

! general loop index
     integer :: i,j 

     integer :: nk1
     integer :: nk2

! kpoint loop index
     integer :: ikp

     real(dp) :: k1(2)
     real(dp) :: k2(2)
     real(dp) :: k3(2)
     real(dp) :: k4(2)
     real(dp) :: k(2)

     real(dp) :: k1min, k1max, k2min, k2max, omega
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

     nk1 = Nk 
     nk2 = Nk

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
           k12(:, ikp)=K2D_start1+ (i-1)*(K2D_end1- K2D_start1)/dble(nk1-1) &
                      +K2D_start2+ (j-1)*(K2D_end2- K2D_start2)/dble(nk2-1)
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
        if (cpuid==0) write(*, *) 'Arc', ikp, nk1*nk2
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
        open(unit=101, file='arc_l.gnu')
        write(101, '(a)')'#set terminal  postscript enhanced color'
        write(101, '(a)')"#set output 'arc_l.eps'"
        write(101, '(3a)')'set terminal  pngcairo truecolor enhanced', &
           ' font ",60" size 1920, 1680'
        write(101, '(a)')"set output 'arc_l.png'"
        write(101,'(2a)') 'set palette defined ( -10 "#194eff", ', &
           '0 "white", 10 "red" )'
        write(101, '(a)')'#set palette rgbformulae 33,13,10'
        write(101, '(a)')'unset ztics'
        write(101, '(a)')'unset key'
        write(101, '(a)')'set pm3d'
        write(101, '(a)')'set border lw 6'
        write(101, '(a)')'set size ratio -1'
        write(101, '(a)')'set view map'
        write(101, '(a)')'set xtics'
        write(101, '(a)')'set ytics'
        write(101, '(a)')'set xlabel "K1"'
        write(101, '(a)')'set ylabel "K2"'
        write(101, '(a)')'set colorbox'
        write(101, '(a, f8.5, a, f8.5, a)')'set xrange [', k1min_shape, ':', k1max_shape, ']'
        write(101, '(a, f8.5, a, f8.5, a)')'set yrange [', k2min_shape, ':', k2max_shape, ']'
        write(101, '(a)')'set pm3d interpolate 2,2'
        write(101, '(2a)')"splot 'arc.dat_l' u 1:2:3 w pm3d"

        close(101)
     endif

     !> write script for gnuplot
     if (cpuid==0) then
        open(unit=101, file='arc_r.gnu')
        write(101, '(a)')'#set terminal  postscript enhanced color'
        write(101, '(a)')"#set output 'arc_r.eps'"
        write(101, '(3a)')'set terminal  pngcairo truecolor enhanced', &
           ' font ",60" size 1920, 1680'
        write(101, '(a)')"set output 'arc_r.png'"
        write(101,'(2a)') 'set palette defined ( -10 "#194eff", ', &
           '0 "white", 10 "red" )'
        write(101, '(a)')'#set palette rgbformulae 33,13,10'
        write(101, '(a)')'unset ztics'
        write(101, '(a)')'unset key'
        write(101, '(a)')'set pm3d'
        write(101, '(a)')'set border lw 6'
        write(101, '(a)')'set size ratio -1'
        write(101, '(a)')'set view map'
        write(101, '(a)')'set xtics'
        write(101, '(a)')'set ytics'
        write(101, '(a)')'set xlabel "K1"'
        write(101, '(a)')'set ylabel "K2"'
        write(101, '(a)')'set colorbox'
        write(101, '(a, f8.5, a, f8.5, a)')'set xrange [', k1min_shape, ':', k1max_shape, ']'
        write(101, '(a, f8.5, a, f8.5, a)')'set yrange [', k2min_shape, ':', k2max_shape, ']'
        write(101, '(a)')'set pm3d interpolate 2,2'
        write(101, '(2a)')"splot 'arc.dat_r' u 1:2:3 w pm3d"

        close(101)
     endif


  return   
  end subroutine
