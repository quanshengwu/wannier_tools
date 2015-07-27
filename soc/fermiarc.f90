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
     
     real(dp), allocatable :: dos(:)
     real(dp), allocatable :: dos_mpi(:)

     complex(dp), allocatable :: ones(:,:)
     complex(dp), allocatable :: GLL(:,:), GRR(:,:)
     complex(dp), allocatable :: H00(:,:), H01(:,:)

     nk1 = Nk 
     nk2 = Nk

     allocate( k12(2, nk1*nk2))
     allocate( k12_shape(2, nk1*nk2))
     allocate( dos(nk1*nk2))
     allocate( dos_mpi(nk1*nk2))
     allocate( GLL(ndim,ndim), GRR(ndim,ndim))
     k12=0d0
     k12_shape=0d0
     dos=0d0
     dos_mpi=0d0

     k1min= 0.35d0/1d0
     k1max= 0.65d0/1d0
     k2min=-0.15d0/1d0
     k2max= 0.15d0/1d0
     ikp=0
     do i= 1, nk1
     do j= 1, nk2
        ikp=ikp+1
        k12(1, ikp)=k1min+ (i-1)*(k1max-k1min)/dble(nk1-1)
        k12(2, ikp)=k2min+ (j-1)*(k2max-k2min)/dble(nk2-1)
        k12_shape(:, ikp)= k12(1, ikp)* Ka2+ k12(2, ikp)* Kb2
     enddo
     enddo

     k1= k1min*Ka2+ k2min*Kb2
     k2= k1max*ka2+ k2max*kb2
     k3= k1min*Ka2+ k2max*Kb2
     k4= k1max*ka2+ k2min*kb2
     k1min_shape= min(k1(1), k2(1), k3(1), k4(1))
     k2min_shape= min(k1(2), k2(2), k3(2), k4(2))
     k1max_shape= max(k1(1), k2(1), k3(1), k4(1))
     k2max_shape= max(k1(2), k2(2), k3(2), k4(2))

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
           dos(ikp)=dos(ikp)- aimag(GLL(i,i))
        enddo
     enddo

     !> we don't have to do allreduce operation
     call mpi_reduce(dos, dos_mpi, size(dos),mpi_double_precision,&
                     mpi_sum, 0, mpi_comm_world, ierr)

     if (cpuid.eq.0)then
        open (unit=12, file='arc.dat_l')
        do ikp=1, nk1*nk2
           write(12, '(3f16.8)')k12_shape(:, ikp), log(dos_mpi(ikp))
           if (mod(ikp, nk2)==0) write(12, *)' '
        enddo
        close(12)
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
           ' font ",60" size 3680, 3360'
        write(101, '(a)')"set output 'arc_l.png'"
        write(101,'(2a)') 'set palette defined ( -10 "white", ', &
           '0 "yellow", 10 "red" )'
        write(101, '(a)')'#set palette rgbformulae 33,13,10'
        write(101, '(a)')'unset ztics'
        write(101, '(a)')'unset key'
        write(101, '(a)')'set pm3d'
        write(101, '(a)')'set border lw 6'
        write(101, '(a)')'set size ratio -1'
        write(101, '(a)')'set view map'
        write(101, '(a)')'set xtics'
        write(101, '(a)')'set ytics'
        write(101, '(a)')'set xlabel "Kx"'
        write(101, '(a)')'set ylabel "Ky"'
        write(101, '(a)')'set colorbox'
        write(101, '(a, f8.5, a, f8.5, a)')'set xrange [', k1min_shape, ':', k1max_shape, ']'
        write(101, '(a, f8.5, a, f8.5, a)')'set yrange [', k2min_shape, ':', k2max_shape, ']'
        write(101, '(a)')'set pm3d interpolate 2,2'
        write(101, '(2a)')"splot 'arc.dat_l' u 1:2:3 w pm3d"

     endif

  return   
  end subroutine
