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

     integer :: nkx
     integer :: nkz

! kpoint loop index
     integer :: ikp

     real(dp) :: k1(2)
     real(dp) :: k2(2)
     real(dp) :: k(2)

     real(dp) :: kxmin, kxmax, kzmin, kzmax, omega

     real(dp), allocatable :: kxz(:,:)
     real(dp), allocatable :: kxz_shape(:,:)
     
     real(dp), allocatable :: dos(:)
     real(dp), allocatable :: dos_mpi(:)

     complex(dp), allocatable :: ones(:,:)
     complex(dp), allocatable :: GLL(:,:), GRR(:,:)
     complex(dp), allocatable :: H00(:,:), H01(:,:)

     nkx = Nk 
     nkz = Nk

     allocate( kxz(2, nkx*nkz))
     allocate( kxz_shape(2, nkx*nkz))
     allocate( dos(nkx*nkz))
     allocate( dos_mpi(nkx*nkz))
     allocate( GLL(ndim,ndim), GRR(ndim,ndim))
     kxz=0d0
     kxz_shape=0d0
     dos=0d0
     dos_mpi=0d0

     kzmin=-0.15d0/1d0
     kzmax= 0.15d0/1d0
     kxmin=-0.15d0/1d0
     kxmax= 0.15d0/1d0
     ikp=0
     do i= 1, nkx
     do j= 1, nkz
        ikp=ikp+1
        kxz(1, ikp)=kxmin+ (i-1)*(kxmax-kxmin)/dble(nkx-1)
        kxz(2, ikp)=kzmin+ (j-1)*(kzmax-kzmin)/dble(nkz-1)
        kxz_shape(:, ikp)= kxz(1, ikp)* Ka2+ kxz(2, ikp)* Kb2
     enddo
     enddo

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


     omega = omegamin 

     do ikp= 1+cpuid, nkx*nkz, num_cpu
        if (cpuid==0) print *, 'Arc', ikp, nkx*nkz
        k(1)= kxz(1, ikp)
        k(2)= kxz(2, ikp)

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
        do ikp=1, nkx*nkz
           write(12, '(3f16.8)')kxz_shape(:, ikp), log(dos_mpi(ikp))
           if (mod(ikp, nkz)==0) write(12, *)' '
        enddo
        close(12)
        write(*,*)'ndim',ndim
        write(*,*) 'Nkx,Nkz,eta',Nkx, Nkz, eta
        write(*,*)'calculate density of state successfully'    
     endif

     !> write script for gnuplot
     if (cpuid==0) then
        open(unit=101, file='arc_l.gnu')
        write(101, '(a)') 'set terminal  postscript enhanced color'
        write(101,'(2a)') '#set palette defined ( -10 "green", ', &
           '0 "yellow", 10 "red" )'
        write(101, '(a)')'set palette rgbformulae 33,13,10'
        write(101, '(a)')"set output 'arc_l.eps'"
        write(101, '(a)')'unset ztics'
        write(101, '(a)')'unset key'
        write(101, '(a)')'set pm3d'
        write(101, '(a)')'set view map'
        write(101, '(a)')'set xtics font ",24"'
        write(101, '(a)')'set ytics font ",24"'
        write(101, '(a, f8.5, a, f8.5, a)')'set xrange [', kxmin, ':', kxmax, ']'
        write(101, '(a, f8.5, a, f8.5, a)')'set yrange [', kzmin, ':', kzmax, ']'
        write(101, '(a)')'set pm3d interpolate 2,2'
        write(101, '(2a)')"splot 'arc.dat_l' u 1:2:3 w pm3d"

     endif

  return   
  end subroutine
