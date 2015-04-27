! calculate bulk's energy band using wannier TB method
  subroutine ek_bulk2D

     use mpi
     use para

     implicit none

     integer :: ik, i, j
     integer :: nkx, nky
	  integer :: knv3
     integer :: ierr
     real(Dp) :: k(3)
     real(dp) :: kx
     real(dp) :: ky
     real(dp) :: kz
     real(Dp) :: W(Num_wann)
     real(dp) :: kxmin, kxmax, kymin, kymax
     real(dp), allocatable :: kxy(:,:)
     
     ! Hamiltonian of bulk system
     complex(Dp) :: Hamk_bulk(Num_wann,Num_wann) 

     ! eigen value of H
	  real(dp), allocatable :: gap(:)
	  real(dp), allocatable :: gap_mpi(:)
	  real(dp), allocatable :: eigv(:,:)
	  real(dp), allocatable :: eigv_mpi(:,:)

     nkx= Nk
     nky= Nk
     knv3= nkx*nky
     allocate( kxy(2, knv3))
     allocate( gap     (   knv3))
     allocate( gap_mpi (   knv3))
     allocate( eigv    (4, knv3))
     allocate( eigv_mpi(4, knv3))
     kxy = 0d0
     gap    = 0d0
     gap_mpi= 0d0
	  eigv    = 0d0
	  eigv_mpi= 0d0

     kxmin= 0.150d0/1d0
     kxmax= 0.170d0/1d0
     kymin= 0.06d0/1d0
     kymax= 0.08d0/1d0
     kz= 0.255
     ik =0
     do i= 1, nkx
     do j= 1, nky
        ik =ik +1
        kxy(1, ik)=kxmin+ (i-1)*(kxmax-kxmin)/dble(nkx-1)
        kxy(2, ik)=kymin+ (j-1)*(kymax-kymin)/dble(nky-1)
     enddo
     enddo

     do ik= 1+cpuid, knv3, num_cpu
	     if (cpuid==0) print * , ik

        k(1:2) = kxy(:, ik)
        k(3)= kz

        ! calculation bulk hamiltonian
        Hamk_bulk= 0d0
        call ham_bulk(k, Hamk_bulk)

        !> diagonalization by call zheev in lapack
        W= 0d0
        call eigensystem_c( 'N', 'U', Num_wann ,Hamk_bulk, W)

        !> for 88 bands
       !eigv(:, ik)= W(55:58)
       !gap (   ik)= W(57)- W(56)

        !> for 44 bands
        eigv(:, ik)= W(27:30)
        gap (   ik)= W(29)- W(28)
     enddo

     call mpi_allreduce(eigv,eigv_mpi,size(eigv),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
     call mpi_allreduce(gap ,gap_mpi,size(gap ),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)

     if (cpuid==0)then
        open(unit=14, file='bulkek2D.dat')
        do ik=1, knv3
           write(14, '(1000f19.9)')kxy(:, ik), eigv_mpi(:, ik)
        enddo
        close(14)
        open(unit=15, file='gap2D.dat')
        do ik=1, knv3
           write(15, '(1000f19.9)')kxy(:, ik), gap_mpi(ik)
        enddo
        close(15)
     endif

   return
   end subroutine ek_bulk2D
