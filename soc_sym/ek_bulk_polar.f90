! calculate bulk's energy band using wannier TB method
  subroutine ek_bulk_polar

     use mpi
     use para

     implicit none

     integer :: ik, i, j
     integer :: nkx, nky
	  integer :: knv3
     integer :: ierr
     integer :: nktheta
     integer :: nkr
     real(Dp) :: k(3)
     real(dp) :: kx
     real(dp) :: ky
     real(dp) :: kz
     real(dp) :: ktheta
     real(dp) :: kr
     real(Dp) :: W(Num_wann)
     real(dp) :: krmax
     real(dp), allocatable :: kxy(:,:)
     real(dp), allocatable :: krtheta(:,:)
     
     ! Hamiltonian of bulk system
     complex(Dp) :: Hamk_bulk(Num_wann,Num_wann) 

     ! eigen value of H
	  real(dp), allocatable :: eigv(:,:)
	  real(dp), allocatable :: eigv_mpi(:,:)

     nktheta= 50
     nkr= Nk
     knv3= nktheta*nkr
     allocate( kxy(2, knv3))
     allocate( krtheta(2, knv3))
     allocate( eigv    (4, knv3))
     allocate( eigv_mpi(4, knv3))
     kxy = 0d0
     krtheta = 0d0
	  eigv    = 0d0
	  eigv_mpi= 0d0

     krmax= 0.2d0
     kz= 0d0
     ik =0

     do i= 1, nktheta
        ktheta= 2d0*pi*(i-1)/nktheta
        do j= 1, nkr
           kr= krmax*(j-1)/nkr
           ik =ik +1
           kxy(1, ik)= kr*cos(ktheta)
           kxy(2, ik)= kr*sin(ktheta)
           krtheta(1, ik)= ktheta
           krtheta(2, ik)= kr
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

        eigv(:, ik)= W(55:58)

     enddo

     call mpi_allreduce(eigv,eigv_mpi,size(eigv),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)

     if (cpuid==0)then
        open(unit=14, file='bulkek_polar.dat')
        do ik=1, knv3
           write(14, '(1000f19.9)')krtheta(:, ik), eigv_mpi(:, ik)
        enddo
        close(14)
     endif

   return
   end subroutine ek_bulk_polar
