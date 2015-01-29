! calculate bulk's energy band using wannier TB method
  subroutine ek_bulk

     use para

     implicit none

     integer :: ik, i, j
	  integer :: knv3
     integer :: ierr
     real(Dp) :: k(3)
     real(Dp) :: W(Num_wann)
     
     ! Hamiltonian of bulk system
     complex(Dp) :: Hamk_bulk(Num_wann,Num_wann) 

     ! eigen value of H
	  real(dp), allocatable :: eigv(:,:)
	  real(dp), allocatable :: eigv_mpi(:,:)

     knv3= nk3_band
     allocate( eigv    (Num_wann, knv3))
     allocate( eigv_mpi(Num_wann, knv3))
	  eigv    = 0d0
	  eigv_mpi= 0d0

     do ik= 1+cpuid, knv3, num_cpu
	     if (cpuid==0) print * , ik

        k = k3points(:, ik)

        ! calculation bulk hamiltonian
        Hamk_bulk= 0d0
        call ham_bulk(k, Hamk_bulk)

        W= 0d0
        call eigensystem_c( 'N', 'U', Num_wann ,Hamk_bulk, W)

        eigv(:, ik)= W

     enddo

     call mpi_allreduce(eigv,eigv_mpi,size(eigv),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)

     if (cpuid==0)then
        open(unit=14, file='ek-tb.dat')
   
        do ik=1, knv3
           write(14, '(1000f19.9)')k3len(ik),eigv_mpi(:, ik)
        enddo
   
        close(14)
     endif

   return
   end subroutine ek_bulk
