  !> calculate landau levels in 3D system in special k line
  !> fix B field
  !> construct on Dec 8 2015
  !> By QuanSheng Wu at ETH Zurich Honggerberg
  subroutine landau_level_k
     use para
     implicit none


     !> magnetic supercell size, perpendicular to the magnetic field
     integer :: Nq
     integer :: iq

     !> Ndim= Nq* Num_wann
     integer :: Ndim

     integer :: knv3
     integer :: NN
     integer :: nlines

     real(dp) :: Bmag
     ! wave vector 
     real(dp) :: k (3)
     real(dp) :: k1(2)
     real(dp) :: k2(2)
     real(Dp) :: k(2), kstart(2), kend(2)
     real(dp) :: kp(16,2), ke(16,2), kpath_stop(16)
     character(4) :: kpath_name(17)

     !> dim= Ndim, knv3
     real(dp), allocatable :: W(:)
     real(dp), allocatable :: eigv(:, :)
     real(dp), allocatable :: eigv_mpi(:, :)

     !> dim= Ndim*Ndim
     complex(dp), allocatable :: ham_landau(:, :)


     Nq= 100
     NN=Nk
     nlines= 4
     knv3=NN*nlines
     Ndim= Num_wann* Nq
     allocate( k_len (knv3))
     allocate( kpoint(knv3, 3))
     allocate( ham_landau(Ndim, Ndim))
     allocate( W( Ndim))
     allocate( eigv( Ndim, knv3))
     allocate( eigv_mpi( Ndim, knv3))
     eigv_mpi= 0d0
     eigv    = 0d0
     kpoint  = 0d0
     ham_landau= 0d0

     kp = 0d0
     ke = 0d0
     kp(1,:)=(/0.5d0, 0.0d0/)  ; kpath_name(1)= 'X'
     ke(1,:)=(/0.0d0, 0.0d0/)  
     kp(2,:)=(/0.0d0, 0.0d0/)  ; kpath_name(2)= 'G'
     ke(2,:)=(/0.0d0, 0.50d0/)  ! K
     kp(3,:)=(/0.0d0, 0.50d0/) ; kpath_name(3)= 'Y'     
     ke(3,:)=(/0.5d0, 0.5d0/)  ! K
     kp(4,:)=(/0.5d0, 0.5d0/)  ; kpath_name(4)= 'M'     
     ke(4,:)=(/0.0d0, 0.0d0/)  ; kpath_name(5)= 'G'  


     t1=0d0
     k_len=0d0
     kpath_stop= 0d0
     do j=1, nlines 
        do i=1, NN
           k = kp(j,:)
           kstart=k
           k1= kstart(1)*Ka2+ kstart(2)*Kb2

           k = ke(j,:)
           kend=k
           k2= kend(1)*Ka2+ kend(2)*Kb2

           kpoint(i+(j-1)*NN,:)= kstart+ (kend-kstart)*dble(i-1)/dble(NN-1)
           
           temp= dsqrt((k2(1)- k1(1))**2 &
                 +(k2(2)- k1(2))**2)/dble(NN-1) 

           if (i.gt.1) then
              t1=t1+temp
           endif
           k_len(i+(j-1)*NN)= t1
        enddo
        kpath_stop(j+1)= t1
     enddo


     !> calculate the landau levels along special k line
     do ik=1+ cpuid, knv3, num_cpu

        k = kpoint(ik, :)

        call ham_3Dlandau(Ndim, Nq, Bmag, k, ham_landau)
        
        !> diagonalization by call zheev in lapack
        W= 0d0
        call eigensystem_c( 'N', 'U', Ndim ,ham_landau, W)

        eigv(:, ik)= W
     enddo !ik

     call mpi_allreduce(eigv,eigv_mpi,size(eigv),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)


     return
  end subroutine landau_level_k

  !> calculate landau levels in 3D system for different B
  !> fix k point, usually Gamma point
  subroutine landau_level_B
     implicit none


  end subroutine landau_level_B

  !> calculate hamiltonian for landau levels
  subroutine ham_3Dlandau(Ndim, Nq, Bmag, k, ham_landau)
     use para
     implicit none

     integer, intent(in) :: Ndim
     integer, intent(in) :: Nq
     real(dp), intent(in) :: bmag
     real(dp), intent(in) :: k(3)
     complex(dp), intent(out) :: ham_landau(Ndim, Ndim)


     return
  end subroutine ham_3Dlandau

