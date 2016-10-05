! calculate bulk's energy band using wannier TB method
! Copyright (c) 2010 QuanSheng Wu. All rights reserved.
  subroutine ek_bulk2D

     use wmpi
     use para

     implicit none

     integer :: ik, i, j, ib
     integer :: nkx, nky
	  integer :: knv3
     integer :: ierr
     integer :: nwann
     real(Dp) :: k(3)
     real(dp) :: sx
     real(dp) :: sy
     real(dp) :: sz
     real(Dp) :: W(Num_wann)
     real(dp) :: kxmin, kxmax, kymin, kymax
     real(dp), allocatable :: kxy(:,:)
     real(dp), allocatable :: kxy_shape(:,:)
     
     ! Hamiltonian of bulk system
     complex(Dp) :: Hamk_bulk(Num_wann,Num_wann) 

     ! eigen value of H
	  real(dp), allocatable :: gap(:)
	  real(dp), allocatable :: gap_mpi(:)
	  real(dp), allocatable :: eigv(:,:)
	  real(dp), allocatable :: eigv_mpi(:,:)
	  real(dp), allocatable :: spin(:, :, :)
	  real(dp), allocatable :: spin_mpi(:, :, :)


     complex(dp), allocatable :: sigmax(:, :)
     complex(dp), allocatable :: sigmay(:, :)
     complex(dp), allocatable :: sigmaz(:, :)
     complex(dp), allocatable :: psi(:)

     nkx= Nk
     nky= Nk
     knv3= nkx*nky
     allocate( psi( Num_wann))
     allocate( sigmax( Num_wann, Num_wann))
     allocate( sigmay( Num_wann, Num_wann))
     allocate( sigmaz( Num_wann, Num_wann))
     allocate( spin    (3, Num_wann, knv3))
     allocate( spin_mpi(3, Num_wann, knv3))
     spin = 0d0
     spin_mpi = 0d0
     sigmax= 0d0
     sigmay= 0d0
     sigmaz= 0d0
 
     allocate( kxy(3, knv3))
     allocate( kxy_shape(3, knv3))
     allocate( gap     (   knv3))
     allocate( gap_mpi (   knv3))
     allocate( eigv    (4, knv3))
     allocate( eigv_mpi(4, knv3))
     kxy = 0d0
     gap    = 0d0
     gap_mpi= 0d0
	  eigv    = 0d0
	  eigv_mpi= 0d0

     nwann= Num_wann
     if (soc>0) nwann= Num_wann/2
     print *, 'nwann', nwann

     if (SOC==0) stop 'you should set soc=0 in the input file'
     !> construct spin matrix
     do i= 1, nwann
        sigmax(i, i+ nwann)= 1d0
        sigmax(i+ nwann, i)= 1d0
        sigmay(i, i+ nwann)= -zi
        sigmay(i+ nwann, i)= zi
        sigmaz(i, i)= 1d0
        sigmaz(i+ nwann, i+ nwann)= -1d0
     enddo

     ik =0
     do i= 1, nkx
     do j= 1, nky
        ik =ik +1
        kxy(:, ik)= K3D_start+ K3D_vec1*(i-1)/dble(nkx-1)+ K3D_vec2*(j-1)/dble(nky-1)
        kxy_shape(:, ik)= kxy(1, ik)* Kua+ kxy(2, ik)* Kub+ kxy(3, ik)* Kuc 
     enddo
     enddo

     do ik= 1+cpuid, knv3, num_cpu
        if (cpuid==0) print * , ik, knv3

        k = kxy(:, ik)

        ! calculation bulk hamiltonian
        Hamk_bulk= 0d0
        call ham_bulk(k, Hamk_bulk)

        !> diagonalization by call zheev in lapack
        W= 0d0
        call eigensystem_c( 'V', 'U', Num_wann ,Hamk_bulk, W)

        do ib= 1, Num_wann
           psi(:)= Hamk_bulk(:, ib)
          !print *, sum(conjg(psi(:))*psi(:))
           sx= 0d0
           sy= 0d0
           sz= 0d0
           do i= 1, Num_wann
              do j= 1, Num_wann
                 sx= sx+ conjg(psi(i))* sigmax(i, j)* psi(j)
                 sy= sy+ conjg(psi(i))* sigmay(i, j)* psi(j)
                 sz= sz+ conjg(psi(i))* sigmaz(i, j)* psi(j)
              enddo ! j
           enddo ! i
           spin(1, ib, ik)= sx
           spin(2, ib, ik)= sy
           spin(3, ib, ik)= sz
          !print *, W(ib),sx*sx+sy*sy+sz*sz
        enddo ! ib

        !> for 44 bands
        eigv(:, ik)= W(Numoccupied-1:Numoccupied+ 2)
        gap (   ik)= W(Numoccupied+1)- W(Numoccupied)
     enddo

#if defined (MPI)
     call mpi_allreduce(eigv,eigv_mpi,size(eigv),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
     call mpi_allreduce(gap ,gap_mpi,size(gap ),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
     call mpi_allreduce(spin ,spin_mpi,size(spin ),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
     eigv_mpi= eigv
     gap_mpi= gap
     spin_mpi= spin
#endif


     if (cpuid==0)then
        open(unit=14, file='bulkek2D.dat')
        do ik=1, knv3
           write(14, '(1000f19.9)')kxy_shape(:, ik), eigv_mpi(:, ik), (spin_mpi(:, ib, ik), ib=Numoccupied-1, Numoccupied+2)
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


   !> using green's function
  subroutine ek_bulk2D_spin

     use wmpi
     use para

     implicit none

     integer :: ik, i, j
     integer :: nkx, nky
	  integer :: nwann
	  integer :: knv3
     integer :: ierr
     real(Dp) :: k(3)
     real(Dp) :: k11(3), k12(3)
     real(Dp) :: k21(3), k22(3)
     real(Dp) :: W(Num_wann)
     real(dp) :: kxmin, kxmax, kymin, kymax
     real(dp), allocatable :: kxy(:,:)
     real(dp), allocatable :: kxy_shape(:,:)
     
     ! Hamiltonian of bulk system
     complex(Dp) :: Hamk_bulk(Num_wann,Num_wann) 

     ! eigen value of H
	  real(dp), allocatable :: gap(:)
	  real(dp), allocatable :: gap_mpi(:)
	  real(dp), allocatable :: eigv(:,:)
	  real(dp), allocatable :: eigv_mpi(:,:)
     real(dp), allocatable :: dos(:)
     real(dp), allocatable :: dos_mpi(:)
	  real(dp), allocatable :: spin(:, :)
	  real(dp), allocatable :: spin_mpi(:, :)


     complex(dp), allocatable :: sigmax(:, :)
     complex(dp), allocatable :: sigmay(:, :)
     complex(dp), allocatable :: sigmaz(:, :)
     complex(dp), allocatable :: psi(:)
     complex(dp), allocatable :: mat1(:, :)

     nkx= Nk
     nky= Nk
     knv3= nkx*nky
     allocate( psi( Num_wann))
     allocate( sigmax( Num_wann, Num_wann))
     allocate( sigmay( Num_wann, Num_wann))
     allocate( sigmaz( Num_wann, Num_wann))
     allocate( mat1  ( Num_wann, Num_wann))
     allocate( dos     (knv3))
     allocate( dos_mpi (knv3))
     allocate( spin    (3, knv3))
     allocate( spin_mpi(3, knv3))
 
     allocate( kxy(3, knv3))
     allocate( kxy_shape(3, knv3))
     allocate( gap     (   knv3))
     allocate( gap_mpi (   knv3))
     allocate( eigv    (4, knv3))
     allocate( eigv_mpi(4, knv3))
     kxy     = 0d0
     gap     = 0d0
     gap_mpi = 0d0
	  eigv    = 0d0
	  eigv_mpi= 0d0
     spin    = 0d0
     spin_mpi= 0d0
     dos    = 0d0
     dos_mpi= 0d0
     sigmax= 0d0
     sigmay= 0d0
     sigmaz= 0d0

     nwann= Num_wann/2

     if (SOC==0) stop 'you should set soc=0 in the input file'
     !> construct spin matrix
     do i= 1, nwann
        sigmax(i, i+ nwann)= 1d0
        sigmax(i+ nwann, i)= 1d0
        sigmay(i, i+ nwann)= -zi
        sigmay(i+ nwann, i)= zi
        sigmaz(i, i)= 1d0
        sigmaz(i+ nwann, i+ nwann)= -1d0
     enddo


     k11=(/-0.0d0,  0.0d0, -0.0d0/) ! 
     k12=(/ 1.0d0,  0.0d0,  1.0d0/) ! X
     k21=(/ 0.0d0,  0.0d0,  0.0d0/) ! 
     k22=(/ 1.0d0,  1.0d0,  0.0d0/) ! Z


     ik =0
     do i= 1, nkx
     do j= 1, nky
        ik =ik +1
        kxy(:, ik)= k11+(k12-k11)*(i-1)/dble(nkx-1)+  k21+ (k22-k21)*(j-1)/dble(nky-1)
        kxy_shape(:, ik)= kxy(1, ik)* Kua+ kxy(2, ik)* Kub+ kxy(3, ik)* Kuc 
     enddo
     enddo

     do ik= 1+cpuid, nkx*nky, num_cpu
	     if (cpuid==0) print * , ik

        k = kxy(:, ik)

        ! calculation bulk hamiltonian
        Hamk_bulk= 0d0
        call ham_bulk(k, Hamk_bulk)

        !> diagonalization by call zheev in lapack
        Hamk_bulk= Hamk_bulk+ zi*eta_arc
        call inv(Num_wann, Hamk_bulk)
        do i=1, Num_wann
           dos(ik)= dos(ik) -aimag(Hamk_bulk(i, i))
        enddo

        call mat_mul(Num_wann, Hamk_bulk, sigmax, mat1)

        do i=1, Num_wann
           spin(1, ik)= spin(1, ik) -aimag(mat1(i, i))
        enddo

        call mat_mul(Num_wann, Hamk_bulk, sigmay, mat1)
        do i=1, Num_wann
           spin(2, ik)= spin(2, ik) -aimag(mat1(i, i))
        enddo

        call mat_mul(Num_wann, Hamk_bulk, sigmaz, mat1)
        do i=1, Num_wann
           spin(3, ik)= spin(3, ik) -aimag(mat1(i, i))
        enddo

        !> for 44 bands
        eigv(:, ik)= W(Numoccupied-1:Numoccupied+ 2)
        gap (   ik)= W(Numoccupied+1)- W(Numoccupied)
     enddo !ik

#if defined (MPI)
     call mpi_allreduce(eigv,eigv_mpi,size(eigv),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
     call mpi_allreduce(gap ,gap_mpi,size(gap ),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
     call mpi_allreduce(spin ,spin_mpi,size(spin ),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
     call mpi_allreduce(dos ,dos_mpi,size(dos ),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
     gap_mpi= gap
     dos_mpi= dos
     eigv_mpi= eigv
     spin_mpi= spin
#endif



     if (cpuid==0)then
        open(unit=14, file='bulkek2D.dat')
        do ik=1, knv3
           write(14, '(1000f20.6)')kxy_shape(:, ik), eigv_mpi(:, ik) 
        enddo
        close(14)

        open(unit=15, file='gap2D.dat')
        do ik=1, knv3
           write(15, '(1000f20.6)')kxy_shape(:, ik), gap_mpi(ik)
        enddo
        close(15)

        open(unit=16, file='bulkspin2D.dat')
        ik= 0
        do i= 1, nkx
           do j= 1, nky
              ik= ik+ 1

              if (ik==1 .or. ik==knv3.or. dos_mpi(ik)<4d0) cycle
              if ((dos_mpi(ik)> dos_mpi(ik-1)) .and. (dos_mpi(ik)> dos_mpi(ik+1))) then
                !print *, ik, dos_mpi(ik)
                 write(16, '(10000f20.5)')kxy_shape(:, ik), dos_mpi(ik) , spin_mpi(:, ik) 

              endif
           enddo
           write(16, '(a)')' '
        enddo
        close(16)

     endif

   return
   end subroutine ek_bulk2D_spin
