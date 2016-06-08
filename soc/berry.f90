!--------+--------+--------+--------+--------+--------+--------+------!
!> Subroutine for calculating Berry phase for a giving path
! Comments:
!          At present, you have to define the k path that you want 
!          in kpoints
! Author : QuanSheng Wu (wuquansheng@gmail.com)
! 31 Mar 2016
! Copyright (c) 2010 QuanSheng Wu. All rights reserved.
!--------+--------+--------+--------+--------+--------+--------+------!

   subroutine  berryphase
      use para
      use mpi
      implicit none

      integer :: i
      integer :: j
      integer :: nfill

      integer :: ik

      !> k points in kx-ky plane
      real(dp), allocatable :: kpoints(:, :)

      !> hamiltonian for each k point
      !> and also the eigenvector of hamiltonian after eigensystem_c
      complex(dp), allocatable :: Hamk(:, :)
      complex(dp), allocatable :: Hamk1(:, :)
      complex(dp), allocatable :: Hamk_dag(:, :)

      !> eigenvector for each kx
      complex(dp), allocatable :: Eigenvector(:, :, :)

      !> Mmnkb=<u_n(k)|u_m(k+b)>
      !> |u_n(k)> is the periodic part of wave function
      complex(dp), allocatable :: Mmnkb(:, :)
      complex(dp), allocatable :: Mmnkb_com(:, :)
      complex(dp), allocatable :: Mmnkb_full(:, :)

      !> 
      complex(dp), allocatable :: Lambda_eig(:)
      complex(dp), allocatable :: Lambda(:, :)
      complex(dp), allocatable :: Lambda0(:, :)

      !> three matrix for SVD 
      !> M= U.Sigma.V^\dag
      !> VT= V^\dag
      complex(dp), allocatable :: U(:, :)
      real   (dp), allocatable :: Sigma(:, :)
      complex(dp), allocatable :: VT(:, :)
   
      complex(dp), allocatable :: phase(:)
      !> eigenvalue
      real(dp), allocatable :: eigenvalue(:)

      complex(dp), allocatable :: mat1(:, :)
      complex(dp), allocatable :: mat2(:, :)
      real(dp) :: kx
      real(dp) :: ky
      real(dp) :: kz
      real(dp) :: r1, r2
      real(dp) :: k(3)
      real(dp) :: k1(3)
      real(dp) :: phi
      complex(dp) :: overlap

      nfill= Numoccupied

      allocate(kpoints(3, Nk))
      kpoints= 0d0

      allocate(Lambda_eig(nfill))
      allocate(Lambda(nfill, nfill))
      allocate(Lambda0(nfill, nfill))
      allocate(Mmnkb(nfill, nfill))
      allocate(Mmnkb_com(nfill, nfill))
      allocate(Mmnkb_full(Num_wann, Num_wann))
      allocate(mat1(Num_wann, Num_wann))
      allocate(mat2(Num_wann, Num_wann))
      allocate(hamk(Num_wann, Num_wann))
      allocate(hamk1(Num_wann, Num_wann))
      allocate(hamk_dag(Num_wann, Num_wann))
      allocate(Eigenvector(Num_wann, Num_wann, Nk))
      allocate(eigenvalue(Num_wann))
      allocate(U(nfill, nfill))
      allocate(Sigma(nfill, nfill))
      allocate(VT(nfill, nfill))
      allocate(phase(Num_wann))
      hamk=0d0
      hamk1=0d0
      eigenvalue=0d0
      Eigenvector=0d0
      Mmnkb_full=0d0
      Mmnkb=0d0
      Mmnkb_com=0d0
      Lambda =0d0
      Lambda0=0d0
      U= 0d0
      Sigma= 0d0
      VT= 0d0

      r1=0.020d0
      r2=0.020d0

      !> define a circle k path
      do ik=1, Nk
         phi= (ik-1d0)*2d0*pi/dble(Nk)
         kx= 0.400000d0+ cos(phi)*r1
         ky= 0.330000d0
         kz= 0.0d0+ sin(phi)*r2
         kpoints(1, ik)= kx
         kpoints(2, ik)= ky
         kpoints(3, ik)= kz
      enddo

      !> for each ky, we can get wanniercenter
      do ik=1, Nk
         k(1)= kpoints(1, ik)
         k(2)= kpoints(2, ik)
         k(3)= kpoints(3, ik)

         k1= k
        !call cart_direct(k, k1)
        !if (cpuid.eq.0)write(stdout, '(a, 3f7.3, a, 3f7.3)')'k',k,'  k1', k1

        !call ham_bulk(k1,hamk1)

        !k = kpoints(:, ik)
        !k(3)= -k(3)
        !Hamk= 0d0
         call ham_bulk    (k, Hamk)
        
         !> symmetrization
        !call mat_mul(Num_wann, mirror_z, hamk, mat1)
        !call mat_mul(Num_wann, mat1, mirror_z, mat2)
        !hamk= (Hamk1+ mat2)/2.d0


         !> diagonal hamk
         call eigensystem_c('V', 'U', Num_wann, hamk, eigenvalue)

         Eigenvector(:, :, ik)= hamk
      enddo

      !> sum over k to get berry phase
      phase= 1d0
      do ik= 1, Nk
         hamk= Eigenvector(:, :, ik)
         hamk_dag= conjg(transpose(hamk))
         if (ik==Nk) then
            hamk= Eigenvector(:, :, 1)
         else
            hamk= Eigenvector(:, :, ik+ 1)
         endif

         !> <u_k|u_k+1>
         do i=1, Num_wann
            overlap= 0d0
            do j=1, Num_wann
               overlap= overlap+ hamk_dag(i, j)* hamk(j, i)
            enddo
           !phase(i)= phase(i)- aimag(log(overlap))/pi
            phase(i)= overlap*phase(i)
         enddo
      enddo  !< ik
      if (cpuid==0)write(stdout, *) 'berry phase'
      do i=1, Num_wann
        !if (cpuid==0)  write(*, '(10f10.5)')aimag(log(phase(i)))/pi
        !write(*, '(10f10.5)')phase(i)
      enddo
      if (cpuid==0) write(stdout, '(f18.6)')sum(aimag(log(phase(1:nfill)))/pi)
     !print *, sum(phase(1:nfill))

      return
   end subroutine berryphase
