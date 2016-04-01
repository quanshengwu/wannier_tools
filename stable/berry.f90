!--------+--------+--------+--------+--------+--------+--------+------!
!> Subroutine for calculating Berry phase for a giving path
! Comments:
!          At present, you have to define the k path that you want 
!          in kpoints
! Author : QuanSheng Wu (wuquansheng@gmail.com)
! 31 Mar 2016
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

      real(dp) :: kx
      real(dp) :: ky
      real(dp) :: kz
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
      allocate(hamk(Num_wann, Num_wann))
      allocate(hamk_dag(Num_wann, Num_wann))
      allocate(Eigenvector(Num_wann, Num_wann, Nk))
      allocate(eigenvalue(Num_wann))
      allocate(U(nfill, nfill))
      allocate(Sigma(nfill, nfill))
      allocate(VT(nfill, nfill))
      allocate(phase(Num_wann))
      hamk=0d0
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

      !> define a k path
      do ik=1, Nk
         phi= (ik-1d0)*2d0*pi/dble(Nk)
         kx= 0.00d0+ cos(phi)*0.100d0!pi
         ky= 0.20d0+ sin(phi)*0.100d0!pi
         kz= 0.00d0
         kpoints(1, ik)= kx
         kpoints(2, ik)= ky
         kpoints(3, ik)= kz
      enddo

      !> for each ky, we can get wanniercenter
      do ik=1, Nk
         k(1)= kpoints(1, ik)
         k(2)= kpoints(2, ik)
         k(3)= kpoints(3, ik)

         call cart_direct(k, k1)
         write(*, '(a, 3f7.3, a, 3f7.3)')'k',k,'  k1', k1

         call ham_bulk(k1,hamk)

         !> diagonal hamk
         call eigensystem_c('V', 'U', Num_wann, hamk, eigenvalue)

         Eigenvector(:, :, ik)= hamk
      enddo

      !> sum over k to get berry phase
      phase= 1d0
      do ik= 1, Nk
         hamk= Eigenvector(:, :, ik)
         hamk_dag= conjg(transpose(hamk))
         if (ik==nk) then
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
            phase(i)= phase(i)*overlap
         enddo
      enddo  !< ik
      print *, 'berry phase'
      do i=1, Num_wann
         write(*, '(10f10.5)')aimag(log(phase(i)))/pi
      enddo

      return
   end subroutine berryphase
