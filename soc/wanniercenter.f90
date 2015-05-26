! this suboutine is used for wannier center calculation for slab system

   subroutine  wannier_center2D
      use para
      use mpi
      implicit none

      integer :: Nkx
      integer :: Nky
      integer :: knv2

      integer :: i
      integer :: j
      integer :: l 
      integer :: nfill

      integer :: ikx
      integer :: iky

      integer :: ierr

      !> k points in kx-ky plane
      complex(dp), allocatable :: kpoints(:, :, :)

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
   
      !> wannier centers for each ky, bands
      real(dp), allocatable :: WannierCenterKy(:, :)
      real(dp), allocatable :: WannierCenterKy_mpi(:, :)

      !> eigenvalue
      real(dp), allocatable :: eigenvalue(:)

      real(dp) :: kx
      real(dp) :: ky
      real(dp) :: k(2)

      Nkx= Nk
      Nky= 40
      knv2= Nkx*Nky

      nfill= 6*Nslab

      allocate(kpoints(2, Nkx, Nky))
      kpoints= 0d0

      allocate(Lambda_eig(nfill))
      allocate(Lambda(nfill, nfill))
      allocate(Lambda0(nfill, nfill))
      allocate(Mmnkb(nfill, nfill))
      allocate(Mmnkb_com(nfill, nfill))
      allocate(Mmnkb_full(Num_wann*Nslab, Num_wann*Nslab))
      allocate(hamk(Num_wann*Nslab, Num_wann*Nslab))
      allocate(hamk_dag(Num_wann*Nslab, Num_wann*Nslab))
      allocate(Eigenvector(Num_wann*Nslab, Num_wann*Nslab, Nkx))
      allocate(eigenvalue(Num_wann*Nslab))
      allocate(U(nfill, nfill))
      allocate(Sigma(nfill, nfill))
      allocate(VT(nfill, nfill))
      allocate(WannierCenterKy(nfill, Nky))
      allocate(WannierCenterKy_mpi(nfill, Nky))
      WannierCenterKy= 0d0
      WannierCenterKy_mpi= 0d0
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

      do iky=1, Nky
         do ikx=1, Nkx
            kx= (ikx-1)/real(Nkx)
            ky= (iky-1)/real(Nky)
            kpoints(1, ikx, iky)= kx
            kpoints(2, ikx, iky)= ky
         enddo
      enddo

      !> for each ky, we can get wanniercenter
      do iky=1+ cpuid, nky, num_cpu
         Lambda0=0d0
         do i=1, nfill
            Lambda0(i, i)= 1d0
         enddo

         if (cpuid==0) print *, iky, nky
         !> for each kx, we get the eigenvectors
         do ikx=1, nkx
            k(1)= kpoints(1, ikx, iky)
            k(2)= kpoints(2, ikx, iky)

            call ham_slab(k,hamk)

            !> diagonal hamk
            call eigensystem_c('V', 'U', Num_wann*Nslab, hamk, eigenvalue)

            Eigenvector(:, :, ikx)= hamk
         enddo

         !> sum over kx to get wanniercenters
         do ikx=1, nkx
            hamk= Eigenvector(:, :, ikx)
            hamk_dag= conjg(transpose(hamk))
            if (ikx==nkx) then
               hamk= Eigenvector(:, :, 1)
            else
               hamk= Eigenvector(:, :, ikx+ 1)
            endif

            !> <u_k|u_k+1>
            call mat_mul(Num_wann*Nslab, hamk_dag, hamk, Mmnkb_full)
            Mmnkb= Mmnkb_full(1:nfill, 1:nfill)
           !Mmnkb_com= 0d0
           !hamk_dag= Eigenvector(:, :, ikx)
           !hamk= Eigenvector(:, :, ikx+1)
           !do i=1, nfill
           !   do j=1, nfill
           !      do l= 1, Num_wann*Nslab
           !         Mmnkb_com(i, j)=  Mmnkb_com(i, j)+ conjg(hamk_dag(l, i))* hamk(l, j)
           !      enddo
           !   enddo
           !enddo

           !print *, maxval(real(Mmnkb-Mmnkb_com))
           !stop


            !> perform Singluar Value Decomposed of Mmnkb
            call zgesvd_pack(nfill, Mmnkb, U, Sigma, VT)

            !> after the calling of zgesvd_pack, Mmnkb becomes a temporal matrix
            U= conjg(transpose(U))
            VT= conjg(transpose(VT))
            call mat_mul(nfill, VT, U, Mmnkb)

            !> check hermicity
           !do i=1, nfill
           !   do j=i, nfill
           !      if (abs(Mmnkb(i, j)-conjg(Mmnkb(j, i)))>0.0001d0)then
           !         print *, 'Mmnkb is not Hermitian'
           !         print*, i, j, Mmnkb(i, j), Mmnkb(j, i)

           !      endif
           !   enddo
           !enddo

           !stop


            call mat_mul(nfill, Mmnkb, Lambda0, Lambda)
            Lambda0 = Lambda
         enddo  !< ikx

         !> diagonalize Lambda to get the eigenvalue 
         call zgeev_pack(nfill, Lambda, Lambda_eig)
         do i=1, nfill
            WannierCenterKy(i, iky)= aimag(log(Lambda_eig(i)))/2d0/pi
         enddo

      enddo !< iky

      call mpi_allreduce(WannierCenterKy, WannierCenterKy_mpi, &
           size(WannierCenterKy), mpi_dp, mpi_sum, mpi_cmw, ierr)

      if (cpuid==0) then
         open(unit=101, file='wanniercenter.dat')

         do iky=1, Nky
            write(101, '(10000f16.8)') kpoints(2, 1, iky), &
               dmod(sum(WannierCenterKy_mpi(:, iky)), 1d0), & 
               WannierCenterKy_mpi(:, iky)
         enddo
         close(101)
      endif

      return
   end subroutine  wannier_center2D


! this suboutine is used for wannier center calculation for slab system

   subroutine  wannier_center2D_alt
      use para
      use mpi
      implicit none

      integer :: Nkx
      integer :: Nky
      integer :: knv2

      integer :: i
      integer :: j
      integer :: l 
      integer :: nfill

      integer :: ikx
      integer :: iky

      integer :: ierr

      !> k points in kx-ky plane
      complex(dp), allocatable :: kpoints(:, :, :)

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
   
      !> wannier centers for each ky, bands
      real(dp), allocatable :: WannierCenterKy(:, :)
      real(dp), allocatable :: WannierCenterKy_mpi(:, :)

      !> eigenvalue
      real(dp), allocatable :: eigenvalue(:)

      real(dp) :: kx
      real(dp) :: ky
      real(dp) :: k(2)

      Nkx= Nk
      Nky= 40
      knv2= Nkx*Nky

      nfill= 6*Nslab

      allocate(kpoints(2, Nkx, Nky))
      kpoints= 0d0

      allocate(Lambda_eig(nfill))
      allocate(Lambda(nfill, nfill))
      allocate(Lambda0(nfill, nfill))
      allocate(Mmnkb(nfill, nfill))
      allocate(Mmnkb_com(nfill, nfill))
      allocate(Mmnkb_full(Num_wann*Nslab, Num_wann*Nslab))
      allocate(hamk(Num_wann*Nslab, Num_wann*Nslab))
      allocate(hamk_dag(Num_wann*Nslab, Num_wann*Nslab))
      allocate(Eigenvector(Num_wann*Nslab, Num_wann*Nslab, Nkx))
      allocate(eigenvalue(Num_wann*Nslab))
      allocate(U(nfill, nfill))
      allocate(Sigma(nfill, nfill))
      allocate(VT(nfill, nfill))
      allocate(WannierCenterKy(nfill, Nky))
      allocate(WannierCenterKy_mpi(nfill, Nky))
      WannierCenterKy= 0d0
      WannierCenterKy_mpi= 0d0
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

      do iky=1, Nky
         do ikx=1, Nkx
            kx= (ikx-1)/real(Nkx)
            ky= (iky-1)/real(Nky)
            kpoints(1, ikx, iky)= kx
            kpoints(2, ikx, iky)= ky
         enddo
      enddo

      !> for each ky, we can get wanniercenter
      do iky=1+ cpuid, nky, num_cpu
         Lambda0=0d0
         do i=1, nfill
            Lambda0(i, i)= 1d0
         enddo

         if (cpuid==0) print *, iky, nky
         !> for each kx, we get the eigenvectors
         do ikx=1, nkx
            k(1)= kpoints(1, ikx, iky)
            k(2)= kpoints(2, ikx, iky)

            call ham_slab(k,hamk)

            !> diagonal hamk
            call eigensystem_c('V', 'U', Num_wann*Nslab, hamk, eigenvalue)

            Eigenvector(:, :, ikx)= hamk
         enddo

         !> sum over kx to get wanniercenters
         do ikx=1, nkx
            hamk= Eigenvector(:, :, ikx)
            hamk_dag= conjg(transpose(hamk))
            if (ikx==nkx) then
               hamk= Eigenvector(:, :, 1)
            else
               hamk= Eigenvector(:, :, ikx+ 1)
            endif

            !> <u_k|u_k+1>
            call mat_mul(Num_wann*Nslab, hamk_dag, hamk, Mmnkb_full)
            Mmnkb= Mmnkb_full(1:nfill, 1:nfill)
           !Mmnkb_com= 0d0
           !hamk_dag= Eigenvector(:, :, ikx)
           !hamk= Eigenvector(:, :, ikx+1)
           !do i=1, nfill
           !   do j=1, nfill
           !      do l= 1, Num_wann*Nslab
           !         Mmnkb_com(i, j)=  Mmnkb_com(i, j)+ conjg(hamk_dag(l, i))* hamk(l, j)
           !      enddo
           !   enddo
           !enddo

           !print *, maxval(real(Mmnkb-Mmnkb_com))
           !stop


            !> perform Singluar Value Decomposed of Mmnkb
            call zgesvd_pack(nfill, Mmnkb, U, Sigma, VT)

            !> after the calling of zgesvd_pack, Mmnkb becomes a temporal matrix
            U= conjg(transpose(U))
            VT= conjg(transpose(VT))
            call mat_mul(nfill, VT, U, Mmnkb)

            !> check hermicity
           !do i=1, nfill
           !   do j=i, nfill
           !      if (abs(Mmnkb(i, j)-conjg(Mmnkb(j, i)))>0.0001d0)then
           !         print *, 'Mmnkb is not Hermitian'
           !         print*, i, j, Mmnkb(i, j), Mmnkb(j, i)

           !      endif
           !   enddo
           !enddo

           !stop


            call mat_mul(nfill, Mmnkb, Lambda0, Lambda)
            Lambda0 = Lambda
         enddo  !< ikx

         !> diagonalize Lambda to get the eigenvalue 
         call zgeev_pack(nfill, Lambda, Lambda_eig)
         do i=1, nfill
            WannierCenterKy(i, iky)= aimag(log(Lambda_eig(i)))/2d0/pi
         enddo

      enddo !< iky

      call mpi_allreduce(WannierCenterKy, WannierCenterKy_mpi, &
           size(WannierCenterKy), mpi_dp, mpi_sum, mpi_cmw, ierr)

      if (cpuid==0) then
         open(unit=101, file='wanniercenter.dat')

         do iky=1, Nky
            write(101, '(10000f16.8)') kpoints(2, 1, iky), &
               sum(WannierCenterKy_mpi(:, iky))/2d0/pi, & 
               WannierCenterKy_mpi(:, iky)
         enddo
         close(101)
      endif

      return
   end subroutine  wannier_center2D_alt


