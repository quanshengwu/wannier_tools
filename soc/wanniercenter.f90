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
      real(dp), allocatable :: kpoints(:, :, :)

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

      nfill= Numoccupied*Nslab

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
!> calculate z2 

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
      integer :: m 
      integer :: ia
      integer :: ia1
      integer :: nfill

      integer :: ikx
      integer :: iky

      integer :: ierr

      !> k points in kx-ky plane
      real(dp), allocatable :: kpoints(:, :, :)

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

      !> atom position in the unit cell
      !> for slab system, dim=Nslab*Num_atoms
      real(dp), allocatable :: AtomsPosition_unitcell(:, :)
      real(dp), allocatable :: AtomsPosition_supercell(:, :)

      !> for each orbital, it correspond to an atom
      !> dim= Num_wann*Nslab
      integer, allocatable :: AtomIndex_orbital(:)

      real(dp) :: Umatrix_t(3, 3)

      integer :: imax
      real(dp) :: maxgap
      real(dp) :: maxgap0
      !> b.R
      real(dp) :: br

      real(dp) :: kx
      real(dp) :: ky
      real(dp) :: k(2)
      real(dp) :: b(2)

      !> larges gap of each two wannier centers for a given k point
      !> dim= Nky
      real(dp), allocatable :: largestgap(:)
      real(dp), allocatable :: largestgap_mpi(:)

      !> exp(-i*b.R)
      complex(dp) :: ratio

      !> Z2 calculation for time reversal invariant system
      integer :: Z2

      !> Chern number 
      real(dp) :: Chern

      Nkx= Nk
      Nky= 40
      knv2= Nkx*Nky

      nfill= Numoccupied*Nslab

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
      allocate(AtomsPosition_unitcell(3, Num_atoms))
      allocate(AtomsPosition_supercell(3, Nslab*Num_atoms))
      allocate(AtomIndex_orbital(Num_wann*Nslab))
      allocate(largestgap(Nky))
      allocate(largestgap_mpi(Nky))
      largestgap= 0d0
      largestgap_mpi= 0d0
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

      !> setup kpoints
      do iky=1, Nky
         do ikx=1, Nkx
            kx= (ikx-1d0)/real(Nkx)
            ky= (iky-1d0)/real(Nky)
            kpoints(1, ikx, iky)= kx
            kpoints(2, ikx, iky)= ky
            b(1)= 1.d0/real(Nkx)
            b(2)= 0.d0
         enddo
      enddo

      !> set up atom index for each orbitals in the basis
      if (soc>0) then  !> with spin orbital coupling
         l= 0
         do i=1, Nslab
            do ia=1, Num_atoms  !> spin up
               do j=1, nprojs(ia)
                  l= l+ 1
                  AtomIndex_orbital(l)= ia+ (i-1)*Num_atoms
               enddo ! l
            enddo ! ia
            do ia=1, Num_atoms  !> spin down
               do j=1, nprojs(ia)
                  l= l+ 1
                  AtomIndex_orbital(l)= ia+ (i-1)*Num_atoms
               enddo ! l
            enddo ! ia
         enddo ! i
      else  !> without spin orbital coupling
         l= 0
         do i=1, Nslab
            do ia=1, Num_atoms  !> spin down
               do j=1, nprojs(ia)
                  l= l+ 1
                  AtomIndex_orbital(l)= ia+ (i-1)*Num_atoms
               enddo ! l
            enddo ! ia
         enddo ! i

      endif

      Umatrix_t= transpose(Umatrix)
      call inv_r(3, Umatrix_t)

      !> set up atoms' position in the unit cell in the new basis
      !> only for 2D slab system
      do ia=1, Num_atoms
         do i=1, 3
         do j=1, 3
            AtomsPosition_unitcell(i, ia)= AtomsPosition_unitcell(i, ia)+ &
               Umatrix_t(i, j)*Atom_position(j, ia)
         enddo ! j
         enddo ! i
      enddo ! ia
     
      !> set up atoms' position in the supercell
      !> actually, we only need the first two coordinates
      ia1= 0
      do i=1, Nslab
         do ia=1, Num_atoms
            ia1= ia1+ 1
            AtomsPosition_supercell(1, ia1)= AtomsPosition_unitcell(1, ia)
            AtomsPosition_supercell(2, ia1)= AtomsPosition_unitcell(2, ia)
         enddo ! ia
      enddo ! i 

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
           !hamk= Eigenvector(:, :, ikx)
           !hamk_dag= conjg(transpose(hamk))
           !if (ikx==nkx) then
           !   hamk= Eigenvector(:, :, 1)
           !else
           !   hamk= Eigenvector(:, :, ikx+ 1)
           !endif

            !> <u_k|u_k+1>
           !call mat_mul(Num_wann*Nslab, hamk_dag, hamk, Mmnkb_full)
           !Mmnkb= Mmnkb_full(1:nfill, 1:nfill)
            Mmnkb= 0d0
            hamk_dag= Eigenvector(:, :, ikx)
            if (ikx==nkx) then
               hamk= Eigenvector(:, :, 1)
            else
               hamk= Eigenvector(:, :, ikx+ 1)
            endif
            do l=1, Nslab
               do m=1, Num_wann
                  ia= AtomIndex_orbital(m+(l-1)*Num_wann)
                  br= b(1)*AtomsPosition_supercell(1, ia)+ &
                      b(2)*AtomsPosition_supercell(2, ia)
                  ratio= cos(br)- zi* sin(br)
            
                  do i=1, nfill
                     do j=1, nfill
                        Mmnkb(i, j)=  Mmnkb(i, j)+ &
                           conjg(hamk_dag((l-1)*Num_wann+m, i))* &
                           hamk((l-1)*Num_wann+m, j)* ratio
                     enddo ! m
                  enddo ! l
               enddo ! j
            enddo ! i

            !> perform Singluar Value Decomposed of Mmnkb
            call zgesvd_pack(nfill, Mmnkb, U, Sigma, VT)

            !> after the calling of zgesvd_pack, Mmnkb becomes a temporal matrix
            U= conjg(transpose(U))
            VT= conjg(transpose(VT))
            call mat_mul(nfill, VT, U, Mmnkb)

            call mat_mul(nfill, Mmnkb, Lambda0, Lambda)
            Lambda0 = Lambda
         enddo  !< ikx

         !> diagonalize Lambda to get the eigenvalue 
         call zgeev_pack(nfill, Lambda, Lambda_eig)
         do i=1, nfill
            WannierCenterKy(i, iky)= aimag(log(Lambda_eig(i)))/2d0/pi
         enddo

         call sortheap(nfill, WannierCenterKy(:, iky))

         maxgap0= -99999d0
         do i=1, nfill
            if (i/=nfill) then
               maxgap= WannierCenterKy(i+1, iky)- WannierCenterKy(i, iky)
            else
               maxgap=1d0+ WannierCenterKy(1, iky)- WannierCenterKy(nfill, iky)
            endif

            if (maxgap>maxgap0) then
               maxgap0= maxgap
               imax= i
            endif

         enddo

         if (imax==nfill) then
            largestgap(iky)= (WannierCenterKy(1, iky)+ &
               WannierCenterKy(nfill, iky) -1d0)/2d0
         else
            largestgap(iky)= (WannierCenterKy(imax+1, iky)+ &
               WannierCenterKy(imax, iky))/2d0
         endif
      enddo !< iky

      call mpi_allreduce(WannierCenterKy, WannierCenterKy_mpi, &
           size(WannierCenterKy), mpi_dp, mpi_sum, mpi_cmw, ierr)
      call mpi_allreduce(largestgap, largestgap_mpi, &
           size(largestgap), mpi_dp, mpi_sum, mpi_cmw, ierr)



      if (cpuid==0) then
         open(unit=101, file='wanniercenter.dat')
         open(unit=102, file='largestgap.dat')

         do iky=1, Nky
            write(101, '(10000f16.8)') kpoints(2, 1, iky), &
               dmod(sum(WannierCenterKy_mpi(:, iky)), 1d0), & 
               WannierCenterKy_mpi(:, iky)
            write(102, '(10000f16.8)') kpoints(2, 1, iky), &
               largestgap_mpi(iky)
         enddo
         close(101)
         close(102)
      endif

      return
   end subroutine  wannier_center2D_alt




   !> this suboutine is used for wannier center calculation for 3D system
   !> only for one plane, calculate mirror chern number. only choose the bands
   !> which have the same mirror eigenvalue
   subroutine  wannier_center3D_plane_mirror
      use para
      use mpi
      implicit none

      integer :: Nk1
      integer :: Nk2
      integer :: knv2

      integer :: i
      integer :: j
      integer :: l 
      integer :: m 
      integer :: n
      integer :: ia
      integer :: ia1
      integer :: nfill
      integer :: nfill_half
      integer :: imax

      integer :: i1

      integer :: i2
      integer :: ik1
      integer :: ik2

      integer :: ierr

      !> k points in kx-ky plane
      real(dp), allocatable :: kpoints(:, :, :)

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

      complex(dp), allocatable :: mat1(:, :)
      complex(dp), allocatable :: mat2(:, :)

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

      !> larges gap of each two wannier centers for a given k point
      !> dim= Nky
      real(dp), allocatable :: largestgap(:)
      real(dp), allocatable :: largestgap_mpi(:)

      !> eigenvalue
      real(dp), allocatable :: eigenvalue(:)

      !> atom position in the unit cell
      !> for slab system, dim=Num_atoms
      real(dp), allocatable :: AtomsPosition_unitcell(:, :)
      real(dp), allocatable :: AtomsPosition_supercell(:, :)

      !> for each orbital, it correspond to an atom
      !> dim= Num_wann
      integer, allocatable :: AtomIndex_orbital(:)

      real(dp) :: Umatrix_t(3, 3)

      !> b.R
      real(dp) :: br

      !> exp(-i*b.R)
      complex(dp) :: ratio

      real(dp) :: kx
      real(dp) :: ky
      real(dp) :: k(3)
      real(dp) :: b(3)

      real(dp) :: maxgap
      real(dp) :: maxgap0

      !> Z2 calculation for time reversal invariant system
      integer :: Z2
      integer :: Delta

      real(dp) :: g
      real(dp) :: phi1
      real(dp) :: phi2
      real(dp) :: phi3
      real(dp) :: zm
      real(dp) :: zm1
      real(dp) :: xnm1
      real(dp) :: Deltam
      real(dp), allocatable :: xnm(:)
      real(dp) :: k11(3), k12(3)
      real(dp) :: k21(3), k22(3)

      !> mirror eigenvalue
      complex(dp), allocatable :: mirror_z_eig(:, :)
      !> the band index that has plus mirror number
      integer, allocatable :: mirror_plus(:, :)
      integer, allocatable :: mirror_minus(:, :)

      Nk1= Nk
      Nk2= Nk
      knv2= Nk1*Nk2

      nfill= Numoccupied
     !nfill_half= Numoccupied
      nfill_half= Numoccupied/2

      allocate(kpoints(3, Nk1, Nk2))
      kpoints= 0d0

      allocate(Lambda_eig(nfill_half))
      allocate(Lambda(nfill_half, nfill_half))
      allocate(Lambda0(nfill_half, nfill_half))
      allocate(Mmnkb(nfill_half, nfill_half))
      allocate(Mmnkb_com(nfill_half, nfill_half))
      allocate(Mmnkb_full(Num_wann, Num_wann))
      allocate(hamk(Num_wann, Num_wann))
      allocate(mat1(Num_wann, Num_wann))
      allocate(mat2(Num_wann, Num_wann))
      allocate(hamk_dag(Num_wann, Num_wann))
      allocate(Eigenvector(Num_wann, Num_wann, Nk1))
      allocate(eigenvalue(Num_wann))
      allocate(mirror_z_eig(nfill, Nk1))
      allocate(mirror_plus(nfill, Nk1))
      allocate(mirror_minus(nfill, Nk1))
      allocate(U(nfill_half, nfill_half))
      allocate(Sigma(nfill_half, nfill_half))
      allocate(VT(nfill_half, nfill_half))
      allocate(WannierCenterKy(nfill_half, Nk2))
      allocate(WannierCenterKy_mpi(nfill_half, Nk2))
      allocate(AtomIndex_orbital(Num_wann))
      allocate(xnm(nfill_half))
      allocate(largestgap(Nk2))
      allocate(largestgap_mpi(Nk2))
      mirror_minus= 0
      mirror_plus= 0
      largestgap= 0d0
      largestgap_mpi= 0d0
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

      !> set k plane
      !> the first dimension should be in one primitive plane
      k11=(/ 0.0d0,  0.0d0,  0.0d0/) ! 
      k12=(/ 1.0d0,  0.0d0,  0.0d0/) ! X
      k21=(/ 0.0d0,  0.0d0,  0.0d0/) ! 
      k22=(/ 0.0d0,  0.5d0,  0.0d0/) ! Y

      do ik2=1, Nk2
         do ik1=1, Nk1
            kpoints(:, ik1, ik2)= k11+(k12-k11)*(ik1-1)/dble(nk1)+  k21+ (k22-k21)*(ik2-1)/dble(nk2)
            b= (k12-k11)/dble(nk1)
         enddo
      enddo


      !> set up atom index for each orbitals in the basis
      if (soc>0) then  !> with spin orbital coupling
         l= 0
         do ia=1, Num_atoms  !> spin up
            do j=1, nprojs(ia)
               l= l+ 1
               AtomIndex_orbital(l)= ia
            enddo ! l
         enddo ! ia
         do ia=1, Num_atoms  !> spin down
            do j=1, nprojs(ia)
               l= l+ 1
               AtomIndex_orbital(l)= ia
            enddo ! l
         enddo ! ia
      else  !> without spin orbital coupling
         l= 0
         do ia=1, Num_atoms  !> spin down
            do j=1, nprojs(ia)
               l= l+ 1
               AtomIndex_orbital(l)= ia
            enddo ! l
         enddo ! ia

      endif

      Umatrix_t= transpose(Umatrix)
      call inv_r(3, Umatrix_t)

      !>> Get wannier center for ky=0 plane
      !> for each ky, we can get wanniercenter
      do ik2=1+ cpuid, Nk2, num_cpu
         if (cpuid.eq.0) print *,  'ik', ik2
         Lambda0=0d0
         do i=1, nfill_half
            Lambda0(i, i)= 1d0
         enddo

         mirror_plus= 0
         mirror_minus= 0
         !> for each k1, we get the eigenvectors
         do ik1=1, Nk1
            !if (cpuid==0) print *, ik1, ik2, Nk1, Nk2
            k= kpoints(:, ik1, ik2)

            call ham_bulk(k,hamk)

           !call mat_mul(Num_wann, mirror_z, hamk, mat1)
           !call mat_mul(Num_wann, mat1, mirror_z, mat2)
           !hamk= (hamk+ mat2)/2.d0

            !> diagonal hamk
            call eigensystem_c('V', 'U', Num_wann, hamk, eigenvalue)

            Eigenvector(:, :, ik1)= hamk

            mat2= conjg(transpose(hamk))

            !> calculate mirror eigenvalue
            call mat_mul(Num_wann, mat2, mirror_z, mat1)
            call mat_mul(Num_wann, mat1, hamk, mat2)
            
            !> get mirror_plus and mirror_minus
            do i=1, nfill
               if (abs(real(mat2(i, i))-1d0)< 1e-5) then
                  mirror_plus(i, ik1)= 1
               else
                  mirror_minus(i, ik1)= 1
               endif
            enddo

           !if (cpuid.eq.0)write(*, '(3f6.2,1000i3)')kpoints(:, ik1, ik2), &
           !   (mirror_plus(i, ik1), i=1, nfill) 
           !if (cpuid.eq.0)write(*, '(3f6.2,1000i3)')kpoints(:, ik1, ik2), &
           !   (mirror_minus(i, ik1), i=1, nfill) 
           !if (cpuid.eq.0)write(*, '(1000f6.2)')kpoints(:, ik1, ik2), &
           !   (real(mat2(i, i)), i=1, Num_wann) 
           !if (cpuid.eq.0)write(*, '(1000f5.1)') &
           !   (real(mirror_z(i, i)), i=1, Num_wann/2) 
           !if (cpuid.eq.0)write(*, '(1000f5.1)') &
           !   (real(mirror_z(i, i)), i=1+ Num_wann/2, Num_wann) 

           !write( *, *)' '
           !pause

         enddo

         !> sum over k1 to get wanniercenters
         do ik1=1, Nk1
            !> <u_k|u_k+1>
            Mmnkb= 0d0
            hamk_dag= Eigenvector(:, :, ik1)
            if (ik1==Nk1) then
               hamk= Eigenvector(:, :, 1)
            else
               hamk= Eigenvector(:, :, ik1+ 1)
            endif
            do m=1, Num_wann
               ia= AtomIndex_orbital(m)
               br= b(1)*Atom_position(1, ia)+ &
                   b(2)*Atom_position(2, ia)+ &
                   b(3)*Atom_position(3, ia)
               ratio= cos(br)- zi* sin(br)
              !ratio= 1d0
        
               i1= 0
               do j=1, nfill
                  if (mirror_plus(j, ik1)==1) cycle
                  i1= i1+ 1
                  i2= 0
                  do i=1, nfill
                     if (mirror_plus(i, ik1)==1) cycle
                     i2= i2+ 1
                     print *, ik1, i1, i2, i, j
                     Mmnkb(i2, i1)=  Mmnkb(i2, i1)+ &
                    !Mmnkb(i, j )=  Mmnkb(i , j )+ &
                        conjg(hamk_dag(m, i))* hamk(m, j)* ratio
                  enddo ! i
               enddo ! j
            enddo ! m

            !> perform Singluar Value Decomposed of Mmnkb
            call zgesvd_pack(nfill_half, Mmnkb, U, Sigma, VT)

            !> after the calling of zgesvd_pack, Mmnkb becomes a temporal matrix
            U= conjg(transpose(U))
            VT= conjg(transpose(VT))
            call mat_mul(nfill_half, VT, U, Mmnkb)

            call mat_mul(nfill_half, Mmnkb, Lambda0, Lambda)
            Lambda0 = Lambda
         enddo  !< ik1

         !> diagonalize Lambda to get the eigenvalue 
         call zgeev_pack(nfill_half, Lambda, Lambda_eig)
         do i=1, nfill_half
            WannierCenterKy(i, ik2)= aimag(log(Lambda_eig(i)))/2d0/pi
            WannierCenterKy(i, ik2)= mod(WannierCenterKy(i, ik2)+10d0, 1d0) 
         enddo

         call sortheap(nfill_half, WannierCenterKy(:, ik2))

         maxgap0= -99999d0
         imax= nfill_half
         do i=1, nfill_half
            if (i/=nfill_half) then
               maxgap= WannierCenterKy(i+1, ik2)- WannierCenterKy(i, ik2)
            else
               maxgap=1d0+ WannierCenterKy(1, ik2)- WannierCenterKy(nfill_half, ik2)
            endif

            if (maxgap>maxgap0) then
               maxgap0= maxgap
               imax= i
            endif

         enddo

         if (imax==nfill_half) then
            largestgap(ik2)= (WannierCenterKy(1, ik2)+ &
               WannierCenterKy(nfill_half, ik2) +1d0)/2d0
            largestgap(ik2)= mod(largestgap(ik2), 1d0)
         else
            largestgap(ik2)= (WannierCenterKy(imax+1, ik2)+ &
               WannierCenterKy(imax, ik2))/2d0
         endif


      enddo !< ik2

      WannierCenterKy_mpi= 0d0
      largestgap_mpi= 0d0
      call mpi_allreduce(WannierCenterKy, WannierCenterKy_mpi, &
           size(WannierCenterKy), mpi_dp, mpi_sum, mpi_cmw, ierr)
      call mpi_allreduce(largestgap, largestgap_mpi, &
           size(largestgap), mpi_dp, mpi_sum, mpi_cmw, ierr)

      if (cpuid==0) then
         open(unit=101, file='wanniercenterky0.dat')

         do ik2=1, Nk2
            write(101, '(10000f16.8)') dble(ik2-1)/dble(Nk2), &
               largestgap_mpi(ik2), dmod(sum(WannierCenterKy_mpi(:, ik2)), 1d0), & 
               WannierCenterKy_mpi(:, ik2)
         enddo
         close(101)
      endif


      !> Z2 calculation Alexey Soluyanov arXiv:1102.5600

      Delta= 0
      !> for each iky, we get a Deltam
      do ik2=1, nk2-1
      
         !> largestgap position
         zm= largestgap_mpi(ik2)
        !if (ik2==nk2) then
        !   zm1= largestgap_mpi(1)
        !   xnm= WannierCenterKy_mpi(1:nfill_half, 1)
        !else
            zm1= largestgap_mpi(ik2+1)
            xnm= WannierCenterKy_mpi(1:nfill_half, ik2+1)
        !endif
         Deltam= 1
         do i=1, nfill_half
            xnm1= xnm(i)
            phi1= 2d0*pi*zm
            phi2= 2d0*pi*zm1
            phi3= 2d0*pi*xnm1
            
            g= sin(phi2-phi1)+ sin(phi3-phi2)+  sin(phi1-phi3) 
            Deltam= Deltam* sign(1d0, g)
         enddo !i 
         if (Deltam<0) then
            Delta= Delta+ 1
         endif
      enddo !ik2

      Z2= mod(Delta, 2)

      if (cpuid==0) print*,'Z2 for ky=0: ', Z2

      return
   end subroutine  wannier_center3D_plane_mirror




   !> this suboutine is used for wannier center calculation for 3D system
   !> only for one plane

   subroutine  wannier_center3D_plane
      use para
      use mpi
      implicit none

      integer :: Nk1
      integer :: Nk2
      integer :: knv2

      integer :: i
      integer :: j
      integer :: l 
      integer :: m 
      integer :: ia
      integer :: ia1
      integer :: nfill
      integer :: imax

      integer :: ik1
      integer :: ik2

      integer :: ierr

      !> k points in kx-ky plane
      real(dp), allocatable :: kpoints(:, :, :)

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

      !> larges gap of each two wannier centers for a given k point
      !> dim= Nky
      real(dp), allocatable :: largestgap(:)
      real(dp), allocatable :: largestgap_mpi(:)

      !> eigenvalue
      real(dp), allocatable :: eigenvalue(:)

      !> atom position in the unit cell
      !> for slab system, dim=Num_atoms
      real(dp), allocatable :: AtomsPosition_unitcell(:, :)
      real(dp), allocatable :: AtomsPosition_supercell(:, :)

      !> for each orbital, it correspond to an atom
      !> dim= Num_wann
      integer, allocatable :: AtomIndex_orbital(:)

      real(dp) :: Umatrix_t(3, 3)

      !> b.R
      real(dp) :: br

      !> exp(-i*b.R)
      complex(dp) :: ratio

      real(dp) :: kx
      real(dp) :: ky
      real(dp) :: k(3)
      real(dp) :: b(3)

      real(dp) :: maxgap
      real(dp) :: maxgap0

      !> Z2 calculation for time reversal invariant system
      integer :: Z2
      integer :: Delta

      real(dp) :: g
      real(dp) :: phi1
      real(dp) :: phi2
      real(dp) :: phi3
      real(dp) :: zm
      real(dp) :: zm1
      real(dp) :: xnm1
      real(dp) :: Deltam
      real(dp), allocatable :: xnm(:)
      real(dp) :: k11(3), k12(3)
      real(dp) :: k21(3), k22(3)



      Nk1= Nk
      Nk2= Nk
      knv2= Nk1*Nk2

      nfill= Numoccupied

      allocate(kpoints(3, Nk1, Nk2))
      kpoints= 0d0

      allocate(Lambda_eig(nfill))
      allocate(Lambda(nfill, nfill))
      allocate(Lambda0(nfill, nfill))
      allocate(Mmnkb(nfill, nfill))
      allocate(Mmnkb_com(nfill, nfill))
      allocate(Mmnkb_full(Num_wann, Num_wann))
      allocate(hamk(Num_wann, Num_wann))
      allocate(hamk_dag(Num_wann, Num_wann))
      allocate(Eigenvector(Num_wann, Num_wann, Nk1))
      allocate(eigenvalue(Num_wann))
      allocate(U(nfill, nfill))
      allocate(Sigma(nfill, nfill))
      allocate(VT(nfill, nfill))
      allocate(WannierCenterKy(nfill, Nk2))
      allocate(WannierCenterKy_mpi(nfill, Nk2))
      allocate(AtomIndex_orbital(Num_wann))
      allocate(xnm(nfill))
      allocate(largestgap(Nk2))
      allocate(largestgap_mpi(Nk2))
      largestgap= 0d0
      largestgap_mpi= 0d0
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

      !> set k plane
      !> the first dimension should be in one primitive plane
      k11=(/ 0.0d0,  0.0d0,  0.0d0/) ! 
      k12=(/ 1.0d0,  0.0d0,  0.0d0/) ! X
      k21=(/ 0.0d0,  0.0d0,  0.0d0/) ! 
      k22=(/ 0.0d0,  0.5d0,  0.0d0/) ! Y

      do ik2=1, Nk2
         do ik1=1, Nk1
            kpoints(:, ik1, ik2)= k11+(k12-k11)*(ik1-1)/dble(nk1)+  k21+ (k22-k21)*(ik2-1)/dble(nk2)
            b= (k12-k11)/dble(nk1)
         enddo
      enddo


      !> set up atom index for each orbitals in the basis
      if (soc>0) then  !> with spin orbital coupling
         l= 0
         do ia=1, Num_atoms  !> spin up
            do j=1, nprojs(ia)
               l= l+ 1
               AtomIndex_orbital(l)= ia
            enddo ! l
         enddo ! ia
         do ia=1, Num_atoms  !> spin down
            do j=1, nprojs(ia)
               l= l+ 1
               AtomIndex_orbital(l)= ia
            enddo ! l
         enddo ! ia
      else  !> without spin orbital coupling
         l= 0
         do ia=1, Num_atoms  !> spin down
            do j=1, nprojs(ia)
               l= l+ 1
               AtomIndex_orbital(l)= ia
            enddo ! l
         enddo ! ia

      endif

      Umatrix_t= transpose(Umatrix)
      call inv_r(3, Umatrix_t)

      !>> Get wannier center for ky=0 plane
      !> for each ky, we can get wanniercenter
      do ik2=1+ cpuid, Nk2, num_cpu
         if (cpuid.eq.0) print *,  'ik', ik2
         Lambda0=0d0
         do i=1, nfill
            Lambda0(i, i)= 1d0
         enddo

         !> for each k1, we get the eigenvectors
         do ik1=1, Nk1
            !if (cpuid==0) print *, ik1, ik2, Nk1, Nk2
            k= kpoints(:, ik1, ik2)

            call ham_bulk(k,hamk)

            !> diagonal hamk
            call eigensystem_c('V', 'U', Num_wann, hamk, eigenvalue)

            Eigenvector(:, :, ik1)= hamk
         enddo

         !> sum over k1 to get wanniercenters
         do ik1=1, Nk1
            !> <u_k|u_k+1>
            Mmnkb= 0d0
            hamk_dag= Eigenvector(:, :, ik1)
            if (ik1==Nk1) then
               hamk= Eigenvector(:, :, 1)
            else
               hamk= Eigenvector(:, :, ik1+ 1)
            endif
            do m=1, Num_wann
               ia= AtomIndex_orbital(m)
               br= b(1)*Atom_position(1, ia)+ &
                   b(2)*Atom_position(2, ia)+ &
                   b(3)*Atom_position(3, ia)
               ratio= cos(br)- zi* sin(br)
              !ratio= 1d0
         
               do j=1, nfill
                  do i=1, nfill
                     Mmnkb(i, j)=  Mmnkb(i, j)+ &
                        conjg(hamk_dag(m, i))* hamk(m, j)* ratio
                  enddo ! i
               enddo ! j
            enddo ! m

            !> perform Singluar Value Decomposed of Mmnkb
            call zgesvd_pack(nfill, Mmnkb, U, Sigma, VT)

            !> after the calling of zgesvd_pack, Mmnkb becomes a temporal matrix
            U= conjg(transpose(U))
            VT= conjg(transpose(VT))
            call mat_mul(nfill, VT, U, Mmnkb)

            call mat_mul(nfill, Mmnkb, Lambda0, Lambda)
            Lambda0 = Lambda
         enddo  !< ik1

         !> diagonalize Lambda to get the eigenvalue 
         call zgeev_pack(nfill, Lambda, Lambda_eig)
         do i=1, nfill
            WannierCenterKy(i, ik2)= aimag(log(Lambda_eig(i)))/2d0/pi
            WannierCenterKy(i, ik2)= mod(WannierCenterKy(i, ik2)+10d0, 1d0) 
         enddo

         call sortheap(nfill, WannierCenterKy(:, ik2))

         maxgap0= -99999d0
         imax= nfill
         do i=1, nfill
            if (i/=nfill) then
               maxgap= WannierCenterKy(i+1, ik2)- WannierCenterKy(i, ik2)
            else
               maxgap=1d0+ WannierCenterKy(1, ik2)- WannierCenterKy(nfill, ik2)
            endif

            if (maxgap>maxgap0) then
               maxgap0= maxgap
               imax= i
            endif

         enddo

         if (imax==nfill) then
            largestgap(ik2)= (WannierCenterKy(1, ik2)+ &
               WannierCenterKy(nfill, ik2) +1d0)/2d0
            largestgap(ik2)= mod(largestgap(ik2), 1d0)
         else
            largestgap(ik2)= (WannierCenterKy(imax+1, ik2)+ &
               WannierCenterKy(imax, ik2))/2d0
         endif


      enddo !< ik2

      WannierCenterKy_mpi= 0d0
      largestgap_mpi= 0d0
      call mpi_allreduce(WannierCenterKy, WannierCenterKy_mpi, &
           size(WannierCenterKy), mpi_dp, mpi_sum, mpi_cmw, ierr)
      call mpi_allreduce(largestgap, largestgap_mpi, &
           size(largestgap), mpi_dp, mpi_sum, mpi_cmw, ierr)

      if (cpuid==0) then
         open(unit=101, file='wanniercenterky0.dat')

         do ik2=1, Nk2
            write(101, '(10000f16.8)') dble(ik2-1)/dble(Nk2), &
               largestgap_mpi(ik2), dmod(sum(WannierCenterKy_mpi(:, ik2)), 1d0), & 
               WannierCenterKy_mpi(:, ik2)
         enddo
         close(101)
      endif


      !> Z2 calculation Alexey Soluyanov arXiv:1102.5600

      Delta= 0
      !> for each iky, we get a Deltam
      do ik2=1, nk2-1
      
         !> largestgap position
         zm= largestgap_mpi(ik2)
        !if (ik2==nk2) then
        !   zm1= largestgap_mpi(1)
        !   xnm= WannierCenterKy_mpi(1:nfill, 1)
        !else
            zm1= largestgap_mpi(ik2+1)
            xnm= WannierCenterKy_mpi(1:nfill, ik2+1)
        !endif
         Deltam= 1
         do i=1, nfill
            xnm1= xnm(i)
            phi1= 2d0*pi*zm
            phi2= 2d0*pi*zm1
            phi3= 2d0*pi*xnm1
            
            g= sin(phi2-phi1)+ sin(phi3-phi2)+  sin(phi1-phi3) 
            Deltam= Deltam* sign(1d0, g)
         enddo !i 
         if (Deltam<0) then
            Delta= Delta+ 1
         endif
      enddo !ik2

      Z2= mod(Delta, 2)

      if (cpuid==0) print*,'Z2 for ky=0: ', Z2

      return
   end subroutine  wannier_center3D_plane



! this suboutine is used for wannier center calculation for 3D system

   subroutine  wannier_center3D
      use para
      use mpi
      implicit none

      integer :: Nk1
      integer :: Nk2
      integer :: knv2

      integer :: i
      integer :: j
      integer :: l 
      integer :: m 
      integer :: ia
      integer :: ia1
      integer :: nfill
      integer :: imax

      integer :: ik1
      integer :: ik2

      integer :: ierr

      !> k points in kx-ky plane
      real(dp), allocatable :: kpoints(:, :, :)

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

      !> larges gap of each two wannier centers for a given k point
      !> dim= Nky
      real(dp), allocatable :: largestgap(:)
      real(dp), allocatable :: largestgap_mpi(:)

      !> eigenvalue
      real(dp), allocatable :: eigenvalue(:)

      !> atom position in the unit cell
      !> for slab system, dim=Num_atoms
      real(dp), allocatable :: AtomsPosition_unitcell(:, :)
      real(dp), allocatable :: AtomsPosition_supercell(:, :)

      !> for each orbital, it correspond to an atom
      !> dim= Num_wann
      integer, allocatable :: AtomIndex_orbital(:)

      real(dp) :: Umatrix_t(3, 3)

      !> b.R
      real(dp) :: br

      !> exp(-i*b.R)
      complex(dp) :: ratio

      real(dp) :: kx
      real(dp) :: ky
      real(dp) :: k(3)
      real(dp) :: b(3)

      real(dp) :: maxgap
      real(dp) :: maxgap0

      !> Z2 calculation for time reversal invariant system
      integer :: Z2
      integer :: Delta

      real(dp) :: g
      real(dp) :: phi1
      real(dp) :: phi2
      real(dp) :: phi3
      real(dp) :: zm
      real(dp) :: zm1
      real(dp) :: xnm1
      real(dp) :: Deltam
      real(dp), allocatable :: xnm(:)


      Nk1= Nk
      Nk2= Nk
      knv2= Nk1*Nk2

      nfill= Numoccupied

      allocate(kpoints(2, Nk1, Nk2))
      kpoints= 0d0

      allocate(Lambda_eig(nfill))
      allocate(Lambda(nfill, nfill))
      allocate(Lambda0(nfill, nfill))
      allocate(Mmnkb(nfill, nfill))
      allocate(Mmnkb_com(nfill, nfill))
      allocate(Mmnkb_full(Num_wann, Num_wann))
      allocate(hamk(Num_wann, Num_wann))
      allocate(hamk_dag(Num_wann, Num_wann))
      allocate(Eigenvector(Num_wann, Num_wann, Nk1))
      allocate(eigenvalue(Num_wann))
      allocate(U(nfill, nfill))
      allocate(Sigma(nfill, nfill))
      allocate(VT(nfill, nfill))
      allocate(WannierCenterKy(nfill, Nk2))
      allocate(WannierCenterKy_mpi(nfill, Nk2))
      allocate(AtomIndex_orbital(Num_wann))
      allocate(xnm(nfill))
      allocate(largestgap(Nk2))
      allocate(largestgap_mpi(Nk2))
      largestgap= 0d0
      largestgap_mpi= 0d0
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

      !> setup kpoints
      do ik2=1, Nk2
         do ik1=1, Nk1
            kx= (ik1-1d0)/real(Nk1)
            ky= (ik2-1d0)/real(Nk2-1)/2d0
            kpoints(1, ik1, ik2)= kx
            kpoints(2, ik1, ik2)= ky
            b(1)= 1.d0/real(Nk1)
            b(2)= 0.d0
            b(3)= 0.d0
         enddo
      enddo

      !> set up atom index for each orbitals in the basis
      if (soc>0) then  !> with spin orbital coupling
         l= 0
         do ia=1, Num_atoms  !> spin up
            do j=1, nprojs(ia)
               l= l+ 1
               AtomIndex_orbital(l)= ia
            enddo ! l
         enddo ! ia
         do ia=1, Num_atoms  !> spin down
            do j=1, nprojs(ia)
               l= l+ 1
               AtomIndex_orbital(l)= ia
            enddo ! l
         enddo ! ia
      else  !> without spin orbital coupling
         l= 0
         do ia=1, Num_atoms  !> spin down
            do j=1, nprojs(ia)
               l= l+ 1
               AtomIndex_orbital(l)= ia
            enddo ! l
         enddo ! ia

      endif

      Umatrix_t= transpose(Umatrix)
      call inv_r(3, Umatrix_t)

      !>> Get wannier center for ky=0 plane
      !> for each ky, we can get wanniercenter
      do ik2=1+ cpuid, Nk2, num_cpu
         Lambda0=0d0
         do i=1, nfill
            Lambda0(i, i)= 1d0
         enddo

         !> for each k1, we get the eigenvectors
         do ik1=1, Nk1
            if (cpuid==0) print *, ik1, ik2, Nk1, Nk2
            k(1)= kpoints(1, ik1, ik2)
            k(2)= 0d0
            k(3)= kpoints(2, ik1, ik2)

            call ham_bulk(k,hamk)

            !> diagonal hamk
            call eigensystem_c('V', 'U', Num_wann, hamk, eigenvalue)

            Eigenvector(:, :, ik1)= hamk
         enddo

         !> sum over k1 to get wanniercenters
         do ik1=1, Nk1
            !> <u_k|u_k+1>
            Mmnkb= 0d0
            hamk_dag= Eigenvector(:, :, ik1)
            if (ik1==Nk1) then
               hamk= Eigenvector(:, :, 1)
            else
               hamk= Eigenvector(:, :, ik1+ 1)
            endif
            do m=1, Num_wann
               ia= AtomIndex_orbital(m)
               br= b(1)*Atom_position(1, ia)+ &
                   b(2)*Atom_position(2, ia)+ &
                   b(3)*Atom_position(3, ia)
              !ratio= cos(br)- zi* sin(br)
               ratio= 1d0
         
               do j=1, nfill
                  do i=1, nfill
                     Mmnkb(i, j)=  Mmnkb(i, j)+ &
                        conjg(hamk_dag(m, i))* hamk(m, j)* ratio
                  enddo ! i
               enddo ! j
            enddo ! m

            !> perform Singluar Value Decomposed of Mmnkb
            call zgesvd_pack(nfill, Mmnkb, U, Sigma, VT)

            !> after the calling of zgesvd_pack, Mmnkb becomes a temporal matrix
            U= conjg(transpose(U))
            VT= conjg(transpose(VT))
            call mat_mul(nfill, VT, U, Mmnkb)

            call mat_mul(nfill, Mmnkb, Lambda0, Lambda)
            Lambda0 = Lambda
         enddo  !< ik1

         !> diagonalize Lambda to get the eigenvalue 
         call zgeev_pack(nfill, Lambda, Lambda_eig)
         do i=1, nfill
            WannierCenterKy(i, ik2)= aimag(log(Lambda_eig(i)))/2d0/pi
            WannierCenterKy(i, ik2)= mod(WannierCenterKy(i, ik2)+10d0, 1d0) 
         enddo

         call sortheap(nfill, WannierCenterKy(:, ik2))

         maxgap0= -99999d0
         imax= nfill
         do i=1, nfill
            if (i/=nfill) then
               maxgap= WannierCenterKy(i+1, ik2)- WannierCenterKy(i, ik2)
            else
               maxgap=1d0+ WannierCenterKy(1, ik2)- WannierCenterKy(nfill, ik2)
            endif

            if (maxgap>maxgap0) then
               maxgap0= maxgap
               imax= i
            endif

         enddo

         if (imax==nfill) then
            largestgap(ik2)= (WannierCenterKy(1, ik2)+ &
               WannierCenterKy(nfill, ik2) +1d0)/2d0
            largestgap(ik2)= mod(largestgap(ik2), 1d0)
         else
            largestgap(ik2)= (WannierCenterKy(imax+1, ik2)+ &
               WannierCenterKy(imax, ik2))/2d0
         endif


      enddo !< ik2

      WannierCenterKy_mpi= 0d0
      largestgap_mpi= 0d0
      call mpi_allreduce(WannierCenterKy, WannierCenterKy_mpi, &
           size(WannierCenterKy), mpi_dp, mpi_sum, mpi_cmw, ierr)
      call mpi_allreduce(largestgap, largestgap_mpi, &
           size(largestgap), mpi_dp, mpi_sum, mpi_cmw, ierr)

      if (cpuid==0) then
         open(unit=101, file='wanniercenterky0.dat')

         do ik2=1, Nk2
            write(101, '(10000f16.8)') kpoints(2, 1, ik2), &
               largestgap_mpi(ik2), dmod(sum(WannierCenterKy_mpi(:, ik2)), 1d0), & 
               WannierCenterKy_mpi(:, ik2)
         enddo
         close(101)
      endif


      !> Z2 calculation Alexey Soluyanov arXiv:1102.5600

      Delta= 0
      !> for each iky, we get a Deltam
      do ik2=1, nk2-1
      
         !> largestgap position
         zm= largestgap_mpi(ik2)
        !if (ik2==nk2) then
        !   zm1= largestgap_mpi(1)
        !   xnm= WannierCenterKy_mpi(1:nfill, 1)
        !else
            zm1= largestgap_mpi(ik2+1)
            xnm= WannierCenterKy_mpi(1:nfill, ik2+1)
        !endif
         Deltam= 1
         do i=1, nfill
            xnm1= xnm(i)
            phi1= 2d0*pi*zm
            phi2= 2d0*pi*zm1
            phi3= 2d0*pi*xnm1
            
            g= sin(phi2-phi1)+ sin(phi3-phi2)+  sin(phi1-phi3) 
            Deltam= Deltam* sign(1d0, g)
         enddo !i 
         if (Deltam<0) then
            Delta= Delta+ 1
         endif
      enddo !ik2

      Z2= mod(Delta, 2)

      if (cpuid==0) print*,'Z2 for ky=0: ', Z2, Delta


      !>> Get wannier center for ky=0.5 plane
      !> for each ky, we can get wanniercenter
      do ik2=1+ cpuid, Nk2, num_cpu
         Lambda0=0d0
         do i=1, nfill
            Lambda0(i, i)= 1d0
         enddo

         !> for each k1, we get the eigenvectors
         do ik1=1, Nk1
            if (cpuid==0) print *, ik1, ik2, Nk1, Nk2
            k(1)= kpoints(1, ik1, ik2)
            k(2)= 0.5d0
            k(3)= kpoints(2, ik1, ik2)

            call ham_bulk(k,hamk)

            !> diagonal hamk
            call eigensystem_c('V', 'U', Num_wann, hamk, eigenvalue)

            Eigenvector(:, :, ik1)= hamk
         enddo

         !> sum over k1 to get wanniercenters
         do ik1=1, Nk1
            !> <u_k|u_k+1>
            Mmnkb= 0d0
            hamk_dag= Eigenvector(:, :, ik1)
            if (ik1==Nk1) then
               hamk= Eigenvector(:, :, 1)
            else
               hamk= Eigenvector(:, :, ik1+ 1)
            endif
            do m=1, Num_wann
               ia= AtomIndex_orbital(m)
               br= b(1)*Atom_position(1, ia)+ &
                   b(2)*Atom_position(2, ia)+ &
                   b(3)*Atom_position(3, ia)
              !ratio= cos(br)- zi* sin(br)
               ratio= 1d0
         
               do j=1, nfill
                  do i=1, nfill
                     Mmnkb(i, j)=  Mmnkb(i, j)+ &
                        conjg(hamk_dag(m, i))* hamk(m, j)* ratio
                  enddo ! i
               enddo ! j
            enddo ! m

            !> perform Singluar Value Decomposed of Mmnkb
            call zgesvd_pack(nfill, Mmnkb, U, Sigma, VT)

            !> after the calling of zgesvd_pack, Mmnkb becomes a temporal matrix
            U= conjg(transpose(U))
            VT= conjg(transpose(VT))
            call mat_mul(nfill, VT, U, Mmnkb)

            call mat_mul(nfill, Mmnkb, Lambda0, Lambda)
            Lambda0 = Lambda
         enddo  !< ik1

         !> diagonalize Lambda to get the eigenvalue 
         call zgeev_pack(nfill, Lambda, Lambda_eig)
         do i=1, nfill
            WannierCenterKy(i, ik2)= aimag(log(Lambda_eig(i)))/2d0/pi
            WannierCenterKy(i, ik2)= mod(WannierCenterKy(i, ik2)+10d0, 1d0) 
         enddo

         call sortheap(nfill, WannierCenterKy(:, ik2))

         maxgap0= -99999d0
         imax= nfill
         do i=1, nfill
            if (i/=nfill) then
               maxgap= WannierCenterKy(i+1, ik2)- WannierCenterKy(i, ik2)
            else
               maxgap=1d0+ WannierCenterKy(1, ik2)- WannierCenterKy(nfill, ik2)
            endif

            if (maxgap>maxgap0) then
               maxgap0= maxgap
               imax= i
            endif

         enddo

         if (imax==nfill) then
            largestgap(ik2)= (WannierCenterKy(1, ik2)+ &
               WannierCenterKy(nfill, ik2) +1d0)/2d0
            largestgap(ik2)= mod(largestgap(ik2), 1d0)
         else
            largestgap(ik2)= (WannierCenterKy(imax+1, ik2)+ &
               WannierCenterKy(imax, ik2))/2d0
         endif

      enddo !< ik2

      WannierCenterKy_mpi= 0d0
      largestgap_mpi= 0d0
      call mpi_allreduce(WannierCenterKy, WannierCenterKy_mpi, &
           size(WannierCenterKy), mpi_dp, mpi_sum, mpi_cmw, ierr)
      call mpi_allreduce(largestgap, largestgap_mpi, &
           size(largestgap), mpi_dp, mpi_sum, mpi_cmw, ierr)

      if (cpuid==0) then
         open(unit=101, file='wanniercenterky05.dat')

         do ik2=1, Nk2
            write(101, '(10000f16.8)') kpoints(2, 1, ik2), &
               largestgap_mpi(ik2), &
               dmod(sum(WannierCenterKy_mpi(:, ik2)), 1d0), & 
               WannierCenterKy_mpi(:, ik2)
         enddo
         close(101)
      endif


      !> Z2 calculation Alexey Soluyanov arXiv:1102.5600

      Delta= 0
      !> for each iky, we get a Deltam
      do ik2=1, nk2- 1
      
         !> largestgap position
         zm= largestgap_mpi(ik2)
         zm1= largestgap_mpi(ik2+1)
         xnm= WannierCenterKy_mpi(1:nfill, ik2+1)
         Deltam= 1
         do i=1, nfill
            xnm1= xnm(i)
            phi1= 2*pi*zm
            phi2= 2*pi*zm1
            phi3= 2*pi*xnm1
            
            g= sin(phi2-phi1)+ sin(phi3-phi2)+  sin(phi1-phi3) 
            Deltam= Deltam* sign(1d0, g)
         enddo !i 
         if (Deltam<0) then
            Delta= Delta+ 1
         endif
      enddo !ik2

      Z2= mod(Delta, 2)

      if (cpuid==0) print*,'Z2 for ky=0.5: ', Z2

      return
   end subroutine  wannier_center3D

   subroutine sortheap(n, arr)
      implicit none
      integer, intent(in) :: n
      real(8), intent(inout) :: arr(n)

      !* local variables
      integer :: i

      do i=n/2, 1, -1
         call sift_down(i, n)
      enddo

      do i=n, 2, -1
         call swap(arr(1), arr(i))
         call sift_down(1, i-1)
      enddo
      contains
      subroutine sift_down(l, r)
         integer, intent(in) :: l, r
         integer :: j, jold
         real(8) :: a
         a= arr(l)
         jold= l
         j= l+ l

         do while (j<=r)
            if (j<r) then
               if (arr(j)<arr(j+1))j=j+1
            endif
            if (a>= arr(j))exit
            arr(jold)= arr(j)
            jold= j
            j= j+ j
         enddo
         arr(jold)= a
         return
      end subroutine sift_down
   end subroutine sortheap

   !>> swap two real numbers
   subroutine swap(a, b)
      real(8), intent(inout) :: a
      real(8), intent(inout) :: b
      real(8) :: c
      c=a
      a=b
      b=c
      return
   end subroutine swap


