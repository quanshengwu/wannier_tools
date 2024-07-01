! this suboutine is used for wannier center calculation for slab system
! Copyright (c) 2010 QuanSheng Wu. All rights reserved.

 !> this suboutine is used for wannier center calculation for 3D system
 !> only for one plane, calculate mirror chern number. only choose the bands
 !> which have the same mirror eigenvalue
 !> when calling this, please make sure that the mirror operator matrix is
 !> properly set in the symmetry.f90.
 !> You can check it with ek_bulk_mirror_z subroutine in ek_bulk.f90
subroutine  wannier_center3D_plane_mirror
   use para
   use wmpi
   implicit none


   integer :: i, j, l, m, ia, nfill, nfill_half

   integer :: i1, i2, ik1, ik2, ikp, ierr

   !> k points in kx-ky plane
   real(dp), allocatable :: kpoints(:, :, :)

   !> hamiltonian for each k point
   !> and also the eigenvector of hamiltonian after eigensystem_c
   complex(dp), allocatable :: Hamk(:, :), hamk_dag(:, :)

   !> eigenvector for each kx
   complex(dp), allocatable :: Eigenvector(:, :, :)

   !> Mmnkb=<u_n(k)|u_m(k+b)>
   !> |u_n(k)> is the periodic part of wave function
   complex(dp), allocatable :: Mmnkb(:, :), Mmnkb_com(:, :)

   complex(dp), allocatable :: mat1(:, :), mat2(:, :)

   complex(dp), allocatable :: Lambda_eig(:), Lambda(:, :), Lambda0(:, :)

   !> three matrix for SVD
   !> M= U.Sigma.V^\dag
   !> VT= V^\dag
   complex(dp), allocatable :: U(:, :), VT(:, :)
   real   (dp), allocatable :: Sigma(:, :)

   !> wannier centers for each ky, bands
   real(dp), allocatable :: WannierCenterKy_minus(:, :),WannierCenterKy_minus_mpi(:, :)
   real(dp), allocatable :: WannierCenterKy_plus(:, :),WannierCenterKy_plus_mpi(:, :)

   !> sumation for Wannier charge center
   real(dp) :: gap_sum, gap_step
   real(dp), allocatable :: wcc_sum(:)

   !> eigenvalue
   real(dp), allocatable :: eigenvalue(:)

   !> for each orbital, it correspond to an atom
   !> dim= Num_wann
   integer, allocatable :: AtomIndex_orbital(:)

   real(dp) :: Umatrix_t(3, 3)

   !> b.R
   real(dp) :: br

   !> exp(-i*b.R)
   complex(dp) :: ratio

   real(dp) :: k(3), b(3)

   real(dp), allocatable :: xnm(:)
   real(dp) :: k0(3), k1(3), k2(3)

   !> mirror eigenvalue
   complex(dp), allocatable :: mirror_z_eig(:, :)

   !> the band index that has plus mirror number
   logical, allocatable :: mirror_plus(:, :), mirror_minus(:, :)

   nfill= NumberofSelectedOccupiedBands
   nfill_half= NumberofSelectedOccupiedBands/2

   allocate(kpoints(3, Nk1, Nk2))
   kpoints= 0d0

   allocate(Lambda_eig(nfill_half))
   allocate(Lambda(nfill_half, nfill_half))
   allocate(Lambda0(nfill_half, nfill_half))
   allocate(Mmnkb(nfill_half, nfill_half))
   allocate(Mmnkb_com(nfill_half, nfill_half))
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
   allocate(WannierCenterKy_minus(nfill_half, Nk2))
   allocate(WannierCenterKy_minus_mpi(nfill_half, Nk2))
   allocate(WannierCenterKy_plus(nfill_half, Nk2))
   allocate(WannierCenterKy_plus_mpi(nfill_half, Nk2))
   allocate(wcc_sum(Nk2))
   allocate(AtomIndex_orbital(Num_wann))
   allocate(xnm(nfill_half))
   mirror_minus= .False.
   mirror_plus= .False.
   WannierCenterKy_minus= 0d0
   WannierCenterKy_minus_mpi= 0d0
   WannierCenterKy_plus= 0d0
   WannierCenterKy_plus_mpi= 0d0
   hamk=0d0
   eigenvalue=0d0
   Eigenvector=0d0
   Mmnkb=0d0
   Mmnkb_com=0d0
   Lambda =0d0
   Lambda0=0d0
   U= 0d0
   Sigma= 0d0
   VT= 0d0

   !> set k plane
   !> the first dimension should be in one primitive cell, [0, 2*pi]
   k0= K3D_start
   k1= K3D_vec1
   k2= K3D_vec2

   do ik2=1, Nk2
      do ik1=1, Nk1
         kpoints(:, ik1, ik2)= k0+k1*dble(ik1-1.d0)/dble(nk1)+ k2*(ik2-1)/dble(nk2-1d0)
      enddo
   enddo
   b= dble(k1)/dble(nk1)
   b= b(1)*Origin_cell%kua+b(2)*Origin_cell%kub+b(3)*Origin_cell%kuc


   !> set up atom index for each orbitals in the basis
   if (soc>0) then  !> with spin orbital coupling
      l= 0
      do ia=1, Origin_cell%Num_atoms  !> spin up
         do j=1, Origin_cell%nprojs(ia)
            l= l+ 1
            AtomIndex_orbital(l)= ia
         enddo ! l
      enddo ! ia
      do ia=1, Origin_cell%Num_atoms  !> spin down
         do j=1, Origin_cell%nprojs(ia)
            l= l+ 1
            AtomIndex_orbital(l)= ia
         enddo ! l
      enddo ! ia
   else  !> without spin orbital coupling
      l= 0
      do ia=1, Origin_cell%Num_atoms  !> spin down
         do j=1, Origin_cell%nprojs(ia)
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
      if (cpuid.eq.0) write(stdout, *) 'ik', ik2

      mirror_plus= .False.
      mirror_minus= .False.
      !> for each k1, we get the eigenvectors
      do ik1=1, Nk1
         k= kpoints(:, ik1, ik2)


         ! generate bulk Hamiltonian
         if (index(KPorTB, 'KP')/=0)then
            call ham_bulk_kp(k, Hamk)
         else
            !> deal with phonon system
            if (index(Particle,'phonon')/=0.and.LOTO_correction) then
               call ham_bulk_LOTO(k, Hamk)
            else
               call ham_bulk_latticegauge    (k, Hamk)
              !call ham_bulk_atomicgauge    (k, Hamk)
            endif
         endif


         !> symmetrization
         call mat_mul(Num_wann, mirror_z, hamk, mat1)
         call mat_mul(Num_wann, mat1, mirror_z, mat2)
         mat1= (hamk+ mat2)/2.d0
         hamk= mat1

         !> diagonal hamk
         call eigensystem_c('V', 'U', Num_wann, hamk, eigenvalue)

         Eigenvector(:, :, ik1)= hamk

         mat2= conjg(transpose(hamk))

         !> calculate mirror eigenvalue
         call mat_mul(Num_wann, mat2, mirror_z, mat1)
         call mat_mul(Num_wann, mat1, hamk, mat2)

         !> get mirror_plus and mirror_minus
         do i=1, nfill
            if (abs(real(mat2(i, i))-1d0)< 1e-3) then
               mirror_plus(i, ik1)= .true.
            else
               mirror_minus(i, ik1)= .true.
            endif
         enddo

      enddo ! ik1

      !> sum over k1 to get wanniercenters for mirror plus
      Lambda0=0d0
      do i=1, nfill_half
         Lambda0(i, i)= 1d0
      enddo
      do ik1=1, Nk1
         !> <u_k|u_k+1>
         Mmnkb= 0d0
         hamk_dag= Eigenvector(:, :, ik1)
         if (ik1==Nk1) then
            hamk= Eigenvector(:, :, 1)
            ikp= 1
         else
            hamk= Eigenvector(:, :, ik1+ 1)
            ikp= ik1+ 1
         endif
         do m=1, Num_wann
            !ia= AtomIndex_orbital(m)
            !br= b(1)*Atom_position(1, ia)+ &
            !    b(2)*Atom_position(2, ia)+ &
            !    b(3)*Atom_position(3, ia)
            br= b(1)*Origin_cell%wannier_centers_cart(1, m )+ &
               b(2)*Origin_cell%wannier_centers_cart(2, m )+ &
               b(3)*Origin_cell%wannier_centers_cart(3, m )
            ratio= cos(br)- zi* sin(br)

            i1= 0
           !do j=1, nfill
            do j=Selected_Occupiedband_index(1), Selected_Occupiedband_index(NumberofSelectedOccupiedBands)
               if (mirror_minus(j, ikp)) cycle
               i1= i1+ 1
               i2= 0
              !do i=1, nfill
               do i=Selected_Occupiedband_index(1), Selected_Occupiedband_index(NumberofSelectedOccupiedBands)
                  if (mirror_minus(i, ik1)) cycle
                  i2= i2+ 1
                  Mmnkb(i2, i1)=  Mmnkb(i2, i1)+ &
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
         WannierCenterKy_plus(i, ik2)= aimag(log(Lambda_eig(i)))/2d0/pi
         WannierCenterKy_plus(i, ik2)= mod(WannierCenterKy_plus(i, ik2)+10d0, 1d0)
      enddo

      call sortheap(nfill_half, WannierCenterKy_plus(:, ik2))


      !> sum over k1 to get wanniercenters for mirror minus
      Lambda0=0d0
      do i=1, nfill_half
         Lambda0(i, i)= 1d0
      enddo
      do ik1=1, Nk1
         !> <u_k|u_k+1>
         Mmnkb= 0d0
         hamk_dag= Eigenvector(:, :, ik1)
         if (ik1==Nk1) then
            hamk= Eigenvector(:, :, 1)
            ikp= 1
         else
            hamk= Eigenvector(:, :, ik1+ 1)
            ikp= ik1+ 1
         endif
         do m=1, Num_wann
            !ia= AtomIndex_orbital(m)
            !br= b(1)*Atom_position(1, ia)+ &
            !    b(2)*Atom_position(2, ia)+ &
            !    b(3)*Atom_position(3, ia)
            br= b(1)*Origin_cell%wannier_centers_cart(1, m )+ &
               b(2)*Origin_cell%wannier_centers_cart(2, m )+ &
               b(3)*Origin_cell%wannier_centers_cart(3, m )
            ratio= cos(br)- zi* sin(br)

            i1= 0
           !do j=1, nfill
            do j=Selected_Occupiedband_index(1), Selected_Occupiedband_index(NumberofSelectedOccupiedBands)
               if (mirror_plus(j, ikp)) cycle
               i1= i1+ 1
               i2= 0
              !do i=1, nfill
               do i=Selected_Occupiedband_index(1), Selected_Occupiedband_index(NumberofSelectedOccupiedBands)
                  if (mirror_plus(i, ik1)) cycle
                  i2= i2+ 1
                  Mmnkb(i2, i1)=  Mmnkb(i2, i1)+ &
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
         WannierCenterKy_minus(i, ik2)= aimag(log(Lambda_eig(i)))/2d0/pi
         WannierCenterKy_minus(i, ik2)= mod(WannierCenterKy_minus(i, ik2)+10d0, 1d0)
      enddo

      call sortheap(nfill_half, WannierCenterKy_minus(:, ik2))

   enddo !< ik2

   WannierCenterKy_minus_mpi= 0d0
   WannierCenterKy_plus_mpi= 0d0
#if defined (MPI)
   call mpi_allreduce(WannierCenterKy_minus, WannierCenterKy_minus_mpi, &
      size(WannierCenterKy_minus), mpi_dp, mpi_sum, mpi_cmw, ierr)
   call mpi_allreduce(WannierCenterKy_plus, WannierCenterKy_plus_mpi, &
      size(WannierCenterKy_plus), mpi_dp, mpi_sum, mpi_cmw, ierr)
#else
   WannierCenterKy_minus_mpi= WannierCenterKy_minus
   WannierCenterKy_plus_mpi= WannierCenterKy_plus
#endif

   outfileindex= outfileindex+ 1
   if (cpuid==0) then
      open(unit=outfileindex, file='wcc-mirrorplus.dat')

      write(outfileindex, '(10000A16)')'#      k',  'sum(wcc(:,ik))', &
         'wcc(:, ik)'
      do ik2=1, Nk2
         write(outfileindex, '(10000f16.8)') dble(ik2-1)/dble(Nk2-1), &
            dmod(sum(WannierCenterKy_plus_mpi(:, ik2)), 1d0), &
            WannierCenterKy_plus_mpi(:, ik2)
      enddo
      close(outfileindex)
   endif

   outfileindex= outfileindex+ 1
   if (cpuid==0) then
      open(unit=outfileindex, file='wcc-mirrorminus.dat')

      write(outfileindex, '(10000A16)')'#      k',  'sum(wcc(:,ik))', &
         'wcc(:, ik)'
      do ik2=1, Nk2
         write(outfileindex, '(10000f16.8)') dble(ik2-1)/dble(Nk2-1), &
            dmod(sum(WannierCenterKy_minus_mpi(:, ik2)), 1d0), &
            WannierCenterKy_minus_mpi(:, ik2)
      enddo
      close(outfileindex)
   endif

   !> generate gnu script for wannier charge center plots
   outfileindex= outfileindex+ 1
   if (cpuid==0) then
      open(unit=outfileindex, file='wcc-mirrorchernnumber.gnu')
      write(outfileindex, '(a)')"set encoding iso_8859_1"
      write(outfileindex, '(a)')'set terminal  postscript enhanced color font ",30"'
      write(outfileindex, '(a)')"set output 'wcc-mirrorchernnumber.eps'"
      write(outfileindex, '(a)')'set key '
      write(outfileindex, '(a)')'set border lw 3 '
      write(outfileindex, '(a)')'set xtics offset 0, 0.2'
      write(outfileindex, '(a)')'set xtics format "%4.1f" nomirror out '
      write(outfileindex, '(a)')'set xlabel "k" '
      write(outfileindex, '(a)')'set xlabel offset 0, 0.7 '
      write(outfileindex, '(a)')'set ytics 0.5 '
      write(outfileindex, '(a)')'set ytics format "%4.1f" nomirror out'
      write(outfileindex, '(a)')'set title "Mirror WCC"'
      write(outfileindex, '(a)')'#set ylabel offset 2, 0.0 '
      write(outfileindex, '(a)')'set ylabel "WCC"'
      write(outfileindex, '(a)')'set xrange [0: 1]'
      write(outfileindex, '(a)')'set yrange [0:1]'
      write(outfileindex, '(a)')"plot 'wcc-mirrorminus.dat' u 1:2 w p pt 7 ps 2  lc 'blue' title 'M=-i', \"
      write(outfileindex, '(a)')"     'wcc-mirrorplus.dat'  u 1:2 w p pt 7 ps 2  lc 'red'  title 'M=+i'"
      write(outfileindex, '(a)')'#unset key '
      write(outfileindex, '(a)')"#plot  \"
      write(outfileindex, '(a, i5, a)')"# for [i=3: ", NumberofSelectedOccupiedBands/2+2, &
         "]  'wcc-mirrorplus.dat' u 1:i w p  pt 7  ps 1.1 lc 'red', \"
      write(outfileindex, '(a, i5, a)')"# for [i=3: ", NumberofSelectedOccupiedBands/2+2, &
         "]  'wcc-mirrorminus.dat' u 1:i w p  pt 7  ps 1.1 lc 'blue'"
      close(outfileindex)
   endif

   wcc_sum= dmod(sum(WannierCenterKy_plus_mpi, dim=1), 1d0)
   !> determine the chirality
   gap_sum= 0d0
   do ik2=1, Nk2-1
      gap_step= wcc_sum(ik2+1)- wcc_sum(ik2)
      if (abs(gap_step+1)<abs(gap_step)) then
         gap_step= gap_step+ 1
      elseif (abs(gap_step-1)<abs(gap_step))then
         gap_step= gap_step- 1
      endif
      gap_sum= gap_sum+ gap_step
   enddo

   if (cpuid==0) write(stdout, '(1X, a, i5)')'MCN for ky=0 mirror +i : ', nint(gap_sum)


   wcc_sum= dmod(sum(WannierCenterKy_minus_mpi, dim=1), 1d0)
   !> determine the chirality
   gap_sum= 0d0
   do ik2=1, Nk2-1
      gap_step= wcc_sum(ik2+1)- wcc_sum(ik2)
      if (abs(gap_step+1)<abs(gap_step)) then
         gap_step= gap_step+ 1
      elseif (abs(gap_step-1)<abs(gap_step))then
         gap_step= gap_step- 1
      endif
      gap_sum= gap_sum+ gap_step
   enddo

   if (cpuid==0) write(stdout, '(1X, a, i5)')'MCN for ky=0 mirror -i : ', nint(gap_sum)

   return
end subroutine  wannier_center3D_plane_mirror


! !> this suboutine is used for wannier center calculation for 3D system
! !> only for one plane, calculate valley chern number. only choose the bands
! !> which have the same valley eigenvalue
! !> when calling this, please make sure that the valley operator matrix is
! !> properly set in the symmetry.f90.
! !> You can check it with ek_bulk_valley_z subroutine in ek_bulk.f90
!subroutine  wannier_center3D_plane_valley
!   use para
!   use wmpi
!   implicit none
!
!
!   integer :: i, j, l, m, ia, nfill, nfill_half
!
!   integer :: i1, i2, ik1, ik2, ikp, ierr
!
!      !> for sparse hamiltonian
!      !> dim= Num_wann*Num_wann
!      integer :: nnzmax, nnz
!      integer, allocatable :: jcoo(:), icoo(:)
!      complex(dp), allocatable :: acoo(:)
!      complex(dp), allocatable :: zeigv(:, :), zeigv_valley(:, :)
!
!      !number of ARPACK eigenvalues
!      integer :: neval
!   
!      ! number of Arnoldi vectors
!      integer :: nvecs
!   
!      !> calculate eigenvector or not
!      logical :: ritzvec
!   
!      !shift-invert shiftsigma
!      complex(dp) :: shiftsigma=(0d0,0d0)
!
!
!   !> k points in kx-ky plane
!   real(dp), allocatable :: kpoints(:, :, :)
!
!   !> hamiltonian for each k point
!   !> and also the eigenvector of hamiltonian after eigensystem_c
!   complex(dp), allocatable :: Hamk(:, :), hamk_dag(:, :)
!
!   !> eigenvector for each kx
!   complex(dp), allocatable :: Eigenvector(:, :, :)
!
!   !> Mmnkb=<u_n(k)|u_m(k+b)>
!   !> |u_n(k)> is the periodic part of wave function
!   complex(dp), allocatable :: Mmnkb(:, :), Mmnkb_com(:, :)
!
!   complex(dp), allocatable :: mat1(:, :), mat2(:, :)
!
!   complex(dp), allocatable :: Lambda_eig(:), Lambda(:, :), Lambda0(:, :)
!
!   !> three matrix for SVD
!   !> M= U.Sigma.V^\dag
!   !> VT= V^\dag
!   complex(dp), allocatable :: U(:, :), VT(:, :)
!   real   (dp), allocatable :: Sigma(:, :)
!
!   !> wannier centers for each ky, bands
!   real(dp), allocatable :: WannierCenterKy_minus(:, :),WannierCenterKy_minus_mpi(:, :)
!   real(dp), allocatable :: WannierCenterKy_plus(:, :),WannierCenterKy_plus_mpi(:, :)
!
!   !> sumation for Wannier charge center
!   real(dp) :: gap_sum, gap_step
!   real(dp), allocatable :: wcc_sum(:)
!
!   !> eigenvalue
!   real(dp), allocatable :: eigenvalue(:)
!
!   !> for each orbital, it correspond to an atom
!   !> dim= Num_wann
!   integer, allocatable :: AtomIndex_orbital(:)
!
!   real(dp) :: Umatrix_t(3, 3)
!
!   !> b.R
!   real(dp) :: br
!
!   !> exp(-i*b.R)
!   complex(dp) :: ratio
!
!   real(dp) :: k(3), b(3)
!
!   real(dp), allocatable :: xnm(:)
!   real(dp) :: k0(3), k1(3), k2(3)
!
!   !> valley eigenvalue
!   complex(dp), allocatable :: valley_z_eig(:, :)
!
!   !> the band index that has plus valley number
!   logical, allocatable :: valley_plus(:, :), valley_minus(:, :)
!
!   nfill= NumberofSelectedOccupiedBands
!   nfill_half= NumberofSelectedOccupiedBands/2
!
!   allocate(kpoints(3, Nk1, Nk2))
!   kpoints= 0d0
!
!   allocate(Lambda_eig(nfill_half))
!   allocate(Lambda(nfill_half, nfill_half))
!   allocate(Lambda0(nfill_half, nfill_half))
!   allocate(Mmnkb(nfill_half, nfill_half))
!   allocate(Mmnkb_com(nfill_half, nfill_half))
!   allocate(mat1(Num_wann, Num_wann))
!   allocate(mat2(Num_wann, Num_wann))
!   if (Is_Sparse) then
!      !> here we only store the wave function that used for calculating the Wilson loop
!      allocate(Eigenvector(Num_wann, nfill, Nk1))
!   else
!      allocate(Eigenvector(Num_wann, Num_wann, Nk1))
!   endif
!   allocate(eigenvalue(Num_wann))
!   allocate(valley_z_eig(nfill, Nk1))
!   allocate(valley_plus(nfill, Nk1))
!   allocate(valley_minus(nfill, Nk1))
!   allocate(U(nfill_half, nfill_half))
!   allocate(Sigma(nfill_half, nfill_half))
!   allocate(VT(nfill_half, nfill_half))
!   allocate(WannierCenterKy_minus(nfill_half, Nk2))
!   allocate(WannierCenterKy_minus_mpi(nfill_half, Nk2))
!   allocate(WannierCenterKy_plus(nfill_half, Nk2))
!   allocate(WannierCenterKy_plus_mpi(nfill_half, Nk2))
!   allocate(wcc_sum(Nk2))
!   allocate(AtomIndex_orbital(Num_wann))
!   allocate(xnm(nfill_half))
!   valley_minus= .False.
!   valley_plus= .False.
!   WannierCenterKy_minus= 0d0
!   WannierCenterKy_minus_mpi= 0d0
!   WannierCenterKy_plus= 0d0
!   WannierCenterKy_plus_mpi= 0d0
!   eigenvalue=0d0
!   Eigenvector=0d0
!   Mmnkb=0d0
!   Mmnkb_com=0d0
!   Lambda =0d0
!   Lambda0=0d0
!   U= 0d0
!   Sigma= 0d0
!   VT= 0d0
!   if (Is_Sparse) then
!      allocate(hamk(Num_wann, NumberofSelectedOccupiedBands))
!      allocate(hamk_dag(Num_wann, NumberofSelectedOccupiedBands))
!   else
!      allocate(hamk(Num_wann, Num_wann))
!      allocate(hamk_dag(Num_wann, Num_wann))
!   endif
!
!
!   if (Is_Sparse) then
!      neval=OmegaNum
!      if (neval>Num_wann-2) neval= Num_wann- 2
!   
!      !> ncv
!      nvecs=int(2*neval)
!      if (nvecs<50) nvecs= int(6*neval)
!   
!      if (nvecs>Num_wann) nvecs= Num_wann
!
!      shiftsigma=(1d0,0d0)*iso_energy
!      nnzmax=splen+Num_wann
!      nnz=splen
!      allocate( acoo(nnzmax))
!      allocate( jcoo(nnzmax))
!      allocate( icoo(nnzmax))
!      allocate( zeigv(Num_wann,nvecs))
!      allocate( zeigv_valley(Num_wann,nvecs))
!   endif
!
!
!   !> set k plane
!   !> the first dimension should be in one primitive cell, [0, 2*pi]
!   k0= K3D_start
!   k1= K3D_vec1
!   k2= K3D_vec2
!
!   do ik2=1, Nk2
!      do ik1=1, Nk1
!         kpoints(:, ik1, ik2)= k0+k1*dble(ik1-1.d0)/dble(nk1)+ k2*(ik2-1)/dble(nk2-1d0)
!      enddo
!   enddo
!   b= dble(k1)/dble(nk1)
!   b= b(1)*Origin_cell%kua+b(2)*Origin_cell%kub+b(3)*Origin_cell%kuc
!
!
!   !> set up atom index for each orbitals in the basis
!   if (soc>0) then  !> with spin orbital coupling
!      l= 0
!      do ia=1, Origin_cell%Num_atoms  !> spin up
!         do j=1, Origin_cell%nprojs(ia)
!            l= l+ 1
!            AtomIndex_orbital(l)= ia
!         enddo ! l
!      enddo ! ia
!      do ia=1, Origin_cell%Num_atoms  !> spin down
!         do j=1, Origin_cell%nprojs(ia)
!            l= l+ 1
!            AtomIndex_orbital(l)= ia
!         enddo ! l
!      enddo ! ia
!   else  !> without spin orbital coupling
!      l= 0
!      do ia=1, Origin_cell%Num_atoms  !> spin down
!         do j=1, Origin_cell%nprojs(ia)
!            l= l+ 1
!            AtomIndex_orbital(l)= ia
!         enddo ! l
!      enddo ! ia
!
!   endif
!
!   Umatrix_t= transpose(Umatrix)
!   call inv_r(3, Umatrix_t)
!
!   !>> Get wannier center for ky=0 plane
!   !> for each ky, we can get wanniercenter
!   do ik2=1+ cpuid, Nk2, num_cpu
!      if (cpuid.eq.0) write(stdout, *) 'ik', ik2
!
!      valley_plus= .False.
!      valley_minus= .False.
!      !> for each k1, we get the eigenvectors
!      do ik1=1, Nk1
!         k= kpoints(:, ik1, ik2)
!
!         if (Is_Sparse) then
!            call ham_bulk_coo_sparsehr_latticegauge(k,acoo,icoo,jcoo)
!            nnz= splen
!           
!            ritzvec= .true.
!            call arpack_sparse_coo_eigs(Num_wann,nnzmax,nnz,acoo,jcoo,icoo,neval,nvecs,eigenvalue,shiftsigma, zeigv, ritzvec)
!         else
!            !> get the TB hamiltonian in k space
!            call ham_bulk_latticegauge(k,hamk)
!            !> diagonal hamk
!            call eigensystem_c('V', 'U', Num_wann, hamk, eigenvalue)
!              !hamk(1:Num_wann, Selected_Occupiedband_index(1):Selected_Occupiedband_index(NumberofSelectedOccupiedBands))
!         endif
! 
!         !> get valley_plus and valley_minus
!         !> get valley operator in coo format
!         call valley_k_coo_sparsehr(nnzmax_valley, k, acoo_valley, icoo_valley, jcoo_valley)
!   
!         !> check the degeneracy of each band
!         IE=1
!         do while (ie.le.neval)
!            ND=1
!            if (ie+ND.le.neval) then
!               do while ((ie+ND).le.neval .and. (W(ie+ND)- W(ie)).lt.tolde)
!                  if (ie+ND .ge. neval) exit
!                  ND= ND+1
!               enddo
!            endif
!   
!   
!            !> if the degeneracy is larger than 1, we need to calculate the matrix
!            valley_k_nd= 0d0
!            do ie1= 1, ND
!               psi1= zeigv(:, IE+ie1-1)
!   
!               do ie2= 1, ND
!                  psi2= zeigv(:, IE+ie2-1)
!                  vpsi=0d0
!                  call mkl_zcoogemv('N', Num_wann, acoo_valley, icoo_valley, jcoo_valley, nnzmax_valley, psi2, vpsi)
!                  valley_k_nd(ie1, ie2)= zdotc(Num_wann, psi1, 1, vpsi, 1)
!               enddo
!            enddo
!   
!            !> diagonalize the matrix
!            VL = 0d0; VR= 0d0
!            if (ND>1) then
!               call zgeev_sys(ND, valley_k_nd(1:ND, 1:ND), valley_eig(1:ND),'N',VL(1:ND, 1:ND),"V",VR(1:ND, 1:ND) )
!               do ie1= 1, ND
!                  weight_valley(IE+ie1-1, ik) = real(valley_eig(ie1))
!               enddo
!
!               !> get new eigenvector
!               do ie1= 1, ND
!                  psi= 0d0
!                  do ie2= 1, ND
!                     psi= psi+ VR(ie2,  ie1)* zeigv(:, IE+ie2-1)
!                  enddo
!                  zeigv_valley(:, ie+ ie1- 1)= psi 
!               enddo
!            else
!               valley_eig(1)= valley_k_nd(1, 1)
!               weight_valley(ie) = real(valley_k_nd(1, 1))
!               zeigv_valley(:, ie)= zeigv_valley(:, ie)
!            endif
!   
!            IE= IE+ ND
!         enddo
!
!         do ie=1, nfill
!            if (weight_valley(ie)>0) valley_plus(ie, ik)=.true.
!         enddo
!
!      enddo ! ik1
!
!      !> sum over k1 to get wanniercenters for valley plus
!      Lambda0=0d0
!      do i=1, nfill_half
!         Lambda0(i, i)= 1d0
!      enddo
!      do ik1=1, Nk1
!         !> <u_k|u_k+1>
!         Mmnkb= 0d0
!         hamk_dag= Eigenvector(:, :, ik1)
!         if (ik1==Nk1) then
!            hamk= Eigenvector(:, :, 1)
!            ikp= 1
!         else
!            hamk= Eigenvector(:, :, ik1+ 1)
!            ikp= ik1+ 1
!         endif
!         do m=1, Num_wann
!            !ia= AtomIndex_orbital(m)
!            !br= b(1)*Atom_position(1, ia)+ &
!            !    b(2)*Atom_position(2, ia)+ &
!            !    b(3)*Atom_position(3, ia)
!            br= b(1)*Origin_cell%wannier_centers_cart(1, m )+ &
!               b(2)*Origin_cell%wannier_centers_cart(2, m )+ &
!               b(3)*Origin_cell%wannier_centers_cart(3, m )
!            ratio= cos(br)- zi* sin(br)
!
!            i1= 0
!           !do j=1, nfill
!            do j=Selected_Occupiedband_index(1), Selected_Occupiedband_index(NumberofSelectedOccupiedBands)
!               if (valley_minus(j, ikp)) cycle
!               i1= i1+ 1
!               i2= 0
!              !do i=1, nfill
!               do i=Selected_Occupiedband_index(1), Selected_Occupiedband_index(NumberofSelectedOccupiedBands)
!                  if (valley_minus(i, ik1)) cycle
!                  i2= i2+ 1
!                  Mmnkb(i2, i1)=  Mmnkb(i2, i1)+ &
!                     conjg(hamk_dag(m, i))* hamk(m, j)* ratio
!               enddo ! i
!            enddo ! j
!         enddo ! m
!
!         !> perform Singluar Value Decomposed of Mmnkb
!         call zgesvd_pack(nfill_half, Mmnkb, U, Sigma, VT)
!
!         !> after the calling of zgesvd_pack, Mmnkb becomes a temporal matrix
!         U= conjg(transpose(U))
!         VT= conjg(transpose(VT))
!         call mat_mul(nfill_half, VT, U, Mmnkb)
!
!         call mat_mul(nfill_half, Mmnkb, Lambda0, Lambda)
!         Lambda0 = Lambda
!      enddo  !< ik1
!
!      !> diagonalize Lambda to get the eigenvalue
!      call zgeev_pack(nfill_half, Lambda, Lambda_eig)
!      do i=1, nfill_half
!         WannierCenterKy_plus(i, ik2)= aimag(log(Lambda_eig(i)))/2d0/pi
!         WannierCenterKy_plus(i, ik2)= mod(WannierCenterKy_plus(i, ik2)+10d0, 1d0)
!      enddo
!
!      call sortheap(nfill_half, WannierCenterKy_plus(:, ik2))
!
!
!      !> sum over k1 to get wanniercenters for mirror minus
!      Lambda0=0d0
!      do i=1, nfill_half
!         Lambda0(i, i)= 1d0
!      enddo
!      do ik1=1, Nk1
!         !> <u_k|u_k+1>
!         Mmnkb= 0d0
!         hamk_dag= Eigenvector(:, :, ik1)
!         if (ik1==Nk1) then
!            hamk= Eigenvector(:, :, 1)
!            ikp= 1
!         else
!            hamk= Eigenvector(:, :, ik1+ 1)
!            ikp= ik1+ 1
!         endif
!         do m=1, Num_wann
!            !ia= AtomIndex_orbital(m)
!            !br= b(1)*Atom_position(1, ia)+ &
!            !    b(2)*Atom_position(2, ia)+ &
!            !    b(3)*Atom_position(3, ia)
!            br= b(1)*Origin_cell%wannier_centers_cart(1, m )+ &
!               b(2)*Origin_cell%wannier_centers_cart(2, m )+ &
!               b(3)*Origin_cell%wannier_centers_cart(3, m )
!            ratio= cos(br)- zi* sin(br)
!
!            i1= 0
!           !do j=1, nfill
!            do j=Selected_Occupiedband_index(1), Selected_Occupiedband_index(NumberofSelectedOccupiedBands)
!               if (mirror_plus(j, ikp)) cycle
!               i1= i1+ 1
!               i2= 0
!              !do i=1, nfill
!               do i=Selected_Occupiedband_index(1), Selected_Occupiedband_index(NumberofSelectedOccupiedBands)
!                  if (mirror_plus(i, ik1)) cycle
!                  i2= i2+ 1
!                  Mmnkb(i2, i1)=  Mmnkb(i2, i1)+ &
!                     conjg(hamk_dag(m, i))* hamk(m, j)* ratio
!               enddo ! i
!            enddo ! j
!         enddo ! m
!
!         !> perform Singluar Value Decomposed of Mmnkb
!         call zgesvd_pack(nfill_half, Mmnkb, U, Sigma, VT)
!
!         !> after the calling of zgesvd_pack, Mmnkb becomes a temporal matrix
!         U= conjg(transpose(U))
!         VT= conjg(transpose(VT))
!         call mat_mul(nfill_half, VT, U, Mmnkb)
!
!         call mat_mul(nfill_half, Mmnkb, Lambda0, Lambda)
!         Lambda0 = Lambda
!      enddo  !< ik1
!
!      !> diagonalize Lambda to get the eigenvalue
!      call zgeev_pack(nfill_half, Lambda, Lambda_eig)
!      do i=1, nfill_half
!         WannierCenterKy_minus(i, ik2)= aimag(log(Lambda_eig(i)))/2d0/pi
!         WannierCenterKy_minus(i, ik2)= mod(WannierCenterKy_minus(i, ik2)+10d0, 1d0)
!      enddo
!
!      call sortheap(nfill_half, WannierCenterKy_minus(:, ik2))
!
!   enddo !< ik2
!
!   WannierCenterKy_minus_mpi= 0d0
!   WannierCenterKy_plus_mpi= 0d0
!#if defined (MPI)
!   call mpi_allreduce(WannierCenterKy_minus, WannierCenterKy_minus_mpi, &
!      size(WannierCenterKy_minus), mpi_dp, mpi_sum, mpi_cmw, ierr)
!   call mpi_allreduce(WannierCenterKy_plus, WannierCenterKy_plus_mpi, &
!      size(WannierCenterKy_plus), mpi_dp, mpi_sum, mpi_cmw, ierr)
!#else
!   WannierCenterKy_minus_mpi= WannierCenterKy_minus
!   WannierCenterKy_plus_mpi= WannierCenterKy_plus
!#endif
!
!   outfileindex= outfileindex+ 1
!   if (cpuid==0) then
!      open(unit=outfileindex, file='wcc-mirrorplus.dat')
!
!      write(outfileindex, '(10000A16)')'#      k',  'sum(wcc(:,ik))', &
!         'wcc(:, ik)'
!      do ik2=1, Nk2
!         write(outfileindex, '(10000f16.8)') dble(ik2-1)/dble(Nk2-1), &
!            dmod(sum(WannierCenterKy_plus_mpi(:, ik2)), 1d0), &
!            WannierCenterKy_plus_mpi(:, ik2)
!      enddo
!      close(outfileindex)
!   endif
!
!   outfileindex= outfileindex+ 1
!   if (cpuid==0) then
!      open(unit=outfileindex, file='wcc-mirrorminus.dat')
!
!      write(outfileindex, '(10000A16)')'#      k',  'sum(wcc(:,ik))', &
!         'wcc(:, ik)'
!      do ik2=1, Nk2
!         write(outfileindex, '(10000f16.8)') dble(ik2-1)/dble(Nk2-1), &
!            dmod(sum(WannierCenterKy_minus_mpi(:, ik2)), 1d0), &
!            WannierCenterKy_minus_mpi(:, ik2)
!      enddo
!      close(outfileindex)
!   endif
!
!   !> generate gnu script for wannier charge center plots
!   outfileindex= outfileindex+ 1
!   if (cpuid==0) then
!      open(unit=outfileindex, file='wcc-mirrorchernnumber.gnu')
!      write(outfileindex, '(a)')"set encoding iso_8859_1"
!      write(outfileindex, '(a)')'set terminal  postscript enhanced color font ",30"'
!      write(outfileindex, '(a)')"set output 'wcc-mirrorchernnumber.eps'"
!      write(outfileindex, '(a)')'set key '
!      write(outfileindex, '(a)')'set border lw 3 '
!      write(outfileindex, '(a)')'set xtics offset 0, 0.2'
!      write(outfileindex, '(a)')'set xtics format "%4.1f" nomirror out '
!      write(outfileindex, '(a)')'set xlabel "k" '
!      write(outfileindex, '(a)')'set xlabel offset 0, 0.7 '
!      write(outfileindex, '(a)')'set ytics 0.5 '
!      write(outfileindex, '(a)')'set ytics format "%4.1f" nomirror out'
!      write(outfileindex, '(a)')'set title "Mirror WCC"'
!      write(outfileindex, '(a)')'#set ylabel offset 2, 0.0 '
!      write(outfileindex, '(a)')'set ylabel "WCC"'
!      write(outfileindex, '(a)')'set xrange [0: 1]'
!      write(outfileindex, '(a)')'set yrange [0:1]'
!      write(outfileindex, '(a)')"plot 'wcc-mirrorminus.dat' u 1:2 w p pt 7 ps 2  lc 'blue' title 'M=-i', \"
!      write(outfileindex, '(a)')"     'wcc-mirrorplus.dat'  u 1:2 w p pt 7 ps 2  lc 'red'  title 'M=+i'"
!      write(outfileindex, '(a)')'#unset key '
!      write(outfileindex, '(a)')"#plot  \"
!      write(outfileindex, '(a, i5, a)')"# for [i=3: ", NumberofSelectedOccupiedBands/2+2, &
!         "]  'wcc-mirrorplus.dat' u 1:i w p  pt 7  ps 1.1 lc 'red', \"
!      write(outfileindex, '(a, i5, a)')"# for [i=3: ", NumberofSelectedOccupiedBands/2+2, &
!         "]  'wcc-mirrorminus.dat' u 1:i w p  pt 7  ps 1.1 lc 'blue'"
!      close(outfileindex)
!   endif
!
!   wcc_sum= dmod(sum(WannierCenterKy_plus_mpi, dim=1), 1d0)
!   !> determine the chirality
!   gap_sum= 0d0
!   do ik2=1, Nk2-1
!      gap_step= wcc_sum(ik2+1)- wcc_sum(ik2)
!      if (abs(gap_step+1)<abs(gap_step)) then
!         gap_step= gap_step+ 1
!      elseif (abs(gap_step-1)<abs(gap_step))then
!         gap_step= gap_step- 1
!      endif
!      gap_sum= gap_sum+ gap_step
!   enddo
!
!   if (cpuid==0) write(stdout, '(1X, a, i5)')'MCN for ky=0 mirror +i : ', nint(gap_sum)
!
!
!   wcc_sum= dmod(sum(WannierCenterKy_minus_mpi, dim=1), 1d0)
!   !> determine the chirality
!   gap_sum= 0d0
!   do ik2=1, Nk2-1
!      gap_step= wcc_sum(ik2+1)- wcc_sum(ik2)
!      if (abs(gap_step+1)<abs(gap_step)) then
!         gap_step= gap_step+ 1
!      elseif (abs(gap_step-1)<abs(gap_step))then
!         gap_step= gap_step- 1
!      endif
!      gap_sum= gap_sum+ gap_step
!   enddo
!
!   if (cpuid==0) write(stdout, '(1X, a, i5)')'MCN for ky=0 mirror -i : ', nint(gap_sum)
!
!   return
!end subroutine  wannier_center3D_plane_valley


subroutine  wannier_center2D
   ! This suboutine is used for wannier center calculation for slab system
   !
   ! Copyright (c) 2010 QuanSheng Wu. All rights reserved.

   use para
   use wmpi
   implicit none

   integer :: Nkx
   integer :: Nky

   integer :: i
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

   nfill= NumberofSelectedOccupiedBands*Nslab

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
         Lambda0(i, i)= 1d0!lam0=I
      enddo

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
         Mmnkb= Mmnkb_full((Selected_Occupiedband_index(1)-1)*Nslab+1:(Selected_Occupiedband_index(1)-1)*Nslab+nfill, &
            (Selected_Occupiedband_index(1)-1)*Nslab+1:(Selected_Occupiedband_index(1)-1)*Nslab+nfill)
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
         WannierCenterKy(i, iky)=(mod(WannierCenterKy(i, iky)+10d0, 1d0))
      enddo

   enddo !< iky

#if defined (MPI)
   call mpi_allreduce(WannierCenterKy, WannierCenterKy_mpi, &
      size(WannierCenterKy), mpi_dp, mpi_sum, mpi_cmw, ierr)
#else
   WannierCenterKy_mpi= WannierCenterKy
#endif


   outfileindex= outfileindex+ 1
   if (cpuid==0) then
      open(unit=outfileindex, file='wcc_slab.dat')

      do iky=1, Nky
         write(outfileindex, '(10000f16.8)') kpoints(2, 1, iky), &
            dmod(sum(WannierCenterKy_mpi(:, iky)), 1d0), &
            WannierCenterKy_mpi(:, iky)
      enddo
      close(outfileindex)
   endif

   outfileindex= outfileindex+ 1
   if (cpuid==0) then
      open(unit=outfileindex, file='wcc_slab.gnu')

      write(outfileindex,*) 'set encoding iso_8859_1'
      write(outfileindex,*) 'set terminal  postscript enhanced color font ",30" '
      write(outfileindex,*) "set output 'wcc_slab.eps'"
      write(outfileindex,*) 'unset key '
      write(outfileindex,*) 'set border lw 3 '
      write(outfileindex,*) 'set xtics offset 0, 0.2'
      write(outfileindex,*) 'set xtics format "%4.1f" nomirror out '
      write(outfileindex,*) 'set xlabel "k" '
      write(outfileindex,*) 'set xlabel offset 0, 0.7 '
      write(outfileindex,*) 'set ytics 0.5 '
      write(outfileindex,*) 'set ytics format "%4.1f" nomirror out'
      write(outfileindex,*) 'set ylabel "WCC-slab"'
      write(outfileindex,*) 'set ylabel offset 2, 0.0 '
      write(outfileindex,*) 'set xrange [0: 1]'
      write(outfileindex,*) 'set yrange [-0.0:1.0]'
      write(outfileindex, '(a, i5, a)')"plot for [i=3: ", nfill+2, "] 'wcc_slab.dat' u 1:i w p  pt 7  ps 1.1 lc 'red'"
      close(outfileindex)
   endif

   return
end subroutine  wannier_center2D


subroutine  wannier_center2D_alt
   ! This suboutine is used for wannier center calculation for slab system
   ! calculate z2

   use para
   use wmpi
   implicit none

   integer :: Nkx
   integer :: Nky

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
   !> for slab system, dim=Nslab*Origin_cell%Num_atoms
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

   Nkx= Nk
   Nky= 40

   nfill= NumberofSelectedOccupiedBands*Nslab

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
   allocate(AtomsPosition_unitcell(3, Origin_cell%Num_atoms))
   allocate(AtomsPosition_supercell(3, Nslab*Origin_cell%Num_atoms))
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
         do ia=1, Origin_cell%Num_atoms  !> spin up
            do j=1, Origin_cell%nprojs(ia)
               l= l+ 1
               AtomIndex_orbital(l)= ia+ (i-1)*Origin_cell%Num_atoms
            enddo ! l
         enddo ! ia
         do ia=1, Origin_cell%Num_atoms  !> spin down
            do j=1, Origin_cell%nprojs(ia)
               l= l+ 1
               AtomIndex_orbital(l)= ia+ (i-1)*Origin_cell%Num_atoms
            enddo ! l
         enddo ! ia
      enddo ! i
   else  !> without spin orbital coupling
      l= 0
      do i=1, Nslab
         do ia=1, Origin_cell%Num_atoms  !> spin down
            do j=1, Origin_cell%nprojs(ia)
               l= l+ 1
               AtomIndex_orbital(l)= ia+ (i-1)*Origin_cell%Num_atoms
            enddo ! l
         enddo ! ia
      enddo ! i

   endif

   Umatrix_t= transpose(Umatrix)
   call inv_r(3, Umatrix_t)

   !> set up atoms' position in the unit cell in the new basis
   !> only for 2D slab system
   do ia=1, Origin_cell%Num_atoms
      do i=1, 3
         do j=1, 3
            AtomsPosition_unitcell(i, ia)= AtomsPosition_unitcell(i, ia)+ &
               Umatrix_t(i, j)*Origin_cell%Atom_position_cart(j, ia)
         enddo ! j
      enddo ! i
   enddo ! ia

   !> set up atoms' position in the supercell
   !> actually, we only need the first two coordinates
   ia1= 0
   do i=1, Nslab
      do ia=1, Origin_cell%Num_atoms
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

               do j=1, nfill
                  do i=1, nfill
                     Mmnkb(i, j)=  Mmnkb(i, j)+ &
                        conjg(hamk_dag((l-1)*Num_wann+m, (Selected_Occupiedband_index(1)-1)*Nslab+i))* &
                        hamk((l-1)*Num_wann+m, (Selected_Occupiedband_index(1)-1)*Nslab+j)* ratio
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

#if defined (MPI)
   call mpi_allreduce(WannierCenterKy, WannierCenterKy_mpi, &
      size(WannierCenterKy), mpi_dp, mpi_sum, mpi_cmw, ierr)
   call mpi_allreduce(largestgap, largestgap_mpi, &
      size(largestgap), mpi_dp, mpi_sum, mpi_cmw, ierr)
#else
   WannierCenterKy_mpi= WannierCenterKy
   largestgap_mpi= largestgap
#endif




   outfileindex= outfileindex+ 1
   if (cpuid==0) then
      open(unit=outfileindex, file='wanniercenter.dat')

      do iky=1, Nky
         write(outfileindex, '(10000f16.8)') kpoints(2, 1, iky), &
            dmod(sum(WannierCenterKy_mpi(:, iky)), 1d0), &
            WannierCenterKy_mpi(:, iky)
      enddo
      close(outfileindex)
   endif

   outfileindex= outfileindex+ 1
   if (cpuid==0) then
      open(unit=outfileindex, file='largestgap.dat')
      do iky=1, Nky
         write(outfileindex, '(10000f16.8)') kpoints(2, 1, iky), &
            largestgap_mpi(iky)
      enddo
      close(outfileindex)
   endif
   return
end subroutine  wannier_center2D_alt


subroutine  wannier_center3D_plane_mirror_minus
   !> this suboutine is used for wannier center calculation for 3D system
   !> only for one plane, calculate mirror chern number. only choose the bands
   !> which have the same mirror eigenvalue
   use para
   use wmpi
   implicit none


   integer :: i
   integer :: j
   integer :: l
   integer :: m
   integer :: ia
   integer :: nfill
   integer :: nfill_half
   integer :: imax

   integer :: i1

   integer :: i2
   integer :: ik1
   integer :: ik2
   integer :: ikp

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

   !> for each orbital, it correspond to an atom
   !> dim= Num_wann
   integer, allocatable :: AtomIndex_orbital(:)

   real(dp) :: Umatrix_t(3, 3)

   !> b.R
   real(dp) :: br

   !> exp(-i*b.R)
   complex(dp) :: ratio

   real(dp) :: k(3)
   real(dp) :: b(3)

   real(dp) :: maxgap
   real(dp) :: maxgap0

   !> Z2 calculation for time reversal invariant system
   real(dp) :: Z2
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
   real(dp) :: k0(3), k1(3), k2(3)

   !> mirror eigenvalue
   complex(dp), allocatable :: mirror_z_eig(:, :)
   !> the band index that has plus mirror number
   logical, allocatable :: mirror_plus(:, :)
   logical, allocatable :: mirror_minus(:, :)


   nfill= NumberofSelectedOccupiedBands
   nfill_half= NumberofSelectedOccupiedBands/2

   allocate(kpoints(3, Nk1, Nk2))
   kpoints= 0d0

   allocate(Lambda_eig(nfill_half))
   allocate(Lambda(nfill_half, nfill_half))
   allocate(Lambda0(nfill_half, nfill_half))
   allocate(Mmnkb(nfill_half, nfill_half))
   allocate(Mmnkb_com(nfill_half, nfill_half))
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
   mirror_minus= .False.
   mirror_plus= .False.
   largestgap= 0d0
   largestgap_mpi= 0d0
   WannierCenterKy= 0d0
   WannierCenterKy_mpi= 0d0
   hamk=0d0
   eigenvalue=0d0
   Eigenvector=0d0
   Mmnkb=0d0
   Mmnkb_com=0d0
   Lambda =0d0
   Lambda0=0d0
   U= 0d0
   Sigma= 0d0
   VT= 0d0

   !> set k plane
   !> the first dimension should be in one primitive cell, [0, 2*pi]
   k0= K3D_start !
   k1= K3D_vec1   !
   k2= K3D_vec2   !

   do ik2=1, Nk2
      do ik1=1, Nk1
         kpoints(:, ik1, ik2)= k0+k1*(ik1-1)/dble(nk1)+ k2*(ik2-1)/dble(nk2-1)
      enddo
   enddo
   b= k1/dble(nk1)
   b= b(1)*Origin_cell%Kua+b(2)*Origin_cell%Kub+b(3)*Origin_cell%Kuc


   !> set up atom index for each orbitals in the basis
   if (soc>0) then  !> with spin orbital coupling
      l= 0
      do ia=1, Origin_cell%Num_atoms  !> spin up
         do j=1, Origin_cell%nprojs(ia)
            l= l+ 1
            AtomIndex_orbital(l)= ia
         enddo ! l
      enddo ! ia
      do ia=1, Origin_cell%Num_atoms  !> spin down
         do j=1, Origin_cell%nprojs(ia)
            l= l+ 1
            AtomIndex_orbital(l)= ia
         enddo ! l
      enddo ! ia
   else  !> without spin orbital coupling
      l= 0
      do ia=1, Origin_cell%Num_atoms  !> spin down
         do j=1, Origin_cell%nprojs(ia)
            l= l+ 1
            AtomIndex_orbital(l)= ia
         enddo ! l
      enddo ! ia

   endif

   Umatrix_t= transpose(Umatrix)
   call inv_r(3, Umatrix_t)

   mirror_z= mirror_z
   !>> Get wannier center for ky=0 plane
   !> for each ky, we can get wanniercenter
   do ik2=1+ cpuid, Nk2, num_cpu
      Lambda0=0d0
      do i=1, nfill_half
         Lambda0(i, i)= 1d0
      enddo

      mirror_plus= .False.
      mirror_minus= .False.
      !> for each k1, we get the eigenvectors
      do ik1=1, Nk1
         k= kpoints(:, ik1, ik2)


         ! generate bulk Hamiltonian
         if (index(KPorTB, 'KP')/=0)then
            call ham_bulk_kp(k, Hamk)
         else
            !> deal with phonon system
            if (index(Particle,'phonon')/=0.and.LOTO_correction) then
               call ham_bulk_LOTO(k, Hamk)
            else
               call ham_bulk_latticegauge    (k, Hamk)
            endif
         endif

         !> symmetrization
         call mat_mul(Num_wann, mirror_z, hamk, mat1)
         call mat_mul(Num_wann, mat1, mirror_z, mat2)
         mat1= (hamk+ mat2)/2.d0
         hamk= mat1

         !> diagonal hamk
         call eigensystem_c('V', 'U', Num_wann, hamk, eigenvalue)

         Eigenvector(:, :, ik1)= hamk

         mat2= conjg(transpose(hamk))

         !> calculate mirror eigenvalue
         call mat_mul(Num_wann, mat2, mirror_z, mat1)
         call mat_mul(Num_wann, mat1, hamk, mat2)

         !> get mirror_plus and mirror_minus
         do i=1, nfill
         if (abs(real(mat2(Selected_Occupiedband_index(1)+i-1, Selected_Occupiedband_index(1)+i-1))-1d0)< 1e-3) then
               mirror_plus(i, ik1)= .true.
            else
               mirror_minus(i, ik1)= .true.
            endif
         enddo

      enddo

      !> sum over k1 to get wanniercenters
      do ik1=1, Nk1
         !> <u_k|u_k+1>
         Mmnkb= 0d0
         hamk_dag= Eigenvector(:, :, ik1)
         if (ik1==Nk1) then
            hamk= Eigenvector(:, :, 1)
            ikp= 1
         else
            hamk= Eigenvector(:, :, ik1+ 1)
            ikp= ik1+ 1
         endif
         do m=1, Num_wann
            !ia= AtomIndex_orbital(m)
            !br= b(1)*Origin_cell%Atom_position_cart(1, ia)+ &
            !    b(2)*Origin_cell%Atom_position_cart(2, ia)+ &
            !    b(3)*Origin_cell%Atom_position_cart(3, ia)
            br= b(1)*Origin_cell%wannier_centers_cart(1, m )+ &
                b(2)*Origin_cell%wannier_centers_cart(2, m )+ &
                b(3)*Origin_cell%wannier_centers_cart(3, m )
            ratio= cos(br)- zi* sin(br)
            !ratio= 1d0

            i1= 0
            do j=1, nfill
               if (mirror_minus(j, ikp)) cycle
               i1= i1+ 1
               i2= 0
               do i=1, nfill
                  if (mirror_minus(i, ik1)) cycle
                  i2= i2+ 1
                  Mmnkb(i2, i1)=  Mmnkb(i2, i1)+ &
                     conjg(hamk_dag(m, Selected_Occupiedband_index(1)-1+i))* &
                     hamk(m, Selected_Occupiedband_index(1)-1+j)* ratio
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
            maxgap= 1d0+ WannierCenterKy(1, ik2)- WannierCenterKy(nfill_half, ik2)
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
#if defined (MPI)
   call mpi_allreduce(WannierCenterKy, WannierCenterKy_mpi, &
      size(WannierCenterKy), mpi_dp, mpi_sum, mpi_cmw, ierr)
   call mpi_allreduce(largestgap, largestgap_mpi, &
      size(largestgap), mpi_dp, mpi_sum, mpi_cmw, ierr)
#else
   WannierCenterKy_mpi= WannierCenterKy
   largestgap_mpi= largestgap
#endif


   outfileindex= outfileindex+ 1
   if (cpuid==0) then
      open(unit=outfileindex, file='wcc-mirrorminus.dat')

      do ik2=1, Nk2
         write(outfileindex, '(10000f16.8)') dble(ik2-1)/dble(Nk2-1)/2d0, &
            largestgap_mpi(ik2), dmod(sum(WannierCenterKy_mpi(:, ik2)), 1d0), &
            WannierCenterKy_mpi(:, ik2)
      enddo
      close(outfileindex)
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

   if (cpuid==0) write(stdout, *)'Z2 for ky=0 mirror -i : ', Z2

   return
end subroutine  wannier_center3D_plane_mirror_minus


subroutine  wannier_center3D_plane_mirror_plus
   !> this suboutine is used for wannier center calculation for 3D system
   !> only for one plane, calculate mirror chern number. only choose the bands
   !> which have the same mirror eigenvalue
   use para
   use wmpi
   implicit none


   integer :: i, j, m , nfill, nfill_half, imax
   integer :: i1, i2, ik1, ik2, ikp
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

   complex(dp), allocatable :: mat1(:, :)
   complex(dp), allocatable :: mat2(:, :)

   !>
   complex(dp), allocatable :: Lambda_eig(:)
   complex(dp), allocatable :: Lambda(:, :)
   complex(dp), allocatable :: Lambda0(:, :)


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

   real(dp) :: Umatrix_t(3, 3)

   !> b.R
   real(dp) :: br

   !> exp(-i*b.R)
   complex(dp) :: ratio

   real(dp) :: k(3)
   real(dp) :: b(3)

   real(dp) :: maxgap
   real(dp) :: maxgap0

   !> Z2 calculation for time reversal invariant system
   real(dp) :: Z2
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
   real(dp) :: k0(3), k1(3), k2(3)

   !> mirror eigenvalue
   complex(dp), allocatable :: mirror_z_eig(:, :)
   !> the band index that has plus mirror number
   logical, allocatable :: mirror_plus(:, :)
   logical, allocatable :: mirror_minus(:, :)


   nfill= NumberofSelectedOccupiedBands
   !nfill_half= NumberofSelectedOccupiedBands
   nfill_half= NumberofSelectedOccupiedBands/2

   allocate(kpoints(3, Nk1, Nk2))
   kpoints= 0d0

   allocate(Lambda_eig(nfill_half))
   allocate(Lambda(nfill_half, nfill_half))
   allocate(Lambda0(nfill_half, nfill_half))
   allocate(Mmnkb(nfill_half, nfill_half))
   allocate(Mmnkb_com(nfill_half, nfill_half))
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
   allocate(xnm(nfill_half))
   allocate(largestgap(Nk2))
   allocate(largestgap_mpi(Nk2))
   mirror_minus= .False.
   mirror_plus= .False.
   largestgap= 0d0
   largestgap_mpi= 0d0
   WannierCenterKy= 0d0
   WannierCenterKy_mpi= 0d0
   hamk=0d0
   eigenvalue=0d0
   Eigenvector=0d0
   Mmnkb=0d0
   Mmnkb_com=0d0
   Lambda =0d0
   Lambda0=0d0
   U= 0d0
   Sigma= 0d0
   VT= 0d0

   !> set k plane
   !> the first dimension should be in one primitive cell, [0, 2*pi]
   k0= K3D_start !
   k1= K3D_vec1   !
   k2= K3D_vec2   !

   do ik2=1, Nk2
      do ik1=1, Nk1
         kpoints(:, ik1, ik2)= k0+k1*(ik1-1)/dble(nk1)+ k2*(ik2-1)/dble(nk2-1)
      enddo
   enddo
   b= k1/dble(nk1)
   b= b(1)*Origin_cell%Kua+b(2)*Origin_cell%Kub+b(3)*Origin_cell%Kuc

   Umatrix_t= transpose(Umatrix)
   call inv_r(3, Umatrix_t)
   mirror_z= mirror_z

   !>> Get wannier center for ky=0 plane
   !> for each ky, we can get wanniercenter
   do ik2=1+ cpuid, Nk2, num_cpu
      Lambda0=0d0
      do i=1, nfill_half
         Lambda0(i, i)= 1d0
      enddo

      mirror_plus= .False.
      mirror_minus= .False.
      !> for each k1, we get the eigenvectors
      do ik1=1, Nk1
         k= kpoints(:, ik1, ik2)


         ! generate bulk Hamiltonian
         if (index(KPorTB, 'KP')/=0)then
            call ham_bulk_kp(k, Hamk)
         else
            !> deal with phonon system
            if (index(Particle,'phonon')/=0.and.LOTO_correction) then
               call ham_bulk_LOTO(k, Hamk)
            else
               call ham_bulk_latticegauge    (k, Hamk)
            endif
         endif


         !> symmetrization
         call mat_mul(Num_wann, mirror_z, hamk, mat1)
         call mat_mul(Num_wann, mat1, mirror_z, mat2)
         mat1= (hamk+ mat2)/2.d0
         hamk= mat1

         !> diagonal hamk
         call eigensystem_c('V', 'U', Num_wann, hamk, eigenvalue)

         Eigenvector(:, :, ik1)= hamk

         mat2= conjg(transpose(hamk))

         !> calculate mirror eigenvalue
         call mat_mul(Num_wann, mat2, mirror_z, mat1)
         call mat_mul(Num_wann, mat1, hamk, mat2)

         !> get mirror_plus and mirror_minus
         do i=1, nfill
         if (abs(real(mat2(Selected_Occupiedband_index(1)-1+i, Selected_Occupiedband_index(1)-1+i))-1d0)< 1e-3) then
               mirror_plus(i, ik1)= .true.
            else
               mirror_minus(i, ik1)= .true.
            endif
         enddo

      enddo

      !> sum over k1 to get wanniercenters
      do ik1=1, Nk1
         !> <u_k|u_k+1>
         Mmnkb= 0d0
         hamk_dag= Eigenvector(:, :, ik1)
         if (ik1==Nk1) then
            hamk= Eigenvector(:, :, 1)
            ikp= 1
         else
            hamk= Eigenvector(:, :, ik1+ 1)
            ikp= ik1+ 1
         endif
         do m=1, Num_wann
            br= b(1)*Origin_cell%wannier_centers_cart(1, m )+ &
               b(2)*Origin_cell%wannier_centers_cart(2, m )+ &
               b(3)*Origin_cell%wannier_centers_cart(3, m )
            ratio= cos(br)- zi* sin(br)

            i1= 0
            do j=1, nfill
               if (mirror_plus(j, ikp)) cycle
               i1= i1+ 1
               i2= 0
               do i=1, nfill
                  if (mirror_plus(i, ik1)) cycle
                  i2= i2+ 1
                  Mmnkb(i2, i1)=  Mmnkb(i2, i1)+ &
                     conjg(hamk_dag(m, Selected_Occupiedband_index(1)-1+i))* hamk(m, Selected_Occupiedband_index(1)-1+j)* ratio
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
#if defined (MPI)
   call mpi_allreduce(WannierCenterKy, WannierCenterKy_mpi, &
      size(WannierCenterKy), mpi_dp, mpi_sum, mpi_cmw, ierr)
   call mpi_allreduce(largestgap, largestgap_mpi, &
      size(largestgap), mpi_dp, mpi_sum, mpi_cmw, ierr)
#else
   WannierCenterKy_mpi= WannierCenterKy
   largestgap_mpi= largestgap
#endif


   outfileindex= outfileindex+ 1
   if (cpuid==0) then
      open(unit=outfileindex, file='wcc-mirrorplus.dat')

      do ik2=1, Nk2
         write(outfileindex, '(10000f16.8)') dble(ik2-1)/dble(Nk2-1)/2, &
            largestgap_mpi(ik2), dmod(sum(WannierCenterKy_mpi(:, ik2)), 1d0), &
            WannierCenterKy_mpi(:, ik2)
      enddo
      close(outfileindex)
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

   if (cpuid==0) write(stdout, *)'Z2 for ky=0: ', Z2

   return
end subroutine  wannier_center3D_plane_mirror_plus


subroutine  wannier_center3D_plane0
   !> this suboutine is used for wannier center calculation for 3D system
   !> only for one plane

   use para
   use wmpi
   implicit none
   integer :: ik2

   real(dp), allocatable :: wcc(:, :)

   !> larges gap of each two wannier centers for a given k point
   !> dim= Nky
   real(dp), allocatable :: largestgap(:)

   !> Z2 calculation for time reversal invariant system
   real(dp) :: Z2

   allocate(wcc(NumberofSelectedOccupiedBands, Nk2))
   allocate(largestgap(Nk2))
   largestgap= 0d0
   wcc= 0d0

   call  wannier_center3D_plane_func(K3D_start, K3D_vec1, K3D_vec2, &
      largestgap, wcc, Z2)
!~    call  w3d_shell(K3D_start, K3D_vec1, K3D_vec2, Num_wann, wcc)


   outfileindex= outfileindex+ 1
   if (cpuid==0) then
      open(unit=outfileindex, file='wcc.dat')
      write(outfileindex, '(10000A16)')'k', 'largestgap', 'sum(wcc(:,ik))', &
         'wcc(:, ik)'
      do ik2=1, Nk2
         write(outfileindex, '(10000f16.8)') dble(ik2-1)/dble(Nk2-1), &
            largestgap(ik2), dmod(sum(wcc(:, ik2)), 1d0), &
            wcc(:, ik2)
      enddo
      close(outfileindex)
   endif

   if (cpuid==0) write(stdout, *)'Z2 for the plane you choose: ', Z2

   !> generate gnu script for wannier charge center plots
   outfileindex= outfileindex+ 1
   if (cpuid==0) then
      open(unit=outfileindex, file='wcc.gnu')
      write(outfileindex, '(a)')"# gnuplot version > 5.4"
      write(outfileindex, '(a)')"set encoding iso_8859_1"
      write(outfileindex, '(a)')'set terminal pdf enhanced color font ",24" size 5, 4'
      write(outfileindex, '(a)')"set output 'wcc.pdf'"
      write(outfileindex, '(a)')'unset key '
      write(outfileindex, '(a)')'set border lw 3 '
      write(outfileindex, '(a)')'set xtics offset 0, 0.2'
      write(outfileindex, '(a)')'set xtics format "%4.1f" nomirror out '
      write(outfileindex, '(a)')'set xlabel "k in unit of line 3 of KPLANE_BULK" font ",16"'
      write(outfileindex, '(a)')'set xlabel offset 0, 0.7 '
      write(outfileindex, '(a)')'set ytics 0.5 '
      write(outfileindex, '(a)')'set ytics format "%4.1f" nomirror out'
      write(outfileindex, '(a)')'set ylabel "WCC"'
      write(outfileindex, '(a)')'set ylabel offset 2, 0.0 '
      write(outfileindex, '(a)')'set xrange [0: 1.0]'
      write(outfileindex, '(a)')'set yrange [0:1]'
      write(outfileindex, '(a, i5, a)')"plot for [i=4: ", NumberofSelectedOccupiedBands+3, &
         "] 'wcc.dat' u 1:i w p  pt 7  ps 0.5 lc 'black'"
      close(outfileindex)
   endif

   return
end subroutine  wannier_center3D_plane0

   subroutine  wannier_center3D_plane_func(kstart, kvec1, kvec2, largest_gap, wcc, Z2)
      !> this suboutine is used for wannier center calculation for 3D system
      !> only for one plane

      use para
      use wmpi
      implicit none

      integer :: i, j, l , m , ia, imax
      integer :: ik1, ik2, ierr

      real(dp), intent(in) :: kstart(3)
      real(dp), intent(in) :: kvec1(3)  ! the integration direction
      real(dp), intent(in) :: kvec2(3)
      real(dp), intent(out) :: largest_gap(Nk2)
      real(dp), intent(out) :: wcc(Numoccupied, Nk2)

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

      real(dp) :: Umatrix_t(3, 3)

      !> b.r
      real(dp) :: br

      !> exp(-i*b.R)
      complex(dp) :: ratio

      real(dp) :: k(3)
      real(dp) :: b(3)

      real(dp) :: maxgap
      real(dp) :: maxgap0

      !> Z2 calculation for time reversal invariant system
      real(dp) :: Z2
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
      real(dp) :: k0(3), k1(3), k2(3)

      allocate(kpoints(3, Nk1, Nk2))
      kpoints= 0d0

      allocate(Lambda_eig(Numoccupied))
      allocate(Lambda(Numoccupied, Numoccupied))
      allocate(Lambda0(Numoccupied, Numoccupied))
      allocate(Mmnkb(Numoccupied, Numoccupied))
      allocate(Mmnkb_com(Numoccupied, Numoccupied))
      allocate(Mmnkb_full(Num_wann, Num_wann))
      allocate(hamk(Num_wann, Num_wann))
      allocate(hamk_dag(Num_wann, Num_wann))
      allocate(Eigenvector(Num_wann, Num_wann, Nk1))
      allocate(eigenvalue(Num_wann))
      allocate(U(Numoccupied, Numoccupied))
      allocate(Sigma(Numoccupied, Numoccupied))
      allocate(VT(Numoccupied, Numoccupied))
      allocate(WannierCenterKy(Numoccupied, Nk2))
      allocate(WannierCenterKy_mpi(Numoccupied, Nk2))
      allocate(xnm(Numoccupied))
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
      !> the first dimension should be in one primitive cell, [0, 1] 
      !> the first dimension is the integration direction
      !> the WCCs are calculated along the second k line
      !> For chern number calculation, the second k line should be in one primitive
      !> vector. 
      !> For  Z2 number clculation, the second k line should 
      !> from one TRIM to another TRIM, usually, we study half of the
      !> reciprocal lattice vector
      k0= kstart ! 
      k1= kvec1   !  
      k2= kvec2   ! 

      do ik2=1, Nk2
         do ik1=1, Nk1
            kpoints(:, ik1, ik2)= k0+ k1*(ik1-1d0)/dble(Nk1)+ k2*(ik2-1d0)/dble(Nk2-1)
         enddo
      enddo
      b= k1/dble(Nk1)
      b= b(1)*Origin_cell%Kua+b(2)*Origin_cell%Kub+b(3)*Origin_cell%Kuc

      Umatrix_t= transpose(Umatrix)
      call inv_r(3, Umatrix_t)

      !>> Get wannier center for ky=0 plane
      !> for each ky, we can get wanniercenter
      do ik2=1+ cpuid, Nk2, num_cpu
         if (cpuid.eq.0) write(stdout, *)' Wilson loop ',  'ik, Nk', ik2, nk2
         Lambda0=0d0
         do i=1, Numoccupied
            Lambda0(i, i)= 1d0
         enddo

         !> for each k1, we get the eigenvectors
         do ik1=1, Nk1
            k= kpoints(:, ik1, ik2)


            ! generate bulk Hamiltonian
            if (index(KPorTB, 'KP')/=0)then
               call ham_bulk_kp(k, Hamk)
            else
              !> deal with phonon system
              if (index(Particle,'phonon')/=0.and.LOTO_correction) then
                 call ham_bulk_LOTO(k, Hamk)
              else
                 call ham_bulk_latticegauge    (k, Hamk)
              endif
            endif


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
               br= b(1)*Origin_cell%wannier_centers_cart(1, m)+ &
                   b(2)*Origin_cell%wannier_centers_cart(2, m)+ &
                   b(3)*Origin_cell%wannier_centers_cart(3, m)
               ratio= cos(br)- zi* sin(br)
         
               do j=1, Numoccupied
                  do i=1, Numoccupied
                     Mmnkb(i, j)=  Mmnkb(i, j)+ &
                        conjg(hamk_dag(m, i))* hamk(m, j)* ratio
                  enddo ! i
               enddo ! j
            enddo ! m

            !> perform Singluar Value Decomposed of Mmnkb
            call zgesvd_pack(Numoccupied, Mmnkb, U, Sigma, VT)

            !> after the calling of zgesvd_pack, Mmnkb becomes a temporal matrix
            U = conjg(transpose(U))
            VT= conjg(transpose(VT))
            call mat_mul(Numoccupied, VT, U, Mmnkb)

            call mat_mul(Numoccupied, Mmnkb, Lambda0, Lambda)
            Lambda0 = Lambda
         enddo  !< ik1

         !> diagonalize Lambda to get the eigenvalue 
         call zgeev_pack(Numoccupied, Lambda, Lambda_eig)
         do i=1, Numoccupied
            WannierCenterKy(i, ik2)= aimag(log(Lambda_eig(i)))/2d0/pi
            WannierCenterKy(i, ik2)= mod(WannierCenterKy(i, ik2)+10d0, 1d0) 
         enddo

         call sortheap(Numoccupied, WannierCenterKy(:, ik2))

         maxgap0= -99999d0
         imax= Numoccupied
         do i=1, Numoccupied
            if (i/=Numoccupied) then
               maxgap= WannierCenterKy(i+1, ik2)- WannierCenterKy(i, ik2)
            else
               maxgap=1d0+ WannierCenterKy(1, ik2)- WannierCenterKy(Numoccupied, ik2)
            endif

            if (maxgap>maxgap0) then
               maxgap0= maxgap
               imax= i
            endif

         enddo

         if (imax==Numoccupied) then
            largestgap(ik2)= (WannierCenterKy(1, ik2)+ &
               WannierCenterKy(Numoccupied, ik2) +1d0)/2d0
            largestgap(ik2)= mod(largestgap(ik2), 1d0)
         else
            largestgap(ik2)= (WannierCenterKy(imax+1, ik2)+ &
               WannierCenterKy(imax, ik2))/2d0
         endif


      enddo !< ik2

      WannierCenterKy_mpi= 0d0
      largestgap_mpi= 0d0
#if defined (MPI)
      call mpi_allreduce(WannierCenterKy, WannierCenterKy_mpi, &
           size(WannierCenterKy), mpi_dp, mpi_sum, mpi_cmw, ierr)
      call mpi_allreduce(largestgap, largestgap_mpi, &
           size(largestgap), mpi_dp, mpi_sum, mpi_cmw, ierr)
#else
      WannierCenterKy_mpi= WannierCenterKy
      largestgap_mpi= largestgap
#endif

      !> Z2 calculation Alexey Soluyanov arXiv:1102.5600

      Delta= 0
      !> for each iky, we get a Deltam
      do ik2=1, nk2-1
      
         !> largestgap position
         zm= largestgap_mpi(ik2)
        !if (ik2==nk2) then
        !   zm1= largestgap_mpi(1)
        !   xnm= WannierCenterKy_mpi(1:Numoccupied, 1)
        !else
            zm1= largestgap_mpi(ik2+1)
            xnm= WannierCenterKy_mpi(1:Numoccupied, ik2+1)
        !endif
         Deltam= 1
         do i=1, Numoccupied
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

      largest_gap= largestgap
      wcc= WannierCenterKy_mpi

      return
   end subroutine  wannier_center3D_plane_func


subroutine  wannier_center3D_k
   !> this suboutine is used for wannier center calculation for 3D system
   !> only for one plane

   use para
   use wmpi
   implicit none
   integer :: ik2,i,j

   real(dp), allocatable :: wcc(:, :)

   !> larges gap of each two wannier centers for a given k point
   !> dim= Nky
   real(dp), allocatable :: largestgap(:)

   !> Z2 calculation for time reversal invariant system
   real(dp) :: Z2

   real(dp),allocatable :: kpoints(:,:,:)
   integer :: nkp1,nkp2
   real(dp) :: kv3(3)

   nkp1=Nk1
   nkp2=knv2

   allocate(wcc(NumberofSelectedOccupiedBands, Nkp2))
   allocate(largestgap(Nkp2))

   largestgap= 0d0
   wcc= 0d0

   allocate(kpoints(3,nkp1,nkp2))
   kpoints=0d0

   kv3(1)=K3D_vec1(2)*K3D_vec2(3)-K3D_vec1(3)*K3D_vec2(2)
   kv3(2)=K3D_vec1(3)*K3D_vec2(1)-K3D_vec1(1)*K3D_vec2(3)
   kv3(3)=K3D_vec1(1)*K3D_vec2(2)-K3D_vec1(2)*K3D_vec2(1)
   kv3=kv3/sqrt(sum(kv3**2))

   do i=1,nkp1
      do j=1,nkp2
         kpoints(:,i,j)=k2_path(j,1)*K3D_vec1+k2_path(j,2)*K3D_vec2+kv3*dble(i-1)/dble(nkp1)
         kpoints(:,i,j)=matmul(kpoints(:,i,j),Urot)
      enddo
   enddo

!~    call  wannier_center3D_plane_func(K3D_start, K3D_vec1, K3D_vec2, &
!~       largestgap, wcc, Z2)
!~    call  w3d_shell(K3D_start, K3D_vec1, K3D_vec2, Num_wann, wcc)
   call wannier_center3D_kpath(kpoints,nkp1,nkp2,largestgap,wcc,Z2)


   outfileindex= outfileindex+ 1
   if (cpuid==0) then
      open(unit=outfileindex, file='wcc.dat')
      write(outfileindex, '(10000A16)')'k', 'largestgap', 'sum(wcc(:,ik))', &
         'wcc(:, ik)'
      do ik2=1, Nkp2
         write(outfileindex, '(10000f16.8)') k2len(ik2), &
            largestgap(ik2), dmod(sum(wcc(:, ik2)), 1d0), &
            wcc(:, ik2)
      enddo
      close(outfileindex)
   endif

   if (cpuid==0) write(stdout, *)'Z2 for the plane you choose: ', Z2

   !> generate gnu script for wannier charge center plots
   outfileindex= outfileindex+ 1
   if (cpuid==0) then
      open(unit=outfileindex, file='wcc.gnu')
      write(outfileindex, '(a)')"set encoding iso_8859_1"
      write(outfileindex, '(a)')'set terminal  postscript enhanced color font ",30"'
      write(outfileindex, '(a)')"set output 'wcc.eps'"
      write(outfileindex, '(a)')'unset key '
      write(outfileindex, '(a)')'set border lw 3 '
      write(outfileindex, '(a)')'set xtics offset 0, 0.2'
      write(outfileindex, '(a)')'set xtics format "%4.1f" nomirror out '
      write(outfileindex, '(a)')'set xlabel "k" '
      write(outfileindex, '(a)')'set xlabel offset 0, 0.7 '
      write(outfileindex, '(a)')'set ytics 0.5 '
      write(outfileindex, '(a)')'set ytics format "%4.1f" nomirror out'
      write(outfileindex, '(a)')'set ylabel "WCC"'
      write(outfileindex, '(a)')'set ylabel offset 2, 0.0 '
      write(outfileindex, '(a, f10.5, a)')'set xrange [0: ', maxval(k2len), ']'
      write(outfileindex, '(a)')'set yrange [0:1]'

      write(outfileindex, 202, advance="no") (trim(k2line_name(i)), k2line_stop(i), i=1, nk2lines)
      write(outfileindex, 203)trim(k2line_name(nk2lines+1)), k2line_stop(nk2lines+1)

      do i=1, nk2lines-1
         write(outfileindex, 204)k2line_stop(i+1), 0.0, k2line_stop(i+1), 1.0

      enddo

      write(outfileindex, '(a)')"plot 'wcc.dat' u 1:2 w l lw 2  lc 'blue', \"
      write(outfileindex, '(a, i5, a)')" for [i=4: ", NumberofSelectedOccupiedBands+3, "] 'wcc.dat' u 1:i w p  pt 7  ps 1.1 lc 'red'"
      close(outfileindex)
   endif

202 format('set xtics (',:20('"',A3,'" ',F10.5,','))
203 format(A3,'" ',F10.5,')')
204 format('set arrow from ',F10.5,',',F10.5, &
      ' to ',F10.5,',',F10.5, ' nohead')

   return
end subroutine  wannier_center3D_k



subroutine  wannier_center3D_kpath(kpoints, nkp1, nkp2, largest_gap, wcc, Z2)
   !> this suboutine is used for wannier center calculation for 3D system
   !> only for one plane

   use para
   use wmpi
   implicit none

   integer :: i, j, m, imax
   integer :: ik1, ik2, ierr

   real(dp), intent(in) :: kpoints(3,nkp1,nkp2)
   integer, intent(in) :: nkp1  ! the integration direction
   integer, intent(in) :: nkp2
   real(dp), intent(out) :: largest_gap(Nkp2)
   real(dp), intent(out) :: wcc(NumberofSelectedOccupiedBands, Nkp2)

   !> k points in kx-ky plane
!~    real(dp), allocatable ::

   !> hamiltonian for each k point
   !> and also the eigenvector of hamiltonian after eigensystem_c
   complex(dp), allocatable :: Hamk(:, :)
!~    complex(dp), allocatable :: Hamk_dag(:, :)

   !> eigenvector for each kx
   complex(dp), allocatable :: Eigenvector(:, :, :)


   complex(dp), allocatable :: Lambda(:, :)
!~    complex(dp), allocatable :: Lambda0(:, :)

   !> three matrix for SVD
   !> M= U.Sigma.V^\dag
   !> VT= V^\dag
!~    complex(dp), allocatable :: U(:, :)
!~    real   (dp), allocatable :: Sigma(:, :)
!~    complex(dp), allocatable :: VT(:, :)

   !> wannier centers for each ky, bands
   real(dp), allocatable :: WannierCenterKy(:, :)
   real(dp), allocatable :: WannierCenterKy_mpi(:, :)

   !> larges gap of each two wannier centers for a given k point
   !> dim= Nky
   real(dp), allocatable :: largestgap(:)
   real(dp), allocatable :: largestgap_mpi(:)

   !> eigenvalue
   real(dp), allocatable :: eigenvalue(:)

   real(dp) :: Umatrix_t(3, 3)

   !> b.r
   real(dp) :: br

   !> exp(-i*b.R)
   complex(dp) :: ratio

   real(dp) :: k(3)
   real(dp) :: b(3)

   real(dp) :: maxgap
   real(dp) :: maxgap0

   !> Z2 calculation for time reversal invariant system
   real(dp) :: Z2
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
   real(dp) :: k0(3), k1(3), k2(3)

   complex(dp),allocatable :: gauge_shift(:,:)



   allocate(Lambda(NumberofSelectedOccupiedBands, NumberofSelectedOccupiedBands))
   allocate(hamk(Num_wann, Num_wann))
   allocate(Eigenvector(Num_wann, Num_wann, Nkp1))
   allocate(eigenvalue(Num_wann))

   allocate(WannierCenterKy(NumberofSelectedOccupiedBands, Nkp2))
   allocate(WannierCenterKy_mpi(NumberofSelectedOccupiedBands, Nkp2))
   allocate(xnm(NumberofSelectedOccupiedBands))
   allocate(largestgap(Nkp2))
   allocate(largestgap_mpi(Nkp2))

   allocate(gauge_shift(Num_wann,Num_wann))
   !~ the gauge shift matrix to guarantee the gauge consistency

   largestgap= 0d0
   largestgap_mpi= 0d0
   WannierCenterKy= 0d0
   WannierCenterKy_mpi= 0d0
   hamk=0d0
   Eigenvector=0d0
   Lambda =0d0
   gauge_shift=0d0

   !> set k plane
   !> the first dimension should be in one primitive cell, [0, 1]
   !> the first dimension is the integration direction
   !> the WCCs are calculated along the second k line
   !> For chern number calculation, the second k line should be in one primitive
   !> vector.
   !> For  Z2 number clculation, the second k line should
   !> from one TRIM to another TRIM, usually, we study half of the
   !> reciprocal lattice vector


   do i=1, Num_wann
      br= sum((kpoints(:, Nkp1, 1)-kpoints(:, 1, 1))* Origin_cell%wannier_centers_direct(:, i))
      gauge_shift(i,i)=cos(twopi*br)- zi*sin(twopi*br)
   enddo

   Umatrix_t= transpose(Umatrix)
   call inv_r(3, Umatrix_t)

   !>> Get wannier center for ky=0 plane
   !> for each ky, we can get wanniercenter
   do ik2=1+ cpuid, Nkp2, num_cpu
      if (cpuid.eq.0) write(stdout, *)' Wilson loop ',  'ik, Nk', ik2, nkp2
      do i=1, NumberofSelectedOccupiedBands
      enddo

      !> for each k1, we get the eigenvectors
      do ik1=1, Nkp1
         k= kpoints(:, ik1, ik2)

         ! generate bulk Hamiltonian
         if (index(KPorTB, 'KP')/=0)then
            call ham_bulk_kp(k, Hamk)
         else
            !> deal with phonon system
            if (index(Particle,'phonon')/=0.and.LOTO_correction) then
               call ham_bulk_LOTO(k, Hamk)
            else
               call ham_bulk_atomicgauge(k, Hamk)
!~                call ham_bulk_latticegauge    (k, Hamk)
            endif
         endif

         call eigensystem_c('V', 'U', Num_wann, hamk, eigenvalue)

         Eigenvector(:, :, ik1)= hamk
      enddo

!~       Eigenvector(:, :, Nk1)= Eigenvector(:, :, 1)
      Eigenvector(:, :, Nkp1)= matmul(gauge_shift,Eigenvector(:, :, 1))

      !~ shift the gauge such that k1 and k1+2pi are in same wavefunction gauge
      call wilson_loop(Eigenvector,Nkp1,Num_wann,Lambda,NumberofSelectedOccupiedBands,WannierCenterKy(:, ik2))

      call sortheap(NumberofSelectedOccupiedBands, WannierCenterKy(:, ik2))

      maxgap0= -99999d0
      imax= NumberofSelectedOccupiedBands
      do i=1, NumberofSelectedOccupiedBands
         if (i/=NumberofSelectedOccupiedBands) then
            maxgap= WannierCenterKy(i+1, ik2)- WannierCenterKy(i, ik2)
         else
            maxgap=1d0+ WannierCenterKy(1, ik2)- WannierCenterKy(NumberofSelectedOccupiedBands, ik2)
         endif

         if (maxgap>maxgap0) then
            maxgap0= maxgap
            imax= i
         endif

      enddo

      if (imax==NumberofSelectedOccupiedBands) then
         largestgap(ik2)= (WannierCenterKy(1, ik2)+ &
            WannierCenterKy(NumberofSelectedOccupiedBands, ik2) +1d0)/2d0
         largestgap(ik2)= mod(largestgap(ik2), 1d0)
      else
         largestgap(ik2)= (WannierCenterKy(imax+1, ik2)+ &
            WannierCenterKy(imax, ik2))/2d0
      endif


   enddo !< ik2

   WannierCenterKy_mpi= 0d0
   largestgap_mpi= 0d0
#if defined (MPI)
   call mpi_allreduce(WannierCenterKy, WannierCenterKy_mpi, &
      size(WannierCenterKy), mpi_dp, mpi_sum, mpi_cmw, ierr)
   call mpi_allreduce(largestgap, largestgap_mpi, &
      size(largestgap), mpi_dp, mpi_sum, mpi_cmw, ierr)
#else
   WannierCenterKy_mpi= WannierCenterKy
   largestgap_mpi= largestgap
#endif

   !> Z2 calculation Alexey Soluyanov arXiv:1102.5600

   Delta= 0
   !> for each iky, we get a Deltam
   do ik2=1, nkp2-1

      !> largestgap position
      zm= largestgap_mpi(ik2)
      !if (ik2==nk2) then
      !   zm1= largestgap_mpi(1)
      !   xnm= WannierCenterKy_mpi(1:NumberofSelectedOccupiedBands, 1)
      !else
      zm1= largestgap_mpi(ik2+1)
      xnm= WannierCenterKy_mpi(1:NumberofSelectedOccupiedBands, ik2+1)
      !endif
      Deltam= 1
      do i=1, NumberofSelectedOccupiedBands
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

   largest_gap= largestgap
   wcc= WannierCenterKy_mpi

   return
end subroutine  wannier_center3D_kpath





subroutine  wannier_center3D_nodalline
   !> this suboutine is used for wannier center calculation for 3D system
   !> for all NL points specified in input.dat

   use para
   use wmpi
   implicit none

   integer :: i, ik2, j
   integer :: chirality

   character(40) :: epsfilename, wccfilename

   real(dp) :: k0(3)
   integer, allocatable :: chirality_all(:)
   real(dp), allocatable :: wcc(:, :)
   real(dp), allocatable :: wcc_all(:, :, :)
   real(dp), allocatable :: wcc_sum_all(:, :)


   allocate(chirality_all(Num_NLs))
   allocate(wcc(NumberofSelectedOccupiedBands, nk2))
   allocate(wcc_all(NumberofSelectedOccupiedBands, nk2, Num_NLs))
   allocate(wcc_sum_all(nk2, Num_NLs))

   wcc= 0d0
   wcc_all= 0d0
   wcc_sum_all= 0d0

   do i= 1, Num_NLs
      k0= NL_center_position_cart(:, i)
      call wannier_center3D_weyl_func_halftorus(k0, Rbig_NL, rsmall_a_NL, wcc, chirality)
      chirality_all(i)= chirality
      wcc_all(:, :, i) = wcc
      wcc_sum_all(:, :)= dmod(sum(wcc_all, dim=1), 1d0)
   enddo

   if (cpuid==0) write(stdout, '(a)')'Chiralities'
   if (cpuid==0) write(stdout, '(a, a8, 5a10,a15)')'#', 'k1', 'k2', 'k3', 'kx', 'ky', 'kz', 'Chirality'
   if (cpuid==0) then
      do i=1, Num_NLs
         write(stdout, '(6f10.5,i8)')NL_center_position_direct(:, i), NL_center_position_cart(:,i)*Angstrom2atomic,chirality_all(i)
      enddo
   endif

   outfileindex= outfileindex+ 1
   if (cpuid==0) then
      open(unit=outfileindex, file='wanniercenter3D_NL.dat')
      write(outfileindex, '(a16, 10000i16)')'# Chirality', chirality_all
      write(outfileindex, '(10000a16)')'# k ', ('phase', j=1, Num_NLs)
      do ik2=1, Nk2
         write(outfileindex, '(10000f16.8)') (ik2-1d0)/(Nk2-1), &
            (wcc_sum_all(ik2, j), j=1, Num_NLs), ((wcc_all(i, ik2, j), i=1, NumberofSelectedOccupiedBands), j=1, Num_NLs)
      enddo
      close(outfileindex)
   endif

   !> generate gnu script for wannier charge center plots
   do i=1, Num_NLs

      if (i<10) then
         write(wccfilename, '(a,i1,a)')'wanniercenter3D_NL_', i, '.gnu'
      elseif (i>=10.and.i<100) then
         write(wccfilename, '(a,i2,a)')'wanniercenter3D_NL_', i, '.gnu'
      elseif (i>=100.and.i<1000) then
         write(wccfilename, '(a,i3,a)')'wanniercenter3D_NL_', i, '.gnu'
      elseif (i>=1000.and.i<10000) then
         write(wccfilename, '(a,i4,a)')'wanniercenter3D_NL_', i, '.gnu'
      endif

      if (i<10) then
         write(epsfilename, '(a,i1,a)')'wanniercenter3D_NL_', i, '.eps'
      elseif (i>=10.and.i<100) then
         write(epsfilename, '(a,i2,a)')'wanniercenter3D_NL_', i, '.eps'
      elseif (i>=100.and.i<1000) then
         write(epsfilename, '(a,i3,a)')'wanniercenter3D_NL_', i, '.eps'
      elseif (i>=1000.and.i<10000) then
         write(epsfilename, '(a,i4,a)')'wanniercenter3D_NL_', i, '.eps'
      endif

      outfileindex= outfileindex+ 1
      if (cpuid==0) then
         open(unit=outfileindex, file=wccfilename)
         write(outfileindex, '(a)')"set encoding iso_8859_1"
         write(outfileindex, '(a)')'set terminal  postscript enhanced color font ",24"'
         write(outfileindex, '(3a)')"set output '", trim(epsfilename), "'"
         write(outfileindex, '(a)')'unset key '
         write(outfileindex, '(a)')'set border lw 3 '
         write(outfileindex, '(a)')'set xtics offset 0, 0.2'
         write(outfileindex, '(a)')'set xtics ("South" 0, "North" 1)'
         write(outfileindex, '(a)')'set xtics format "%4.1f" nomirror in'
         write(outfileindex, '(a)')'set xlabel "k" '
         write(outfileindex, '(a)')'set xlabel offset 0, 0.7 '
         write(outfileindex, '(a)')'set ytics 1.0 '
         write(outfileindex, '(a)')'set ytics format "%4.1f" nomirror in'
         write(outfileindex, '(a)')'set ylabel "WCC"'
         write(outfileindex, '(a, f8.5, a, f8.5, a, f8.5, a)') &
            'set title "NL (', NL_center_position_direct(1, i), ',', &
            NL_center_position_direct(2, i), ',', &
            NL_center_position_direct(3, i), ') "'
         write(outfileindex, '(a)')'set ylabel offset 2, 0.0 '
         write(outfileindex, '(a)')'set xrange [0: 1.0]'
         write(outfileindex, '(a)')'set yrange [0:1]'
         write(outfileindex, '(a,i5,a)')"plot 'wanniercenter3D_NL.dat' u 1:",i+1, &
            " w p ps 1.5 pt 7 lc rgb '#696969'"
         close(outfileindex)
      endif
   enddo


   return
end subroutine wannier_center3D_nodalline


subroutine  wannier_center3D_weyl
   !> this suboutine is used for wannier center calculation for 3D system
   !> for all weyl points specified in input.dat

   use para
   use wmpi
   implicit none

   integer :: i, ik2, j
   integer :: chirality

   character(40) :: epsfilename, wccfilename

   real(dp) :: k0(3)
   integer, allocatable :: chirality_all(:)
   real(dp), allocatable :: wcc(:, :)
   real(dp), allocatable :: wcc_all(:, :, :)
   real(dp), allocatable :: wcc_sum_all(:, :)


   allocate(chirality_all(Num_Weyls))
   allocate(wcc(NumberofSelectedOccupiedBands, nk2))
   allocate(wcc_all(NumberofSelectedOccupiedBands, nk2, Num_Weyls))
   allocate(wcc_sum_all(nk2, Num_Weyls))

   wcc= 0d0
   wcc_all= 0d0
   wcc_sum_all= 0d0

   do i= 1, Num_Weyls
      if (cpuid==0) write(stdout, '(a)')' '
      if (cpuid==0) write(stdout, '(a, i3, a)')'>> We are calculation the chirality of the ', i, "'th Weyl point"
      k0= weyl_position_cart(:, i)
      if (cpuid==0) write(stdout, '(a, 3f12.6)')'>> k (Cartesian coordinates): ', k0
      call wannier_center3D_weyl_func_sphere(k0, kr0, wcc, chirality)
      chirality_all(i)= chirality
      wcc_all(:, :, i) = wcc
      wcc_sum_all(:, :)= dmod(sum(wcc_all, dim=1), 1d0)
   enddo

   if (cpuid==0) write(stdout, '(a)')'Chiralities'
   if (cpuid==0) write(stdout, '(a, a8, 5a10,a15)')'#', 'k1', 'k2', 'k3', 'kx', 'ky','kz', 'Chirality'
   if (cpuid==0) then
      do i=1, Num_Weyls
         write(stdout, '(6f10.5,i8)')weyl_position_direct(:, i), weyl_position_cart(:,i)*Angstrom2atomic, chirality_all(i)
      enddo
   endif

   outfileindex= outfileindex+ 1
   if (cpuid==0) then
      open(unit=outfileindex, file='wanniercenter3D_Weyl.dat')
      write(outfileindex, '(a16, 10000i16)')'# Chirality', chirality_all
      write(outfileindex, '(10000a16)')'# k ', ('phase', j=1, Num_Weyls)
      do ik2=1, Nk2
         write(outfileindex, '(10000f16.8)') (ik2-1d0)/(Nk2-1), &
            (wcc_sum_all(ik2, j), j=1, Num_Weyls), ((wcc_all(i, ik2, j), i=1, NumberofSelectedOccupiedBands), j=1, Num_Weyls)
      enddo
      close(outfileindex)
   endif

   !> generate gnu script for wannier charge center plots
   do i=1, Num_Weyls

      if (i<10) then
         write(wccfilename, '(a,i1,a)')'wanniercenter3D_Weyl_', i, '.gnu'
      elseif (i>=10.and.i<100) then
         write(wccfilename, '(a,i2,a)')'wanniercenter3D_Weyl_', i, '.gnu'
      elseif (i>=100.and.i<1000) then
         write(wccfilename, '(a,i3,a)')'wanniercenter3D_Weyl_', i, '.gnu'
      elseif (i>=1000.and.i<10000) then
         write(wccfilename, '(a,i4,a)')'wanniercenter3D_Weyl_', i, '.gnu'
      endif

      if (i<10) then
         write(epsfilename, '(a,i1,a)')'wanniercenter3D_Weyl_', i, '.eps'
      elseif (i>=10.and.i<100) then
         write(epsfilename, '(a,i2,a)')'wanniercenter3D_Weyl_', i, '.eps'
      elseif (i>=100.and.i<1000) then
         write(epsfilename, '(a,i3,a)')'wanniercenter3D_Weyl_', i, '.eps'
      elseif (i>=1000.and.i<10000) then
         write(epsfilename, '(a,i4,a)')'wanniercenter3D_Weyl_', i, '.eps'
      endif

      outfileindex= outfileindex+ 1
      if (cpuid==0) then
         open(unit=outfileindex, file=wccfilename)
         write(outfileindex, '(a)')"set encoding iso_8859_1"
         write(outfileindex, '(a)')'set terminal  postscript enhanced color font ",24"'
         write(outfileindex, '(3a)')"set output '", trim(epsfilename), "'"
         write(outfileindex, '(a)')'unset key '
         write(outfileindex, '(a)')'set border lw 3 '
         write(outfileindex, '(a)')'set xtics offset 0, 0.2'
         write(outfileindex, '(a)')'set xtics ("South" 0, "North" 1)'
         write(outfileindex, '(a)')'set xtics format "%4.1f" nomirror in'
         write(outfileindex, '(a)')'set xlabel "k" '
         write(outfileindex, '(a)')'set xlabel offset 0, 0.7 '
         write(outfileindex, '(a)')'set ytics 1.0 '
         write(outfileindex, '(a)')'set ytics format "%4.1f" nomirror in'
         write(outfileindex, '(a)')'set ylabel "WCC"'
         write(outfileindex, '(a, f8.5, a, f8.5, a, f8.5, a)') &
            'set title "Weyl (', weyl_position_direct(1, i), ',', &
            weyl_position_direct(2, i), ',', &
            weyl_position_direct(3, i), ') "'
         write(outfileindex, '(a)')'set ylabel offset 2, 0.0 '
         write(outfileindex, '(a)')'set xrange [0: 1.0]'
         write(outfileindex, '(a)')'set yrange [0:1]'
         write(outfileindex, '(a,i5,a)')"plot 'wanniercenter3D_Weyl.dat' u 1:",i+1, &
            " w p ps 1.5 pt 7 lc rgb '#696969'"
         close(outfileindex)
      endif
   enddo


   return
end subroutine wannier_center3D_weyl


subroutine  wannier_center3D_weyl_func_kpoints(kpoints, wcc, chirality)
   !> this suboutine is used for wannier center calculation for 3D system
   !> only for one weyl point, with given kpoints

   use para
   use wmpi
   implicit none

   integer :: i, j, m, ik1, ik2, ierr

   !> inout variables
   !> the first dimension should be in one primitive cell, [0, 2*pi]
   !> the first dimension is the integration direction
   !> the WCCs are calculated along the second k line
   real(dp), intent(in) :: kpoints(3, Nk1, Nk2)
   real(dp), intent(out) :: wcc(NumberofSelectedOccupiedBands, Nk2)
   integer , intent(out) :: chirality

   !> hamiltonian for each k point
   !> and also the eigenvector of hamiltonian after eigensystem_c
   complex(dp), allocatable :: Hamk(:, :), Hamk_dag(:, :)

   !> eigenvalue
   real(dp), allocatable :: eigenvalue(:)

   !> eigenvector for each kx
   complex(dp), allocatable :: Eigenvector(:, :, :)

   !> Mmnkb=<u_n(k)|u_m(k+b)>
   !> |u_n(k)> is the periodic part of wave function
   complex(dp), allocatable :: Mmnkb(:, :), Mmnkb_com(:, :), Mmnkb_full(:, :)
   complex(dp), allocatable :: Lambda_eig(:), Lambda(:, :), Lambda0(:, :)

   !> three matrix for SVD
   !> M= U.Sigma.V^\dag
   !> VT= V^\dag
   complex(dp), allocatable :: U(:, :)
   real   (dp), allocatable :: Sigma(:, :)
   complex(dp), allocatable :: VT(:, :)

   !> wannier centers for each ky, bands
   real(dp), allocatable :: WannierCenterKy(:, :), WannierCenterKy_mpi(:, :)

   !> sumation for chirality number
   real(dp), allocatable :: wcc_sum(:)

   real(dp) :: Umatrix_t(3, 3)

   !> b.r
   real(dp) :: br, k(3), b(3)

   !> exp(-i*b.R)
   complex(dp) :: ratio


   real(dp) :: gap_sum, gap_step
   real(dp), allocatable :: xnm(:)


   allocate(Lambda_eig(NumberofSelectedOccupiedBands),Lambda(NumberofSelectedOccupiedBands, NumberofSelectedOccupiedBands),Lambda0(NumberofSelectedOccupiedBands, NumberofSelectedOccupiedBands))
   allocate(Mmnkb(NumberofSelectedOccupiedBands, NumberofSelectedOccupiedBands),Mmnkb_com(NumberofSelectedOccupiedBands, NumberofSelectedOccupiedBands))
   allocate(Mmnkb_full(Num_wann, Num_wann),hamk(Num_wann, Num_wann))
   allocate(hamk_dag(Num_wann, Num_wann),Eigenvector(Num_wann, Num_wann, Nk1))
   allocate(eigenvalue(Num_wann),U(NumberofSelectedOccupiedBands, NumberofSelectedOccupiedBands))
   allocate(Sigma(NumberofSelectedOccupiedBands, NumberofSelectedOccupiedBands),VT(NumberofSelectedOccupiedBands, NumberofSelectedOccupiedBands))
   allocate(WannierCenterKy(NumberofSelectedOccupiedBands, Nk2),WannierCenterKy_mpi(NumberofSelectedOccupiedBands, Nk2))
   allocate(wcc_sum(Nk2),xnm(NumberofSelectedOccupiedBands))
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
   wcc_sum= 0d0

   Umatrix_t= transpose(Umatrix)
   call inv_r(3, Umatrix_t)

   !> for each ky, we can get wanniercenter
   do ik2=1+ cpuid, Nk2, num_cpu
      if (cpuid.eq.0) write(stdout, *)' Wilson loop ',  'ik, Nk', ik2, nk2
      Lambda0=0d0
      do i=1, NumberofSelectedOccupiedBands
         Lambda0(i, i)= 1d0
      enddo

      !> for each k1, we get the eigenvectors
      do ik1=1, Nk1
         k= kpoints(:, ik1, ik2)


         ! generate bulk Hamiltonian
         if (index(KPorTB, 'KP')/=0)then
            call ham_bulk_kp(k, Hamk)
         else
            !> deal with phonon system
            if (index(Particle,'phonon')/=0.and.LOTO_correction) then
               call ham_bulk_LOTO(k, Hamk)
            else
               call ham_bulk_latticegauge    (k, Hamk)
            endif
         endif


         !> diagonal hamk
         call eigensystem_c('V', 'U', Num_wann, hamk, eigenvalue)

         Eigenvector(:, :, ik1)= hamk
      enddo

      !> sum over k1 to get wanniercenters
      do ik1=1, Nk1
         if (ik1<Nk1) then
            b= kpoints(:, ik1+1, ik2)- kpoints(:, ik1, ik2)
         else
            b= kpoints(:, 1, ik2)- kpoints(:, 2, ik2)
         endif
         b= b(1)*Origin_cell%Kua+b(2)*Origin_cell%Kub+b(3)*Origin_cell%Kuc

         !> <u_k|u_k+1>
         Mmnkb= 0d0
         hamk_dag= Eigenvector(:, :, ik1)
         if (ik1==Nk1) then
            hamk= Eigenvector(:, :, 1)
         else
            hamk= Eigenvector(:, :, ik1+ 1)
         endif
         do m=1, Num_wann
            br= b(1)*Origin_cell%wannier_centers_cart(1, m)+ &
               b(2)*Origin_cell%wannier_centers_cart(2, m)+ &
               b(3)*Origin_cell%wannier_centers_cart(3, m)
            ratio= cos(br)- zi* sin(br)

            do j=1, NumberofSelectedOccupiedBands
               do i=1, NumberofSelectedOccupiedBands
                  Mmnkb(i, j)=  Mmnkb(i, j)+ &
                     conjg(hamk_dag(m, Selected_Occupiedband_index(i)))* hamk(m, Selected_Occupiedband_index(j))* ratio
               enddo ! i
            enddo ! j
         enddo ! m

         !> perform Singluar Value Decomposed of Mmnkb
         call zgesvd_pack(NumberofSelectedOccupiedBands, Mmnkb, U, Sigma, VT)

         !> after the calling of zgesvd_pack, Mmnkb becomes a temporal matrix
         U = conjg(transpose(U))
         VT= conjg(transpose(VT))
         call mat_mul(NumberofSelectedOccupiedBands, VT, U, Mmnkb)

         call mat_mul(NumberofSelectedOccupiedBands, Mmnkb, Lambda0, Lambda)
         Lambda0 = Lambda
      enddo  !< ik1

      !> diagonalize Lambda to get the eigenvalue
      call zgeev_pack(NumberofSelectedOccupiedBands, Lambda, Lambda_eig)
      do i=1, NumberofSelectedOccupiedBands
         WannierCenterKy(i, ik2)= aimag(log(Lambda_eig(i)))/2d0/pi
         WannierCenterKy(i, ik2)= mod(WannierCenterKy(i, ik2)+10d0, 1d0)
      enddo

      call sortheap(NumberofSelectedOccupiedBands, WannierCenterKy(:, ik2))

   enddo !< ik2

   WannierCenterKy_mpi= 0d0
#if defined (MPI)
   call mpi_allreduce(WannierCenterKy, WannierCenterKy_mpi, &
      size(WannierCenterKy), mpi_dp, mpi_sum, mpi_cmw, ierr)
#else
   WannierCenterKy_mpi= WannierCenterKy
#endif

   wcc= WannierCenterKy_mpi

   wcc_sum= dmod(sum(wcc, dim=1), 1d0)

   !> determine the chirality
   gap_sum= 0d0
   do ik2=1, Nk2-1
      gap_step= wcc_sum(ik2+1)- wcc_sum(ik2)
      if (abs(gap_step+1)<abs(gap_step)) then
         gap_step= gap_step+ 1
      elseif (abs(gap_step-1)<abs(gap_step))then
         gap_step= gap_step- 1
      endif
      gap_sum= gap_sum+ gap_step
   enddo

   chirality= nint(gap_sum)

   return
end subroutine  wannier_center3D_weyl_func_kpoints


subroutine  wannier_center3D_weyl_func_sphere(k0, r0, wcc, chirality)
   !> this suboutine is used for wannier center calculation for 3D system
   !> only for one weyl point
   !> k0 must be in Cartesian coordinates

   use para
   use wmpi
   implicit none

   integer :: ik1, ik2
   real(dp) :: theta, phi, r_para
   real(dp) :: k_cart(3), k_direct(3)

   !> inout variables
   real(dp), intent(in) :: k0(3)
   real(dp), intent(in) :: r0
   integer, intent(out) :: chirality
   real(dp), intent(out) :: wcc(NumberofSelectedOccupiedBands, Nk2)

   !> k points in kx-ky plane
   real(dp), allocatable :: kpoints(:, :, :)


   allocate(kpoints(3, Nk1, Nk2))
   kpoints= 0d0

   !> set k plane
   !> the first dimension should be in one primitive cell, [0, 2*pi]
   !> the first dimension is the integration direction
   !> the WCCs are calculated along the second k line
   do ik2=1, Nk2 ! the mesh for k line of that wcc are calculated
      theta= (1d0-(ik2- 1d0)/(Nk2- 1d0))* pi
      if (ik2== 1) theta= (1d0-(ik2- 1d0+ 0.10)/(Nk2- 1d0))* pi  ! avoid the south pole
      if (ik2== Nk2) theta= (1d0-(ik2- 1d0- 0.10)/(Nk2- 1d0))* pi  ! avoid the north pole
      do ik1=1, Nk1  ! the mesh along the integration direction
         phi= (ik1- 1d0)/Nk1* 2d0* pi
         r_para= r0* sin(theta)
         k_cart(1)= k0(1)+ r_para* cos(phi)
         k_cart(2)= k0(2)+ r_para* sin(phi)
         k_cart(3)= k0(3)+ r0* cos(theta)
         call cart_direct_rec(k_cart, k_direct)
         kpoints(:, ik1, ik2)= k_direct
      enddo
   enddo

   !> write out all kpoints in to weylball_kpoints.txt
   outfileindex= outfileindex+ 1
   if (cpuid.eq.0) then
      open (unit=outfileindex, file='weylball_kpoints.txt')

      do ik2=1, Nk2
         do ik1=1, NK1
            k_direct= kpoints(:, ik1, ik2)
            call direct_cart_rec(k_direct, k_cart)
            write(outfileindex, '(6f10.4)') k_cart*Angstrom2atomic, k_direct
         enddo
      enddo
      close(outfileindex)

   endif

   call wannier_center3D_weyl_func_kpoints(kpoints, wcc, chirality)

   return
end subroutine  wannier_center3D_weyl_func_sphere


subroutine  wannier_center3D_weyl_func_halftorus(k0, Rbig, rsmall, wcc, chirality)
   !> this suboutine is used for wannier center calculation for 3D system
   !> only for one weyl point
   !> The half torus is defined like this
   !> kx= (R+r cos\theta)cos\phi, ky=(R+r cos\theta) sin\phi, kz=r sin\theta
   !> theta is in [0, pi], phi is in [0, 2pi]
   !> R is Rbig, r is rsmall

   use para
   use wmpi
   implicit none

   integer :: ik1, ik2
   real(dp) :: theta, phi, r_para

   !> inout variables
   real(dp), intent(in) :: k0(3)
   real(dp), intent(in) :: Rbig
   real(dp), intent(in) :: rsmall
   integer, intent(out) :: chirality
   real(dp), intent(out) :: wcc(NumberofSelectedOccupiedBands, Nk2)

   real(dp) :: k_cart(3), k_direct(3)

   !> k points in kx-ky plane
   real(dp), allocatable :: kpoints(:, :, :)


   allocate(kpoints(3, Nk1, Nk2))
   kpoints= 0d0

   !> set k plane
   !> the first dimension should be in one primitive cell, phi=[0, 2*pi]
   !> the first dimension is the integration direction
   !> the WCCs are calculated along the second k line theta=[0, pi]
   do ik2=1, Nk2 ! the mesh for k line of that wcc are calculated ->theta
      theta= (1d0-(ik2- 1d0)/(Nk2- 1d0))* 1d0*pi
      if (ik2== 1) theta= (1d0-(ik2- 1d0+ 0.01)/(Nk2- 1d0))* pi  ! avoid the North pole
      if (ik2== Nk2) theta= (1d0-(ik2- 1d0- 0.01)/(Nk2- 1d0))* pi  ! avoid the south pole
      do ik1=1, Nk1  ! the mesh along the integration direction -> phi
         phi= (ik1- 1d0)/Nk1* 2d0* pi
         r_para= Rbig+ rsmall* cos(theta)
         k_cart(1)= k0(1)+ r_para* cos(phi)
         k_cart(2)= k0(2)+ r_para* sin(phi)
         k_cart(3)= k0(3)+ rsmall* sin(theta)
         call cart_direct_rec(k_cart, k_direct)
         kpoints(:, ik1, ik2)= k_direct
      enddo
   enddo

   call wannier_center3D_weyl_func_kpoints(kpoints, wcc, chirality)

   return
end subroutine  wannier_center3D_weyl_func_halftorus


subroutine  Chern_3D
   ! this suboutine is used for wannier center calculation for 3D system
   use para
   use wmpi
   implicit none

   integer :: ik2, i, j

   real(dp) :: kstart(3)
   real(dp) :: kvec1(3)
   real(dp) :: kvec2(3)

   real(dp) :: Chern
   real(dp) :: Chern_all(6)

   !> wannier centers for each ky, bands
   real(dp), allocatable :: wcc(:, :)
   real(dp), allocatable :: wcc_all(:, :, :)

   !> larges gap of each two wannier centers for a given k point
   !> dim= Nky
   real(dp), allocatable :: largestgap(:)
   real(dp), allocatable :: largestgap_all(:,:)

   allocate(wcc(NumberofSelectedOccupiedBands, Nk2))
   allocate(wcc_all(NumberofSelectedOccupiedBands, Nk2, 6))
   allocate(largestgap(Nk2))
   allocate(largestgap_all(Nk2,6))
   largestgap= 0d0
   wcc= 0d0
   wcc_all= 0d0

   !> integration over kc (\bar{c}}, wcc along kb, fixed k1=0
   kstart= (/0.0d0, 0.0d0, 0.0d0/)
   kvec1 = (/0.0d0, 0.0d0, 1.0d0/)
   kvec2 = (/0.0d0, 1.0d0, 0.0d0/)

   call  wannier_center3D_plane_func(kstart, kvec1, kvec2, &
      largestgap, wcc, Chern)
   call get_chern_number_from_wcc(NumberofSelectedOccupiedBands, Nk2, wcc, Chern)
   Chern_all(1)= Chern
   wcc_all(:, :, 1)= wcc
   largestgap_all(:, 1)= largestgap

   !> integration over kc (\bar{c}}, wcc along kb, fixed k1=ka/2
   kstart= (/0.5d0, 0.0d0, 0.0d0/)
   kvec1 = (/0.0d0, 0.0d0, 1.0d0/)
   kvec2 = (/0.0d0, 1.0d0, 0.0d0/)

   call  wannier_center3D_plane_func(kstart, kvec1, kvec2, &
      largestgap, wcc, Chern)
   call get_chern_number_from_wcc(NumberofSelectedOccupiedBands, Nk2, wcc, Chern)
   Chern_all(2)= Chern
   wcc_all(:, :, 2)= wcc
   largestgap_all(:, 2)= largestgap


   !> integration over kc (\bar{c}}, wcc along ka, fixed k2=0
   kstart= (/0.0d0, 0.0d0, 0.0d0/)
   kvec1 = (/0.0d0, 0.0d0, 1.0d0/)
   kvec2 = (/1.0d0, 0.0d0, 0.0d0/)

   call  wannier_center3D_plane_func(kstart, kvec1, kvec2, &
      largestgap, wcc, Chern)
   call get_chern_number_from_wcc(NumberofSelectedOccupiedBands, Nk2, wcc, Chern)
   Chern_all(3)= Chern
   wcc_all(:, :, 3)= wcc
   largestgap_all(:, 3)= largestgap

   !> integration over kc (\bar{c}}, wcc along ka, fixed k2=kb/2
   kstart= (/0.0d0, 0.5d0, 0.0d0/)
   kvec1 = (/0.0d0, 0.0d0, 1.0d0/)
   kvec2 = (/1.0d0, 0.0d0, 0.0d0/)

   call  wannier_center3D_plane_func(kstart, kvec1, kvec2, &
      largestgap, wcc, Chern)
   call get_chern_number_from_wcc(NumberofSelectedOccupiedBands, Nk2, wcc, Chern)
   Chern_all(4)= Chern
   wcc_all(:, :, 4)= wcc
   largestgap_all(:, 4)= largestgap

   !> integration over ka (\bar{a}}, wcc along kb, fixed k3=0
   kstart= (/0.0d0, 0.0d0, 0.0d0/)
   kvec1 = (/1.0d0, 0.0d0, 0.0d0/)
   kvec2 = (/0.0d0, 1.0d0, 0.0d0/)

   call  wannier_center3D_plane_func(kstart, kvec1, kvec2, &
      largestgap, wcc, Chern)
   call get_chern_number_from_wcc(NumberofSelectedOccupiedBands, Nk2, wcc, Chern)
   Chern_all(5)= Chern
   wcc_all(:, :, 5)= wcc
   largestgap_all(:, 5)= largestgap

   !> integration over ka (\bar{a}}, wcc along kb, fixed k3=kc/2
   kstart= (/0.0d0, 0.0d0, 0.5d0/)
   kvec1 = (/1.0d0, 0.0d0, 0.0d0/)
   kvec2 = (/0.0d0, 1.0d0, 0.0d0/)

   call  wannier_center3D_plane_func(kstart, kvec1, kvec2, &
      largestgap, wcc, Chern)
   call get_chern_number_from_wcc(NumberofSelectedOccupiedBands, Nk2, wcc, Chern)
   Chern_all(6)= Chern
   wcc_all(:, :, 6)= wcc
   largestgap_all(:, 6)= largestgap

   outfileindex= outfileindex+ 1
   if (cpuid==0) then
      open(unit=outfileindex, file='wanniercenter3D_Chern.dat')
      do ik2=1, Nk2
         write(outfileindex, '(10000f16.8)') (ik2-1d0)/(Nk2-1), &
            (dmod(sum((wcc_all(:, ik2, j))), 1d0), j=1, 6)
      enddo
      close(outfileindex)
   endif


   outfileindex= outfileindex+ 1
   if (cpuid==0) then
      open(unit=outfileindex, file='wanniercenter3D_Chern.gnu')

      write(outfileindex,'(a)')' set encoding iso_8859_1 '
      write(outfileindex,'(a)')' set terminal  postscript enhanced color font ",18"'
      write(outfileindex,'(a)')' set output "wanniercenter3D_Chern.eps"'
      write(outfileindex,'(a)')' set size 0.6,1.0'
      write(outfileindex,'(a)')' set multiplot '
      write(outfileindex,'(a)')' unset key '
      write(outfileindex,'(a)')' set border lw 1 '
      write(outfileindex,'(3a)')' NOXTICS = "set format x ', "''" ,&
         '; unset xtics; unset xlabel"'
      write(outfileindex,'(3a)')' XTICS = "set xtics format ', "'%1.0f'", &
         '; set xtics 1.0   nomirror in offset 0, 0.3; set mxtics 5;"'
      write(outfileindex,'(3a)')' NOYTICS = "set format y '," '';", 'unset ylabel"'
      write(outfileindex,'(3a)')' YTICS = "set ytics format '," '%1.0f' ",&
         '1.0  nomirror in offset 0.7,0; set mytics 2;"'
      write(outfileindex,'(a)')' TMARGIN = "set tmargin at screen 0.96; set bmargin at screen 0.71"'
      write(outfileindex,'(a)')' MMARGIN = "set tmargin at screen 0.63; set bmargin at screen 0.38"'
      write(outfileindex,'(a)')' BMARGIN = "set tmargin at screen 0.30; set bmargin at screen 0.05"'
      write(outfileindex,'(a)')' LMARGIN = "set lmargin at screen 0.20; set rmargin at screen 0.45"'
      write(outfileindex,'(a)')' RMARGIN = "set lmargin at screen 0.50; set rmargin at screen 0.75"'
      write(outfileindex,'(a)')' TITLE = "offset 0, -0.7"'
      write(outfileindex,'(3a)')' LCOLOR = "rgb '," '#696969'",'"'
      write(outfileindex,'(a)')' '
      write(outfileindex,'(3a)')' POS = "at graph -0.23,1.0 font', " ',18' ",'"'
      write(outfileindex,'(3a)')' POS2 = "at graph -0.15,1.0 font'," ',18' ",'"'
      write(outfileindex,'(a)')' '
      write(outfileindex,'(a)')' set xrange [0:1]'
      write(outfileindex,'(a)')' set yrange [0:1]'
      write(outfileindex,'(a)')' @TMARGIN; @LMARGIN'
      write(outfileindex,'(a)')' @XTICS; @YTICS'
      write(outfileindex,'(a)')' #set title "k_1=0.0" @TITLE'
      write(outfileindex,'(a)')' set xlabel "k_2" offset 0,1.6'
      write(outfileindex,'(a)')' set ylabel "c" rotate by  0 offset 1.8,0'
      write(outfileindex,'(a)')' set label 1 "(a)"  @POS front '
      write(outfileindex,'(a)')' plot "wanniercenter3D_Chern.dat" u 1:2 w p  pt 7  ps 0.6 lc @LCOLOR'
      write(outfileindex,'(a)')' '
      write(outfileindex,'(a)')' @TMARGIN; @RMARGIN'
      write(outfileindex,'(a)')' @NOYTICS; @XTICS '
      write(outfileindex,'(a)')' #set title "k_1=0.5" @TITLE'
      write(outfileindex,'(a)')' set label 1 "(b)" @POS2 front'
      write(outfileindex,'(a)')' set xlabel "k_2" offset 0,1.6'
      write(outfileindex,'(a)')' unset ylabel'
      write(outfileindex,'(a)')' plot "wanniercenter3D_Chern.dat" u 1:3 w p  pt 7  ps 0.6 lc @LCOLOR'
      write(outfileindex,'(a)')' '
      write(outfileindex,'(a)')' @MMARGIN; @LMARGIN'
      write(outfileindex,'(a)')' @YTICS; @XTICS '
      write(outfileindex,'(a)')' #set title "k_2=0.0" @TITLE'
      write(outfileindex,'(a)')' set label 1 "(c)" @POS front'
      write(outfileindex,'(a)')' set xlabel "k_1" offset 0,1.6'
      write(outfileindex,'(a)')' set ylabel "c" rotate by  0 offset 1.8,0'
      write(outfileindex,'(a)')' plot "wanniercenter3D_Chern.dat" u 1:4 w p  pt 7  ps 0.6 lc @LCOLOR'
      write(outfileindex,'(a)')' '
      write(outfileindex,'(a)')' @MMARGIN; @RMARGIN'
      write(outfileindex,'(a)')' @NOYTICS; @XTICS '
      write(outfileindex,'(a)')' set label 1 "(d)" @POS2 front'
      write(outfileindex,'(a)')' #set title "k_2=0.5" @TITLE'
      write(outfileindex,'(a)')' set xlabel "k_1" offset 0,1.6'
      write(outfileindex,'(a)')' plot "wanniercenter3D_Chern.dat" u 1:5 w p  pt 7  ps 0.6 lc @LCOLOR'
      write(outfileindex,'(a)')''
      write(outfileindex,'(a)')' @BMARGIN; @LMARGIN'
      write(outfileindex,'(a)')' @YTICS; @XTICS '
      write(outfileindex,'(a)')' #set title "k_3=0.0" @TITLE'
      write(outfileindex,'(a)')' set label 1 "(e)" @POS front'
      write(outfileindex,'(a)')' set ylabel "a" rotate by  0 offset 1.8,0'
      write(outfileindex,'(a)')' set xlabel "k_2" offset 0,1.6'
      write(outfileindex,'(a)')' plot "wanniercenter3D_Chern.dat" u 1:6 w p  pt 7  ps 0.6 lc @LCOLOR'
      write(outfileindex,'(a)')' '
      write(outfileindex,'(a)')' '
      write(outfileindex,'(a)')' @BMARGIN; @RMARGIN'
      write(outfileindex,'(a)')' @NOYTICS; @XTICS '
      write(outfileindex,'(a)')' #set title "k_3=0.5" @TITLE'
      write(outfileindex,'(a)')' set label 1 "(f)" @POS2 front'
      write(outfileindex,'(a)')' set xlabel "k_2" offset 0,1.6'
      write(outfileindex,'(a)')' plot "wanniercenter3D_Chern.dat" u 1:7 w p  pt 7  ps 0.6 lc @LCOLOR'

      close(outfileindex)
   endif

   if (cpuid==0) then
      write(stdout, *)'# Chern number for 6 planes'
      write(stdout, *)'k1=0.0, k2-k3 plane: ', int(Chern_all(1))
      write(stdout, *)'k1=0.5, k2-k3 plane: ', int(Chern_all(2))
      write(stdout, *)'k2=0.0, k1-k3 plane: ', int(Chern_all(3))
      write(stdout, *)'k2=0.5, k1-k3 plane: ', int(Chern_all(4))
      write(stdout, *)'k3=0.0, k1-k2 plane: ', int(Chern_all(5))
      write(stdout, *)'k3=0.5, k1-k2 plane: ', int(Chern_all(6))
   endif

   return
end subroutine  Chern_3D

subroutine get_chern_number_from_wcc(Noccpied, Numk, wcc, Chern_number)
   !> get chern number from the given wcc
   use para
   implicit none

   integer, intent(in) :: Noccpied
   integer, intent(in) :: Numk
   real(dp), intent(out) :: Chern_number
   real(dp), intent(in) :: wcc(Noccpied, Numk)

   integer :: ik
   real(dp) :: diff
   real(dp), allocatable :: wcc_sum(:)

   allocate(wcc_sum(Numk))
   wcc_sum= 0d0

   wcc_sum= dmod(sum(wcc, dim=1), 1d0)

   Chern_number= 0d0
   do ik=1, Numk-1
      call periodic_diff_1D(wcc_sum(ik+1), wcc_sum(ik), diff )
      Chern_number= Chern_number+ diff
   enddo

   return
end subroutine get_chern_number_from_wcc

subroutine  Z2_3D
   ! this suboutine is used for wannier center calculation for 3D system
   use para
   use wmpi
   implicit none

   integer :: ik2, i, j

   real(dp) :: kstart(3)
   real(dp) :: kvec1(3)
   real(dp) :: kvec2(3)

   real(dp) :: Z2
   real(dp) :: Z2_all(6)

   !> wannier centers for each ky, bands
   real(dp), allocatable :: wcc(:, :)
   real(dp), allocatable :: wcc_all(:, :, :)

   !> larges gap of each two wannier centers for a given k point
   !> dim= Nky
   real(dp), allocatable :: largestgap(:)
   real(dp), allocatable :: largestgap_all(:,:)

   allocate(wcc(NumberofSelectedOccupiedBands, Nk2))
   allocate(wcc_all(NumberofSelectedOccupiedBands, Nk2, 6))
   allocate(largestgap(Nk2))
   allocate(largestgap_all(Nk2,6))
   largestgap= 0d0
   wcc= 0d0
   wcc_all= 0d0

   !> integration over kc (\bar{c}}, wcc along kb, fixed k1=0
   kstart= (/0.0d0, 0.0d0, 0.0d0/)
   kvec1 = (/0.0d0, 0.0d0, 1.0d0/)
   kvec2 = (/0.0d0, 0.5d0, 0.0d0/)

   call  wannier_center3D_plane_func(kstart, kvec1, kvec2, &
      largestgap, wcc, Z2)
   Z2_all(1)= Z2
   wcc_all(:, :, 1)= wcc
   largestgap_all(:, 1)= largestgap

   !> integration over kc (\bar{c}}, wcc along kb, fixed k1=ka/2
   kstart= (/0.5d0, 0.0d0, 0.0d0/)
   kvec1 = (/0.0d0, 0.0d0, 1.0d0/)
   kvec2 = (/0.0d0, 0.5d0, 0.0d0/)

   call  wannier_center3D_plane_func(kstart, kvec1, kvec2, &
      largestgap, wcc, Z2)
   Z2_all(2)= Z2
   wcc_all(:, :, 2)= wcc
   largestgap_all(:, 2)= largestgap


   !> integration over kc (\bar{c}}, wcc along ka, fixed k2=0
   kstart= (/0.0d0, 0.0d0, 0.0d0/)
   kvec1 = (/0.0d0, 0.0d0, 1.0d0/)
   kvec2 = (/0.5d0, 0.0d0, 0.0d0/)

   call  wannier_center3D_plane_func(kstart, kvec1, kvec2, &
      largestgap, wcc, Z2)
   Z2_all(3)= Z2
   wcc_all(:, :, 3)= wcc
   largestgap_all(:, 3)= largestgap

   !> integration over kc (\bar{c}}, wcc along ka, fixed k2=kb/2
   kstart= (/0.0d0, 0.5d0, 0.0d0/)
   kvec1 = (/0.0d0, 0.0d0, 1.0d0/)
   kvec2 = (/0.5d0, 0.0d0, 0.0d0/)

   call  wannier_center3D_plane_func(kstart, kvec1, kvec2, &
      largestgap, wcc, Z2)
   Z2_all(4)= Z2
   wcc_all(:, :, 4)= wcc
   largestgap_all(:, 4)= largestgap

   !> integration over ka (\bar{a}}, wcc along kb, fixed k3=0
   kstart= (/0.0d0, 0.0d0, 0.0d0/)
   kvec1 = (/1.0d0, 0.0d0, 0.0d0/)
   kvec2 = (/0.0d0, 0.5d0, 0.0d0/)

   call  wannier_center3D_plane_func(kstart, kvec1, kvec2, &
      largestgap, wcc, Z2)
   Z2_all(5)= Z2
   wcc_all(:, :, 5)= wcc
   largestgap_all(:, 5)= largestgap

   !> integration over ka (\bar{a}}, wcc along kb, fixed k3=kc/2
   kstart= (/0.0d0, 0.0d0, 0.5d0/)
   kvec1 = (/1.0d0, 0.0d0, 0.0d0/)
   kvec2 = (/0.0d0, 0.5d0, 0.0d0/)

   call  wannier_center3D_plane_func(kstart, kvec1, kvec2, &
      largestgap, wcc, Z2)
   Z2_all(6)= Z2
   wcc_all(:, :, 6)= wcc
   largestgap_all(:, 6)= largestgap

   outfileindex= outfileindex+ 1
   if (cpuid==0) then
      open(unit=outfileindex, file='wanniercenter3D_Z2.dat')
      do i=1, NumberofSelectedOccupiedBands
         do ik2=1, Nk2
            write(outfileindex, '(10000f16.8)') (ik2-1d0)/(Nk2-1)/2d0, &
               (dmod((wcc_all(i, ik2, j)), 1d0), j=1, 6)
         enddo
         write(outfileindex, *)' '
      enddo
      close(outfileindex)
   endif


   outfileindex= outfileindex+ 1
   if (cpuid==0) then
      open(unit=outfileindex, file='wanniercenter3D_Z2.gnu')

      write(outfileindex,'(a)')' #gnuplot version>5.4 '
      write(outfileindex,'(a)')' set encoding iso_8859_1'
      write(outfileindex,'(a)')' set terminal pdf enhanced color font ",12" size 5,4'
      write(outfileindex,'(a)')' set output "wanniercenter3D_Z2.pdf"'
      write(outfileindex,'(a)')' set size 0.6,1.0'
      write(outfileindex,'(a)')' set multiplot '
      write(outfileindex,'(a)')' unset key '
      write(outfileindex,'(a)')' set border lw 1 '
      write(outfileindex,'(3a)')' NOXTICS = "set format x ', "''" ,&
         '; unset xtics; unset xlabel"'
      write(outfileindex,'(3a)')' XTICS = "set xtics format ', "'%4.1f'", &
         '; set xtics 0.5 nomirror in offset 0, 0.3; set mxtics 5;"'
      write(outfileindex,'(3a)')' NOYTICS = "set format y '," '';", 'unset ylabel"'
      write(outfileindex,'(3a)')' YTICS = "set ytics format '," '%1.0f' ",&
         '1.0  nomirror in offset 0.7,0; set mytics 2;"'
      write(outfileindex,'(a)')' TMARGIN = "set tmargin at screen 0.96; set bmargin at screen 0.72"'
      write(outfileindex,'(a)')' MMARGIN = "set tmargin at screen 0.63; set bmargin at screen 0.39"'
      write(outfileindex,'(a)')' BMARGIN = "set tmargin at screen 0.30; set bmargin at screen 0.06"'
      write(outfileindex,'(a)')' LMARGIN = "set lmargin at screen 0.20; set rmargin at screen 0.45"'
      write(outfileindex,'(a)')' RMARGIN = "set lmargin at screen 0.50; set rmargin at screen 0.75"'
      write(outfileindex,'(a)')' TITLE = "offset 0, -0.7"'
      write(outfileindex,'(3a)')' LCOLOR = "rgb '," '#696969'",'"'
      write(outfileindex,'(a)')' '
      write(outfileindex,'(3a)')' POS = "at graph -0.23,1.0 font', " ',12' ",'"'
      write(outfileindex,'(3a)')' POS2 = "at graph -0.15,1.0 font'," ',12' ",'"'
      write(outfileindex,'(a)')' '
      write(outfileindex,'(a)')' set xrange [0: 0.5]'
      write(outfileindex,'(a)')' set yrange [0:1]'
      write(outfileindex,'(a)')' @TMARGIN; @LMARGIN'
      write(outfileindex,'(a)')' @XTICS; @YTICS'
      write(outfileindex,'(a)')' #set title "k_1=0.0" @TITLE'
      write(outfileindex,'(a)')' set xlabel "k_2" offset 0,1.6'
      write(outfileindex,'(a)')' set ylabel "c" rotate by  0 offset 1.8,0'
      write(outfileindex,'(a)')' set label 1 "(a)"  @POS front '
      write(outfileindex,'(a)')' plot "wanniercenter3D_Z2.dat" u 1:2 w p  pt 7  ps 0.2 lc @LCOLOR'
      write(outfileindex,'(a)')' '
      write(outfileindex,'(a)')' @TMARGIN; @RMARGIN'
      write(outfileindex,'(a)')' @NOYTICS; @XTICS '
      write(outfileindex,'(a)')' #set title "k_1=0.5" @TITLE'
      write(outfileindex,'(a)')' set label 1 "(b)" @POS2 front'
      write(outfileindex,'(a)')' set xlabel "k_2" offset 0,1.6'
      write(outfileindex,'(a)')' unset ylabel'
      write(outfileindex,'(a)')' plot "wanniercenter3D_Z2.dat" u 1:3 w p  pt 7  ps 0.2 lc @LCOLOR'
      write(outfileindex,'(a)')' '
      write(outfileindex,'(a)')' @MMARGIN; @LMARGIN'
      write(outfileindex,'(a)')' @YTICS; @XTICS '
      write(outfileindex,'(a)')' #set title "k_2=0.0" @TITLE'
      write(outfileindex,'(a)')' set label 1 "(c)" @POS front'
      write(outfileindex,'(a)')' set xlabel "k_1" offset 0,1.6'
      write(outfileindex,'(a)')' set ylabel "c" rotate by  0 offset 1.8,0'
      write(outfileindex,'(a)')' plot "wanniercenter3D_Z2.dat" u 1:4 w p  pt 7  ps 0.2 lc @LCOLOR'
      write(outfileindex,'(a)')' '
      write(outfileindex,'(a)')' @MMARGIN; @RMARGIN'
      write(outfileindex,'(a)')' @NOYTICS; @XTICS '
      write(outfileindex,'(a)')' set label 1 "(d)" @POS2 front'
      write(outfileindex,'(a)')' #set title "k_2=0.5" @TITLE'
      write(outfileindex,'(a)')' set xlabel "k_1" offset 0,1.6'
      write(outfileindex,'(a)')' plot "wanniercenter3D_Z2.dat" u 1:5 w p  pt 7  ps 0.2 lc @LCOLOR'
      write(outfileindex,'(a)')''
      write(outfileindex,'(a)')' @BMARGIN; @LMARGIN'
      write(outfileindex,'(a)')' @YTICS; @XTICS '
      write(outfileindex,'(a)')' #set title "k_3=0.0" @TITLE'
      write(outfileindex,'(a)')' set label 1 "(e)" @POS front'
      write(outfileindex,'(a)')' set ylabel "a" rotate by  0 offset 1.8,0'
      write(outfileindex,'(a)')' set xlabel "k_2" offset 0,1.6'
      write(outfileindex,'(a)')' plot "wanniercenter3D_Z2.dat" u 1:6 w p  pt 7  ps 0.2 lc @LCOLOR'
      write(outfileindex,'(a)')' '
      write(outfileindex,'(a)')' '
      write(outfileindex,'(a)')' @BMARGIN; @RMARGIN'
      write(outfileindex,'(a)')' @NOYTICS; @XTICS '
      write(outfileindex,'(a)')' #set title "k_3=0.5" @TITLE'
      write(outfileindex,'(a)')' set label 1 "(f)" @POS2 front'
      write(outfileindex,'(a)')' set xlabel "k_2" offset 0,1.6'
      write(outfileindex,'(a)')' plot "wanniercenter3D_Z2.dat" u 1:7 w p  pt 7  ps 0.2 lc @LCOLOR'

      close(outfileindex)
   endif

   if (cpuid==0) then
      write(stdout, *)'# z2 number for 6 planes'
      write(stdout, *)'k1=0.0, k2-k3 plane: ', Z2_all(1)
      write(stdout, *)'k1=0.5, k2-k3 plane: ', Z2_all(2)
      write(stdout, *)'k2=0.0, k1-k3 plane: ', Z2_all(3)
      write(stdout, *)'k2=0.5, k1-k3 plane: ', Z2_all(4)
      write(stdout, *)'k3=0.0, k1-k2 plane: ', Z2_all(5)
      write(stdout, *)'k3=0.5, k1-k2 plane: ', Z2_all(6)
   endif

   return
end subroutine  Z2_3D


subroutine wilson_loop(Eigenvector,lkline,Nev,WProj,Nocc0,WEigen)
   use para
   use wmpi
   implicit none
!~    The subroutine calculating wilson loop by naive projector methods
!~    Input: eigenvectors "Eigenvector" in wannier_center2D
!~    Output: loop
!

   complex(dp),intent(in) :: Eigenvector(Nev,Nev,lkline)
   integer,intent(in) :: Nev,Nocc0,lkline
   complex(dp),intent(out) :: WProj(Nocc0,Nocc0)
   real(dp),intent(out) :: WEigen(Nocc0)

   integer :: ik1,i,j

   complex(dp), allocatable :: Hamk(:, :)
   complex(dp), allocatable :: Hamk_dag(:, :)

   complex(dp),allocatable :: Mmnkb(:,:)
   complex(dp), allocatable :: U(:, :)
   real   (dp), allocatable :: Sigma(:, :)
   complex(dp), allocatable :: VT(:, :)

   complex(dp), allocatable :: Lambda_eig(:)
   complex(dp), allocatable :: Lambda(:, :)
   complex(dp), allocatable :: Lambda0(:, :)

   complex(dp),allocatable :: mVL(:,:),mVR(:,:)
   complex(dp),allocatable :: wilham(:,:)

   allocate(Lambda_eig(Nocc0))
   allocate(Lambda(Nocc0, Nocc0))
   allocate(Lambda0(Nocc0, Nocc0))

!~    MmnKb= U Sigma VT - the SVD matrices
!~    Mmnkb : the multipliy of eigenvectors PI <u_k|u_k+1>
!> Mmnkb=<u_n(k)|u_m(k+b)>
   !> |u_n(k)> is the unitary part of wave function
   allocate(Mmnkb(Nocc0, Nocc0))
   !> three matrix for SVD
   !> M= U.Sigma.V^\dag
   !> VT= V^\dag
   allocate(U(Nocc0, Nocc0))
   allocate(Sigma(Nocc0, Nocc0))
   allocate(VT(Nocc0, Nocc0))

   allocate(hamk(Nev, Nev))
   allocate(hamk_dag(Nev, Nev))
   ! the eigenvalues


   allocate(mVL(Nocc0,Nocc0),mVR(Nocc0,Nocc0))
   ! the matrices to sort the eigenvalues of wilson loop hamiltonian
   allocate(wilham(Nocc0,Nocc0))
   wilham=0


   Lambda0=0d0
   do i=1, Nocc0
      Lambda0(i, i)= 1d0
   enddo
   lambda=0d0
   do ik1=1, lkline-1
      !  <u_k|u_k+1>
      Mmnkb= 0d0
      hamk_dag=Eigenvector(:, :, ik1)
      hamk_dag= conjg(transpose(hamk_dag))
      if (ik1==lkline) then
         hamk= Eigenvector(:, :, 1)
      else
         hamk= Eigenvector(:, :, ik1+ 1)
      endif
      Mmnkb=matmul(hamk_dag(Selected_Occupiedband_index(1):Selected_Occupiedband_index(1)+Nocc0,:), &
         hamk(:,Selected_Occupiedband_index(1):Selected_Occupiedband_index(1)+Nocc0))

      !> perform Singluar Value Decomposed of Mmnkb
      call zgesvd_pack(Nocc0, Mmnkb, U, Sigma, VT)

      !> after the calling of zgesvd_pack, Mmnkb becomes a temporal matrix
      U = conjg(transpose(U))
      VT= conjg(transpose(VT))
      call mat_mul(Nocc0, VT, U, Mmnkb)

!~       call zgeev_sys(Nocc0, Mmnkb, Lambda_eig,'V',mVL,'V',mVR)
!~       Mmnkb=matmul(mVR,conjg(transpose(mVL)))

      call mat_mul(Nocc0, Mmnkb, Lambda0, Lambda)
      Lambda0 = Lambda
   enddo  !< ik1

   !~       Lambda is then constructed as WilsonLoop Matrix
   wilham=Lambda
   !> diagonalize Lambda to get the eigenvalue
   call zgeev_sys(Nocc0, Lambda, Lambda_eig,'V',mVL,'V',mVR)
!~    call zgeev_pack(Nocc0, Lambda, Lambda_eig)
!~ return
   do i=1, Nocc0
      WEigen(i)= aimag(log(Lambda_eig(i)))/2d0/pi
      WEigen(i)= (mod(WEigen(i)+10d0, 1d0))
   enddo

   Lambda=0d0
   do i=1,Nocc0
      Lambda(i,i)=WEigen(i)!*2d0*pi
!~       if(abs(Lambda(i,i))>0.5) Lambda(i,i)=0d0

   enddo
   write(19888,*) (weigen(i),i=1,Nocc0)
   WProj=0d0

   wilham=matmul(mVR,conjg(transpose(mVL)))
!~    wilham=Lambda
!~    write(*,*)
   do i=1,Nocc0

      if(abs(real(wilham(i,i))-1.00)>0.000001) then
         write(*,*) (wilham(i,j),j=1,Nocc0)
         stop
      endif

   enddo

!~    do i=1,Nocc0
!~       WProj=WProj+matmul(mvR(:,i)*WEigen(i),mVL(:,i))!conjg(transpose(mVL)))
!~    enddo


   WProj=matmul(matmul(mVR,Lambda),conjg(transpose(mVR)))
!~    WProj=mVR

   return
end subroutine


