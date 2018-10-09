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

      nfill= Numoccupied
      nfill_half= Numoccupied/2

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
            kpoints(:, ik1, ik2)= k0+k1*(ik1-1)/dble(nk1)+ k2*(ik2-1)/dble(nk2-1)
         enddo
      enddo
      b= k1/dble(nk1)
      b= b(1)*kua+b(2)*kub+b(3)*kuc


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
               br= b(1)*wannier_centers_cart(1, m )+ &
                   b(2)*wannier_centers_cart(2, m )+ &
                   b(3)*wannier_centers_cart(3, m )
               ratio= cos(br)- zi* sin(br)
        
               i1= 0
               do j=1, nfill
                  if (mirror_minus(j, ikp)) cycle
                  i1= i1+ 1
                  i2= 0
                  do i=1, nfill
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
               br= b(1)*wannier_centers_cart(1, m )+ &
                   b(2)*wannier_centers_cart(2, m )+ &
                   b(3)*wannier_centers_cart(3, m )
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
         write(outfileindex, '(a, i5, a)')"# for [i=3: ", Numoccupied/2+2, &
            "]  'wcc-mirrorplus.dat' u 1:i w p  pt 7  ps 1.1 lc 'red', \"  
         write(outfileindex, '(a, i5, a)')"# for [i=3: ", Numoccupied/2+2, &
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

   !> this suboutine is used for wannier center calculation for 3D system
   !> only for one plane
   subroutine  wannier_center3D_plane
      use para
      use wmpi
      implicit none

      integer :: i, j, m, nfill, imax

      integer :: ik1, ik2,  ierr

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

      real(dp) :: k(3), b(3)

      real(dp) :: maxgap
      real(dp) :: maxgap0

      !> Z2 calculation for time reversal invariant system
      integer :: Z2
      integer :: Delta

      real(dp) :: g, phi1, phi2, phi3
      real(dp) :: zm, zm1, xnm1, Deltam
      real(dp), allocatable :: xnm(:)
      real(dp) :: k0(3), k1(3), k2(3)



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
      !> the first dimension should be in one primitive cell, [0, 2*pi]
      !> the first dimension is the integration direction
      !> the WCCs are calculated along the second k line
      k0= K3D_start ! 
      k1= K3D_vec1   !  
      k2= K3D_vec2   ! 

      do ik2=1, Nk2
         do ik1=1, Nk1
            kpoints(:, ik1, ik2)= k0+ k1*(ik1-1d0)/dble(Nk1)+ k2*(ik2-1d0)/dble(Nk2-1)
         enddo
      enddo
      b= k1/dble(Nk1)
      b= b(1)*kua+b(2)*kub+b(3)*kuc

      Umatrix_t= transpose(Umatrix)
      call inv_r(3, Umatrix_t)

      !>> Get wannier center for ky=0 plane
      !> for each ky, we can get wanniercenter
      do ik2=1+ cpuid, Nk2, num_cpu
         if (cpuid.eq.0) write(stdout, *)' Wilson loop ',  'ik, Nk', ik2, nk2
         Lambda0=0d0
         do i=1, nfill
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
               br= b(1)*wannier_centers_cart(1, m)+ &
                   b(2)*wannier_centers_cart(2, m)+ &
                   b(3)*wannier_centers_cart(3, m)
               ratio= cos(br)- zi* sin(br)
         
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
         open(unit=outfileindex, file='wcc.dat')

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
         write(outfileindex, '(a)')'set xrange [0: 0.5]'
         write(outfileindex, '(a)')'set yrange [0:1]'
         write(outfileindex, '(a)')"plot 'wcc.dat' u 1:2 w l lw 2  lc 'blue', \"   
         write(outfileindex, '(a, i5, a)')" for [i=4: ", nfill+3, "] 'wcc.dat' u 1:i w p  pt 7  ps 1.1 lc 'red'"
         close(outfileindex)
      endif

      return
   end subroutine  wannier_center3D_plane


   subroutine  wannier_center3D_plane_func(kstart, kvec1, kvec2, largest_gap, wcc, Z2)
      !> this suboutine is used for wannier center calculation for 3D system
      !> only for one plane

      use para
      use wmpi
      implicit none

      integer :: i, j, m , imax
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
      b= b(1)*kua+b(2)*kub+b(3)*kuc

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
               br= b(1)*wannier_centers_cart(1, m)+ &
                   b(2)*wannier_centers_cart(2, m)+ &
                   b(3)*wannier_centers_cart(3, m)
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



! this suboutine is used for wannier center calculation for 3D system

   subroutine  wannier_center3D
      use para
      use wmpi
      implicit none

      integer :: i, j, m , nfill, imax

      integer :: ik1, ik2, ierr

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


      Umatrix_t= transpose(Umatrix)
      call inv_r(3, Umatrix_t)

      !>> Get wannier center for ky=0 plane
      !> for each ky, we can get wanniercenter
      do ik2=1+ cpuid, Nk2, num_cpu
         if (cpuid.eq.0) write(stdout, *)' Wilson loop',  'ik, Nk', ik2, nk2
         Lambda0=0d0
         do i=1, nfill
            Lambda0(i, i)= 1d0
         enddo

         !> for each k1, we get the eigenvectors
         do ik1=1, Nk1
            k(1)= kpoints(1, ik1, ik2)
            k(2)= 0d0
            k(3)= kpoints(2, ik1, ik2)


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
               br= b(1)*wannier_centers_cart(1, m )+ &
                   b(2)*wannier_centers_cart(2, m )+ &
                   b(3)*wannier_centers_cart(3, m )
               ratio= cos(br)- zi* sin(br)
         
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
         open(unit=outfileindex, file='wanniercenterky0.dat')

         do ik2=1, Nk2
            write(outfileindex, '(10000f16.8)') kpoints(2, 1, ik2), &
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
               br= b(1)*wannier_centers_cart(1, m )+ &
                   b(2)*wannier_centers_cart(2, m )+ &
                   b(3)*wannier_centers_cart(3, m )
               ratio= cos(br)- zi* sin(br)
         
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
         open(unit=outfileindex, file='wanniercenterky05.dat')

         do ik2=1, Nk2
            write(outfileindex, '(10000f16.8)') kpoints(2, 1, ik2), &
               largestgap_mpi(ik2), &
               dmod(sum(WannierCenterKy_mpi(:, ik2)), 1d0), & 
               WannierCenterKy_mpi(:, ik2)
         enddo
         close(outfileindex)
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
      implicit none
      real(8), intent(inout) :: a
      real(8), intent(inout) :: b
      real(8) :: c
      c=a
      a=b
      b=c
      return
   end subroutine swap


   subroutine  wannier_center3D_weyl
   !> this suboutine is used for wannier center calculation for 3D system
   !> for all weyl points specified in input.dat

      use para
      use wmpi
      implicit none

      integer :: i, j, ik2
      integer :: chirality

      character(40) :: epsfilename, wccfilename

      real(dp) :: k0(3)
      integer, allocatable :: chirality_all(:)
      real(dp), allocatable :: wcc(:, :)
      real(dp), allocatable :: wcc_all(:, :, :)
      real(dp), allocatable :: wcc_sum_all(:, :)

   
      allocate(chirality_all(Num_Weyls))
      allocate(wcc(Numoccupied, nk2))
      allocate(wcc_all(Numoccupied, nk2, Num_Weyls))
      allocate(wcc_sum_all(nk2, Num_Weyls))

      wcc= 0d0
      wcc_all= 0d0
      wcc_sum_all= 0d0

      do i= 1, Num_Weyls
         k0= weyl_position_direct(:, i)
         call wannier_center3D_weyl_func(k0, kr0, wcc, chirality)
         chirality_all(i)= chirality
         wcc_all(:, :, i) = wcc
         wcc_sum_all(:, :)= dmod(sum(wcc_all, dim=1), 1d0)
      enddo

      if (cpuid==0) write(stdout, '(a)')'Chiralities'
      if (cpuid==0) write(stdout, '(a, 3a9,a15)')'#', 'k1', 'k2', 'k3', 'Chirality'
      if (cpuid==0) then
         do i=1, Num_Weyls
            write(stdout, '(3f10.5,i8)')weyl_position_direct(:, i), chirality_all(i)
         enddo
      endif

      outfileindex= outfileindex+ 1
      if (cpuid==0) then
         open(unit=outfileindex, file='wanniercenter3D_Weyl.dat')
         write(outfileindex, '(a16, 10000i16)')'# Chirality', chirality_all
         write(outfileindex, '(10000a16)')'# k ', ('phase', j=1, Num_Weyls)
         do ik2=1, Nk2
            write(outfileindex, '(10000f16.8)') (ik2-1d0)/(Nk2-1), &
               (wcc_sum_all(ik2, j), j=1, Num_Weyls), ((wcc_all(i, ik2, j), i=1, Numoccupied), j=1, Num_Weyls)
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


   subroutine  wannier_center3D_weyl_func(k0, r0, wcc, chirality)
   !> this suboutine is used for wannier center calculation for 3D system
   !> only for one weyl point

      use para
      use wmpi
      implicit none

      integer :: i, j, m , ik1, ik2, ierr
      real(dp) :: r_para, theta1, theta2

      !> inout variables
      real(dp), intent(in) :: k0(3)
      real(dp), intent(in) :: r0
      integer, intent(out) :: chirality
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

      !> eigenvalue
      real(dp), allocatable :: eigenvalue(:)

      !> sumation for wannier charge center
      real(dp), allocatable :: wcc_sum(:)

      real(dp) :: Umatrix_t(3, 3)

      !> b.r
      real(dp) :: br

      !> exp(-i*b.R)
      complex(dp) :: ratio

      real(dp) :: k(3), b(3)
      real(dp) :: gap_sum, gap_step
      real(dp), allocatable :: xnm(:)


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
      allocate(wcc_sum(Nk2))
      allocate(xnm(Numoccupied))
      WannierCenterKy= 0d0; WannierCenterKy_mpi= 0d0
      hamk=0d0; eigenvalue=0d0; Eigenvector=0d0
      Mmnkb_full=0d0; Mmnkb=0d0; Mmnkb_com=0d0
      Lambda =0d0; Lambda0=0d0
      U= 0d0; Sigma= 0d0; VT= 0d0; wcc_sum= 0d0

      !> set k plane
      !> the first dimension should be in one primitive cell, [0, 2*pi]
      !> the first dimension is the integration direction
      !> the WCCs are calculated along the second k line
      do ik2=1, Nk2 ! the mesh for k line of that wcc are calculated
         theta2= (ik2- 1d0)/(Nk2- 1d0)* pi
         if (ik2== 1) theta2= (ik2- 1d0+ 0.10)/(Nk2- 1d0)* pi  ! avoid the North pole
         if (ik2== Nk2) theta2= (ik2- 1d0- 0.10)/(Nk2- 1d0)* pi  ! avoid the south pole
         do ik1=1, Nk1  ! the mesh along the integration direction
            theta1= (ik1- 1d0)/Nk1* 2d0* pi
            r_para= r0* sin(theta2)
            kpoints(1, ik1, ik2)= k0(1)+ r_para* cos (theta1)
            kpoints(2, ik1, ik2)= k0(2)+ r_para* sin (theta1)
            kpoints(3, ik1, ik2)= k0(3)+ kr0* cos (theta2)
         enddo
      enddo

      Umatrix_t= transpose(Umatrix)
      call inv_r(3, Umatrix_t)

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
            if (ik1<Nk1) then
               b= kpoints(:, ik1+1, ik2)- kpoints(:, ik1, ik2)
            else
               b= kpoints(:, 1, ik2)- kpoints(:, 2, ik2)
            endif
            b= b(1)*kua+b(2)*kub+b(3)*kuc

            !> <u_k|u_k+1>
            Mmnkb= 0d0
            hamk_dag= Eigenvector(:, :, ik1)
            if (ik1==Nk1) then
               hamk= Eigenvector(:, :, 1)
            else
               hamk= Eigenvector(:, :, ik1+ 1)
            endif
            do m=1, Num_wann
               br= b(1)*wannier_centers_cart(1, m)+ &
                   b(2)*wannier_centers_cart(2, m)+ &
                   b(3)*wannier_centers_cart(3, m)
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
   end subroutine  wannier_center3D_weyl_func


   subroutine  Chern_3D
      ! this suboutine is used for wannier center calculation for 3D system
      use para
      use wmpi
      implicit none

      integer :: ik2, i, j 

      real(dp) :: kstart(3), kvec1(3), kvec2(3)

      integer :: Chern, Chern_all(6)

      !> wannier centers for each ky, bands
      real(dp), allocatable :: wcc(:, :)
      real(dp), allocatable :: wcc_all(:, :, :)

      !> larges gap of each two wannier centers for a given k point
      !> dim= Nky
      real(dp), allocatable :: largestgap(:)
      real(dp), allocatable :: largestgap_all(:,:)

      !> sumation for chirality number
      real(dp), allocatable :: wcc_sum(:)

      real(dp) :: gap_sum, gap_step

      allocate(wcc_sum(Nk2))
      allocate(wcc(Numoccupied, Nk2))
      allocate(wcc_all(Numoccupied, Nk2, 6))
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
      
      Chern= nint(gap_sum)
      
      Chern_all(1)= Chern

      wcc_all(:, :, 1)= wcc
      largestgap_all(:, 1)= largestgap

      !> integration over kc (\bar{c}}, wcc along kb, fixed k1=ka/2
      kstart= (/0.5d0, 0.0d0, 0.0d0/)
      kvec1 = (/0.0d0, 0.0d0, 1.0d0/)
      kvec2 = (/0.0d0, 1.0d0, 0.0d0/)

      call  wannier_center3D_plane_func(kstart, kvec1, kvec2, &
         largestgap, wcc, Chern)

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
      
      Chern= nint(gap_sum)
      

      Chern_all(2)= Chern
      wcc_all(:, :, 2)= wcc
      largestgap_all(:, 2)= largestgap


      !> integration over kc (\bar{c}}, wcc along ka, fixed k2=0
      kstart= (/0.0d0, 0.0d0, 0.0d0/)
      kvec1 = (/0.0d0, 0.0d0, 1.0d0/)
      kvec2 = (/1.0d0, 0.0d0, 0.0d0/)

      call  wannier_center3D_plane_func(kstart, kvec1, kvec2, &
         largestgap, wcc, Chern)

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
      
      Chern= nint(gap_sum)
      
      Chern_all(3)= Chern
      wcc_all(:, :, 3)= wcc
      largestgap_all(:, 3)= largestgap

      !> integration over kc (\bar{c}}, wcc along ka, fixed k2=kb/2
      kstart= (/0.0d0, 0.5d0, 0.0d0/)
      kvec1 = (/0.0d0, 0.0d0, 1.0d0/)
      kvec2 = (/1.0d0, 0.0d0, 0.0d0/)

      call  wannier_center3D_plane_func(kstart, kvec1, kvec2, &
         largestgap, wcc, Chern)

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
      
      Chern= nint(gap_sum)
      
      Chern_all(4)= Chern
      wcc_all(:, :, 4)= wcc
      largestgap_all(:, 4)= largestgap

      !> integration over ka (\bar{a}}, wcc along kb, fixed k3=0
      kstart= (/0.0d0, 0.0d0, 0.0d0/)
      kvec1 = (/1.0d0, 0.0d0, 0.0d0/)
      kvec2 = (/0.0d0, 1.0d0, 0.0d0/)

      call  wannier_center3D_plane_func(kstart, kvec1, kvec2, &
         largestgap, wcc, Chern)

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
      
      Chern= nint(gap_sum)
      
      Chern_all(5)= Chern
      wcc_all(:, :, 5)= wcc
      largestgap_all(:, 5)= largestgap

      !> integration over ka (\bar{a}}, wcc along kb, fixed k3=kc/2
      kstart= (/0.0d0, 0.0d0, 0.5d0/)
      kvec1 = (/1.0d0, 0.0d0, 0.0d0/)
      kvec2 = (/0.0d0, 1.0d0, 0.0d0/)

      call  wannier_center3D_plane_func(kstart, kvec1, kvec2, &
         largestgap, wcc, Chern)

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
      
      Chern= nint(gap_sum)
      
      Chern_all(6)= Chern
      wcc_all(:, :, 6)= wcc
      largestgap_all(:, 6)= largestgap

      outfileindex= outfileindex+ 1
      if (cpuid==0) then
         open(unit=outfileindex, file='wanniercenter3D_Chern.dat')
         do i=1, Numoccupied
            do ik2=1, Nk2
               write(outfileindex, '(10000f16.8)') (ik2-1d0)/(Nk2-1), &
                  (dmod((wcc_all(i, ik2, j)), 1d0), j=1, 6)
            enddo
            write(outfileindex, *)' '
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
         write(stdout, *)'k1=0.0, k2-k3 plane: ', Chern_all(1)
         write(stdout, *)'k1=0.5, k2-k3 plane: ', Chern_all(2)
         write(stdout, *)'k2=0.0, k1-k3 plane: ', Chern_all(3)
         write(stdout, *)'k2=0.5, k1-k3 plane: ', Chern_all(4)
         write(stdout, *)'k3=0.0, k1-k2 plane: ', Chern_all(5)
         write(stdout, *)'k3=0.5, k1-k2 plane: ', Chern_all(6)
      endif

      return
   end subroutine  Chern_3D

