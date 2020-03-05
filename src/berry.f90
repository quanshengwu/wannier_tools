   subroutine  berryphase
      !> Subroutine for calculating Berry phase for a giving path
      !
      ! Comments:
      !
      !          At present, you have to define the k path that you want 
      !          in kpoints
      !
      ! Author : QuanSheng Wu (wuquansheng@gmail.com)
      !
      ! 31 Mar 2016
      !
      ! Copyright (c) 2010 QuanSheng Wu. All rights reserved.

      use para
      use wmpi
      implicit none

      integer :: i, j, it, nfill, ik, Nk_seg, NK_Berry_tot

      !> k points in kx-ky plane
      real(dp), allocatable :: kpoints(:, :)

      !> hamiltonian for each k point
      !> and also the eigenvector of hamiltonian after eigensystem_c
      complex(dp), allocatable :: Hamk(:, :), Hamk1(:, :), Hamk_dag(:, :)

      !> eigenvector for each kx
      complex(dp), allocatable :: Eigenvector(:, :, :)

      !> Mmnkb=<u_n(k)|u_m(k+b)>
      !> |u_n(k)> is the periodic part of wave function
      complex(dp), allocatable :: Mmnkb(:, :), Mmnkb_com(:, :), Mmnkb_full(:, :)

      complex(dp), allocatable :: Lambda_eig(:), Lambda(:, :), Lambda0(:, :)

      !> three matrix for SVD 
      !> M= U.Sigma.V^\dag
      !> VT= V^\dag
      real   (dp), allocatable :: Sigma(:, :), eigenvalue(:)
      complex(dp), allocatable :: U(:, :), VT(:, :)
      complex(dp), allocatable :: phase(:), mat1(:, :), mat2(:, :)
      real(dp) :: kx, ky, kz, r1, r2, phi, br
      real(dp) :: k(3), b(3), k1(3)
      complex(dp) :: overlap, ratio

      nfill= Numoccupied

      ! number of k points in each segment
      Nk_seg= Nk

      ! total number of k points in the loop for Berry phase calculation
      NK_Berry_tot= (NK_Berry-1)*Nk_seg 

      allocate(kpoints(3, NK_Berry_tot))
      kpoints= 0d0

      allocate(Lambda_eig(nfill), Lambda(nfill, nfill), Lambda0(nfill, nfill))
      allocate(Mmnkb(nfill, nfill), Mmnkb_com(nfill, nfill), Mmnkb_full(Num_wann, Num_wann))
      allocate(mat1(Num_wann, Num_wann), mat2(Num_wann, Num_wann))
      allocate(hamk(Num_wann, Num_wann), hamk1(Num_wann, Num_wann), hamk_dag(Num_wann, Num_wann))
      allocate(Eigenvector(Num_wann, Num_wann, NK_Berry_tot), eigenvalue(Num_wann))
      allocate(U(nfill, nfill), Sigma(nfill, nfill), VT(nfill, nfill))
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

      r1=0.050d0
      r2=0.050d0

      !> define a circle k path
     !do ik=1, Nk
     !   phi= (ik-1d0)*2d0*pi/dble(Nk)
     !   kx= 0.400000d0+ cos(phi)*r1
     !   ky= 0.330000d0
     !   kz= 0.0d0+ sin(phi)*r2
     !   kpoints(1, ik)= kx
     !   kpoints(2, ik)= ky
     !   kpoints(3, ik)= kz
     !enddo

      !> define a circle k path
      ! for C8 Bin Wen
     !do ik=1, Nk
     !   phi= (ik-1d0)*2d0*pi/dble(Nk)
     !   kx= 0.200000d0+ cos(phi)*r1
     !   ky= -0.20000d0+ sin(phi)*r2
     !   kz= 0.0d0
     !   kpoints(1, ik)= kx
     !   kpoints(2, ik)= ky
     !   kpoints(3, ik)= kz
     !enddo

     
      !> set k path for berry phase calculation
      !> kpoints, k3points_Berry are in fractional/direct coordinates
      it = 0
      do ik=1, NK_Berry- 1
         do i= 1, Nk_seg
            it= it+ 1
            kpoints(:, it)= k3points_Berry(:, ik)+ &
               (k3points_Berry(:, ik+1)- k3points_Berry(:, ik))*(i-1)/Nk_seg
         enddo ! i
      enddo ! ik

      !> for each ky, we can get wanniercenter
      do ik=1, NK_Berry_tot
         k(1)= kpoints(1, ik)
         k(2)= kpoints(2, ik)
         k(3)= kpoints(3, ik)

         if (index(KPorTB, 'KP')/=0)then
            call ham_bulk_kp (k, Hamk)
         else
            call ham_bulk_latticegauge(k, Hamk)
         endif
        
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
      do ik= 1, NK_Berry_tot
         hamk= Eigenvector(:, :, ik)
         hamk_dag= conjg(transpose(hamk))
         if (ik==NK_Berry_tot) then
            b= kpoints(:, 1)- kpoints(:, NK_Berry_tot)
            hamk= Eigenvector(:, :, 1)
         else
            b= kpoints(:, ik+1)- kpoints(:, ik)
            hamk= Eigenvector(:, :, ik+ 1)
         endif
         b= b(1)*Origin_cell%Kua+b(2)*Origin_cell%Kub+b(3)*Origin_cell%Kuc

         !> <u_k|u_k+1>
         do i=1, Num_wann
            br= b(1)*Origin_cell%wannier_centers_cart(1, i)+ &
                b(2)*Origin_cell%wannier_centers_cart(2, i)+ &
                b(3)*Origin_cell%wannier_centers_cart(3, i)
            ratio= cos(br)- zi* sin(br)

            overlap= 0d0
            do j=1, Num_wann
               overlap= overlap+ hamk_dag(i, j)* hamk(j, i)* ratio
            enddo
           !phase(i)= phase(i)- aimag(log(overlap))/pi
            phase(i)= overlap*phase(i)
         enddo

      enddo  !< ik

      if (cpuid==0)write(stdout, *) 'Berry phase for the loop you chose: in unit of \pi'
      do i=1, Num_wann
        !if (cpuid==0)  write(*, '(10f10.5)')aimag(log(phase(i)))/pi
        !write(*, '(10f10.5)')phase(i)
      enddo
      if (cpuid==0) write(stdout, '(f18.6)')mod(sum(aimag(log(phase(1:nfill)))/pi), 2d0)

      outfileindex= outfileindex+ 1
      if (cpuid==0) then
         open(unit= outfileindex, file="kpath_berry.txt")
         write(outfileindex, '("#",a11, a12, a12, a, f12.6, a)')"kx", "ky", "kz", " Berry phase= ", &
                              mod(sum(aimag(log(phase(1:nfill)))/pi), 2d0), ' pi'
         do ik=1, NK_Berry_tot
            b= kpoints(1, ik)*Origin_cell%Kua+ kpoints(2, ik)*Origin_cell%Kub+ kpoints(3, ik)*Origin_cell%Kuc
            write(outfileindex, '(3f12.6)')b
         enddo
      endif
      
      deallocate(kpoints,Lambda_eig,Lambda,Lambda0)
      deallocate(Mmnkb,Mmnkb_com,Mmnkb_full)
      deallocate(mat1,mat2,hamk,hamk1,hamk_dag)
      deallocate(Eigenvector,eigenvalue)
      deallocate(U, Sigma, VT, phase)
 
      return
   end subroutine berryphase



   subroutine  nonabelian_berryphase
      !> Subroutine for calculating Berry phase for a giving path
      !
      ! Comments:
      !
      !          At present, you have to define the k path that you want 
      !          in kpoints
      !
      ! Author : QuanSheng Wu (wuquansheng@gmail.com)
      !
      ! 31 Mar 2016
      !
      ! Copyright (c) 2010 QuanSheng Wu. All rights reserved.

      use para
      use wmpi
      implicit none

      integer :: i, j, m, it, nfill, ik, Nk_seg, NK_Berry_tot

      !> k points in kx-ky plane
      real(dp), allocatable :: kpoints(:, :)

      !> hamiltonian for each k point
      !> and also the eigenvector of hamiltonian after eigensystem_c
      complex(dp), allocatable :: Hamk(:, :), Hamk1(:, :), Hamk_dag(:, :)

      !> eigenvector for each kx
      complex(dp), allocatable :: Eigenvector(:, :, :)

      !> Mmnkb=<u_n(k)|u_m(k+b)>
      !> |u_n(k)> is the periodic part of wave function
      complex(dp), allocatable :: Mmnkb(:, :), Mmnkb_com(:, :), Mmnkb_full(:, :)

      !> 
      complex(dp), allocatable :: Lambda_eig(:), Lambda(:, :), Lambda0(:, :)

      !> three matrix for SVD 
      !> M= U.Sigma.V^\dag
      !> VT= V^\dag
      real   (dp), allocatable :: Sigma(:, :), eigenvalue(:)
      complex(dp), allocatable :: U(:, :), VT(:, :)
      complex(dp), allocatable :: phase(:), mat1(:, :), mat2(:, :)
      
      real(dp) :: kx, ky, kz, r1, r2, phi, br
      real(dp) :: k(3), b(3), k1(3)
      complex(dp) :: overlap, ratio

      nfill= NumberofSelectedBands
      
      ! number of k points in each segment
      Nk_seg= Nk

      ! total number of k points in the loop for Berry phase calculation
      NK_Berry_tot= (NK_Berry-1)*Nk_seg 

      allocate(kpoints(3, NK_Berry_tot))
      kpoints= 0d0

      allocate(Lambda_eig(nfill), Lambda(nfill, nfill), Lambda0(nfill, nfill))
      allocate(Mmnkb(nfill, nfill), Mmnkb_com(nfill, nfill), Mmnkb_full(Num_wann, Num_wann))
      allocate(mat1(Num_wann, Num_wann), mat2(Num_wann, Num_wann))
      allocate(hamk(Num_wann, Num_wann), hamk1(Num_wann, Num_wann), hamk_dag(Num_wann, Num_wann))
      allocate(Eigenvector(Num_wann, Num_wann, NK_Berry_tot), eigenvalue(Num_wann))
      allocate(U(nfill, nfill), Sigma(nfill, nfill), VT(nfill, nfill))
      allocate(phase(Num_wann))
      hamk=0d0; hamk1=0d0
      eigenvalue=0d0; Eigenvector=0d0
      Mmnkb_full=0d0; Mmnkb= 0d0; Mmnkb_com=0d0
      Lambda =0d0; Lambda0=0d0
      U= 0d0; VT= 0d0; Sigma= 0d0

      r1=0.05d0; r2=0.05d0

      !> set k path for berry phase calculation from the input.dat or wt.in
      !> kpoints, k3points_Berry are in fractional/direct unit
      it = 0
      do ik=1, NK_Berry- 1
         do i= 1, Nk_seg
            it= it+ 1
            kpoints(:, it)= k3points_Berry(:, ik)+ &
               (k3points_Berry(:, ik+1)- k3points_Berry(:, ik))*(i-1)/Nk_seg
         enddo ! i
      enddo ! ik

      Lambda0=0d0
      do i=1, Numoccupied
         Lambda0(i, i)= 1d0
      enddo
      do ik=1, NK_Berry_tot
         k(1)= kpoints(1, ik)
         k(2)= kpoints(2, ik)
         k(3)= kpoints(3, ik)

         if (index(KPorTB, 'KP')/=0)then
            call ham_bulk_kp (k, Hamk)
         else
            call ham_bulk_latticegauge(k, Hamk)
         endif
        
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
      do ik= 1, NK_Berry_tot  
         hamk= Eigenvector(:, :, ik)
         hamk_dag= hamk
         Mmnkb= 0d0
         if (ik==NK_Berry_tot) then
            b= kpoints(:, 1)- kpoints(:, NK_Berry_tot)
            hamk= Eigenvector(:, :, 1)
         else
            b= kpoints(:, ik+1)- kpoints(:, ik)
            hamk= Eigenvector(:, :, ik+ 1)
         endif
         b= b(1)*Origin_cell%Kua+b(2)*Origin_cell%Kub+b(3)*Origin_cell%Kuc

         !> <u_k|u_k+1>
         do m=1, Num_wann
            br= b(1)*Origin_cell%wannier_centers_cart(1, m)+ &
                b(2)*Origin_cell%wannier_centers_cart(2, m)+ &
                b(3)*Origin_cell%wannier_centers_cart(3, m)
            ratio= cos(br)- zi* sin(br)

            do j=1, NumberofSelectedBands
               do i=1, NumberofSelectedBands
                  Mmnkb(i, j)=  Mmnkb(i, j)+ &
                     conjg(hamk_dag(m, Selected_band_index(i)))* hamk(m, Selected_band_index(j))* ratio
               enddo ! i
            enddo ! j
         enddo

         !> perform Singluar Value Decomposed of Mmnkb
        !call zgesvd_pack(NumberofSelectedBands, Mmnkb, U, Sigma, VT)

         !> after the calling of zgesvd_pack, Mmnkb becomes a temporal matrix
        !U = conjg(transpose(U))
        !VT= conjg(transpose(VT))
        !call mat_mul(NumberofSelectedBands, VT, U, Mmnkb)

         call mat_mul(NumberofSelectedBands, Mmnkb, Lambda0, Lambda)
         Lambda0 = Lambda
      enddo
      
      !> diagonalize Lambda to get the eigenvalue 
      call zgeev_pack(NumberofSelectedBands, Lambda, Lambda_eig)
      if (cpuid==0)write(stdout, *) 'Eigenvalues of W for the selected bands'
      if (cpuid==0)then
         do i=1, NumberofSelectedBands
            write(stdout, '(10f10.5)') Lambda_eig(i)
         enddo
      endif

      if (cpuid==0)write(stdout, *) 'Berry phase for the loop you chose: in unit of \pi'
      if (cpuid==0)write(stdout, '(10f10.5)')aimag(log(Lambda_eig(:)))/pi
      if (cpuid==0)write(stdout, *) ' '

      outfileindex= outfileindex+ 1
      if (cpuid==0) then
         open(unit= outfileindex, file="kpath_berry.txt")
         write(outfileindex, '("#",a11, a12, a12, a, f12.6, a)')"kx", "ky", "kz", " Berry phase= ", &
                              mod(sum(aimag(log(phase(1:nfill)))/pi), 2d0), ' pi'
         do ik=1, NK_Berry_tot
            b= kpoints(1, ik)*Origin_cell%Kua+ kpoints(2, ik)*Origin_cell%Kub+ kpoints(3, ik)*Origin_cell%Kuc
            write(outfileindex, '(3f12.6)')b
         enddo
      endif
      
      deallocate(kpoints,Lambda_eig,Lambda,Lambda0)
      deallocate(Mmnkb,Mmnkb_com,Mmnkb_full)
      deallocate(mat1,mat2,hamk,hamk1,hamk_dag)
      deallocate(Eigenvector,eigenvalue)
      deallocate(U, Sigma, VT, phase)
 
      return
   end subroutine nonabelian_berryphase
