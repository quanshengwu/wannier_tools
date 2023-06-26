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

      integer :: i, j, it, ik, Nk_seg, NK_Berry_tot

      !> k points in kx-ky plane
      real(dp), allocatable :: kpoints(:, :)

      !> hamiltonian for each k point
      !> and also the eigenvector of hamiltonian after eigensystem_c
      complex(dp), allocatable :: uk(:, :), uk_dag(:, :)

      !> eigenvector for each kx
      complex(dp), allocatable :: Eigenvector(:, :, :)

      real   (dp), allocatable :: eigenvalue(:)
      complex(dp), allocatable :: phase(:)
      real(dp) :: br
      real(dp) :: k(3), b(3)
      complex(dp) :: overlap, ratio


      ! number of k points in each segment
      Nk_seg= Nk

      ! total number of k points in the loop for Berry phase calculation
      NK_Berry_tot= (NK_Berry-1)*Nk_seg 

      allocate(kpoints(3, NK_Berry_tot))
      kpoints= 0d0

      allocate(uk(Num_wann, Num_wann),  uk_dag(Num_wann, Num_wann))
      allocate(Eigenvector(Num_wann, Num_wann, NK_Berry_tot), eigenvalue(Num_wann))
      allocate(phase(Num_wann))
      
      uk=0d0
      eigenvalue=0d0
      Eigenvector=0d0

      !> set k path for berry phase calculation
      !> kpoints, k3points_Berry are in fractional/direct coordinates
      it = 0
      do ik=1, NK_Berry- 1
         do i= 1, Nk_seg
            it= it+ 1
            kpoints(:, it)= k3points_Berry(:, ik)+ &
               (k3points_Berry(:, ik+1)- k3points_Berry(:, ik))*(i-1d0)/(Nk_seg-1d0)
         enddo ! i
      enddo ! ik

      !> for each ky, we can get wanniercenter
      do ik=1, NK_Berry_tot
         k= kpoints(:, ik)

         if (index(KPorTB, 'KP')/=0)then
            call ham_bulk_kp (k, uk)
         else
            call ham_bulk_latticegauge(k, uk)
         endif
        
         !> diagonal uk
         call eigensystem_c('V', 'U', Num_wann, uk, eigenvalue)

         Eigenvector(:, :, ik)= uk
      enddo

      !> sum over k to get berry phase
      phase= 1d0
      do ik= 1, NK_Berry_tot-1
         uk= Eigenvector(:, :, ik)
         uk_dag= conjg(transpose(uk))
         if (ik==NK_Berry_tot-1) then
            uk= Eigenvector(:, :, 1)
         else
            uk= Eigenvector(:, :, ik+ 1)
         endif
         b= kpoints(:, ik+1)- kpoints(:, ik)

         !> <u_k|u_k+1>
         do i=1, Num_wann
            overlap= 0d0
            do j=1, Num_wann
               br= b(1)*Origin_cell%wannier_centers_direct(1, j)+ &
                   b(2)*Origin_cell%wannier_centers_direct(2, j)+ &
                   b(3)*Origin_cell%wannier_centers_direct(3, j)
               ratio= cos(2d0*pi*br)- zi*sin(2d0*pi*br)

               overlap= overlap+ uk_dag(i, j)* uk(j, i)* ratio
            enddo
            phase(i)= overlap*phase(i)
         enddo

      enddo  !< ik

      if (cpuid==0)write(stdout, *) ">> WARNING: Please increase NK1 until Berry phase is converged!"
      if (cpuid==0)write(stdout, *) ">> WARNING: The starting point and the ending point should be different by a reciprocal lattice vector"
      if (cpuid==0)write(stdout, *) 'Berry phase for the loop you chose: in unit of \pi'
      if (cpuid==0) write(stdout, '(f18.6)') mod(sum(aimag(log(phase(1:NumOccupied)))/pi), 2d0)

      outfileindex= outfileindex+ 1
      if (cpuid==0) then
         open(unit= outfileindex, file="kpath_berry.txt")
         write(outfileindex, '("#",a11, 5a12, a, f12.6, a)')"kx", "ky", "kz", &
                              "k1", "k2", "k3", " Berry phase= ", &
                              mod(sum(aimag(log(phase(1:NumOccupied)))/pi), 2d0), ' pi'
         do ik=1, NK_Berry_tot
            b= kpoints(1, ik)*Origin_cell%Kua+ kpoints(2, ik)*Origin_cell%Kub+ kpoints(3, ik)*Origin_cell%Kuc
            write(outfileindex, '(6f12.6)')b, kpoints(:, ik)
         enddo
      endif
      
      deallocate(kpoints)
      deallocate(uk,uk_dag)
      deallocate(Eigenvector,eigenvalue)
      deallocate(phase)
 
      return
   end subroutine berryphase


   subroutine  berryphase_atomic
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

      integer :: i, j, it, ik, Nk_seg, NK_Berry_tot

      !> k points in kx-ky plane
      real(dp), allocatable :: kpoints(:, :)

      !> hamiltonian for each k point
      !> and also the eigenvector of hamiltonian after eigensystem_c
      complex(dp), allocatable :: uk(:, :), uk_dag(:, :)

      !> eigenvector for each kx
      complex(dp), allocatable :: Eigenvector(:, :, :)

      real   (dp), allocatable :: eigenvalue(:)
      complex(dp), allocatable :: phase(:), mat1(:, :), mat2(:, :)
      real(dp) :: br
      real(dp) :: k(3), b(3)
      complex(dp) :: overlap, ratio


      ! number of k points in each segment
      Nk_seg= Nk

      ! total number of k points in the loop for Berry phase calculation
      NK_Berry_tot= (NK_Berry-1)*Nk_seg 

      allocate(kpoints(3, NK_Berry_tot))
      kpoints= 0d0

      allocate(mat1(Num_wann, Num_wann), mat2(Num_wann, Num_wann))
      allocate(uk(Num_wann, Num_wann),  uk_dag(Num_wann, Num_wann))
      allocate(Eigenvector(Num_wann, Num_wann, NK_Berry_tot), eigenvalue(Num_wann))
      allocate(phase(Num_wann))
      
      uk=0d0
      eigenvalue=0d0
      Eigenvector=0d0

      !> set k path for berry phase calculation
      !> kpoints, k3points_Berry are in fractional/direct coordinates
      it = 0
      do ik=1, NK_Berry- 1
         do i= 1, Nk_seg
            it= it+ 1
            kpoints(:, it)= k3points_Berry(:, ik)+ &
               (k3points_Berry(:, ik+1)- k3points_Berry(:, ik))*(i-1d0)/(Nk_seg-1d0)
         enddo ! i
      enddo ! ik

      !> for each ky, we can get wanniercenter
      do ik=1, NK_Berry_tot
         k= kpoints(:, ik)

         if (index(KPorTB, 'KP')/=0)then
            call ham_bulk_kp (k, uk)
         else
            call ham_bulk_atomicgauge(k, uk)
         endif
        
         !> diagonal uk
         call eigensystem_c('V', 'U', Num_wann, uk, eigenvalue)

         Eigenvector(:, :, ik)= uk
      enddo

      !> sum over k to get berry phase
      phase= 1d0
      do ik= 1, NK_Berry_tot-1
         uk= Eigenvector(:, :, ik)
         uk_dag= conjg(transpose(uk))
         if (ik==NK_Berry_tot-1) then
            b=-kpoints(:, 1)+ kpoints(:, NK_Berry_tot)
            uk= Eigenvector(:, :, 1)
         else
            uk= Eigenvector(:, :, ik+ 1)
         endif

         !> <u_k|u_k+1>
         do i=1, Num_wann
            overlap= 0d0
            do j=1, Num_wann
               if (ik==NK_Berry_tot-1) then
                  br= b(1)*Origin_cell%wannier_centers_direct(1, j)+ &
                      b(2)*Origin_cell%wannier_centers_direct(2, j)+ &
                      b(3)*Origin_cell%wannier_centers_direct(3, j)
                  ratio= cos(2d0*pi*br)- zi*sin(2d0*pi*br)
               else
                  ratio=1d0
               endif

               overlap= overlap+ uk_dag(i, j)* uk(j, i)* ratio
            enddo
            phase(i)= overlap*phase(i)
         enddo

      enddo  !< ik

      if (cpuid==0)write(stdout, *) ">> WARNING: Please increase NK1 until Berry phase is converged!"
      if (cpuid==0)write(stdout, *) ">> WARNING: The starting point and the ending point should be different by a reciprocal lattice vector"
      if (cpuid==0)write(stdout, *) 'Berry phase for the loop you chose: in unit of \pi'
      if (cpuid==0) write(stdout, '(f18.6)') mod(sum(aimag(log(phase(1:NumOccupied)))/pi), 2d0)

      outfileindex= outfileindex+ 1
      if (cpuid==0) then
         open(unit= outfileindex, file="kpath_berry.txt")
         write(outfileindex, '("#",a11, 5a12, a, f12.6, a)')"kx", "ky", "kz", &
                              "k1", "k2", "k3", " Berry phase= ", &
                              mod(sum(aimag(log(phase(1:NumOccupied)))/pi), 2d0), ' pi'
         do ik=1, NK_Berry_tot
            b= kpoints(1, ik)*Origin_cell%Kua+ kpoints(2, ik)*Origin_cell%Kub+ kpoints(3, ik)*Origin_cell%Kuc
            write(outfileindex, '(6f12.6)')b, kpoints(:, ik)
         enddo
      endif
      
      deallocate(kpoints)
      deallocate(mat1,mat2,uk,uk_dag)
      deallocate(Eigenvector,eigenvalue)
      deallocate(phase)
 
      return
   end subroutine berryphase_atomic


