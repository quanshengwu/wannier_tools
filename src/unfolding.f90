
subroutine unfolding_kpath
   !> Unfold the energy bands of the supercell to a specified unit cell.
   !> Calculate unfolded band at (k, \omega)
   !> 
   !> Implemented by QSWU 2019
   use para
   use sparse
   use me_calculate

   implicit none

   real(dp), allocatable :: omega(:), eta_array(:)

   real(dp), allocatable :: spectrum_unfold(:, :, :, :), spectrum_unfold_mpi(:, :, :, :)
   !> unfolded spectrum, dimension (omeganum_unfold, NumberofEta, NumberofSelectedOrbitals_groups, nk3_band)

   real(dp) :: k_PBZ_direct(3), k_PBZ_direct_in_SBZ(3), k_cart(3), k_SBZ_direct(3)

   integer :: nnzmax, nnz
   integer, allocatable :: jcoo(:), icoo(:)
   complex(dp), allocatable :: acoo(:)
   !> Hamiltonian for sparse version

   complex(dp), allocatable :: hamk_bulk(:, :)
   !> Hamiltonian for dense version

   real(dp), allocatable :: W(:)

   complex(dp), allocatable :: psi(:), zeigv(:, :)
   !> eigenvector of the sparse matrix acoo. Dim=(ndim, neval)

   integer :: neval
   !number of ARPACK eigenvalues

   integer :: nvecs
   ! number of Arnoldi vectors

   logical :: ritzvec
   !> calculate eigenvector or not

   complex(dp) :: sigma=(0d0,0d0)
   !> shift-invert sigma

   integer :: n, i, ie, ik, ig, iq, io, ierr, ieta,NumberofEta, Ndimq
  
   !> selected atoms and orbitals
   integer :: NumberofSelectedOrbitals_groups_local
   integer, allocatable :: NumberofSelectedOrbitals_local(:)
   type(int_array1D), allocatable :: Selected_WannierOrbitals_local(:)

   ! time measurement
   real(dp) :: time1, time2, time3

   real(dp), allocatable :: weight(:)
   real(dp), external :: delta
   NumberofEta = 9


   !> set up selected atoms and orbitals
   if (Landaulevel_unfold_line_calc) then
      !> in the magnetic supercell, number of selected orbitals should be Magq times.
      NumberofSelectedOrbitals_groups_local= NumberofSelectedOrbitals_groups
      allocate(NumberofSelectedOrbitals_local(NumberofSelectedOrbitals_groups_local))
      allocate(Selected_WannierOrbitals_local(NumberofSelectedOrbitals_groups_local))
      do ig=1, NumberofSelectedOrbitals_groups_local
         NumberofSelectedOrbitals_local(ig)= NumberofSelectedOrbitals(ig)*Magq
         allocate(Selected_WannierOrbitals_local(ig)%iarray(NumberofSelectedOrbitals_local(ig)))
         do iq=1, Magq
            Selected_WannierOrbitals_local(ig)%iarray( &
               (iq-1)*NumberofSelectedOrbitals(ig)+1:iq*NumberofSelectedOrbitals(ig)) &
               = Selected_WannierOrbitals(ig)%iarray+(iq-1)*NumberofSelectedOrbitals(ig)   
         enddo
      enddo
   else
      NumberofSelectedOrbitals_groups_local= NumberofSelectedOrbitals_groups
      allocate(NumberofSelectedOrbitals_local(NumberofSelectedOrbitals_groups_local))
      allocate(Selected_WannierOrbitals_local(NumberofSelectedOrbitals_groups_local))
      do ig=1, NumberofSelectedOrbitals_groups_local
         NumberofSelectedOrbitals_local(ig)= NumberofSelectedOrbitals(ig)
         allocate(Selected_WannierOrbitals_local(ig)%iarray(NumberofSelectedOrbitals_local(ig)))
         Selected_WannierOrbitals_local(ig)%iarray= Selected_WannierOrbitals(ig)%iarray
      enddo
   endif

  

   sigma=(1d0,0d0)*iso_energy
   if (Is_Sparse_Hr) then
      if (NumSelectedEigenVals==0) NumSelectedEigenVals=OmegaNum
      neval=NumSelectedEigenVals
      if (neval>=Num_wann) neval= Num_wann- 2
 
      !> ncv
      nvecs=int(2*neval)
      if (nvecs<50) nvecs= 50
      if (nvecs>Num_wann) nvecs= Num_wann
   
   
      nnzmax=splen+Num_wann
      nnz=splen
      allocate( acoo(nnzmax))
      allocate( jcoo(nnzmax))
      allocate( icoo(nnzmax))
   else
      if (Landaulevel_unfold_line_calc.and..not.Is_Sparse_Hr) then
         neval= Num_wann*Magq
         nvecs= Num_wann*Magq
         allocate( hamk_bulk(num_wann*Magq, num_wann*Magq))
      else
         neval= Num_wann
         nvecs= Num_wann
         allocate( hamk_bulk(num_wann, num_wann))
      endif
   endif
   allocate( W( neval))
   allocate( psi(Num_wann))
   allocate( zeigv(Num_wann,nvecs))
   allocate( weight(NumberofSelectedOrbitals_groups_local))
   allocate( spectrum_unfold(omeganum_unfold, NumberofEta, NumberofSelectedOrbitals_groups_local, nk3_band)) 
   allocate( spectrum_unfold_mpi(omeganum_unfold, NumberofEta, NumberofSelectedOrbitals_groups_local, nk3_band)) 
   spectrum_unfold= 0d0
   spectrum_unfold_mpi= 0d0

   NumberofEta=9
   allocate(eta_array(NumberofEta))
   eta_array=(/0.1d0, 0.2d0, 0.4d0, 0.8d0, 1.0d0, 2d0, 4d0, 8d0, 10d0/)
   eta_array= eta_array*Fermi_broadening

   allocate(omega(omeganum_unfold))
   omega= 0d0
   do i= 1, omeganum_unfold
      omega(i)=omegamin+(i-1)*(omegamax-omegamin)/dble(omeganum_unfold)
   enddo

   !> first unfold the kpoints from the kpath of the supercell
   do ik= 1+cpuid, nk3_band, num_cpu
      if (cpuid==0) write(stdout, '(a, i10," /", i10)') 'BulkBand unfolding at :', ik, nk3_band
      k_PBZ_direct= kpath_3d(:, ik)
      call direct_cart_rec_unfold(k_PBZ_direct, k_cart)
      if (Landaulevel_unfold_line_calc) then 
         call cart_direct_rec_magneticcell(k_cart, k_PBZ_direct_in_SBZ)
      else
         call cart_direct_rec(k_cart, k_PBZ_direct_in_SBZ)
      endif

      k_SBZ_direct= k_PBZ_direct_in_SBZ- floor(k_PBZ_direct_in_SBZ)

      zeigv=0d0
      if (Landaulevel_unfold_line_calc) then 
         Bx=-2d0*pi*Magp/Magq; By=0d0
         Ndimq= Num_wann*Magq
         if (NumSelectedEigenVals==0) NumSelectedEigenVals=OmegaNum
         neval=NumSelectedEigenVals
         if (neval>=Ndimq) neval= Ndimq- 2
         nvecs=int(2*neval)
         if (nvecs<50) nvecs= 50
         if (nvecs>Ndimq) nvecs= Ndimq
         
         sigma=(1d0,0d0)*iso_energy

         if (allocated(zeigv)) deallocate(zeigv)
         if (allocated(psi)) deallocate(psi)
         if (allocated(W)) deallocate(W)
         allocate( zeigv(Ndimq,nvecs))
         allocate( psi(Ndimq))
         allocate( W(neval))
         if (Is_Sparse_Hr) then
            nnzmax=splen*Magq+Ndimq
            if (allocated(acoo)) deallocate(acoo)
            if (allocated(icoo)) deallocate(icoo)
            if (allocated(jcoo)) deallocate(jcoo)
            allocate( acoo(nnzmax))
            allocate( jcoo(nnzmax))
            allocate( icoo(nnzmax))
            nnz=nnzmax
            call ham_3Dlandau_sparseHR(nnz, Ndimq, Magq, k_SBZ_direct, acoo,jcoo,icoo)
            
            !> diagonalization by call zheev in lapack
            W= 0d0
            !> after arpack_sparse_coo_eigs, nnz will be updated.
            ritzvec= .true.
            call arpack_sparse_coo_eigs(Ndimq,nnzmax,nnz,acoo,jcoo,icoo,neval,nvecs,W,sigma, zeigv, ritzvec)
         else
            call ham_3Dlandau(Ndimq, Magq, k_SBZ_direct, hamk_bulk)
            zeigv=hamk_bulk
            call eigensystem_c('V', 'U', Ndimq ,zeigv, W)
         endif

         call now(time3)
      else
 
         if (Is_Sparse_Hr) then
            ! sparse hr
            call ham_bulk_coo_sparsehr(k_SBZ_direct,acoo,icoo,jcoo)
            nnz= splen
         
            !> diagonalization by call zheev in lapack
            W= 0d0
            !> after arpack_sparse_coo_eigs, nnz will be updated.
            ritzvec= .true.
            call arpack_sparse_coo_eigs(Num_wann,nnzmax,nnz,acoo,jcoo,icoo,neval,nvecs,W,sigma, zeigv, ritzvec)
            call now(time3)
         else
            ! dense hr
            call ham_bulk_atomicgauge(k_SBZ_direct, zeigv)
            
            ! diagonalize the Hamiltonian, zeigv is the Hamiltonian matrix before eigensystem_c calling.
            ! It will be replaced by the eigenvectors after eigensystem_c calling. 
            call eigensystem_c('V', 'U', Num_wann, zeigv, W)
         endif
      endif  ! unfold landaulevel or not



      do n= 1, neval
         psi= zeigv(:, n)
         k_cart_abs = sqrt(W(n) + photon_energy_arpes)
         if (Landaulevel_unfold_line_calc) then
            call get_projection_weight_bulk_unfold(Ndimq, k_SBZ_direct, k_PBZ_direct, psi, weight, Magnetic_cell, &
              NumberofSelectedOrbitals_groups_local, NumberofSelectedOrbitals_local, Selected_WannierOrbitals_local)
         else
            call get_projection_weight_bulk_unfold(num_wann, k_SBZ_direct, k_PBZ_direct, psi, weight, Origin_cell, &
              NumberofSelectedOrbitals_groups_local, NumberofSelectedOrbitals_local, Selected_WannierOrbitals_local)
         endif
         do ig=1, NumberofSelectedOrbitals_groups_local
            do ieta= 1, NumberofEta
               do ie=1, omeganum_unfold
                  spectrum_unfold(ie, ieta, ig, ik)= spectrum_unfold(ie, ieta, ig, ik) + &
                     weight(ig)*delta(eta_array(ieta), W(n)-omega(ie))
               enddo ! ie
            enddo ! ieta
         enddo ! ig
      enddo ! sum over n
   enddo !ik

#if defined (MPI)
   call mpi_allreduce(spectrum_unfold, spectrum_unfold_mpi,size(spectrum_unfold),&
      mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
   spectrum_unfold_mpi= spectrum_unfold
#endif

   outfileindex= outfileindex+ 1
   if (cpuid.eq.0)then
      open (unit=outfileindex, file='spectrum_unfold_kpath.dat')
      write(outfileindex, '("# Column", I5, 100I16)')(i, i=1, 11)
      write(outfileindex, '("#", a, 6X, 300f16.2)')'Brodening \eta (meV): ', Eta_array(:)*1000d0/eV2Hartree
      write(outfileindex, '("# ", a12, 3a16)')'k', ' E(eV)', 'A(k,E)'
      do ik=1, nk3_band
         do ie=1, omeganum_unfold
            write(outfileindex, '(300f16.8)')k3len_unfold(ik)*Angstrom2atomic, omega(ie)/eV2Hartree, &
               ((spectrum_unfold_mpi(ie, ieta, ig, ik), ieta=1, NumberofEta), ig=1, NumberofSelectedOrbitals_groups_local)
         enddo
         write(outfileindex, *) ' '
      enddo
      close(outfileindex)
      write(stdout,*)'<<< Unfold bands successfully'    
   endif
     
#if defined (MPI)
     call mpi_barrier(mpi_cmw, ierr)
#endif

   outfileindex= outfileindex+ 1
   if (cpuid==0) then
      open(unit=outfileindex, file='spectrum_unfold_kpath.gnu')
      write(outfileindex, '(a)') '#set terminal  postscript enhanced color font ",30"'
      write(outfileindex, '(a)')"#set output 'spectrum_unfold.eps'"
      write(outfileindex, '(a)') 'set terminal pngcairo enhanced color font ",60" size 1920,1680'
      write(outfileindex, '(a)') 'set palette defined ( 0 "white", 1  "#D72F01" )'
      write(outfileindex, '(a)')"set output 'spectrum_unfold_kpath.png'"
      write(outfileindex, '(a)')'set style data linespoints'
      write(outfileindex, '(a)')'set size 0.9, 1'
      write(outfileindex, '(a)')'set origin 0.05,0'
      write(outfileindex, '(a)')'unset key'
      write(outfileindex, '(a)')'set border lw 3'
      write(outfileindex, '(a)')'set view map'
      write(outfileindex, '(a)')'#set xtics font ",24"'
      write(outfileindex, '(a)')'#set ytics font ",24"'
      write(outfileindex, '(a)')'#set ylabel font ",24"'
      write(outfileindex, '(a)')'#set ylabel offset 1.5,0'
      write(outfileindex, '(a,f12.6)')'emin=', omegamin/eV2Hartree
      write(outfileindex, '(a,f12.6)')'emax=', omegamax/eV2Hartree
      write(outfileindex, '(a, f10.5, a)')'set xrange [0: ', maxval(k3len_unfold*Angstrom2atomic), ']'
      if (index(Particle,'phonon')/=0) then
         write(outfileindex, '(a, f10.5, a)')'set yrange [0: emax ]'
         write(outfileindex, '(a)')'set ylabel "Frequency (THz)"'
      else
         write(outfileindex, '(a)')'set ylabel "Energy (eV)"'
         write(outfileindex, '(a)')'set yrange [ emin : emax ]'
      endif
      write(outfileindex, 202, advance="no") (k3line_name(i), k3line_unfold_stop(i)*Angstrom2atomic, i=1, nk3lines)
      write(outfileindex, 203)k3line_name(nk3lines+1), k3line_unfold_stop(nk3lines+1)*Angstrom2atomic

      do i=1, nk3lines-1
         if (index(Particle,'phonon')/=0) then
            write(outfileindex, 204)k3line_unfold_stop(i+1)*Angstrom2atomic, '0.0', k3line_unfold_stop(i+1)*Angstrom2atomic, 'emax'
         else
            write(outfileindex, 204)k3line_unfold_stop(i+1)*Angstrom2atomic, 'emin', k3line_unfold_stop(i+1)*Angstrom2atomic, 'emax'
         endif
      enddo

      write(outfileindex, '(a)')"set colorbox"
      write(outfileindex, '(a)')'set pm3d interpolate 2,2'
      write(outfileindex, '(2a)')"splot 'spectrum_unfold_kpath.dat' u 1:2:(($6 )) ",  &
         " w pm3d"
      close(outfileindex)
   endif

202 format('set xtics (',20('"',A3,'" ',F10.5,','))
203 format(A3,'" ',F10.5,')')
204 format('set arrow from ',F10.5,',',A5,' to ',F10.5,',',A5, ' nohead front lw 3')

     
#if defined (MPI)
     call mpi_barrier(mpi_cmw, ierr)
#endif

   return
end subroutine unfolding_kpath


subroutine unfolding_kplane
   !> Unfold the energy bands of the supercell to a specified unit cell.
   !> Calculate unfolded band at (k1, k2) at a given energy iso_energy.
   !> Implemented by QSWU 2019
   use para
   use sparse
   use me_calculate
   implicit none

   real(dp), allocatable :: eta_array(:)

   real(dp), allocatable :: spectrum_unfold(:, :, :, :), spectrum_unfold_mpi(:, :, :, :)
   !> unfolded spectrum, dimension (omeganum_unfold, NumberofEta, NumberofSelectedOrbitals_groups, nk3_band)

   real(dp), allocatable :: qpi_unfold(:, :, :), qpi_unfold_mpi(:, :, :)
   !> For the QPI calculations, we only support the 2D systems

   real(dp) :: k_PBZ_direct(3), k_PBZ_direct_in_SBZ(3), k_cart(3), k_SBZ_direct(3)


   integer, allocatable :: jcoo(:), icoo(:)
   complex(dp), allocatable :: acoo(:)
   !> Hamiltonian for sparse version

   complex(dp), allocatable :: hamk_bulk(:, :)
   !> Hamiltonian for dense version

   real(dp), allocatable :: W(:)

   complex(dp), allocatable :: psi(:), zeigv(:, :)
   !> eigenvector of the sparse matrix acoo. Dim=(Num_wann, neval)

   integer :: neval
   !number of ARPACK eigenvalues

   integer :: nvecs
   ! number of Arnoldi vectors

   logical :: ritzvec
   !> calculate eigenvector or not

   complex(dp) :: sigma=(0d0,0d0)
   !> shift-invert sigma

   integer , allocatable :: ik12(:,:)
   real(dp), allocatable :: kxy(:,:), kxy_shape(:,:), kxy_plane(:,:)

   integer :: nnzmax, nnz, knv3, Ndimq, nkx, nky, Nk1_half, Nk2_half 
   integer :: iq, iq1, iq2, imin1, imin2, imax1, imax2
   integer :: n, i, j, ie, ik, io, ik1, ik2, ig, ierr, ieta,NumberofEta
    
   
   ! time measurement
   real(dp) :: time1, time2, time3

   real(dp), allocatable :: weight(:)
   real(dp), external :: delta
   NumberofEta = 9

   !> Nk1 and Nk2 should be odd number so that the center of the kslice is (0,0)
   !> if you want to calculate the QPI
   nkx=Nk1; nky=Nk2
   if (mod(Nk1, 2)==0) nkx= Nk1+1
   if (mod(Nk2, 2)==0) nky= Nk2+1

   Nk1_half= (nkx-1)/2
   Nk2_half= (nky-1)/2
 
   knv3= Nkx*Nky
   allocate( kxy(3, knv3))
   allocate( ik12(2, nkx*nky))
   allocate( kxy_shape(3, knv3))
   allocate( kxy_plane(3, knv3))
   kxy=0d0
   kxy_shape=0d0
   kxy_plane=0d0
   
   ik =0
   do i= 1, nkx
      do j= 1, nky
         ik =ik +1
         ik12(1, ik)= i
         ik12(2, ik)= j
         kxy(:, ik)= K3D_start+ K3D_vec1*(i-1)/dble(nkx-1)+ K3D_vec2*(j-1)/dble(nky-1) &
            -(K3D_vec1+K3D_vec2)/2d0
         kxy_shape(:, ik)= kxy(1, ik)* Folded_cell%Kua+ kxy(2, ik)* Folded_cell%Kub+ kxy(3, ik)* Folded_cell%Kuc 
         call rotate_k3_to_kplane(kxy_shape(:, ik), kxy_plane(:, ik))
      enddo
   enddo

   k_cart_abs = sqrt(iso_energy+ photon_energy_arpes)

   sigma=(1d0,0d0)*iso_energy
   if (Is_Sparse_Hr) then

      if (NumSelectedEigenVals==0) then
         if (OmegaNum==0) then
            NumSelectedEigenVals=Num_wann
         else
            NumSelectedEigenVals=OmegaNum
         endif
      endif
   
      neval=NumSelectedEigenVals
      if (neval>Num_wann-2) neval= Num_wann- 2

      !> ncv
      nvecs=int(2*neval)
      if (nvecs<50) nvecs= 50
      if (nvecs>Num_wann) nvecs= Num_wann
   
   
      nnzmax=splen+Num_wann
      nnz=splen
      allocate( acoo(nnzmax))
      allocate( jcoo(nnzmax))
      allocate( icoo(nnzmax))
   else
      neval= Num_wann
      nvecs= Num_wann
      allocate( hamk_bulk(num_wann, num_wann))
   endif
   allocate( W( neval))
   allocate( psi(Num_wann))
   allocate( zeigv(Num_wann,nvecs))
   allocate( weight(NumberofSelectedOrbitals_groups))
   allocate( qpi_unfold(NumberofEta, NumberofSelectedOrbitals_groups, knv3)) 
   allocate( qpi_unfold_mpi(NumberofEta, NumberofSelectedOrbitals_groups, knv3)) 
   allocate( spectrum_unfold(NumberofEta, NumberofSelectedOrbitals_groups, nkx, nky)) 
   allocate( spectrum_unfold_mpi(NumberofEta, NumberofSelectedOrbitals_groups, nkx, nky)) 
   psi= 0d0; zeigv= 0d0; weight= 0d0
   qpi_unfold= 0d0; qpi_unfold_mpi= 0d0
   spectrum_unfold= 0d0; spectrum_unfold_mpi= 0d0

   NumberofEta=9
   allocate(eta_array(NumberofEta))
   eta_array=(/0.1d0, 0.2d0, 0.4d0, 0.8d0, 1.0d0, 2d0, 4d0, 8d0, 10d0/)
   eta_array= eta_array*Fermi_broadening

   !> first unfold the kpoints from the kpath of the supercell
   do ik= 1+cpuid, knv3, num_cpu
      ik1= ik12(1, ik)
      ik2= ik12(2, ik)
      if (cpuid==0.and.mod(ik, 40)==1) write(stdout, '(a, i10," /", i10)') 'BulkBand unfolding at :', ik, knv3
      k_PBZ_direct= kxy(:, ik)
      call direct_cart_rec_unfold(k_PBZ_direct, k_cart)
      call cart_direct_rec(k_cart, k_PBZ_direct_in_SBZ)

      k_SBZ_direct= k_PBZ_direct_in_SBZ- floor(k_PBZ_direct_in_SBZ)

      if (Landaulevel_unfold_line_calc) then 
         stop " we don't support Landaulevel_unfold_line_calc for LandauLevel_kplane_calc calc"
         Bx=-2d0*pi*Magp/Magq; By=0d0
         Ndimq= Num_wann*Magq
         nnzmax= Num_wann*(2*ijmax+1)*Ndimq+Ndimq
         if (allocated(acoo)) deallocate(acoo)
         if (allocated(icoo)) deallocate(icoo)
         if (allocated(jcoo)) deallocate(jcoo)
         allocate( acoo(nnzmax))
         allocate( jcoo(nnzmax))
         allocate( icoo(nnzmax))
         nnz=nnzmax
         if (Is_Sparse_Hr) then
            call ham_3Dlandau_sparseHR(nnz, Ndimq, Magq, k_SBZ_direct, acoo,jcoo,icoo)
         else
            call ham_3Dlandau_sparse1(nnz, Ndimq, Magq, k_SBZ_direct, acoo,jcoo,icoo)
         endif

         !> diagonalization by call zheev in lapack
         W= 0d0
         !> after arpack_sparse_coo_eigs, nnz will be updated.
         ritzvec= .true.
         call arpack_sparse_coo_eigs(Num_wann,nnzmax,nnz,acoo,jcoo,icoo,neval,nvecs,W,sigma, zeigv, ritzvec)
         call now(time3)
      else
 
         if (Is_Sparse_Hr) then
            ! sparse hr
            call ham_bulk_coo_sparsehr(k_SBZ_direct,acoo,icoo,jcoo)
            nnz= splen
         
            !> diagonalization by call zheev in lapack
            W= 0d0
            !> after arpack_sparse_coo_eigs, nnz will be updated.
            ritzvec= .true.
            call arpack_sparse_coo_eigs(Num_wann,nnzmax,nnz,acoo,jcoo,icoo,neval,nvecs,W,sigma, zeigv, ritzvec)
            call now(time3)
         else
            ! dense hr
            call ham_bulk_atomicgauge(k_SBZ_direct, zeigv)
            
            ! diagonalize the Hamiltonian, zeigv is the Hamiltonian matrix before eigensystem_c calling.
            ! It will be replaced by the eigenvectors after eigensystem_c calling. 
            call eigensystem_c('V', 'U', Num_wann ,zeigv, W)
         endif
      endif  ! landaulevel or not


      do n= 1, neval
         psi= zeigv(:, n)
         
         call get_projection_weight_bulk_unfold(Num_wann, k_SBZ_direct, k_PBZ_direct, psi, weight, Origin_cell, &
              NumberofSelectedOrbitals_groups, NumberofSelectedOrbitals, Selected_WannierOrbitals)
         do ig=1, NumberofSelectedOrbitals_groups
            do ieta= 1, NumberofEta
               spectrum_unfold(ieta, ig, ik1, ik2)= spectrum_unfold(ieta, ig, ik1, ik2) + &
                  weight(ig)*delta(eta_array(ieta), W(n)-iso_energy)
            enddo ! ieta
         enddo ! ig
      enddo ! sum over n
   enddo !ik

#if defined (MPI)
   call mpi_allreduce(spectrum_unfold, spectrum_unfold_mpi,size(spectrum_unfold),&
      mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
   spectrum_unfold_mpi= spectrum_unfold
#endif

   outfileindex= outfileindex+ 1
   if (cpuid.eq.0)then
      open (unit=outfileindex, file='spectrum_unfold_kplane.dat')
      write(outfileindex, "('#', a8, 5a16, 3X, '| A(k,E)', a6, 100(8X,'group ', i2))")&
         'kx', 'ky', 'kz', 'kp1', 'kp2', 'kp3', 'total',&
         (i, i=1, NumberofSelectedOrbitals_groups)
      write(outfileindex, "('#column', i5, 3000i16)")(i, i=1, 6+NumberofSelectedOrbitals_groups*NumberofEta)
      write(outfileindex, '("#", a, 70X, 300f16.2)')'Brodening \eta (meV): ', Eta_array(:)*1000d0/eV2Hartree
      do ik=1, knv3
         ik1= ik12(1, ik)
         ik2= ik12(2, ik)
         write(outfileindex, '(3000E16.5)')kxy_shape(:, ik)*Angstrom2atomic, kxy_plane(:, ik)*Angstrom2atomic, &
              ((spectrum_unfold_mpi(ieta, ig, ik1, ik2), ieta=1, NumberofEta), ig=1, NumberofSelectedOrbitals_groups)
         if (mod(ik, nky)==0) write(outfileindex, *)' '
      enddo
      close(outfileindex)
      write(stdout,*)'<<< Unfold bands successfully'    
   endif
     
#if defined (MPI)
     call mpi_barrier(mpi_cmw, ierr)
#endif

   outfileindex= outfileindex+ 1
   if (cpuid==0) then
      open(unit=outfileindex, file='spectrum_unfold_kplane.gnu')
      write(outfileindex, '(a)') 'set terminal pngcairo enhanced color font ",60" size 1920,1680'
      write(outfileindex, '(a)') 'set palette defined ( 0 "white", 1  "#D72F01" )'
      write(outfileindex, '(a)')"set output 'spectrum_unfold_kplane.png'"
      write(outfileindex, '(a)')'set size 0.9, 1'
      write(outfileindex, '(a)')'set origin 0.05,0'
      write(outfileindex, '(a)')'set border lw 3'
      write(outfileindex, '(a)')'set pm3d'
      write(outfileindex, '(a)')'unset key'
      write(outfileindex, '(a)')'set view map'
      write(outfileindex, '(a)')'#set xtics font ",24"'
      write(outfileindex, '(a)')'#set ytics font ",24"'
      write(outfileindex, '(a)')'#set ylabel font ",24"'
      write(outfileindex, '(a)')'#set ylabel offset 1.5,0'
      write(outfileindex, '(a)')'set size ratio -1'
      write(outfileindex, '(a)')"set colorbox"
      write(outfileindex, '(a)')'set pm3d interpolate 2,2'
      write(outfileindex, '(2a)')"splot 'spectrum_unfold_kplane.dat' u 4:5:($11) ",  &
         " w pm3d"
      close(outfileindex)
   endif

   IF (QPI_unfold_plane_calc) then

     !> calculate QPI (jdos)
     do iq= 1+ cpuid, nkx*nky, num_cpu
        iq1= ik12(1, iq)- Nk1_half
        iq2= ik12(2, iq)- Nk2_half
        if (cpuid==0.and. mod(iq/num_cpu, 100)==0) &
           write(stdout, *) 'JDOS, iq ', iq, 'Nq',nkx*nky, 'time left', &
           (nkx*nky-iq)*(time3- time1)/num_cpu, ' s'
        call now(time1)
        imin1= max(-Nk1_half-iq1, -Nk1_half)+ Nk1_half+ 1
        imax1= min(Nk1_half-iq1, Nk1_half)+ Nk1_half+ 1
        imin2= max(-Nk2_half-iq2, -Nk2_half)+ Nk2_half+ 1
        imax2= min(Nk2_half-iq2, Nk2_half)+ Nk2_half+ 1
        do ik2= imin2, imax2
           do ik1= imin1, imax1
              do ig=1, NumberofSelectedOrbitals_groups
                 do ieta= 1, NumberofEta
                    qpi_unfold(ieta, ig, iq)= qpi_unfold(ieta, ig, iq)+ &
                       spectrum_unfold_mpi(ieta, ig, ik1, ik2)* spectrum_unfold_mpi(ieta, ig, ik1+iq1, ik2+iq2)
                 enddo !ieta
              enddo !ig
           enddo !ik1
        enddo !ik2

        call now(time3)
     enddo !iq

     qpi_unfold_mpi= 1D-12
#if defined (MPI)
     call mpi_reduce(qpi_unfold, qpi_unfold_mpi, size(qpi_unfold),mpi_double_precision,&
                     mpi_sum, 0, mpi_comm_world, ierr)
#else
     qpi_unfold_mpi= qpi_unfold
#endif
     qpi_unfold_mpi= qpi_unfold_mpi/nkx/nky

     outfileindex= outfileindex+ 1
     if (cpuid.eq.0.and.QPI_unfold_plane_calc)then
        write(stdout,*)'>>The calculation of joint density of state in the unfold model was done.'
        open (unit=outfileindex, file='qpi_unfold.dat')
        write(outfileindex,'(a)')'# Bulk joint density of states in the unfold mode'
        write(outfileindex,'(a)')'# The coordinate of k is redefined according to the KPLANE_BULK'
        write(outfileindex,'(a)')"# x axis is parallel to K1'"
        write(outfileindex,'(a)')"# z axis is parallel to K1'xK2'"
        write(outfileindex,'(a)')"# y axis is parallel to z x x"
        write(outfileindex,'(30a16)')'#kx', 'ky', 'log(dos)'
        write(outfileindex, "('#column', i5, 3000i12)")(i, i=1, 7+NumberofSelectedOrbitals_groups*NumberofEta)
        write(outfileindex, '("#", a, 6X, 300f16.2)')'Brodening \eta (meV): ', Eta_array(:)*1000d0/eV2Hartree
        write(outfileindex, "('#', a6, 6a12, 3X, '| A(k,E)', a6, 100(8X,'group ', i2))") &
           'kx', 'ky', 'kz', 'kp1', 'kp2', 'kp3' 
        do ik=1, nkx*nky
           write(outfileindex, '(3000E12.5)')kxy_shape(:, ik)*Angstrom2atomic, kxy_plane(:, ik)*Angstrom2atomic, &
              ((qpi_unfold_mpi(ieta, ig, ik), ieta=1, NumberofEta), ig=1, NumberofSelectedOrbitals_groups)
           if (mod(ik, nky)==0) write(outfileindex, *)' '
        enddo
        close(outfileindex)
        write(stdout,*)'<<< qpi_unfold finished!'    
     endif

     outfileindex= outfileindex+ 1
     if (cpuid==0) then
        open(unit=outfileindex, file='qpi_unfold.gnu')
        write(outfileindex, '(a)') 'set terminal pngcairo enhanced color font ",60" size 1920,1680'
        write(outfileindex, '(a)') 'set palette defined ( 0 "white", 1  "#D72F01" )'
        write(outfileindex, '(a)')"set output 'qpi_unfold.png'"
        write(outfileindex, '(a)')'set size 0.9, 1'
        write(outfileindex, '(a)')'set origin 0.05,0'
        write(outfileindex, '(a)')'set border lw 3'
        write(outfileindex, '(a)')'set pm3d'
        write(outfileindex, '(a)')'unset key'
        write(outfileindex, '(a)')'set view map'
        write(outfileindex, '(a)')'#set xtics font ",24"'
        write(outfileindex, '(a)')'#set ytics font ",24"'
        write(outfileindex, '(a)')'set xlabel "k1"'
        write(outfileindex, '(a)')'set ylabel "k2"'
        write(outfileindex, '(a)')'#set ylabel offset 1.5,0'
        write(outfileindex, '(a)')'set size ratio -1'
        write(outfileindex, '(a)')"set colorbox"
        write(outfileindex, '(a)')'set pm3d interpolate 2,2'
        write(outfileindex, '(2a)')"splot 'qpi_unfold.dat' u 4:5:($11) ",  &
           " w pm3d"
        close(outfileindex)
     endif
   ENDIF ! qpi_unfold
#if defined (MPI)
     call mpi_barrier(mpi_cmw, ierr)
#endif

   return
end subroutine unfolding_kplane


subroutine get_projection_weight_bulk_unfold(ndim, k_SBZ_direct, k_PBZ_direct, psi, weight, origincell, &
      NumberofSelectedOrbitals_groups, NumberofSelectedOrbitals, Selected_WannierOrbitals)
   !> Calculate the weights of given selected orbitals and mode for given wavefunction psi
   !> There are two modes. One is project on the selected orbitals. 
   !> The other one is projected on a special mode of another lattice which is usually a folded lattice.
   !> For the first mode, numk should be one.
   !> SBZ : Supercell Brillouin zone, usually, this is the Origin_cell.
   !> PBZ : primitive cell Brillouin zone, usually, this is the Folded_cell.
   !> Implemented by QSWU 2019
   use para, only : dp, projection_weight_mode, &
      cpuid, stdout, Nrpts, irvec, global_shift_SC_to_PC_cart, &
      Folded_cell, eps3, cell_type, int_array1D, Landaulevel_unfold_line_calc, twopi, zi, Matrix_Element_calc, SOC!@

   use me_calculate !@ use "!@" to mark the place modified by ycshen
   implicit none

   integer, intent(in) :: ndim
   !> a k vector associated with the wave function psi in the supercell Brillouin zone
   real(dp), intent(in) :: k_SBZ_direct(3)
   real(dp), intent(in) :: k_PBZ_direct(3)

   !> psi is a given vector for a k point and one band
   complex(dp), intent(in) :: psi(ndim)
   type(cell_type) :: origincell

   real(dp), intent(out) :: weight(NumberofSelectedOrbitals_groups)
   real(dp), allocatable :: weight_matix_element(:)

   integer, intent(in) :: NumberofSelectedOrbitals_groups
   integer, intent(in) :: NumberofSelectedOrbitals(NumberofSelectedOrbitals_groups)
   type(int_array1D) :: Selected_WannierOrbitals(NumberofSelectedOrbitals_groups)

   !> local variables
   integer :: io, io_SC, io_PC, projector_SC, projector_PC
   integer :: ia, ig, ir, ik, i, j, icount
   real(dp) :: posi_cart(3), posi_direct(3), posi_direct_unfold(3)
   real(dp) :: tau_i_tilde(3), tau_j_tilde(3), dij_tilde_cart(3), dij_tilde_direct(3)
   real(dp) :: k_cart(3), k_SBZ_direct_in_PBZ(3), k_t(3), k_t2(3), k_PBZ_direct_in_SBZ(3)
   real(dp) :: kdotr, k_abs         
   complex(dp) :: overlp, overlp_matrix_element
   character(10) :: atom_name_PC, atom_name_SC, projector_name_PC
   complex(dp) :: me !@ matrix element
   complex(dp), allocatable :: me_values(:)!@ matrix element values

   !> delta function
   real(dp), external :: delta, norm

   allocate(weight_matix_element(NumberofSelectedOrbitals_groups))
   weight_matix_element= 0d0


   !> k_PBZ_direct and k_SBZ_direct should be different by an reciprocal lattice vector of the Origin_cell (SBZ)
   call direct_cart_rec_unfold(k_PBZ_direct, k_cart)
   if (Landaulevel_unfold_line_calc) then
      call cart_direct_rec_magneticcell(k_cart, k_PBZ_direct_in_SBZ)
   else
      call cart_direct_rec(k_cart, k_PBZ_direct_in_SBZ)
   endif
   call periodic_diff(k_PBZ_direct_in_SBZ, k_SBZ_direct, k_t)

   !> check if k_t is the integer times of the reciprocal lattice vector of Origin_cell
   if (norm(k_t)>eps3) then
      weight= 0d0
      return
   endif

   !> use Folded_cell as a reference cell 
   if (Landaulevel_unfold_line_calc) then
      call direct_cart_rec_magneticcell(k_SBZ_direct, k_cart)
   else
      call direct_cart_rec(k_SBZ_direct, k_cart)
   endif
   call cart_direct_rec_unfold(k_cart, k_SBZ_direct_in_PBZ)

   !@ let k_cart be the transformation of k_PBZ, not the k_SPZ in the first brillouin zone
   call direct_cart_rec_unfold(k_PBZ_direct, k_cart)

   weight= 0d0

   !> k_t is in unit of the reciprocal lattice vector of the primitive unit cell (Folded_cell).
   k_t=k_PBZ_direct-k_SBZ_direct_in_PBZ

   allocate(me_values(Folded_cell%NumberofSpinOrbitals))

   if ( Matrix_Element_calc .eqv. .True. ) then
      !@ k_abs is the k_f(3) considering the photon energy
      if ( (k_cart_abs**2 - k_cart(1)**2 - k_cart(2)**2 ) .le. 0 ) then
         weight = 0
         return
      end if
      k_abs = sqrt(k_cart_abs**2 - k_cart(1)**2 - k_cart(2)**2)
      
      
   
     
    
         !@ calculate the matrix element for each spinorbital
         do io_PC=1, Folded_cell%NumberofSpinOrbitals
      
            atom_name_PC= adjustl(trim(Folded_cell%atom_name(Folded_cell%spinorbital_to_atom_index(io_PC))))
            projector_PC= Folded_cell%spinorbital_to_projector_index(io_PC)
   
            projector_name_PC= adjustl(trim(Folded_cell%proj_name(Folded_cell%spinorbital_to_projector_index(io_PC), Folded_cell%spinorbital_to_atom_index(io_PC))))
            me = 0d0 !! initialize the variable
            call get_matrix_element(atom_name_PC, projector_name_PC, k_cart, Folded_cell%wannier_centers_direct(:, io_PC), me)
            me_values(io_PC) = me
            
         enddo ! io_PC
   
      
   end if

   
   do ig=1, NumberofSelectedOrbitals_groups
      overlp_matrix_element=0d0
      do io_PC=1, Folded_cell%NumberofSpinOrbitals
         icount = 0
         overlp=0d0
         atom_name_PC= adjustl(trim(Folded_cell%atom_name(Folded_cell%spinorbital_to_atom_index(io_PC))))
         projector_PC= Folded_cell%spinorbital_to_projector_index(io_PC)
         tau_j_tilde= Folded_cell%wannier_centers_direct(:, io_PC)

         do io=1, NumberofSelectedOrbitals(ig)
            io_SC = Selected_WannierOrbitals(ig)%iarray(io)
            atom_name_SC= adjustl(trim(origincell%atom_name(origincell%spinorbital_to_atom_index(io_SC))))
            projector_SC= origincell%spinorbital_to_projector_index(io_SC)

            !> The atom name and the orbital should be the same between SC and PC
            if (atom_name_SC/=atom_name_PC .or. projector_SC/=projector_PC)cycle

            !> the atom position in the SuperCell. global_shift_SC_to_PC_cart is defined in readinput.f90
            posi_cart=origincell%wannier_centers_cart(:, io_SC)+ global_shift_SC_to_PC_cart
            call cart_direct_real_unfold(posi_cart, posi_direct_unfold)
            
            tau_i_tilde= posi_direct_unfold- floor(posi_direct_unfold)

            !> here we only take the lattice part
            kdotr=dot_product(posi_direct_unfold, k_t)

            call periodic_diff(tau_j_tilde, tau_i_tilde, dij_tilde_direct)
            call direct_cart_real_unfold(dij_tilde_direct, dij_tilde_cart)

            if (delta(0.1d0, norm(dij_tilde_cart))>0.01d0) then
               icount=icount+ 1
            endif

            !@ make sure the spin part can get match

            if (SOC > 0) then
               if ((io_PC-Folded_cell%NumberofSpinOrbitals/2d0 > 0) .neqv. (io_SC-origincell%NumberofSpinOrbitals/2d0 > 0) ) then
                  cycle
               endif
            end if

            
            !> brodening is 0.2 Bohr
            overlp= overlp+ delta(0.2d0, norm(dij_tilde_cart))*(cos(twopi*kdotr)-zi*sin(twopi*kdotr))*psi(io_SC)/delta(0.2d0, 0d0)
            
            overlp_matrix_element= overlp_matrix_element+ delta(0.2d0, norm(dij_tilde_cart))*(cos(twopi*kdotr)-zi*sin(twopi*kdotr))*psi(io_SC)/delta(0.2d0, 0d0)*me_values(io_PC)&
            *(cos((k_abs-k_cart(3))*posi_cart(3))+zi*sin((k_abs-k_cart(3))*posi_cart(3)))&
            *exp(-(posi_cart(3)/penetration_lambda_arpes))
            
            
         enddo ! io
         weight(ig)= weight(ig)+ abs(overlp)**2/NumberofSelectedOrbitals(ig)
      enddo ! io_PC
      weight_matix_element(ig)= abs(overlp_matrix_element)**2/NumberofSelectedOrbitals(ig)
   enddo ! ig
   weight= weight/origincell%CellVolume*Folded_cell%CellVolume
   weight_matix_element= weight_matix_element/origincell%CellVolume*Folded_cell%CellVolume

   if ( Matrix_Element_calc .eqv. .True. ) then
   weight = weight_matix_element
   end if 
   

   return
end subroutine get_projection_weight_bulk_unfold
