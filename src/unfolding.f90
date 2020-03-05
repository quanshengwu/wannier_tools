
subroutine unfolding_kpath
   !> Unfold the energy bands of the supercell to a specified unit cell.
   !> Calculate unfolded band at (k, \omega)
   !> 
   !> Implemented by QSWU 2019
   use para
  !use sparse
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
   !> eigenvector of the sparse matrix acoo. Dim=(Num_wann, neval)

   integer :: neval
   !number of ARPACK eigenvalues

   integer :: nvecs
   ! number of Arnoldi vectors

   logical :: ritzvec
   !> calculate eigenvector or not

   complex(dp) :: sigma=(0d0,0d0)
   !> shift-invert sigma

   integer :: n, i, ie, ik, ig, ierr, ieta,NumberofEta
   
   ! time measurement
   real(dp) :: time1, time2, time3

   real(dp), allocatable :: weight(:)
   real(dp), external :: delta
   NumberofEta = 9

   sigma=(1d0,0d0)*E_arc
   if (Is_Sparse_Hr) then

      neval=OmegaNum
      if (neval>Num_wann-2) then
         neval= Num_wann- 2
      endif
   
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
   allocate( spectrum_unfold(omeganum_unfold, NumberofEta, NumberofSelectedOrbitals_groups, nk3_band)) 
   allocate( spectrum_unfold_mpi(omeganum_unfold, NumberofEta, NumberofSelectedOrbitals_groups, nk3_band)) 
   spectrum_unfold= 0d0
   spectrum_unfold_mpi= 0d0

   NumberofEta=9
   allocate(eta_array(NumberofEta))
   eta_array=(/0.1d0, 0.2d0, 0.4d0, 0.8d0, 1.0d0, 2d0, 4d0, 8d0, 10d0/)
   eta_array= eta_array*Eta_Arc

   allocate(omega(omeganum_unfold))
   omega= 0d0
   do i= 1, omeganum_unfold
      omega(i)=omegamin+(i-1)*(omegamax-omegamin)/dble(omeganum_unfold)
   enddo

   !> first unfold the kpoints from the kpath of the supercell
   do ik= 1+cpuid, nk3_band, num_cpu
      if (cpuid==0.and.mod(ik, 40)==1) write(stdout, '(a, i8," /", i8)') 'BulkBand unfolding at :', ik, nk3_band
      k_PBZ_direct= k3points(:, ik)
      call direct_cart_rec_unfold(k_PBZ_direct, k_cart)
      call cart_direct_rec(k_cart, k_PBZ_direct_in_SBZ)

      k_SBZ_direct= k_PBZ_direct_in_SBZ- floor(k_PBZ_direct_in_SBZ)

      zeigv=0d0
      if (Is_Sparse_Hr) then
         stop " we don't support sparse hr.dat yet"
         ! sparse hr
         call ham_bulk_coo_sparsehr(k_SBZ_direct,acoo,icoo,jcoo)
         nnz= splen
 
         !> diagonalization by call zheev in lapack
         W= 0d0
         !> after arpack_sparse_coo_eigs, nnz will be updated.
         ritzvec= .true.
        !call arpack_sparse_coo_eigs(Num_wann,nnzmax,nnz,acoo,jcoo,icoo,neval,nvecs,W,sigma, zeigv, ritzvec)
         call now(time3)
      else
         ! dense hr
         call ham_bulk_atomicgauge(k_SBZ_direct, zeigv)
         
         ! diagonalize the Hamiltonian, zeigv is the Hamiltonian matrix before eigensystem_c calling.
         ! It will be replaced by the eigenvectors after eigensystem_c calling. 
         call eigensystem_c('V', 'U', Num_wann ,zeigv, W)
      endif
      do n= 1, neval
         psi= zeigv(:, n)
         call get_projection_weight_bulk_unfold(k_SBZ_direct, k_PBZ_direct, psi, weight)
         do ig=1, NumberofSelectedOrbitals_groups
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
      write(outfileindex, "('#column', i5, 3000i12)")(i, i=1, 2+NumberofSelectedOrbitals_groups*NumberofEta)
      do ik=1, nk3_band
         do ie=1, omeganum_unfold
            write(outfileindex, '(300f12.6)')k3len_folded(ik), omega(ie), &
               ((spectrum_unfold_mpi(ie, ieta, ig, ik), ieta=1, NumberofEta), ig=1, NumberofSelectedOrbitals_groups)
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
      write(outfileindex, '(a)')"set output 'spectrum_unfoldz_kpath.png'"
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
      write(outfileindex, '(a,f12.6)')'emin=', omegamin
      write(outfileindex, '(a,f12.6)')'emax=', omegamax
      write(outfileindex, '(a, f10.5, a)')'set xrange [0: ', maxval(k3len_folded), ']'
      if (index(Particle,'phonon')/=0) then
         write(outfileindex, '(a, f10.5, a)')'set yrange [0: emax ]'
         write(outfileindex, '(a)')'set ylabel "Frequency (THz)"'
      else
         write(outfileindex, '(a)')'set ylabel "Energy (eV)"'
         write(outfileindex, '(a)')'set yrange [ emin : emax ]'
      endif
      write(outfileindex, 202, advance="no") (k3line_name(i), k3line_folded_stop(i), i=1, nk3lines)
      write(outfileindex, 203)k3line_name(nk3lines+1), k3line_folded_stop(nk3lines+1)

      do i=1, nk3lines-1
         if (index(Particle,'phonon')/=0) then
            write(outfileindex, 204)k3line_folded_stop(i+1), '0.0', k3line_folded_stop(i+1), 'emax'
         else
            write(outfileindex, 204)k3line_folded_stop(i+1), 'emin', k3line_folded_stop(i+1), 'emax'
         endif
      enddo

      write(outfileindex, '(a)')"set colorbox"
      write(outfileindex, '(a)')'set pm3d interpolate 2,2'
      write(outfileindex, '(2a)')"splot 'spectrum_unfold_kpath.dat' u 1:2:(($11)) ",  &
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
   !> Calculate unfolded band at (k1, k2) at a given energy E_arc.
   !> Implemented by QSWU 2019
   use para
  !use sparse
   implicit none

   real(dp), allocatable :: eta_array(:)

   real(dp), allocatable :: spectrum_unfold(:, :, :), spectrum_unfold_mpi(:, :, :)
   !> unfolded spectrum, dimension (omeganum_unfold, NumberofEta, NumberofSelectedOrbitals_groups, nk3_band)

   real(dp) :: k_PBZ_direct(3), k_PBZ_direct_in_SBZ(3), k_cart(3), k_SBZ_direct(3)


   integer :: nnzmax, nnz, knv3
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

   real(dp), allocatable :: kxy(:,:), kxy_shape(:,:), kxy_plane(:,:)

   integer :: n, i, j, ie, ik, ig, ierr, ieta,NumberofEta
    
   
   ! time measurement
   real(dp) :: time1, time2, time3

   real(dp), allocatable :: weight(:)
   real(dp), external :: delta
   NumberofEta = 9

   knv3= Nk1*Nk2
   allocate( kxy(3, knv3))
   allocate( kxy_shape(3, knv3))
   allocate( kxy_plane(3, knv3))
   kxy=0d0
   kxy_shape=0d0
   kxy_plane=0d0
   
   ik =0
   do i= 1, Nk1
      do j= 1, Nk2
         ik =ik +1
         kxy(:, ik)= K3D_start+ K3D_vec1*(i-1)/dble(Nk1-1)+ K3D_vec2*(j-1)/dble(Nk2-1) &
            -(K3D_vec1+K3D_vec2)/2d0
         kxy_shape(:, ik)= kxy(1, ik)* Folded_cell%Kua+ kxy(2, ik)* Folded_cell%Kub+ kxy(3, ik)* Folded_cell%Kuc 
         call rotate_k3_to_kplane(kxy_shape(:, ik), kxy_plane(:, ik))
      enddo
   enddo


   sigma=(1d0,0d0)*E_arc
   if (Is_Sparse_Hr) then

      neval=OmegaNum
      if (neval>Num_wann-2) then
         neval= Num_wann- 2
      endif
   
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
   allocate( spectrum_unfold(NumberofEta, NumberofSelectedOrbitals_groups, knv3)) 
   allocate( spectrum_unfold_mpi(NumberofEta, NumberofSelectedOrbitals_groups, knv3)) 
   spectrum_unfold= 0d0
   spectrum_unfold_mpi= 0d0

   NumberofEta=9
   allocate(eta_array(NumberofEta))
   eta_array=(/0.1d0, 0.2d0, 0.4d0, 0.8d0, 1.0d0, 2d0, 4d0, 8d0, 10d0/)
   eta_array= eta_array*Eta_Arc

   !> first unfold the kpoints from the kpath of the supercell
   do ik= 1+cpuid, knv3, num_cpu
      if (cpuid==0.and.mod(ik, 40)==1) write(stdout, '(a, i10," /", i10)') 'BulkBand unfolding at :', ik, knv3
      k_PBZ_direct= kxy(:, ik)
      call direct_cart_rec_unfold(k_PBZ_direct, k_cart)
      call cart_direct_rec(k_cart, k_PBZ_direct_in_SBZ)

      k_SBZ_direct= k_PBZ_direct_in_SBZ- floor(k_PBZ_direct_in_SBZ)

      if (Is_Sparse_Hr) then
         stop " we don't support sparse hr.dat yet"
         ! sparse hr
         call ham_bulk_coo_sparsehr(k_SBZ_direct,acoo,icoo,jcoo)
         nnz= splen
 
         !> diagonalization by call zheev in lapack
         W= 0d0
         !> after arpack_sparse_coo_eigs, nnz will be updated.
         ritzvec= .true.
        !call arpack_sparse_coo_eigs(Num_wann,nnzmax,nnz,acoo,jcoo,icoo,neval,nvecs,W,sigma, zeigv, ritzvec)
         call now(time3)
      else
         ! dense hr
         call ham_bulk_atomicgauge(k_SBZ_direct, zeigv)
         
         ! diagonalize the Hamiltonian, zeigv is the Hamiltonian matrix before eigensystem_c calling.
         ! It will be replaced by the eigenvectors after eigensystem_c calling. 
         call eigensystem_c('V', 'U', Num_wann ,zeigv, W)
      endif
      do n= 1, neval
         psi= zeigv(:, n)
         call get_projection_weight_bulk_unfold(k_SBZ_direct, k_PBZ_direct, psi, weight)
         do ig=1, NumberofSelectedOrbitals_groups
            do ieta= 1, NumberofEta
               spectrum_unfold(ieta, ig, ik)= spectrum_unfold(ieta, ig, ik) + &
                  weight(ig)*delta(eta_array(ieta), W(n)-E_arc)
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
      write(outfileindex, "('#', a8, 5a12, 3X, '| A(k,E)', a6, 100(8X,'group ', i2))")&
         'kx', 'ky', 'kz', 'kp1', 'kp2', 'kp3', 'total',&
         (i, i=1, NumberofSelectedOrbitals_groups)
      write(outfileindex, "('#column', i5, 3000i12)")(i, i=1, 7+NumberofSelectedOrbitals_groups*NumberofEta)
      do ik=1, knv3
         write(outfileindex, '(3000f12.5)')kxy_shape(:, ik), kxy_plane(:, ik), &
              ((spectrum_unfold_mpi(ieta, ig, ik), ieta=1, NumberofEta), ig=1, NumberofSelectedOrbitals_groups)
         if (mod(ik, nk2)==0) write(outfileindex, *)' '
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
      write(outfileindex, '(a)')"set output 'spectrum_unfoldz_kplane.png'"
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

#if defined (MPI)
     call mpi_barrier(mpi_cmw, ierr)
#endif

   return
end subroutine unfolding_kplane


subroutine get_projection_weight_bulk_unfold(k_SBZ_direct, k_PBZ_direct, psi, weight)
   !> Calculate the weights of given selected orbitals and mode for given wavefunction psi
   !> There are two modes. One is project on the selected orbitals. 
   !> The other one is projected on a special mode of another lattice which is usually a folded lattice.
   !> For the first mode, numk should be one.
   !> SBZ : Supercell Brillouin zone, usually, this is the Origin_cell.
   !> PBZ : primitive cell Brillouin zone, usually, this is the Folded_cell.
   !> Implemented by QSWU 2019
   use para, only : dp, num_wann, NumberofSelectedOrbitals_groups, projection_weight_mode, &
      NumberofSelectedOrbitals, Selected_WannierOrbitals, cpuid, stdout, Nrpts, irvec, &
      Origin_cell, Folded_cell, pi2zi, eps3
   implicit none

   !> a k vector associated with the wave function psi in the supercell Brillouin zone
   real(dp), intent(in) :: k_SBZ_direct(3)
   real(dp), intent(in) :: k_PBZ_direct(3)

   !> psi is a given vector for a k point and one band
   complex(dp), intent(in) :: psi(num_wann)

   real(dp), intent(out) :: weight(NumberofSelectedOrbitals_groups)

   !> local variables
   integer :: io, io_SC, io_PC, projector_SC, projector_PC
   integer :: ia, ig, ir, ik, i, j, icount
   real(dp) :: posi_cart(3), posi_direct(3), posi_direct_unfold(3)
   real(dp) :: tau_i_tilde(3), tau_j_tilde(3), dij_tilde_cart(3), dij_tilde_direct(3)
   real(dp) :: k_cart(3), k_SBZ_direct_in_PBZ(3), k_t(3), k_t2(3), k_PBZ_direct_in_SBZ(3)
   real(dp) :: kdotr        
   complex(dp) :: overlp
   character(10) :: atom_name_PC, atom_name_SC

   !> delta function
   real(dp), external :: delta, norm


   !> k_PBZ_direct and k_SBZ_direct should be different by an reciprocal lattice vector of the Origin_cell (SBZ)
   call direct_cart_rec_unfold(k_PBZ_direct, k_cart)
   call cart_direct_rec(k_cart, k_PBZ_direct_in_SBZ)
   call periodic_diff(k_PBZ_direct_in_SBZ, k_SBZ_direct, k_t)

   !> check if k_t is the integer times of the reciprocal lattice vector of Origin_cell
   if (norm(k_t)>eps3) then
      weight= 0d0
      return
   endif

   !> use Folded_cell as a reference cell 
   call direct_cart_rec(k_SBZ_direct, k_cart)
   call cart_direct_rec_unfold(k_cart, k_SBZ_direct_in_PBZ)

   weight= 0d0

   !> k_t is in unit of the reciprocal lattice vector of the primitive unit cell (Folded_cell).
   k_t=k_PBZ_direct-k_SBZ_direct_in_PBZ


   do ig=1, NumberofSelectedOrbitals_groups
      do io_PC=1, Folded_cell%NumberofSpinOrbitals
         icount = 0
         overlp=0d0
         atom_name_PC= adjustl(trim(Folded_cell%atom_name(Folded_cell%spinorbital_to_atom_index(io_PC))))
         projector_PC= Folded_cell%spinorbital_to_projector_index(io_PC)
         tau_j_tilde= Folded_cell%wannier_centers_direct(:, io_PC)

         do io=1, NumberofSelectedOrbitals(ig)
            io_SC = Selected_WannierOrbitals(ig)%iarray(io)
            atom_name_SC= adjustl(trim(Origin_cell%atom_name(Origin_cell%spinorbital_to_atom_index(io_SC))))
            projector_SC= Origin_cell%spinorbital_to_projector_index(io_SC)

            !> The atom name and the orbital should be the same between SC and PC
            if (atom_name_SC/=atom_name_PC .or. projector_SC/=projector_PC)cycle

            !> the atom position in the SuperCell.
            posi_direct=Origin_cell%wannier_centers_direct(:, io_SC)
            call direct_cart_real(posi_direct, posi_cart)
            call cart_direct_real_unfold(posi_cart, posi_direct_unfold)
            
            tau_i_tilde= posi_direct_unfold- floor(posi_direct_unfold)

            !> here we only take the lattice part
            kdotr=dot_product(posi_direct_unfold, k_t)

            call periodic_diff(tau_j_tilde, tau_i_tilde, dij_tilde_direct)
            call direct_cart_real_unfold(dij_tilde_direct, dij_tilde_cart)

            if (delta(0.1d0, norm(dij_tilde_cart))>0.01d0) then
               icount=icount+ 1
            endif

            !> brodening is 0.1 Angstrom
            overlp= overlp+ delta(0.1d0, norm(dij_tilde_cart))*exp(-pi2zi*(kdotr))*psi(io_SC)/delta(0.1d0, 0d0)

         enddo ! io
         weight(ig)= weight(ig)+ abs(overlp)**2
      enddo ! io_PC
   enddo ! ig
   weight= weight/Origin_cell%CellVolume*Folded_cell%CellVolume

   return
end subroutine get_projection_weight_bulk_unfold
