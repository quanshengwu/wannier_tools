
subroutine dos_sparse
!> calculate density of state for 3D bulk system
!
!> DOS(\omega)= \sum_k \delta(\omega- E(k))
   use sparse
   use wmpi
   use para
   implicit none

   !> the integration k space
   real(dp) :: emin, emax

   integer :: ik, ie, ib, ikx, iky, ikz
   integer :: knv3, NE, ierr, ieta
   integer :: NumberofEta

   real(dp) :: x, dk3, eta0

   real(dp) :: k(3)
   real(dp) :: time_start, time_end

   real(dp), allocatable :: eigval(:)
   real(dp), allocatable :: W(:)
   real(dp), allocatable :: omega(:)
   real(dp), allocatable :: dos(:, :), dos_mpi(:, :)
   real(dp), external :: delta
   real(dp), allocatable :: eta_array(:)

   logical :: ritzvec


   complex(dp), allocatable :: acoo(:),zeigv(:,:)
   integer,allocatable :: icoo(:),jcoo(:)
   integer :: neval,Ndimq,nnz,nnzmax,nvecs
   complex(dp) :: sigma=(0d0,0d0)

   ndimq=Num_wann
   ritzvec= .false.

   if (NumSelectedEigenVals==0) NumSelectedEigenVals=Num_wann

   neval=NumSelectedEigenVals
   if (neval>Ndimq-2) neval= Ndimq- 2

   !> ncv
   nvecs=int(1.5*neval)

   if (nvecs<50) nvecs= 50
   if (nvecs>Ndimq) nvecs= Ndimq
   !> delta function

   knv3= Nk1*Nk2*Nk3

   if (OmegaNum<2) OmegaNum=2
   NE= OmegaNum
   sigma=(1d0,0d0)*iso_energy
   nnzmax=splen+ Ndimq
   nnz=splen

   NumberofEta=9

   allocate(W(Num_wann))
   allocate(eigval(neval))
   allocate(eta_array(NumberofEta))
   allocate(dos(NE, NumberofEta))
   allocate(dos_mpi(NE, NumberofEta))
   allocate(omega(NE))
   allocate( acoo(nnzmax))
   allocate( jcoo(nnzmax))
   allocate( icoo(nnzmax))
   allocate( zeigv(ndimq,nvecs))
   dos=0d0
   dos_mpi=0d0
   eigval= 0d0


   emin= OmegaMin
   emax= OmegaMax

   eta_array=(/0.1d0, 0.2d0, 0.4d0, 0.8d0, 1.0d0, 2d0, 4d0, 8d0, 10d0/)
   eta_array= eta_array*Fermi_broadening

   !eta= (emax- emin)/ dble(NE)*3d0

   !> energy
   do ie=1, NE
      omega(ie)= emin+ (emax-emin)* (ie-1d0)/dble(NE-1)
   enddo ! ie

   !dk3= kCubeVolume/dble(knv3)
   dk3= 1d0/dble(knv3)

   !> get eigenvalue
   time_start= 0d0
   time_end= 0d0
   do ik=1+cpuid, knv3, num_cpu

      if (cpuid.eq.0.and. mod(ik/num_cpu, 100).eq.0) &
         write(stdout, '(a, i18, "/", i18, a, f10.3, "s")') 'ik/knv3', &
         ik, knv3, ' time left', (knv3-ik)*(time_end-time_start)/num_cpu


      call now(time_start)
      ikx= (ik-1)/(nk2*nk3)+1
      iky= ((ik-1-(ikx-1)*Nk2*Nk3)/nk3)+1
      ikz= (ik-(iky-1)*Nk3- (ikx-1)*Nk2*Nk3)
      k= K3D_start_cube+ K3D_vec1_cube*(ikx-1)/dble(nk1)  &
         + K3D_vec2_cube*(iky-1)/dble(nk2)  &
         + K3D_vec3_cube*(ikz-1)/dble(nk3)

      call ham_bulk_coo_sparsehr(k,acoo,icoo,jcoo)
      W= 0d0
      call arpack_sparse_coo_eigs(Ndimq,nnzmax,nnz,acoo,jcoo,icoo,neval,nvecs,W,sigma, zeigv, ritzvec)

      eigval(:)= W(1:neval)

      !> get density of state
      !> dos(e)= \sum_nk \delta(e-e_nk)
      do ie= 1, NE
         do ib= 1, neval
            x= omega(ie)- eigval(ib)
            do ieta= 1, NumberofEta
               eta0= eta_array(ieta)
               dos_mpi(ie, ieta) = dos_mpi(ie, ieta)+ delta(eta0, x)
            enddo
         enddo ! ib
      enddo ! ie
      call now(time_end)

   enddo  ! ik

#if defined (MPI)
   call mpi_allreduce(dos_mpi,dos,size(dos),&
      mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
   dos= dos_mpi
#endif
   dos= dos*dk3

   outfileindex= outfileindex+ 1
   if (cpuid.eq.0) then
      open(unit=outfileindex, file='dos.dat')
      write(outfileindex, *)'# Density of state of bulk system'
      write(outfileindex, '(a16, a)')'# E(eV)', 'DOS(E) (states/unit cell/eV)'
      write(outfileindex, '("#", a, f6.2, 300f16.2)')'Broadening \eta (meV): ', Eta_array(:)*1000d0/eV2Hartree
      do ie=1, NE
         write(outfileindex, '(90f16.6)')omega(ie)/eV2Hartree, dos(ie, :)*eV2Hartree
      enddo ! ie
      close(outfileindex)
   endif

   outfileindex= outfileindex+ 1
   !> write script for gnuplot
   if (cpuid==0) then
      open(unit=outfileindex, file='dos.gnu')
      write(outfileindex, '(a)')"set encoding iso_8859_1"
      write(outfileindex, '(a)')'set terminal pdf enhanced color font ",16" size 5,4'
      write(outfileindex, '(a)')"set output 'dos.pdf'"
      write(outfileindex, '(a)')'set border lw 2'
      write(outfileindex, '(a)')'set autoscale fix'
      write(outfileindex, '(a, f16.6,a)')'set yrange [0:', maxval(dos)*eV2Hartree+0.5, '1]'
      write(outfileindex, '(a)')'set key samplen 0.8 spacing 1 font ",12"'
      write(outfileindex, '(a)')'set xlabel "Energy (eV)"'
      write(outfileindex, '(a)')'set ylabel "DOS (states/eV/unit cell)"'
      write(outfileindex, '(a, f6.1, a)')"plot 'dos.dat' u 1:2 w l lw 2 title '",&
         Eta_array(1)*1000/eV2Hartree, "meV', \"
      do ieta= 2, NumberofEta-1
         write(outfileindex, 201)" '' u 1:", ieta, " w l lw 2 title '", &
            Eta_array(ieta)*1000/eV2Hartree, "meV', \"
      enddo
      write(outfileindex, '(a, f6.1, a)')" '' u 1:10 w l lw 2 title '",&
         Eta_array(NumberofEta)*1000/eV2Hartree, "meV'"
      close(outfileindex)
   endif
201 format(a, i3, a, f6.1, a)

#if defined (MPI)
   call mpi_barrier(mpi_cmw, ierr)
#endif
   deallocate(W)
!~    deallocate(Hk)
   deallocate(eigval)
   deallocate(dos)
   deallocate(dos_mpi)
   deallocate(omega)

   return
end subroutine dos_sparse

subroutine charge_density_sparse
!> calculate charge density
!
!> rho(r)= \sum_k |\psi_k(r)|^2
   use sparse
   use wmpi
   use para
   implicit none

   integer :: ik, ie, ib, ikx, iky, ikz, i, j, ia, n
   integer :: knv3, NE, ierr, ieta
   integer :: NumberofEta

   real(dp) :: x, dk3, eta0

   real(dp) :: k(3)
   real(dp) :: time_start, time_end

   real(dp), allocatable :: eigval(:)
   real(dp), allocatable :: W(:)
   real(dp), allocatable :: omega(:)
   real(dp), allocatable :: chargedensity(:), chargedensity_mpi(:)
   real(dp), external :: delta

   logical :: ritzvec


   complex(dp), allocatable :: acoo(:),zeigv(:,:)
   integer,allocatable :: icoo(:),jcoo(:)
   integer :: neval,Ndimq,nnz,nnzmax,nvecs
   complex(dp) :: sigma=(0d0,0d0)

   ndimq=Num_wann
   ritzvec= .true.

   if (NumSelectedEigenVals==0) then
      if (OmegaNum==0) then
         NumSelectedEigenVals=Num_wann
      else
         NumSelectedEigenVals=OmegaNum
      endif
   endif

   neval=NumSelectedEigenVals
   if (neval>Ndimq-2) neval= Ndimq- 2

   !> ncv
   nvecs=int(1.5*neval)

   if (nvecs<50) nvecs= 50
   if (nvecs>Ndimq) nvecs= Ndimq
   !> delta function

   knv3= Nk1*Nk2*Nk3

   if (OmegaNum<2) OmegaNum=2
   NE= OmegaNum
   sigma=(1d0,0d0)*iso_energy
   nnzmax=splen+ Ndimq
   nnz=splen


   allocate(W(Num_wann))
   allocate(eigval(neval))
   allocate(chargedensity(Origin_cell%Num_atoms))
   allocate(chargedensity_mpi(Origin_cell%Num_atoms))
   allocate(omega(NE))
   allocate( acoo(nnzmax))
   allocate( jcoo(nnzmax))
   allocate( icoo(nnzmax))
   allocate( zeigv(ndimq,nvecs))
   chargedensity=0d0
   chargedensity_mpi=0d0
   eigval= 0d0


   !> get eigenvalue
   time_start= 0d0
   time_end= 0d0
   do ik=1+cpuid, knv3, num_cpu

      if (cpuid.eq.0.and. mod((ik-1)/num_cpu, 100).eq.0) &
         write(stdout, '(a, i18, "/", i18, a, f10.3, "s")') 'ik/knv3', &
         ik, knv3, ' time left', (knv3-ik)*(time_end-time_start)/num_cpu


      call now(time_start)
      ikx= (ik-1)/(nk2*nk3)+1
      iky= ((ik-1-(ikx-1)*Nk2*Nk3)/nk3)+1
      ikz= (ik-(iky-1)*Nk3- (ikx-1)*Nk2*Nk3)
      k= K3D_start_cube+ K3D_vec1_cube*(ikx-1)/dble(nk1)  &
         + K3D_vec2_cube*(iky-1)/dble(nk2)  &
         + K3D_vec3_cube*(ikz-1)/dble(nk3)

      call ham_bulk_coo_sparsehr(k,acoo,icoo,jcoo)
      W= 0d0
      call arpack_sparse_coo_eigs(Ndimq,nnzmax,nnz,acoo,jcoo,icoo,neval,nvecs,W,sigma, zeigv, ritzvec)

      !> get charge density
      do i=1, neval
         if (W(i)>OmegaMin .and.W(i)<OmegaMax) then
            n=0
            do ia= 1, Origin_cell%Num_atoms
               do j= 1, Origin_cell%nprojs(ia)
                  n=n+1
                  chargedensity_mpi(ia)= chargedensity_mpi(ia)+ abs(zeigv(n, i))**2
                  if (SOC>0) then
                     chargedensity_mpi(ia)= chargedensity_mpi(ia)+ abs(zeigv(n+Num_wann/2, i))**2
                  endif
               enddo
            enddo
         endif
      enddo ! i
      call now(time_end)
   enddo  ! ik

#if defined (MPI)
   call mpi_allreduce(chargedensity_mpi,chargedensity,size(chargedensity),&
      mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
   chargedensity= chargedensity_mpi
#endif
   chargedensity= chargedensity/knv3

   outfileindex= outfileindex+ 1
   if (cpuid.eq.0) then
      open(unit=outfileindex, file='chargedensity.dat')
      write(outfileindex, *)'# Density of state of bulk system'
      write(outfileindex, '(a16, a)')'# E(eV)', 'chargedensity(E) (states/unit cell/eV)'
      do ia=1, Origin_cell%Num_atoms
         write(outfileindex, '(90f16.6)')Origin_cell%Atom_position_cart(:, ia), chargedensity(ia)
      enddo ! ia
      close(outfileindex)
   endif

   outfileindex= outfileindex+ 1
   !> write script for gnuplot
   if (cpuid==0) then
      open(unit=outfileindex, file='chargedensity.gnu')
      write(outfileindex, '(a)')"set encoding iso_8859_1"
      write(outfileindex, '(a)')'set terminal  postscript enhanced color font ",24" '
      write(outfileindex, '(a)')"set output 'chargedensity.eps'"
      write(outfileindex, '(a)')'set border lw 2'
      write(outfileindex, '(a)')'set autoscale fix'
      write(outfileindex, '(a)')'set yrange [0:1]'
      write(outfileindex, '(a)')'unset key'
      write(outfileindex, '(a)')'set xlabel "Energy (eV)"'
      write(outfileindex, '(a)')'set ylabel "chargedensity (states/eV/unit cell)"'
      write(outfileindex, '(a)')"splot 'chargedensity.dat' u 1:2:3 w lp pt 6 ps 0.2 lw 1.0 lc rgb 'black'  "
      close(outfileindex)
   endif

#if defined (MPI)
   call mpi_barrier(mpi_cmw, ierr)
#endif
   deallocate(W)
   deallocate(eigval)
   deallocate(chargedensity)
   deallocate(chargedensity_mpi)

   return
end subroutine charge_density_sparse


subroutine dos_sub
!> calculate density of state for 3D bulk system
!
!> DOS(\omega)= \sum_k \delta(\omega- E(k))

   use wmpi
   use para
   implicit none

   !> the integration k space
   real(dp) :: emin, emax

   integer :: ik,ie,ib,ikx,iky,ikz,ieta
   integer :: knv3,NE,ierr
   integer :: NumberofEta

   !> integration for band
   integer :: iband_low,iband_high,iband_tot

   real(dp) :: x, dk3, eta0

   real(dp) :: k(3)
   real(dp) :: time_start, time_end

   real(dp), allocatable :: eigval(:)
   real(dp), allocatable :: W(:)
   real(dp), allocatable :: omega(:)
   real(dp), allocatable :: dos(:, :)
   real(dp), allocatable :: dos_mpi(:, :)
   real(dp), allocatable :: eta_array(:)
   complex(dp), allocatable :: Hk(:, :)

   !> delta function
   real(dp), external :: delta

   knv3= Nk1*Nk2*Nk3

   if (OmegaNum<2) OmegaNum=2
   NE= OmegaNum

   iband_low= Numoccupied- 10000
   iband_high= Numoccupied+ 10000

   if (iband_low <1) iband_low = 1
   if (iband_high >Num_wann) iband_high = Num_wann

   iband_tot= iband_high- iband_low+ 1

   NumberofEta = 9

   allocate(W(Num_wann))
   allocate(Hk(Num_wann, Num_wann))
   allocate(eigval(iband_tot))
   allocate(eta_array(NumberofEta))
   allocate(dos(NE, NumberofEta))
   allocate(dos_mpi(NE, NumberofEta))
   allocate(omega(NE))
   dos=0d0
   dos_mpi=0d0
   eigval= 0d0


   emin= OmegaMin
   emax= OmegaMax
   eta_array=(/0.1d0, 0.2d0, 0.4d0, 0.8d0, 1.0d0, 2d0, 4d0, 8d0, 10d0/)
   eta_array= eta_array*Fermi_broadening


   !> energy
   do ie=1, NE
      omega(ie)= emin+ (emax-emin)* (ie-1d0)/dble(NE-1)
   enddo ! ie

   !dk3= kCubeVolume/dble(knv3)
   dk3= 1d0/dble(knv3)

   !> get eigenvalue
   time_start= 0d0
   time_end= 0d0
   do ik=1+cpuid, knv3, num_cpu

      if (cpuid.eq.0.and. mod(ik/num_cpu, 100).eq.0) &
         write(stdout, '(a, i18, "/", i18, a, f10.3, "s")') 'ik/knv3', &
         ik, knv3, ' time left', (knv3-ik)*(time_end-time_start)/num_cpu

      call now(time_start)
      ikx= (ik-1)/(nk2*nk3)+1
      iky= ((ik-1-(ikx-1)*Nk2*Nk3)/nk3)+1
      ikz= (ik-(iky-1)*Nk3- (ikx-1)*Nk2*Nk3)
      k= K3D_start_cube+ K3D_vec1_cube*(ikx-1)/dble(nk1)  &
         + K3D_vec2_cube*(iky-1)/dble(nk2)  &
         + K3D_vec3_cube*(ikz-1)/dble(nk3)

      if (index(KPorTB, 'KP')/=0)then
         call ham_bulk_kp_abcb_graphene(k, Hk)
      else
         call ham_bulk_atomicgauge(k, Hk)
      endif

      W= 0d0
      call eigensystem_c( 'N', 'U', Num_wann ,Hk, W)
      eigval(:)= W(iband_low:iband_high)

      !> get density of state
      do ie= 1, NE
         do ib= 1, iband_tot
            x= omega(ie)- eigval(ib)
            do ieta= 1, NumberofEta
               eta0= eta_array(ieta)
               dos_mpi(ie, ieta) = dos_mpi(ie, ieta)+ delta(eta0, x)
            enddo
         enddo ! ib
      enddo ! ie
      call now(time_end)


      call now(time_end)

   enddo  ! ik

#if defined (MPI)
   call mpi_allreduce(dos_mpi,dos,size(dos),&
      mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
   dos= dos_mpi
#endif
   dos= dos*dk3

   !> include the spin degeneracy if there is no SOC in the tight binding Hamiltonian.
   if (SOC<=0) dos=dos*2d0

   outfileindex= outfileindex+ 1
   if (cpuid.eq.0) then
      open(unit=outfileindex, file='dos.dat')
      write(outfileindex, *)'# Density of state of bulk system'
      write(outfileindex, '(2a16)')'# E(eV)', 'DOS(E) (1/eV)'
      write(outfileindex, '("#", a, f6.2, 300f16.2)')'Broadening \eta (meV): ', Eta_array(:)*1000d0/eV2Hartree
      do ie=1, NE
         write(outfileindex, '(90f16.6)')omega(ie)/eV2Hartree, dos(ie, :)*eV2Hartree
      enddo ! ie
      close(outfileindex)
   endif

   outfileindex= outfileindex+ 1
   !> write script for gnuplot
   if (cpuid==0) then
      open(unit=outfileindex, file='dos.gnu')
      write(outfileindex, '(a)')"set encoding iso_8859_1"
      write(outfileindex, '(a)')'set terminal pdf enhanced color font ",16" size 5,4 '
      write(outfileindex, '(a)')"set output 'dos.pdf'"
      write(outfileindex, '(a)')'set border lw 2'
      write(outfileindex, '(a)')'set autoscale fix'
      write(outfileindex, '(a, f16.6,a)')'set yrange [0:', maxval(dos)*eV2Hartree+0.5, '1]'
      write(outfileindex, '(a)')'set key samplen 0.8 spacing 1 font ",12"'
      write(outfileindex, '(a)')'set xlabel "Energy (eV)"'
      write(outfileindex, '(a)')'set key title  "Broadening"'
      write(outfileindex, '(a)')'set title "DOS with different broadenings"'
      write(outfileindex, '(a)')'set ylabel "DOS (states/eV/unit cell)"'
      write(outfileindex, '(a, f6.1, a)')"plot 'dos.dat' u 1:2 w l lw 2 title '",&
         Eta_array(1)*1000/eV2Hartree, "meV', \"
      do ieta= 2, NumberofEta-1
         write(outfileindex, 202)" '' u 1:", ieta, " w l lw 2 title '", &
            Eta_array(ieta)*1000/eV2Hartree, "meV', \"
      enddo
      write(outfileindex, '(a, f6.1, a)')" '' u 1:10 w l lw 2 title '",&
         Eta_array(NumberofEta)*1000/eV2Hartree, "meV'"
      close(outfileindex)
   endif
202 format(a, i3, a, f6.1, a)

#if defined (MPI)
   call mpi_barrier(mpi_cmw, ierr)
#endif
   deallocate(W)
   deallocate(Hk)
   deallocate(eigval)
   deallocate(dos)
   deallocate(dos_mpi)
   deallocate(omega)

   return
end subroutine dos_sub



subroutine dos_slab
!> calculate density of state for slab system
!
!> DOS(\omega)= \sum_k \delta(\omega- E(k))

   use wmpi
   use para
   implicit none

   !> the integration k space
   real(dp) :: emin, emax

   integer :: ik,ie,ib,ikx,iky,ikz,ieta
   integer :: knv2_slab,NE,ierr, ndim_slab
   integer :: NumberofEta

   !> integration for band
   integer :: iband_low,iband_high,iband_tot

   real(dp) :: x, dk2, eta0

   real(dp) :: k(2)
   real(dp) :: time_start, time_end

   real(dp), allocatable :: eigval(:)
   real(dp), allocatable :: W(:)
   real(dp), allocatable :: omega(:)
   real(dp), allocatable :: dos(:, :)
   real(dp), allocatable :: dos_mpi(:, :)
   real(dp), allocatable :: eta_array(:)
   complex(dp), allocatable :: Hk(:, :)

   !> delta function
   real(dp), external :: delta

   knv2_slab= Nk1*Nk2
   ndim_slab= Num_wann*Nslab

   if (OmegaNum<2) OmegaNum=2
   NE= OmegaNum

   iband_low= Numoccupied- 10000
   iband_high= Numoccupied+ 10000

   if (iband_low <1) iband_low = 1
   if (iband_high >ndim_slab) iband_high = ndim_slab

   iband_tot= iband_high- iband_low+ 1

   NumberofEta = 9

   allocate(W(ndim_slab))
   allocate(Hk(ndim_slab, ndim_slab))
   allocate(eigval(iband_tot))
   allocate(eta_array(NumberofEta))
   allocate(dos(NE, NumberofEta))
   allocate(dos_mpi(NE, NumberofEta))
   allocate(omega(NE))
   dos=0d0
   dos_mpi=0d0
   eigval= 0d0


   emin= OmegaMin
   emax= OmegaMax
   eta_array=(/0.1d0, 0.2d0, 0.4d0, 0.8d0, 1.0d0, 2d0, 4d0, 8d0, 10d0/)
   eta_array= eta_array*Fermi_broadening


   !> energy
   do ie=1, NE
      omega(ie)= emin+ (emax-emin)* (ie-1d0)/dble(NE-1)
   enddo ! ie

   !dk2= kCubeVolume/dble(knv2_slab)
   dk2= 1d0/dble(knv2_slab)

   !> get eigenvalue
   time_start= 0d0
   time_end= 0d0
   do ik=1+cpuid, knv2_slab, num_cpu

      if (cpuid.eq.0.and. mod(ik/num_cpu, 100).eq.0) &
         write(stdout, '(a, i18, "/", i18, a, f10.3, "s")') 'ik/knv2_slab', &
         ik, knv2_slab, ' time left', (knv2_slab-ik)*(time_end-time_start)/num_cpu

      call now(time_start)
      ikx= (ik-1)/(Nk2)+1
      iky= ((ik-1-(ikx-1)*Nk2))+1
      k= K2D_start+ K2D_vec1*(ikx-1)/dble(Nk1)  &
         + K2D_vec2*(iky-1)/dble(Nk2) 

      !> get the Hamlitonian
      call ham_slab(k, Hk)

      W= 0d0
      call eigensystem_c('N', 'U', ndim_slab ,Hk, W)
      eigval(:)= W(iband_low:iband_high)

      !> get density of state
      do ie= 1, NE
         do ib= 1, iband_tot
            x= omega(ie)- eigval(ib)
            do ieta= 1, NumberofEta
               eta0= eta_array(ieta)
               dos_mpi(ie, ieta) = dos_mpi(ie, ieta)+ delta(eta0, x)
            enddo
         enddo ! ib
      enddo ! ie
      call now(time_end)

      call now(time_end)

   enddo  ! ik

#if defined (MPI)
   call mpi_allreduce(dos_mpi,dos,size(dos),&
      mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
   dos= dos_mpi
#endif
   dos= dos*dk2

   !> include the spin degeneracy if there is no SOC in the tight binding Hamiltonian.
   if (SOC<=0) dos=dos*2d0

   outfileindex= outfileindex+ 1
   if (cpuid.eq.0) then
      open(unit=outfileindex, file='dos_slab.dat')
      write(outfileindex, "(a, i10)")'# Density of state of slab system, Nslab= ', Nslab
      write(outfileindex, '(2a16)')'# E(eV)', 'DOS(E) (1/eV)'
      write(outfileindex, '("#", a, f6.2, 300f16.2)')'Broadening \eta (meV): ', Eta_array(:)*1000d0/eV2Hartree
      do ie=1, NE
         write(outfileindex, '(90f16.6)')omega(ie)/eV2Hartree, dos(ie, :)*eV2Hartree
      enddo ! ie
      close(outfileindex)
   endif

   outfileindex= outfileindex+ 1
   !> write script for gnuplot
   if (cpuid==0) then
      open(unit=outfileindex, file='dos_slab.gnu')
      write(outfileindex, '(a)')"set encoding iso_8859_1"
      write(outfileindex, '(a)')'set terminal pdf enhanced color font ",16" size 5,4 '
      write(outfileindex, '(a)')"set output 'dos_slab.pdf'"
      write(outfileindex, '(a)')'set border lw 2'
      write(outfileindex, '(a)')'set autoscale fix'
      write(outfileindex, '(a, f16.6,a)')'set yrange [0:', maxval(dos)*eV2Hartree+0.5, '1]'
      write(outfileindex, '(a)')'set key samplen 0.8 spacing 1 font ",12"'
      write(outfileindex, '(a)')'set xlabel "Energy (eV)"'
      write(outfileindex, '(a)')'set key title  "Broadening"'
      write(outfileindex, '(a)')'set title "DOS with different broadenings"'
      write(outfileindex, '(a)')'set ylabel "DOS (states/eV/unit cell)"'
      write(outfileindex, '(a, f6.1, a)')"plot 'dos_slab.dat' u 1:2 w l lw 2 title '",&
         Eta_array(1)*1000/eV2Hartree, "meV', \"
      do ieta= 2, NumberofEta-1
         write(outfileindex, 202)" '' u 1:", ieta, " w l lw 2 title '", &
            Eta_array(ieta)*1000/eV2Hartree, "meV', \"
      enddo
      write(outfileindex, '(a, f6.1, a)')" '' u 1:10 w l lw 2 title '",&
         Eta_array(NumberofEta)*1000/eV2Hartree, "meV'"
      close(outfileindex)
   endif
202 format(a, i3, a, f6.1, a)

#if defined (MPI)
   call mpi_barrier(mpi_cmw, ierr)
#endif
   deallocate(W)
   deallocate(Hk)
   deallocate(eigval)
   deallocate(dos)
   deallocate(dos_mpi)
   deallocate(omega)

   return
end subroutine dos_slab

subroutine joint_dos
! calculate joint density of state for 3D bulk system
!
! JDOS(\omega)= \sum_k (f_c(k)-f_v(k) \delta(\omega- Ec(k)+ Ev(k))

   use wmpi
   use para
   implicit none

   !> the integration k space
   real(dp) :: emin
   real(dp) :: emax
   real(dp) :: eta_broadening

   integer :: ik, ie, ib, ib1, ib2
   integer :: ikx, iky, ikz, knv3, NE, ierr

   !> integration for band
   integer :: iband_low, iband_high, iband_tot

   real(dp) :: x, dk3

   real(dp) :: k(3)

   real(dp), allocatable :: kpoints(:, :), eigval(:, :), eigval_mpi(:, :)
   real(dp), allocatable :: W(:), omega(:), jdos(:), jdos_mpi(:)
   complex(dp), allocatable :: Hk(:, :)

   !> fermi distribution
   real(dp), allocatable :: fermi_dis(:, :)

   !> delta function
   real(dp), external :: delta

   knv3= Nk1*Nk2*Nk3

   NE= OmegaNum
   iband_low= Numoccupied- 10
   iband_high= Numoccupied+ 10

   if (iband_low <1) iband_low = 1
   if (iband_high >Num_wann) iband_high = Num_wann

   iband_tot= iband_high- iband_low+ 1


   allocate(jdos(NE))
   allocate(jdos_mpi(NE))
   allocate(omega(NE))
   allocate(W(Num_wann))
   allocate(kpoints(3, knv3))
   allocate(Hk(Num_wann, Num_wann))
   allocate(eigval(iband_tot, knv3))
   allocate(eigval_mpi(iband_tot, knv3))
   allocate(fermi_dis(iband_tot, knv3))
   W= 0d0
   Hk= 0d0
   eigval= 0d0
   eigval_mpi= 0d0
   fermi_dis= 0d0
   kpoints= 0d0
   jdos= 0d0
   jdos_mpi= 0d0
   omega= 0d0

   ik =0

   do ikx= 1, nk1
      do iky= 1, nk2
         do ikz= 1, nk3
            ik= ik+ 1
            kpoints(:, ik)= K3D_start_cube+ K3D_vec1_cube*(ikx-1)/dble(nk1)  &
               + K3D_vec2_cube*(iky-1)/dble(nk2)  &
               + K3D_vec3_cube*(ikz-1)/dble(nk3)
         enddo
      enddo
   enddo

   dk3= kCubeVolume/dble(knv3)

   !> get eigenvalue
   do ik=1+cpuid, knv3, num_cpu
      if (cpuid.eq.0) write(stdout, *) 'ik, knv3', ik, knv3
      k= kpoints(:, ik)
      call ham_bulk_atomicgauge(k, Hk)
      W= 0d0
      call eigensystem_c( 'N', 'U', Num_wann ,Hk, W)
      eigval_mpi(:, ik)= W(iband_low:iband_high)
   enddo

#if defined (MPI)
   call mpi_allreduce(eigval_mpi,eigval,size(eigval),&
      mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
   eigval= eigval_mpi
#endif


   !> calculate fermi-dirac distribution
   do ik=1, knv3
      do ib=1, iband_tot
         if (eigval(ib, ik)<0) then
            fermi_dis(ib, ik)= 1d0
         else
            fermi_dis(ib, ik)= 0d0
         endif
      enddo
   enddo

   emin= 0d0
   emax= OmegaMax
   eta_broadening= (emax- emin)/ dble(NE)*3d0


   !> energy
   do ie=1, NE
      omega(ie)= emin+ (emax-emin)* (ie-1d0)/dble(NE-1)
   enddo ! ie

   !> get density of state
   jdos_mpi= 0d0
   do ie= 1, NE
      if (cpuid.eq.0) write(stdout, *)'ie, NE', ie, NE
      !> intergrate with k
      do ik= 1+cpuid, knv3, num_cpu
         do ib1= 1, iband_tot-1
            do ib2= ib1+1, iband_tot
               x= omega(ie)- eigval(ib2, ik) + eigval(ib1, ik)
               jdos_mpi(ie) = jdos_mpi(ie)+ delta(eta_broadening, x)* (fermi_dis(ib1, ik)- fermi_dis(ib2, ik))
            enddo ! ib2
         enddo ! ib1
      enddo ! ik
      jdos_mpi(ie)= jdos_mpi(ie)*dk3
   enddo ! ie

   jdos = 0d0
#if defined (MPI)
   call mpi_allreduce(jdos_mpi,jdos,size(jdos),&
      mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
   jdos= jdos_mpi
#endif

   outfileindex= outfileindex+ 1
   if (cpuid.eq.0) then
      open(unit=outfileindex, file='jdos.dat')
      do ie=1, NE
         write(outfileindex, *)omega(ie)/eV2Hartree, jdos(ie)
      enddo ! ie
      close(outfileindex)
   endif

#if defined (MPI)
   call mpi_barrier(mpi_cmw, ierr)
#endif
   deallocate(jdos)
   deallocate(jdos_mpi)
   deallocate(omega)
   deallocate(W)
   deallocate(kpoints)
   deallocate(Hk)
   deallocate(eigval)
   deallocate(eigval_mpi)
   deallocate(fermi_dis)


   return
end subroutine joint_dos


subroutine dos_joint_dos
!  calculate density of state and joint density of state for 3D bulk system
!
!  JDOS(\omega)= \sum_k (f_c(k)-f_v(k) \delta(\omega- Ec(k)+ Ev(k))

   use wmpi
   use para
   implicit none

   !> the integration k space
   real(dp) :: emin
   real(dp) :: emax
   real(dp) :: eta_broadening

   integer :: ik, ie, ib, ib1, ib2
   integer :: ikx, iky, ikz, knv3, NE, ierr

   !> integration for band
   integer :: iband_low, iband_high, iband_tot

   real(dp) :: x, dk3

   real(dp) :: k(3)
   real(dp), allocatable :: W(:), omega_dos(:), omega_jdos(:)
   real(dp), allocatable :: dos(:), dos_mpi(:), jdos(:), jdos_mpi(:)
   complex(dp), allocatable :: Hk(:, :)

   !> fermi distribution
   real(dp), allocatable :: fermi_dis(:)

   !> delta function
   real(dp), external :: delta

   knv3= Nk1*Nk2*Nk3

   NE= OmegaNum
   iband_low= Numoccupied- 40
   iband_high= Numoccupied+ 40

   if (iband_low <1) iband_low = 1
   if (iband_high >Num_wann) iband_high = Num_wann

   iband_tot= iband_high- iband_low+ 1


   allocate(dos(NE))
   allocate(dos_mpi(NE))
   allocate(jdos(NE))
   allocate(jdos_mpi(NE))
   allocate(omega_dos(NE))
   allocate(omega_jdos(NE))
   allocate(W(Num_wann))
   allocate(Hk(Num_wann, Num_wann))
   allocate(fermi_dis(Num_wann))
   W= 0d0
   Hk= 0d0
   fermi_dis= 0d0
   jdos= 0d0
   jdos_mpi= 0d0
   dos= 0d0
   dos_mpi= 0d0
   omega_dos= 0d0
   omega_jdos= 0d0


   dk3= kCubeVolume/dble(knv3)

   emin= 0d0
   emax= OmegaMax
   eta_broadening= (emax- emin)/ dble(NE)*5d0

   !> energy
   do ie=1, NE
      omega_jdos(ie)= emin+ (emax-emin)* (ie-1d0)/dble(NE-1)
   enddo ! ie

   emin= OmegaMin
   emax= OmegaMax

   !> energy
   do ie=1, NE
      omega_dos(ie)= emin+ (emax-emin)* (ie-1d0)/dble(NE-1)
   enddo ! ie


   !> get eigenvalue
   dos_mpi= 0d0
   jdos_mpi= 0d0
   do ik=1+cpuid, knv3, num_cpu
      if (cpuid.eq.0) write(stdout, *) 'ik, knv3', ik, knv3
      ikx= (ik- 1)/(Nk2*Nk3)+ 1
      iky= (ik- (ikx-1)*Nk2*Nk3- 1)/Nk3+ 1
      ikz= ik- (ikx-1)*Nk2*Nk3- (iky-1)*Nk3

      k= K3D_start_cube+ K3D_vec1_cube*(ikx-1)/dble(nk1-1)  &
         + K3D_vec2_cube*(iky-1)/dble(nk2-1)  &
         + K3D_vec3_cube*(ikz-1)/dble(nk3-1)
      call ham_bulk_atomicgauge(k, Hk)
      W= 0d0
      call eigensystem_c( 'N', 'U', Num_wann ,Hk, W)

      !> calculate fermi-dirac distribution
      do ib=iband_low, iband_high
         if (W(ib)<0) then
            fermi_dis(ib)= 1d0
         else
            fermi_dis(ib)= 0d0
         endif
      enddo !ib

      !> get density of state
      do ie= 1, NE
         do ib1= iband_low, iband_high-1
            do ib2= ib1+1, iband_high
               x= omega_jdos(ie)- W(ib2) + W(ib1)
               jdos_mpi(ie)= jdos_mpi(ie)+ delta(eta_broadening, x)* (fermi_dis(ib1)- fermi_dis(ib2))
            enddo ! ib2
         enddo ! ib1
      enddo ! ie

      !> get density of state
      do ie= 1, NE
         !> intergrate with k
         do ib= iband_low, iband_high-1
            x= omega_dos(ie)- W(ib)
            dos_mpi(ie) = dos_mpi(ie)+ delta(eta_broadening, x)
         enddo ! ib
      enddo ! ie

   enddo ! ik

#if defined (MPI)
   call mpi_allreduce(dos_mpi,dos,size(dos),&
      mpi_dp,mpi_sum,mpi_cmw,ierr)
   call mpi_allreduce(jdos_mpi,jdos,size(jdos),&
      mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
   dos= dos_mpi
   jdos= jdos_mpi
#endif


   outfileindex= outfileindex+ 1
   if (cpuid.eq.0) then
      open(unit=outfileindex, file='jdos.dat')
      do ie=1, NE
         write(outfileindex, *)omega_jdos(ie)/eV2Hartree, jdos(ie)*dk3*eV2Hartree
      enddo ! ie
      close(outfileindex)
   endif

   outfileindex= outfileindex+ 1
   if (cpuid.eq.0) then
      open(unit=outfileindex, file='dos.dat')
      do ie=1, NE
         write(outfileindex, *)omega_dos(ie)/eV2Hartree, dos(ie)*dk3*eV2Hartree
      enddo ! ie
      close(outfileindex)
   endif
#if defined (MPI)
   call mpi_barrier(mpi_cmw, ierr)
#endif

   deallocate(dos)
   deallocate(dos_mpi)
   deallocate(jdos)
   deallocate(jdos_mpi)
   deallocate(omega_dos)
   deallocate(omega_jdos)
   deallocate(W)
   deallocate(Hk)
   deallocate(fermi_dis)

   return
end subroutine dos_joint_dos


function delta(eta, x)
   !>  Lorentz or gaussian expansion of the Delta function
   use para, only : dp, pi
   implicit none
   real(dp), intent(in) :: eta
   real(dp), intent(in) :: x
   real(dp) :: delta, y

   !> Lorentz expansion
   !delta= 1d0/pi*eta/(eta*eta+x*x)

   y= x*x/eta/eta/2d0

   !> Gaussian broadening
   !> exp(-60) = 8.75651076269652e-27
   if (y>60d0) then
      delta = 0d0
   else
      delta= exp(-y)/sqrt(2d0*pi)/eta
   endif

   return
end function delta
