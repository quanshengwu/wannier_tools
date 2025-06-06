!> There are three modes to calculate the energy bands in the BZ
!> 1. Line mode : k points are defined in the KPATH_BULK, subroutine ek_bulk_line
!> 2. Plane mode: k points are defined in the KPLANE_BULK, subroutine ek_bulk_plane
!> 3. Cube mode : k points are defined in the KCUBE_BULK, subroutine ek_bulk_cube

subroutine ek_bulk_line
   ! Calculate bulk's energy bands using wannier TB method
   ! Line mode
   ! Copyright (c) 2010 QuanSheng Wu. All rights reserved.

   use wmpi
   use para

   implicit none

   integer :: ik, il, ig, io, i, j, knv3, ierr
   real(dp) :: emin,  emax,  k(3)
   character*40 :: filename

   !> eigenvalues of H
   real(Dp), allocatable :: W(:)

   ! Hamiltonian of bulk system
   complex(Dp), allocatable :: Hamk_bulk(:, :), Sk_bulk(:, :)

   ! eigenectors of H
   real(dp), allocatable :: eigv(:,:), eigv_mpi(:,:)
   real(dp), allocatable :: weight(:,:,:), weight_mpi(:,:,:), weight_sum(:,:)

   knv3= nk3_band
   allocate(W(Num_wann))
   allocate(Hamk_bulk(Num_wann, Num_wann))
   allocate( eigv    (Num_wann, knv3))
   allocate( eigv_mpi(Num_wann, knv3))
   allocate( weight    (NumberofSelectedOrbitals_groups,Num_wann, knv3))
   allocate( weight_mpi(NumberofSelectedOrbitals_groups,Num_wann, knv3))
   allocate( weight_sum(Num_wann, knv3))
   W       = 0d0; Hamk_bulk = 0d0
   eigv    = 0d0; eigv_mpi= 0d0
   weight  = 0d0; weight_sum = 0d0; weight_mpi = 0d0

   if (.not. Orthogonal_Basis) then
      allocate(Sk_bulk(Num_wann, Num_wann))
      Sk_bulk= 0d0
   endif

   do ik= 1+cpuid, knv3, num_cpu

      k = kpath_3d(:, ik)
      
      ! calculation bulk hamiltonian
      Hamk_bulk= 0d0

      ! generate bulk Hamiltonian
      if (index(KPorTB, 'KP')/=0)then
         call ham_bulk_kp_abcb_graphene(k, Hamk_bulk)
      else
         !> deal with phonon system
         if (index(Particle,'phonon')/=0.and.LOTO_correction) then
            call ham_bulk_LOTO(k, Hamk_bulk)
         else
            if (Is_Sparse) then
               stop 'ERROR: Sparse mode is not supported in ek_bulk_line'
            else
               if (Orthogonal_Basis) then
                  call ham_bulk_latticegauge(k, Hamk_bulk)
               else
                  
                  Sk_bulk= 0d0
                  call S_bulk_latticegauge(k, Sk_bulk)
                  call ham_bulk_latticegauge(k, Hamk_bulk)
                  call orthogonalize_hamiltonian(Hamk_bulk, Sk_bulk, Num_wann)
               endif

            endif
           !call ham_bulk_atomicgauge(k, Hamk_bulk)
         endif
      endif

      !> diagonalization by call zheev in lapack
      W= 0d0
      call eigensystem_c('V', 'U', Num_wann ,Hamk_bulk, W)
      eigv(:, ik)= W

      do j= 1, Num_wann  !> band
           do ig= 1, NumberofSelectedOrbitals_groups
              do i= 1, NumberofSelectedOrbitals(ig)
                 io= Selected_WannierOrbitals(ig)%iarray(i)
                 weight(ig, j, ik)= weight(ig, j, ik)+ abs(Hamk_bulk(io, j))**2 
              enddo !i
           enddo !ig
      enddo !j
   enddo !ik

#if defined (MPI)
   call mpi_allreduce(eigv,eigv_mpi,size(eigv),&
      mpi_dp,mpi_sum,mpi_cmw,ierr)
   call mpi_allreduce(weight, weight_mpi,size(weight),&
      mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
   eigv_mpi= eigv
   weight_mpi= weight
#endif

   if (index(Particle,'phonon')/=0) then
      eigv_mpi = eigv_mpi - MINVAL(eigv_mpi)
   endif

   !> deal with phonon system
   if (index(Particle,'phonon')/=0) then
      do ik=1, knv3
         do j=1, Num_wann
            eigv_mpi(j, ik)= sqrt(abs(eigv_mpi(j, ik)))*sign(1d0, eigv_mpi(j, ik))
            !eigv_mpi(j, ik)=     (abs(eigv_mpi(j, ik)))*sign(1d0, eigv_mpi(j, ik))
         enddo
      enddo
   endif
   eigv_mpi= eigv_mpi/eV2Hartree

   outfileindex= outfileindex+ 1
   if (cpuid==0)then
      do il= 1, nk3lines
         if (il<10) then
            write(filename, '(a,i1)')'bulkek.dat-segment', il
         elseif (il>=10.and.il<100) then
            write(filename, '(a,i2)')'bulkek.dat-segment', il
         elseif (il>=100.and.il<1000) then
            write(filename, '(a,i3)')'bulkek.dat-segment', il
         endif
         write(filename, '(5a)')'bulkek.dat-', trim(adjustl(k3line_name(il))), '-', trim(adjustl(k3line_name(il+1)))
         outfileindex= outfileindex+ 1
         open(unit=outfileindex, file=filename)
         write(outfileindex, '(a, a6, a, a6)')'# segment between', k3line_name(il), ' and ',  k3line_name(il+1)
         write(outfileindex, '(a, 3f10.6)')'# kstart (fractional units)', k3line_start(:, il)!*Angstrom2atomic   !fractional units
         write(outfileindex, '(a, 3f10.6)')'# kend   (fractional units)', k3line_end(:, il)!*Angstrom2atomic
   
         do i=1, Num_wann
            do ik=1+(il-1)*Nk1, il*Nk1
               write(outfileindex, '(2f19.9, 10000i5)')(k3len(ik)*Angstrom2atomic-k3len((il-1)*Nk1+1)*Angstrom2atomic),eigv_mpi(i, ik)
            enddo
            write(outfileindex, *)' '
         enddo ! i
         close(outfileindex)
      enddo ! il
   endif

   outfileindex= outfileindex+ nk3lines+1
   if (cpuid==0)then
      open(unit=outfileindex, file='bulkek.dat')

      write(outfileindex, &
         "('#', a12, a14, 1X, '| projection', 100(3X,'|group', i2, ': A '))")&
         'klen', 'E', (i, i=1, NumberofSelectedOrbitals_groups)
      write(outfileindex, "('#column', i5, 200i16)")(i, i=1, 2+NumberofSelectedOrbitals_groups)
      do i=1, Num_wann
         do ik=1, knv3
            write(outfileindex, '(200f16.9)')k3len(ik)*Angstrom2atomic,eigv_mpi(i, ik), &
               weight_mpi(:, i, ik)
         enddo
         write(outfileindex, *)' '
      enddo
      close(outfileindex)
   endif

   outfileindex= outfileindex+ nk3lines+1
   if (cpuid==0)then
      open(unit=outfileindex, file='bulkek.dat-matlab')
      write(outfileindex, '(2a19)')'% klen', 'E(n)'

      do ik=1, knv3
         write(outfileindex, '(2000f19.9)')k3len(ik)*Angstrom2atomic,eigv_mpi(:, ik)
      enddo
      close(outfileindex)
   endif


   !> minimum and maximum value of energy bands
   emin=  minval(eigv_mpi)-0.5d0
   emax=  maxval(eigv_mpi)+0.5d0

   call generate_ek_kpath_gnu('bulkek.dat', 'bulkek.gnu', 'bulkek.pdf', &
                                 emin, emax, knv3, Nk3lines, &
                                 k3line_name, k3line_stop, k3len)

   deallocate(W)
   deallocate(Hamk_bulk)
   deallocate( eigv    )
   deallocate( eigv_mpi)
   deallocate( weight    )
   deallocate( weight_mpi)
   deallocate( weight_sum)

   return
end subroutine ek_bulk_line


subroutine ek_bulk_line_valley
   ! Calculate bulk's energy bands using wannier TB method
   ! Line mode
   ! Copyright (c) 2010 QuanSheng Wu. All rights reserved.

   use wmpi
   use para

   implicit none

   integer :: ik, il, ig, io, i, j, knv3, ierr
   integer :: ie, ND, NDMAX, ie1, ie2
   real(dp) :: emin,  emax,  k(3)
   character*40 :: filename
   complex(dp), external :: zdotc

   !> eigenvalues of H
   real(Dp), allocatable :: W(:)

   ! Hamiltonian of bulk system
   complex(Dp), allocatable :: psi(:), vpsi(:), valley_k_nd(:, :), psi1(:), psi2(:)
   complex(Dp), allocatable :: Hamk_bulk(:, :), valley_k(:, :)
   complex(Dp), allocatable :: VL(:, :), VR(:, :), valley_eig(:)
   logical, allocatable :: valley_plus(:, :), valley_plus_mpi(:, :)

   ! eigenectors of H
   real(dp), allocatable :: eigv(:,:), eigv_mpi(:,:)
   real(dp), allocatable :: weight(:,:), weight_mpi(:,:)
   real(dp) :: tolde= 1E-4

   NDMAX= 12
   knv3= nk3_band
   allocate(W(Num_wann))
   allocate(Hamk_bulk(Num_wann, Num_wann))
   allocate( eigv    (Num_wann, knv3))
   allocate( eigv_mpi(Num_wann, knv3))
   allocate( weight    (Num_wann, knv3))
   allocate( weight_mpi(Num_wann, knv3))
   allocate( psi(Num_wann), vpsi(Num_wann))
   allocate( psi1(Num_wann), psi2(Num_wann))
   allocate(valley_k(Num_wann, Num_wann))
   allocate(valley_plus(Num_wann, knv3))
   allocate(valley_plus_mpi(Num_wann, knv3))
   allocate( valley_k_nd(NDMAX, NDMAX), valley_eig(NDMAX))
   allocate( VL(NDMAX, NDMAX), VR(NDMAX, NDMAX))
   W       = 0d0; Hamk_bulk = 0d0
   eigv    = 0d0; eigv_mpi= 0d0
   weight  = 0d0; weight_mpi = 0d0
   psi     = 0d0; vpsi= 0d0
   valley_k_nd= 0d0
   valley_plus= .False.; valley_plus_mpi= .False.


   do ik= 1+cpuid, knv3, num_cpu

      k = kpath_3d(:, ik)
      
      ! calculation bulk hamiltonian
      Hamk_bulk= 0d0

      call ham_bulk_atomicgauge(k, Hamk_bulk)

      !> diagonalization by call zheev in lapack
      W= 0d0
      call eigensystem_c('V', 'U', Num_wann ,Hamk_bulk, W)

      !> get the fourier transform of a valley operator
      call valley_k_atomicgauge(k, valley_k)

      eigv(:, ik)= W

      !> check the degeneracy of each band
      IE=1
      do while (ie.le.Num_wann)
         ND=1
         if (ie+ND.le.Num_wann) then
            do while ((ie+ND).le.Num_wann .and. (W(ie+ND)- W(ie)).lt.tolde)
               if (ie+ND .ge. Num_wann) exit
               ND= ND+1
            enddo
         endif

         !> if the degeneracy is larger than 1, we need to calculate the matrix
         valley_k_nd= 0d0
         do ie1= 1, ND
            psi1= hamk_bulk(:, IE+ie1-1)
            do ie2= 1, ND
               psi2= hamk_bulk(:, IE+ie2-1)
               vpsi=0d0
               call zgemv('N', Num_wann, Num_wann, One_complex,  valley_k, Num_wann, psi2, 1, zzero, vpsi, 1)
               valley_k_nd(ie1, ie2)= zdotc(Num_wann, psi1, 1, vpsi, 1)
            enddo
           !if (abs(W(IE+ie1-1)/eV2Hartree)<2d0) &
           !write( *, '(a, 2i5, 200f8.3)') "ND, ie1: ", ND, ie1, valley_k_nd(ie1, 1:ND)
         enddo

         !> diagonalize the matrix
         VL = 0d0; VR= 0d0
         if (ND>1) then
            call zgeev_sys(ND, valley_k_nd(1:ND, 1:ND), valley_eig(1:ND),'N',VL(1:ND, 1:ND),"V",VR(1:ND, 1:ND) )
            do ie1= 1, ND
               weight(IE+ie1-1, ik) = real(valley_eig(ie1))
            enddo
         else
            valley_eig(1)= valley_k_nd(1, 1)
            weight(IE, ik) = real(valley_k_nd(1, 1))
         endif
         !  if (abs(W(IE+ie1-1)/eV2Hartree)<2d0) &
         ! write(*, '(a, 2i5, a, 200f8.3)'), 'ND, ie', ND, ie, '  valley_eig', real(valley_eig(1:ND))
         !  if (abs(W(IE+ie1-1)/eV2Hartree)<2d0) &
         ! print *, ' '

         !> get new eigenvector
        !do ie1= 1, ND
        !   psi= 0d0
        !   do ie2= 1, ND
        !      psi= psi+ VR(ie2,  ie1)* ham_bulk(:, IE+ie2-1)
        !   enddo
        !enddo

         !> if the degeneracy is larger than 1, we need to calculate the matrix
         IE= IE+ ND
      enddo
      do ie=1, Num_wann
         if (weight(ie, ik)>0) valley_plus(ie, ik)=.true.
      enddo

   enddo ! ik


#if defined (MPI)
   call mpi_allreduce(eigv,eigv_mpi,size(eigv),&
      mpi_dp,mpi_sum,mpi_cmw,ierr)
   call mpi_allreduce(weight, weight_mpi,size(weight),&
      mpi_dp,mpi_sum,mpi_cmw,ierr)
   call mpi_allreduce(valley_plus, valley_plus_mpi,size(valley_plus),&
      mpi_logical,mpi_lor,mpi_cmw,ierr)
#else
   eigv_mpi= eigv
   weight_mpi= weight
   valley_plus_mpi= valley_plus
#endif
   eigv_mpi= eigv_mpi/eV2Hartree

   outfileindex= outfileindex+ 1
   if (cpuid==0)then
      open(unit=outfileindex, file='bulkek.dat')
      do i=1, Num_wann
         do ik=1, knv3
            write(outfileindex, '(200F16.8)')k3len(ik)*Angstrom2atomic,eigv_mpi(i, ik), &
               weight_mpi(i, ik)
         enddo
         write(outfileindex, *)' '
      enddo
      close(outfileindex)
   endif
   outfileindex= outfileindex+ 1
   if (cpuid==0)then
      open(unit=outfileindex, file='bulkek_valley_plus.dat')
      do i=1, Num_wann
         do ik=1, knv3
            if (valley_plus_mpi(i, ik)) then
               write(outfileindex, '(200F16.8)')k3len(ik)*Angstrom2atomic,eigv_mpi(i, ik), &
                  weight_mpi(i, ik)
            endif
         enddo
         write(outfileindex, *)' '
      enddo
      close(outfileindex)
   endif

   outfileindex= outfileindex+ 1
   if (cpuid==0)then
      open(unit=outfileindex, file='bulkek_valley_minus.dat')
      do i=1, Num_wann
         do ik=1, knv3
            if (.not.valley_plus_mpi(i, ik)) then
               write(outfileindex, '(200F16.8)')k3len(ik)*Angstrom2atomic,eigv_mpi(i, ik), &
                  weight_mpi(i, ik)
            endif
         enddo
         write(outfileindex, *)' '
      enddo
      close(outfileindex)
   endif



   !> minimum and maximum value of energy bands
   emin=  minval(eigv_mpi)-0.5d0
   emax=  maxval(eigv_mpi)+0.5d0

   call generate_ek_kpath_gnu('bulkek.dat', 'bulkek.gnu', 'bulkek.pdf', &
                                 emin, emax, knv3, Nk3lines, &
                                 k3line_name, k3line_stop, k3len)

   deallocate(W)
   deallocate(Hamk_bulk)
   deallocate( eigv    )
   deallocate( eigv_mpi)
   deallocate( weight    )
   deallocate( weight_mpi)

   return
end subroutine ek_bulk_line_valley

  !>> calculate energy band levels at given kpoints
  subroutine ek_bulk_point_mode

     use wmpi
     use para

     implicit none

     integer :: ik, i, j, knv3, ierr
     real(dp) :: k(3)
     real(Dp), allocatable :: W(:)
     
     ! Hamiltonian of bulk system
     complex(Dp), allocatable :: Hamk_bulk(:, :)

     ! eigen value of H
	  real(dp), allocatable :: eigv(:,:)
	  real(dp), allocatable :: eigv_mpi(:,:)
	  real(dp), allocatable :: weight(:,:,:)
	  real(dp), allocatable :: weight_mpi(:,:,:)

     knv3= Nk3_point_mode
     allocate( W(Num_wann))
     allocate( eigv    (Num_wann, knv3))
     allocate( eigv_mpi(Num_wann, knv3))
     allocate( Hamk_bulk (Num_wann,Num_wann))
     allocate( weight    (Num_wann,Num_wann, knv3))
     allocate( weight_mpi(Num_wann,Num_wann, knv3))
     eigv    = 0d0
     weight  = 0d0
     eigv_mpi= 0d0
     weight_mpi = 0d0

     do ik= 1+cpuid, Nk3_point_mode, num_cpu
        if (cpuid==0) write(stdout, *)'# BulkBand in point mode, ik, knv3 ', ik, knv3

        k = k3points_pointmode_direct(:, ik)

        !> calculation bulk hamiltonian
        Hamk_bulk= 0d0
        call ham_bulk_latticegauge(k, Hamk_bulk)
       !call ham_bulk    (k, Hamk_bulk)

        !> diagonalization by call zheev in lapack
        W= 0d0
        call eigensystem_c( 'V', 'U', Num_wann ,Hamk_bulk, W)

        eigv(:, ik)= W
        do i=1, Num_wann  !> band 
           if (SOC==0) then
              do j=1, Num_wann  !> projector
                 weight(j, i, ik)= (abs(Hamk_bulk(j, i))**2)
              enddo ! j
           else
              do j=1, Num_wann  !> projector
                 weight(j, i, ik)= abs(Hamk_bulk(j, i))**2     
              enddo ! j
           endif
        enddo ! i
 
     enddo ! ik

#if defined (MPI)
     call mpi_allreduce(eigv,eigv_mpi,size(eigv),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
     call mpi_allreduce(weight, weight_mpi,size(weight),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
     eigv_mpi= eigv
     weight_mpi= weight
#endif

     eigv_mpi= eigv_mpi/eV2Hartree
     weight= weight_mpi/maxval(weight_mpi)*255d0

     outfileindex= outfileindex+ 1
     if (cpuid==0)then
        open(unit=outfileindex, file='bulkek-pointsmode.dat')
   
        do ik=1, Nk3_point_mode
           write(outfileindex, '(a, i10)')'# No. of k point', ik
           write(outfileindex, '("#",a9,100a10)')'k1', 'k2', 'k3', 'kx', 'ky', 'kz'  
           write(outfileindex, '(100f10.6)')k3points_pointmode_direct(:, ik), k3points_pointmode_cart(:, ik)*Angstrom2atomic
           write(outfileindex, '("#", a11, a19, a)')'band index', 'Eigenvalue', '     orbital weights (0-255)'
           do i=1, Num_wann
              write(outfileindex, '(i12, f19.10, 1000i5)')i, eigv_mpi(i, ik), &
                 int(weight(:, i, ik))
           enddo
           write(outfileindex, *)' '
        enddo
        close(outfileindex)
     endif

     return
   end subroutine ek_bulk_point_mode



subroutine ek_bulk_plane_C2yT
   ! Calculate bulk's energy band using wannier TB method for a given k plane
   ! Copyright (c) 2017 QuanSheng Wu. All rights reserved.
   ! Revised by QS.Wu on Sep 19 2017 @EPFL, Switzerland

   use wmpi
   use para

   implicit none

   integer :: ik, i, j, ib, ik1, ik2
   integer :: knv3
   integer :: ierr
   integer :: nwann

   integer :: nband_min
   integer :: nband_max
   integer :: nband_store

   real(dp) :: time_start, time_end

   real(Dp) :: k(3), kxyz(3)
   real(Dp), allocatable :: W(:), D(:)
   real(dp) :: kxmin, kxmax, kymin, kymax
   real(dp), allocatable :: kxy(:,:)
   real(dp), allocatable :: kxy_shape(:,:)
   real(dp), allocatable :: kxy_plane(:,:)
   
   ! Hamiltonian of bulk system
   complex(Dp), allocatable :: Hamk_bulk(:, :)
   real(Dp), allocatable :: eigenvectors(:, :, :)
   real(Dp), allocatable :: eigenvectors_mpi(:, :, :)

   ! eigen value of H
   real(dp), allocatable :: gap(:)
   real(dp), allocatable :: gap_mpi(:)
   real(dp), allocatable :: eigv(:,:)
   real(dp), allocatable :: eigv_mpi(:,:)

   complex(dp), allocatable :: U(:, :)
   complex(dp), allocatable :: sqrt_C2yT(:, :)
   complex(dp), allocatable :: psi_tilde(:, :)

   allocate(D(Num_wann), U(Num_wann, Num_wann))
   D=0d0; U=0d0


  !call TakagiFactor(Num_wann, C2yT, num_wann, D, U, num_wann, 0, 1)
  !write(*, *)'real(C2yT)'
  !do i=1, num_wann
  !   write(*, '(100f6.2)') real(C2yT(i, :))
  !enddo

  !write(*, *)'imag(C2yT)'
  !do i=1, num_wann
  !   write(*, '(100f6.2)') aimag(C2yT(i, :))
  !enddo

  !write(*, *)'D'
  !do i=1, num_wann
  !   write(*, '(100f6.2)') D(i)
  !enddo

  !write(*, *)'real(U)'
  !do i=1, num_wann
  !   write(*, '(100f6.2)') real(U(i, :))
  !enddo

  !write(*, *)'imag(U)'
  !do i=1, num_wann
  !   write(*, '(100f6.2)') aimag(U(i, :))
  !enddo
  !stop


   if (.not.Symmetry_Import_calc) stop "Please set Symmetry_Import_calc=T in the input.dat"

   allocate(sqrt_C2yT(Num_wann, Num_wann))
   allocate(psi_tilde(Num_wann, Num_wann))
   psi_tilde= 0d0
   sqrt_C2yT= 0d0

   do i=1, Num_wann
      sqrt_C2yT(i, i)= sqrt(C2yT(i, i))
   enddo

   nband_min= Numoccupied- 1
   nband_max= Numoccupied+ 2
   if (nband_min< 1) nband_min= 1
   if (nband_max> Num_wann) nband_max= Num_wann
   nband_store= nband_max- nband_min+ 1


   knv3= nk1*Nk2
   allocate(W(Num_wann))
   allocate(Hamk_bulk(Num_wann, Num_wann))
   allocate(eigenvectors(Num_wann, 2, knv3))
   allocate(eigenvectors_mpi(Num_wann, 2, knv3))
   allocate( kxy(3, nk1*Nk2))
   allocate( kxy_shape(3, nk1*Nk2))
   allocate( kxy_plane(3, nk1*Nk2))
   allocate( eigv    (nband_store, knv3))
   allocate( eigv_mpi(nband_store, knv3))

   eigv    = 0d0
   eigv_mpi= 0d0
   eigenvectors= 0d0
   eigenvectors_mpi= 0d0

   kxy=0d0
   kxy_shape=0d0
   kxy_plane=0d0
  
   if (Nk1<2 .or. Nk2<2) stop 'ERROR: I refuse to do this job because you give me so small Nk1 and Nk2'      
   if (Numoccupied> Num_wann) stop 'ERROR: please set correct Numoccupied, it should be small than num_wann'
   
   ik =0
   do i= 1, Nk1
      do j= 1, Nk2
         ik =ik +1
         kxy(:, ik)= K3D_start+ K3D_vec1*(i-1)/dble(Nk1-1)+ K3D_vec2*(j-1)/dble(Nk2-1) &
            -(K3D_vec1+K3D_vec2)/2d0
         kxy_shape(:, ik)= kxy(1, ik)* Origin_cell%Kua+ kxy(2, ik)* Origin_cell%Kub+ kxy(3, ik)* Origin_cell%Kuc 
         call rotate_k3_to_kplane(kxy_shape(:, ik), kxy_plane(:, ik))
      enddo
   enddo

   time_start= 0d0
   time_end= 0d0
   do ik= 1+cpuid, knv3, num_cpu
      if (cpuid==0.and. mod(ik/num_cpu, 500)==0) &
         write(stdout, '(a, i12, a, i12, a, f10.2, a)') &
         'ek_bulk_plane, ik ', ik, ' knv3',knv3, ' time left', &
         (knv3-ik)*(time_end- time_start)/num_cpu, ' s'
      call now(time_start)

      k = kxy(:, ik)
      call direct_cart_rec(k, kxyz)

      ! calculation bulk hamiltonian
      Hamk_bulk= 0d0
      call ham_bulk_atomicgauge(k, Hamk_bulk)

      !> diagonalization by call zheev in lapack
      W= 0d0
      call eigensystem_c( 'V', 'U', Num_wann ,Hamk_bulk, W)
      call mat_mul(Num_wann, sqrt_C2yT, Hamk_bulk, psi_tilde)
      eigenvectors(:, 1:2, ik)= real(psi_tilde(:, Numoccupied:Numoccupied+1))

     !write(*, *)'psi real'
     !do i=1, num_wann
     !   write(*, '(100f6.2)') real(Hamk_bulk(i, :))
     !enddo
   
     !write(*, *)'psi imag'
     !do i=1, num_wann
     !   write(*, '(100f6.2)') aimag(Hamk_bulk(i, :))
     !enddo
    

     !write(*, *)'psi_tilde real'
     !do i=1, num_wann
     !   write(*, '(100f6.2)') real(psi_tilde(i, :))
     !enddo
   
     !write(*, *)'psi_tilde imag'
     !do i=1, num_wann
     !   write(*, '(100f6.2)') aimag(psi_tilde(i, :))
     !enddo
     !stop


      eigv(:, ik)= W(nband_min:nband_max)

      call now(time_end)
   enddo ! ik

#if defined (MPI)
   call mpi_allreduce(eigv,eigv_mpi,size(eigv),&
                     mpi_dp,mpi_sum,mpi_cmw,ierr)
   call mpi_allreduce(eigenvectors, eigenvectors_mpi,size(eigenvectors),&
                     mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
   eigv_mpi= eigv
#endif
   eigv_mpi= eigv_mpi/eV2Hartree

   !> write out the data in gnuplot format
   outfileindex= outfileindex+ 1
   if (cpuid==0)then
      open(unit=outfileindex, file='eigenvectors_plane_band1.txt')
      do i=1, num_wann
         ik=0
         do ik1=1, Nk1
            do ik2=1, Nk2-1
               ik=ik+1
               write(outfileindex, '(2000f19.9)', advance='no')eigenvectors_mpi(i, 1, ik)
            enddo
            ik=ik+1
            write(outfileindex, '(2000f19.9)', advance='yes')eigenvectors_mpi(i, 1, ik)
         enddo
         write(outfileindex, *)' '
      enddo
      close(outfileindex)
   endif

   !> write out the data in gnuplot format
   outfileindex= outfileindex+ 1
   if (cpuid==0)then
      open(unit=outfileindex, file='eigenvectors_plane_band2.txt')
      do i=1, num_wann
         ik=0
         do ik1=1, Nk1
            do ik2=1, Nk2-1
               ik=ik+1
               write(outfileindex, '(2000f19.9)', advance='no')eigenvectors_mpi(i, 2, ik)
            enddo
            ik=ik+1
            write(outfileindex, '(2000f19.9)', advance='yes')eigenvectors_mpi(i, 2, ik)
         enddo
         write(outfileindex, *)' '
      enddo
      close(outfileindex)
   endif


   !> write out the data in gnuplot format
   outfileindex= outfileindex+ 1
   if (cpuid==0)then
      open(unit=outfileindex, file='bulkek_plane.dat')
      write(outfileindex, '(1000a19)')'# kx', 'ky', 'kz', 'k1', 'k2', 'k3', &
         'E(Numoccupied-1)', 'E(Numoccupied)' , 'E(Numoccupied+1)', 'E(Numoccupied+2)'
      do ik=1, knv3
         write(outfileindex, '(1000f19.9)')kxy_shape(:, ik)*Angstrom2atomic, &
            kxy_plane(:, ik)*Angstrom2atomic, eigv_mpi(:, ik)  
         if (mod(ik, nk2)==0) write(outfileindex, *)' '
      enddo
      close(outfileindex)
   endif

   !> write out the data in gnuplot format
   outfileindex= outfileindex+ 1
   if (cpuid==0)then
      open(unit=outfileindex, file='bulkek_plane-matlab.dat')
      write(outfileindex, '(1000a19)')'% kx', 'ky', 'kz', 'k1', 'k2', &
         'E(Numoccupied-1)', 'E(Numoccupied)' , 'E(Numoccupied+1)', 'E(Numoccupied+2)'
      do ik=1, knv3
         write(outfileindex, '(1000f19.9)')kxy_shape(:, ik)*Angstrom2atomic, &
            kxy_plane(:, ik)*Angstrom2atomic, eigv_mpi(:, ik)  
      enddo
      close(outfileindex)
   endif

   !> write out a script that can be used for gnuplot
   outfileindex= outfileindex+ 1
   if (cpuid==0)then
      open(unit=outfileindex, file='bulkek_plane.gnu')
      write(outfileindex, '(a)')"set encoding iso_8859_1"
      write(outfileindex, '(a)')'#set terminal  postscript enhanced color'
      write(outfileindex, '(a)')"#set output 'bulkek_plane.eps'"
      write(outfileindex, '(3a)')'set terminal  png truecolor enhanced', &
         ' size 1920, 1680 font ",36"'
      write(outfileindex, '(a)')"set output 'bulkek_plane.png'"
      write(outfileindex, '(a)')'set palette rgbformulae 33,13,10'
      write(outfileindex, '(a)')'unset key'
      write(outfileindex, '(a)')'set pm3d'
      write(outfileindex, '(a)')'set origin 0.2, 0'
      write(outfileindex, '(a)')'set size 0.8, 1'
      write(outfileindex, '(a)')'set border lw 3'
      write(outfileindex, '(a)')'#set xtics font ",24"'
      write(outfileindex, '(a)')'#set ytics font ",24"'
      write(outfileindex, '(a)')'set size ratio -1'
      write(outfileindex, '(a)')'set xtics'
      write(outfileindex, '(a)')'set ytics'
      write(outfileindex, '(a)')'set view 80,60'
      write(outfileindex, '(a)')'set xlabel "k_1"'
      write(outfileindex, '(a)')'set ylabel "k_2"'
      write(outfileindex, '(a)')'set zlabel "Energy (eV)" rotate by 90'
      write(outfileindex, '(a)')'unset colorbox'
      write(outfileindex, '(a)')'set autoscale fix'
      write(outfileindex, '(a)')'set pm3d interpolate 4,4'
      write(outfileindex, '(2a)')"splot 'bulkek_plane.dat' u 4:5:8 w pm3d, \"
      write(outfileindex, '(2a)')"      'bulkek_plane.dat' u 4:5:9 w pm3d"

      close(outfileindex)

   endif ! cpuid

   deallocate(W)
   deallocate(Hamk_bulk)
   deallocate( kxy)
   deallocate( kxy_shape)
   deallocate( kxy_plane)
   deallocate( eigv    )
   deallocate( eigv_mpi)

   return
end subroutine ek_bulk_plane_C2yT


subroutine ek_bulk_plane
   ! Calculate bulk's energy band using wannier TB method for a given k plane
   ! Copyright (c) 2017 QuanSheng Wu. All rights reserved.
   ! Revised by QS.Wu on Sep 19 2017 @EPFL, Switzerland
   ! This subroutine is usefull for Dirac cone plot

   use wmpi
   use para

   implicit none

   integer :: ik, i, j, ib
   integer :: knv3
   integer :: ierr
   integer :: nwann

   integer :: nband_min
   integer :: nband_max
   integer :: nband_store

   real(dp) :: time_start, time_end

   real(Dp) :: k(3)
   real(Dp), allocatable :: W(:)
   real(dp) :: kxmin, kxmax, kymin, kymax
   real(dp), allocatable :: kxy(:,:)
   real(dp), allocatable :: kxy_shape(:,:)
   real(dp), allocatable :: kxy_plane(:,:)
   
   ! Hamiltonian of bulk system
   complex(Dp), allocatable :: Hamk_bulk(:, :)

   ! eigen value of H
   real(dp), allocatable :: gap(:)
   real(dp), allocatable :: gap_mpi(:)
   real(dp), allocatable :: eigv(:,:)
   real(dp), allocatable :: eigv_mpi(:,:)

   nband_min= Numoccupied- 1
   nband_max= Numoccupied+ 2
   if (nband_min< 1) nband_min= 1
   if (nband_max> Num_wann) nband_max= Num_wann
   nband_store= nband_max- nband_min+ 1


   knv3= nk1*Nk2
   allocate(W(Num_wann))
   allocate(Hamk_bulk(Num_wann, Num_wann))
   allocate( kxy(3, nk1*Nk2))
   allocate( kxy_shape(3, nk1*Nk2))
   allocate( kxy_plane(3, nk1*Nk2))
   allocate( eigv    (nband_store, knv3))
   allocate( eigv_mpi(nband_store, knv3))

   eigv    = 0d0
   eigv_mpi= 0d0

   kxy=0d0
   kxy_shape=0d0
   kxy_plane=0d0
  
   if (Nk1<2 .or. Nk2<2) stop 'ERROR: I refuse to do this job because you give me so small Nk1 and Nk2'      
   if (Numoccupied> Num_wann) stop 'ERROR: please set correct Numoccupied, it should be small than num_wann'
   
   ik =0
   do i= 1, Nk1
      do j= 1, Nk2
         ik =ik +1
         kxy(:, ik)= K3D_start+ K3D_vec1*(i-1)/dble(Nk1-1)+ K3D_vec2*(j-1)/dble(Nk2-1) &
            -(K3D_vec1+K3D_vec2)/2d0
         kxy_shape(:, ik)= kxy(1, ik)* Origin_cell%Kua+ kxy(2, ik)* Origin_cell%Kub+ kxy(3, ik)* Origin_cell%Kuc 
         call rotate_k3_to_kplane(kxy_shape(:, ik), kxy_plane(:, ik))
      enddo
   enddo

   time_start= 0d0
   time_end= 0d0
   do ik= 1+cpuid, knv3, num_cpu
      if (cpuid==0.and. mod(ik/num_cpu, 500)==0) &
         write(stdout, '(a, i12, a, i12, a, f10.2, a)') &
         'ek_bulk_plane, ik ', ik, ' knv3',knv3, ' time left', &
         (knv3-ik)*(time_end- time_start)/num_cpu, ' s'
      call now(time_start)

      k = kxy(:, ik)

      ! generate bulk Hamiltonian
      Hamk_bulk= 0d0
      if (index(KPorTB, 'KP')/=0)then
         call ham_bulk_kp_abcb_graphene(k, Hamk_bulk)
      else
         !> deal with phonon system
         if (index(Particle,'phonon')/=0.and.LOTO_correction) then
            call ham_bulk_LOTO(k, Hamk_bulk)
         else
            call ham_bulk_latticegauge(k, Hamk_bulk)
           !call ham_bulk_atomicgauge(k, Hamk_bulk)
         endif
      endif

      !> diagonalization by call zheev in lapack
      W= 0d0
      call eigensystem_c( 'N', 'U', Num_wann ,Hamk_bulk, W)
      eigv(:, ik)= W(nband_min:nband_max)

      call now(time_end)
   enddo ! ik

#if defined (MPI)
   call mpi_allreduce(eigv,eigv_mpi,size(eigv),&
                     mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
   eigv_mpi= eigv
#endif

   if (index(Particle,'phonon')/=0) then
      eigv_mpi = eigv_mpi - MINVAL(eigv_mpi)
   endif

   !> deal with phonon system
   if (index(Particle,'phonon')/=0) then
      do ik=1, knv3
         do j=1, Num_wann
            eigv_mpi(j, ik)= sqrt(abs(eigv_mpi(j, ik)))*sign(1d0, eigv_mpi(j, ik))
            !eigv_mpi(j, ik)=     (abs(eigv_mpi(j, ik)))*sign(1d0, eigv_mpi(j, ik))
         enddo
      enddo
   endif
   eigv_mpi= eigv_mpi/eV2Hartree


   !> write out the data in gnuplot format
   outfileindex= outfileindex+ 1
   if (cpuid==0)then
      open(unit=outfileindex, file='bulkek_plane.dat')
      write(outfileindex, '(1000a19)')'# kx', 'ky', 'kz', 'k1', 'k2', 'k3', &
         'E(Numoccupied-1)', 'E(Numoccupied)' , 'E(Numoccupied+1)', 'E(Numoccupied+2)'
      do ik=1, knv3
         write(outfileindex, '(1000f19.9)')kxy_shape(:, ik)*Angstrom2atomic, &
            kxy_plane(:, ik)*Angstrom2atomic, eigv_mpi(:, ik)  
         if (mod(ik, nk2)==0) write(outfileindex, *)' '
      enddo
      close(outfileindex)
   endif

   !> write out the data in gnuplot format
   outfileindex= outfileindex+ 1
   if (cpuid==0)then
      open(unit=outfileindex, file='bulkek_plane-matlab.dat')
      write(outfileindex, '(1000a19)')'% kx', 'ky', 'kz', 'k1', 'k2', &
         'E(Numoccupied-1)', 'E(Numoccupied)' , 'E(Numoccupied+1)', 'E(Numoccupied+2)'
      do ik=1, knv3
         write(outfileindex, '(1000f19.9)')kxy_shape(:, ik)*Angstrom2atomic, &
            kxy_plane(:, ik)*Angstrom2atomic, eigv_mpi(:, ik)  
      enddo
      close(outfileindex)
   endif

   !> write out a script that can be used for gnuplot
   outfileindex= outfileindex+ 1
   if (cpuid==0)then
      open(unit=outfileindex, file='bulkek_plane.gnu')
      write(outfileindex, '(a)')"set encoding iso_8859_1"
      write(outfileindex, '(a)')'#set terminal  postscript enhanced color'
      write(outfileindex, '(a)')"#set output 'bulkek_plane.eps'"
      write(outfileindex, '(3a)')'set terminal  png truecolor enhanced', &
         ' size 1920, 1680 font ",36"'
      write(outfileindex, '(a)')"set output 'bulkek_plane.png'"
      write(outfileindex, '(a)')'set palette rgbformulae 33,13,10'
      write(outfileindex, '(a)')'unset key'
      write(outfileindex, '(a)')'set pm3d'
      write(outfileindex, '(a)')'set origin 0.2, 0'
      write(outfileindex, '(a)')'set size 0.8, 1'
      write(outfileindex, '(a)')'set border lw 3'
      write(outfileindex, '(a)')'#set xtics font ",24"'
      write(outfileindex, '(a)')'#set ytics font ",24"'
      write(outfileindex, '(a)')'set size ratio -1'
      write(outfileindex, '(a)')'set ticslevel 0'
      write(outfileindex, '(a)')'unset xtics'
      write(outfileindex, '(a)')'unset ytics'
      write(outfileindex, '(a)')'set view 80,60'
      write(outfileindex, '(a)')'set xlabel "k_1"'
      write(outfileindex, '(a)')'set ylabel "k_2"'
      write(outfileindex, '(a)')'set zlabel "Energy (eV)" rotate by 90'
      write(outfileindex, '(a)')'unset colorbox'
      write(outfileindex, '(a)')'set autoscale fix'
      write(outfileindex, '(a)')'set pm3d interpolate 4,4'
      write(outfileindex, '(2a)')"splot 'bulkek_plane.dat' u 4:5:8 w pm3d, \"
      write(outfileindex, '(2a)')"      'bulkek_plane.dat' u 4:5:9 w pm3d"

      close(outfileindex)

   endif ! cpuid

   deallocate(W)
   deallocate(Hamk_bulk)
   deallocate( kxy)
   deallocate( kxy_shape)
   deallocate( kxy_plane)
   deallocate( eigv    )
   deallocate( eigv_mpi)

   return
end subroutine ek_bulk_plane

subroutine ek_bulk_cube
   ! Calculate bulk's energy band using wannier TB method for a given k cube
   ! Revised by QS.Wu on Feb 18 2019 @EPFL, Switzerland
   ! Under the licence GPL V3

   use wmpi
   use para

   implicit none

   integer :: ik, i, j, l, ib, knv3, ierr, nwann
   integer :: nband_min, nband_max, nband_store

   real(dp) :: time_start, time_end

   real(Dp) :: k(3)
   real(Dp), allocatable :: W(:)
   real(dp), allocatable :: kabc(:,:)
   real(dp), allocatable :: kxyz(:,:)
   
   ! Hamiltonian of bulk system
   complex(Dp), allocatable :: Hamk_bulk(:, :)

   ! eigen value of H
   real(dp), allocatable :: gap(:)
   real(dp), allocatable :: gap_mpi(:)
   real(dp), allocatable :: eigv(:,:)
   real(dp), allocatable :: eigv_mpi(:,:)

   nband_min= Numoccupied- 1
   nband_max= Numoccupied+ 2
   if (nband_min< 1) nband_min= 1
   if (nband_max> Num_wann) nband_max= Num_wann
   nband_store= nband_max- nband_min+ 1

   knv3= Nk1*Nk2*Nk3
   allocate(W(Num_wann))
   allocate(Hamk_bulk(Num_wann, Num_wann))
   allocate( kabc(3, knv3), kxyz(3, knv3))
   allocate( eigv    (nband_store, knv3))
   allocate( eigv_mpi(nband_store, knv3))

   eigv    = 0d0
   eigv_mpi= 0d0
   kabc= 0d0; kxyz= 0d0

   if (Nk1<2 .or. Nk2<2 .or. Nk3<2)  &
      stop 'ERROR: I refuse to do this job because you give me so small Nk1, Nk2 and Nk3'      
   if (Numoccupied> Num_wann) &
      stop 'ERROR: please set correct Numoccupied, it should be small than num_wann'
  
   !> The center of the cube is K3D_start_cube.
   ik =0
   do i= 1, Nk1
      do j= 1, Nk2
         do l= 1, Nk3
            ik =ik +1
            kabc(:, ik)= K3D_start_cube+ K3D_vec1_cube*(i-1)/dble(nk1-1)  &
                         + K3D_vec2_cube*(j-1)/dble(nk2-1)  &
                         + K3D_vec3_cube*(l-1)/dble(nk3-1) &
                         - (K3D_vec1_cube+K3D_vec2_cube+K3D_vec3_cube)/2d0
            kxyz(:, ik)= kabc(1, ik)* Origin_cell%Kua+ kabc(2, ik)* Origin_cell%Kub+ kabc(3, ik)* Origin_cell%Kuc 
         enddo
      enddo
   enddo

   time_start= 0d0
   time_end= 0d0
   do ik= 1+cpuid, knv3, num_cpu
      if (cpuid==0.and. mod(ik/num_cpu, 500)==0) &
         write(stdout, '(a, i12, a, i12, a, f10.2, a)') &
         'ek_bulk_cube, ik ', ik, ' knv3',knv3, ' time left', &
         (knv3-ik)*(time_end- time_start)/num_cpu, ' s'
      call now(time_start)

      k = kabc(:, ik)

      ! calculation bulk hamiltonian
      Hamk_bulk= 0d0
      call ham_bulk_atomicgauge(k, Hamk_bulk)
      Hamk_bulk= Hamk_bulk/eV2Hartree

      !> diagonalization by call zheev in lapack
      W= 0d0
      call eigensystem_c( 'N', 'U', Num_wann ,Hamk_bulk, W)
      eigv(:, ik)= W(nband_min:nband_max)

      call now(time_end)
   enddo ! ik

#if defined (MPI)
     call mpi_allreduce(eigv,eigv_mpi,size(eigv),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
     eigv_mpi= eigv
#endif

     !> write out the data in gnuplot format
     outfileindex= outfileindex+ 1
     if (cpuid==0)then
        open(unit=outfileindex, file='bulkek_cube.txt')
        write(outfileindex, '(3(a6, i6))')'# Nk1:', Nk1, ' Nk2:', Nk2, ' Nk3:', Nk3 
        write(outfileindex, '(1000a19)')'# kx', 'ky', 'kz', 'k1', 'k2', 'k3', &
           'E(Numoccupied-1)', 'E(Numoccupied)' , 'E(Numoccupied+1)', 'E(Numoccupied+2)'
        do ik=1, knv3
           write(outfileindex, '(1000f19.9)')kxyz(:, ik)*Angstrom2atomic, &
              kabc(:, ik), eigv_mpi(:, ik)  
        enddo
        close(outfileindex)
     endif

     deallocate(W)
     deallocate(Hamk_bulk)
     deallocate( kabc, kxyz)
     deallocate( eigv    )
     deallocate( eigv_mpi)

   return
end subroutine ek_bulk_cube



subroutine sparse_ekbulk
   use sparse
   use para
   implicit none


   !> some temporary integers
   integer :: ik, i, j, ierr, ib, ig, numk

   ! wave vector
   real(dp) :: k3(3)

   !> dim= Num_wann, knv3
   real(dp), allocatable :: W(:)
   real(dp), allocatable :: eigv(:, :)
   real(dp), allocatable :: eigv_mpi(:, :)

   real(dp) :: emin, emax

   !> dim= Num_wann*Num_wann
   integer :: nnzmax, nnz
   complex(dp), allocatable :: acoo(:)
   integer, allocatable :: jcoo(:)
   integer, allocatable :: icoo(:)

   !> storage for overlap matrix
   integer :: snnzmax, snnz
   complex(dp), allocatable :: sacoo_k(:)
   integer, allocatable :: sjcoo_k(:)
   integer, allocatable :: sicoo_k(:)


   !> eigenvector of the sparse matrix acoo. Dim=(Num_wann, neval)
   complex(dp), allocatable :: psi(:)
   complex(dp), allocatable :: psi_project(:)
   complex(dp), allocatable :: zeigv(:, :)

   !> print the weight for the Selected_WannierOrbitals
   real(dp), allocatable :: weight(:, :)
   real(dp), allocatable :: dos_selected(:, :, :, :)
   real(dp), allocatable :: dos_selected_mpi(:, :, :, :)

   !number of ARPACK eigenvalues
   integer :: neval

   ! number of Arnoldi vectors
   integer :: nvecs

   !> calculate eigenvector or not
   logical :: ritzvec

   !shift-invert sigma
   complex(dp) :: sigma=(0d0,0d0)

   !> time measurement
   real(dp) :: time1, time2, time3


   if (OmegaNum==0) OmegaNum= Num_wann
   if (NumSelectedEigenVals==0) NumSelectedEigenVals=OmegaNum

   !> first use NumSelectedEigenVals, if NumSelectedEigenVals is not set, 
   !> then use OmegaNum; if OmegaNum is also not set, 
   !> then use Num_wann
   if (NumSelectedEigenVals>0) then
      neval= NumSelectedEigenVals
   else if (OmegaNum>0) then
      neval= OmegaNum
   else
      neval = Num_wann
   endif

   if (neval>Num_wann-2) neval= Num_wann- 2

   !> ncv
   if (NumLCZVecs.ne.0) then
      nvecs= NumLCZVecs
   else
      nvecs=int(2*neval)
   endif

   if (neval+2>=nvecs) neval= nvecs-2

   if (nvecs<20) nvecs= 20
   if (nvecs>Num_wann) nvecs= Num_wann

   if (trim(adjustl(projection_weight_mode))=='FOLDEDKPOITNS') then
      numk= Nk3_unfold_point_mode
   else
      numk=1
   endif

   nnzmax=splen+Num_wann
   if (.not.Orthogonal_Basis) then
      nnzmax=splen+Num_wann+splen_overlap_input
   else
      nnzmax=splen+Num_wann
   endif
   nnz=splen
   snnzmax=splen_overlap_input
   allocate( acoo(nnzmax))
   allocate( jcoo(nnzmax))
   allocate( icoo(nnzmax))
   if (.not.Orthogonal_Basis) then
      allocate( sacoo_k(snnzmax))
      allocate( sjcoo_k(snnzmax))
      allocate( sicoo_k(snnzmax))
   endif
   allocate( W( neval))
   allocate( eigv( neval, nk3_band))
   allocate( eigv_mpi( neval, nk3_band))
   allocate( psi(Num_wann))
   allocate( psi_project(Num_wann))
   allocate( zeigv(Num_wann,nvecs))
   allocate( weight(NumberofSelectedOrbitals_groups, numk))
   allocate( dos_selected     (neval,   nk3_band, NumberofSelectedOrbitals_groups, numk))
   allocate( dos_selected_mpi (neval,   nk3_band, NumberofSelectedOrbitals_groups, numk))
   psi=0d0; psi_project= 0d0; zeigv= 0d0
   dos_selected= 0d0; dos_selected_mpi= 0d0

   eigv_mpi= 0d0;  eigv    = 0d0
   acoo= 0d0; icoo=0; jcoo=0
   sacoo_k= 0d0; sicoo_k=0; sjcoo_k=0

   ritzvec= BulkFatBand_calc

   if (cpuid.eq.0)then
      write(stdout, *)'# parameters for ARPACK in ek_bulk'
      write(stdout, *)'NEV=', neval
      write(stdout, *)'NCV=', nvecs
      write(stdout, *)'Dimension=', Num_wann
   endif

   !> calculate the energy bands along special k line
   k3= 0

   !> change the energy unit from Hatree to eV
   sigma=(1d0,0d0)*iso_energy/eV2Hartree

   do ik=1+ cpuid, nk3_band, num_cpu
      if (cpuid==0) write(stdout, '(a, 2i10)') 'BulkBand_calc in sparse mode:', ik,nk3_band
      k3 = kpath_3d(:, ik)
      call now(time1)
      !> use atomic gauge here
      call ham_bulk_coo_sparsehr(k3,acoo,icoo,jcoo)
      if (.not.Orthogonal_Basis) then
         call overlap_bulk_coo_sparse(k3, sacoo_k, sicoo_k, sjcoo_k)
      endif
      !> change the energy unit from Hatree to eV
      acoo= acoo/eV2Hartree
      nnz= splen
      snnz=splen_overlap_input
      call now(time2)
      
      !> diagonalization by call zheev in lapack
      W= 0d0
      !> after arpack_sparse_coo_eigs, nnz will be updated.
      if (.not.Orthogonal_Basis) then
         !> non-orthogonal basis like, openmx
         call arpack_sparse_coo_eigs_nonorth(Num_wann, nnzmax, nnz, acoo, jcoo, icoo, &
             snnzmax, snnz, sacoo_k, sjcoo_k, sicoo_k, neval,nvecs,W,sigma,zeigv, ritzvec)
      else
         call arpack_sparse_coo_eigs(Num_wann,nnzmax,nnz,acoo,jcoo,icoo,neval,nvecs,W,sigma, zeigv, ritzvec)
      endif
      call now(time3)
      eigv(1:neval, ik)= W(1:neval)

      do ib= 1, neval
         psi(:)= zeigv(:, ib)  !> the eigenvector of ib'th band
         call get_projection_weight_bulk(k3, numk, psi, weight)
         dos_selected(ib, ik, :, :)= weight(:, :)
      enddo
      !> calculate the weight on the selected orbitals
     !do ig=1, NumberofSelectedOrbitals_groups
     !   do ib= 1, neval
     !      psi(:)= zeigv(:, ib)  !> the eigenvector of ib'th band
     !      do i= 1, NumberofSelectedOrbitals(ig)
     !         j= Selected_WannierOrbitals(ig)%iarray(i)
     !         dos_selected(ib, ik, ig)= dos_selected(ib, ik, ig)+ abs(psi(j))**2
     !      enddo ! sweep the selected orbitals
     !   enddo ! ib sweep the eigenvalue
     !enddo

      if (cpuid==0)write(stdout, '(a, f20.2, a)')'  >> Time cost for constructing H: ', time2-time1, ' s'
      if (cpuid==0)write(stdout, '(a, f20.2, a)')'  >> Time cost for diagonalize H: ', time3-time2, ' s'
   enddo !ik

#if defined (MPI)
   call mpi_allreduce(eigv,eigv_mpi,size(eigv),&
      mpi_dp,mpi_sum,mpi_cmw,ierr)
   call mpi_allreduce(dos_selected, dos_selected_mpi,size(dos_selected),&
      mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
   eigv_mpi= eigv
   dos_selected_mpi= dos_selected
#endif

   !> minimum and maximum value of energy bands
   emin= minval(eigv_mpi)+0.05*(maxval(eigv_mpi)-minval(eigv_mpi))
   emax= maxval(eigv_mpi)-0.05*(maxval(eigv_mpi)-minval(eigv_mpi))
   


   outfileindex= outfileindex+ 1
   if (cpuid==0)then
      open(unit=outfileindex, file='bulkek.dat')

      !> 
      if (BulkFatBand_calc) then
         write(outfileindex, &
            "('#', a12, a14, 1X, '| projection', 100(3X,'|group', i2, ': A '))")&
            'klen', 'E', (i, i=1, NumberofSelectedOrbitals_groups)
         write(outfileindex, "('#column', i5, 200i16)")(i, i=1, 2+NumberofSelectedOrbitals_groups)
      else
         write(outfileindex, '(2a19)') '# klen', 'E(eV)'
      endif

      do i=1, neval
         do ik=1, nk3_band
            if (BulkFatBand_calc) then
               write(outfileindex, '(300f16.9)')k3len(ik)*Angstrom2atomic,eigv_mpi(i, ik), &
                  (dos_selected_mpi(i, ik, ig, :), ig=1, NumberofSelectedOrbitals_groups)
            else
               write(outfileindex, '(2f16.9)')k3len(ik)*Angstrom2atomic,eigv_mpi(i, ik)
            endif
         enddo
         write(outfileindex, *)' '
      enddo
      close(outfileindex)
   endif

   call generate_ek_kpath_gnu('bulkek.dat', 'bulkek.gnu', 'bulkek.pdf', &
                                 emin, emax, nk3_band, Nk3lines, &
                                 k3line_name, k3line_stop, k3len)

#if defined (MPI)
   call mpi_barrier(mpi_cmw, ierr)
#endif

   deallocate( acoo)
   deallocate( jcoo)
   deallocate( icoo)
   deallocate( W)
   deallocate( eigv)
   deallocate( eigv_mpi)
   deallocate( zeigv)
   deallocate( dos_selected)
   deallocate( dos_selected_mpi)

   return
end subroutine sparse_ekbulk


!> only test for valley projection
subroutine sparse_ekbulk_valley
   use sparse
   use para
   implicit none


   !> some temporary integers
   integer :: ik, i, j, ierr, ib, ig
   integer :: ie, ND, NDMAX, ie1, ie2

   ! wave vector
   real(dp) :: k3(3)

   !> dim= Num_wann, knv3
   real(dp), allocatable :: W(:)
   real(dp), allocatable :: eigv(:, :)
   real(dp), allocatable :: eigv_mpi(:, :)

   real(dp) :: emin, emax

   complex(dp) :: alpha_mv, beta_mv

   !> dim= Num_wann*Num_wann
   integer :: nnzmax, nnz, nnzmax_valley, nnz_valley
   complex(dp), allocatable :: acoo(:), acoo_valley(:)
   integer, allocatable :: jcoo(:), jcoo_valley(:)
   integer, allocatable :: icoo(:), icoo_valley(:)

   !> storage for overlap matrix
   integer :: snnzmax, snnz
   complex(dp), allocatable :: sacoo_k(:)
   integer, allocatable :: sjcoo_k(:)
   integer, allocatable :: sicoo_k(:)


   !> eigenvector of the sparse matrix acoo. Dim=(Num_wann, neval)
   complex(dp), allocatable :: psi_project(:)
   complex(dp), allocatable :: zeigv(:, :)

   ! Hamiltonian of bulk system
   complex(Dp), allocatable :: psi(:), vpsi(:), valley_k_nd(:, :), psi1(:), psi2(:)
   complex(Dp), allocatable :: valley_k(:, :)
   complex(Dp), allocatable :: VL(:, :), VR(:, :), valley_eig(:)
   logical, allocatable :: valley_plus(:, :), valley_plus_mpi(:, :)


   !> print the weight for the Selected_WannierOrbitals
   real(dp), allocatable :: weight_valley(:, :), weight_valley_mpi(:, :)

   !number of ARPACK eigenvalues
   integer :: neval

   ! number of Arnoldi vectors
   integer :: nvecs

   !> calculate eigenvector or not
   logical :: ritzvec

   !shift-invert sigma
   complex(dp) :: sigma=(0d0,0d0)

   !> time measurement
   real(dp) :: time1, time2, time3

   !> 

   ! vector multiply vector v1*v2
   complex(dp), external :: zdotc
   real(dp) :: tolde= 1E-4


   !if (OmegaNum==0) OmegaNum= Num_wann
   !if (NumSelectedEigenVals==0) NumSelectedEigenVals=OmegaNum

   !> first use NumSelectedEigenVals, if NumSelectedEigenVals is not set, 
   !> then use OmegaNum; if OmegaNum is also not set, 
   !> then use Num_wann
   if (NumSelectedEigenVals>0) then
      neval= NumSelectedEigenVals
   else if (OmegaNum>0) then
      neval= OmegaNum
   else
      neval = Num_wann
   endif

   if (neval>Num_wann-2) neval= Num_wann- 2

   !> ncv
   nvecs=int(2*neval)

   if (nvecs<20) nvecs= 20
   if (nvecs>Num_wann) nvecs= Num_wann

   nnzmax=splen+Num_wann
   if (.not.Orthogonal_Basis) then
      nnzmax=splen+Num_wann+splen_overlap_input
   else
      nnzmax=splen+Num_wann
   endif
   nnz=splen
   snnzmax=splen_overlap_input
 
   NDMAX= 12
   allocate( acoo(nnzmax))
   allocate( jcoo(nnzmax))
   allocate( icoo(nnzmax))
   if (.not.Orthogonal_Basis) then
      allocate( sacoo_k(snnzmax))
      allocate( sjcoo_k(snnzmax))
      allocate( sicoo_k(snnzmax))
   endif
 
   nnzmax_valley = splen_valley_input
   allocate( acoo_valley(nnzmax_valley))
   allocate( jcoo_valley(nnzmax_valley))
   allocate( icoo_valley(nnzmax_valley))

   allocate( W( neval))
   allocate( eigv( neval, nk3_band))
   allocate( eigv_mpi( neval, nk3_band))
   allocate( psi_project(Num_wann))
   allocate( zeigv(Num_wann,nvecs))
   allocate( weight_valley(neval, nk3_band), weight_valley_mpi(neval, nk3_band))
   allocate( psi(Num_wann), vpsi(Num_wann))
   allocate( psi1(Num_wann), psi2(Num_wann))
   allocate(valley_k(Num_wann, Num_wann))
   allocate(valley_plus(Num_wann, nk3_band))
   allocate(valley_plus_mpi(Num_wann, nk3_band))
   allocate( valley_k_nd(NDMAX, NDMAX), valley_eig(NDMAX))
   allocate( VL(NDMAX, NDMAX), VR(NDMAX, NDMAX))
   valley_k_nd= 0d0
   valley_plus= .False.; valley_plus_mpi= .False.
   psi=0d0; vpsi=0d0; psi_project= 0d0; zeigv= 0d0
   weight_valley_mpi=0d0; weight_valley=0d0

   eigv_mpi= 0d0;  eigv    = 0d0
   acoo= 0d0; icoo=0; jcoo=0

   ritzvec= .True.

   !> change the energy unit from Hatree to eV
   sigma=(1d0,0d0)*iso_energy/eV2Hartree

   !> calculate the energy bands along special k line
   k3= 0
   do ik=1+ cpuid, nk3_band, num_cpu
      if (cpuid==0) write(stdout, '(a, 2i10)') 'BulkBand_calc in sparse mode:', ik,nk3_band
      k3 = kpath_3d(:, ik)
      call now(time1)
      call ham_bulk_coo_sparsehr(k3,acoo,icoo,jcoo)
      !> change the energy unit from Hatree to eV
      acoo= acoo/eV2Hartree
      nnz= splen
      if (.not.Orthogonal_Basis) then
         call overlap_bulk_coo_sparse(k3, sacoo_k, sicoo_k, sjcoo_k)
      endif
      call now(time2)
 
      !> diagonalization by call zheev in lapack
      W= 0d0
      !> after arpack_sparse_coo_eigs, nnz will be updated.
      if (.not.Orthogonal_Basis) then
         !> non-orthogonal basis like, openmx
         call arpack_sparse_coo_eigs_nonorth(Num_wann, nnzmax, nnz, acoo, jcoo, icoo, &
             snnzmax, snnz, sacoo_k, sjcoo_k, sicoo_k, neval,nvecs,W,sigma,zeigv, ritzvec)
      else
         call arpack_sparse_coo_eigs(Num_wann,nnzmax,nnz,acoo,jcoo,icoo,neval,nvecs,W,sigma, zeigv, ritzvec)
      endif
      call now(time3)
      eigv(1:neval, ik)= W(1:neval)


      !> get valley operator in coo format
      call valley_k_coo_sparsehr(nnzmax_valley, k3, acoo_valley, icoo_valley, jcoo_valley)

      !> get the valley projection
!     do ib= 1, neval
!        psi(:)= zeigv(:, ib)  !> the eigenvector of ib'th band

!        !> weight_valley= <psi|vz|psi>
!        call mkl_zcoogemv('N', Num_wann, acoo_valley, icoo_valley, jcoo_valley, nnzmax_valley, psi, vpsi)

!        weight_valley(ib, ik)= real(zdotc(Num_wann, psi, 1, vpsi, 1))

!     enddo

      !> check the degeneracy of each band
      IE=1
      do while (ie.le.neval)
         ND=1
         if (ie+ND.le.neval) then
            do while ((ie+ND).le.neval .and. (W(ie+ND)- W(ie)).lt.tolde)
               if (ie+ND .ge. neval) exit
               ND= ND+1
            enddo
         endif


         !> if the degeneracy is larger than 1, we need to calculate the matrix
         valley_k_nd= 0d0
         do ie1= 1, ND
            psi1= zeigv(:, IE+ie1-1)

            do ie2= 1, ND
               psi2= zeigv(:, IE+ie2-1)
               vpsi=0d0
               !call mkl_zcoogemv('N', Num_wann, acoo_valley, icoo_valley, jcoo_valley, nnzmax_valley, psi2, vpsi)
               call coomv_z(Num_wann, nnzmax_valley, acoo_valley, icoo_valley, jcoo_valley, psi2, vpsi)
               valley_k_nd(ie1, ie2)= zdotc(Num_wann, psi1, 1, vpsi, 1)
            enddo
           !if (abs(W(IE+ie1-1)/eV2Hartree)<2d0) &
           !write( *, '(a, 2i5, 200f8.3)') "ND, ie1: ", ND, ie1, valley_k_nd(ie1, 1:ND)
         enddo

         !> diagonalize the matrix
         VL = 0d0; VR= 0d0
         if (ND>1) then
            call zgeev_sys(ND, valley_k_nd(1:ND, 1:ND), valley_eig(1:ND),'N',VL(1:ND, 1:ND),"V",VR(1:ND, 1:ND) )
            do ie1= 1, ND
               weight_valley(IE+ie1-1, ik) = real(valley_eig(ie1))
            enddo
         else
            valley_eig(1)= valley_k_nd(1, 1)
            weight_valley(IE, ik) = real(valley_k_nd(1, 1))
         endif
         !  if (abs(W(IE+ie1-1)/eV2Hartree)<2d0) &
         ! write(*, '(a, 2i5, a, 200f8.3)'), 'ND, ie', ND, ie, '  valley_eig', real(valley_eig(1:ND))
         !  if (abs(W(IE+ie1-1)/eV2Hartree)<2d0) &
         ! print *, ' '

         !> get new eigenvector
        !do ie1= 1, ND
        !   psi= 0d0
        !   do ie2= 1, ND
        !      psi= psi+ VR(ie2,  ie1)* zeigv(:, IE+ie2-1)
        !   enddo
        !enddo

         !> if the degeneracy is larger than 1, we need to calculate the matrix
         IE= IE+ ND
      enddo
      do ie=1, neval
         if (weight_valley(ie, ik)>0) valley_plus(ie, ik)=.true.
      enddo


      if (cpuid==0)write(stdout, '(a, f20.2, a)')'  >> Time cost for constructing H: ', time2-time1, ' s'
      if (cpuid==0)write(stdout, '(a, f20.2, a)')'  >> Time cost for diagonalize H: ', time3-time2, ' s'
   enddo !ik

#if defined (MPI)
   call mpi_allreduce(eigv,eigv_mpi,size(eigv),&
      mpi_dp,mpi_sum,mpi_cmw,ierr)
   call mpi_allreduce(weight_valley, weight_valley_mpi,size(weight_valley),&
      mpi_dp,mpi_sum,mpi_cmw,ierr)
   call mpi_allreduce(valley_plus, valley_plus_mpi,size(valley_plus),&
      mpi_logical,mpi_lor,mpi_cmw,ierr)
#else
   eigv_mpi= eigv
   weight_valley_mpi= weight_valley
   valley_plus_mpi= valley_plus
#endif

   !> minimum and maximum value of energy bands
   emin= minval(eigv_mpi)+0.05*(maxval(eigv_mpi)-minval(eigv_mpi))
   emax= maxval(eigv_mpi)-0.05*(maxval(eigv_mpi)-minval(eigv_mpi))
   


   outfileindex= outfileindex+ 1
   if (cpuid==0)then
      open(unit=outfileindex, file='bulkek.dat')
      do i=1, neval
         do ik=1, nk3_band
            write(outfileindex, '(300f16.9)')k3len(ik)*Angstrom2atomic,eigv_mpi(i, ik), &
               weight_valley_mpi(i, ik)
         enddo
         write(outfileindex, *)' '
      enddo
      close(outfileindex)
   endif

   outfileindex= outfileindex+ 1
   if (cpuid==0)then
      open(unit=outfileindex, file='bulkek_valley_plus.dat')
      do i=1, neval
         do ik=1, nk3_band
            if (valley_plus_mpi(i, ik)) then
               write(outfileindex, '(200F16.8)')k3len(ik)*Angstrom2atomic,eigv_mpi(i, ik), &
                  weight_valley_mpi(i, ik)
            endif
         enddo
         write(outfileindex, *)' '
      enddo
      close(outfileindex)
   endif

   outfileindex= outfileindex+ 1
   if (cpuid==0)then
      open(unit=outfileindex, file='bulkek_valley_minus.dat')
      do i=1, neval
         do ik=1, nk3_band
            if (.not.valley_plus_mpi(i, ik)) then
               write(outfileindex, '(200F16.8)')k3len(ik)*Angstrom2atomic,eigv_mpi(i, ik), &
                  weight_valley_mpi(i, ik)
            endif
         enddo
         write(outfileindex, *)' '
      enddo
      close(outfileindex)
   endif


   call generate_ek_kpath_gnu('bulkek.dat', 'bulkek.gnu', 'bulkek.pdf', &
                                 emin, emax, nk3_band, Nk3lines, &
                                 k3line_name, k3line_stop, k3len)

#if defined (MPI)
   call mpi_barrier(mpi_cmw, ierr)
#endif

   deallocate( acoo)
   deallocate( jcoo)
   deallocate( icoo)
   deallocate( W)
   deallocate( eigv)
   deallocate( eigv_mpi)
   deallocate( zeigv)
   deallocate( weight_valley, weight_valley_mpi)

   return
end subroutine sparse_ekbulk_valley


subroutine sparse_ekbulk_plane
   use sparse
   use para
   implicit none


   !> some temporary integers
   integer :: ik, i, j, ierr, ib, ig, knv3

   ! wave vector
   real(dp) :: k3(3)

   !> dim= Num_wann, knv3
   real(dp), allocatable :: W(:)
   real(dp), allocatable :: eigv(:, :)
   real(dp), allocatable :: eigv_mpi(:, :)

   real(dp), allocatable :: kxy(:,:)
   real(dp), allocatable :: kxy_shape(:,:)
   real(dp), allocatable :: kxy_plane(:,:)
   real(dp) :: emin, emax

   !> dim= Num_wann*Num_wann
   integer :: nnzmax, nnz
   complex(dp), allocatable :: acoo(:)
   integer, allocatable :: jcoo(:)
   integer, allocatable :: icoo(:)

   !> eigenvector of the sparse matrix acoo. Dim=(Num_wann, neval)
   complex(dp), allocatable :: psi(:)
   complex(dp), allocatable :: psi_project(:)
   complex(dp), allocatable :: zeigv(:, :)

   !> print the weight for the Selected_WannierOrbitals
   real(dp), allocatable :: weight(:, :)
   real(dp), allocatable :: dos_selected(:, :, :, :)
   real(dp), allocatable :: dos_selected_mpi(:, :, :, :)

   !number of ARPACK eigenvalues
   integer :: neval

   ! number of Arnoldi vectors
   integer :: nvecs

   !> calculate eigenvector or not
   logical :: ritzvec

   !shift-invert sigma
   complex(dp) :: sigma=(0d0,0d0)

   !> time measurement
   real(dp) :: time1, time2, time3, time_start, time_end


   !> first use NumSelectedEigenVals, if NumSelectedEigenVals is not set, 
   !> then use OmegaNum; if OmegaNum is also not set, 
   !> then use Num_wann
   if (NumSelectedEigenVals>0) then
      neval= NumSelectedEigenVals
   else if (OmegaNum>0) then
      neval= OmegaNum
   else
      neval = Num_wann
   endif

   if (neval>Num_wann-2) neval= Num_wann- 2

   !> ncv
   nvecs=int(2*neval)

   if (nvecs<20) nvecs= 20
   if (nvecs>Num_wann) nvecs= Num_wann

   knv3= nk1*Nk2
   allocate( kxy(3, knv3))
   allocate( kxy_shape(3, knv3))
   allocate( kxy_plane(3, knv3))

   eigv    = 0d0
   eigv_mpi= 0d0

   kxy=0d0
   kxy_shape=0d0
   kxy_plane=0d0
  
   if (Nk1<2 .or. Nk2<2) stop 'ERROR: I refuse to do this job because you give me so small Nk1 and Nk2'      
   if (Numoccupied> Num_wann) stop 'ERROR: please set correct Numoccupied, it should be small than num_wann'
   
   ik =0
   do i= 1, Nk1
      do j= 1, Nk2
         ik =ik +1
         kxy(:, ik)= K3D_start+ K3D_vec1*(i-1)/dble(Nk1-1)+ K3D_vec2*(j-1)/dble(Nk2-1) &
            -(K3D_vec1+K3D_vec2)/2d0
         kxy_shape(:, ik)= kxy(1, ik)* Origin_cell%Kua+ kxy(2, ik)* Origin_cell%Kub+ kxy(3, ik)* Origin_cell%Kuc 
         call rotate_k3_to_kplane(kxy_shape(:, ik), kxy_plane(:, ik))
      enddo
   enddo

   sigma=(1d0,0d0)*iso_energy
   nnzmax=splen+Num_wann
   nnz=splen
   allocate( acoo(nnzmax))
   allocate( jcoo(nnzmax))
   allocate( icoo(nnzmax))
   allocate( W( neval))
   allocate( eigv( neval, knv3))
   allocate( eigv_mpi( neval, knv3))
   allocate( zeigv(Num_wann,nvecs))
   zeigv= 0d0

   eigv_mpi= 0d0;  eigv    = 0d0
   acoo= 0d0; icoo=0; jcoo=0

   ritzvec= .False.

   !> calculate the energy bands along special k line
   k3= 0
   time_start= 0d0
   time_end= 0d0
   do ik=1+ cpuid, knv3, num_cpu
      if (cpuid==0.and. mod(ik/num_cpu, 20 )==0) &
         write(stdout, '(a, i12, a, i12, a, f10.2, a)') &
         'BulkBand_plane_calc sparse, ik ', ik, ' knv3',knv3, ' time left', &
         (knv3-ik)*(time_end- time_start)/num_cpu, ' s'
      call now(time_start)

      k3 = kxy(:, ik)
      call ham_bulk_coo_sparsehr(k3,acoo,icoo,jcoo)
      acoo= acoo/eV2Hartree
      nnz= splen
      call now(time2)
      
      !> diagonalization by call zheev in lapack
      W= 0d0
      !> after arpack_sparse_coo_eigs, nnz will be updated.
      call arpack_sparse_coo_eigs(Num_wann,nnzmax,nnz,acoo,jcoo,icoo,neval,nvecs,W,sigma, zeigv, ritzvec)
      call now(time_end)
      eigv(1:neval, ik)= W(1:neval)

      if (cpuid==0)write(stdout, '(a, f20.2, a)')'  >> Time cost for constructing H: ', time2-time_start, ' s'
      if (cpuid==0)write(stdout, '(a, f20.2, a)')'  >> Time cost for diagonalize H: ', time_end-time2, ' s'
   enddo !ik

#if defined (MPI)
   call mpi_allreduce(eigv,eigv_mpi,size(eigv),&
      mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
   eigv_mpi= eigv
#endif

   !> minimum and maximum value of energy bands
   emin= minval(eigv_mpi)+0.05*(maxval(eigv_mpi)-minval(eigv_mpi))
   emax= maxval(eigv_mpi)-0.05*(maxval(eigv_mpi)-minval(eigv_mpi))
   

   !> write out the data in gnuplot format
   outfileindex= outfileindex+ 1
   if (cpuid==0)then
      open(unit=outfileindex, file='bulkek_plane.dat')
      write(outfileindex, '(1000a19)')'# kx', 'ky', 'kz', 'k1', 'k2', 'k3', &
         'E(Numoccupied-1)', 'E(Numoccupied)' , 'E(Numoccupied+1)', 'E(Numoccupied+2)'
      do ik=1, knv3
         write(outfileindex, '(1000f19.9)')kxy_shape(:, ik)*Angstrom2atomic, &
            kxy_plane(:, ik)*Angstrom2atomic, eigv_mpi(:, ik)  
         if (mod(ik, nk2)==0) write(outfileindex, *)' '
      enddo
      close(outfileindex)
   endif

   !> write out the data in gnuplot format
   outfileindex= outfileindex+ 1
   if (cpuid==0)then
      open(unit=outfileindex, file='bulkek_plane-matlab.dat')
      write(outfileindex, '(1000a19)')'% kx', 'ky', 'kz', 'k1', 'k2', &
         'E(Numoccupied-1)', 'E(Numoccupied)' , 'E(Numoccupied+1)', 'E(Numoccupied+2)'
      do ik=1, knv3
         write(outfileindex, '(1000f19.9)')kxy_shape(:, ik)*Angstrom2atomic, &
            kxy_plane(:, ik)*Angstrom2atomic, eigv_mpi(:, ik)  
      enddo
      close(outfileindex)
   endif

   !> write out a script that can be used for gnuplot
   outfileindex= outfileindex+ 1
   if (cpuid==0)then
      open(unit=outfileindex, file='bulkek_plane.gnu')
      write(outfileindex, '(a)')"set encoding iso_8859_1"
      write(outfileindex, '(a)')'#set terminal  postscript enhanced color'
      write(outfileindex, '(a)')"#set output 'bulkek_plane.eps'"
      write(outfileindex, '(3a)')'set terminal  png truecolor enhanced', &
         ' size 1920, 1680 font ",36"'
      write(outfileindex, '(a)')"set output 'bulkek_plane.png'"
      write(outfileindex, '(a)')'set palette rgbformulae 33,13,10'
      write(outfileindex, '(a)')'unset key'
      write(outfileindex, '(a)')'set pm3d'
      write(outfileindex, '(a)')'set origin 0.2, 0'
      write(outfileindex, '(a)')'set size 0.8, 1'
      write(outfileindex, '(a)')'set border lw 3'
      write(outfileindex, '(a)')'#set xtics font ",24"'
      write(outfileindex, '(a)')'#set ytics font ",24"'
      write(outfileindex, '(a)')'set size ratio -1'
      write(outfileindex, '(a)')'set xtics'
      write(outfileindex, '(a)')'set ytics'
      write(outfileindex, '(a)')'set view 80,60'
      write(outfileindex, '(a)')'set xlabel "k_1"'
      write(outfileindex, '(a)')'set ylabel "k_2"'
      write(outfileindex, '(a)')'set zlabel "Energy (eV)" rotate by 90'
      write(outfileindex, '(a)')'unset colorbox'
      write(outfileindex, '(a)')'set autoscale fix'
      write(outfileindex, '(a)')'set pm3d interpolate 4,4'
      write(outfileindex, '(2a)')"splot 'bulkek_plane.dat' u 4:5:8 w pm3d, \"
      write(outfileindex, '(2a)')"      'bulkek_plane.dat' u 4:5:9 w pm3d"

      close(outfileindex)

   endif ! cpuid



#if defined (MPI)
   call mpi_barrier(mpi_cmw, ierr)
#endif

   deallocate( acoo)
   deallocate( jcoo)
   deallocate( icoo)
   deallocate( W)
   deallocate( eigv)
   deallocate( eigv_mpi)
   deallocate( zeigv)

   return
end subroutine sparse_ekbulk_plane


#if defined (ELPA)
subroutine ekbulk_elpa
   use para
   use elpa1
   use elpa2
   implicit none

   !> magnetic supercell size, perpendicular to the magnetic field
   !> Ndimq= Nq* Num_wann

   integer :: ik, ierr,i,j

   ! wave vector
   real(dp) :: k3(3)

   !> dim= Ndimq, knv3
!~    real(dp), allocatable :: W(:)
   real(dp), allocatable :: eigv(:, :)
   real(dp), allocatable :: eigv_mpi(:, :)
   real(dp), allocatable :: eigv_mpi2(:, :)

   !> wave function



   !> dim= Ndimq*Ndimq

   real(dp),external :: PhasePdGauge

   integer :: np_cols,np_rows
   !ProcessorCore_
   integer :: my_blacs_ctxt
   integer ::  nprow, npcol, my_prow, my_pcol
   integer :: switch001
   integer :: info_mpi,nblk=1,na_rows,na_cols

   complex(dp),allocatable :: haml_block(:,:) ,z(:,:)

   integer, external :: numroc
   integer, external :: indxl2g,indxg2p,indxg2l

   integer :: sc_desc(9)
   integer :: mpi_comm_rows, mpi_comm_cols

   integer :: neval


   logical success


!~    allocate( ham_landau(Ndimq, Ndimq))
   allocate( eigv(Num_wann,nk3_band))
   allocate( eigv_mpi(Num_wann,nk3_band))
   allocate( eigv_mpi2(Num_wann,nk3_band))
!~    allocate( psi(Ndimq))




   eigv_mpi= 0d0
   eigv    = 0d0



   neval=OmegaNum






   do np_cols = NINT(SQRT(REAL(num_cpu))),2,-1
      if(mod(num_cpu,np_cols) == 0 ) exit
   enddo
   ! search for uniform distribution near SQRT

! at the end of the above loop, nprocs is always divisible by np_cols
   np_rows = num_cpu/np_cols

   ! Blacs preparation
   my_blacs_ctxt = mpi_cmw
   call BLACS_Gridinit( my_blacs_ctxt, 'C', np_rows, np_cols )
   call BLACS_Gridinfo( my_blacs_ctxt, nprow, npcol, my_prow, my_pcol )

   ! All ELPA routines need MPI communicators for communicating within
! rows or columns of processes, these are set in get_elpa_row_col_comms.

   switch001= get_elpa_row_col_comms(mpi_comm_world, my_prow, my_pcol, mpi_comm_rows, mpi_comm_cols)

!~    switch001 = get_elpa_communicators(mpi_cmw, my_prow, my_pcol, mpi_comm_rows, mpi_comm_cols)

   na_rows = numroc(Num_wann, nblk, my_prow, 0, np_rows)
   na_cols = numroc(Num_wann, nblk, my_pcol, 0, np_cols)

   call descinit( sc_desc, Num_wann, Num_wann, nblk, nblk, 0, 0, my_blacs_ctxt, na_rows, info_mpi )

   !< NO eigenvectors>

   allocate(haml_block(na_rows,na_cols),z(na_rows,na_cols))
   haml_block=0.0d0
   z=0.0d0
   outfileindex= outfileindex+ 1
   if(cpuid==0) then
      open(outfileindex,file='bulk.ek',status='replace',form='unformatted')
      write(outfileindex) nk3_band,neval
   end if
   do ik=1,nk3_band

      haml_block=0.0d0
      z=0.0d0
      k3  = kpath_3d(:, ik)
      call ham_bulk_elpa( k3, haml_block,&
         np_rows,np_cols,&
         na_rows,na_cols,&
         nblk,&
         my_prow,my_pcol)

      call mpi_barrier(mpi_cmw, ierr)
!~       do i=1,na_rows
!~          do j=1,na_cols
!~             write(233,*) haml_block(i,j)
!~          end do
!~       end do
      success = solve_evp_complex_2stage&
         (Num_wann, neval,&
         haml_block, na_rows,eigv_mpi(:, ik),&
         z, na_rows,&
         nblk,na_cols,&
         mpi_comm_rows, mpi_comm_cols, mpi_cmw)


      if (.not.(success)) then
         if (cpuid==0) write(*,*) "solve_evp_complex_2stage produced an error! Aborting..."
         call MPI_ABORT(mpi_comm_world, ierr)
      endif
      if (cpuid==0) write (*,*) "#Diagonalization done"



   enddo !ib
   !,status='replace',form='unformatted')

   if(cpuid==0) then
      do ik=1,nk3_band
         write(outfileindex) k3len(ik)*Angstrom2atomic,eigv_mpi(1:neval,ik)
!~     write(239,*) k3len(ik),eigv_mpi(1:neval,ik)
      enddo
   endif



end subroutine ekbulk_elpa
#endif




   subroutine orbitaltexture3D
     ! Calculate bulk's energy band and orbital texture using wannier TB method for a given k plane
     ! Copyright (c) 2018 QuanSheng Wu. All rights reserved.
     ! at EPFL Feb 4th 2018
     ! wuquansheng@gmail.com

     use wmpi
     use para

     implicit none

     integer :: ik, i, j, ib, ig
     integer :: nkx, nky
     integer :: knv3
     integer :: ierr
     character(40) :: filename
     
     real(Dp) :: k(3)
     real(Dp), allocatable :: W(:)
     real(dp), allocatable :: kxy(:,:)
     real(dp), allocatable :: kxy_shape(:,:)
     real(dp), allocatable :: kxy_plane(:,:)
     
     ! Hamiltonian of bulk system
     complex(Dp), allocatable :: Hamk_bulk(:, :)

     ! eigen value of H
	  real(dp), allocatable :: eigv(:,:)
	  real(dp), allocatable :: eigv_mpi(:,:)
     real(dp), allocatable :: dos_selected(:, :, :)
     real(dp), allocatable :: dos_selected_mpi(:, :, :)

     !> wave function
     complex(dp), allocatable :: psi(:)

     nkx= Nk1
     nky= Nk2
     knv3= nkx*nky
     allocate( W( Num_wann))
     allocate( psi( Num_wann))
     allocate( Hamk_bulk( Num_wann, Num_wann))
     allocate( kxy(3, knv3))
     allocate( kxy_shape(3, knv3))
     allocate( kxy_plane(3, knv3))
     allocate( dos_selected     (NumberofSelectedBands,   knv3, NumberofSelectedOrbitals_groups))
     allocate( dos_selected_mpi (NumberofSelectedBands,   knv3, NumberofSelectedOrbitals_groups))
     allocate( eigv    (NumberofSelectedBands, knv3))
     allocate( eigv_mpi(NumberofSelectedBands, knv3))
     W= 0d0
     psi= 0d0
     Hamk_bulk= 0d0
     kxy = 0d0
     kxy_shape= 0d0
     eigv    = 0d0
     eigv_mpi= 0d0
     dos_selected= 0d0
     dos_selected_mpi= 0d0
 
     !> set up k slice
     ik =0
     do i= 1, nkx
        do j= 1, nky
           ik =ik +1
           kxy(:, ik)= K3D_start+ K3D_vec1*(i-1)/dble(Nk1-1)+ K3D_vec2*(j-1)/dble(Nk2-1) &
              -(K3D_vec1+K3D_vec2)/2d0
           kxy_shape(:, ik)= kxy(1, ik)* Origin_cell%Kua+ kxy(2, ik)* Origin_cell%Kub+ kxy(3, ik)* Origin_cell%Kuc 
           call rotate_k3_to_kplane(kxy_shape(:, ik), kxy_plane(:, ik))
        enddo
     enddo

     do ik= 1+cpuid, knv3, num_cpu
        if (cpuid==0) write(stdout, '(a, 2i10)')' In orbitaltexture3D ', ik, knv3

        k = kxy(:, ik)

        ! calculation bulk hamiltonian
        Hamk_bulk= 0d0

        call ham_bulk_atomicgauge(k, Hamk_bulk)
        Hamk_bulk= Hamk_bulk/eV2Hartree

        !> diagonalization by call zheev in lapack
        W= 0d0
        call eigensystem_c( 'V', 'U', Num_wann ,Hamk_bulk, W )

        do ig=1, NumberofSelectedOrbitals_groups
           do ib= 1, NumberofSelectedBands
              psi(:)= Hamk_bulk(:, Selected_band_index(ib))
              do i=1, NumberofSelectedOrbitals(ig)
                 dos_selected(ib, ik, ig)= dos_selected(ib, ik, ig)+ &
                    abs(psi(Selected_WannierOrbitals(ig)%iarray(i)))**2
              enddo ! i
           enddo ! ib
        enddo ! ig

        !> for 44 bands
        do ib=1, NumberofSelectedBands
           eigv(ib, ik)= W(Selected_band_index(ib))
        enddo ! i
     enddo ! ik

#if defined (MPI)
     call mpi_allreduce(eigv,eigv_mpi,size(eigv),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
     call mpi_allreduce(dos_selected ,dos_selected_mpi,size(dos_selected ),&
        mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
     eigv_mpi= eigv
     dos_selected_mpi= dos_selected
#endif

     do ig=1, NumberofSelectedOrbitals_groups
        write(filename, '(a,i1)')'orbitaltexture3D.dat-group',ig
        outfileindex= outfileindex+ 1
        if (cpuid==0)then
           open(unit=outfileindex, file=filename)
           write(outfileindex, '(6a16, a)')'# kx', ' ky', ' kz', ' k1', ' k2', ' k3',  &
              ' (E(ib), dos(ib)), ib=1, NumberofSelectedBands'
           write(outfileindex, '(a, 2i10)')'# Nk1, Nk2=', Nk1, Nk2
           do ik=1, knv3
              write(outfileindex, '(2000f16.8)') kxy_shape(:, ik)*Angstrom2atomic, kxy_plane(:, ik)*Angstrom2atomic, &
                 (  eigv_mpi(ib, ik), dos_selected_mpi(ib, ik, ig) , ib=1, NumberofSelectedBands)
              if (mod(ik, nky)==0) write(outfileindex, *)' '
           enddo
           close(outfileindex)
        endif
     enddo !ig

#if defined (MPI)
     call mpi_barrier(mpi_cmw, ierr)
#endif

     outfileindex= outfileindex+ 1
     if (cpuid==0)then
        open(unit=outfileindex, file='orbitaltexture3D.dat-matlab')
        write(outfileindex, '(6a16, a)')'% here we only write out the first group, please modify code'
        write(outfileindex, '(6a16, a)')'% kx', ' ky', ' kz', ' k1', ' k2', ' k3', &
           ' (E(ib), dos(ib)), ib=1, NumberofSelectedOrbitals'
        write(outfileindex, '(a, 2i10)')'% Nk1, Nk2=', Nk1, Nk2
        do ik=1, knv3
           write(outfileindex, '(1000f16.8)')kxy_shape(:, ik)*Angstrom2atomic, kxy_plane(:, ik)*Angstrom2atomic, & 
              (eigv_mpi(ib, ik), dos_selected_mpi(ib, ik, 1), ib=1, NumberofSelectedBands)
        enddo
        close(outfileindex)
     endif

#if defined (MPI)
     call mpi_barrier(mpi_cmw, ierr)
#endif


     deallocate( W)
     deallocate( psi)
     deallocate( Hamk_bulk)
     deallocate( kxy)
     deallocate( kxy_shape)
     deallocate( eigv    )
     deallocate( eigv_mpi)
     deallocate( dos_selected)
     deallocate( dos_selected_mpi)
 
   return
   end subroutine orbitaltexture3D



   subroutine ek_bulk2D_spin
     ! Calculate bulk's energy band using wannier TB method for a given k plane
     ! Copyright (c) 2010 QuanSheng Wu. All rights reserved.

     use wmpi
     use para

     implicit none

     integer :: ik, i, j, ib
     integer :: nkx, nky
     integer :: knv3
     integer :: ierr
     integer :: nwann
     real(Dp) :: k(3)
     real(dp) :: sx
     real(dp) :: sy
     real(dp) :: sz
     real(Dp), allocatable :: W(:)
     real(dp), allocatable :: kxy(:,:)
     real(dp), allocatable :: kxy_shape(:,:)
     
     ! Hamiltonian of bulk system
     complex(Dp), allocatable :: Hamk_bulk(:, :)

     ! eigen value of H
     real(dp), allocatable :: gap(:)
     real(dp), allocatable :: gap_mpi(:)
     real(dp), allocatable :: eigv(:,:)
     real(dp), allocatable :: eigv_mpi(:,:)
     real(dp), allocatable :: spin(:, :, :)
     real(dp), allocatable :: spin_mpi(:, :, :)


     complex(dp), allocatable :: sigmax(:, :)
     complex(dp), allocatable :: sigmay(:, :)
     complex(dp), allocatable :: sigmaz(:, :)
     complex(dp), allocatable :: psi(:)

     nkx= Nk
     nky= Nk
     knv3= nkx*nky
     allocate( W( Num_wann))
     allocate( psi( Num_wann))
     allocate( Hamk_bulk( Num_wann, Num_wann))
     allocate( sigmax( Num_wann, Num_wann))
     allocate( sigmay( Num_wann, Num_wann))
     allocate( sigmaz( Num_wann, Num_wann))
     allocate( spin    (3, Num_wann, knv3))
     allocate( spin_mpi(3, Num_wann, knv3))
     allocate( kxy(3, knv3))
     allocate( kxy_shape(3, knv3))
     allocate( gap     (   knv3))
     allocate( gap_mpi (   knv3))
     allocate( eigv    (4, knv3))
     allocate( eigv_mpi(4, knv3))
     kxy = 0d0
     gap    = 0d0
     gap_mpi= 0d0
     eigv    = 0d0
     eigv_mpi= 0d0
     spin = 0d0
     spin_mpi = 0d0
     sigmax= 0d0
     sigmay= 0d0
     sigmaz= 0d0
 
 
     nwann= Num_wann
     if (soc>0) nwann= Num_wann/2
    !print *, 'nwann', nwann

     if (SOC==0) stop 'you should set soc=0 in the input file'
     !> construct spin matrix
     do i= 1, nwann
        sigmax(i, i+ nwann)= 1d0
        sigmax(i+ nwann, i)= 1d0
        sigmay(i, i+ nwann)= -zi
        sigmay(i+ nwann, i)= zi
        sigmaz(i, i)= 1d0
        sigmaz(i+ nwann, i+ nwann)= -1d0
     enddo

     ik =0
     do i= 1, nkx
     do j= 1, nky
        ik =ik +1
        kxy(:, ik)= K3D_start+ K3D_vec1*(i-1)/dble(nkx-1)+ K3D_vec2*(j-1)/dble(nky-1)
        kxy_shape(:, ik)= kxy(1, ik)* Origin_cell%Kua+ kxy(2, ik)* Origin_cell%Kub+ kxy(3, ik)* Origin_cell%Kuc 
     enddo
     enddo

     do ik= 1+cpuid, knv3, num_cpu
        if (cpuid==0) print * , ik, knv3

        k = kxy(:, ik)

        ! calculation bulk hamiltonian
        Hamk_bulk= 0d0
        call ham_bulk_atomicgauge(k, Hamk_bulk)
        Hamk_bulk= Hamk_bulk/eV2Hartree

        !> diagonalization by call zheev in lapack
        W= 0d0
        call eigensystem_c( 'V', 'U', Num_wann ,Hamk_bulk, W)

        do ib= 1, Num_wann
           psi(:)= Hamk_bulk(:, ib)
          !print *, sum(conjg(psi(:))*psi(:))
           sx= 0d0
           sy= 0d0
           sz= 0d0
           do i= 1, Num_wann
              do j= 1, Num_wann
                 sx= sx+ conjg(psi(i))* sigmax(i, j)* psi(j)
                 sy= sy+ conjg(psi(i))* sigmay(i, j)* psi(j)
                 sz= sz+ conjg(psi(i))* sigmaz(i, j)* psi(j)
              enddo ! j
           enddo ! i
           spin(1, ib, ik)= sx
           spin(2, ib, ik)= sy
           spin(3, ib, ik)= sz
          !print *, W(ib),sx*sx+sy*sy+sz*sz
        enddo ! ib

        !> for 44 bands
        eigv(:, ik)= W(Numoccupied-1:Numoccupied+ 2)
        gap (   ik)= W(Numoccupied+1)- W(Numoccupied)
     enddo

#if defined (MPI)
     call mpi_allreduce(eigv,eigv_mpi,size(eigv),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
     call mpi_allreduce(gap ,gap_mpi,size(gap ),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
     call mpi_allreduce(spin ,spin_mpi,size(spin ),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
     eigv_mpi= eigv
     gap_mpi= gap
     spin_mpi= spin
#endif

     outfileindex= outfileindex+ 1
     if (cpuid==0)then
        open(unit=outfileindex, file='bulkek2D.dat')
        do ik=1, knv3
           write(outfileindex, '(1000f19.9)')kxy_shape(:, ik)*Angstrom2atomic, eigv_mpi(:, ik), &
              (spin_mpi(:, ib, ik), ib=Numoccupied-1, Numoccupied+2)
        enddo
        close(outfileindex)
     endif

     outfileindex= outfileindex+ 1
     if (cpuid==0)then
        open(unit=outfileindex, file='gap2D.dat')
        do ik=1, knv3
           write(outfileindex, '(1000f19.9)')kxy(:, ik)*Angstrom2atomic, gap_mpi(ik)
        enddo
        close(outfileindex)
     endif

     deallocate( W)
     deallocate( psi)
     deallocate( Hamk_bulk)
     deallocate( sigmax)
     deallocate( sigmay)
     deallocate( sigmaz)
     deallocate( spin    )
     deallocate( spin_mpi)
     deallocate( kxy)
     deallocate( kxy_shape)
     deallocate( gap     )
     deallocate( gap_mpi )
     deallocate( eigv    )
     deallocate( eigv_mpi)
 
   return
   end subroutine ek_bulk2D_spin


  subroutine ek_bulk2D_spin_green
     !> using green's function

     use wmpi
     use para

     implicit none

     integer :: ik, i, j
     integer :: nkx, nky
     integer :: nwann
     integer :: knv3
     integer :: ierr
     real(Dp) :: k(3)
     real(Dp) :: k11(3), k12(3)
     real(Dp) :: k21(3), k22(3)
     real(dp) :: kxmin, kxmax, kymin, kymax
     real(dp), allocatable :: kxy(:,:)
     real(dp), allocatable :: kxy_shape(:,:)
     
     ! Hamiltonian of bulk system
     real(Dp), allocatable :: W(:)
     complex(Dp), allocatable :: Hamk_bulk(:, :) 

     ! eigen value of H
	  real(dp), allocatable :: gap(:)
	  real(dp), allocatable :: gap_mpi(:)
	  real(dp), allocatable :: eigv(:,:)
	  real(dp), allocatable :: eigv_mpi(:,:)
     real(dp), allocatable :: dos(:)
     real(dp), allocatable :: dos_mpi(:)
	  real(dp), allocatable :: spin(:, :)
	  real(dp), allocatable :: spin_mpi(:, :)


     complex(dp), allocatable :: sigmax(:, :)
     complex(dp), allocatable :: sigmay(:, :)
     complex(dp), allocatable :: sigmaz(:, :)
     complex(dp), allocatable :: psi(:)
     complex(dp), allocatable :: mat1(:, :)

     nkx= Nk
     nky= Nk
     knv3= nkx*nky
     allocate( W( Num_wann))
     allocate( psi( Num_wann))
     allocate( Hamk_bulk( Num_wann, Num_wann))
     allocate( sigmax( Num_wann, Num_wann))
     allocate( sigmay( Num_wann, Num_wann))
     allocate( sigmaz( Num_wann, Num_wann))
     allocate( mat1  ( Num_wann, Num_wann))
     allocate( dos     (knv3))
     allocate( dos_mpi (knv3))
     allocate( spin    (3, knv3))
     allocate( spin_mpi(3, knv3))
     allocate( kxy(3, knv3))
     allocate( kxy_shape(3, knv3))
     allocate( gap     (   knv3))
     allocate( gap_mpi (   knv3))
     allocate( eigv    (4, knv3))
     allocate( eigv_mpi(4, knv3))
     kxy     = 0d0
     gap     = 0d0
     gap_mpi = 0d0
     eigv    = 0d0
     eigv_mpi= 0d0
     spin    = 0d0
     spin_mpi= 0d0
     dos    = 0d0
     dos_mpi= 0d0
     sigmax= 0d0
     sigmay= 0d0
     sigmaz= 0d0

     nwann= Num_wann/2

     if (SOC==0) stop 'you should set soc=0 in the input file'
     !> construct spin matrix
     do i= 1, nwann
        sigmax(i, i+ nwann)= 1d0
        sigmax(i+ nwann, i)= 1d0
        sigmay(i, i+ nwann)= -zi
        sigmay(i+ nwann, i)= zi
        sigmaz(i, i)= 1d0
        sigmaz(i+ nwann, i+ nwann)= -1d0
     enddo


     k11=(/-0.0d0,  0.0d0, -0.0d0/) ! 
     k12=(/ 1.0d0,  0.0d0,  1.0d0/) ! X
     k21=(/ 0.0d0,  0.0d0,  0.0d0/) ! 
     k22=(/ 1.0d0,  1.0d0,  0.0d0/) ! Z


     ik =0
     do i= 1, nkx
     do j= 1, nky
        ik =ik +1
        kxy(:, ik)= k11+(k12-k11)*(i-1)/dble(nkx-1)+  k21+ (k22-k21)*(j-1)/dble(nky-1)
        kxy_shape(:, ik)= kxy(1, ik)* Origin_cell%Kua+ kxy(2, ik)* Origin_cell%Kub+ kxy(3, ik)* Origin_cell%Kuc 
     enddo
     enddo

     do ik= 1+cpuid, nkx*nky, num_cpu
        if (cpuid==0) print * , ik

        k = kxy(:, ik)

        ! calculation bulk hamiltonian
        Hamk_bulk= 0d0
        call ham_bulk_atomicgauge(k, Hamk_bulk)
        Hamk_bulk= Hamk_bulk/eV2Hartree

        !> diagonalization by call zheev in lapack
        Hamk_bulk= Hamk_bulk+ zi*Fermi_broadening
        call inv(Num_wann, Hamk_bulk)
        do i=1, Num_wann
           dos(ik)= dos(ik) -aimag(Hamk_bulk(i, i))
        enddo

        call mat_mul(Num_wann, Hamk_bulk, sigmax, mat1)

        do i=1, Num_wann
           spin(1, ik)= spin(1, ik) -aimag(mat1(i, i))
        enddo

        call mat_mul(Num_wann, Hamk_bulk, sigmay, mat1)
        do i=1, Num_wann
           spin(2, ik)= spin(2, ik) -aimag(mat1(i, i))
        enddo

        call mat_mul(Num_wann, Hamk_bulk, sigmaz, mat1)
        do i=1, Num_wann
           spin(3, ik)= spin(3, ik) -aimag(mat1(i, i))
        enddo

        !> for 44 bands
        eigv(:, ik)= W(Numoccupied-1:Numoccupied+ 2)
        gap (   ik)= W(Numoccupied+1)- W(Numoccupied)
     enddo !ik

#if defined (MPI)
     call mpi_allreduce(eigv,eigv_mpi,size(eigv),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
     call mpi_allreduce(gap ,gap_mpi,size(gap ),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
     call mpi_allreduce(spin ,spin_mpi,size(spin ),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
     call mpi_allreduce(dos ,dos_mpi,size(dos ),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
     gap_mpi= gap
     dos_mpi= dos
     eigv_mpi= eigv
     spin_mpi= spin
#endif


     outfileindex= outfileindex+ 1
     if (cpuid==0)then
        open(unit=outfileindex, file='bulkek2D.dat')
        do ik=1, knv3
           write(outfileindex, '(1000f20.6)')kxy_shape(:, ik)*Angstrom2atomic, eigv_mpi(:, ik) 
        enddo
        close(outfileindex)
     endif

     outfileindex= outfileindex+ 1
     if (cpuid==0)then
        open(unit=outfileindex, file='gap2D.dat')
        do ik=1, knv3
           write(outfileindex, '(1000f20.6)')kxy_shape(:, ik)*Angstrom2atomic, gap_mpi(ik)
        enddo
        close(outfileindex)
     endif

     outfileindex= outfileindex+ 1
     if (cpuid==0)then
        open(unit=outfileindex, file='bulkspin2D.dat')
        ik= 0
        do i= 1, nkx
           do j= 1, nky
              ik= ik+ 1

              if (ik==1 .or. ik==knv3.or. dos_mpi(ik)<4d0) cycle
              if ((dos_mpi(ik)> dos_mpi(ik-1)) .and. (dos_mpi(ik)> dos_mpi(ik+1))) then
                 write(outfileindex, '(10000f20.5)')kxy_shape(:, ik)*Angstrom2atomic, dos_mpi(ik) , spin_mpi(:, ik) 
              endif
           enddo
           write(outfileindex, '(a)')' '
        enddo
        close(outfileindex)
     endif

     deallocate( W)
     deallocate( psi)
     deallocate( Hamk_bulk)
     deallocate( sigmax)
     deallocate( sigmay)
     deallocate( sigmaz)
     deallocate( mat1  )
     deallocate( dos     )
     deallocate( dos_mpi )
     deallocate( spin    )
     deallocate( spin_mpi)
     deallocate( kxy)
     deallocate( kxy_shape)
     deallocate( gap     )
     deallocate( gap_mpi )
     deallocate( eigv    )
     deallocate( eigv_mpi)
 
   return
   end subroutine ek_bulk2D_spin_green



subroutine ek_bulk_mengyu
   ! Calculate bulk's energy bands using wannier TB method
   !
   ! Copyright (c) 2010 QuanSheng Wu. All rights reserved.

   use wmpi
   use para

   implicit none

   integer :: ik, i, j, knv3, ierr
   real(dp) :: emin, emax, k(3), kx, ky, kz, kxy_plane(3), kxy_cart(3), k1(3), k2(3)
   real(Dp), allocatable :: W(:)

   real(dp) :: kzbroden

   ! Hamiltonian of bulk system
   complex(Dp), allocatable :: Hamk_bulk(:, :)

   ! eigen value of H
   real(dp), allocatable :: eigv(:,:)
   real(dp), allocatable :: eigv_mpi(:,:)
   real(dp), allocatable :: weight(:,:,:)
   real(dp), allocatable :: weight_mpi(:,:,:)
   real(dp), allocatable :: weight_sum(:,:)

   kzbroden= omegamin

   knv3= nk3_band
   allocate(W(Num_wann))
   allocate(Hamk_bulk(Num_wann, Num_wann))
   allocate( eigv    (Num_wann, knv3))
   allocate( eigv_mpi(Num_wann, knv3))
   allocate( weight    (Num_wann,Num_wann, knv3))
   allocate( weight_mpi(Num_wann,Num_wann, knv3))
   allocate( weight_sum(Num_wann, knv3))
   eigv    = 0d0
   eigv_mpi= 0d0
   weight = 0d0
   weight_sum = 0d0
   weight_mpi = 0d0

   do ik= 1+cpuid, knv3, num_cpu
      !if (cpuid==0) write(stdout, *)'BulkBand, ik, knv3 ', ik, knv3

      !> in fractional coordinates
      k = kpath_3d(:, ik)

      kxy_cart(:)= k(1)* Origin_cell%Kua+ k(2)* Origin_cell%Kub+ k(3)* Origin_cell%Kuc
      call rotate_k3_to_kplane(kxy_cart(:), kxy_plane(:))
      kx=kxy_plane(1)
      ky=kxy_plane(2)

      !> due to the bending effect in the ARPES measurement
      kz=kxy_plane(3)+ kzbroden*dsqrt(kx*kx+ky*ky)
      k1(1)= kx
      k1(2)= ky
      k1(3)= kz
      call rotate_kplane_to_k3(k1, k2)
      call cart_direct_rec(k2, k(:))

      ! calculation bulk hamiltonian
      Hamk_bulk= 0d0

      ! generate bulk Hamiltonian
      if (index(KPorTB, 'KP')/=0)then
         call ham_bulk_kp(k, Hamk_bulk)
      else
         !> deal with phonon system
         if (index(Particle,'phonon')/=0.and.LOTO_correction) then
            call ham_bulk_LOTO(k, Hamk_bulk)
         else
            call ham_bulk_latticegauge    (k, Hamk_bulk)
            Hamk_bulk= Hamk_bulk/eV2Hartree
         endif
      endif

      !> diagonalization by call zheev in lapack
      W= 0d0
      call eigensystem_c('V', 'U', Num_wann ,Hamk_bulk, W)

      !if (sum(abs(k))<1e-5) call  orbital_momenta(k, Hamk_bulk)

      eigv(:, ik)= W
      do i=1, Num_wann  !> band
         if (soc==0) then
            do j=1, Num_wann  !> projector
               weight(j, i, ik)= abs(Hamk_bulk(j, i))**2
            enddo ! j
         else
            do j=1, Num_wann/2  !> projector
               weight(j, i, ik)= (abs(Hamk_bulk(j, i))**2+ &
                  abs(Hamk_bulk(j+ Num_wann/2, i))**2)
            enddo ! j
         endif
      enddo ! i
   enddo ! ik

#if defined (MPI)
   call mpi_allreduce(eigv,eigv_mpi,size(eigv),&
      mpi_dp,mpi_sum,mpi_cmw,ierr)
   call mpi_allreduce(weight, weight_mpi,size(weight),&
      mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
   eigv_mpi= eigv
   weight_mpi= weight
#endif

   if (index(Particle,'phonon')/=0) then
      eigv_mpi = eigv_mpi - MINVAL(eigv_mpi)
   endif

   !> deal with phonon system
   if (index(Particle,'phonon')/=0) then
      do ik=1, knv3
         do j=1, Num_wann
            eigv_mpi(j, ik)= sqrt(abs(eigv_mpi(j, ik)))*sign(1d0, eigv_mpi(j, ik))
            !eigv_mpi(j, ik)=     (abs(eigv_mpi(j, ik)))*sign(1d0, eigv_mpi(j, ik))
         enddo
      enddo
   endif

   do i=1, Num_wann
      do ik=1, knv3
         weight_sum(i, ik)= sum(weight_mpi(:, i, ik))
      enddo
   enddo

   weight= weight_mpi/maxval(weight_sum)*255d0

   outfileindex= outfileindex+ 1
   if (cpuid==0)then
      open(unit=outfileindex, file='bulkek.dat')

      do i=1, Num_wann
         do ik=1, knv3
            write(outfileindex, '(2f19.9, 10000i5)')k3len(ik)*Angstrom2atomic,eigv_mpi(i, ik), &
               int(weight(:, i, ik))
         enddo
         write(outfileindex, *)' '
      enddo
      close(outfileindex)
   endif

   !> minimum and maximum value of energy bands
   emin=  minval(eigv_mpi)-0.5d0
   emax=  maxval(eigv_mpi)+0.5d0

   call generate_ek_kpath_gnu('bulkek.dat', 'bulkek.gnu', 'bulkek.pdf', &
                                 emin, emax, knv3, Nk3lines, &
                                 k3line_name, k3line_stop, k3len)

   deallocate(W)
   deallocate(Hamk_bulk)
   deallocate( eigv    )
   deallocate( eigv_mpi)
   deallocate( weight    )
   deallocate( weight_mpi)
   deallocate( weight_sum)

   return
end subroutine ek_bulk_mengyu


subroutine ek_bulk_spin
   ! Calculate bulk's energy band using wannier TB method
   ! also calculate spin direction for each band and each kpoint

   use wmpi
   use para

   implicit none

   integer :: ik, i, j, ib
   integer :: knv3
   integer :: ierr
   integer :: nwann
   real(dp) :: sx, sy, sz
   real(dp) :: emin
   real(dp) :: emax
   real(dp) :: k(3)
   real(dp), allocatable :: W(:)

   ! Hamiltonian of bulk system
   complex(dp), allocatable :: Hamk_bulk(:, :)
   complex(dp), allocatable :: sigmax(:, :)
   complex(dp), allocatable :: sigmay(:, :)
   complex(dp), allocatable :: sigmaz(:, :)
   complex(dp), allocatable :: psi(:)

   ! eigen value of H
   real(dp), allocatable :: eigv(:,:)
   real(dp), allocatable :: eigv_mpi(:,:)
   real(dp), allocatable :: spin(:, :, :)
   real(dp), allocatable :: spin_mpi(:, :, :)

   knv3= nk3_band
   allocate( W( Num_wann))
   allocate( psi( Num_wann))
   allocate( Hamk_bulk( Num_wann, Num_wann))
   allocate( sigmax( Num_wann, Num_wann))
   allocate( sigmay( Num_wann, Num_wann))
   allocate( sigmaz( Num_wann, Num_wann))
   allocate( eigv    (Num_wann, knv3))
   allocate( eigv_mpi(Num_wann, knv3))
   allocate( spin    (3, Num_wann, knv3))
   allocate( spin_mpi(3, Num_wann, knv3))
   psi = 0d0
   eigv = 0d0
   eigv_mpi= 0d0
   spin= 0d0
   spin_mpi= 0d0
   sigmax= 0d0
   sigmay= 0d0
   sigmaz= 0d0

   nwann= Num_wann/2

   if (SOC==0) return !stop 'you should set soc=0 in the input file'
   !> construct spin matrix
   do i= 1, nwann
      sigmax(i, i+ nwann)= 1d0
      sigmax(i+ nwann, i)= 1d0
      sigmay(i, i+ nwann)= -zi
      sigmay(i+ nwann, i)= zi
      sigmaz(i, i)= 1d0
      sigmaz(i+ nwann, i+ nwann)= -1d0
   enddo

   do ik= 1+cpuid, knv3, num_cpu
      if (cpuid==0) write(stdout, *)'BulkBandSpin, ik, knv3 ', ik, knv3

      k = kpath_3d(:, ik)

      ! calculation bulk hamiltonian
      Hamk_bulk= 0d0
      call ham_bulk_latticegauge(k, Hamk_bulk)
      Hamk_bulk= Hamk_bulk/eV2Hartree

      !> diagonalization by call zheev in lapack
      W= 0d0
      call eigensystem_c( 'V', 'U', Num_wann ,Hamk_bulk, W)

      eigv(:, ik)= W

      do ib= 1, Num_wann
         psi(:)= Hamk_bulk(:, ib)
         sx= 0d0
         sy= 0d0
         sz= 0d0
         do i= 1, Num_wann
            do j= 1, Num_wann
               sx= sx+ real(conjg(psi(i))* sigmax(i, j)* psi(j))
               sy= sy+ real(conjg(psi(i))* sigmay(i, j)* psi(j))
               sz= sz+ real(conjg(psi(i))* sigmaz(i, j)* psi(j))
            enddo ! j
         enddo ! i
         spin(1, ib, ik)= sx
         spin(2, ib, ik)= sy
         spin(3, ib, ik)= sz
      enddo ! ib

   enddo ! ik

#if defined (MPI)
   call mpi_allreduce(eigv,eigv_mpi,size(eigv),&
      mpi_dp,mpi_sum,mpi_cmw,ierr)
   call mpi_allreduce(spin,spin_mpi,size(spin),&
      mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
   eigv_mpi= eigv
   spin_mpi= spin
#endif


   outfileindex= outfileindex+ 1
   if (cpuid==0)then
      open(unit=outfileindex, file='bulkek.dat')

      do i=1, Num_wann
         do ik=1, knv3
            write(outfileindex, '(100000f19.9)')k3len(ik)*Angstrom2atomic,eigv_mpi(i, ik), spin(:, i, ik)
         enddo
         write(outfileindex, *)' '
      enddo
      close(outfileindex)
   endif

   !> minimum and maximum value of energy bands
   emin= minval(eigv_mpi)-0.5d0
   emax= maxval(eigv_mpi)+0.5d0

   call generate_ek_kpath_gnu('bulkek.dat', 'bulkek.gnu', 'bulkek.pdf', &
                                 emin, emax, knv3, Nk3lines, &
                                 k3line_name, k3line_stop, k3len)

   deallocate( W)
   deallocate( psi)
   deallocate( Hamk_bulk)
   deallocate( sigmax)
   deallocate( sigmay)
   deallocate( sigmaz)
   deallocate( eigv    )
   deallocate( eigv_mpi)
   deallocate( spin    )
   deallocate( spin_mpi)

   return
end subroutine ek_bulk_spin

subroutine ek_bulk_mirror_z
   ! Use the eigenvalue of mirror z to label the energy bands

   use wmpi
   use para

   implicit none

   integer :: ik, i
   integer :: knv3
   integer :: ierr
   real(dp) :: emin
   real(dp) :: emax
   real(Dp) :: k(3)
   real(Dp), allocatable :: W(:)

   ! Hamiltonian of bulk system
   complex(Dp), allocatable :: Hamk_bulk(:, :)
   complex(Dp), allocatable :: Hamk     (:, :)
   complex(Dp), allocatable :: mat1     (:, :)
   complex(Dp), allocatable :: mat2     (:, :)

   ! eigen value of H
   real(dp), allocatable :: eigv(:,:)
   real(dp), allocatable :: eigv_mpi(:,:)
   logical, allocatable :: mirror_plus(:, :)
   logical, allocatable :: mirror_minus(:, :)

   knv3= nk3_band
   allocate(W(Num_wann))
   allocate(Hamk_bulk(Num_wann, Num_wann))
   allocate(Hamk     (Num_wann, Num_wann))
   allocate(mat1     (Num_wann, Num_wann))
   allocate(mat2     (Num_wann, Num_wann))
   allocate( eigv    (Num_wann, knv3))
   allocate( eigv_mpi(Num_wann, knv3))
   allocate( mirror_plus(Num_wann, knv3))
   allocate( mirror_minus(Num_wann, knv3))
   mirror_plus= .False.
   mirror_minus= .False.
   eigv    = 0d0
   eigv_mpi= 0d0

   do ik= 1+cpuid, knv3, num_cpu
      if (cpuid==0) write(stdout, *)'BulkBandmirrorz, ik, knv3 ', ik, knv3

      k = kpath_3d(:, ik)

      ! calculation bulk hamiltonian
      Hamk_bulk= 0d0
      call ham_bulk_atomicgauge    (k, Hamk_bulk)
      Hamk_bulk= Hamk_bulk/eV2Hartree

      k = kpath_3d(:, ik)
      !k(2)= -k(2)
      Hamk= 0d0
      call ham_bulk_latticegauge    (k, Hamk)
      Hamk= Hamk/eV2Hartree

      !> symmetrization
     !call mat_mul(Num_wann, mirror_z, hamk, mat1)
     !call mat_mul(Num_wann, mat1, mirror_z, mat2)
     !hamk= (Hamk_bulk+ mat2)/2.d0

      !> diagonal hamk
      call eigensystem_c('V', 'U', Num_wann, hamk, W)

      mat2= conjg(transpose(hamk))

      !> calculate mirror eigenvalue
      call mat_mul(Num_wann, mat2, mirror_z, mat1)
      call mat_mul(Num_wann, mat1, hamk, mat2)

      !> get mirror_plus and mirror_minus
      do i=1, Num_wann
         print *, real(mat2(i, i))
         if (abs(real(mat2(i, i))-1d0)< 1e-1) then
            mirror_plus(i, ik)= .true.
         else
            mirror_minus(i, ik)= .true.
         endif
      enddo
     !pause

      !if (cpuid.eq.0)write(*, *)ik,&
      !   (mirror_plus(i, ik), i=1, Num_wann)
      eigv(:, ik)= W

   enddo  ! ik

#if defined (MPI)
   call mpi_allreduce(eigv,eigv_mpi,size(eigv),&
      mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
   eigv_mpi= eigv
#endif


   outfileindex= outfileindex+ 1
   if (cpuid==0)then
      open(unit=outfileindex, file='bulkek.dat')
      do i=1, Num_wann
         do ik=1, knv3
            write(outfileindex, '(1000f19.9)')k3len(ik)*Angstrom2atomic,eigv_mpi(i, ik)
         enddo
         write(outfileindex, *)' '
      enddo
      close(outfileindex)
   endif

   outfileindex= outfileindex+ 1
   if (cpuid==0)then
      open(unit=outfileindex, file='bulkekmirrorplus.dat')
      do i=1, Num_wann
         do ik=1, knv3
            if (mirror_plus(i, ik)) then
               write(outfileindex, '(1000f19.9)')k3len(ik)*Angstrom2atomic,eigv_mpi(i, ik)
            endif
         enddo
         write(outfileindex, *)' '
      enddo
      close(outfileindex)
   endif


   outfileindex= outfileindex+ 1
   if (cpuid==0)then
      open(unit=outfileindex, file='bulkekmirrorminus.dat')
      do i=1, Num_wann
         do ik=1, knv3
            if (.not.mirror_plus(i, ik)) then
               write(outfileindex, '(1000f19.9)')k3len(ik)*Angstrom2atomic,eigv_mpi(i, ik)
            endif
         enddo
         write(outfileindex, *)' '
      enddo
      close(outfileindex)
   endif

   !> minimum and maximum value of energy bands
   emin= minval(eigv_mpi)-0.5d0
   emax= maxval(eigv_mpi)+0.5d0

   call generate_ek_kpath_gnu('bulkek.dat', 'bulkek.gnu', 'bulkek.pdf', &
                                 emin, emax, knv3, Nk3lines, &
                                 k3line_name, k3line_stop, k3len)

   deallocate(W)
   deallocate(Hamk_bulk)
   deallocate(Hamk     )
   deallocate(mat1     )
   deallocate(mat2     )
   deallocate( eigv    )
   deallocate( eigv_mpi)
   deallocate( mirror_plus)
   deallocate( mirror_minus)

   return
end subroutine ek_bulk_mirror_z




subroutine ek_bulk_mirror_x
   !> calculate bulk band for a given mirror symmetry

   use wmpi
   use para

   implicit none

   integer :: ik, i
   integer :: knv3
   integer :: ierr
   real(dp) :: emin
   real(dp) :: emax
   real(Dp) :: k(3)

   ! Hamiltonian of bulk system
   real(Dp), allocatable :: W(:)

   ! Hamiltonian of bulk system
   complex(Dp), allocatable :: Hamk_bulk(:, :)
   complex(Dp), allocatable :: Hamk     (:, :)
   complex(Dp), allocatable :: mat1     (:, :)
   complex(Dp), allocatable :: mat2     (:, :)
   complex(dp), allocatable :: inv_mirror_x(:, :)

   ! eigen value of H
   real(dp), allocatable :: eigv(:,:)
   real(dp), allocatable :: eigv_mpi(:,:)
   logical, allocatable :: mirror_plus(:, :)
   logical, allocatable :: mirror_minus(:, :)


   knv3= nk3_band
   allocate(W(Num_wann))
   allocate(Hamk_bulk(Num_wann, Num_wann))
   allocate(Hamk     (Num_wann, Num_wann))
   allocate(mat1     (Num_wann, Num_wann))
   allocate(mat2     (Num_wann, Num_wann))
   allocate(inv_mirror_x     (Num_wann, Num_wann))
   allocate( eigv    (Num_wann, knv3))
   allocate( eigv_mpi(Num_wann, knv3))
   allocate( mirror_plus(Num_wann, knv3))
   allocate( mirror_minus(Num_wann, knv3))
   mirror_plus= .False.
   mirror_minus= .False.
   eigv    = 0d0
   eigv_mpi= 0d0

   do ik= 1+cpuid, knv3, num_cpu
      if (cpuid==0) write(stdout, *)'EkBulk_mirror, ik, knv3 ', ik, knv3

      k = kpath_3d(:, ik)

      ! calculation bulk hamiltonian
      Hamk_bulk= 0d0
      call ham_bulk_latticegauge    (k, Hamk_bulk)
      Hamk_bulk= Hamk_bulk/eV2Hartree

      !k = kpath_3d(:, ik)
      !k(1)= -k(1)
      !Hamk= 0d0
      !call ham_bulk_latticegauge    (k, Hamk)


      hamk=  Hamk_bulk
      !> symmetrization
      inv_mirror_x=mirror_x

      call inv(Num_wann, inv_mirror_x)

      !print *, maxval(abs(inv_mirror_x-mirror_x))
      !stop

      call mat_mul(Num_wann, inv_mirror_x, hamk, mat1)
      call mat_mul(Num_wann, mat1, mirror_x, mat2)
      hamk= (Hamk_bulk+ mat2)/2.d0
      !hamk=  mat2

      !> diagonal hamk
      call eigensystem_c('V', 'U', Num_wann, hamk, W)

      mat2= conjg(transpose(hamk))

      !> calculate mirror eigenvalue
      call mat_mul(Num_wann, mat2, mirror_x, mat1)
      call mat_mul(Num_wann, mat1, hamk, mat2)

      !> get mirror_plus and mirror_minus
      do i=1, Num_wann
         if (abs(real(mat2(i, i))-1d0)< 1e-3) then
            mirror_plus(i, ik)= .true.
         else
            mirror_minus(i, ik)= .true.
         endif
      enddo

      !if (cpuid.eq.0)write(*, *)ik,&
      !   (mirror_plus(i, ik), i=1, Num_wann)
      eigv(:, ik)= W

   enddo  ! ik

#if defined (MPI)
   call mpi_allreduce(eigv,eigv_mpi,size(eigv),&
      mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
   eigv_mpi= eigv
#endif

   outfileindex= outfileindex+ 1
   if (cpuid==0)then
      open(unit=outfileindex, file='bulkek.dat')
      do i=1, Num_wann
         do ik=1, knv3
            write(outfileindex, '(1000f19.9)')k3len(ik)*Angstrom2atomic,eigv_mpi(i, ik)
         enddo
         write(outfileindex, *)' '
      enddo
      close(outfileindex)
   endif

   outfileindex= outfileindex+ 1
   if (cpuid==0)then
      open(unit=outfileindex, file='bulkekmirrorplus.dat')
      do i=1, Num_wann
         do ik=1, knv3
            if (mirror_plus(i, ik)) then
               write(outfileindex, '(1000f19.9)')k3len(ik)*Angstrom2atomic,eigv_mpi(i, ik)
            endif
         enddo
         write(outfileindex, *)' '
      enddo
      close(outfileindex)
   endif

   outfileindex= outfileindex+ 1
   if (cpuid==0)then
      open(unit=outfileindex, file='bulkekmirrorminus.dat')
      do i=1, Num_wann
         do ik=1, knv3
            if (.not.mirror_plus(i, ik)) then
               write(outfileindex, '(1000f19.9)')k3len(ik)*Angstrom2atomic,eigv_mpi(i, ik)
            endif
         enddo
         write(outfileindex, *)' '
      enddo
      close(outfileindex)
   endif

   !> minimum and maximum value of energy bands
   emin= minval(eigv_mpi)-0.5d0
   emax= maxval(eigv_mpi)+0.5d0

   call generate_ek_kpath_gnu('bulkek.dat', 'bulkek.gnu', 'bulkek.pdf', &
                                 emin, emax, knv3, Nk3lines, &
                                 k3line_name, k3line_stop, k3len)

   deallocate(W)
   deallocate(Hamk_bulk)
   deallocate(Hamk     )
   deallocate(mat1     )
   deallocate(mat2     )
   deallocate(inv_mirror_x )
   deallocate( eigv    )
   deallocate( eigv_mpi)
   deallocate( mirror_plus)
   deallocate( mirror_minus)

   return
end subroutine ek_bulk_mirror_x

subroutine get_projection_weight_bulk(k_SBZ_direct, numk, psi, weight)
   !> calculate the weights of given selected orbitals and mode for given wavefunction psi
   !> There are two modes. One is project on the selected orbitals. 
   !> The other one is projected on a special mode of another lattice which is usually a folded lattice.
   !> For the first mode, numk should be one.
   use para, only : dp, num_wann, NumberofSelectedOrbitals_groups, projection_weight_mode, &
      NumberofSelectedOrbitals, Selected_WannierOrbitals, cpuid, stdout, Nrpts, irvec, &
      Origin_cell, k3points_unfold_pointmode_direct, pi2zi
   implicit none

   !> psi is a given vector for a k point and one band
   integer, intent(in) :: numk

   !> a k vector associated with the wave function psi in the supercell Brillouin zone
   real(dp), intent(in) :: k_SBZ_direct(3)
   complex(dp), intent(in) :: psi(num_wann)

   real(dp), intent(out) :: weight(NumberofSelectedOrbitals_groups, numK)

   !> local variables
   integer :: ig, ik, i, j
   real(dp) :: k_PBZ_direct(3)
   real(dp), allocatable :: weight_t(:)
   complex(dp) :: overlp

   allocate(weight_t(NumberofSelectedOrbitals_groups))
   weight= 0d0
   weight_t= 0d0

   if (trim(adjustl(projection_weight_mode))=='NORMAL'.and.numk/=1) then
      if (cpuid==0) write(stdout, *) '  Warning: There is something wrong with get_projection_weight_bulk calling.' 
   endif

   if (trim(adjustl(projection_weight_mode))=='NORMAL') then
      do ig=1, NumberofSelectedOrbitals_groups
         do i= 1, NumberofSelectedOrbitals(ig)
            j= Selected_WannierOrbitals(ig)%iarray(i)
            weight(ig, 1)= weight(ig, 1)+ abs(psi(j))**2
         enddo ! sweep the selected orbitals
      enddo ! ig
   elseif (trim(adjustl(projection_weight_mode))=='FOLDEDKPOITNS') then
      do ik=1, numk
         k_PBZ_direct= k3points_unfold_pointmode_direct(:, ik)
         call get_projection_weight_bulk_unfold(k_SBZ_direct, k_PBZ_direct, psi, weight_t)
         weight(:, ik)= weight_t(:)
      enddo ! ik
   endif


   return
end subroutine get_projection_weight_bulk

subroutine generate_ek_kpath_gnu(datafilename, gnufilename, gnuoutfilename, &
                                 emin, emax, num_k, nklines, &
                                 kname, k_stop, klen)
   use para
   implicit none

   character(*), intent(in) :: gnufilename, datafilename, gnuoutfilename

   !> in unit of eV
   real(dp), intent(in) :: emin, emax

   !> number of segments of the kpath
   integer, intent(in) :: nklines, num_k
   character(*), intent(in) :: kname(nklines+1)
   real(dp), intent(in) :: k_stop(nklines+1)
   real(dp), intent(in) :: klen(num_k)

   integer :: i

    !> write script for gnuplot
   outfileindex= outfileindex+ 1
   if (cpuid==0) then
      open(unit=outfileindex, file=gnufilename)
      write(outfileindex, '(a)') '# requirement: gnuplot version>5.4'
      write(outfileindex, '(2a)') '# Please open the data file to check the data:   ', trim(adjustl(datafilename))
      write(outfileindex, '(a)') 'set terminal pdf enhanced color font ",24"'
      write(outfileindex,'(2a)') 'set palette defined ( 0  "green", ', &
         '5 "yellow", 10 "red" )'
      write(outfileindex, '(3a)')"set output '", trim(adjustl(gnuoutfilename)), "' "
      write(outfileindex, '(a)')'set style data linespoints'
      write(outfileindex, '(a)')'unset key'
      write(outfileindex, '(a)')'set pointsize 0.8'
      write(outfileindex, '(a)')'#set xtics font ",24"'
      write(outfileindex, '(a)')'#set ytics font ",24"'
      write(outfileindex, '(a)')'#set ylabel font ",24"'
      write(outfileindex, '(a)')'set ylabel offset 0.5,0'
      write(outfileindex, '(a)')'set border lw 2'
      write(outfileindex, '(a, f10.5, a)')'set xrange [0: ', maxval(klen*Angstrom2atomic), ']'
      write(outfileindex, '(a,f12.6)')'emin=', emin
      write(outfileindex, '(a,f12.6)')'emax=', emax
      if (index(Particle,'phonon')/=0) then
         write(outfileindex, '(a)')'set yrange [0: emax]'
         write(outfileindex, '(a)')'set ylabel "Frequency (THz)"'
      else
         write(outfileindex, '(a)')'set ylabel "Energy (eV)"'
         write(outfileindex, '(a)')'set yrange [ emin : emax ]'
      endif
      write(outfileindex, 202, advance="no") (kname(i), k_stop(i)*Angstrom2atomic, i=1, nklines)
      write(outfileindex, 203)kname(nklines+1), k_stop(nklines+1)*Angstrom2atomic

      do i=1, nklines-1
         if (index(Particle,'phonon')/=0) then
            write(outfileindex, 204)k_stop(i+1)*Angstrom2atomic, '0.0', k_stop(i+1)*Angstrom2atomic, 'emax'
         else
            write(outfileindex, 204)k_stop(i+1)*Angstrom2atomic, 'emin', k_stop(i+1)*Angstrom2atomic, 'emax'
         endif
      enddo
      write(outfileindex, '(2a)')"# please comment the following lines to plot the fatband "
      if (Is_Sparse) then
         write(outfileindex, '(4a)')"plot '", trim(adjustl(datafilename)), "' u 1:2 ",  &
            " w p pt 7  ps 0.2 lc rgb 'black', 0 w l lw 2 dt 2"
      else
         write(outfileindex, '(4a)')"plot '", trim(adjustl(datafilename)), "' u 1:2 ",  &
            " w lp lw 2 pt 7  ps 0.2 lc rgb 'black', 0 w l lw 2 dt 2"
      endif
      write(outfileindex, '(2a)')" " 
      write(outfileindex, '(2a)')"# uncomment the following lines to plot the fatband "
      write(outfileindex, '(4a)')"#plot '", trim(adjustl(datafilename)), "' u 1:2:3 ",  &
         " w lp lw 2 pt 7  ps 0.2 lc palette, 0 w l lw 2 dt 2"
      write(outfileindex, '(4a)')"# uncomment the following lines to plot the spin if necessary"
      write(outfileindex, '(4a)')"#plot '", trim(adjustl(datafilename)), "' u 1:2 ",  &
         "w lp lw 2 pt 7  ps 0.2, \"
      write(outfileindex, '(4a)')"#plot '", trim(adjustl(datafilename)), "' u 1:2:($3/6):($4/6) ",  &
         "w vec"
      close(outfileindex)
   endif

202 format('set xtics (',20('"',A3,'" ',F10.5,','))
203 format(A3,'" ',F10.5,')')
204 format('set arrow from ',F10.5,',',A5,' to ',F10.5,',',A5, ' nohead lw 2')

    return
end subroutine generate_ek_kpath_gnu
