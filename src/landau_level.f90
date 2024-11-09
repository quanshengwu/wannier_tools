!> calculate landau levels in 3D system in special k line
!> fix B field
!> the magnetic field is in the a1 a2 plane
!> B= (B cos\theta, B sin\theta, 0)
!>Periodic Landau gauge by Guan Yifei. contact: guan.yifei@outlook.com
!> construct on Dec 8 2015
!> By QuanSheng Wu at ETH Zurich Honggerberg
! Copyright (c) 2010 QuanSheng Wu. All rights reserved.
!
!> At present, only \theta=0 is well tested
!include "sparse.f90"
subroutine landau_level_k
   use para
   implicit none


   !> magnetic supercell size, perpendicular to the magnetic field, Ndimq= Nq*Num_wann
   integer :: Nq, Ndimq

   integer :: ik, i, j, ib, ierr, iq, ig

   !> lowest and highest band index to be calculated
   integer :: il, iu

   !> mdimq= iu-il+1, number of LLs to be calculated
   integer :: mdimq


   real(dp) :: B0, B0Tesla, B0Tesla_quantumflux_magsupcell
   real(dp) :: time_start, time_end, time_start0

   real(dp) :: theta
   real(dp) :: emin, emax

   ! wave vector
   real(dp) :: k3(3)

   !> dim= Ndimq, knv3
   real(dp), allocatable :: W(:)
   real(dp), allocatable :: W_full(:)
   real(dp), allocatable :: eigv(:, :)
   real(dp), allocatable :: eigv_mpi(:, :)
   complex(dp), allocatable :: psi(:)

   !> print the weight for the Selected_WannierOrbitals
   real(dp), allocatable :: dos_selected(:, :, :)
   real(dp), allocatable :: dos_selected_mpi(:, :, :)

   real(dp), allocatable :: dos_l_selected(:, :, :)
   real(dp), allocatable :: dos_l_selected_mpi(:, :, :)
   real(dp), allocatable :: dos_r_selected(:, :, :)
   real(dp), allocatable :: dos_r_selected_mpi(:, :, :)

   !> dim= Ndimq*Ndimq
   complex(dp), allocatable :: eigvec(:, :)
   complex(dp), allocatable :: ham_landau(:, :)

   Nq= Magq
   Ndimq= Num_wann* Nq

   il= Numoccupied*Nq- 2000
   iu= Numoccupied*Nq+ 2000
   il=1
   iu=Ndimq
   if (il< 1) il= 1
   if (iu> Ndimq) iu= Ndimq
   mdimq= iu-il+1

   if (cpuid==0) then
      write(stdout, '(a,i8)')' >> The magnetic-supercell size is ', Nq
      write(stdout, '(a,i8)')' >> The dimension of magnetic-supercell is ', Ndimq
      write(stdout, '(a,i8,a,i8)')' >> calculated LLs from il', il, ' to', iu, ' bands'
      write(stdout, '(a,i8,a,i8)')' >> in total, there are  ', mdimq, ' bands'
   endif

   allocate( ham_landau(Ndimq, Ndimq))
   allocate( W( mdimq))
   allocate( W_full( Ndimq))
   allocate( eigv( mdimq, nk3_band))
   allocate( eigv_mpi( mdimq, nk3_band))
   allocate( eigvec( Ndimq, mdimq))
   allocate( psi(Ndimq))
   allocate( dos_selected     (mdimq,   nk3_band, NumberofSelectedOrbitals_groups))
   allocate( dos_selected_mpi (mdimq,   nk3_band, NumberofSelectedOrbitals_groups))
   if (landau_chern_calc) then
      allocate( dos_l_selected     (mdimq,   nk3_band, NumberofSelectedOrbitals_groups))
      allocate( dos_l_selected_mpi (mdimq,   nk3_band, NumberofSelectedOrbitals_groups))
      allocate( dos_r_selected     (mdimq,   nk3_band, NumberofSelectedOrbitals_groups))
      allocate( dos_r_selected_mpi (mdimq,   nk3_band, NumberofSelectedOrbitals_groups))
      dos_l_selected= 0d0
      dos_l_selected_mpi= 0d0
      dos_r_selected= 0d0
      dos_r_selected_mpi= 0d0
   endif
   dos_selected= 0d0
   dos_selected_mpi= 0d0
   eigv_mpi   = 0d0
   eigv       = 0d0
   eigvec     = 0d0
   ham_landau = 0d0


   !> deal with the magnetic field parallel with R1' direction which 
   !> is the same as the first vector in SURFACE card.
   !> the flux in the unit cell not the magnetic supercell
   theta=0d0
   B0= 2d0*pi/dble(Nq)*Magp
   Bx= B0* dcos(theta)
   By= B0* dsin(theta)


   !> transform it into Tesla
   !> B=2*pi*\phi_0/S0, where \phi_0 is the quantum flux, S0 is the projected area
   !> \phi_0 = h/2e  h=6.62607004*1E-34, e= 1.6*1E-19
   !> B0Tesla_quantumflux_magsupcell is the magnetic field that makes the flux through the magnetic
   !> supercell to be one quantum flux
   B0Tesla_quantumflux_magsupcell= 6.62607004*1E-34/1.6021766208*1E19/MagneticSuperProjectedArea*1E20
   B0Tesla= 6.62607004*1E-34/1.6021766208*1E19/MagneticSuperProjectedArea*1E20*Magp

   if (cpuid==0) then
      write(stdout, '(a, 2E18.8)')' Magnetic field B in Tesla that makes the flux through '
      write(stdout, '(a, 2E18.8)')' the magnetic supercell to be one quantum flux : ', B0Tesla_quantumflux_magsupcell
      write(stdout, '(a, 2E18.8)')' Magnetic field B in Tesla that fits to Magp : ', B0Tesla
      write(stdout, '(a, 2f18.8)')' Magnetic field Bx, By= ', Bx, By
   endif

   !> calculate the landau levels along special k line
   k3= 0
   time_start= 0d0
   time_start0= 0d0
   call now(time_start0)
   time_start= time_start0
   time_end  = time_start0
   do ik=1+ cpuid, nk3_band, num_cpu
      if (cpuid.eq.0) &
         write(stdout, '(2a, i9, "  /", i10, a, f10.1, "s", a, f10.1, "s")') &
         'In LandauLevel_kband ', ' ik/NK ', ik, nk3_band, &
         ' time elapsed: ', time_end-time_start0, &
         ' time left: ', ((nk3_band-ik)/num_cpu)*(time_end-time_start)

      call now(time_start)

      k3= kpath_3d(:, ik)
      call ham_3Dlandau(Ndimq, Nq, k3, ham_landau)
      ham_landau= ham_landau/eV2Hartree

      !> diagonalization by call zheev in lapack
      !W_full= 0d0
      if (mdimq<2000) then
         if (LandauLevel_wavefunction_calc) then
            call eigensystem_c( 'V', 'U', Ndimq ,ham_landau, W_full)
            eigv(:, ik)= W_full(il:iu)
            eigvec(1:Ndimq, 1:mdimq)= ham_landau(1:Ndimq, 1:mdimq)
         else
            call eigensystem_c( 'N', 'U', Ndimq ,ham_landau, W_full)
            eigv(:, ik)= W_full(il:iu)
            eigvec(1:Ndimq, 1:mdimq)= 0d0
         endif
      else
         W= 0d0
         if (LandauLevel_wavefunction_calc) then
            call zheevx_pack('V', 'U', Ndimq, il, iu, ham_landau, W, eigvec)
         else
            call zheevx_pack('N', 'U', Ndimq, il, iu, ham_landau, W, eigvec)
         endif
         eigv(:, ik)= W
      endif

      !> calculate the weight on the selected orbitals
      do ib= 1, mdimq
         psi(:)= eigvec(:, ib)  !> the eigenvector of ib'th band
         do ig=1, NumberofSelectedOrbitals_groups
            do iq=1, Nq
               do i= 1, NumberofSelectedOrbitals(ig)
                  j= Num_wann*(iq-1)+ Selected_WannierOrbitals(ig)%iarray(i)
                  dos_selected(ib, ik, ig)= dos_selected(ib, ik, ig)+ abs(psi(j))**2
               enddo ! sweep the selected orbitals
            enddo ! iq sweep the magnetic supercell
            if (landau_chern_calc) then
               do iq=1, 1  ! edge states
                  if (iq>Nq) cycle
                  do i= 1, NumberofSelectedOrbitals(ig)
                     j= Num_wann*(iq-1)+ Selected_WannierOrbitals(ig)%iarray(i)
                     dos_l_selected(ib, ik, ig)= dos_l_selected(ib, ik, ig)+ abs(psi(j))**2
                  enddo ! sweep the selected orbitals
               enddo ! iq sweep the magnetic supercell
               do iq=Nq, Nq ! edge states
                  if (iq<1) cycle
                  do i= 1, NumberofSelectedOrbitals(ig)
                     j= Num_wann*(iq-1)+ Selected_WannierOrbitals(ig)%iarray(i)
                     dos_r_selected(ib, ik, ig)= dos_r_selected(ib, ik, ig)+ abs(psi(j))**2
                  enddo ! sweep the selected orbitals
               enddo ! iq sweep the magnetic supercell
            endif
         enddo ! ig
      enddo ! ib sweep the eigenvalue

      call now(time_end)
   enddo !ik

#if defined (MPI)
   call mpi_allreduce(eigv, eigv_mpi, size(eigv), &
      mpi_dp, mpi_sum, mpi_cmw, ierr)
   call mpi_allreduce(dos_selected, dos_selected_mpi,size(dos_selected),&
      mpi_dp,mpi_sum,mpi_cmw,ierr)

   if (landau_chern_calc) then
      call mpi_allreduce(dos_l_selected, dos_l_selected_mpi,size(dos_l_selected),&
         mpi_dp,mpi_sum,mpi_cmw,ierr)
      call mpi_allreduce(dos_r_selected, dos_r_selected_mpi,size(dos_r_selected),&
         mpi_dp,mpi_sum,mpi_cmw,ierr)
   endif

#else
   eigv_mpi= eigv
   dos_selected_mpi= dos_selected
   if (landau_chern_calc) then
      dos_l_selected_mpi= dos_l_selected
      dos_r_selected_mpi= dos_r_selected
   endif
#endif

   !> minimum and maximum value of energy bands
   emin= minval(eigv_mpi)*1.05-maxval(eigv_mpi)*0.05
   emax= -minval(eigv_mpi)*0.05+maxval(eigv_mpi)*1.05

   outfileindex= outfileindex+ 1
   if (cpuid.eq.0) then
      open(unit=outfileindex, file='landaulevel_k.dat')
      if (LandauLevel_wavefunction_calc) then
         write(outfileindex, '("#", a14, f16.8)')'Magnetic field in Tesla: ', B0Tesla
         if (landau_chern_calc) then
            write(outfileindex, '("#", a14, a15, 2a16)')'k ', ' Eig of LL', ' left side', ' right side'
         else
            write(outfileindex, '("#", a14, a15, a)')'k ', ' Eig of LL', ' Weight on the selected orbitals'
         endif
         do j=1, mdimq
            do i=1, nk3_band
               if (landau_chern_calc) then
                  write(outfileindex,'(20f16.8)')k3len_mag(i)*Angstrom2atomic, eigv_mpi(j, i), &
                     (dos_l_selected_mpi(j, i, ig), dos_r_selected_mpi(j, i, ig), &
                     ig=1, NumberofSelectedOrbitals_groups)
               else
                  write(outfileindex,'(20f16.8)')k3len_mag(i)*Angstrom2atomic, eigv_mpi(j, i), &
                     (dos_selected_mpi(j, i, ig), ig=1, NumberofSelectedOrbitals_groups)
               endif
            enddo
            write(outfileindex , *)' '
         enddo
      else
         write(outfileindex, '("#", a14, f16.8)')'Magnetic field in Tesla: ', B0Tesla
         write(outfileindex, '("#", a14, a15, a)')'k ', ' Eig of LL', ' Weight on the selected orbitals'
         do j=1, mdimq
            do i=1, nk3_band
               write(outfileindex,'(20f16.8)')k3len_mag(i)*Angstrom2atomic, eigv_mpi(j, i)
            enddo
            write(outfileindex , *)' '
         enddo
      endif

      close(outfileindex)
      write(stdout,*) 'Calculate landau level done'
   endif

   outfileindex= outfileindex+ 1
   if (cpuid==0) then
      open(unit=outfileindex, file='landaulevel_k.gnu')
      write(outfileindex, '(a)')"set encoding iso_8859_1"
      write(outfileindex, '(a)')'#set terminal  postscript enhanced color'
      write(outfileindex, '(a)')"#set output 'landaulevel_k.eps'"
      write(outfileindex, '(3a)')'#set terminal  pngcairo truecolor enhanced', &
         '  font ",60" size 1920, 1680'
      write(outfileindex, '(3a)')'set terminal  png truecolor enhanced', &
         '  font ",60" size 1920, 1680'
      write(outfileindex, '(a)')"set output 'landaulevel_k.png'"
      write(outfileindex, '(a)')'unset key'
      write(outfileindex, '(a)')'set pointsize 0.8'
      write(outfileindex, '(a)')'set border lw 3 '
      write(outfileindex, '(a)')'#set xtics font ",36"'
      write(outfileindex, '(a)')'#set ytics font ",36"'
      write(outfileindex, '(a)')'#set xlabel font ",36"'
      write(outfileindex, '(a)')'set ylabel "Energy (eV)"'
      write(outfileindex, '(a, f10.5, a)')'set xrange [0: ', maxval(k3len_mag*Angstrom2atomic), ']'
      write(outfileindex, '(a, i6, a, i6, a)')'set title "Landau level with p/q=', magp, "/", Nq,  '"'
      write(outfileindex, '(a)')'#set xtics offset 0, -1'
      write(outfileindex, '(a)')'#set ylabel offset -1, 0 '
      write(outfileindex, '(a)')'rgb(r,g,b) = int(r)*65536 + int(g)*256 + int(b)'
      write(outfileindex, 202, advance="no") (k3line_name(i), k3line_mag_stop(i)*Angstrom2atomic, i=1, nk3lines)
      write(outfileindex, 203)k3line_name(nk3lines+1), k3line_mag_stop(nk3lines+1)*Angstrom2atomic
      write(outfileindex, '(a, f10.5, a, f10.5, a)')'set yrange [', emin, ':', emax, ']'

      do i=1, nk3lines-1
         write(outfileindex, 204)k3line_mag_stop(i+1)*Angstrom2atomic, emin, &
            k3line_mag_stop(i+1)*Angstrom2atomic, emax
      enddo

202   format('set xtics (',:20('"',A3,'" ',F10.5,','))
203   format(A3,'" ',F10.5,')')
204   format('set arrow from ',F10.5,',',F10.5,' to ',F10.5,',',F10.5, ' nohead')

      if (LandauLevel_wavefunction_calc) then
         if (landau_chern_calc) then
            write(outfileindex,'(2a)') 'set palette defined ( -1  "blue", ', &
               '0 "grey", 1 "red" )'
            write(outfileindex, '(2a)')"plot 'landaulevel_k.dat' u 1:2:($4-$3)",  &
               "w p pt 7  ps 2 lc palette"
         else
            write(outfileindex, '(2a)')"plot 'landaulevel_k.dat' u 1:2:(rgb(255,255-255*$3, 3)) ",  &
               " w p  pt 7  ps 2 lc rgb variable"
         endif
      else
         write(outfileindex, '(2a)')"plot 'landaulevel_k.dat' u 1:2",  &
            "w p pt 7  ps 2"
      endif
      close(outfileindex)
   endif


   deallocate( ham_landau)
   deallocate( W)
   deallocate( W_full)
   deallocate( eigv)
   deallocate( eigv_mpi)
   deallocate( eigvec)

   return
end subroutine landau_level_k

!> calculate landau levels in 3D system in kplane
!> fix B field
!> the magnetic field is in the a1 a2 plane
!> B= (B cos\theta, B sin\theta, 0)
!> Periodic Landau gauge by Guan Yifei. contact: guan.yifei@outlook.com
!> construct on Dec 8 2015
!> By QuanSheng Wu at ETH Zurich Honggerberg
!> License : GPL V3 
!>
!include "sparse.f90"
subroutine landau_level_kplane
   use para
   implicit none


   !> magnetic supercell size, perpendicular to the magnetic field, Ndimq= Nq*Num_wann
   integer :: Nq, Ndimq

   integer :: ik, i, j, ib, ierr, iq, ig

   !> lowest and highest band index to be calculated
   integer :: il, iu

   !> mdimq= iu-il+1, number of LLs to be calculated
   integer :: mdimq

   integer :: nkx, nky

   real(dp) :: B0, B0Tesla, B0Tesla_quantumflux_magsupcell
   real(dp) :: time_start, time_end, time_start0

   real(dp) :: theta
   real(dp) :: emin, emax

   ! wave vector
   real(dp) :: k3(3)
   real(dp), allocatable :: kxy(:,:), kxy_shape(:,:), kxy_plane(:,:)

   !> dim= Ndimq, knv3
   real(dp), allocatable :: W(:)
   real(dp), allocatable :: W_full(:)
   real(dp), allocatable :: eigv(:, :)
   real(dp), allocatable :: eigv_mpi(:, :)
   complex(dp), allocatable :: psi(:)

   !> print the weight for the Selected_WannierOrbitals
   real(dp), allocatable :: dos_selected(:, :, :)
   real(dp), allocatable :: dos_selected_mpi(:, :, :)

   !> dim= Ndimq*Ndimq
   complex(dp), allocatable :: eigvec(:, :)
   complex(dp), allocatable :: ham_landau(:, :)

   Nq= Magq
   Ndimq= Num_wann* Nq
   nkx= Nk
   nky= Nk

   if (cpuid==0) then
      write(stdout, '(a,i8)')' >> The magnetic-supercell size is ', Nq
      write(stdout, '(a,i8)')' >> The dimension of magnetic-supercell is ', Ndimq
      write(stdout, '(a,1000i5)')' >> calculated LLs for selected bands',   Selected_band_index(:)
      write(stdout, '(a,i8,a,i8)')' >> in total, there are  ', NumberofSelectedBands, ' bands'
   endif

   allocate( ham_landau(Ndimq, Ndimq))
   allocate( W( NumberofSelectedBands))
   allocate( W_full( Ndimq))
   allocate( eigv( NumberofSelectedBands, nkx*nky))
   allocate( eigv_mpi( NumberofSelectedBands, nkx*nky))
   allocate( eigvec( Ndimq, NumberofSelectedBands))
   allocate( psi(Ndimq))
   allocate( dos_selected     (NumberofSelectedBands,   nkx*nky, NumberofSelectedOrbitals_groups))
   allocate( dos_selected_mpi (NumberofSelectedBands,   nkx*nky, NumberofSelectedOrbitals_groups))
   dos_selected= 0d0
   dos_selected_mpi= 0d0
   eigv_mpi   = 0d0
   eigv       = 0d0
   eigvec     = 0d0
   ham_landau = 0d0


   !> deal with the magnetic field
   !> first transform the Bx By into B*Cos\theta, B*Sin\theta
   if (abs(By)<1e-8) then
      if (Bx<0) then
         theta= pi
      else
         theta= 0d0
      endif
   elseif (By>0) then
      theta = atan(Bx/By)
   else
      theta = atan(Bx/By)+ pi
   endif


   !> The flux in the supercell should be 2*pi

   !> the flux in the unit cell not the magnetic supercell
   theta=0d0
   B0= 2d0*pi/dble(Nq)*Magp
   Bx= B0* dcos(theta)
   By= B0* dsin(theta)


   !> transform it into Tesla
   !> B=2*pi*\phi_0/S0, where \phi_0 is the quantum flux, S0 is the projected area
   !> \phi_0 = h/2e  h=6.62607004*1E-34, e= 1.6*1E-19
   !> B0Tesla_quantumflux_magsupcell is the magnetic field that makes the flux through the magnetic
   !> supercell to be one quantum flux
   B0Tesla_quantumflux_magsupcell= 6.62607004*1E-34/1.6021766208*1E19/MagneticSuperProjectedArea*1E20
   B0Tesla= 6.62607004*1E-34/1.6021766208*1E19/MagneticSuperProjectedArea*1E20*Magp

   if (cpuid==0) then
      write(stdout, '(a, 2E18.8)')' Magnetic field B in Tesla that makes the flux through '
      write(stdout, '(a, 2E18.8)')' the magnetic supercell to be one quantum flux : ', B0Tesla_quantumflux_magsupcell
      write(stdout, '(a, 2E18.8)')' Magnetic field B in Tesla that fits to Magp : ', B0Tesla
      write(stdout, '(a, 2f18.8)')' Magnetic field Bx, By= ', Bx, By
   endif

   !> set k slice
   allocate( kxy(3, nkx*nky))
   allocate( kxy_shape(3, nkx*nky))
   allocate( kxy_plane(3, nkx*nky))
   kxy=0d0
   kxy_shape=0d0
   kxy_plane=0d0


   ik =0
   do i= 1, nkx
      do j= 1, nky
         ik =ik +1
         kxy(:, ik)= K3D_start+ K3D_vec1*(i-1)/dble(nkx-1)+ K3D_vec2*(j-1)/dble(nky-1) &
            -(K3D_vec1+K3D_vec2)/2d0
         kxy_shape(:, ik)= kxy(1, ik)* Kua_mag+ kxy(2, ik)* Kub_mag+ kxy(3, ik)* Kuc_mag
         call rotate_k3_to_kplane_mag(kxy_shape(:, ik), kxy_plane(:, ik))
      enddo
   enddo

   !> calculate the landau levels along special k line
   k3= 0
   time_start= 0d0
   time_start0= 0d0
   call now(time_start0)
   time_start= time_start0
   time_end  = time_start0
   do ik=1+ cpuid, nkx*nky, num_cpu
      if (cpuid.eq.0) &
         write(stdout, '(2a, i9, "  /", i10, a, f10.1, "s", a, f10.1, "s")') &
         'In LandauLevel_k_dos_Lanczos ', ' ik/NK ', ik, Nkx*nky, &
         ' time elapsed: ', time_end-time_start0, &
         ' time left: ', ((nkx*nky-ik)/num_cpu)*(time_end-time_start)

      call now(time_start)

      k3= kxy(:, ik)
      call ham_3Dlandau(Ndimq, Nq, k3, ham_landau)
      ham_landau= ham_landau/eV2Hartree

      !> diagonalization by call zheev in lapack
      !W_full= 0d0
      if (Ndimq<2000) then
         if (LandauLevel_wavefunction_calc) then
            call eigensystem_c( 'V', 'U', Ndimq ,ham_landau, W_full)
            do ib=1, NumberofSelectedBands
               eigv(ib, ik)= W_full(Selected_band_index(ib))
               eigvec(1:Ndimq, ib)= ham_landau(1:Ndimq, Selected_band_index(ib))
            enddo
         else
            call eigensystem_c( 'N', 'U', Ndimq ,ham_landau, W_full)
            do ib=1, NumberofSelectedBands
               eigv(ib, ik)= W_full(Selected_band_index(ib))
               eigvec(1:Ndimq, ib)= 0d0
            enddo
         endif
      else
         do i =1, NumberofSelectedBands
            il= Selected_band_index(i)
            iu= Selected_band_index(i)
            W= 0d0
            if (LandauLevel_wavefunction_calc) then
               call zheevx_pack('V', 'U', Ndimq, il, iu, ham_landau, W(1), eigvec(:, i))
            else
               call zheevx_pack('N', 'U', Ndimq, il, iu, ham_landau, W(1), eigvec(:, i))
            endif
            eigv(i, ik)= W(1)
         enddo
      endif

      if (LandauLevel_wavefunction_calc) then
         !> calculate the weight on the selected orbitals
         do ib= 1, NumberofSelectedBands
            psi(:)= eigvec(:, ib)  !> the eigenvector of ib'th band
            do ig=1, NumberofSelectedOrbitals_groups
               do iq=1, Nq
                  do i= 1, NumberofSelectedOrbitals(i)
                     j= Num_wann*(iq-1)+ Selected_WannierOrbitals(ig)%iarray(i)
                     dos_selected(ib, ik, ig)= dos_selected(ib, ik, ig)+ abs(psi(j))**2
                  enddo ! sweep the selected orbitals
               enddo ! iq sweep the magnetic supercell
            enddo ! ig
         enddo ! ib sweep the eigenvalue
      endif

      call now(time_end)
   enddo !ik

#if defined (MPI)
   call mpi_allreduce(eigv, eigv_mpi, size(eigv), &
      mpi_dp, mpi_sum, mpi_cmw, ierr)
   call mpi_allreduce(dos_selected, dos_selected_mpi,size(dos_selected),&
      mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
   eigv_mpi= eigv
   dos_selected_mpi= dos_selected
#endif

   !> minimum and maximum value of energy bands
   emin= minval(eigv_mpi)-0.5d0
   emax= maxval(eigv_mpi)+0.5d0

   !> write out the data in gnuplot format
   outfileindex= outfileindex+ 1
   if (cpuid==0)then
      open(unit=outfileindex, file='landaulevel_kplane.dat')
      write(outfileindex, "(6a19, ' E', i17, 1000i19 )")'# kx', 'ky', 'kz', 'k1', 'k2', 'k3', &
         Selected_band_index(:)
      do ik=1, nkx*nky
         write(outfileindex, '(10000f19.9)')kxy_shape(:, ik)*Angstrom2atomic, &
            kxy_plane(:, ik)*Angstrom2atomic, eigv_mpi(:, ik)
         if (mod(ik, nk2)==0) write(outfileindex, *)' '
      enddo
      close(outfileindex)
   endif

   outfileindex= outfileindex+ 1
   if (cpuid==0)then
      open(unit=outfileindex, file='landaulevel_kplane-matlab.txt')
      write(outfileindex, "(6a19, ' E', i17, 1000i19 )")'% kx', 'ky', 'kz', 'k1', 'k2', 'k3', &
         Selected_band_index(:)
      do ik=1, nkx*nky
         write(outfileindex, '(10000f19.9)')kxy_shape(:, ik)*Angstrom2atomic, &
            kxy_plane(:, ik)*Angstrom2atomic, eigv_mpi(:, ik)
      enddo
      close(outfileindex)
   endif




   !> write out a script that can be used for gnuplot
   outfileindex= outfileindex+ 1
   if (cpuid==0)then
      open(unit=outfileindex, file='landaulevel_kplane.gnu')
      write(outfileindex, '(a)')"set encoding iso_8859_1"
      write(outfileindex, '(a)')'#set terminal  postscript enhanced color'
      write(outfileindex, '(a)')"#set output 'bulkek_plane.eps'"
      write(outfileindex, '(3a)')'set terminal  png truecolor enhanced', &
         ' size 1920, 1680 font ",36"'
      write(outfileindex, '(a)')"set output 'landaulevel_kplane.png'"
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
      write(outfileindex, '(2a)')"splot 'landaulevel_kplane.dat' u 4:5:7 w pm3d, \"
      write(outfileindex, '(2a)')"      'landaulevel_kplane.dat' u 4:5:8 w pm3d"

      close(outfileindex)

   endif ! cpuid


   deallocate( ham_landau)
   deallocate( W)
   deallocate( W_full)
   deallocate( eigv)
   deallocate( eigv_mpi)
   deallocate( eigvec)

   return
end subroutine landau_level_kplane



!> calculate landau levels in 3D system for different B
!> fix k point, usually Gamma point
subroutine landau_level_B
   use para
   implicit none

   !> magnetic supercell size, perpendicular to the magnetic field
   !> Ndimq= Nq* Num_wann
   integer :: Nq, Ndimq
   integer :: Nmag

   integer :: ib, i, j, ie, ierr, iq, ig, ieta

   real(dp) :: B0, B0Tesla, B0Tesla_quantumflux_magsupcell, x, eta
   real(dp) :: time_start, time_end, time_start0

   ! wave vector
   real(dp) :: k3(3)
   real(dp), external :: delta

   !> dim= Ndimq, knv3
   real(dp), allocatable :: mag(:), mag_Tesla(:)
   real(dp), allocatable :: W(:)
   real(dp), allocatable :: eigv(:, :)
   real(dp), allocatable :: eigv_mpi(:, :)
   real(dp), allocatable :: eigv_mpi2(:, :)

   !> wave function
   complex(dp), allocatable :: psi(:)

   !> print the weight for the Selected_WannierOrbitals
   real(dp), allocatable :: dos_selected(:, :, :)
   real(dp), allocatable :: dos_selected_mpi(:, :, :)

   !> energy interval
   !> OmegaNum is defined in the module.f90 and read from the input.dat or wt.in
   real, allocatable :: omega(:)

   !> spectrum calculated
   real(dp), allocatable :: dos_B_omega(:, :, :), dos_B_omega_mpi(:, :, :)
   real(dp), allocatable ::  n_int(:, :)

   integer :: NumberofEta, ie_Earc
   real(dp), allocatable :: eta_array(:), n_Earc(:)



   !> dim= Ndimq*Ndimq
   complex(dp), allocatable :: ham_landau(:, :)

   Nq= nslab
   Nmag= Magp_max-Magp_min-1
   if (Nmag<=0) Nmag=Magp
   Ndimq= Num_wann* Nq
   allocate( ham_landau(Ndimq, Ndimq))
   allocate( W( Ndimq))
   allocate( eigv( Ndimq, Nmag+1))
   allocate( eigv_mpi( Ndimq, Nmag+1))
   allocate( eigv_mpi2( Ndimq, Nmag+1))
   allocate( mag(Nmag+1), mag_Tesla(Nmag+1))
   allocate( psi(Ndimq))
   allocate( dos_selected     (Ndimq,   Nmag+ 1, NumberofSelectedOrbitals_groups))
   allocate( dos_selected_mpi (Ndimq,   Nmag+ 1, NumberofSelectedOrbitals_groups))
   dos_selected= 0d0
   dos_selected_mpi= 0d0
   mag= 0d0
   eigv_mpi= 0d0
   eigv    = 0d0
   ham_landau= 0d0

   NumberofEta=9
   allocate(eta_array(NumberofEta))
   allocate(n_Earc(NumberofEta))
   eta_array=(/0.1d0, 0.2d0, 0.4d0, 0.8d0, 1.0d0, 2d0, 4d0, 8d0, 10d0/)
   eta_array= eta_array*Fermi_broadening

   allocate(omega(OmegaNum))
   allocate(n_int(0:Omeganum, NumberofEta))
   allocate(dos_B_omega(Nmag, OmegaNum, NumberofEta), dos_B_omega_mpi(Nmag, OmegaNum, NumberofEta))
   dos_B_omega= 0d0; dos_B_omega_mpi= 0d0

   !> energy
   do ie=1, OmegaNum
      omega(ie)= OmegaMin+ (OmegaMax-OmegaMin)* (ie-1d0)/dble(OmegaNum-1)
   enddo ! ie

   !> deal with the magnetic field
   !> first transform the Bx By into B*Cos\theta, B*Sin\theta
!~    if (abs(By)<1e-8) then
!~       if (Bx<0) then
!~          theta= pi
!~       else
!~          theta= 0d0
!~       endif
!~    elseif (By>0) then
!~       theta = atan(Bx/By)
!~    else
!~       theta = atan(Bx/By)+ pi
!~    endif

   !> The flux in the supercell should be the integer of 2*pi
   !    if (dis1< 1e-9) stop 'something wrong with the atom position'
   B0=2d0*pi/(Nq)
   B0=abs(B0)

   !> transform it into Tesla
   !> B=2*pi*\phi_0/S0, where \phi_0 is the quantum flux, S0 is the projected area
   !> \phi_0 = h/2e  h=6.62607004*1E-34, e= 1.6*1E-19
   !> B0Tesla_quantumflux_magsupcell is the magnetic field that makes the flux through the magnetic
   !> supercell to be one quantum flux
   B0Tesla_quantumflux_magsupcell= 6.62607004*1E-34/1.6021766208*1E19/MagneticSuperProjectedArea*1E20
   B0Tesla= 6.62607004*1E-34/1.6021766208*1E19/MagneticSuperProjectedArea*1E20*Magp

   if (cpuid==0) then
      write(stdout, '(a, 2E18.8)')' Magnetic field B in Tesla that makes the flux through '
      write(stdout, '(a, 2E18.8)')' the magnetic supercell to be one quantum flux : ', B0Tesla_quantumflux_magsupcell
      write(stdout, '(a, 2E18.8)')' Magnetic field B in Tesla that fits to Magp : ', B0Tesla
      write(stdout, '(a, 2f18.8)')' Magnetic field Bx, By= ', Bx, By
   endif

   k3=Single_KPOINT_3D_DIRECT
   
   if (cpuid==0) then
      write(stdout, *) k3
   endif

   !> calculate the landau levels along special k line
   eigv_mpi2=0
   do ib=1, Nmag+1
      mag(ib)= B0* (ib+Magp_min-1) 
      mag_Tesla(ib)= B0Tesla_quantumflux_magsupcell* (ib-1)
   enddo

   time_start= 0d0
   time_start0= 0d0
   call now(time_start0)
   time_start= time_start0
   time_end  = time_start0
   do ib=1+ cpuid, Nmag+1, num_cpu

      !Bx= mag(ib)* Cos(theta)
      !By= mag(ib)* Sin(theta)
      Bx= -mag(ib)
      By=0d0
      if (cpuid.eq.0) &
         write(stdout, '(2a, i9, "  /", i10, a, f10.1, "s", a, f10.1, "s")') &
         'In LandauLevel_B ', ' ib/Nmag ', ib, Nmag+1, &
         ' time elapsed: ', time_end-time_start0, &
         ' time left: ', ((Nmag-ib)/num_cpu)*(time_end-time_start)

      call now(time_start)


      call ham_3Dlandau(Ndimq, Nq, k3, ham_landau)

      !> diagonalization by call zheev in lapack
      W= 0d0

      if (LandauLevel_wavefunction_calc) then
         call eigensystem_c( 'V', 'U', Ndimq ,ham_landau, W)
         !> calculate the weight on the selected orbitals
         do ie= 1, Ndimq
            psi(:)= ham_landau(:, ie)  !> the eigenvector of ib'th band
            do ig=1, NumberofSelectedOrbitals_groups
               do iq=1, Nq
                  do i= 1, NumberofSelectedOrbitals(ig)
                     j= Num_wann*(iq-1)+ Selected_WannierOrbitals(ig)%iarray(i)
                     dos_selected(ie, ib, ig)= dos_selected(ie, ib, ig)+ abs(psi(j))**2
                  enddo ! sweep the selected orbitals
               enddo ! iq sweep the magnetic supercell
            enddo ! ig
         enddo ! ib sweep the eigenvalue

      else
         call eigensystem_c( 'N', 'U', Ndimq ,ham_landau, W)
      endif

      eigv(:, ib)= W

      !> calculate density of states using the eigenvalues
      do ieta=1, NumberofEta
         eta= eta_array(ieta)
         do ie=1, OmegaNum
            do i = 1, Ndimq
               x= omega(ie)- W(i)
               dos_B_omega(ib, ie, ieta)= dos_B_omega(ib, ie, ieta)+ &
                  delta(eta, x)
            enddo
         enddo
      enddo
    
      call now(time_end)
   enddo !ib

#if defined (MPI)
   call mpi_allreduce(eigv,eigv_mpi,size(eigv),&
      mpi_dp,mpi_sum,mpi_cmw,ierr)
   call mpi_allreduce(dos_selected, dos_selected_mpi,size(dos_selected),&
      mpi_dp,mpi_sum,mpi_cmw,ierr)
   call mpi_allreduce(dos_B_omega, dos_B_omega_mpi, size(dos_B_omega), &
      mpi_dp, mpi_sum, mpi_cmw, ierr)
#else
   eigv_mpi= eigv
   dos_selected_mpi= dos_selected
   dos_B_omega_mpi= dos_B_omega
#endif
   eigv_mpi = eigv_mpi/eV2Hartree

   outfileindex= outfileindex+ 1
   if (cpuid.eq.0) then
      open(unit=outfileindex, file='landaulevel_B.dat')
      if (LandauLevel_wavefunction_calc) then
         write(outfileindex, '("#", a, 3f12.6, a)')'k points in fractional coordinates (', k3, ')'
         write(outfileindex, '("#", a14, 2a15, a)')'Phi per cell', 'B (Tesla)', ' Eig of LL', ' Weight on the selected orbitals'
         do j=1, Ndimq
            do ib=1, Nmag+1
               write(outfileindex, '(30f19.8)')mag(ib)/2d0/pi, mag_Tesla(ib), eigv_mpi(j, ib),  &
                  (dos_selected_mpi(j, ib, ig), ig=1, NumberofSelectedOrbitals_groups)
            enddo
            write(outfileindex , *)' '
         enddo
      else
         write(outfileindex, '("#", a, 3f12.6, a)')'k points in fractional coordinates (', k3, ')'
         write(outfileindex, '("#", a14, 2a15)')'Phi per cell', 'B (Tesla)', ' Eig of LL'
         do j=1, Ndimq
            do ib=1, Nmag+1
               write(outfileindex, '(30f19.8)')mag(ib)/2d0/pi, mag_Tesla(ib), eigv_mpi(j, ib)
            enddo
            write(outfileindex , *)' '
         enddo

      endif
      close(outfileindex)
      write(stdout,*) 'calculate landau level done'
   endif

   outfileindex= outfileindex+ 1
   if (cpuid==0) then
      open(unit=outfileindex, file='landaulevel_B.gnu')
      write(outfileindex, '(a)')"set encoding iso_8859_1"
      write(outfileindex, '(a)')'#set terminal  postscript enhanced color'
      write(outfileindex, '(a)')"#set output 'landaulevel_B.eps'"
      write(outfileindex, '(3a)')'#set terminal  pngcairo truecolor enhanced', &
         '  font ",60" size 1920, 1680'
      write(outfileindex, '(3a)')'set terminal  png truecolor enhanced', &
         '  font ",60" size 1920, 1680'
      write(outfileindex, '(a)')"set output 'landaulevel_B.png'"
      write(outfileindex, '(a)')'unset key'
      write(outfileindex, '(a)')'set pointsize 0.8'
      write(outfileindex, '(a)')'set border lw 3 '
      write(outfileindex, '(a)')'#set xtics font ",36"'
      write(outfileindex, '(a)')'#set ytics font ",36"'
      write(outfileindex, '(a)')'#set xlabel font ",36"'
      write(outfileindex, '(a)')'set xlabel "Phi per cell"'
      write(outfileindex,*) '#set xlabel "{/Symbol F}/{/Symbol F}_0"'
      write(outfileindex, '(a)')'set ylabel "Energy (eV)"'
      write(outfileindex, '(a, i6,a)')'set title "Hofstadter butterfly with Nq=', Nq, '" font ",40"'
      write(outfileindex, '(a)')'#set xtics offset 0, -1'
      write(outfileindex, '(a)')'#set ylabel offset -1, 0 '
      write(outfileindex, '(a)')'rgb(r,g,b) = int(r)*65536 + int(g)*256 + int(b)'
      if (LandauLevel_wavefunction_calc) then
         write(outfileindex, '(2a)')"plot 'landaulevel_B.dat' u 1:3:(rgb(255,255-255*$4, 4)) ",  &
            "w p  pt 7  ps 1 lc rgb variable"
      else
         write(outfileindex, '(2a)')"plot 'landaulevel_B.dat' u 1:3",  &
            " w p  pt 7  ps 1"
      endif
      close(outfileindex)
   endif

   outfileindex= outfileindex+ 1
   if (cpuid.eq.0)then
      open (unit=outfileindex, file='LandauLevel_B_dos.dat')
      open (unit=outfileindex+1,file='wannierdiagram.dat')
      write(outfileindex, '("#", a)')'Hofstadter butterfly (Landau level in (E,B) axis) '
      write(outfileindex, '("#", a, i6)')'Magnetic supercell size : ', Magq
      write(outfileindex+1, '("#", a)')'Wannier diagram '
      write(outfileindex+1, '("#", a, i6)')'Magnetic supercell size : ', Magq
      write(outfileindex, '("#", a, 3f12.6, a)')'k points in fractional coordinates (', k3, ')'
      write(outfileindex, '("# Column", I5, 100I16)')(i, i=1, 12)
      write(outfileindex, '("#", a15, 5a16)', advance='NO')'Flux', 'B(Tesla)', 'E(eV)'
      write(outfileindex+1, '("# Column", I5, 100I16)')(i, i=1, 20)
      write(outfileindex+1, '("#", a15, 2a16)', advance='NO')'Flux', 'B(Tesla)'
      do ieta=1, NumberofEta-1
         write(outfileindex, '(a16)', advance='NO') 'A(B, E)'
         write(outfileindex+1, '(2a16)', advance='NO') 'n   ', 'A(n, E)'
      enddo
      write(outfileindex, '(a16)') 'A(B, E)'
      write(outfileindex+1, '(2a16)') 'n   ', 'A(n, E)'

      write(outfileindex, '("#", a, 20X, 300f16.2)')'Broadening \eta (meV): ', Eta_array(:)*1000d0/eV2Hartree
      write(outfileindex+1, '("#", a, 5X, f26.2,300f32.2)')'Broadening \eta (meV): ', Eta_array(:)*1000d0/eV2Hartree

      !> find n_int at EF
      do ie=1, omeganum
         if (omega(ie)>iso_energy) then
            ie_Earc= ie- 1
            exit
         endif
      enddo

      do ib=1, Nmag
         n_int=0
         do ie=1, omeganum
            write(outfileindex, '(300E16.4)')mag(ib)/2d0/pi, mag_Tesla(ib), omega(ie)/eV2Hartree, dos_B_omega_mpi(ib, ie, :)
         enddo

         !> get number of electrons between the lowest energy level and iso_energy
         n_Earc= 0d0
         do ie=1, ie_Earc
            n_Earc(:)= n_Earc(:)+ dos_B_omega_mpi(ib, ie, :)
         enddo

         !> get number of electrons between the lowest energy level and omega(ie)
         do ie=1, omeganum
            n_int(ie, :)=n_int(ie-1, :)+ dos_B_omega_mpi(ib, ie, :)
         enddo

         !> set n(E)= n_int- n_Earc= \int_Earc^E \rho(\epsilon)d\epsilon
         do ie=1, omeganum
            n_int(ie, :)= n_int(ie, :)-n_Earc(:)
         enddo

        !do ieta=1, NumberofEta
        !   n_int(:, ieta)=n_int(:, ieta)/n_int(omeganum, ieta)
        !enddo

         do ie=1, omeganum
            write(outfileindex+1,'(300E16.4)')  mag(ib)/2d0/pi, mag_Tesla(ib), &
               (n_int(ie, ieta), dos_B_omega_mpi(ib, ie, ieta), ieta=1, NumberofEta)
         enddo

         write(outfileindex+1, *) ' '
         write(outfileindex, *) ' '
      enddo
      close(outfileindex)
      write(stdout,*)'calculate Landau level spectrum in B-E mode successfully'
   endif

   outfileindex= outfileindex+ 2
   if (cpuid.eq.0)then
      open (unit=outfileindex, file='LandauLevel_B_dos.gnu')
      write(outfileindex, '(a)')'#set terminal  postscript enhanced color font ",24"'
      write(outfileindex, '(a)')'set terminal pngcairo enhanced color font ",60" size 1920, 1680'
      write(outfileindex, '(a)')"set output 'LandauLevel_B_dos.png'"
      write(outfileindex, '(a)')'set pm3d'
      write(outfileindex, '(a)')'set palette rgb 21,22,23'
      write(outfileindex, '(a)')'#set isosamples 50,50'
      write(outfileindex, '(a)')'set size 0.9,1'
      write(outfileindex, '(a)')'set origin 0.05,0'
      write(outfileindex, '(a)')'set view map'
      write(outfileindex, '(a)')'unset ztics'
      write(outfileindex, '(a)')'unset surface'
      write(outfileindex, '(a)')'unset key'
      write(outfileindex, '(a)')'#set xtics font ",24"'
      write(outfileindex, '(a)')'#set ytics font ",24"'
      write(outfileindex, '(a)')'#set xlabel font ",24"'
      write(outfileindex, '(a)')'set ylabel "Energy (eV)"'
      write(outfileindex, '(a, i6,a)')'set title "Hofstadter butterfly with Nq=', Nq, '" font ",40"'
      write(outfileindex, '(a, f12.6, a, f12.6, a)')'set yrange [',OmegaMin/eV2Hartree,':',OmegaMax/eV2Hartree,']'
      write(outfileindex, '(a)')'set xlabel "Phi per unit cell"'
      write(outfileindex,*) '#set xlabel "{/Symbol F}/{/Symbol F}_0"'
      write(outfileindex, '(a, f12.6, a)')     'set xrange [ 0.0000 :', maxval(mag/2d0/pi) ,']'
      write(outfileindex, '(a)')     "splot 'LandauLevel_B_dos.dat' u 1:3:(log($4)) w pm3d"
      write(outfileindex, '(a)')'set xlabel "B (Tesla)"'
      write(outfileindex, '(a, f12.6, a)')     '#set xrange [ 0.0000 :', maxval(mag_Tesla) ,']'
      write(outfileindex, '(a)')     "#splot 'LandauLevel_B_dos.dat' u 2:3:(log($4)) w pm3d"
   endif

   outfileindex=outfileindex+1
   if(cpuid == 0) then
      open(unit=outfileindex, file='wannierdiagram.gnu')
      write(outfileindex,*) 'set terminal pngcairo enhanced color font ",60" size 1920, 1680'
      write(outfileindex,*) "set output 'wannierdiagram.png'"
      write(outfileindex,*) 'set pm3d'
      write(outfileindex, '(a)')'set size 0.9,1'
      write(outfileindex, '(a)')'set origin 0.05,0'
      write(outfileindex, '(a)')'set palette rgb 21,22,23'
      write(outfileindex,*) '#set isosamples 50,50'
      write(outfileindex,*) 'set view map'
      write(outfileindex,*) 'unset ztics'
      write(outfileindex,*) 'unset surface'
      write(outfileindex,*) 'unset key'

      write(outfileindex,*) 'set ylabel "n"'
      write(outfileindex, '(a, i6, a)') 'set title "Wannier diagram with Nq=', Nq, '"'
      write(outfileindex,*) '#set yrange [   ] noextend'
      write(outfileindex,*) '#set xlabel "{/Symbol F}/{/Symbol F}_0"'
      write(outfileindex,*) 'set xlabel "Phi/Phi_0 per unit cell"'
      write(outfileindex,*) '#set xlabel "{/Symbol F}/{/Symbol F}_0"'
      write(outfileindex,*) '#set xrange [ ] noextend'

      write(outfileindex,*) "splot 'wannierdiagram.dat' u 1:3:(log($4)) w pm3d #lc palette"
   end if


#if defined (MPI)
   call mpi_barrier(mpi_cmw, ierr)
#endif

   deallocate( ham_landau)
   deallocate( W)
   deallocate( eigv)
   deallocate( eigv_mpi)
   deallocate( eigv_mpi2)
   deallocate( mag)

   return
end subroutine landau_level_B



subroutine landau_sf
   use para
   implicit none


   !> magnetic supercell size, perpendicular to the magnetic field, Ndimq= Nq*Num_wann
   integer :: Nq, Ndimq

   integer :: ik, i, j, ib, ierr, iq, ig

   !> lowest and highest band index to be calculated
   integer :: il, iu

   !> mdimq= iu-il+1, number of LLs to be calculated
   integer :: mdimq


   real(dp) :: B0, B0Tesla, B0Tesla_quantumflux_magsupcell
   real(dp) :: time_start, time_end, time_start0

   real(dp) :: theta
   real(dp) :: emin, emax

   ! wave vector
   real(dp) :: k3(3)

   !> dim= Ndimq, knv3
   real(dp), allocatable :: W(:)
   real(dp), allocatable :: W_full(:)
   real(dp), allocatable :: eigv(:, :)
   real(dp), allocatable :: eigv_mpi(:, :)
   complex(dp), allocatable :: psi(:)

   !> print the weight for the Selected_WannierOrbitals
   real(dp), allocatable :: dos_selected(:, :, :)
   real(dp), allocatable :: dos_selected_mpi(:, :, :)

   !> dim= Ndimq*Ndimq
   complex(dp), allocatable :: eigvec(:, :)
   complex(dp), allocatable :: ham_landau(:, :)

   Nq= Magq
   Ndimq= Num_wann* Nq

   il= Numoccupied*Nq- 2000
   iu= Numoccupied*Nq+ 2000
   il=1
   iu=Ndimq
   if (il< 1) il= 1
   if (iu> Ndimq) iu= Ndimq
   mdimq= iu-il+1

   if (cpuid==0) then
      write(stdout, '(a,i8)')' >> The magnetic-supercell size is ', Nq
      write(stdout, '(a,i8)')' >> The dimension of magnetic-supercell is ', Ndimq
      write(stdout, '(a,i8,a,i8)')' >> calculated LLs from il', il, ' to', iu, ' bands'
      write(stdout, '(a,i8,a,i8)')' >> in total, there are  ', mdimq, ' bands'
   endif

   allocate( ham_landau(Ndimq, Ndimq))
   allocate( W( mdimq))
   allocate( W_full( Ndimq))
   allocate( eigv( mdimq, knv2))
   allocate( eigv_mpi( mdimq, knv2))
   allocate( eigvec( Ndimq, mdimq))
   allocate( psi(Ndimq))
   allocate( dos_selected     (mdimq,   knv2, NumberofSelectedOrbitals_groups))
   allocate( dos_selected_mpi (mdimq,   knv2, NumberofSelectedOrbitals_groups))
   dos_selected= 0d0
   dos_selected_mpi= 0d0
   eigv_mpi   = 0d0
   eigv       = 0d0
   eigvec     = 0d0
   ham_landau = 0d0


   !> deal with the magnetic field
   !> first transform the Bx By into B*Cos\theta, B*Sin\theta
   if (abs(By)<1e-8) then
      if (Bx<0) then
         theta= pi
      else
         theta= 0d0
      endif
   elseif (By>0) then
      theta = atan(Bx/By)
   else
      theta = atan(Bx/By)+ pi
   endif


   !> The flux in the supercell should be 2*pi
   !    if (dis1< 1e-9) stop 'something wrong with the atom position'

   !> the flux in the unit cell not the magnetic supercell
   B0= 2d0*pi/dble(Nq)*Magp
   Bx= B0* dcos(theta)
   By= B0* dsin(theta)


   !> transform it into Tesla
   !> B=2*pi*\phi_0/S0, where \phi_0 is the quantum flux, S0 is the projected area
   !> \phi_0 = h/2e  h=6.62607004*1E-34, e= 1.6*1E-19
   !> B0Tesla_quantumflux_magsupcell is the magnetic field that makes the flux through the magnetic
   !> supercell to be one quantum flux
   B0Tesla_quantumflux_magsupcell= 6.62607004*1E-34/1.6021766208*1E19/MagneticSuperProjectedArea*1E20
   B0Tesla= 6.62607004*1E-34/1.6021766208*1E19/MagneticSuperProjectedArea*1E20*Magp

   if (cpuid==0) then
      write(stdout, '(a, 2E18.8)')' Magnetic field B in Tesla that makes the flux through '
      write(stdout, '(a, 2E18.8)')' the magnetic supercell to be one quantum flux : ', B0Tesla_quantumflux_magsupcell
      write(stdout, '(a, 2E18.8)')' Magnetic field B in Tesla that fits to Magp : ', B0Tesla
      write(stdout, '(a, 2f18.8)')' Magnetic field Bx, By= ', Bx, By
   endif

   !> calculate the landau levels along special k line
   k3= 0
   time_start= 0d0
   time_start0= 0d0
   call now(time_start0)
   time_start= time_start0
   time_end  = time_start0
   do ik=1+cpuid,knv2,num_cpu
      if (cpuid.eq.0) &
         write(stdout, '(2a, i9, "  /", i10, a, f10.1, "s", a, f10.1, "s")') &
         'In LandauLevel_k_dos_Lanczos ', ' ik/NK ', ik, knv2, &
         ' time elapsed: ', time_end-time_start0, &
         ' time left: ', ((knv2-ik)/num_cpu)*(time_end-time_start)

      call now(time_start)

      k3(2:3)= k2_path(ik,:)
      k3(1)=0d0
      call ham_3Dlandau(Ndimq, Nq, k3, ham_landau)
      ham_landau= ham_landau/eV2Hartree

      !> diagonalization by call zheev in lapack
      !W_full= 0d0
      if (mdimq<2000) then
         if (LandauLevel_wavefunction_calc) then
            call eigensystem_c( 'V', 'U', Ndimq ,ham_landau, W_full)
            eigv(:, ik)= W_full(il:iu)
            eigvec(1:Ndimq, 1:mdimq)= ham_landau(1:Ndimq, 1:mdimq)
         else
            call eigensystem_c( 'N', 'U', Ndimq ,ham_landau, W_full)
            eigv(:, ik)= W_full(il:iu)
            eigvec(1:Ndimq, 1:mdimq)= 0d0
         endif
      else
         W= 0d0
         if (LandauLevel_wavefunction_calc) then
            call zheevx_pack('V', 'U', Ndimq, il, iu, ham_landau, W, eigvec)
         else
            call zheevx_pack('N', 'U', Ndimq, il, iu, ham_landau, W, eigvec)
         endif
         eigv(:, ik)= W
      endif

      !> calculate the weight on the selected orbitals
      do ib= 1, mdimq
         psi(:)= eigvec(:, ib)  !> the eigenvector of ib'th band
         do ig=1, NumberofSelectedOrbitals_groups
            do iq=1, Nq
               do i= 1, NumberofSelectedOrbitals(ig)
                  j= Num_wann*(iq-1)+ Selected_WannierOrbitals(ig)%iarray(i)
                  dos_selected(ib, ik, ig)= dos_selected(ib, ik, ig)+ abs(psi(j))**2
               enddo ! sweep the selected orbitals
            enddo ! iq sweep the magnetic supercell
         enddo ! ig
      enddo ! ib sweep the eigenvalue

      call now(time_end)
   enddo !ik

#if defined (MPI)
   call mpi_allreduce(eigv, eigv_mpi, size(eigv), &
      mpi_dp, mpi_sum, mpi_cmw, ierr)
   call mpi_allreduce(dos_selected, dos_selected_mpi,size(dos_selected),&
      mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
   eigv_mpi= eigv
   dos_selected_mpi= dos_selected
#endif

   !> minimum and maximum value of energy bands
   emin= minval(eigv_mpi)-0.5d0
   emax= maxval(eigv_mpi)+0.5d0

   outfileindex= outfileindex+ 1
   if (cpuid.eq.0) then
      open(unit=outfileindex, file='landaulevel_k.dat')
      if (LandauLevel_wavefunction_calc) then
         write(outfileindex, '("#", a14, f16.8)')'Magnetic field in Tesla: ', B0Tesla
         write(outfileindex, '("#", a14, a15, a)')'k ', ' Eig of LL', ' Weight on the selected orbitals'
         do j=1, mdimq
            do i=1, knv2
               write(outfileindex,'(200f16.8)')k2len(i), eigv_mpi(j, i), &
                  (dos_selected_mpi(j, i, ig), ig=1, NumberofSelectedOrbitals_groups)
            enddo
            write(outfileindex , *)' '
         enddo
      else
         write(outfileindex, '("#", a14, f16.8)')'Magnetic field in Tesla: ', B0Tesla
         write(outfileindex, '("#", a14, a15, a)')'k ', ' Eig of LL', ' Weight on the selected orbitals'
         do j=1, mdimq
            do i=1, knv2
               write(outfileindex,'(20f16.8)')k2len(i), eigv_mpi(j, i)
            enddo
            write(outfileindex , *)' '
         enddo
      endif

      close(outfileindex)
      write(stdout,*) 'Calculate landau level done'
   endif

   outfileindex= outfileindex+ 1
   if (cpuid==0) then
      open(unit=outfileindex, file='landaulevel_k.gnu')
      write(outfileindex, '(a)')"set encoding iso_8859_1"
      write(outfileindex, '(a)')'#set terminal  postscript enhanced color'
      write(outfileindex, '(a)')"#set output 'landaulevel_k.eps'"
      write(outfileindex, '(3a)')'#set terminal  pngcairo truecolor enhanced', &
         '  font ",60" size 1920, 1680'
      write(outfileindex, '(3a)')'set terminal  png truecolor enhanced', &
         '  font ",60" size 1920, 1680'
      write(outfileindex, '(a)')"set output 'landaulevel_k.png'"
      write(outfileindex, '(a)')'unset key'
      write(outfileindex, '(a)')'set pointsize 0.8'
      write(outfileindex, '(a)')'set border lw 3 '
      write(outfileindex, '(a)')'#set xtics font ",36"'
      write(outfileindex, '(a)')'#set ytics font ",36"'
      write(outfileindex, '(a)')'#set xlabel font ",36"'
      write(outfileindex, '(a)')'set ylabel "Energy (eV)"'
      write(outfileindex, '(a, f10.5, a)')'set xrange [0: ', maxval(k2len), ']'
      write(outfileindex, '(a, i6, a, i6, a)')'set title "Landau level with p/q=', magp, "/", Nq,  '"'
      write(outfileindex, '(a)')'#set xtics offset 0, -1'
      write(outfileindex, '(a)')'#set ylabel offset -1, 0 '
      write(outfileindex, '(a)')'rgb(r,g,b) = int(r)*65536 + int(g)*256 + int(b)'
      write(outfileindex, 202, advance="no") (trim(k2line_name(i)), k2line_stop(i), i=1, nk2lines)
      write(outfileindex, 203)trim(k2line_name(nk2lines+1)), k2line_stop(nk2lines+1)

      write(outfileindex, '(a, f10.5, a, f10.5, a)')'set yrange [', emin, ':', emax, ']'

      do i=1, nk2lines-1
         write(outfileindex, 204)k2line_stop(i+1), emin, k2line_stop(i+1), emax
      enddo

202   format('set xtics (',:20('"',A3,'" ',F10.5,','))
203   format(A3,'" ',F10.5,')')
204   format('set arrow from ',F10.5,',',F10.5,' to ',F10.5,',',F10.5, ' nohead')

      if (LandauLevel_wavefunction_calc) then
         write(outfileindex, '(2a)')"plot 'landaulevel_k.dat' u 1:2:(rgb(255,255-255*$3, 3)) ",  &
            " w p  pt 7  ps 1 lc rgb variable"
      else
         write(outfileindex, '(2a)')"plot 'landaulevel_k.dat' u 1:2",  &
            " w p  pt 7  ps 1"
      endif
      close(outfileindex)
   endif


   deallocate( ham_landau)
   deallocate( W)
   deallocate( W_full)
   deallocate( eigv)
   deallocate( eigv_mpi)
   deallocate( eigvec)

   return
end subroutine landau_sf



!> get hybridization function
!> summation the g(i*\omega_n)=\sum_k G(k, i*\omega_n)
!>  G(k, i*\omega_n)= 1/(i*\omega_n+ mu- H(k))
!>  g^{-1}(i*\omega_n)= i*\omega_n+ mu - \Delta(i*\omega_n)
subroutine get_hybridization(omegan, mu, Nq, Ndimq, Delta)
   use para
   implicit none

   integer, intent(in) :: Ndimq, Nq
   real(dp), intent(in) :: omegan, mu
   complex(dp), intent(out) :: Delta(Ndimq, Ndimq)

   integer :: knv3, ik1, ik2, ik3, ierr, ik, i
   real(dp) :: k3(3)
   complex(dp), allocatable :: g0(:, :)
   complex(dp), allocatable :: g0_mpi(:, :)
   complex(dp), allocatable :: G_latt(:, :)
   complex(dp), allocatable :: ham_landau(:, :)
   complex(dp), allocatable :: z_unit(:, :)

   allocate(G_latt(Ndimq, Ndimq))
   allocate(z_unit(Ndimq, Ndimq))
   allocate(ham_landau(Ndimq, Ndimq))
   allocate(g0(Ndimq, Ndimq))
   allocate(g0_mpi(Ndimq, Ndimq))
   G_latt= 0d0; z_unit= 0d0; g0= 0d0; ham_landau= 0d0; g0_mpi= 0d0

   !> an unit matrix
   do i=1, Ndimq
      z_unit(i, i)= 1d0
   enddo

   knv3= Nk1*Nk2*Nk3
   do ik=1+cpuid, knv3, num_cpu
      ik1= (ik-1)/(Nk2*Nk3)+1
      ik2= ((ik-1-(ik1-1)*Nk2*Nk3)/Nk3)+1
      ik3= (ik-(ik2-1)*Nk3- (ik1-1)*Nk2*Nk3)
      k3(1)=(ik1-1)/dble(Nk1)
      k3(2)=(ik2-1)/dble(Nk2)
      k3(3)=(ik3-1)/dble(Nk3)

      !> Get the hamiltonian with given k point and size of magnetic supercell.
      !> The magnetic field is determined by the size of magnetic supercell Nq.
      call ham_3Dlandau(Ndimq, Nq, k3, ham_landau)

      !> Get the lattice green's function
      G_latt= (zi*omegan + mu)*z_unit - ham_landau

      !> Do a summation over G_latt
      g0_mpi= g0_mpi+ G_latt
   enddo ! ik
#if defined (MPI)
   call mpi_allreduce(g0_mpi, g0, size(g0), &
      mpi_dc, mpi_sum, mpi_cmw, ierr)
#else
   g0= g0_mpi
#endif

   !> inverse of g0
   call inv(Ndimq, g0)

   !> calculate the hybridization function
   Delta= (zi*omegan+ mu)*z_unit- g0

   !nwann= Num_wann/2
   !allocate( orbital_start(Origin_cell%Num_atoms+ 1))
   !orbital_start= 0
   !orbital_start(1)= 1
   !do ia=1, Origin_cell%Num_atoms
   !   orbital_start(ia+1)= orbital_start(ia)+ Origin_cell%nprojs(ia)
   !enddo


!   !> change the basis order from (o1, up) (o2, up) (o1, dn) (o2, dn) to
!   !> (o1, up) (o1, dn) (o2, up), (o2, dn)
!   !> old basis the first index
!   do iq1=1, Nq
!      do ia1=1, Origin_cell%Num_atoms
!         !> old basis the second index
!         do iq2=1, Nq
!            do ia2=1, Origin_cell%Num_atoms
!               do jq1=1, Nq
!                  do ja1=1, Origin_cell%Num_atoms
!                     !> old basis the second index
!                     do jq2=1, Nq
!                        do ja2=1, Origin_cell%Num_atoms
!                           istart1= (iq1-1)*Num_wann+ orbital_start(ia1)
!                           iend1= (iq1-1)*Num_wann+ orbital_start(ia1+1)- 1
!                           istart2= (iq2-1)*Num_wann+ orbital_start(ia2)
!                           iend2= (iq2-1)*Num_wann+ orbital_start(ia2+1)- 1
!                           jstart1= (jq1-1)*Num_wann+ orbital_start(ja1)
!                           jend1= (jq1-1)*Num_wann+ orbital_start(ja1+1)- 1
!                           jstart2= (jq2-1)*Num_wann+ orbital_start(ja2)
!                           jend2= (jq2-1)*Num_wann+ orbital_start(ja2+1)- 1
!                           ham_landau_newbasis(jstart1:jend1, jstart2:jend2)=  &
!                           ham_landau_newbasis(istart1:iend1, istart2:iend2)
!                           if (soc>0) then
!                              istart1= (iq1-1)*Num_wann+ orbital_start(ia1)+ Nwann
!                              iend1= (iq1-1)*Num_wann+ orbital_start(ia1+1)- 1+ Nwann
!                              istart2= (iq2-1)*Num_wann+ orbital_start(ia2)
!                              iend2= (iq2-1)*Num_wann+ orbital_start(ia2+1)- 1
!                              jstart1= (jq1-1)*Num_wann+ orbital_start(ja1)
!                              jend1= (jq1-1)*Num_wann+ orbital_start(ja1+1)- 1
!                              jstart2= (jq2-1)*Num_wann+ orbital_start(ja2)
!                              jend2= (jq2-1)*Num_wann+ orbital_start(ja2+1)- 1
!                              ham_landau_newbasis(jstart1:jend1, jstart2:jend2)=  &
!                              ham_landau_newbasis(istart1:iend1, istart2:iend2)
!                           endif
!                        enddo
!                     enddo
!                  enddo
!               enddo
!            enddo
!         enddo
!      enddo
!   enddo


   return

end subroutine get_hybridization


!> calculate hamiltonian for landau levels
!> consider the internal atom's position
subroutine ham_3Dlandau(Ndimq, Nq, k, ham_landau)
   use para
   implicit none

   integer, intent(in) :: Ndimq
   integer, intent(in) :: Nq
   real(dp), intent(in) :: k(3)
   complex(dp), intent(out) :: ham_landau(Ndimq, Ndimq)

   !> inta-hopping for the supercell
   complex(dp), allocatable :: H00(:, :)

   !> inter-hopping for the supercell
   complex(dp), allocatable :: H01(:, :)

   ! loop index
   integer :: i1, i2, iR

   ! index used to sign irvec
   integer :: ia1, ia2
   real(dp) :: ia, ib, ic

   integer :: istart1, istart2, iend1, iend2

   integer :: inew_ic

   !> nwann= Num_wann/2
   integer :: nwann

   integer, allocatable :: orbital_start(:)

   ! new index used to sign irvec
   real(dp) :: new_ia,new_ib,new_ic

   ! wave vector k times lattice vector R
   real(dp) :: kdotr, phase
   complex(dp) :: ratio, fac

   real(dp) :: Rp1(3), Rp2(3), R1(3), R2(3)
   real(dp) :: Ri(3), Rj(3), tau1(3), tau2(3)

   real(dp) :: hct,ndq
   real(dp),external :: phase2d,phase2

   allocate( H00( Ndimq, Ndimq))
   allocate( H01( Ndimq, Ndimq))

   nwann= Num_wann/2
   allocate( orbital_start(Origin_cell%Num_atoms+ 1))
   orbital_start= 0
   orbital_start(1)= 1
   do ia1=1, Origin_cell%Num_atoms
      orbital_start(ia1+1)= orbital_start(ia1)+ Origin_cell%nprojs(ia1)
   enddo
   !> calculate intra-hopping
   H00=0.0d0
   ! i1 column index
   do i1=1, Nq
      ! i2 row index
      do i2=1, Nq
         if (abs(i2-i1)> ijmax) cycle
         !> sum over R points to get H(k1, k2)
         do iR=1, Nrpts !Rth in lattice vector
            ia=dble(irvec(1,iR))
            ib=dble(irvec(2,iR))
            ic=dble(irvec(3,iR))

            !> new lattice
            call latticetransform(ia, ib, ic, new_ia, new_ib, new_ic)

            !> the magnetic supercell is along the third direction
            inew_ic= int(new_ic)
            if (inew_ic /= (i2-i1)) cycle

            !> exp(i k.R)
            kdotr= k(1)*new_ia+ k(2)*new_ib
            ratio= cos(2d0*pi*kdotr)+ zi*sin(2d0*pi*kdotr)

            R1= (i1-1)*Ruc_new
            R2= new_ia*Rua_new+ new_ib*Rub_new+ (i2-1)*Ruc_new
            Rp1= R1
            Rp2= R2

            do ia1=1, Origin_cell%Num_atoms
               do ia2=1, Origin_cell%Num_atoms
                  R1= Origin_cell%Atom_position_cart(:, ia1)
                  R2= Origin_cell%Atom_position_cart(:, ia2)
                  call rotate(R1, tau1)
                  call rotate(R2, tau2)

                  Ri= Rp1+ tau1
                  Rj= Rp2+ tau2

                  phase=phase2(Ri,Rj)
                  fac= cos(phase)+ zi*sin(phase)

                  istart1= (i1-1)*Num_wann+ orbital_start(ia1)
                  iend1= (i1-1)*Num_wann+ orbital_start(ia1+1)- 1
                  istart2= (i2-1)*Num_wann+ orbital_start(ia2)
                  iend2= (i2-1)*Num_wann+ orbital_start(ia2+1)- 1

                  H00( istart1:iend1, istart2:iend2) &
                     = H00( istart1:iend1, istart2:iend2) &
                     + HmnR( istart1- (i1-1)*Num_wann:iend1- (i1-1)*Num_wann, &
                     istart2- (i2-1)*Num_wann:iend2- (i2-1)*Num_wann, iR)*ratio/ndegen(iR)* fac

                  !> there is soc term in the hr file
                  if (soc>0) then
                     istart1= (i1-1)*Num_wann+ orbital_start(ia1) + Nwann
                     iend1= (i1-1)*Num_wann+ orbital_start(ia1+1)- 1 + Nwann
                     istart2= (i2-1)*Num_wann+ orbital_start(ia2)
                     iend2= (i2-1)*Num_wann+ orbital_start(ia2+1)- 1

                     H00( istart1:iend1, istart2:iend2) &
                        = H00( istart1:iend1, istart2:iend2) &
                        + HmnR( istart1- (i1-1)*Num_wann:iend1- (i1-1)*Num_wann, &
                        istart2- (i2-1)*Num_wann:iend2- (i2-1)*Num_wann, iR)*ratio/ndegen(iR)* fac

                     istart1= (i1-1)*Num_wann+ orbital_start(ia1)
                     iend1= (i1-1)*Num_wann+ orbital_start(ia1+1)- 1
                     istart2= (i2-1)*Num_wann+ orbital_start(ia2) + Nwann
                     iend2= (i2-1)*Num_wann+ orbital_start(ia2+1)- 1 + Nwann

                     H00( istart1:iend1, istart2:iend2) &
                        = H00( istart1:iend1, istart2:iend2) &
                        + HmnR( istart1- (i1-1)*Num_wann:iend1- (i1-1)*Num_wann, &
                        istart2- (i2-1)*Num_wann:iend2- (i2-1)*Num_wann, iR)*ratio/ndegen(iR)* fac

                     istart1= (i1-1)*Num_wann+ orbital_start(ia1) + Nwann
                     iend1= (i1-1)*Num_wann+ orbital_start(ia1+1)- 1 + Nwann
                     istart2= (i2-1)*Num_wann+ orbital_start(ia2) + Nwann
                     iend2= (i2-1)*Num_wann+ orbital_start(ia2+1)- 1 + Nwann

                     H00( istart1:iend1, istart2:iend2) &
                        = H00( istart1:iend1, istart2:iend2) &
                        + HmnR( istart1- (i1-1)*Num_wann:iend1- (i1-1)*Num_wann, &
                        istart2- (i2-1)*Num_wann:iend2- (i2-1)*Num_wann, iR)*ratio/ndegen(iR)* fac
                  endif ! soc
               enddo ! ia2
            enddo ! ia1
         enddo ! iR
      enddo ! i2
   enddo ! i1



   !>> calculate inter-hopping
   H01=0.0d0
   ! i1 column index
   do i1=1, Nq
      ! i2 row index
      do i2=1, Nq
         !> sum over R points to get H(k1, k2)
         do iR=1, Nrpts
            ia=dble(irvec(1,iR))
            ib=dble(irvec(2,iR))
            ic=dble(irvec(3,iR))

            !> new lattice
            call latticetransform(ia, ib, ic, new_ia, new_ib, new_ic)
            !                write(*,*) 'i',ia,ib,ic
            !                write(*,*) 'n',new_ia, new_ib, new_ic
            inew_ic= int(new_ic)
            if (inew_ic /= (i2+ Nq -i1)) cycle

            !> exp(i k.R)
            kdotr= k(1)*new_ia+ k(2)*new_ib
            ratio= cos(2d0*pi*kdotr)+ zi*sin(2d0*pi*kdotr)

            R1= (i1-1)*Ruc_new
            R2= new_ia*Rua_new+ new_ib*Rub_new+ (i2-1+ Nq)*Ruc_new
            Rp1= R1
            Rp2= R2
            !call rotate(R1, Rp1)
            !call rotate(R2, Rp2)

            do ia1=1, Origin_cell%Num_atoms
               do ia2=1, Origin_cell%Num_atoms
                  R1= Origin_cell%Atom_position_cart(:, ia1)
                  R2= Origin_cell%Atom_position_cart(:, ia2)
                  call rotate(R1, tau1)
                  call rotate(R2, tau2)

                  Ri= Rp1+ tau1
                  Rj= Rp2+ tau2
                  phase=phase2(Ri,Rj)

                  fac= cos(phase)+ zi*sin(phase)

                  istart1= (i1-1)*Num_wann+ orbital_start(ia1)
                  iend1= (i1-1)*Num_wann+ orbital_start(ia1+1)- 1
                  istart2= (i2-1)*Num_wann+ orbital_start(ia2)
                  iend2= (i2-1)*Num_wann+ orbital_start(ia2+1)- 1

                  H01( istart1:iend1, istart2:iend2) &
                     = H01( istart1:iend1, istart2:iend2) &
                     + HmnR( istart1- (i1-1)*Num_wann:iend1- (i1-1)*Num_wann, &
                     istart2- (i2-1)*Num_wann:iend2- (i2-1)*Num_wann, iR)*ratio/ndegen(iR)* fac

                  !> there is soc term in the hr file
                  if (soc>0) then
                     istart1= (i1-1)*Num_wann+ orbital_start(ia1) + Nwann
                     iend1= (i1-1)*Num_wann+ orbital_start(ia1+1)- 1 + Nwann
                     istart2= (i2-1)*Num_wann+ orbital_start(ia2)
                     iend2= (i2-1)*Num_wann+ orbital_start(ia2+1)- 1

                     H01( istart1:iend1, istart2:iend2) &
                        = H01( istart1:iend1, istart2:iend2) &
                        + HmnR( istart1- (i1-1)*Num_wann:iend1- (i1-1)*Num_wann, &
                        istart2- (i2-1)*Num_wann:iend2- (i2-1)*Num_wann, iR)*ratio/ndegen(iR)* fac

                     istart1= (i1-1)*Num_wann+ orbital_start(ia1)
                     iend1= (i1-1)*Num_wann+ orbital_start(ia1+1)- 1
                     istart2= (i2-1)*Num_wann+ orbital_start(ia2) + Nwann
                     iend2= (i2-1)*Num_wann+ orbital_start(ia2+1)- 1 + Nwann

                     H01( istart1:iend1, istart2:iend2) &
                        = H01( istart1:iend1, istart2:iend2) &
                        + HmnR( istart1- (i1-1)*Num_wann:iend1- (i1-1)*Num_wann, &
                        istart2- (i2-1)*Num_wann:iend2- (i2-1)*Num_wann, iR)*ratio/ndegen(iR)* fac

                     istart1= (i1-1)*Num_wann+ orbital_start(ia1) + Nwann
                     iend1= (i1-1)*Num_wann+ orbital_start(ia1+1)- 1 + Nwann
                     istart2= (i2-1)*Num_wann+ orbital_start(ia2) + Nwann
                     iend2= (i2-1)*Num_wann+ orbital_start(ia2+1)- 1 + Nwann

                     H01( istart1:iend1, istart2:iend2) &
                        = H01( istart1:iend1, istart2:iend2) &
                        + HmnR( istart1- (i1-1)*Num_wann:iend1- (i1-1)*Num_wann, &
                        istart2- (i2-1)*Num_wann:iend2- (i2-1)*Num_wann, iR)*ratio/ndegen(iR)* fac
                  endif ! soc
               enddo ! ia2
            enddo ! ia1
         enddo ! iR
      enddo ! i2
   enddo ! i1

   !> periodic boundary
   ham_landau= H00
   if(.not. landau_chern_calc) ham_landau=ham_landau+ H01* exp(zi*k(3)*2d0*pi) + conjg(transpose(H01))* exp(-zi*k(3)*2d0*pi)
   hct=0d0
   ndq=Ndimq
   !> check hermitcity
   do i1=1,Nq*Num_wann
      do i2=1,Nq*Num_wann

         if(abs(ham_landau(i1,i2)) > 1e-6) hct=hct+1d0
         if(abs(H00(i1,i2)-conjg(H00(i2,i1))).ge.1e-6)then
            write(*,*)'there is something wrong with ham_landau'
            write(*,*) i1,i2
            stop
         endif
!            if(abs(H00(i1,i2))>1e-6)
      enddo
   enddo
!    write(*,*) 'sparse%=',hct/(ndq*ndq),hct/ndq

   deallocate(H00, H01)

   return
end subroutine ham_3Dlandau





subroutine LTrans(Ri,Ru1,Ru2,Ru3,Ro)
   use para
   implicit none

   real(dp) :: Ri(3)
   real(dp) :: Ro(3)
   real(dp) :: Ru1(3),Ru2(3),Ru3(3)
   real(dp), allocatable :: Uinv(:, :)

   allocate(Uinv(3, 3))
   Uinv= 0d0
   Uinv(:,1)=Ru1(:)
   Uinv(:,2)=Ru2(:)
   Uinv(:,3)=Ru3(:)

   call inv_r(3, Uinv)

   Ro=matmul(Uinv,Ri)

   return
end subroutine


!> Calculate the magnetic phase between two positions
!> Author : YiFei Guan
!> Here we use the new coordinates
!> set R1 to the new x direction ex'
!> set R1\cross R2 to the new z direction ez'
!> set ey'= ez'\cross ex'
!> in this subroutine, only the magnetic field along x direction which is along 
!> the first vector in SURFACE card is well tested.
!> References :
!> 1. Hofstadter butterflies of bilayer graphene.
!>    Norbert Nemec and Gianaurelio Cuniberti, Phys. Rev. B 75, 201404(R) (2007)
!> 2. Periodic Landau gauge and quantum Hall effect in twisted bilayer graphene
!>    Yasumasa Hasegawa and Mahito Kohmoto Phys. Rev. B 88, 125426 (2013)
function phase2(Ra,Rb) result(phase)!-+--+-++
   use para
   implicit none
   integer :: i,k
   integer :: i1,i2
   real(dp) :: phase
   real(dp) :: Ri(3),Rj(3)
   real(dp) :: Ra(3),Rb(3)
   real(dp) :: frazx,frazy
   !array for counting the pieces of jump
   real(dp) :: kk1
   real(dp) :: z1,z2
   real(dp),external :: frac
   real(dp), allocatable :: piece(:)
   real(dp),parameter :: edge=1e-4,frc=1e-4

   allocate(piece(999))


   call ltrans(Ra,Rua_new,Rub_new,Ruc_new,Ri)
   call ltrans(Rb,Rua_new,Rub_new,Ruc_new,Rj)
   phase=0
   i1=floor(Ri(1)+edge)
   i2=floor(Rj(1)+edge)
   if(i1==i2) then
      frazx=frac(Ri(1))+frac(Rj(1))
      phase=phase- By*(frazx)*(Rj(3)-Ri(3))/2d0
   elseif(i1<i2) then
      piece(1)=Ri(1)
      k=2
      do i=i1+1,i2
         piece(k)=dble(i)
         k=k+1
      end do
      piece(k)=Rj(1)
      kk1=(Rj(3)-Ri(3))/(Rj(1)-Ri(1))
      do i=1,k-1
         frazx=frac(piece(i)+frc)+frac(piece(i+1)-frc)
         z1=(Ri(3)+(piece(i)-Ri(1))*kk1)
         z2=(Ri(3)+(piece(i+1)-Ri(1))*kk1)
         phase=phase-By*frazx*(z2-z1)/2d0!y \delta dz
         if(i>1) phase=phase+z1*By
      end do
   elseif(i1>i2) then
      piece(1)=Ri(1)
      k=2
      do i=i1,i2+1,-1
         piece(k)=dble(i)
         k=k+1
      end do
      piece(k)=Rj(1)
      kk1=(Rj(3)-Ri(3))/(Rj(1)-Ri(1))
      do i=1,k-1
         frazx=frac(piece(i)-frc)+frac(piece(i+1)+frc)
         z1=(Ri(3)+(piece(i)-Ri(1))*kk1)
         z2=(Ri(3)+(piece(i+1)-Ri(1))*kk1)
         phase=phase-By*frazx*(z2-z1)/2d0!y \delta dz
         if(i>1) phase=phase-z1*By
      end do
   end if

!    phase=0
   k=0
   i1=floor(Ri(2)+edge)
   i2=floor(Rj(2)+edge)
   if(i1==i2) then
      frazy=frac(Ri(2)+frc)+frac(Rj(2)+frc)
      phase=phase+Bx*(frazy)*(Rj(3)-Ri(3))/2d0
      k=0

   elseif(i1<i2) then
      piece(1)=Ri(2)
      k=2
      do i=i1+1,i2
         piece(k)=dble(i)
         k=k+1
      end do
      piece(k)=Rj(2)
      kk1=(Rj(3)-Ri(3))/(Rj(2)-Ri(2))
!        write(233,*) 'what',Ri(2),Rj(2),Ri(3),Rj(3),Bx
      do i=1,k-1
         frazy=frac(piece(i)+frc)+frac(piece(i+1)-frc)
         z1=(Ri(3)+(piece(i)-Ri(2))*kk1)
         z2=(Ri(3)+(piece(i+1)-Ri(2))*kk1)
         phase=phase+Bx*frazy*(z2-z1)/2d0!y \delta dz
         if(i>1) phase=phase-z1*Bx
      end do

   elseif(i1>i2) then
      piece(1)=Ri(2)
      k=2
      do i=i1,i2+1,-1
         piece(k)=dble(i)
         k=k+1
      end do
      piece(k)=Rj(2)
      !        write(*,*) k
      kk1=(Rj(3)-Ri(3))/(Rj(2)-Ri(2))
!        write(233,*) 'what',Ri(2),Rj(2),Ri(3),Rj(3),Bx
      do i=1,k-1
         frazy=frac(piece(i)-frc)+frac(piece(i+1)+frc)
         z1=(Ri(3)+(piece(i)-Ri(2))*kk1)
         z2=(Ri(3)+(piece(i+1)-Ri(2))*kk1)
         phase=phase+Bx*frazy*(z2-z1)/2d0!y \delta dz
         if(i>1) phase=phase+z1*Bx

      end do
   end if
   phase=-phase
!    phase=phase+pi
   return
end function phase2

function frac(x) result(fr)
   use para
   implicit none
   real(dp) :: x,fr
   fr=x-dble(floor(x))
   return
end function frac


!Mistaken and abandoned
function phase2d(Ra,Rb) result(phase)
   use para
   implicit none
   real(dp) :: phase
   real(dp) :: Ri(3),Rj(3)
   real(dp) :: Ra(3),Rb(3)
   real(dp) :: frazx,frazy
   real(dp) :: t1,t2,kk1,kk2
   real(dp),external :: frac
   integer :: i,i1,i2,itp

   call ltrans(Ra,Rua_new,Rub_new,Ruc_new,Ri)
   call ltrans(Rb,Rua_new,Rub_new,Ruc_new,Rj)
   frazx=frac(Ri(1))+frac(Rj(1))!Ri(1)-floor(Ri(1)+1e-4)+Rj(1)-floor(Rj(1)+1e-4)
   frazy=Ri(2)-floor(Ri(2)+1e-4)+Rj(2)-floor(Rj(2)+1e-4)
   phase=0
   phase= Bx*(frazy)*(Rj(3)-Ri(3))/2d0- By*(frazx)*(Rj(3)-Ri(3))/2d0
   t1=0

   i1=floor(Ri(1)+0.05)
   i2=floor(Rj(1)+0.05)
   if(i1<i2) then
      !if(abs(Rj(1)-Ri(1))<1e-8) write(*,*) '1',i,j
      kk1=(Rj(3)-Ri(3))/(Rj(1)-Ri(1))
      itp=0
      do i=i1+1,i2
         t1=t1-(Ri(3)+(dble(i)-Ri(1))*kk1)
      end do
   elseif(i1>i2)then
      !if(abs(Rj(1)-Ri(1))<1e-8) write(*,*) '1',i,j
      kk1=(Rj(3)-Ri(3))/(Rj(1)-Ri(1))
      itp=0
      do i=i2+1,i1
         t1=t1+(Ri(3)+(dble(i)-Ri(1))*kk1)
      end do
   elseif(i1==i2) then
      t1=0
   end if


   i1=floor(Ri(2)+0.05)
   i2=floor(Rj(2)+0.05)
   t2=0
   if(i1<i2) then
      !if(abs(Rj(2)-Ri(2))<1e-6) write(*,*) '2',i,j
      kk2=(Rj(3)-Ri(3))/(Rj(2)-Ri(2))
      itp=0
      do i=i1+1,i2
         t2=t2+(Ri(3)+(dble(i)-Ri(2))*kk2)!-0.2)

      end do
   elseif(i1>i2)then
      !if(abs(Rj(2)-Ri(2))<1e-6) write(*,*) '2',i,j
      kk2=(Rj(3)-Ri(3))/(Rj(2)-Ri(2))
      itp=0
      do i=i2+1,i1
         t2=t2-(Ri(3)+(dble(i)-Ri(2))*kk2)!-0.2)
      end do
   elseif(i1==i2) then
      t2=0
   end if

   !    phase= By*(Rj(3)+Ri(3))*(Rj(1)-Ri(1))  &
   !                    - Bx*(Rj(3)+Ri(3))*(Rj(2)-Ri(2))
   phase=phase-Bx*t2+By*t1
   !    itp=int(phase/(2d0*pi))
   !    phase=phase-itp*2d0*pi
   !                phase= By*(Rj(3)+Ri(3))*(Rj(1)-Ri(1))- Bx*(Rj(3)+Ri(3))*(Rj(2)-Ri(2))
   !    phase=phase+0.2*(Rj(3)-Ri(3))*Bx
!    write(233,*) Ri(2),',',Ri(3),',',Rj(2),',',Rj(3),',',phase/Bx,',',phase,',',Bx

   !    write(*,*) 'i',Ri
   !    write(*,*) 'a',Ra

   return
end function

!> calculate hamiltonian for landau levels
!> don't consider the internal atom's position
subroutine ham_3Dlandau2(Ndimq, Nq, k, ham_landau)
   use para
   implicit none

   integer, intent(in) :: Ndimq
   integer, intent(in) :: Nq
   real(dp), intent(in) :: k(3)
   complex(dp), intent(out) :: ham_landau(Ndimq, Ndimq)

   !> inta-hopping for the supercell
   complex(dp), allocatable :: H00(:, :)

   !> inter-hopping for the supercell
   complex(dp), allocatable :: H01(:, :)

   ! loop index
   integer :: i1, i2, iR, ia1

   ! index used to sign irvec
   real(dp) :: ia,ib,ic

   integer :: istart1, istart2
   integer :: iend1, iend2

   integer :: inew_ic

   !> nwann= Num_wann/2
   integer :: nwann

   integer, allocatable :: orbital_start(:)

   ! new index used to sign irvec
   real(dp) :: new_ia,new_ib,new_ic

   ! wave vector k times lattice vector R
   real(Dp) :: kdotr
   real(dp) :: phase
   complex(dp) :: ratio
   complex(dp) :: fac

   real(dp) :: Ri(3)
   real(dp) :: Rj(3)

   allocate( H00( Ndimq, Ndimq))
   allocate( H01( Ndimq, Ndimq))

   nwann= Num_wann/2
   allocate( orbital_start(Origin_cell%Num_atoms+ 1))
   orbital_start= 0
   orbital_start(1)= 1
   do ia1=1, Origin_cell%Num_atoms
      orbital_start(ia1+1)= orbital_start(ia1)+ Origin_cell%nprojs(ia1)
   enddo

   !> calculate intra-hopping
   H00=0.0d0
   ! i1 column index
   do i1=1, Nq
      ! i2 row index
      do i2=1, Nq
         if (abs(i2-i1)> ijmax) cycle
         !> sum over R points to get H(k1, k2)
         do iR=1, Nrpts
            ia=dble(irvec(1,iR))
            ib=dble(irvec(2,iR))
            ic=dble(irvec(3,iR))

            !> new lattice
            call latticetransform(ia, ib, ic, new_ia, new_ib, new_ic)

            inew_ic= int(new_ic)
            if (inew_ic /= (i2-i1)) cycle

            !> exp(i k.R)
            kdotr= k(1)*new_ia+ k(2)*new_ib
            ratio= cos(2d0*pi*kdotr)+ zi*sin(2d0*pi*kdotr)

            Ri= (i1-1)*Ruc_new
            Rj= new_ia*Rua_new+ new_ib*Rub_new+ (i2-1)*Ruc_new

            phase= By*(Rj(3)+Ri(3))*(Rj(1)-Ri(1))- Bx*(Rj(3)+Ri(3))*(Rj(2)-Ri(2))
            fac= cos(phase)+ zi*sin(phase)

            istart1= (i2-1)*Num_wann+ 1
            iend1= (i2-1)*Num_wann+ Num_wann
            istart2= (i1-1)*Num_wann+ 1
            iend2= (i1-1)*Num_wann+ Num_wann

            H00(istart1:iend1, istart2:iend2)= &
               H00(istart1:iend1, istart2:iend2)+ &
               HmnR(:, :, iR)/dble(ndegen(iR))* ratio*fac

         enddo ! iR
         !pause
      enddo ! i2
   enddo ! i1



   !> check hermitcity
   do i1=1,Nq*Num_wann
      do i2=1,Nq*Num_wann
         if(abs(H00(i1,i2)-conjg(H00(i2,i1))).ge.1e-6)then
            write(stdout,*)'there is something wrong with H00'
            write(*,*)'there is something wrong with H00'
            print *, i1, i2, H00(i1, i2)
            stop
         endif
      enddo
   enddo

   !> calculate inter-hopping
   H01=0.0d0
   ! i1 column index
   do i1=1, Nq
      ! i2 row index
      do i2=1, Nq
         if (abs(i2+Nq-i1)> ijmax) cycle
         !> sum over R points to get H(k1, k2)
         do iR=1, Nrpts
            ia=dble(irvec(1,iR))
            ib=dble(irvec(2,iR))
            ic=dble(irvec(3,iR))

            !> new lattice
            call latticetransform(ia, ib, ic, new_ia, new_ib, new_ic)

            inew_ic= int(new_ic)
            if (inew_ic /= (i2+Nq-i1)) cycle

            !> exp(i k.R)
            kdotr= k(1)*new_ia+ k(2)*new_ib
            ratio= cos(2d0*pi*kdotr)+ zi*sin(2d0*pi*kdotr)

            Ri= (i1-1)*Ruc_new
            Rj= new_ia*Rua_new+ new_ib*Rub_new+ (i2+Nq-1)*Ruc_new

            phase= By*(Rj(3)+Ri(3))*(Rj(1)-Ri(1))  &
               - Bx*(Rj(3)+Ri(3))*(Rj(2)-Ri(2))

            fac= cos(phase)+ zi*sin(phase)

            istart1= (i2-1)*Num_wann+ 1
            iend1= (i2-1)*Num_wann+ Num_wann
            istart2= (i1-1)*Num_wann+ 1
            iend2= (i1-1)*Num_wann+ Num_wann

            H01(istart1:iend1, istart2:iend2)= &
               H01(istart1:iend1, istart2:iend2)+ &
               HmnR(:, :, iR)/dble(ndegen(ir))* ratio*fac

         enddo ! iR
      enddo ! i2
   enddo ! i1


   !> periodic boundary
   ham_landau= H00 + H01* exp(zi*k(3)*2d0*pi) + conjg(transpose(H01))* exp(-zi*k(3)*2d0*pi)

   !> check hermitcity
   do i1=1,Nq*Num_wann
      do i2=1,Nq*Num_wann
         if(abs(H00(i1,i2)-conjg(H00(i2,i1))).ge.1e-6)then
            write(stdout,*)'there are something wrong with ham_landau'
            stop
         endif
      enddo
   enddo

   deallocate(H00, H01)
   return
end subroutine ham_3Dlandau2
