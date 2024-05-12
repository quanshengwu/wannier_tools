  subroutine fermisurface3D
     ! This subroutine calculates 3D fermi surface in the 1st BZ using wannier TB method
     ! 
     ! Copyright (c) 2010 QuanSheng Wu. All rights reserved.

     use wmpi
     use para

     implicit none

     integer :: ik, i, knv3, ikx, iky, ikz, ierr

     integer :: nband_min, nband_max, nband_store

     character(40) :: fsfile

     real(dp) :: k(3)
     real(dp) :: time_start, time_end
     
     ! Hamiltonian of bulk system
     complex(Dp), allocatable :: Hamk_bulk(:, :)

     real(dp) :: kxmin, kxmax, kymin, kymax, kzmin, kzmax

     real(dp), allocatable :: W(:)
     real(dp), allocatable :: eigval(:,:)
     real(dp), allocatable :: eigval_mpi(:,:)

     ! only for output the FS3D.bxsf, we don't have to output all the bands,
     ! only consider the bands close to the Fermi level 
     if (SOC == 0) then
        nband_min= Numoccupied- 5
        nband_max= Numoccupied+ 6
     else
        nband_min= Numoccupied- 7
        nband_max= Numoccupied+ 8
     endif

     if (nband_min< 1) then 
        nband_min= 1
     endif
     if (nband_max> Num_wann) then
        nband_max= Num_wann
     endif
     if (nband_min>nband_max) then
        nband_min= max(1, nband_max-4)
        if (SOC>0) nband_min= max(1, nband_max-8)
     endif

     nband_store= nband_max- nband_min+ 1
     if (cpuid.eq.0) then
        write(stdout, *)">> In fermisurface3D"
        write(stdout, '(a, i8, a, i8)')">> We are going to write out the bands from band ", nband_min, ' to', nband_max
     endif

     kxmin= 0.00d0/1d0
     kxmax= 1.00d0/1d0
     kymin= 0.00d0/1d0
     kymax= 1.00d0/1d0
     kzmin= 0.00d0/1d0
     kzmax= 1.00d0/1d0
     ik =0

     knv3= nk1*nk2*nk3
     allocate(W(Num_wann))
     allocate(Hamk_bulk(Num_wann, Num_wann))
     allocate(eigval(nband_store, knv3), stat=ierr)
     if (ierr>0) stop 'not enough memory'
     allocate(eigval_mpi(nband_store, knv3), stat=ierr)
     if (ierr>0) stop 'not enough memory'

     eigval_mpi= 0d0
     eigval= 0d0
     time_start= 0d0
     time_end= 0d0
     do ik= 1+cpuid, knv3, num_cpu
        if (cpuid==0.and. mod(ik/num_cpu, 500)==0) &
           write(stdout, *) '3DFS, ik ', ik, 'knv3',knv3, 'time left', &
           (knv3-ik)*(time_end- time_start)/num_cpu, ' s'
        call now(time_start)

        ikx= (ik-1)/(nk2*nk3)+1
        iky= ((ik-1-(ikx-1)*Nk2*Nk3)/nk3)+1
        ikz= (ik-(iky-1)*Nk3- (ikx-1)*Nk2*Nk3)
        k= K3D_start_cube+ K3D_vec1_cube*(ikx-1)/dble(nk1-1)  &
         + K3D_vec2_cube*(iky-1)/dble(nk2-1)  &
         + K3D_vec3_cube*(ikz-1)/dble(nk3-1)


        ! calculation bulk hamiltonian
        Hamk_bulk= 0d0
       !call ham_bulk_atomicgauge    (k, Hamk_bulk)
        call ham_bulk_latticegauge    (k, Hamk_bulk)
        call eigensystem_c( 'N', 'U', Num_wann, Hamk_bulk, W)
        eigval_mpi(:, ik)= W(nband_min:nband_max)
        call now(time_end)
     enddo

#if defined (MPI)
     call mpi_allreduce(eigval_mpi, eigval,size(eigval),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
     if (ierr>0) then
        print *, 'Something wrong in mpi_allreduce in fermisurface3D', ierr
        stop
     endif
#else
     eigval= eigval_mpi 
#endif
     if (cpuid==0) then
        write(stdout, *)'>> All processors finished their job for fermisurface3D'
     endif

     !> writeout eigenvalues data for each band into separate files for matlab isosurface function use.
     do i=1, nband_store
        outfileindex= outfileindex+ 1
        if (cpuid==0) then
           if (i>0 .and. i<10) then
              write(fsfile, '(a, i1, a)')'FS3D_matlab_band_',i, '.txt'
           else if (i>=10 .and. i<100) then
              write(fsfile, '(a, i2, a)')'FS3D_matlab_band_',i, '.txt'
           else if (i>=100 .and. i<1000) then
              write(fsfile, '(a, i3, a)')'FS3D_matlab_band_',i, '.txt'
           else if (i>=1000 .and. i<10000) then
              write(fsfile, '(a, i4, a)')'FS3D_matlab_band_',i, '.txt'
           endif
           open(unit=outfileindex,FILE=fsfile,STATUS='UNKNOWN',FORM='FORMATTED')
           write(outfileindex,'(a, i10)') '%BAND: ',i
           write(outfileindex,'(a, 4i16)') '%Nk1, Nk2, Nk3, total', Nk1, Nk2, Nk3, knv3
           do ik=1, knv3
              write(outfileindex,'(f16.8)') eigval(i, ik)/eV2Hartree
           enddo
           close(outfileindex)
        endif
     enddo

     outfileindex= outfileindex+ 1
     if (cpuid==0)then
        open(unit=outfileindex,FILE='FS3D.bxsf',STATUS='UNKNOWN',FORM='FORMATTED')
        write(outfileindex,'(a)') ' BEGIN_INFO'
        write(outfileindex,'(a)') '      #'
        write(outfileindex,'(a)') '      # this is a Band-XCRYSDEN-Structure-File'
        write(outfileindex,'(a)') '      # for Fermi Surface Visualisation'
        write(outfileindex,'(a)') '      #'
        write(outfileindex,'(a)') '      #  Launch as: xcrysden --bxsf FS3D.bxsf'
        write(outfileindex,'(a)') '       Fermi Energy: 0'
        write(outfileindex,'(a)') ' END_INFO'
        write(outfileindex,'(a)') 
        write(outfileindex,'(a)') ' BEGIN_BLOCK_BANDGRID_3D'
        write(outfileindex,'(a)') 'from_wannier_code'
        write(outfileindex,'(a)') ' BEGIN_BANDGRID_3D_fermi'
        write(outfileindex,'(i10)') nband_store
        write(outfileindex,'(3i10)') nk1, nk2, nk3
        write(outfileindex,'(a)') '0.0 0.0 0.0'
        write(outfileindex,'(3f16.8)') (Origin_cell%Kua(i)*Angstrom2atomic, i=1,3)
        write(outfileindex,'(3f16.8)') (Origin_cell%Kub(i)*Angstrom2atomic, i=1,3)
        write(outfileindex,'(3f16.8)') (Origin_cell%Kuc(i)*Angstrom2atomic, i=1,3)
        do i=1,nband_store
           write(outfileindex,'(a, i10)') 'BAND: ',i
           do ik=1, knv3
              write(outfileindex,'(E16.8)') eigval(i, ik)/eV2Hartree
           enddo
        enddo
        write(outfileindex,'(a)') 'END_BANDGRID_3D'
        write(outfileindex,'(a)') ' END_BLOCK_BANDGRID_3D'
        close(outfileindex)
     endif

#if defined (MPI)
     call mpi_barrier(mpi_cmw, ierr)
#endif

     deallocate(W)
     deallocate(Hamk_bulk)
     deallocate(eigval)
     deallocate(eigval_mpi)
     return
   end subroutine fermisurface3D


  subroutine orbitaltexture
     ! This subroutine calculates orbital texture in a k plane
     ! using wannier TB method

     use wmpi
     use para

     implicit none

     integer :: ia, ig
     integer :: ik, i, j, i1, i2
     integer :: knv3, nkx, nky

     integer :: ierr
     real(dp) :: kz
     real(Dp) :: k(3)
     
     real(Dp) :: dos_max
     ! Hamiltonian of bulk system
     complex(Dp), allocatable :: Hamk_bulk(:, :)

     real(dp) :: time_start, time_end
     real(dp) :: kxmin, kxmax, kymin, kymax
     real(dp) :: zmin, zmax
     real(dp) :: kxmin_shape, kxmax_shape, kymin_shape, kymax_shape

     real(dp), allocatable :: kxy(:,:)
     real(dp), allocatable :: kxy_shape(:,:)
     real(dp), allocatable :: kxy_plane(:,:)
    
     !> atomic and k resolved spectral function
     real(dp), allocatable :: dos(:, :)
     real(dp), allocatable :: dos_mpi(:, :)
     real(dp), allocatable :: dos_total(:)
     real(dp), allocatable :: dos_selected(:,:)
     real(dp), allocatable :: dos_unselected(:,:)

     complex(dp), allocatable :: ones(:,:)


     allocate(Hamk_bulk(Num_wann, Num_wann))

     nkx= Nk
     nky= Nk
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
           kxy(:, ik)= K3D_start+ K3D_vec1*(i-1)/dble(nkx-1)+ K3D_vec2*(j-1)/dble(nky-1)
           kxy_shape(:, ik)= kxy(1, ik)* Origin_cell%Kua+ kxy(2, ik)* Origin_cell%Kub+ kxy(3, ik)* Origin_cell%Kuc 
           call rotate_k3_to_kplane(kxy_shape(:, ik), kxy_plane(:, ik))
        enddo
     enddo

     i1=1
     i2=2
     kxmin_shape=minval(kxy_shape(i1,:))
     kxmax_shape=maxval(kxy_shape(i1,:))
     kymin_shape=minval(kxy_shape(i2,:))
     kymax_shape=maxval(kxy_shape(i2,:))
      
     knv3= nkx*nky
     allocate( dos    (knv3, Num_wann))
     allocate( dos_mpi(knv3, Num_wann))
     allocate( dos_total(knv3))
     allocate( dos_selected(knv3, NumberofSelectedOrbitals_groups))
     allocate( dos_unselected(knv3, NumberofSelectedOrbitals_groups))
     dos    = 0d0
     dos_mpi= 0d0

     allocate(ones(Num_wann, Num_wann))
     ones= 0d0
     do i=1, Num_wann
        ones(i, i)= 1d0
     enddo

     time_start= 0d0
     time_end= 0d0
     do ik= 1+cpuid, knv3, num_cpu
        if (cpuid==0.and. mod(ik/num_cpu, 500)==0) &
           write(stdout, *) 'FS_Plane, ik ', ik, 'knv3',knv3, 'time left', &
           (knv3-ik)*(time_end- time_start)/num_cpu, ' s'
        call now(time_start)

        k = kxy(:, ik)

        ! calculation bulk hamiltonian
        Hamk_bulk= 0d0
        if (index(KPorTB, 'KP')/=0)then
           call ham_bulk_kp(k, Hamk_bulk)
        else
          !> deal with phonon system
          if (index(Particle,'phonon')/=0.and.LOTO_correction) then
             call ham_bulk_LOTO(k, Hamk_bulk)
          else
             call ham_bulk_latticegauge    (k, Hamk_bulk)
          endif
        endif


        Hamk_bulk= (iso_energy -zi* Fermi_broadening)* ones - Hamk_bulk
        call inv(Num_wann, Hamk_bulk)
        do i=1, Num_wann
           dos(ik, i)= aimag(Hamk_bulk(i, i))/pi
        enddo
        call now(time_end)

     enddo

#if defined (MPI)
     call mpi_allreduce(dos,dos_mpi,size(dos),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
     dos_mpi= dos
#endif

     dos_total= eps9
     do ik=1, knv3
        do i=1, Num_wann
           dos_total(ik)= dos_total(ik)+ dos_mpi(ik, i)
        enddo
     enddo

     dos_selected= 0d0
     do ik=1, knv3
        do ig=1, NumberofSelectedOrbitals_groups
           do i=1, NumberofSelectedOrbitals(ig)
              dos_selected(ik, ig)= dos_selected(ik, ig)+ dos_mpi(ik, Selected_WannierOrbitals(ig)%iarray(i))
           enddo
        enddo
        dos_unselected(ik, ig)= dos_total(ik)- dos_selected(ik, ig)
     enddo


     dos_max= maxval(dos_total)/4d0

     outfileindex= outfileindex+ 1
     if (cpuid==0)then
        open(unit=outfileindex, file='fs.dat')
        write(outfileindex, '(80a16, a)')'# kx', 'ky', 'kz', 'kp1', 'kp2', 'kp3', 'total A', 'A1(k,E)', 'A2(k,E)'
        do ik=1, knv3
           if (dos_total(ik)>dos_max) then
              write(outfileindex, '(3000f16.8)')kxy_shape(:, ik)*Angstrom2atomic, kxy_plane(:, ik)*Angstrom2atomic, &
                dos_total(ik), ((dos_selected(ik, ig))/dos_total(ik)*500d0+100, &
                (dos_unselected(ik, ig))/dos_total(ik)*500d0+100, ig=1, NumberofSelectedOrbitals_groups)
           else
              write(outfileindex, '(3000f16.8)')kxy_shape(:, ik)*Angstrom2atomic, kxy_plane(:, ik)*Angstrom2atomic, &
                 dos_total(ik), ((dos_selected(ik, ig)), (dos_unselected(ik, ig)), &
                 ig=1, NumberofSelectedOrbitals_groups)
           endif
           if (mod(ik, nky)==0) write(outfileindex, *)' '
        enddo
        close(outfileindex)
     endif
#if defined (MPI)
     call mpi_barrier(mpi_cmw, ierr)
#endif


     outfileindex= outfileindex+ 1
     if (cpuid==0)then
        open(unit=outfileindex, file='fs.dat-matlab')
        write(outfileindex, '(80a16, a)')'% kx', 'ky', 'kz', 'kp1', 'kp2', 'kp3', 'total A', 'A1(k,E)', 'A2(k,E)'
        write(outfileindex, '(a16, 2i16)')'% Nk1, Nk2=', Nk1, Nk2
        do ik=1, knv3
           write(outfileindex, '(3000f16.8)')kxy_shape(:, ik)*Angstrom2atomic, kxy_plane(:, ik)*Angstrom2atomic, &
                dos_total(ik), ((dos_selected(ik, ig)), (dos_unselected(ik, ig)), &
                ig=1, NumberofSelectedOrbitals_groups)
        enddo
        close(outfileindex)
     endif
#if defined (MPI)
     call mpi_barrier(mpi_cmw, ierr)
#endif


     outfileindex= outfileindex+ 1
     if (cpuid==0)then
        open(unit=outfileindex, file='orbitaltexture.dat')
        write(outfileindex, '(8a16, a)')'# kx', 'ky', 'kz', 'kp1', 'kp2', 'kp3', '(A1(k,E))', 'A2(k,E)'
        do ik=1, knv3
           if (dos_total(ik)>dos_max) then
              write(outfileindex, '(300f16.8)')kxy_shape(:, ik)*Angstrom2atomic, kxy_plane(:, ik)*Angstrom2atomic, &
                   (dos_selected(ik, ig), dos_unselected(ik, ig), ig=1, NumberofSelectedOrbitals_groups)
           endif

           if (mod(ik, nky)==0) write(outfileindex, *)' '
        enddo
        close(outfileindex)
     endif
#if defined (MPI)
     call mpi_barrier(mpi_cmw, ierr)
#endif



     zmax= maxval(log(dos_total))
     zmin= minval(log(dos_total))

     !> minimum and maximum value of energy bands

     outfileindex= outfileindex+ 1
     !> write script for gnuplot
     if (cpuid==0) then
        open(unit=outfileindex, file='fs.gnu')
        write(outfileindex, '(a)')"set encoding iso_8859_1"
        write(outfileindex, '(a)')'#set terminal  postscript enhanced color'
        write(outfileindex, '(a)')"#set output 'fs.eps'"
        write(outfileindex, '(3a)')'set terminal  pngcairo truecolor enhanced', &
           ' size 1920, 1680 font ",36"'
        write(outfileindex, '(a)')"set output 'fs.png'"
        write(outfileindex,'(a, f10.4, 2a, f10.4, a)') &
           '#set palette defined ( ', zmin, ' "white", ', &
          '0 "black", ', zmax,'  "red" )'
        write(outfileindex,'(a)') &
           'set palette defined ( 0 "white", 100 "white", 101 "green", 500 "red" )'
        write(outfileindex, '(a)')'#set palette rgbformulae 33,13,10'
        write(outfileindex, '(a)')'unset ztics'
        write(outfileindex, '(a)')'unset key'
        write(outfileindex, '(a)')'set pm3d'
        write(outfileindex, '(a)')'#set view equal xyz'
        write(outfileindex, '(a)')'set view map'
        write(outfileindex, '(a)')'set border lw 3'
        write(outfileindex, '(a)')'#set xtics font ",24"'
        write(outfileindex, '(a)')'#set ytics font ",24"'
        write(outfileindex, '(a)')'set size ratio -1'
        write(outfileindex, '(a)')'unset xtics'
        write(outfileindex, '(a)')'unset ytics'
        write(outfileindex, '(a)')'set colorbox'
        write(outfileindex, '(a)')'set autoscale fix'
        write(outfileindex, '(a)')'set pm3d interpolate 2,2'
        write(outfileindex, '(2a)')"splot 'fs.dat' u 4:5:(($8)) w pm3d , \"
        write(outfileindex, '(a)')"     'orbitaltexture.dat' u 4:5:(0):($7/90000):($8/90000):(0) \"
        write(outfileindex, '(a)')" w vec  head lw 5 lc rgb 'orange' front"

        close(outfileindex)
     endif


   return
   end subroutine orbitaltexture


  subroutine fermisurface_kplane
     ! This subroutine calculates 3D fermi surface in the a fixed k plane
     ! using wannier TB method

     use wmpi
     use para

     implicit none

     integer :: ia, ik, i, j, i1, i2, ig, io
     integer :: knv3, nkx, nky, ierr, nwann
     real(dp) :: kz, k(3), s0(3), s0t(3), s0_len
     
     ! Hamiltonian of bulk system
     complex(Dp), allocatable :: Hamk_bulk(:, :)

     real(dp) :: time_start, time_end
     real(dp) :: kxmin, kxmax, kymin, kymax
     real(dp) :: zmin, zmax, dos_selected_max
     real(dp) :: kxmin_shape, kxmax_shape, kymin_shape, kymax_shape

     real(dp), allocatable :: kxy(:,:), kxy_shape(:,:), kxy_plane(:,:)
    
     !> Wannier orbital projected and k-resolved spectral function
     real(dp), allocatable :: dos_selected(:, :), dos_selected_mpi(:, :)
     real(dp), allocatable :: dos_total(:), dos_total_mpi(:), s1(:, :)
     real(dp), allocatable :: sx_selected(:, :), sy_selected(:, :), sz_selected(:, :)
     real(dp), allocatable :: sx_selected_mpi(:, :), sy_selected_mpi(:, :), sz_selected_mpi(:, :)

     !> spin 1/2 Pauli matrix in Wannier basis
     complex(Dp),allocatable :: spin_sigma_x(:,:), spin_sigma_y(:,:), spin_sigma_z(:,:)
     
     complex(dp), allocatable :: ones(:, :), ctemp(:, :)

     logical, external :: if_in_the_WS
     real(dp), external :: norm


     nkx= Nk
     nky= Nk
     knv3= nkx*nky
     allocate( kxy(3, knv3))
     allocate( kxy_shape(3, knv3))
     allocate( kxy_plane(3, knv3))
     kxy=0d0
     kxy_shape=0d0
     kxy_plane=0d0
     
     
     ik =0
     do i= 1, nkx
        do j= 1, nky
           ik =ik +1
           kxy(:, ik)= K3D_start+ K3D_vec1*(i-1)/dble(nkx-1)+ K3D_vec2*(j-1)/dble(nky-1) &
              -(K3D_vec1+K3D_vec2)/2d0
           kxy_shape(:, ik)= kxy(1, ik)* Origin_cell%Kua+ kxy(2, ik)* Origin_cell%Kub+ kxy(3, ik)* Origin_cell%Kuc 
           call rotate_k3_to_kplane(kxy_shape(:, ik), kxy_plane(:, ik))
        enddo
     enddo

     i1=1
     i2=2
     kxmin_shape=minval(kxy_shape(i1,:))
     kxmax_shape=maxval(kxy_shape(i1,:))
     kymin_shape=minval(kxy_shape(i2,:))
     kymax_shape=maxval(kxy_shape(i2,:))
      
     allocate(Hamk_bulk(Num_wann, Num_wann))
     allocate( dos_total(knv3), dos_total_mpi(knv3))
     allocate( dos_selected    (knv3, NumberofSelectedOrbitals_groups))
     allocate( dos_selected_mpi(knv3, NumberofSelectedOrbitals_groups))
     dos_total= eps12; dos_total_mpi= eps12
     dos_selected= eps12; dos_selected_mpi= eps12

     if (SOC>0 .and. BulkSpintexture_calc) then
        allocate( sx_selected(knv3, NumberofSelectedOrbitals_groups))
        allocate( sy_selected(knv3, NumberofSelectedOrbitals_groups))
        allocate( sz_selected(knv3, NumberofSelectedOrbitals_groups))
        allocate( sx_selected_mpi(knv3, NumberofSelectedOrbitals_groups))
        allocate( sy_selected_mpi(knv3, NumberofSelectedOrbitals_groups))
        allocate( sz_selected_mpi(knv3, NumberofSelectedOrbitals_groups))
        allocate(spin_sigma_x(Num_wann,Num_wann), spin_sigma_y(Num_wann,Num_wann), spin_sigma_z(Num_wann,Num_wann))
        sx_selected= 0d0; sx_selected_mpi= 0d0
        sy_selected= 0d0; sy_selected_mpi= 0d0
        sz_selected= 0d0; sz_selected_mpi= 0d0
        spin_sigma_x= 0d0; spin_sigma_y= 0d0; spin_sigma_z= 0d0
     endif

     if (SOC>0 .and. BulkSpintexture_calc) then
        if (index(Particle,'phonon')/=0) then
           stop "ERROR: we don't support spintexture calculation for phonon system"
        endif
        Nwann= Num_wann/2
        !> spin operator matrix
        !> this part is package dependent. 
       !if (index( Package, 'VASP')/=0.or. index( Package, 'Wien2k')/=0 &
       !   .or. index( Package, 'Abinit')/=0.or. index( Package, 'openmx')/=0) then
           do j=1, Nwann
              spin_sigma_x(j, Nwann+j)=1.0d0
              spin_sigma_x(j+Nwann, j)=1.0d0
              spin_sigma_y(j, Nwann+j)=-zi
              spin_sigma_y(j+Nwann, j)=zi
              spin_sigma_z(j, j)= 1d0
              spin_sigma_z(j+Nwann, j+Nwann)=-1d0
           enddo
       !else
       !   if (cpuid.eq.0) write(stdout, *)'Error: please report your software generating tight binding and wannier90.wout to me'
       !   if (cpuid.eq.0) write(stdout, *)'wuquansheng@gmail.com'
       !   stop 'Error: please report your software and wannier90.wout to wuquansheng@gmail.com'
       !endif
     endif

     allocate(ones(Num_wann, Num_wann), ctemp(Num_wann, Num_wann))
     ones= 0d0
     do i=1, Num_wann
        ones(i, i)= 1d0
     enddo

     time_start= 0d0
     time_end= 0d0
     do ik= 1+cpuid, knv3, num_cpu
        if (cpuid==0.and. mod(ik/num_cpu, 500)==0) &
           write(stdout, *) 'FS_Plane, ik ', ik, 'knv3',knv3, 'time left', &
           (knv3-ik)*(time_end- time_start)/num_cpu, ' s'
        call now(time_start)

        k = kxy(:, ik)

        !> if the k points is not in the first Wigner-Seitz cell, then ignore it
        if (Translate_to_WS_calc)then
           if (.not.if_in_the_WS(k)) then
              cycle
           endif
        endif

        ! calculation bulk hamiltonian
        Hamk_bulk= 0d0
        call ham_bulk_latticegauge(k, Hamk_bulk)

        Hamk_bulk= (iso_energy -zi* Fermi_broadening)* ones - Hamk_bulk

        !> Hamk_bulk is the Green's function after subroutine inv.
        call inv(Num_wann, Hamk_bulk)

        do i=1, Num_wann
           dos_total(ik)= dos_total(ik)+ aimag(Hamk_bulk(i, i))/pi
        enddo

        do ig=1, NumberofSelectedOrbitals_groups
           do i=1, NumberofSelectedOrbitals(ig)
              io=Selected_WannierOrbitals(ig)%iarray(i)
              dos_selected(ik, ig)= dos_selected(ik, ig)+ aimag(Hamk_bulk(io, io))/pi
           enddo
        enddo
        call now(time_end)

        if (SOC>0.and.BulkSpintexture_calc) then
           !>> calculate spin-resolved bulk energy spectrum
           !> spin x
           call mat_mul(Num_wann, Hamk_bulk, spin_sigma_x, ctemp)
           do ig=1, NumberofSelectedOrbitals_groups
              do i=1, NumberofSelectedOrbitals(ig)
                 io=Selected_WannierOrbitals(ig)%iarray(i)
                 sx_selected(ik, ig)= sx_selected(ik, ig)+ aimag(ctemp(io, io))/pi
              enddo
           enddo
   
           !> spin y
           call mat_mul(Num_wann, Hamk_bulk, spin_sigma_y, ctemp)
           do ig=1, NumberofSelectedOrbitals_groups
              do i=1, NumberofSelectedOrbitals(ig)
                 io=Selected_WannierOrbitals(ig)%iarray(i)
                 sy_selected(ik, ig)= sy_selected(ik, ig)+ aimag(ctemp(io, io))/pi
              enddo
           enddo
   
           !> spin z
           call mat_mul(Num_wann, Hamk_bulk, spin_sigma_z, ctemp)
           do ig=1, NumberofSelectedOrbitals_groups
              do i=1, NumberofSelectedOrbitals(ig)
                 io=Selected_WannierOrbitals(ig)%iarray(i)
                 sz_selected(ik, ig)= sz_selected(ik, ig)+ aimag(ctemp(io, io))/pi
              enddo
           enddo
        endif
     enddo

#if defined (MPI)
     call mpi_allreduce(dos_total,dos_total_mpi,size(dos_total),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
     call mpi_allreduce(dos_selected,dos_selected_mpi,size(dos_selected),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)

     if (SOC>0.and.BulkSpintexture_calc) then
        call mpi_allreduce(sx_selected,sx_selected_mpi,size(sx_selected),&
                          mpi_dp,mpi_sum,mpi_cmw,ierr)
        call mpi_allreduce(sy_selected,sy_selected_mpi,size(sy_selected),&
                          mpi_dp,mpi_sum,mpi_cmw,ierr)
        call mpi_allreduce(sz_selected,sz_selected_mpi,size(sz_selected),&
                          mpi_dp,mpi_sum,mpi_cmw,ierr)
     endif
#else
     dos_total_mpi= dos_total
     dos_selected_mpi= dos_selected
     if (SOC>0.and.BulkSpintexture_calc) then
        sx_selected_mpi= sx_selected
        sy_selected_mpi= sy_selected
        sz_selected_mpi= sz_selected
     endif
#endif

     if (SOC>0.and.BulkSpintexture_calc) then
        do ig=1, NumberofSelectedOrbitals_groups
           dos_selected_max= maxval(dos_selected_mpi(:, ig))
           do ik=1, knv3
              if (dos_selected_mpi(ik, ig)<dos_selected_max/20d0) then
                 sx_selected_mpi(ik, ig)= eps12
                 sy_selected_mpi(ik, ig)= eps12
                 sz_selected_mpi(ik, ig)= eps12
              endif
           enddo
        enddo
   
        allocate(s1(3, NumberofSelectedOrbitals_groups))
        s1= 0d0
     endif

     outfileindex= outfileindex+ 1
     if (SOC>0.and.BulkSpintexture_calc) then
        if (cpuid==0)then
           open(unit=outfileindex, file='bulkspintext.dat')
           write(outfileindex, &
              "('#', a12, 5a16, 3X, '| A(k,E)', a6, 100(4X,'|group ', i2, ': A ', 13X, 'sx', 14X, 'sy', 14X, 'sz'))")&
              'kx', 'ky', 'kz', 'kp1', 'kp2', 'kp3', 'total',&
              (i, i=1, NumberofSelectedOrbitals_groups)
           write(outfileindex, "('#column', i5, 3000i16)")(i, i=1, 7+4*NumberofSelectedOrbitals_groups)
           do ik=1, knv3
              do ig=1, NumberofSelectedOrbitals_groups
                 s0(1)= sx_selected_mpi(ik, ig)   /dos_selected_mpi(ik, ig)
                 s0(2)= sy_selected_mpi(ik, ig)   /dos_selected_mpi(ik, ig)
                 s0(3)= sz_selected_mpi(ik, ig)   /dos_selected_mpi(ik, ig)

                !s0_len= norm(s0t)
                !if (s0_len>eps6) then
                !   s0= s0t/s0_len
                !endif
                 call rotate_k3_to_kplane(s0, s1(:, ig))
              enddo
              write(outfileindex, '(3000E16.8)')kxy_shape(:, ik)*Angstrom2atomic, kxy_plane(:, ik)*Angstrom2atomic, &
                 dos_total_mpi(ik), (dos_selected_mpi(ik, ig), s1(:, ig), ig=1, NumberofSelectedOrbitals_groups)
              if (mod(ik, nky)==0) write(outfileindex, *)' '
           enddo
           close(outfileindex)
        endif
     else
        if (cpuid==0)then
           open(unit=outfileindex, file='fs_kplane.dat')
           write(outfileindex, "('#', a12, 5a16, 3X, '| A(k,E)', a6, 100(8X,'group ', i2))")&
              'kx', 'ky', 'kz', 'kp1', 'kp2', 'kp3', 'total',&
              (i, i=1, NumberofSelectedOrbitals_groups)
           write(outfileindex, "('#column', i5, 3000i16)")(i, i=1, 7+NumberofSelectedOrbitals_groups)
           do ik=1, knv3
              write(outfileindex, '(3000E16.8)')kxy_shape(:, ik)*Angstrom2atomic, kxy_plane(:, ik)*Angstrom2atomic, &
                 dos_total_mpi(ik), (dos_selected_mpi(ik, ig), ig=1, NumberofSelectedOrbitals_groups)
              if (mod(ik, nky)==0) write(outfileindex, *)' '
           enddo
           close(outfileindex)
        endif
     endif
     
#if defined (MPI)
     call mpi_barrier(mpi_cmw, ierr)
#endif

     !> minimum and maximum value of dos_total
     zmax= maxval(log(dos_total_mpi))
     zmin= minval(log(dos_total_mpi))

     outfileindex= outfileindex+ 1
     !> write script for gnuplot
     if (cpuid==0) then
        if (SOC>0.and.BulkSpintexture_calc) then
           open(unit=outfileindex, file='bulkspintext.gnu')
           write(outfileindex, '(a)')"set encoding iso_8859_1"
           write(outfileindex, '(3a)')'set terminal  pngcairo truecolor enhanced', &
              ' size 1920, 1680 font ",36"'
           write(outfileindex, '(a)')"set output 'bulkspintext.png'"
           write(outfileindex,'(a, f10.4, a, f10.4, a, f10.4, a)') &
              'set palette defined ( ', zmin, ' "white", ', &
             (zmin+zmax)/2d0,' "black", ', zmax,'  "red" )'
           write(outfileindex, '(a)')'#set palette rgbformulae 33,13,10'
           write(outfileindex, '(a)')'unset ztics'
           write(outfileindex, '(a)')'unset key'
           write(outfileindex, '(a)')'set pm3d'
           write(outfileindex, '(a)')'#set view equal xyz'
           write(outfileindex, '(a)')'set view map'
           write(outfileindex, '(a)')'set border lw 3'
           write(outfileindex, '(a)')'#set xtics font ",24"'
           write(outfileindex, '(a)')'#set ytics font ",24"'
           write(outfileindex, '(a)')'set size ratio -1'
           write(outfileindex, '(a)')'unset xtics'
           write(outfileindex, '(a)')'unset ytics'
           write(outfileindex, '(a)')'set colorbox'
           write(outfileindex, '(a)')'set autoscale fix'
           write(outfileindex, '(a)')'set pm3d interpolate 2,2'
           write(outfileindex, '(2a)')"splot 'bulkspintext.dat' u 4:5:(log($8+0.1)) w pm3d, \"
           write(outfileindex, '(a)')"    'bulkspintext.dat' u 4:5:(0):($9/5.00):($10/5.00):(0)  w vec  head lw 5 lc rgb 'orange' front"
           close(outfileindex)
       
        else
           open(unit=outfileindex, file='fs_kplane.gnu')
           write(outfileindex, '(a)')"set encoding iso_8859_1"
           write(outfileindex, '(3a)')'set terminal  pngcairo truecolor enhanced', &
              ' size 1920, 1680 font ",36"'
           write(outfileindex, '(a)')"set output 'fs_kplane.png'"
           write(outfileindex,'(a, f10.4, a, f10.4, a, f10.4, a)') &
              'set palette defined ( ', zmin, ' "white", ', &
             (zmin+zmax)/2d0,' "black", ', zmax,'  "red" )'
           write(outfileindex, '(a)')'#set palette rgbformulae 33,13,10'
           write(outfileindex, '(a)')'unset ztics'
           write(outfileindex, '(a)')'unset key'
           write(outfileindex, '(a)')'set pm3d'
           write(outfileindex, '(a)')'#set view equal xyz'
           write(outfileindex, '(a)')'set view map'
           write(outfileindex, '(a)')'set border lw 3'
           write(outfileindex, '(a)')'#set xtics font ",24"'
           write(outfileindex, '(a)')'#set ytics font ",24"'
           write(outfileindex, '(a)')'set size ratio -1'
           write(outfileindex, '(a)')'unset xtics'
           write(outfileindex, '(a)')'unset ytics'
           write(outfileindex, '(a)')'set colorbox'
           write(outfileindex, '(a)')'set autoscale fix'
           write(outfileindex, '(a)')'set pm3d interpolate 2,2'
           write(outfileindex, '(2a)')"splot 'fs_kplane.dat' u 4:5:(log($7+0.1)) w pm3d"
       
           close(outfileindex)
        endif
     endif


   return
end subroutine fermisurface_kplane

subroutine fermisurface_stack
     !> This subroutine calculates 3D fermi surfaces in the a fixed k plane
     !> using wannier TB method
     !> With the given k-plane set by KPLANE_BULK, we calculate a set of spectrum with the kplane
     !> parallel to the KPLANE_BULK

     use wmpi
     use para

     implicit none

     integer :: ia, ik, i, j, i1, i2, ig
     integer :: knv3, nkx, nky, ierr, ikz
     real(dp) :: kz(3), k(3), k1(3), k2(3)
     
     ! Hamiltonian of bulk system
     complex(Dp), allocatable :: Hamk_bulk(:, :)

     real(dp) :: time_start, time_end
     real(dp) :: kxmin, kxmax, kymin, kymax
     real(dp) :: zmin, zmax
     real(dp) :: kxmin_shape, kxmax_shape, kymin_shape, kymax_shape

     real(dp), allocatable :: kxy(:,:), kxy_shape(:,:), kxy_plane(:,:)
    
     !> atomic and k resolved spectral function
     real(dp), allocatable :: dos(:, :), dos_mpi(:, :), dos_atom(:, :)
     real(dp), allocatable :: dos_atom_mpi(:, :), dos_total(:, :)

     complex(dp), allocatable :: ones(:,:)

     integer, allocatable :: orbitals_end(:)

     logical, external :: if_in_the_WS
     real(dp), external :: norm

     allocate(Hamk_bulk(Num_wann, Num_wann))
     allocate(orbitals_end(Origin_cell%Num_atoms))

     orbitals_end(1)= Origin_cell%nprojs(1)
     do ia=2, Origin_cell%Num_atoms
        orbitals_end(ia)= orbitals_start(ia)+ Origin_cell%nprojs(ia)- 1
     enddo

     nkx= Nk
     nky= Nk
     allocate( kxy(3, nkx*nky))
     allocate( kxy_shape(3, nkx*nky))
     allocate( kxy_plane(3, nkx*nky))
     kxy=0d0
     kxy_shape=0d0
     kxy_plane=0d0
    
     !> a k-plane defined by KPLANE_BULK
     ik =0
     do i= 1, nkx
        do j= 1, nky
           ik =ik +1
           kxy(:, ik)= K3D_start+ K3D_vec1*(i-1)/dble(nkx-1)+ K3D_vec2*(j-1)/dble(nky-1) &
              -(K3D_vec1+K3D_vec2)/2d0
           kxy_shape(:, ik)= kxy(1, ik)* Origin_cell%Kua+ kxy(2, ik)* Origin_cell%Kub+ kxy(3, ik)* Origin_cell%Kuc 
           call rotate_k3_to_kplane(kxy_shape(:, ik), kxy_plane(:, ik))
        enddo
     enddo

     !> find a vector that perpendicular to the k-plane
     call direct_cart_rec(K3D_vec1, k1)
     call direct_cart_rec(K3D_vec2, k2)
    
     call cross_product(k1, k2, kz)
     kz=kz/norm(kz)*2d0
     k2=kz
     call cart_direct_rec(k2, kz)

     i1=1
     i2=2
     kxmin_shape=minval(kxy_shape(i1,:))
     kxmax_shape=maxval(kxy_shape(i1,:))
     kymin_shape=minval(kxy_shape(i2,:))
     kymax_shape=maxval(kxy_shape(i2,:))
      
     knv3= nkx*nky
     allocate( dos_atom    (knv3, Origin_cell%Num_atoms))
     allocate( dos_atom_mpi(knv3, Origin_cell%Num_atoms))
     allocate( dos    (knv3, Num_wann))
     allocate( dos_mpi(knv3, Num_wann))
     allocate( dos_total(knv3, NumberofSelectedOrbitals_groups))
     dos    = 0d0
     dos_mpi= 0d0
     dos_atom    = 0d0
     dos_atom_mpi= 0d0

     allocate(ones(Num_wann, Num_wann))
     ones= 0d0
     do i=1, Num_wann
        ones(i, i)= 1d0
     enddo

     time_start= 0d0
     time_end= 0d0
     do ik= 1+cpuid, knv3, num_cpu
        if (cpuid==0.and. mod(ik/num_cpu, 500)==0) &
           write(stdout, *) 'FS_Plane, ik ', ik, 'knv3',knv3, 'time left', &
           (knv3-ik)*(time_end- time_start)/num_cpu, ' s'
        call now(time_start)

        k = kxy(:, ik)

        if (Translate_to_WS_calc)then
           if (.not.if_in_the_WS(k)) then
              dos(ik, :)= eps6
              dos_atom(ik, :)= eps6
              cycle
           endif
        endif

        do ikz=1, Nk3
           k = kxy(:, ik)+ kz*(ikz-1)/Nk3

           ! calculation bulk hamiltonian
           Hamk_bulk= 0d0
           call ham_bulk_latticegauge(k, Hamk_bulk)
   
           Hamk_bulk= (iso_energy -zi* Fermi_broadening)* ones - Hamk_bulk
           call inv(Num_wann, Hamk_bulk)
   
           do i=1, Num_wann
              dos(ik, i)=dos(ik, i)+ aimag(Hamk_bulk(i, i))/pi
           enddo
   
           do ia=1, Origin_cell%Num_atoms
              do i=orbitals_start(ia), orbitals_end(ia)
                 dos_atom(ik, ia)= dos_atom(ik, ia)+ aimag(Hamk_bulk(i, i))/pi
              enddo
              if (SOC>0) then
                 do i=orbitals_start(ia)+num_wann/2, orbitals_end(ia)+num_wann/2
                    dos_atom(ik, ia)= dos_atom(ik, ia)+ aimag(Hamk_bulk(i, i))/pi
                 enddo
              endif
           enddo ! ia
        enddo ! ikz

        call now(time_end)
     enddo

#if defined (MPI)
     call mpi_allreduce(dos,dos_mpi,size(dos),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
     call mpi_allreduce(dos_atom,dos_atom_mpi,size(dos_atom),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
     dos_mpi= dos
     dos_atom_mpi= dos_atom
#endif

     dos_total= eps9
     do ik=1, knv3
        do ig=1, NumberofSelectedOrbitals_groups
           do i=1, NumberofSelectedOrbitals(ig)
              dos_total(ik, ig)= dos_total(ik, ig)+ dos_mpi(ik, Selected_WannierOrbitals(ig)%iarray(i))
           enddo
        enddo
     enddo


     outfileindex= outfileindex+ 1
     if (cpuid==0)then
        open(unit=outfileindex, file='fs_stack.dat')
        write(outfileindex, '(7a16, a)')'# kx', 'ky', 'kz', 'kp1', 'kp2', 'kp3', '(A(k,E))', 'A(i,k,E),i=1,numatoms'
        do ik=1, knv3
           write(outfileindex, '(3000f16.8)')kxy_shape(:, ik)*Angstrom2atomic, kxy_plane(:, ik)*Angstrom2atomic, &
             (dos_total(ik, ig), ig=1, NumberofSelectedOrbitals_groups), (dos_atom_mpi(ik, :))
           if (mod(ik, nky)==0) write(outfileindex, *)' '
        enddo
        close(outfileindex)
     endif
     zmax= maxval(log(dos_total))
     zmin= minval(log(dos_total))

     !> minimum and maximum value of energy bands

     outfileindex= outfileindex+ 1
     !> write script for gnuplot
     if (cpuid==0) then
        open(unit=outfileindex, file='fs_stack.gnu')
        write(outfileindex, '(a)')"set encoding iso_8859_1"
        write(outfileindex, '(a)')'#set terminal  postscript enhanced color'
        write(outfileindex, '(a)')"#set output 'fs_stack.eps'"
        write(outfileindex, '(3a)')'set terminal  pngcairo truecolor enhanced', &
           ' size 1920, 1680 font ",36"'
        write(outfileindex, '(a)')"set output 'fs_stack.png'"
        write(outfileindex,'(a, f10.4, 2a, f10.4, a)') &
           'set palette defined ( ', zmin, ' "white", ', &
          '0 "black", ', zmax,'  "red" )'
        write(outfileindex, '(a)')'#set palette rgbformulae 33,13,10'
        write(outfileindex, '(a)')'unset ztics'
        write(outfileindex, '(a)')'unset key'
        write(outfileindex, '(a)')'set pm3d'
        write(outfileindex, '(a)')'#set view equal xyz'
        write(outfileindex, '(a)')'set view map'
        write(outfileindex, '(a)')'set border lw 3'
        write(outfileindex, '(a)')'#set xtics font ",24"'
        write(outfileindex, '(a)')'#set ytics font ",24"'
        write(outfileindex, '(a)')'set size ratio -1'
        write(outfileindex, '(a)')'unset xtics'
        write(outfileindex, '(a)')'unset ytics'
        write(outfileindex, '(a)')'set colorbox'
       !write(outfileindex, '(a, f10.5, a, f10.5, a)')'set xrange [', kxmin, ':', kxmax, ']'
       !write(outfileindex, '(a, f10.5, a, f10.5, a)')'set yrange [', kymin, ':', kymax, ']'
        write(outfileindex, '(a)')'set autoscale fix'
       !write(outfileindex, '(a, f10.5, a, f10.5, a)')'set xrange [', kxmin_shape, ':', kxmax_shape, ']'
       !write(outfileindex, '(a, f10.5, a, f10.5, a)')'set yrange [', kymin_shape, ':', kymax_shape, ']'
        write(outfileindex, '(a)')'set pm3d interpolate 2,2'
        write(outfileindex, '(2a)')"splot 'fs_stack.dat' u 4:5:(log($7+0.1)) w pm3d"

        close(outfileindex)
     endif


   return
   end subroutine fermisurface_stack

   subroutine gapshape3D
      ! This subroutine get the k points in the BZ at which the gap is smaller then
      ! Gap_threshold

      use wmpi
      use para
      
      implicit none
      
      integer :: ik, i, j, l
      integer :: knv3
      integer :: nkx
      integer :: nky
      integer :: nkz
      
      integer :: ierr
      real(Dp) :: k(3)
      
      ! Hamiltonian of bulk system
      complex(Dp), allocatable :: Hamk_bulk(:, :) 
      
      real(dp) :: kxmin_shape, kxmax_shape, kymin_shape
      real(dp) :: kymax_shape, kzmin_shape, kzmax_shape
      
      real(dp), allocatable :: kxy(:,:)
      real(dp), allocatable :: kxy_shape(:,:)
      
      real(dp), allocatable :: gap(:, :)
      real(dp), allocatable :: gap_mpi(:, :)
      real(dp), allocatable :: W(:)
      
      complex(dp), allocatable :: ones(:,:)
      
      nkx= Nk1
      nky= Nk2
      nkz= Nk3
      allocate( kxy(3, nkx*nky*nkz))
      allocate( kxy_shape(3, nkx*nky*nkz))
      kxy=0d0
      kxy_shape=0d0
      
      ik =0

      do i= 1, nkx
         do j= 1, nky
            do l= 1, nkz
               ik= ik+ 1
               kxy(:, ik)= K3D_start_cube+ K3D_vec1_cube*(i-1)/dble(nkx-1)  &
                         + K3D_vec2_cube*(j-1)/dble(nky-1)  &
                         + K3D_vec3_cube*(l-1)/dble(nkz-1)
               kxy_shape(:, ik)= kxy(1, ik)* Origin_cell%Kua+ kxy(2, ik)* Origin_cell%Kub+ kxy(3, ik)* Origin_cell%Kuc 
            enddo
         enddo
      enddo

      kxmin_shape=minval(kxy_shape(1,:))
      kxmax_shape=maxval(kxy_shape(1,:))
      kymin_shape=minval(kxy_shape(2,:))
      kymax_shape=maxval(kxy_shape(2,:))
      kzmin_shape=minval(kxy_shape(3,:))
      kzmax_shape=maxval(kxy_shape(3,:))
      
      
      knv3= nkx*nky*nkz
      allocate( gap    (3, knv3))
      allocate( gap_mpi(3, knv3))
      gap    = 0d0
      gap_mpi= 0d0
      
      allocate(W(Num_wann))
      allocate(ones(Num_wann, Num_wann))
      allocate(Hamk_bulk(Num_wann, Num_wann))
      W= 0d0
      ones= 0d0
      do i=1, Num_wann
         ones(i, i)= 1d0
      enddo
      
      if (Numoccupied> Num_wann) then
         stop 'Numoccupied should less than Num_wann'
      endif
      
      do ik= 1+cpuid, knv3, num_cpu
         if (cpuid==0) write(stdout, *) 'Gap3D, ik, knv3', ik, knv3
      
         k(1) = kxy(1, ik)
         k(2) = kxy(2, ik)
         k(3) = kxy(3, ik)
      
         ! calculation bulk hamiltonian
         Hamk_bulk= 0d0
         call ham_bulk_latticegauge(k, Hamk_bulk)
      
         call eigensystem_c( 'N', 'U', Num_wann ,Hamk_bulk, W)
         gap(1, ik)= W(Numoccupied+1)- W(Numoccupied)
         gap(2, ik)= W(Numoccupied)
         gap(3, ik)= W(Numoccupied+1)
      
      enddo
     
      gap_mpi = 0d0
#if defined (MPI)
      call mpi_allreduce(gap,gap_mpi,size(gap),&
                        mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
     gap_mpi= gap
#endif
      
      outfileindex= outfileindex+ 1
      if (cpuid==0)then
         open(unit=outfileindex, file='GapCube.dat')
         write(outfileindex, '(80a16)')'kx (1/A)', 'ky (1/A)', 'kz (1/A)', &
            'Energy gap', 'Ev', 'Ec', &
            'k1 (2pi/a)', 'k2 (2pi/b)', 'k3 (2pi/c)'
         do ik=1, knv3
            if (abs(gap_mpi(1, ik))< Gap_threshold) then
               write(outfileindex, '(80f16.8)') kxy_shape(:, ik)*Angstrom2atomic, &
                  (gap_mpi(:, ik)/eV2Hartree), kxy(:, ik)*Angstrom2atomic
            endif
         enddo
         close(outfileindex)
      endif
      
     !!> minimum and maximum value of energy bands
     !
     !zmax= maxval((gap_mpi))
     !zmin= minval((gap_mpi))
      
      !> write script for gnuplot
      outfileindex= outfileindex+ 1
      if (cpuid==0) then
         open(unit=outfileindex, file='GapCube.gnu')
         write(outfileindex, '(a)')"set encoding iso_8859_1"
         write(outfileindex, '(a)')'#set terminal  postscript enhanced color'
         write(outfileindex, '(a)')"#set output 'gap.eps'"
         write(outfileindex, '(3a)')'set terminal  png      truecolor enhanced', &
            ' size 1920, 1680 font ",36"'
         write(outfileindex, '(a)')"set output 'GapCube.png'"
         write(outfileindex, '(a)')'unset ztics'
         write(outfileindex, '(a)')'unset key'
         write(outfileindex, '(a)')'set ticslevel 0'
         write(outfileindex, '(a)')'#set view equal xyz'
         write(outfileindex, '(a)')'set border lw 3'
         write(outfileindex, '(a)')'set xlabel "k_x"'
         write(outfileindex, '(a)')'set ylabel "k_y"'
         write(outfileindex, '(a)')'set zlabel "k_z"'
         write(outfileindex, '(a)')'set xtics nomirror scale 0.5'
         write(outfileindex, '(a)')'set ytics nomirror scale 0.5'
         write(outfileindex, '(a)')'set ztics nomirror scale 0.5'
         write(outfileindex, '(a)')'set xtics offset  0.0,-1.0 , 0'
         write(outfileindex, '(a)')'set ytics offset -2.5,   0 , 0'
         write(outfileindex, '(a)')'set size ratio -1'
         write(outfileindex, '(a)')'set view 60, 140, 1, 1'
         write(outfileindex, '(a, f10.5, a, f10.5, a)')'set xrange [', kxmin_shape*Angstrom2atomic, ':', kxmax_shape*Angstrom2atomic, ']'
         write(outfileindex, '(a, f10.5, a, f10.5, a)')'set yrange [', kymin_shape*Angstrom2atomic, ':', kymax_shape*Angstrom2atomic, ']'
         write(outfileindex, '(a, f10.5, a, f10.5, a)')'set zrange [', kzmin_shape*Angstrom2atomic, ':', kzmax_shape*Angstrom2atomic, ']'
         write(outfileindex, '(2a)')"splot 'GapCube.dat' u 1:2:3 w p pt 7 ps 2"
         close(outfileindex)
     
      endif
      
      
      return
   end subroutine gapshape3D


   subroutine gapshape
      ! This subroutine gets the gap at each k 
      ! points on the k plane specified by user

      use wmpi
      use para
      
      implicit none
      
      integer :: ik, i, j
      integer :: knv3
      integer :: nkx
      integer :: nky
      
      integer :: ierr, i1, i2
      real(Dp) :: k(3)
      
      ! Hamiltonian of bulk system
      complex(Dp), allocatable :: Hamk_bulk(:, :)
      
      real(dp) :: zmin, zmax
      real(dp) :: kxmin_shape, kxmax_shape, kymin_shape, kymax_shape
      
      real(dp), allocatable :: kxy(:,:)
      real(dp), allocatable :: kxy_shape(:,:)
      real(dp), allocatable :: kxy_plane(:,:)
      
      real(dp), allocatable :: gap(:, :)
      real(dp), allocatable :: gap_mpi(:, :)
      real(dp), allocatable :: W(:)
      
      complex(dp), allocatable :: ones(:,:)
     
      allocate(Hamk_bulk(Num_wann, Num_wann))

      nkx= Nk1
      nky= Nk2
      allocate( kxy(3, nkx*nky))
      allocate( kxy_shape(3, nkx*nky))
      allocate( kxy_plane(3, nk1*Nk2))
      kxy=0d0
      kxy_shape=0d0
      
     
      ik =0
      do i= 1, nkx
         do j= 1, nky
            ik =ik +1
            kxy(:, ik)= K3D_start+ K3D_vec1*(i-1)/dble(nkx-1)+ K3D_vec2*(j-1)/dble(nky-1) &
               - (K3D_vec1+ K3D_vec2)/2d0
            kxy_shape(:, ik)= kxy(1, ik)* Origin_cell%Kua+ kxy(2, ik)* Origin_cell%Kub+ kxy(3, ik)* Origin_cell%Kuc 
            call rotate_k3_to_kplane(kxy_shape(:, ik), kxy_plane(:, ik))
         enddo
      enddo

      i1=2
      i2=3
      kymin_shape=minval(kxy_shape(i2,:))
      kymax_shape=maxval(kxy_shape(i2,:))
      kxmin_shape=minval(kxy_shape(i1,:))
      kxmax_shape=maxval(kxy_shape(i1,:))
      
      
      knv3= nkx*nky
      allocate( gap    ( 9, knv3))
      allocate( gap_mpi( 9, knv3))
      gap    = 0d0
      gap_mpi= 0d0
      
      allocate(W(Num_wann))
      allocate(ones(Num_wann, Num_wann))
      W= 0d0
      ones= 0d0
      do i=1, Num_wann
         ones(i, i)= 1d0
      enddo
      
      if (Numoccupied> Num_wann) then
         stop 'Numoccupied should less than Num_wann'
      endif
      
      do ik= 1+cpuid, knv3, num_cpu
         if (cpuid==0)write(stdout, *)'Gap plane', ik, knv3
      
         k(1) = kxy(1, ik)
         k(2) = kxy(2, ik)
         k(3) = kxy(3, ik)
      
         !> calculation bulk hamiltonian
         Hamk_bulk= 0d0
         if (index(KPorTB, 'KP')/=0)then
           call ham_bulk_kp(k, Hamk_bulk)
        else
          !> deal with phonon system
          if (index(Particle,'phonon')/=0.and.LOTO_correction) then
             call ham_bulk_LOTO(k, Hamk_bulk)
          else
             call ham_bulk_latticegauge    (k, Hamk_bulk)
          endif
        endif

     
         !> diagonalization
         call eigensystem_c( 'N', 'U', Num_wann ,Hamk_bulk, W)
         gap(1, ik)= W(Numoccupied+1)- W(Numoccupied)
         gap(2:9 , ik)= W(Numoccupied-3:Numoccupied+4)
      
      enddo
      
#if defined (MPI)
      call mpi_allreduce(gap,gap_mpi,size(gap),&
                        mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
     gap_mpi= gap
#endif
      
      outfileindex= outfileindex+ 1 
      if (cpuid==0)then
         open(unit=outfileindex, file='GapPlane.dat')
     
         write(outfileindex, '(100a16)')'# kx', 'ky', 'kz', "kx'", "ky'", "kz'", 'gap', 'Ev4', 'Ev3', &
            'Ev2', 'Ev1', 'Ec1', 'Ec2', 'Ec3', 'Ec4', 'k1', 'k2', 'k3'
         do ik=1, knv3
            write(outfileindex, '(300f16.8)')kxy_shape(:, ik)*Angstrom2atomic, &
               kxy_plane(:, ik)*Angstrom2atomic, (gap_mpi(:, ik)/eV2Hartree), kxy(:, ik)*Angstrom2atomic
            if (mod(ik, nky)==0) write(outfileindex, *)' '
         enddo
         close(outfileindex)
      endif

      outfileindex= outfileindex+ 1 
      if (cpuid==0)then
         open(unit=outfileindex, file='gap2d.dat')
         write(outfileindex, '(100a16)')'% kx', 'ky', 'kz', "kx'", "ky'", "kz'", &
            'gap', 'Ev2', 'Ev1', 'Ec1', &
            'Ec2', 'k1', 'k2', 'k3'
         do ik=1, knv3
            if (abs(gap_mpi(1, ik))< Gap_threshold) then
               write(outfileindex, '(800f16.8)')kxy_shape(:, ik)*Angstrom2atomic, &
               kxy_plane(:, ik)*Angstrom2atomic, (gap_mpi(:, ik)/eV2Hartree), kxy(:, ik)*Angstrom2atomic
            endif
         enddo
         close(outfileindex)
      endif
 
      outfileindex= outfileindex+ 1 
      if (cpuid==0)then
         open(unit=outfileindex, file='GapPlane_matlab.dat')
         write(outfileindex, '(100a16)')'% kx', 'ky', 'kz', 'gap', 'Ev2', 'Ev1', 'Ec1', &
            'Ec2', 'k1', 'k2', 'k3'
         do ik=1, knv3
            write(outfileindex, '(300f16.8)')kxy_shape(:, ik)*Angstrom2atomic, &
               kxy_plane(:, ik)*Angstrom2atomic, (gap_mpi(:, ik)/eV2Hartree), kxy(:, ik)*Angstrom2atomic
         enddo
         close(outfileindex)
      endif
      
      !> minimum and maximum value of energy bands
      
      zmax= maxval(gap_mpi(1, :))
      zmin= minval(gap_mpi(1, :))
      
      !> write script for gnuplot
      outfileindex= outfileindex+ 1 
      if (cpuid==0) then
         open(unit=outfileindex, file='GapPlane.gnu')
         write(outfileindex, '(a)')"set encoding iso_8859_1"
         write(outfileindex, '(a)')'#set terminal  postscript enhanced color'
         write(outfileindex, '(a)')"#set output 'GapPlane.eps'"
         write(outfileindex, '(3a)')'#set terminal  pngcairo   truecolor enhanced', &
            '  size 1920, 1680 font ",60"'
         write(outfileindex, '(3a)')'set terminal  png   truecolor enhanced', &
            ' size 1920, 1680 font ",60"'
         write(outfileindex, '(a)')"set output 'GapPlane.png'"
         write(outfileindex,'(a, f10.4, a, f10.4, a, f10.4, a)') &
            'set palette defined ( ', zmin, ' "black", ', &
            (zmin+zmax)/20d0,' "orange", ',zmax,'  "white" )'
         write(outfileindex, '(a)')"set origin 0.10, 0.0"
         write(outfileindex, '(a)')"set size 0.85, 1.0"
         write(outfileindex, '(a)')'unset ztics'
         write(outfileindex, '(a)')'unset key'
         write(outfileindex, '(a)')'set pm3d'
         write(outfileindex, '(a)')'set view map'
         write(outfileindex, '(a)')'set border lw 3'
         write(outfileindex, '(a)')'#set size ratio -1'
         write(outfileindex, '(a)')'set title "Gap in k plane"'
         write(outfileindex, '(a)')'set xtics nomirror scale 0.5'
         write(outfileindex, '(a)')'set ytics nomirror scale 0.5'
         write(outfileindex, '(a)')"set xlabel 'k (1/{\305})'"
         write(outfileindex, '(a)')"set ylabel 'k (1/{\305})'"
         write(outfileindex, '(a)')'set colorbox'
         write(outfileindex, '(a)')'set xrange [ ] noextend'
         write(outfileindex, '(a)')'set yrange [ ] noextend'
         write(outfileindex, '(a)')'set pm3d interpolate 2,2'
         write(outfileindex, '(2a)')"splot 'GapPlane.dat' u 4:5:7 w pm3d"
     
         close(outfileindex)
      endif
      
      return
   end subroutine gapshape


   subroutine get_fermilevel
      !> Calculate fermilevel for the given hamiltonian
      use wmpi
      use para
      implicit none

      integer :: ikx, iky, ikz, io, ik, ik_first, ik_last

      !> number of k points
      integer :: knv3, ierr, iter, itermax, ibeta

      !> fermi level
      real(dp) :: EF, k(3)

      real(dp) ::  Beta_fake, lmin0, lmax0, lmin, lmax, tot, tot_mpi, lmin_mpi, lmax_mpi

      !> fermi-dirac distribution function
      real(dp), external :: fermi

      !> eigen value for each kpoint
      real(dp), allocatable :: W(:)
      real(dp), allocatable :: eigvals(:, :)

      !> we calculate Fermi level at different temperature
      integer :: Beta_num
      real(dp), allocatable :: Beta_array(:), EF_array(:)

      complex(dp), allocatable :: ham(:, :)

      knv3= Nk1*Nk2*Nk3
      call WTGenerateLocalPartition(knv3, num_cpu, cpuid, ik_first, ik_last)

      allocate(W(Num_wann))
      allocate(eigvals(Num_wann, ik_first:ik_last))
      allocate(ham(Num_wann, Num_wann))
      eigvals= 0d0
      ham= 0d0
      Beta_fake = Beta
      Beta_num=30
      allocate(Beta_array(Beta_num), EF_array(Beta_num))
      Beta_array= 0d0
      EF_array= 0d0

      do ibeta=1, Beta_num
         Beta_array(ibeta)= 11600d0/(10d0+ (300d0-10d0)*(ibeta-1)/(Beta_num-1))
      enddo

      do ik=ik_first, ik_last
 
         ikx= (ik-1)/(nk2*nk3)+1
         iky= ((ik-1-(ikx-1)*Nk2*Nk3)/nk3)+1
         ikz= (ik-(iky-1)*Nk3- (ikx-1)*Nk2*Nk3)
         k= K3D_start_cube+ K3D_vec1_cube*(ikx-1)/dble(nk1-1)  &
          + K3D_vec2_cube*(iky-1)/dble(nk2-1)  &
          + K3D_vec3_cube*(ikz-1)/dble(nk3-1)

         ham= 0d0
         call ham_bulk_latticegauge(k, ham)
         call eigensystem_c( 'N', 'U', num_wann, ham, W)
         eigvals(:, ik)= W
      enddo ! ik

      ! using bisection algorithm to search the fermi level
      iter= 0 
      itermax= 100
      tot= 9999d0
      lmin_mpi= minval(eigvals)
      lmax_mpi= maxval(eigvals)

#if defined (MPI)
         call mpi_allreduce(lmin_mpi, lmin, 1, &
                            mpi_dp, mpi_min, mpi_cmw, ierr)
         call mpi_allreduce(lmax_mpi, lmax, 1, &
                            mpi_dp, mpi_max, mpi_cmw, ierr)
#else
         lmin= lmin_mpi
         lmax= lmax_mpi
#endif
 

      if (cpuid==0) write(stdout, *)' Lowest energy level in the whold energy bands', lmin
      if (cpuid==0) write(stdout, *)' highest energy level in the whold energy bands', lmax
      lmin0= lmin
      lmax0= lmax

      do ibeta= 1, Beta_num
         Beta_fake= Beta_array(ibeta)/eV2Hartree
         lmin= lmin0
         lmax= lmax0
         tot= 9999d0
         iter = 0
         
         if (cpuid==0) then
             write(stdout, '(a,f12.6,a,f12.6,a)') ' Beta :' , Beta_fake*eV2Hartree, ' T: ', 11600d0/Beta_fake/eV2Hartree, ' Kelvin'
         endif
         do while( abs(tot- Ntotch).gt. eps6 .and. iter.lt.itermax)
         
            iter= iter+ 1
         
            EF= (lmin+ lmax)* half
         
            tot_mpi= 0d0
            do ik=ik_first, ik_last
               do io=1, Num_wann
                  tot_mpi= tot_mpi+ fermi(eigvals(io, ik)- EF, Beta_fake)
               enddo ! io
            enddo ! ik
         
            tot = 0d0
#if defined (MPI)
            call mpi_allreduce(tot_mpi, tot, 1, &
                               mpi_dp, mpi_sum, mpi_cmw, ierr)
#else   
            tot= tot_mpi
#endif   
         
            tot= tot/dble(knv3)
         
            if (SOC==0) then
               tot= tot*2
            endif
         
            !> bisection
            if (tot > Ntotch)then
               lmax= EF
            else
               lmin= EF
            endif
         
            if (cpuid==0) then
                write(stdout, 100)iter, tot-Ntotch, EF/eV2Hartree, '  Charge: ', tot
            endif
         100   format(2x,">iter",i4,2x,"diff:",f12.6,2x,"EF: ",f12.6,a,f12.6)
         
         enddo ! bisection
         EF_array( ibeta) = EF
      enddo ! ibeta

      E_fermi= EF_array(1)

      if (cpuid==0) write(stdout, '(a,f16.6)')" >>Fermi level we found by bisection method  at different temperature: "

      if (cpuid==0) then
         write(stdout, '(a, f12.6)') 'Number of electrons : ', Ntotch
         write(stdout, '(3a12)') 'Beta    ', ' T (Kelvin)  ', ' E_F (eV)   '
         do ibeta=1, Beta_num
            write(stdout, '(3f12.6)') Beta_array(ibeta), 11600d0/Beta_array(ibeta), EF_array(ibeta)/eV2Hartree
         enddo
         write(stdout, '(3a12)') ' '
      endif

      deallocate(W)
      deallocate(eigvals)
      deallocate(ham)
      return
   end subroutine get_fermilevel

   subroutine get_density
      !> Calculate fermilevel for the given hamiltonian
      use wmpi
      use para
      implicit none

      integer :: ikx, iky, ikz, ik, iT, ie, iwan, ierr, ieta

      !> number of k points
      integer :: knv3
      integer :: NumberofEta

      !> fermi level
      real(dp) :: k(3)

      real(dp), allocatable :: occupation_fermi(:), occupation_mu(:, :)
      real(dp), allocatable :: occupation_fermi_mpi(:), occupation_mu_mpi(:, :)
      real(dp), allocatable :: density(:)
      real(dp) :: time_start, time_end

      !> fermi-dirac distribution function
      real(dp), external :: fermi

      !> eigen value for each kpoint
      real(dp), allocatable :: W(:)
      real(dp), allocatable :: eigvals(:, :), eigvals_mpi(:, :)

      !> we calculate Fermi level at different temperature
      integer :: Beta_num
      real(dp), allocatable :: mu_array(:), KBT_array(:), eta_array(:)

      complex(dp), allocatable :: ham(:, :)
      character(40) :: etaname
      character(40), allocatable :: etanamelist(:)

      knv3= Nk1*Nk2*Nk3
      ! call WTGenerateLocalPartition(knv3, num_cpu, cpuid, ik_first, ik_last)
      NumberofEta = 9 

      allocate(eta_array(NumberofEta))
      allocate(W(Num_wann))
      ! allocate(eigvals(Num_wann, ik_first:ik_last))
      ! allocate(eigvals_mpi(Num_wann, ik_first:ik_last))
      allocate(ham(Num_wann, Num_wann))
      allocate(KBT_array(NumT))
      allocate(mu_array(OmegaNum))
      allocate(occupation_fermi(NumberofEta))
      allocate(occupation_fermi_mpi(NumberofEta))
      allocate(occupation_mu(OmegaNum, NumT))
      allocate(occupation_mu_mpi(OmegaNum, NumT))
      allocate(density(NumberofEta))
      allocate(etanamelist(NumberofEta))
      occupation_fermi= 0d0
      occupation_mu= 0d0
      occupation_fermi_mpi= 0d0
      occupation_mu_mpi= 0d0
      mu_array= 0d0
      KBT_array= 0d0
      ! eigvals= 0d0
      ham= 0d0
      W = 0d0
      density = 0d0
      
      eta_array=(/0.1d0, 0.2d0, 0.4d0, 0.8d0, 1.0d0, 2d0, 4d0, 8d0, 10d0/)
      eta_array= eta_array*Fermi_broadening

      if (NumT>1) then
         do iT=1, NumT
            KBT_array(iT)= Tmin+ (iT-1.0d0)/(NumT-1.0d0)*(Tmax-Tmin)
         enddo
      else
         KBT_array = Tmin
      endif

      if (cpuid.eq.0) then
         write(stdout, *) ' '
         write(stdout, *)' KBT array in the calculation in unit of Kelvin'
         write(stdout, '(10f8.2)') KBT_array
         write(stdout, *) ' '
      endif
   
      !> transform from Kelvin to eV
      !> The SI unit of temperature is the kelvin (K), but using the above relation the electron temperature is often expressed in
      !> terms of the energy unit electronvolt (eV). Each kelvin (1 K) corresponds to 8.6173324(78)×10−5 eV; this factor is the ratio
      !> of the Boltzmann constant to the elementary charge. After version 2.6, we 
      !> adopt the atomic unit
      KBT_array= KBT_array*8.6173324E-5*eV2Hartree

      if (OmegaNum>1) then
         do ie=1, OmegaNum
            mu_array(ie)= OmegaMin+ (ie-1.0d0)/(OmegaNum-1.0d0)*(OmegaMax-OmegaMin)
         enddo
      else
         mu_array = OmegaMin
      endif

      time_start= 0d0
      time_end= 0d0
      do ik= 1+cpuid, knv3, num_cpu
      ! do ik= ik_first, ik_last
         if (cpuid==0.and. mod(ik/num_cpu, 500)==0) &
           write(stdout, '(a,I16,a,I16,a,f16.3,a)') 'occupation, ik ', ik, '/',knv3, 'time left', &
           (knv3-ik)*(time_end- time_start)/num_cpu, ' s'
        call now(time_start)
        
         ikx= (ik-1)/(nk2*nk3)+1
         iky= ((ik-1-(ikx-1)*Nk2*Nk3)/nk3)+1
         ikz= (ik-(iky-1)*Nk3- (ikx-1)*Nk2*Nk3)
         k= K3D_start_cube+ K3D_vec1_cube*(ikx-1)/dble(nk1-1)  &
          + K3D_vec2_cube*(iky-1)/dble(nk2-1)  &
          + K3D_vec3_cube*(ikz-1)/dble(nk3-1)
 
 
         ! calculation bulk hamiltonian
         ham= 0d0
        !call ham_bulk_atomicgauge    (k, ham)
         call ham_bulk_latticegauge    (k, ham)
         call eigensystem_c( 'N', 'U', Num_wann, ham, W)
         ! eigvals_mpi(:, ik)= W

         do iT = 1, NumT
            do iwan = 1, Num_wann
               do ie = 1, OmegaNum
                  occupation_mu_mpi(ie, iT)= occupation_mu_mpi(ie, iT)+ &
                     fermi(W(iwan)-mu_array(ie), 1/KBT_array(iT))
               enddo
            enddo
         enddo
         
         do iwan = 1, Num_wann
            do ieta = 1, NumberofEta
               occupation_fermi_mpi(ieta)= occupation_fermi_mpi(ieta)+ &
                  fermi(W(iwan), 1/eta_array(ieta))
            enddo
         enddo
         call now(time_end)
      enddo

 
#if defined (MPI)
   ! call mpi_allreduce(eigvals_mpi, eigvals,size(eigvals),&
                     ! mpi_dp,mpi_sum,mpi_cmw,ierr)
   call mpi_allreduce(occupation_fermi_mpi, occupation_fermi, &
                     size(occupation_fermi),mpi_dp,mpi_sum,mpi_cmw,ierr)
   call mpi_allreduce(occupation_mu_mpi, occupation_mu, &
                     size(occupation_mu),mpi_dp,mpi_sum,mpi_cmw,ierr)
   if (ierr>0) then
      print *, 'Something wrong in mpi_allreduce in get_occupation', ierr
      stop
   endif
#else
   ! eigvals= eigvals_mpi 
   occupation_fermi= occupation_fermi_mpi
   occupation_mu = occupation_mu_mpi
#endif
   if (cpuid==0) then
      write(stdout, *)'>> All processors finished their job for get_density' 
   endif

   !> calculate the occupation in the unit of cm^-3
      occupation_fermi= occupation_fermi/knv3/Origin_cell%CellVolume/kCubeVolume &
                        *Origin_cell%ReciprocalCellVolume/(Bohr_radius*100)**3
      occupation_mu= occupation_mu/knv3/Origin_cell%CellVolume/kCubeVolume &
                        *Origin_cell%ReciprocalCellVolume/(Bohr_radius*100)**3

      do ieta = 1, NumberofEta
         write(etaname,'(f12.2)') eta_array(ieta)*1000d0/eV2Hartree
         write(etanamelist(ieta),'(3a)') 'occ@eta_',trim(adjustl(etaname)),'meV'
      enddo
      outfileindex = outfileindex+1
      open(unit=outfileindex, file='density.dat', status='unknown')
      write(outfileindex,'(a)') '# density in the unit of cm^-3 at different temperatures and chemical potentials' 
      write(outfileindex, '(a,1000f8.3)')'# Tlist  =  ', KBT_array(:)/8.6173324E-5/eV2Hartree
      write(outfileindex, '(a,1000f8.3)')'# mulist  =  ', mu_array(:)/eV2Hartree
      do iT = 1, NumT
         write(outfileindex, '(a, f9.3, a)')'# T  =  ', KBT_array(iT)/8.6173324E-5/eV2Hartree, ' K'
         write(outfileindex, '(a9,9a20)')'# mu (eV)', etanamelist
         do ie = 1, OmegaNum
            density = 0d0
            do ieta = 1, NumberofEta
               density(ieta) = occupation_mu(ie, iT)-occupation_fermi(ieta)
            enddo
            write(outfileindex, '(f9.3, 9E20.6)') mu_array(ie)/eV2Hartree, density
         enddo
         write(outfileindex, '(a)')''
      enddo

      return
   end subroutine get_density
   
   function fermi(omega, Beta_fake) result(value)
      ! This function sets the Fermi-Dirac distribution

      use para
      implicit none

      ! >> inout variables
      real(dp), intent(in) :: omega
      real(dp), intent(in) :: Beta_fake

      ! return value
      real(dp) :: value
    
      ! avoid numerical instability 
      if (Beta_fake*omega .ge. 20d0) then
         value = zero
      elseif (Beta_fake*omega.le. -20d0)then
         value = one
      else
         value= one/(one+exp(Beta_fake*omega))
      endif

      return
   end function fermi


   function dfde(omega, Beta_fake) result(value)
      ! This function sets the Fermi-Dirac distribution

      use para
      implicit none

      ! >> inout variables
      real(dp), intent(in) :: omega
      real(dp), intent(in) :: Beta_fake

      ! return value
      real(dp) :: value
    
      ! avoid numerical instability 
      if (Beta_fake*omega .ge. 20d0) then
         value = zero
      elseif (Beta_fake*omega.le. -20d0)then
         value = one
      else
         value= one/(one+exp(Beta_fake*omega))
      endif

      return
   end function dfde

   subroutine rotate_k3_to_kplane_mag(k3, kplane)
      use para, only : dp, K3D_vec1, K3D_vec2, Kua_mag, Kub_mag, Kuc_mag
      implicit none

      real(dp), intent(in) :: k3(3)
      real(dp), intent(out) :: kplane(3)
      !> three new unit vectors
      real(dp), allocatable :: Urot_t(:, :)

      real(dp) :: kvec1(3), kvec2(3)

      real(dp), external :: norm

      allocate(Urot_t(3, 3))

      kvec1= K3D_vec1(1)*Kua_mag+ K3D_vec1(2)*Kub_mag+ K3D_vec1(3)*Kuc_mag
      kvec2= K3D_vec2(1)*Kua_mag+ K3D_vec2(2)*Kub_mag+ K3D_vec2(3)*Kuc_mag

      Urot_t(1, :)= kvec1/norm(kvec1)

      !> e_z'
      Urot_t(3, 1)= (kvec1(2)*kvec2(3)- kvec1(3)*kvec2(2))
      Urot_t(3, 2)= (kvec1(3)*kvec2(1)- kvec1(1)*kvec2(3))
      Urot_t(3, 3)= (kvec1(1)*kvec2(2)- kvec1(2)*kvec2(1))
      Urot_t(3, :)= Urot_t(3, :)/norm(Urot_t(3, :))
 
      !> e_y'= e_z'\cross e_x'
      Urot_t(2, 1)= (Urot_t(3, 2)*Urot_t(1, 3)- Urot_t(3, 3)*Urot_t(1, 2))
      Urot_t(2, 2)= (Urot_t(3, 3)*Urot_t(1, 1)- Urot_t(3, 1)*Urot_t(1, 3))
      Urot_t(2, 3)= (Urot_t(3, 1)*Urot_t(1, 2)- Urot_t(3, 2)*Urot_t(1, 1))
      Urot_t(2, :)= Urot_t(2, :)/norm(Urot_t(2, :))

      kplane(1)= Urot_t(1, 1)*k3(1)+ Urot_t(1, 2)*k3(2)+ Urot_t(1, 3)*k3(3)
      kplane(2)= Urot_t(2, 1)*k3(1)+ Urot_t(2, 2)*k3(2)+ Urot_t(2, 3)*k3(3)
      kplane(3)= Urot_t(3, 1)*k3(1)+ Urot_t(3, 2)*k3(2)+ Urot_t(3, 3)*k3(3)

      return
   end subroutine  rotate_k3_to_kplane_mag


   subroutine rotate_k3_to_kplane(k3, kplane)
      use para, only : dp, K3D_vec1, K3D_vec2, Origin_cell
      implicit none

      real(dp), intent(in) :: k3(3)
      real(dp), intent(out) :: kplane(3)
      !> three new unit vectors
      real(dp), allocatable :: Urot_t(:, :)

      real(dp) :: kvec1(3), kvec2(3)

      real(dp), external :: norm

      allocate(Urot_t(3, 3))

      kvec1= K3D_vec1(1)*Origin_cell%Kua+ K3D_vec1(2)*Origin_cell%Kub+ K3D_vec1(3)*Origin_cell%Kuc
      kvec2= K3D_vec2(1)*Origin_cell%Kua+ K3D_vec2(2)*Origin_cell%Kub+ K3D_vec2(3)*Origin_cell%Kuc

      Urot_t(1, :)= kvec1/norm(kvec1)

      !> e_z'
      Urot_t(3, 1)= (kvec1(2)*kvec2(3)- kvec1(3)*kvec2(2))
      Urot_t(3, 2)= (kvec1(3)*kvec2(1)- kvec1(1)*kvec2(3))
      Urot_t(3, 3)= (kvec1(1)*kvec2(2)- kvec1(2)*kvec2(1))
      Urot_t(3, :)= Urot_t(3, :)/norm(Urot_t(3, :))
 
      !> e_y'= e_z'\cross e_x'
      Urot_t(2, 1)= (Urot_t(3, 2)*Urot_t(1, 3)- Urot_t(3, 3)*Urot_t(1, 2))
      Urot_t(2, 2)= (Urot_t(3, 3)*Urot_t(1, 1)- Urot_t(3, 1)*Urot_t(1, 3))
      Urot_t(2, 3)= (Urot_t(3, 1)*Urot_t(1, 2)- Urot_t(3, 2)*Urot_t(1, 1))
      Urot_t(2, :)= Urot_t(2, :)/norm(Urot_t(2, :))

      kplane(1)= Urot_t(1, 1)*k3(1)+ Urot_t(1, 2)*k3(2)+ Urot_t(1, 3)*k3(3)
      kplane(2)= Urot_t(2, 1)*k3(1)+ Urot_t(2, 2)*k3(2)+ Urot_t(2, 3)*k3(3)
      kplane(3)= Urot_t(3, 1)*k3(1)+ Urot_t(3, 2)*k3(2)+ Urot_t(3, 3)*k3(3)

      return
   end subroutine  rotate_k3_to_kplane


   subroutine rotate_kplane_to_k3(kplane, k3)
      use para, only : dp, K3D_vec1, K3D_vec2, Origin_cell
      implicit none

      real(dp), intent(in ) :: kplane(3)
      real(dp), intent(out) :: k3(3)
      
      !> three new unit vectors
      real(dp), allocatable :: Urot_t(:, :)

      real(dp) :: kvec1(3), kvec2(3)

      real(dp), external :: norm

      allocate(Urot_t(3, 3))

      kvec1= K3D_vec1(1)*Origin_cell%Kua+ K3D_vec1(2)*Origin_cell%Kub+ K3D_vec1(3)*Origin_cell%Kuc
      kvec2= K3D_vec2(1)*Origin_cell%Kua+ K3D_vec2(2)*Origin_cell%Kub+ K3D_vec2(3)*Origin_cell%Kuc

      Urot_t(1, :)= kvec1/norm(kvec1)

      !> e_z'
      Urot_t(3, 1)= (kvec1(2)*kvec2(3)- kvec1(3)*kvec2(2))
      Urot_t(3, 2)= (kvec1(3)*kvec2(1)- kvec1(1)*kvec2(3))
      Urot_t(3, 3)= (kvec1(1)*kvec2(2)- kvec1(2)*kvec2(1))
      Urot_t(3, :)= Urot_t(3, :)/norm(Urot_t(3, :))
 
      !> e_y'= e_z'\cross e_x'
      Urot_t(2, 1)= (Urot_t(3, 2)*Urot_t(1, 3)- Urot_t(3, 3)*Urot_t(1, 2))
      Urot_t(2, 2)= (Urot_t(3, 3)*Urot_t(1, 1)- Urot_t(3, 1)*Urot_t(1, 3))
      Urot_t(2, 3)= (Urot_t(3, 1)*Urot_t(1, 2)- Urot_t(3, 2)*Urot_t(1, 1))
      Urot_t(2, :)= Urot_t(2, :)/norm(Urot_t(2, :))

      call inv_r(3, Urot_t)

      k3(1)= Urot_t(1, 1)*kplane(1)+ Urot_t(1, 2)*kplane(2)+ Urot_t(1, 3)*kplane(3)
      k3(2)= Urot_t(2, 1)*kplane(1)+ Urot_t(2, 2)*kplane(2)+ Urot_t(2, 3)*kplane(3)
      k3(3)= Urot_t(3, 1)*kplane(1)+ Urot_t(3, 2)*kplane(2)+ Urot_t(3, 3)*kplane(3)

      return
   end subroutine  rotate_kplane_to_k3

   subroutine project_k3_to_kplane_defined_by_direction(k3, direction, kplane)
      !> In this subroutine, we define a new coordinate system. z axis is parallel to direction
      !> x axis is generated by the cross product of direction and [0 1 0]
      !> y axis is generated by the cross product of direction and x
      use para, only : dp, eps6
      implicit none

      real(dp), intent(in) :: k3(3)
      real(dp), intent(in) :: direction(3)
      real(dp), intent(out) :: kplane(3)

      !> three new unit vectors
      real(dp), allocatable :: Urot_t(:, :)

      real(dp) :: kvec1(3), kvec2(3), kvec3(3)

      real(dp), external :: norm

      allocate(Urot_t(3, 3))

      kvec1= direction/norm(direction)
      kvec2=(/0.d0, 1.d0, 0.d0/) 
      call cross_product(kvec1, kvec2, kvec3)
      if (norm(kvec3)<eps6) then
         kvec2=(/1.d0, 0.d0, 0.d0/) 
         call cross_product(kvec1, kvec2, kvec3)
      endif

      !> e_z'
      Urot_t(3, :)= kvec1/norm(kvec1) 

      !> e_x'
      Urot_t(1, :)= kvec3/norm(kvec3) 

      !> e_y'=e_z' x e_x'
      call cross_product(kvec1, kvec3, kvec2)
      Urot_t(2, :)= kvec2/norm(kvec2) 

      kplane(1)= Urot_t(1, 1)*k3(1)+ Urot_t(1, 2)*k3(2)+ Urot_t(1, 3)*k3(3)
      kplane(2)= Urot_t(2, 1)*k3(1)+ Urot_t(2, 2)*k3(2)+ Urot_t(2, 3)*k3(3)
      kplane(3)= Urot_t(3, 1)*k3(1)+ Urot_t(3, 2)*k3(2)+ Urot_t(3, 3)*k3(3)

      return
   end subroutine  project_k3_to_kplane_defined_by_direction


