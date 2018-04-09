  subroutine fermisurface3D
     ! This subroutine calculates 3D fermi surface in the 1st BZ using wannier TB method
     ! 
     ! Copyright (c) 2010 QuanSheng Wu. All rights reserved.

     use wmpi
     use para

     implicit none

     integer :: ik, i, j, l
     integer :: knv3
     integer :: ikx, iky, ikz

     integer :: ierr
     real(dp) :: kz
     real(Dp) :: k(3)

     integer :: nband_min
     integer :: nband_max
     integer :: nband_store

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
        nband_min= Numoccupied- 3
        nband_max= Numoccupied+ 4
     else
        nband_min= Numoccupied- 7
        nband_max= Numoccupied+ 8
     endif

     if (nband_min< 1) then
        nband_min= 1
        nband_max= 4
        if (SOC>0) nband_max = 8
     endif

     if (nband_max> Num_wann) then
        nband_max= Num_wann
        nband_min= Num_wann-4
        if (SOC>0) nband_min= max(1, Num_wann-8)
     endif

     nband_store= nband_max- nband_min+ 1

     kxmin= 0.00d0/1d0
     kxmax= 1.00d0/1d0
     kymin= 0.00d0/1d0
     kymax= 1.00d0/1d0
     kzmin= 0.00d0/1d0
     kzmax= 1.00d0/1d0
     kz= 0.0d0
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
        k= K3D_start_cube+ K3D_vec1_cube*(ikx-1)/dble(nk1)  &
         + K3D_vec2_cube*(iky-1)/dble(nk2)  &
         + K3D_vec3_cube*(ikz-1)/dble(nk3)


        ! calculation bulk hamiltonian
        Hamk_bulk= 0d0
        call ham_bulk_old(k, Hamk_bulk)
        call eigensystem_c( 'N', 'U', Num_wann ,Hamk_bulk, W)
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
        write(outfileindex,'(3f16.8)') (Kua(i), i=1,3)
        write(outfileindex,'(3f16.8)') (Kub(i), i=1,3)
        write(outfileindex,'(3f16.8)') (Kuc(i), i=1,3)
        do i=1,nband_store
           write(outfileindex,'(a, i10)') 'BAND: ',i
           do ik=1, knv3
              write(outfileindex,'(E16.8)') eigval(i, ik)
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

     integer :: ia
     integer :: ik, i, j, i1, i2
     integer :: knv3
     integer :: nkx
     integer :: nky

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
     real(dp), allocatable :: dos_selected(:)
     real(dp), allocatable :: dos_unselected(:)

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
           kxy_shape(:, ik)= kxy(1, ik)* Kua+ kxy(2, ik)* Kub+ kxy(3, ik)* Kuc 
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
     allocate( dos_selected(knv3))
     allocate( dos_unselected(knv3))
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
        call ham_bulk_old(k, Hamk_bulk)

        Hamk_bulk= (E_arc -zi* eta_arc)* ones - Hamk_bulk
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

     dos_total= 0d0
     do ik=1, knv3
        do i=1, Num_wann
           dos_total(ik)= dos_total(ik)+ dos_mpi(ik, i)
        enddo
     enddo

     dos_selected= 0d0
     do ik=1, knv3
        do i=1, NumberofSelectedOrbitals
           dos_selected(ik)= dos_selected(ik)+ dos_mpi(ik, Selected_Orbitals(i))
        enddo
        dos_unselected(ik)= dos_total(ik)- dos_selected(ik)
     enddo


     dos_max= maxval(dos_total)/4d0

     outfileindex= outfileindex+ 1
     if (cpuid==0)then
        open(unit=outfileindex, file='fs.dat')
        write(outfileindex, '(80a16, a)')'# kx', 'ky', 'kz', 'kp1', 'kp2', 'kp3', 'total A', 'A1(k,E)', 'A2(k,E)'
        do ik=1, knv3
           if (dos_total(ik)>dos_max) then
              write(outfileindex, '(3000f16.8)')kxy_shape(:, ik), kxy_plane(:, ik), &
                dos_total(ik), (dos_selected(ik))/dos_total(ik)*500d0+100, (dos_unselected(ik))/dos_total(ik)*500d0+100
           else
              write(outfileindex, '(3000f16.8)')kxy_shape(:, ik), kxy_plane(:, ik), &
                 dos_total(ik), (dos_selected(ik)), (dos_unselected(ik))
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
           write(outfileindex, '(3000f16.8)')kxy_shape(:, ik), kxy_plane(:, ik), &
                dos_total(ik), (dos_selected(ik)), (dos_unselected(ik))
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
              write(outfileindex, '(300f16.8)')kxy_shape(:, ik), kxy_plane(:, ik), &
                   (dos_selected(ik)), (dos_unselected(ik))
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
        write(outfileindex, '(a)')"     'orbitaltexture.dat' u 4:5:(0):($7/90000):($8/90000):(0)  w vec  head lw 5 lc rgb 'orange' front"

        close(outfileindex)
     endif


   return
   end subroutine orbitaltexture



  subroutine fermisurface
     ! This subroutine calculates 3D fermi surface in the a fixed k plane
     ! using wannier TB method

     use wmpi
     use para

     implicit none

     integer :: ia
     integer :: ik, i, j, i1, i2
     integer :: knv3
     integer :: nkx
     integer :: nky

     integer :: ierr
     real(dp) :: kz
     real(Dp) :: k(3)
     
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
     real(dp), allocatable :: dos_atom(:, :)
     real(dp), allocatable :: dos_atom_mpi(:, :)
     real(dp), allocatable :: dos_total(:)

     complex(dp), allocatable :: ones(:,:)

     integer, allocatable :: orbitals_end(:)

     allocate(Hamk_bulk(Num_wann, Num_wann))
     allocate(orbitals_end(Num_atoms))

     orbitals_end(1)= nprojs(1)
     do ia=2, Num_atoms
        orbitals_end(ia)= orbitals_start(ia)+ nprojs(ia)- 1
     enddo

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
           kxy_shape(:, ik)= kxy(1, ik)* Kua+ kxy(2, ik)* Kub+ kxy(3, ik)* Kuc 
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
     allocate( dos_atom    (knv3, Num_atoms))
     allocate( dos_atom_mpi(knv3, Num_atoms))
     allocate( dos    (knv3, Num_wann))
     allocate( dos_mpi(knv3, Num_wann))
     allocate( dos_total(knv3))
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

        ! calculation bulk hamiltonian
        Hamk_bulk= 0d0
        call ham_bulk_old(k, Hamk_bulk)

        Hamk_bulk= (E_arc -zi* eta_arc)* ones - Hamk_bulk
        call inv(Num_wann, Hamk_bulk)
        do i=1, Num_wann
           dos(ik, i)= aimag(Hamk_bulk(i, i))/pi
        enddo
        do ia=1, Num_atoms
           do i=orbitals_start(ia), orbitals_end(ia)
              dos_atom(ik, ia)= dos_atom(ik, ia)+ aimag(Hamk_bulk(i, i))/pi
           enddo
           if (SOC>0) then
              do i=orbitals_start(ia)+num_wann/2, orbitals_end(ia)+num_wann/2
                 dos_atom(ik, ia)= dos_atom(ik, ia)+ aimag(Hamk_bulk(i, i))/pi
              enddo
           endif
        enddo ! ia
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

     dos_total= 0d0
     do ik=1, knv3
        do i=1, NumberofSelectedOrbitals
           dos_total(ik)= dos_total(ik)+ dos_mpi(ik, Selected_Orbitals(i))
        enddo
     enddo


     outfileindex= outfileindex+ 1
     if (cpuid==0)then
        open(unit=outfileindex, file='fs.dat')
        write(outfileindex, '(7a16, a)')'# kx', 'ky', 'kz', 'kp1', 'kp2', 'kp3', '(A(k,E))', 'A(i,k,E),i=1,numatoms'
        do ik=1, knv3
           write(outfileindex, '(3000f16.8)')kxy_shape(:, ik), kxy_plane(:, ik), &
             (dos_total(ik)), (dos_atom_mpi(ik, :))
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
        open(unit=outfileindex, file='fs.gnu')
        write(outfileindex, '(a)')"set encoding iso_8859_1"
        write(outfileindex, '(a)')'#set terminal  postscript enhanced color'
        write(outfileindex, '(a)')"#set output 'fs.eps'"
        write(outfileindex, '(3a)')'set terminal  pngcairo truecolor enhanced', &
           ' size 1920, 1680 font ",36"'
        write(outfileindex, '(a)')"set output 'fs.png'"
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
        write(outfileindex, '(2a)')"splot 'fs.dat' u 4:5:(log($7)) w pm3d"

        close(outfileindex)
     endif


   return
   end subroutine fermisurface

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
               kxy_shape(:, ik)= kxy(1, ik)* Kua+ kxy(2, ik)* Kub+ kxy(3, ik)* Kuc 
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
         call ham_bulk_old(k, Hamk_bulk)
      
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
      
      if (cpuid==0)then
         open(unit=15, file='GapCube.dat')
         write(15, '(80a16)')'kx (1/A)', 'ky (1/A)', 'kz (1/A)', &
            'Energy gap', 'Ev', 'Ec', &
            'k1 (2pi/a)', 'k2 (2pi/b)', 'k3 (2pi/c)'
         do ik=1, knv3
            if (abs(gap_mpi(1, ik))< Gap_threshold) then
               write(15, '(80f16.8)') kxy_shape(:, ik), (gap_mpi(:, ik)), kxy(:, ik)
            endif
         enddo
         close(15)
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
         write(outfileindex, '(a, f10.5, a, f10.5, a)')'set xrange [', kxmin_shape, ':', kxmax_shape, ']'
         write(outfileindex, '(a, f10.5, a, f10.5, a)')'set yrange [', kymin_shape, ':', kymax_shape, ']'
         write(outfileindex, '(a, f10.5, a, f10.5, a)')'set zrange [', kzmin_shape, ':', kzmax_shape, ']'
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
      
      real(dp), allocatable :: gap(:, :)
      real(dp), allocatable :: gap_mpi(:, :)
      real(dp), allocatable :: W(:)
      
      complex(dp), allocatable :: ones(:,:)
     
      allocate(Hamk_bulk(Num_wann, Num_wann))

      nkx= Nk
      nky= Nk
      allocate( kxy(3, nkx*nky))
      allocate( kxy_shape(3, nkx*nky))
      kxy=0d0
      kxy_shape=0d0
      
     
      ik =0
      do i= 1, nkx
         do j= 1, nky
            ik =ik +1
            kxy(:, ik)= K3D_start+ K3D_vec1*(i-1)/dble(nkx-1)+ K3D_vec2*(j-1)/dble(nky-1)
            kxy_shape(:, ik)= kxy(1, ik)* Kua+ kxy(2, ik)* Kub+ kxy(3, ik)* Kuc 
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
         call ham_bulk_old(k, Hamk_bulk)
     
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
     
         write(outfileindex, '(100a16)')'% kx', 'ky', 'kz', 'gap', 'Ev4', 'Ev3', &
            'Ev2', 'Ev1', 'Ec1', 'Ec2', 'Ec3', 'Ec4', 'k1', 'k2', 'k3'
         do ik=1, knv3
            write(outfileindex, '(30f16.8)')kxy_shape(:, ik), (gap_mpi(:, ik)), kxy(:, ik)
            if (mod(ik, nky)==0) write(outfileindex, *)' '
         enddo
         close(outfileindex)
      endif

      outfileindex= outfileindex+ 1 
      if (cpuid==0)then
         open(unit=outfileindex, file='gap2d.dat')
         write(outfileindex, '(100a16)')'% kx', 'ky', 'kz', 'gap', 'Ev2', 'Ev1', 'Ec1', &
            'Ec2', 'k1', 'k2', 'k3'
         do ik=1, knv3
            if (abs(gap_mpi(1, ik))< Gap_threshold) then
               write(outfileindex, '(80f16.8)')kxy_shape(:, ik), (gap_mpi(:, ik)), kxy(:, ik)
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
            write(outfileindex, '(30f16.8)')kxy_shape(:, ik), (gap_mpi(:, ik)), kxy(:, ik)
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
         write(outfileindex, '(2a)')"splot 'GapPlane.dat' u 1:2:4 w pm3d"
     
         close(outfileindex)
      endif
      
      return
   end subroutine gapshape


   subroutine get_fermilevel
      !> Calculate fermilevel for the given hamiltonian
      use wmpi
      use para
      implicit none

      integer :: i1
      integer :: i2
      integer :: i3
      integer :: io
      integer :: ik

      !> number of k points
      integer :: knv3

      integer :: ierr
      integer :: iter
      integer :: itermax

      !> fermi level
      real(dp) :: EF

      real(dp) :: k(3)

      real(dp) ::  Beta 

      real(dp) :: lmin
      real(dp) :: lmax
      real(dp) :: tot


      !> fermi-dirac distribution function
      real(dp), external :: fermi

      !> kpoint coordinates
      real(dp), allocatable :: kpoints(:, :)

      !> eigen value for each kpoint
      real(dp), allocatable :: W(:)
      real(dp), allocatable :: eigvals(:, :)
      real(dp), allocatable :: eigvals_mpi(:, :)

      complex(dp), allocatable :: ham(:, :)

      knv3= Nk*Nk*Nk

      allocate(W(Num_wann))
      allocate(eigvals(Num_wann, knv3))
      allocate(eigvals_mpi(Num_wann, knv3))
      allocate(ham(Num_wann, Num_wann))
      allocate(kpoints(3, knv3))
      eigvals= 0d0
      eigvals_mpi= 0d0
      ham= 0d0
      kpoints= 0d0
      Beta= 200d0

      ik= 0
      do i1=1, Nk
      do i2=1, Nk
      do i3=1, Nk
         ik= ik+ 1
         kpoints(1, ik)= (i1-1d0)/dble(Nk)
         kpoints(2, ik)= (i2-1d0)/dble(Nk)
         kpoints(3, ik)= (i3-1d0)/dble(Nk)
      enddo
      enddo
      enddo

      do ik=1+ cpuid, knv3, num_cpu
 
         ham= 0d0
         k= kpoints(:, ik)
         call ham_bulk_old(k, ham)
         call eigensystem_c( 'N', 'U', num_wann, ham, W)
         eigvals_mpi(:, ik)= W
      enddo ! ik

#if defined (MPI)
      call mpi_allreduce(eigvals_mpi, eigvals, size(eigvals), &
                         mpi_dp, mpi_sum, mpi_cmw, ierr)
#else
     eigvals= eigvals_mpi
#endif

      ! using bisection algorithm to search the fermi level
      iter= 0 
      itermax= 100
      tot= 9999d0
      lmin= minval(eigvals)
      lmax= maxval(eigvals)
      if (cpuid==0) print *, 'Emin= ', lmin
      if (cpuid==0) print *, 'Emax= ', lmax
      do while( abs(tot- Ntotch).gt. eps6 .and. iter.lt.itermax)

         iter= iter+ 1

         EF= (lmin+ lmax)* half

         tot= 0d0
         do ik=1, knv3
            do io=1, Num_wann
               tot= tot+ fermi(eigvals(io, ik)- EF, Beta)
            enddo ! io
         enddo ! ik
         tot= tot/dble(knv3)

         if (tot > Ntotch)then
            lmax= EF
         else
            lmin= EF
         endif

         if (cpuid==0) then
             write(stdout, 100)iter, tot-Ntotch, EF, '  Charge: ', tot
         endif
      100   format(2x,">iter",i4,2x,"diff:",f12.6,2x,"EF: ",f12.6,a,f12.6)

      enddo ! bisection

      E_fermi= EF

      deallocate(W)
      deallocate(eigvals)
      deallocate(eigvals_mpi)
      deallocate(ham)
      deallocate(kpoints)
      return
   end subroutine get_fermilevel

   function fermi(omega, Beta) result(value)
      ! This function sets the Fermi-Dirac distribution

      use para
      implicit none

      ! >> inout variables
      real(dp), intent(in) :: omega
      real(dp), intent(in) :: Beta

      ! return value
      real(dp) :: value
    
      ! avoid numerical instability 
      if (beta*omega .ge. 20d0) then
         value = zero
      elseif (beta*omega.le. -20d0)then
         value = one
      else
         value= one/(one+exp(beta*omega))
      endif

      return
   end function fermi

   subroutine rotate_k3_to_kplane(k3, kplane)
      use para, only : dp, K3D_vec1, K3D_vec2, Kua, Kub, Kuc
      implicit none

      real(dp), intent(in) :: k3(3)
      real(dp), intent(out) :: kplane(3)
      !> three new unit vectors
      real(dp), allocatable :: Urot_t(:, :)

      real(dp) :: kvec1(3), kvec2(3)

      real(dp), external :: norm

      allocate(Urot_t(3, 3))

      kvec1= K3D_vec1(1)*Kua+ K3D_vec1(2)*Kub+ K3D_vec1(3)*Kuc
      kvec2= K3D_vec2(1)*Kua+ K3D_vec2(2)*Kub+ K3D_vec2(3)*Kuc

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
   end subroutine 


