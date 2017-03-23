  subroutine fermisurface3D
     ! This subroutine calculates 3D fermi surface in the 1st BZ using wannier TB method
     ! 
     ! Copyright (c) 2010 QuanSheng Wu. All rights reserved.

     use wmpi
     use para

     implicit none

     integer :: ik, i, j, l
     integer :: knv3
     integer :: nkx
     integer :: nky
     integer :: nkz

     integer :: ierr
     real(dp) :: kz
     real(Dp) :: k(3)

     integer :: nband_min
     integer :: nband_max
     integer :: nband_store
     
     ! Hamiltonian of bulk system
     complex(Dp) :: Hamk_bulk(Num_wann,Num_wann) 

     real(dp) :: kxmin, kxmax, kymin, kymax, kzmin, kzmax

     real(dp), allocatable :: kxyz(:,:)
     real(dp), allocatable :: W(:)
     real(dp), allocatable :: eigval(:,:)
     real(dp), allocatable :: eigval_mpi(:,:)

     ! only for output the FS3D.bxsf, we don't have to output all the bands,
     ! only consider the bands close to the Fermi level 
     nband_min= Numoccupied- 7
     nband_max= Numoccupied+ 8

     if (nband_min< 1) nband_min= 1
     if (nband_max> Num_wann) nband_max= Num_wann

     nband_store= nband_max- nband_min+ 1

     nkx= Nk1
     nky= Nk2
     nkz= Nk3
     allocate(kxyz(3, nkx*nky*nkz))
     kxyz=0d0

     kxmin= 0.00d0/1d0
     kxmax= 1.00d0/1d0
     kymin= 0.00d0/1d0
     kymax= 1.00d0/1d0
     kzmin= 0.00d0/1d0
     kzmax= 1.00d0/1d0
     kz= 0.0d0
     ik =0
     do i= 1, nkx
        do j= 1, nky
        do l= 1, nkz
           ik =ik +1
           kxyz(1, ik)=kxmin+ (i-1)*(kxmax-kxmin)/dble(nkx-1)
           kxyz(2, ik)=kymin+ (j-1)*(kymax-kymin)/dble(nky-1)
           kxyz(3, ik)=kzmin+ (l-1)*(kzmax-kzmin)/dble(nkz-1)
        enddo
        enddo
     enddo


     knv3= nkx*nky*nkz
     allocate(     W(Num_wann))
     allocate(eigval(nband_store, knv3))
     allocate(eigval_mpi(nband_store, knv3))
     eigval_mpi= 0d0
     eigval= 0d0
     do ik= 1+cpuid, knv3, num_cpu
        if (cpuid==0.and. mod(ik,100)<3)write(stdout, *)'3DFS, ik, knv3', ik, knv3

        k(1) = kxyz(1, ik)
        k(2) = kxyz(2, ik)
        k(3) = kxyz(3, ik)

        ! calculation bulk hamiltonian
        Hamk_bulk= 0d0
        call ham_bulk_old(k, Hamk_bulk)
        call eigensystem_c( 'N', 'U', Num_wann ,Hamk_bulk, W)
        eigval_mpi(:, ik)= W(nband_min:nband_max)
     enddo

#if defined (MPI)
     call mpi_allreduce(eigval_mpi, eigval,size(eigval),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
     eigval= eigval_mpi 
#endif

     outfileindex= outfileindex+ 1
     if (cpuid==0)then
        open(unit=outfileindex,FILE='FS3D.bxsf',STATUS='UNKNOWN',FORM='FORMATTED')
        write(outfileindex,*) ' BEGIN_INFO'
        write(outfileindex,*) '      #'
        write(outfileindex,*) '      # this is a Band-XCRYSDEN-Structure-File'
        write(outfileindex,*) '      # for Fermi Surface Visualisation'
        write(outfileindex,*) '      #'
        write(outfileindex,*) '      #  Launch as: xcrysden --bxsf FS3D.bxsf'
        write(outfileindex,*) '       Fermi Energy: 0'
        write(outfileindex,*) ' END_INFO'
        write(outfileindex,*) 
        write(outfileindex,*) ' BEGIN_BLOCK_BANDGRID_3D'
        write(outfileindex,*) 'from_wannier_code'
        write(outfileindex,*) ' BEGIN_BANDGRID_3D_fermi'
        write(outfileindex,*) nband_store
        write(outfileindex,*) nkx, nky, nkz
        write(outfileindex,*) '0.0 0.0 0.0'
        write(outfileindex,*) (Kua(i), i=1,3)
        write(outfileindex,*) (Kub(i), i=1,3)
        write(outfileindex,*) (Kuc(i), i=1,3)
        do i=1,nband_store
           write(outfileindex,*) 'BAND: ',i
           do ik=1, knv3
              write(outfileindex,'(2E16.8)') eigval(i, ik)
           enddo
        enddo
        write(outfileindex,*) 'END_BANDGRID_3D'
        write(outfileindex,*) ' END_BLOCK_BANDGRID_3D'
        close(outfileindex)
    
     endif

   return
   end subroutine fermisurface3D



  subroutine fermisurface
     ! This subroutine calculates 3D fermi surface in the a fixed k plane
     ! using wannier TB method

     use wmpi
     use para

     implicit none

     integer :: ik, i, j, i1, i2
     integer :: knv3
     integer :: nkx
     integer :: nky

     integer :: ierr
     real(dp) :: kz
     real(Dp) :: k(3)
     
     ! Hamiltonian of bulk system
     complex(Dp) :: Hamk_bulk(Num_wann,Num_wann) 

     real(dp) :: kxmin, kxmax, kymin, kymax
     real(dp) :: zmin, zmax
     real(dp) :: kxmin_shape, kxmax_shape, kymin_shape, kymax_shape

     real(dp), allocatable :: kxy(:,:)
     real(dp), allocatable :: kxy_shape(:,:)
     
     real(dp), allocatable :: dos(:)
     real(dp), allocatable :: dos_mpi(:)

     complex(dp), allocatable :: ones(:,:)

     nkx= Nk1
     nky= Nk2
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

     i1=1
     i2=2
     kxmin_shape=minval(kxy_shape(i1,:))
     kxmax_shape=maxval(kxy_shape(i1,:))
     kymin_shape=minval(kxy_shape(i2,:))
     kymax_shape=maxval(kxy_shape(i2,:))


     knv3= nkx*nky
     allocate( dos    (knv3))
     allocate( dos_mpi(knv3))
     dos    = 0d0
     dos_mpi= 0d0

     allocate(ones(Num_wann, Num_wann))
     ones= 0d0
     do i=1, Num_wann
        ones(i, i)= 1d0
     enddo

     do ik= 1+cpuid, knv3, num_cpu
        if (cpuid==0) write(stdout, *)'FS, ik, knv3' , ik, knv3

        k = kxy(:, ik)

        ! calculation bulk hamiltonian
        Hamk_bulk= 0d0
        call ham_bulk_old(k, Hamk_bulk)

        Hamk_bulk= (E_arc -zi* eta_arc)* ones - Hamk_bulk
        call inv(Num_wann, Hamk_bulk)
        do i=1, Num_wann
           dos(ik)= dos(ik)+ aimag(Hamk_bulk(i, i))/pi
        enddo

     enddo

#if defined (MPI)
     call mpi_allreduce(dos,dos_mpi,size(dos),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
     dos_mpi= dos
#endif

     outfileindex= outfileindex+ 1
     if (cpuid==0)then
        open(unit=outfileindex, file='fs.dat')
   
        do ik=1, knv3
           write(outfileindex, '(30f16.8)')kxy_shape(:, ik), log(dos_mpi(ik))
           if (mod(ik, nky)==0) write(outfileindex, *)' '
        enddo
        close(outfileindex)
     endif
     zmax= maxval(log(dos_mpi))
     zmin= minval(log(dos_mpi))

     !> minimum and maximum value of energy bands

     !> write script for gnuplot
     if (cpuid==0) then
        open(unit=101, file='fs.gnu')
        write(101, '(a)')"set encoding iso_8859_1"
        write(101, '(a)')'#set terminal  postscript enhanced color'
        write(101, '(a)')"#set output 'fs.eps'"
        write(101, '(3a)')'set terminal  pngcairo truecolor enhanced', &
           ' size 1920, 1680 font ",36"'
        write(101, '(a)')"set output 'fs.png'"
        write(101,'(a, f10.4, a, f10.4, a, f10.4, a)') &
           'set palette defined ( ', zmin, ' "white", ', &
         zmin+(zmin+zmax)/6d0, '"black", ', zmax,'  "red" )'
        write(101, '(a)')'#set palette rgbformulae 33,13,10'
        write(101, '(a)')'unset ztics'
        write(101, '(a)')'unset key'
        write(101, '(a)')'set pm3d'
        write(101, '(a)')'#set view equal xyz'
        write(101, '(a)')'set view map'
        write(101, '(a)')'set border lw 3'
        write(101, '(a)')'#set xtics font ",24"'
        write(101, '(a)')'#set ytics font ",24"'
        write(101, '(a)')'set size ratio -1'
        write(101, '(a)')'unset xtics'
        write(101, '(a)')'unset ytics'
        write(101, '(a)')'set colorbox'
       !write(101, '(a, f10.5, a, f10.5, a)')'set xrange [', kxmin, ':', kxmax, ']'
       !write(101, '(a, f10.5, a, f10.5, a)')'set yrange [', kymin, ':', kymax, ']'
       !write(101, '(a, f10.5, a, f10.5, a)')'set xrange [', kxmin_shape, ':', kxmax_shape, ']'
       !write(101, '(a, f10.5, a, f10.5, a)')'set yrange [', kymin_shape, ':', kymax_shape, ']'
        write(101, '(a)')'set pm3d interpolate 2,2'
        write(101, '(2a)')"splot 'fs.dat' u 1:2:4 w pm3d"

        close(101)
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
      complex(Dp) :: Hamk_bulk(Num_wann,Num_wann) 
      
      real(dp) :: zmin, zmax
      real(dp) :: kxmin_shape, kxmax_shape, kymin_shape, kymax_shape
      
      real(dp), allocatable :: kxy(:,:)
      real(dp), allocatable :: kxy_shape(:,:)
      
      real(dp), allocatable :: gap(:, :)
      real(dp), allocatable :: gap_mpi(:, :)
      real(dp), allocatable :: W(:)
      
      complex(dp), allocatable :: ones(:,:)
      
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
         call ham_bulk(k, Hamk_bulk)
     
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
            if (abs(gap_mpi(1, ik))< 0.10d0) then
               write(outfileindex, '(8f16.8)')kxy_shape(:, ik), (gap_mpi(:, ik)), kxy(:, ik)
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
         call ham_bulk(k, ham)
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



