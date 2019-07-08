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

   integer :: ik, il, i, j, knv3, ierr
   real(dp) :: emin,  emax,  k(3)
   character*40 :: filename
   real(Dp), allocatable :: W(:)

   ! Hamiltonian of bulk system
   complex(Dp), allocatable :: Hamk_bulk(:, :)

   ! eigen value of H
   real(dp), allocatable :: eigv(:,:), eigv_mpi(:,:)
   real(dp), allocatable :: weight(:,:,:), weight_mpi(:,:,:), weight_sum(:,:)

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

      k = k3points(:, ik)

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
            call ham_bulk_latticegauge(k, Hamk_bulk)
         endif
      endif

      !> diagonalization by call zheev in lapack
      W= 0d0
      call eigensystem_c('V', 'U', Num_wann ,Hamk_bulk, W)

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
      do il= 1, nk3lines
         if (il<10) then
            write(filename, '(a,i1)')'bulkek.dat-segment', il
         elseif (il>=10.and.il<100) then
            write(filename, '(a,i2)')'bulkek.dat-segment', il
         elseif (il>=100.and.il<1000) then
            write(filename, '(a,i3)')'bulkek.dat-segment', il
         endif
         outfileindex= outfileindex+ 1
         open(unit=outfileindex, file=filename)
         write(outfileindex, '(a, a6, a, a6)')'# segment between', k3line_name(il), ' and ',  k3line_name(il+1)
         write(outfileindex, '(a, 3f10.6)')'# kstart (fractional units)', k3line_start(:, il)
         write(outfileindex, '(a, 3f10.6)')'# kend   (fractional units)', k3line_end(:, il)
   
         do i=1, Num_wann
            do ik=1+(il-1)*Nk1, il*Nk1
               write(outfileindex, '(2f19.9, 10000i5)')(k3len(ik)-k3len((il-1)*Nk1+1)),eigv_mpi(i, ik)
            enddo
            write(outfileindex, *)' '
         enddo ! i
         close(outfileindex)
      enddo ! il
   endif

   outfileindex= outfileindex+ nk3lines+1
   if (cpuid==0)then
      open(unit=outfileindex, file='bulkek.dat')

      do i=1, Num_wann
         do ik=1, knv3
            write(outfileindex, '(2f19.9, 10000i5)')k3len(ik),eigv_mpi(i, ik), &
               int(weight(:, i, ik))
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
         write(outfileindex, '(2000f19.9)')k3len(ik),eigv_mpi(:, ik)
      enddo
      close(outfileindex)
   endif


   !> minimum and maximum value of energy bands
   emin=  minval(eigv_mpi)-0.5d0
   emax=  maxval(eigv_mpi)+0.5d0

   !> write script for gnuplot
   outfileindex= outfileindex+ 1
   if (cpuid==0) then
      open(unit=outfileindex, file='bulkek.gnu')
      write(outfileindex, '(a)') 'set terminal  postscript enhanced color font ",30"'
      write(outfileindex,'(2a)') '#set palette defined ( 0  "green", ', &
         '5 "yellow", 10 "red" )'
      write(outfileindex, '(a)')"set output 'bulkek.eps'"
      write(outfileindex, '(a)')'set style data linespoints'
      write(outfileindex, '(a)')'unset ztics'
      write(outfileindex, '(a)')'unset key'
      write(outfileindex, '(a)')'set pointsize 0.8'
      write(outfileindex, '(a)')'set view 0,0'
      write(outfileindex, '(a)')'set xtics font ",24"'
      write(outfileindex, '(a)')'set ytics font ",24"'
      write(outfileindex, '(a)')'set ylabel font ",24"'
      write(outfileindex, '(a)')'set ylabel offset 1.5,0'
      write(outfileindex, '(a, f10.5, a)')'set xrange [0: ', maxval(k3len), ']'
      write(outfileindex, '(a,f12.6)')'emin=', emin
      write(outfileindex, '(a,f12.6)')'emax=', emax
      if (index(Particle,'phonon')/=0) then
         write(outfileindex, '(a)')'set yrange [0: emax]'
         write(outfileindex, '(a)')'set ylabel "Frequency (THz)"'
      else
         write(outfileindex, '(a)')'set ylabel "Energy (eV)"'
         write(outfileindex, '(a)')'set yrange [ emin : emax ]'
      endif
      write(outfileindex, 202, advance="no") (k3line_name(i), k3line_stop(i), i=1, nk3lines)
      write(outfileindex, 203)k3line_name(nk3lines+1), k3line_stop(nk3lines+1)

      do i=1, nk3lines-1
         if (index(Particle,'phonon')/=0) then
            write(outfileindex, 204)k3line_stop(i+1), '0.0', k3line_stop(i+1), 'emax'
         else
            write(outfileindex, 204)k3line_stop(i+1), 'emin', k3line_stop(i+1), 'emax'
         endif
      enddo
      write(outfileindex, '(2a)')"plot 'bulkek.dat' u 1:2 ",  &
         " w lp lw 2 pt 7  ps 0.2 lc rgb 'black', 0 w l lw 2"
      close(outfileindex)
   endif

202 format('set xtics (',20('"',A3,'" ',F10.5,','))
203 format(A3,'" ',F10.5,')')
204 format('set arrow from ',F10.5,',',A5,' to ',F10.5,',',A5, ' nohead')

   deallocate(W)
   deallocate(Hamk_bulk)
   deallocate( eigv    )
   deallocate( eigv_mpi)
   deallocate( weight    )
   deallocate( weight_mpi)
   deallocate( weight_sum)

   return
end subroutine ek_bulk_line

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

     weight= weight_mpi/maxval(weight_mpi)*255d0

     outfileindex= outfileindex+ 1
     if (cpuid==0)then
        open(unit=outfileindex, file='bulkek-pointsmode.dat')
   
        do ik=1, Nk3_point_mode
           write(outfileindex, '(a, i10)')'# No. of k point', ik
           write(outfileindex, '("#",a9,100a10)')'k1', 'k2', 'k3', 'kx', 'ky', 'kz'  
           write(outfileindex, '(100f10.6)')k3points_pointmode_direct(:, ik), k3points_pointmode_cart(:, ik)
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

      ! calculation bulk hamiltonian
      Hamk_bulk= 0d0
      call ham_bulk_atomicgauge(k, Hamk_bulk)

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
      open(unit=outfileindex, file='bulkek_plane.dat')
      write(outfileindex, '(1000a19)')'# kx', 'ky', 'kz', 'k1', 'k2', 'k3', &
         'E(Numoccupied-1)', 'E(Numoccupied)' , 'E(Numoccupied+1)', 'E(Numoccupied+2)'
      do ik=1, knv3
         write(outfileindex, '(1000f19.9)')kxy_shape(:, ik), &
            kxy_plane(:, ik), eigv_mpi(:, ik)  
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
         write(outfileindex, '(1000f19.9)')kxy_shape(:, ik), &
            kxy_plane(:, ik), eigv_mpi(:, ik)  
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
           write(outfileindex, '(1000f19.9)')kxyz(:, ik), &
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
              write(outfileindex, '(2000f16.8)') kxy_shape(:, ik), kxy_plane(:, ik), &
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
           write(outfileindex, '(1000f16.8)')kxy_shape(:, ik), kxy_plane(:, ik), & 
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
     real(dp) :: kxmin, kxmax, kymin, kymax
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
     print *, 'nwann', nwann

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
           write(outfileindex, '(1000f19.9)')kxy_shape(:, ik), eigv_mpi(:, ik), &
              (spin_mpi(:, ib, ik), ib=Numoccupied-1, Numoccupied+2)
        enddo
        close(outfileindex)
     endif

     outfileindex= outfileindex+ 1
     if (cpuid==0)then
        open(unit=outfileindex, file='gap2D.dat')
        do ik=1, knv3
           write(outfileindex, '(1000f19.9)')kxy(:, ik), gap_mpi(ik)
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

        !> diagonalization by call zheev in lapack
        Hamk_bulk= Hamk_bulk+ zi*eta_arc
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
           write(outfileindex, '(1000f20.6)')kxy_shape(:, ik), eigv_mpi(:, ik) 
        enddo
        close(outfileindex)
     endif

     outfileindex= outfileindex+ 1
     if (cpuid==0)then
        open(unit=outfileindex, file='gap2D.dat')
        do ik=1, knv3
           write(outfileindex, '(1000f20.6)')kxy_shape(:, ik), gap_mpi(ik)
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
                 write(outfileindex, '(10000f20.5)')kxy_shape(:, ik), dos_mpi(ik) , spin_mpi(:, ik) 
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

   do ik= 1+cpuid, knv3, num_cpu
      if (cpuid==0) write(stdout, *)'BulkBandSpin, ik, knv3 ', ik, knv3

      k = k3points(:, ik)

      ! calculation bulk hamiltonian
      Hamk_bulk= 0d0
      call ham_bulk_latticegauge(k, Hamk_bulk)

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
            write(outfileindex, '(100000f19.9)')k3len(ik),eigv_mpi(i, ik), spin(:, i, ik)
         enddo
         write(outfileindex, *)' '
      enddo
      close(outfileindex)
   endif

   !> minimum and maximum value of energy bands
   emin= minval(eigv_mpi)-0.5d0
   emax= maxval(eigv_mpi)+0.5d0

   !> write script for gnuplot
   outfileindex= outfileindex+ 1
   if (cpuid==0) then
      open(unit=outfileindex, file='bulkek.gnu')
      write(outfileindex, '(a)') 'set terminal  postscript enhanced color font 24'
      write(outfileindex,'(2a)') '#set palette defined ( 0  "green", ', &
         '5 "yellow", 10 "red" )'
      write(outfileindex, '(a)')"set output 'bulkek.eps'"
      write(outfileindex, '(a)')'set style data linespoints'
      write(outfileindex, '(a)')'unset ztics'
      write(outfileindex, '(a)')'unset key'
      write(outfileindex, '(a)')'set pointsize 0.8'
      write(outfileindex, '(a)')'set view 0,0'
      write(outfileindex, '(a, f10.5, a)')'set xrange [0: ', maxval(k3len), ']'
      if (index(Particle,'phonon')/=0) then
         write(outfileindex, '(a, f10.5, a)')'set yrange [0:', emax, ']'
         write(outfileindex, '(a)')'set ylabel "Frequency (THz)"'
      else
         write(outfileindex, '(a)')'set ylabel "Energy (eV)"'
         write(outfileindex, '(a, f10.5, a, f10.5, a)')'set yrange [', emin, ':', emax, ']'
      endif
      write(outfileindex, 202, advance="no") (k3line_name(i), k3line_stop(i), i=1, nk3lines)
      write(outfileindex, 203)k3line_name(nk3lines+1), k3line_stop(nk3lines+1)

      do i=1, nk3lines-1
         if (index(Particle,'phonon')/=0) then
            write(outfileindex, 204)k3line_stop(i+1), 0.0, k3line_stop(i+1), emax
         else
            write(outfileindex, 204)k3line_stop(i+1), emin, k3line_stop(i+1), emax
         endif
      enddo
      write(outfileindex, '(2a)')"plot 'bulkek.dat' u 1:2 ",  &
         "w lp lw 2 pt 7  ps 0.2, \"
      write(outfileindex, '(2a)')"     'bulkek.dat' u 1:2:($3/6):($4/6) ",  &
         "w vec"
      close(outfileindex)
   endif

202 format('set xtics (',20('"',A3,'" ',F10.5,','))
203 format(A3,'" ',F10.5,')')
204 format('set arrow from ',F10.5,',',F10.5,' to ',F10.5,',',F10.5, ' nohead')

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

      k = k3points(:, ik)

      ! calculation bulk hamiltonian
      Hamk_bulk= 0d0
      call ham_bulk_latticegauge    (k, Hamk_bulk)

      k = k3points(:, ik)
      !k(2)= -k(2)
      Hamk= 0d0
      call ham_bulk_latticegauge    (k, Hamk)

      !> symmetrization
      call mat_mul(Num_wann, mirror_z, hamk, mat1)
      call mat_mul(Num_wann, mat1, mirror_z, mat2)
      hamk= (Hamk_bulk+ mat2)/2.d0

      !> diagonal hamk
      call eigensystem_c('V', 'U', Num_wann, hamk, W)

      mat2= conjg(transpose(hamk))

      !> calculate mirror eigenvalue
      call mat_mul(Num_wann, mat2, mirror_z, mat1)
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
            write(outfileindex, '(1000f19.9)')k3len(ik),eigv_mpi(i, ik)
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
               write(outfileindex, '(1000f19.9)')k3len(ik),eigv_mpi(i, ik)
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
               write(outfileindex, '(1000f19.9)')k3len(ik),eigv_mpi(i, ik)
            endif
         enddo
         write(outfileindex, *)' '
      enddo
      close(outfileindex)
   endif

   !> minimum and maximum value of energy bands
   emin= minval(eigv_mpi)-0.5d0
   emax= maxval(eigv_mpi)+0.5d0

   !> write script for gnuplot
   outfileindex= outfileindex+ 1
   if (cpuid==0) then
      open(unit=outfileindex, file='bulkek.gnu')
      write(outfileindex, '(a)') 'set terminal  postscript enhanced color'
      write(outfileindex,'(2a)') '#set palette defined ( 0  "green", ', &
         '5 "yellow", 10 "red" )'
      write(outfileindex, '(a)')"set output 'bulkek.eps'"
      write(outfileindex, '(a)')'set style data linespoints'
      write(outfileindex, '(a)')'unset ztics'
      write(outfileindex, '(a)')'unset key'
      write(outfileindex, '(a)')'set pointsize 0.8'
      write(outfileindex, '(a)')'set view 0,0'
      write(outfileindex, '(a)')'set xtics font ",24"'
      write(outfileindex, '(a)')'set ytics font ",24"'
      write(outfileindex, '(a)')'set ylabel font ",24"'
      write(outfileindex, '(a, f10.5, a)')'set xrange [0: ', maxval(k3len), ']'
      if (index(Particle,'phonon')/=0) then
         write(outfileindex, '(a, f10.5, a)')'set yrange [0:', emax, ']'
         write(outfileindex, '(a)')'set ylabel "Frequency (THz)"'
      else
         write(outfileindex, '(a, f10.5, a, f10.5, a)')'set yrange [', emin, ':', emax, ']'
         write(outfileindex, '(a)')'set ylabel "Energy (eV)"'
      endif
      write(outfileindex, 202, advance="no") (k3line_name(i), k3line_stop(i), i=1, nk3lines)
      write(outfileindex, 203)k3line_name(nk3lines+1), k3line_stop(nk3lines+1)

      do i=1, nk3lines-1
         if (index(Particle,'phonon')/=0) then
            write(outfileindex, 204)k3line_stop(i+1), 0.0, k3line_stop(i+1), emax
         else
            write(outfileindex, 204)k3line_stop(i+1), emin, k3line_stop(i+1), emax
         endif
      enddo
      write(outfileindex, '(2a)')"plot 'bulkek.dat' u 1:2 ",  &
         "w lp lw 2 pt 7  ps 0.2, 0 w l"
      close(outfileindex)
   endif

202 format('set xtics (',20('"',A3,'" ',F10.5,','))
203 format(A3,'" ',F10.5,')')
204 format('set arrow from ',F10.5,',',F10.5,' to ',F10.5,',',F10.5, ' nohead')

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

      k = k3points(:, ik)

      ! calculation bulk hamiltonian
      Hamk_bulk= 0d0
      call ham_bulk_latticegauge    (k, Hamk_bulk)

      !k = k3points(:, ik)
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
            write(outfileindex, '(1000f19.9)')k3len(ik),eigv_mpi(i, ik)
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
               write(outfileindex, '(1000f19.9)')k3len(ik),eigv_mpi(i, ik)
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
               write(outfileindex, '(1000f19.9)')k3len(ik),eigv_mpi(i, ik)
            endif
         enddo
         write(outfileindex, *)' '
      enddo
      close(outfileindex)
   endif

   !> minimum and maximum value of energy bands
   emin= minval(eigv_mpi)-0.5d0
   emax= maxval(eigv_mpi)+0.5d0

   !> write script for gnuplot
   outfileindex= outfileindex+ 1
   if (cpuid==0) then
      open(unit=outfileindex, file='bulkek.gnu')
      write(outfileindex, '(a)') 'set terminal  postscript enhanced color'
      write(outfileindex,'(2a)') '#set palette defined ( 0  "green", ', &
         '5 "yellow", 10 "red" )'
      write(outfileindex, '(a)')"set output 'bulkek.eps'"
      write(outfileindex, '(a)')'set style data linespoints'
      write(outfileindex, '(a)')'unset ztics'
      write(outfileindex, '(a)')'unset key'
      write(outfileindex, '(a)')'set pointsize 0.8'
      write(outfileindex, '(a)')'set view 0,0'
      write(outfileindex, '(a)')'set xtics font ",24"'
      write(outfileindex, '(a)')'set ytics font ",24"'
      write(outfileindex, '(a)')'set ylabel font ",24"'
      write(outfileindex, '(a, f10.5, a)')'set xrange [0: ', maxval(k3len), ']'
      if (index(Particle,'phonon')/=0) then
         write(outfileindex, '(a, f10.5, a)')'set yrange [0:', emax, ']'
         write(outfileindex, '(a)')'set ylabel "Frequency (THz)"'
      else
         write(outfileindex, '(a)')'set ylabel "Energy (eV)"'
         write(outfileindex, '(a, f10.5, a, f10.5, a)')'set yrange [', emin, ':', emax, ']'
      endif
      write(outfileindex, 202, advance="no") (k3line_name(i), k3line_stop(i), i=1, nk3lines)
      write(outfileindex, 203)k3line_name(nk3lines+1), k3line_stop(nk3lines+1)

      do i=1, nk3lines-1
         if (index(Particle,'phonon')/=0) then
            write(outfileindex, 204)k3line_stop(i+1), 0.0, k3line_stop(i+1), emax
         else
            write(outfileindex, 204)k3line_stop(i+1), emin, k3line_stop(i+1), emax
         endif
      enddo
      write(outfileindex, '(2a)')"plot 'bulkek.dat' u 1:2 ",  &
         "w lp lw 2 pt 7  ps 0.2, 0 w l"
      close(outfileindex)
   endif

202 format('set xtics (',20('"',A3,'" ',F10.5,','))
203 format(A3,'" ',F10.5,')')
204 format('set arrow from ',F10.5,',',F10.5,' to ',F10.5,',',F10.5, ' nohead')

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


subroutine ek_bulk_out
   ! Calculate bulk's energy bands using wannier TB method
   !
   ! Copyright (c) 2010 QuanSheng Wu. All rights reserved.

   use wmpi
   use para

   implicit none

   integer :: ik, i, j
   integer :: knv3
!~    integer :: Nwann
   integer :: ierr
   real(Dp) :: k(3)
   real(Dp), allocatable :: W(:)

   ! Hamiltonian of bulk system
   complex(Dp), allocatable :: Hamk_bulk(:, :)

   ! eigen value of H
   real(dp), allocatable :: eigv(:,:)
   real(dp), allocatable :: eigv_mpi(:,:)
   real(dp), allocatable :: weight(:,:,:)
   real(dp), allocatable :: weight_mpi(:,:,:)
   real(dp), allocatable :: weight_sum(:,:)

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

      k = k3points(:, ik)

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
            call ham_bulk_atomicgauge        (k, Hamk_bulk)
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
      open(unit=outfileindex, file='bulkek.txt')
      write(outfileindex, "('#', a18, 3a19, 1000(9x,' band', i5))")&
         'klen (1/A)', 'k1 (Fractional)', 'k2 (Fractional)', 'k3 (fractional)', &
         (i, i=1, Num_wann) 

      do ik=1, knv3
        !write(outfileindex, '(200f19.9)')k3len(ik),(eigv_mpi(i, ik), i=Numoccupied-1, Numoccupied+2)
         write(outfileindex, '(200f19.9)')k3len(ik),k3points(:, ik), (eigv_mpi(i, ik), i=1, Num_wann)
      enddo
      close(outfileindex)
   endif

   deallocate(W)
   deallocate(Hamk_bulk)
   deallocate( eigv    )
   deallocate( eigv_mpi)
   deallocate( weight    )
   deallocate( weight_mpi)
   deallocate( weight_sum)

   return
end subroutine ek_bulk_out


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
   integer :: ia, ig, ir, ik, io, i, j
   real(dp) :: posij_cart(3), posij_direct(3), posij_direct_unfold(3), kdotr
   real(dp) :: k_cart(3), k_PBZ_direct(3), k_t(3)
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

