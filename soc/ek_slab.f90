  subroutine ek_slab
     !> This subroutine is used for calculating energy 
     !> dispersion with wannier functions for 2D slab system
     !
     ! Copyright (c) 2010 QuanSheng Wu. All rights reserved.
    
     use wmpi
     use para
     implicit none 

     ! loop index
     integer :: i, j, l, lwork, ierr, io

     real(Dp) :: k(2), emin, emax, maxweight
 
     ! time measurement
     real(dp) :: time_start, time_end, time_start0

     ! parameters for zheev
     real(Dp), allocatable ::  rwork(:)
     complex(Dp), allocatable :: work(:)

     ! eigenvalue 
     real(Dp), allocatable :: eigenvalue(:)

     ! energy dispersion
     real(Dp),allocatable :: ekslab(:,:), ekslab_mpi(:,:)

     !> color for plot, surface state weight
     real(dp), allocatable :: surf_l_weight(:, :), surf_l_weight_mpi(:, :)
     real(dp), allocatable :: surf_r_weight(:, :), surf_r_weight_mpi(:, :)

     ! hamiltonian slab
     complex(Dp),allocatable ::CHamk(:,:)

     lwork= 16*Nslab*Num_wann
     ierr = 0


     allocate(eigenvalue(nslab*Num_wann))
     allocate( surf_l_weight (Nslab* Num_wann, knv2))
     allocate( surf_l_weight_mpi (Nslab* Num_wann, knv2))
     allocate( surf_r_weight (Nslab* Num_wann, knv2))
     allocate( surf_r_weight_mpi (Nslab* Num_wann, knv2))
     allocate(ekslab(Nslab*Num_wann,knv2))
     allocate(ekslab_mpi(Nslab*Num_wann,knv2))
     allocate(CHamk(nslab*Num_wann,nslab*Num_wann))
     allocate(work(lwork))
     allocate(rwork(lwork))
 
     surf_l_weight= 0d0
     surf_l_weight_mpi= 0d0
     surf_r_weight= 0d0
     surf_r_weight_mpi= 0d0

     ! sweep k
     ekslab=0.0d0
     ekslab_mpi=0.0d0
     time_start= 0d0
     time_start0= 0d0
     call now(time_start0)
     time_start= time_start0
     time_end  = time_start0
     do i= 1+cpuid, knv2, num_cpu
        if (cpuid==0.and. mod(i/num_cpu, 4)==0) &
           write(stdout, '(a, i9, "  /", i10, a, f10.1, "s", a, f10.1, "s")') &
           ' Slabek: ik', i, knv2, ' time left', &
           (knv2-i)*(time_end- time_start)/num_cpu, &
           ' time elapsed: ', time_end-time_start0 

        call now(time_start)

        k= k2_path(i, :)
        chamk=0.0d0 

        !> no magnetgic field
        if (abs(Bx)<eps9.and. abs(By)<eps9.and. abs(Bz)<eps9)then
           call ham_slab(k,Chamk)
        !> in-plane magnetic field
        elseif (abs(Bx)>eps9 .or. abs(By)>eps9)then
           write(stdout, *)'>> magnetic field is larger than zero'
           call ham_slab_parallel_B(k,Chamk)
        !> vertical magnetic field
        else
           print *, 'Error: we only support in-plane magnetic field at present'
           stop 'please set Bz= 0'
        endif


        eigenvalue=0.0d0

        ! diagonal Chamk
        call eigensystem_c('V', 'U', Num_wann*Nslab, CHamk, eigenvalue)
       
        ekslab(:,i)=eigenvalue

        do j=1, Nslab* Num_wann
           !> left is the bottom surface
           do l= 1, NBottomOrbitals
              io= BottomOrbitals(l) 
              surf_l_weight(j, i)= surf_l_weight(j, i) &
                 + abs(CHamk(io, j))**2  ! first slab
           enddo ! l sweeps the selected orbitals

           !> right is the top surface
           do l= 1, NTopOrbitals
              io= Num_wann*(Nslab-1)+ TopOrbitals(l)   
              surf_r_weight(j, i)= surf_r_weight(j, i) &
                 + abs(CHamk(io, j))**2  ! first slab
           enddo ! l sweeps the selected orbitals

          !do l=1, Num_wann
          !   surf_l_weight(j, i)= surf_l_weight(j, i) &
          !      + abs(CHamk(l, j))**2  ! first slab
          !     !+ abs(CHamk(Num_wann+ l, j))**2 & ! the second slab
          !   surf_r_weight(j, i)= surf_r_weight(j, i) &
          !      + abs(CHamk(Num_wann*Nslab- l+ 1, j))**2 !& ! last slab
          !     !+ abs(CHamk(Num_wann*(Nslab-1)- l, j))**2 ! last second slab
          !enddo ! l
        enddo ! j 
        call now(time_end)
     enddo ! i

#if defined (MPI)
     call mpi_allreduce(ekslab,ekslab_mpi,size(ekslab),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
     call mpi_allreduce(surf_l_weight, surf_l_weight_mpi,size(surf_l_weight),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
     call mpi_allreduce(surf_r_weight, surf_r_weight_mpi,size(surf_r_weight),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
     ekslab_mpi= ekslab
     surf_l_weight_mpi= surf_l_weight
     surf_r_weight_mpi= surf_r_weight
#endif

 
     !> deal with phonon system
     if (index(Particle,'phonon')/=0) then
        do i=1, knv2
           do j=1, Num_wann*Nslab
              ekslab_mpi(j, i)= sqrt(abs(ekslab_mpi(j, i)))*sign(1d0, ekslab_mpi(j, i))
           enddo
        enddo
     endif

     ekslab=ekslab_mpi

     maxweight=maxval(surf_r_weight_mpi+ surf_l_weight_mpi)
     surf_l_weight= surf_l_weight_mpi/ maxweight
     surf_r_weight= surf_r_weight_mpi/ maxweight
     
     outfileindex= outfileindex+ 1
     if(cpuid==0)then
        open(unit=outfileindex, file='slabek.dat')
        write(outfileindex, "('#', a10, a15, 5X, 2a16 )")'# k', ' E', 'BS weight', 'TS weight'
        do j=1, Num_wann*Nslab
           do i=1, knv2
             !write(outfileindex,'(3f15.7, i8)')k2len(i), ekslab(j,i), &
             !   (surf_weight(j, i))
              write(outfileindex,'(2f15.7, 2f16.7)')k2len(i), ekslab(j,i), &
                 (surf_l_weight(j, i)), &
                 (surf_r_weight(j, i))
           enddo
           write(outfileindex , *)''
        enddo
        close(outfileindex)
        write(stdout,*) 'calculate energy band  done'
     endif

     emin= minval(ekslab)-0.5d0
     emax= maxval(ekslab)+0.5d0
     !> write script for gnuplot
     outfileindex= outfileindex+ 1
     if (cpuid==0) then
        open(unit=outfileindex, file='slabek.gnu')
        write(outfileindex, '(a)')"set encoding iso_8859_1"
        write(outfileindex, '(a)')'#set terminal  postscript enhanced color'
        write(outfileindex, '(a)')"#set output 'slabek.eps'"
        write(outfileindex, '(3a)')'#set terminal  pngcairo truecolor enhanced', &
           '  font ",60" size 1920, 1680'
        write(outfileindex, '(3a)')'set terminal  png truecolor enhanced', &
           '  font ",60" size 1920, 1680'
        write(outfileindex, '(a)')"set output 'slabek.png'"
        write(outfileindex,'(2a)') 'set palette defined ( 0  "green", ', &
           '5 "yellow", 10 "red" )'
        write(outfileindex, '(a)')'set style data linespoints'
        write(outfileindex, '(a)')'unset ztics'
        write(outfileindex, '(a)')'unset key'
        write(outfileindex, '(a)')'set pointsize 0.8'
        write(outfileindex, '(a)')'set border lw 3 '
        write(outfileindex, '(a)')'set view 0,0'
        write(outfileindex, '(a)')'#set xtics font ",36"'
        write(outfileindex, '(a)')'#set ytics font ",36"'
        write(outfileindex, '(a)')'#set ylabel font ",36"'
        write(outfileindex, '(a)')'#set xtics offset 0, -1'
        write(outfileindex, '(a)')'set ylabel offset -1, 0 '
        write(outfileindex, '(a, f10.5, a)')'set xrange [0: ', maxval(k2len), ']'
        if (index(Particle,'phonon')/=0) then
           write(outfileindex, '(a, f10.5, a)')'set yrange [0:', emax, ']'
           write(outfileindex, '(a)')'set ylabel "Frequency (THz)"'
        else
           write(outfileindex, '(a)')'set ylabel "Energy (eV)"'
           write(outfileindex, '(a, f10.5, a, f10.5, a)')'set yrange [', emin, ':', emax, ']'
        endif
        write(outfileindex, 202, advance="no") (trim(k2line_name(i)), k2line_stop(i), i=1, nk2lines)
        write(outfileindex, 203)trim(k2line_name(nk2lines+1)), k2line_stop(nk2lines+1)

        do i=1, nk2lines-1
           if (index(Particle,'phonon')/=0) then
              write(outfileindex, 204)k2line_stop(i+1), 0.0, k2line_stop(i+1), emax
           else
              write(outfileindex, 204)k2line_stop(i+1), emin, k2line_stop(i+1), emax
           endif
        enddo
        write(outfileindex, '(a)')'#rgb(r,g,b) = int(r)*65536 + int(g)*256 + int(b)'
        write(outfileindex, '(2a)')"#plot 'slabek.dat' u 1:2:(rgb(255,$3, 3)) ",  &
            "w lp lw 2 pt 7  ps 1 lc rgb variable"
        write(outfileindex, '(2a)')"# (a) "
        write(outfileindex, '(2a)')"# plot the top and bottom surface's weight together"
        write(outfileindex, '(2a)')"#plot 'slabek.dat' u 1:2:($3+$4) ",  &
            "w lp lw 2 pt 7  ps 1 lc palette"
        write(outfileindex, '(2a)')"# (b) "
        write(outfileindex, '(2a)') &
           "# plot top and bottom surface's weight with red and blue respectively"
        write(outfileindex,'(2a)') 'set palette defined ( -1  "blue", ', &
           '0 "grey", 1 "red" )'
        write(outfileindex, '(2a)')"plot 'slabek.dat' u 1:2:($4-$3) ",  &
            "w lp lw 2 pt 7  ps 1 lc palette"

        !write(outfileindex, '(2a)')"splot 'slabek.dat' u 1:2:3 ",  &
        !   "w lp lw 2 pt 13 palette"
        close(outfileindex)
     endif

     202 format('set xtics (',:20('"',A3,'" ',F10.5,','))
     203 format(A3,'" ',F10.5,')')
     204 format('set arrow from ',F10.5,',',F10.5, &
        ' to ',F10.5,',',F10.5, ' nohead')
   
     deallocate(eigenvalue)
     deallocate( surf_l_weight )
     deallocate( surf_l_weight_mpi )
     deallocate( surf_r_weight )
     deallocate( surf_r_weight_mpi )
     deallocate(ekslab)
     deallocate(ekslab_mpi)
     deallocate(CHamk)
     deallocate(work)
     deallocate(rwork)
   
  return
  end subroutine ek_slab

  subroutine ek_slab_kplane
     !> This subroutine is used for calculating energy 
     !> dispersion with wannier functions for 2D slab system
     !
     ! Copyright (c) 2010 QuanSheng Wu. All rights reserved.
    
     use wmpi
     use para
     implicit none 


     ! loop index
     integer :: i, j, l, lwork, ierr, kn12, ik, istart, iend

     ! wave vector 
     real(Dp) :: k(2)

     real(Dp) :: time_start, time_end

     real(Dp), allocatable :: eigenvalue(:)
   
     ! energy dispersion
     real(Dp),allocatable :: ekslab(:,:)
     real(Dp),allocatable :: ekslab_mpi(:,:)

     real(dp), allocatable :: k12(:,:)
     real(dp), allocatable :: k12_shape(:,:)

     !> color for plot, surface state weight
     real(dp), allocatable :: surf_weight(:, :)
     real(dp), allocatable :: surf_weight_mpi(:, :)

     complex(Dp),allocatable ::CHamk(:,:)

     kn12= nk1*nk2
     lwork= 16*Nslab*Num_wann
     ierr = 0

     allocate(eigenvalue(nslab*Num_wann))
     allocate( surf_weight (Nslab* Num_wann, kn12))
     allocate( surf_weight_mpi (Nslab* Num_wann, kn12))
     allocate(ekslab(Nslab*Num_wann,kn12))
     allocate(ekslab_mpi(Nslab*Num_wann,kn12))
     allocate(CHamk(nslab*Num_wann,nslab*Num_wann))
     allocate(k12(2,kn12))
     allocate(k12_shape(2,kn12))
 
     surf_weight= 0d0
     surf_weight_mpi= 0d0

     !> set up k slice
     ik =0
     do i= 1, nk1
        do j= 1, nk2
           ik =ik +1
           k12(:, ik)=K2D_start+ (i-1)*K2D_vec1/dble(nk1-1) &
                      + (j-1)*K2D_vec2/dble(nk2-1)
           k12_shape(:, ik)= k12(1, ik)* Ka2+ k12(2, ik)* Kb2
        enddo
     enddo

     ! sweep k
     ekslab=0.0d0
     ekslab_mpi=0.0d0
     time_start= 0d0
     time_end= 0d0
     do i=1+cpuid, kn12, num_cpu
        if (cpuid==0.and. mod(i/num_cpu, 100)==0) &
           write(stdout, *) 'SlabBand_plane, ik ', i, 'Nk',nk1*nk2, 'time left', &
           (nk1*nk2-i)*(time_end- time_start)/num_cpu, ' s'
        call now(time_start)
        
        k= k12(:, i)
        chamk=0.0d0 

        call ham_slab(k,Chamk)

        eigenvalue=0.0d0

        ! diagonal Chamk
        call eigensystem_c('V', 'U', Num_wann*Nslab, CHamk, eigenvalue)
       
        ekslab(:,i)=eigenvalue

        do j=1, Nslab* Num_wann
           do l=1, Num_wann
              surf_weight(j, i)= surf_weight(j, i) &
                 + abs(CHamk(l, j))**2 & ! first slab
                 + abs(CHamk(Num_wann*Nslab- l+ 1, j))**2 !& ! last slab
                !+ abs(CHamk(Num_wann+ l, j))**2 & ! the second slab
                !+ abs(CHamk(Num_wann*(Nslab-1)- l, j))**2 ! last second slab
           enddo ! l
          !surf_weight(j, i)= (surf_weight(j, i))
        enddo ! j 
        call now(time_end)
     enddo ! i

#if defined (MPI)
     call mpi_allreduce(ekslab,ekslab_mpi,size(ekslab),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
     call mpi_allreduce(surf_weight, surf_weight_mpi,size(surf_weight),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
     ekslab_mpi= ekslab
     surf_weight_mpi= surf_weight
#endif

 
     !> deal with phonon system
     if (index(Particle,'phonon')/=0) then
        do i=1, kn12
           do j=1, Num_wann*Nslab
              ekslab_mpi(j, i)= sqrt(abs(ekslab_mpi(j, i)))*sign(1d0, ekslab_mpi(j, i))
           enddo
        enddo
     endif

     ekslab=ekslab_mpi
     if (maxval(surf_weight_mpi)<0.00001d0)surf_weight_mpi=1d0
     surf_weight= surf_weight_mpi/ maxval(surf_weight_mpi)
     

     istart= Numoccupied*Nslab-1
     iend= Numoccupied*Nslab+2
     outfileindex= outfileindex+ 1
     if (cpuid==0) then
        open(unit=outfileindex, file='slabek_plane.dat')
        write(outfileindex, '(4a16, a)')'# kx', ' ky', ' k1', ' k2', ' (E(ib), dos(ib)), ib=1, NumberofSelectedOrbitals'
        write(outfileindex, '(a, 2i10)')'# Nk1, Nk2=', Nk1, Nk2
        do i=1, kn12
           write(outfileindex,'(2000f16.7)')k12_shape(:,i), k12(:,i), &
              (ekslab(j,i), (255-surf_weight(j, i)*255d0), j=istart, iend)
           if (mod(i, nk1)==0) write (outfileindex, *)' '
        enddo
        close(outfileindex)
        write(stdout,*) 'calculate energy band  done'
     endif

     outfileindex= outfileindex+ 1
     if (cpuid==0) then
        open(unit=outfileindex, file='slabek_plane-matlab.dat')
        write(outfileindex, '(4a16, a)')'% kx', ' ky', ' k1', ' k2', ' (E(ib), dos(ib)), ib=1, NumberofSelectedOrbitals'
        write(outfileindex, '(a, 2i10)')'% Nk1, Nk2=', Nk1, Nk2
        do i=1, kn12
           write(outfileindex,'(2000f16.7)')k12_shape(:,i), k12(:,i), &
              (ekslab(j,i), (255-surf_weight(j, i)*255d0), j=istart, iend)
        enddo
        close(outfileindex)
     endif


     !> write out a script that can be used for gnuplot
     outfileindex= outfileindex+ 1
     if (cpuid==0)then
        open(unit=outfileindex, file='slabek_plane.gnu')
        write(outfileindex, '(a)')"set encoding iso_8859_1"
        write(outfileindex, '(a)')'#set terminal  postscript enhanced color'
        write(outfileindex, '(a)')"#set output 'slabek_plane.eps'"
        write(outfileindex, '(3a)')'set terminal  png truecolor enhanced', &
           ' size 1920, 1680 font ",36"'
        write(outfileindex, '(a)')"set output 'slabek_plane.png'"
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
        write(outfileindex, '(2a)')"splot 'slabek_plane.dat' u 1:2:7 w pm3d, \"
        write(outfileindex, '(2a)')"      'slabek_plane.dat' u 1:2:9 w pm3d"

        close(outfileindex)

     endif ! cpuid

#if defined (MPI)
     call mpi_barrier(mpi_cmw, ierr)
#endif

     deallocate(eigenvalue)
     deallocate( surf_weight )
     deallocate( surf_weight_mpi )
     deallocate(ekslab)
     deallocate(ekslab_mpi)
     deallocate(CHamk)
   
     return
  end subroutine ek_slab_kplane

  subroutine ek_slab_b
  !> This subroutine is used for calculating energy dispersion
  !> with wannier functions in plane magnetic field 
  !> 
    
     use wmpi
     use para
     implicit none 

! loop index
     integer :: i     
     integer :: j
     integer :: l

     integer :: lwork

     real(Dp) :: k(2)

     real(dp) :: emin, emax
 
! parameters for zheev
     integer :: ierr
     real(Dp), allocatable ::  rwork(:)
     complex(Dp), allocatable :: work(:)
      
! eigenvalue 
     real(Dp), allocatable :: eigenvalue(:)
   
! energy dispersion
     real(Dp),allocatable :: ekslab(:,:)
     real(Dp),allocatable :: ekslab_mpi(:,:)

     !> color for plot, surface state weight
     real(dp), allocatable :: surf_weight(:, :)
     real(dp), allocatable :: surf_weight_mpi(:, :)

! hamiltonian slab
     complex(Dp),allocatable ::CHamk(:,:)

     lwork= 16*Nslab*Num_wann
     ierr = 0

     allocate(eigenvalue(nslab*Num_wann))
     allocate(surf_weight (Nslab* Num_wann, knv2))
     allocate(surf_weight_mpi (Nslab* Num_wann, knv2))
     allocate(ekslab(Nslab*Num_wann,knv2))
     allocate(ekslab_mpi(Nslab*Num_wann,knv2))
     allocate(CHamk(nslab*Num_wann,nslab*Num_wann))
     allocate(work(lwork))
     allocate(rwork(lwork))
 
     surf_weight= 0d0
     surf_weight_mpi= 0d0

     ! sweep k
     ekslab=0.0d0
     ekslab_mpi=0.0d0
     do i=1+cpuid,knv2,num_cpu
        if (cpuid==0)write(stdout, *) 'SlabEkB, ik ',  i, knv2
        k= k2_path(i, :)
        chamk=0.0d0 
        eigenvalue=0.0d0
        ! diagonal Chamk
        call eigensystem_c('V', 'U', Num_wann*Nslab, CHamk, eigenvalue)
       
        ekslab(:,i)=eigenvalue

        do j=1, Nslab* Num_wann
           do l=1, Num_wann
              surf_weight(j, i)= surf_weight(j, i) &
                 + abs(CHamk(l, j))**2 & ! first slab
                 + abs(CHamk(Num_wann*Nslab- l+ 1, j))**2 !& ! last slab
                !+ abs(CHamk(Num_wann+ l, j))**2 & ! the second slab
                !+ abs(CHamk(Num_wann*(Nslab-1)- l+ 1, j))**2 ! last second slab
           enddo ! l
           surf_weight(j, i)= (surf_weight(j, i))
        enddo ! j 
        if (cpuid==0) write(stdout, *)'SlabEk,k', i, knv2
     enddo ! i

#if defined (MPI)
     call mpi_allreduce(ekslab,ekslab_mpi,size(ekslab),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
     call mpi_allreduce(surf_weight, surf_weight_mpi,size(surf_weight),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
     ekslab_mpi= ekslab
     surf_weight_mpi= surf_weight
#endif

      !> deal with phonon system
     if (index(Particle,'phonon')/=0) then
        do i=1, knv2
           do j=1, Num_wann*Nslab
              ekslab_mpi(j, i)= sqrt(abs(ekslab_mpi(j, i)))*sign(1d0, ekslab_mpi(j, i))
           enddo
        enddo
     endif


     ekslab=ekslab_mpi
     surf_weight= surf_weight_mpi/ maxval(surf_weight_mpi)
     
     outfileindex= outfileindex+ 1
     if(cpuid==0)then
        open(unit=outfileindex, file='slabek.dat')
        do j=1, Num_wann*Nslab
           do i=1, knv2
             !write(outfileindex,'(3f15.7, i8)')k2len(i), ekslab(j,i), &
             !   (surf_weight(j, i))
              write(outfileindex,'(2f15.7, i8)')k2len(i), ekslab(j,i), &
                 int(255-surf_weight(j, i)*255d0)
           enddo
           write(outfileindex , *)''
        enddo
        close(outfileindex)
        write(stdout,*) 'calculate energy band  done'
     endif

     emin= minval(ekslab)-0.5d0
     emax= maxval(ekslab)+0.5d0
     !> write script for gnuplot
     outfileindex= outfileindex+ 1
     if (cpuid==0) then
        open(unit=outfileindex, file='slabek.gnu')
        write(outfileindex, '(a)')"set encoding iso_8859_1"
        write(outfileindex, '(a)')'#set terminal  postscript enhanced color'
        write(outfileindex, '(a)')"#set output 'slabek.eps'"
        write(outfileindex, '(3a)')'#set terminal  pngcairo truecolor enhanced', &
           '  font ",60" size 1920, 1680'
        write(outfileindex, '(3a)')'set terminal  png truecolor enhanced', &
           '  font ",60" size 1920, 1680'
        write(outfileindex, '(a)')"set output 'slabek.png'"
        write(outfileindex,'(2a)') 'set palette defined ( 0  "green", ', &
           '5 "yellow", 10 "red" )'
        write(outfileindex, '(a)')'set style data linespoints'
        write(outfileindex, '(a)')'unset ztics'
        write(outfileindex, '(a)')'unset key'
        write(outfileindex, '(a)')'set pointsize 0.8'
        write(outfileindex, '(a)')'set border lw 3 '
        write(outfileindex, '(a)')'set view 0,0'
        write(outfileindex, '(a)')'#set xtics font ",36"'
        write(outfileindex, '(a)')'#set ytics font ",36"'
        write(outfileindex, '(a)')'#set ylabel font ",36"'
        write(outfileindex, '(a)')'#set xtics offset 0, -1'
        write(outfileindex, '(a)')'set ylabel offset -1, 0 '
        write(outfileindex, '(a, f10.5, a)')'set xrange [0: ', maxval(k2len), ']'
        if (index(Particle,'phonon')/=0) then
           write(outfileindex, '(a)')'set ylabel "Frequency (THz)"'
           write(outfileindex, '(a, f10.5, a)')'set yrange [0:', emax, ']'
        else
           write(outfileindex, '(a)')'set ylabel "Energy (eV)"'
           write(outfileindex, '(a, f10.5, a, f10.5, a)')'set yrange [', emin, ':', emax, ']'
        endif
        write(outfileindex, 202, advance="no") (trim(k2line_name(i)), k2line_stop(i), i=1, nk2lines)
        write(outfileindex, 203)trim(k2line_name(nk2lines+1)), k2line_stop(nk2lines+1)

        do i=1, nk2lines-1
           if (index(Particle,'phonon')/=0) then
              write(outfileindex, 204)k2line_stop(i+1), 0.0, k2line_stop(i+1), emax
           else
              write(outfileindex, 204)k2line_stop(i+1), emin, k2line_stop(i+1), emax
           endif
        enddo
        write(outfileindex, '(a)')'rgb(r,g,b) = int(r)*65536 + int(g)*256 + int(b)'
        write(outfileindex, '(2a)')"plot 'slabek.dat' u 1:2:(rgb(255,$3, 3)) ",  &
            "w lp lw 2 pt 7  ps 1 lc rgb variable"

        !write(outfileindex, '(2a)')"splot 'slabek.dat' u 1:2:3 ",  &
        !   "w lp lw 2 pt 13 palette"
        close(outfileindex)
     endif

     202 format('set xtics (',:20('"',A3,'" ',F10.5,','))
     203 format(A3,'" ',F10.5,')')
     204 format('set arrow from ',F10.5,',',F10.5,' to ',F10.5,',',F10.5, ' nohead')
 
     deallocate(surf_weight )
     deallocate(surf_weight_mpi )
     deallocate(ekslab)
     deallocate(ekslab_mpi)
     deallocate(CHamk)
     deallocate(work)
     deallocate(rwork)
     deallocate(eigenvalue)
     
  return
  end subroutine ek_slab_b
