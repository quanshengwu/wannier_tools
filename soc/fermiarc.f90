  subroutine fermiarc
     ! This subroutine calculates surface states using
     ! iterative Green's function method  ---  J.Phys.F.Met.Phys.15(1985)851-858
     !
     ! History:
     !
     !         by Quan Sheng Wu on 4/20/2010
     !
     !            mpi version      4/21/2010


     use wmpi
     use para
     implicit none

     integer :: ierr

     integer :: arclfile, arcrfile, arcbulkfile

     ! general loop index
     integer :: i,j,io

     ! kpoint loop index
     integer :: ikp


     real(dp) :: k(2)

     real(dp) :: omega
     real(dp) :: k1min_shape, k1max_shape, k2min_shape, k2max_shape

     real(dp) :: time_start, time_end

     real(dp), allocatable :: k12(:,:)
     real(dp), allocatable :: k12_shape(:,:)

     real(dp), allocatable :: dos_l(:)
     real(dp), allocatable :: dos_l_mpi(:)
     real(dp), allocatable :: dos_r(:)
     real(dp), allocatable :: dos_r_mpi(:)
     real(dp), allocatable :: dos_l_only(:)
     real(dp), allocatable :: dos_r_only(:)
     real(dp), allocatable :: dos_bulk(:)
     real(dp), allocatable :: dos_bulk_mpi(:)

     complex(dp), allocatable :: ones(:,:)
     complex(dp), allocatable :: GLL(:,:), GRR(:,:), GB(:,:)
     complex(dp), allocatable :: H00(:,:), H01(:,:)


     allocate( k12(2, nk1*nk2))
     allocate( k12_shape(2, nk1*nk2))
     allocate( dos_l(nk1*nk2))
     allocate( dos_l_mpi(nk1*nk2))
     allocate( dos_r(nk1*nk2))
     allocate( dos_l_only(nk1*nk2))
     allocate( dos_r_only(nk1*nk2))
     allocate( dos_r_mpi(nk1*nk2))
     allocate( dos_bulk(nk1*nk2))
     allocate( dos_bulk_mpi(nk1*nk2))
     allocate( GLL(ndim,ndim), GRR(ndim,ndim), GB(ndim, ndim))
     k12=0d0
     k12_shape=0d0
     dos_l=0d0
     dos_l_mpi=0d0
     dos_r=0d0
     dos_r_mpi=0d0
     dos_bulk=0d0
     dos_bulk_mpi=0d0

     ikp=0
     do i= 1, nk1
        do j= 1, nk2
           ikp=ikp+1
           k12(:, ikp)=K2D_start+ (i-1)*K2D_vec1/dble(nk1-1) &
                      + (j-1)*K2D_vec2/dble(nk2-1)
           k12_shape(:, ikp)= k12(1, ikp)* Ka2+ k12(2, ikp)* Kb2
        enddo
     enddo

     k1min_shape= minval(k12_shape(1, :))
     k2min_shape= minval(k12_shape(2, :))
     k1max_shape= maxval(k12_shape(1, :))
     k2max_shape= maxval(k12_shape(2, :))

     allocate(H00(Ndim, Ndim))
     allocate(H01(Ndim, Ndim))
     allocate(ones(Ndim, Ndim))
     GLL= 0d0
     GRR= 0d0
     H00= 0d0
     H01= 0d0
     ones= 0d0

     do i=1,Ndim
        ones(i,i)=1.0d0
     enddo


     omega = E_arc
     eta= eta_arc

     !> deal with phonon system
     !> for phonon system, omega should be changed to omega^2
     if (index(Particle,'phonon')/=0) then
        omega= omega*omega
     endif

     time_start= 0d0
     time_end= 0d0
     do ikp= 1+cpuid, nk1*nk2, num_cpu
        if (cpuid==0.and. mod(ikp/num_cpu, 100)==0) &
           write(stdout, *) 'Arc, ik ', ikp, 'Nk',Nk1*Nk2, 'time left', &
           (nk1*nk2-ikp)/num_cpu*(time_end- time_start), ' s'
        call now(time_start)
        k(1)= k12(1, ikp)
        k(2)= k12(2, ikp)


        if (index(Particle,'phonon')/=0.and.LOTO_correction) then
           call ham_qlayer2qlayer_LOTO(k,H00,H01)
        else
           call ham_qlayer2qlayer(k,H00,H01)
        endif


        !> calculate surface green function
        ! there are two method to calculate surface green's function
        ! the method in 1985 is better, you can find the ref in the
        ! subroutine
        call surfgreen_1985(omega,GLL,GRR,GB,H00,H01,ones)
        ! call surfgreen_1984(omega,GLL,GRR,H00,H01,ones)

        ! calculate spectral function
        do i= 1, NtopOrbitals
           io= TopOrbitals(i)
           dos_l(ikp)=dos_l(ikp)- aimag(GLL(io,io))
        enddo ! i
        do i= 1, NBottomOrbitals
           io= Ndim- Num_wann+ BottomOrbitals(i)
           dos_r(ikp)=dos_r(ikp)- aimag(GRR(io,io))
        enddo ! i
        do i= 1, Ndim
           dos_bulk(ikp)=dos_bulk(ikp)- aimag(GB (i,i))
        enddo ! i
        call now(time_end)

     enddo

     !> we don't have to do allreduce operation
#if defined (MPI)
     call mpi_reduce(dos_l, dos_l_mpi, size(dos_l),mpi_double_precision,&
                     mpi_sum, 0, mpi_comm_world, ierr)
     call mpi_reduce(dos_r, dos_r_mpi, size(dos_r),mpi_double_precision,&
                     mpi_sum, 0, mpi_comm_world, ierr)
     call mpi_reduce(dos_bulk, dos_bulk_mpi, size(dos_bulk), mpi_double_precision,&
                     mpi_sum, 0, mpi_comm_world, ierr)
#else
     dos_l_mpi= dos_l
     dos_r_mpi= dos_r
     dos_bulk_mpi= dos_bulk
#endif

     do ikp=1, Nk1*Nk2
        dos_l_only(ikp)= dos_l_mpi(ikp)- dos_bulk_mpi(ikp)
        if (dos_l_only(ikp)<0) dos_l_only(ikp)=eps9
        dos_r_only(ikp)= dos_r_mpi(ikp)- dos_bulk_mpi(ikp)
        if (dos_r_only(ikp)<0) dos_r_only(ikp)=eps9
     enddo

     outfileindex= outfileindex+ 1
     arclfile= outfileindex
     outfileindex= outfileindex+ 1
     arcrfile= outfileindex
     outfileindex= outfileindex+ 1
     arcbulkfile= outfileindex
     if (cpuid.eq.0)then
        open (unit=arclfile, file='arc.dat_l')
        open (unit=arcrfile, file='arc.dat_r')
        open (unit=arcbulkfile, file='arc.dat_bulk')
        do ikp=1, nk1*nk2
           write(arclfile, '(30f16.8)')k12_shape(:, ikp), log(dos_l_mpi(ikp)), log(dos_l_only(ikp))
           if (mod(ikp, nk2)==0) write(arclfile, *)' '
           write(arcrfile, '(30f16.8)')k12_shape(:, ikp), log(dos_r_mpi(ikp)), log(dos_r_only(ikp))
           if (mod(ikp, nk2)==0) write(arcrfile, *)' '
           write(arcbulkfile, '(30f16.8)')k12_shape(:, ikp), log(abs(dos_bulk_mpi(ikp)))
           if (mod(ikp, nk2)==0) write(arcbulkfile, *)' '
        enddo
        close(arclfile)
        close(arcrfile)
        close(arcbulkfile)

        write(stdout,*)'ndim',ndim
        write(stdout,*)'Nk1,Nk2,eta', Nk1, Nk2, eta
        write(stdout,*)'calculate density of state successfully'
     endif

     !> write script for gnuplot for bulk green's function
     outfileindex= outfileindex+ 1
     if (cpuid==0) then
        open(unit=outfileindex, file='arc_bulk.gnu')
        write(outfileindex, '(a)')"set encoding iso_8859_1"
        write(outfileindex, '(a)')'#set terminal  postscript enhanced color'
        write(outfileindex, '(a)')"#set output 'arc_l.eps'"
        write(outfileindex, '(3a)')'#set terminal pngcairo truecolor enhanced', &
           '  font ",50" size 1920, 1680'
        write(outfileindex, '(3a)')'set terminal png truecolor enhanced', &
           ' font ",50" size 1920, 1680'
        write(outfileindex, '(a)')"set output 'arc_bulk.png'"
        write(outfileindex,'(2a)') 'set palette defined ( -10 "#194eff", ', &
           '0 "white", 10 "red" )'
        write(outfileindex, '(a)')'#set palette rgbformulae 33,13,10'
        write(outfileindex, '(a)')'unset ztics'
        write(outfileindex, '(a)')'unset key'
        write(outfileindex, '(a)')'set pm3d'
        write(outfileindex, '(a)')'set border lw 6'
        write(outfileindex, '(a)')'set size ratio -1'
        write(outfileindex, '(a)')'set view map'
        write(outfileindex, '(a)')'set xtics'
        write(outfileindex, '(a)')'set ytics'
        write(outfileindex, '(a)')'set xlabel "K_1 (1/{\305})"'
        write(outfileindex, '(a)')'set ylabel "K_2 (1/{\305})"'
        write(outfileindex, '(a)')'set ylabel offset 1, 0'
        write(outfileindex, '(a)')'set colorbox'
        write(outfileindex, '(a, f8.5, a, f8.5, a)')'set xrange [', k1min_shape, ':', k1max_shape, ']'
        write(outfileindex, '(a, f8.5, a, f8.5, a)')'set yrange [', k2min_shape, ':', k2max_shape, ']'
        write(outfileindex, '(a)')'set pm3d interpolate 2,2'
        write(outfileindex, '(2a)')"splot 'arc.dat_bulk' u 1:2:3 w pm3d"
        close(outfileindex)
     endif

     !> write script for gnuplot
     outfileindex= outfileindex+ 1
     if (cpuid==0) then
        open(unit=outfileindex, file='arc_l_only.gnu')
        write(outfileindex, '(a)')"set encoding iso_8859_1"
        write(outfileindex, '(a)')'#set terminal  postscript enhanced color'
        write(outfileindex, '(a)')"#set output 'arc_l.eps'"
        write(outfileindex, '(3a)')'#set terminal pngcairo truecolor enhanced', &
           '  font ",50" size 1920, 1680'
        write(outfileindex, '(3a)')'set terminal png truecolor enhanced', &
           ' font ",50" size 1920, 1680'
        write(outfileindex, '(a)')"set output 'arc_l_only.png'"
        write(outfileindex,'(2a)') 'set palette defined (0  "white", ', &
           '6 "red", 20 "black" )'
        write(outfileindex, '(a)')'#set palette rgbformulae 33,13,10'
        write(outfileindex, '(a)')'unset ztics'
        write(outfileindex, '(a)')'unset key'
        write(outfileindex, '(a)')'set pm3d'
        write(outfileindex, '(a)')'set border lw 6'
        write(outfileindex, '(a)')'set size ratio -1'
        write(outfileindex, '(a)')'set view map'
        write(outfileindex, '(a)')'set xtics'
        write(outfileindex, '(a)')'set ytics'
        write(outfileindex, '(a)')'set xlabel "K_1 (1/{\305})"'
        write(outfileindex, '(a)')'set ylabel "K_2 (1/{\305})"'
        write(outfileindex, '(a)')'set ylabel offset 1, 0'
        write(outfileindex, '(a)')'set colorbox'
        write(outfileindex, '(a)')'unset cbtics'
        write(outfileindex, '(a, f8.5, a, f8.5, a)')'set xrange [', k1min_shape, ':', k1max_shape, ']'
        write(outfileindex, '(a, f8.5, a, f8.5, a)')'set yrange [', k2min_shape, ':', k2max_shape, ']'
        write(outfileindex, '(a)')'set pm3d interpolate 2,2'
        write(outfileindex, '(2a)')"splot 'arc.dat_l' u 1:2:(exp($4)) w pm3d"

        close(outfileindex)
     endif



     !> write script for gnuplot
     outfileindex= outfileindex+ 1
     if (cpuid==0) then
        open(unit=outfileindex, file='arc_l.gnu')
        write(outfileindex, '(a)')"set encoding iso_8859_1"
        write(outfileindex, '(a)')'#set terminal  postscript enhanced color'
        write(outfileindex, '(a)')"#set output 'arc_l.eps'"
        write(outfileindex, '(3a)')'#set terminal pngcairo truecolor enhanced', &
           '  font ",50" size 1920, 1680'
        write(outfileindex, '(3a)')'set terminal png truecolor enhanced', &
           ' font ",50" size 1920, 1680'
        write(outfileindex, '(a)')"set output 'arc_l.png'"
        write(outfileindex,'(2a)') 'set palette defined ( -10 "#194eff", ', &
           '0 "white", 10 "red" )'
        write(outfileindex, '(a)')'#set palette rgbformulae 33,13,10'
        write(outfileindex, '(a)')'unset ztics'
        write(outfileindex, '(a)')'unset key'
        write(outfileindex, '(a)')'set pm3d'
        write(outfileindex, '(a)')'set border lw 6'
        write(outfileindex, '(a)')'set size ratio -1'
        write(outfileindex, '(a)')'set view map'
        write(outfileindex, '(a)')'set xtics'
        write(outfileindex, '(a)')'set ytics'
        write(outfileindex, '(a)')'set xlabel "K_1 (1/{\305})"'
        write(outfileindex, '(a)')'set ylabel "K_2 (1/{\305})"'
        write(outfileindex, '(a)')'set ylabel offset 1, 0'
        write(outfileindex, '(a)')'set colorbox'
        write(outfileindex, '(a, f8.5, a, f8.5, a)')'set xrange [', k1min_shape, ':', k1max_shape, ']'
        write(outfileindex, '(a, f8.5, a, f8.5, a)')'set yrange [', k2min_shape, ':', k2max_shape, ']'
        write(outfileindex, '(a)')'set pm3d interpolate 2,2'
        write(outfileindex, '(2a)')"splot 'arc.dat_l' u 1:2:3 w pm3d"

        close(outfileindex)
     endif

     !> write script for gnuplot
     outfileindex= outfileindex+ 1
     if (cpuid==0) then
        open(unit=outfileindex, file='arc_r_only.gnu')
        write(outfileindex, '(a)')"set encoding iso_8859_1"
        write(outfileindex, '(a)')'#set terminal  postscript enhanced color'
        write(outfileindex, '(a)')"#set output 'arc_r.eps'"
        write(outfileindex, '(3a)')'#set terminal pngcairo truecolor enhanced', &
           '  font ",50" size 1920, 1680'
        write(outfileindex, '(3a)')'set terminal png truecolor enhanced', &
           ' font ",50" size 1920, 1680'
        write(outfileindex, '(a)')"set output 'arc_r_only.png'"
        write(outfileindex,'(2a)') 'set palette defined (0  "white", ', &
           '6 "red", 20 "black" )'
        write(outfileindex, '(a)')'#set palette rgbformulae 33,13,10'
        write(outfileindex, '(a)')'unset ztics'
        write(outfileindex, '(a)')'unset key'
        write(outfileindex, '(a)')'set pm3d'
        write(outfileindex, '(a)')'set border lw 6'
        write(outfileindex, '(a)')'set size ratio -1'
        write(outfileindex, '(a)')'set view map'
        write(outfileindex, '(a)')'set xtics'
        write(outfileindex, '(a)')'set ytics'
        write(outfileindex, '(a)')'set xlabel "K_1 (1/{\305})"'
        write(outfileindex, '(a)')'set ylabel "K_2 (1/{\305})"'
        write(outfileindex, '(a)')'set ylabel offset 1, 0'
        write(outfileindex, '(a)')'set colorbox'
        write(outfileindex, '(a)')'unset cbtics'
        write(outfileindex, '(a, f8.5, a, f8.5, a)')'set xrange [', k1min_shape, ':', k1max_shape, ']'
        write(outfileindex, '(a, f8.5, a, f8.5, a)')'set yrange [', k2min_shape, ':', k2max_shape, ']'
        write(outfileindex, '(a)')'set pm3d interpolate 2,2'
        write(outfileindex, '(2a)')"splot 'arc.dat_r' u 1:2:(exp($4)) w pm3d"

        close(outfileindex)
     endif


     !> write script for gnuplot
     outfileindex= outfileindex+ 1
     if (cpuid==0) then
        open(unit=outfileindex, file='arc_r.gnu')
        write(outfileindex, '(a)')"set encoding iso_8859_1"
        write(outfileindex, '(a)')'#set terminal  postscript enhanced color'
        write(outfileindex, '(a)')"#set output 'arc_r.eps'"
        write(outfileindex, '(3a)')'#set terminal pngcairo truecolor enhanced', &
           '  font ",50" size 1920, 1680'
        write(outfileindex, '(3a)')'set terminal png truecolor enhanced', &
           ' font ",50" size 1920, 1680'
        write(outfileindex, '(a)')"set output 'arc_r.png'"
        write(outfileindex,'(2a)') 'set palette defined ( -10 "#194eff", ', &
           '0 "white", 10 "red" )'
        write(outfileindex, '(a)')'#set palette rgbformulae 33,13,10'
        write(outfileindex, '(a)')'unset ztics'
        write(outfileindex, '(a)')'unset key'
        write(outfileindex, '(a)')'set pm3d'
        write(outfileindex, '(a)')'set border lw 6'
        write(outfileindex, '(a)')'set size ratio -1'
        write(outfileindex, '(a)')'set view map'
        write(outfileindex, '(a)')'set xtics'
        write(outfileindex, '(a)')'set ytics'
        write(outfileindex, '(a)')'set xlabel "K_1 (1/{\305})"'
        write(outfileindex, '(a)')'set ylabel "K_2 (1/{\305})"'
        write(outfileindex, '(a)')'set ylabel offset 1, 0'
        write(outfileindex, '(a)')'set colorbox'
        write(outfileindex, '(a, f8.5, a, f8.5, a)')'set xrange [', k1min_shape, ':', k1max_shape, ']'
        write(outfileindex, '(a, f8.5, a, f8.5, a)')'set yrange [', k2min_shape, ':', k2max_shape, ']'
        write(outfileindex, '(a)')'set pm3d interpolate 2,2'
        write(outfileindex, '(2a)')"splot 'arc.dat_r' u 1:2:3 w pm3d"

        close(outfileindex)
     endif

#if defined (MPI)
     call mpi_barrier(mpi_cmw, ierr)
#endif
     deallocate( k12)
     deallocate( k12_shape)
     deallocate( dos_l)
     deallocate( dos_l_mpi)
     deallocate( dos_r)
     deallocate( dos_l_only)
     deallocate( dos_r_only)
     deallocate( dos_r_mpi)
     deallocate( dos_bulk)
     deallocate( dos_bulk_mpi)
     deallocate( GLL, GRR, GB)
     deallocate(H00)
     deallocate(H01)
     deallocate(ones)


     return
  end subroutine fermiarc


  subroutine fermiarc_jdos
     !> Calculation joint density of state
     !> refs :http://science.sciencemag.org/content/sci/suppl/2016/03/09/351.6278.1184.DC1/Inoue.SM.pdf
     !> in this version, we also calculate the spin density jdos


     use wmpi
     use para
     implicit none

     integer :: ierr

     integer :: arclfile, arcrfile, arcbulkfile
     integer :: arcljfile, arcrjfile
     integer :: arcljsfile, arcrjsfile
     integer :: spindosrfile, spindoslfile

     ! general loop index
     integer :: i,j,io

     integer :: Nwann
     ! kpoint loop index
     integer :: ikp, ik1, ik2, iq, Nk1_half, Nk2_half
     integer :: imin1, imax1, imin2, imax2, iq1, iq2

     real(dp) :: k(2), dis

     real(dp) :: omega
     real(dp) :: k1min_shape, k1max_shape, k2min_shape, k2max_shape

     real(dp) :: dos_l_max, dos_r_max
     real(dp) :: time_start, time_end

     real(dp) :: sx_bulk, sy_bulk, sz_bulk

     integer , allocatable :: ik12(:,:)
     real(dp), allocatable :: k12(:,:)
     real(dp), allocatable :: k12_shape(:,:)

     real(dp), allocatable :: dos_l(:,:)
     real(dp), allocatable :: dos_l_mpi(:,:)
     real(dp), allocatable :: dos_r(:,:)
     real(dp), allocatable :: dos_r_mpi(:,:)
     real(dp), allocatable :: dos_bulk(:,:)
     real(dp), allocatable :: dos_bulk_mpi(:,:)
     real(dp), allocatable :: jdos_l(:)
     real(dp), allocatable :: jdos_l_mpi(:)
     real(dp), allocatable :: jdos_r(:)
     real(dp), allocatable :: jdos_r_mpi(:)
     real(dp), allocatable :: jsdos_l(:)
     real(dp), allocatable :: jsdos_l_mpi(:)
     real(dp), allocatable :: jsdos_r(:)
     real(dp), allocatable :: jsdos_r_mpi(:)
     real(dp), allocatable :: dos_l_only(:, :)
     real(dp), allocatable :: dos_r_only(:, :)
     real(dp), allocatable :: sx_l(:, :)
     real(dp), allocatable :: sy_l(:, :)
     real(dp), allocatable :: sz_l(:, :)
     real(dp), allocatable :: sx_l_mpi(:, :)
     real(dp), allocatable :: sy_l_mpi(:, :)
     real(dp), allocatable :: sz_l_mpi(:, :)
     real(dp), allocatable :: sx_r(:, :)
     real(dp), allocatable :: sy_r(:, :)
     real(dp), allocatable :: sz_r(:, :)
     real(dp), allocatable :: sx_r_mpi(:, :)
     real(dp), allocatable :: sy_r_mpi(:, :)
     real(dp), allocatable :: sz_r_mpi(:, :)

     ! spin operator matrix sigma_x,sigma_y in sigma_z representation
     complex(Dp),allocatable :: sigma_x(:,:)
     complex(Dp),allocatable :: sigma_y(:,:)
     complex(Dp),allocatable :: sigma_z(:,:)
     complex(Dp),allocatable :: ctemp(:,:)
     complex(dp), allocatable :: ones(:,:)
     complex(dp), allocatable :: GLL(:,:), GRR(:,:), GB(:,:)
     complex(dp), allocatable :: H00(:,:), H01(:,:)

     Nk1_half= (Nk1-1)/2
     Nk2_half= (Nk2-1)/2
     Nwann= Num_wann/2

     allocate( ik12(2, nk1*nk2))
     allocate( k12(2, nk1*nk2))
     allocate( k12_shape(2, nk1*nk2))
     allocate( dos_l(nk1,nk2))
     allocate( dos_l_mpi(nk1,nk2))
     allocate( dos_r(nk1,nk2))
     allocate( dos_r_mpi(nk1,nk2))
     allocate( dos_r_only(nk1,nk2))
     allocate( dos_l_only(nk1,nk2))
     allocate( dos_bulk(nk1,nk2))
     allocate( dos_bulk_mpi(nk1,nk2))
     allocate( jdos_l(nk1*nk2))
     allocate( jdos_l_mpi(nk1*nk2))
     allocate( jdos_r(nk1*nk2))
     allocate( jdos_r_mpi(nk1*nk2))
     allocate( jsdos_l(nk1*nk2))
     allocate( jsdos_l_mpi(nk1*nk2))
     allocate( jsdos_r(nk1*nk2))
     allocate( jsdos_r_mpi(nk1*nk2))
     allocate( sx_l(nk1,nk2))
     allocate( sx_l_mpi(nk1,nk2))
     allocate( sy_l(nk1,nk2))
     allocate( sy_l_mpi(nk1,nk2))
     allocate( sz_l(nk1,nk2))
     allocate( sz_l_mpi(nk1,nk2))
     allocate( sx_r(nk1,nk2))
     allocate( sx_r_mpi(nk1,nk2))
     allocate( sy_r(nk1,nk2))
     allocate( sy_r_mpi(nk1,nk2))
     allocate( sz_r(nk1,nk2))
     allocate( sz_r_mpi(nk1,nk2))
     allocate( GLL(ndim,ndim), GRR(ndim,ndim), GB(ndim, ndim))
     allocate(sigma_x(ndim,ndim))
     allocate(sigma_y(ndim,ndim))
     allocate(sigma_z(ndim,ndim))
     allocate(ctemp(ndim,ndim))

     sigma_x   =0.0d0
     sigma_y   =0.0d0
     sigma_z   =0.0d0
     ik12=0
     k12=0d0
     k12_shape=0d0
     dos_l=0d0
     dos_l_mpi=0d0
     dos_r=0d0
     dos_r_mpi=0d0
     jdos_l=0d0
     jdos_l_mpi=1d-12
     jdos_r=0d0
     jdos_r_mpi=1d-12
     jsdos_l=0d0
     jsdos_l_mpi=1d-12
     jsdos_r=0d0
     jsdos_r_mpi=1d-12
     sx_l=0d0
     sy_l=0d0
     sz_l=0d0
     sx_l_mpi=0d0
     sy_l_mpi=0d0
     sz_l_mpi=0d0
     sx_r=0d0
     sy_r=0d0
     sz_r=0d0
     sx_r_mpi=0d0
     sy_r_mpi=0d0
     sz_r_mpi=0d0

     ikp=0
     do i= 1, nk1
        do j= 1, nk2
           ikp=ikp+1
           ik12(1, ikp)= i
           ik12(2, ikp)= j
           k12(:, ikp)=K2D_start+ (i-1)*K2D_vec1/dble(nk1-1) &
                      + (j-1)*K2D_vec2/dble(nk2-1)
           k12_shape(:, ikp)= k12(1, ikp)* Ka2+ k12(2, ikp)* Kb2
        enddo
     enddo

     k1min_shape= minval(k12_shape(1, :))
     k2min_shape= minval(k12_shape(2, :))
     k1max_shape= maxval(k12_shape(1, :))
     k2max_shape= maxval(k12_shape(2, :))

     allocate(H00(Ndim, Ndim))
     allocate(H01(Ndim, Ndim))
     allocate(ones(Ndim, Ndim))
     GLL= 0d0
     GRR= 0d0
     H00= 0d0
     H01= 0d0
     ones= 0d0

     do i=1,Ndim
        ones(i,i)=1.0d0
     enddo

     Nwann= Num_wann/2

     !> spin operator matrix
     do i=1, Np
        do j=1, Nwann
           sigma_x(Num_wann*(i-1)+j, Num_wann*(i-1)+Nwann+j)=1.0d0
           sigma_x(Num_wann*(i-1)+j+Nwann, Num_wann*(i-1)+j)=1.0d0
           sigma_y(Num_wann*(i-1)+j, Num_wann*(i-1)+Nwann+j)=-zi
           sigma_y(Num_wann*(i-1)+j+Nwann, Num_wann*(i-1)+j)=zi
           sigma_z(Num_wann*(i-1)+j, Num_wann*(i-1)+j)= 1d0
           sigma_z(Num_wann*(i-1)+j+Nwann, Num_wann*(i-1)+j+Nwann)=-1d0
        enddo
     enddo

     omega = E_arc
     eta= eta_arc

     !> deal with phonon system
     !> for phonon system, omega should be changed to omega^2
     if (index(Particle,'phonon')/=0) then
        omega= omega*omega
     endif


     time_start= 0d0
     time_end= 0d0
     do ikp= 1+cpuid, nk1*nk2, num_cpu
        if (cpuid==0.and. mod(ikp/num_cpu, 100)==0) &
           write(stdout, *) 'Arc, ik ', ikp, 'Nk',Nk1*Nk2, 'time left', &
           (nk1*nk2-ikp)*(time_end- time_start)/num_cpu, ' s'
        call now(time_start)
        k(1)= k12(1, ikp)
        k(2)= k12(2, ikp)

        if (index(Particle,'phonon')/=0.and.LOTO_correction) then
           call ham_qlayer2qlayer_LOTO(k,H00,H01)
        else
           call ham_qlayer2qlayer(k,H00,H01)
        endif


        !> calculate surface green function
        ! there are two method to calculate surface green's function
        ! the method in 1985 is better, you can find the ref in the
        ! subroutine
        call surfgreen_1985(omega,GLL,GRR,GB,H00,H01,ones)
        ! call surfgreen_1984(omega,GLL,GRR,H00,H01,ones)

        ik1= ik12(1, ikp)
        ik2= ik12(2, ikp)

        ! calculate spectral function
        do i= 1, NtopOrbitals
           io= TopOrbitals(i)
           dos_l(ik1, ik2)=dos_l(ik1, ik2)- aimag(GLL(io,io))
        enddo ! i
        do i= 1, NBottomOrbitals
           io= Ndim- Num_wann+ BottomOrbitals(i)
           dos_r(ik1, ik2)=dos_r(ik1, ik2)- aimag(GRR(io,io))
        enddo ! i

        do i= 1, Ndim
        !  dos_l(ik1, ik2)=dos_l(ik1, ik2)- aimag(GLL(i,i))
        enddo ! i
        do i= 1, Ndim
        !  dos_r(ik1, ik2)=dos_r(ik1, ik2)- aimag(GRR(i,i))
        enddo ! i
        DO i= 1, Ndim
           dos_bulk(ik1, ik2)=dos_bulk(ik1, ik2)- AIMAG(GB(i,i))
        ENDDO ! i

        !>> calculate spin-resolved bulk spectrum
        sx_bulk= 0d0
        call mat_mul(ndim,GB ,sigma_x,ctemp)
        do i= 1, Ndim
           sx_bulk=sx_bulk- aimag(ctemp(i,i))
        enddo ! i

        sy_bulk= 0d0
        call mat_mul(ndim,GB ,sigma_y,ctemp)
        do i= 1, Ndim
           sy_bulk=sy_bulk- aimag(ctemp(i,i))
        enddo ! i

        sz_bulk= 0d0
        call mat_mul(ndim,GB ,sigma_z,ctemp)
        do i= 1, Ndim
           sz_bulk=sz_bulk- aimag(ctemp(i,i))
        enddo ! i

        !>> calculate spin-resolved surface spectrum

        !ctemp=matmul(surfgreen,sigma_x)
        call mat_mul(ndim,GLL,sigma_x,ctemp)
        do i= 1, NtopOrbitals
           io= TopOrbitals(i)
           sx_l(ik1, ik2)=sx_l(ik1, ik2)- aimag(ctemp(io,io))
        enddo ! i
       !sx_l(ik1, ik2)= sx_l(ik1, ik2)- sx_bulk
       !if (sx_l(ik1, ik2)<0) sx_l(ik1, ik2)= eps9

        !ctemp=matmul(surfgreen,sigma_y)
        call mat_mul(ndim,GLL,sigma_y,ctemp)
        do i= 1, NtopOrbitals
           io= TopOrbitals(i)
           sy_l(ik1, ik2)=sy_l(ik1, ik2)- aimag(ctemp(io,io))
        enddo ! i
       !sy_l(ik1, ik2)= sy_l(ik1, ik2)- sy_bulk
       !if (sy_l(ik1, ik2)<0) sy_l(ik1, ik2)= eps9


        !ctemp=matmul(surfgreen,sigma_z)
        call mat_mul(ndim,GLL,sigma_z,ctemp)
        do i= 1, NtopOrbitals
           io= TopOrbitals(i)
           sz_l(ik1, ik2)=sz_l(ik1, ik2)- aimag(ctemp(io,io))
        enddo ! i
       !sz_l(ik1, ik2)= sz_l(ik1, ik2)- sz_bulk
       !if (sz_l(ik1, ik2)<0) sz_l(ik1, ik2)= eps9


        !ctemp=matmul(surfgreen,sigma_x)
        call mat_mul(ndim,GRR,sigma_x,ctemp)
        do i= 1, NBottomOrbitals
           io= Ndim- Num_wann+ BottomOrbitals(i)
           sx_r(ik1, ik2)=sx_r(ik1, ik2)- aimag(ctemp(io,io))
        enddo ! i
       !sx_r(ik1, ik2)= sx_r(ik1, ik2)- sz_bulk
       !if (sx_r(ik1, ik2)<0) sx_r(ik1, ik2)= eps9


        !ctemp=matmul(surfgreen,sigma_y)
        call mat_mul(ndim,GRR,sigma_y,ctemp)
        do i= 1, NBottomOrbitals
           io= Ndim- Num_wann+ BottomOrbitals(i)
           sy_r(ik1, ik2)=sy_r(ik1, ik2)- aimag(ctemp(io,io))
        enddo ! i
       !sy_r(ik1, ik2)= sy_r(ik1, ik2)- sz_bulk
       !if (sy_r(ik1, ik2)<0) sy_r(ik1, ik2)= eps9


        !ctemp=matmul(surfgreen,sigma_z)
        call mat_mul(ndim,GRR,sigma_z,ctemp)
        do i= 1, NBottomOrbitals
           io= Ndim- Num_wann+ BottomOrbitals(i)
           sz_r(ik1, ik2)=sz_r(ik1, ik2)- aimag(ctemp(io,io))
        enddo ! i
       !sz_r(ik1, ik2)= sz_r(ik1, ik2)- sz_bulk
       !if (sz_r(ik1, ik2)<0) sz_r(ik1, ik2)= eps9

        call now(time_end)

     enddo

#if defined (MPI)
     !> we don't have to do allreduce operation
     call mpi_allreduce(dos_l, dos_l_mpi, size(dos_l),mpi_double_precision,&
                     mpi_sum, mpi_comm_world, ierr)
     call mpi_allreduce(dos_r, dos_r_mpi, size(dos_r),mpi_double_precision,&
                     mpi_sum, mpi_comm_world, ierr)
     call mpi_allreduce(dos_bulk, dos_bulk_mpi, size(dos_bulk),mpi_double_precision,&
                     mpi_sum, mpi_comm_world, ierr)
     call mpi_allreduce(sx_l, sx_l_mpi, size(sx_l),mpi_double_precision,&
                     mpi_sum, mpi_comm_world, ierr)
     call mpi_allreduce(sy_l, sy_l_mpi, size(sy_l),mpi_double_precision,&
                     mpi_sum, mpi_comm_world, ierr)
     call mpi_allreduce(sz_l, sz_l_mpi, size(sz_l),mpi_double_precision,&
                     mpi_sum, mpi_comm_world, ierr)
     call mpi_allreduce(sx_r, sx_r_mpi, size(sx_r),mpi_double_precision,&
                     mpi_sum, mpi_comm_world, ierr)
     call mpi_allreduce(sy_r, sy_r_mpi, size(sy_r),mpi_double_precision,&
                     mpi_sum, mpi_comm_world, ierr)
     call mpi_allreduce(sz_r, sz_r_mpi, size(sz_r),mpi_double_precision,&
                     mpi_sum, mpi_comm_world, ierr)
#else
     dos_l_mpi= dos_l
     dos_r_mpi= dos_r
     dos_bulk_mpi= dos_bulk
     sx_l_mpi= sx_l
     sy_l_mpi= sy_l
     sz_l_mpi= sz_l
     sx_r_mpi= sx_r
     sy_r_mpi= sy_r
     sz_r_mpi= sz_r
#endif

!    do ikp=1, Nk1*Nk2
!       ik1= ik12(1, ikp)
!       ik2= ik12(2, ikp)

!       dos_l_only(ik1, ik2)= dos_l_mpi(ik1, ik2)- dos_bulk_mpi(ik1, ik2)
!       if (dos_l_only(ik1, ik2)<0) then
!          dos_l_only(ik1, ik2)=eps9
!          sx_l_mpi(ik1, ik2)=eps9
!          sy_l_mpi(ik1, ik2)=eps9
!          sz_l_mpi(ik1, ik2)=eps9
!       endif
!       dos_r_only(ik1, ik2)= dos_r_mpi(ik1, ik2)- dos_bulk_mpi(ik1, ik2)
!       if (dos_r_only(ik1, ik2)<0) then
!          dos_r_only(ik1, ik2)=eps9
!          sx_r_mpi(ik1, ik2)=eps9
!          sy_r_mpi(ik1, ik2)=eps9
!          sz_r_mpi(ik1, ik2)=eps9
!       endif
!    enddo

     dos_l_max= maxval(dos_l_mpi)
     dos_r_max= maxval(dos_r_mpi)

     do ikp=1, Nk1*Nk2
        ik1= ik12(1, ikp)
        ik2= ik12(2, ikp)

        dos_l_only(ik1, ik2)= dos_l_mpi(ik1, ik2)
        if (dos_l_only(ik1, ik2)<dos_l_max/10d0) then
           dos_l_only(ik1, ik2)=eps9
           sx_l_mpi(ik1, ik2)=eps9
           sy_l_mpi(ik1, ik2)=eps9
           sz_l_mpi(ik1, ik2)=eps9
        endif
        dos_r_only(ik1, ik2)= dos_r_mpi(ik1, ik2)
        if (dos_r_only(ik1, ik2)<dos_r_max/10d0) then
           dos_r_only(ik1, ik2)=eps9
           sx_r_mpi(ik1, ik2)=eps9
           sy_r_mpi(ik1, ik2)=eps9
           sz_r_mpi(ik1, ik2)=eps9
        endif
     enddo

     outfileindex= outfileindex+ 1
     arclfile= outfileindex
     outfileindex= outfileindex+ 1
     arcrfile= outfileindex
     outfileindex= outfileindex+ 1
     spindoslfile= outfileindex
     outfileindex= outfileindex+ 1
     spindosrfile= outfileindex
     if (cpuid.eq.0)then
        open (unit=arclfile, file='arc.dat_l')
        open (unit=arcrfile, file='arc.dat_r')
        open(spindoslfile,file='spindos_l.dat')
        open(spindoslfile,file='spindos_r.dat')
        do ikp=1, nk1*nk2
           write(arclfile, '(30f16.8)')k12_shape(:, ikp), log(dos_l_mpi(ik12(1, ikp), ik12(2, ikp))) &
                                     , log(dos_l_only(ik12(1, ikp), ik12(2, ikp)))
           if (mod(ikp, nk2)==0) write(arclfile, *)' '
           write(arcrfile, '(30f16.8)')k12_shape(:, ikp), log(dos_r_mpi(ik12(1, ikp), ik12(2, ikp))) &
                                     , log(dos_r_only(ik12(1, ikp), ik12(2, ikp)))
           if (mod(ikp, nk2)==0) write(arcrfile, *)' '
           write(spindoslfile, '(30f16.8)')k12_shape(:, ikp), &
                                     log(dos_l_only(ik12(1, ikp), ik12(2, ikp))), &
                                     (sx_l_mpi(ik12(1, ikp), ik12(2, ikp)))/dos_l_mpi(ik12(1, ikp), ik12(2, ikp)), &
                                     (sy_l_mpi(ik12(1, ikp), ik12(2, ikp)))/dos_l_mpi(ik12(1, ikp), ik12(2, ikp)), &
                                     (sz_l_mpi(ik12(1, ikp), ik12(2, ikp)))/dos_l_mpi(ik12(1, ikp), ik12(2, ikp))
           write(spindoslfile, '(30f16.8)')k12_shape(:, ikp), &
                                     log(dos_r_only(ik12(1, ikp), ik12(2, ikp))), &
                                     (sx_r_mpi(ik12(1, ikp), ik12(2, ikp)))/dos_r_mpi(ik12(1, ikp), ik12(2, ikp)), &
                                     (sy_r_mpi(ik12(1, ikp), ik12(2, ikp)))/dos_r_mpi(ik12(1, ikp), ik12(2, ikp)), &
                                     (sz_r_mpi(ik12(1, ikp), ik12(2, ikp)))/dos_r_mpi(ik12(1, ikp), ik12(2, ikp))
           if (mod(ikp, nk2)==0) write(spindoslfile, *)' '
           if (mod(ikp, nk2)==0) write(spindoslfile, *)' '
        enddo
        close(arclfile)
        close(arcrfile)
        close(spindoslfile)
        close(spindoslfile)
        write(stdout,*)'ndim',ndim
        write(stdout,*) 'Nk1,Nk2,eta',Nk1, Nk2, eta
        write(stdout,*)'calculate density of state successfully'
     endif

     !> calculate jdos
     do iq= 1+ cpuid, Nk1*Nk2, num_cpu
        iq1= ik12(1, iq)- Nk1_half
        iq2= ik12(2, iq)- Nk2_half
        if (cpuid==0.and. mod(iq/num_cpu, 100)==0) &
           write(stdout, *) 'JDOS, iq ', iq, 'Nq',Nk1*Nk2, 'time left', &
           (nk1*nk2-iq)*(time_end- time_start)/num_cpu, ' s'
        call now(time_start)
        imin1= max(-Nk1_half-iq1, -Nk1_half)+ Nk1_half+ 1
        imax1= min(Nk1_half-iq1, Nk1_half)+ Nk1_half+ 1
        imin2= max(-Nk2_half-iq2, -Nk2_half)+ Nk2_half+ 1
        imax2= min(Nk2_half-iq2, Nk2_half)+ Nk2_half+ 1
        do ik2= imin2, imax2
           do ik1= imin1, imax1
              jdos_l(iq)= jdos_l(iq)+ dos_l_only(ik1, ik2)* dos_l_only(ik1+iq1, ik2+iq2)
              jdos_r(iq)= jdos_r(iq)+ dos_r_only(ik1, ik2)* dos_r_only(ik1+iq1, ik2+iq2)
              jsdos_l(iq)= jsdos_l(iq)+ dos_l_only(ik1, ik2)* dos_l_only(ik1+iq1, ik2+iq2) &
                                      + sx_l_mpi(ik1, ik2)* sx_l_mpi(ik1+iq1, ik2+iq2) &
                                      + sy_l_mpi(ik1, ik2)* sy_l_mpi(ik1+iq1, ik2+iq2) &
                                      + sz_l_mpi(ik1, ik2)* sz_l_mpi(ik1+iq1, ik2+iq2)
              jsdos_r(iq)= jsdos_r(iq)+ dos_r_only(ik1, ik2)* dos_r_only(ik1+iq1, ik2+iq2) &
                                      + sx_r_mpi(ik1, ik2)* sx_r_mpi(ik1+iq1, ik2+iq2) &
                                      + sy_r_mpi(ik1, ik2)* sy_r_mpi(ik1+iq1, ik2+iq2) &
                                      + sz_r_mpi(ik1, ik2)* sz_r_mpi(ik1+iq1, ik2+iq2)
           enddo !ik1
        enddo !ik2
        call now(time_end)
     enddo !iq

     jdos_l_mpi=1d-12
     jdos_r_mpi=1d-12
     jsdos_l_mpi=1d-12
     jsdos_r_mpi=1d-12

#if defined (MPI)
     call mpi_reduce(jdos_l, jdos_l_mpi, size(jdos_l),mpi_double_precision,&
                     mpi_sum, 0, mpi_comm_world, ierr)
     call mpi_reduce(jdos_r, jdos_r_mpi, size(jdos_r),mpi_double_precision,&
                     mpi_sum, 0, mpi_comm_world, ierr)
     call mpi_reduce(jsdos_l, jsdos_l_mpi, size(jsdos_l),mpi_double_precision,&
                     mpi_sum, 0, mpi_comm_world, ierr)
     call mpi_reduce(jsdos_r, jsdos_r_mpi, size(jsdos_r),mpi_double_precision,&
                     mpi_sum, 0, mpi_comm_world, ierr)
#else
     jdos_l_mpi= jdos_l
     jdos_r_mpi= jdos_r
     jsdos_l_mpi= jsdos_l
     jsdos_r_mpi= jsdos_r
#endif

     outfileindex= outfileindex+ 1
     arcljfile= outfileindex
     outfileindex= outfileindex+ 1
     arcrjfile= outfileindex
     outfileindex= outfileindex+ 3
     arcljsfile= outfileindex
     outfileindex= outfileindex+ 3
     arcrjsfile= outfileindex
     if (cpuid.eq.0)then
        write(stdout,*)'The calculation of joint density of state was done.'
        write(stdout,*)'Now it is ready to write out.'
        open (unit=arcljfile, file='arc.jdat_l')
        do ikp=1, nk1*nk2
           write(arcljfile, '(30f16.8)')k12_shape(:, ikp),  log(abs(jdos_l_mpi(ikp)))
           if (mod(ikp, nk2)==0) write(arcljfile, *)' '
        enddo
        close(arcljfile)

        open (unit=arcrjfile, file='arc.jdat_r')
        do ikp=1, nk1*nk2
           write(arcrjfile, '(30f16.8)')k12_shape(:, ikp),  log(abs(jdos_r_mpi(ikp)))
           if (mod(ikp, nk2)==0) write(arcrjfile, *)' '
        enddo
        close(arcrjfile)

        open (unit=arcljsfile, file='arc.jsdat_l')
        do ikp=1, nk1*nk2
           write(arcljsfile, '(30f16.8)')k12_shape(:, ikp), log(abs(jsdos_l_mpi(ikp)))
           if (mod(ikp, nk2)==0) write(arcljsfile, *)' '
        enddo
        close(arcljsfile)

        open (unit=arcrjsfile, file='arc.jsdat_r')
        do ikp=1, nk1*nk2
           write(arcrjsfile, '(30f16.8)')k12_shape(:, ikp), log(abs(jsdos_r_mpi(ikp)))
           if (mod(ikp, nk2)==0) write(arcrjsfile, *)' '
        enddo
        close(arcrjsfile)
        write(stdout,*)'calculate joint density of state successfully'
     endif ! cpuid==0

     !> write script for gnuplot
     outfileindex= outfileindex+ 1
     if (cpuid==0) then
        open(unit=outfileindex, file='arc_l_jsdos.gnu')
        write(outfileindex, '(a)')"set encoding iso_8859_1"
        write(outfileindex, '(a)')'#set terminal  postscript enhanced color'
        write(outfileindex, '(a)')"#set output 'arc_l_jsdos.eps'"
        write(outfileindex, '(3a)')'#set terminal pngcairo truecolor enhanced', &
           '  font ",50" size 1920, 1680'
        write(outfileindex, '(3a)')'set terminal png truecolor enhanced', &
           ' font ",50" size 1920, 1680'
        write(outfileindex, '(a)')"set output 'arc_l_jsdos.png'"
        write(outfileindex,'(2a)') 'set palette defined (0  "white", ', &
           '6 "red", 20 "black" )'
        write(outfileindex, '(a)')'#set palette rgbformulae 33,13,10'
        write(outfileindex, '(a)')'unset ztics'
        write(outfileindex, '(a)')'unset key'
        write(outfileindex, '(a)')'set pm3d'
        write(outfileindex, '(a)')'set border lw 6'
        write(outfileindex, '(a)')'set size ratio -1'
        write(outfileindex, '(a)')'set view map'
        write(outfileindex, '(a)')'set xtics'
        write(outfileindex, '(a)')'set ytics'
        write(outfileindex, '(a)')'set xlabel "K_1 (1/{\305})"'
        write(outfileindex, '(a)')'set ylabel "K_2 (1/{\305})"'
        write(outfileindex, '(a)')'set ylabel offset 1, 0'
        write(outfileindex, '(a)')'set colorbox'
        write(outfileindex, '(a)')'unset cbtics'
        write(outfileindex, '(a, f8.5, a, f8.5, a)')'set xrange [', k1min_shape, ':', k1max_shape, ']'
        write(outfileindex, '(a, f8.5, a, f8.5, a)')'set yrange [', k2min_shape, ':', k2max_shape, ']'
        write(outfileindex, '(a)')'set pm3d interpolate 2,2'
        write(outfileindex, '(2a)')"splot 'arc.jsdat_l' u 1:2:(exp($3)) w pm3d"
        close(outfileindex)
     endif


     !> write script for gnuplot
     outfileindex= outfileindex+ 1
     if (cpuid==0) then
        open(unit=outfileindex, file='arc_l_jdos.gnu')
        write(outfileindex, '(a)')"set encoding iso_8859_1"
        write(outfileindex, '(a)')'#set terminal  postscript enhanced color'
        write(outfileindex, '(a)')"#set output 'arc_l_jdos.eps'"
        write(outfileindex, '(3a)')'#set terminal pngcairo truecolor enhanced', &
           '  font ",50" size 1920, 1680'
        write(outfileindex, '(3a)')'set terminal png truecolor enhanced', &
           ' font ",50" size 1920, 1680'
        write(outfileindex, '(a)')"set output 'arc_l_jdos.png'"
        write(outfileindex,'(2a)') 'set palette defined (0  "white", ', &
           '6 "red", 20 "black" )'
        write(outfileindex, '(a)')'#set palette rgbformulae 33,13,10'
        write(outfileindex, '(a)')'unset ztics'
        write(outfileindex, '(a)')'unset key'
        write(outfileindex, '(a)')'set pm3d'
        write(outfileindex, '(a)')'set border lw 6'
        write(outfileindex, '(a)')'set size ratio -1'
        write(outfileindex, '(a)')'set view map'
        write(outfileindex, '(a)')'set xtics'
        write(outfileindex, '(a)')'set ytics'
        write(outfileindex, '(a)')'set xlabel "K_1 (1/{\305})"'
        write(outfileindex, '(a)')'set ylabel "K_2 (1/{\305})"'
        write(outfileindex, '(a)')'set ylabel offset 1, 0'
        write(outfileindex, '(a)')'set colorbox'
        write(outfileindex, '(a)')'unset cbtics'
        write(outfileindex, '(a, f8.5, a, f8.5, a)')'set xrange [', k1min_shape, ':', k1max_shape, ']'
        write(outfileindex, '(a, f8.5, a, f8.5, a)')'set yrange [', k2min_shape, ':', k2max_shape, ']'
        write(outfileindex, '(a)')'set pm3d interpolate 2,2'
        write(outfileindex, '(2a)')"splot 'arc.jdat_l' u 1:2:(exp($3)) w pm3d"
        close(outfileindex)
     endif

     !> write script for gnuplot
     outfileindex= outfileindex+ 1
     if (cpuid==0) then
        open(unit=outfileindex, file='arc_l_only.gnu')
        write(outfileindex, '(a)')"set encoding iso_8859_1"
        write(outfileindex, '(a)')'#set terminal  postscript enhanced color'
        write(outfileindex, '(a)')"#set output 'arc_l_only.eps'"
        write(outfileindex, '(3a)')'#set terminal pngcairo truecolor enhanced', &
           '  font ",50" size 1920, 1680'
        write(outfileindex, '(3a)')'set terminal png truecolor enhanced', &
           ' font ",50" size 1920, 1680'
        write(outfileindex, '(a)')"set output 'arc_l_only.png'"
        write(outfileindex,'(2a)') 'set palette defined (0  "white", ', &
           '6 "red", 20 "black" )'
        write(outfileindex, '(a)')'#set palette rgbformulae 33,13,10'
        write(outfileindex, '(a)')'unset ztics'
        write(outfileindex, '(a)')'unset key'
        write(outfileindex, '(a)')'set pm3d'
        write(outfileindex, '(a)')'set border lw 6'
        write(outfileindex, '(a)')'set size ratio -1'
        write(outfileindex, '(a)')'set view map'
        write(outfileindex, '(a)')'set xtics'
        write(outfileindex, '(a)')'set ytics'
        write(outfileindex, '(a)')'set xlabel "K_1 (1/{\305})"'
        write(outfileindex, '(a)')'set ylabel "K_2 (1/{\305})"'
        write(outfileindex, '(a)')'set ylabel offset 1, 0'
        write(outfileindex, '(a)')'set colorbox'
        write(outfileindex, '(a)')'unset cbtics'
        write(outfileindex, '(a, f8.5, a, f8.5, a)')'set xrange [', k1min_shape, ':', k1max_shape, ']'
        write(outfileindex, '(a, f8.5, a, f8.5, a)')'set yrange [', k2min_shape, ':', k2max_shape, ']'
        write(outfileindex, '(a)')'set pm3d interpolate 2,2'
        write(outfileindex, '(2a)')"splot 'arc.dat_l' u 1:2:(exp($4)) w pm3d"
        close(outfileindex)
     endif

     !> write script for gnuplot
     outfileindex= outfileindex+ 1
     if (cpuid==0) then
        open(unit=outfileindex, file='arc_r_jsdos.gnu')
        write(outfileindex, '(a)')"set encoding iso_8859_1"
        write(outfileindex, '(a)')'#set terminal  postscript enhanced color'
        write(outfileindex, '(a)')"#set output 'arc_r_jsdos.eps'"
        write(outfileindex, '(3a)')'#set terminal pngcairo truecolor enhanced', &
           '  font ",50" size 1920, 1680'
        write(outfileindex, '(3a)')'set terminal png truecolor enhanced', &
           ' font ",50" size 1920, 1680'
        write(outfileindex, '(a)')"set output 'arc_r_jsdos.png'"
        write(outfileindex,'(2a)') 'set palette defined (0  "white", ', &
           '6 "red", 20 "black" )'
        write(outfileindex, '(a)')'#set palette rgbformulae 33,13,10'
        write(outfileindex, '(a)')'unset ztics'
        write(outfileindex, '(a)')'unset key'
        write(outfileindex, '(a)')'set pm3d'
        write(outfileindex, '(a)')'set border lw 6'
        write(outfileindex, '(a)')'set size ratio -1'
        write(outfileindex, '(a)')'set view map'
        write(outfileindex, '(a)')'set xtics'
        write(outfileindex, '(a)')'set ytics'
        write(outfileindex, '(a)')'set xlabel "K_1 (1/{\305})"'
        write(outfileindex, '(a)')'set ylabel "K_2 (1/{\305})"'
        write(outfileindex, '(a)')'set ylabel offset 1, 0'
        write(outfileindex, '(a)')'set colorbox'
        write(outfileindex, '(a)')'unset cbtics'
        write(outfileindex, '(a, f8.5, a, f8.5, a)')'set xrange [', k1min_shape, ':', k1max_shape, ']'
        write(outfileindex, '(a, f8.5, a, f8.5, a)')'set yrange [', k2min_shape, ':', k2max_shape, ']'
        write(outfileindex, '(a)')'set pm3d interpolate 2,2'
        write(outfileindex, '(2a)')"splot 'arc.jsdat_r' u 1:2:(exp($3)) w pm3d"
        close(outfileindex)
     endif


     !> write script for gnuplot
     outfileindex= outfileindex+ 1
     if (cpuid==0) then
        open(unit=outfileindex, file='arc_r_jdos.gnu')
        write(outfileindex, '(a)')"set encoding iso_8859_1"
        write(outfileindex, '(a)')'#set terminal  postscript enhanced color'
        write(outfileindex, '(a)')"#set output 'arc_r_jdos.eps'"
        write(outfileindex, '(3a)')'#set terminal pngcairo truecolor enhanced', &
           '  font ",50" size 1920, 1680'
        write(outfileindex, '(3a)')'set terminal png truecolor enhanced', &
           ' font ",50" size 1920, 1680'
        write(outfileindex, '(a)')"set output 'arc_r_jdos.png'"
        write(outfileindex,'(2a)') 'set palette defined (0  "white", ', &
           '6 "red", 20 "black" )'
        write(outfileindex, '(a)')'#set palette rgbformulae 33,13,10'
        write(outfileindex, '(a)')'unset ztics'
        write(outfileindex, '(a)')'unset key'
        write(outfileindex, '(a)')'set pm3d'
        write(outfileindex, '(a)')'set border lw 6'
        write(outfileindex, '(a)')'set size ratio -1'
        write(outfileindex, '(a)')'set view map'
        write(outfileindex, '(a)')'set xtics'
        write(outfileindex, '(a)')'set ytics'
        write(outfileindex, '(a)')'set xlabel "K_1 (1/{\305})"'
        write(outfileindex, '(a)')'set ylabel "K_2 (1/{\305})"'
        write(outfileindex, '(a)')'set ylabel offset 1, 0'
        write(outfileindex, '(a)')'set colorbox'
        write(outfileindex, '(a)')'unset cbtics'
        write(outfileindex, '(a, f8.5, a, f8.5, a)')'set xrange [', k1min_shape, ':', k1max_shape, ']'
        write(outfileindex, '(a, f8.5, a, f8.5, a)')'set yrange [', k2min_shape, ':', k2max_shape, ']'
        write(outfileindex, '(a)')'set pm3d interpolate 2,2'
        write(outfileindex, '(2a)')"splot 'arc.jdat_r' u 1:2:(exp($3)) w pm3d"
        close(outfileindex)
     endif


     !> write script for gnuplot
     outfileindex= outfileindex+ 1
     if (cpuid==0) then
        open(unit=outfileindex, file='arc_r_only.gnu')
        write(outfileindex, '(a)')"set encoding iso_8859_1"
        write(outfileindex, '(a)')'#set terminal  postscript enhanced color'
        write(outfileindex, '(a)')"#set output 'arc_r_only.eps'"
        write(outfileindex, '(3a)')'#set terminal pngcairo truecolor enhanced', &
           '  font ",50" size 1920, 1680'
        write(outfileindex, '(3a)')'set terminal png truecolor enhanced', &
           ' font ",50" size 1920, 1680'
        write(outfileindex, '(a)')"set output 'arc_r_only.png'"
        write(outfileindex,'(2a)') 'set palette defined (0  "white", ', &
           '6 "red", 20 "black" )'
        write(outfileindex, '(a)')'#set palette rgbformulae 33,13,10'
        write(outfileindex, '(a)')'unset ztics'
        write(outfileindex, '(a)')'unset key'
        write(outfileindex, '(a)')'set pm3d'
        write(outfileindex, '(a)')'set border lw 6'
        write(outfileindex, '(a)')'set size ratio -1'
        write(outfileindex, '(a)')'set view map'
        write(outfileindex, '(a)')'set xtics'
        write(outfileindex, '(a)')'set ytics'
        write(outfileindex, '(a)')'set xlabel "K_1 (1/{\305})"'
        write(outfileindex, '(a)')'set ylabel "K_2 (1/{\305})"'
        write(outfileindex, '(a)')'set ylabel offset 1, 0'
        write(outfileindex, '(a)')'set colorbox'
        write(outfileindex, '(a)')'unset cbtics'
        write(outfileindex, '(a, f8.5, a, f8.5, a)')'set xrange [', k1min_shape, ':', k1max_shape, ']'
        write(outfileindex, '(a, f8.5, a, f8.5, a)')'set yrange [', k2min_shape, ':', k2max_shape, ']'
        write(outfileindex, '(a)')'set pm3d interpolate 2,2'
        write(outfileindex, '(2a)')"splot 'arc.dat_r' u 1:2:(exp($4)) w pm3d"

        close(outfileindex)
     endif

#if defined (MPI)
     call mpi_barrier(mpi_cmw, ierr)
#endif

     deallocate( ik12)
     deallocate( k12)
     deallocate( k12_shape)
     deallocate( dos_l)
     deallocate( dos_l_mpi)
     deallocate( dos_r)
     deallocate( dos_r_mpi)
     deallocate( dos_r_only)
     deallocate( dos_l_only)
     deallocate( dos_bulk)
     deallocate( dos_bulk_mpi)
     deallocate( jdos_l)
     deallocate( jdos_l_mpi)
     deallocate( jdos_r)
     deallocate( jdos_r_mpi)
     deallocate( jsdos_l)
     deallocate( jsdos_l_mpi)
     deallocate( jsdos_r)
     deallocate( jsdos_r_mpi)
     deallocate( sx_l)
     deallocate( sx_l_mpi)
     deallocate( sy_l)
     deallocate( sy_l_mpi)
     deallocate( sz_l)
     deallocate( sz_l_mpi)
     deallocate( sx_r)
     deallocate( sx_r_mpi)
     deallocate( sy_r)
     deallocate( sy_r_mpi)
     deallocate( sz_r)
     deallocate( sz_r_mpi)
     deallocate( GLL, GRR, GB)
     deallocate(sigma_x)
     deallocate(sigma_y)
     deallocate(sigma_z)
     deallocate(ctemp)
     deallocate(H00)
     deallocate(H01)
     deallocate(ones)



  return
  end subroutine fermiarc_jdos

