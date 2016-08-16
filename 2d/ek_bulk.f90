! calculate bulk's energy band using wannier TB method
! Copyright (c) 2010 QuanSheng Wu. All rights reserved.

  subroutine ek_bulk

     use wmpi
     use para

     implicit none

     integer :: ik, i, j
     integer :: Nwann
     integer :: ierr
     real(dp) :: emin
     real(dp) :: emax
     real(Dp) :: k(2)
     real(Dp) :: W(Num_wann)
     
     ! Hamiltonian of bulk system
     complex(Dp) :: Hamk_bulk(Num_wann,Num_wann) 

     ! eigen value of H
	  real(dp), allocatable :: eigv(:,:)
	  real(dp), allocatable :: eigv_mpi(:,:)
	  real(dp), allocatable :: weight(:,:,:)
	  real(dp), allocatable :: weight_mpi(:,:,:)

     allocate( eigv    (Num_wann, knv2))
     allocate( eigv_mpi(Num_wann, knv2))
     allocate( weight    (Num_wann,Num_wann, knv2))
     allocate( weight_mpi(Num_wann,Num_wann, knv2))
     eigv    = 0d0
     eigv_mpi= 0d0
     weight = 0d0
     weight_mpi = 0d0

     do ik= 1+cpuid, knv2, num_cpu
        if (cpuid==0) write(stdout, *)'BulkBand, ik, knv2 ', ik, knv2

        k = k2_path(ik, :)

        ! calculation bulk hamiltonian
        Hamk_bulk= 0d0
        call ham_bulk_old(k, Hamk_bulk)
       !call ham_bulk    (k, Hamk_bulk)

        !> diagonalization by call zheev in lapack
        W= 0d0
        call eigensystem_c( 'V', 'U', Num_wann ,Hamk_bulk, W)

        eigv(:, ik)= W
        do i=1, Num_wann  !> band 
           do j=1, Num_wann/2  !> projector
              weight(j, i, ik)= sqrt(abs(Hamk_bulk(j, i))**2+ &
                 abs(Hamk_bulk(j+ Num_wann/2, i))**2)
           enddo ! j
        enddo ! i
     enddo ! ik

     call mpi_allreduce(eigv,eigv_mpi,size(eigv),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
     call mpi_allreduce(weight, weight_mpi,size(weight),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)

     weight= weight_mpi/maxval(weight_mpi)*255d0

     if (cpuid==0)then
        open(unit=14, file='bulkek.dat')
   
        do i=1, Num_wann
           do ik=1, knv2
              write(14, '(2f19.9, 1000i5)')k2len(ik),eigv_mpi(i, ik), &
                 int(weight(:, i, ik))
           enddo
           write(14, *)' '
        enddo
        close(14)
     endif

     !> minimum and maximum value of energy bands
    !emin=  minval(eigv_mpi)-0.5d0
    !emax=  maxval(eigv_mpi)+0.5d0
     emin=  -5d0
     emax=   5d0

     !> write script for gnuplot
     if (cpuid==0) then
        open(unit=101, file='bulkek.gnu')
        write(101, '(a)') 'set terminal  postscript enhanced color font ",30"'
        write(101,'(2a)') '#set palette defined ( 0  "green", ', &
           '5 "yellow", 10 "red" )'
        write(101, '(a)')"set output 'bulkek.eps'"
        write(101, '(a)')'set style data linespoints'
        write(101, '(a)')'unset ztics'
        write(101, '(a)')'unset key'
        write(101, '(a)')'set pointsize 0.8'
        write(101, '(a)')'set view 0,0'
        write(101, '(a)')'set xtics font ",24"'
        write(101, '(a)')'set ytics font ",24"'
        write(101, '(a)')'set ylabel font ",24"'
        write(101, '(a)')'set ylabel "Energy (eV)"'
        write(101, '(a)')'set ylabel offset 1.5,0'
        write(101, '(a, f10.5, a)')'set xrange [0: ', maxval(k2len), ']'
        write(101, '(a, f10.5, a, f10.5, a)')'set yrange [', emin, ':', emax, ']'
        write(101, 202, advance="no") (k2line_name(i), k2line_stop(i), i=1, nk2lines)
        write(101, 203)k2line_name(nk2lines+1), k2line_stop(nk2lines+1)

        do i=1, nk2lines-1
           write(101, 204)k2line_stop(i+1), emin, k2line_stop(i+1), emax
        enddo
        write(101, '(2a)')"plot 'bulkek.dat' u 1:2 ",  &
            "w lp lw 2 pt 7  ps 0.2 lc rgb 'black', 0 w l lw 2"
        close(101)
     endif

     202 format('set xtics (',:20('"',A3,'" ',F10.5,','))
     203 format(A3,'" ',F10.5,')')
     204 format('set arrow from ',F10.5,',',F10.5,' to ',F10.5,',',F10.5, ' nohead')
     
   return
   end subroutine ek_bulk


  ! calculate bulk's energy band using wannier TB method
  !> calculate spin direction for each band and each kpoint
  subroutine ek_bulk_spin

     use wmpi
     use para

     implicit none

     integer :: ik, i, j, ib
     integer :: ierr
     integer :: nwann
     real(dp) :: sx, sy, sz
     real(dp) :: emin
     real(dp) :: emax
     real(dp) :: k(2)
     real(dp) :: W(Num_wann)
     
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

     allocate( psi( Num_wann))
     allocate( Hamk_bulk( Num_wann, Num_wann))
     allocate( sigmax( Num_wann, Num_wann))
     allocate( sigmay( Num_wann, Num_wann))
     allocate( sigmaz( Num_wann, Num_wann))
     allocate( eigv    (Num_wann, knv2))
     allocate( eigv_mpi(Num_wann, knv2))
     allocate( spin    (2, Num_wann, knv2))
     allocate( spin_mpi(2, Num_wann, knv2))
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

     do ik= 1+cpuid, knv2, num_cpu
        if (cpuid==0) write(stdout, *)'BulkBandSpin, ik, knv2 ', ik, knv2

        k = k2_path(ik, :)

        ! calculation bulk hamiltonian
        Hamk_bulk= 0d0
        call ham_bulk(k, Hamk_bulk)

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
                 sx= sx+ conjg(psi(i))* sigmax(i, j)* psi(j)
                 sy= sy+ conjg(psi(i))* sigmay(i, j)* psi(j)
              enddo ! j
           enddo ! i
           spin(1, ib, ik)= sx
           spin(2, ib, ik)= sy
        enddo ! ib

     enddo ! ik

     call mpi_allreduce(eigv,eigv_mpi,size(eigv),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
     call mpi_allreduce(spin,spin_mpi,size(spin),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)

     if (cpuid==0)then
        open(unit=14, file='bulkek.dat')
   
        do i=1, Num_wann
           do ik=1, knv2
              write(14, '(1000f19.9)')k2len(ik),eigv_mpi(i, ik), spin(:, i, ik)
           enddo
           write(14, *)' '
        enddo
        close(14)
     endif

     !> minimum and maximum value of energy bands
    !emin= minval(eigv_mpi)-0.5d0
    !emax= maxval(eigv_mpi)+0.5d0
     emin= -5d0
     emax= 5d0 

     !> write script for gnuplot
     if (cpuid==0) then
        open(unit=101, file='bulkek.gnu')
        write(101, '(a)') 'set terminal  postscript enhanced color font 24'
        write(101,'(2a)') '#set palette defined ( 0  "green", ', &
           '5 "yellow", 10 "red" )'
        write(101, '(a)')"set output 'bulkek.eps'"
        write(101, '(a)')'set style data linespoints'
        write(101, '(a)')'unset ztics'
        write(101, '(a)')'unset key'
        write(101, '(a)')'set pointsize 0.8'
        write(101, '(a)')'set view 0,0'
        write(101, '(a)')'set ylabel "Energy (eV)"'
        write(101, '(a, f10.5, a)')'set xrange [0: ', maxval(k2len), ']'
        write(101, '(a, f10.5, a, f10.5, a)')'set yrange [', emin, ':', emax, ']'
        write(101, 202, advance="no") (k2line_name(i), k2line_stop(i), i=1, nk2lines)
        write(101, 203)k2line_name(nk2lines+1), k2line_stop(nk2lines+1)

        do i=1, nk2lines-1
           write(101, 204)k2line_stop(i+1), emin, k2line_stop(i+1), emax
        enddo
        write(101, '(2a)')"plot 'bulkek.dat' u 1:2 ",  &
            "w lp lw 2 pt 7  ps 0.2, \"
        write(101, '(2a)')"     'bulkek.dat' u 1:2:($3/6):($4/6) ",  &
            "w vec"
        close(101)
     endif

     202 format('set xtics (',:20('"',A3,'" ',F10.5,','))
     203 format(A3,'" ',F10.5,')')
     204 format('set arrow from ',F10.5,',',F10.5,' to ',F10.5,',',F10.5, ' nohead')
     
   return
   end subroutine ek_bulk_spin

