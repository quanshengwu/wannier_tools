! calculate bulk's energy band using wannier TB method
  subroutine ek_bulk

     use mpi
     use para

     implicit none

     integer :: ik, i, j
	  integer :: knv3
     integer :: ierr
     real(dp) :: emin
     real(dp) :: emax
     real(Dp) :: k(3)
     real(Dp) :: W(Num_wann)
     
     ! Hamiltonian of bulk system
     complex(Dp) :: Hamk_bulk(Num_wann,Num_wann) 

     ! eigen value of H
	  real(dp), allocatable :: eigv(:,:)
	  real(dp), allocatable :: eigv_mpi(:,:)

     knv3= nk3_band
     allocate( eigv    (Num_wann, knv3))
     allocate( eigv_mpi(Num_wann, knv3))
	  eigv    = 0d0
	  eigv_mpi= 0d0

     do ik= 1+cpuid, knv3, num_cpu
	     if (cpuid==0) write(stdout, *)'ik, knv3 ', ik, knv3

        k = k3points(:, ik)

        ! calculation bulk hamiltonian
        Hamk_bulk= 0d0
        call ham_bulk(k, Hamk_bulk)

        !> diagonalization by call zheev in lapack
        W= 0d0
        call eigensystem_c( 'N', 'U', Num_wann ,Hamk_bulk, W)

        eigv(:, ik)= W

     enddo

     call mpi_allreduce(eigv,eigv_mpi,size(eigv),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)

     if (cpuid==0)then
        open(unit=14, file='bulkek.dat')
   
        do i=1, Num_wann
           do ik=1, knv3
              write(14, '(1000f19.9)')k3len(ik),eigv_mpi(i, ik)
           enddo
           write(14, *)' '
        enddo
        close(14)
     endif

     !> minimum and maximum value of energy bands
     emin= minval(eigv_mpi)-0.5d0
     emax= maxval(eigv_mpi)+0.5d0

     !> write script for gnuplot
     if (cpuid==0) then
        open(unit=101, file='bulkek.gnu')
        write(101, '(a)') 'set terminal  postscript enhanced color'
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
        write(101, '(a, f8.5, a)')'set xrange [0: ', maxval(k3len), ']'
        write(101, '(a, f8.5, a, f8.5, a)')'set yrange [', emin, ':', emax, ']'
        write(101, 202, advance="no") (k3line_name(i), k3line_stop(i), i=1, nk3lines)
        write(101, 203)k3line_name(nk3lines+1), k3line_stop(nk3lines+1)

        do i=1, nk3lines-1
           write(101, 204)k3line_stop(i+1), emin, k3line_stop(i+1), emax
        enddo
        write(101, '(2a)')"plot 'bulkek.dat' u 1:2 ",  &
            "w lp lw 2 pt 7  ps 0.2"
        close(101)
     endif

     202 format('set xtics (',:20('"',A3,'" ',F8.5,','))
     203 format(A3,'" ',F8.5,')')
     204 format('set arrow from ',F8.5,',',F10.5,' to ',F8.5,',',F10.5, ' nohead')
     
   return
   end subroutine ek_bulk
