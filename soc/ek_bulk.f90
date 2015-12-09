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
	     if (cpuid==0) write(*, *)'ik, knv3 ', ik, knv3
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
        write(101, '(a, f10.5, a)')'set xrange [0: ', maxval(k3len), ']'
        write(101, '(a, f10.5, a, f10.5, a)')'set yrange [', emin, ':', emax, ']'
        write(101, 202, advance="no") (k3line_name(i), k3line_stop(i), i=1, nk3lines)
        write(101, 203)k3line_name(nk3lines+1), k3line_stop(nk3lines+1)

        do i=1, nk3lines-1
           write(101, 204)k3line_stop(i+1), emin, k3line_stop(i+1), emax
        enddo
        write(101, '(2a)')"plot 'bulkek.dat' u 1:2 ",  &
            "w lp lw 2 pt 7  ps 0.2, 0 w l"
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

     use mpi
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

     knv3= nk3_band
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
	     if (cpuid==0) write(*, *)'ik, knv3 ', ik, knv3
	     if (cpuid==0) write(stdout, *)'ik, knv3 ', ik, knv3

        k = k3points(:, ik)

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
                 sz= sz+ conjg(psi(i))* sigmaz(i, j)* psi(j)
              enddo ! j
           enddo ! i
           spin(1, ib, ik)= sx
           spin(2, ib, ik)= sy
           spin(3, ib, ik)= sz
        enddo ! ib

     enddo ! ik

     call mpi_allreduce(eigv,eigv_mpi,size(eigv),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
     call mpi_allreduce(spin,spin_mpi,size(spin),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)

     if (cpuid==0)then
        open(unit=14, file='bulkek.dat')
   
        do i=1, Num_wann
           do ik=1, knv3
              write(14, '(1000f19.9)')k3len(ik),eigv_mpi(i, ik), spin(:, i, ik)
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
        write(101, '(a, f10.5, a)')'set xrange [0: ', maxval(k3len), ']'
        write(101, '(a, f10.5, a, f10.5, a)')'set yrange [', emin, ':', emax, ']'
        write(101, 202, advance="no") (k3line_name(i), k3line_stop(i), i=1, nk3lines)
        write(101, 203)k3line_name(nk3lines+1), k3line_stop(nk3lines+1)

        do i=1, nk3lines-1
           write(101, 204)k3line_stop(i+1), emin, k3line_stop(i+1), emax
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

   subroutine ek_bulk_fortomas

     use mpi
     use para

     implicit none

     integer :: ik1, ik2, ik3
	  integer :: Nk_t
     integer :: ierr
     real(Dp) :: k(3), k1, k2, k3
     real(Dp) :: W(Num_wann)
     
     ! Hamiltonian of bulk system
     complex(Dp) :: Hamk_bulk(Num_wann,Num_wann) 

     if (cpuid==0)then
        open(unit=14, file='bulkek-fortomas.dat')
        open(unit=15, file='bulkek-math.dat')
     endif

     Nk_t= 10
     if (cpuid==0) write (15, '(A)', advance="NO")'{'
     do ik1=0, Nk_t
     if (cpuid==0) write (15, '(A)', advance="NO")'{'
     do ik2=0, Nk_t
     if (cpuid==0) write (15, '(A)', advance="NO")'{'
     do ik3=0, Nk_t
        k1= dble(ik1)/Nk_t
        k2= dble(ik2)/Nk_t
        k3= dble(ik3)/Nk_t

        k(1)= 0.5d0*k2+ 0.5d0*k3
        k(2)= 0.5d0*k3+ 0.5d0*k1
        k(3)= 0.5d0*k1+ 0.5d0*k2

        ! calculation bulk hamiltonian
        Hamk_bulk= 0d0
        call ham_bulk(k, Hamk_bulk)

        !> diagonalization by call zheev in lapack
        W= 0d0
        call eigensystem_c( 'N', 'U', Num_wann ,Hamk_bulk, W)

        if (cpuid==0) then
           write(14, '(3i5, 4f12.8)')ik1, ik2, ik3, W(9:12)
           if (ik3/=Nk_t) then
              write(15, '(a, 4(f12.7, a))', advance="No")"{", W(9), ",", W(10), ",", W(11), ",", W(12), "},"
           else
              write(15, '(a, 4(f12.7, a))', advance="No")"{", W(9), ",", W(10), ",", W(11), ",", W(12), "}"
           endif
        endif
     enddo !ik3
     if (ik2/=Nk_t) then
     if (cpuid==0) write (15, '(A)', advance="NO")'},'
     else
     if (cpuid==0) write (15, '(A)', advance="NO")'}'
     endif
     enddo !ik2
     if (ik1/=Nk_t) then
     if (cpuid==0) write (15, '(A)', advance="NO")'},'
     else
     if (cpuid==0) write (15, '(A)', advance="NO")'}'
     endif
     enddo !ik1
     if (cpuid==0) write (15, '(A)', ADVANCE="NO")'}'
     if (cpuid==0) close(14)
   
   return
   end subroutine ek_bulk_fortomas


