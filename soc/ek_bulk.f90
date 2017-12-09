! calculate bulk's energy band using wannier TB method
! Copyright (c) 2010 QuanSheng Wu. All rights reserved.

  subroutine ek_bulk

     use wmpi
     use para

     implicit none

     integer :: ik, i, j
     integer :: knv3
     integer :: Nwann
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
	  real(dp), allocatable :: weight(:,:,:)
	  real(dp), allocatable :: weight_mpi(:,:,:)

     knv3= nk3_band
     allocate( eigv    (Num_wann, knv3))
     allocate( eigv_mpi(Num_wann, knv3))
     allocate( weight    (Num_wann,Num_wann, knv3))
     allocate( weight_mpi(Num_wann,Num_wann, knv3))
     eigv    = 0d0
     eigv_mpi= 0d0
     weight = 0d0
     weight_mpi = 0d0

     do ik= 1+cpuid, knv3, num_cpu
        if (cpuid==0) write(stdout, *)'BulkBand, ik, knv3 ', ik, knv3

        k = k3points(:, ik)

        ! calculation bulk hamiltonian
        Hamk_bulk= 0d0
       !call ham_bulk_old(k, Hamk_bulk)
        call ham_bulk    (k, Hamk_bulk)

        !> diagonalization by call zheev in lapack
        W= 0d0
        call eigensystem_c( 'V', 'U', Num_wann ,Hamk_bulk, W)

        if (sum(abs(k))<1e-5) call  orbital_momenta(k, Hamk_bulk)

        eigv(:, ik)= W
        do i=1, Num_wann  !> band 
           if (SOC==0) then
              do j=1, Num_wann  !> projector
                 weight(j, i, ik)= (abs(Hamk_bulk(j, i))**2)
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

     weight= weight_mpi/maxval(weight_mpi)*255d0

     outfileindex= outfileindex+ 1
     if (cpuid==0)then
        open(unit=outfileindex, file='bulkek.dat')
   
        do i=1, Num_wann
           do ik=1, knv3
              write(outfileindex, '(2f19.9, 1000i5)')k3len(ik),eigv_mpi(i, ik), &
                 int(weight(:, i, ik))
           enddo
           write(outfileindex, *)' '
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
        write(outfileindex, '(a)')'set ylabel "Energy (eV)"'
        write(outfileindex, '(a)')'set ylabel offset 1.5,0'
        write(outfileindex, '(a, f10.5, a)')'set xrange [0: ', maxval(k3len), ']'
        write(outfileindex, '(a, f10.5, a, f10.5, a)')'set yrange [', emin, ':', emax, ']'
        write(outfileindex, 202, advance="no") (k3line_name(i), k3line_stop(i), i=1, nk3lines)
        write(outfileindex, 203)k3line_name(nk3lines+1), k3line_stop(nk3lines+1)

        do i=1, nk3lines-1
           write(outfileindex, 204)k3line_stop(i+1), emin, k3line_stop(i+1), emax
        enddo
        write(outfileindex, '(2a)')"plot 'bulkek.dat' u 1:2 ",  &
            "w lp lw 2 pt 7  ps 0.2 lc rgb 'black', 0 w l lw 2"
        close(outfileindex)
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
        if (cpuid==0) write(stdout, *)'BulkBandSpin, ik, knv3 ', ik, knv3

        k = k3points(:, ik)

        ! calculation bulk hamiltonian
        Hamk_bulk= 0d0
        call ham_bulk_old(k, Hamk_bulk)

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
              write(outfileindex, '(1000f19.9)')k3len(ik),eigv_mpi(i, ik), spin(:, i, ik)
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
        write(outfileindex, '(a)')'set ylabel "Energy (eV)"'
        write(outfileindex, '(a, f10.5, a)')'set xrange [0: ', maxval(k3len), ']'
        write(outfileindex, '(a, f10.5, a, f10.5, a)')'set yrange [', emin, ':', emax, ']'
        write(outfileindex, 202, advance="no") (k3line_name(i), k3line_stop(i), i=1, nk3lines)
        write(outfileindex, 203)k3line_name(nk3lines+1), k3line_stop(nk3lines+1)

        do i=1, nk3lines-1
           write(outfileindex, 204)k3line_stop(i+1), emin, k3line_stop(i+1), emax
        enddo
        write(outfileindex, '(2a)')"plot 'bulkek.dat' u 1:2 ",  &
            "w lp lw 2 pt 7  ps 0.2, \"
        write(outfileindex, '(2a)')"     'bulkek.dat' u 1:2:($3/6):($4/6) ",  &
            "w vec"
        close(outfileindex)
     endif

     202 format('set xtics (',:20('"',A3,'" ',F10.5,','))
     203 format(A3,'" ',F10.5,')')
     204 format('set arrow from ',F10.5,',',F10.5,' to ',F10.5,',',F10.5, ' nohead')
     
   return
   end subroutine ek_bulk_spin

   subroutine ek_bulk_fortomas

     use wmpi
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
        call ham_bulk_old(k, Hamk_bulk)

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


  subroutine ek_bulk_mirror_z

     use wmpi
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
     complex(Dp) :: Hamk     (Num_wann,Num_wann) 
     complex(Dp) :: mat1     (Num_wann,Num_wann) 
     complex(Dp) :: mat2     (Num_wann,Num_wann) 

     ! eigen value of H
	  real(dp), allocatable :: eigv(:,:)
	  real(dp), allocatable :: eigv_mpi(:,:)
     logical, allocatable :: mirror_plus(:, :)
     logical, allocatable :: mirror_minus(:, :)

     knv3= nk3_band
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
        call ham_bulk_old    (k, Hamk_bulk)

        k = k3points(:, ik)
        k(3)= -k(3)
        Hamk= 0d0
        call ham_bulk_old    (k, Hamk)

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


     if (cpuid==0)then
        open(unit=14, file='bulkek.dat')
        do i=1, Num_wann
           do ik=1, knv3
              write(14, '(1000f19.9)')k3len(ik),eigv_mpi(i, ik)
           enddo
           write(14, *)' '
        enddo
        close(14)

        open(unit=15, file='bulkekmirrorplus.dat')
        open(unit=16, file='bulkekmirrorminus.dat')
        do i=1, Num_wann
           do ik=1, knv3
              if (mirror_plus(i, ik)) then
                 write(15, '(1000f19.9)')k3len(ik),eigv_mpi(i, ik)
              else
                 write(16, '(1000f19.9)')k3len(ik),eigv_mpi(i, ik)
              endif
           enddo
           write(15, *)' '
           write(16, *)' '
        enddo
        close(15)
        close(16)
 

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
        write(outfileindex, '(a)')'set ylabel "Energy (eV)"'
        write(outfileindex, '(a, f10.5, a)')'set xrange [0: ', maxval(k3len), ']'
        write(outfileindex, '(a, f10.5, a, f10.5, a)')'set yrange [', emin, ':', emax, ']'
        write(outfileindex, 202, advance="no") (k3line_name(i), k3line_stop(i), i=1, nk3lines)
        write(outfileindex, 203)k3line_name(nk3lines+1), k3line_stop(nk3lines+1)

        do i=1, nk3lines-1
           write(outfileindex, 204)k3line_stop(i+1), emin, k3line_stop(i+1), emax
        enddo
        write(outfileindex, '(2a)')"plot 'bulkek.dat' u 1:2 ",  &
            "w lp lw 2 pt 7  ps 0.2, 0 w l"
        close(outfileindex)
     endif

     202 format('set xtics (',:20('"',A3,'" ',F10.5,','))
     203 format(A3,'" ',F10.5,')')
     204 format('set arrow from ',F10.5,',',F10.5,' to ',F10.5,',',F10.5, ' nohead')
     
   return
   end subroutine ek_bulk_mirror_z




  !> calculate bulk band for a given mirror symmetry
  subroutine ek_bulk_mirror_x

     use wmpi
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
     complex(Dp) :: Hamk     (Num_wann,Num_wann) 
     complex(Dp) :: mat1     (Num_wann,Num_wann) 
     complex(Dp) :: mat2     (Num_wann,Num_wann) 
     complex(dp) :: inv_mirror_x(Num_wann, Num_wann)

     ! eigen value of H
	  real(dp), allocatable :: eigv(:,:)
	  real(dp), allocatable :: eigv_mpi(:,:)
     logical, allocatable :: mirror_plus(:, :)
     logical, allocatable :: mirror_minus(:, :)

     knv3= nk3_band
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
        call ham_bulk_old    (k, Hamk_bulk)

       !k = k3points(:, ik)
       !k(1)= -k(1)
       !Hamk= 0d0
       !call ham_bulk_old    (k, Hamk)


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

     if (cpuid==0)then
        open(unit=14, file='bulkek.dat')
        do i=1, Num_wann
           do ik=1, knv3
              write(14, '(1000f19.9)')k3len(ik),eigv_mpi(i, ik)
           enddo
           write(14, *)' '
        enddo
        close(14)

        open(unit=15, file='bulkekmirrorplus.dat')
        open(unit=16, file='bulkekmirrorminus.dat')
        do i=1, Num_wann
           do ik=1, knv3
              if (mirror_plus(i, ik)) then
                 write(15, '(1000f19.9)')k3len(ik),eigv_mpi(i, ik)
              else
                 write(16, '(1000f19.9)')k3len(ik),eigv_mpi(i, ik)
              endif
           enddo
           write(15, *)' '
           write(16, *)' '
        enddo
        close(15)
        close(16)
 

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
        write(outfileindex, '(a)')'set ylabel "Energy (eV)"'
        write(outfileindex, '(a, f10.5, a)')'set xrange [0: ', maxval(k3len), ']'
        write(outfileindex, '(a, f10.5, a, f10.5, a)')'set yrange [', emin, ':', emax, ']'
        write(outfileindex, 202, advance="no") (k3line_name(i), k3line_stop(i), i=1, nk3lines)
        write(outfileindex, 203)k3line_name(nk3lines+1), k3line_stop(nk3lines+1)

        do i=1, nk3lines-1
           write(outfileindex, 204)k3line_stop(i+1), emin, k3line_stop(i+1), emax
        enddo
        write(outfileindex, '(2a)')"plot 'bulkek.dat' u 1:2 ",  &
            "w lp lw 2 pt 7  ps 0.2, 0 w l"
        close(outfileindex)
     endif

     202 format('set xtics (',:20('"',A3,'" ',F10.5,','))
     203 format(A3,'" ',F10.5,')')
     204 format('set arrow from ',F10.5,',',F10.5,' to ',F10.5,',',F10.5, ' nohead')
     
   return
   end subroutine ek_bulk_mirror_x



