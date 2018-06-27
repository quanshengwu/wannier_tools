  subroutine berry_curvarture_singlek_EF(k, Omega_x, Omega_y, Omega_z)
     !> Calculate Berry curvature for a sigle k point
     !> The Fermi distribution is determined by the Fermi level E_arc
     !> ref : Physical Review B 74, 195118(2006)
     !
     !> eqn (34)
     !
     !> Jun. 11 2018 by Quansheng Wu @ ETHZ
     !> Try to write some simple subroutines. 
     !> on the output, only the real part is meaningfull.
     ! Copyright (c) 2018 QuanSheng Wu. All rights reserved.

     use wmpi
     use para
     implicit none

     !> input parameter, k is in unit of reciprocal lattice vectors
     real(dp), intent(in) :: k(3)
     complex(dp), intent(out) :: Omega_x(Num_wann)
     complex(dp), intent(out) :: Omega_y(Num_wann)
     complex(dp), intent(out) :: Omega_z(Num_wann)

     integer :: iR, m, n, i, j
     real(dp), allocatable :: W(:)
     real(dp) :: kdotr, Beta_fake
     complex(dp), allocatable :: Amat(:, :), DHDk(:, :, :), DHDkdag(:, :, :)
     complex(dp), allocatable :: UU(:, :), UU_dag(:, :), Hamk_bulk(:, :)
     complex(dp), allocatable :: vx(:, :), vy(:, :), vz(:, :)

     !> Fermi-Dirac distribution
     real(dp), external :: fermi

     allocate(W(Num_wann))
     allocate(vx(Num_wann, Num_wann), vy(Num_wann, Num_wann), vz(Num_wann, Num_wann))
     allocate(UU(Num_wann, Num_wann), UU_dag(Num_wann, Num_wann), Hamk_bulk(Num_wann, Num_wann))
     allocate(Amat(Num_wann, Num_wann), DHDk(Num_wann, Num_wann, 3), DHDkdag(Num_wann, Num_wann, 3))
     W=0d0; vx= 0d0; vy= 0d0; vz= 0d0; UU= 0d0; UU_dag= 0d0

     ! calculation bulk hamiltonian by a direct Fourier transformation of HmnR
     ! generate bulk Hamiltonian
     if (index(KPorTB, 'KP')/=0)then
        call ham_bulk_kp(k, Hamk_bulk)
     else
       !> deal with phonon system
       if (index(Particle,'phonon')/=0.and.LOTO_correction) then
          call ham_bulk_LOTO(k, Hamk_bulk)
       else
          call ham_bulk_old    (k, Hamk_bulk)
       endif
     endif

     !> diagonalization by call zheev in lapack
     UU=Hamk_bulk
     call eigensystem_c( 'V', 'U', Num_wann, UU, W)
    !call zhpevx_pack(hamk_bulk,Num_wann, W, UU)

     UU_dag= conjg(transpose(UU))
     do iR= 1, Nrpts
        kdotr= k(1)*irvec(1,iR) + k(2)*irvec(2,iR) + k(3)*irvec(3,iR)
        vx= vx+ zi*crvec(1, iR)*HmnR(:,:,iR)*Exp(pi2zi*kdotr)/ndegen(iR)
        vy= vy+ zi*crvec(2, iR)*HmnR(:,:,iR)*Exp(pi2zi*kdotr)/ndegen(iR)
        vz= vz+ zi*crvec(3, iR)*HmnR(:,:,iR)*Exp(pi2zi*kdotr)/ndegen(iR)
     enddo ! iR

     !> unitility rotate velocity
     call mat_mul(Num_wann, vx, UU, Amat) 
     call mat_mul(Num_wann, UU_dag, Amat, vx) 
     call mat_mul(Num_wann, vy, UU, Amat) 
     call mat_mul(Num_wann, UU_dag, Amat, vy) 
     call mat_mul(Num_wann, vz, UU, Amat) 
     call mat_mul(Num_wann, UU_dag, Amat, vz) 

     Omega_x=0d0;Omega_y=0d0; Omega_z=0d0
     do m= 1, Num_wann
        do n= 1, Num_wann
           if (m==n) cycle
           Omega_x(m)= Omega_x(m)+ vy(n, m)*vz(m, n)/((W(m)-W(n))**2)
           Omega_y(m)= Omega_y(m)+ vz(n, m)*vx(m, n)/((W(m)-W(n))**2)
           Omega_z(m)= Omega_z(m)+ vx(n, m)*vy(m, n)/((W(m)-W(n))**2)
        enddo ! m
     enddo ! n

     Omega_x= -Omega_x*2d0*zi
     Omega_y= -Omega_y*2d0*zi
     Omega_z= -Omega_z*2d0*zi

     !> consider the Fermi-distribution according to the brodening Earc_eta
     Beta_fake= 1d0/Eta_Arc
     do m= 1, Num_wann
        Omega_x(m)= Omega_x(m)*fermi(W(m), Beta_fake)
        Omega_y(m)= Omega_y(m)*fermi(W(m), Beta_fake)
        Omega_z(m)= Omega_z(m)*fermi(W(m), Beta_fake)
     enddo

     return
  end subroutine berry_curvarture_singlek_EF

  subroutine berry_curvarture_singlek_numoccupied(k, Omega_x, Omega_y, Omega_z)
     !> Calculate Berry curvature for a sigle k point
     !> The Fermi distribution is determined by the NumOccupied bands, not the Fermi level
     !> ref : Physical Review B 74, 195118(2006)
     !
     !> eqn (30), we only calculate yz, zx, xy
     !
     !> Jun. 25 2018 by Quansheng Wu @ airplane CA781 from Beijing to Zurich
     !> Try to write some simple subroutines. 
     ! Copyright (c) 2018 QuanSheng Wu. All rights reserved.

     use wmpi
     use para
     implicit none

     !> input parameter, k is in unit of reciprocal lattice vectors
     real(dp), intent(in) :: k(3)
     complex(dp), intent(out) :: Omega_x(Num_wann)
     complex(dp), intent(out) :: Omega_y(Num_wann)
     complex(dp), intent(out) :: Omega_z(Num_wann)

     integer :: iR, m, n, i, j
     real(dp), allocatable :: W(:)
     real(dp) :: kdotr
     complex(dp), allocatable :: Amat(:, :), DHDk(:, :, :), DHDkdag(:, :, :)
     complex(dp), allocatable :: UU(:, :), UU_dag(:, :), Hamk_bulk(:, :)
     complex(dp), allocatable :: vx(:, :), vy(:, :), vz(:, :)

     allocate(W(Num_wann))
     allocate(vx(Num_wann, Num_wann), vy(Num_wann, Num_wann), vz(Num_wann, Num_wann))
     allocate(UU(Num_wann, Num_wann), UU_dag(Num_wann, Num_wann), Hamk_bulk(Num_wann, Num_wann))
     allocate(Amat(Num_wann, Num_wann), DHDk(Num_wann, Num_wann, 3), DHDkdag(Num_wann, Num_wann, 3))
     W=0d0; vx= 0d0; vy= 0d0; vz= 0d0; UU= 0d0; UU_dag= 0d0

     ! calculation bulk hamiltonian by a direct Fourier transformation of HmnR
     if (index(KPorTB, 'KP')/=0)then
        call ham_bulk_kp(k, Hamk_bulk)
     else
       !> deal with phonon system
       if (index(Particle,'phonon')/=0.and.LOTO_correction) then
          call ham_bulk_LOTO(k, Hamk_bulk)
       else
          call ham_bulk_old    (k, Hamk_bulk)
       endif
     endif

     !> diagonalization by call zheev in lapack
     UU=Hamk_bulk
     call eigensystem_c( 'V', 'U', Num_wann, UU, W)
    !call zhpevx_pack(hamk_bulk,Num_wann, W, UU)

     UU_dag= conjg(transpose(UU))
     do iR= 1, Nrpts
        kdotr= k(1)*irvec(1,iR) + k(2)*irvec(2,iR) + k(3)*irvec(3,iR)
        vx= vx+ zi*crvec(1, iR)*HmnR(:,:,iR)*Exp(pi2zi*kdotr)/ndegen(iR)
        vy= vy+ zi*crvec(2, iR)*HmnR(:,:,iR)*Exp(pi2zi*kdotr)/ndegen(iR)
        vz= vz+ zi*crvec(3, iR)*HmnR(:,:,iR)*Exp(pi2zi*kdotr)/ndegen(iR)
     enddo ! iR

     !> unitility rotate velocity
     call mat_mul(Num_wann, vx, UU, Amat) 
     call mat_mul(Num_wann, UU_dag, Amat, vx) 
     call mat_mul(Num_wann, vy, UU, Amat) 
     call mat_mul(Num_wann, UU_dag, Amat, vy) 
     call mat_mul(Num_wann, vz, UU, Amat) 
     call mat_mul(Num_wann, UU_dag, Amat, vz) 

     Omega_x=0d0;Omega_y=0d0; Omega_z=0d0
     do m= 1, NumOccupied
        do n= 1, Num_wann
           if (m==n) cycle
           Omega_x(m)= Omega_x(m)+ vy(n, m)*vz(m, n)/((W(m)-W(n))**2)
           Omega_y(m)= Omega_y(m)+ vz(n, m)*vx(m, n)/((W(m)-W(n))**2)
           Omega_z(m)= Omega_z(m)+ vx(n, m)*vy(m, n)/((W(m)-W(n))**2)
        enddo ! m
     enddo ! n

     Omega_x= -Omega_x*2d0*zi
     Omega_y= -Omega_y*2d0*zi
     Omega_z= -Omega_z*2d0*zi

     return
  end subroutine berry_curvarture_singlek_numoccupied

  subroutine berry_curvarture
     !> Calculate Berry curvature 
     !
     !> ref : Physical Review B 74, 195118(2006)
     !
     !> eqn (34)
     !
     !> Sep. 22 2015 by Quansheng Wu @ ETHZ
     !
     ! Copyright (c) 2010 QuanSheng Wu. All rights reserved.

     use wmpi
     use para
     implicit none
    
     integer :: ik, m, n, i, j, ierr

     real(dp) :: kdotr, k(3), o1(3)

     !> k points slice
     real(dp), allocatable :: kslice(:, :), kslice_shape(:, :)
   
     real(dp), external :: norm

     !> velocities
     real(dp), allocatable :: vx(:), vy(:), vz(:)
    
     !> Berry curvature  (3, bands, k)
     complex(dp), allocatable :: Omega_x(:), Omega_y(:), Omega_z(:)
     complex(dp), allocatable :: Omega(:, :), Omega_mpi(:, :)

     allocate( kslice(3, Nk1*Nk2))
     allocate( kslice_shape(3, Nk1*Nk2))
     allocate( Omega_x(Num_wann))
     allocate( Omega_y(Num_wann))
     allocate( Omega_z(Num_wann))
     allocate( Omega    (3, Nk1*Nk2))
     allocate( Omega_mpi(3, Nk1*Nk2))
     allocate( vx      (Num_wann))
     allocate( vy      (Num_wann))
     allocate( vz      (Num_wann))
     kslice=0d0
     kslice_shape=0d0
     omega= 0d0
     omega_mpi= 0d0
     vx=0d0
     vy=0d0
     vz=0d0
    
     !> kslice is centered at K3d_start
     ik =0
     do i= 1, nk1
        do j= 1, nk2
           ik =ik +1
           kslice(:, ik)= K3D_start+ K3D_vec1*(i-1)/dble(nk1-1)  &
                     + K3D_vec2*(j-1)/dble(nk2-1) - (K3D_vec1+ K3D_vec2)/2d0
           kslice_shape(:, ik)= kslice(1, ik)* Kua+ kslice(2, ik)* Kub+ kslice(3, ik)* Kuc 
        enddo
     enddo

     do ik= 1+ cpuid, Nk1*Nk2, num_cpu
        if (cpuid==0) write(stdout, *)'Berry curvature ik, nk1*nk2 ', ik, Nk1*Nk2

        !> diagonalize hamiltonian
        k= kslice(:, ik)

        Omega_x= 0d0
        Omega_y= 0d0
        Omega_z= 0d0

        call berry_curvarture_singlek_numoccupied(k, Omega_x, Omega_y, Omega_z)
       !call berry_curvarture_singlek_EF(k, Omega_x, Omega_y, Omega_z)
        Omega(1, ik) = sum(Omega_x)
        Omega(2, ik) = sum(Omega_y)
        Omega(3, ik) = sum(Omega_z)

     enddo ! ik

     Omega_mpi= 0d0

#if defined (MPI)
     call mpi_allreduce(Omega,Omega_mpi,size(Omega_mpi),&
                       mpi_dc,mpi_sum,mpi_cmw,ierr)
#else
     Omega_mpi= Omega
#endif

     !> output the Berry curvature to file
     outfileindex= outfileindex+ 1
     if (cpuid==0) then
        open(unit=outfileindex, file='Berrycurvature.dat')
        write(outfileindex, '(20a28)')'# Please take the real part for your use'
        write(outfileindex, '(20a28)')'# kx (1/A)', 'ky (1/A)', 'kz (1/A)', &
           'real(Omega_x)', 'imag(Omega_x)', &
           'real(Omega_y)', 'imag(Omega_y)', &
           'real(Omega_z)', 'imag(Omega_z)'
        ik= 0
        do i= 1, nk1
           do j= 1, nk2
              ik= ik+ 1
              write(outfileindex, '(20E28.10)')kslice_shape(:, ik), Omega_mpi(:, ik)
           enddo
           write(outfileindex, *) ' '
        enddo

        close(outfileindex)

     endif

     !> output the Berry curvature to file
     outfileindex= outfileindex+ 1
     if (cpuid==0) then
        open(unit=outfileindex, file='Berrycurvature-normalized.dat')
        write(outfileindex, '(20a28)')'# Please take the real part for your use'
        write(outfileindex, '(20a28)')'# kx (1/A)', 'ky (1/A)', 'kz (1/A)', &
           'real(Omega_x)',  &
           'real(Omega_y)',  &
           'real(Omega_z)'
        ik= 0
        do i= 1, nk1
           do j= 1, nk2
              ik= ik+ 1
              o1=real(Omega_mpi(:,ik))
              if (norm(o1)>eps9) o1= o1/norm(o1)
              write(outfileindex, '(20f28.10)')kslice_shape(:, ik), o1
           enddo
           write(outfileindex, *) ' '
        enddo

        close(outfileindex)

     endif

     !> generate gnuplot script to plot the Berry curvature
     outfileindex= outfileindex+ 1
     if (cpuid==0) then

        open(unit=outfileindex, file='Berrycurvature.gnu')
        write(outfileindex, '(a)')"set encoding iso_8859_1"
        write(outfileindex, '(a)')'#set terminal  pngcairo  truecolor enhanced size 1920, 1680 font ",40"'
        write(outfileindex, '(a)')'set terminal  png       truecolor enhanced size 1920, 1680 font ",40"'
        write(outfileindex, '(a)')"set output 'Berrycurvature.png'"
        write(outfileindex, '(a)')'if (!exists("MP_LEFT"))   MP_LEFT = .12'
        write(outfileindex, '(a)')'if (!exists("MP_RIGHT"))  MP_RIGHT = .92'
        write(outfileindex, '(a)')'if (!exists("MP_BOTTOM")) MP_BOTTOM = .12'
        write(outfileindex, '(a)')'if (!exists("MP_TOP"))    MP_TOP = .88'
        write(outfileindex, '(a)')'if (!exists("MP_GAP"))    MP_GAP = 0.08'
        write(outfileindex, '(a)')'set multiplot layout 2,3 rowsfirst title "\{ Berry Curvature\}" \'
        write(outfileindex, '(a)')"              margins screen MP_LEFT, MP_RIGHT, MP_BOTTOM, MP_TOP spacing screen MP_GAP"
        write(outfileindex, '(a)')" "
        write(outfileindex, '(a)')"set palette rgbformulae 33,13,10"
        write(outfileindex, '(a)')"unset ztics"
        write(outfileindex, '(a)')"unset key"
        write(outfileindex, '(a)')"set pm3d"
        write(outfileindex, '(a)')"set view map"
        write(outfileindex, '(a)')"set border lw 3"
        write(outfileindex, '(a)')"set xlabel 'k (1/{\305})'"
        write(outfileindex, '(a)')"set ylabel 'k (1/{\305})'"
        write(outfileindex, '(a)')"unset colorbox"
        write(outfileindex, '(a)')"unset xtics"
        write(outfileindex, '(a)')"unset xlabel"
        write(outfileindex, '(a)')"set xrange [] noextend"
        write(outfileindex, '(a)')"set yrange [] noextend"
        write(outfileindex, '(a)')"set ytics 0.5 nomirror scale 0.5"
        write(outfileindex, '(a)')"set pm3d interpolate 2,2"
        write(outfileindex, '(a)')"set title 'Omega_x real'"
        write(outfileindex, '(a)')"splot 'Berrycurvature.dat' u 1:2:4 w pm3d"
        write(outfileindex, '(a)')"unset ylabel"
        write(outfileindex, '(a)')"unset ytics"
        write(outfileindex, '(a)')"set title 'Omega_y real'"
        write(outfileindex, '(a)')"splot 'Berrycurvature.dat' u 1:2:6 w pm3d"
        write(outfileindex, '(a)')"set title 'Omega_z real'"
        write(outfileindex, '(a)')"set colorbox"
        write(outfileindex, '(a)')"splot 'Berrycurvature.dat' u 1:2:8 w pm3d"
 
        write(outfileindex, '(a)')"set xtics 0.5 nomirror scale 0.5"
        write(outfileindex, '(a)')"set ytics nomirror scale 0.5"
        write(outfileindex, '(a)')"set xlabel 'k (1/{\305})'"
        write(outfileindex, '(a)')"set ylabel 'k (1/{\305})'"
        write(outfileindex, '(a)')"unset colorbox"
        write(outfileindex, '(a)')"set title 'Omega_x imag'"
        write(outfileindex, '(a)')"splot 'Berrycurvature.dat' u 1:2:5 w pm3d"
        write(outfileindex, '(a)')"unset ytics"
        write(outfileindex, '(a)')"set title 'Omega_y imag'"
        write(outfileindex, '(a)')"splot 'Berrycurvature.dat' u 1:2:7 w pm3d"
        write(outfileindex, '(a)')"set title 'Omega_z imag'"
        write(outfileindex, '(a)')"set colorbox"
        write(outfileindex, '(a)')"splot 'Berrycurvature.dat' u 1:2:9 w pm3d"
        close(outfileindex)
     endif


     !> generate gnuplot script to plot the Berry curvature
     outfileindex= outfileindex+ 1
     if (cpuid==0) then

        open(unit=outfileindex, file='Berrycurvature-normalized.gnu')
        write(outfileindex, '(a)')"set encoding iso_8859_1"
        write(outfileindex, '(a)')'#set terminal  pngcairo  truecolor enhanced size 1920, 1680 font ",40"'
        write(outfileindex, '(a)')'set terminal  png       truecolor enhanced size 1920, 1680 font ",40"'
        write(outfileindex, '(a)')"set output 'Berrycurvature-normalized.png'"
        write(outfileindex, '(a)')"unset ztics"
        write(outfileindex, '(a)')"unset key"
        write(outfileindex, '(a)')"set border lw 3"
        write(outfileindex, '(a)')"set xlabel 'k (1/{\305})'"
        write(outfileindex, '(a)')"set ylabel 'k (1/{\305})'"
        write(outfileindex, '(a)')"unset colorbox"
        write(outfileindex, '(a)')"unset xtics"
        write(outfileindex, '(a)')"unset xlabel"
        write(outfileindex, '(a)')"set xrange [] noextend"
        write(outfileindex, '(a)')"set yrange [] noextend"
        write(outfileindex, '(a)')"set ytics 0.5 nomirror scale 0.5"
        write(outfileindex, '(a)')"set pm3d interpolate 2,2"
        write(outfileindex, '(a)')"set title '(Omega_x, Omega_y) real'"
        write(outfileindex, '(a)')"plot 'Berrycurvature-normalized.dat' u 1:2:($4/500):($5/500) w vec"
        close(outfileindex)
     endif

#if defined (MPI)
     call mpi_barrier(mpi_cmw, ierr)
#endif

     deallocate( kslice)
     deallocate( kslice_shape)
     deallocate( Omega_x, Omega_y, Omega_z)
     deallocate( Omega, Omega_mpi)
 
     return

  end subroutine berry_curvarture

  subroutine Fourier_R_to_k(k, ham)
     !> Fourier transform the Hamiltonian from R space to k space
     use para, only: irvec, HmnR, Nrpts, ndegen, pi, zi, Num_wann, dp
     implicit none

     real(dp), intent(in) :: k(3)
     complex(dp), intent(out) :: ham(Num_wann, Num_wann)
     integer :: iR
    !real(dp) :: R(3)
     real(dp) :: kdotr

     ham= 0d0
     do iR= 1, Nrpts
       !R= Rua*irvec(1,iR) + Rub*irvec(2,iR) + Ruc*irvec(3,iR)
        kdotr= k(1)*irvec(1,iR) + k(2)*irvec(2,iR) + k(3)*irvec(3,iR)
        Ham= Ham+ HmnR(:,:,iR)*Exp(2d0*pi*zi*kdotr)/ndegen(iR)
     enddo

     return
  end subroutine Fourier_R_to_k

  subroutine Im_trace(ndim, A, tr)
     !> Calculate trace only with the imaginary part of a matrix A with dimension ndim
     use para, only : dp
     implicit none
     integer :: ndim
     complex(dp), intent(out) :: tr
     complex(dp), intent(in) :: A(ndim, ndim)

     integer :: i

     tr = 0d0
     do i=1, ndim
        tr= tr+ aimag(A(i, i))
     enddo

     return
  end subroutine Im_trace

  subroutine trace(ndim, A, tr)
     !> Calculate trace of a matrix A with dimension ndim
     use para, only : dp
     use para, only : dp
     implicit none
     integer :: ndim
     complex(dp), intent(out) :: tr
     complex(dp), intent(in) :: A(ndim, ndim)

     integer :: i

     tr = 0d0
     do i=1, ndim
        tr= tr+ (A(i, i))
     enddo

     return
  end subroutine trace
