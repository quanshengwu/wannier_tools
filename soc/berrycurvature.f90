  subroutine berry_curvarture_singlek_EF_old(k, Omega_x, Omega_y, Omega_z)
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
     call ham_bulk_latticegauge(k, Hamk_bulk)

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

     DHDk= 0d0
     DHDkdag= 0d0
     do m= 1, Num_wann
        do n= 1, Num_wann
           if (W(n) > E_arc .and. W(m)<E_arc) then
              DHDk(n, m, 1)= zi*vx(n, m)/(W(m)-W(n))
              DHDk(n, m, 2)= zi*vy(n, m)/(W(m)-W(n))
              DHDk(n, m, 3)= zi*vz(n, m)/(W(m)-W(n))
           else
              DHDk(n, m, 1)= 0d0
              DHDk(n, m, 2)= 0d0
              DHDk(n, m, 3)= 0d0
           endif
        enddo ! m
     enddo ! n

     do m= 1, Num_wann
        do n= 1, Num_wann
           if (W(m) > E_arc .and. W(n)< E_arc) then
              DHDkdag(n, m, 1)= zi*vx(n, m)/(W(m)-W(n))
              DHDkdag(n, m, 2)= zi*vy(n, m)/(W(m)-W(n))
              DHDkdag(n, m, 3)= zi*vz(n, m)/(W(m)-W(n))
           else              
              DHDkdag(n, m, 1)= 0d0
              DHDkdag(n, m, 2)= 0d0
              DHDkdag(n, m, 3)= 0d0
           endif
        enddo ! m
     enddo ! n

     !> rotate DHDk and DHDkdag to diagonal basis
     do i=1, 3
        call mat_mul(Num_wann, DHDk(:, :, i), UU_dag, Amat) 
        call mat_mul(Num_wann, UU, Amat, DHDk(:, :, i)) 
        call mat_mul(Num_wann, DHDkdag(:, :, i), UU_dag, Amat) 
        call mat_mul(Num_wann, UU, Amat, DHDkdag(:, :, i)) 
     enddo

     call mat_mul(Num_wann, DHDk(:, :, 1), DHDkdag(:, :, 2), vz)
     call mat_mul(Num_wann, DHDk(:, :, 2), DHDkdag(:, :, 3), vx)
     call mat_mul(Num_wann, DHDk(:, :, 3), DHDkdag(:, :, 1), vy)

     do i=1, Num_wann
        Omega_x(i)= -2d0*zi*vx(i, i)
        Omega_y(i)= -2d0*zi*vy(i, i)
        Omega_z(i)= -2d0*zi*vz(i, i)
     enddo

     !> add a factor of 2
     Omega_x= Omega_x*2d0
     Omega_y= Omega_y*2d0
     Omega_z= Omega_z*2d0

     return
  end subroutine berry_curvarture_singlek_EF_old

  subroutine berry_curvarture_singlek_EF(k, mu, Omega_x, Omega_y, Omega_z)
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
     !> mu: chemical potential which can be tuned by the back gate
     real(dp), intent(in) :: k(3), mu
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
     call ham_bulk_latticegauge(k, Hamk_bulk)

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
        Omega_x(m)= Omega_x(m)*fermi(W(m)-mu, Beta_fake)
        Omega_y(m)= Omega_y(m)*fermi(W(m)-mu, Beta_fake)
        Omega_z(m)= Omega_z(m)*fermi(W(m)-mu, Beta_fake)
     enddo

     return
  end subroutine berry_curvarture_singlek_EF



  subroutine berry_curvarture_singlek_numoccupied_old(k, Omega_x, Omega_y, Omega_z)
     !> Calculate Berry curvature for a sigle k point
     !> The Fermi distribution is determined by the NumOccupied bands, not the Fermi level
     !> ref : Physical Review B 74, 195118(2006)
     !
     !> eqn (34)
     !
     !> Jun. 11 2018 by Quansheng Wu @ ETHZ
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
     call ham_bulk_latticegauge(k, Hamk_bulk)

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

     DHDk= 0d0
     DHDkdag= 0d0
     do m= 1, Num_wann
        do n= 1, Num_wann
           if (n> Numoccupied .and. m<= Numoccupied) then 
              DHDk(n, m, 1)= zi*vx(n, m)/(W(m)-W(n))
              DHDk(n, m, 2)= zi*vy(n, m)/(W(m)-W(n))
              DHDk(n, m, 3)= zi*vz(n, m)/(W(m)-W(n))
           else
              DHDk(n, m, 1)= 0d0
              DHDk(n, m, 2)= 0d0
              DHDk(n, m, 3)= 0d0
           endif
        enddo ! m
     enddo ! n

     do m= 1, Num_wann
        do n= 1, Num_wann
           if (m>Numoccupied .and. n<=Numoccupied) then
              DHDkdag(n, m, 1)= zi*vx(n, m)/(W(m)-W(n))
              DHDkdag(n, m, 2)= zi*vy(n, m)/(W(m)-W(n))
              DHDkdag(n, m, 3)= zi*vz(n, m)/(W(m)-W(n))
           else              
              DHDkdag(n, m, 1)= 0d0
              DHDkdag(n, m, 2)= 0d0
              DHDkdag(n, m, 3)= 0d0
           endif
        enddo ! m
     enddo ! n

     !> rotate DHDk and DHDkdag to diagonal basis
     do i=1, 3
        call mat_mul(Num_wann, DHDk(:, :, i), UU_dag, Amat) 
        call mat_mul(Num_wann, UU, Amat, DHDk(:, :, i)) 
        call mat_mul(Num_wann, DHDkdag(:, :, i), UU_dag, Amat) 
        call mat_mul(Num_wann, UU, Amat, DHDkdag(:, :, i)) 
     enddo

     call mat_mul(Num_wann, DHDk(:, :, 1), DHDkdag(:, :, 2), vz)
     call mat_mul(Num_wann, DHDk(:, :, 2), DHDkdag(:, :, 3), vx)
     call mat_mul(Num_wann, DHDk(:, :, 3), DHDkdag(:, :, 1), vy)

     do i=1, Num_wann
        Omega_x(i)= -2d0*zi*vx(i, i)
        Omega_y(i)= -2d0*zi*vy(i, i)
        Omega_z(i)= -2d0*zi*vz(i, i)
     enddo

     Omega_x= real(Omega_x*2d0)
     Omega_y= real(Omega_y*2d0)
     Omega_z= real(Omega_z*2d0)

     return
  end subroutine berry_curvarture_singlek_numoccupied_old

  subroutine berry_curvarture_singlek_numoccupied_slab_total(k, Omega_z)
     !> Calculate Berry curvature for a sigle k point
     !> The Fermi distribution is determined by the NumOccupied bands, not the Fermi level
     !> ref : Physical Review B 74, 195118(2006)
     !
     !> eqn (11), we only calculate xy
     !
     !> Aug. 06 2018 by Quansheng Wu 
     !> Try to write some simple subroutines. 
     ! Copyright (c) 2018 QuanSheng Wu. All rights reserved.

     use wmpi
     use para
     implicit none

     !> input parameter, k is in unit of reciprocal lattice vectors
     real(dp), intent(in) :: k(2)
     complex(dp), intent(out) :: Omega_z(1)

     integer :: iR, m, n, i, j, Mdim, i1, i2
     real(dp), allocatable :: W(:)
     real(dp) :: kdotr
     complex(dp), allocatable :: Amat(:, :), DHDk(:, :, :), DHDkdag(:, :, :)
     complex(dp), allocatable :: UU(:, :), UU_dag(:, :), Hamk_slab(:, :)
     complex(dp), allocatable :: vx(:, :), vy(:, :)
     complex(dp), allocatable :: Vij_x(:, :, :), Vij_y(:, :, :)

     !> leading order of the hamiltonian matrix for slab system specified with Nslab
     Mdim = Num_wann*Nslab

     allocate(W(Mdim))
     allocate(vx(Mdim, Mdim), vy(Mdim, Mdim))
     allocate(UU(Mdim, Mdim), UU_dag(Mdim, Mdim), Hamk_slab(Mdim, Mdim))
     allocate(Amat(Mdim, Mdim), DHDk(Mdim, Mdim, 3), DHDkdag(Mdim, Mdim, 3))
     allocate(Vij_x(-ijmax:ijmax, Num_wann, Num_wann))
     allocate(Vij_y(-ijmax:ijmax, Num_wann, Num_wann))
     W=0d0; vx= 0d0; vy= 0d0; UU= 0d0; UU_dag= 0d0

     call ham_qlayer2qlayer_velocity(k, Vij_x, Vij_y) 
     do i1=1, nslab
        ! i2 row index
        do i2=1, nslab
          if (abs(i2-i1).le.ijmax)then
            vx((i2-1)*Num_wann+1:(i2-1)*Num_wann+Num_wann,&
                      (i1-1)*Num_wann+1:(i1-1)*Num_wann+Num_wann )&
            = Vij_x(i2-i1,1:Num_wann,1:Num_wann)
            vy((i2-1)*Num_wann+1:(i2-1)*Num_wann+Num_wann,&
                      (i1-1)*Num_wann+1:(i1-1)*Num_wann+Num_wann )&
            = Vij_y(i2-i1,1:Num_wann,1:Num_wann)
          endif 
        enddo ! i2
     enddo ! i1

     ! calculation slab hamiltonian by a direct Fourier transformation of HmnR
     call ham_slab(k, Hamk_slab)

     !> diagonalization by call zheev in lapack
     UU=Hamk_slab
     call eigensystem_c( 'V', 'U', Mdim, UU, W)
    !call zhpevx_pack(hamk_slab,Mdim, W, UU)

     UU_dag= conjg(transpose(UU))
 
     !> unitility rotate velocity
     call mat_mul(Mdim, vx, UU, Amat) 
     call mat_mul(Mdim, UU_dag, Amat, vx) 
     call mat_mul(Mdim, vy, UU, Amat) 
     call mat_mul(Mdim, UU_dag, Amat, vy) 

     Omega_z=0d0
     do m= 1, NumOccupied*Nslab
        do n= NumOccupied*Nslab+1, Mdim
           Omega_z(1)= Omega_z(1)+ vx(n, m)*vy(m, n)/((W(m)-W(n))**2)
        enddo ! m
     enddo ! n

     Omega_z= -aimag(Omega_z*2d0)

     return
  end subroutine berry_curvarture_singlek_numoccupied_slab_total

  subroutine berry_curvarture_singlek_numoccupied_total(k, Omega_x, Omega_y, Omega_z)
     !> Calculate Berry curvature for a sigle k point using Kubo-formula
     !> The Fermi distribution is determined by the NumOccupied bands, not the Fermi level
     !> ref : Physical Review B 74, 195118(2006)
     !
     !> eqn (11), we only calculate yz, zx, xy
     !
     !> Jun. 25 2018 by Quansheng Wu @ airplane CA781 from Beijing to Zurich
     !> Try to write some simple subroutines. 
     ! Copyright (c) 2018 QuanSheng Wu. All rights reserved.

     use wmpi
     use para
     implicit none

     !> input parameter, k is in unit of reciprocal lattice vectors
     real(dp), intent(in) :: k(3)
     complex(dp), intent(out) :: Omega_x
     complex(dp), intent(out) :: Omega_y
     complex(dp), intent(out) :: Omega_z

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
     !call ham_bulk_latticegauge(k, Hamk_bulk)
     call ham_bulk    (k, Hamk_bulk)

     !> diagonalization by call zheev in lapack
     UU=Hamk_bulk
     call eigensystem_c( 'V', 'U', Num_wann, UU, W)
    !call zhpevx_pack(hamk_bulk,Num_wann, W, UU)

     UU_dag= conjg(transpose(UU))
    !do iR= 1, Nrpts
    !   kdotr= k(1)*irvec(1,iR) + k(2)*irvec(2,iR) + k(3)*irvec(3,iR)
    !   vx= vx+ zi*crvec(1, iR)*HmnR(:,:,iR)*Exp(pi2zi*kdotr)/ndegen(iR)
    !   vy= vy+ zi*crvec(2, iR)*HmnR(:,:,iR)*Exp(pi2zi*kdotr)/ndegen(iR)
    !   vz= vz+ zi*crvec(3, iR)*HmnR(:,:,iR)*Exp(pi2zi*kdotr)/ndegen(iR)
    !enddo ! iR

     !call dHdk_latticegauge(k, vx, vy, vz)
     call dHdk_atomicgauge(k, vx, vy, vz)

     !> unitility rotate velocity
     call mat_mul(Num_wann, vx, UU, Amat) 
     call mat_mul(Num_wann, UU_dag, Amat, vx) 
     call mat_mul(Num_wann, vy, UU, Amat) 
     call mat_mul(Num_wann, UU_dag, Amat, vy) 
     call mat_mul(Num_wann, vz, UU, Amat) 
     call mat_mul(Num_wann, UU_dag, Amat, vz) 

     Omega_x=0d0; Omega_y=0d0; Omega_z=0d0
     do m= 1, NumOccupied !> sum over valence band
        do n= NumOccupied+1, Num_wann !> sum over conduction band
           Omega_x= Omega_x+ vy(n, m)*vz(m, n)/((W(m)-W(n))**2)
           Omega_y= Omega_y+ vz(n, m)*vx(m, n)/((W(m)-W(n))**2)
           Omega_z= Omega_z+ vx(n, m)*vy(m, n)/((W(m)-W(n))**2)
        enddo ! m
     enddo ! n

     Omega_x= -aimag(Omega_x*2d0)  
     Omega_y= -aimag(Omega_y*2d0)  
     Omega_z= -aimag(Omega_z*2d0)  

     return
  end subroutine berry_curvarture_singlek_numoccupied_total

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
     call ham_bulk    (k, Hamk_bulk)
     !call ham_bulk_latticegauge(k, Hamk_bulk)

     !> diagonalization by call zheev in lapack
     UU=Hamk_bulk
     call eigensystem_c( 'V', 'U', Num_wann, UU, W)
    !call zhpevx_pack(hamk_bulk,Num_wann, W, UU)

     UU_dag= conjg(transpose(UU))
    !do iR= 1, Nrpts
    !   kdotr= k(1)*irvec(1,iR) + k(2)*irvec(2,iR) + k(3)*irvec(3,iR)
    !   vx= vx+ zi*crvec(1, iR)*HmnR(:,:,iR)*Exp(pi2zi*kdotr)/ndegen(iR)
    !   vy= vy+ zi*crvec(2, iR)*HmnR(:,:,iR)*Exp(pi2zi*kdotr)/ndegen(iR)
    !   vz= vz+ zi*crvec(3, iR)*HmnR(:,:,iR)*Exp(pi2zi*kdotr)/ndegen(iR)
    !enddo ! iR

     !call dHdk_latticegauge(k, vx, vy, vz)
     call dHdk_atomicgauge(k, vx, vy, vz)

     !> unitility rotate velocity
     call mat_mul(Num_wann, vx, UU, Amat) 
     call mat_mul(Num_wann, UU_dag, Amat, vx) 
     call mat_mul(Num_wann, vy, UU, Amat) 
     call mat_mul(Num_wann, UU_dag, Amat, vy) 
     call mat_mul(Num_wann, vz, UU, Amat) 
     call mat_mul(Num_wann, UU_dag, Amat, vz) 

     Omega_x=0d0;Omega_y=0d0; Omega_z=0d0
     do m= NumOccupied, NumOccupied
        do n= NumOccupied+1, NumOccupied+1
    !do m= 1, NumOccupied
    !   do n= 1, Num_wann
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

  subroutine berry_curvarture_slab
     !> Calculate Berry curvature 
     !
     !> ref : Physical Review B 74, 195118(2006)
     !
     !> Aug. 06 2018 by Quansheng Wu @ EPFL
     !
     ! Copyright (c) 2018 QuanSheng Wu. All rights reserved.

     use wmpi
     use para
     implicit none
    
     integer :: ik, m, n, i, j, ierr, Mdim

     real(dp) :: kdotr, k(2), o1(3)

     !> k points slice
     real(dp), allocatable :: k12(:, :)
     real(dp), allocatable :: k12_shape(:, :)
   
     real(dp), external :: norm

     !> Berry curvature  (k)
     complex(dp), allocatable :: Omega_z(:)
     complex(dp), allocatable :: Omega(:)
     complex(dp), allocatable :: Omega_mpi(:)

     Mdim = Num_wann*Nslab

     allocate( k12(2, Nk1*Nk2))
     allocate( k12_shape(2, Nk1*Nk2))
     allocate( Omega_z(Mdim))
     allocate( Omega    (Nk1*Nk2))
     allocate( Omega_mpi(Nk1*Nk2))
     k12=0d0
     k12_shape=0d0
     omega= 0d0
     omega_mpi= 0d0
    
     !> k12 is centered at K2d_start
     ik=0
     do i= 1, nk1
        do j= 1, nk2
           ik=ik+1
           k12(:, ik)=K2D_start+ (i-1)*K2D_vec1/dble(nk1-1) &
                      + (j-1)*K2D_vec2/dble(nk2-1)- (K2D_vec1+K2D_vec2)/2d0
           k12_shape(:, ik)= k12(1, ik)* Ka2+ k12(2, ik)* Kb2
        enddo
     enddo

     do ik= 1+ cpuid, Nk1*Nk2, num_cpu
        if (cpuid==0) write(stdout, *)'Berry curvature ik, nk1*nk2 ', ik, Nk1*Nk2

        !> diagonalize hamiltonian
        k= k12(:, ik)

        Omega_z= 0d0

        call berry_curvarture_singlek_numoccupied_slab_total(k, Omega_z(1))
        Omega(ik) = sum(Omega_z)

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
        open(unit=outfileindex, file='Berrycurvature_slab.dat')
        write(outfileindex, '(20a28)')'# Please take the real part for your use'
        write(outfileindex, '(20a28)')'# kx (1/A)', 'ky (1/A)', &
           'Omega_z'
        ik= 0
        do i= 1, nk1
           do j= 1, nk2
              ik= ik+ 1
              write(outfileindex, '(20E28.10)')k12_shape(:, ik), real(Omega_mpi(ik))
           enddo
           write(outfileindex, *) ' '
        enddo

        close(outfileindex)

     endif

     !> generate gnuplot script to plot the Berry curvature
     outfileindex= outfileindex+ 1
     if (cpuid==0) then

        open(unit=outfileindex, file='Berrycurvature_slab.gnu')
        write(outfileindex, '(a)')"set encoding iso_8859_1"
        write(outfileindex, '(a)')'#set terminal  pngcairo  truecolor enhanced size 1920, 1680 font ",60"'
        write(outfileindex, '(a)')'set terminal  png       truecolor enhanced size 1920, 1680 font ",60"'
        write(outfileindex, '(a)')"set output 'Berrycurvature_slab.png'"
        write(outfileindex, '(a)')"set palette rgbformulae 33,13,10"
        write(outfileindex, '(a)')"unset ztics"
        write(outfileindex, '(a)')"set size ratio -1"
        write(outfileindex, '(a)')"set size 0.9, 0.95"
        write(outfileindex, '(a)')"set origin 0.05, 0.02"
        write(outfileindex, '(a)')"unset key"
        write(outfileindex, '(a)')"set pm3d"
        write(outfileindex, '(a)')"#set zbrange [ -10: 10] "
        write(outfileindex, '(a)')"#set cbrange [ -100: 100] "
        write(outfileindex, '(a)')"set view map"
        write(outfileindex, '(a)')"set border lw 3"
        write(outfileindex, '(a)')"set xlabel 'k (1/{\305})'"
        write(outfileindex, '(a)')"set ylabel 'k (1/{\305})'"
        write(outfileindex, '(a)')"set xrange [] noextend"
        write(outfileindex, '(a)')"set yrange [] noextend"
        write(outfileindex, '(a)')"set ytics 0.5 nomirror scale 0.5"
        write(outfileindex, '(a)')"set pm3d interpolate 2,2"
        write(outfileindex, '(a)')"set title 'Berry Curvature {/Symbol W}_z'"
        write(outfileindex, '(a)')"set colorbox"
        write(outfileindex, '(a)')"splot 'Berrycurvature_slab.dat' u 1:2:3 w pm3d"
        close(outfileindex)
     endif


#if defined (MPI)
     call mpi_barrier(mpi_cmw, ierr)
#endif

     deallocate( k12)
     deallocate( k12_shape)
     deallocate( Omega_z)
     deallocate( Omega, Omega_mpi)
 
     return

  end subroutine berry_curvarture_slab


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

     real(dp) :: kdotr, k(3), o1(3), vmin, vmax

     real(dp) :: time_start, time_end, time_start0

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

     time_start= 0d0
     time_start0= 0d0
     call now(time_start0)
     time_start= time_start0
     time_end  = time_start0
     do ik= 1+ cpuid, Nk1*Nk2, num_cpu
        if (cpuid==0.and. mod(ik/num_cpu, 100)==0) &
           write(stdout, '(a, i9, "  /", i10, a, f10.1, "s", a, f10.1, "s")') &
           ' Berry curvature: ik', ik, Nk1*Nk2, ' time left', &
           (nk1*nk2-ik)*(time_end- time_start)/num_cpu, &
           ' time elapsed: ', time_end-time_start0 

        call now(time_start)
        !> diagonalize hamiltonian
        k= kslice(:, ik)

        Omega_x= 0d0
        Omega_y= 0d0
        Omega_z= 0d0

        !call berry_curvarture_singlek_numoccupied_old(k, Omega_x, Omega_y, Omega_z)
        !call berry_curvarture_singlek_numoccupied(k, Omega_x(1), Omega_y(1), Omega_z(1))
        call berry_curvarture_singlek_numoccupied_total(k, Omega_x(1), Omega_y(1), Omega_z(1))
        ! call berry_curvarture_singlek_EF(k, 0d0, Omega_x, Omega_y, Omega_z)
        Omega(1, ik) = sum(Omega_x)
        Omega(2, ik) = sum(Omega_y)
        Omega(3, ik) = sum(Omega_z)

        call now(time_end)
     enddo ! ik

     Omega_mpi= 0d0

#if defined (MPI)
     call mpi_allreduce(Omega,Omega_mpi,size(Omega_mpi),&
                       mpi_dc,mpi_sum,mpi_cmw,ierr)
#else
     Omega_mpi= Omega
#endif
     vmax=maxval(real(Omega_mpi))
     vmin=minval(real(Omega_mpi))

     !> output the Berry curvature to file
     outfileindex= outfileindex+ 1
     if (cpuid==0) then
        open(unit=outfileindex, file='Berrycurvature.dat')
        write(outfileindex, '(20a28)')'# Please take the real part for your use'
        write(outfileindex, '(20a28)')'# kx (1/A)', 'ky (1/A)', 'kz (1/A)', &
           'real(Omega_x)', 'real(Omega_y)', 'real(Omega_z)' 

        ik= 0
        do i= 1, nk1
           do j= 1, nk2
              ik= ik+ 1
              write(outfileindex, '(20E28.10)')kslice_shape(:, ik), real(Omega_mpi(:, ik))
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
           'real(Omega_x)', 'real(Omega_y)', 'real(Omega_z)'

        ik= 0
        do i= 1, nk1
           do j= 1, nk2
              ik= ik+ 1
              o1= real(Omega_mpi(:,ik))
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
        write(outfileindex, '(a)')'#set terminal  pngcairo  truecolor enhanced size 3680, 1920 font ",40"'
        write(outfileindex, '(a)')'set terminal  png       truecolor enhanced size 3680, 1920 font ",40"'
        write(outfileindex, '(a)')"set output 'Berrycurvature.png'"
        write(outfileindex, '(a)')'if (!exists("MP_LEFT"))   MP_LEFT = .12'
        write(outfileindex, '(a)')'if (!exists("MP_RIGHT"))  MP_RIGHT = .92'
        write(outfileindex, '(a)')'if (!exists("MP_BOTTOM")) MP_BOTTOM = .12'
        write(outfileindex, '(a)')'if (!exists("MP_TOP"))    MP_TOP = .88'
        write(outfileindex, '(a)')'if (!exists("MP_GAP"))    MP_GAP = 0.08'
        write(outfileindex, '(a)')'set multiplot layout 1,3 rowsfirst \'
        write(outfileindex, '(a)')"              margins screen MP_LEFT, MP_RIGHT, MP_BOTTOM, MP_TOP spacing screen MP_GAP"
        write(outfileindex, '(a)')" "
        write(outfileindex, '(a)')"set palette rgbformulae 33,13,10"
        write(outfileindex, '(a)')"unset ztics"
        write(outfileindex, '(a)')"unset key"
        write(outfileindex, '(a)')"set pm3d"
        write(outfileindex, '(a)')"#set zbrange [ -10: 10] "
        write(outfileindex, '(a, f10.3, a, f10.3, a)')"set cbrange [ ", vmin, ':', vmax, " ] "
        write(outfileindex, '(a)')"set view map"
        write(outfileindex, '(a)')"set size ratio -1"
        write(outfileindex, '(a)')"set border lw 3"
        write(outfileindex, '(a)')"set xlabel 'k (1/{\305})'"
        write(outfileindex, '(a)')"set ylabel 'k (1/{\305})'"
        write(outfileindex, '(a)')"unset colorbox"
        write(outfileindex, '(a)')"#unset xtics"
        write(outfileindex, '(a)')"#unset xlabel"
        write(outfileindex, '(a)')"set xrange [] noextend"
        write(outfileindex, '(a)')"set yrange [] noextend"
        write(outfileindex, '(a)')"set ytics 0.5 nomirror scale 0.5"
        write(outfileindex, '(a)')"set pm3d interpolate 2,2"
        write(outfileindex, '(a)')"set title 'Berry Curvature {/Symbol W}_x'"
        write(outfileindex, '(a)')"splot 'Berrycurvature.dat' u 1:2:4 w pm3d"
        write(outfileindex, '(a)')"unset ylabel"
        write(outfileindex, '(a)')"unset ytics"
        write(outfileindex, '(a)')"set title 'Berry Curvature {/Symbol W}_y'"
        write(outfileindex, '(a)')"splot 'Berrycurvature.dat' u 1:2:5 w pm3d"
        write(outfileindex, '(a)')"set title 'Berry Curvature {/Symbol W}_z'"
        write(outfileindex, '(a)')"set colorbox"
        write(outfileindex, '(a)')"splot 'Berrycurvature.dat' u 1:2:6 w pm3d"
 
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
        write(outfileindex, '(a)')"set title '({/Symbol W}_x, {/Symbol W}_y) real'"
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


  subroutine Chern_sphere_single(k0, r0, Chern)
     !> calculate Chern number on a sphere by integration over a sphere
     !> QuanSheng Wu at EPFL June 12 2018
     !> wuquansheng@gmail.com
     !> The sphere is defined by a center k0 and a radius in WEYL_CHIRALITY
     !> C= \int d\theta d\phi 
     !> k0 must be in Cartesian coordinates

     use wmpi
     use para
     implicit none
     
     !> inout variables
     real(dp), intent(in) :: k0(3)
     real(dp), intent(in) :: r0
     real(dp), intent(out) :: Chern


     integer :: ik, ik1, ik2, nkmesh2, ierr
     real(dp) :: theta, phi, r_para, dtheta, dphi, Chern_mpi
     real(dp) :: st, ct, sp, cp, O_x, O_y, O_z

     !> k points in the sphere
     real(dp) :: k_cart(3), k_direct(3)
     real(dp), allocatable :: kpoints(:, :)
     real(dp), allocatable :: thetas(:)
     real(dp), allocatable :: phis(:)

     !> Berry curvature  (3, bands, k)
     complex(dp), allocatable :: Omega_x(:)
     complex(dp), allocatable :: Omega_y(:)
     complex(dp), allocatable :: Omega_z(:)

     nkmesh2= Nk1*Nk2
     allocate(Omega_x(Num_wann), Omega_y(Num_wann), Omega_z(Num_wann))
     allocate(thetas(nkmesh2), phis(nkmesh2), kpoints(3, nkmesh2))

     !> set the sphere defined in the Cartesian coordinates
     !> the first dimension should be in one primitive cell, [0, 2*pi]
     !> the first dimension is the integration direction
     !> the WCCs are calculated along the second k line
     ik= 0
     do ik2=1, Nk2 ! the angle respect to z axis
        theta= (ik2- 1d0)/(Nk2- 1d0)* pi
        if (ik2== 1) theta= (ik2- 1d0+ 0.10)/(Nk2- 1d0)* pi  ! avoid the North pole
        if (ik2== Nk2) theta= (ik2- 1d0- 0.10)/(Nk2- 1d0)* pi  ! avoid the south pole
        do ik1=1, Nk1  ! polar angle in kx-ky plane
           ik = ik+ 1
           phi= (ik1- 1d0)/Nk1* 2d0* pi
           r_para= r0* sin(theta)
           k_cart(1)= k0(1)+ r_para* cos(phi)
           k_cart(2)= k0(2)+ r_para* sin(phi)
           k_cart(3)= k0(3)+ r0* cos(theta)
           call cart_direct_rec(k_cart, k_direct)
           kpoints(:, ik)= k_direct
           thetas(ik)= theta
           phis(ik)= phi
        enddo
      enddo
      dtheta= 1.0d0/(Nk2-1d0)*pi
      dphi  = 1.0d0/Nk1*2d0*pi

      Chern_mpi = 0d0
      do ik=1+cpuid, Nk1*Nk2, num_cpu
         k_direct= kpoints(:, ik)
         theta= thetas(ik)
         phi= phis(ik)
         st=sin(theta)
         ct=cos(theta)
         sp=sin(phi)
         cp=cos(phi)
        !call berry_curvarture_singlek_numoccupied_old(k_direct, Omega_x, Omega_y, Omega_z)
         call berry_curvarture_singlek_numoccupied    (k_direct, Omega_x, Omega_y, Omega_z)
         O_x= real(sum(Omega_x(1:Numoccupied)))
         O_y= real(sum(Omega_y(1:Numoccupied)))
         O_z= real(sum(Omega_z(1:Numoccupied)))
         O_x= real((Omega_x(Numoccupied)))
         O_y= real((Omega_y(Numoccupied)))
         O_z= real((Omega_z(Numoccupied)))
         Chern_mpi= Chern_mpi+ st*st*cp*O_x+ st*st*sp*O_y+ st*ct*O_z
      enddo

#if defined (MPI)
     call mpi_allreduce(Chern_mpi, Chern, 1, &
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
     Chern= Chern_mpi 
#endif


      Chern = Chern* r0* r0* dtheta* dphi/pi/2d0

     return
  end subroutine Chern_sphere_single

  subroutine Chern_sphere
     use wmpi
     use para
     implicit none

     integer :: i
     real(dp) :: k0(3), Chern
     real(dp), allocatable :: Chern_array(:)

     allocate(Chern_array(Num_Weyls))
     Chern_array= 0d0

     do i=1, Num_Weyls
        k0= weyl_position_cart(:, i)
        call Chern_sphere_single(k0, kr0, Chern)
        Chern_array(i)= Chern
     enddo

     if (cpuid==0)then
        write(stdout, *)'# Chern number for the Weyl points'
        write(stdout, '("#",a8,2a9, a16)')'kx', 'ky', 'kz', 'Chern'
        do i=1, Num_Weyls
           write(stdout, '(3f9.5, f16.8)')weyl_position_cart(:, i), Chern_array(i)
        enddo
     endif

     return
  end subroutine Chern_sphere


