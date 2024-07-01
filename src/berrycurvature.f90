  subroutine Berry_curvature_singlek_EF(k, mu, Omega_x, Omega_y, Omega_z)
     !> Calculate Berry curvature for a sigle k point
     !> The Fermi distribution is determined by the Fermi level iso_energy
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

     integer :: m, n, i
     real(dp), allocatable :: W(:)
     real(dp) :: beta_fake
     complex(dp), allocatable :: Amat(:, :), DHDk(:, :, :), DHDkdag(:, :, :)
     complex(dp), allocatable :: UU(:, :), UU_dag(:, :), Hamk_bulk(:, :)
     complex(dp), allocatable :: velocity_wann(:, :, :), velocity_Ham(:, :, :)

     !> Fermi-Dirac distribution
     real(dp), external :: fermi

     allocate(W(Num_wann))
     allocate(velocity_wann(Num_wann, Num_wann, 3), velocity_Ham(Num_wann, Num_wann, 3))
     allocate(UU(Num_wann, Num_wann), UU_dag(Num_wann, Num_wann), Hamk_bulk(Num_wann, Num_wann))
     allocate(Amat(Num_wann, Num_wann), DHDk(Num_wann, Num_wann, 3), DHDkdag(Num_wann, Num_wann, 3))
     W=0d0; velocity_wann= 0d0; UU= 0d0; UU_dag= 0d0; velocity_Ham= 0d0

     ! calculation bulk hamiltonian by a direct Fourier transformation of HmnR
     call ham_bulk_atomicgauge(k, Hamk_bulk)

     !> diagonalization by call zheev in lapack
     UU=Hamk_bulk
     call eigensystem_c( 'V', 'U', Num_wann, UU, W)
    !call zhpevx_pack(hamk_bulk,Num_wann, W, UU)

     UU_dag= conjg(transpose(UU))
     call dHdk_atomicgauge(k, velocity_wann)

     !> unitility rotate velocity
     do i=1, 3
        call mat_mul(Num_wann, velocity_wann(:, :, i), UU, Amat) 
        call mat_mul(Num_wann, UU_dag, Amat, velocity_Ham(:, :, i)) 
     enddo

     Omega_x=0d0;Omega_y=0d0; Omega_z=0d0
     do m= 1, Num_wann
        do n= 1, Num_wann
           if (abs(W(m)-W(n))<eps6) cycle
           Omega_x(m)= Omega_x(m)+ velocity_Ham(n, m, 2)*velocity_Ham(m, n, 3)/((W(m)-W(n))**2)
           Omega_y(m)= Omega_y(m)+ velocity_Ham(n, m, 3)*velocity_Ham(m, n, 1)/((W(m)-W(n))**2)
           Omega_z(m)= Omega_z(m)+ velocity_Ham(n, m, 1)*velocity_Ham(m, n, 2)/((W(m)-W(n))**2)
        enddo ! m
     enddo ! n

     Omega_x= -Omega_x*2d0*zi
     Omega_y= -Omega_y*2d0*zi
     Omega_z= -Omega_z*2d0*zi

     !> consider the Fermi-distribution according to the broadening Earc_eta
     if (Fermi_broadening<eps6) Fermi_broadening= eps6
     beta_fake= 1d0/Fermi_broadening
     do m= 1, Num_wann
        Omega_x(m)= Omega_x(m)*fermi(W(m)-mu, beta_fake)
        Omega_y(m)= Omega_y(m)*fermi(W(m)-mu, beta_fake)
        Omega_z(m)= Omega_z(m)*fermi(W(m)-mu, Beta_fake)
     enddo

     return
  end subroutine Berry_curvature_singlek_EF

  subroutine Berry_curvature_singlek_numoccupied_slab_total(k, Omega_z)
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

     integer :: m, n, Mdim, i1, i2
     real(dp), allocatable :: W(:)
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

     Omega_z= -aimag(Omega_z*2d0)/Angstrom2atomic**2

     return
  end subroutine Berry_curvature_singlek_numoccupied_slab_total

  subroutine Berry_curvature_singlek_numoccupied_total(k, Omega_x, Omega_y, Omega_z)
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

     integer :: m, n, i
     real(dp), allocatable :: W(:)
     complex(dp), allocatable :: Amat(:, :)
     complex(dp), allocatable :: UU(:, :), UU_dag(:, :), Hamk_bulk(:, :)
     complex(dp), allocatable :: velocity_wann(:, :, :), velocity_Ham(:, :, :)

     allocate(W(Num_wann))
     allocate(velocity_wann(Num_wann, Num_wann, 3), velocity_Ham(Num_wann, Num_wann, 3))
     allocate(UU(Num_wann, Num_wann), UU_dag(Num_wann, Num_wann), Hamk_bulk(Num_wann, Num_wann))
     allocate(Amat(Num_wann, Num_wann))
     W=0d0; velocity_wann= 0d0; UU= 0d0; UU_dag= 0d0; velocity_Ham= 0d0

     ! calculation bulk hamiltonian by a direct Fourier transformation of HmnR
     !call ham_bulk_latticegauge(k, Hamk_bulk)
     call ham_bulk_atomicgauge    (k, Hamk_bulk)

     !> diagonalization by call zheev in lapack
     UU=Hamk_bulk
     call eigensystem_c( 'V', 'U', Num_wann, UU, W)
    !call zhpevx_pack(hamk_bulk,Num_wann, W, UU)

     UU_dag= conjg(transpose(UU))

     call dHdk_atomicgauge(k, velocity_wann)

     !> unitility rotate velocity
     do i=1, 3
        call mat_mul(Num_wann, velocity_wann(:, :, i), UU, Amat) 
        call mat_mul(Num_wann, UU_dag, Amat, velocity_Ham(:, :, i)) 
     enddo

     Omega_x=0d0; Omega_y=0d0; Omega_z=0d0
     do m= 1, NumOccupied !> sum over valence band
        do n= NumOccupied+1, Num_wann !> sum over conduction band
           Omega_x= Omega_x+ velocity_Ham(n, m, 2)*velocity_Ham(m, n, 3)/((W(m)-W(n))**2)
           Omega_y= Omega_y+ velocity_Ham(n, m, 3)*velocity_Ham(m, n, 1)/((W(m)-W(n))**2)
           Omega_z= Omega_z+ velocity_Ham(n, m, 1)*velocity_Ham(m, n, 2)/((W(m)-W(n))**2)
        enddo ! m
     enddo ! n

     Omega_x=  aimag(Omega_x*2d0)  
     Omega_y=  aimag(Omega_y*2d0)  
     Omega_z=  aimag(Omega_z*2d0)  

     deallocate(W, UU, UU_dag, hamk_bulk)
     deallocate(Amat, velocity_wann, velocity_Ham)

     return
  end subroutine Berry_curvature_singlek_numoccupied_total

  subroutine Berry_curvature_singlek_numoccupied(k, Omega_x, Omega_y, Omega_z)
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

     integer :: m, n, i
     real(dp), allocatable :: W(:)
     complex(dp), allocatable :: Amat(:, :), DHDk(:, :, :), DHDkdag(:, :, :)
     complex(dp), allocatable :: UU(:, :), UU_dag(:, :), Hamk_bulk(:, :)
     complex(dp), allocatable :: velocity_wann(:, :, :), velocity_Ham(:, :, :)

     allocate(W(Num_wann))
     allocate(velocity_wann(Num_wann, Num_wann, 3), velocity_Ham(Num_wann, Num_wann, 3))
     allocate(UU(Num_wann, Num_wann), UU_dag(Num_wann, Num_wann), Hamk_bulk(Num_wann, Num_wann))
     allocate(Amat(Num_wann, Num_wann), DHDk(Num_wann, Num_wann, 3), DHDkdag(Num_wann, Num_wann, 3))
     W=0d0; velocity_wann= 0d0; UU= 0d0; UU_dag= 0d0; velocity_Ham= 0d0

     ! calculation bulk hamiltonian by a direct Fourier transformation of HmnR
     call ham_bulk_atomicgauge    (k, Hamk_bulk)

     !> diagonalization by call zheev in lapack
     UU=Hamk_bulk
     call eigensystem_c( 'V', 'U', Num_wann, UU, W)
    !call zhpevx_pack(hamk_bulk,Num_wann, W, UU)

     UU_dag= conjg(transpose(UU))

     call dHdk_atomicgauge(k, velocity_wann)

     !> unitility rotate velocity
     do i=1, 3
        call mat_mul(Num_wann, velocity_wann(:, :, i), UU, Amat) 
        call mat_mul(Num_wann, UU_dag, Amat, velocity_Ham(:, :, i)) 
     enddo

     Omega_x=0d0;Omega_y=0d0; Omega_z=0d0
    !do m= NumOccupied  , NumOccupied
    !   do n= NumOccupied+1, NumOccupied+1
     do m= 1, NumOccupied
        do n= 1, Num_wann
           if (m==n) cycle
           Omega_x(m)= Omega_x(m)+ velocity_Ham(n, m, 2)*velocity_Ham(m, n, 3)/((W(m)-W(n))**2)
           Omega_y(m)= Omega_y(m)+ velocity_Ham(n, m, 3)*velocity_Ham(m, n, 1)/((W(m)-W(n))**2)
           Omega_z(m)= Omega_z(m)+ velocity_Ham(n, m, 1)*velocity_Ham(m, n, 2)/((W(m)-W(n))**2)
        enddo ! m
     enddo ! n

     Omega_x= -Omega_x*2d0*zi
     Omega_y= -Omega_y*2d0*zi
     Omega_z= -Omega_z*2d0*zi

     return
  end subroutine Berry_curvature_singlek_numoccupied

  subroutine Berry_curvature_singlek_allbands(Dmn_Ham, Omega_BerryCurv)
     !> Calculate Berry curvature for a sigle k point and all bands
     !> ref : eqn (30) Physical Review B 74, 195118(2006)
     !> \Omega_n^{\gamma}(k)=i\sum_{\alpha\beta}\epsilon_{\gamma\alpha\beta}(D^{\alpha\dag}D^{\beta})_{nn}
     !> Dec. 30 2019 by Quansheng Wu
     !> Copyright (c) 2019 QuanSheng Wu. All rights reserved.

     use wmpi
     use para
     implicit none

     !> D_mn^H=V_mn/(En-Em) for m!=n
     !> D_nn^H=0 
     complex(dp), intent(in) :: Dmn_Ham(Num_wann, Num_wann, 3)

     !> Berry curvature vectors for all bands
     real(dp), intent(out) :: Omega_BerryCurv(Num_wann, 3)

     integer :: m, n, i, alphai(3), betai(3)
     complex(dp), allocatable :: Amat(:, :)
     complex(dp), allocatable :: Dmn_Ham_dag(:, :, :)
     complex(dp) :: zz(3)

     allocate(Amat(Num_wann, Num_wann))
     allocate(Dmn_Ham_dag(Num_wann, Num_wann, 3))
     Amat= 0d0
     do i=1, 3
        Dmn_Ham_dag(:, :, i)=conjg(transpose(Dmn_Ham(:, :, i)))
     !  Dmn_Ham_dag(:, :, i)=-(Dmn_Ham(:, :, i))
     enddo
    !print *, Dmn_Ham_dag(1, 2, 2)
    !print *, Dmn_Ham(2, 1, 3)

     Omega_BerryCurv= 0d0
     alphai(1)=2;alphai(2)=3;alphai(3)=1;
     betai(1) =3; betai(2)=1; betai(3)=2;

     zz=0d0
     do i=1, 3
        call mat_mul(Num_wann, Dmn_Ham_dag(:, :, alphai(i)), Dmn_Ham(:, :, betai(i)), Amat)
        do m=1, Num_wann
          !Omega_BerryCurv(m, i)= Omega_BerryCurv(m, i)- aimag(Amat(m, m))
           do n=1, Num_wann
              Omega_BerryCurv(m, i)= Omega_BerryCurv(m, i)- &
                 aimag(Dmn_Ham_dag(m, n, alphai(i))*Dmn_Ham(n, m, betai(i)))
           enddo
        enddo
        call mat_mul(Num_wann, Dmn_Ham_dag(:, :, betai(i)), Dmn_Ham(:, :, alphai(i)), Amat)
        do m=1, Num_wann
          !Omega_BerryCurv(m, i)= Omega_BerryCurv(m, i)+ aimag(Amat(m, m))
           do n=1, Num_wann
              Omega_BerryCurv(m, i)= Omega_BerryCurv(m, i)+ &
                 aimag(Dmn_Ham_dag(m, n, betai(i))*Dmn_Ham(n, m, alphai(i)))
           enddo
        enddo
     enddo

     return
  end subroutine Berry_curvature_singlek_allbands


  subroutine orbital_magnetization_singlek_allbands(Dmn_Ham, Vmn_Ham_nondiag, m_OrbMag)
     !> Calculate orbital magnetization for a sigle k point and all bands
     !> m_n^{\gamma}(k)=-i*e/2/\hbar\sum_{\alpha\beta}\epsilon_{\gamma\alpha\beta}(D^{\alpha\dag}V^{\beta})_{nn}
     !> Dec. 30 2019 by Quansheng Wu
     !> Copyright (c) 2019 QuanSheng Wu. All rights reserved.

     use wmpi
     use para
     implicit none

     !> D_mn^H=V_mn/(En-Em) for m!=n
     !> D_nn^H=0 
     complex(dp), intent(in) :: Dmn_Ham(Num_wann, Num_wann, 3)

     !> Vmn_Ham_nondiag(m, n, i)= Vmn_Ham(m, n, i) for m!=n
     !> Vmn_Ham_nondiag(n, n, i)= 0
     complex(dp), intent(in) :: Vmn_Ham_nondiag(Num_wann, Num_wann, 3)

     !> Berry curvature vectors for all bands
     real(dp), intent(out) :: m_OrbMag(Num_wann, 3)

     integer :: m, i, ialpha(3), ibeta(3)
     complex(dp), allocatable :: Amat(:, :)
     complex(dp), allocatable :: Dmn_Ham_dag(:, :, :)

     allocate(Amat(Num_wann, Num_wann))
     allocate(Dmn_Ham_dag(Num_wann, Num_wann, 3))
     Amat=0d0
     Dmn_Ham_dag = -Dmn_Ham

     m_OrbMag= 0d0
     ialpha(1)=2;ialpha(2)=3;ialpha(3)=1;
     ibeta(1)=3;ibeta(2)=1;ibeta(3)=2;

     do i=1, 3
        call mat_mul(Num_wann, Dmn_Ham_dag(:, :, ialpha(i)), Vmn_Ham_nondiag(:, :, ibeta(i)), Amat)
        do m=1, Num_wann
           m_OrbMag(m, i)= m_OrbMag(m, i)+ aimag(Amat(m, m))
        enddo
        call mat_mul(Num_wann, Dmn_Ham_dag(:, :, ibeta(i)), Vmn_Ham_nondiag(:, :, ialpha(i)), Amat)
        do m=1, Num_wann
           m_OrbMag(m, i)= m_OrbMag(m, i)- aimag(Amat(m, m))
        enddo
     enddo

     !> 1/2d0 comes from e/\hbar/2d0, here we use atomic unit
     m_OrbMag= m_OrbMag/2d0

     !> if we want to use Bohr magneton \mu_B as unit, we need multiply it by two.

     return
  end subroutine orbital_magnetization_singlek_allbands

  subroutine get_Vmn_Ham_nondiag(velocity_Ham, velocity_Ham_nondiag)
     !> the diagonal of velocity_Ham_nondiag is zero
     !> and velocity_Ham_nondiag(m, n)= -velocity_Ham(m, n)
     use para, only : dp, Num_wann
     implicit none

     complex(dp), intent(in) :: velocity_Ham(Num_wann, Num_wann, 3)
     complex(dp), intent(out) :: velocity_Ham_nondiag(Num_wann, Num_wann, 3)
     integer :: n

     velocity_Ham_nondiag= -velocity_Ham
     do n=1, Num_wann
        velocity_Ham_nondiag(n, n, :)=0d0
     enddo

     return
  end subroutine get_Vmn_Ham_nondiag

  subroutine get_Dmn_Ham(W, velocity_Ham, Dmn_Ham)
     !> define D matrix 
     !> Eq (24) Phys.Rev.B 74, 195118(2006)
     !> Dmn_Ham_dag=-Dmn_Ham
     use para, only : dp, Num_wann, eps12, zzero
     implicit none

     real(dp), intent(in) :: W(Num_wann)
     complex(dp), intent(in) :: velocity_Ham(Num_wann, Num_wann, 3)
     complex(dp), intent(out) :: Dmn_Ham(Num_wann, Num_wann, 3)

     integer :: m, n, i

     Dmn_Ham = zzero 
     do i=1, 3
        do n=1, Num_wann
           do m=1, Num_wann
              if (abs(W(n)-W(m))>eps12) then
                 Dmn_Ham(m,n,i)=velocity_Ham(m,n,i)/(W(n)-W(m))
              endif
           enddo
           Dmn_Ham(n, n, i)=zzero
        enddo
     enddo

     return
  end subroutine get_Dmn_Ham

  subroutine Berry_curvature_slab
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
    
     integer :: ik, ierr, Mdim, i, j

     real(dp) :: k(2)  

     !> k points slice
     real(dp), allocatable :: k12(:, :)
     real(dp), allocatable :: k12_xyz(:, :)
   
     real(dp), external :: norm

     !> Berry curvature  (k)
     complex(dp), allocatable :: Omega_z(:)
     complex(dp), allocatable :: Omega(:)
     complex(dp), allocatable :: Omega_mpi(:)

     Mdim = Num_wann*Nslab

     allocate( k12(2, Nk1*Nk2))
     allocate( k12_xyz(2, Nk1*Nk2))
     allocate( Omega_z(Mdim))
     allocate( Omega    (Nk1*Nk2))
     allocate( Omega_mpi(Nk1*Nk2))
     k12=0d0
     k12_xyz=0d0
     omega= 0d0
     omega_mpi= 0d0
    
     !> k12 is centered at K2d_start
     ik=0
     do i= 1, nk1
        do j= 1, nk2
           ik=ik+1
           k12(:, ik)=K2D_start+ (i-1)*K2D_vec1/dble(nk1-1) &
                      + (j-1)*K2D_vec2/dble(nk2-1)- (K2D_vec1+K2D_vec2)/2d0
           k12_xyz(:, ik)= k12(1, ik)* Ka2+ k12(2, ik)* Kb2
        enddo
     enddo

     do ik= 1+ cpuid, Nk1*Nk2, num_cpu
        if (cpuid==0) write(stdout, *)'Berry curvature ik, nk1*nk2 ', ik, Nk1*Nk2

        !> diagonalize hamiltonian
        k= k12(:, ik)

        Omega_z= 0d0

        call Berry_curvature_singlek_numoccupied_slab_total(k, Omega_z(1))
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
        write(outfileindex, '(20a28)')'# Unit of Berry curvature is Angstrom^2'
        write(outfileindex, '(20a28)')'# kx (1/A)', 'ky (1/A)', &
           'Omega_z'
        ik= 0
        do i= 1, nk1
           do j= 1, nk2
              ik= ik+ 1
              write(outfileindex, '(20E28.10)')k12_xyz(:, ik)*Angstrom2atomic, real(Omega_mpi(ik))/Angstrom2atomic**2
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
        write(outfileindex, '(a)')"set title 'Berry Curvature {/Symbol W}_z (Ang^2)'"
        write(outfileindex, '(a)')"set colorbox"
        write(outfileindex, '(a)')"splot 'Berrycurvature_slab.dat' u 1:2:3 w pm3d"
        close(outfileindex)
     endif


#if defined (MPI)
     call mpi_barrier(mpi_cmw, ierr)
#endif

     deallocate( k12)
     deallocate( k12_xyz)
     deallocate( Omega_z)
     deallocate( Omega, Omega_mpi)
 
     return

  end subroutine Berry_curvature_slab

  subroutine Berry_curvature_line_occupied
     !> Calculate Berry curvature for a k line defined bu KPATH_BULK
     !
     !> ref : Physical Review B 74, 195118(2006)
     !
     !> eqn (9), Eqn (34)
     !
     !> Sep. 28 2018 by Quansheng Wu @ EPFL
     !
     ! Copyright (c) 2018 QuanSheng Wu. All rights reserved.

     use wmpi
     use para
     implicit none
    
     integer :: ik, ierr, i

     real(dp) :: k(3), o1(3), k_cart(3)
     real(dp) :: time_start, time_end, time_start0, ybound_min, ybound_max

     real(dp), external :: norm

     complex(dp), allocatable :: Omega_x(:), Omega_y(:), Omega_z(:)
     !> Berry curvature  (3, k)
     real(dp), allocatable :: Omega(:, :), Omega_mpi(:, :)

     !> Berry curvature  (3, bands, k)
     real(dp), allocatable :: Omega_sep_bk(:, :), Omega_sep_bk_mpi(:, :)

     !> energy bands
     real(dp), allocatable :: eigv(:,:)
     real(dp), allocatable :: eigv_mpi(:,:)

     allocate( Omega_x(Num_wann))
     allocate( Omega_y(Num_wann))
     allocate( Omega_z(Num_wann))
     allocate( eigv    (Num_wann, nk3_band))
     allocate( eigv_mpi(Num_wann, nk3_band))
     allocate( Omega    (3, nk3_band))
     allocate( Omega_mpi(3, nk3_band))
     omega= 0d0
     omega_mpi= 0d0
    
     time_start= 0d0
     time_start0= 0d0
     call now(time_start0)
     time_start= time_start0
     time_end  = time_start0
     do ik= 1+ cpuid, nk3_band, num_cpu
        if (cpuid==0.and. mod(ik/num_cpu, 100)==0) &
           write(stdout, '(a, i9, "  /", i10, a, f10.1, "s", a, f10.1, "s")') &
           ' Berry curvature: ik', ik, nk3_band, ' time left', &
           (nk3_band-ik)*(time_end- time_start)/num_cpu, &
           ' time elapsed: ', time_end-time_start0 

        !> a k point in fractional coordinates
        k= kpath_3d(:, ik)

        call now(time_start)
 
        Omega_x= 0d0
        Omega_y= 0d0
        Omega_z= 0d0

        !call Berry_curvature_singlek_numoccupied_old(k, Omega_x, Omega_y, Omega_z)
        if (Berrycurvature_kpath_EF_calc) then
           call Berry_curvature_singlek_EF(k, iso_energy, Omega_x, Omega_y, Omega_z)
        else if (BerryCurvature_kpath_Occupied_calc) then
           call Berry_curvature_singlek_numoccupied_total(k, Omega_x(1), Omega_y(1), Omega_z(1))
        else
           write(*, *) 'ERROR: In subroutine Berry_curvature_line, we only support BerryCurvature_kpath_Occupied_calc '
           write(*, *) ' and Berrycurvature_kpath_EF_calc '
           stop
        endif
 
        Omega(1, ik) = real(sum(Omega_x))
        Omega(2, ik) = real(sum(Omega_y))
        Omega(3, ik) = real(sum(Omega_z))
        call now(time_end)
     enddo ! ik

     Omega_mpi= 0d0

#if defined (MPI)
     call mpi_allreduce(Omega,Omega_mpi,size(Omega_mpi),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
     Omega_mpi= Omega
#endif

     !> output the Berry curvature to file
     outfileindex= outfileindex+ 1
     if (cpuid==0) then
        open(unit=outfileindex, file='Berrycurvature_line.dat')
        write(outfileindex, '(20a18)')'# Column 1 kpath with accumulated length in the kpath'
        write(outfileindex, '(20a18)')'# Column 2-4 Berry curvature (Ang^2)'
        write(outfileindex, '(20a18)')'# k (1/A)', &
           'real(Omega_x)', 'real(Omega_y)', 'real(Omega_z)'

        do ik= 1, nk3_band
           k=kpath_3d(:, ik)
           write(outfileindex, '(20E18.8)')k3len(ik)*Angstrom2atomic, real(Omega_mpi(:, ik))/Angstrom2atomic**2
        enddo

        close(outfileindex)
     endif

     ybound_min= minval(real(Omega_mpi))-2
     ybound_max= maxval(real(Omega_mpi))+5 
     !> write script for gnuplot
     outfileindex= outfileindex+ 1
     if (cpuid==0) then
        open(unit=outfileindex, file='Berrycurvature_line.gnu')
        write(outfileindex, '(a)')"set encoding iso_8859_1"
        write(outfileindex, '(a)') 'set terminal pdf enhanced color font ",16"'
        write(outfileindex, '(a)')"set output 'Berrycurvature_line.pdf'"
        write(outfileindex, '(a)')'set style data linespoints'
        write(outfileindex, '(a)')'unset key'
        write(outfileindex, '(a)')'set pointsize 0.8'
        write(outfileindex, '(a)')'set view 0,0'
        write(outfileindex, '(a, f10.5, a)')'set xrange [0: ', maxval(k3len)*Angstrom2atomic, ']'
        if (index(Particle,'phonon')/=0) then
           write(outfileindex, '(a, f10.5, a)')'set yrange [0:', ybound_max, ']'
           write(outfileindex, '(a)')'set ylabel "Frequency (THz)"'
        else
           write(outfileindex, '(a)')'set ylabel "Energy (eV)"'
           write(outfileindex, '(a, f10.5, a, f10.5, a)')'set yrange [', ybound_min, ':', ybound_max, ']'
        endif
        write(outfileindex, 202, advance="no") (k3line_name(i), k3line_stop(i)*Angstrom2atomic, i=1, nk3lines)
        write(outfileindex, 203)k3line_name(nk3lines+1), k3line_stop(nk3lines+1)*Angstrom2atomic
  
        do i=1, nk3lines-1
           if (index(Particle,'phonon')/=0) then
              write(outfileindex, 204)k3line_stop(i+1)*Angstrom2atomic, 0.0, k3line_stop(i+1)*Angstrom2atomic, ybound_max
           else
              write(outfileindex, 204)k3line_stop(i+1)*Angstrom2atomic, ybound_min, k3line_stop(i+1)*Angstrom2atomic, ybound_max
           endif
        enddo
        write(outfileindex, '(a)')"plot 'Berrycurvature_line.dat' \"  
        write(outfileindex, '(a)')"u 1:2 w lp lc rgb 'red'   lw 2 pt 7 ps 0.2 title '{/Symbol W}_x', \" 
        write(outfileindex, '(a)')"'' u 1:3 w lp lc rgb 'green' lw 2 pt 7 ps 0.2 title '{/Symbol W}_y', \" 
        write(outfileindex, '(a)')"'' u 1:4 w lp lc rgb 'blue'  lw 2 pt 7 ps 0.2 title '{/Symbol W}_z' "
        close(outfileindex)
     endif

202 format('set xtics (',20('"',A3,'" ',F10.5,','))
203 format(A3,'" ',F10.5,')')
204 format('set arrow from ',F10.5,',',F10.5,' to ',F10.5,',',F10.5, ' nohead')


#if defined (MPI)
     call mpi_barrier(mpi_cmw, ierr)
#endif

     deallocate( Omega_x, Omega_y, Omega_z)
     deallocate( Omega, Omega_mpi)
 
     return

  end subroutine Berry_curvature_line_occupied


  subroutine berry_curvature_cube
     !> Calculate Berry curvature  in a cube defined in KCUBE_BULK
     !
     !> ref : Physical Review B 74, 195118(2006)
     !
     !> eqn (34)
     !
     !> Sep. 28 2018 by Quansheng Wu @ EPFL
     !
     ! Copyright (c) 2018 QuanSheng Wu. All rights reserved.

     use wmpi
     use para
     implicit none
    
     integer :: ik, ierr, ikx, iky, ikz, n_kpoints, i, m, n

     real(dp) :: k(3), o1(3), k_cart(3), emin, emax
     real(dp) :: time_start, time_end, time_start0

     real(dp), external :: norm

     !> Berry curvature and orbital magnetization (3, bands, k)
     real(dp), allocatable :: Omega_allk(:, :, :), Omega_allk_mpi(:, :, :)
     real(dp), allocatable :: m_OrbMag_allk(:, :, :), m_OrbMag_allk_mpi(:, :, :)
     real(dp), allocatable :: Omega_BerryCurv(:, :), m_OrbMag(:, :)

     !> velocity matrix in Wannier basis 
     complex(dp), allocatable :: Vmn_wann(:, :, :)

     !> velocity matrix in Hamiltonian basis 
     complex(dp), allocatable :: Vmn_Ham(:, :, :), Vmn_Ham_nondiag(:, :, :)

     complex(dp), allocatable :: Dmn_Ham(:, :, :)

     !> Hamiltonian, eigenvalue and eigenvectors
     complex(dp), allocatable :: UU(:, :)
     real(dp), allocatable :: W(:)
     real(dp), allocatable :: eigval_allk(:, :), eigval_allk_mpi(:, :)

     if (BerryCurvature_Cube_calc) then
        n_kpoints=Nk1*Nk2*Nk3
     else
        n_kpoints= nk3_band 
     endif

     allocate(Vmn_wann(Num_wann, Num_wann, 3), Vmn_Ham(Num_wann, Num_wann, 3))
     allocate(Dmn_Ham(Num_wann,Num_wann,3), Vmn_Ham_nondiag(Num_wann, Num_wann, 3))
     allocate(W(Num_wann))
     allocate(UU(Num_wann, Num_wann))
     allocate(eigval_allk(Num_wann, n_kpoints))
     allocate(eigval_allk_mpi(Num_wann, n_kpoints))
     allocate( Omega_BerryCurv(Num_wann, 3), m_OrbMag(Num_wann, 3))
     allocate( Omega_allk    (Num_wann, 3, n_kpoints))
     allocate( Omega_allk_mpi(Num_wann, 3, n_kpoints))
     allocate( m_OrbMag_allk    (Num_wann, 3, n_kpoints))
     allocate( m_OrbMag_allk_mpi(Num_wann, 3, n_kpoints))
     Omega_BerryCurv= 0d0
     m_OrbMag=0d0
     
     Omega_allk= 0d0
     eigval_allk= 0d0
     m_OrbMag_allk=0d0
     Omega_allk_mpi= 0d0
     eigval_allk_mpi= 0d0
     m_OrbMag_allk_mpi=0d0
    
     time_start= 0d0
     time_start0= 0d0
     call now(time_start0)
     time_start= time_start0
     time_end  = time_start0
     do ik= 1+ cpuid, n_kpoints, num_cpu
        if (cpuid==0.and. mod(ik/num_cpu, 100)==0) &
           write(stdout, '(a, i9, "  /", i10, a, f10.1, "s", a, f10.1, "s")') &
           ' Berry curvature: ik', ik, n_kpoints, ' time left', &
           (n_kpoints-ik)*(time_end- time_start)/num_cpu, &
           ' time elapsed: ', time_end-time_start0 

        !> if we calculate BC in the BZ, we generate kpoints in the BZ
        if (BerryCurvature_Cube_calc) then
           !> kbulk mode
           ikx= (ik-1)/(nk2*nk3)+1
           iky= ((ik-1-(ikx-1)*Nk2*Nk3)/nk3)+1
           ikz= (ik-(iky-1)*Nk3- (ikx-1)*Nk2*Nk3)
           k= K3D_start_cube+ K3D_vec1_cube*(ikx-1)/dble(nk1)  &
              + K3D_vec2_cube*(iky-1)/dble(nk2)  &
              + K3D_vec3_cube*(ikz-1)/dble(nk3) 
 
        elseif (BerryCurvature_kpath_sepband_calc) then
           !> kpath mode
           k= kpath_3d(:, ik)
        endif

        call now(time_start)
        
        !> calculation bulk hamiltonian by a direct Fourier transformation of HmnR
        call ham_bulk_atomicgauge(k, UU)
   
        !> diagonalization by call zheev in lapack
        call eigensystem_c( 'V', 'U', Num_wann, UU, W)
        eigval_allk(:, ik) = W

        !> get velocity operator in Hamiltonian basis
        call dHdk_atomicgauge_Ham(k, UU, Vmn_Ham)

        call get_Dmn_Ham(W, Vmn_Ham, Dmn_Ham)
        call get_Vmn_Ham_nondiag(Vmn_Ham, Vmn_Ham_nondiag)

        call Berry_curvature_singlek_allbands(Dmn_Ham, Omega_BerryCurv)
        call orbital_magnetization_singlek_allbands(Dmn_Ham, Vmn_Ham_nondiag, m_OrbMag)
        Omega_allk(:, :, ik) = Omega_BerryCurv
        m_OrbMag_allk(:, :, ik) = m_OrbMag

        call now(time_end)
     enddo ! ik

#if defined (MPI)
     call mpi_allreduce(Omega_allk,Omega_allk_mpi,size(Omega_allk_mpi),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
     call mpi_allreduce(eigval_allk,eigval_allk_mpi,size(eigval_allk_mpi),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
     call mpi_allreduce(m_OrbMag_allk,m_OrbMag_allk_mpi,size(m_OrbMag_allk_mpi),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
     Omega_allk_mpi= Omega_allk
     eigval_allk_mpi= eigval_allk
     m_OrbMag_allk_mpi= m_OrbMag_allk
#endif

     IF (BerryCurvature_Cube_calc) THEN
     !> write out Berry curvature and orbital magnetization to a file which
     !> can be open by software Fermisurfer.
     outfileindex= outfileindex+ 1
     if (cpuid==0) then
        open(unit=outfileindex, file='Berrycurvature_norm.frmsf')
        write(outfileindex, *) Nk1, Nk2, Nk3
        write(outfileindex, *) 1
        write(outfileindex, *) Num_wann
        write(outfileindex, '(3f12.6)') Origin_cell%Kua
        write(outfileindex, '(3f12.6)') Origin_cell%Kub
        write(outfileindex, '(3f12.6)') Origin_cell%Kuc
        do m=1, Num_wann
           do ik= 1, n_kpoints
              write(outfileindex, '(E18.10)') eigval_allk_mpi(m, ik)-iso_energy
           enddo
        enddo
        do m=1, Num_wann
           do ik= 1, n_kpoints
              o1= Omega_allk_mpi(m, :, ik)/Angstrom2atomic**2
              write(outfileindex, '(E18.10)') norm(o1)
           enddo
        enddo
        close(outfileindex)
     endif


     !> write out Berry curvature and orbital magnetization to a file which
     !> can be open by software Fermisurfer.
     outfileindex= outfileindex+ 1
     if (cpuid==0) then
        open(unit=outfileindex, file='Berrycurvature_x.frmsf')
        write(outfileindex, *) Nk1, Nk2, Nk3
        write(outfileindex, *) 1
        write(outfileindex, *) Num_wann
        write(outfileindex, '(3f12.6)') Origin_cell%Kua
        write(outfileindex, '(3f12.6)') Origin_cell%Kub
        write(outfileindex, '(3f12.6)') Origin_cell%Kuc
        do m=1, Num_wann
           do ik= 1, n_kpoints
              write(outfileindex, '(E18.10)') eigval_allk_mpi(m, ik)-iso_energy
           enddo
        enddo
        do m=1, Num_wann
           do ik= 1, n_kpoints
              o1= Omega_allk_mpi(m, :, ik)/Angstrom2atomic**2
              write(outfileindex, '(E18.10)') o1(1)
           enddo
        enddo
        close(outfileindex)
     endif


     !> write out Berry curvature and orbital magnetization to a file which
     !> can be open by software Fermisurfer.
     outfileindex= outfileindex+ 1
     if (cpuid==0) then
        open(unit=outfileindex, file='Berrycurvature_z.frmsf')
        write(outfileindex, *) Nk1, Nk2, Nk3
        write(outfileindex, *) 1
        write(outfileindex, *) Num_wann
        write(outfileindex, '(3f12.6)') Origin_cell%Kua
        write(outfileindex, '(3f12.6)') Origin_cell%Kub
        write(outfileindex, '(3f12.6)') Origin_cell%Kuc
        do m=1, Num_wann
           do ik= 1, n_kpoints
              write(outfileindex, '(E18.10)') eigval_allk_mpi(m, ik)-iso_energy
           enddo
        enddo
        do m=1, Num_wann
           do ik= 1, n_kpoints
              o1= Omega_allk_mpi(m, :, ik)/Angstrom2atomic**2
              write(outfileindex, '(E18.10)') o1(3)
           enddo
        enddo
        close(outfileindex)
     endif

     !> orbital_magnetization
     outfileindex= outfileindex+ 1
     if (cpuid==0) then
        open(unit=outfileindex, file='Orbital_magnetization_norm.frmsf')
        write(outfileindex, *) Nk1, Nk2, Nk3
        write(outfileindex, *) 1
        write(outfileindex, *) Num_wann
        write(outfileindex, '(3f12.6)') Origin_cell%Kua
        write(outfileindex, '(3f12.6)') Origin_cell%Kub
        write(outfileindex, '(3f12.6)') Origin_cell%Kuc
        do m=1, Num_wann
           do ik= 1, n_kpoints
              write(outfileindex, '(E18.10)') eigval_allk_mpi(m, ik)-iso_energy
           enddo
        enddo
        do m=1, Num_wann
           do ik= 1, n_kpoints
              o1= m_OrbMag_allk_mpi(m, :, ik)
              write(outfileindex, '(E18.10)') norm(o1)
           enddo
        enddo
        close(outfileindex)
     endif

     !> orbital_magnetization
     outfileindex= outfileindex+ 1
     if (cpuid==0) then
        open(unit=outfileindex, file='Orbital_magnetization_z.frmsf')
        write(outfileindex, *) Nk1, Nk2, Nk3
        write(outfileindex, *) 1
        write(outfileindex, *) Num_wann
        write(outfileindex, '(3f12.6)') Origin_cell%Kua
        write(outfileindex, '(3f12.6)') Origin_cell%Kub
        write(outfileindex, '(3f12.6)') Origin_cell%Kuc
        do m=1, Num_wann
           do ik= 1, n_kpoints
              write(outfileindex, '(E18.10)') eigval_allk_mpi(m, ik)-iso_energy
           enddo
        enddo
        do m=1, Num_wann
           do ik= 1, n_kpoints
              o1= m_OrbMag_allk_mpi(m, :, ik)
              write(outfileindex, '(E18.10)') o1(3)
           enddo
        enddo
        close(outfileindex)
     endif

     !> orbital_magnetization
     outfileindex= outfileindex+ 1
     if (cpuid==0) then
        open(unit=outfileindex, file='Orbital_magnetization_x.frmsf')
        write(outfileindex, *) Nk1, Nk2, Nk3
        write(outfileindex, *) 1
        write(outfileindex, *) Num_wann
        write(outfileindex, '(3f12.6)') Origin_cell%Kua
        write(outfileindex, '(3f12.6)') Origin_cell%Kub
        write(outfileindex, '(3f12.6)') Origin_cell%Kuc
        do m=1, Num_wann
           do ik= 1, n_kpoints
              write(outfileindex, '(E18.10)') eigval_allk_mpi(m, ik)-iso_energy
           enddo
        enddo
        do m=1, Num_wann
           do ik= 1, n_kpoints
              o1= m_OrbMag_allk_mpi(m, :, ik)
              write(outfileindex, '(E18.10)') o1(1)
           enddo
        enddo
        close(outfileindex)
     endif

     ELSE


     !> output the Berry curvature to file
     outfileindex= outfileindex+ 1
     if (cpuid==0) then
        open(unit=outfileindex, file='Berrycurvature_line_sepband.txt')
        write(outfileindex, '(a)')'# Column 1 kpath with accumulated length in the kpath, Coloum 2: energy'
        write(outfileindex, '(a)')'# Column 2-4 Berry curvature (Ang^2)'
        write(outfileindex, "('#column', i5, 20i16)")(i, i=1, 8)
        write(outfileindex, '(20a16)')'# k (1/A)', " eig", &
           'Omega_x', 'Omega_y', 'Omega_z', &
           'm_x', 'm_y', 'm_z'
        do i=1, Num_wann
           do ik=1, n_kpoints
              write(outfileindex, '(200E16.6)')k3len(ik)*Angstrom2atomic, eigval_allk_mpi(i, ik)/eV2Hartree, &
                 Omega_allk_mpi(i, :, ik), m_OrbMag_allk_mpi(i, :, ik)
           enddo
           write(outfileindex, *)' '
        enddo
        close(outfileindex)
 

        close(outfileindex)
     endif

     !> minimum and maximum value of energy bands
     emin=  minval(eigval_allk_mpi)/eV2Hartree-0.5d0
     emax=  maxval(eigval_allk_mpi)/eV2Hartree+0.5d0

      !> write script for gnuplot
     outfileindex= outfileindex+ 1
     if (cpuid==0) then
        open(unit=outfileindex, file='Berrycurvature_line_sepband.gnu')
        write(outfileindex, '(a)') '# requirement: gnuplot version>5.4'
        write(outfileindex, '(2a)') '# Please open the data file to check the data: Berrycurvature_line_sepband.txt  '
        write(outfileindex, '(a)') 'set terminal pdf enhanced color font ",24"'
        write(outfileindex,'(2a)') 'set palette defined ( 0  "green", ', &
           '5 "yellow", 10 "red" )'
        write(outfileindex, '(3a)')"set output 'Berrycurvature_line_sepband.pdf'"
        write(outfileindex, '(a)')'set style data linespoints'
        write(outfileindex, '(a)')'unset key'
        write(outfileindex, '(a)')'set pointsize 0.8'
        write(outfileindex, '(a)')'#set xtics font ",24"'
        write(outfileindex, '(a)')'#set ytics font ",24"'
        write(outfileindex, '(a)')'#set ylabel font ",24"'
        write(outfileindex, '(a)')'set ylabel offset 0.5,0'
        write(outfileindex, '(a)')'set border lw 2'
        write(outfileindex, '(a, f10.5, a)')'set xrange [0: ', maxval(k3len*Angstrom2atomic), ']'
        write(outfileindex, '(a,f12.6)')'emin=', emin
        write(outfileindex, '(a,f12.6)')'emax=', emax
        write(outfileindex, '(a)')' '
        write(outfileindex, '(a)')'#set cbrange if the number is too large meeting some band crossings'
        write(outfileindex, '(a)')'set cbrange [-100:100]'
        write(outfileindex, '(a)')' '
        if (index(Particle,'phonon')/=0) then
           write(outfileindex, '(a)')'set yrange [0: emax]'
           write(outfileindex, '(a)')'set ylabel "Frequency (THz)"'
        else
           write(outfileindex, '(a)')'set ylabel "Energy (eV)"'
           write(outfileindex, '(a)')'set yrange [ emin : emax ]'
        endif
        write(outfileindex, 202, advance="no") (k3line_name(i), k3line_stop(i)*Angstrom2atomic, i=1, Nk3lines)
        write(outfileindex, 203)k3line_name(Nk3lines+1), k3line_stop(Nk3lines+1)*Angstrom2atomic
  
        do i=1, Nk3lines-1
           if (index(Particle,'phonon')/=0) then
              write(outfileindex, 204)k3line_stop(i+1)*Angstrom2atomic, '0.0', k3line_stop(i+1)*Angstrom2atomic, 'emax'
           else
              write(outfileindex, 204)k3line_stop(i+1)*Angstrom2atomic, 'emin', k3line_stop(i+1)*Angstrom2atomic, 'emax'
           endif
        enddo
        write(outfileindex, '(4a)')"plot 'Berrycurvature_line_sepband.txt' u 1:2:5 ",  &
           " w lp lw 2 pt 7  ps 0.2 lc palette, 0 w l lw 2 dt 2"
        close(outfileindex)
     endif

202 format('set xtics (',20('"',A3,'" ',F10.5,','))
203 format(A3,'" ',F10.5,')')
204 format('set arrow from ',F10.5,',',A5,' to ',F10.5,',',A5, ' nohead lw 2')

     ENDIF


#if defined (MPI)
     call mpi_barrier(mpi_cmw, ierr)
#endif

     return

  end subroutine berry_curvature_cube

  subroutine Berry_curvature_plane_full
     !> Calculate Berry curvature and orbital magnetization for the selected bands
     !> ref : Physical Review B 74, 195118(2006)
     !> eqn (34)
     !
     !> August 20 2020 by Quansheng Wu @ EPFL
     !
     ! Copyright (c) 2020 QuanSheng Wu. All rights reserved.

     use wmpi
     use para
     implicit none
    
     integer :: ik, i, j, m, n, ierr, nkmesh2

     real(dp) :: k(3), o1(3), vmin, vmax, kxy_plane(3)

     !> k points slice
     real(dp), allocatable :: kslice(:, :), kslice_xyz(:, :)
   
     real(dp) :: time_start, time_end, time_start0

     real(dp), allocatable :: W(:)
     !> Berry curvature and orbital magnetization (3, bands, k)
     real(dp), allocatable :: Omega_allk_Occ(:, :, :), Omega_allk_Occ_mpi(:, :, :)
     real(dp), allocatable :: m_OrbMag_allk_Occ(:, :, :), m_OrbMag_allk_Occ_mpi(:, :, :)
     real(dp), allocatable :: Omega_allk_EF(:, :, :), Omega_allk_EF_mpi(:, :, :)
     real(dp), allocatable :: m_OrbMag_allk_EF(:, :, :), m_OrbMag_allk_EF_mpi(:, :, :)
     real(dp), allocatable :: Omega_BerryCurv(:, :), m_OrbMag(:, :)
     real(dp) :: beta_fake, fermi, delta, norm

     nkmesh2= Nk1*Nk2
     allocate( Omega_BerryCurv(Num_wann, 3), m_OrbMag(Num_wann, 3))
     allocate( Omega_allk_Occ    (2, 3, nkmesh2))
     allocate( Omega_allk_Occ_mpi(2, 3, nkmesh2))
     allocate( m_OrbMag_allk_Occ    (2, 3, nkmesh2))
     allocate( m_OrbMag_allk_Occ_mpi(2, 3, nkmesh2))
     allocate( Omega_allk_EF    (2, 3, nkmesh2))
     allocate( Omega_allk_EF_mpi(2, 3, nkmesh2))
     allocate( m_OrbMag_allk_EF    (2, 3, nkmesh2))
     allocate( m_OrbMag_allk_EF_mpi(2, 3, nkmesh2))
     allocate( kslice(3, nkmesh2))
     allocate( kslice_xyz(3, nkmesh2))
     allocate( W(Num_wann))
     W=0d0
     m_OrbMag=0d0
     Omega_BerryCurv= 0d0
     Omega_allk_Occ= 0d0
     m_OrbMag_allk_Occ=0d0
     Omega_allk_Occ_mpi= 0d0
     m_OrbMag_allk_Occ_mpi=0d0
     Omega_allk_EF= 0d0
     m_OrbMag_allk_EF=0d0
     Omega_allk_EF_mpi= 0d0
     m_OrbMag_allk_EF_mpi=0d0

     kslice=0d0
     kslice_xyz=0d0
    
     !> kslice is centered at K3d_start
     ik =0
     do i= 1, nk1
        do j= 1, nk2
           ik =ik +1
           kslice(:, ik)= K3D_start+ K3D_vec1*(i-1)/dble(nk1-1)  &
                     + K3D_vec2*(j-1)/dble(nk2-1) - (K3D_vec1+ K3D_vec2)/2d0
           kslice_xyz(:, ik)= kslice(1, ik)* Origin_cell%Kua+ kslice(2, ik)* Origin_cell%Kub+ kslice(3, ik)* Origin_cell%Kuc 
        enddo
     enddo

     time_start= 0d0
     time_start0= 0d0
     call now(time_start0)
     time_start= time_start0
     time_end  = time_start0
     do ik= 1+ cpuid, nkmesh2, num_cpu
        if (cpuid==0.and. mod((ik-1)/num_cpu, 100)==0) &
           write(stdout, '(a, i9, "  /", i10, a, f10.1, "s", a, f10.1, "s")') &
           ' Berry curvature: ik', ik, nkmesh2, ' time left', &
           (nkmesh2-ik)*(time_end- time_start)/num_cpu, &
           ' time elapsed: ', time_end-time_start0 


        call now(time_start)
 
        !> diagonalize hamiltonian
        k= kslice(:, ik)

        call now(time_start)
        call Berry_curvature_orb_mag_singlek_allbands_pack(k, Omega_BerryCurv, m_OrbMag, W)
       
        do i=1, 3
           Omega_allk_Occ(1, i, ik) = sum(Omega_BerryCurv(1:NumOccupied, i))
           m_OrbMag_allk_Occ(1, i, ik) = sum(m_OrbMag(1:NumOccupied, i))
           Omega_allk_Occ(2, i, ik) = Omega_BerryCurv(NumOccupied, i)
           m_OrbMag_allk_Occ(2, i, ik) = m_OrbMag(NumOccupied, i)
        enddo

        !> consider the Fermi-distribution according to the broadening Earc_eta
        if (Fermi_broadening<eps6) Fermi_broadening= eps6
        beta_fake= 1d0/Fermi_broadening
        do i=1, 3
           do m= 1, Num_wann
              Omega_allk_EF(1, i, ik)= Omega_allk_EF(1, i, ik)+ Omega_BerryCurv(m, i)*fermi(W(m)-iso_energy, beta_fake)
              m_OrbMag_allk_EF(1, i, ik)= m_OrbMag_allk_EF(1, i, ik)+ m_OrbMag(m, i)*fermi(W(m)-iso_energy, beta_fake)
              Omega_allk_EF(2, i, ik)= Omega_allk_EF(2, i, ik)+ Omega_BerryCurv(m, i)*delta(Fermi_broadening, W(m)-iso_energy)
              m_OrbMag_allk_EF(2, i, ik)= m_OrbMag_allk_EF(2, i, ik)+ m_OrbMag(m, i)*delta(Fermi_broadening, W(m)-iso_energy)
           enddo
        enddo


        call now(time_end)
     enddo ! ik

#if defined (MPI)
     call mpi_allreduce(Omega_allk_Occ,Omega_allk_Occ_mpi,size(Omega_allk_Occ_mpi),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
     call mpi_allreduce(m_OrbMag_allk_Occ,m_OrbMag_allk_Occ_mpi,size(m_OrbMag_allk_Occ_mpi),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
     call mpi_allreduce(Omega_allk_EF,Omega_allk_EF_mpi,size(Omega_allk_EF_mpi),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
     call mpi_allreduce(m_OrbMag_allk_EF,m_OrbMag_allk_EF_mpi,size(m_OrbMag_allk_EF_mpi),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
     Omega_allk_Occ_mpi= Omega_allk_Occ
     m_OrbMag_allk_Occ_mpi= m_OrbMag_allk_Occ
     Omega_allk_EF_mpi= Omega_allk_EF
     m_OrbMag_allk_EF_mpi= m_OrbMag_allk_EF
#endif

     !> write the Berry curvature to file
     outfileindex= outfileindex+ 1
     if (cpuid==0) then
        open(unit=outfileindex, file='Berrycurvature.dat')
        write(outfileindex, '(a)')'# Berry curvature in unit of (Angstrom^2)  '
        write(outfileindex, '("#", 30X,4a36)')'Sum over 1-NumOccupied bands', &
           "values at NumOccupied'th band", &
           ' Sum below Fermi level','values at Fermi level'
        write(outfileindex, '(a10,2000a12)')'# kx (1/A)', 'ky (1/A)', 'kz (1/A)', &
                                            "kx' (1/A)", "ky' (1/A)", "kz' (1/A)", &
            'Omega_x', 'Omega_y', 'Omega_z', &
            'Omega_x', 'Omega_y', 'Omega_z', &
            'Omega_x', 'Omega_y', 'Omega_z', &
            'Omega_x', 'Omega_y', 'Omega_z'
        write(outfileindex, '("#col ", i5, 200i12)')(i, i=1,18)

        ik= 0
        do i= 1, nk1
           do j= 1, nk2
              ik= ik+ 1
              call rotate_k3_to_kplane(kslice_xyz(:, ik), kxy_plane)
              write(outfileindex, '(6f12.6,2000E12.4)')kslice_xyz(:, ik)*Angstrom2atomic, kxy_plane*Angstrom2atomic, &
                 Omega_allk_Occ_mpi(1, :, ik)/Angstrom2atomic**2, &
                 Omega_allk_Occ_mpi(2, :, ik)/Angstrom2atomic**2, &
                 Omega_allk_EF_mpi(1, :, ik)/Angstrom2atomic**2, &
                 Omega_allk_EF_mpi(2, :, ik)/Angstrom2atomic**2
           enddo
           write(outfileindex, *) ' '
        enddo

        close(outfileindex)

     endif

     !> Convert to unit of Bohr magneton
     m_OrbMag_allk_EF_mpi= m_OrbMag_allk_EF_mpi*eV2Hartree*Ang2Bohr**2*2
     m_OrbMag_allk_Occ_mpi= m_OrbMag_allk_Occ_mpi*eV2Hartree*Ang2Bohr**2*2


     !> write the orbital magnetization to file
     outfileindex= outfileindex+ 1
     if (cpuid==0) then
        open(unit=outfileindex, file='Orbitalmagnetization.dat')
        write(outfileindex, '(a)')'# Orbital magnetization in unit of Bohr magneton \mu_B'
        write(outfileindex, '("#col ", i5, 200i12)')(i, i=1, 15)
        write(outfileindex, '("#", 30X,4a36)')'Sum over 1~NumOccupied bands', &
           "values at NumOccupied'th band", &
           ' Sum below Fermi level','values at Fermi level'
        write(outfileindex, '(a10,2000a12)')'# kx (1/A)', 'ky (1/A)', 'kz (1/A)', &
                                            "kx' (1/A)", "ky' (1/A)", "kz' (1/A)", &
              'm_x', 'm_y', 'm_z', 'm_x', 'm_y', 'm_z', &
              'm_x', 'm_y', 'm_z', 'm_x', 'm_y', 'm_z'

        ik= 0
        do i= 1, nk1
           do j= 1, nk2
              ik= ik+ 1
              write(outfileindex, '(6f12.6,2000E12.4)')kslice_xyz(:, ik)*Angstrom2atomic, kxy_plane*Angstrom2atomic, &
                 m_OrbMag_allk_Occ_mpi(1, :, ik),  m_OrbMag_allk_Occ_mpi(2, :, ik), &
                 m_OrbMag_allk_EF_mpi(1, :, ik), m_OrbMag_allk_EF_mpi(2, :, ik)
           enddo
           write(outfileindex, *) ' '
        enddo

        close(outfileindex)

     endif

     !> generate gnuplot script to plot the Berry curvature
     outfileindex= outfileindex+ 1
     if (cpuid==0) then
        vmin=minval(Omega_allk_Occ_mpi(1, :, :))
        vmax=maxval(Omega_allk_Occ_mpi(1, :, :))
        open(unit=outfileindex, file='Berrycurvature.gnu')
        write(outfileindex, '(a)')"set encoding iso_8859_1"
        write(outfileindex, '(a)')'set terminal  pngcairo  truecolor enhanced size 3680, 1920 font ",40"'
        write(outfileindex, '(a)')'#set terminal  png       truecolor enhanced size 3680, 1920 font ",40"'
        write(outfileindex, '(a)')"set output 'Berrycurvature.png'"
        write(outfileindex, '(a)')'if (!exists("MP_LEFT"))   MP_LEFT = .12'
        write(outfileindex, '(a)')'if (!exists("MP_RIGHT"))  MP_RIGHT = .92'
        write(outfileindex, '(a)')'if (!exists("MP_BOTTOM")) MP_BOTTOM = .12'
        write(outfileindex, '(a)')'if (!exists("MP_TOP"))    MP_TOP = .88'
        write(outfileindex, '(a)')'if (!exists("MP_GAP"))    MP_GAP = 0.08'
        write(outfileindex, '(3a)')'set label 1 "Sum over all bands below Occupied', &
          "'th band", ' " at screen 0.5 ,0.98 center'
        write(outfileindex, '(a)')'set multiplot layout 1,3 rowsfirst \'
        write(outfileindex, '(a)')"              margins screen MP_LEFT, MP_RIGHT, MP_BOTTOM, MP_TOP spacing screen MP_GAP"
        write(outfileindex, '(a)')" "
        write(outfileindex, '(a)')"set palette rgbformulae 33,13,10"
        write(outfileindex, '(a)')"unset ztics"
        write(outfileindex, '(a)')"unset key"
        write(outfileindex, '(a)')"set pm3d"
        write(outfileindex, '(a)')"#set zbrange [ -10: 10] "
        write(outfileindex, '(a, E10.3, a, E10.3, a)')"set cbrange [ ", vmin, ':', vmax, " ] "
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
        write(outfileindex, '(a)')"set title 'Berry Curvature {/Symbol W}_x ({\305}^2)'"
        write(outfileindex, '(a)')"splot 'Berrycurvature.dat' u 4:5:7 w pm3d"
        write(outfileindex, '(a)')"unset ylabel"
        write(outfileindex, '(a)')"unset ytics"
        write(outfileindex, '(a)')"set title 'Berry Curvature {/Symbol W}_y ({\305}^2)'"
        write(outfileindex, '(a)')"splot 'Berrycurvature.dat' u 4:5:8 w pm3d"
        write(outfileindex, '(a)')"set title 'Berry Curvature {/Symbol W}_z ({\305}^2)'"
        write(outfileindex, '(a)')"set colorbox"
        write(outfileindex, '(a)')"splot 'Berrycurvature.dat' u 4:5:9 w pm3d"

        write(outfileindex, '(a)')"unset multiplot"
        write(outfileindex, '(a)')"unset label 1"
        write(outfileindex, '(3a)')'set label 2 "Sum over all bands below Fermi level', &
          ' E\\\_arc" at screen 0.5 ,0.98 center'
        write(outfileindex, '(a)')"set output 'Berrycurvature_EF.png'"
        write(outfileindex, '(a)')'set multiplot layout 1,3 rowsfirst \'
        write(outfileindex, '(a)')"              margins screen MP_LEFT, MP_RIGHT, MP_BOTTOM, MP_TOP spacing screen MP_GAP"
        vmin=minval(Omega_allk_EF_mpi(1, :, :))
        vmax=maxval(Omega_allk_EF_mpi(1, :, :))
        open(unit=outfileindex, file='Berrycurvature.gnu')
        write(outfileindex, '(a, E10.3, a, E10.3, a)')"set cbrange [ ", vmin, ':', vmax, " ] "
        write(outfileindex, '(a)')"set title 'Berry Curvature {/Symbol W}_x ({\305}^2)'"
        write(outfileindex, '(a)')"splot 'Berrycurvature.dat' u 4:5:10 w pm3d"
        write(outfileindex, '(a)')"unset ylabel"
        write(outfileindex, '(a)')"unset ytics"
        write(outfileindex, '(a)')"set title 'Berry Curvature {/Symbol W}_y ({\305}^2)'"
        write(outfileindex, '(a)')"splot 'Berrycurvature.dat' u 4:5:11 w pm3d"
        write(outfileindex, '(a)')"set title 'Berry Curvature {/Symbol W}_z ({\305}^2)'"
        write(outfileindex, '(a)')"set colorbox"
        write(outfileindex, '(a)')"splot 'Berrycurvature.dat' u 4:5:12 w pm3d"


        close(outfileindex)
     endif

     !> generate gnuplot script to plot the Berry curvature
     outfileindex= outfileindex+ 1
     if (cpuid==0) then

        open(unit=outfileindex, file='Orbitalmagnetization.gnu')
        write(outfileindex, '(a)')"set encoding iso_8859_1"
        write(outfileindex, '(a)')'set terminal  pngcairo  truecolor enhanced size 3680, 1920 font ",40"'
        write(outfileindex, '(a)')'#set terminal  png       truecolor enhanced size 3680, 1920 font ",40"'
        write(outfileindex, '(a)')"set output 'Orbitalmagnetization.png'"
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
        write(outfileindex, '(a)')"#set cbrange [  -100 : 100 ] "
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
        write(outfileindex, '(a)')"set title 'Orbital magnetization m_x ({/Symbol m}_B)'"
        write(outfileindex, '(a)')"splot 'Orbitalmagnetization.dat' u 4:5:7 w pm3d"
        write(outfileindex, '(a)')"unset ylabel"
        write(outfileindex, '(a)')"unset ytics"
        write(outfileindex, '(a)')"set title 'Orbital magnetization m_y ({/Symbol m}_B)'"
        write(outfileindex, '(a)')"splot 'Orbitalmagnetization.dat' u 4:5:8 w pm3d"
        write(outfileindex, '(a)')"set title 'Orbital magnetization m_{z} ({/Symbol m}_B)'"
        write(outfileindex, '(a)')"set colorbox"
        write(outfileindex, '(a)')"splot 'Orbitalmagnetization.dat' u 4:5:9 w pm3d"
        close(outfileindex)
     endif


     !> write the normalized Berry curvature to file in order to generate vector plot
     outfileindex= outfileindex+ 1
     if (cpuid==0) then
        open(unit=outfileindex, file='Berrycurvature-normalized.dat')
        write(outfileindex, '(a)')'# Omega(k)/|Omega(k)|, Omega(k) comes from the 4,5,6 columns of Berrycurvature.dat.'
        write(outfileindex, '("#col ", i5, 200i12)')(i, i=1, 6)
        write(outfileindex, '(20a12)')'# kx (1/A)', 'ky (1/A)', 'kz (1/A)', &
           'Omega_x', 'Omega_y', 'Omega_z'

        ik= 0
        do i= 1, nk1
           do j= 1, nk2
              ik= ik+ 1
              o1= Omega_allk_Occ_mpi(1, :, ik)
              if (norm(o1)>eps4) o1= o1/norm(o1)
              write(outfileindex, '(20f12.6)')kslice_xyz(:, ik)*Angstrom2atomic, o1
           enddo
           write(outfileindex, *) ' '
        enddo
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
        write(outfileindex, '(a)')"#unset xtics"
        write(outfileindex, '(a)')"#unset xlabel"
        write(outfileindex, '(a)')"set xrange [] noextend"
        write(outfileindex, '(a)')"set yrange [] noextend"
        write(outfileindex, '(a)')"set ytics 0.5 nomirror scale 0.5"
        write(outfileindex, '(a)')"set pm3d interpolate 2,2"
        write(outfileindex, '(a)')"set title '({/Symbol W}_x, {/Symbol W}_y)'"
        write(outfileindex, '(a)')"plot 'Berrycurvature-normalized.dat' u 4:5:($7/10):($8/10) w vec"
        close(outfileindex)
     endif




#if defined (MPI)
     call mpi_barrier(mpi_cmw, ierr)
#endif

     deallocate( kslice)
     deallocate( kslice_xyz)

     return

  end subroutine Berry_curvature_plane_full

  subroutine Berry_curvature_orb_mag_singlek_allbands_pack(k, Omega_BerryCurv, m_OrbMag, W)
     !> Calculate Berry curvature and orbital magnetization for the selected bands
     !> ref : Physical Review B 74, 195118(2006)
     !> eqn (34)
     !
     !> Sep 29 2020 by Quansheng Wu @ EPFL
     !
     ! Copyright (c) 2020 QuanSheng Wu. All rights reserved.

     use para, only : Num_wann, dp, zi, NumOccupied
     implicit none

     !> Input:  the k point coordinate in fractional coordinate
     real(dp), intent(in) :: k(3)
     
     !> Output: Berry curvature and orbital magnetization (3, bands)
     real(dp), intent(out) :: m_OrbMag(Num_wann, 3)
     real(dp), intent(out) :: Omega_BerryCurv(Num_wann, 3)

     !> Ouput: eigenvalue
     real(dp), intent(out) :: W(Num_wann)
    
     integer :: i

     !> velocity matrix in Wannier basis 
     complex(dp), allocatable :: Vmn_wann(:, :, :)

     !> velocity matrix in Hamiltonian basis 
     complex(dp), allocatable :: Vmn_Ham(:, :, :), Vmn_Ham_nondiag(:, :, :)
     complex(dp), allocatable :: Dmn_Ham(:, :, :)

     !> Hamiltonian, eigenvalue and eigenvectors
     complex(dp), allocatable :: UU(:, :)
    
     allocate(Vmn_wann(Num_wann, Num_wann, 3), Vmn_Ham(Num_wann, Num_wann, 3))
     allocate(Dmn_Ham(Num_wann,Num_wann,3), Vmn_Ham_nondiag(Num_wann, Num_wann, 3))
     allocate(UU(Num_wann, Num_wann))
     W=0d0
     m_OrbMag=0d0
     Omega_BerryCurv= 0d0

     !> calculation bulk hamiltonian by a direct Fourier transformation of HmnR
     call ham_bulk_atomicgauge(k, UU)

     !> diagonalization by call zheev in lapack
     call eigensystem_c( 'V', 'U', Num_wann, UU, W)

     !> get velocity operator in Wannier basis
     call dHdk_atomicgauge(k, Vmn_wann)
     
     !> Rotate Vmn_wann from Wannier basis to Hamiltonian basis
     do i=1, 3
        call rotation_to_Ham_basis(UU, Vmn_wann(:, :, i), Vmn_Ham(:, :, i))
     enddo

     call get_Dmn_Ham(W, Vmn_Ham, Dmn_Ham)
     call get_Vmn_Ham_nondiag(Vmn_Ham, Vmn_Ham_nondiag)

     call Berry_curvature_singlek_allbands(Dmn_Ham, Omega_BerryCurv)
     call orbital_magnetization_singlek_allbands(Dmn_Ham, Vmn_Ham_nondiag, m_OrbMag)

     deallocate(Vmn_wann, Vmn_Ham, Vmn_Ham_nondiag)
     deallocate(UU, Dmn_Ham)
     return

  end subroutine Berry_curvature_orb_mag_singlek_allbands_pack




  subroutine Berry_curvature_plane_EF
     !> Calculate Berry curvature and orbital magnetization for the selected bands
     !> ref : Physical Review B 74, 195118(2006)
     !> eqn (34)
     !
     !> August 20 2020 by Quansheng Wu @ EPFL
     !
     ! Copyright (c) 2020 QuanSheng Wu. All rights reserved.

     use wmpi
     use para
     implicit none
    
     integer :: ik, i, j, n, ierr, nkmesh2

     real(dp) :: k(3), o1(3), vmin, vmax

     !> k points slice
     real(dp), allocatable :: kslice(:, :), kslice_xyz(:, :)
   
     real(dp), external :: norm
     
     real(dp) :: time_start, time_end, time_start0

     !> Berry curvature and orbital magnetization (3, bands, k)
     real(dp), allocatable :: Omega_allk(:, :, :), Omega_allk_mpi(:, :, :)
     real(dp), allocatable :: m_OrbMag_allk(:, :, :), m_OrbMag_allk_mpi(:, :, :)
     real(dp), allocatable :: Omega_BerryCurv(:, :), m_OrbMag(:, :)

     !> velocity matrix in Wannier basis 
     complex(dp), allocatable :: Vmn_wann(:, :, :)

     !> velocity matrix in Hamiltonian basis 
     complex(dp), allocatable :: Vmn_Ham(:, :, :), Vmn_Ham_nondiag(:, :, :)

     complex(dp), allocatable :: Dmn_Ham(:, :, :)

     !> Hamiltonian, eigenvalue and eigenvectors
     complex(dp), allocatable :: UU(:, :)
     real(dp), allocatable :: W(:)

     nkmesh2= Nk1*Nk2
     allocate(Vmn_wann(Num_wann, Num_wann, 3), Vmn_Ham(Num_wann, Num_wann, 3))
     allocate(Dmn_Ham(Num_wann,Num_wann,3), Vmn_Ham_nondiag(Num_wann, Num_wann, 3))
     allocate(W(Num_wann))
     allocate(UU(Num_wann, Num_wann))
     allocate( Omega_BerryCurv(Num_wann, 3), m_OrbMag(Num_wann, 3))
     allocate( Omega_allk    (Num_wann, 3, nkmesh2))
     allocate( Omega_allk_mpi(Num_wann, 3, nkmesh2))
     allocate( m_OrbMag_allk    (Num_wann, 3, nkmesh2))
     allocate( m_OrbMag_allk_mpi(Num_wann, 3, nkmesh2))
     allocate( kslice(3, nkmesh2))
     allocate( kslice_xyz(3, nkmesh2))
     m_OrbMag=0d0
     Omega_allk= 0d0
     m_OrbMag_allk=0d0
     Omega_allk_mpi= 0d0
     Omega_BerryCurv= 0d0
     m_OrbMag_allk_mpi=0d0

     kslice=0d0
     kslice_xyz=0d0
    
     !> kslice is centered at K3d_start
     ik =0
     do i= 1, nk1
        do j= 1, nk2
           ik =ik +1
           kslice(:, ik)= K3D_start+ K3D_vec1*(i-1)/dble(nk1-1)  &
                     + K3D_vec2*(j-1)/dble(nk2-1) - (K3D_vec1+ K3D_vec2)/2d0
           kslice_xyz(:, ik)= kslice(1, ik)* Origin_cell%Kua+ kslice(2, ik)* Origin_cell%Kub+ kslice(3, ik)* Origin_cell%Kuc 
        enddo
     enddo

     time_start= 0d0
     time_start0= 0d0
     call now(time_start0)
     time_start= time_start0
     time_end  = time_start0
     do ik= 1+ cpuid, nkmesh2, num_cpu
        if (cpuid==0.and. mod(ik/num_cpu, 100)==0) &
           write(stdout, '(a, i9, "  /", i10, a, f10.1, "s", a, f10.1, "s")') &
           ' Berry curvature: ik', ik, nkmesh2, ' time left', &
           (nkmesh2-ik)*(time_end- time_start)/num_cpu, &
           ' time elapsed: ', time_end-time_start0 


        call now(time_start)
 
        !> diagonalize hamiltonian
        k= kslice(:, ik)

        call now(time_start)
        
        !> calculation bulk hamiltonian by a direct Fourier transformation of HmnR
        call ham_bulk_atomicgauge(k, UU)
   
        !> diagonalization by call zheev in lapack
        call eigensystem_c( 'V', 'U', Num_wann, UU, W)

        !> get velocity operator in Wannier basis
        call dHdk_atomicgauge(k, Vmn_wann)
        
        !> Rotate Vmn_wann from Wannier basis to Hamiltonian basis
        do i=1, 3
           call rotation_to_Ham_basis(UU, Vmn_wann(:, :, i), Vmn_Ham(:, :, i))
        enddo

        call get_Dmn_Ham(W, Vmn_Ham, Dmn_Ham)
        call get_Vmn_Ham_nondiag(Vmn_Ham, Vmn_Ham_nondiag)

        call Berry_curvature_singlek_allbands(Dmn_Ham, Omega_BerryCurv)
        call orbital_magnetization_singlek_allbands(Dmn_Ham, Vmn_Ham_nondiag, m_OrbMag)
        Omega_allk(:, :, ik) = Omega_BerryCurv
        m_OrbMag_allk(:, :, ik) = m_OrbMag

        call now(time_end)
     enddo ! ik

#if defined (MPI)
     call mpi_allreduce(Omega_allk,Omega_allk_mpi,size(Omega_allk_mpi),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
     call mpi_allreduce(m_OrbMag_allk,m_OrbMag_allk_mpi,size(m_OrbMag_allk_mpi),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
     Omega_allk_mpi= Omega_allk
     m_OrbMag_allk_mpi= m_OrbMag_allk
#endif

     !> output the Berry curvature to file
     outfileindex= outfileindex+ 1
     if (cpuid==0) then
        open(unit=outfileindex, file='Berrycurvature_Orbitalmagnetization.dat')
        write(outfileindex, '("#col ", i5, 200i12)')(i, i=1, NumberofSelectedBands*6+3)
        write(outfileindex, '(a10,2000a12)')'# kx (1/A)', 'ky (1/A)', 'kz (1/A)', &
           'Omega_x(A^2)', 'Omega_y(A^2)', 'Omega_z(A^2)' , 'm_x', 'm_y', 'm_z'

        ik= 0
        do i= 1, nk1
           do j= 1, nk2
              ik= ik+ 1
              write(outfileindex, '(3f12.6,2000E12.4)')kslice_xyz(:, ik)*Angstrom2atomic, &
                 (Omega_allk_mpi(Selected_band_index(n), :, ik)/Angstrom2atomic**2, &
                 m_OrbMag_allk_mpi(Selected_band_index(n), :, ik), n=1, NumberofSelectedBands)   
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
        write(outfileindex, '(a)')'set terminal  pngcairo  truecolor enhanced size 3680, 1920 font ",40"'
        write(outfileindex, '(a)')'#set terminal  png       truecolor enhanced size 3680, 1920 font ",40"'
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
        write(outfileindex, '(a, f10.3, a, f10.3, a)')"#set cbrange [ ", vmin, ':', vmax, " ] "
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
        write(outfileindex, '(a)')"splot 'Berrycurvature_Orbitalmagnetization.dat' u 1:2:4 w pm3d"
        write(outfileindex, '(a)')"unset ylabel"
        write(outfileindex, '(a)')"unset ytics"
        write(outfileindex, '(a)')"set title 'Berry Curvature {/Symbol W}_y'"
        write(outfileindex, '(a)')"splot 'Berrycurvature_Orbitalmagnetization.dat' u 1:2:5 w pm3d"
        write(outfileindex, '(a)')"set title 'Berry Curvature {/Symbol W}_z'"
        write(outfileindex, '(a)')"set colorbox"
        write(outfileindex, '(a)')"splot 'Berrycurvature_Orbitalmagnetization.dat' u 1:2:6 w pm3d"
        close(outfileindex)
     endif


#if defined (MPI)
     call mpi_barrier(mpi_cmw, ierr)
#endif

     deallocate( kslice)
     deallocate( kslice_xyz)
 
     return

  end subroutine Berry_curvature_plane_EF



  subroutine Berry_curvature_plane_selectedbands
     !> Calculate Berry curvature and orbital magnetization for the selected bands
     !> ref : Physical Review B 74, 195118(2006)
     !> eqn (34)
     !
     !> August 20 2020 by Quansheng Wu @ EPFL
     !
     ! Copyright (c) 2020 QuanSheng Wu. All rights reserved.

     use wmpi
     use para
     implicit none
    
     integer :: ik, i, j, n, ierr, nkmesh2

     real(dp) :: k(3), o1(3), vmin, vmax

     !> k points slice
     real(dp), allocatable :: kslice(:, :), kslice_xyz(:, :)
   
     real(dp), external :: norm
     
     real(dp) :: time_start, time_end, time_start0

     !> Berry curvature and orbital magnetization (3, bands, k)
     real(dp), allocatable :: Omega_allk(:, :, :), Omega_allk_mpi(:, :, :)
     real(dp), allocatable :: m_OrbMag_allk(:, :, :), m_OrbMag_allk_mpi(:, :, :)
     real(dp), allocatable :: Omega_BerryCurv(:, :), m_OrbMag(:, :)

     !> velocity matrix in Wannier basis 
     complex(dp), allocatable :: Vmn_wann(:, :, :)

     !> velocity matrix in Hamiltonian basis 
     complex(dp), allocatable :: Vmn_Ham(:, :, :), Vmn_Ham_nondiag(:, :, :)

     complex(dp), allocatable :: Dmn_Ham(:, :, :)

     !> Hamiltonian, eigenvalue and eigenvectors
     complex(dp), allocatable :: UU(:, :)
     real(dp), allocatable :: W(:)

     nkmesh2= Nk1*Nk2
     allocate(Vmn_wann(Num_wann, Num_wann, 3), Vmn_Ham(Num_wann, Num_wann, 3))
     allocate(Dmn_Ham(Num_wann,Num_wann,3), Vmn_Ham_nondiag(Num_wann, Num_wann, 3))
     allocate(W(Num_wann))
     allocate(UU(Num_wann, Num_wann))
     allocate( Omega_BerryCurv(Num_wann, 3), m_OrbMag(Num_wann, 3))
     allocate( Omega_allk    (Num_wann, 3, nkmesh2))
     allocate( Omega_allk_mpi(Num_wann, 3, nkmesh2))
     allocate( m_OrbMag_allk    (Num_wann, 3, nkmesh2))
     allocate( m_OrbMag_allk_mpi(Num_wann, 3, nkmesh2))
     allocate( kslice(3, nkmesh2))
     allocate( kslice_xyz(3, nkmesh2))
     m_OrbMag=0d0
     Omega_allk= 0d0
     m_OrbMag_allk=0d0
     Omega_allk_mpi= 0d0
     Omega_BerryCurv= 0d0
     m_OrbMag_allk_mpi=0d0

     kslice=0d0
     kslice_xyz=0d0
    
     !> kslice is centered at K3d_start
     ik =0
     do i= 1, nk1
        do j= 1, nk2
           ik =ik +1
           kslice(:, ik)= K3D_start+ K3D_vec1*(i-1)/dble(nk1-1)  &
                     + K3D_vec2*(j-1)/dble(nk2-1) - (K3D_vec1+ K3D_vec2)/2d0
           kslice_xyz(:, ik)= kslice(1, ik)* Origin_cell%Kua+ kslice(2, ik)* Origin_cell%Kub+ kslice(3, ik)* Origin_cell%Kuc 
        enddo
     enddo

     time_start= 0d0
     time_start0= 0d0
     call now(time_start0)
     time_start= time_start0
     time_end  = time_start0
     do ik= 1+ cpuid, nkmesh2, num_cpu
        if (cpuid==0.and. mod(ik/num_cpu, 100)==0) &
           write(stdout, '(a, i9, "  /", i10, a, f10.1, "s", a, f10.1, "s")') &
           ' Berry curvature: ik', ik, nkmesh2, ' time left', &
           (nkmesh2-ik)*(time_end- time_start)/num_cpu, &
           ' time elapsed: ', time_end-time_start0 


        call now(time_start)
 
        !> diagonalize hamiltonian
        k= kslice(:, ik)

        call now(time_start)
        
        !> calculation bulk hamiltonian by a direct Fourier transformation of HmnR
        call ham_bulk_atomicgauge(k, UU)
   
        !> diagonalization by call zheev in lapack
        call eigensystem_c( 'V', 'U', Num_wann, UU, W)

        !> get velocity operator in Wannier basis
        call dHdk_atomicgauge(k, Vmn_wann)
        
        !> Rotate Vmn_wann from Wannier basis to Hamiltonian basis
        do i=1, 3
           call rotation_to_Ham_basis(UU, Vmn_wann(:, :, i), Vmn_Ham(:, :, i))
        enddo

        call get_Dmn_Ham(W, Vmn_Ham, Dmn_Ham)
        call get_Vmn_Ham_nondiag(Vmn_Ham, Vmn_Ham_nondiag)

        call Berry_curvature_singlek_allbands(Dmn_Ham, Omega_BerryCurv)
        call orbital_magnetization_singlek_allbands(Dmn_Ham, Vmn_Ham_nondiag, m_OrbMag)
        Omega_allk(:, :, ik) = Omega_BerryCurv
        m_OrbMag_allk(:, :, ik) = m_OrbMag

        call now(time_end)
     enddo ! ik

#if defined (MPI)
     call mpi_allreduce(Omega_allk,Omega_allk_mpi,size(Omega_allk_mpi),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
     call mpi_allreduce(m_OrbMag_allk,m_OrbMag_allk_mpi,size(m_OrbMag_allk_mpi),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
     Omega_allk_mpi= Omega_allk
     m_OrbMag_allk_mpi= m_OrbMag_allk
#endif

     !> output the Berry curvature to file
     outfileindex= outfileindex+ 1
     if (cpuid==0) then
        open(unit=outfileindex, file='Berrycurvature_Orbitalmagnetization.dat')
        write(outfileindex, '("#col ", i5, 200i12)')(i, i=1, NumberofSelectedBands*6+3)
        write(outfileindex, '(a10,2000a12)')'# kx (1/A)', 'ky (1/A)', 'kz (1/A)', &
           'Omega_x', 'Omega_y', 'Omega_z' , 'm_x', 'm_y', 'm_z'

        ik= 0
        do i= 1, nk1
           do j= 1, nk2
              ik= ik+ 1
              write(outfileindex, '(3f12.6,2000E12.4)')kslice_xyz(:, ik)*Angstrom2atomic, &
                 (Omega_allk_mpi(Selected_band_index(n), :, ik)/Angstrom2atomic**2, &
                 m_OrbMag_allk_mpi(Selected_band_index(n), :, ik), n=1, NumberofSelectedBands)   
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
        write(outfileindex, '(a)')'set terminal  pngcairo  truecolor enhanced size 3680, 1920 font ",40"'
        write(outfileindex, '(a)')'#set terminal  png       truecolor enhanced size 3680, 1920 font ",40"'
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
        write(outfileindex, '(a, f10.3, a, f10.3, a)')"#set cbrange [ ", vmin, ':', vmax, " ] "
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
        write(outfileindex, '(a)')"splot 'Berrycurvature_Orbitalmagnetization.dat' u 1:2:4 w pm3d"
        write(outfileindex, '(a)')"unset ylabel"
        write(outfileindex, '(a)')"unset ytics"
        write(outfileindex, '(a)')"set title 'Berry Curvature {/Symbol W}_y'"
        write(outfileindex, '(a)')"splot 'Berrycurvature_Orbitalmagnetization.dat' u 1:2:5 w pm3d"
        write(outfileindex, '(a)')"set title 'Berry Curvature {/Symbol W}_z'"
        write(outfileindex, '(a)')"set colorbox"
        write(outfileindex, '(a)')"splot 'Berrycurvature_Orbitalmagnetization.dat' u 1:2:6 w pm3d"
        close(outfileindex)
     endif


#if defined (MPI)
     call mpi_barrier(mpi_cmw, ierr)
#endif

     deallocate( kslice)
     deallocate( kslice_xyz)
 
     return

  end subroutine Berry_curvature_plane_selectedbands



  subroutine Berry_curvature_plane
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
    
     integer :: ik, i, j, ierr

     real(dp) :: k(3), o1(3), vmin, vmax

     !> k points slice
     real(dp), allocatable :: kslice(:, :), kslice_xyz(:, :)
   
     real(dp), external :: norm
     
     real(dp) :: time_start, time_end, time_start0

     !> Berry curvature  (3, bands, k)
     complex(dp), allocatable :: Omega_x(:), Omega_y(:), Omega_z(:)
     complex(dp), allocatable :: Omega(:, :), Omega_mpi(:, :)

     allocate( kslice(3, Nk1*Nk2))
     allocate( kslice_xyz(3, Nk1*Nk2))
     allocate( Omega_x(Num_wann))
     allocate( Omega_y(Num_wann))
     allocate( Omega_z(Num_wann))
     allocate( Omega    (3, Nk1*Nk2))
     allocate( Omega_mpi(3, Nk1*Nk2))
     kslice=0d0
     kslice_xyz=0d0
     omega= 0d0
     omega_mpi= 0d0
    
     !> kslice is centered at K3d_start
     ik =0
     do i= 1, nk1
        do j= 1, nk2
           ik =ik +1
           kslice(:, ik)= K3D_start+ K3D_vec1*(i-1)/dble(nk1-1)  &
                     + K3D_vec2*(j-1)/dble(nk2-1) - (K3D_vec1+ K3D_vec2)/2d0
           kslice_xyz(:, ik)= kslice(1, ik)* Origin_cell%Kua+ kslice(2, ik)* Origin_cell%Kub+ kslice(3, ik)* Origin_cell%Kuc 
        enddo
     enddo

     time_start= 0d0
     time_start0= 0d0
     call now(time_start0)
     time_start= time_start0
     time_end  = time_start0
     do ik= 1+ cpuid, Nk1*Nk2, num_cpu
        if (cpuid==0.and. mod((ik-1)/num_cpu, 100)==0) &
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

        !call Berry_curvature_singlek_numoccupied_old(k, Omega_x, Omega_y, Omega_z)
        if (Berrycurvature_EF_calc) then
           call Berry_curvature_singlek_EF(k, iso_energy, Omega_x, Omega_y, Omega_z)
        else
           call Berry_curvature_singlek_numoccupied_total(k, Omega_x(1), Omega_y(1), Omega_z(1))
        endif
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
           'Omega_x', 'Omega_y', 'Omega_z' 

        ik= 0
        do i= 1, nk1
           do j= 1, nk2
              ik= ik+ 1
              write(outfileindex, '(20E28.10)')kslice_xyz(:, ik)*Angstrom2atomic, &
                 real(Omega_mpi(:, ik))/Angstrom2atomic**2
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
           'Omega_x(A^2)', 'Omega_y(A^2)', 'Omega_z(A^2)'

        ik= 0
        do i= 1, nk1
           do j= 1, nk2
              ik= ik+ 1
              o1= real(Omega_mpi(:,ik))
              if (norm(o1)>eps9) o1= o1/norm(o1)
              write(outfileindex, '(20f28.10)')kslice_xyz(:, ik)*Angstrom2atomic, o1
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
        write(outfileindex, '(a)')'set terminal  pngcairo  truecolor enhanced size 3680, 1920 font ",40"'
        write(outfileindex, '(a)')'#set terminal  png       truecolor enhanced size 3680, 1920 font ",40"'
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
     deallocate( kslice_xyz)
     deallocate( Omega_x, Omega_y, Omega_z)
     deallocate( Omega, Omega_mpi)
 
     return

  end subroutine Berry_curvature_plane

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
        !call Berry_curvature_singlek_numoccupied_old(k_direct, Omega_x, Omega_y, Omega_z)
         call Berry_curvature_singlek_numoccupied    (k_direct, Omega_x, Omega_y, Omega_z)
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



  subroutine Chern_halftorus_single(k0, Rbig, rsmall_a, rsmall_b, Chern)
     !> calculate Chern number on a halftorus by integration over a halftorus
     !> QuanSheng Wu at EPFL June 12 2018
     !> wuquansheng@gmail.com
     !> The halftorus is defined by a center k0 and a radius in WEYL_CHIRALITY
     !> C= \int d\theta d\phi 
     !> k0 must be in Cartesian coordinates

     use wmpi
     use para
     implicit none
     
     !> inout variables
     real(dp), intent(in) :: k0(3)
     real(dp), intent(in) :: Rbig
     real(dp), intent(in) :: rsmall_a
     real(dp), intent(in) :: rsmall_b
     real(dp), intent(out) :: Chern
     real(dp), external :: norm


     integer :: ik, ik1, ik2, nkmesh2, ierr
     real(dp) :: theta, phi, r_para, dtheta, dphi, Chern_mpi
     real(dp) :: time_start, time_end
     real(dp) :: st, ct, sp, cp, O_x, O_y, O_z, rt

     !> k points in the halftorus
     real(dp) :: k_cart(3), k_direct(3), o1(3)
     real(dp), allocatable :: kpoints(:, :)
     real(dp), allocatable :: thetas(:)
     real(dp), allocatable :: phis(:)
     real(dp), allocatable :: Omega_k(:, :)
     real(dp), allocatable :: Omega_k_mpi(:, :)

     !> Berry curvature  (bands)
     complex(dp), allocatable :: Omega_x(:)
     complex(dp), allocatable :: Omega_y(:)
     complex(dp), allocatable :: Omega_z(:)

     nkmesh2= Nk1*Nk2
     allocate(Omega_x(Num_wann), Omega_y(Num_wann), Omega_z(Num_wann))
     allocate(thetas(nkmesh2), phis(nkmesh2), kpoints(3, nkmesh2))
     allocate(Omega_k(3, nkmesh2), Omega_k_mpi(3, nkmesh2))
     Omega_k= 0d0; Omega_k_mpi= 0d0

     !> set k plane
     !> the first dimension should be in one primitive cell, phi=[0, 2*pi]
     !> the first dimension is the integration direction
     !> the WCCs are calculated along the second k line theta=[0, pi]
     ik= 0
     do ik2=1, Nk2 ! the mesh for k line of that wcc are calculated ->theta
        theta= (ik2- 1d0)/(Nk2- 1d0)* 1d0*pi
        if (ik2== 1) theta= (ik2- 1d0+ 0.01)/(Nk2- 1d0)* pi  ! avoid the North pole
        if (ik2== Nk2) theta= (ik2- 1d0- 0.01)/(Nk2- 1d0)* pi  ! avoid the south pole
        do ik1=1, Nk1  ! the mesh along the integration direction -> phi
           ik = ik+ 1
           phi= (ik1- 1d0)/Nk1* 2.0d0* pi
           r_para= Rbig+ rsmall_a* cos(theta)
           k_cart(1)= k0(1)+ r_para* cos(phi)
           k_cart(2)= k0(2)+ r_para* sin(phi)
           k_cart(3)= k0(3)+ rsmall_b* sin(theta)
           call cart_direct_rec(k_cart, k_direct)
           kpoints(:, ik)= k_direct
           thetas(ik)= theta
           phis(ik) = phi
        enddo
     enddo
     dtheta= 1d0/(Nk2- 1d0)*pi
     dphi= 1d0/Nk1*2d0*pi

     Chern_mpi = 0d0
     time_start= 0d0
     time_end= 0d0
     call now(time_start)
     time_end= time_start
     do ik=cpuid+1, Nk1*Nk2, num_cpu
        if (cpuid==0.and. mod(ik/num_cpu, 100)==0) &
           write(stdout, '(a, i8, a, i7, a, f16.1, a)') &
           'Chern_halftorus_single, ik ', ik, ' nkmesh2 ',nkmesh2, ' time left ', &
           (nkmesh2-ik)*(time_end- time_start)/num_cpu, ' s'
        call now(time_start)
 
        k_direct= kpoints(:, ik)
        theta= thetas(ik)
        phi= phis(ik)
        st=sin(theta)
        ct=cos(theta)
        sp=sin(phi)
        cp=cos(phi)
        rt= Rbig+ rsmall_a* ct
        call Berry_curvature_singlek_numoccupied(k_direct, Omega_x, Omega_y, Omega_z)
       !O_x= real(Omega_x(Numoccupied))
       !O_y= real(Omega_y(Numoccupied))
       !O_z= real(Omega_z(Numoccupied))
        O_x= real(sum(Omega_x))
        O_y= real(sum(Omega_y))
        O_z= real(sum(Omega_z))
        Chern_mpi= Chern_mpi-( rsmall_b*rt*ct*cp*O_x+ rsmall_b*rt*ct*sp*O_y+rsmall_a*rt*st*O_z)
        call now(time_end)
        Omega_k_mpi(1, ik)= O_x
        Omega_k_mpi(2, ik)= O_y
        Omega_k_mpi(3, ik)= O_z
     enddo

#if defined (MPI)
     Chern= 0d0
     call mpi_allreduce(Chern_mpi, Chern, 1, &
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
     call mpi_allreduce(Omega_k_mpi, Omega_k, size(Omega_k), &
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
     Chern= Chern_mpi 
     Omega_k= Omega_k_mpi
#endif

     outfileindex= outfileindex+ 1
     if (cpuid==0) then
        open(unit=outfileindex, file='halftorus.txt')
        write(outfileindex, '(a, 3f12.5)')'# Torus centered at : ', k0
        write(outfileindex, '(a, f12.5)') '# Rbig (1/A) ', Rbig
        write(outfileindex, '(a, f12.5)') '# rsmall_a (1/A) ', rsmall_a
        write(outfileindex, '(a, f12.5)') '# rsmall_b (1/A) ', rsmall_b
        do ik=1, Nk1*Nk2
           k_direct= kpoints(:, ik)
           call direct_cart_rec(k_direct, k_cart)
           o1= Omega_k(:, ik)
           if (norm(o1)>eps9) o1= o1/norm(o1)
           write(outfileindex, '(6f16.5)') k_cart, o1
        enddo
        close(outfileindex)
     endif



      Chern = Chern* dtheta* dphi/pi/2d0

     return
  end subroutine Chern_halftorus_single

  subroutine Chern_halftorus
     use wmpi
     use para
     implicit none

     integer :: i
     real(dp) :: k0(3), Chern
     real(dp), allocatable :: Chern_array(:)

     allocate(Chern_array(Num_NLs))
     Chern_array= 0d0

     do i=1, Num_NLs
        k0= NL_center_position_cart(:, i)
        call Chern_halftorus_single(k0, Rbig_NL, rsmall_a_NL, rsmall_b_NL, Chern)
        Chern_array(i)= Chern
     enddo

     if (cpuid==0)then
        write(stdout, *)'# Chern number for the Weyl points'
        write(stdout, '("#",a8,5a9, a16)')'kx', 'ky', 'kz', 'k1', 'k2', 'k3', 'Chern'
        do i=1, Num_NLs
           write(stdout, '(6f9.5, f16.8)')NL_center_position_cart(:, i), NL_center_position_direct(:, i), Chern_array(i)
        enddo
     endif

     return
  end subroutine Chern_halftorus


