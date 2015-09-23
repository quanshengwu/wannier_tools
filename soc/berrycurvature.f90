!> calculate Berry curvature 
!> ref : Physical Review B 74, 195118(2006)
!> eqn (34)
!> Sep. 22 2015 by Quansheng Wu @ ETHZ
  subroutine berry_curvarture

     use para
     implicit none
    
     integer :: iR
     integer :: ia
     integer :: ib
     integer :: ic
     integer :: ik

     integer :: knv2

     real(dp) :: k1min, k1max, k2min, k2max, kz
   
     !> R points coordinates (3, nrpts)
     real(dp) :: crvec(:, :)
     real(dp) :: R(3)

     !> k points slice
     real(dp) :: kslice(:, :)
   
     ! eigen value of H
	  real(dp), allocatable :: eigv(:,:)
	  real(dp), allocatable :: eigv_mpi(:,:)

     !> velocities
     complex(dp), allocatable :: vx(:, :)
     complex(dp), allocatable :: vy(:, :)
     complex(dp), allocatable :: vz(:, :)
     complex(dp), allocatable :: DHDk(:, :, :)
     complex(dp), allocatable :: DHDkdag(:, :, :)
    
     !> Berry curvature  (3, bands, k)
     complex(dp), allocatable :: Omega(:, :)
     complex(dp), allocatable :: Omega_mpi(:, :)


     nk1= Nk
     nk2= Nk
     knv2= nk1*nk2
     allocate( kslice(2, knv2))
     allocate( kslice_shape(2, knv2))
     allocate( eigv    (Num_wann, knv3))
     allocate( eigv_mpi(Num_wann, knv3))
     allocate( Omega    (3, knv3))
     allocate( Omega_mpi(3, knv3))
     allocate( DHDk    (Num_wann, Num_wann, 3))
     allocate( DHDkdag (Num_wann, Num_wann, 3))
     kslice=0d0
     kslice_shape=0d0

     k1min= 0.10d0/1d0
     k1max= 0.30d0/1d0
     k2min=-0.05d0/1d0
     k2max= 0.05d0/1d0
     kz= 0.0d0
     ik =0
     do i= 1, nk1
        do j= 1, nk2
           ik =ik +1
           kslice(1, ik)=k1min+ (i-1)*(k1max-k1min)/dble(nk1-1)
           kslice(2, ik)=k2min+ (j-1)*(k2max-k2min)/dble(nk2-1)
           kslice_shape(1, ik)= kslice(1, ik)* Kua(1)+ kslice(2, ik)* Kub(1)
           kslice_shape(2, ik)= kslice(1, ik)* Kua(2)+ kslice(2, ik)* Kub(2)
        enddo
     enddo
     k2min_shape=minval(kslice_shape(2,:))
     k2max_shape=maxval(kslice_shape(2,:))
     k1min_shape=minval(kslice_shape(1,:))
     k1max_shape=maxval(kslice_shape(1,:))

     !> sweep over k points
     do ik= 1+ cpuid, knv2, num_cpu
        k(1)= kslice(1, ik)
        k(2)= kslice(2, ik)
        k(3)= kz

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

     !> get R coordinates 
     do iR=1, Nrpts
        crvec(:, iR)= Rua*irvec(1,iR) + Rub*irvec(2,iR) + Ruc*irvec(3,iR)
     enddo

     do ik= 1+ cpuid, knv2, num_cpu
        ham= 0d0
        vx= 0d0
        vy= 0d0
        vz= 0d0
        do iR= 1, Nrpts
           kdotr= k(1)*irvec(1,iR) + k(2)*irvec(2,iR) + k(3)*irvec(3,iR)
           vx= vx+ crvec(1, iR)*HmnR(:,:,iR)*Exp(pi2zi*kdotr)/ndegen(iR)
           vy= vy+ crvec(2, iR)*HmnR(:,:,iR)*Exp(pi2zi*kdotr)/ndegen(iR)
           vz= vz+ crvec(3, iR)*HmnR(:,:,iR)*Exp(pi2zi*kdotr)/ndegen(iR)
        enddo ! iR

        do n= 1, Num_wann
           do m= 1, Num_wann
              if (eigv(n, ik) > E_fermi .and. eigv(m, ik)< E_fermi) then
                 DHDk(n, m, 1)= zi*vx(n, m)/(eigv(m, ik)-eigv(n, ik))
                 DHDk(n, m, 2)= zi*vy(n, m)/(eigv(m, ik)-eigv(n, ik))
                 DHDk(n, m, 3)= zi*vz(n, m)/(eigv(m, ik)-eigv(n, ik))
              else
                 DHDk(n, m, 1)= 0d0
                 DHDk(n, m, 2)= 0d0
                 DHDk(n, m, 3)= 0d0
              endif
           enddo ! m
        enddo ! n

        do n= 1, Num_wann
           do m= 1, Num_wann
              if (eigv(m, ik) > E_fermi .and. eigv(n, ik)< E_fermi) then
                 DHDkdag(n, m, 1)= zi*vx(n, m)/(eigv(m, ik)-eigv(n, ik))
                 DHDkdag(n, m, 2)= zi*vy(n, m)/(eigv(m, ik)-eigv(n, ik))
                 DHDkdag(n, m, 3)= zi*vz(n, m)/(eigv(m, ik)-eigv(n, ik))
              else              
                 DHDkdag(n, m, 1)= 0d0
                 DHDkdag(n, m, 2)= 0d0
                 DHDkdag(n, m, 3)= 0d0
              endif
           enddo ! m
        enddo ! n

        call mat_mul(Num_wann, DHDk(:, :, 1), DHDkdag(:, :, 2), vz)
        call mat_mul(Num_wann, DHDk(:, :, 2), DHDkdag(:, :, 3), vx)
        call mat_mul(Num_wann, DHDk(:, :, 3), DHDkdag(:, :, 1), vy)

        call trace(Num_wann, vx, Omega(1, ik))
        call trace(Num_wann, vy, Omega(1, ik))
        call trace(Num_wann, vz, Omega(1, ik))

     enddo ! ik

     return

  end subroutine berry_curvarture

  subroutine Fourier_R_to_k(k, ham)
     use para
     implicit none

     real(dp), intent(in) :: k(3)
     complex(dp), intent(out) :: ham(Num_wann, Num_wann)
     integer :: iR
     real(dp) :: R(3)

     ham= 0d0
     do iR= 1, Nrpts
        R= Rua*irvec(1,iR) + Rub*irvec(2,iR) + Ruc*irvec(3,iR)
        kdotr= k(1)*irvec(1,iR) + k(2)*irvec(2,iR) + k(3)*irvec(3,iR)
        Ham= Ham+ HmnR(:,:,iR)*Exp(2d0*pi*zi*kdotr)/ndegen(iR)
     enddo

     return
  end subroutine Fourier_R_to_k

  subroutine trace(ndim, A, tr)
     use para, only : dp
     implicit none
     integer :: ndim
     complex(dp), intent(out) :: tr
     complex(dp), intent(in) :: A(ndim, ndim)

     integer :: i

     tr = 0d0
     do i=1, ndim
        tr= tr+ A(i, i)
     enddo

     return
  end subroutine trace
