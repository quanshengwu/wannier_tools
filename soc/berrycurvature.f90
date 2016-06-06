!> calculate Berry curvature 
!> ref : Physical Review B 74, 195118(2006)
!> eqn (34)
!> Sep. 22 2015 by Quansheng Wu @ ETHZ
! Copyright (c) 2010 QuanSheng Wu. All rights reserved.
  subroutine berry_curvarture

     use para
     implicit none
    
     integer :: iR
     integer :: ik

     integer :: m, n, i, j

     integer :: ierr


     real(dp) :: kdotr
     real(dp) :: k(3)

     !> R points coordinates (3, nrpts)
     real(dp), allocatable :: crvec(:, :)

     !> k points slice
     real(dp), allocatable :: kslice(:, :)
     real(dp), allocatable :: kslice_shape(:, :)
   
     ! eigen value of H
	  real(dp), allocatable :: W(:)
     complex(dp), allocatable :: Hamk_bulk(:, :)
     complex(dp), allocatable :: Amat(:, :)
     complex(dp), allocatable :: UU(:, :)
     complex(dp), allocatable :: UU_dag(:, :)

     !> velocities
     complex(dp), allocatable :: vx(:, :)
     complex(dp), allocatable :: vy(:, :)
     complex(dp), allocatable :: vz(:, :)
     complex(dp), allocatable :: DHDk(:, :, :)
     complex(dp), allocatable :: DHDkdag(:, :, :)
    
     !> Berry curvature  (3, bands, k)
     complex(dp), allocatable :: Omega(:, :)
     complex(dp), allocatable :: Omega_mpi(:, :)


     knv2= nk1*nk2
     allocate( kslice(3, knv2))
     allocate( kslice_shape(3, knv2))
     allocate( W       (Num_wann))
     allocate( crvec    (3, nrpts))
     allocate( Omega    (3, knv2))
     allocate( Omega_mpi(3, knv2))
     allocate( vx      (Num_wann, Num_wann))
     allocate( vy      (Num_wann, Num_wann))
     allocate( vz      (Num_wann, Num_wann))
     allocate( Hamk_bulk(Num_wann, Num_wann))
     allocate( Amat(Num_wann, Num_wann))
     allocate( UU(Num_wann, Num_wann))
     allocate( UU_dag(Num_wann, Num_wann))
     allocate( DHDk    (Num_wann, Num_wann, 3))
     allocate( DHDkdag (Num_wann, Num_wann, 3))
     kslice=0d0
     kslice_shape=0d0


     ik =0
     do i= 1, nk1
        do j= 1, nk2
           ik =ik +1
           kslice(:, ik)= K3D_start+ K3D_vec1*(i-1)/dble(nk1-1)  &
                     + K3D_vec2*(j-1)/dble(nk2-1)  
           kslice_shape(:, ik)= kslice(1, ik)* Kua+ kslice(2, ik)* Kub+ kslice(3, ik)* Kuc 
        enddo
     enddo

     !> get R coordinates 
     do iR=1, Nrpts
        crvec(:, iR)= Rua*irvec(1,iR) + Rub*irvec(2,iR) + Ruc*irvec(3,iR)
     enddo

     do ik= 1+ cpuid, knv2, num_cpu
        if (cpuid==0) print *, ik, knv2

        !> diagonalize hamiltonian
        k= kslice(:, ik)

        ! calculation bulk hamiltonian
        UU= 0d0
        call ham_bulk(k, Hamk_bulk)
        

        !> diagonalization by call zheev in lapack
        W= 0d0
       !call eigensystem_c( 'V', 'U', Num_wann, UU, W)
        call utility_diagonalize(hamk_bulk,Num_wann, W, UU)

        UU_dag= conjg(transpose(UU))

        vx= 0d0
        vy= 0d0
        vz= 0d0
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
              if (W(n) > 0d0 .and. W(m)<0d0) then
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
              if (W(m) > 0d0 .and. W(n)< 0d0) then
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

        call Im_trace(Num_wann, vx, Omega(1, ik))
        call Im_trace(Num_wann, vy, Omega(2, ik))
        call Im_trace(Num_wann, vz, Omega(3, ik))

     enddo ! ik

     Omega_mpi= 0d0
     call mpi_allreduce(Omega,Omega_mpi,size(Omega_mpi),&
                       mpi_dc,mpi_sum,mpi_cmw,ierr)

     if (cpuid==0) then
        open(unit=18, file='berry.dat')
        ik= 0
        do i= 1, nk1
           do j= 1, nk2
              ik= ik+ 1
              write(18, '(20f18.10)')kslice_shape(:, ik), Omega_mpi(:, ik)
           enddo
           write(18, *) ' '
        enddo

        close(18)

     endif


     return

  end subroutine berry_curvarture

  subroutine Fourier_R_to_k(k, ham)
     use para
     implicit none

     real(dp), intent(in) :: k(3)
     complex(dp), intent(out) :: ham(Num_wann, Num_wann)
     integer :: iR
     real(dp) :: R(3)
     real(dp) :: kdotr

     ham= 0d0
     do iR= 1, Nrpts
        R= Rua*irvec(1,iR) + Rub*irvec(2,iR) + Ruc*irvec(3,iR)
        kdotr= k(1)*irvec(1,iR) + k(2)*irvec(2,iR) + k(3)*irvec(3,iR)
        Ham= Ham+ HmnR(:,:,iR)*Exp(2d0*pi*zi*kdotr)/ndegen(iR)
     enddo

     return
  end subroutine Fourier_R_to_k

  subroutine Im_trace(ndim, A, tr)
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
