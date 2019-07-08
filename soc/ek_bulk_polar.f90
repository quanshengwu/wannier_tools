  subroutine ek_bulk_polar
     ! Calculate bulk's energy band using wannier TB method. 
     ! This subroutine is very useful for nodal-line materials. 
     ! You have to define the center, k plane and radius
     ! Copyright (c) 2010 QuanSheng Wu. All rights reserved.

     use wmpi
     use para

     implicit none

     integer :: ik, i, j, ia, knv3
     integer :: ierr, nktheta, nkr, nwan
     real(Dp) :: k (3), k1(3)
     real(dp) :: kz, ktheta, kr, krmax, mingap
     real(Dp) :: W(Num_wann)
     real(dp), allocatable :: k12(:,:), k123(:,:), krtheta(:,:)

     !> center
     real(dp) :: center_direct(3)
     real(dp) :: center_cart(3)
     
     ! Hamiltonian of bulk system
     complex(Dp) :: Hamk_bulk(Num_wann,Num_wann) 

     ! eigen value of H
	  real(dp), allocatable :: gap(:,:)
	  real(dp), allocatable :: eigv(:,:)
	  real(dp), allocatable :: eigv_mpi(:,:)

     !> pauli matrix in basis
     complex(dp), allocatable :: sigmax(:, :)
     complex(dp), allocatable :: sigmay(:, :)
     complex(dp), allocatable :: sigmaz(:, :)
     complex(dp), allocatable :: mat1(:, :)
     complex(dp), allocatable :: eigvec_conj(:, :)

     !> spin expect value for each band and each k point
     real(dp), allocatable :: sx(:, :)
     real(dp), allocatable :: sy(:, :)
     real(dp), allocatable :: sz(:, :)


     nktheta= 100 
     nkr= Nk
     knv3= nktheta*nkr
     allocate( k12(2, knv3))
     allocate( k123(3, knv3))
     allocate( krtheta(2, knv3))
     allocate( gap    (4, nktheta))
     allocate( eigv    (4, knv3))
     allocate( eigv_mpi(4, knv3))
     allocate( eigvec_conj(Num_wann, Num_wann))
     allocate( sigmax(Num_wann, Num_wann))
     allocate( sigmay(Num_wann, Num_wann))
     allocate( sigmaz(Num_wann, Num_wann))
     allocate( mat1(Num_wann, Num_wann))
     allocate( sx(Num_wann, knv3))
     allocate( sy(Num_wann, knv3))
     allocate( sz(Num_wann, knv3))
     k12 = 0d0
     k123= 0d0
     krtheta = 0d0
     eigv    = 0d0
     eigv_mpi= 0d0

     krmax= 0.6d0
     kz= 0d0
     ik =0

     !> for IrF  DFT
     !> xz plane, center is (0.0, 0.5, 0.5)
     !> for IrF  TB
     !> xz plane, center is (0.5, 0.0, 0.5)
    !center_direct= (/ 0.5d0,  0.0d0,  0.5d0/)
    !call direct_cart_rec(center_direct, center_cart)
    !ik = 0
    !do i= 1, nktheta
    !   ktheta= 2d0*pi*(i-1)/dble(nktheta)
    !   do j= 1, nkr
    !      kr= krmax*(j-1)/dble(nkr)
    !      ik =ik +1
    !      k12(1, ik)= kr*cos(ktheta)+ center_cart(1)
    !      k12(2, ik)= kr*sin(ktheta)+ center_cart(3)
    !      krtheta(1, ik)= ktheta
    !      krtheta(2, ik)= kr
    !      k123(1, ik) = k12(1, ik)
    !      k123(2, ik)= center_cart(2)
    !      k123(3, ik) = k12(2, ik)
    !   enddo
    !enddo

     !> For DFT yz plane, center is (0.5, 0.0, 0.5)
    !center_direct= (/ 0.5d0, 0.0d0, 0.5d0/)
    !call direct_cart_rec(center_direct, center_cart)
    !ik = 0
    !do i= 1, nktheta
    !   ktheta= 2d0*pi*(i-1)/dble(nktheta)
    !   do j= 1, nkr
    !      kr= krmax*(j-1)/dble(nkr)
    !      ik =ik +1
    !      k12(1, ik)= kr*cos(ktheta)+ center_cart(2)
    !      k12(2, ik)= kr*sin(ktheta)+ center_cart(3)
    !      krtheta(1, ik)= ktheta
    !      krtheta(2, ik)= kr
    !      k123(1, ik)= center_cart(1)
    !      k123(2, ik) = k12(1, ik)
    !      k123(3, ik) = k12(2, ik)
    !   enddo
    !enddo

     !>> for TB
     !> xy plane, center is (0.0, 0.0, 0.0)
    !center_direct= (/ 0.0d0,  0.0d0, -0.0d0/)
    !call direct_cart_rec(center_direct, center_cart)
    !ik = 0
    !do i= 1, nktheta
    !   ktheta= 2d0*pi*(i-1)/dble(nktheta)
    !   do j= 1, nkr
    !      kr= krmax*(j-1)/dble(nkr)
    !      ik =ik +1
    !      k12(1, ik)= kr*cos(ktheta)+ center_cart(1)
    !      k12(2, ik)= kr*sin(ktheta)+ center_cart(2)
    !      krtheta(1, ik)= ktheta
    !      krtheta(2, ik)= kr
    !      k123(1, ik) = k12(1, ik)
    !      k123(2, ik) = k12(2, ik)
    !      k123(3, ik)= center_cart(3)
    !   enddo
    !enddo


     !> xz plane, center is (0.5, 0.0, 0.5)
    !center_direct= (/ 0.5d0,  0.0d0,  0.5d0/)
    !center_direct= (/ 0.5d0,  0.0d0, -0.5d0/)
    !call direct_cart_rec(center_direct, center_cart)
    !ik = 0
    !do i= 1, nktheta
    !   ktheta= 2d0*pi*(i-1)/dble(nktheta)
    !   do j= 1, nkr
    !      kr= krmax*(j-1)/dble(nkr)
    !      ik =ik +1
    !      k12(1, ik)= kr*cos(ktheta)+ center_cart(1)
    !      k12(2, ik)= kr*sin(ktheta)+ center_cart(3)
    !      krtheta(1, ik)= ktheta
    !      krtheta(2, ik)= kr
    !      k123(1, ik) = k12(1, ik)
    !      k123(2, ik)= center_cart(2)
    !      k123(3, ik) = k12(2, ik)
    !   enddo
    !enddo


     !> For yz plane, center is (0.0, 0.5, 0.5)
    !center_direct= (/ 0.0d0,-0.5d0,-0.5d0/)
    !center_direct= (/ 1.0d0, 0.5d0, 0.5d0/)
    !call direct_cart_rec(center_direct, center_cart)
    !ik = 0
    !do i= 1, nktheta
    !   ktheta= 2d0*pi*(i-1)/dble(nktheta)
    !   do j= 1, nkr
    !      kr= krmax*(j-1)/dble(nkr)
    !      ik =ik +1
    !      k12(1, ik)= kr*cos(ktheta)+ center_cart(2)
    !      k12(2, ik)= kr*sin(ktheta)+ center_cart(3)
    !      krtheta(1, ik)= ktheta
    !      krtheta(2, ik)= kr
    !      k123(1, ik)= center_cart(1)
    !      k123(2, ik) = k12(1, ik)
    !      k123(3, ik) = k12(2, ik)
    !   enddo
    !enddo


     !>> for CuTlTe2
     !> 110 plane, center is X (0.0, 0.0, 0.5)
     center_direct= (/ 0.0d0,  0.0d0,  0.5d0/)
     call direct_cart_rec(center_direct, center_cart)

     if (cpuid.eq.0) write(*, '(10f8.5)') center_direct
     if (cpuid.eq.0) write(*, '(10f8.5)') center_cart

     ik = 0
     krmax= 0.3d0
     do i= 1, nktheta
        ktheta= 2d0*pi*(i-1)/dble(nktheta)
        do j= 1, nkr
           kr= krmax*(j-1)/dble(nkr)
           ik =ik +1
           k123(1, ik)= kr*cos(ktheta)+ center_cart(1)
           k123(2, ik)= kr*cos(ktheta)+ center_cart(2)
           k123(3, ik)= kr*sin(ktheta)+ center_cart(3)
        enddo
     enddo

     !>> for CuTlTe2
     !> 1-10 plane, center is U (0.5, 0.5, 0.0)
    !center_direct= (/ 0.5d0,  0.5d0,  0.0d0/)
    !call direct_cart_rec(center_direct, center_cart)

    !if (cpuid.eq.0) write(*, '(10f8.5)') center_direct
    !if (cpuid.eq.0) write(*, '(10f8.5)') center_cart

    !ik = 0
    !krmax= 0.3d0
    !do i= 1, nktheta
    !   ktheta= 2d0*pi*(i-1)/dble(nktheta)
    !   do j= 1, nkr
    !      kr= krmax*(j-1)/dble(nkr)
    !      ik =ik +1
    !      k123(1, ik)= kr*cos(ktheta)+ center_cart(1)
    !      k123(2, ik)= -kr*cos(ktheta)+ center_cart(2)
    !      k123(3, ik)= kr*sin(ktheta)+ center_cart(3)
    !   enddo
    !enddo




     !> set Pauli matrix
     sigmax= 0d0
     sigmay= 0d0
     sigmaz= 0d0
     if (soc>0) then
        nwan= Num_wann/2
     else
        stop 'calculate spin texture need soc'
     endif

     i= 0
     do ia=1, Origin_cell%Num_atoms
        do j=1, Origin_cell%nprojs(ia)
           i= i+ 1
           sigmax(i, i+nwan) = 1d0
           sigmax(i+nwan, i) = 1d0
           sigmay(i, i+nwan) = -zi
           sigmay(i+nwan, i) = zi
           sigmaz(i, i) = 1d0
           sigmaz(i+nwan, i+nwan) = -1d0
        enddo  !> number of projectors
     enddo  ! ia  number of atoms


     if (Numoccupied> Num_wann ) stop ' Numoccupied should less than Num_wann'

     do ik= 1+cpuid, knv3, num_cpu
        if (cpuid==0) print * , 'ik' , ik

        k1= k123(:, ik)
        !> from cartisen coordinate to direct coordinate
        call cart_direct_rec(k1, k)

        ! calculation bulk hamiltonian
        Hamk_bulk= 0d0
        call ham_bulk_atomicgauge(k, Hamk_bulk)

        !> diagonalization by call zheev in lapack
        W= 0d0
        call eigensystem_c( 'V', 'U', Num_wann ,Hamk_bulk, W)
        eigvec_conj= conjg(transpose(Hamk_bulk))

        eigv(:, ik)= W(Numoccupied-1:Numoccupied+2)

        call mat_mul(Num_wann,eigvec_conj, sigmax, mat1)
        call mat_mul(Num_wann,mat1, Hamk_bulk, sigmax)

        call mat_mul(Num_wann,eigvec_conj, sigmay, mat1)
        call mat_mul(Num_wann,mat1, Hamk_bulk, sigmay)

        call mat_mul(Num_wann,eigvec_conj, sigmaz, mat1)
        call mat_mul(Num_wann,mat1, Hamk_bulk, sigmaz)

        do i=1, Num_wann
           sx(i, ik)= sigmax(i, i)
           sy(i, ik)= sigmay(i, i)
           sz(i, ik)= sigmaz(i, i)
        enddo
     enddo ! ik

#if defined (MPI)
     call mpi_allreduce(eigv,eigv_mpi,size(eigv),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
     eigv_mpi= eigv
#endif

     eigv= eigv_mpi
     ik= 0
     do i=1, nktheta
        mingap= 999d0
        do j= 1, nkr
           ik= ik+ 1
           if ((eigv(3, ik)- eigv(2, ik))< mingap) then
              mingap= eigv(3, ik)- eigv(2, ik)
              gap(1, i)= mingap
              gap(2:4, i)= k123(:, ik)
           endif
        enddo
     enddo

     outfileindex= outfileindex+ 1
     if (cpuid==0)then
        open(unit=outfileindex, file='bulkek_polar.dat')
        do i=1, nktheta
           write(outfileindex, '(1000f19.9)')gap(:, i)
        enddo
        do i=1, knv3
        !  write(outfileindex, '(1000f19.9)')k123(:, i)
        enddo
        close(outfileindex)
     endif

   return
   end subroutine ek_bulk_polar

