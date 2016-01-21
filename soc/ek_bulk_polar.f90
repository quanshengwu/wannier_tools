! calculate bulk's energy band using wannier TB method
  subroutine ek_bulk_polar

     use mpi
     use para

     implicit none

     integer :: ik, i, j
     integer :: nkx, nky
	  integer :: knv3
     integer :: ierr
     integer :: nktheta
     integer :: nkr
     real(Dp) :: k (3)
     real(Dp) :: k1(3)
     real(Dp) :: k2(3)
     real(dp) :: kx
     real(dp) :: ky
     real(dp) :: kz
     real(dp) :: ktheta
     real(dp) :: kr
     real(Dp) :: W(Num_wann)
     real(dp) :: krmax
     real(dp), allocatable :: k12(:,:)
     real(dp), allocatable :: k123(:,:)
     real(dp), allocatable :: krtheta(:,:)

     real(dp) :: mingap

     !> center
     real(dp) :: center_direct(3)
     real(dp) :: center_cart(3)
     
     ! Hamiltonian of bulk system
     complex(Dp) :: Hamk_bulk(Num_wann,Num_wann) 

     ! eigen value of H
	  real(dp), allocatable :: gap(:,:)
	  real(dp), allocatable :: eigv(:,:)
	  real(dp), allocatable :: eigv_mpi(:,:)


     nktheta= 100 
     nkr= Nk
     knv3= nktheta*nkr
     allocate( k12(2, knv3))
     allocate( k123(3, knv3))
     allocate( krtheta(2, knv3))
     allocate( gap    (4, nktheta))
     allocate( eigv    (4, knv3))
     allocate( eigv_mpi(4, knv3))
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
    !call direct_cart(center_direct, center_cart)
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
    !call direct_cart(center_direct, center_cart)
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
     center_direct= (/ 0.0d0,  0.0d0, -0.0d0/)
     call direct_cart(center_direct, center_cart)
     ik = 0
     do i= 1, nktheta
        ktheta= 2d0*pi*(i-1)/dble(nktheta)
        do j= 1, nkr
           kr= krmax*(j-1)/dble(nkr)
           ik =ik +1
           k12(1, ik)= kr*cos(ktheta)+ center_cart(1)
           k12(2, ik)= kr*sin(ktheta)+ center_cart(2)
           krtheta(1, ik)= ktheta
           krtheta(2, ik)= kr
           k123(1, ik) = k12(1, ik)
           k123(2, ik) = k12(2, ik)
           k123(3, ik)= center_cart(3)
        enddo
     enddo


     !> xz plane, center is (0.5, 0.0, 0.5)
    !center_direct= (/ 0.5d0,  0.0d0,  0.5d0/)
    !center_direct= (/ 0.5d0,  0.0d0, -0.5d0/)
    !call direct_cart(center_direct, center_cart)
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
    !call direct_cart(center_direct, center_cart)
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



     if (Numoccupied> Num_wann ) stop ' Numoccupied should less than Num_wann'

     do ik= 1+cpuid, knv3, num_cpu
	     if (cpuid==0) print * , ik

        k1= k123(:, ik)
        !> from cartisen coordinate to direct coordinate
        call cart_direct(k1, k)

        ! calculation bulk hamiltonian
        Hamk_bulk= 0d0
        call ham_bulk(k, Hamk_bulk)

        !> diagonalization by call zheev in lapack
        W= 0d0
        call eigensystem_c( 'N', 'U', Num_wann ,Hamk_bulk, W)

        eigv(:, ik)= W(Numoccupied-1:Numoccupied+2)

     enddo

     call mpi_allreduce(eigv,eigv_mpi,size(eigv),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)

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

     if (cpuid==0)then
        open(unit=14, file='bulkek_polar.dat')
        do i=1, nktheta
           write(14, '(1000f19.9)')gap(:, i)
        enddo
        do i=1, knv3
        !  write(14, '(1000f19.9)')k123(:, i)
        enddo
        close(14)
     endif

   return
   end subroutine ek_bulk_polar

   subroutine cart_direct(k1, k2)
      use para
      implicit none
      real(dp), intent(in) :: k1(3)
      real(dp), intent(inout) :: k2(3)
      real(dp) :: mata(3, 3)

      mata(1, :)= Kua
      mata(2, :)= Kub
      mata(3, :)= Kuc

      call inv_r(3, mata)
      K2= k1(1)*mata(1, :)+ k1(2)*mata(2, :)+ k1(3)*mata(3, :)

      return
   end subroutine cart_direct

   subroutine direct_cart(k1, k2)
      use para
      implicit none
      real(dp), intent(in) :: k1(3)
      real(dp), intent(inout) :: k2(3)

      K2= k1(1)*Kua+ k1(2)*Kub+ k1(3)*Kuc

      return
   end subroutine direct_cart
