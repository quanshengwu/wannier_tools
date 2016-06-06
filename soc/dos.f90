subroutine dos_calc

   use wmpi
   use para
   implicit none

   !> the integration k space
   real(dp) :: kxmin
   real(dp) :: kxmax
   real(dp) :: kymin
   real(dp) :: kymax
   real(dp) :: kzmin
   real(dp) :: kzmax
   real(dp) :: emin
   real(dp) :: emax

   integer :: ik
   integer :: ie
   integer :: ib
   integer :: ikx
   integer :: iky
   integer :: ikz
   integer :: knv3
   integer :: NE
   integer :: ierr

   !> integration for band
   integer :: iband_low
   integer :: iband_high

   real(dp) :: x
   real(dp) :: dkx
   real(dp) :: dky
   real(dp) :: dkz

   real(dp) :: k(3)

   real(dp), allocatable :: kpoints(:, :)
   real(dp), allocatable :: eigval(:, :)
   real(dp), allocatable :: eigval_mpi(:, :)
   real(dp), allocatable :: W(:)
   real(dp), allocatable :: omega(:)
   real(dp), allocatable :: dos(:)
   complex(dp), allocatable :: Hk(:, :)

   !> delta function
   real(dp), external :: delta

   knv3= Nk*Nk*Nk
   NE= 100

   allocate(W(Num_wann))
   allocate(Hk(Num_wann, Num_wann))
   allocate(eigval(Num_wann, knv3))
   allocate(eigval_mpi(Num_wann, knv3))
   allocate(kpoints(3, knv3))
   allocate(dos(NE))
   allocate(omega(NE))
   kpoints= 0d0
   eigval= 0d0
   eigval_mpi= 0d0

   iband_low= 56
   iband_high= 57
   !> around (0.1225, 0.012)
   kxmin= 0.12d0
   kxmax= 0.125d0
   kymin= 0.005d0
   kymax= 0.018d0
   kzmin= -0.010d0
   kzmax=  0.01d0

   !> around (0.1140, 0.038)
  !kxmin= 0.11d0
  !kxmax= 0.120d0
  !kymin= 0.032d0
  !kymax= 0.045d0
  !kzmin= -0.010d0
  !kzmax=  0.01d0

   !> around (0.1140, 0.038)
  !kxmin= 0.11d0
  !kxmax= 0.125d0
  !kymin= 0.005d0
  !kymax= 0.045d0
  !kzmin= -0.010d0
  !kzmax=  0.01d0

   ik= 0
   do ikx= 1, Nk
   do iky= 1, Nk
   do ikz= 1, Nk
      ik= ik+ 1
      kpoints(1, ik)= kxmin+ (kxmax-kxmin)*(ikx-1d0)/dble(Nk)
      kpoints(2, ik)= kymin+ (kymax-kymin)*(iky-1d0)/dble(Nk)
      kpoints(3, ik)= kzmin+ (kzmax-kzmin)*(ikz-1d0)/dble(Nk)
   enddo
   enddo
   enddo
   dkx= (kxmax-kxmin)/dble(Nk)
   dky= (kymax-kymin)/dble(Nk)
   dkz= (kzmax-kzmin)/dble(Nk)

   !> get eigenvalue
   do ik=1+cpuid, knv3, num_cpu
      if (cpuid.eq.0) print *, 'ik, knv3', ik, knv3
      k= kpoints(:, ik)
      call ham_bulk(k, Hk)
      W= 0d0
      call eigensystem_c( 'N', 'U', Num_wann ,Hk, W)
      eigval_mpi(:, ik)= W
   enddo
   call mpi_allreduce(eigval_mpi,eigval,size(eigval),&
                      mpi_dp,mpi_sum,mpi_cmw,ierr)

   emin= minval(eigval(iband_low,:))
   emax= maxval(eigval(iband_high,:))
   eta= (emax- emin)/ dble(NE)*3d0


   !> energy
   do ie=1, NE
      omega(ie)= emin+ (emax-emin)* (ie-1d0)/dble(NE-1)
   enddo ! ie

   !> get density of state
   dos= 0d0
   do ie= 1, NE
      if (cpuid.eq.0) print *, 'ie, NE', ie, NE
      !> intergrate with k
      do ik= 1, knv3
         do ib= iband_low, iband_high
            x= omega(ie)- eigval(ib, ik)
            dos(ie) = dos(ie)+ delta(eta, x)
         enddo ! ib
      enddo ! ik
      dos(ie)= dos(ie)*dkx*dky*dkz
   enddo ! ie

   if (cpuid.eq.0) then
      open(unit=100, file='dos.dat')
      do ie=1, NE
         write(100, *)omega(ie), dos(ie)
      enddo ! ie 
      close(100)
   endif

 
   return
end subroutine dos_calc

function delta(eta, x)
   implicit none
   integer, parameter :: dp=kind(1d0)
   integer, parameter :: pi= 3.1415926535d0
   real(dp), intent(in) :: eta
   real(dp), intent(in) :: x
   real(dp) :: delta

   delta= 1d0/pi*eta/(eta*eta+x*x)
  !delta= exp(-x*x/eta/eta/2d0)/sqrt(2d0*pi)/eta

end function


