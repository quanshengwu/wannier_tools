
!> calculate density of state and joint density of state for 3D bulk system
!> JDOS(\omega)= \sum_k (f_c(k)-f_v(k) \delta(\omega- Ec(k)+ Ev(k))
subroutine dos_joint_dos

   use wmpi
   use para
   implicit none

   !> the integration k space
   real(dp) :: kxmin
   real(dp) :: kxmax
   real(dp) :: kymin
   real(dp) :: kymax
   real(dp) :: emin
   real(dp) :: emax

   integer :: ik
   integer :: ik_start
   integer :: ik_end
   integer :: Nk_local

   integer :: ie
   integer :: ib, ib1, ib2
   integer :: ikx
   integer :: iky
   integer :: NE
   integer :: ierr

   !> integration for band
   integer :: iband_low
   integer :: iband_high
   integer :: iband_tot

   real(dp) :: x
   real(dp) :: dk2

   real(dp) :: k(2)

   real(dp), allocatable :: W(:)
   real(dp), allocatable :: omega_dos(:)
   real(dp), allocatable :: omega_jdos(:)
   real(dp), allocatable :: dos(:)
   real(dp), allocatable :: dos_mpi(:)
   real(dp), allocatable :: jdos(:)
   real(dp), allocatable :: jdos_mpi(:)
   complex(dp), allocatable :: Hk(:, :)

   !> fermi distribution
   real(dp), allocatable :: fermi_dis(:)

   !> delta function
   real(dp), external :: delta

   NE= OmegaNum
   iband_low= Numoccupied- 40
   iband_high= Numoccupied+ 40

   if (iband_low <1) iband_low = 1
   if (iband_high >Num_wann) iband_high = Num_wann

   iband_tot= iband_high- iband_low+ 1


   allocate(dos(NE))
   allocate(dos_mpi(NE))
   allocate(jdos(NE))
   allocate(jdos_mpi(NE))
   allocate(omega_dos(NE))
   allocate(omega_jdos(NE))
   allocate(W(Num_wann))
   allocate(Hk(Num_wann, Num_wann))
   allocate(fermi_dis(Num_wann))
   W= 0d0
   Hk= 0d0
   fermi_dis= 0d0
   jdos= 0d0
   jdos_mpi= 0d0
   dos= 0d0
   dos_mpi= 0d0
   omega_dos= 0d0
   omega_jdos= 0d0
 

   dk2= 1d0/dble(Nk1*Nk2)

   emin= 0d0
   emax= OmegaMax
   eta= (emax- emin)/ dble(NE)*5d0

   !> energy
   do ie=1, NE
      omega_jdos(ie)= emin+ (emax-emin)* (ie-1d0)/dble(NE-1)
   enddo ! ie

   emin= OmegaMin
   emax= OmegaMax

   !> energy
   do ie=1, NE
      omega_dos(ie)= emin+ (emax-emin)* (ie-1d0)/dble(NE-1)
   enddo ! ie


   !> get eigenvalue
   dos_mpi= 0d0
   jdos_mpi= 0d0
   do ik=1+cpuid, Nk1*Nk2, num_cpu
      if (cpuid.eq.0) write(stdout, *) 'ik, knv2', ik, Nk1*Nk2 
      ikx= (ik- 1)/(Nk2)+ 1
      iky=  ik- (ikx-1)*Nk2 

      k= K2D_start+ K2D_vec1*(ikx-1)/dble(nk1-1)  &
                + K2D_vec2*(iky-1)/dble(nk2-1)  
      call ham_bulk_old(k, Hk)
      W= 0d0
      call eigensystem_c('N', 'U', Num_wann ,Hk, W)

      !> calculate fermi-dirac distribution
      do ib=iband_low, iband_high
         if (W(ib)<0) then
            fermi_dis(ib)= 1d0
         else
            fermi_dis(ib)= 0d0
         endif
      enddo !ib

      !> get density of state
      do ie= 1, NE
         do ib1= iband_low, iband_high-1
            do ib2= ib1+1, iband_high
               x= omega_jdos(ie)- W(ib2) + W(ib1)
               jdos_mpi(ie)= jdos_mpi(ie)+ delta(eta, x)* (fermi_dis(ib1)- fermi_dis(ib2))
            enddo ! ib2
         enddo ! ib1
      enddo ! ie

      !> get density of state
      do ie= 1, NE
         !> intergrate with k
         do ib= iband_low, iband_high-1
            x= omega_dos(ie)- W(ib)
            dos_mpi(ie) = dos_mpi(ie)+ delta(eta, x)
         enddo ! ib
      enddo ! ie

   enddo ! ik

   call mpi_allreduce(dos_mpi,dos,size(dos),&
                      mpi_dp,mpi_sum,mpi_cmw,ierr)
   call mpi_allreduce(jdos_mpi,jdos,size(jdos),&
                      mpi_dp,mpi_sum,mpi_cmw,ierr)

   if (cpuid.eq.0) then
      open(unit=100, file='jdos.dat')
      do ie=1, NE
         write(100, *)omega_jdos(ie), jdos(ie)*dk2
      enddo ! ie 
      close(100)
   endif

   if (cpuid.eq.0) then
      open(unit=100, file='dos.dat')
      do ie=1, NE
         write(100, *)omega_dos(ie), dos(ie)*dk2
      enddo ! ie 
      close(100)
   endif
 
   return
end subroutine dos_joint_dos


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


