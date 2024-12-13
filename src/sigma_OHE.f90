!---------------------------------------------------------!
!> Calculate magnetoresistance with R.G.Chambers's formula!
!> based on Boltzmann transport                           !
!> Written By QuanSheng Wu (wuquansheng@gmail.com)        !
!> Modified by Hanqi Pi (hqpi1999@gmail.com)              !       
!> Thanks Yi Liu for helpful discussions                  !
!>                                                        !
!> References :                                           !
!> [1] Electrons in metals and semiconductors,            !
!>     R.G. Chambers,                                     !
!> [2] Ab initio investigation of magnetic transport      !
!>     properties by Wannier interpolation, PHYSICAL      !
!>     REVIEW B 79, 245123 (2009), Yi Liu, Hai-Jun Zhang, !
!>     and Yugui Yao                                      !
!> [3] Magnetoresistance from Fermi surface topology,     !
!>     ShengNan Zhang, QuanSheng Wu, Yi Liu, and          !
!>     Oleg V. Yazyev,Phys. Rev. B 99, 035142 (2019)      !
!>                                                        !
!> Implemented on Oct. 07 2017                            !
!> uploaded on Sep. 05. 2019                              !
!---------------------------------------------------------!

subroutine sigma_ohe_calc_symm(mu_array, KBT_array, BTau_array, Nband_Fermi_Level, bands_fermi_level, sigma_ohe_tensor)
!-----------------------------------------------------------!
!> This is the main routine to calculate the                !  
!> magnetoresistance/conductivity based on Boltzmann        !
!> transport equation.                                      !
!>                                                          !
!> We calculate the conductivity and resistivity under      !
!> the band-resolved constantly relaxation time. We give    !
!> the simga/tau of Btau and rho*tau of Btau instead of     !
!> sigma and rho.                                           !
!>                                                          !
!> This version impose symmetry constrain to get the        !
!> transport coefficients under magnetic field.             !
!>                                                          !
!> It produces 4 types of output files,                     !
!> every type has OmegaNum files at different chemical      !
!>                                                          ! 
!> 1. sigma_bands_mu_***eV.dat : the band-resolved sigma/tau!
!> 2. rho_bands_mu_***eV.dat   : the band-resolved rho*tau  !
!> 3. sigma_total_mu_***eV.dat : the total sigma/tau        !
!> 4. rho_total_mu_***eV.dat   : the total rho*tau          !
!>                                                          !
!> The first two files are output on a grid of              !
!> (Nband_Fermi_Level, NumT, NBTau), while the last two     !
!> files are output on a grid of (NumT, NBTau). NumT is the !
!> number of temperature points, NBTau is the number of Btau! 
!> and Nband_Fermi_Level is the number of bands.            !
!>                                                          !
!> The output sigma/tau is in the unit of 1/(Ohm*m*s),      !
!> while the output rho*tau is in the unit of Ohm*m*s.      !
!-----------------------------------------------------------!
      use wmpi
      use para
      implicit none


      integer, intent(inout) :: Nband_Fermi_Level
      integer, intent(inout) :: bands_fermi_level(Nband_Fermi_Level)
      real(dp), intent(in) :: KBT_array(NumT) ! K_Boltzmann*Temperature in eV
      real(dp), intent(in) :: mu_array(OmegaNum) ! chemical potential relative to Fermi level in eV
      real(dp), intent(in) :: BTau_array(NBTau) ! omega*tau without units
      real(dp), intent(inout) :: sigma_ohe_tensor(9, NBTau, OmegaNum, NumT, Nband_Fermi_Level) 

      real(dp) :: coeff,  mu, BTau, KBT, exponent_max
      integer :: ie, ibtau, ikt

      integer :: Nk_total, Nk_current, Nk_start, Nk_end
      integer  :: knv3_left, knv3_left_mod, knv3, knv3_mod
      integer :: ik, iband, ik1, ik2, ik3, ik_temp
      integer :: ierr, it, i, ix, j1, j2, j
      integer :: nrecevs

      real(dp) :: v_t(3), vo_k, v_k(3)
      real(dp) :: k(3), k_start(3), magnetic_field(3)
      real(dp) :: sigma_symm_t(9), rho(3, 3)
      
      real(dp) :: time_start, time_end
      real(dp) :: time_start0, time_end0
      real(dp) :: time_start1, time_end1
      real(dp) :: time_rkf45, time_velocity, time_integral
      integer :: NSlice_Btau_inuse

      !> Btau slices for Runge-Kutta integration
      real(dp) :: Btau_start, Btau_final, DeltaBtau
      logical :: fail
      integer :: icycle, ishift


      !> energy bands 
      real(dp) :: EE
      real(dp), allocatable :: Ek(:)  ! in eV
      real(dp), allocatable :: Enk(:, :)

      !> minus fermi derivation
      real(dp) :: minusdfde

      !> 3-component velocity for each band and each k point 
      real(dp), allocatable :: velocity_k(:, :)
      real(dp), allocatable :: velocity_bar_k(:)
      real(dp), allocatable :: velocity(:, :, :)

      !> fileindex list for sigma_bands and rho_bands
      integer, allocatable :: sigmafileindex(:)
      integer, allocatable :: rhofileindex(:)
      
      !> 3-component velocity for each band and each k point 
      !> Eq.(3) in PRB 79, 245123(2009)
      real(dp), allocatable :: velocity_bar(:, :, :)

      type(kcube_type) :: KCube3D_total
      type(kcube_type) :: KCube3D_left(Nband_Fermi_Level)

      !> some arrays for mpi
      integer, allocatable :: info(:, :) !> howmany integers to be sent on each cpu 
      integer, allocatable :: Displs(:, :)

      !> number of steps used in the Runge-Kutta integration
      integer :: NSlice_Btau
      integer :: NSlice_Btau_local
      real(dp), allocatable :: kout(:, :), kout_all(:, :)
      real(dp), allocatable :: weights(:)

      !> define some arrays for different bands. Since there are different number
      !> of k points left for different bands.
      type klist_iband_type
         !> dim=3*NSlice_Btau
         real(dp), allocatable :: klist_rkfs(:, :)
         !> dim=3*NSlice_Btau
         real(dp), allocatable :: velocity_k(:, :)

         !> calculate -df(e)/de, where f(e) is the Fermi-Dirac distribution
         !> dim=(NumT, OmegaNum))
         real(dp), allocatable :: minusdfde(:, :)
      end type klist_iband_type
      type(klist_iband_type) :: klist_iband(Nband_Fermi_Level)

      type sigma_iband_type
           ! sigma_ohe_tensor_k(9, NBTau, OmegaNum, NumT) 
           real(dp), allocatable :: sigma_ohe_tensor_k_mpi(:, :, :, :)
           real(dp), allocatable :: sigma_ohe_tensor_k(:, :, :, :)
           ! sigma_ohe_tensor_kz(9, NBTau, OmegaNum, NumT, Nk3) 
         !   real(dp), allocatable :: sigma_ohe_tensor_kz(:, :, :, :, :)

           !> time cost for k point
           real(dp), allocatable :: time_cost(:)
           real(dp), allocatable :: time_cost_mpi(:)
      end type sigma_iband_type
      type(sigma_iband_type) :: sigma_iband_k(Nband_Fermi_Level)

      !> file index
      !integer, allocatable  :: myfileindex(:)
      character(80) :: sigmafilename, bandname, tname, muname, ibname, ikname, filename


      !> inverse of group operator
      real(dp) :: Tmat(3, 3)
      real(dp), allocatable :: pgop_cart_inverse(:, :, :)

      real(dp), external :: det3

      ! default value
      icycle = 1
      fail = .False.

      allocate(pgop_cart_inverse(3, 3, number_group_operators))
      do i= 1, number_group_operators
         Tmat= pgop_cart(:, :, i) 
         call inv_r(3, Tmat)
         pgop_cart_inverse(:, :, i)= Tmat
      enddo

      !> set up 

      !> Distribute k kpoints in the 3DCube into different MPI threads
      knv3= KCube3D_symm%Nk_total_symm
      KCube3D_total%Nk_total= knv3
      knv3_mod= mod(knv3, num_cpu)
      if (knv3_mod==0) then  !> perfect divided
         KCube3D_total%Nk_current= knv3/num_cpu
         KCube3D_total%Nk_start=1+ knv3*cpuid/num_cpu
         KCube3D_total%Nk_end  =(1+cpuid)*knv3/num_cpu
      else if (knv3/num_cpu==0) then    !> Number of MPI threads is large than knv3
         KCube3D_total%Nk_current= 1 !> one k piont per MPI thread
         KCube3D_total%Nk_start= cpuid+ 1 !> one k piont per MPI thread
         KCube3D_total%Nk_end  = cpuid+ 1
         if (cpuid+1 > knv3) then
            KCube3D_total%Nk_start= 1
            KCube3D_total%Nk_end  = 0
         endif
      else
         KCube3D_total%Nk_current= knv3/num_cpu+ 1
         if (cpuid< knv3_mod) then
            KCube3D_total%Nk_start= 1+ cpuid*KCube3D_total%Nk_current
            KCube3D_total%Nk_end  = (1+cpuid)*KCube3D_total%Nk_current
         else
            KCube3D_total%Nk_start= knv3_mod*KCube3D_total%Nk_current+ &
               (cpuid-knv3_mod)*(KCube3D_total%Nk_current-1)+1
            KCube3D_total%Nk_end  = knv3_mod*KCube3D_total%Nk_current+ &
               (cpuid-knv3_mod+1)*(KCube3D_total%Nk_current-1)
         endif
      endif

      !> calculate the volume of the k cube

      allocate(KCube3D_total%k_direct(3, KCube3D_total%Nk_start:KCube3D_total%Nk_end))
      allocate(KCube3D_total%weight_k(KCube3D_total%Nk_start:KCube3D_total%Nk_end))

      do ik= KCube3D_total%Nk_start, KCube3D_total%Nk_end
         j1= KCube3D_symm%ik_array_symm(ik)
         ik1= (j1-1)/(Nk2*Nk3)+1
         ik2= ((j1-1-(ik1-1)*Nk2*Nk3)/Nk3)+1
         ik3= (j1-(ik2-1)*Nk3- (ik1-1)*Nk2*Nk3)
         KCube3D_total%k_direct(1, ik)=  (ik1-1d0)/dble(Nk1)  
         KCube3D_total%k_direct(2, ik)=  (ik2-1d0)/dble(Nk2)  
         KCube3D_total%k_direct(3, ik)=  (ik3-1d0)/dble(Nk3)  
         KCube3D_total%weight_k(ik)= KCube3D_symm%weight_k(ik)
      enddo

#if defined (MPI)
     call mpi_barrier(mpi_cmw, ierr)
#endif

      !> setup NSlice_Btau
      !> NSlice_Btau should be the integer times of NBTau
      if (NBTau>1) then
         NSlice_Btau= (NBTau-1)*(Nslice_BTau_Max/(NBTau-1))
      else
         NSlice_Btau= 1
      endif

      if (cpuid.eq.0) write(stdout, *) ' NSlice_Btau :', NSlice_Btau

      Nk_total= KCube3D_total%Nk_total
      Nk_current= KCube3D_total%Nk_current
      Nk_start= KCube3D_total%Nk_start
      Nk_end= KCube3D_total%Nk_end

      allocate( weights(Nslice_BTau_Max))
      allocate( Ek(Nband_Fermi_Level))
      allocate( Enk(Nk_start:Nk_end, Nband_Fermi_Level))
      allocate( velocity(3, Nk_start:Nk_end, Nband_Fermi_Level))
      allocate( velocity_k(3, Nband_Fermi_Level))
      allocate( velocity_bar(3, Nk_start:Nk_end, Nband_Fermi_Level))
      allocate( velocity_bar_k(3))
      Ek= 0d0
      Enk= 0d0
      velocity= 0d0
      velocity_k= 0d0
      velocity_bar= 0d0
      velocity_bar_k= 0d0
      weights = 0d0

      time_start= 0d0
      time_end= 0d0
      do ik= Nk_start, Nk_end
         if (cpuid.eq.0.and. mod(ik, 100).eq.0) &
            write(stdout, '(a, i18, "   /", i18, a, f10.3, "s")') 'ik/NK', &
            ik-Nk_start,Nk_current, 'time left', &
            (Nk_current+Nk_start-ik)*(time_end-time_start)/num_cpu

         call now(time_start)
         k= KCube3D_total%k_direct(:, ik)

         call velocity_calc(Nband_Fermi_Level, bands_fermi_level, k, velocity_k, Ek)
         velocity(:, ik, :)= velocity_k
         Enk(ik, :)= Ek
         call now(time_end)
      enddo

      !> we get the kpath by Btau_final=-exponent_max*BTauMax, but we only use half of them
      !>  means that we can reach the accuracy as to exp(-exponent_max)
      ! exponent_max= 30
      exponent_max= 15

      !> exclude all kpoints with zero velocity x B and large energy away from Fermi level
      do iband= 1, Nband_Fermi_Level 

         !> first check howmany k points left
         it= 0
         do ik= Nk_start, Nk_end
            !> check whether v x B=0, which means the magnetic field is parallel with velocity
            v_t= velocity(:, ik, iband)
           !vcrossB(1)= -v_t(2)*Bdirection(3)+ v_t(3)*Bdirection(2)
           !vcrossB(2)= -v_t(3)*Bdirection(1)+ v_t(1)*Bdirection(3)
           !vcrossB(3)= -v_t(1)*Bdirection(2)+ v_t(2)*Bdirection(1)
           !if (abs(Enk(ik, iband))<0.05d0.and.dsqrt(sum((abs(vcrossB)**2)))>eps3) then
            if (abs(Enk(ik, iband))/eV2Hartree<EF_integral_range) then
               it = it+ 1
            endif
         enddo ! ik

         KCube3D_left(iband)%Nk_current= it
         if (it>0) then
            allocate(KCube3D_left(iband)%ik_array(it))
            allocate(KCube3D_left(iband)%Ek_local(it))
            allocate(KCube3D_left(iband)%vx_local(it))
            allocate(KCube3D_left(iband)%vy_local(it))
            allocate(KCube3D_left(iband)%vz_local(it))
            allocate(KCube3D_left(iband)%weight_k_local(it))
         else 
            allocate(KCube3D_left(iband)%ik_array(1))  !> only useful for mpi_allgatherv
            allocate(KCube3D_left(iband)%Ek_local(1))
            allocate(KCube3D_left(iband)%vx_local(1))
            allocate(KCube3D_left(iband)%vy_local(1))
            allocate(KCube3D_left(iband)%vz_local(1))
            allocate(KCube3D_left(iband)%weight_k_local(1))
         endif
         KCube3D_left(iband)%ik_array= 0
         KCube3D_left(iband)%Ek_local= 0d0
         KCube3D_left(iband)%vx_local= 0d0
         KCube3D_left(iband)%vy_local= 0d0
         KCube3D_left(iband)%vz_local= 0d0
         KCube3D_left(iband)%weight_k_local= 0d0


         it= 0
         do ik= Nk_start, Nk_end
            !> check whether v x B=0, which means the magnetic field is parallel with velocity
            v_t= velocity(:, ik, iband)
           !vcrossB(1)= -v_t(2)*Bdirection(3)+ v_t(3)*Bdirection(2)
           !vcrossB(2)= -v_t(3)*Bdirection(1)+ v_t(1)*Bdirection(3)
           !vcrossB(3)= -v_t(1)*Bdirection(2)+ v_t(2)*Bdirection(1)
           !if (abs(Enk(ik, iband))<0.05d0.and.dsqrt(sum((abs(vcrossB)**2)))>eps3) then
            if (abs(Enk(ik, iband))/eV2Hartree<EF_integral_range) then
               it = it+ 1
               KCube3D_left(iband)%weight_k_local(it) = KCube3D_symm%weight_k(ik)
               KCube3D_left(iband)%ik_array(it) = KCube3D_symm%ik_array_symm(ik)
               KCube3D_left(iband)%Ek_local(it) = Enk(ik, iband)
               KCube3D_left(iband)%vx_local(it) = v_t(1)
               KCube3D_left(iband)%vy_local(it) = v_t(2)
               KCube3D_left(iband)%vz_local(it) = v_t(3)
            endif
         enddo ! ik
      enddo ! iband

      !> try to get the total number of k points left for each band 
      do iband=1, Nband_Fermi_Level
#if defined (MPI)
         call mpi_allreduce(KCube3D_left(iband)%Nk_current,KCube3D_left(iband)%Nk_total,1,&
                      mpi_in,mpi_sum,mpi_cmw,ierr)
#else
         KCube3D_left(iband)%Nk_total= KCube3D_left(iband)%Nk_current
#endif
      enddo

      !> gather the number of k points left into a array info
      allocate(info(num_cpu, Nband_Fermi_Level))
      info= -1
      do iband= 1, Nband_Fermi_Level
         nrecevs= KCube3D_left(iband)%Nk_current
         if (nrecevs<0) nrecevs= 0
#if defined (MPI)
         call mpi_allgather(nrecevs, 1, mpi_in, info(:, iband), 1, mpi_in, mpi_cmw, ierr)
#else
         info(1, iband)= KCube3D_left(iband)%Nk_current
#endif
      enddo


      !> An array for mpi_allgatherv
      allocate(Displs(num_cpu+1, Nband_Fermi_Level))
      Displs= 0
      do iband=1, Nband_Fermi_Level
         Displs(1, iband)=0
         do i=2, num_cpu+1
            Displs(i, iband)=Displs(i-1, iband)+ info(i-1, iband)
         enddo
      enddo ! iband
#if defined (MPI)
     call mpi_barrier(mpi_cmw, ierr)
#endif

      !> put all the kpoints left together
      do iband=1, Nband_Fermi_Level
         allocate(KCube3D_left(iband)%IKleft_array(KCube3D_left(iband)%Nk_total))
         allocate(KCube3D_left(iband)%Ek_total(KCube3D_left(iband)%Nk_total))
         allocate(KCube3D_left(iband)%vx_total(KCube3D_left(iband)%Nk_total))
         allocate(KCube3D_left(iband)%vy_total(KCube3D_left(iband)%Nk_total))
         allocate(KCube3D_left(iband)%vz_total(KCube3D_left(iband)%Nk_total))
         allocate(KCube3D_left(iband)%weight_k(KCube3D_left(iband)%Nk_total))
         KCube3D_left(iband)%IKleft_array = 0
         KCube3D_left(iband)%Ek_total= 0d0
         KCube3D_left(iband)%vx_total= 0d0
         KCube3D_left(iband)%vy_total= 0d0
         KCube3D_left(iband)%vz_total= 0d0
         KCube3D_left(iband)%weight_k= 0d0
      enddo  ! iband
     !> gather Enk and velocity 
#if defined (MPI)
      do iband=1, Nband_Fermi_Level
        !nrecevs= KCube3D_left(iband)%Nk_end-KCube3D_left(iband)%Nk_start+ 1
         nrecevs= KCube3D_left(iband)%Nk_current
         if (nrecevs<0) nrecevs= 0
         call mpi_allgatherv(KCube3D_left(iband)%ik_array, nrecevs, &
                             mpi_in, KCube3D_left(iband)%IKleft_array, &
                             info(:, iband), Displs(:, iband), mpi_in, mpi_cmw, ierr)
         call mpi_allgatherv(KCube3D_left(iband)%Ek_local, nrecevs, &
                             mpi_dp, KCube3D_left(iband)%Ek_total, &
                             info(:, iband), Displs(:, iband), mpi_dp, mpi_cmw, ierr)
         call mpi_allgatherv(KCube3D_left(iband)%vx_local, nrecevs, &
                             mpi_dp, KCube3D_left(iband)%vx_total, &
                             info(:, iband), Displs(:, iband), mpi_dp, mpi_cmw, ierr)
         call mpi_allgatherv(KCube3D_left(iband)%vy_local, nrecevs, &
                             mpi_dp, KCube3D_left(iband)%vy_total, &
                             info(:, iband), Displs(:, iband), mpi_dp, mpi_cmw, ierr)
         call mpi_allgatherv(KCube3D_left(iband)%vz_local, nrecevs, &
                             mpi_dp, KCube3D_left(iband)%vz_total, &
                             info(:, iband), Displs(:, iband), mpi_dp, mpi_cmw, ierr)
         call mpi_allgatherv(KCube3D_left(iband)%weight_k_local, nrecevs, &
                             mpi_dp, KCube3D_left(iband)%weight_k, &
                             info(:, iband), Displs(:, iband), mpi_dp, mpi_cmw, ierr)
      enddo ! iband
#else
      do iband=1, Nband_Fermi_Level
         KCube3D_left(iband)%IKleft_array= KCube3D_left(iband)%ik_array
         KCube3D_left(iband)%Ek_total= KCube3D_left(iband)%Ek_local
         KCube3D_left(iband)%vx_total= KCube3D_left(iband)%vx_local
         KCube3D_left(iband)%vy_total= KCube3D_left(iband)%vy_local
         KCube3D_left(iband)%vz_total= KCube3D_left(iband)%vz_local
         KCube3D_left(iband)%weight_k= KCube3D_left(iband)%weight_k_local
      enddo ! iband
#endif

      !> redistribute all those k points into different cpus
      if (cpuid.eq.0) then
         write(stdout, '(a)')' '
         write(stdout, '(a, i10, a)')' There are ', KCube3D_total%Nk_total, ' k points generated by the input files' 
         do iband= 1, Nband_Fermi_Level
            write(stdout, '(a, i10, a, i10)')' However there are only ', KCube3D_left(iband)%Nk_total, &
               ' k points contribute to the conductance calculations for band ', bands_fermi_level(iband)
         enddo
      endif


      !> redistribute the left kpoints into different processors
      !> for different bands, the number of left kpoints is different.
      do iband= 1, Nband_Fermi_Level
         knv3_left= KCube3D_left(iband)%Nk_total
         knv3_left_mod= mod(knv3_left, num_cpu)
         if (knv3_left_mod==0) then  !> perfect divided
            KCube3D_left(iband)%Nk_current= knv3_left/num_cpu
            KCube3D_left(iband)%Nk_start=1+ knv3_left*cpuid/num_cpu
            KCube3D_left(iband)%Nk_end  =(1+cpuid)*knv3_left/num_cpu
         else if (knv3_left/num_cpu==0) then    !> Number of MPI threads is large than knv3_left
            KCube3D_left(iband)%Nk_current= 1 !> one k piont per MPI thread
            KCube3D_left(iband)%Nk_start= cpuid+ 1 !> one k piont per MPI thread
            KCube3D_left(iband)%Nk_end  = cpuid+ 1
            if (cpuid+1 > knv3_left) then
               KCube3D_left(iband)%Nk_start= 1
               KCube3D_left(iband)%Nk_end  = 0
            endif
         else
            KCube3D_left(iband)%Nk_current= knv3_left/num_cpu+ 1
            if (cpuid< knv3_left_mod) then
               KCube3D_left(iband)%Nk_start= 1+ cpuid*KCube3D_left(iband)%Nk_current
               KCube3D_left(iband)%Nk_end  = (1+cpuid)*KCube3D_left(iband)%Nk_current
            else
               KCube3D_left(iband)%Nk_start= knv3_left_mod*KCube3D_left(iband)%Nk_current+ &
                  (cpuid-knv3_left_mod)*(KCube3D_left(iband)%Nk_current-1)+1
               KCube3D_left(iband)%Nk_end  = knv3_left_mod*KCube3D_left(iband)%Nk_current+ &
                  (cpuid-knv3_left_mod+1)*(KCube3D_left(iband)%Nk_current-1)
            endif
         endif


         allocate(KCube3D_left(iband)%k_direct(3, KCube3D_left(iband)%Nk_start:KCube3D_left(iband)%Nk_end))
   
         do ik= KCube3D_left(iband)%Nk_start, KCube3D_left(iband)%Nk_end
            i= KCube3D_left(iband)%IKleft_array(ik)
            ik1= (i-1)/(Nk2*Nk3)+1
            ik2= ((i-1-(ik1-1)*Nk2*Nk3)/Nk3)+1
            ik3= (i-(ik2-1)*Nk3- (ik1-1)*Nk2*Nk3)
            KCube3D_left(iband)%k_direct(1, ik)= (ik1-1)/dble(Nk1) 
            KCube3D_left(iband)%k_direct(2, ik)= (ik2-1)/dble(Nk2) 
            KCube3D_left(iband)%k_direct(3, ik)= (ik3-1)/dble(Nk3) 
         enddo ! ik
         
         !> allocate array to store the conductivity for eack k point  
         allocate( sigma_iband_k(iband)%sigma_ohe_tensor_k(9, NBTau, OmegaNum, NumT), stat= ierr)
         if (ierr>0) stop ' Error : no enough memory'
         allocate( sigma_iband_k(iband)%sigma_ohe_tensor_k_mpi(9, NBTau, OmegaNum, NumT), stat= ierr)
         if (ierr>0) stop ' Error : no enough memory'
         ! allocate( sigma_iband_k(iband)%sigma_ohe_tensor_kz(9, NBTau, OmegaNum, NumT, Nk3))
         allocate( sigma_iband_k(iband)%time_cost(KCube3D_left(iband)%Nk_total))
         allocate( sigma_iband_k(iband)%time_cost_mpi(KCube3D_left(iband)%Nk_total))
         sigma_iband_k(iband)%time_cost= 0d0
         sigma_iband_k(iband)%time_cost_mpi= 0d0
         sigma_iband_k(iband)%sigma_ohe_tensor_k= 0d0
         ! sigma_iband_k(iband)%sigma_ohe_tensor_kz= 0d0
         sigma_iband_k(iband)%sigma_ohe_tensor_k_mpi= 0d0
      enddo  ! iband=1, Nband_Fermi_Level


      !> gather the number of receive buffs left into a array info
      info= -1
      do iband= 1, Nband_Fermi_Level
         nrecevs= (KCube3D_left(iband)%Nk_end-KCube3D_left(iband)%Nk_start+1)*9*NBTau*OmegaNum*NumT
         if (nrecevs<0) nrecevs=0
#if defined (MPI)
         call mpi_allgather(nrecevs, 1, mpi_in, info(:, iband), 1, mpi_in, mpi_cmw, ierr)
#else
         info(1, iband)= nrecevs
#endif
      enddo

      if (cpuid.eq.0) write(stdout, '(a, 10i8)') ' '
      do iband=1, Nband_Fermi_Level
         if (cpuid.eq.0) write(stdout, '(a, i7)')'>> Number of k points at different CPUs at band', iband
         if (cpuid.eq.0) write(stdout, '(10i8)')info(:, iband)/9/NBTau/OmegaNum/NumT
      enddo
      if (cpuid.eq.0) write(stdout, '(a, 10i8)') ' '


      !> An array for mpi_allgatherv
      Displs= 0
      do iband=1, Nband_Fermi_Level
         Displs(1, iband)=0
         do i=2, num_cpu+1
            Displs(i, iband)=Displs(i-1, iband)+ info(i-1, iband)
         enddo
      enddo ! iband

      do iband=1, Nband_Fermi_Level
         allocate(klist_iband(iband)%klist_rkfs(3, NSlice_Btau))
         allocate(klist_iband(iband)%velocity_k(3, NSlice_Btau))
         klist_iband(iband)%klist_rkfs= 0d0
         klist_iband(iband)%velocity_k= 0d0
      enddo

      !> a temp array used in RKFS
      allocate(kout(3, NSlice_Btau))
      allocate(kout_all(3, NSlice_Btau))
      kout= 0d0
      kout_all= 0d0

      !> allocate array for sigmafileindex and rhofileindex
      allocate(sigmafileindex(OmegaNum))
      allocate(rhofileindex(OmegaNum))
      sigmafileindex= (/(outfileindex + ie, ie = 1, OmegaNum)/)
      rhofileindex= (/(outfileindex + ie +OmegaNum, ie = 1, OmegaNum)/)
      outfileindex = outfileindex + 2*OmegaNum

      if (cpuid.eq.0) then
         do ie = 1, OmegaNum
            write(muname, '(f12.1)')mu_array(ie)/eV2Hartree*1000

            !> open file for conductivity at mu and write header
            write(sigmafilename, '(3a)')'sigma_bands_mu_',trim(adjustl(muname)),'meV.dat'
            open(unit=sigmafileindex(ie), file=sigmafilename)
            write(sigmafileindex(ie), '(4a)')'# Conductivity tensor/tau (in unit of (\Omega*m*s)^-1) for every contributing band' , &
               ', at mu= ', trim(adjustl(muname)), ' meV'
            write(sigmafileindex(ie), '(a,100I5)')'# SELECTEDBANDS  =  ', bands_fermi_level(:)
            write(sigmafileindex(ie), '(a,1000f8.3)')'# Tlist  =  ', KBT_array(:)/8.6173324E-5/eV2Hartree
            write(sigmafileindex(ie), '(a, 1000E16.6)') '# Btaulist  =  ', BTau_array(:)*Magneticfluxdensity_atomic/Relaxation_Time_Tau

            !> open file for resistivity at mu and write header
            write(sigmafilename, '(3a)')'rho_bands_mu_', trim(adjustl(muname)),'meV.dat'
            open(unit=rhofileindex(ie), file=sigmafilename)
            write(rhofileindex(ie), '(4a)')'# Resistivity \tau*\rho (in unit of \Omega*m*s) for every contributing band ', &
               ', at mu= ', trim(adjustl(muname)), ' meV'
            write(rhofileindex(ie), '(a,100I5)')'# SELECTEDBANDS  =  ', bands_fermi_level(:)
            write(rhofileindex(ie), '(a,1000f8.3)')'# Tlist  =  ', KBT_array(:)/8.6173324E-5/eV2Hartree
            write(rhofileindex(ie), '(a, 1000E16.6)') '# Btaulist  =  ', BTau_array(:)*Magneticfluxdensity_atomic/Relaxation_Time_Tau
         
         enddo ! ie
      endif ! cpuid.eq.0

      !> now we turn to use Runge-Kutta method to get all the kpoints from (0, BTauMax)
      !> and we calculate the conductivity/Tau over different bands and different k points
      time_start= 0d0
      time_end  = 0d0

      do iband= 1, Nband_Fermi_Level
         call now(time_start0)
      
         !> dim=(Nk_start: Nk_end, NumT, OmegaNum))
         allocate(klist_iband(iband)%minusdfde(OmegaNum, NumT))
         call cal_sigma_iband_k

#if defined (MPI)
         call mpi_allreduce(sigma_iband_k(iband)%time_cost_mpi, sigma_iband_k(iband)%time_cost, &
                            size(sigma_iband_k(iband)%time_cost), &
                            mpi_dp,mpi_sum,mpi_cmw,ierr)
         call mpi_allreduce(sigma_iband_k(iband)%sigma_ohe_tensor_k_mpi, sigma_iband_k(iband)%sigma_ohe_tensor_k, &
                            size(sigma_iband_k(iband)%sigma_ohe_tensor_k), mpi_dp,mpi_sum,mpi_cmw,ierr)
         if (ierr>0) then
            write(stdout, *)'>>>Error happends in mpi_allgatherv at cpuid', cpuid, ' ierr=', ierr
            stop
         endif

#else
         sigma_iband_k(iband)%time_cost= sigma_iband_k(iband)%time_cost_mpi
         sigma_iband_k(iband)%sigma_ohe_tensor_k= sigma_iband_k(iband)%sigma_ohe_tensor_k_mpi
#endif

        !call mpi_barrier(mpi_cmw, ierr)
        !call now(time_2)
        !if (cpuid.eq.0) write(stdout, *) 'after mpi_allgatherv in sigma_ohe_calc', time_2-time_1, 's'

         if (cpuid.eq.0) then
            write(stdout, '(a)')' '
            write(stdout, '(a, i10)')'>> Time cost for each k point at iband ', iband
            write(stdout, '(10f16.2)')(sigma_iband_k(iband)%time_cost(ik), ik= 1, KCube3D_left(iband)%Nk_total)
            write(stdout, '(a)')' '
         endif
   
   
         !> sum over all the kpoints to get the band dependent conductivity/tau 
         do ikt=1, NumT
            do ie=1, OmegaNum
               do ibtau=1, NBTau
                  do i=1, 9
                     sigma_ohe_tensor(i, ibtau, ie, ikt, iband)=  &
                     sigma_iband_k(iband)%sigma_ohe_tensor_k(i, ibtau, ie, ikt)
                !print *, i, ibtau, ie, ikt, sigma_ohe_tensor(i, ibtau, ie, ikt, iband)
                !pause
                  enddo ! i
               enddo ! ibtau
            enddo ! ie
         enddo ! ikt
   
         call now(time_end0)
         if (cpuid.eq.0) write(stdout, '(a, i6, a, f16.2, a)')'>> Time cost for calculate MR at iband=', iband, &
         'is ', time_end0- time_start0, ' s'
   
         !> In the end, we start to care about the units of the conductivity/tau
         !> change from the Hatree Atomic units to SI units
         !> the conductivity/tau is in units of Ohm^-1*m^-1*s^-1
         !> 
         coeff= Echarge**2/hbar/Bohr_radius/Origin_cell%CellVolume/kCubeVolume*Origin_cell%ReciprocalCellVolume
         coeff= coeff/Time_atomic 
         sigma_ohe_tensor(:, :, :, :, iband)= sigma_ohe_tensor(:, :, :, :, iband)*coeff
         
         if (cpuid.eq.0) then
            do ie = 1, OmegaNum

               !> open file for conductivity at mu and write data
               write(sigmafileindex(ie), '(2a, i5)')'# ',' iband = ', bands_fermi_level(iband)
               write(sigmafileindex(ie),'(a)') ''
               do ikt = 1, NumT
                  KBT= KBT_array(ikt)/8.6173324E-5/eV2Hartree
                  write(sigmafileindex(ie), '(2a, f16.4, a)') '#', ' T = ', KBT, ' K'
                  write(sigmafileindex(ie), '("# Column", i5, 100i16)')(i, i=1, 10)
                  write(sigmafileindex(ie), '("#",19a16)')'BTau (T.ps)', 'xx', 'xy', 'xz', 'yx', 'yy', 'yz', 'zx', 'zy','zz'
                  !> write out the conductivity/tau into file
                  do i=1, NBTau
                     write(sigmafileindex(ie), '(19E16.6)') &
                        BTau_array(i)*Magneticfluxdensity_atomic/Relaxation_Time_Tau, &
                        sigma_ohe_tensor(:, i, ie, ikt, iband)
                  enddo ! i, NBTau
                  write(sigmafileindex(ie),'(a)') ''
               enddo
               write(sigmafileindex(ie),'(a)') ''

               !> open file for resistivity at mu and write data
               write(rhofileindex(ie), '(2a, i5)')'# ',' iband = ', bands_fermi_level(iband)
                  write(rhofileindex(ie),'(a)') ''
                  do ikt = 1, NumT
                     KBT= KBT_array(ikt)/8.6173324E-5/eV2Hartree
                     write(rhofileindex(ie), '(2a, f16.4, a)') '#', ' T = ', KBT, ' K'
                     write(rhofileindex(ie), '("# Column", i5, 100i16)')(i, i=1, 10)
                     write(rhofileindex(ie), '("#",19a16)')'BTau (T.ps)', 'xx', 'xy', 'xz', 'yx', 'yy', 'yz', 'zx', 'zy','zz'               
                     !> write out the inverse of conductivity/tau into file
                     do i=1, NBTau
                        rho(1, 1:3)=sigma_ohe_tensor(1:3, i, ie, ikt, iband)
                        rho(2, 1:3)=sigma_ohe_tensor(4:6, i, ie, ikt, iband)
                        rho(3, 1:3)=sigma_ohe_tensor(7:9, i, ie, ikt, iband)
                        if (abs(det3(rho))>eps6) then
                           call inv_r(3, rho)
                           write(rhofileindex(ie), '(19E16.6)')& 
                              BTau_array(i)*Magneticfluxdensity_atomic/Relaxation_Time_Tau, &
                              rho(1, 1:3), rho(2, 1:3), rho(3, 1:3)
                        else
                           write(rhofileindex(ie), '(a)')'# error: sigma is zero since no k points contribute to the calculations of MR'
                        endif
                     enddo ! i, NBTau
                     write(rhofileindex(ie),'(a)')''
                  enddo
               write(rhofileindex(ie),'(a)')''

            enddo ! ie, OmegaNum
         endif 

         !---------------------------------------------------------------------------------------------------------------
         !> This part is to calculate and store sigma_kz and was commentted out by pihanqi on Apr. 29. 2023 
         ! !calculate sigma_kz
         ! do ik3=1, Nk3
         !    do ik= 1, KCube3D_left(iband)%Nk_total
         !       i= KCube3D_left(iband)%IKleft_array(ik)
         !       ik1= (i-1)/(Nk2*Nk3)+1
         !       ik2= ((i-1-(ik1-1)*Nk2*Nk3)/Nk3)+1
         !       if (ik3== (i-(ik2-1)*Nk3- (ik1-1)*Nk2*Nk3)) then
         !          do  ikt=1, NumT
         !             do ie=1, OmegaNum
         !                do ibtau=1, NBTau
         !                   do ix=1, 9
         !                      sigma_iband_k(iband)%sigma_ohe_tensor_kz(ix, ibtau, ie, ikt, ik3)= &
         !                      sigma_iband_k(iband)%sigma_ohe_tensor_kz(ix, ibtau, ie, ikt, ik3)+ &
         !                      sigma_iband_k(iband)%sigma_ohe_tensor_k (ix, ibtau, ie, ikt, ik )   
         !                   enddo ! ix
         !                enddo ! ibtau
         !             enddo ! ie
         !          enddo ! ikt
         !       endif
         !    enddo ! ik
         ! enddo ! ik3

         ! sigma_iband_k(iband)%sigma_ohe_tensor_kz= &
         ! sigma_iband_k(iband)%sigma_ohe_tensor_kz*coeff
   
         !allocate(myfileindex(Nband_Fermi_Level))
         !!> file index for different bands
         !do iband=1, Nband_Fermi_Level
         !   outfileindex= outfileindex+ 1
         !   myfileindex(iband)= outfileindex
         !enddo
 
         ! do ikt=1, NumT
         !    do ie=1, OmegaNum
         !    outfileindex= outfileindex+ 1
         !    KBT= KBT_array(ikt)/8.6173324E-5/eV2Hartree
         !    write(tname, '(f12.1)')KBT
         !    write(muname, '(f12.1)')mu_array(ie)/eV2Hartree
         !    if (cpuid.eq.0) then
         !       write(bandname, '(i10)')bands_fermi_level(iband)
         !       write(sigmafilename, '(7a)')'sigma_kz_band_', trim(adjustl(bandname)),'_mu_',&
         !          trim(adjustl(muname)),'meV_T_', trim(adjustl(tname)), 'K.dat'
         !       open(unit=outfileindex, file=sigmafilename)
         !       write(outfileindex, '(a20, i5, 2(a, f16.4, a))')'# sigma_k at band ', & 
         !          bands_fermi_level(iband), ' temperature at ', KBT, ' K', ' chemical potential at ', mu_array(ie), ' eV'
         !       write(outfileindex, '("#",a12,20a16)')'ik3', 'xx', 'xy', 'xz', 'yx', 'yy', 'yz', 'zx', 'zy','zz'
         !    endif
   
         !    if (cpuid.eq.0) then
         !       do ibtau=1, NBTau
         !          !> write out the conductivity/tau into file
         !          write(outfileindex, '(a6, f16.8, a12, f16.8)')'# Btau', &
         !             BTau_array(ibtau)*Magneticfluxdensity_atomic/Relaxation_Time_Tau, &
         !             ' omega*tau', BTau_array(ibtau)*Magneticfluxdensity_atomic/Relaxation_Time_Tau*0.175874356d0
         !          write(outfileindex, '("#",a12,20a16)')'ik3', 'xx', 'xy', 'xz', 'yx', 'yy', 'yz', 'zx', 'zy','zz'
         !          do ik3=1, Nk3
         !             write(outfileindex, '(i6,20E16.6)')ik3, sigma_iband_k(iband)%sigma_ohe_tensor_kz(:, ibtau, ie, ikt, ik3)
         !          enddo
         !       enddo ! ibtau
         !    endif ! cpuid=0
            
         !    if (cpuid.eq.0) then
         !       close(outfileindex)
         !    endif ! cpuid=0
         !---------------------------------------------------------------------------------------------------------------
      
      enddo ! iband =1, Nband_Fermi_Level

      if (cpuid.eq.0) then
         do ie = 1, OmegaNum

            !> close file for conductivity at mu
            close(sigmafileindex(ie))
            !> close file for resistivity at mu
            close(rhofileindex(ie))
            
         enddo ! ie
      endif ! cpuid.eq.0

#if defined (MPI)
            call mpi_barrier(mpi_cmw, ierr)
#endif
      !> write data into file, modified by pihanqi on Apr. 29. 2023
      !> write out the conductivity/tau and resistivity*tau with the same relaxation time for every contributing bands
      if (cpuid.eq.0) then
         do ie=1, OmegaNum
            write(muname, '(f12.1)')mu_array(ie)/eV2Hartree*1000
            
            !> write out the total conductivity with the same relaxation time for all bands
            outfileindex= outfileindex+ 1
            write(sigmafilename, '(3a)')'sigma_total_mu_',trim(adjustl(muname)),'meV.dat'
            open(unit=outfileindex, file=sigmafilename)
            write(outfileindex, '(4a)')'# \sigma/\tau with unit (Ohm*m*s)^-1 is the summation of all bands \sum_n\sigma_n/\tau_n ', &
               ', at mu= ', trim(adjustl(muname)), ' meV'
            write(outfileindex, '(a)')'# relaxation time \tau_n=\tau is the same for all bands. '
            write(outfileindex, '(a,2I6)')'# NumT  NumBtau  =  ', NumT, NBTau 
            write(outfileindex, '(a,1000f8.3)')'# Tlist  =  ', KBT_array(:)/8.6173324E-5/eV2Hartree
            do ikt = 1, NumT
               KBT= KBT_array(ikt)/8.6173324E-5/eV2Hartree
               write(outfileindex, '(2a, f16.4, a)')'# ',' T = ', KBT, ' K '
               write(outfileindex, '("# Column", i5, 100i16)')(i, i=1, 10)
               write(outfileindex, '("#",19a16)')'BTau (T.ps)', 'xx', 'xy', 'xz', 'yx', 'yy', 'yz', 'zx', 'zy','zz'
               !> write out the inverse of conductivity/tau into file
               !> the name of rho is meaningless here, just a temp variable
               do i=1, NBTau
                  rho(1, 1:3)=sum(sigma_ohe_tensor(1:3, i, ie, ikt, :), dim=2)
                  rho(2, 1:3)=sum(sigma_ohe_tensor(4:6, i, ie, ikt, :), dim=2)
                  rho(3, 1:3)=sum(sigma_ohe_tensor(7:9, i, ie, ikt, :), dim=2)
                  write(outfileindex, '(19E16.6)')&
                     BTau_array(i)*Magneticfluxdensity_atomic/Relaxation_Time_Tau, &
                     rho(1, 1:3), rho(2, 1:3), rho(3, 1:3)
               enddo ! i, NBTau
               write(outfileindex, '(a, 10i8)')''
            enddo
            close(outfileindex)


            !> write out the total resistivity with the same relaxation time for all bands
            outfileindex= outfileindex+ 1
            write(sigmafilename, '(3a)')'rho_total_mu_',trim(adjustl(muname)),'meV.dat'
            open(unit=outfileindex, file=sigmafilename)
            write(outfileindex, '(4a)')'# \tau*\rho with unit (Ohm*m*s) is the inverse of Conductivity tensor \sum_n\sigma_n/\tau ', &
               ', at mu= ', trim(adjustl(muname)), ' meV'
            write(outfileindex, '(a)')'# relaxation time \tau_n=\tau is the same for all bands. '
            write(outfileindex, '(a,2I6)')'# NumT  NumBtau  =  ', NumT, NBTau 
            write(outfileindex, '(a,1000f8.3)')'# Tlist  =  ', KBT_array(:)/8.6173324E-5/eV2Hartree
            do ikt =1, NumT
               KBT= KBT_array(ikt)/8.6173324E-5/eV2Hartree
               write(outfileindex, '(2a, f16.4, a)')'# ', ' T = ', KBT, ' K '
               write(outfileindex, '("# Column", i5, 100i16)')(i, i=1, 10)
               write(outfileindex, '("#",19a16)')'BTau (T.ps)', 'xx', 'xy', 'xz', 'yx', 'yy', 'yz', 'zx', 'zy','zz'
               !> write out the inverse of conductivity/tau into file
               do i=1, NBTau
                  rho(1, 1:3)=sum(sigma_ohe_tensor(1:3, i, ie, ikt, :), dim=2)
                  rho(2, 1:3)=sum(sigma_ohe_tensor(4:6, i, ie, ikt, :), dim=2)
                  rho(3, 1:3)=sum(sigma_ohe_tensor(7:9, i, ie, ikt, :), dim=2)
                  if (abs(det3(rho))>eps6) then
                     call inv_r(3, rho)
                     write(outfileindex, '(19E16.6)')&
                        BTau_array(i)*Magneticfluxdensity_atomic/Relaxation_Time_Tau, &
                        rho(1, 1:3), rho(2, 1:3), rho(3, 1:3)
                  else
                     write(outfileindex, '(a)')'# error: sigma is zero since no k points contribute to the calculations of MR'
                  endif
               enddo ! i, NBTau
               write(outfileindex, '(a, 10i8)')''
            enddo
            close(outfileindex)
            
         enddo  ! ie=1, OmegaNum
      endif ! cpuid=0
       

#if defined (MPI)
      call mpi_barrier(mpi_cmw, ierr)
#endif

      contains
   
      subroutine cal_sigma_iband_k
  
         do ik= KCube3D_left(iband)%Nk_start, KCube3D_left(iband)%Nk_end
            if (cpuid.eq.0) &
               write(stdout, '(a, i8, a, i18, "   /", i18, a, f10.3, "s", a, f10.3, "s")') &
               'In sigma_OHE iband', iband, ' ik/NK', &
               ik-KCube3D_left(iband)%Nk_start+1,KCube3D_left(iband)%Nk_current, &
               ' time cost', time_end-time_start, &
               ' time left', &
               (KCube3D_left(iband)%Nk_current+KCube3D_left(iband)%Nk_start-ik)*(time_end-time_start)
   
            call now(time_start)
            EE= KCube3D_left(iband)%Ek_total(ik)
            v_k(1)= KCube3D_left(iband)%vx_total(ik)
            v_k(2)= KCube3D_left(iband)%vy_total(ik)
            v_k(3)= KCube3D_left(iband)%vz_total(ik)
            
            !> calculate df/de for each k point and each band
            do ikt=1, NumT
               KBT= KBT_array(ikt)
               do ie=1, OmegaNum
                  mu= mu_array(ie)
                  call minusdfde_calc_single(EE, KBT, mu,  minusdfde)
                  klist_iband(iband)%minusdfde(ie, ikt)= minusdfde
               enddo ! ie
            enddo ! ikt
   
   
            !> start to get the evolution of k points under magnetic field using Runge-Kutta method
            k_start= KCube3D_left(iband)%k_direct(:, ik)
            kout= 0d0
            Btau_start= 0d0
   
            !> we get the kpath by Btau_final=-30*BTauMax, but we only use half of them
            Btau_final= -exponent_max*BTauMax   !< -15 means that we can reach the accuracy as to exp(-15d0)
            !Btau_final= -10d0*BTauMax   !< test for arXiv:2201.03292 (2022)
   
            !> Runge-Kutta only applied with BTauMax>0
            !> if the magnetic field is zero. 
            call now(time_start1)
            if (BTauMax>eps3) then
               NSlice_Btau_inuse= NSlice_Btau
               call RKF45_pack(magnetic_field, bands_fermi_level(iband),  &
                    NSlice_Btau_inuse, k_start, Btau_start, Btau_final, kout, icycle, fail)
            else
               icycle= 1
               do ibtau=1, NSlice_Btau
                  kout(:, ibtau)= k_start(:)
               enddo
               NSlice_Btau_inuse = NSlice_Btau
            endif
            call now(time_end1)
            time_rkf45= time_end1-time_start1
   
            if (NSlice_Btau_inuse==1) then
               write(stdout, '(a, i6, a, i4, a, i6, a, 3f12.6)')&
                  '>>> NSlice_Btau_inuse=1 at cpuid=', cpuid, ' iband=', iband, ' ik', ik, ' k', k_start
            endif
   
            !> we omit the 
            if (fail) then
               write(stdout, '(a, i6, a, i4, a, i6, a, 3f12.6)')&
                  '>>> Runge-Kutta integration fails at cpuid=', cpuid, ' iband=', iband, ' ik', ik, ' k', k_start
               write(stdout, *)' '
               cycle
            endif
  
            call now(time_start1)
            if (NSlice_Btau_inuse==1)then  ! vxB=0
               call velocity_calc_iband(bands_fermi_level(iband), k_start, v_t)
               do ibtau=1, NSlice_Btau
                  kout(:, ibtau)= k_start(:)
                  kout_all(:, ibtau)= k_start(:)
                  klist_iband(iband)%velocity_k(:, ibtau)= v_t
               enddo
               NSlice_Btau_inuse = NSlice_Btau
            else
               do it= 1, icycle
                  k= kout(:, it) 
                  call velocity_calc_iband(bands_fermi_level(iband), k, v_t)
                  klist_iband(iband)%velocity_k(:, it)= v_t
                  kout_all(:, it) = k
               enddo ! integrate over time step
   
               !> periodic kpath in the BZ can be reused
               do i=2, NSlice_Btau_inuse/icycle
                  do j=1, icycle
                     klist_iband(iband)%velocity_k(:, j+(i-1)*icycle)= klist_iband(iband)%velocity_k(:, j)
                     kout_all(:, j+(i-1)*icycle) = kout(:, j)
                  enddo
               enddo
               do i=(NSlice_Btau_inuse/icycle)*icycle+1, NSlice_Btau
                  klist_iband(iband)%velocity_k(:, i)= klist_iband(iband)%velocity_k(:, i-(NSlice_Btau_inuse/icycle)*icycle)
                  kout_all(:, i) = kout(:, i-(NSlice_Btau_inuse/icycle)*icycle)
               enddo
            endif
            call now(time_end1)
            time_velocity=  time_end1- time_start1

            if (cpuid==0.and.iprint_level==3) then
               write(ibname, '(i6)')iband
               write(ikname, '(i6)')ik
               write(filename, '(5a)')'klist_ib_', trim(Adjustl(ibname)), '_ik', trim(Adjustl(ikname)), '.txt'
               open(unit=324232, file=filename)

               write(324232, *)'#icycle= ', icycle, ', NSlice_Btau= ', NSlice_Btau
               do i=1, NSlice_Btau
                  write(324232, "(i9, 6f16.6)") i, kout_all(:, i), klist_iband(iband)%velocity_k(:, i)
               enddo
               close(324232)

            endif

            call now(time_start1)
            !> calculate the conductivity/tau
            do ikt=1, NumT
               KBT= KBT_array(ikt)
               do ie=1, OmegaNum
                  mu= mu_array(ie)
                  minusdfde= klist_iband(iband)%minusdfde(ie, ikt)
   
                  do ibtau=1, NBTau
                     BTau= BTau_array(ibtau)

                     if (NBTau==1)then
                        ! NSlice_Btau_local= 2
                        NSlice_Btau_local= 1
                     else
                        ! NSlice_Btau_local= (ibtau-1)*NSlice_Btau_inuse/(NBTau-1)/2
                        NSlice_Btau_local= (ibtau-1)*NSlice_Btau_inuse/(NBTau-1)
                        if (NSlice_Btau_local==0)then
                           ! NSlice_Btau_local= 2
                           NSlice_Btau_local= 1
                        else
                           ! DeltaBtau= exponent_max/2d0/NSlice_Btau_local
                           DeltaBtau= exponent_max/NSlice_Btau_local
                        endif
                     endif
   
                    !do ishift=0, NSlice_Btau_local-1
                     do ishift=0, 0
                        !> here, velocity is in the cartesian coordinate
                        !> The core of Chamber formular is to get the average of velocity 
                        !> in the relaxation time approximation
                        v_k= klist_iband(iband)%velocity_k(:, 1+ishift)
                        ! if (BTau>eps3.and.NSlice_Btau_local>2) then
                        if (BTau>eps3.and.NSlice_Btau_local>1) then
                           velocity_bar_k= 0d0
                           !> five point integration(n=4)
                           do it=1, NSlice_Btau_local, 4 !(hj, j=it-1, j mod 4==0, Aj=28/45)
                             velocity_bar_k= velocity_bar_k+ &
                                 28.0d0/45.0d0*DeltaBtau*exp(-(it-1d0)*DeltaBtau)*klist_iband(iband)%velocity_k(:, it+ishift)
                           enddo
                           do it=2, NSlice_Btau_local, 4 !(hj, j=it-1, j mod 4==1, Aj=64/45)
                              velocity_bar_k= velocity_bar_k+ &
                                 64.0d0/45.0d0*DeltaBtau*exp(-(it-1d0)*DeltaBtau)*klist_iband(iband)%velocity_k(:, it+ishift)
                           enddo
                           do it=3, NSlice_Btau_local, 4 !(hj, j=it-1, j mod 4==2, Aj=24/45)
                              velocity_bar_k= velocity_bar_k+ &
                                 24.0d0/45.0d0*DeltaBtau*exp(-(it-1d0)*DeltaBtau)*klist_iband(iband)%velocity_k(:, it+ishift)
                           enddo
                           do it=4, NSlice_Btau_local, 4 !(hj, j=it-1, j mod 4==3, Aj=64/45)
                              velocity_bar_k= velocity_bar_k+ &
                                 64.0d0/45.0d0*DeltaBtau*exp(-(it-1d0)*DeltaBtau)*klist_iband(iband)%velocity_k(:, it+ishift)
                           enddo
                           velocity_bar_k= velocity_bar_k-14.0d0/45.0d0*DeltaBtau*klist_iband(iband)%velocity_k(:, 1+ishift)
                           ! set weight 
                           ! do it=1, NSlice_Btau_local
                           !    weights(it) = exp(-DeltaBtau*(it-1d0))
                           ! enddo
                           ! !> trapezoidal integral
                           ! velocity_bar_k= 0d0
                           ! do it=1, NSlice_Btau_local
                           !   velocity_bar_k= velocity_bar_k+ &
                           !      DeltaBtau*exp(-(it-1d0)*DeltaBtau)*klist_iband(iband)%velocity_k(:, it+ishift)
                           ! enddo
                           ! velocity_bar_k= velocity_bar_k &
                           ! - 0.5d0*DeltaBtau*(exp(-(NSlice_Btau_local-1d0)*DeltaBtau)&
                           ! * klist_iband(iband)%velocity_k(:, NSlice_Btau_local+ishift)  &
                           ! + klist_iband(iband)%velocity_k(:, 1+ishift))

                          !> Simpson's integral
                          !> add the first and the last term
                         ! velocity_bar_k=  weights(1)*klist_iband(iband)%velocity_k(:, 1)
                         ! velocity_bar_k= velocity_bar_k+ weights(NSlice_Btau_local)*klist_iband(iband)%velocity_k(:, NSlice_Btau_local)
                         !
                         ! !> add the middle part
                         ! do it = 2, NSlice_Btau_local-1, 2
                         !    velocity_bar_k = velocity_bar_k + 4d0 * weights(it) * klist_iband(iband)%velocity_k(:, it)
                         ! end do
                         !
                         ! do it = 3, NSlice_Btau_local-2, 2
                         !    velocity_bar_k = velocity_bar_k + 2d0 * weights(it) * klist_iband(iband)%velocity_k(:, it)
                         ! end do
                         !
                         ! velocity_bar_k = velocity_bar_k * (DeltaBtau / 3d0)

                        else
                           velocity_bar_k= v_k
                        endif
       
                        !> calculate the conductivity now
                        ik_temp= ik - KCube3D_left(iband)%Nk_start+ 1
                        sigma_symm_t= 0d0
       
                        !> Apply point group operations to the velocities, and average them
                        do j1=1, 3
                        do j2=1, 3
                        do j=1, number_group_operators
                        !> sigma_xx
                        sigma_symm_t(1)= sigma_symm_t(1)+ &
                           pgop_cart_inverse(1, j1, j)* pgop_cart_inverse(1, j2, j) &
                           *v_k(j1)*velocity_bar_k(j2)
                        !> sigma_xy
                        sigma_symm_t(2)= sigma_symm_t(2)+ &
                           pgop_cart_inverse(1, j1, j)* pgop_cart_inverse(2, j2, j) &
                           *v_k(j1)*velocity_bar_k(j2)
                        !> sigma_xz
                        sigma_symm_t(3)= sigma_symm_t(3)+ &
                           pgop_cart_inverse(1, j1, j)* pgop_cart_inverse(3, j2, j) &
                           *v_k(j1)*velocity_bar_k(j2)
                        !> sigma_yx
                        sigma_symm_t(4)= sigma_symm_t(4)+ &
                           pgop_cart_inverse(2, j1, j)* pgop_cart_inverse(1, j2, j) &
                           *v_k(j1)*velocity_bar_k(j2)
                        !> sigma_yy
                        sigma_symm_t(5)= sigma_symm_t(5)+ &
                           pgop_cart_inverse(2, j1, j)* pgop_cart_inverse(2, j2, j) &
                           *v_k(j1)*velocity_bar_k(j2)
                        !> sigma_yz
                        sigma_symm_t(6)= sigma_symm_t(6)+ &
                           pgop_cart_inverse(2, j1, j)* pgop_cart_inverse(3, j2, j) &
                           *v_k(j1)*velocity_bar_k(j2)
                        !> sigma_zx
                        sigma_symm_t(7)= sigma_symm_t(7)+ &
                           pgop_cart_inverse(3, j1, j)* pgop_cart_inverse(1, j2, j) &
                           *v_k(j1)*velocity_bar_k(j2)
                        !> sigma_zy
                        sigma_symm_t(8)= sigma_symm_t(8)+ &
                           pgop_cart_inverse(3, j1, j)* pgop_cart_inverse(2, j2, j) &
                           *v_k(j1)*velocity_bar_k(j2)
                        !> sigma_zz
                        sigma_symm_t(9)= sigma_symm_t(9)+ &
                           pgop_cart_inverse(3, j1, j)* pgop_cart_inverse(3, j2, j) &
                           *v_k(j1)*velocity_bar_k(j2)
                        enddo
                        enddo
                        enddo
                        sigma_iband_k(iband)%sigma_ohe_tensor_k_mpi(:, ibtau, ie, ikt) = &
                        sigma_iband_k(iband)%sigma_ohe_tensor_k_mpi(:, ibtau, ie, ikt) + &
                           sigma_symm_t/dble(number_group_operators)*minusdfde*KCube3D_left(iband)%weight_k(ik)
                          !sigma_symm_t/dble(number_group_operators)*minusdfde*KCube3D_left(iband)%weight_k(ik)/NSlice_Btau_local
   
                     enddo ! ishift
                  enddo ! ibtau  Btau
               enddo ! ie  mu
            enddo ! ikt KBT
            call now(time_end1)
            time_integral= time_end1- time_start1
            if (cpuid.eq.0) write(stdout, '(a, i9, a, i10)')'>> icycle = ', icycle, ' , NSlice_Btau', NSlice_Btau
            if (cpuid.eq.0) write(stdout, '(a, f16.2, a)')'>> time cost for RKF45_pack is    ', time_rkf45, ' s'
            if (cpuid.eq.0) write(stdout, '(a, f16.2, a)')'>> time cost for velocity_calc is ', time_velocity, ' s'
            if (cpuid.eq.0) write(stdout, '(a, f16.2, a)')'>> time cost for integral is      ', time_integral, ' s'
            call now(time_end)
            sigma_iband_k(iband)%time_cost_mpi(ik)= time_end- time_start
            if (cpuid.eq.0) write(stdout, '(a, f16.2, a)')'>> time cost for this loop is     ', time_end- time_start, ' s'
         enddo ! ik  kpoints
   
      end subroutine cal_sigma_iband_k

   end subroutine sigma_ohe_calc_symm

   !> calculate -df(e)/de, where f(e) is the Fermi-Dirac distribution
   !> KBT, mu, and e is in unit of Hartree
   subroutine minusdfde_calc_single(e, KBT, mu,  minusdfde)
      use wmpi
      use para, only : dp
      implicit none

      real(dp), intent(in) :: KBT  ! K_Boltzmann*Temperature in Hartree
      real(dp), intent(in) :: mu   ! Chemical potential related to the Fermi level in Hartree
      real(dp), intent(in) :: e
      real(dp), intent(out) :: minusdfde

      real(dp) :: factor
      real(dp), parameter :: maxexp= 37d0

      factor= (e- mu)/KBT
      if (abs(factor)> maxexp) then
         minusdfde = 0d0
      else
         minusdfde = 1d0/KBT* exp(factor)/ ((exp(factor)+1d0)**2)
      endif

      return
   end subroutine minusdfde_calc_single


   !> calculate -df(e)/de, where f(e) is the Fermi-Dirac distribution
   subroutine minusdfde_calc(Nband_Fermi_Level, Nk_start, Nk_end, Enk, minusdfde, KBT, mu)
      use wmpi
      use para, only : dp
      implicit none

      integer, intent(in) :: Nband_Fermi_Level
      integer, intent(in) :: Nk_start
      integer, intent(in) :: Nk_end
      real(dp), intent(in) :: KBT  ! K_Boltzmann*Temperature in eV
      real(dp), intent(in) :: mu   ! Chemical potential related to the Fermi level in eV
      real(dp), intent(in) :: Enk(Nk_start:Nk_end, Nband_Fermi_Level)
      real(dp), intent(out) :: minusdfde(Nk_start:Nk_end, Nband_Fermi_Level)

      integer :: iband, ik
      real(dp) :: factor
      real(dp), parameter :: maxexp= 37d0

      do iband=1, Nband_Fermi_Level
         do ik= Nk_start, Nk_end
            factor= (Enk(ik, iband)- mu)/KBT
            if (abs(factor)> maxexp) then
               minusdfde(ik, iband) = 0d0
            else
               minusdfde(ik, iband) = 1d0/KBT* exp(factor)/ ((exp(factor)+1d0)**2)
            endif
         enddo !ik
      enddo !iband

      return
   end subroutine minusdfde_calc

   !> calculate velocity for given k points for all bands in bands_fermi_level
   !> Eq.(2) in PRB 79, 245123(2009)
   subroutine velocity_calc(Nband_Fermi_Level, bands_fermi_level, k, velocity_k, Ek)

      use wmpi
      use para
      implicit none

      !> inout parameters
      integer, intent(in) :: Nband_Fermi_Level
      integer, intent(in) :: bands_fermi_level(Nband_Fermi_Level)

      !> k must be in fractional coordinates
      real(dp), intent(in) :: k(3) 
      real(dp), intent(inout) :: velocity_k(3, Nband_Fermi_Level)
      real(dp), intent(out) :: Ek(Nband_Fermi_Level) 

      integer :: iR, iband
      real(dp) :: kdotr
      complex(dp) :: factor
      !> velocity operator
      real(dp), allocatable :: W(:)
      complex(dp), allocatable :: vx(:, :)
      complex(dp), allocatable :: vy(:, :)
      complex(dp), allocatable :: vz(:, :)
      complex(dp), allocatable :: UU(:, :)
      complex(dp), allocatable :: Hamk_bulk(:, :)

      allocate( W(Num_wann))
      allocate( vx(Num_wann, Num_wann))
      allocate( vy(Num_wann, Num_wann))
      allocate( vz(Num_wann, Num_wann))
      allocate( UU(Num_wann, Num_wann))
      allocate( Hamk_bulk(Num_wann, Num_wann))
 
      ! calculation bulk hamiltonian
      call ham_bulk_latticegauge(k, Hamk_bulk)

      !> diagonalization by call zheev in lapack
      W= 0d0
      UU=Hamk_bulk
      call eigensystem_c( 'V', 'U', Num_wann, UU, W)
      do iband=1, Nband_Fermi_Level
         Ek(iband)= W(bands_fermi_level(iband))
      enddo

      vx= 0d0
      vy= 0d0
      vz= 0d0
      do iR= 1, Nrpts
         kdotr= k(1)*irvec(1,iR) + k(2)*irvec(2,iR) + k(3)*irvec(3,iR)
         factor= (cos(twopi*kdotr)+zi*sin(twopi*kdotr))/ndegen(iR)
         vx= vx+ zi*crvec(1, iR)*HmnR(:,:,iR)*factor
         vy= vy+ zi*crvec(2, iR)*HmnR(:,:,iR)*factor
         vz= vz+ zi*crvec(3, iR)*HmnR(:,:,iR)*factor
      enddo ! iR

      do iband= 1, Nband_Fermi_Level
         velocity_k(1, iband )= dot_product(UU(:, bands_fermi_level(iband)), matmul(vx, UU(:, bands_fermi_level(iband)))) 
         velocity_k(2, iband )= dot_product(UU(:, bands_fermi_level(iband)), matmul(vy, UU(:, bands_fermi_level(iband)))) 
         velocity_k(3, iband )= dot_product(UU(:, bands_fermi_level(iband)), matmul(vz, UU(:, bands_fermi_level(iband)))) 
      enddo

      if (allocated(W))deallocate(W)
      if (allocated(vx))deallocate(vx)
      if (allocated(vy))deallocate(vy)
      if (allocated(vz))deallocate(vz)
      if (allocated(UU))deallocate(UU)
      if (allocated(Hamk_bulk))deallocate(Hamk_bulk)

      return
   end subroutine velocity_calc


   !> calculate averaged velocity for given k points and a given band index
   !> Eq.(3) in PRB 79, 245123(2009)
   subroutine velocity_bar_calc(magnetic_field, Nband_Fermi_Level, bands_fermi_level, Btau, k0, velocity_bar_k)

      use wmpi
      use para
      implicit none

      !> inout parameters
      integer, intent(in) :: Nband_Fermi_Level
      integer, intent(in) :: bands_fermi_level(Nband_Fermi_Level)

      !> B*tau in Tesla*ps
      real(dp), intent(in) :: BTau
      real(dp), intent(in) :: magnetic_field(3)

      !> k must be in fractional coordinates
      real(dp), intent(in) :: k0(3) 
      real(dp), intent(inout) :: velocity_bar_k(3, Nband_Fermi_Level)

      integer :: iband, it
      logical :: fail
      integer :: icycle


      !> velocity operator
	   real(dp), allocatable :: W(:)
      complex(dp), allocatable :: Hamk_bulk(:, :)

      !> number of steps used in the Runge-Kutta integration
      integer :: NSlice_Btau= 2000
      integer :: NSlice_Btau_inuse

      !> Btau slices for Runge-Kutta integration
      real(dp) :: Btau_start
      real(dp) :: Btau_final
      real(dp) :: DeltaBtau

      real(dp) :: k_start(3), k(3)
      real(dp) :: velocity_k(3)
      real(dp) :: velocity_k0(3)
      real(dp), allocatable :: kout(:, :)
      real(dp), allocatable :: velocity_init(:, :, :)

      !> the transform from \omega*\tau to B*\tau is
      !> \omega*\tau=e/m_e*B*\tau
      !> so \omega*\tau = 1/0.178*B*\tau with out any units

      Btau_start= 0d0
      Btau_final= -5.6d0*Btau*2d0
      DeltaBtau= 5.6d0/NSlice_Btau*2d0

      allocate(W(Num_wann))
      allocate(Hamk_bulk(Num_wann, Num_wann))
      allocate(kout(3, NSlice_Btau))
      allocate(velocity_init(Nband_Fermi_Level, 3, NSlice_Btau))

      call ham_bulk_latticegauge(k0, Hamk_bulk)
      call eigensystem_c( 'N', 'U', Num_wann, Hamk_bulk, W)

      velocity_bar_k= 0d0
      !> calculate velocity for a given k point
      do iband=1, Nband_Fermi_Level
         !> skip the values far away from the fermi surface
         !if (abs(W(bands_fermi_level(iband)))/eV2Hartree>0.05d0)then
         if (abs(W(bands_fermi_level(iband)))/eV2Hartree>EF_integral_range)then
            cycle
         endif

         k_start= k0
         NSlice_Btau_inuse= NSlice_Btau
         call RKF45_pack(magnetic_field, bands_fermi_level(iband),  &
              NSlice_Btau, k_start, Btau_start, Btau_final, kout, icycle, fail)

         k= kout(:, 1) 
         if (NSlice_Btau_inuse==1) cycle
         !call velocity_calc_iband(iband+bands_fermi_level(1)-1, k, velocity_k0)
         call velocity_calc_iband(bands_fermi_level(iband), k, velocity_k0)
         do it= 1, NSlice_Btau_inuse
            k= kout(:, it) 
            !call velocity_calc_iband(iband+bands_fermi_level(1)-1, k, velocity_k)
            call velocity_calc_iband(bands_fermi_level(iband), k, velocity_k)
            velocity_bar_k(:, iband)= velocity_bar_k(:, iband)+ &
               DeltaBtau*exp(-(it-1d0)*DeltaBtau)*velocity_k
         enddo ! integrate over time step
         velocity_bar_k(:, iband)= velocity_bar_k(:, iband) &
            - 0.5d0*DeltaBtau*(exp(-(NSlice_Btau_inuse-1d0)*DeltaBtau)*velocity_k &
            + velocity_k0)
      enddo ! iband

      return
   end subroutine velocity_bar_calc

   !> calculate the derivation of kn(t) for each band crossing Fermi level
   !> Eq.(4) in PRB 79, 245123(2009)
   !> Btau is the time, basically, Btau is not used here, 
   !> since dk/dBt is not explicitly dependent on Btau
   !> kt is the fractional coordinates of a k point
   !> kdot is dk/d(Bt)
   !> iband is the band index
   !> \hbar*dkn(t)/d(Bt)= − e \vec{vn(kn(t))} \times \vec{Bdirection}  
   !> this subroutine would be called by RKFS as an external subroutine
   subroutine dkdt(magnetic_field, Btau, kt, kdot, iband)

      use wmpi
      use para, only : dp, Bdirection, Echarge, hbar, eps3
      implicit none

      !> inout variables
      integer, intent(in) :: iband
      real(dp), intent(in) :: magnetic_field(3)
      real(dp), intent(in) :: Btau !> not used in this subroutine
      real(dp), intent(in) :: kt(3)  !> 3 means the three axias x, y and z
      real(dp), intent(inout) :: kdot(3)

      !> v_x, v_y, v_z
      real(dp) :: velocity_k(3)
      real(dp) :: kdot_cart(3)

      kdot = 0d0
      !> calculate velocity for a given k point
      call velocity_calc_iband(iband, kt, velocity_k)
      if (sum(abs(velocity_k))<eps3) return

      !> \hbar*dkn(t)/d(Bt)= − e \vec{vn(kn(t))} \times \vec{Bdirection}  
      kdot_cart(1)= -velocity_k(2)*Bdirection(3)+ velocity_k(3)*Bdirection(2)
      kdot_cart(2)= -velocity_k(3)*Bdirection(1)+ velocity_k(1)*Bdirection(3)
      kdot_cart(3)= -velocity_k(1)*Bdirection(2)+ velocity_k(2)*Bdirection(1)

      !> coordinates transformation
      !> transform from cartesian to fractional coordinates
      call cart_direct_rec(kdot_cart, kdot)

      !> 1E-12 comes from the relaxation time units ps
      !> 1E-10 comes from Angstrom
      !> \hbar dk/d(Bt)=-e\vec{v(k)}\times \vec{e_B}
      !> \vec{v(k)}= 1/\hbar*\partitial Hmn(k)/\partitial k
      !> In the velocity subroutine, we don't consider the 1/\hbar
     !kdot= kdot*(Echarge/hbar)**2*1E-10*1E-12   ! 1/m
      !> 1E-10 comes from Angstrom
     !kdot= kdot*1E-10   ! 1/Angstrom

      return
   end subroutine dkdt

   !> calculate velocity for given k points for all bands in bands_fermi_level
   !> Eq.(2) in PRB 79, 245123(2009)
   subroutine velocity_calc_iband2(iband, k, velocity_k, E_iband)

      use wmpi
      use para
      implicit none

      !> inout parameters
      integer, intent(in) :: iband

      !> k must be in fractional coordinates
      real(dp), intent(in) :: k(3) 
      real(dp), intent(out) :: E_iband
      real(dp), intent(inout) :: velocity_k(3)

      integer :: iR 
      real(dp) :: kdotr
      complex(dp) :: factor
      !> velocity operator
	   real(dp), allocatable :: W(:)
      complex(dp), allocatable :: vx(:, :)
      complex(dp), allocatable :: vy(:, :)
      complex(dp), allocatable :: vz(:, :)
      complex(dp), allocatable :: UU(:)
      complex(dp), allocatable :: Hamk_bulk(:, :)

      allocate( W (Num_wann))
      allocate( vx(Num_wann, Num_wann))
      allocate( vy(Num_wann, Num_wann))
      allocate( vz(Num_wann, Num_wann))
      allocate( UU(Num_wann))
      allocate( Hamk_bulk(Num_wann, Num_wann))
      Hamk_bulk= 0d0
      W= 0d0
      UU= 0d0


      ! calculation bulk hamiltonian
      call ham_bulk_latticegauge(k, Hamk_bulk)

      call zheevx_pack('V', 'U', Num_wann, iband, iband, Hamk_bulk, W, UU)

      E_iband= W(1)
      !> Only the energy levels close to the Fermi level contribute to the conductivity
      if (abs(W(1))/eV2Hartree>EF_integral_range) then
         velocity_k= 0d0
         return
      endif
 
      vx= 0d0
      vy= 0d0
      vz= 0d0
      do iR= 1, Nrpts
         kdotr= k(1)*irvec(1,iR) + k(2)*irvec(2,iR) + k(3)*irvec(3,iR)
         factor= (cos(twopi*kdotr)+zi*sin(twopi*kdotr))/ndegen(iR)
         vx= vx+ zi*crvec(1, iR)*HmnR(:,:,iR)*factor
         vy= vy+ zi*crvec(2, iR)*HmnR(:,:,iR)*factor
         vz= vz+ zi*crvec(3, iR)*HmnR(:,:,iR)*factor
      enddo ! iR

      !> velocity is in cartesian coordinate
      velocity_k(1)= dot_product(UU(:), matmul(vx, UU(:))) 
      velocity_k(2)= dot_product(UU(:), matmul(vy, UU(:))) 
      velocity_k(3)= dot_product(UU(:), matmul(vz, UU(:))) 

      !write(stdout, '(a, i6, a, 3f12.8, a, f10.5, a, 3f14.8)')"iband", iband, ' k', k, ' eig', W(1),  &
      !   ' velocity_k', velocity_k

      deallocate(W, vx, vy, vz, Hamk_bulk, UU)

      return
   end subroutine velocity_calc_iband2



   !> calculate velocity for given k points for all bands in bands_fermi_level
   !> Eq.(2) in PRB 79, 245123(2009)
   subroutine velocity_calc_iband(iband, k, velocity_k)

      use wmpi
      use para
      implicit none

      !> inout parameters
      integer, intent(in) :: iband

      !> k must be in fractional coordinates
      real(dp), intent(in) :: k(3) 
      real(dp), intent(inout) :: velocity_k(3)

      integer :: iR
      real(dp) :: kdotr
      complex(dp) :: factor
      !> velocity operator
      real(dp), allocatable :: W(:)
      complex(dp), allocatable :: vx(:, :)
      complex(dp), allocatable :: vy(:, :)
      complex(dp), allocatable :: vz(:, :)
      complex(dp), allocatable :: UU(:)
      complex(dp), allocatable :: Hamk_bulk(:, :)

      allocate( W (Num_wann))
      allocate( vx(Num_wann, Num_wann))
      allocate( vy(Num_wann, Num_wann))
      allocate( vz(Num_wann, Num_wann))
      allocate( UU(Num_wann))
      allocate( Hamk_bulk(Num_wann, Num_wann))
      Hamk_bulk= 0d0
      W= 0d0
      UU= 0d0

      ! calculation bulk hamiltonian
      call ham_bulk_latticegauge(k, Hamk_bulk)

!     call eigensystem_c( 'V', 'U', Num_wann, Hamk_bulk, W)

!     !> Only the energy levels close to the Fermi level contribute to the conductivity
!     if (abs(W(iband))/eV2Hartree>EF_integral_range) then
!        velocity_k= 0d0
!        return
!     endif

!     UU= Hamk_bulk(:, iband)

      call zheevx_pack('V', 'U', Num_wann, iband, iband, Hamk_bulk, W, UU)

      !> Only the energy levels close to the Fermi level contribute to the conductivity
      if (abs(W(1))/eV2Hartree>EF_integral_range) then
         velocity_k= 0d0
         return
      endif


      vx= 0d0
      vy= 0d0
      vz= 0d0
      do iR= 1, Nrpts
         kdotr= k(1)*irvec(1,iR) + k(2)*irvec(2,iR) + k(3)*irvec(3,iR)
         factor= (cos(twopi*kdotr)+zi*sin(twopi*kdotr))/ndegen(iR)
         vx= vx+ zi*crvec(1, iR)*HmnR(:,:,iR)*factor
         vy= vy+ zi*crvec(2, iR)*HmnR(:,:,iR)*factor
         vz= vz+ zi*crvec(3, iR)*HmnR(:,:,iR)*factor
      enddo ! iR

      !> velocity is in cartesian coordinate
      velocity_k(1)= dot_product(UU(:), matmul(vx, UU(:))) 
      velocity_k(2)= dot_product(UU(:), matmul(vy, UU(:))) 
      velocity_k(3)= dot_product(UU(:), matmul(vz, UU(:))) 

      !write(stdout, '(a, i6, a, 3f12.8, a, f10.5, a, 3f14.8)')"iband", iband, ' k', k, ' eig', W(1),  &
      !   ' velocity_k', velocity_k

      deallocate(W, vx, vy, vz, Hamk_bulk, UU)

      return
   end subroutine velocity_calc_iband


   !> a pack for rkf45
   subroutine RKF45_pack(magnetic_field, iband, NSlice_Btau, k_start, Btau_start, Btau_final, kout, icycle, fail)
      use wmpi
      use para, only : dp, eps9, eps8, cpuid, stdout, eps6, Bdirection, eps3, RKF45_PERIODIC_LEVEL
      implicit none

      !> inout variable 
      integer , intent(in) :: iband
      real(dp), intent(in) :: magnetic_field(3)
      integer , intent(inout) :: NSlice_Btau
      real(dp), intent(in) :: k_start(3)
      real(dp), intent(in) :: Btau_start
      real(dp), intent(in) :: Btau_final
      real(dp), intent(out) :: kout(3, NSlice_Btau)
      integer, intent(out) :: icycle
      logical, intent(out) :: fail

      !> three axes: x, y and z
      integer, parameter :: neqn= 3
      external :: dkdt
      real(dp), external :: norm


      !> parameters for RKFS
      integer  :: it, i, j
      real(dp) :: t, tout
      real(dp) :: kt(neqn), kdot(neqn), kdiff(3)
      real(dp) :: DeltaBtau, relerr, abserr 
      integer  :: iflag, iter

      real(dp) :: time_start, dis, dis_smallest

      real(dp) :: velocity_k0(3), vcrossB(3)


      icycle=NSlice_Btau
      fail= .false.

      DeltaBtau = Btau_final/dble(NSlice_Btau)
      iflag = 1
      relerr= 1d-10
      abserr= 0d0

      !> starting point
      kt= k_start
      t = Btau_start
      tout= Btau_start

      !> check whether the velocity is zero, if so, we just abandom those points
      call velocity_calc_iband(iband, k_start, velocity_k0)

      !> Vk x B
      call cross_product(-velocity_k0, Bdirection, vcrossB)

      if (sum(abs(vcrossB))<eps6) then
         NSlice_Btau= 1
         do it= 1, NSlice_Btau
            kout(:, it)= k_start
         enddo
         return
      endif

      it = 0
      iter= 0
      10 continue
         !> integration from t to tout starts from kt, for the first call of rkfs, t=tout
         !> after rkfs, t=tout
         !> 
         call now(time_start)

        !> if Runge-Kutta integration is not converged, then this calculation is failed.
         if (iter>3) then 
            fail=.true.
            goto 80
         endif

         call r8_rkf45(dkdt, neqn, kt, kdot, t, tout, relerr, abserr, iflag, iband, magnetic_field)
        !if (iter>0) print*, 'iter', iter, 'cpuid', cpuid, kt, iflag
        !call now(time_end)

        !if (cpuid==0)write(stdout, '(a, i5, a,2i6,a,f6.2,2a, 2f12.6, a, 3f8.4, a, 3f8.4)')&
        !            'iband', iband, &
        !            'In RKFS_PACK', it, iflag, ' time for rkfs ', time_end- time_start, ' s', &
        !            ' t, tout', t, tout, ' kt', kt, ' kdot', kdot

         !> kdot = 0 means the magnetic field don't change the k anymore.
         if (sum(abs(kdot))<eps6)then
            it= it+ 1
            kout(:, it)= kt
            NSlice_Btau= it
            goto 80
         endif
         goto (80, 20, 30, 40, 50, 60, 70, 80), iflag

      !iflag = 2 -- integration reached tout. indicates successful retur
      !             and is the normal mode for continuing integration.
      20 tout= t+ DeltaBtau  !> goto the next step
         it= it+ 1
         kout(:, it)= kt
         iter=0
         if (it+1> NSlice_Btau) goto 80
         

         call periodic_diff(kout(:,2), kout(:,1), kdiff)
         if (it>2)dis_smallest= norm(kdiff)/RKF45_PERIODIC_LEVEL

         !> check if kout(:, it)==kout(:, 1)
         !> if it's a close orbit, we don't have to calculate all of them.
         if (it>10) then
            call periodic_diff(kout(:,it), kout(:,1), kdiff)
            dis= norm(kdiff)
            if (dis<dis_smallest) then
               icycle= it- 1
               do i=2, NSlice_Btau/icycle
                  do j=1, icycle
                     kout(:, j+(i-1)*icycle)=kout(:, j)
                  enddo
               enddo
               do i=(NSlice_Btau/icycle)*icycle+1, NSlice_Btau
                  kout(:, i)=kout(:, i-(NSlice_Btau/icycle)*icycle)
               enddo
               goto 80
            endif
         endif


         if (dabs(t-Btau_final) >eps8 .and. t>Btau_final) then
            goto 10
         else
            goto 80  ! exit
         endif

      !iflag= 3 -- integration was not completed because relative error
      !       tolerance was too small. relerr has been increased
      !       appropriately for continuing.
      30 continue
         iter= iter+ 1
         goto 10

      !iflag= 4 -- integration was not completed because more than
      !       3000 derivative evaluations were needed. this
      !       is approximately 500 steps.
      40 continue 
         iter= iter+ 1
         goto 10

      !iflag= 5 -- integration was not completed because solution
      !        vanished making a pure relative error test
      !        impossible. must use non-zero abserr to continue.
      !        using the one-step integration mode for one step
      !        is a good way to proceed.
      50 abserr=eps9
         iter= iter+ 1
         continue
         goto 10

      !iflag= 6 -- integration was not completed because requested
      !       accuracy could not be achieved using smallest
      !       allowable stepsize. user must increase the error
      !       tolerance before continued integration can be
      !       attempted.
      60 relerr=10d0*relerr
         iter= iter+ 1
         continue
         iflag= 2
         goto 10

      !iflag = 7 -- it is likely that rkf45 is inefficient for solving
      !        this problem. too much output is restricting the
      !        natural stepsize choice. use the one-step integrator
      !        mode.
      70 continue
         iter= iter+ 1
         iflag= 2
         goto 10

      !iflag= 8 -- invalid input parameters
      !       this indicator occurs if any of the following is
      !       satisfied -   neqn <= 0
      !                     t=tout  and  iflag /= +1 or -1
      !                     relerr or abserr < 0.
      !                     iflag == 0  or  < -2  or  > 8
      80 continue

      return
   end subroutine RKF45_pack

!> Get the k points evolution under the magnetic field. Using Runge-Kutta method
!> Rewrote By QuanSheng Wu (wuquansheng@gmail.com)
!> References : 
!> Ref1: Electrons in metals and semiconductors, R.G. Chambers,
!> Ref2: Ab initio investigation of magnetic transport properties by Wannier interpolation, 
!> PHYSICAL REVIEW B 79, 245123  2009 , Yi Liu, Hai-Jun Zhang, and Yugui Yao
!> see Eq.(4) in Ref2
   subroutine evolve_k_ohe
      use wmpi
      use para
      implicit none
     
      !real(dp) :: OmegaTau  !> e*B/m*Tau

      integer :: ik, ib
      integer :: ierr, it, i, ibtau

      real(dp) :: v_t(3), v_t2(3)

      real(dp) :: k(3), k_start(3)
      real(dp) :: magnetic_field(3)
      
      !> Btau slices for Runge-Kutta integration
      real(dp) :: Btau, Btau_start, Btau_final
      logical :: fail
      integer :: icycle

      !> energy bands 
      real(dp) :: E_iband
      real(dp), allocatable :: Ek(:)  ! in eV

      !> 3-component velocity for each band and each k point 
      real(dp), allocatable :: velocity_k(:, :)
      
      integer :: NSlice_Btau 
      integer, allocatable :: NSlice_Btau_all(:), NSlice_Btau_all_mpi(:)
      real(dp), allocatable :: kout(:, :), kout_all(:, :, :), kout_all_mpi(:, :, :)

      !> Bands crossing Fermi level
      integer :: Nband_Fermi_Level
      integer, allocatable :: bands_fermi_level_temp(:)
      integer, allocatable :: bands_fermi_level(:)

      !> file index
      integer, allocatable  :: myfileindex(:)
      character(40) :: evolvefilename
      character(40) :: bandname


      !> First we calculate howmany and which bands cross the Fermi level
      !> Nband_Fermi_Level and bands_fermi_level will be updated after
      !> get_bands_cross_fermilevel, we hope 1000 is quite enough
      Nband_Fermi_Level= 1000
      allocate(bands_fermi_level_temp(Nband_Fermi_Level))
      allocate(myfileindex(Nband_Fermi_Level))

      if (NumberofSelectedBands/=0) then
         Nband_Fermi_Level= NumberofSelectedBands
         allocate(bands_fermi_level(Nband_Fermi_Level))
         do i=1, Nband_Fermi_Level
            bands_fermi_level(i)= Selected_band_index(i) 
         enddo
      else
         !> First we calculate howmany and which bands cross the Fermi level
         call get_bands_cross_fermilevel(Nband_Fermi_Level, bands_fermi_level_temp)
         allocate(bands_fermi_level(Nband_Fermi_Level))
         bands_fermi_level= bands_fermi_level_temp(1:Nband_Fermi_Level)
      endif

      allocate(Ek(Nband_Fermi_Level))
      Ek= 0d0

      !> setup NSlice_Btau
      !> NSlice_Btau should be the integer times of NBTau
      NSlice_Btau= Nslice_BTau_Max
      if (cpuid.eq.0) write(stdout, *) ' NSlice_Btau :', NSlice_Btau
      allocate(NSlice_Btau_all(nk3_band))
      allocate(NSlice_Btau_all_mpi(nk3_band))
      allocate(kout(3, NSlice_Btau))
      allocate(kout_all(3, NSlice_Btau, nk3_band))
      allocate(kout_all_mpi(3, NSlice_Btau, nk3_band))
      kout_all_mpi= 0d0; kout_all=0d0; kout=0d0
      NSlice_Btau_all= 0; NSlice_Btau_all_mpi= 0

      allocate( velocity_k(3, Nband_Fermi_Level))
      velocity_k= 0d0
   
   
      k= Single_KPOINT_3D_DIRECT
      call velocity_calc(Nband_Fermi_Level, bands_fermi_level, k, velocity_k, Ek)
 
      !> file index for different bands
      do ib=1, Nband_Fermi_Level
         outfileindex= outfileindex+ 1
         myfileindex(ib)= outfileindex
      enddo

      do ib=1, Nband_Fermi_Level
         if (cpuid.eq.0) then
            write(bandname, '(i10)')bands_fermi_level(ib)
            write(evolvefilename, '(3a)')'evolve_band_', trim(adjustl(bandname)), '.txt'
            open(unit=myfileindex(ib), file=evolvefilename)
            write(myfileindex(ib), '(a10, i5, a16, f16.8)')'# evolve k ', bands_fermi_level(ib), 'energy level', Ek(ib)
            write(myfileindex(ib), '("#", a13, 24a16)')'BTau (T.ps)', &
                'kx', 'ky', 'kz', 'vx', 'vy', 'vz', 'k1', 'k2','k3', 'Energy(ev)', "vx'", "vy'", "vz'"
         endif
      enddo
#if defined (MPI)
      call mpi_barrier(mpi_cmw, ierr)
#endif

     magnetic_field(1)=Bx
     magnetic_field(2)=By
     magnetic_field(3)=Bz
  
      !> exclude all kpoints with zero velocity x B and large energy away from Fermi level
     do ib= 1, Nband_Fermi_Level 

        NSlice_Btau_all_mpi = 0; NSlice_Btau_all=0
        kout_all= 0d0; kout_all_mpi= 0d0
        do ik= 1+cpuid, nk3_band, num_cpu
           k_start = kpath_3d(:, ik)
           if (cpuid.eq.0) then
              write( stdout, *)'ib, ik', ib, ik
           endif
 
           !> start to get the evolution of k points under magnetic field using Runge-Kutta method
           !k_start= Single_KPOINT_3D_DIRECT
           kout= 0d0
           Btau_start= 0d0
           Btau_final= -15d0*BTauMax
           NSlice_Btau = Nslice_BTau_Max
       
           !> Runge-Kutta only applied with BTauMax>0
           !> if the magnetic field is zero. 
           icycle = 0
           if (BTauMax>eps3) then
              call RKF45_pack(magnetic_field, bands_fermi_level(ib),  &
                 NSlice_Btau, k_start, Btau_start, Btau_final, kout, icycle, fail)
           else
              stop "ERROR: please set finite BtauMax in the file wt.in"
              do ibtau=1, NSlice_Btau
                 kout(:, ibtau)= k_start(:)
              enddo
           endif
           NSlice_Btau_all_mpi(ik)= NSlice_Btau
           kout_all_mpi(:, :, ik)= kout
        enddo

#if defined (MPI)
        call mpi_allreduce(NSlice_Btau_all_mpi,NSlice_Btau_all,size(NSlice_Btau_all),&
                         mpi_in,mpi_sum,mpi_cmw,ierr)
        call mpi_allreduce(kout_all_mpi,kout_all,size(kout_all),&
                         mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
        NSlice_Btau_all= NSlice_Btau_all_mpi
        kout_all= kout_all_mpi
#endif
 

        do ik= 1, nk3_band
           k_start = kpath_3d(:, ik)
           call velocity_calc(Nband_Fermi_Level, bands_fermi_level, k_start, velocity_k, Ek)
           if (cpuid.eq.0) then
           endif
       
           if (cpuid.eq.0) then
              if (NSlice_Btau_all(ik)==1) then
                 write(myfileindex(ib), '(1X, a, 3f12.5, a, i8)')"# Starting k point (fractional coordinates) :", k_start, ' ith-band ', bands_fermi_level(ib)
                 write(myfileindex(ib), '(a, f16.8, a)')'# At energy level', Ek(ib)/eV2Hartree, ' eV'
                 write(myfileindex(ib), '(a)')">> No data for output since VxB=0 or far away from Fermi level |E-E_F|>0.05eV" 
                 write(myfileindex(ib), '(a)')" "
              else
                 write(myfileindex(ib), '(a)')" "
                 write(myfileindex(ib), '(a)')"# Quasi-particle's trajectory under magnetic field"
                 write(myfileindex(ib), '(a, 3f12.5)')'# Magnetic field is along (cartesian coordinates)', Bdirection
                 write(myfileindex(ib), '(1X, a, 3f12.5, a, i8)')"# Starting k point (fractional coordinates) :", k_start, ' ith-band ', bands_fermi_level(ib)
                 write(myfileindex(ib), '(a, f16.8, a)')'# At energy level', Ek(ib)/eV2Hartree, ' eV'
                 write(myfileindex(ib), '("#", a13, 20a16)')'BTau (T.ps)', &
                     'kx', 'ky', 'kz', 'vx', 'vy', 'vz', 'k1', 'k2','k3', "vx'", "vy'", "vz'"
                 do it=1, NSlice_Btau_all(ik)
                    k= kout_all(:, it, ik) 
                    !call velocity_calc_iband(ib+bands_fermi_level(1)-1, k, v_t)
                    call velocity_calc_iband2(bands_fermi_level(ib), k, v_t, E_iband)
                    BTau = -(it-1.0)/NSlice_Btau_all(ik)*15d0*BTauMax
                    call direct_cart_rec(kout_all(:, it, ik), k)
                    call project_k3_to_kplane_defined_by_direction(v_t, Bdirection, v_t2)
                    write(myfileindex(ib), '(100f16.8)')Btau*Magneticfluxdensity_atomic/Relaxation_Time_Tau, &
                        k, v_t, kout_all(:, it, ik), v_t2
                 enddo
              endif
           endif

        enddo ! ik
        if (cpuid.eq.0) then
           close(myfileindex(ib))
        endif
#if defined (MPI)
        call mpi_barrier(mpi_cmw, ierr)
#endif
     enddo ! ib
    
     return
  end subroutine evolve_k_ohe


!> Calculate band resolved magnetoconductance with R.G.Chambers's formula
!> Rewrote By QuanSheng Wu (wuquansheng@gmail.com)
!> References : Electrons in metals and semiconductors, R.G. Chambers,
!> Ab initio investigation of magnetic transport properties by Wannier interpolation, 
!> PHYSICAL REVIEW B 79, 245123  2009 , Yi Liu, Hai-Jun Zhang, and Yugui Yao
!> This subroutine will only give the conductivity/tau(kz) instead of conductivity.
   subroutine sigma_k_ohe
      use wmpi
      use para
      implicit none
     
      real(dp), allocatable :: sigma_k_ohe_tensor(:, :) 
      real(dp), allocatable :: sigma_k_ohe_tensor_mpi(:, :)

      !real(dp) :: OmegaTau  !> e*B/m*Tau
      real(dp) :: coeff

      real(dp) :: mu, BTau, KBT
      integer :: ibtau

      integer :: Nk_total, Nk_current, Nk_start, Nk_end, knv3
      integer :: knv3_left, knv3_left_mod, knv3_mod
      integer :: ik, ib, ik1, ik2, ik3, ierr, it, i
      logical :: fail
      integer :: icycle

      real(dp) :: v_t(3), v_k(3), k(3), k_start(3), magnetic_field(3)
      
      real(dp) :: time_start, time_end, time_start0, time_end0
      integer :: NSlice_Btau_inuse

      !> Btau slices for Runge-Kutta integration
      real(dp) :: Btau_start, Btau_final, DeltaBtau


      !> energy bands 
      real(dp) :: EE
      real(dp), allocatable :: Ek(:)  ! in eV
      real(dp), allocatable :: Enk(:, :)

      !> minus fermi derivation
      real(dp) :: minusdfde

      !> 3-component velocity for each band and each k point 
      real(dp), allocatable :: velocity_k(:, :)
      real(dp), allocatable :: velocity_bar_k(:)
      real(dp), allocatable :: velocity(:, :, :)
      
      !> 3-component velocity for each band and each k point 
      !> Eq.(3) in PRB 79, 245123(2009)
      real(dp), allocatable :: velocity_bar(:, :, :)

      type(kcube_type), allocatable :: KCube2D_left(:)

      !> some arrays for mpi
      integer, allocatable :: info(:, :) !> howmany integers to be sent on each cpu 
      integer, allocatable :: Displs(:, :)

      !> number of steps used in the Runge-Kutta integration
      integer :: NSlice_Btau
      integer :: NSlice_Btau_local
      real(dp), allocatable :: kout(:, :)

      !> define some arrays for different bands. Since there are different number
      !> of k points left for different bands.
      type klist_iband_type
         !> dim=3*NSlice_Btau
         real(dp), allocatable :: klist_rkfs(:, :)
         !> dim=3*NSlice_Btau
         real(dp), allocatable :: velocity_k(:, :)

         !> calculate -df(e)/de, where f(e) is the Fermi-Dirac distribution
         !> dim=(NumT, OmegaNum))
         real(dp), allocatable :: minusdfde(:, :)
      end type klist_iband_type
      type(klist_iband_type), allocatable :: klist_iband(:)

      !> Bands crossing Fermi level
      integer :: Nband_Fermi_Level
      integer, allocatable :: bands_fermi_level_temp(:)
      integer, allocatable :: bands_fermi_level(:)

      !> file index
      integer, allocatable  :: myfileindex(:)
      character(40) :: sigmafilename, bandname


      type(kcube_type) :: KCube2D


      !> Nband_Fermi_Level and bands_fermi_level will be updated after
      !> get_bands_cross_fermilevel, we hope 1000 is quite enough
      Nband_Fermi_Level= 1000
      allocate(bands_fermi_level_temp(Nband_Fermi_Level))
      allocate(myfileindex(Nband_Fermi_Level))
      allocate(KCube2D_left(Nband_Fermi_Level))
      allocate(klist_iband(Nband_Fermi_Level))

      !> First we calculate howmany and which bands cross the Fermi level
      call get_bands_cross_fermilevel(Nband_Fermi_Level, bands_fermi_level_temp)
      allocate(bands_fermi_level(Nband_Fermi_Level))
      bands_fermi_level= bands_fermi_level_temp(1:Nband_Fermi_Level)

      !> setup NSlice_Btau
      !> NSlice_Btau should be the integer times of NBTau
      NSlice_Btau= Nslice_BTau_Max
      if (cpuid.eq.0) write(stdout, *) ' NSlice_Btau :', NSlice_Btau

      allocate( sigma_k_ohe_tensor    (9, Nband_Fermi_Level))
      allocate( sigma_k_ohe_tensor_mpi(9, Nband_Fermi_Level))
      sigma_k_ohe_tensor    = 0d0
      sigma_k_ohe_tensor_mpi= 0d0

      !> file index for different bands
      do ib=1, Nband_Fermi_Level
         outfileindex= outfileindex+ 1
         myfileindex(ib)= outfileindex
      enddo

      do ib=1, Nband_Fermi_Level
         if (cpuid.eq.0) then
            write(bandname, '(i10)')bands_fermi_level(ib)
            write(sigmafilename, '(3a)')'sigma_k_band_', trim(adjustl(bandname)), '.dat'
            open(unit=myfileindex(ib), file=sigmafilename)
            write(myfileindex(ib), '(a20, i5, a16, f16.8)')'# sigma_k at band ', bands_fermi_level(ib)
            write(myfileindex(ib), '("#",a12,20a16)')'ik3', 'xx', 'xy', 'xz', 'yx', 'yy', 'yz', 'zx', 'zy','zz'
         endif
      enddo


      ! for each k3 point the kmesh is different
      do ik3= 1, Nk3
         sigma_k_ohe_tensor    = 0d0
         sigma_k_ohe_tensor_mpi= 0d0

         !> Distribute k kpoints in the 2DCube into different MPI threads
         knv3= Nk1*Nk2
         KCube2D%Nk_total= knv3
         knv3_mod= mod(knv3, num_cpu)
         if (knv3_mod==0) then  !> perfect divided
            KCube2D%Nk_current= knv3/num_cpu
            KCube2D%Nk_start=1+ knv3*cpuid/num_cpu
            KCube2D%Nk_end  =(1+cpuid)*knv3/num_cpu
         else if (knv3/num_cpu==0) then    !> Number of MPI threads is large than knv3
            KCube2D%Nk_current= 1 !> one k piont per MPI thread
            KCube2D%Nk_start= cpuid+ 1 !> one k piont per MPI thread
            KCube2D%Nk_end  = cpuid+ 1
            if (cpuid+1 > knv3) then
               KCube2D%Nk_start= 1
               KCube2D%Nk_end  = 0
            endif
         else
            KCube2D%Nk_current= knv3/num_cpu+ 1
            if (cpuid< knv3_mod) then
               KCube2D%Nk_start= 1+ cpuid*KCube2D%Nk_current
               KCube2D%Nk_end  = (1+cpuid)*KCube2D%Nk_current
            else
               KCube2D%Nk_start= knv3_mod*KCube2D%Nk_current+ &
                  (cpuid-knv3_mod)*(KCube2D%Nk_current-1)+1
               KCube2D%Nk_end  = knv3_mod*KCube2D%Nk_current+ &
                  (cpuid-knv3_mod+1)*(KCube2D%Nk_current-1)
            endif
         endif
    
         !> calculate the volume of the k cube
         allocate(KCube2D%k_direct(3, KCube2D%Nk_start:KCube2D%Nk_end))
    
         do ik= KCube2D%Nk_start, KCube2D%Nk_end
            ik1= (ik-1)/(Nk2)+1
            ik2= (ik- (ik1-1)*Nk2)
            KCube2D%k_direct(:, ik)= K3D_start_cube+ K3D_vec1_cube*(ik1-1)/dble(Nk1)  &
             + K3D_vec2_cube*(ik2-1)/dble(Nk2)  &
             + K3D_vec3_cube*(ik3-1)/dble(Nk3)
         enddo

         Nk_end= KCube2D%Nk_end
         Nk_total= KCube2D%Nk_total
         Nk_start= KCube2D%Nk_start
         Nk_current= KCube2D%Nk_current
   
         allocate( Ek(Nband_Fermi_Level))
         allocate( Enk(Nk_start:Nk_end, Nband_Fermi_Level))
         allocate( velocity(3, Nk_start:Nk_end, Nband_Fermi_Level))
         allocate( velocity_k(3, Nband_Fermi_Level))
         allocate( velocity_bar(3, Nk_start:Nk_end, Nband_Fermi_Level))
         allocate( velocity_bar_k(3))
         Ek= 0d0
         Enk= 0d0
         velocity= 0d0
         velocity_k= 0d0
         velocity_bar= 0d0
         velocity_bar_k= 0d0
   
         time_start= 0d0
         time_end= 0d0
         do ik= Nk_start, Nk_end
            if (cpuid.eq.0.and. mod(ik, 100).eq.0) &
               write(stdout, '(a, i18, "   /", i18, a, f10.3, "s")') 'ik/NK', &
               ik-Nk_start,Nk_current, 'time left', &
               (Nk_current+Nk_start-ik)*(time_end-time_start)/num_cpu
   
            call now(time_start)
            k= KCube2D%k_direct(:, ik)
   
            call velocity_calc(Nband_Fermi_Level, bands_fermi_level, k, velocity_k, Ek)
            velocity(:, ik, :)= velocity_k
            Enk(ik, :)= Ek
            call now(time_end)
         enddo
   
         !> exclude all kpoints with zero velocity x B and large energy away from Fermi level
         do ib= 1, Nband_Fermi_Level 
   
            !> first check howmany k points left
            it= 0
            do ik= Nk_start, Nk_end
               !> check whether v x B=0, which means the magnetic field is parallel with velocity
               v_t= velocity(:, ik, ib)
              !vcrossB(1)= -v_t(2)*Bdirection(3)+ v_t(3)*Bdirection(2)
              !vcrossB(2)= -v_t(3)*Bdirection(1)+ v_t(1)*Bdirection(3)
              !vcrossB(3)= -v_t(1)*Bdirection(2)+ v_t(2)*Bdirection(1)
              !if (abs(Enk(ik, ib))/eV2Hartree<0.05d0.and.dsqrt(sum((abs(vcrossB)**2)))>eps3) then
               if (abs(Enk(ik, ib))/eV2Hartree<EF_integral_range) then
                  it = it+ 1
               endif
            enddo ! ik
   
            KCube2D_left(ib)%Nk_current= it
            if (it>0) then
               allocate(KCube2D_left(ib)%ik_array(it))
               allocate(KCube2D_left(ib)%Ek_local(it))
               allocate(KCube2D_left(ib)%vx_local(it))
               allocate(KCube2D_left(ib)%vy_local(it))
               allocate(KCube2D_left(ib)%vz_local(it))
            else 
               allocate(KCube2D_left(ib)%ik_array(1))  !> only useful for mpi_allgatherv
               allocate(KCube2D_left(ib)%Ek_local(1))
               allocate(KCube2D_left(ib)%vx_local(1))
               allocate(KCube2D_left(ib)%vy_local(1))
               allocate(KCube2D_left(ib)%vz_local(1))
            endif
            KCube2D_left(ib)%ik_array= 0
            KCube2D_left(ib)%Ek_local= 0d0
            KCube2D_left(ib)%vx_local= 0d0
            KCube2D_left(ib)%vy_local= 0d0
            KCube2D_left(ib)%vz_local= 0d0
   
   
            it= 0
            do ik= Nk_start, Nk_end
               !> check whether v x B=0, which means the magnetic field is parallel with velocity
               v_t= velocity(:, ik, ib)
              !vcrossB(1)= -v_t(2)*Bdirection(3)+ v_t(3)*Bdirection(2)
              !vcrossB(2)= -v_t(3)*Bdirection(1)+ v_t(1)*Bdirection(3)
              !vcrossB(3)= -v_t(1)*Bdirection(2)+ v_t(2)*Bdirection(1)
              !if (abs(Enk(ik, ib))/eV2Hartree<0.05d0.and.dsqrt(sum((abs(vcrossB)**2)))>eps3) then
               if (abs(Enk(ik, ib))/eV2Hartree<EF_integral_range) then
                  it = it+ 1
                  KCube2D_left(ib)%ik_array(it) = ik
                  KCube2D_left(ib)%Ek_local(it) = Enk(ik, ib)
                  KCube2D_left(ib)%vx_local(it) = v_t(1)
                  KCube2D_left(ib)%vy_local(it) = v_t(2)
                  KCube2D_left(ib)%vz_local(it) = v_t(3)
               endif
            enddo ! ik
         enddo ! ib
   
         !> try to get the total number of k points left for each band 
         do ib=1, Nband_Fermi_Level
#if defined (MPI)
            call mpi_allreduce(KCube2D_left(ib)%Nk_current,KCube2D_left(ib)%Nk_total,1,&
                         mpi_in,mpi_sum,mpi_cmw,ierr)
#else
            KCube2D_left(ib)%Nk_total= KCube2D_left(ib)%Nk_current
#endif
         enddo
   
         !> gather the number of k points left into a array info
         allocate(info(num_cpu, Nband_Fermi_Level))
         info= -1
         do ib= 1, Nband_Fermi_Level
#if defined (MPI)
            call mpi_allgather(KCube2D_left(ib)%Nk_current, 1, mpi_in, info(:, ib), 1, mpi_in, mpi_cmw, ierr)
#else     
            info(1, ib)= KCube2D_left(ib)%Nk_current
#endif
         enddo
     
         !> An array for mpi_allgatherv
         allocate(Displs(num_cpu+1, Nband_Fermi_Level))
         Displs= 0
         do ib=1, Nband_Fermi_Level
            Displs(1, ib)=0
            do i=2, num_cpu+1
               Displs(i, ib)=Displs(i-1, ib)+ info(i-1, ib)
            enddo
         enddo ! ib
     
     
         !> put all the kpoints left together
         do ib=1, Nband_Fermi_Level
            allocate(KCube2D_left(ib)%IKleft_array(KCube2D_left(ib)%Nk_total))
            allocate(KCube2D_left(ib)%Ek_total(KCube2D_left(ib)%Nk_total))
            allocate(KCube2D_left(ib)%vx_total(KCube2D_left(ib)%Nk_total))
            allocate(KCube2D_left(ib)%vy_total(KCube2D_left(ib)%Nk_total))
            allocate(KCube2D_left(ib)%vz_total(KCube2D_left(ib)%Nk_total))
            KCube2D_left(ib)%IKleft_array = 0
            KCube2D_left(ib)%Ek_total= 0d0
            KCube2D_left(ib)%vx_total= 0d0
            KCube2D_left(ib)%vy_total= 0d0
            KCube2D_left(ib)%vz_total= 0d0
         enddo  ! ib
        !> gather Enk and velocity 
#if defined (MPI)
         do ib=1, Nband_Fermi_Level
            call mpi_allgatherv(KCube2D_left(ib)%ik_array, KCube2D_left(ib)%Nk_current, &
                                mpi_in, KCube2D_left(ib)%IKleft_array, &
                                info(:, ib), Displs(:, ib), mpi_in, mpi_cmw, ierr)
     
            call mpi_allgatherv(KCube2D_left(ib)%Ek_local, KCube2D_left(ib)%Nk_current, &
                                mpi_dp, KCube2D_left(ib)%Ek_total, &
                                info(:, ib), Displs(:, ib), mpi_dp, mpi_cmw, ierr)
            call mpi_allgatherv(KCube2D_left(ib)%vx_local, KCube2D_left(ib)%Nk_current, &
                                mpi_dp, KCube2D_left(ib)%vx_total, &
                                info(:, ib), Displs(:, ib), mpi_dp, mpi_cmw, ierr)
            call mpi_allgatherv(KCube2D_left(ib)%vy_local, KCube2D_left(ib)%Nk_current, &
                                mpi_dp, KCube2D_left(ib)%vy_total, &
                                info(:, ib), Displs(:, ib), mpi_dp, mpi_cmw, ierr)
            call mpi_allgatherv(KCube2D_left(ib)%vz_local, KCube2D_left(ib)%Nk_current, &
                                mpi_dp, KCube2D_left(ib)%vz_total, &
                                info(:, ib), Displs(:, ib), mpi_dp, mpi_cmw, ierr)
         enddo ! ib
#else
#endif
     
         !> redistribute all those k points into different cpus
         if (cpuid.eq.0) then
            write(stdout, '(a)')' '
            write(stdout, '(a10, i10, a20, i10, a5, i10)')' There are ', KCube2D%Nk_total, ' k points at k3= ' , ik3, ' in' , Nk3
            do ib= 1, Nband_Fermi_Level
               write(stdout, '(a, i10, a, i10)')' However there are only ', KCube2D_left(ib)%Nk_total, &
                  ' k points contribute to the conductance calculations for band ', bands_fermi_level(ib)
            enddo
         endif
     
         !> distribute the left kpoints into different processors
         !> for different bands, the number of left kpoints is different.
         do ib= 1, Nband_Fermi_Level
            knv3_left= KCube2D_left(ib)%Nk_total
            knv3_left_mod= mod(knv3_left, num_cpu)
            if (knv3_left_mod==0) then  !> perfect divided
               KCube2D_left(ib)%Nk_current= knv3_left/num_cpu
               KCube2D_left(ib)%Nk_start=1+ knv3_left*cpuid/num_cpu
               KCube2D_left(ib)%Nk_end  =(1+cpuid)*knv3_left/num_cpu
            else if (knv3_left/num_cpu==0) then    !> Number of MPI threads is large than knv3_left
               KCube2D_left(ib)%Nk_current= 1 !> one k piont per MPI thread
               KCube2D_left(ib)%Nk_start= cpuid+ 1 !> one k piont per MPI thread
               KCube2D_left(ib)%Nk_end  = cpuid+ 1
               if (cpuid+1 > knv3_left) then
                  KCube2D_left(ib)%Nk_start= 1
                  KCube2D_left(ib)%Nk_end  = 0
               endif
            else
               KCube2D_left(ib)%Nk_current= knv3_left/num_cpu+ 1
               if (cpuid< knv3_left_mod) then
                  KCube2D_left(ib)%Nk_start= 1+ cpuid*KCube2D_left(ib)%Nk_current
                  KCube2D_left(ib)%Nk_end  = (1+cpuid)*KCube2D_left(ib)%Nk_current
               else
                  KCube2D_left(ib)%Nk_start= knv3_left_mod*KCube2D_left(ib)%Nk_current+ &
                     (cpuid-knv3_left_mod)*(KCube2D_left(ib)%Nk_current-1)+1
                  KCube2D_left(ib)%Nk_end  = knv3_left_mod*KCube2D_left(ib)%Nk_current+ &
                     (cpuid-knv3_left_mod+1)*(KCube2D_left(ib)%Nk_current-1)
               endif
            endif
     
     
            allocate(KCube2D_left(ib)%k_direct(3, KCube2D_left(ib)%Nk_start:KCube2D_left(ib)%Nk_end))
       
            do ik= KCube2D_left(ib)%Nk_start, KCube2D_left(ib)%Nk_end
               i= KCube2D_left(ib)%IKleft_array(ik)
               ik1= (i-1)/(Nk2)+1
               ik2= (i- (ik1-1)*Nk2)
               KCube2D_left(ib)%k_direct(:, ik)= K3D_start_cube+ K3D_vec1_cube*(ik1-1)/dble(Nk1)  &
                + K3D_vec2_cube*(ik2-1)/dble(Nk2)  &
                + K3D_vec3_cube*(ik3-1)/dble(Nk3)
            enddo
         enddo
     
         do ib=1, Nband_Fermi_Level
            allocate(klist_iband(ib)%klist_rkfs(3, NSlice_Btau))
            allocate(klist_iband(ib)%velocity_k(3, NSlice_Btau))
            klist_iband(ib)%klist_rkfs= 0d0
            klist_iband(ib)%velocity_k= 0d0
         enddo
     
         !> a temp array used in RKFS
         allocate(kout(3, NSlice_Btau))
         kout= 0d0
     
         !> now we turn to use Runge-Kutta method to get all the kpoints from (0, BTauMax)
         !> and we calculate the conductivity/Tau over different bands and different k points
         time_start= 0d0
         time_end  = 0d0
    
         magnetic_field(1)=Bx
         magnetic_field(2)=By
         magnetic_field(3)=Bz
         call now(time_start0)
         do ib= 1, Nband_Fermi_Level
!           print *, 'cpuid, ib, nk', cpuid, ib, KCube2D_left(ib)%Nk_current
            !> dim=(Nk_start: Nk_end, NumT, OmegaNum))
            do ik= KCube2D_left(ib)%Nk_start, KCube2D_left(ib)%Nk_end
               if (cpuid.eq.0) &
                  write(stdout, '(a, i8, a, i18, "   /", i18, a, f10.3, "s")') 'In sigma_OHE iband', ib, ' ik/NK', &
                  ik-KCube2D_left(ib)%Nk_start,KCube2D_left(ib)%Nk_current, ' time left', &
                  (KCube2D_left(ib)%Nk_current+KCube2D_left(ib)%Nk_start-ik)*(time_end-time_start)
     
               call now(time_start)
               EE= KCube2D_left(ib)%Ek_total(ik)
               v_k(1)= KCube2D_left(ib)%vx_total(ik)
               v_k(2)= KCube2D_left(ib)%vy_total(ik)
               v_k(3)= KCube2D_left(ib)%vz_total(ik)
              
               !> Kelvin to Hartree
               KBT= Tmin*8.6173324E-5*eV2Hartree
               mu= OmegaMin
               call minusdfde_calc_single(EE, KBT, mu,  minusdfde)
     
               !> start to get the evolution of k points under magnetic field using Runge-Kutta method
               k_start= KCube2D_left(ib)%k_direct(:, ik)
               kout= 0d0
               Btau_start= 0d0
               Btau_final= -15d0*BTauMax
     
               !> Runge-Kutta only applied with BTauMax>0
               !> if the magnetic field is zero. 
               NSlice_Btau_inuse= NSlice_Btau
               if (BTauMax>eps3) then
                  call RKF45_pack(magnetic_field, bands_fermi_level(ib),  &
                       NSlice_Btau, k_start, Btau_start, Btau_final, kout, icycle, fail)
               else
                  NSlice_Btau_inuse = 1
                  do ibtau=1, NSlice_Btau
                     kout(:, ibtau)= k_start(:)
                  enddo
               endif
     
      
               klist_iband(ib)%klist_rkfs= kout
     
              !if (NSlice_Btau_inuse==1) cycle
               do it= 1, NSlice_Btau_inuse
                  k= kout(:, it) 
                 !call velocity_calc_iband(ib+bands_fermi_level(1)-1, k, v_t)
                  call velocity_calc_iband(bands_fermi_level(ib), k, v_t)
                  klist_iband(ib)%velocity_k(:, it)= v_t
               enddo ! integrate over time step
     
     
               !> calculate the conductivity/tau
               DeltaBtau= 15d0/NSlice_Btau_local
     
               if (BTau>eps3) then
                  velocity_bar_k= 0d0
                  do it=1, NSlice_Btau_local
                     velocity_bar_k= velocity_bar_k+ &
                        DeltaBtau*exp(-(it-1d0)*DeltaBtau)*klist_iband(ib)%velocity_k(:, it)
                  enddo
                  velocity_bar_k= velocity_bar_k &
                  - 0.5d0*DeltaBtau*(exp(-(NSlice_Btau_local-1d0)*DeltaBtau)&
                  * klist_iband(ib)%velocity_k(:, NSlice_Btau_local)  &
                  + klist_iband(ib)%velocity_k(:, 1))
               else
                  velocity_bar_k= v_k
               endif
     
               !> calculate the conductivity now
               !> sigma_xx
               sigma_k_ohe_tensor_mpi(1, ib) = &
                  sigma_k_ohe_tensor_mpi(1, ib) + &
                  v_k(1)*velocity_bar_k(1)*minusdfde
               !> sigma_xy
               sigma_k_ohe_tensor_mpi(2, ib) = &
                  sigma_k_ohe_tensor_mpi(2, ib) + &
                  v_k(1)*velocity_bar_k(2)*minusdfde
               !> sigma_xz
               sigma_k_ohe_tensor_mpi(3, ib) = &
                  sigma_k_ohe_tensor_mpi(3, ib) + &
                  v_k(1)*velocity_bar_k(3)*minusdfde
               !> sigma_yx
               sigma_k_ohe_tensor_mpi(4, ib) = &
                  sigma_k_ohe_tensor_mpi(4, ib) + &
                  v_k(2)*velocity_bar_k(1)*minusdfde
               !> sigma_yy
               sigma_k_ohe_tensor_mpi(5, ib) = &
                  sigma_k_ohe_tensor_mpi(5, ib) + &
                  v_k(2)*velocity_bar_k(2)*minusdfde
               !> sigma_yz
               sigma_k_ohe_tensor_mpi(6, ib) = &
                  sigma_k_ohe_tensor_mpi(6, ib) + &
                  v_k(2)*velocity_bar_k(3)*minusdfde
               !> sigma_zx
               sigma_k_ohe_tensor_mpi(7, ib) = &
                  sigma_k_ohe_tensor_mpi(7, ib) + &
                  v_k(3)*velocity_bar_k(1)*minusdfde
               !> sigma_zy
               sigma_k_ohe_tensor_mpi(8, ib) = &
                  sigma_k_ohe_tensor_mpi(8, ib) + &
                  v_k(3)*velocity_bar_k(2)*minusdfde
               !> sigma_zz
               sigma_k_ohe_tensor_mpi(9, ib) = &
                  sigma_k_ohe_tensor_mpi(9, ib) + &
                  v_k(3)*velocity_bar_k(3)*minusdfde
     
               call now(time_end)
            enddo ! ik  kpoints
         enddo ! ib bands
     
         call now(time_end0)

#if defined (MPI)
         call mpi_allreduce(sigma_k_ohe_tensor_mpi,sigma_k_ohe_tensor,size(sigma_k_ohe_tensor),&
                      mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
         sigma_k_ohe_tensor= sigma_k_ohe_tensor_mpi
#endif
         sigma_k_ohe_tensor= sigma_k_ohe_tensor/Nk_total

         deallocate( Ek)
         deallocate( Enk)
         deallocate( kout)
         deallocate( info)
         deallocate( Displs)
         deallocate( velocity)
         deallocate( velocity_k)
         deallocate( velocity_bar)
         deallocate( velocity_bar_k)
         deallocate(KCube2D%k_direct)

         do ib= 1, Nband_Fermi_Level
            deallocate(KCube2D_left(ib)%ik_array)
            deallocate(KCube2D_left(ib)%Ek_local)
            deallocate(KCube2D_left(ib)%vx_local)
            deallocate(KCube2D_left(ib)%vy_local)
            deallocate(KCube2D_left(ib)%vz_local)
            deallocate(KCube2D_left(ib)%k_direct)
            deallocate(klist_iband(ib)%klist_rkfs)
            deallocate(klist_iband(ib)%velocity_k)
            deallocate(KCube2D_left(ib)%IKleft_array)
            deallocate(KCube2D_left(ib)%Ek_total)
            deallocate(KCube2D_left(ib)%vx_total)
            deallocate(KCube2D_left(ib)%vy_total)
            deallocate(KCube2D_left(ib)%vz_total)
         enddo

         if (cpuid.eq.0) then
            !> write out the conductivity/tau into file
            do ib=1, Nband_Fermi_Level
               coeff= Echarge**2/hbar/Bohr_radius/Origin_cell%CellVolume/kCubeVolume*Origin_cell%ReciprocalCellVolume*1E-12
               coeff= coeff/Time_atomic 
               write(myfileindex(ib), '(i6,20E16.6)')ik3, sigma_k_ohe_tensor(:, ib)*coeff
            enddo ! ib, band
         endif ! cpuid=0
 
      enddo ! ik3

      if (cpuid.eq.0) then
         do ib=1, Nband_Fermi_Level
            close(myfileindex(ib))
         enddo ! ib, band
      endif ! cpuid=0

      !> In the end, we start to care about the units of the conductivity/tau
      !> the conductivity/tau is in units of Ohm^-1*m^-1*s^-1

      return
   end subroutine sigma_k_ohe

subroutine r8_rkf45 (f, neqn, y, yp, t, tout, relerr, abserr, flag, iband, magnetic_field)

!*****************************************************************************80
!
!! R8_RKF45 carries out the Runge-Kutta-Fehlberg method (double precision).
!
!  Discussion:
!
!    This routine is primarily designed to solve non-stiff and mildly stiff
!    differential equations when derivative evaluations are inexpensive.
!    It should generally not be used when the user is demanding
!    high accuracy.
!
!    This routine integrates a system of NEQN first-order ordinary differential
!    equations of the form:
!
!      dY(i)/dT = F(T,Y(1),Y(2),...,Y(NEQN))
!
!    where the Y(1:NEQN) are given at T.
!
!    Typically the subroutine is used to integrate from T to TOUT but it
!    can be used as a one-step integrator to advance the solution a
!    single step in the direction of TOUT.  On return, the parameters in
!    the call list are set for continuing the integration.  The user has
!    only to call again (and perhaps define a new value for TOUT).
!
!    Before the first call, the user must
!
!    * supply the subroutine F(T,Y,YP) to evaluate the right hand side;
!      and declare F in an EXTERNAL statement;
!
!    * initialize the parameters:
!      NEQN, Y(1:NEQN), T, TOUT, RELERR, ABSERR, FLAG.
!      In particular, T should initially be the starting point for integration,
!      Y should be the value of the initial conditions, and FLAG should
!      normally be +1.
!
!    Normally, the user only sets the value of FLAG before the first call, and
!    thereafter, the program manages the value.  On the first call, FLAG should
!    normally be +1 (or -1 for single step mode.)  On normal return, FLAG will
!    have been reset by the program to the value of 2 (or -2 in single
!    step mode), and the user can continue to call the routine with that
!    value of FLAG.
!
!    (When the input magnitude of FLAG is 1, this indicates to the program
!    that it is necessary to do some initialization work.  An input magnitude
!    of 2 lets the program know that that initialization can be skipped,
!    and that useful information was computed earlier.)
!
!    The routine returns with all the information needed to continue
!    the integration.  If the integration reached TOUT, the user need only
!    define a new TOUT and call again.  In the one-step integrator
!    mode, returning with FLAG = -2, the user must keep in mind that
!    each step taken is in the direction of the current TOUT.  Upon
!    reaching TOUT, indicated by the output value of FLAG switching to 2,
!    the user must define a new TOUT and reset FLAG to -2 to continue
!    in the one-step integrator mode.
!
!    In some cases, an error or difficulty occurs during a call.  In that case,
!    the output value of FLAG is used to indicate that there is a problem
!    that the user must address.  These values include:
!
!    * 3, integration was not completed because the input value of RELERR, the
!      relative error tolerance, was too small.  RELERR has been increased
!      appropriately for continuing.  If the user accepts the output value of
!      RELERR, then simply reset FLAG to 2 and continue.
!
!    * 4, integration was not completed because more than MAXNFE derivative
!      evaluations were needed.  This is approximately (MAXNFE/6) steps.
!      The user may continue by simply calling again.  The function counter
!      will be reset to 0, and another MAXNFE function evaluations are allowed.
!
!    * 5, integration was not completed because the solution vanished,
!      making a pure relative error test impossible.  The user must use
!      a non-zero ABSERR to continue.  Using the one-step integration mode
!      for one step is a good way to proceed.
!
!    * 6, integration was not completed because the requested accuracy
!      could not be achieved, even using the smallest allowable stepsize.
!      The user must increase the error tolerances ABSERR or RELERR before
!      continuing.  It is also necessary to reset FLAG to 2 (or -2 when
!      the one-step integration mode is being used).  The occurrence of
!      FLAG = 6 indicates a trouble spot.  The solution is changing
!      rapidly, or a singularity may be present.  It often is inadvisable
!      to continue.
!
!    * 7, it is likely that this routine is inefficient for solving
!      this problem.  Too much output is restricting the natural stepsize
!      choice.  The user should use the one-step integration mode with
!      the stepsize determined by the code.  If the user insists upon
!      continuing the integration, reset FLAG to 2 before calling
!      again.  Otherwise, execution will be terminated.
!
!    * 8, invalid input parameters, indicates one of the following:
!      NEQN <= 0;
!      T = TOUT and |FLAG| /= 1;
!      RELERR < 0 or ABSERR < 0;
!      FLAG == 0  or FLAG < -2 or 8 < FLAG.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2011
!
!  Author:
!
!    Original FORTRAN77 version by Herman Watts, Lawrence Shampine.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Erwin Fehlberg,
!    Low-order Classical Runge-Kutta Formulas with Stepsize Control,
!    NASA Technical Report R-315, 1969.
!
!    Lawrence Shampine, Herman Watts, S Davenport,
!    Solving Non-stiff Ordinary Differential Equations - The State of the Art,
!    SIAM Review,
!    Volume 18, pages 376-411, 1976.
!
!  Parameters:
!
!    Input, external F, a user-supplied subroutine to evaluate the
!    derivatives Y'(T), of the form:
!
!      subroutine f ( t, y, yp )
!      real ( kind = 8 ) t
!      real ( kind = 8 ) y(*)
!      real ( kind = 8 ) yp(*)
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations to be integrated.
!
!    Input/output, real ( kind = 8 ) Y(NEQN), the current solution vector at T.
!
!    Input/output, real ( kind = 8 ) YP(NEQN), the current value of the
!    derivative of the dependent variable.  The user should not set or alter
!    this information!
!
!    Input/output, real ( kind = 8 ) T, the current value of the independent
!    variable.
!
!    Input, real ( kind = 8 ) TOUT, the output point at which solution is
!    desired.  TOUT = T is allowed on the first call only, in which case
!    the routine returns with FLAG = 2 if continuation is possible.
!
!    Input, real ( kind = 8 ) RELERR, ABSERR, the relative and absolute
!    error tolerances for the local error test.  At each step the code
!    requires:
!      abs ( local error ) <= RELERR * abs ( Y ) + ABSERR
!    for each component of the local error and the solution vector Y.
!    RELERR cannot be "too small".  If the routine believes RELERR has been
!    set too small, it will reset RELERR to an acceptable value and return
!    immediately for user action.
!
!    Input/output, integer ( kind = 4 ) FLAG, indicator for status of 
!    integration.  On the first call, set FLAG to +1 for normal use, or to -1 
!    for single step mode.  On return, a value of 2 or -2 indicates normal 
!    progress, while any other value indicates a problem that should 
!    be addressed.
!
  implicit none

  integer ( kind = 4 ) neqn
  integer ( kind = 4 ) iband
  real ( kind=8) ::  magnetic_field(3)

  real ( kind = 8 ) abserr
  real ( kind = 8 ), save :: abserr_save = -1.0D+00
  real ( kind = 8 ) ae
  real ( kind = 8 ) dt
  real ( kind = 8 ) ee
  real ( kind = 8 ) eeoet
  real ( kind = 8 ) eps
  real ( kind = 8 ) esttol
  real ( kind = 8 ) et
  external f
  real ( kind = 8 ) f1(neqn)
  real ( kind = 8 ) f2(neqn)
  real ( kind = 8 ) f3(neqn)
  real ( kind = 8 ) f4(neqn)
  real ( kind = 8 ) f5(neqn)
  real ( kind = 8 ), save :: h = -1.0D+00
  logical hfaild
  real ( kind = 8 ) hmin
  integer ( kind = 4 ) flag
  integer ( kind = 4 ), save :: flag_save = -1000
  integer ( kind = 4 ), save :: init = -1000
  integer ( kind = 4 ) k
  integer ( kind = 4 ), save :: kflag = -1000
  integer ( kind = 4 ), save :: kop = -1
  integer ( kind = 4 ), parameter :: maxnfe = 3000
  integer ( kind = 4 ) mflag
  integer ( kind = 4 ), save :: nfe = -1
  logical output
  real ( kind = 8 ) relerr
  real ( kind = 8 ) relerr_min
  real ( kind = 8 ), save :: relerr_save = -1.0D+00
  real ( kind = 8 ), parameter :: remin = 1.0E-12
  real ( kind = 8 ) s
  real ( kind = 8 ) scale
  real ( kind = 8 ) t
  real ( kind = 8 ) tol
  real ( kind = 8 ) toln
  real ( kind = 8 ) tout
  real ( kind = 8 ) y(neqn)
  real ( kind = 8 ) yp(neqn)
  real ( kind = 8 ) ypk
!
!  Check the input parameters.
!
  eps = epsilon ( eps )

  if ( neqn < 1 ) then
    flag = 8
    return
  end if

  if ( relerr < 0.0D+00 ) then
    flag = 8
    return
  end if

  if ( abserr < 0.0D+00 ) then
    flag = 8
    return
  end if

  if ( flag == 0 .or. 8 < flag .or. flag < -2 ) then
    flag = 8
    return
  end if

  mflag = abs ( flag )
!
!  Is this a continuation call?
!
  if ( mflag /= 1 ) then

    if ( abs(t- tout)<eps .and. kflag /= 3 ) then
      flag = 8
      return
    end if

    if ( mflag == 2 ) then

      if ( kflag == 3 ) then

        flag = flag_save
        mflag = abs ( flag )

      else if ( init == 0 ) then

        flag = flag_save

      else if ( kflag == 4 ) then

        nfe = 0

      else if ( kflag == 5 .and. abserr <= eps ) then

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_RKF45 - Fatal error!'
        write ( *, '(a)' ) '  KFLAG = 5 and ABSERR = 0.0'
        stop

      else if ( &
        kflag == 6 .and. &
        relerr <= relerr_save .and. &
        abserr <= abserr_save ) then

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_RKF45 - Fatal error!'
        write ( *, '(a)' ) '  KFLAG = 6 and'
        write ( *, '(a)' ) '  RELERR <= RELERR_SAVE and'
        write ( *, '(a)' ) '  ABSERR <= ABSERR_SAVE'

        stop

      end if
!
!  FLAG = 3, 4, 5, 6, 7 or 8.
!
    else

      if ( flag == 3 ) then

        flag = flag_save
        if ( kflag == 3 ) then
          mflag = abs ( flag )
        end if

      else if ( flag == 4 ) then

        nfe = 0
        flag = flag_save
        if ( kflag == 3 ) then
          mflag = abs ( flag )
        end if

      else if ( flag == 5 .and. 0.0D+00 < abserr ) then

        flag = flag_save
        if ( kflag == 3 ) then
          mflag = abs ( flag )
        end if
!
!  Integration cannot be continued because the user did not respond to
!  the instructions pertaining to FLAG = 5, 6, 7 or 8.
!
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_RKF45 - Fatal error!'
        write ( *, '(a)' ) '  Integration cannot be continued.'
        write ( *, '(a)' ) '  The user did not respond to the output'
        write ( *, '(a)' ) '  value FLAG = 5, 6, 7, or 8.'
        stop
      end if

    end if

  end if
!
!  Save the input value of FLAG.
!  Set the continuation flag KFLAG for subsequent input checking.
!
  flag_save = flag
  kflag = 0
!
!  Save RELERR and ABSERR for checking input on subsequent calls.
!
  relerr_save = relerr
  abserr_save = abserr
!
!  Restrict the relative error tolerance to be at least
!
!    2 * EPS + REMIN
!
!  to avoid limiting precision difficulties arising from impossible
!  accuracy requests.
!
  relerr_min = 2.0D+00 * epsilon ( relerr_min ) + remin
!
!  Is the relative error tolerance too small?
!
  if ( relerr < relerr_min ) then
    relerr = relerr_min
    flag = 3
    kflag = 3
    return
  end if

  dt = tout - t
!
!  Initialization:
!
!  Set the initialization completion indicator, INIT;
!  set the indicator for too many output points, KOP;
!  evaluate the initial derivatives;
!  set the counter for function evaluations, NFE;
!  estimate the starting stepsize.
!
  if ( mflag == 1 ) then

    init = 0
    kop = 0
    call f (magnetic_field, t, y, yp, iband )
    nfe = 1

    if ( abs(t- tout)<eps ) then
      flag = 2
      return
    end if

  end if

  if ( init == 0 ) then

    init = 1
    h = abs ( dt )
    toln = 0.0D+00

    do k = 1, neqn
      tol = relerr * abs ( y(k) ) + abserr
      if ( 0.0D+00 < tol ) then
        toln = tol
        ypk = abs ( yp(k) )
        if ( tol < ypk * h**5 ) then
          h = ( tol / ypk )**0.2D+00
        end if
      end if
    end do

    if ( toln <= 0.0D+00 ) then
      h = 0.0D+00
    end if

    h = max ( h, 26.0D+00 * eps * max ( abs ( t ), abs ( dt ) ) )
    flag_save = sign ( 2, flag )

  end if
!
!  Set the stepsize for integration in the direction from T to TOUT.
!
  h = sign ( h, dt )
!
!  Test to see if too may output points are being requested.
!
  if ( 2.0D+00 * abs ( dt ) <= abs ( h ) ) then
    kop = kop + 1
  end if
!
!  Unnecessary frequency of output.
!  Seems like an error I'm willing to tolerate!
!
! if ( kop == 100 ) then
  if ( kop == 10000 ) then
    kop = 0
    flag = 7
    return
  end if
!
!  If we are too close to the output point, then simply extrapolate and return.
!
  if ( abs ( dt ) <= 26.0D+00 * eps * abs ( t ) ) then
    t = tout
    y(1:neqn) = y(1:neqn) + dt * yp(1:neqn)
    call f (magnetic_field, t, y, yp, iband )
    nfe = nfe + 1
    flag = 2
    return
  end if
!
!  Initialize the output point indicator.
!
  output = .false.
!
!  To avoid premature underflow in the error tolerance function,
!  scale the error tolerances.
!
  scale = 2.0D+00 / relerr
  ae = scale * abserr
!
!  Step by step integration.
!
  do

    hfaild = .false.
!
!  Set the smallest allowable stepsize.
!
    hmin = 26.0D+00 * eps * abs ( t )
!
!  Adjust the stepsize if necessary to hit the output point.
!
!  Look ahead two steps to avoid drastic changes in the stepsize and
!  thus lessen the impact of output points on the code.
!
    dt = tout - t

    if ( 2.0D+00 * abs ( h ) <= abs ( dt ) ) then

    else
!
!  Will the next successful step complete the integration to the output point?
!
      if ( abs ( dt ) <= abs ( h ) ) then
        output = .true.
        h = dt
      else
        h = 0.5D+00 * dt
      end if

    end if
!
!  Here begins the core integrator for taking a single step.
!
!  The tolerances have been scaled to avoid premature underflow in
!  computing the error tolerance function ET.
!  To avoid problems with zero crossings, relative error is measured
!  using the average of the magnitudes of the solution at the
!  beginning and end of a step.
!  The error estimate formula has been grouped to control loss of
!  significance.
!
!  To distinguish the various arguments, H is not permitted
!  to become smaller than 26 units of roundoff in T.
!  Practical limits on the change in the stepsize are enforced to
!  smooth the stepsize selection process and to avoid excessive
!  chattering on problems having discontinuities.
!  To prevent unnecessary failures, the code uses 9/10 the stepsize
!  it estimates will succeed.
!
!  After a step failure, the stepsize is not allowed to increase for
!  the next attempted step.  This makes the code more efficient on
!  problems having discontinuities and more effective in general
!  since local extrapolation is being used and extra caution seems
!  warranted.
!
!  Test the number of derivative function evaluations.
!  If okay, try to advance the integration from T to T+H.
!
    do
!
!  Have we done too much work?
!
      if ( maxnfe < nfe ) then
        flag = 4
        kflag = 4
        return
      end if
!
!  Advance an approximate solution over one step of length H.
!
      call r8_fehl ( f, neqn, y, t, h, yp, f1, f2, f3, f4, f5, f1, iband, magnetic_field)
      nfe = nfe + 5
!
!  Compute and test allowable tolerances versus local error estimates
!  and remove scaling of tolerances.  The relative error is
!  measured with respect to the average of the magnitudes of the
!  solution at the beginning and end of the step.
!
      eeoet = 0.0D+00

      do k = 1, neqn

        et = abs ( y(k) ) + abs ( f1(k) ) + ae

        if ( et <= 0.0D+00 ) then
          flag = 5
          return
        end if

        ee = abs &
        ( ( -2090.0D+00 * yp(k) &
          + ( 21970.0D+00 * f3(k) - 15048.0D+00 * f4(k) ) &
          ) &
        + ( 22528.0D+00 * f2(k) - 27360.0D+00 * f5(k) ) &
        )

        eeoet = max ( eeoet, ee / et )

      end do

      esttol = abs ( h ) * eeoet * scale / 752400.0D+00

      if ( esttol <= 1.0D+00 ) then
        exit
      end if
!
!  Unsuccessful step.  Reduce the stepsize, try again.
!  The decrease is limited to a factor of 1/10.
!
      hfaild = .true.
      output = .false.

      if ( esttol < 59049.0D+00 ) then
        s = 0.9D+00 / esttol**0.2D+00
      else
        s = 0.1D+00
      end if

      h = s * h

      if ( abs ( h ) < hmin ) then
        flag = 6
        kflag = 6
        return
      end if

    end do
!
!  We exited the loop because we took a successful step.
!  Store the solution for T+H, and evaluate the derivative there.
!
    t = t + h
    y(1:neqn) = f1(1:neqn)
    call f (magnetic_field, t, y, yp, iband )
    nfe = nfe + 1
!
!  Choose the next stepsize.  The increase is limited to a factor of 5.
!  If the step failed, the next stepsize is not allowed to increase.
!
    if ( 0.0001889568D+00 < esttol ) then
      s = 0.9D+00 / esttol**0.2D+00
    else
      s = 5.0D+00
    end if

    if ( hfaild ) then
      s = min ( s, 1.0D+00 )
    end if

    h = sign ( max ( s * abs ( h ), hmin ), h )
!
!  End of core integrator
!
!  Should we take another step?
!
    if ( output ) then
      t = tout
      flag = 2
      return
    end if

    if ( flag <= 0 ) then
      exit
    end if

  end do
!
!  One step integration mode.
!
  flag = -2

  return
end subroutine r8_rkf45

subroutine r8_fehl ( f, neqn, y, t, h, yp, f1, f2, f3, f4, f5, s, iband, magnetic_field)

!*****************************************************************************80
!
!! R8_FEHL takes one Fehlberg fourth-fifth order step (double precision).
!
!  Discussion:
!
!    This routine integrates a system of NEQN first order ordinary differential
!    equations of the form
!      dY(i)/dT = F(T,Y(1:NEQN))
!    where the initial values Y and the initial derivatives
!    YP are specified at the starting point T.
!
!    The routine advances the solution over the fixed step H and returns
!    the fifth order (sixth order accurate locally) solution
!    approximation at T+H in array S.
!
!    The formulas have been grouped to control loss of significance.
!    The routine should be called with an H not smaller than 13 units of
!    roundoff in T so that the various independent arguments can be
!    distinguished.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 March 2004
!
!  Author:
!
!    Original FORTRAN77 version by Herman Watts, Lawrence Shampine.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Erwin Fehlberg,
!    Low-order Classical Runge-Kutta Formulas with Stepsize Control,
!    NASA Technical Report R-315, 1969.
!
!    Lawrence Shampine, Herman Watts, S Davenport,
!    Solving Non-stiff Ordinary Differential Equations - The State of the Art,
!    SIAM Review,
!    Volume 18, pages 376-411, 1976.
!
!  Parameters:
!
!    Input, external F, a user-supplied subroutine to evaluate the
!    derivatives Y'(T), of the form:
!
!      subroutine f ( t, y, yp )
!      real ( kind = 8 ) t
!      real ( kind = 8 ) y(*)
!      real ( kind = 8 ) yp(*)
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations to be integrated.
!
!    Input, real ( kind = 8 ) Y(NEQN), the current value of the
!    dependent variable.
!
!    Input, real ( kind = 8 ) T, the current value of the independent
!    variable.
!
!    Input, real ( kind = 8 ) H, the step size to take.
!
!    Input, real ( kind = 8 ) YP(NEQN), the current value of the
!    derivative of the dependent variable.
!
!    Output, real ( kind = 8 ) F1(NEQN), F2(NEQN), F3(NEQN), F4(NEQN),
!    F5(NEQN), derivative values needed for the computation.
!
!    Output, real ( kind = 8 ) S(NEQN), the estimate of the solution at T+H.
!
  implicit none

  integer ( kind = 4 ) neqn
  integer ( kind = 4 ) iband

  real ( kind = 8 ) ch
  external f
  real ( kind = 8 ) f1(neqn)
  real ( kind = 8 ) f2(neqn)
  real ( kind = 8 ) f3(neqn)
  real ( kind = 8 ) f4(neqn)
  real ( kind = 8 ) f5(neqn)
  real ( kind = 8 ) h
  real ( kind = 8 ) s(neqn)
  real ( kind = 8 ) t
  real ( kind = 8 ) y(neqn)
  real ( kind = 8 ) yp(neqn)
  real ( kind = 8 ) magnetic_field(3)


  ch = h / 4.0D+00

  f5(1:neqn) = y(1:neqn) + ch * yp(1:neqn)

  call f (magnetic_field, t + ch, f5, f1, iband)

  ch = 3.0D+00 * h / 32.0D+00

  f5(1:neqn) = y(1:neqn) + ch * ( yp(1:neqn) + 3.0D+00 * f1(1:neqn) )

  call f (magnetic_field, t + 3.0D+00 * h / 8.0D+00, f5, f2, iband)

  ch = h / 2197.0D+00

  f5(1:neqn) = y(1:neqn) + ch * &
  ( 1932.0D+00 * yp(1:neqn) &
  + ( 7296.0D+00 * f2(1:neqn) - 7200.0D+00 * f1(1:neqn) ) &
  )

  call f (magnetic_field, t + 12.0D+00 * h / 13.0D+00, f5, f3, iband)

  ch = h / 4104.0D+00

  f5(1:neqn) = y(1:neqn) + ch * &
  ( &
    ( 8341.0D+00 * yp(1:neqn) - 845.0D+00 * f3(1:neqn) ) &
  + ( 29440.0D+00 * f2(1:neqn) - 32832.0D+00 * f1(1:neqn) ) &
  )

  call f (magnetic_field, t + h, f5, f4, iband)

  ch = h / 20520.0D+00

  f1(1:neqn) = y(1:neqn) + ch * &
  ( &
    ( -6080.0D+00 * yp(1:neqn) &
    + ( 9295.0D+00 * f3(1:neqn) - 5643.0D+00 * f4(1:neqn) ) &
    ) &
  + ( 41040.0D+00 * f1(1:neqn) - 28352.0D+00 * f2(1:neqn) ) &
  )

  call f (magnetic_field, t + h / 2.0D+00, f1, f5, iband)
!
!  Ready to compute the approximate solution at T+H.
!
  ch = h / 7618050.0D+00

  s(1:neqn) = y(1:neqn) + ch * &
  ( &
    ( 902880.0D+00 * yp(1:neqn) &
    + ( 3855735.0D+00 * f3(1:neqn) - 1371249.0D+00 * f4(1:neqn) ) ) &
  + ( 3953664.0D+00 * f2(1:neqn) + 277020.0D+00 * f5(1:neqn) ) &
  )

  return
end  subroutine r8_fehl

