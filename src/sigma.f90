!> Calculate magnetoresistance with R.G.Chambers's formula based on Boltzmann transport
!> Written By QuanSheng Wu (wuquansheng@gmail.com)
!> Thanks Yi Liu for help discussions
!> References : 
!> [1] Electrons in metals and semiconductors, R.G. Chambers,
!> [2] Ab initio investigation of magnetic transport properties by Wannier interpolation, 
!> PHYSICAL REVIEW B 79, 245123 (2009), Yi Liu, Hai-Jun Zhang, and Yugui Yao
!> [3] Magnetoresistance from Fermi surface topology, ShengNan Zhang, QuanSheng Wu, Yi Liu, and Oleg V. Yazyev
!> Phys. Rev. B 99, 035142 (2019)
!> In this subroutine, we can calculate the total conductivity and
!> resistivity under the band-resolved constantly relaxation time
!> approximation.
!> Implemented on Oct. 07 2017
!> uploaded on Sep. 05. 2017
   subroutine sigma_resistivity
      use wmpi
      use para
      implicit none

      integer :: i, ib

      !> Bands crossing Fermi level
      integer :: Nband_Fermi_Level
      integer, allocatable :: bands_fermi_level_temp(:)
      integer, allocatable :: bands_fermi_level(:)

      !> file index
      integer, allocatable  :: myfileindex(:)

      !real(dp) :: KBT ! K_Boltzmann*Temperature in eV
      !real(dp) :: mu ! chemical potential relative to Fermi level in eV

      !> In this method, we can't treat the magnetic field and relaxation time individually
      !> They always come together as Btau. Omega= eB/m*
      !> BTau is in units of Tesla*ps where ps means 10^-12 second.
      !> BTau is in order of 1
      !> the relaxation time tau is in order of 1
      !> For Si, at zero temperature, tau=1ps
      !> For Ge, tau= 0.26ps
      !> For GaAs, tau= 0.48ps
      !> For InAs, tau= 0.08ps 
      !> reference  http://www.iue.tuwien.ac.at/phd/palankovski/node51.html
      !real(dp) :: BTau
      !real(dp) :: BTau_max
      real(dp), allocatable :: BTau_array(:)
      real(dp), allocatable :: mu_array(:)
      real(dp), allocatable :: KBT_array(:)

      real(dp) :: time_start, time_end

      !> conductivity tensor(9, Btau, iband, mu, KBT)
      real(dp), allocatable :: sigma_ohe_tensor(:, :, :, :, :)

      !> plasma frequencies
      real(dp), allocatable :: Plasma_Frequencies(:, :)

      !> file name
      character(40) :: bandname
      character(40) :: sigmafilename
    

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


      !> set Btau array
      allocate(BTau_array(NBTau))
      BTau_array= 0d0
      if (NBTau>1) then
         do i=1, NBTau
            BTau_array(i)= (i-1.0d0)/(NBTau-1)*BTauMax
         enddo
      else
         BTau_array = BTauMax
      endif

      !> set chemical potential array
      allocate(mu_array(OmegaNum))
      mu_array= 0d0
      if (OmegaNum>1) then
         do i=1, OmegaNum
            mu_array(i)= OmegaMin+ (i-1.0d0)/(OmegaNum-1.0d0)*(OmegaMax-OmegaMin)
         enddo
      else
         mu_array = OmegaMin
      endif

      !> set temperature array
      allocate(KBT_array(NumT))
      KBT_array= 0d0
      if (NumT>1) then
         do i=1, NumT
            KBT_array(i)= Tmin+ (i-1.0d0)/(NumT-1.0d0)*(Tmax-Tmin)
         enddo
      else
         KBT_array = Tmin
      endif

      if (cpuid.eq.0) then
         write(stdout, *) ' '
         write(stdout, *)' KBT array in the calculation in unit of Kelvin'
         write(stdout, '(10f8.2)') KBT_array
         write(stdout, *) ' '
      endif
      !> transform from Kelvin to eV
      !> The SI unit of temperature is the kelvin (K), but using the above relation the electron temperature is often expressed in
      !> terms of the energy unit electronvolt (eV). Each kelvin (1 K) corresponds to 8.6173324(78)×10−5 eV; this factor is the ratio
      !> of the Boltzmann constant to the elementary charge. After version 2.6, we 
      !> adopt the atomic unit
      KBT_array= KBT_array*8.6173324E-5*eV2Hartree

      !>> calculate the band resolved conductivity tensor
      !> The tensor is like
      !> xx  xy  xz
      !> yx  yy  yz
      !> zx  zy  zz
      !> sigma1=xx, sigma2=xy, sigma3=xz
      !> sigma4=yx, sigma5=yy, sigma6=yz
      !> sigma7=zx, sigma8=zy, sigma9=zz

      allocate(Plasma_Frequencies(3, Nband_Fermi_Level))
      allocate(sigma_ohe_tensor(9, NBTau, OmegaNum, NumT, Nband_Fermi_Level))
      Plasma_Frequencies= 0d0
      sigma_ohe_tensor= 0d0


      !> file index for different bands
     !do ib=1, Nband_Fermi_Level
     !   outfileindex= outfileindex+ 1
     !   myfileindex(ib)= outfileindex
     !enddo

     !do ib=1, Nband_Fermi_Level
     !   if (cpuid.eq.0) then
     !      write(bandname, '(i10)')bands_fermi_level(ib)
     !      write(sigmafilename, '(3a)')'sigma_band_', trim(adjustl(bandname)), '.dat'
     !      open(unit=myfileindex(ib), file=sigmafilename)
     !      write(myfileindex(ib), '(a,i5)')'# Conductivity tensor for band ', bands_fermi_level(ib)
     !      write(myfileindex(ib), '("#",20a16)')'BTau (T.ps)', 'OmegaTau', 'xx', 'xy', 'xz', 'yx', 'yy', 'yz', 'zx', 'zy','zz'
     !   endif
     !enddo

      time_start= 0d0
      time_end= 0d0

      call sigma_ohe_calc_symm(mu_array, KBT_array, BTau_array, Nband_Fermi_Level, bands_fermi_level, sigma_ohe_tensor)
     
     !if (cpuid.eq.0) then
     !   do ib=1, Nband_Fermi_Level
     !      close(myfileindex(ib))
     !   enddo
     !endif
      !> get the plasma frequency at zero magnetic field
      !> omega_n^2=1/epsilon0*sigma_n/tau_n
      !> \omega_n(x)
      Plasma_Frequencies(1, :)= sqrt(abs(sigma_ohe_tensor(1, 1, 1, 1, :))/epsilon0)&
         *hbar/Echarge
      !> \omega_n(y)
      Plasma_Frequencies(2, :)= sqrt(abs(sigma_ohe_tensor(5, 1, 1, 1, :))/epsilon0)&
         *hbar/Echarge
      !> \omega_n(z)
      Plasma_Frequencies(3, :)= sqrt(abs(sigma_ohe_tensor(9, 1, 1, 1, :))/epsilon0)&
         *hbar/Echarge

    
      if (cpuid.eq.0) then
         write(stdout, *)' The calculation of sigma_OHE is finished. Now, we print out some results:'
         write(stdout, *)' The plasma frequencies at the zero field case for different bands  in units of eV'
         write(stdout, '(a6, 3a10)')' # nband', 'omega_x', 'omega_y', 'omega_z'  
         do i=1, Nband_Fermi_Level
            write(stdout, '(i6, 3E16.5)') i, Plasma_Frequencies(:, i)
         enddo
      endif
      

    
      outfileindex= outfileindex+ 1
      !> write script for gnuplot
      if (cpuid==0) then
         open(unit=outfileindex, file='sigma.gnu')
         write(outfileindex, '(a)')"set encoding iso_8859_1"
         write(outfileindex, '(a)')'set terminal  postscript enhanced color font ",24" '
         write(outfileindex, '(a)')"set output 'sigma.eps'"
         write(outfileindex, '(a)')'set border lw 2'
         write(outfileindex, '(a)')'set autoscale fix'
         write(outfileindex, '(a)')'unset key'
         write(outfileindex, '(a)')"set ylabel {/Symbol s}/{/Symbol t} (1/{/Symbol W}/m/s)"
         write(outfileindex, '(a)')"set xlabel B{/Symbol t} (T.ps)"
         do ib=1, Nband_Fermi_Level-1
            write(bandname, '(i10)')bands_fermi_level(ib)
            write(sigmafilename, '(3a)')'sigma_band_',trim(adjustl(bandname)), '.dat'
            write(outfileindex, '(3a)')"plot '",trim(adjustl(sigmafilename)), "' u 1:3 w lp lw 3 lc rgb 'black', \"
         enddo
         write(bandname, '(i10)')bands_fermi_level(Nband_Fermi_Level)
         write(sigmafilename, '(3a)')'sigma_band_',trim(adjustl(bandname)), '.dat'
         write(outfileindex, '(3a)')"     '",trim(adjustl(sigmafilename)), "' u 1:3 w lp lw 3 lc rgb 'black'"
         close(outfileindex)
      endif

      if (cpuid.eq.0) write(stdout, '(a)') ' <<< The conductivity calculation finished'
      if (cpuid.eq.0) write(stdout, '(a)') ' '

      return
   end subroutine sigma_resistivity

   !> We calculate howmany and which bands cross the Fermi level
   subroutine get_bands_cross_fermilevel(Nband_Fermi_Level, bands_fermi_level)
      use wmpi
      use para
      implicit none

      integer, intent(inout) :: Nband_Fermi_Level
      integer, intent(inout) :: bands_fermi_level(Nband_Fermi_Level)

      logical :: lower, higher
      integer :: ik1, nkx, ik
      integer :: ik2, nky, knv3
      integer :: ik3, nkz, ib
      integer :: ierr

      real(dp) :: k(3)
      real(dp), allocatable :: W(:)
      real(dp), allocatable :: eigval(:, :)
      real(dp), allocatable :: eigval_mpi(:, :)
      complex(dp), allocatable :: Hk(:, :)

      nkx= 10
      nky= 10
      nkz= 10
      knv3= nkx*nky*nkz

      allocate(W(Num_wann))
      allocate(Hk(Num_wann, Num_wann))
      allocate(eigval(Num_wann, knv3))
      allocate(eigval_mpi(Num_wann, knv3))
      W=0d0
      Hk= 0d0
      eigval= 0d0
      eigval_mpi= 0d0

      do ik=1+cpuid, knv3, num_cpu
         if (cpuid.eq.0) write(stdout, '(a, i18, "      /", i18)') 'ik/knv3', ik, knv3
   
         ik1= (ik-1)/(nky*nkz)+1
         ik2= ((ik-1-(ik1-1)*nky*nkz)/nkz)+1
         ik3= (ik-(ik2-1)*nkz- (ik1-1)*nky*nkz)
         k(1)= (ik1-1)/dble(nkx)  
         k(2)= (ik2-1)/dble(nky) 
         k(3)= (ik3-1)/dble(nkz)
   
         W= 0d0
         call ham_bulk_atomicgauge(k, Hk)
         call eigensystem_c( 'N', 'U', Num_wann ,Hk, W)
         eigval_mpi(:, ik)= W
      enddo  ! ik

#if defined (MPI)
     call mpi_allreduce(eigval_mpi, eigval,size(eigval),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
     eigval= eigval_mpi 
#endif

      !> find the lowest band which cross the fermi level
      Nband_Fermi_Level= 0
      do ib=1, Num_wann
         lower= .FALSE.
         do ik=1, knv3
            if (eigval(ib, ik)< 0d0) then
               lower =.TRUE.
               exit
            endif
         enddo ! ik

         higher= .FALSE.
         do ik=1, knv3
            if (eigval(ib, ik)> 0d0) then
               higher =.TRUE.
               exit
            endif
         enddo ! ik

         !> check if there are bands lower and higher than EF
         if (lower .and. higher) then
            Nband_Fermi_Level= Nband_Fermi_Level+ 1
            bands_fermi_level(Nband_Fermi_Level)= ib
         endif
      enddo

      if (cpuid.eq.0) then
         write(stdout, '(a)') ' '
         write(stdout, '(a)') ' >> Runing in subroutine get_bands_cross_fermilevel'
         write(stdout, '(a, i5, a)')' There are ', Nband_Fermi_Level, ' bands crossing the fermi level'
         write(stdout, '(a, 100i5)')' Those bands are ', bands_fermi_level(1:Nband_Fermi_Level)
         write(stdout, '(a)') ' '
      endif

      return
   end subroutine get_bands_cross_fermilevel
