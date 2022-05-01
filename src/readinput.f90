subroutine readinput
   ! Read in the control paramters from wt.in,
   ! and set default values if not specified in the wt.in
   !
   ! Constructed on 4/22/2010 by QS.Wu

   use wmpi
   use para
   implicit none

   character*12 :: fname='wt.in'
   character*25 :: char_temp
   character*256 :: inline
   logical ::  exists, lfound
   real(dp) :: cell_volume, cell_volume2

   integer  :: i, j, ia, ik, iq, io, it, idummy, ig, NN, stat, istart, iend
   integer  :: NumberOfspinorbitals, NumberOfspinorbitals_unfold

   integer, allocatable :: iarray_temp(:)

   real(dp) :: t1, temp
   real(dp) :: pos(3), k1(3), k2(3), k(3), kstart(3), kend(3)
   real(dp) :: R1(3), R2(3), R3(3), Rt(3), Rt2(3)
   real(dp), external :: norm, angle

   real(dp), allocatable :: mass_temp(:)
   real(dp), allocatable :: Born_Charge_temp(:, :, :)

   inquire(file=fname,exist=exists)
   if (exists)then
      if(cpuid==0)write(stdout,*) '  '
      if(cpuid==0)write(stdout,*) '>>>Read some paramters from wt.in'
      open(unit=1001,file=fname,status='old')
   else
      if(cpuid==0)write(stdout,*)'file' ,fname, 'dosnot exist'
      stop
   endif

!===============================================================================================================!
!  TB_FILE namelist
!===============================================================================================================!

   Particle='electron'
   Package= 'VASP'
   KPorTB = 'TB'
   Is_HrFile= .TRUE.
   Is_Sparse_Hr= .FALSE.
   Is_Sparse   = .FALSE.
   Use_ELPA= .FALSE.
   vef=0d0
   read(1001, TB_FILE, iostat= stat)
   if (stat/=0) then
      Hrfile='wannier90_hr.dat'
      Particle='electron'
      inquire(file='wannier90_hr.dat',exist=exists)

      backspace(1001)
      read(1001,fmt='(A)') inline
      write(*,'(A)') &
         '>>> ERROR : Invalid line in namelist TB_FILE : '//trim(inline)

      if (.not.exists) stop "ERROR>> TB_FIlE namelist should be given or wannier90_hr.dat should exist"

   endif
   if(cpuid==0)write(stdout,'(1x, a, a6, a)')"You are using : ", KPorTB, " model"
   Is_Sparse_Hr= (Is_Sparse_Hr.or.Is_Sparse)
   if(cpuid==0) then
      if(Is_HrFile) then
         write(stdout,'(1x, a)')"Tight-binding Hamiltonian FROM: Hr File"
         write(stdout,'(1x, a, L2)')"Is_Sparse_Hr= ", Is_Sparse_Hr
         if(Is_Sparse_Hr) write(stdout,'(1x, a)')"Tight-binding Hamiltonian FROM: Sparse Hr File"
      else
         write(stdout,'(1x, a)')"Tight-binding Hamiltonian FROM: fitting"
      end if
   end if
   if(cpuid==0)write(stdout,'(1x, a, a25)')"Tight-binding Hamiltonian filename: ",Hrfile
   if(cpuid==0)write(stdout,'(1x, a, a25)')"System of particle: ", Particle
   if(cpuid==0)write(stdout,'(1x, a, a25)')"Tight-binding Hamiltonian obtained from : ",Package

   if (index(Particle, 'electron')==0 .and. index(Particle, 'phonon')==0 &
      .and. index(Particle, 'photon')==0) then
      write(stdout, *)' ERROR: Particle shoule equal either "electron"' , &
         '"phonon", or "photon"'
      stop
   endif

!===============================================================================================================!
!> CONTROL namelist
!===============================================================================================================!

   BulkBand_calc         = .FALSE.
   BulkBand_line_calc    = .FALSE.
   BulkBand_unfold_line_calc    = .FALSE.
   Landaulevel_unfold_line_calc    = .FALSE.
   BulkBand_unfold_plane_calc    = .FALSE.
   QPI_unfold_plane_calc    = .FALSE.
   BulkFatBand_calc      = .FALSE.
   BulkBand_plane_calc   = .FALSE.
   BulkBand_cube_calc    = .FALSE.
   BulkFS_calc           = .FALSE.
   BulkFS_Plane_calc     = .FALSE.
   BulkFS_plane_stack_calc     = .FALSE.
   BulkGap_cube_calc     = .FALSE.
   BulkGap_plane_calc    = .FALSE.
   SlabBand_calc         = .FALSE.
   SlabBandWaveFunc_calc = .FALSE.
   WireBand_calc         = .FALSE.
   SlabSS_calc           = .FALSE.
   SlabArc_calc          = .FALSE.
   SlabQPI_calc          = .FALSE.
   SlabQPI_kpath_calc    = .FALSE.
   SlabQPI_kplane_calc   = .FALSE.
   SlabSpintexture_calc  = .FALSE.
   BulkSpintexture_calc  = .FALSE.
   wanniercenter_calc    = .FALSE.
   Z2_3D_calc            = .FALSE.
   Chern_3D_calc         = .FALSE.
   WeylChirality_calc    = .FALSE.
   NLChirality_calc    = .FALSE.
   BerryPhase_calc       = .FALSE.
   BerryCurvature_calc   = .FALSE.
   BerryCurvature_EF_calc   = .FALSE.
   BerryCurvature_plane_selectedbands_calc = .FALSE.
   BerryCurvature_Cube_calc   = .FALSE.
   BerryCurvature_slab_calc = .FALSE.
   Berrycurvature_kpath_EF_calc = .FALSE.
   BerryCurvature_kpath_Occupied_calc = .FALSE.
   MirrorChern_calc      = .FALSE.
   Dos_calc              = .FALSE.
   JDos_calc             = .FALSE.
   EffectiveMass_calc    = .FALSE.
   FindNodes_calc        = .FALSE.
   LOTO_correction       = .FALSE.
   Boltz_OHE_calc        = .FALSE.
   Boltz_Berry_correction= .FALSE.
   AHC_calc              = .FALSE.
   Hof_Butt_calc         = .FALSE.
   LandauLevel_k_calc    = .FALSE.
   LandauLevel_B_calc    = .FALSE.
   LandauLevel_wavefunction_calc    = .FALSE.
   OrbitalTexture_calc   = .FALSE.
   OrbitalTexture_3D_calc  = .FALSE.
   Fit_kp_calc             = .FALSE.
   DMFT_MAG_calc         = .FALSE.
   Symmetry_Import_calc  = .FALSE.
   LanczosSeqDOS_calc    = .FALSE.
   LandauLevel_kplane_calc = .FALSE.
   LandauLevel_k_dos_calc = .FALSE.
   LandauLevel_B_dos_calc = .FALSE.
   Translate_to_WS_calc    = .FALSE.
   FermiLevel_calc    = .FALSE.
   w3d_nested_calc =.false.
   ChargeDensity_selected_bands_calc= .FALSE.
   ChargeDensity_selected_energies_calc= .FALSE.

   read(1001, CONTROL, iostat=stat)
   SlabQPI_kplane_calc= SlabQPI_kplane_calc.or.SlabQPI_calc

   if (stat/=0) then
      write(*, *)"ERROR: namelist CONTROL should be set"
      write(*, *)"You should set one of these functions to be T."
      write(*, *)"And please make sure that the spelling are correct."
      write(*, *)"BulkBand_calc, BulkBand_plane_calc, BulkFS_calc"
      write(*, *)"BulkBand_line_calc, BulkBand_cube_calc"
      write(*, *)"Landaulevel_unfold_line_calc, "
      write(*, *)"BulkBand_unfold_line_calc, "
      write(*, *)"BulkBand_unfold_plane_calc, "
      write(*, *)"QPI_unfold_plane_calc, "
      write(*, *)"BulkFatBand_calc, "
      write(*, *)"BulkGap_cube_calc,BulkGap_plane_calc"
      write(*, *)"SlabBand_calc,SlabBandWaveFunc_calc"
      write(*, *)"WireBand_calc,SlabSS_calc,SlabArc_calc "
      write(*, *)"SlabQPI_calc"
      write(*, *)"SlabQPI_kpath_calc"
      write(*, *)"SlabQPI_kplane_calc"
      write(*, *)"SlabSpintexture,wanniercenter_calc"
      write(*, *)"BerryPhase_calc,BerryCurvature_calc, BerryCurvature_EF_calc"
      write(*, *)"Berrycurvature_kpath_EF_calc, BerryCurvature_kpath_Occupied_calc"
      write(*, *)"BerryCurvature_slab_calc, BerryCurvature_Cube_calc"
      write(*, *)"Dos_calc, JDos_calc, FindNodes_calc"
      write(*, *)"BulkFS_plane_calc"
      write(*, *)"BulkFS_plane_stack_calc"
      write(*, *)"Z2_3D_calc"
      write(*, *)"Chern_3D_calc"
      write(*, *)"MirrorChern_calc"
      write(*, *)"WeylChirality_calc"
      write(*, *)"NLChirality_calc"
      write(*, *)"LOTO_correction"
      write(*, *)"AHC_calc"
      write(*, *)"Hof_Butt_calc"
      write(*, *)"Boltz_OHE_calc"
      write(*, *)"Boltz_Berry_correction"
      write(*, *)"DMFT_MAG_calc"
      write(*, *)"Fit_kp_calc"
      write(*, *)"OrbitalTexture_calc"
      write(*, *)"OrbitalTexture_3D_calc"
      write(*, *)"LandauLevel_wavefunction_calc"
      write(*, *)"LandauLevel_k_calc"
      write(*, *)"LandauLevel_kplane_calc"
      write(*, *)"LandauLevel_B_calc"
      write(*, *)"Hof_Butt_calc"
      write(*, *)"LandauLevel_k_dos_calc"
      write(*, *)"LandauLevel_B_dos_calc"
      write(*, *)"Translate_to_WS_calc"
      write(*, *)"FermiLevel_calc"
      write(*, *)"ChargeDensity_selected_energies_calc"
      write(*, *)"ChargeDensity_selected_bands_calc"
      write(*, *)"The default Vaule is F"

      backspace(1001)
      read(1001,fmt='(A)') inline
      write(*,'(A)') &
         '>>> ERROR : Invalid line in namelist CONTROL : '//trim(inline)


      stop
   endif

   !> In order to be compatiable with the old version, we keep the bulkband_calc.
   BulkBand_line_calc= BulkBand_line_calc.or.BulkBand_calc
   BulkBand_unfold_line_calc= BulkBand_unfold_line_calc.or.Landaulevel_unfold_line_calc

   if (MirrorChern_calc) Symmetry_Import_calc = .true.

   !> control parameters
   if (cpuid==0) then
      write(stdout, *) "  "
      write(stdout, *) ">>>Control parameters: "
      write(stdout, *) "BulkBand_line_calc                : ",  BulkBand_line_calc
      write(stdout, *) "BulkBand_plane_calc               : ",  BulkBand_plane_calc
      write(stdout, *) "Landaulevel_unfold_line_calc      : ",  Landaulevel_unfold_line_calc
      write(stdout, *) "BulkBand_unfold_line_calc         : ",  BulkBand_unfold_line_calc
      write(stdout, *) "BulkBand_unfold_plane_calc        : ",  BulkBand_unfold_plane_calc
      write(stdout, *) "QPI_unfold_plane_calc             : ",  QPI_unfold_plane_calc
      write(stdout, *) "BulkFatBand_calc                  : ",  BulkFatBand_calc
      write(stdout, *) "BulkBand_cube_calc                : ",  BulkBand_cube_calc
      write(stdout, *) "BulkFS_calc                       : ",  BulkFS_calc
      write(stdout, *) "BulkFS_Plane_calc                 : ",  BulkFS_Plane_calc
      write(stdout, *) "BulkFS_plane_stack_calc           : ",  BulkFS_plane_stack_calc
      write(stdout, *) "BulkGap_cube_calc                 : ",  BulkGap_cube_calc
      write(stdout, *) "BulkGap_plane_calc                : ",  BulkGap_plane_calc
      write(stdout, *) "SlabBand_calc                     : ",  SlabBand_calc
      write(stdout, *) "SlabBandWaveFunc_calc             : ",  SlabBandWaveFunc_calc
      write(stdout, *) "SlabSS_calc                       : ",  SlabSS_calc
      write(stdout, *) "SlabArc_calc                      : ",  SlabArc_calc
      write(stdout, *) "SlabSpintexture_calc              : ",  SlabSpintexture_calc
      write(stdout, *) "wanniercenter_calc                : ", wanniercenter_calc
      write(stdout, *) "BerryPhase_calc                   : ", BerryPhase_calc
      write(stdout, *) "BerryCurvature_calc               : ", BerryCurvature_calc
      write(stdout, *) "BerryCurvature_EF_calc            : ", BerryCurvature_EF_calc
      write(stdout, *) "BerryCurvature_kpath_EF_calc      : ", BerryCurvature_kpath_EF_calc
      write(stdout, *) "BerryCurvature_kpath_Occupied_calc: ", BerryCurvature_kpath_Occupied_calc
      write(stdout, *) "BerryCurvature_Cube_calc          : ", BerryCurvature_Cube_calc
      write(stdout, *) "BerryCurvature_slab_calc          : ", BerryCurvature_slab_calc
      write(stdout, *) "Dos_calc                          : ",  DOS_calc
      write(stdout, *) "Z2_3D_calc                        : ",  Z2_3D_calc
      write(stdout, *) "WeylChirality_calc                : ",  WeylChirality_calc
      write(stdout, *) "NLChirality_calc                  : ",  NLChirality_calc
      write(stdout, *) "Chern_3D_calc                     : ",  Chern_3D_calc
      write(stdout, *) "MirrorChern_calc                  : ",  MirrorChern_calc
      write(stdout, *) "JDos_calc                         : ",  JDOS_calc
      write(stdout, *) "FindNodes_calc                    : ",  FindNodes_calc
      write(stdout, *) "EffectiveMass_calc                : ", EffectiveMass_calc
      write(stdout, *) "AHC_calc                          : ", AHC_calc
      write(stdout, *) "Boltz_OHE_calc                    : ", Boltz_OHE_calc
      write(stdout, *) "Boltz_Berry_correction            : ", Boltz_Berry_correction
      write(stdout, *) "LOTO_correction                   : ", LOTO_correction
      write(stdout, *) "OrbitalTexture_calc               : ", OrbitalTexture_calc
      write(stdout, *) "OrbitalTexture_3D_calc            : ", OrbitalTexture_3D_calc
      write(stdout, *) "LandauLevel_k_calc                : ", LandauLevel_k_calc
      write(stdout, *) "LandauLevel_B_calc                : ", LandauLevel_B_calc
      write(stdout, *) "LandauLevel_wavefunction_calc     : ", LandauLevel_wavefunction_calc
      write(stdout, *) "Fit_kp_calc                       : ", Fit_kp_calc
      write(stdout, *) "DMFT_MAG_calc                     : ", DMFT_MAG_calc
      write(stdout, *) "Translate_to_WS_calc              : ", Translate_to_WS_calc
      write(stdout, *) "LandauLevel_kplane_calc           : ", LandauLevel_kplane_calc
      write(stdout, *) "LandauLevel_k_dos_calc            : ", LandauLevel_k_dos_calc
      write(stdout, *) "LandauLevel_B_dos_calc            : ", LandauLevel_B_dos_calc
      write(stdout, *) "FermiLevel_calc                   : ", FermiLevel_calc
      write(stdout, *) "Symmetry_Import_calc              : ", Symmetry_Import_calc
      write(stdout, *) "ChargeDensity_selected_bands_calc : ", ChargeDensity_selected_bands_calc
      write(stdout, *) "ChargeDensity_selected_energies_calc : ", ChargeDensity_selected_energies_calc
   endif

!===============================================================================================================!
!> SYSTEM namelist
!===============================================================================================================!

   !> set system parameters by default
   Nslab= 10
   Nslab1= 1
   Nslab2= 1
   Numoccupied = 0
   Ntotch = 0
   SOC = 0
   SOC_in = 0
   E_FERMI = 0d0

   !> By default magnetic field is zero
   Bx = 0d0
   By = 0d0
   Bz = 0d0
   
   Bmagnitude = 0d0
   Btheta = -99999d0
   Bphi = -99999d0
   surf_onsite = 0d0

   !> By default, we don't add zeeman field
   Add_Zeeman_Field = .FALSE.

   !> by default, g-factor is 2
   Effective_gfactor = 2d0
   Zeeman_energy_in_eV = 0d0

   !> by default, Electric_field_in_eVpA=0
   Electric_field_in_eVpA= 0d0
   Symmetrical_Electric_field_in_eVpA= 0d0
   Inner_symmetrical_Electric_Field= .False.

   !> by default, Vacuum_thickness_in_Angstrom= 20 Angstrom
   Vacuum_thickness_in_Angstrom = 20d0

   !> read system parameters from file
   read(1001, SYSTEM, iostat=stat)
   if (stat/=0) then
      write(*, *)"ERROR: namelist SYSTEM is wrong and should be set correctly"

      backspace(1001)
      read(1001,fmt='(A)') inline
      write(*,'(A)') &
         '>>> ERROR : Invalid line in namelist SYSTEM : '//trim(inline)

      stop
   endif
   SOC_in=SOC

   if (SOC == -1) then
      write(*, *)"ERROR: you should set SOC in namelist SYSTEM correctly"
      stop
   endif

   if (Numoccupied == 0) then
      if (Z2_3D_calc.or.Chern_3D_calc.or.BulkFS_calc.or.BulkFS_Plane_calc &
      .or.BulkFS_plane_stack_calc.or.BulkGap_plane_calc.or.WannierCenter_calc.or.&
      BerryPhase_calc.or.BerryCurvature_EF_calc.or.BerryCurvature_calc.or.&
      BerryCurvature_plane_selectedbands_calc.or.BerryCurvature_slab_calc.or.&
      MirrorChern_calc.or.WeylChirality_calc.or.NLChirality_calc.or.&
      FindNodes_calc) then
         write(*, *)"ERROR: you should set Numoccupied in namelist SYSTEM correctly"
         stop
      else 
         Numoccupied = 1
      endif
   endif


   if (abs(Ntotch) <eps3) then
      if (SOC>0) then
         Ntotch = Numoccupied
      else
         Ntotch = Numoccupied*2
      endif
   endif

   if (.not.Add_Zeeman_Field) then
      Zeeman_energy_in_eV = 0d0
      Effective_gfactor = 0d0
   endif

   if (cpuid==0) then
      write(stdout, *) "  "
      write(stdout, *) ">>>System parameters: "
      write(stdout, '(1x, a, i6 )')"NumSlabs :", Nslab
      write(stdout, '(1x, a, i6)')"Nslab1 for ribbon  :", Nslab1
      write(stdout, '(1x, a, i6)')"Nslab2 for ribbon  :", Nslab2
      write(stdout, '(1x, a, i6)')"Number of Occupied bands:", NumOccupied
      write(stdout, '(1x, a, f12.6)')"Number of total electrons:", Ntotch
      write(stdout, '(1x, a, i6)')"With SOC or not in Hrfile:", SOC
      write(stdout, '(1x, a, 3f16.6)')"Fermi energy (eV) :", E_FERMI
      write(stdout, '(1x, a, 3f16.6)')"surf_onsite (eV): ", surf_onsite
      write(stdout, '(1x, a, L)')"Add_Zeeman_Field: ", Add_Zeeman_Field
      write(stdout, '(1x, a, 3f16.6)')"Zeeman_energy_in_eV (eV): ",  Zeeman_energy_in_eV
      write(stdout, '(1x, a, 3f16.6)')"Electric_field_in_eVpA (eV/Angstrom): ",  Electric_field_in_eVpA
      write(stdout, '(1x, a, 3f16.6)')"Symmetrical_Electric_field_in_eVpA (eV/Angstrom): ",  Symmetrical_Electric_field_in_eVpA
      write(stdout, '(1x, a, L)')"Inner_symmetrical_Electric_Field: ",  Inner_symmetrical_Electric_Field
      write(stdout, '(1x, a, i6 )')"ijmax :", ijmax
   endif

   if (cpuid==0) then
      write(stdout, *) "  "
      write(stdout, "(1x,a)") "**Notes**: There are two ways to specify magnetic field."
      write(stdout, "(1x,a)") "1. specify Bmagnitude, Btheta, Bphi"
      write(stdout, "(1x,a)") "2. specify Bx, By, Bz"
      write(stdout, "(1x,a)") "Bmagnitude, Bx, By, Bz are real numbers in unit of Tesla. and Bx, By, Bz "
      write(stdout, "(1x,a)") "are in the cartesian coordinates. Btheta is the angle with z axial "
      write(stdout, "(1x,a)") " and Bphi is the angle with respect to x axial in the x-y plane"
      write(stdout, "(1x,a)") " Btheta is in [0, 180], Bphi is in [0, 360)."
      write(stdout, "(1x,a)") "If you specify both of them together, we will choose the first one."
      write(stdout, "(1x,a)") "If choose the first one, but not specify Btheta, Bphi, then "
      write(stdout, "(1x,a)") "by default we set Btheta=0, Bphi=0 which means B is along z direction."
   endif

   !> check if Bmagnitude is specified in the input.dat/wt.in
   if (abs(Bmagnitude)>eps6 .and. (abs(Bx)+abs(By)+abs(Bz))>eps6) then
      if (cpuid==0) then
         write(stdout, *) " "
         write(stdout, *) " Warning: You specify Bmagnitude and Bx, By, Bz at the same time "
         write(stdout, *) " in the SYSTEM namelist in the wt.in/input.dat."
         write(stdout, *) " However, we will only take Bmagnitude and Btheta, Bphi but discard Bx,By,Bz. "
         write(stdout, *) " Bx,By,Bz will be calculated as Bx=Bmagnitude*sin(Btheta*pi/180)*cos(bphi/180*pi). "
         write(stdout, *) "  By=Bmagnitude*sin(Btheta*pi/180)*sin(bphi/180*pi), "
         write(stdout, *) "  Bz=cos(btheta/180*pi). "
      endif
   endif

   if (abs(Bmagnitude)<eps6 .and. (abs(Bx)+abs(By)+abs(Bz))<eps6) then
      if (cpuid==0) then
         write(stdout, *) " "
         write(stdout, *) " Warning: You didn't specify the magnitude of magnetic field."
         write(stdout, *) " "
      endif
   endif


   if (abs(Bmagnitude)>eps6) then
      if (abs(Btheta+99999d0)<eps6) then
         Btheta=0d0
      endif
      if (abs(Bphi+99999d0)<eps6) then
         Bphi=0d0
      endif

      Bx=Bmagnitude*sin(Btheta/180d0*pi)*cos(Bphi/180d0*pi)
      By=Bmagnitude*sin(Btheta/180d0*pi)*sin(Bphi/180d0*pi)
      Bz=Bmagnitude*cos(Btheta/180d0*pi)

      Bdirection(1)=Bx/Bmagnitude
      Bdirection(2)=By/Bmagnitude
      Bdirection(3)=Bz/Bmagnitude

      if (cpuid==0)then
         write(stdout, '(1x, a,  f16.6)')">>You specified the magnitude of magnetic field, "
         write(stdout, '(1x, a,  f16.6)')"So we will take this option but discard the setting of Bx, By, Bz"
         write(stdout, '(1x, a,  f16.6)')"Bmagnitude (in Tesla) :", Bmagnitude
         write(stdout, '(1x, a, 3f16.6)')"Btheta, Bphi :", Btheta, Bphi
         write(stdout, '(1x, a, 3f16.6)')"Bx, By, Bz (in Tesla) :", Bx, By, Bz
         write(stdout, '(1x, a, 3f16.6)')"B direction cosines :", Bdirection
      endif
   else
      if (abs(Bx)+abs(By)+abs(Bz)>eps6) then
         if (cpuid==0)then
            write(stdout, '(1x, a,  f16.6)')">>You did not specified the magnitude of magnetic field, "
            write(stdout, '(1x, a,  f16.6)')" but set Bx, By, Bz. So we will set the magnetic direction"
            write(stdout, '(1x, a,  f16.6)')" from Bx, By, Bz but discard the settings of Btheta, Bphi."
         endif
         Bmagnitude= sqrt(Bx*Bx+By*By+Bz*Bz)  !> in Tesla
         Bdirection(1)=Bx/Bmagnitude
         Bdirection(2)=By/Bmagnitude
         Bdirection(3)=Bz/Bmagnitude

         Btheta= acos(Bdirection(3))*180d0/pi
         if (abs(Bx)+abs(By)<eps6) then
            Bphi=0d0
         else
            Bphi= atan2(Bdirection(2),Bdirection(1))*180d0/pi
         endif
      else   
         if (cpuid==0)then
            write(stdout, '(1x, a,  f16.6)')">>You did not specified the magnitude of magnetic field, "
            write(stdout, '(1x, a,  f16.6)')" and also didn't set Bx, By, Bz. So we will set the magnetic direction"
            write(stdout, '(1x, a,  f16.6)')" from the settings of Btheta, Bphi."
            write(stdout, '(1x, a,  f16.6)')" If you even didn't set Btheta, Bphi, we will take the default value, "
            write(stdout, '(1x, a,  f16.6)')" Btheta=0, Bphi=0 which means B is along z direction"
         endif

         if (abs(Btheta+99999d0)<eps6) then
            Btheta=0d0
            if (cpuid==0)then
               write(stdout, '(1x, a,  f16.6)')">>You didn't set Btheta, and we set it to 0."
            endif
         else
            if (cpuid==0)then
               write(stdout, '(1x, a,  f16.6)')">>You set Btheta, and we will take it."
            endif
         endif
         if (abs(Bphi+99999d0)<eps6) then
            Bphi=0d0
            if (cpuid==0)then
               write(stdout, '(1x, a,  f16.6)')">>You didn't set Bphi, and we set it to 0."
            endif
         else
            if (cpuid==0)then
               write(stdout, '(1x, a,  f16.6)')">>You set Bphi, and we will take it."
            endif
         endif
         Bdirection(1)=sin(Btheta/180d0*pi)*cos(Bphi/180d0*pi)
         Bdirection(2)=sin(Btheta/180d0*pi)*sin(Bphi/180d0*pi)
         Bdirection(3)=cos(Btheta/180d0*pi)
         Bmagnitude= 1d0
         Bx=Bmagnitude*sin(Btheta/180d0*pi)*cos(Bphi/180d0*pi)
         By=Bmagnitude*sin(Btheta/180d0*pi)*sin(Bphi/180d0*pi)
         Bz=Bmagnitude*cos(Btheta/180d0*pi)
      endif

      if (cpuid==0)then
         write(stdout, '(1x, a,  f16.6)')"Bmagnitude (in Tesla) :", Bmagnitude
         write(stdout, '(1x, a, 3f16.6)')"Btheta, Bphi :", Btheta, Bphi
         write(stdout, '(1x, a, 3f16.6)')"Bx, By, Bz (in Tesla) :", Bx, By, Bz
         write(stdout, '(1x, a, 3f16.6)')"B direction cosines :", Bdirection
      endif
   endif

   !> change the unit of magnetic field from Tesla to atomic unit
   !> 1 atomic unit = hbar/(e*a0^2) Tesla
   !> 1 Tesla= (e*a0^2)/hbar a.u.
   Bmagnitude= Bmagnitude*Echarge*bohr_radius**2/hbar
   Bx=Bmagnitude*sin(Btheta/180d0*pi)*cos(Bphi/180d0*pi)
   By=Bmagnitude*sin(Btheta/180d0*pi)*sin(Bphi/180d0*pi)
   Bz=Bmagnitude*cos(Btheta/180d0*pi)

!===============================================================================================================!
!> PARAMETERS namelist
!===============================================================================================================!


   !> set up parameters for calculation
   E_arc = 0.0d0
   Eta_Arc= 0.001d0
   OmegaNum = 100
   OmegaNum_unfold = 0
   OmegaMin = -1d0
   OmegaMax =  1d0
   Nk1 = 10
   Nk2 = 10
   Nk3 = 1 
   NP  = 2
   Gap_threshold= 0.01d0
   Tmin = 100.  ! in Kelvin
   Tmax = 100.  ! in Kelvin
   NumT= 1
   NBTau = 1
   BTauNum = 1
   Nslice_BTau_Max = 5000
   BTauMax = 0d0
   Rcut = 999999d0
   Magp= 1
   Magq= 0
   Magp_min=0
   Magp_max=0   
   wcc_calc_tol= 0.08
   wcc_neighbour_tol= 0.3
   NumLCZVecs=400
   NumRandomConfs=1
   NumSelectedEigenVals=0 
   Beta= 100
   Relaxation_Time_Tau= 1d0  ! in ps
   topsurface_atom_index= 0


   !> by default, we only project on atoms for a given wave function
   projection_weight_mode = "NORMAL"


   read(1001, PARAMETERS, iostat= stat)
   if (Magp<1) Magp= 0
   if (Magp_max<1) Magp_max= Magp
   if (Magq==0) Magq= Nslab
   if (Is_Sparse_Hr) then
      if (OmegaNum_unfold==0) OmegaNum_unfold= 4*OmegaNum
   else
      if (OmegaNum_unfold==0) OmegaNum_unfold= 200
   endif

   if (stat>0) then

      backspace(1001)
      read(1001,fmt='(A)') inline
      write(*,'(A)') &
         '>>> ERROR : Invalid line in namelist PARAMETERS : '//trim(inline)

   endif

   NBTau= max(NBTau, BTauNum)
  
   projection_weight_mode= upper(projection_weight_mode)
   if (cpuid==0) then
      write(stdout, *) "  "
      write(stdout, *) ">>>calculation parameters : "
      write(stdout, '(1x, a, f16.5)')'E_arc : ', E_arc
      write(stdout, '(1x, a, f16.5)')'Eta_arc : ', Eta_arc
      write(stdout, '(1x, a, f16.5)')'Gap_threshold', Gap_threshold
      write(stdout, '(1x, a, f16.5)')'OmegaMin : ', OmegaMin
      write(stdout, '(1x, a, f16.5)')'OmegaMax : ', OmegaMax
      write(stdout, '(1x, a, i6   )')'OmegaNum : ', OmegaNum
      write(stdout, '(1x, a, i6   )')'OmegaNum_unfold : ', OmegaNum_unfold
      write(stdout, '(1x, a, i6   )')'Nk1 : ', Nk1
      write(stdout, '(1x, a, i6   )')'Nk2 : ', Nk2
      write(stdout, '(1x, a, i6   )')'Nk3 : ', Nk3
      write(stdout, '(1x, a, i6   )')'NP number of principle layers  : ', Np
      write(stdout, '(1x, a, f16.5)')'Tmin(Kelvin)  : ', Tmin
      write(stdout, '(1x, a, f16.5)')'Tmax(Kelvin)  : ', Tmax
      write(stdout, '(1x, a, i6   )')'NumT  : ', NumT
      write(stdout, '(1x, a, i6   )')'NBTau  : ', NBTau
      write(stdout, '(1x, a, f16.5)')'Beta  : ', Beta
      write(stdout, '(1x, a, i6   )')'Nslice_BTau_Max  : ', Nslice_BTau_Max
      write(stdout, '(1x, a, f16.5)')'BTauMax(Tesla.ps)', BTauMax
      write(stdout, '(1x, a, f16.5)')'Relaxation_Time_Tau (ps)', Relaxation_Time_Tau
      write(stdout, '(1x, a, f16.5)')'Rcut', Rcut
      write(stdout, '(1x, a, i16  )')'Magp', Magp
      write(stdout, '(1x, a, i16  )')'Magp_min', Magp_min
      write(stdout, '(1x, a, i16  )')'Magp_max', Magp_max
      write(stdout, '(1x, a, f16.5)')'wcc_calc_tol', wcc_calc_tol
      write(stdout, '(1x, a, f16.5)')'wcc_neighbour_tol', wcc_neighbour_tol
      write(stdout, '(1x, a, i6   )')'NumLCZVecs', NumLCZVecs
      write(stdout, '(1x, a, i6   )')'NumSelectedEigenVals', NumSelectedEigenVals
      write(stdout, '(1x, a, i6   )')'NumRandomConfs:', NumRandomConfs
      write(stdout, '(1x, a, a    )')'Projection weight mode:', projection_weight_mode
      write(stdout, '(1x, a, i8   )')'The size of magnetic supercell is Magq= :', Magq
   endif

   !> changed to atomic units
   E_arc= E_arc*eV2Hartree
   Eta_Arc = Eta_Arc*eV2Hartree
   OmegaMin= OmegaMin*eV2Hartree
   OmegaMax= OmegaMax*eV2Hartree
   Gap_threshold= Gap_threshold*eV2Hartree
   Rcut= Rcut*Ang2Bohr

   !> change the unit of relaxtion time from ps to atomic unit
   Relaxation_Time_Tau= Relaxation_Time_Tau*1E-12/Time_atomic

   !> change the unit of B*Tau from T*ps to atomic unit
   BTauMax= BTauMax/Magneticfluxdensity_atomic*Relaxation_Time_Tau

   !> 
   allocate(Omega_array(OmegaNum))
   if (OmegaNum==1) then
      Omega_array(1)= OmegaMin
   else
      do i=1, OmegaNum
         Omega_array(i)= OmegaMin+ (OmegaMax-OmegaMin)* (i-1d0)/dble(OmegaNum-1)
      enddo ! i
   endif

!===============================================================================================================!
!> LATTICE card
!===============================================================================================================!

   NK = Nk1

   !> Read the cell information include the lattice and atom's position
   rewind(1001)
   lfound = .false.
   do while (.true.)
      read(1001, *, end= 100)inline
      inline= upper(inline)
      if (trim(adjustl(inline))=='LATTICE') then
         lfound= .true.
         if (cpuid==0) write(stdout, *)' '
         if (cpuid==0) write(stdout, *)'We found LATTICE card'
         exit
      endif
   enddo
100 continue

   if (lfound) then
      read(1001, *)inline   ! The unit of lattice vector
      inline= upper(inline)
      AngOrBohr=trim(adjustl(inline))
      read(1001, *)Origin_cell%Rua
      read(1001, *)Origin_cell%Rub
      read(1001, *)Origin_cell%Ruc
   else
      stop 'ERROR: please set lattice information'
   endif

   if (index(AngOrBohr, 'ANG')>0) then
      !> the input unit is Angstrom
      Origin_cell%Rua= Origin_cell%Rua*Angstrom2atomic
      Origin_cell%Rub= Origin_cell%Rub*Angstrom2atomic
      Origin_cell%Ruc= Origin_cell%Ruc*Angstrom2atomic
   endif
   Origin_cell%lattice(:, 1)= Origin_cell%Rua
   Origin_cell%lattice(:, 2)= Origin_cell%Rub
   Origin_cell%lattice(:, 3)= Origin_cell%Ruc

   !> cell parameters
   Origin_cell%cell_parameters(1)= norm(Origin_cell%Rua)
   Origin_cell%cell_parameters(2)= norm(Origin_cell%Rub)
   Origin_cell%cell_parameters(3)= norm(Origin_cell%Ruc)
   Origin_cell%cell_parameters(4)= angle(Origin_cell%Rub, Origin_cell%Ruc)
   Origin_cell%cell_parameters(5)= angle(Origin_cell%Rua, Origin_cell%Ruc)
   Origin_cell%cell_parameters(6)= angle(Origin_cell%Rua, Origin_cell%Rub)

   !> transform lattice from direct space to reciprocal space

   Origin_cell%Kua= 0d0
   Origin_cell%Kub= 0d0
   Origin_cell%Kuc= 0d0

   call get_volume(Origin_cell%Rua, Origin_cell%Rub, Origin_cell%Ruc, Origin_cell%CellVolume)
   Origin_cell%ReciprocalCellVolume= (2d0*pi)**3/Origin_cell%CellVolume

   call get_reciprocal_lattice(Origin_cell%Rua, Origin_cell%Rub, Origin_cell%Ruc, &
                               Origin_cell%Kua, Origin_cell%Kub, Origin_cell%Kuc)
   
   Origin_cell%reciprocal_lattice(:, 1)= Origin_cell%Kua
   Origin_cell%reciprocal_lattice(:, 2)= Origin_cell%Kub
   Origin_cell%reciprocal_lattice(:, 3)= Origin_cell%Kuc

   !> cell parameters
   Origin_cell%reciprocal_cell_parameters(1)= norm(Origin_cell%Kua)
   Origin_cell%reciprocal_cell_parameters(2)= norm(Origin_cell%Kub)
   Origin_cell%reciprocal_cell_parameters(3)= norm(Origin_cell%Kuc)
   Origin_cell%reciprocal_cell_parameters(4)= angle(Origin_cell%Kub, Origin_cell%Kuc)
   Origin_cell%reciprocal_cell_parameters(5)= angle(Origin_cell%Kua, Origin_cell%Kuc)
   Origin_cell%reciprocal_cell_parameters(6)= angle(Origin_cell%Kua, Origin_cell%Kub)

   if(cpuid==0)write(stdout, '(a)') '>> lattice information (Angstrom)'
   if(cpuid==0)write(stdout, '(6a12)')" a", " b", " c", 'alpha', 'beta', 'gamma'
   if(cpuid==0)write(stdout, '(6f12.6)')Origin_cell%cell_parameters/Angstrom2atomic
   if(cpuid==0)write(stdout, '(a)')" Three Lattice vectors of the unfolded cell: "
   if(cpuid==0)write(stdout, '(3f12.6)')Origin_cell%Rua/Angstrom2atomic
   if(cpuid==0)write(stdout, '(3f12.6)')Origin_cell%Rub/Angstrom2atomic
   if(cpuid==0)write(stdout, '(3f12.6)')Origin_cell%Ruc/Angstrom2atomic

   if(cpuid==0)write(stdout, '(a)') '>> Reciprocal lattice information (1/Angstrom)'
   if(cpuid==0)write(stdout, '(6a12)')" a", " b", " c", 'alpha', 'beta', 'gamma'
   if(cpuid==0)write(stdout, '(6f12.6)')Origin_cell%reciprocal_cell_parameters*Angstrom2atomic
   if(cpuid==0)write(stdout, '(a)')" Three reciprocal lattice vectors of the primitive cell: "
   if(cpuid==0)write(stdout, '(3f12.6)')Origin_cell%Kua*Angstrom2atomic
   if(cpuid==0)write(stdout, '(3f12.6)')Origin_cell%Kub*Angstrom2atomic
   if(cpuid==0)write(stdout, '(3f12.6)')Origin_cell%Kuc*Angstrom2atomic

!===============================================================================================================!
!> ATOM_POSITIONS card
!===============================================================================================================!

   !> Read atom positions information
   rewind(1001)
   lfound = .false.
   do while (.true.)
      read(1001, *, end= 101)inline
      inline=upper(inline)
      if (trim(adjustl(inline))=='ATOM_POSITIONS') then
         lfound= .true.
         if (cpuid==0) write(stdout, *)' '
         if (cpuid==0) write(stdout, *)'We found ATOM_POSITIONS card'
         exit
      endif
   enddo
101 continue

   if (lfound) then
      read(1001, *)Origin_cell%Num_atoms   ! total number of atoms
      if(cpuid==0)write(stdout, '(a, i5)')'Origin_cell%Num_atoms', Origin_cell%Num_atoms
      allocate(Origin_cell%atom_name(Origin_cell%Num_atoms))
      allocate(Origin_cell%Atom_position_cart(3, Origin_cell%Num_atoms))
      allocate(Origin_cell%Atom_position_direct(3, Origin_cell%Num_atoms))
      allocate(Atom_position_cart_newcell(3, Origin_cell%Num_atoms))
      allocate(Atom_position_direct_newcell(3, Origin_cell%Num_atoms))
      allocate(Origin_cell%Atom_magnetic_moment(3, Origin_cell%Num_atoms))
      Origin_cell%Atom_magnetic_moment= 0d0
      read(1001, *)inline   ! The unit of lattice vector
      DirectOrCart= trim(adjustl(inline))

      !> check whether we have the magnetic moment in the POSITION card
      do i=1, Origin_cell%Num_atoms
         read(1001, *, err=132) Origin_cell%atom_name(i), Origin_cell%Atom_position_cart(:, i), Origin_cell%Atom_magnetic_moment(:, i)
         if(cpuid==0)write(stdout, '(a4,3f12.6)')Origin_cell%atom_name(i), Origin_cell%Atom_position_cart(:, i)
         if (index(DirectOrCart, "D")>0)then
            pos= Origin_cell%Atom_position_cart(:, i)
            Origin_cell%Atom_position_cart(:, i)= pos(1)*Origin_cell%Rua+ pos(2)*Origin_cell%Rub+ pos(3)*Origin_cell%Ruc
         endif
      enddo
      go to 133

132   continue
      !> if the code comes to here, it means there is no atom's magnetic moment in the POSITION card
      if (cpuid==0) write(stdout, *) ' '
      if (cpuid==0) write(stdout, *) &
         "Warning: You didn't specify the atom magnetic moment in the ATOMIC_POSITION card", &
         " Or the format is wrong. ", &
         "So we set all the Atom-magnetic-moments to zero."
      Origin_cell%Atom_magnetic_moment= 0d0
      rewind(1001)
      do while (.true.)
         read(1001, *)inline
         inline=upper(inline)
         if (trim(adjustl(inline))=='ATOM_POSITIONS') then
            exit
         endif
      enddo
      !> skip two lines
      read(1001, *)
      read(1001, *)

      do i=1, Origin_cell%Num_atoms
         read(1001, *, err=134) Origin_cell%atom_name(i), Origin_cell%Atom_position_cart(:, i)
         !> Origin_cell%Atom_position_cart is in the cartesian coordinate.
         if (index(DirectOrCart, "D")>0)then
            pos= Origin_cell%Atom_position_cart(:, i)
            Origin_cell%Atom_position_cart(:, i)= pos(1)*Origin_cell%Rua+ pos(2)*Origin_cell%Rub+ pos(3)*Origin_cell%Ruc
         endif
      enddo
      go to 133
134   continue
      write(*, *)"ERROR happens in the ATOMIC_POSITION card"
      write(*, *)"This is a free format card, between lines there should be any comments"
      write(*, *)"The number in the second line should be the same as the number of lines of the atom positions."
      stop "ERROR: please set ATOMIC_POSITION card correctly, see manual on www.wanniertools.com"

133   continue

      do ia=1, Origin_cell%Num_atoms
         call cart_direct_real(Origin_cell%Atom_position_cart(:, ia), &
            Origin_cell%Atom_position_direct(:, ia), Origin_cell%lattice)
      enddo

      if(cpuid==0)write(stdout,'(a)')' '
      if(cpuid==0)write(stdout,'(a)')'>> Atom position and magnetic moment of Original lattice'
      if(cpuid==0)write(stdout,'(13X, 2a36, a24)')' Catesian', 'Fractional(Direct)', 'Magnetic moment'
      if(cpuid==0)write(stdout,'(a)')'------------------------------------------------------------------------------------------------------------------'
      if(cpuid==0)write(stdout,'(a6, 2X, a10, 6a12, 3a8)')'index', 'Atom Name ', ' x', ' y', ' z', 'a', 'b', 'c', 'Mx', 'My', 'Mz'
      if(cpuid==0)write(stdout,'(a)')'------------------------------------------------------------------------------------------------------------------'
      do i=1, Origin_cell%Num_atoms
         if(cpuid==0)write(stdout, '(i6,2X, a10,6f12.6,3f8.3)')i, Origin_cell%atom_name(i), &
            Origin_cell%Atom_position_cart(:, i), Origin_cell%Atom_position_direct(:,i), Origin_cell%Atom_magnetic_moment(:, i)
      enddo

   else
      stop "ERROR: please set atom's positions information correctly"
   endif


   !> setup atom type
   if (allocated(iarray_temp))deallocate(iarray_temp)
   allocate(iarray_temp(Origin_cell%Num_atoms))
   iarray_temp= 1
   do ia=1, Origin_cell%Num_atoms
      char_temp= Origin_cell%atom_name(ia)
      do i=ia+1, Origin_cell%Num_atoms
         if (char_temp==Origin_cell%atom_name(i).and.iarray_temp(i)/=0)then
            iarray_temp(i)=0
         endif
      enddo
   enddo
   Origin_cell%Num_atom_type= sum(iarray_temp)

   allocate(Origin_cell%Name_of_atomtype(Origin_cell%Num_atom_type))
   allocate(Origin_cell%Num_atoms_eachtype(Origin_cell%Num_atom_type))
   allocate(Origin_cell%itype_atom(Origin_cell%Num_atoms))
   it = 0
   do ia=1, Origin_cell%Num_atoms
      if (iarray_temp(ia)/=0) then
         it= it+ 1
         Origin_cell%Name_of_atomtype(it)= Origin_cell%atom_name(ia)
      endif
   enddo

   !> find the type of atoms and label them
   do ia=1, Origin_cell%Num_atoms
      do i=1, Origin_cell%Num_atom_type
         if (Origin_cell%atom_name(ia)==Origin_cell%Name_of_atomtype(i))then
            Origin_cell%itype_atom(ia)= i
         endif
      enddo
   enddo

   do i=1, Origin_cell%Num_atom_type
      it = 0
      do ia=1, Origin_cell%Num_atoms
         if (Origin_cell%atom_name(ia)==Origin_cell%Name_of_atomtype(i))then
            it = it+ 1
         endif
      enddo
      Origin_cell%Num_atoms_eachtype(i)= it
   enddo


!===============================================================================================================!
!> PROJECTORS card
!===============================================================================================================!

   !> Read projectors information
   rewind(1001)
   lfound = .false.
   do while (.true.)
      read(1001, *, end= 102)inline
      inline=upper(inline)
      if (trim(adjustl(inline))=='PROJECTORS'.or.&
         trim(adjustl(inline))=='PROJECTOR') then
         lfound= .true.
         if (cpuid==0) write(stdout, *)' '
         if (cpuid==0) write(stdout, *)'We found PROJECTORS card'
         exit
      endif
   enddo
102 continue

   if (lfound) then
      allocate(Origin_cell%nprojs(Origin_cell%Num_atoms))
      Origin_cell%nprojs= 0
      read(1001, *)Origin_cell%nprojs
      if(cpuid==0)write(stdout, '(a)')' >>Number of projectors per atoms:'
      if(cpuid==0)write(stdout, '(10i6)')Origin_cell%nprojs

      Origin_cell%max_projs= maxval(Origin_cell%nprojs)
      allocate(Origin_cell%proj_name(Origin_cell%max_projs, Origin_cell%Num_atoms))
      Origin_cell%proj_name= ' '
      do i=1, Origin_cell%Num_atoms
         read(1001, *)char_temp, Origin_cell%proj_name(1:Origin_cell%nprojs(i), i)
         if(cpuid==0)write(stdout, '(40a8)') &
            char_temp, Origin_cell%proj_name(1:Origin_cell%nprojs(i), i)
      enddo
   else
      stop "ERROR: please set projectors for Wannier functions information"
   endif

   !> set up orbitals_start
   allocate(orbitals_start(Origin_cell%Num_atoms))
   orbitals_start= 1
   do i=1, Origin_cell%Num_atoms-1
      orbitals_start(i+1)= orbitals_start(i)+ Origin_cell%nprojs(i)
   enddo

   !> orbital index order
   allocate(index_start(Origin_cell%Num_atoms))
   allocate(index_end  (Origin_cell%Num_atoms))
   index_start= 0
   index_end= 0
   index_start(1)= 1
   index_end(1)= Origin_cell%nprojs(1)
   do i=2, Origin_cell%Num_atoms
      index_start(i)= index_start(i-1)+ Origin_cell%nprojs(i-1)
      index_end(i)= index_end(i-1)+ Origin_cell%nprojs(i)
   enddo



   !> read Wannier centres
   NumberOfspinorbitals= sum(Origin_cell%nprojs)
   if (SOC>0.or.Add_Zeeman_Field) NumberOfspinorbitals= 2*NumberOfspinorbitals
   Origin_cell%NumberOfspinorbitals= NumberOfspinorbitals
   allocate(Origin_cell%spinorbital_to_atom_index(NumberOfspinorbitals))
   allocate(Origin_cell%spinorbital_to_projector_index(NumberOfspinorbitals))
   allocate(Origin_cell%wannier_centers_cart(3, NumberOfspinorbitals))
   allocate(Origin_cell%wannier_centers_direct(3, NumberOfspinorbitals))
   Origin_cell%wannier_centers_direct= 0d0
   Origin_cell%wannier_centers_cart= 0d0
   !> default wannier centers
   i= 0
   do ia= 1, Origin_cell%Num_atoms
      do j= 1, Origin_cell%nprojs(ia)
         i= i+ 1
         Origin_cell%spinorbital_to_atom_index(i)= ia
         Origin_cell%spinorbital_to_projector_index(i)= j
         Origin_cell%wannier_centers_cart(:, i)=  Origin_cell%Atom_position_cart(:, ia)
         call cart_direct_real(Origin_cell%wannier_centers_cart(:, i),  &
            Origin_cell%wannier_centers_direct(:, i), &
            Origin_cell%lattice)
         if (SOC>0.or.Add_Zeeman_Field) then
            Origin_cell%spinorbital_to_atom_index(i+NumberOfspinorbitals/2)= ia
            Origin_cell%spinorbital_to_projector_index(i+NumberOfspinorbitals/2)= j
            Origin_cell%wannier_centers_cart(:, i+NumberOfspinorbitals/2)=  Origin_cell%Atom_position_cart(:, ia)
            call cart_direct_real(Origin_cell%wannier_centers_cart(:, i+NumberOfspinorbitals/2),  &
               Origin_cell%wannier_centers_direct(:, i+NumberOfspinorbitals/2), &
               Origin_cell%lattice)
         endif
      enddo ! j
   enddo ! ia

   rewind(1001)
   lfound = .false.
   do while (.true.)
      read(1001, *, end= 110)inline
      inline=upper(inline)
      if (trim(adjustl(inline))=='Origin_cell%wannier_ceters' &
         .or. trim(adjustl(inline))=='WANNIER_CENTRES') then
         lfound= .true.
         if (cpuid==0) write(stdout, *)' '
         if (cpuid==0) write(stdout, *)'We found Origin_cell%wannier_ceters card'
         exit
      endif
   enddo

   if (lfound) then
      read(1001, *)inline   ! Direct or Cartesian
      inline= upper(inline)
      DirectOrCart= trim(adjustl(inline))

      it= 0
      if (index(DirectOrCart, "D")>0)then
         if (SOC==0.and.Add_Zeeman_Field) then
            do i=1, NumberOfspinorbitals/2
               read(1001, *, end=207, err=207) Origin_cell%wannier_centers_direct(:, i)
               it= it+ 2
               call direct_cart_real(Origin_cell%wannier_centers_direct(:, i), &
                  Origin_cell%wannier_centers_cart(:, i), Origin_cell%lattice)
               Origin_cell%wannier_centers_cart(:, i+NumberOfspinorbitals/2)= &
               Origin_cell%wannier_centers_cart(:, i)
               Origin_cell%wannier_centers_direct(:, i+NumberOfspinorbitals/2)= &
               Origin_cell%wannier_centers_direct(:, i)
            enddo
         else
            do i=1, NumberOfspinorbitals
               read(1001, *, end=207, err=207) Origin_cell%wannier_centers_direct(:, i)
               it= it+ 1
               call direct_cart_real(Origin_cell%wannier_centers_direct(:, i), &
                  Origin_cell%wannier_centers_cart(:, i), Origin_cell%lattice)
            enddo
         endif

      else
         if (SOC==0.and.Add_Zeeman_Field) then
            do i=1, NumberOfspinorbitals/2
               read(1001, *, end=207, err=207) Origin_cell%wannier_centers_cart(:, i)
               it= it+ 2
               call cart_direct_real(Origin_cell%wannier_centers_cart(:, i), &
                  Origin_cell%wannier_centers_direct(:, i), Origin_cell%lattice)
               Origin_cell%wannier_centers_cart(:, i+NumberOfspinorbitals/2)= &
               Origin_cell%wannier_centers_cart(:, i)
               Origin_cell%wannier_centers_direct(:, i+NumberOfspinorbitals/2)= &
               Origin_cell%wannier_centers_direct(:, i)
            enddo
         else
            do i=1, NumberOfspinorbitals
               read(1001, *, end=207, err=207) Origin_cell%wannier_centers_cart(:, i)
               it= it+ 1
               call cart_direct_real(Origin_cell%wannier_centers_cart(:, i), &
                  Origin_cell%wannier_centers_direct(:, i), Origin_cell%lattice)
            enddo
         endif
      endif
   endif ! found Origin_cell%wannier_ceters card
207 continue
   if (it< NumberOfspinorbitals.and.cpuid==0) then
      write(stdout, *)' '
      write(stdout, *)' >>>> Error happens in Wannier_centres card'
      write(stdout, *)' Error: the number of Origin_cell%wannier_ceters lines should '
      write(stdout, *)' equal to the number wannier functions (include spin)'
      write(stdout, *)' Num_wann', NumberOfspinorbitals, ' the centres lines you given ', it
      write(stdout, *)' Otherwise, if you do not know the meaning of this,'
      write(stdout, *)' please delete this card'
      stop
   endif



110 continue

   if (lfound) then
      if (cpuid==0) then
         write(stdout, *)" "
         write(stdout, *)">> Wannier centers from wt.in, in unit of unit lattice vectors"
         write(stdout, '(a6, 6a10)')'iwann', 'R1', 'R2', 'R3', 'x', 'y', 'z'
         do i=1, NumberOfspinorbitals
            write(stdout, '(i6, 6f10.6)')i, Origin_cell%wannier_centers_direct(:, i), Origin_cell%wannier_centers_cart(:, i)
         enddo
      endif
   else
      if (cpuid==0) then
         write(stdout, *)" "
         write(stdout, *)">> Wannier centers by default, in unit of unit lattice vectors"
         write(stdout, '(a6, 6a10)')'iwann', 'R1', 'R2', 'R3', 'x', 'y', 'z'
         do i=1, NumberOfspinorbitals
            write(stdout, '(i6, 6f10.6)')i, Origin_cell%wannier_centers_direct(:, i), Origin_cell%wannier_centers_cart(:, i)
         enddo
      endif
   endif

  
!===============================================================================================================!
!> LATTICE_UNFOLD card
!===============================================================================================================!

   !>> This segment is for folded lattice which is smaller than the original lattice usually
   !> read folded lattice information
   rewind(1001)
   lfound = .false.
   do while (.true.)
      read(1001, *, end= 1080)inline
      inline= upper(inline)
      if (trim(adjustl(inline))=='LATTICE_UNFOLD') then
         lfound= .true.
         if (cpuid==0) write(stdout, *)' '
         if (cpuid==0) write(stdout, *)'We found LATTICE_UNFOLD card'
         exit
      endif
   enddo
1080 continue

   if (lfound) then
      read(1001, *)inline   ! The unit of lattice vector
      inline= upper(inline)
      AngOrBohr=trim(adjustl(inline))
      read(1001, *)Folded_cell%Rua
      read(1001, *)Folded_cell%Rub
      read(1001, *)Folded_cell%Ruc
      if (index(AngOrBohr, 'ANG')>0) then
         Folded_cell%Rua= Folded_cell%Rua*Angstrom2atomic
         Folded_cell%Rub= Folded_cell%Rub*Angstrom2atomic
         Folded_cell%Ruc= Folded_cell%Ruc*Angstrom2atomic
      endif

   else
      if (cpuid==0) write(stdout, *)"We didn't found LATTICE_UNFOLD card, so it is the same as the unit cell."
      Folded_cell%Rua=Origin_cell%Rua
      Folded_cell%Rub=Origin_cell%Rub
      Folded_cell%Ruc=Origin_cell%Ruc
   endif

   !> cell parameters
   Folded_cell%cell_parameters(1)= norm(Folded_cell%Rua)
   Folded_cell%cell_parameters(2)= norm(Folded_cell%Rub)
   Folded_cell%cell_parameters(3)= norm(Folded_cell%Ruc)
   Folded_cell%cell_parameters(4)= angle(Folded_cell%Rub, Folded_cell%Ruc)
   Folded_cell%cell_parameters(5)= angle(Folded_cell%Rua, Folded_cell%Ruc)
   Folded_cell%cell_parameters(6)= angle(Folded_cell%Rua, Folded_cell%Rub)


   !> transform lattice from direct space to reciprocal space

   Folded_cell%Kua= 0d0
   Folded_cell%Kub= 0d0
   Folded_cell%Kuc= 0d0
   call get_volume(Folded_cell%Rua, Folded_cell%Rub, Folded_cell%Ruc, Folded_cell%CellVolume )
   Folded_cell%ReciprocalCellVolume= (2d0*3.1415926535d0)**3/Folded_cell%CellVolume

   call get_reciprocal_lattice(Folded_cell%Rua, Folded_cell%Rub, Folded_cell%Ruc, & 
                               Folded_cell%Kua, Folded_cell%Kub, Folded_cell%Kuc)
   
   !> reciprlcal cell parameters
   Folded_cell%reciprocal_cell_parameters(1)= norm(Folded_cell%Kua)
   Folded_cell%reciprocal_cell_parameters(2)= norm(Folded_cell%Kub)
   Folded_cell%reciprocal_cell_parameters(3)= norm(Folded_cell%Kuc)
   Folded_cell%reciprocal_cell_parameters(4)= angle(Folded_cell%Kub, Folded_cell%Kuc)
   Folded_cell%reciprocal_cell_parameters(5)= angle(Folded_cell%Kua, Folded_cell%Kuc)
   Folded_cell%reciprocal_cell_parameters(6)= angle(Folded_cell%Kua, Folded_cell%Kub)

   if(cpuid==0)write(stdout, '(a)') '>> Folded lattice information (Angstrom)'
   if(cpuid==0)write(stdout, '(6a12)')" a", " b", " c", 'alpha', 'beta', 'gamma'
   if(cpuid==0)write(stdout, '(6f12.6)')Folded_cell%cell_parameters/Angstrom2atomic
   if(cpuid==0)write(stdout, '(a)')" Three Lattice vectors of unfolded cell: "
   if(cpuid==0)write(stdout, '(3f12.6)')Folded_cell%Rua/Angstrom2atomic
   if(cpuid==0)write(stdout, '(3f12.6)')Folded_cell%Rub/Angstrom2atomic
   if(cpuid==0)write(stdout, '(3f12.6)')Folded_cell%Ruc/Angstrom2atomic

   if(cpuid==0)write(stdout, '(a)') '>> Folded Reciprocal lattice information (1/Angstrom)'
   if(cpuid==0)write(stdout, '(6a12)')" a", " b", " c", 'alpha', 'beta', 'gamma'
   if(cpuid==0)write(stdout, '(6f12.6)')Folded_cell%reciprocal_cell_parameters*Angstrom2atomic
   if(cpuid==0)write(stdout, '(a)')" Three reciprocal lattice vectors of unfolded cell: "
   if(cpuid==0)write(stdout, '(3f12.6)')Folded_cell%Kua*Angstrom2atomic
   if(cpuid==0)write(stdout, '(3f12.6)')Folded_cell%Kub*Angstrom2atomic
   if(cpuid==0)write(stdout, '(3f12.6)')Folded_cell%Kuc*Angstrom2atomic

   !> calculate the coordinates of Origin_cell in the unit of Folded_cell
   call cart_direct_real_unfold(Origin_cell%Rua, R1)
   call cart_direct_real_unfold(Origin_cell%Rub, R2)
   call cart_direct_real_unfold(Origin_cell%Ruc, R3)

   if(cpuid==0) then
      write(stdout, '(a)') ' '
      write(stdout, '(a)') '>> The relation between the original lattice and the unfolded lattice'
      write(stdout, '(a, f12.6)') '>> The cell volume ratio is: ', Origin_cell%CellVolume/Folded_cell%CellVolume
      write(stdout, '(a, f12.6)') '>> The lattice constant a ratio is: ', Origin_cell%cell_parameters(1)/Folded_cell%cell_parameters(1)
      write(stdout, '(a, f12.6)') '>> The lattice constant b ratio is: ', Origin_cell%cell_parameters(2)/Folded_cell%cell_parameters(2)
      write(stdout, '(a, f12.6)') '>> The lattice constant c ratio is: ', Origin_cell%cell_parameters(3)/Folded_cell%cell_parameters(3)
      write(stdout, '(a, 3f12.6)') '>> Origin_cell in unit of Folded_cell: R1  ', R1 
      write(stdout, '(a, 3f12.6)') '>> Origin_cell in unit of Folded_cell: R2  ', R2 
      write(stdout, '(a, 3f12.6)') '>> Origin_cell in unit of Folded_cell: R3  ', R3 
      write(stdout, '(a)') ' '
   endif


!===============================================================================================================!
!> ATOM_POSITIONS_UNFOLD card
!===============================================================================================================!

   !> Read atom positions information
   rewind(1001)
   lfound = .false.
   do while (.true.)
      read(1001, *, end= 1013)inline
      inline=upper(inline)
      if (trim(adjustl(inline))=='ATOM_POSITIONS_UNFOLD') then
         lfound= .true.
         if (cpuid==0) write(stdout, *)' '
         if (cpuid==0) write(stdout, *)'We found ATOM_POSITIONS_UNFOLD card'
         exit
      endif
   enddo
1013 continue

   if (lfound) then
      read(1001, *)Folded_cell%Num_atoms   ! total number of atoms
      if(cpuid==0)write(stdout, '(a, i5)')'Folded_cell%Num_atoms', Folded_cell%Num_atoms
      allocate(Folded_cell%atom_name(Folded_cell%Num_atoms))
      allocate(Folded_cell%Atom_position_cart(3, Folded_cell%Num_atoms))
      allocate(Folded_cell%Atom_position_direct(3, Folded_cell%Num_atoms))
      allocate(Folded_cell%Atom_magnetic_moment(3, Folded_cell%Num_atoms))
      Folded_cell%Atom_magnetic_moment= 0d0
      read(1001, *)inline   ! The unit of lattice vector
      DirectOrCart= trim(adjustl(inline))

      !> check whether we have the magnetic moment in the POSITION card
      do i=1, Folded_cell%Num_atoms
         read(1001, *, err=1320) Folded_cell%atom_name(i), Folded_cell%Atom_position_cart(:, i), Folded_cell%Atom_magnetic_moment(:, i)
         if(cpuid==0)write(stdout, '(a4,3f12.6)')Folded_cell%atom_name(i), Folded_cell%Atom_position_cart(:, i)
         if (index(DirectOrCart, "D")>0)then
            pos= Folded_cell%Atom_position_cart(:, i)
            Folded_cell%Atom_position_cart(:, i)= pos(1)*Folded_cell%Rua+ pos(2)*Folded_cell%Rub+ pos(3)*Folded_cell%Ruc
         endif
      enddo
      go to 1330

1320  continue
      !> if the code comes to here, it means there is no atom's magnetic moment in the POSITION card
      if (cpuid==0) write(stdout, *) ' '
      if (cpuid==0) write(stdout, *) &
         "Warning: You didn't specify the atom magnetic moment in the ATOMIC_POSITION card", &
         " Or the format is wrong. ", &
         "So we set all the Atom-magnetic-moments to zero."
      Folded_cell%Atom_magnetic_moment= 0d0
      rewind(1001)
      do while (.true.)
         read(1001, *)inline
         inline=upper(inline)
         if (trim(adjustl(inline))=='ATOM_POSITIONS_UNFOLD') then
            exit
         endif
      enddo
      !> skip two lines
      read(1001, *)
      read(1001, *)

      do i=1, Folded_cell%Num_atoms
         read(1001, *, err=1340) Folded_cell%atom_name(i), Folded_cell%Atom_position_cart(:, i)
         !> Folded_cell%Atom_position_cart is in the cartesian coordinate.
         if (index(DirectOrCart, "D")>0)then
            pos= Folded_cell%Atom_position_cart(:, i)
            Folded_cell%Atom_position_cart(:, i)= pos(1)*Folded_cell%Rua+ pos(2)*Folded_cell%Rub+ pos(3)*Folded_cell%Ruc
         endif
      enddo
      go to 1330
1340  continue
      write(*, *)"ERROR happens in the ATOM_POSITION_UNFOLD card"
      write(*, *)"This is a free format card, between lines there should be any comments"
      write(*, *)"The number in the second line should be the same as the number of lines of the atom positions."
      stop "ERROR: please set ATOM_POSITION_UNFOLD card correctly, see manual on www.wanniertools.com"

1330  continue

      if(cpuid==0)write(stdout,'(a)')' '
      do ia=1, Folded_cell%Num_atoms
         call cart_direct_real_unfold(Folded_cell%Atom_position_cart(:, ia), Folded_cell%Atom_position_direct(:, ia))
         if(cpuid==0)write(stdout, '(a4,3f12.6)')Folded_cell%atom_name(ia), Folded_cell%Atom_position_direct(:, ia)
      enddo
   else
      if (abs(Folded_cell%CellVolume-Origin_cell%CellVolume)>eps6) then
         call printerrormsg("ERROR: please set ATOM_POSITIONS_UNFOLD since you set LATTICE_UNFOLD")
      else
         Folded_cell%Num_atoms= Origin_cell%Num_atoms
         allocate(Folded_cell%atom_name(Folded_cell%Num_atoms))
         allocate(Folded_cell%Atom_position_cart(3, Folded_cell%Num_atoms))
         allocate(Folded_cell%Atom_position_direct(3, Folded_cell%Num_atoms))
         allocate(Folded_cell%Atom_magnetic_moment(3, Folded_cell%Num_atoms))
         Folded_cell%atom_name= Origin_cell%atom_name
         Folded_cell%Atom_position_cart= Origin_cell%Atom_position_cart
         Folded_cell%Atom_position_direct= Origin_cell%Atom_position_direct
         Folded_cell%Atom_magnetic_moment= Origin_cell%Atom_magnetic_moment
      endif
   endif

   if(cpuid==0)write(stdout,'(a)')' '
   if(cpuid==0)write(stdout,'(a)')'>>> Atom position and magnetic moment of Unfolded lattice'
   if(cpuid==0)write(stdout,'(13X, 2a36, a24)')' Catesian(Ang)', 'Fractional(Direct)', 'Magnetic moment'
   if(cpuid==0)write(stdout,'(a)')'------------------------------------------------------------------------------------------------------------------'
   if(cpuid==0)write(stdout,'(a6, 2X, a10, 6a12, 3a8)')'index', 'Atom Name ', ' x', ' y', ' z', 'a', 'b', 'c', 'Mx', 'My', 'Mz'
   if(cpuid==0)write(stdout,'(a)')'------------------------------------------------------------------------------------------------------------------'
   do i=1, Folded_cell%Num_atoms
      if(cpuid==0)write(stdout, '(i6,2X, a10,6f12.6,3f8.3)')i, Folded_cell%atom_name(i), &
         Folded_cell%Atom_position_cart(:, i)/Angstrom2atomic, Folded_cell%Atom_position_direct(:,i), Folded_cell%Atom_magnetic_moment(:, i)
   enddo


   !> setup atom type
   if (allocated(iarray_temp))deallocate(iarray_temp)
   allocate(iarray_temp(Folded_cell%Num_atoms))
   iarray_temp= 1
   do ia=1, Folded_cell%Num_atoms
      char_temp= Folded_cell%atom_name(ia)
      do i=ia+1, Folded_cell%Num_atoms
         if (char_temp==Folded_cell%atom_name(i).and.iarray_temp(i)/=0)then
            iarray_temp(i)=0
         endif
      enddo
   enddo
   Folded_cell%Num_atom_type= sum(iarray_temp)

   allocate(Folded_cell%Name_of_atomtype(Folded_cell%Num_atom_type))
   allocate(Folded_cell%Num_atoms_eachtype(Folded_cell%Num_atom_type))
   allocate(Folded_cell%itype_atom(Folded_cell%Num_atoms))
   it = 0
   do ia=1, Folded_cell%Num_atoms
      if (iarray_temp(ia)/=0) then
         it= it+ 1
         Folded_cell%Name_of_atomtype(it)= Folded_cell%atom_name(ia)
      endif
   enddo

   !> find the type of atoms and label them
   do ia=1, Folded_cell%Num_atoms
      do i=1, Folded_cell%Num_atom_type
         if (Folded_cell%atom_name(ia)==Folded_cell%Name_of_atomtype(i))then
            Folded_cell%itype_atom(ia)= i
         endif
      enddo
   enddo

   do i=1, Folded_cell%Num_atom_type
      it = 0
      do ia=1, Folded_cell%Num_atoms
         if (Folded_cell%atom_name(ia)==Folded_cell%Name_of_atomtype(i))then
            it = it+ 1
         endif
      enddo
      Folded_cell%Num_atoms_eachtype(i)= it
   enddo

   call writeout_poscar(Folded_cell, "POSCAR-Folded")
!===============================================================================================================!
!> PROJECTORS_UNFOLD card
!===============================================================================================================!

   !> Read projectors information for unfolded lattice
   rewind(1001)
   lfound = .false.
   do while (.true.)
      read(1001, *, end= 1020)inline
      inline=upper(inline)
      if (trim(adjustl(inline))=='PROJECTORS_UNFOLD'.or.&
         trim(adjustl(inline))=='PROJECTOR_UNFOLD') then
         lfound= .true.
         if (cpuid==0) write(stdout, *)' '
         if (cpuid==0) write(stdout, *)'We found PROJECTORS_UNFOLD card'
         exit
      endif
   enddo
1020 continue

   if (lfound) then
      allocate(Folded_cell%nprojs(Folded_cell%Num_atoms))
      Folded_cell%nprojs= 0
      read(1001, *)Folded_cell%nprojs
      if(cpuid==0)write(stdout, '(a)')' >>Number of projectors per atoms:'
      if(cpuid==0)write(stdout, '(10i6)')Folded_cell%nprojs

      Folded_cell%max_projs= maxval(Folded_cell%nprojs)
      allocate(Folded_cell%proj_name(Folded_cell%max_projs, Folded_cell%Num_atoms))
      Folded_cell%proj_name= ' '
      do i=1, Folded_cell%Num_atoms
         read(1001, *)char_temp, Folded_cell%proj_name(1:Folded_cell%nprojs(i), i)
         if(cpuid==0)write(stdout, '(40a8)') &
            char_temp, Folded_cell%proj_name(1:Folded_cell%nprojs(i), i)
      enddo
   else
      if (abs(Folded_cell%CellVolume-Origin_cell%CellVolume)>eps6) then
         call printerrormsg("ERROR: please set PROJECTORS_UNFOLD since you set LATTICE_UNFOLD")
      else
         allocate(Folded_cell%nprojs(Folded_cell%Num_atoms))
         Folded_cell%nprojs= Origin_cell%nprojs
         Folded_cell%max_projs= maxval(Folded_cell%nprojs)
         allocate(Folded_cell%proj_name(Folded_cell%max_projs, Folded_cell%Num_atoms))
         Folded_cell%proj_name= Origin_cell%proj_name
      endif
   endif

   !> Wannier centres for unfoled lattice
   NumberOfspinorbitals_unfold= sum(Folded_cell%nprojs)
   if (SOC>0) NumberOfspinorbitals_unfold= 2*NumberOfspinorbitals_unfold
   Folded_cell%NumberOfspinorbitals= NumberOfspinorbitals_unfold
   allocate(Folded_cell%spinorbital_to_atom_index(NumberOfspinorbitals_unfold))
   allocate(Folded_cell%spinorbital_to_projector_index(NumberOfspinorbitals_unfold))
   allocate(Folded_cell%wannier_centers_cart(3, NumberOfspinorbitals_unfold))
   allocate(Folded_cell%wannier_centers_direct(3, NumberOfspinorbitals_unfold))
   Folded_cell%wannier_centers_direct= 0d0
   Folded_cell%wannier_centers_cart= 0d0
   !> default wannier centers
   i= 0
   do ia= 1, Folded_cell%Num_atoms
      do j= 1, Folded_cell%nprojs(ia)
         i= i+ 1
         Folded_cell%spinorbital_to_atom_index(i)= ia
         Folded_cell%spinorbital_to_projector_index(i)= j
         Folded_cell%wannier_centers_cart(:, i)=  Folded_cell%Atom_position_cart(:, ia)
         call cart_direct_real_unfold(Folded_cell%wannier_centers_cart(:, i),  &
            Folded_cell%wannier_centers_direct(:, i))
         if (SOC>0) then
            Folded_cell%spinorbital_to_atom_index(i+NumberOfspinorbitals_unfold/2)= ia
            Folded_cell%spinorbital_to_projector_index(i+NumberOfspinorbitals_unfold/2)= j
            Folded_cell%wannier_centers_cart(:, i+NumberOfspinorbitals_unfold/2)=  Folded_cell%Atom_position_cart(:, ia)
            call cart_direct_real_unfold(Folded_cell%wannier_centers_cart(:, i+NumberOfspinorbitals_unfold/2),  &
               Folded_cell%wannier_centers_direct(:, i+NumberOfspinorbitals_unfold/2))
         endif
      enddo ! j
   enddo ! ia

!===============================================================================================================!
!> MILLER_INDICES card
!===============================================================================================================!

   !> read surface information by Miller indices
   rewind(1001)
   lfound = .false.
   do while (.true.)
      read(1001, *, end= 224)inline
      inline=upper(inline)
      if (trim(adjustl(inline))=='MILLER_INDICES'  &
         .or. trim(adjustl(inline))=='MILLER_INDEX') then
         lfound= .true.
         if (cpuid==0) write(stdout, *)' '
         if (cpuid==0) write(stdout, *)'We found MILLER_INDEX card for slab calculations'
         exit
      endif
   enddo
224 continue

   MillerIndices= 0
   if (lfound) then
      read(1001, *, err=225, iostat=stat) MillerIndices(:)
225   continue
      if (stat/=0) stop "Something wrong with setting of MillerIndices, they should be like this 1 0 0"
      if (cpuid.eq.0) then
         write(stdout, '(a, 3i6)')'  Miller indices are :', MillerIndices
      endif
      call MillerIndicestoumatrix()
   endif

!===============================================================================================================!
!> SURFACE card
!===============================================================================================================!

   !> read surface information
   rewind(1001)
   lfound = .false.
   do while (.true.)
      read(1001, *, end= 103)inline
      inline=upper(inline)
      if (trim(adjustl(inline))=='SURFACE') then
         lfound= .true.
         if (cpuid==0) write(stdout, *)' '
         if (cpuid==0) write(stdout, *)'We found SURFACE card'
         exit
      endif
   enddo
103 continue

   if (.not.lfound.and.sum(abs(MillerIndices))==0) then
      Umatrix(1, :)=(/1d0, 0d0, 0d0/)
      Umatrix(2, :)=(/0d0, 1d0, 0d0/)
      Umatrix(3, :)=(/0d0, 0d0, 1d0/)
      if (cpuid==0) then
         write(stdout, *) "Warnning: You didn't set SURFACE card, by default, it's (001) surface."
         write(stdout, '(a, 3f12.6)')' The 1st vector on surface     :', Umatrix(1, :)
         write(stdout, '(a, 3f12.6)')' The 2nd vector on surface     :', Umatrix(2, :)
         write(stdout, '(a, 3f12.6)')' The 3rd vector out of surface :', Umatrix(3, :)
      endif
   endif

   if (lfound) then
      !> read information for new lattice
      !> in order to get different surface state
      !> R1'=U11*R1+U12*R2+U13*R3
      !> R2'=U21*R1+U22*R2+U23*R3
      !> R3'=U31*R1+U32*R2+U33*R3
      read(1001, *)Umatrix(1, :)
      read(1001, *)Umatrix(2, :)

      Umatrix(3,:)=(/0.0,0.0,1.0/)
      read(1001, *, err=260, iostat=stat)Umatrix(3, :)
260   continue

      if (cpuid==0) then
         write(stdout, '(a)')' '
         write(stdout, '(a)')'>> new vectors to define the surface (in unit of lattice vector)'
         write(stdout, '(a, 3f12.6)')' The 1st vector on surface     :', Umatrix(1, :)
         write(stdout, '(a, 3f12.6)')' The 2nd vector on surface     :', Umatrix(2, :)
         write(stdout, '(a, 3f12.6)')' The 3rd vector out of surface :', Umatrix(3, :)
      endif

      if (sum(abs(MillerIndices))>0 .and. cpuid==0) then
         write(stdout, '(a, a, a)')' Attention: you specified surface information twice.' , &
            ' However, we only take the information from SURFACE card, ', &
            ' and omit the settings by MILLER_INDEX card'
      endif
   endif

   !> check whether Umatrix is right
   !> the volume of the new cell should be the same as the old ones
   !> Here R1, R2, R3 are vectors defined by SURFACE CARD in original cartesian coordinates
   R1= Umatrix(1, 1)*Origin_cell%Rua+ Umatrix(1, 2)*Origin_cell%Rub+ Umatrix(1, 3)*Origin_cell%Ruc
   R2= Umatrix(2, 1)*Origin_cell%Rua+ Umatrix(2, 2)*Origin_cell%Rub+ Umatrix(2, 3)*Origin_cell%Ruc
   R3= Umatrix(3, 1)*Origin_cell%Rua+ Umatrix(3, 2)*Origin_cell%Rub+ Umatrix(3, 3)*Origin_cell%Ruc

   cell_volume2= R1(1)*(R2(2)*R3(3)- R2(3)*R3(2)) &
      +R1(2)*(R2(3)*R3(1)- R2(1)*R3(3)) &
      +R1(3)*(R2(1)*R3(2)- R2(2)*R3(1))

   if (cell_volume2<0) then
      R3=-R3
      Umatrix(3, :)= -Umatrix(3, :)
   endif

   if (abs(abs(cell_volume2)-abs(Origin_cell%CellVolume))> 0.001d0.and.cpuid==0) then
      write(stdout, '(a)')' '
      write(stdout, '(2a)')' Warnning: The Umatrix is wrongly set, the new cell', &
         'volume should be the same as the old ones. '
      write(stdout, '(a,2f10.4)')' cell_volume vs cell_volume-new', Origin_cell%CellVolume, cell_volume2
      write(stdout, '(a)')" However, don't worry, WannierTools will help you to find a suitable rotation matrix."
      write(stdout, '(a)')" I am looking for new unit cell atuomatically: "
   endif
   if (abs(abs(cell_volume2)-abs(Origin_cell%CellVolume))> 0.001d0) then
      call FindTheThirdLatticeVector()
      R1= Umatrix(1, 1)*Origin_cell%Rua+ Umatrix(1, 2)*Origin_cell%Rub+ Umatrix(1, 3)*Origin_cell%Ruc
      R2= Umatrix(2, 1)*Origin_cell%Rua+ Umatrix(2, 2)*Origin_cell%Rub+ Umatrix(2, 3)*Origin_cell%Ruc
      R3= Umatrix(3, 1)*Origin_cell%Rua+ Umatrix(3, 2)*Origin_cell%Rub+ Umatrix(3, 3)*Origin_cell%Ruc
      if (cpuid==0) then
         write(stdout, '(a)')' '
         write(stdout, '(a)')'>> New SURFACE CARD:'
         write(stdout, '(a, 3f12.6)')' The 1st vector on surface     :', Umatrix(1, :)
         write(stdout, '(a, 3f12.6)')' The 2nd vector on surface     :', Umatrix(2, :)
         write(stdout, '(a, 3f12.6)')' The 3rd vector out of surface :', Umatrix(3, :)
      endif
   endif

   !> print out the new basis
   if (cpuid.eq.0) then
      write(stdout, *)" "
      write(stdout, *)"The rotated new unit cell in cartesian coordinates : "
      write(stdout, '(3f12.6)') R1
      write(stdout, '(3f12.6)') R2
      write(stdout, '(3f12.6)') R3

      call get_volume(R1, R2, R3, cell_volume2)
      write(stdout, '(a, f18.5, a)')"New cell's Volume is ", cell_volume2/(Angstrom2atomic**3), 'Ang^3'
      write(stdout, *)" "
   endif

   if (cpuid.eq.0) then
      write(stdout, *)"Fractional coordinates of atoms in units of new lattice vectors : "
      do ia=1, Origin_cell%Num_atoms
         call rotate_newlattice(Origin_cell%Atom_position_direct(:, ia), Rt)
         call transformtohomecell(Rt)
         if(cpuid==0)write(stdout, '(a4,3f12.6)')Origin_cell%atom_name(ia), Rt
      enddo
      write(stdout, *)" "
   endif

   !> print out the new basis
   outfileindex= outfileindex+ 1
   if (cpuid.eq.0) then
      open(outfileindex, file="POSCAR-rotated")
      write(outfileindex, '(a)')"Rotated POSCAR by SURFACE card in wt.in by WannierTools"
      write(outfileindex, '(a)')"1.0"
      write(outfileindex, '(3f12.6)') R1/Angstrom2atomic
      write(outfileindex, '(3f12.6)') R2/Angstrom2atomic
      write(outfileindex, '(3f12.6)') R3/Angstrom2atomic
      write(outfileindex, '(300A6)') Origin_cell%Name_of_atomtype
      write(outfileindex, '(300i6)') Origin_cell%Num_atoms_eachtype
      write(outfileindex, '(a)')"Direct"
      do ia=1, Origin_cell%Num_atoms
         call rotate_newlattice(Origin_cell%Atom_position_direct(:, ia), Rt)
         call transformtohomecell(Rt)
         if(cpuid==0)write(outfileindex, '(3f12.6, a9)')Rt, trim(adjustl(Origin_cell%atom_name(ia)))
      enddo
      close(outfileindex)
   endif

   !> three lattice vectors in old cartesian coordinates
   Rua_newcell= R1
   Rub_newcell= R2
   Ruc_newcell= R3


   !> Setting the new cell defined by SURFACE card
   Cell_defined_by_surface%Rua = R1
   Cell_defined_by_surface%Rub = R2
   Cell_defined_by_surface%Ruc = R3
   Cell_defined_by_surface%lattice(:, 1)= R1
   Cell_defined_by_surface%lattice(:, 2)= R2
   Cell_defined_by_surface%lattice(:, 3)= R3

   Cell_defined_by_surface%Num_atoms = Origin_cell%Num_atoms
   Cell_defined_by_surface%max_projs = Origin_cell%max_projs
   Cell_defined_by_surface%NumberOfspinorbitals = Origin_cell%NumberOfspinorbitals
   Cell_defined_by_surface%Num_atom_type= Origin_cell%Num_atom_type
   allocate(Cell_defined_by_surface%Num_atoms_eachtype(Origin_cell%Num_atom_type))
   allocate(Cell_defined_by_surface%Name_of_atomtype(Origin_cell%Num_atom_type))
   allocate(Cell_defined_by_surface%itype_atom(Origin_cell%Num_atoms))
   allocate(Cell_defined_by_surface%Atom_name(Origin_cell%Num_atoms))
   allocate(Cell_defined_by_surface%Atom_position_cart  (3, Origin_cell%Num_atoms))
   allocate(Cell_defined_by_surface%Atom_position_direct(3, Origin_cell%Num_atoms))
   allocate(Cell_defined_by_surface%nprojs(Origin_cell%Num_atoms))
   allocate(Cell_defined_by_surface%spinorbital_to_atom_index(Origin_cell%NumberOfspinorbitals))
   allocate(Cell_defined_by_surface%spinorbital_to_projector_index(Origin_cell%NumberOfspinorbitals))
   Cell_defined_by_surface%Num_atoms_eachtype= Origin_cell%Num_atoms_eachtype
   Cell_defined_by_surface%Name_of_atomtype= Origin_cell%Name_of_atomtype
   Cell_defined_by_surface%itype_atom= Origin_cell%itype_atom
   Cell_defined_by_surface%Atom_name= Origin_cell%Atom_name
   Cell_defined_by_surface%nprojs= Origin_cell%nprojs
   do ia=1, Origin_cell%Num_atoms
      call rotate_newlattice(Origin_cell%Atom_position_direct(:, ia), Rt)
      call transformtohomecell(Rt)
      Atom_position_direct_newcell(:, ia)= Rt
      Cell_defined_by_surface%Atom_position_direct(:, ia)= Rt
      call direct_cart_real_newcell(Rt, Atom_position_cart_newcell(:, ia))
      Cell_defined_by_surface%Atom_position_cart(:, ia)= Atom_position_cart_newcell(:, ia)
   enddo

   !> try to find  one atom on the top surface according to the third coordinate of Atom_position_direct
   idummy= 0
   temp=Cell_defined_by_surface%Atom_position_direct(3, 1)
   if (topsurface_atom_index==0) then
      topsurface_atom_index= 1
      do ia=2, Origin_cell%Num_atoms
         if (Cell_defined_by_surface%Atom_position_direct(3, ia)>temp) then
            temp= Cell_defined_by_surface%Atom_position_direct(3, ia)
            topsurface_atom_index= ia
         endif
      enddo
   endif
   !> shift the top surface atom to 0.99 along R3'
   Rt= (/0d0, 0d0, &
      0.95d0-Cell_defined_by_surface%Atom_position_direct(3, topsurface_atom_index)/)
   call direct_cart_real(Rt, shift_to_topsurface_cart, Cell_defined_by_surface%lattice)

   !> shift all atoms by Rt
   do ia=1, Origin_cell%Num_atoms
      Rt2= Cell_defined_by_surface%Atom_position_direct(:, ia)+ Rt
      call transformtohomecell(Rt2)
      Cell_defined_by_surface%Atom_position_direct(:, ia)= Rt2
      call direct_cart_real(Rt2, Atom_position_cart_newcell(:, ia), Cell_defined_by_surface%lattice)
      Cell_defined_by_surface%Atom_position_cart(:, ia)= Atom_position_cart_newcell(:, ia)
   enddo

   !> write out
   if (cpuid==0) then
      write(stdout, '(a)') ' '
      write(stdout, '(a, i6, 1x, a)') " >> Top surface atom is ", topsurface_atom_index, &
         Cell_defined_by_surface%Atom_name(topsurface_atom_index)
      write(stdout, '(a)') ' '
   endif

   !> default wannier centers
   i= 0
   do ia= 1, Cell_defined_by_surface%Num_atoms
      do j= 1, Cell_defined_by_surface%nprojs(ia)
         i= i+ 1
         Cell_defined_by_surface%spinorbital_to_atom_index(i)= ia
         Cell_defined_by_surface%spinorbital_to_projector_index(i)= j
         if (SOC>0.or.Add_Zeeman_Field) then
            Cell_defined_by_surface%spinorbital_to_atom_index(i+NumberOfspinorbitals/2)= ia
            Cell_defined_by_surface%spinorbital_to_projector_index(i+NumberOfspinorbitals/2)= j
         endif
      enddo ! j
   enddo ! ia

   call writeout_poscar(Cell_defined_by_surface, 'POSCAR-SURFACE')

   !> generate POSCAR for slab system 
   call generate_slab_poscar(Cell_defined_by_surface)

   !> get the surface vector, we should set the new coordinate system
   !> set R1 to the new x direction ex'
   !> set R1\cross R2 to the new z direction ez'
   !> set ey'= ez'\cross ex'
   !> then e_i'= \sum_j U_ij e_j
   Urot= 0d0
   !> e_x'
   Urot(1, :)= R1/norm(R1)

   !> e_z'
   Urot(3, 1)= (R1(2)*R2(3)- R1(3)*R2(2))
   Urot(3, 2)= (R1(3)*R2(1)- R1(1)*R2(3))
   Urot(3, 3)= (R1(1)*R2(2)- R1(2)*R2(1))
   Urot(3, :)= Urot(3, :)/norm(Urot(3, :))

   !> e_y'= e_z'\cross e_x'
   Urot(2, 1)= (Urot(3, 2)*Urot(1, 3)- Urot(3, 3)*Urot(1, 2))
   Urot(2, 2)= (Urot(3, 3)*Urot(1, 1)- Urot(3, 1)*Urot(1, 3))
   Urot(2, 3)= (Urot(3, 1)*Urot(1, 2)- Urot(3, 2)*Urot(1, 1))
   Urot(2, :)= Urot(2, :)/norm(Urot(2, :))
   
   !> Here Rua_new, Origin_cell%Rub_new, Origin_cell%Ruc_new are vectors defined by SURFACE CARD in new coordinates
   call rotate(R1, Rua_new)
   call rotate(R2, Rub_new)
   call rotate(R3, Ruc_new)

   !> then transform R1, R2 to the new coordinates
   !> R1'_j= \sum_i U_ij R_i
   !> because the z direction is perpendicular to R1, R2,
   !> so the z coordinates for R1, R2 in the new axis are zero
   Ra2(1)= Urot(1, 1)*R1(1)+ Urot(1, 2)*R1(2)+ Urot(1, 3)*R1(3)
   Ra2(2)= Urot(2, 1)*R1(1)+ Urot(2, 2)*R1(2)+ Urot(2, 3)*R1(3)
   Rb2(1)= Urot(1, 1)*R2(1)+ Urot(1, 2)*R2(2)+ Urot(1, 3)*R2(3)
   Rb2(2)= Urot(2, 1)*R2(1)+ Urot(2, 2)*R2(2)+ Urot(2, 3)*R2(3)

   !> get the surface reciprocal vector
   cell_volume=Ra2(1)*Rb2(2)- Rb2(1)*Ra2(2)
   cell_volume= abs(cell_volume)

   if (abs(cell_volume)<1e-6) stop 'cell_volume equal zero'

   Ka2(1)= 2d0*pi/cell_volume*Rb2(2)
   Ka2(2)=-2d0*pi/cell_volume*Rb2(1)
   Kb2(1)=-2d0*pi/cell_volume*Ra2(2)
   Kb2(2)= 2d0*pi/cell_volume*Ra2(1)

   if (cpuid==0) then
      write(stdout, *)'2D Primitive Cell_Volume: ', Cell_Volume/Angstrom2atomic/Angstrom2atomic
      write(stdout, *)'Ra2, Rb2'
      write(stdout, '(3f10.4)')Ra2/Angstrom2atomic
      write(stdout, '(3f10.4)')Rb2/Angstrom2atomic
      write(stdout, *)'Ka2, Kb2'
      write(stdout, '(3f10.4)')ka2/Angstrom2atomic
      write(stdout, '(3f10.4)')kb2/Angstrom2atomic
   endif


!===============================================================================================================!
!> Set magnetic super cell
!===============================================================================================================!

   !> magnetic supercell stacks along Origin_cell%Ruc_newcell direction which is defined 
   !> as the third vector in SURFACE card
   !> The size of the supercell is Nslab
   !> Magnetic field is along Rua_newcell
   Rua_mag= Rua_newcell
   Rub_mag= Rub_newcell
   Ruc_mag= Ruc_newcell*Magq
   Magnetic_cell%Rua= Rua_mag
   Magnetic_cell%Rub= Rub_mag
   Magnetic_cell%Ruc= Ruc_mag

   Magnetic_cell%cell_parameters(1)= norm(Magnetic_cell%Rua)
   Magnetic_cell%cell_parameters(2)= norm(Magnetic_cell%Rub)
   Magnetic_cell%cell_parameters(3)= norm(Magnetic_cell%Ruc)
   Magnetic_cell%cell_parameters(4)= angle(Magnetic_cell%Rub, Magnetic_cell%Ruc)
   Magnetic_cell%cell_parameters(5)= angle(Magnetic_cell%Ruc, Magnetic_cell%Rua)
   Magnetic_cell%cell_parameters(6)= angle(Magnetic_cell%Rua, Magnetic_cell%Rub)

   !> transform lattice from direct space to reciprocal space
   call get_volume(Magnetic_cell%Rua, Magnetic_cell%Rub, Magnetic_cell%Ruc, Magnetic_cell%CellVolume)

   !> Volume of reciprocal lattice of magnetic supercell, in unit of (1/Bohr)**3
   Magnetic_cell%ReciprocalCellVolume= (2d0*pi)**3/Magnetic_cell%CellVolume

   !> Reciprocal lattice vectors in unit of 1/Bohr
   call get_reciprocal_lattice(Magnetic_cell%Rua, Magnetic_cell%Rub, Magnetic_cell%Ruc, &
                               Magnetic_cell%Kua, Magnetic_cell%Kub, Magnetic_cell%Kuc)

   if(cpuid==0)write(stdout, '(a)') '>> lattice information of the magnetic supercell (Angstrom)'
   if(cpuid==0)write(stdout, '(3f12.6)')Magnetic_cell%Rua/Angstrom2atomic
   if(cpuid==0)write(stdout, '(3f12.6)')Magnetic_cell%Rub/Angstrom2atomic
   if(cpuid==0)write(stdout, '(3f12.6)')Magnetic_cell%Ruc/Angstrom2atomic

   if(cpuid==0)write(stdout, '(a)') '>> Reciprocal lattice information of the magnetic supercell (1/Angstrom)'
   if(cpuid==0)write(stdout, '(3f12.6)')Magnetic_cell%Kua*Angstrom2atomic
   if(cpuid==0)write(stdout, '(3f12.6)')Magnetic_cell%Kub*Angstrom2atomic
   if(cpuid==0)write(stdout, '(3f12.6)')Magnetic_cell%Kuc*Angstrom2atomic

   MagneticSuperProjectedArea= Magnetic_cell%CellVolume/norm(R1)

   !> get atoms' position in the magnetic supercell
   Magnetic_cell%Num_atoms = Origin_cell%Num_atoms*Magq
   Magnetic_cell%max_projs = Origin_cell%max_projs*Magq
   Magnetic_cell%NumberOfspinorbitals = Origin_cell%NumberOfspinorbitals*Magq
   Magnetic_cell%Num_atom_type= Origin_cell%Num_atom_type
   allocate(Magnetic_cell%Num_atoms_eachtype(Origin_cell%Num_atom_type))
   allocate(Magnetic_cell%Name_of_atomtype(Origin_cell%Num_atom_type))
   allocate(Magnetic_cell%itype_atom(Magnetic_cell%Num_atoms))
   allocate(Magnetic_cell%Atom_name(Magnetic_cell%Num_atoms))
   allocate(Magnetic_cell%nprojs(Magnetic_cell%Num_atoms))
   allocate(Magnetic_cell%Atom_position_cart  (3, Magnetic_cell%Num_atoms))
   allocate(Magnetic_cell%Atom_position_direct(3, Magnetic_cell%Num_atoms))
   allocate(Magnetic_cell%wannier_centers_cart(3, Magnetic_cell%NumberOfspinorbitals))
   allocate(Magnetic_cell%wannier_centers_direct(3, Magnetic_cell%NumberOfspinorbitals))
   Magnetic_cell%Num_atoms_eachtype= Origin_cell%Num_atoms_eachtype*Magq
   Magnetic_cell%Name_of_atomtype= Origin_cell%Name_of_atomtype
   do iq=1, Magq
      do ia=1, Origin_cell%Num_atoms
         Magnetic_cell%itype_atom((iq-1)*Origin_cell%Num_atoms+ia)= Origin_cell%itype_atom(ia)
         Magnetic_cell%Atom_name((iq-1)*Origin_cell%Num_atoms+ia)= Origin_cell%Atom_name(ia)
      enddo
   enddo
   do iq=1, Magq
      do ia=1, Origin_cell%Num_atoms
         call rotate_newlattice(Origin_cell%Atom_position_direct(:, ia), Rt)
         Magnetic_cell%Atom_position_direct(1:2, (iq-1)*Origin_cell%Num_atoms+ia)= Rt(1:2)
         Magnetic_cell%Atom_position_direct(3, (iq-1)*Origin_cell%Num_atoms+ia)= (Rt(3)+ iq-1d0)/Magq
         call direct_cart_real_magneticcell(Magnetic_cell%Atom_position_direct(:, (iq-1)*Origin_cell%Num_atoms+ia), &
         Magnetic_cell%Atom_position_cart(:, (iq-1)*Origin_cell%Num_atoms+ia))
      enddo
      do i=1, NumberOfspinorbitals
         call rotate_newlattice(Origin_cell%wannier_centers_direct(:, i), Rt)
         Magnetic_cell%wannier_centers_direct(1:2, (iq-1)*NumberOfspinorbitals+i)= Rt(1:2)
         Magnetic_cell%wannier_centers_direct(3, (iq-1)*NumberOfspinorbitals+i)= (Rt(3)+ iq-1d0)/Magq
         call direct_cart_real_magneticcell(Magnetic_cell%wannier_centers_direct(:, (iq-1)*NumberOfspinorbitals+i), &
         Magnetic_cell%wannier_centers_cart(:, (iq-1)*NumberOfspinorbitals+i))
      enddo
   enddo

   Magnetic_cell%NumberOfspinorbitals= Origin_cell%NumberOfspinorbitals*Magq
   allocate(Magnetic_cell%spinorbital_to_atom_index(Magnetic_cell%NumberOfspinorbitals))
   allocate(Magnetic_cell%spinorbital_to_projector_index(Magnetic_cell%NumberOfspinorbitals))
   do iq=1, Magq
      istart= (iq-1)*Origin_cell%NumberOfspinorbitals+1
      iend  = iq*Origin_cell%NumberOfspinorbitals
      Magnetic_cell%spinorbital_to_atom_index(istart:iend)= Origin_cell%spinorbital_to_atom_index+istart-1
      Magnetic_cell%spinorbital_to_projector_index(istart:iend)= Origin_cell%spinorbital_to_projector_index+istart-1
   enddo

   do iq=1, Magq
      do ia= 1, Origin_cell%Num_atoms
         Magnetic_cell%nprojs(ia+(iq-1)*Origin_cell%Num_atoms)= Origin_cell%nprojs(ia)
      enddo
   enddo

   i=0
   do ia= 1, Magnetic_cell%Num_atoms
      do j= 1, Magnetic_cell%nprojs(ia)
         i= i+ 1
         Magnetic_cell%spinorbital_to_atom_index(i)= ia
         Magnetic_cell%spinorbital_to_projector_index(i)= j
         if (SOC>0.or.Add_Zeeman_Field) then
            Magnetic_cell%spinorbital_to_atom_index(i+Magq*NumberOfspinorbitals/2)= ia
            Magnetic_cell%spinorbital_to_projector_index(i+Magq*NumberOfspinorbitals/2)= j
         endif
      enddo
   enddo

   call writeout_poscar(Magnetic_cell, "POSCAR-mag")

   if (cpuid==0) then
      write(stdout, *)" "
      write(stdout, '(a, f16.8, a)')'3D Primitive Origin_cell%CellVolume: ', Origin_cell%CellVolume/(Angstrom2atomic**3), ' in Angstrom^3'
      write(stdout, '(a, f16.8, a)')'3D Magnetic supercell Volume: ', Magnetic_cell%CellVolume/(Angstrom2atomic**3), ' in Angstrom^3'
      write(stdout, '(a, f16.8, a)')'Projected area of magnetic supercell normal to the first vector specifed in SURFACE card: ',  &
         MagneticSuperProjectedArea, ' in Angstrom^2'
   endif

!===============================================================================================================!
!> KPATH_BULK card
!===============================================================================================================!


   !> read kpath_bulk information
   rewind(1001)
   lfound = .false.
   do while (.true.)
      read(1001, *, end= 104)inline
      inline=upper(inline)
      if (trim(adjustl(inline))=='KPATH_BULK') then
         lfound= .true.
         if (cpuid==0) write(stdout, *)' '
         if (cpuid==0) write(stdout, *)'We found KPATH_BULK card'
         exit
      endif
   enddo

   !> kline for 3d band stOrigin_cell%Ructure
   !> high symmetry k points
   read(1001, *) nk3lines
   if(cpuid==0)write(stdout, '(a, 40i5)')'Number of K lines : ', nk3lines
   allocate(k3line_start(3, nk3lines))
   allocate(k3line_end(3, nk3lines))
   allocate(k3line_name(nk3lines+1))
   allocate(k3line_stop(nk3lines+1))
   allocate(k3line_mag_stop(nk3lines+1))
   allocate(k3line_unfold_stop(nk3lines+1))
   k3line_mag_stop= 0d0
   k3line_unfold_stop= 0d0
   k3line_stop= 0d0
   k3line_start= 0d0
   k3line_end= 0d0
   k3line_name= ' '
   it=0
   do i=1, nk3lines
      read(1001, *, err=201) k3line_name(i), k3line_start(:, i), &
         char_temp, k3line_end(:, i)
      it= it+ 1
      if(cpuid==0)write(stdout, '(a5, 3f9.4, 2x, a5, 3f9.4)')&
         k3line_name(i), k3line_start(:, i), &
         char_temp, k3line_end(:, i)

   enddo
201 continue
   if (it< nk3lines.and.cpuid==0) then
      write(stdout, *)' '
      write(stdout, *)' >>>> Error happens in KPATH_BULK card'
      write(stdout, *)' Error: the number of kpath lines should consistent'
      write(stdout, *)' Nk3lines ', nk3lines, ' the k lines you given ', it
      write(stdout, *)' Please set nk3lines to be  ', it
      stop
   endif

   k3line_name(nk3lines+1)= char_temp

   NN= Nk
   nk3_band= NN*nk3lines

   allocate(k3len(nk3_band))
   allocate(k3len_mag(nk3_band))
   allocate(k3len_unfold(nk3_band))
   allocate(k3points(3, nk3_band))
   k3len=0d0
   k3len_mag=0d0
   k3len_unfold=0d0
   k3points= 0d0
   t1= 0d0
   do j=1, nk3lines
      do i=1, NN
         kstart= k3line_start(:, j)
         kend  = k3line_end(:, j)
         k1= kstart(1)*Origin_cell%Kua+ kstart(2)*Origin_cell%Kub+ kstart(3)*Origin_cell%Kuc
         k2= kend(1)*Origin_cell%Kua+ kend(2)*Origin_cell%Kub+ kend(3)*Origin_cell%Kuc
         !k1= kstart
         !k2= kend

         k3points(:, i+ (j-1)*NN)= kstart+ (kend- kstart)*dble(i-1)/dble(NN-1)

         temp= dsqrt((k2(1)- k1(1))**2 &
            +(k2(2)- k1(2))**2  &
            +(k2(3)- k1(3))**2)/dble(NN-1)

         if (i.gt.1) then
            t1=t1+temp
         endif
         k3len(i+(j-1)*NN)= t1
      enddo
      k3line_stop(j+1)= t1
   enddo

   !> for magnetic supercell
   t1=0
   do j=1, nk3lines
      do i=1, NN
         kstart= k3line_start(:, j)
         kend  = k3line_end(:, j)
         k1= kstart(1)*Magnetic_cell%Kua+ kstart(2)*Magnetic_cell%Kub+ kstart(3)*Magnetic_cell%Kuc
         k2= kend(1)*Magnetic_cell%Kua+ kend(2)*Magnetic_cell%Kub+ kend(3)*Magnetic_cell%Kuc

         temp= dsqrt((k2(1)- k1(1))**2 &
            +(k2(2)- k1(2))**2  &
            +(k2(3)- k1(3))**2)/dble(NN-1)

         if (i.gt.1) then
            t1=t1+temp
         endif
         k3len_mag(i+(j-1)*NN)= t1
      enddo
      k3line_mag_stop(j+1)= t1
   enddo

   !> for unfolded cell

   !> for magnetic supercell
   t1=0
   do j=1, nk3lines
      do i=1, NN
         kstart= k3line_start(:, j)
         kend  = k3line_end(:, j)
         k1= kstart(1)*Folded_cell%Kua+ kstart(2)*Folded_cell%Kub+ kstart(3)*Folded_cell%Kuc
         k2= kend(1)*Folded_cell%Kua+ kend(2)*Folded_cell%Kub+ kend(3)*Folded_cell%Kuc

         temp= dsqrt((k2(1)- k1(1))**2 &
            +(k2(2)- k1(2))**2  &
            +(k2(3)- k1(3))**2)/dble(NN-1)

         if (i.gt.1) then
            t1=t1+temp
         endif
         k3len_unfold(i+(j-1)*NN)= t1
      enddo
      k3line_unfold_stop(j+1)= t1
   enddo



104 continue
   if (.not.lfound .and. (BulkBand_line_calc.or.LandauLevel_k_calc)) then
      stop 'ERROR: please set KPATH_BULK for bulk band stOrigin_cell%Ructure calculation'
   endif

!===============================================================================================================!
!> KPATH_SLAB card
!===============================================================================================================!

   !> read kpath_slab information
   rewind(1001)
   lfound = .false.
   do while (.true.)
      read(1001, *, end= 105)inline
      inline=upper(inline)
      if (trim(adjustl(inline))=='KPATH_SLAB') then
         lfound= .true.
         if (cpuid==0) write(stdout, *)' '
         if (cpuid==0) write(stdout, *)'We found KPATH_SLAB card'
         exit
      endif
   enddo

   !> read in k lines for 2D system
   k2line_name= ' '
   if (cpuid==0) write(stdout, *)'k lines for 2D system'
   read(1001, *)nk2lines
   if (cpuid==0) write(stdout, *)'No. of k lines for 2D system : ', nk2lines
   it = 0
   do i=1, nk2lines
      read(1001, *, err=202) k2line_name(i), kp(:, i), &
         char_temp, ke(:, i)
      it= it+ 1
      if (cpuid==0) write(stdout, '(a6, 2f9.5, 4x, a6, 2f9.5)')&
         k2line_name(i), kp(:, i), &
         char_temp, ke(:, i)
   enddo
202 continue

   if (it< nk2lines.and.cpuid==0) then
      write(stdout, *)' '
      write(stdout, *)' >>>> Error happens in KPATH_SLAB card'
      write(stdout, *)' Error: the number of kpath lines should consistent'
      write(stdout, *)' Nk2lines ', nk2lines, ' the k lines you given ', it
      write(stdout, *)' Please set nk2lines to be  ', it
      stop
   endif

   k2line_name(nk2lines+1) = char_temp


   NN= Nk
   knv2= NN*nk2lines
   allocate( k2_path(knv2, 2))
   allocate( k2len (knv2))
   k2_path= 0d0
   k2len= 0d0

   t1=0d0
   k2len=0d0
   k2line_stop= 0d0
   do j=1, nk2lines
      do i=1, NN
         kstart(1:2)= kp(:, j)
         kend(1:2)  = ke(:, j)
         k1(1:2)= kstart(1)*Ka2+ kstart(2)*Kb2
         k2(1:2)= kend(1)*Ka2+ kend(2)*Kb2
         k2_path(i+(j-1)*NN,:)= kstart(1:2)+ (kend(1:2)-kstart(1:2))*(i-1)/dble(NN-1)

         temp= dsqrt((k2(1)- k1(1))**2 &
            +(k2(2)- k1(2))**2)/dble(NN-1)

         if (i.gt.1) then
            t1=t1+temp
         endif
         k2len(i+(j-1)*NN)= t1
      enddo
      k2line_stop(j+1)= t1

   enddo


105 continue
   if (.not.lfound .and.(SlabBand_calc .or. SlabSS_calc)) then
      stop 'ERROR: please set KPATH_SLAB for slab band stOrigin_cell%Ructure calculation'
   endif


   !> read kplane_slab information
   !> default value for KPLANE_SLAB
   K2D_start= (/-0.5, -0.5/)
   K2D_vec1 = (/ 1.0,  0.0/)
   K2D_vec2 = (/ 0.0,  1.0/)

   rewind(1001)
   lfound = .false.
   do while (.true.)
      read(1001, *, end= 106)inline
      inline=upper(inline)
      if (trim(adjustl(inline))=='KPLANE_SLAB') then
         lfound= .true.
         if (cpuid==0) write(stdout, *)' '
         if (cpuid==0) write(stdout, *)'We found KPLANE_SLAB card'
         exit
      endif
   enddo

   !> kpoints plane for 2D system--> arcs
   it = 0
   read(1001, *, err=203)K2D_start
   it= it+ 1
   read(1001, *, err=203)K2D_vec1
   it= it+ 1
   read(1001, *, err=203)K2D_vec2
   it= it+ 1
203 continue
   if (it< 3.and.cpuid==0) then
      write(stdout, *)' '
      write(stdout, *)' >>>> Error happens in KPLANE_SLAB card'
      write(stdout, *)' Error: There are three lines in this card to specify the start point'
      write(stdout, *)' of vectors, and two lines to assign two vectors.'
      write(stdout, *)" If you don't know the meaning of this card, please delete this card"
      stop
   endif


106 continue

   if (cpuid==0) write(stdout, *)'>> Kpoints plane for 2D system--> arcs  '
   if (cpuid==0) write(stdout, '((a, 2f8.4))')'K2D_start:', K2D_start
   if (cpuid==0) write(stdout, '((a, 2f8.4))')'The first vector: ', K2D_vec1
   if (cpuid==0) write(stdout, '((a, 2f8.4))')'The second vector: ', K2D_vec2
   if (.not.lfound .and.(SlabArc_calc  .or. SlabSpintexture_calc)) then
      stop 'ERROR: please set KPLANE_SLAB for arc or spintexture calculations'
   endif


!===============================================================================================================!
!> KPLANE_BULK card
!===============================================================================================================!

   !> read kplane_bulk information
   !> default value for KPLANE_BULK
   K3D_start= (/ 0.0,  0.0,   0.0/)
   K3D_vec1 = (/ 1.0,  0.0,   0.0/)
   K3D_vec2 = (/ 0.0,  0.5,   0.0/)

   rewind(1001)
   lfound = .false.
   do while (.true.)
      read(1001, *, end= 107)inline
      inline=upper(inline)
      if (trim(adjustl(inline))=='KPLANE_BULK') then
         lfound= .true.
         if (cpuid==0) write(stdout, *)' '
         if (cpuid==0) write(stdout, *)'We found KPLANE_BULK card'
         exit
      endif
   enddo

   !> check whether we have a line to determine the coordinates
   read(1001, *, end= 107)inline
   inline=trim(adjustl(inline))
   if (index(inline, 'C')==0.and.index(inline, 'D')==0 &
       .and.index(inline, 'c')==0.and.index(inline, 'd')==0)then
      DirectOrCart='D'

      rewind(1001)
      do while (.true.)
         read(1001, *, end= 107)inline
         inline=upper(inline)
         if (trim(adjustl(inline))=='KPLANE_BULK') then
            lfound= .true.
            exit
         endif
      enddo
      goto 2061
   elseif (index(inline, 'c')==1.or.index(inline, 'C')==1)then
      DirectOrCart='C'
   else
      DirectOrCart='D'
   endif

2061 continue

   !> kpoints plane for 3D system--> gapshape
   it= 0
   read(1001, *, err=206)K3D_start
   it= it+ 1
   read(1001, *, err=206)K3D_vec1
   it= it+ 1
   read(1001, *, err=206)K3D_vec2
   it= it+ 1

   if (index(DirectOrCart, 'C')>0) then
      call cart_direct_rec(K3D_start, k1)
      K3D_start= k1
      call cart_direct_rec(K3D_vec1, k1)
      K3D_vec1= k1
      call cart_direct_rec(K3D_vec2, k1)
      K3D_vec2= k1
   endif

206 continue
   if (it< 3.and.cpuid==0) then
      write(stdout, *)' '
      write(stdout, *)' >>>> Error happens in KPLANE_BULK card'
      write(stdout, *)' Error: There are three lines in this card to specify the start point'
      write(stdout, *)' of vectors, and two lines to assign two vectors.'
      write(stdout, *)" If you don't know the meaning of this card, please delete this card"
      stop
   endif


107 continue

   if (cpuid==0) write(stdout, *)'>> Kpoints plane for 3D system--> gapshape  '
   if (cpuid==0) write(stdout, '((a, 3f8.4))')'k3D_start : ', K3D_start
   if (cpuid==0) write(stdout, '((a, 3f8.4))')'The 1st vector: ', K3D_vec1
   if (cpuid==0) write(stdout, '((a, 3f8.4))')'The 2nd vector: ', K3D_vec2
   if (.not.lfound .and.(BulkGap_plane_calc  .or. wanniercenter_calc)) then
      stop 'ERROR: please set KPLANE_bulk for gap or WCC calculations'
   endif

!===============================================================================================================!
!> KCUBE_BULK card
!===============================================================================================================!

   !> read kcube_bulk information
   !> default value for KCUBE_BULK
   K3D_start_cube= (/ 0.0,  0.0,   0.0/)
   K3D_vec1_cube = (/ 1.0,  0.0,   0.0/)
   K3D_vec2_cube = (/ 0.0,  1.0,   0.0/)
   K3D_vec3_cube = (/ 0.0,  0.0,   1.0/)

   rewind(1001)
   lfound = .false.
   do while (.true.)
      read(1001, *, end= 108)inline
      inline= upper(inline)
      if (trim(adjustl(inline))=='KCUBE_BULK') then
         lfound= .true.
         if (cpuid==0) write(stdout, *)' '
         if (cpuid==0) write(stdout, *)'We found KCUBE_BULK card'
         exit
      endif
   enddo

   !> kpoints plane for 3D system--> gapshape
   it= 0
   read(1001, *, err=204)K3D_start_cube
   it= it+ 1
   read(1001, *, err=204)K3D_vec1_cube
   it= it+ 1
   read(1001, *, err=204)K3D_vec2_cube
   it= it+ 1
   read(1001, *, err=204)K3D_vec3_cube
   it= it+ 1
204 continue
   if (it< 3.and.cpuid==0) then
      write(stdout, *)' '
      write(stdout, *)' >>>> Error happens in KCUBE_BULK card'
      write(stdout, *)' Error: There are four lines in this card to specify the start point'
      write(stdout, *)' of vectors, and three lines to assign two vectors.'
      write(stdout, *)" If you don't know the meaning of this card, please delete this card"
      stop
   endif
108 continue

   kCubeVolume= K3D_vec1_cube(1)*(K3D_vec2_cube(2)*K3D_vec3_cube(3) &
      - K3D_vec2_cube(3)*K3D_vec3_cube(2)) &
      + K3D_vec1_cube(2)*(K3D_vec2_cube(3)*K3D_vec3_cube(1) &
      - K3D_vec2_cube(1)*K3D_vec3_cube(3)) &
      + K3D_vec1_cube(3)*(K3D_vec2_cube(1)*K3D_vec3_cube(2) &
      - K3D_vec2_cube(2)*K3D_vec3_cube(1))

   kCubeVolume= kCubeVolume*Origin_cell%ReciprocalCellVolume


   if (cpuid==0) write(stdout, *)'>> Kpoints cube for 3D system--> gapshape3D  '
   if (cpuid==0) write(stdout, '((a, 3f8.4))')'k3D_start :', K3D_start_cube
   if (cpuid==0) write(stdout, '((a, 3f8.4))')'The 1st vector: ', K3D_vec1_cube
   if (cpuid==0) write(stdout, '((a, 3f8.4))')'The 2nd vector: ', K3D_vec2_cube
   if (cpuid==0) write(stdout, '((a, 3f8.4))')'The 3rd vector: ', K3D_vec3_cube
   if (cpuid==0) write(stdout, '((a, 3f8.4))')'kCubeVolume: ', kCubeVolume*Angstrom2atomic**3
   if (cpuid==0) write(stdout, '((a, 3f8.4))')'ReciprocalOrigin_cell%CellVolume: ', &
   Origin_cell%ReciprocalCellVolume*Angstrom2atomic**3
   if (.not.lfound .and.(BulkGap_cube_calc)) then
      stop 'ERROR: please set KCUBE_BULK for gap3D calculations'
   endif

!===============================================================================================================!
!> KPATH_BERRY card
!===============================================================================================================!

   !> set default parameters for Berry phase calculation
   NK_Berry= 2
   allocate(k3points_Berry(3, NK_Berry))
   DirectOrCart_Berry='Direct'
   k3points_Berry= 0d0
   rewind(1001)
   lfound = .false.
   do while (.true.)
      read(1001, *, end= 113)inline
      inline=upper(inline)
      if (trim(adjustl(inline))=='KPATH_BERRY') then
         lfound= .true.
         if (cpuid==0) write(stdout, *)' '
         if (cpuid==0) write(stdout, *)'We found KPATH_BERRY card'
         exit
      endif
   enddo

   read(1001, *, end=208, err=208)NK_Berry
   if (cpuid==0) write(stdout, '(a, i10)')'NK_Berry', NK_Berry
   read(1001, *, end=208, err=208)inline   ! The unit of lattice vector
   DirectOrCart_Berry= trim(adjustl(inline))

   deallocate(k3points_Berry)
   allocate(k3points_Berry(3, NK_Berry))
   k3points_Berry= 0d0

   it= 0
   if (index(DirectOrCart_Berry, "D")>0)then
      do ik=1, NK_Berry
         read(1001, *, end=208, err=208)k3points_Berry(:, ik)   ! The unit of lattice vector
         it = it+ 1
      enddo
   else
      do ik=1, NK_Berry
         read(1001, *, end=208, err=208)k    ! The unit of lattice vector
         call cart_direct_rec(k, k3points_Berry(:, ik))
         it = it+ 1
      enddo
   endif  ! Direct or Cart coordinates
208 continue

   if (it< NK_Berry.and. cpuid==0) then
      write(stdout, *)"ERROR: something wrong in the KPATH_BERRY card"
      write(stdout, *)"No. of kpoints for berry is not consistent with No. of lines"
      write(stdout, *)"I found ", it, " lines"
      write(stdout, *)"while you set NK_Berry to be ", NK_Berry
      stop
   endif

113 continue

   if (cpuid==0) write(stdout, *)' '
   if (.not.lfound.and.cpuid==0.and.BerryPhase_calc)then
      write(stdout, *)'Error : you have to set KPATH_BERRY card with a list of k points'
      stop
   endif

   !< end of Berry phase setting


!===============================================================================================================!
!> KPOINT_BULK card
!===============================================================================================================!


   !> default parameters for KPOINT
   Kpoint_3D_direct = 0
   Kpoint_3D_cart = 0
   rewind(1001)
   lfound = .false.
   do while (.true.)
      read(1001, *, end= 114)inline
      inline= upper(inline)
      if (trim(adjustl(inline))=='KPOINT_BULK') then
         lfound= .true.
         if (cpuid==0) write(stdout, *)' '
         if (cpuid==0) write(stdout, *)'We found KPOINT_BULK card'
         if (cpuid==0) write(stdout, *)'This card should like this:'
         if (cpuid==0) write(stdout, *)'KPOINT_BULK'
         if (cpuid==0) write(stdout, *)'direct'
         if (cpuid==0) write(stdout, *)'0.0 0.0 0.0  ! three real number'
         exit
      endif
   enddo

   read(1001, *)inline   ! The unit of lattice vector
   inline=upper(inline)
   DirectOrCart= trim(adjustl(inline))
   if (index(DirectOrCart, "D")>0)then
      read(1001, *)Kpoint_3D_direct
      call direct_cart_rec(Kpoint_3D_direct, Kpoint_3D_cart)
   else
      read(1001, *)Kpoint_3D_cart
      call cart_direct_rec(Kpoint_3D_cart, Kpoint_3D_direct)
   endif

114 continue
   if (cpuid==0) write(stdout, *)' '
   if (.not.lfound.and.cpuid==0)write(stdout, *)'>> Using default parameters for Kpoint_3D'
   if (cpuid==0) then
      write(stdout, '(a, 3f7.4, a)')'>> k points ', Kpoint_3D_direct, " in unit of reciprocal primitive cell"
      write(stdout, '(a, 3f7.4, a)')'>> k points ', Kpoint_3D_cart*Angstrom2atomic, " in Cartesian coordinates"
   endif

!===============================================================================================================!
!> EFFECTIVE_MASS card
!===============================================================================================================!

   !>> setting up effective mass calculation
   !> default parameters for effective mass calculation
   dk_mass= 0.01  ! in unit of 1/Ang
   iband_mass= NumOccupied
   k_mass= 0
   rewind(1001)
   lfound = .false.
   do while (.true.)
      read(1001, *, end= 109)inline
      inline= upper(inline)
      if (trim(adjustl(inline))=='EFFECTIVE_MASS') then
         lfound= .true.
         if (cpuid==0) write(stdout, *)' '
         if (cpuid==0) write(stdout, *)'We found EFFECTIVE_MASS card'
         exit
      endif
   enddo

   it= 0
   read(1001, *, end=205, err=205)iband_mass
   it= it+ 1
   read(1001, *, end=205, err=205)dk_mass
   it= it+ 1
   read(1001, *, end=205, err=205)k_mass
   it= it+ 1
205 continue
   if (it< 3.and.cpuid==0) then
      write(stdout, *)' '
      write(stdout, *)' >>>> Error happens in EFFECTIVE_MASS card'
      write(stdout, *)' Error: There are three lines in this card to specify the iband_mass'
      write(stdout, *)" , dk_mass and k_mass, like this: "
      write(stdout, *)"EFFECTIVE_MASS"
      write(stdout, *)" 6       ! the 6'th band"
      write(stdout, *)" 0.01    ! in unit of 1/Bohr"
      write(stdout, *)" 0 0 0   ! k point"
      stop
   endif


109 continue
   dk_mass= dk_mass/Angstrom2atomic
   if (cpuid==0) write(stdout, *)' '
   if (.not.lfound.and.cpuid==0)write(stdout, *)'>> Using default parameters for effective mass calculation'
   if (cpuid==0) write(stdout, *)'>> Effective mass calculation parameters  '
   if (cpuid==0) write(stdout, '(a, i5, a)')'>> The ', iband_mass, "'th band"
   if (cpuid==0) write(stdout, '(a, f7.4, a)')'>> k step ', dk_mass*Angstrom2atomic, " in unit of 1/Angstrom"
   if (cpuid==0) write(stdout, '(a, 3f7.4, a)')'>> k points ', k_mass, " in unit of reciprocal primitive cell"
   k1=k_mass
   call direct_cart_rec(k1, k_mass)
   if (cpuid==0) write(stdout, '(a, 3f7.4, a)')'>> k points ', k_mass*Angstrom2atomic, " in unit of 1/Angstrom"

!===============================================================================================================!
!> KPOINTS_3D card
!===============================================================================================================!

     !>> setting up a series of k points in 3D BZ
     rewind(1001)
     lfound = .false.
     do while (.true.)
        read(1001, *, end= 321)inline
        inline=upper(inline)
        if (trim(adjustl(inline))=='KPOINTS_3D') then
           lfound= .true.
           if (cpuid==0) write(stdout, *)' '
           if (cpuid==0) write(stdout, *)'We found KPOINTS_3D card'
           exit
        endif
     enddo

     read(1001, *, end=319, err=319, iostat=stat)Nk3_point_mode   ! The unit of lattice vector
     allocate(k3points_pointmode_cart(3, Nk3_point_mode))
     allocate(k3points_pointmode_direct(3, Nk3_point_mode))
     k3points_pointmode_cart= 0d0
     k3points_pointmode_direct= 0d0

     read(1001, *, end=319, err=319, iostat=stat)inline   ! The unit of lattice vector
     inline= upper(inline)
     if (index(trim(adjustl(inline)), "D")>0)then
        do ik= 1, Nk3_point_mode
           read(1001, *, end=319, err=319, iostat=stat)k3points_pointmode_direct(:, ik)
           call direct_cart_rec(k3points_pointmode_direct(:, ik), &
              k3points_pointmode_cart(:, ik))
        enddo
     else
        do ik= 1, Nk3_point_mode
           read(1001, *, end=319, err=319, iostat=stat)k3points_pointmode_cart(:, ik)
           call cart_direct_rec(k3points_pointmode_cart(:, ik), &
              k3points_pointmode_direct(:, ik))
        enddo
     endif
     
     319 continue
     if (stat/=0 .and. cpuid==0) then
        write(stdout, '(8f10.5)') "ERROR: there are something wrong in KPOINTS_3D card"
        write(stdout, '(8f10.5)') "It should be like this:"
        write(stdout, '(8f10.5)') "The number of lines below 'Direct' "
        write(stdout, '(8f10.5)') "should be the same as the number of k points"
        write(stdout, '(8f10.5)') "KPOINTS_3D"
        write(stdout, '(8f10.5)') "4 ! number of k points"
        write(stdout, '(8f10.5)') "Direct"
        write(stdout, '(8f10.5)') "0.0  0.0  0.0"
        write(stdout, '(8f10.5)') "0.5  0.0  0.0"
        write(stdout, '(8f10.5)') "0.0  0.5  0.0"
        write(stdout, '(8f10.5)') "0.0  0.0  0.5"
        stop
     endif

     !> print out the single kpoint positions
     if (cpuid==0) then
        write(stdout, '(a)')" "
        write(stdout, '(a)')"KPOINTS_3D positions"
        write(stdout, '(8a10)')"index", "kx", 'ky', 'kz', 'k1', 'k2', 'k3'
        do ik=1, Nk3_point_mode
           write(stdout, '(i8,4x,8f10.5)')ik, k3points_pointmode_cart(:, ik)*Angstrom2atomic, k3points_pointmode_direct(:, ik)
        enddo
        write(stdout, '(a)')" "
     endif

     321 continue
     if (cpuid==0) write(stdout, *)' '
     if (.not. lfound) then
        Nk3_point_mode = 1
        allocate(k3points_pointmode_cart(3, Nk3_point_mode))
        allocate(k3points_pointmode_direct(3, Nk3_point_mode))
        k3points_pointmode_cart= 0d0
        k3points_pointmode_direct= 0d0
     endif
     if (.not.lfound.and.cpuid==0)write(stdout, *)'>> We use the default values for k3points_pointmode_direct=[0,0,0]' 
     if (.not.lfound.and.cpuid==0)write(stdout, *)'>> and Nk3_point_mode = 1'

!===============================================================================================================!
!> KPOINTS_FOLD_3D card
!===============================================================================================================!

     !>> setting up a series of k points in the folded 3D BZ
     rewind(1001)
     lfound = .false.
     do while (.true.)
        read(1001, *, end= 3210)inline
        inline=upper(inline)
        if (trim(adjustl(inline))=='KPOINTS_FOLD_3D'.or.trim(adjustl(inline))=='KPOINTS_FOLDED_3D') then
           lfound= .true.
           if (cpuid==0) write(stdout, *)' '
           if (cpuid==0) write(stdout, *)'We found KPOINTS_FOLD_3D card'
           exit
        endif
     enddo

     read(1001, *, end=3190, err=3190, iostat=stat)Nk3_unfold_point_mode   ! The unit of lattice vector
     allocate(k3points_unfold_pointmode_cart(3, Nk3_unfold_point_mode))
     allocate(k3points_unfold_pointmode_direct(3, Nk3_unfold_point_mode))
     k3points_unfold_pointmode_cart= 0d0
     k3points_unfold_pointmode_direct= 0d0

     read(1001, *, end=3190, err=3190, iostat=stat)inline   ! The unit of lattice vector
     inline=upper(inline)
     if (index(trim(adjustl(inline)), "D")>0)then
        do ik= 1, Nk3_unfold_point_mode
           read(1001, *, end=3190, err=3190, iostat=stat)k3points_unfold_pointmode_direct(:, ik)
           call direct_cart_rec_unfold(k3points_unfold_pointmode_direct(:, ik), &
              k3points_unfold_pointmode_cart(:, ik))
        enddo
     else
        do ik= 1, Nk3_point_mode
           read(1001, *, end=3190, err=3190, iostat=stat)k3points_unfold_pointmode_cart(:, ik)
           call cart_direct_rec_unfold(k3points_unfold_pointmode_cart(:, ik), &
              k3points_unfold_pointmode_direct(:, ik))
        enddo
     endif
     
     3190 continue
     if (stat/=0 .and. cpuid==0) then
        write(stdout, '(8f10.5)') "ERROR: there are something wrong in KPOINTS_FOLD_3D card"
        write(stdout, '(8f10.5)') "It should be like this:"
        write(stdout, '(8f10.5)') "The number of lines below 'Direct' "
        write(stdout, '(8f10.5)') "should be the same as the number of k points"
        write(stdout, '(8f10.5)') "KPOINTS_FOLD_3D"
        write(stdout, '(8f10.5)') "4 ! number of k points"
        write(stdout, '(8f10.5)') "Direct"
        write(stdout, '(8f10.5)') "0.0  0.0  0.0"
        write(stdout, '(8f10.5)') "0.5  0.0  0.0"
        write(stdout, '(8f10.5)') "0.0  0.5  0.0"
        write(stdout, '(8f10.5)') "0.0  0.0  0.5"
     endif

     !> print out the set of kpoints' positions
     if (cpuid==0) then
        write(stdout, '(a)')" "
        write(stdout, '(a)')"KPOINTS_FOLD_3D positions"
        write(stdout, '(8a10)')"index", "kx", 'ky', 'kz', 'k1', 'k2', 'k3'
        do ik=1, Nk3_unfold_point_mode
           write(stdout, '(i8, 4x,8f10.5)')ik, k3points_unfold_pointmode_cart(:, ik)*Angstrom2atomic, &
              k3points_unfold_pointmode_direct(:, ik)
        enddo
        write(stdout, '(a)')" "
     endif

     3210 continue
     if (cpuid==0) write(stdout, *)' '
     if (.not. lfound) then
        Nk3_unfold_point_mode = 1
        allocate(k3points_unfold_pointmode_cart(3, Nk3_unfold_point_mode))
        allocate(k3points_unfold_pointmode_direct(3, Nk3_unfold_point_mode))
        k3points_unfold_pointmode_cart= 0d0
        k3points_unfold_pointmode_direct= 0d0
     endif
     if (.not.lfound.and.cpuid==0)write(stdout, *)'>> We use the default values for k3points_unfold_pointmode_direct=[0,0,0]' 
     if (.not.lfound.and.cpuid==0)write(stdout, *)'>> and Nk3_unfold_point_mode = 1'


!===============================================================================================================!
!> SINGLEKPOINT_2D card
!===============================================================================================================!

   !>> setting up a single k points in 2D BZ
   Single_KPOINT_2D_CART= [0.d0, 0d0]
   Single_KPOINT_2D_DIRECT= [0.d0, 0d0]
   rewind(1001)
   lfound = .false.
   do while (.true.)
      read(1001, *, end= 308)inline
      inline=upper(inline)
      if (trim(adjustl(inline))=='SINGLEKPOINT_2D') then
         lfound= .true.
         if (cpuid==0) write(stdout, *)' '
         if (cpuid==0) write(stdout, *)'We found SINGLEKPOINT_2D card'
         exit
      endif
   enddo

   read(1001, *, end=307, err=307, iostat=stat)inline   ! The unit of lattice vector
   inline=upper(inline)
   DirectOrCart_SINGLE= trim(adjustl(inline))
   if (index(DirectOrCart_SINGLE, "D")>0)then
      read(1001, *, end=307, err=307, iostat=stat)Single_KPOINT_2D_DIRECT(1:2)
   else
      stop " for SINGLEKPOINT_2D, we only support Direct coordinates"
   endif

307 continue
   if (stat/=0 .and. cpuid==0) then
      write(stdout, '(8f10.5)') "ERROR: there is something wrong in SINGLEKPOINT_2D card"
      write(stdout, '(8f10.5)') "It should be like this:"
      write(stdout, '(8f10.5)') "SINGLEKPOINT_2D"
      write(stdout, '(8f10.5)') "Direct"
      write(stdout, '(8f10.5)') "0.0  0.0"
   endif

   !> print out the single kpoint positions
   if (cpuid==0) then
      write(stdout, '(a)')" "
      write(stdout, '(a)')"SingleKPOINT_2D positions in fractional coordinates"
      write(stdout, '(8a10)')'k1', 'k2', 'k3'
      write(stdout, '(8f10.5)')Single_KPOINT_2D_DIRECT
      write(stdout, '(a)')" "
   endif

308 continue
   if (cpuid==0) write(stdout, *)' '
   if (.not.lfound.and.cpuid==0)write(stdout, *)'>> We use the default values for Single_KPOINT_2D_DIRECT=[0,0]'



!===============================================================================================================!
!> SINGLEKPOINT_3D card
!===============================================================================================================!


   !>> setting up a single k points in 3D BZ
   Single_KPOINT_3D_CART= [0.d0, 0d0, 0d0]
   Single_KPOINT_3D_DIRECT= [0.d0, 0d0, 0d0]
   rewind(1001)
   lfound = .false.
   do while (.true.)
      read(1001, *, end= 311)inline
      inline=upper(inline)
      if (trim(adjustl(inline))=='SINGLEKPOINT_3D') then
         lfound= .true.
         if (cpuid==0) write(stdout, *)' '
         if (cpuid==0) write(stdout, *)'We found SINGLEKPOINT_3D card'
         exit
      endif
   enddo

   read(1001, *, end=309, err=309, iostat=stat)inline   ! The unit of lattice vector
   inline=upper(inline)
   DirectOrCart_SINGLE= trim(adjustl(inline))
   if (index(DirectOrCart_SINGLE, "D")>0)then
      read(1001, *, end=309, err=309, iostat=stat)Single_KPOINT_3D_DIRECT(1:3)
      call direct_cart_rec(Single_KPOINT_3D_DIRECT, Single_KPOINT_3D_CART)
   else
      read(1001, *, end=309, err=309, iostat=stat)Single_KPOINT_3D_CART(1:3)
      call cart_direct_rec(Single_KPOINT_3D_CART, Single_KPOINT_3D_DIRECT)
   endif

309 continue
   if (stat/=0 .and. cpuid==0) then
      write(stdout, '(8f10.5)') "ERROR: there is something wrong in SINGLEKPOINT_3D card"
      write(stdout, '(8f10.5)') "It should be like this:"
      write(stdout, '(8f10.5)') "SINGLEKPOINT_3D"
      write(stdout, '(8f10.5)') "Direct"
      write(stdout, '(8f10.5)') "0.0  0.0  0.0"
   endif

   !> print out the single kpoint positions
   if (cpuid==0) then
      write(stdout, '(a)')" "
      write(stdout, '(a)')"Single_KPOINT_3D positions"
      write(stdout, '(8a10)')"kx", 'ky', 'kz', 'k1', 'k2', 'k3'
      write(stdout, '(8f10.5)') Single_KPOINT_3D_CART*Angstrom2atomic, Single_KPOINT_3D_DIRECT
      write(stdout, '(a)')" "
   endif

311 continue
   if (cpuid==0) write(stdout, *)' '
   if (.not.lfound.and.cpuid==0)write(stdout, *)'>> We use the default values for Single_KPOINT_3D_DIRECT=[0,0,0]'


!===============================================================================================================!
!> SURFACE_ATOMS card
!===============================================================================================================!

   !> setup the atoms on top and bottom surface that used for output the
   !> surface-state spectrum
   !> by default we output all the atoms' weight
   NtopAtoms   = Origin_cell%Num_atoms
   NbottomAtoms= Origin_cell%Num_atoms
   allocate(TopAtoms(NtopAtoms))
   allocate(BottomAtoms(NbottomAtoms))
   do i=1, NTopAtoms
      TopAtoms(i)= i
   enddo
   do i=1, NBottomAtoms
      BottomAtoms(i)= i
   enddo

   rewind(1001)
   lfound = .false.
   do while (.true.)
      read(1001, *, end= 112, err=116, iostat=stat)inline
      inline=upper(inline)
      if (trim(adjustl(inline))=='SURFACE_ATOMS') then
         lfound= .true.
         if (cpuid==0) write(stdout, *)' '
         if (cpuid==0) write(stdout, *)'We found SURFACE_ATOMS card'
         exit
      endif
   enddo
   if (allocated(TopAtoms))deallocate(TopAtoms)
   allocate(TopAtoms(1))
   read(1001, '(A)', end= 112, err=116, iostat=stat)inline
   !> first count howmany values
   call param_get_range_vector('TOPATOMS',inline,NTopAtoms,.true., TOPATOMS)
   if (allocated(TopAtoms))deallocate(TopAtoms)
   allocate(TopAtoms(NTopAtoms))
   !> then get values
   call param_get_range_vector('TOPATOMS',inline,NTopAtoms,.false., TopAtoms)


   if (allocated(BottomAtoms))deallocate(BottomAtoms)
   allocate(BottomAtoms(1))
   read(1001, '(A)', end= 112, err=116, iostat=stat)inline
   !> first count howmany values
   call param_get_range_vector('BOTTOMATOMS',inline,NBottomAtoms,.true., BottomAtoms)
   if (allocated(BottomAtoms))deallocate(BottomAtoms)
   allocate(BottomAtoms(NBottomAtoms))
   !> then get values
   call param_get_range_vector('BOTTOMATOMS',inline,NBottomAtoms,.false., BottomAtoms)

   !> error happens when reading SURFACE_ATOMS
116 if (stat/=0 .and. cpuid==0) then
       write(stdout, '(a)')'>>> ERROR: There are something wrong with the SURFACE_ATOMS card' 
       write(stdout, '(a)')'    It should like this:'
       write(stdout, '(a)')'SURFACE_ATOMS '
       write(stdout, '(a)')"1 3 ! top surface's atom index '"
       write(stdout, '(a)')"2 4 ! bottom surface's atom index '"
    endif


112 continue
   if (.not.lfound.and.cpuid==0)write(stdout, *)'>> Output all atoms weight for surface state spectrum'
   if (cpuid==0) write(stdout, *)'> NtopAtoms ', NtopAtoms
   if (cpuid==0) write(stdout, '(a)')'> TopAtoms '
   if (cpuid==0) write(stdout, '(10i6)')TopAtoms
   if (cpuid==0) write(stdout, '(a, i10)')'> NbottomAtoms ', NbottomAtoms
   if (cpuid==0) write(stdout, '(a)')'> BottomAtoms '
   if (cpuid==0) write(stdout, '(10i6)')BottomAtoms

   NtopOrbitals=0
   do i=1, NTopAtoms
      NtopOrbitals= NtopOrbitals+ Origin_cell%nprojs(TopAtoms(i))
   enddo
   if (SOC>0) NtopOrbitals= NtopOrbitals*2
   allocate(TopOrbitals(NtopOrbitals))
   TopOrbitals= 1

   !> set up top surface orbitals for output the surface spectrum
   io=0
   do i=1, NTopAtoms
      do j=1, Origin_cell%nprojs(TopAtoms(i))
         io =io+ 1
         TopOrbitals(io)= orbitals_start(TopAtoms(i))+ j- 1
         if (SOC>0)TopOrbitals(io+ NtopOrbitals/2 )= orbitals_start(TopAtoms(i))+ j- 1+ NumberOfspinorbitals/2
      enddo ! j
   enddo ! i

   NBottomOrbitals=0
   do i=1, NBottomAtoms
      NBottomOrbitals= NBottomOrbitals+ Origin_cell%nprojs(BottomAtoms(i))
   enddo
   if (SOC>0) NBottomOrbitals= NBottomOrbitals*2
   allocate(BottomOrbitals(NBottomOrbitals))
   BottomOrbitals= 1

   !> set up Bottom surface orbitals for output the surface spectrum
   io=0
   do i=1, NBottomAtoms
      do j=1, Origin_cell%nprojs(BottomAtoms(i))
         io =io+ 1
         BottomOrbitals(io)= orbitals_start(BottomAtoms(i))+ j- 1
         if (SOC>0)BottomOrbitals(io+ NBottomOrbitals/2)= orbitals_start(BottomAtoms(i))+ j- 1+ NumberOfspinorbitals/2
      enddo ! j
   enddo ! i

   if (cpuid==0) write(stdout, *)'> NtopOrbitals ', NtopOrbitals
   if (cpuid==0) write(stdout, '(a)')'> TopOrbitals '
   if (cpuid==0) write(stdout, '(10i6)')TopOrbitals
   if (cpuid==0) write(stdout, '(a,999i4)')'> NBottomOrbitals ', NBottomOrbitals
   if (cpuid==0) write(stdout, '(a)')'> BottomOrbitals '
   if (cpuid==0) write(stdout, '(10i6)')BottomOrbitals

!===============================================================================================================!
!> NL_CHIRALITY card
!===============================================================================================================!


   !> setup for Weyl points chirality calculation
   !> default
   Num_NLs= 0 ! in unit of 1/Ang
   Rbig_NL= 0d0
   rsmall_a_NL= 0d0
   rsmall_b_NL= 0d0
   rewind(1001)
   lfound = .false.
   do while (.true.)
      read(1001, *, end= 211)inline
      inline=upper(inline)
      if (trim(adjustl(inline))=='NL_CHIRALITY') then
         lfound= .true.
         if (cpuid==0) write(stdout, *)' '
         if (cpuid==0) write(stdout, *)'We found NL_CHIRALITY card'
         exit
      endif
   enddo

   read(1001, *, end=219, err=219)Num_NLs
   allocate(NL_center_position_cart(3, Num_NLs))
   allocate(NL_center_position_direct(3, Num_NLs))
   NL_center_position_cart= 0d0
   NL_center_position_direct= 0d0
   read(1001, *, end=219, err=219)inline   ! The unit of lattice vector
   inline=upper(inline)
   DirectOrCart_NL= trim(adjustl(inline))
   read(1001, *, end=219, err=219)Rbig_NL, rsmall_a_NL, rsmall_b_NL
   it= 0
   do i=1, Num_NLs
      if (index(DirectOrCart_NL, "D")>0)then
         read(1001, *, end=219, err=219)NL_center_position_direct(:, i)
         call direct_cart_rec(NL_center_position_direct(:, i), NL_center_position_cart(:, i))
         it = it+ 1
      else
         read(1001, *, end=219, err=219)NL_center_position_cart(:, i)
         call cart_direct_rec(NL_center_position_cart(:, i), NL_center_position_direct(:, i))
         it = it+ 1
      endif
   enddo

219 continue

   !> print out the Weyl positions
   if (cpuid==0.and.lfound) then
      write(stdout, '(a)')" "
      write(stdout, '(a)')"Nodal Line center positions"
      write(stdout, '(8a10)')"kx", 'ky', 'kz', 'k1', 'k2', 'k3'
      do i=1, Num_NLs
         write(stdout, '(8f10.5)')NL_center_position_cart(:, i)*Angstrom2atomic, NL_center_position_direct(:, i)
      enddo

      write(stdout, '(a)')" "
   endif


   if (it< Num_NLs.and.cpuid==0) then
      write(stdout, *)' Error: Num_NLs should the same as ', &
         ' the number of weyl position lines'
      write(stdout, *)' Num_NLs = ', Num_NLs
      write(stdout, *)' Num of pos lines = ', it
      stop
   endif

211 continue
   if (cpuid==0) write(stdout, *)' '
   if (.not.lfound.and.cpuid==0)write(stdout, *)'>> We do not calculate chirality for nodal lines'
   if (.not.lfound.and.NLChirality_calc.and.cpuid==0) then
      write(stdout, *) 'ERROR: you should specify the NL_CHIRALITY card, see documentation'
   endif


   !> setup for Weyl points chirality calculation
   !> default
   Num_Weyls= 0 ! in unit of 1/Ang
   kr0=0
   rewind(1001)
   lfound = .false.
   do while (.true.)
      read(1001, *, end= 111)inline
      inline=upper(inline)
      if (trim(adjustl(inline))=='WEYL_CHIRALITY') then
         lfound= .true.
         if (cpuid==0) write(stdout, *)' '
         if (cpuid==0) write(stdout, *)'We found WEYL_CHIRALITY card'
         exit
      endif
   enddo

   read(1001, *, end=209, err=209)Num_Weyls
   allocate(weyl_position_direct(3, Num_Weyls))
   allocate(weyl_position_cart(3, Num_Weyls))
   weyl_position_direct= 0d0
   weyl_position_cart= 0d0
   read(1001, *, end=209, err=209)inline   ! The unit of lattice vector
   inline=upper(inline)
   DirectOrCart_Weyl= trim(adjustl(inline))
   read(1001, *, end=209, err=209)kr0
   it= 0
   do i=1, Num_Weyls
      if (index(DirectOrCart_Weyl, "D")>0)then
         read(1001, *, end=209, err=209)weyl_position_direct(:, i)
         call direct_cart_rec(weyl_position_direct(:, i), weyl_position_cart(:, i))
         it = it+ 1
      else
         read(1001, *, end=209, err=209)weyl_position_cart(:, i)
         call cart_direct_rec(weyl_position_cart(:, i), weyl_position_direct(:, i))
         it = it+ 1
      endif
   enddo

209 continue

   !> print out the Weyl positions
   if (cpuid==0.and.lfound) then
      write(stdout, '(a)')" "
      write(stdout, '(a)')"Weyl point positions"
      write(stdout, '(8a10)')"kx", 'ky', 'kz', 'k1', 'k2', 'k3'
      do i=1, Num_Weyls
         write(stdout, '(8f10.5)')weyl_position_cart(:, i)*Angstrom2atomic, weyl_position_direct(:, i)
      enddo

      write(stdout, '(a)')" "
   endif


   if (it< Num_Weyls.and.cpuid==0) then
      write(stdout, *)' Error: Num_Weyls should the same as ', &
         ' the number of weyl position lines'
      write(stdout, *)' Num_Weyls = ', Num_Weyls
      write(stdout, *)' Num of pos lines = ', it
      stop
   endif

111 continue
   if (cpuid==0) write(stdout, *)' '
   if (.not.lfound.and.cpuid==0)write(stdout, *)'>> We do not calculate chirality for weyl points'
   if (.not.lfound.and.WeylChirality_calc.and.cpuid==0) then
      write(stdout, *) 'ERROR: you should specify the WEYL_CHIRALITY card, see documentation'
   endif

!===============================================================================================================!
!> SELECTED_ATOMS card
!===============================================================================================================!

! SELECTED_ATOMS
! 2 ! NumberofSelectedAtoms_groups
! 1-3  8  ! atom indicies of group 1
! 4-6  9  ! atom indicies of group 2


   !>> parameters for SelectedAtoms
   !> this part is useful for surfstat or other slab or bulk band stOrigin_cell%Ructure calculations
   rewind(1001)
   lfound = .false.
   stat=0
   do while (.true.)
      read(1001, *, end= 332)inline
      inline=upper(inline)
      if (trim(adjustl(inline))=='SELECTED_ATOMS'.OR.trim(adjustl(inline))=='SELECTEDATOMS') then
         lfound= .true.
         if (cpuid==0) write(stdout, *)' '
         if (cpuid==0) write(stdout, *)'We found SELECTED_ATOMS card'
         exit
      endif
   enddo

   it= 0
   stat= 0
   NumberofSelectedAtoms_groups= 0
   read(1001, *, err=332, iostat=stat)NumberofSelectedAtoms_groups
   allocate(NumberofSelectedAtoms(NumberofSelectedAtoms_groups))
   allocate(Selected_Atoms(NumberofSelectedAtoms_groups))
   do i=1, NumberofSelectedAtoms_groups
      read(1001, '(A)', err=332, iostat=stat)inline
      idummy=1
      if (allocated(Selected_Atoms(i)%iarray))deallocate(Selected_Atoms(i)%iarray)
      allocate(Selected_Atoms(i)%iarray(1))

      !> first count howmany atoms for each line
      call param_get_range_vector('SelectedAtoms',inline,idummy,.true., Selected_Atoms(i)%iarray)
      NumberofSelectedAtoms(i)=idummy
      Selected_Atoms(i)%length=idummy
      if (idummy>0) then
         if (allocated(Selected_Atoms(i)%iarray))deallocate(Selected_Atoms(i)%iarray)
         allocate(Selected_Atoms(i)%iarray(idummy))
         !> then get values
         call param_get_range_vector('SelectedAtoms',inline,idummy,.false., Selected_Atoms(i)%iarray)
      else
         stop 'NumberofSelectedAtoms should be an integer and larger than zero'
      endif
   enddo

332 continue
   if (stat/=0) then
      if (cpuid==0) write(stdout, *)' '
      if (cpuid==0) write(stdout, *)' '
      if (cpuid==0) write(stdout, '(a)')'Error: Please set the right number of atoms in wt.in like this:'
      if (cpuid==0) write(stdout, *)' '
      if (cpuid==0) write(stdout, '(a)')'SELECTED_ATOMS'
      if (cpuid==0) write(stdout, '(a)')'2 ! number of groups'
      if (cpuid==0) write(stdout, '(a)')'1 2-4 ! atomic indices '
      if (cpuid==0) write(stdout, '(a)')'3 5 ! atomic indices '
      stop 'Errors happen in the WT.in, please check informations in the WT.out'
   endif

   !> setup SelectedAtoms if not specified by wt.in
   if (.not.allocated(Selected_Atoms)) then
      NumberofSelectedAtoms_groups=1
      allocate(NumberofSelectedAtoms(NumberofSelectedAtoms_groups))
      allocate(Selected_Atoms(NumberofSelectedAtoms_groups))
      allocate(Selected_Atoms(1)%iarray(Origin_cell%Num_atoms))
      NumberofSelectedAtoms(1)=Origin_cell%Num_atoms
      Selected_Atoms(1)%length=Origin_cell%Num_atoms
      do ia=1, Origin_cell%Num_atoms
         Selected_Atoms(1)%iarray(ia)= ia
      enddo
   endif

   if (cpuid==0) write(stdout, *)' '
   if (cpuid==0) write(stdout, '(a,i4,a)')'>> There are ', NumberofSelectedAtoms_groups, ' groups of SelectedAtoms'
   do i=1, NumberofSelectedAtoms_groups
      if (cpuid==0) write(stdout, '(a, i3)')'Group : ', i
      if (cpuid==0) write(stdout, '(a, 3i10)')'>> Number of atoms selected', &
         NumberofSelectedAtoms(i)
      if (cpuid==0) write(stdout, '(a)')'>> Selected atoms are'
      if (cpuid==0) write(stdout, '(10(i5, 2X, a))') &
         (Selected_Atoms(i)%iarray(ia), Origin_cell%atom_name(Selected_Atoms(i)%iarray(ia)), ia=1, NumberofSelectedAtoms(i))
   enddo

!===============================================================================================================!
!> SELECTEDWANNIERORBITALS card
!===============================================================================================================!

   !>> parameters for selectedorbitals
   !> this part is useful for surfstat or other slab or bulk band stOrigin_cell%Ructure calculations
   rewind(1001)
   lfound = .false.
   stat=0
   do while (.true.)
      read(1001, *, end= 331)inline
      inline=upper(inline)
      if (trim(adjustl(inline))=='SELECTED_WANNIERORBITALS'.or.trim(adjustl(inline))=='SELECTEDWANNIERORBITALS') then
         lfound= .true.
         if (cpuid==0) write(stdout, *)' '
         if (cpuid==0) write(stdout, *)'We found SELECTED_WANNIERORBITALS card'
         exit
      endif
   enddo

   it= 0
   stat= 0
   NumberofSelectedOrbitals_groups= 0
   read(1001, *, err=332, iostat=stat)NumberofSelectedOrbitals_groups
   allocate(Selected_WannierOrbitals(NumberofSelectedOrbitals_groups))
   allocate(NumberofSelectedOrbitals(NumberofSelectedOrbitals_groups))

   do i=1, NumberofSelectedOrbitals_groups
      read(1001, '(A)', err=332, iostat=stat)inline
      !> first count howmany orbitals selected
      idummy= 1
      if (allocated(Selected_WannierOrbitals(i)%iarray))deallocate(Selected_WannierOrbitals(i)%iarray)
      allocate(Selected_WannierOrbitals(i)%iarray(1))
      
      call param_get_range_vector('SelectedOrbitals',inline,idummy,.true., Selected_WannierOrbitals(i)%iarray)
      NumberofSelectedOrbitals(i)= idummy
      Selected_WannierOrbitals(i)%length= idummy

      if (NumberofSelectedOrbitals(i)>0) then
         if (allocated(Selected_WannierOrbitals(i)%iarray))deallocate(Selected_WannierOrbitals(i)%iarray)
         allocate(Selected_WannierOrbitals(i)%iarray(idummy))
         
         !> then get values
         call param_get_range_vector('SelectedOrbitals',inline,idummy,.false., Selected_WannierOrbitals(i)%iarray)
      else
         stop 'NumberofSelectedOrbitals should be an integer and larger than zero'
      endif
   enddo

331 continue
   if (stat/=0) then
      if (cpuid==0) write(stdout, *)' '
      if (cpuid==0) write(stdout, *)' '
      if (cpuid==0) write(stdout, '(a)')'Error: Please set the right number of orbitals and orbitals in wt.in like this:'
      if (cpuid==0) write(stdout, *)' '
      if (cpuid==0) write(stdout, '(a)')'SELECTED_WANNIERORBITALS'
      if (cpuid==0) write(stdout, '(a)')'2  ! number of groups of selectedorbitals'
      if (cpuid==0) write(stdout, '(a)')'1-3 ! orbitals indices '
      if (cpuid==0) write(stdout, '(a)')'4-7 ! orbitals indices '
      stop 'Errors happen in the WT.in, please check informations in the WT.out'
   endif


   !> setup SelectedOrbitals if not specified by wt.in
   !> by default we take the orbitals associated with Selected_Atoms
   if (.not.allocated(Selected_WannierOrbitals)) then
      NumberofSelectedOrbitals_groups= NumberofSelectedAtoms_groups
      allocate(NumberofSelectedOrbitals(NumberofSelectedAtoms_groups))
      allocate(Selected_WannierOrbitals(NumberofSelectedAtoms_groups))
      NumberofSelectedOrbitals= 0

      do ig=1, NumberofSelectedOrbitals_groups
         do i=1, NumberofSelectedAtoms(ig)
            ia = Selected_Atoms(ig)%iarray(i)
            NumberofSelectedOrbitals(ig)= NumberofSelectedOrbitals(ig)+ Origin_cell%nprojs(ia)
         enddo
         if (SOC>0) NumberofSelectedOrbitals(ig)= NumberofSelectedOrbitals(ig)*2
   
         allocate(Selected_WannierOrbitals(ig)%iarray(NumberofSelectedOrbitals(ig)))
         Selected_WannierOrbitals(ig)%iarray = 0
         io= 0
         do i=1, NumberofSelectedAtoms(ig)
            ia = Selected_Atoms(ig)%iarray(i)
            do j=1, Origin_cell%nprojs(ia)
               io = io+ 1
               Selected_WannierOrbitals(ig)%iarray(io)= index_start(ia)+ j- 1
            enddo
         enddo
         if (SOC>0) then
            do i=1, NumberofSelectedAtoms(ig)
               ia = Selected_Atoms(ig)%iarray(i)
               do j=1, Origin_cell%nprojs(ia)
                  io = io+ 1
                  Selected_WannierOrbitals(ig)%iarray(io)=index_start(ia)+ j+ NumberOfspinorbitals/2- 1
               enddo
            enddo
         endif
      enddo ! groups
   endif

   if (cpuid==0) write(stdout, *)' '
   if (cpuid==0) write(stdout, '(a,i3,a)')'>> There are ', NumberofSelectedOrbitals_groups, ' groups of SelectedOrbitals'
   do ig=1, NumberofSelectedOrbitals_groups
      if (cpuid==0) write(stdout, *)'>> SelectedOrbitals'
      if (cpuid==0) write(stdout, '(a, 3i10)')'>> Number of orbitals selected (including spin degenarcy)', &
         NumberofSelectedOrbitals(ig)
      if (cpuid==0) write(stdout, '(a)')'>> Orbitals are'
      if (cpuid==0) write(stdout, '(12i8)')Selected_WannierOrbitals(ig)%iarray(:)
   enddo

!===============================================================================================================!
!> SELECTED_OCCUPIEDBANDS card
!===============================================================================================================!

   !> parameters for selectedOccupiedBands
   rewind(1001)
   lfound = .false.
   stat=0
   do while (.true.)
      read(1001, *, end= 232)inline
      inline=upper(inline)
      if (trim(adjustl(inline))=='SELECTED_OCCUPIEDBANDS' .or.&
          trim(adjustl(inline))=='SELECTED_OCCUPIED_BANDS'.or.&
          trim(adjustl(inline))=='SELECTEDOCCUPIEDBANDS') then
         lfound= .true.
         if (cpuid==0) write(stdout, *)' '
         if (cpuid==0) write(stdout, *)'We found SELECTED_OCCUPIED_BANDS card'
         exit
      endif
   enddo

   it= 0
   stat= 0
   NumberofSelectedOccupiedBands= 0
   read(1001, '(A)', err=232, iostat=stat)inline
   !> get howmany integer numbers specified in the inline string
   call param_get_range_vector('SelectedOccupiedBands',inline,idummy,.true., Selected_Occupiedband_index)
   NumberofSelectedOccupiedBands= idummy

   if (NumberofSelectedOccupiedBands>0) then
      allocate(Selected_Occupiedband_index(NumberofSelectedOccupiedBands))
      Selected_Occupiedband_index= 0
      call param_get_range_vector('SelectedOccupiedBands',inline,idummy,.false., Selected_Occupiedband_index)
   else
      stop 'NumberofSelectedOccupiedBands should be an integer and larger than zero'
   endif

232 continue
   if (stat/=0) then
      if (cpuid==0) write(stdout, *)' '
      if (cpuid==0) write(stdout, *)' '
      if (cpuid==0) write(stdout, '(a)')'Error: Please set the right number of bands and band indices in wt.in like this:'
      if (cpuid==0) write(stdout, *)' '
      if (cpuid==0) write(stdout, '(a)')'SELECTED_OCCUPIED_BANDS'
      if (cpuid==0) write(stdout, '(a)')'4-7 ! band indices '
      stop 'Errors happen in the WT.in, please check informations in the WT.out'
   endif

   !> setup SELECTEDOccupiedBANDS
   !> if not given SELECTED_OCCUPIED_BANDS section, we will use NumOccupied as the inputs
   if (.not.allocated(Selected_Occupiedband_index))then
      NumberofSelectedOccupiedBands= NumOccupied
      allocate(Selected_Occupiedband_index(NumberofSelectedOccupiedBands))
      do i=1, NumberofSelectedOccupiedBands
         Selected_Occupiedband_index(i)= i
      enddo
   endif

   if (cpuid==0) write(stdout, *)' '
   if (cpuid==0) write(stdout, *)'>> SELECTED_OCCUPIED_BANDS'
   if (cpuid==0) write(stdout, '(a, 3i10)')'>> Number of Occupied bands selected ', &
      NumberofSelectedOccupiedBands
   if (cpuid==0) write(stdout, '(a)')'>> OccupiedBand indices are'
   if (cpuid==0) write(stdout, '(12i6)')Selected_Occupiedband_index(:)
   if (cpuid==0) write(stdout, *) ' '

!===============================================================================================================!
!> SELECTEDBANDS card
!===============================================================================================================!

   !> parameters for selectedbands
   rewind(1001)
   lfound = .false.
   stat=0
   do while (.true.)
      read(1001, *, end= 231)inline
      inline=upper(inline)
      if (trim(adjustl(inline))=='SELECTED_BANDS'.or.trim(adjustl(inline))=='SELECTEDBANDS') then
         lfound= .true.
         if (cpuid==0) write(stdout, *)' '
         if (cpuid==0) write(stdout, *)'We found SELECTED_BANDS card'
         exit
      endif
   enddo

   it= 0
   stat= 0
   NumberofSelectedBands= 0
   read(1001, *, err=231, iostat=stat)NumberofSelectedBands

   if (NumberofSelectedBands>0) then
      allocate(Selected_band_index(NumberofSelectedBands))
      Selected_band_index= 0
      read(1001, *, err=231, iostat=stat) (Selected_band_index(i), i=1, NumberofSelectedBands)
   else
      stop 'NumberofSelectedBands should be an integer and larger than zero'
   endif

231 continue
   if (stat/=0) then
      if (cpuid==0) write(stdout, *)' '
      if (cpuid==0) write(stdout, *)' '
      if (cpuid==0) write(stdout, '(a)')'Error: Please set the right number of bands and band indices in wt.in like this:'
      if (cpuid==0) write(stdout, *)' '
      if (cpuid==0) write(stdout, '(a)')'SELECTED_BANDS'
      if (cpuid==0) write(stdout, '(a)')'4  ! number of selected bands'
      if (cpuid==0) write(stdout, '(a)')'4 5 6 7 ! band indices '
      stop 'Errors happen in the WT.in, please check informations in the WT.out'
   endif

   !> setup SELECTEDBANDS
   if (.not.allocated(Selected_band_index))then
      NumberofSelectedBands= NumberOfspinorbitals
      allocate(Selected_band_index(NumberofSelectedBands))
      do i=1, NumberOfspinorbitals
         Selected_band_index(i)= i
      enddo
   endif

   if (cpuid==0) write(stdout, *)' '
   if (cpuid==0) write(stdout, *)'>> SELECTEDBANDS'
   if (cpuid==0) write(stdout, '(a, 3i10)')'>> Number of bands selected ', &
      NumberofSelectedBands
   if (cpuid==0) write(stdout, '(a)')'>> Band indices are'
   if (cpuid==0) write(stdout, '(12i6)')Selected_band_index(:)
   if (cpuid==0) write(stdout, *) ' '


!===============================================================================================================!
!> TBTOKP card
!===============================================================================================================!


   !> parameters for tbtokp
   Num_selectedbands_tbtokp = 0
   rewind(1001)
   lfound = .false.
   do while (.true.)
      read(1001, *, end= 220)inline
      inline=upper(inline)
      if (trim(adjustl(inline))=='TBTOKP') then
         lfound= .true.
         if (cpuid==0) write(stdout, *)' '
         if (cpuid==0) write(stdout, *)'We found TBTOKP card'
         exit
      endif
   enddo

   it= 0
   stat= 0
   Num_selectedbands_tbtokp= 0
   k_tbtokp= 0d0
   read(1001, *, err=220, iostat=stat)Num_selectedbands_tbtokp
   if (Num_selectedbands_tbtokp==0 .and. TBtoKP_calc) then
      stop 'Num_selectedbands_tbtokp should be an integer which is larger than zero if TBtoKP_calc=T'
   endif

   if (Num_selectedbands_tbtokp>0) then
      allocate(Selected_bands_tbtokp(Num_selectedbands_tbtokp))
      Selected_bands_tbtokp= 0
      read(1001, *, err=220, iostat=stat) (Selected_bands_tbtokp(i), i=1, Num_selectedbands_tbtokp)
   else
      stop 'Num_selectedbands_tbtokp should be larger than zero'
   endif
   read(1001, *, err=220, iostat=stat) k_tbtokp

   if (cpuid==0) write(stdout, *)' '
   if (.not.lfound.and.cpuid==0)write(stdout, *)'>>We donot constOrigin_cell%Ruct kp model.'
   if (cpuid==0) write(stdout, '(a, i10)')'>> Number of bands selected for kp model', &
      Num_selectedbands_tbtokp
   if (cpuid==0) write(stdout, '(a)')'>> Band indices are'
   if (cpuid==0) write(stdout, '(10i5)')Selected_bands_tbtokp(:)
   if (cpuid==0) write(stdout, '(a, 3f10.6)')'k point to constOrigin_cell%Ruct kp model in fractional coordinates', k_tbtokp

220 continue
   if (stat/=0) then
      if (cpuid==0) write(stdout, *)' '
      if (cpuid==0) write(stdout, *)' '
      if (cpuid==0) write(stdout, '(a)')'Error: Please set the right number of bands and band indices in wt.in like this:'
      if (cpuid==0) write(stdout, *)' '
      if (cpuid==0) write(stdout, '(a)')'TBTOKP'
      if (cpuid==0) write(stdout, '(a)')'8  ! number of selected bands to constOrigin_cell%Ruct kp model'
      if (cpuid==0) write(stdout, '(a)')'1 2 3 4 5 6 7 8 ! band indices '
      if (cpuid==0) write(stdout, '(a)')'0 0 0 ! k point in fractional coordinates '
      stop 'Errors happen in the WT.in, please check informations in the WT.out'
   endif

!===============================================================================================================!
!> ATOM_MASS card
!===============================================================================================================!

   !> for phonon system,  LO-TO correction, by T.T Zhang
   !> Atomic MASS in unit of g/mol
   rewind(1001)
   lfound = .false.
   do while (.true.)
      read(1001, *, end= 221)inline
      inline=upper(inline)
      if (trim(adjustl(inline))=='ATOM_MASS') then
         lfound= .true.
         if (cpuid==0) write(stdout, *)' '
         if (cpuid==0) write(stdout, *)'We found ATOM_MASS card'
         exit
      endif
   enddo
221 continue

   if (lfound) then
      read(1001,*)Origin_cell%Num_atom_type
      if(cpuid==0)write(stdout,'(a,i10,a)')'There are', Origin_cell%Num_atom_type, 'kind of atoms'
      if (.not.allocated(Origin_cell%Num_atoms_eachtype))allocate(Origin_cell%Num_atoms_eachtype(Origin_cell%Num_atom_type))
      allocate(mass_temp(Origin_cell%Num_atom_type))
      allocate(ATOM_MASS(Origin_cell%Num_atoms))
      read(1001,*)Origin_cell%Num_atoms_eachtype(1:Origin_cell%Num_atom_type)
      read(1001,*)mass_temp(1:Origin_cell%Num_atom_type)
      do i= 1, Origin_cell%Num_atom_type
         if (cpuid==0)write(stdout,'(a,i10,a)')'Each type have', Origin_cell%Num_atoms_eachtype(i), ' atoms'
         if (cpuid==0)write(stdout,'(a,f12.6)')'And their mass is', mass_temp(i)
      enddo
      it=0
      do i=1, Origin_cell%Num_atom_type
         do j=1, Origin_cell%Num_atoms_eachtype(i)
            it=it+1
            ATOM_MASS(it)=mass_temp(i)
         enddo
      enddo
   else
      if (LOTO_correction)stop "ERROR: please set ATOM_MASS card for LOTO correction of phonon spectrum"
   endif

   !>> setup Dielectric tensor for a given material
   rewind(1001)
   lfound = .false.
   do while (.true.)
      read(1001, *, end= 222)inline
      inline=upper(inline)
      if (trim(adjustl(inline))=='LOTO_DT') then
         lfound= .true.
         if (cpuid==0) write(stdout, *)' '
         if (cpuid==0) write(stdout, *)'We found LOTO_DT card for LOTO correction'
         exit
      endif
   enddo
222 continue

   if (lfound) then
      read(1001, *)Diele_Tensor(1,:)   ! Diele_Tensor is a 3*3 tensor for a material
      if (cpuid==0)write(stdout,'(a,3f12.5)')'Diele_tensor(1,:)',Diele_Tensor(1,:)
      read(1001, *)Diele_Tensor(2,:)
      if (cpuid==0)write(stdout,'(a,3f12.5)')'Diele_tensor(2,:)',Diele_Tensor(2,:)
      read(1001, *)Diele_Tensor(3,:)
      if (cpuid==0)write(stdout,'(a,3f12.5)')'Diele_tensor(3,:)',Diele_Tensor(3,:)
   else
      if (LOTO_correction) then
         if (cpuid==0) then
            write(stdout, *)"ERROR: please set LOTO_DT card for LOTO correction of phonon spectrum"
            write(stdout, *)"ERROR: please set Dielectronic Tensor information"
            stop "ERROR: Check error messages in WT.out"
         endif
      endif
   endif

   !>> setup Born charge for a given material
   rewind(1001)
   lfound = .false.
   do while (.true.)
      read(1001, *, end= 223)inline
      inline=upper(inline)
      if (trim(adjustl(inline))=='LOTO_BC') then
         lfound= .true.
         if (cpuid==0) write(stdout, *)' '
         if (cpuid==0) write(stdout, *)'We found LOTO_BC card for LOTO correction'
         exit
      endif
   enddo
223 continue

   if (lfound) then
      it=0
      allocate(Born_Charge(Origin_cell%Num_atoms,3,3))
      allocate(Born_Charge_temp(Origin_cell%Num_atom_type,3,3))
      do i=1,Origin_cell%Num_atom_type
         read(1001, *)Born_Charge_temp(i,1,:)
         read(1001, *)Born_Charge_temp(i,2,:)
         read(1001, *)Born_Charge_temp(i,3,:)
         do j=1,Origin_cell%Num_atoms_eachtype(i)
            it=it+1
            Born_Charge(it,:,:)=Born_Charge_temp(i,:,:)
            if (cpuid==0) then
               write(stdout,'(a,i3,2X,a6)')'Born_Charge for atom ', it, Origin_cell%atom_name(it)
               write(stdout,'(3f12.5)')Born_Charge(it,1,:)
               write(stdout,'(3f12.5)')Born_Charge(it,2,:)
               write(stdout,'(3f12.5)')Born_Charge(it,3,:)
            endif
         enddo
      enddo
   else
      if (LOTO_correction) then
         if (cpuid==0) then
            write(stdout, *)"ERROR: please set LOTO_BC card for LOTO correction of phonon spectrum"
            write(stdout, *)"ERROR: please set Born charge information"
            stop "ERROR: Check error messages in WT.out"
         endif
      endif
   endif


   !> close wt.in
   close(1001)

   eta=(omegamax- omegamin)/omeganum*2d0

   if(cpuid==0)write(stdout,*)'<<<Read wt.in file successfully'

   contains

   function upper(s1) result (s2)
      character(*)       :: s1
      character(len(s1)) :: s2
      character          :: ch
      integer, parameter :: DUC = ICHAR('A') - ICHAR('a')
      integer            :: i

      do i = 1,LEN(s1)
         ch = s1(i:i)
         if (ch >= 'a'.AND.ch <= 'z') ch = CHAR(ICHAR(ch)+DUC)
         s2(i:i) = ch
      enddo
   end function upper

end subroutine readinput


 !> rotate a vector in unit of  the original lattice vector into the new lattice
 !> vector defined by Umatrix
subroutine rotate_newlattice(R1, R2)
   use para, only : dp, Umatrix
   implicit none
   real(dp), intent(in) :: R1(3)
   real(dp), intent(inout) :: R2(3)
   real(dp), allocatable :: Umatrix_inv(:, :)

   allocate(Umatrix_inv(3, 3))
   Umatrix_inv= Umatrix

   call inv_r(3, Umatrix_inv)

   R2(1)= Umatrix_inv(1, 1)*R1(1)+ Umatrix_inv(2, 1)*R1(2)+ Umatrix_inv(3, 1)*R1(3)
   R2(2)= Umatrix_inv(1, 2)*R1(1)+ Umatrix_inv(2, 2)*R1(2)+ Umatrix_inv(3, 2)*R1(3)
   R2(3)= Umatrix_inv(1, 3)*R1(1)+ Umatrix_inv(2, 3)*R1(2)+ Umatrix_inv(3, 3)*R1(3)

   deallocate(Umatrix_inv)

   return
end subroutine rotate_newlattice


subroutine set_kcube3d
   use para
   implicit none

   integer :: knv3, knv3_mod
   integer  ::  ik, ik1, ik2, ik3

   !> Distribute k kpoints in the 3DCube into different MPI threads
   knv3= Nk1*Nk2*NK3
   KCube3D%Nk_total= knv3
   knv3_mod= mod(knv3, num_cpu)
   if (knv3_mod==0) then  !> perfect divided
      KCube3D%Nk_current= knv3/num_cpu
      KCube3D%Nk_start=1+ knv3*cpuid/num_cpu
      KCube3D%Nk_end  =(1+cpuid)*knv3/num_cpu
   else if (knv3/num_cpu==0) then    !> Number of MPI threads is large than knv3
      KCube3D%Nk_current= 1 !> one k piont per MPI thread
      KCube3D%Nk_start= cpuid+ 1 !> one k piont per MPI thread
      KCube3D%Nk_end  = cpuid+ 1
      if (cpuid+1 > knv3) then
         KCube3D%Nk_start= 1
         KCube3D%Nk_end  = 0
      endif
   else
      KCube3D%Nk_current= knv3/num_cpu+ 1
      if (cpuid< knv3_mod) then
         KCube3D%Nk_start= 1+ cpuid*KCube3D%Nk_current
         KCube3D%Nk_end  = (1+cpuid)*KCube3D%Nk_current
      else
         KCube3D%Nk_start= knv3_mod*KCube3D%Nk_current+ &
            (cpuid-knv3_mod)*(KCube3D%Nk_current-1)+1
         KCube3D%Nk_end  = knv3_mod*KCube3D%Nk_current+ &
            (cpuid-knv3_mod+1)*(KCube3D%Nk_current-1)
      endif
   endif

   !> calculate the volume of the k cube

   if (allocated(KCube3D%k_direct))deallocate(KCube3D%k_direct)
   allocate(KCube3D%k_direct(3, KCube3D%Nk_start:KCube3D%Nk_end))

   do ik= KCube3D%Nk_start, KCube3D%Nk_end
      ik1= (ik-1)/(Nk2*Nk3)+1
      ik2= ((ik-1-(ik1-1)*Nk2*Nk3)/Nk3)+1
      ik3= (ik-(ik2-1)*Nk3- (ik1-1)*Nk2*Nk3)
      KCube3D%k_direct(:, ik)= K3D_start_cube+ K3D_vec1_cube*(ik1-1)/dble(Nk1)  &
         + K3D_vec2_cube*(ik2-1)/dble(Nk2)  &
         + K3D_vec3_cube*(ik3-1)/dble(Nk3)
   enddo

end subroutine set_kcube3d

!> rotate a vector to the new coordinate system
!> rotate the vector from the original coordinate to the new coordinate
!> which defined like this: x is along R1', z is along R1'xR2', y is along z x y
!> Urot is a matrix linking the old coordinate and the new coordinate
subroutine rotate(R1, R2)
   use para, only : dp, Urot
   implicit none
   real(dp), intent(in) :: R1(3)
   real(dp), intent(inout) :: R2(3)

   R2(1)= Urot(1, 1)*R1(1)+ Urot(1, 2)*R1(2)+ Urot(1, 3)*R1(3)
   R2(2)= Urot(2, 1)*R1(1)+ Urot(2, 2)*R1(2)+ Urot(2, 3)*R1(3)
   R2(3)= Urot(3, 1)*R1(1)+ Urot(3, 2)*R1(2)+ Urot(3, 3)*R1(3)

   return
end subroutine rotate


!> transform from Cartesian coordinates to direct lattice vector basis for the magnetic cell
subroutine cart_direct_real_magneticcell(R1, R2)
   use para, only : dp, Magnetic_cell
   implicit none
   real(dp), intent(in) :: R1(3)
   real(dp), intent(inout) :: R2(3)
   real(dp), allocatable :: mata(:, :)

   allocate(mata(3, 3))

   mata(1, :)= Magnetic_cell%Rua
   mata(2, :)= Magnetic_cell%Rub
   mata(3, :)= Magnetic_cell%Ruc

   call inv_r(3, mata)
   R2= R1(1)*mata(1, :)+ R1(2)*mata(2, :)+ R1(3)*mata(3, :)

   deallocate(mata)

   return
end subroutine cart_direct_real_magneticcell



!> transform from Cartesian coordinates to direct lattice vector basis for the newcell defined by SURFACE card
subroutine cart_direct_real_newcell(R1, R2)
   use para
   implicit none
   real(dp), intent(in) :: R1(3)
   real(dp), intent(inout) :: R2(3)
   real(dp), allocatable :: mata(:, :)

   allocate(mata(3, 3))

   mata(1, :)= Rua_newcell
   mata(2, :)= Rub_newcell
   mata(3, :)= Ruc_newcell

   call inv_r(3, mata)
   R2= R1(1)*mata(1, :)+ R1(2)*mata(2, :)+ R1(3)*mata(3, :)

   deallocate(mata)

   return
end subroutine cart_direct_real_newcell


 !> transform from Cartesian coordinates to direct lattice vector basis
subroutine cart_direct_real_unfold(R1, R2)
   use para
   implicit none
   real(dp), intent(in) :: R1(3)
   real(dp), intent(inout) :: R2(3)
   real(dp), allocatable :: mata(:, :)

   allocate(mata(3, 3))

   mata(1, :)= Folded_cell%Rua
   mata(2, :)= Folded_cell%Rub
   mata(3, :)= Folded_cell%Ruc

   call inv_r(3, mata)
   R2= R1(1)*mata(1, :)+ R1(2)*mata(2, :)+ R1(3)*mata(3, :)

   deallocate(mata)

   return
end subroutine cart_direct_real_unfold



 !> transform from Cartesian coordinates to direct lattice vector basis
subroutine cart_direct_real(R1, R2, lattice)
   use para, only : dp
   implicit none
   real(dp), intent(in) :: R1(3)
   real(dp), intent(inout) :: R2(3)
   real(dp), intent(in) :: lattice(3, 3)
   real(dp), allocatable :: mata(:, :)

   allocate(mata(3, 3))

   mata= transpose(lattice)

   call inv_r(3, mata)
   R2= R1(1)*mata(1, :)+ R1(2)*mata(2, :)+ R1(3)*mata(3, :)

   deallocate(mata)

   return
end subroutine cart_direct_real

 !> transform from direct lattice vector basis to Cartesian coordinates for the newcell defined by SURFACE card
subroutine direct_cart_real_newcell(R1, R2)
   use para, only : dp, Rua_newcell, Rub_newcell, Ruc_newcell
   implicit none
   real(dp), intent(in) :: R1(3)
   real(dp), intent(inout) :: R2(3)

   R2= R1(1)*Rua_newcell+ R1(2)*Rub_newcell+ R1(3)*Ruc_newcell

   return
end subroutine direct_cart_real_newcell

!> transform from direct lattice vector basis to Cartesian coordinates for the magnetic cell
subroutine direct_cart_real_magneticcell(R1, R2)
   use para, only : dp, Magnetic_cell
   implicit none
   real(dp), intent(in) :: R1(3)
   real(dp), intent(inout) :: R2(3)

   R2= R1(1)*Magnetic_cell%Rua+ R1(2)*Magnetic_cell%Rub+ R1(3)*Magnetic_cell%Ruc

   return
end subroutine direct_cart_real_magneticcell



 !> transform from direct lattice vector basis to Cartesian coordinates
subroutine direct_cart_real_unfold(R1, R2)
   use para
   implicit none
   real(dp), intent(in) :: R1(3)
   real(dp), intent(inout) :: R2(3)

   R2= R1(1)*Folded_cell%Rua+ R1(2)*Folded_cell%Rub+ R1(3)*Folded_cell%Ruc

   return
end subroutine direct_cart_real_unfold


 !> transform from direct lattice vector basis to Cartesian coordinates
subroutine direct_cart_real(R1, R2, lattice)
   use para, only : dp
   implicit none
   real(dp), intent(in) :: R1(3)
   real(dp), intent(inout) :: R2(3)
   real(dp), intent(in) :: lattice(3, 3)

   R2= R1(1)*lattice(:, 1)+ R1(2)*lattice(:, 2)+ R1(3)*lattice(:, 3)

   return
end subroutine direct_cart_real

 !> transform from Cartesian coordinates to reciprocal lattice vector basis
subroutine cart_direct_rec_newcell(k1, k2)
   use para
   implicit none
   real(dp), intent(in) :: k1(3)
   real(dp), intent(inout) :: k2(3)

   real(dp), allocatable :: mata(:, :)

   allocate(mata(3, 3))

   mata(1, :)= Kua_newcell
   mata(2, :)= Kub_newcell
   mata(3, :)= Kuc_newcell

   call inv_r(3, mata)
   K2= k1(1)*mata(1, :)+ k1(2)*mata(2, :)+ k1(3)*mata(3, :)

   deallocate(mata)

   return
end subroutine cart_direct_rec_newcell


 !> transform from Cartesian coordinates to reciprocal lattice vector basis
subroutine cart_direct_rec_unfold(k1, k2)
   use para
   implicit none
   real(dp), intent(in) :: k1(3)
   real(dp), intent(inout) :: k2(3)

   real(dp), allocatable :: mata(:, :)

   allocate(mata(3, 3))

   mata(1, :)= Folded_cell%Kua
   mata(2, :)= Folded_cell%Kub
   mata(3, :)= Folded_cell%Kuc

   call inv_r(3, mata)
   K2= k1(1)*mata(1, :)+ k1(2)*mata(2, :)+ k1(3)*mata(3, :)

   deallocate(mata)

   return
end subroutine cart_direct_rec_unfold


 !> transform from Cartesian coordinates to reciprocal lattice vector basis for magnetic supercell
subroutine cart_direct_rec_magneticcell(k1, k2)
   use para, only : dp, Magnetic_cell
   implicit none
   real(dp), intent(in) :: k1(3)
   real(dp), intent(inout) :: k2(3)

   real(dp), allocatable :: mata(:, :)

   allocate(mata(3, 3))

   mata(1, :)= Magnetic_cell%Kua
   mata(2, :)= Magnetic_cell%Kub
   mata(3, :)= Magnetic_cell%Kuc

   call inv_r(3, mata)
   K2= k1(1)*mata(1, :)+ k1(2)*mata(2, :)+ k1(3)*mata(3, :)

   deallocate(mata)

   return
end subroutine cart_direct_rec_magneticcell


 !> transform from Cartesian coordinates to reciprocal lattice vector basis
subroutine cart_direct_rec(k1, k2)
   use para
   implicit none
   real(dp), intent(in) :: k1(3)
   real(dp), intent(inout) :: k2(3)

   real(dp), allocatable :: mata(:, :)

   allocate(mata(3, 3))

   mata(1, :)= Origin_cell%Kua
   mata(2, :)= Origin_cell%Kub
   mata(3, :)= Origin_cell%Kuc

   call inv_r(3, mata)
   k2= k1(1)*mata(1, :)+ k1(2)*mata(2, :)+ k1(3)*mata(3, :)

   deallocate(mata)

   return
end subroutine cart_direct_rec

subroutine direct_cart_rec_newcell(k1, k2)
   use para
   implicit none
   real(dp), intent(in) :: k1(3)
   real(dp), intent(inout) :: k2(3)

   K2= k1(1)*Kua_newcell+ k1(2)*Kub_newcell+ k1(3)*Kuc_newcell

   return
end subroutine direct_cart_rec_newcell

subroutine direct_cart_rec_magneticcell(k1, k2)
   use para, only : dp, Magnetic_cell
   implicit none
   real(dp), intent(in) :: k1(3)
   real(dp), intent(inout) :: k2(3)

   K2= k1(1)*Magnetic_cell%Kua+ k1(2)*Magnetic_cell%Kub+ k1(3)*Magnetic_cell%Kuc

   return
end subroutine direct_cart_rec_magneticcell


subroutine direct_cart_rec(k1, k2)
   use para
   implicit none
   real(dp), intent(in) :: k1(3)
   real(dp), intent(inout) :: k2(3)

   K2= k1(1)*Origin_cell%Kua+ k1(2)*Origin_cell%Kub+ k1(3)*Origin_cell%Kuc

   return
end subroutine direct_cart_rec

subroutine direct_cart_rec_unfold(k1, k2)
   use para
   implicit none
   real(dp), intent(in) :: k1(3)
   real(dp), intent(inout) :: k2(3)

   K2= k1(1)*Folded_cell%Kua+ k1(2)*Folded_cell%Kub+ k1(3)*Folded_cell%Kuc

   return
end subroutine direct_cart_rec_unfold


 !> define a new unit cell with the given MillerIndices [hkl]
subroutine MillerIndicestoumatrix()
   use para
   implicit none
   integer :: i1, i2, i3, h, k, l, it
   real(dp) :: R1(3), R2(3), R3(3), Rhkl(3), dot

   integer, allocatable :: vector_on_hkl_surface(:, :)
   integer :: Nvectors_on_hkl_surface
   real(dp) :: smallest_area, area
   real(dp) :: largestangle, angle
   real(dp) :: smallest_volume, cell_volume
   real(dp) :: smallest_length
   real(dp) :: norm_1, norm_2, norm_3

   integer :: iRmax
   iRmax= 6

   allocate(vector_on_hkl_surface(3, (2*iRmax+1)**3))
   vector_on_hkl_surface= 0

   h= MillerIndices(1)
   k= MillerIndices(2)
   l= MillerIndices(3)
   Rhkl= h*Origin_cell%Rua+ k*Origin_cell%Rub+ l*Origin_cell%Ruc

   !> Firstly, find all vectors that are orthorgonal to hkl
   it= 0
   do i1=-iRmax, iRmax
      do i2=-iRmax, iRmax
         do i3=-iRmax, iRmax
            if (i1==0 .and. i2==0 .and. i3==0) cycle
            !R= i1*Origin_cell%Rua+i2*Origin_cell%Rub+i3*Origin_cell%Ruc
            !dot=abs(R(1)*Rhkl(1)+ R(2)*Rhkl(2)+ R(3)*Rhkl(3))
            dot=abs(i1*h+i2*k+i3*l)
            if (dot<eps9) then
               it= it+1
               vector_on_hkl_surface(1, it)= i1
               vector_on_hkl_surface(2, it)= i2
               vector_on_hkl_surface(3, it)= i3
            endif
         enddo
      enddo
   enddo
   Nvectors_on_hkl_surface= it

   !> secondly, find the smallest area and the largest
   !> angle of two vectors
   smallest_area= 99999999d0
   largestangle= 0
   do i1=1, Nvectors_on_hkl_surface
      do i2=i1+1, Nvectors_on_hkl_surface
         R1= vector_on_hkl_surface(1, i1)*Origin_cell%Rua+ &
            vector_on_hkl_surface(2, i1)*Origin_cell%Rub+ &
            vector_on_hkl_surface(3, i1)*Origin_cell%Ruc
         R2= vector_on_hkl_surface(1, i2)*Origin_cell%Rua+ &
            vector_on_hkl_surface(2, i2)*Origin_cell%Rub+ &
            vector_on_hkl_surface(3, i2)*Origin_cell%Ruc
         dot= R1(1)*R2(1)+ R1(2)*R2(2)+ R1(3)*R2(3)
         norm_1= sqrt(R1(1)*R1(1)+ R1(2)*R1(2)+ R1(3)*R1(3))
         norm_2= sqrt(R2(1)*R2(1)+ R2(2)*R2(2)+ R2(3)*R2(3))
         if (norm_1*norm_2<eps9) stop 'norm of vector should larger than zero'
         R3(1)= R1(2)*R2(3)- R1(3)*R2(2)
         R3(2)= R1(3)*R2(1)- R1(1)*R2(3)
         R3(3)= R1(1)*R2(2)- R1(2)*R2(1)
         area= sqrt(R3(1)*R3(1)+ R3(2)*R3(2)+ R3(3)*R3(3))
         if (area<eps6) cycle
         if (area< smallest_area)smallest_area= area
      enddo
   enddo

   !> thirdly, find the largest
   !> angle in those two vectors which have smallest area
   largestangle= 0
   do i1=1, Nvectors_on_hkl_surface
      do i2=i1+1, Nvectors_on_hkl_surface
         R1= vector_on_hkl_surface(1, i1)*Origin_cell%Rua+ &
            vector_on_hkl_surface(2, i1)*Origin_cell%Rub+ &
            vector_on_hkl_surface(3, i1)*Origin_cell%Ruc
         R2= vector_on_hkl_surface(1, i2)*Origin_cell%Rua+ &
            vector_on_hkl_surface(2, i2)*Origin_cell%Rub+ &
            vector_on_hkl_surface(3, i2)*Origin_cell%Ruc

         R3(1)= R1(2)*R2(3)- R1(3)*R2(2)
         R3(2)= R1(3)*R2(1)- R1(1)*R2(3)
         R3(3)= R1(1)*R2(2)- R1(2)*R2(1)
         area= dsqrt(R3(1)*R3(1)+ R3(2)*R3(2)+ R3(3)*R3(3))

         dot= R1(1)*R2(1)+ R1(2)*R2(2)+ R1(3)*R2(3)
         norm_1= dsqrt(R1(1)*R1(1)+ R1(2)*R1(2)+ R1(3)*R1(3))
         norm_2= dsqrt(R2(1)*R2(1)+ R2(2)*R2(2)+ R2(3)*R2(3))
         angle= dacos(dot/norm_1/norm_2)
         if (angle>pi/2d0) angle= abs(angle-pi)

         if (dabs(area- smallest_area)<eps9)then
            if (angle> largestangle)largestangle= angle
         endif
      enddo
   enddo

   !> thirdly, find the two vectors which have smallest area and largest
   !> angle
   l1: do i1=1, Nvectors_on_hkl_surface
      do i2=i1+1, Nvectors_on_hkl_surface
         R1= vector_on_hkl_surface(1, i1)*Origin_cell%Rua+ &
            vector_on_hkl_surface(2, i1)*Origin_cell%Rub+ &
            vector_on_hkl_surface(3, i1)*Origin_cell%Ruc
         R2= vector_on_hkl_surface(1, i2)*Origin_cell%Rua+ &
            vector_on_hkl_surface(2, i2)*Origin_cell%Rub+ &
            vector_on_hkl_surface(3, i2)*Origin_cell%Ruc
         dot= R1(1)*R2(1)+ R1(2)*R2(2)+ R1(3)*R2(3)
         norm_1= dsqrt(R1(1)*R1(1)+ R1(2)*R1(2)+ R1(3)*R1(3))
         norm_2= dsqrt(R2(1)*R2(1)+ R2(2)*R2(2)+ R2(3)*R2(3))
         angle= dacos(dot/norm_1/norm_2)
         R3(1)= R1(2)*R2(3)- R1(3)*R2(2)
         R3(2)= R1(3)*R2(1)- R1(1)*R2(3)
         R3(3)= R1(1)*R2(2)- R1(2)*R2(1)
         area= dsqrt(R3(1)*R3(1)+ R3(2)*R3(2)+ R3(3)*R3(3))
         if (angle>pi/2d0) angle= abs(angle-pi)
         if (dabs(area- smallest_area)<eps9 .and. dabs(angle-largestangle)<eps9)then
            Umatrix(1, :)= vector_on_hkl_surface(:, i1)
            Umatrix(2, :)= vector_on_hkl_surface(:, i2)
            exit l1
         endif
      enddo
   enddo l1

   !> The last step, find the third vector that makes the new unit cell has
   !> the same volume as the old unit cell
   smallest_volume= 9999999d0
   R1= Umatrix(1, 1)*Origin_cell%Rua+  Umatrix(1, 2)*Origin_cell%Rub+  Umatrix(1, 3)*Origin_cell%Ruc
   R2= Umatrix(2, 1)*Origin_cell%Rua+  Umatrix(2, 2)*Origin_cell%Rub+  Umatrix(2, 3)*Origin_cell%Ruc
   do i1=-iRmax, iRmax
      do i2=-iRmax, iRmax
         do i3=-iRmax, iRmax
            if (i1==0 .and. i2==0 .and. i3==0) cycle
            R3= i1*Origin_cell%Rua+i2*Origin_cell%Rub+i3*Origin_cell%Ruc
            cell_volume= R1(1)*(R2(2)*R3(3)- R2(3)*R3(2)) &
               +R1(2)*(R2(3)*R3(1)- R2(1)*R3(3)) &
               +R1(3)*(R2(1)*R3(2)- R2(2)*R3(1))
            cell_volume= dabs(cell_volume)
            if (cell_volume< eps9) cycle
            if (cell_volume< smallest_volume) smallest_volume= cell_volume
         enddo
      enddo
   enddo

   !> find the third vector with the shortest length
   smallest_length= 9999999d0
   do i1=-iRmax, iRmax
      do i2=-iRmax, iRmax
         do i3=-iRmax, iRmax
            if (i1==0 .and. i2==0 .and. i3==0) cycle
            R3= i1*Origin_cell%Rua+i2*Origin_cell%Rub+i3*Origin_cell%Ruc
            cell_volume= R1(1)*(R2(2)*R3(3)- R2(3)*R3(2)) &
               +R1(2)*(R2(3)*R3(1)- R2(1)*R3(3)) &
               +R1(3)*(R2(1)*R3(2)- R2(2)*R3(1))
            cell_volume= dabs(cell_volume)
            if (dabs(cell_volume- smallest_volume)<eps9) then
               norm_3= dsqrt(R3(1)*R3(1)+ R3(2)*R3(2)+ R3(3)*R3(3))
               if (norm_3< smallest_length) smallest_length= norm_3
            endif
         enddo
      enddo
   enddo

   l2: do i1=-iRmax, iRmax
      do i2=-iRmax, iRmax
         do i3=-iRmax, iRmax
            if (i1==0 .and. i2==0 .and. i3==0) cycle
            R3= i1*Origin_cell%Rua+i2*Origin_cell%Rub+i3*Origin_cell%Ruc
            cell_volume= R1(1)*(R2(2)*R3(3)- R2(3)*R3(2)) &
               +R1(2)*(R2(3)*R3(1)- R2(1)*R3(3)) &
               +R1(3)*(R2(1)*R3(2)- R2(2)*R3(1))
            cell_volume= dabs(cell_volume)
            norm_3= dsqrt(R3(1)*R3(1)+ R3(2)*R3(2)+ R3(3)*R3(3))
            if (dabs(cell_volume- smallest_volume)<eps9 .and. &
               dabs(norm_3- smallest_length)<eps9) then
               Umatrix(3, 1)= i1
               Umatrix(3, 2)= i2
               Umatrix(3, 3)= i3
               exit l2
            endif
         enddo
      enddo
   enddo l2

   R1= Umatrix(1, 1)*Origin_cell%Rua+  Umatrix(1, 2)*Origin_cell%Rub+ Umatrix(1, 3)*Origin_cell%Ruc
   R2= Umatrix(2, 1)*Origin_cell%Rua+  Umatrix(2, 2)*Origin_cell%Rub+ Umatrix(2, 3)*Origin_cell%Ruc
   R3= Umatrix(3, 1)*Origin_cell%Rua+  Umatrix(3, 2)*Origin_cell%Rub+ Umatrix(3, 3)*Origin_cell%Ruc
   cell_volume= R1(1)*(R2(2)*R3(3)- R2(3)*R3(2)) &
      + R1(2)*(R2(3)*R3(1)- R2(1)*R3(3)) &
      + R1(3)*(R2(1)*R3(2)- R2(2)*R3(1))
   if (cell_volume<0) Umatrix(3, :)= -Umatrix(3, :)

   R1= Umatrix(1, 1)*Origin_cell%Rua+  Umatrix(1, 2)*Origin_cell%Rub+ Umatrix(1, 3)*Origin_cell%Ruc
   R2= Umatrix(2, 1)*Origin_cell%Rua+  Umatrix(2, 2)*Origin_cell%Rub+ Umatrix(2, 3)*Origin_cell%Ruc
   R3= Umatrix(3, 1)*Origin_cell%Rua+  Umatrix(3, 2)*Origin_cell%Rub+ Umatrix(3, 3)*Origin_cell%Ruc
   cell_volume= R1(1)*(R2(2)*R3(3)- R2(3)*R3(2)) &
      + R1(2)*(R2(3)*R3(1)- R2(1)*R3(3)) &
      + R1(3)*(R2(1)*R3(2)- R2(2)*R3(1))

   if (abs(cell_volume- Origin_cell%CellVolume)< eps9 .and. cpuid==0) then
      write(stdout, *)'  Congratulations, you got a unit cell that has ', &
         ' the same volume as the original unit cell '
      write(stdout, *)' The unitary rotation matrix is : '
      write(stdout, '(3f10.3)')Umatrix(1,:)
      write(stdout, '(3f10.3)')Umatrix(2,:)
      write(stdout, '(3f10.3)')Umatrix(3,:)

      write(stdout, *)' '
      write(stdout, *)'The lattice vectors for new cell are : '
      write(stdout, '(a,3f10.3)')' R1=', R1
      write(stdout, '(a,3f10.3)')' R2=', R2
      write(stdout, '(a,3f10.3)')' R3=', R3
      write(stdout, *)' Where R1 and R2 are in (hkl) plane'
   endif

   return
end subroutine MillerIndicestoumatrix

 !> define a new unit cell with the given two vectors of the SURFACE card
subroutine FindTheThirdLatticeVector()
   use para
   implicit none
   integer :: i1, i2, i3, h, k, l, it
   real(dp) :: R1(3), R2(3), R3(3), cross(3), dot

   integer, allocatable :: vectors_parallel_umatrix1(:, :)
   integer, allocatable :: vectors_parallel_umatrix2(:, :)
   integer :: Nvectors_parallel_umatrix1, Nvectors_parallel_umatrix2
   real(dp) :: smallest_volume, cell_volume
   real(dp) :: smallest_length
   real(dp) :: norm_3

   integer :: iRmax
   iRmax= 6

   allocate(vectors_parallel_umatrix1(3, (2*iRmax+1)**3))
   allocate(vectors_parallel_umatrix2(3, (2*iRmax+1)**3))
   vectors_parallel_umatrix1= 0
   vectors_parallel_umatrix2= 0

   !> Firstly, find all vectors that are parallel to Umatrix(:,1)
   it= 0
   do i1=-iRmax, iRmax
      do i2=-iRmax, iRmax
         do i3=-iRmax, iRmax
            if (i1==0 .and. i2==0 .and. i3==0) cycle
            !R= i1*Origin_cell%Rua+i2*Origin_cell%Rub+i3*Origin_cell%Ruc
            !dot=abs(R(1)*Rhkl(1)+ R(2)*Rhkl(2)+ R(3)*Rhkl(3))
            cross(1)= Umatrix(2, 1)*i3- i2*Umatrix(3, 1)
            cross(2)= Umatrix(3, 1)*i1- i3*Umatrix(1, 1)
            cross(3)= Umatrix(1, 1)*i2- i1*Umatrix(2, 1)
            dot = i1*Umatrix(1, 1)+ i2*Umatrix(2, 1)+ i3*Umatrix(3, 1)
            if ((abs(cross(1))+abs(cross(2))+abs(cross(3)))<eps9.and.dot>0) then
               it= it+1
               vectors_parallel_umatrix1(1, it)= i1
               vectors_parallel_umatrix1(2, it)= i2
               vectors_parallel_umatrix1(3, it)= i3
            endif
         enddo
      enddo
   enddo
   Nvectors_parallel_umatrix1= it

   !> and select the shortest vectors_parallel_umatrix1
   smallest_length= 9999999d0
   do it= 1, Nvectors_parallel_umatrix1
      R1= vectors_parallel_umatrix1(1, it)*Origin_cell%Rua+ &
         vectors_parallel_umatrix1(2, it)*Origin_cell%Rub+ &
         vectors_parallel_umatrix1(3, it)*Origin_cell%Ruc
      norm_3= dsqrt(R1(1)*R1(1)+ R1(2)*R1(2)+ R1(3)*R1(3))
      if (norm_3< smallest_length) then
         smallest_length= norm_3
         Umatrix(:, 1) = vectors_parallel_umatrix1(:, it)
      endif
   enddo


   !> secondly, find all vectors that are parallel to Umatrix(:,2)
   it= 0
   do i1=-iRmax, iRmax
      do i2=-iRmax, iRmax
         do i3=-iRmax, iRmax
            if (i1==0 .and. i2==0 .and. i3==0) cycle
            !R= i1*Origin_cell%Rua+i2*Origin_cell%Rub+i3*Origin_cell%Ruc
            !dot=abs(R(1)*Rhkl(1)+ R(2)*Rhkl(2)+ R(3)*Rhkl(3))
            cross(1)= Umatrix(2, 2)*i3- i2*Umatrix(3, 2)
            cross(2)= Umatrix(3, 2)*i1- i3*Umatrix(1, 2)
            cross(3)= Umatrix(1, 2)*i2- i1*Umatrix(2, 2)
            dot = i1*Umatrix(1, 2)+ i2*Umatrix(2, 2)+ i3*Umatrix(3, 2)
            if ((abs(cross(1))+abs(cross(2))+abs(cross(3)))<eps9.and.dot>0) then
               it= it+1
               vectors_parallel_umatrix2(1, it)= i1
               vectors_parallel_umatrix2(2, it)= i2
               vectors_parallel_umatrix2(3, it)= i3
            endif
         enddo
      enddo
   enddo
   Nvectors_parallel_umatrix2= it

   !> and select the shortest vectors_parallel_umatrix1
   smallest_length= 9999999d0
   do it= 1, Nvectors_parallel_umatrix2
      R2= vectors_parallel_umatrix2(1, it)*Origin_cell%Rua+ &
         vectors_parallel_umatrix2(2, it)*Origin_cell%Rub+ &
         vectors_parallel_umatrix2(3, it)*Origin_cell%Ruc
      norm_3= dsqrt(R2(1)*R2(1)+ R2(2)*R2(2)+ R2(3)*R2(3))
      if (norm_3< smallest_length) then
         smallest_length= norm_3
         Umatrix(:, 2) = vectors_parallel_umatrix2(:, it)
      endif
   enddo

   !> The last step, find the third vector that makes the new unit cell has
   !> the same volume as the old unit cell
   smallest_volume= 9999999d0
   R1= Umatrix(1, 1)*Origin_cell%Rua+  Umatrix(1, 2)*Origin_cell%Rub+  Umatrix(1, 3)*Origin_cell%Ruc
   R2= Umatrix(2, 1)*Origin_cell%Rua+  Umatrix(2, 2)*Origin_cell%Rub+  Umatrix(2, 3)*Origin_cell%Ruc
   do i1=-iRmax, iRmax
      do i2=-iRmax, iRmax
         do i3=-iRmax, iRmax
            if (i1==0 .and. i2==0 .and. i3==0) cycle
            R3= i1*Origin_cell%Rua+i2*Origin_cell%Rub+i3*Origin_cell%Ruc
            cell_volume= R1(1)*(R2(2)*R3(3)- R2(3)*R3(2)) &
               +R1(2)*(R2(3)*R3(1)- R2(1)*R3(3)) &
               +R1(3)*(R2(1)*R3(2)- R2(2)*R3(1))
            cell_volume= dabs(cell_volume)
            if (cell_volume< eps9) cycle
            if (cell_volume< smallest_volume) smallest_volume= cell_volume
         enddo
      enddo
   enddo

   !> find the third vector with the shortest length
   smallest_length= 9999999d0
   do i1=-iRmax, iRmax
      do i2=-iRmax, iRmax
         do i3=-iRmax, iRmax
            if (i1==0 .and. i2==0 .and. i3==0) cycle
            R3= i1*Origin_cell%Rua+i2*Origin_cell%Rub+i3*Origin_cell%Ruc
            cell_volume= R1(1)*(R2(2)*R3(3)- R2(3)*R3(2)) &
               +R1(2)*(R2(3)*R3(1)- R2(1)*R3(3)) &
               +R1(3)*(R2(1)*R3(2)- R2(2)*R3(1))
            cell_volume= dabs(cell_volume)
            if (dabs(cell_volume- smallest_volume)<eps9) then
               norm_3= dsqrt(R3(1)*R3(1)+ R3(2)*R3(2)+ R3(3)*R3(3))
               if (norm_3< smallest_length) smallest_length= norm_3
            endif
         enddo
      enddo
   enddo

   l2: do i1=-iRmax, iRmax
      do i2=-iRmax, iRmax
         do i3=-iRmax, iRmax
            if (i1==0 .and. i2==0 .and. i3==0) cycle
            R3= i1*Origin_cell%Rua+i2*Origin_cell%Rub+i3*Origin_cell%Ruc
            cell_volume= R1(1)*(R2(2)*R3(3)- R2(3)*R3(2)) &
               +R1(2)*(R2(3)*R3(1)- R2(1)*R3(3)) &
               +R1(3)*(R2(1)*R3(2)- R2(2)*R3(1))
            cell_volume= dabs(cell_volume)
            norm_3= dsqrt(R3(1)*R3(1)+ R3(2)*R3(2)+ R3(3)*R3(3))
            if (dabs(cell_volume- smallest_volume)<eps9 .and. &
               dabs(norm_3- smallest_length)<eps9) then
               Umatrix(3, 1)= i1
               Umatrix(3, 2)= i2
               Umatrix(3, 3)= i3
               exit l2
            endif
         enddo
      enddo
   enddo l2

   R1= Umatrix(1, 1)*Origin_cell%Rua+  Umatrix(1, 2)*Origin_cell%Rub+ Umatrix(1, 3)*Origin_cell%Ruc
   R2= Umatrix(2, 1)*Origin_cell%Rua+  Umatrix(2, 2)*Origin_cell%Rub+ Umatrix(2, 3)*Origin_cell%Ruc
   R3= Umatrix(3, 1)*Origin_cell%Rua+  Umatrix(3, 2)*Origin_cell%Rub+ Umatrix(3, 3)*Origin_cell%Ruc
   cell_volume= R1(1)*(R2(2)*R3(3)- R2(3)*R3(2)) &
      + R1(2)*(R2(3)*R3(1)- R2(1)*R3(3)) &
      + R1(3)*(R2(1)*R3(2)- R2(2)*R3(1))
   if (cell_volume<0) Umatrix(3, :)= -Umatrix(3, :)

   R1= Umatrix(1, 1)*Origin_cell%Rua+  Umatrix(1, 2)*Origin_cell%Rub+ Umatrix(1, 3)*Origin_cell%Ruc
   R2= Umatrix(2, 1)*Origin_cell%Rua+  Umatrix(2, 2)*Origin_cell%Rub+ Umatrix(2, 3)*Origin_cell%Ruc
   R3= Umatrix(3, 1)*Origin_cell%Rua+  Umatrix(3, 2)*Origin_cell%Rub+ Umatrix(3, 3)*Origin_cell%Ruc
   cell_volume= R1(1)*(R2(2)*R3(3)- R2(3)*R3(2)) &
      + R1(2)*(R2(3)*R3(1)- R2(1)*R3(3)) &
      + R1(3)*(R2(1)*R3(2)- R2(2)*R3(1))

   if (abs(cell_volume- Origin_cell%CellVolume)< eps9 .and. cpuid==0) then
      write(stdout, *)'  Congratulations, you got a unit cell that has ', &
         ' the same volume as the original unit cell '
      write(stdout, *)' The unitary rotation matrix is : '
      write(stdout, '(3f10.3)')Umatrix(1,:)
      write(stdout, '(3f10.3)')Umatrix(2,:)
      write(stdout, '(3f10.3)')Umatrix(3,:)

      write(stdout, *)' '
      write(stdout, *)'The lattice vectors for new cell are : '
      write(stdout, '(a,3f10.3)')' R1=', R1
      write(stdout, '(a,3f10.3)')' R2=', R2
      write(stdout, '(a,3f10.3)')' R3=', R3
      write(stdout, *)' Where R1, R2, R3 are in cartesian coordinates'
   else
      write(stdout, *) &
         " Warning:  I am sorry that I can't properly find unit cell with the first two vectors", &
         " defined in the SURFACE card that have the same volume as the original one." , &
         " Now, I will use my own method to choose the SURFACE card. But don't worry, ", &
         " The new surface card we found is just for the surface you defined. However, ", &
         " you should notice that the first and the second vectors in the SURFACE card ", &
         " could be changed which would affect the slab reciprocal lattice vectors."

   endif

   !> use MillerIndicestoumatrix
   if (abs(cell_volume- Origin_cell%CellVolume)> eps9 ) then
      !> first find the Miller indices
      h=int(Umatrix(1, 2)*Umatrix(2, 3)-Umatrix(1, 3)*Umatrix(2, 2))
      k=int(Umatrix(1, 3)*Umatrix(2, 1)-Umatrix(1, 1)*Umatrix(2, 3))
      l=int(Umatrix(1, 1)*Umatrix(2, 2)-Umatrix(1, 2)*Umatrix(2, 1))

      call gcd_reduce(h, k, l)
      MillerIndices(1)=h
      MillerIndices(2)=k
      MillerIndices(3)=l

      if (cpuid.eq.0) then
         write(stdout, '(a, 3i5)')'>> Miller indices for SURFACE are : ', h, k, l
      endif

      call MillerIndicestoumatrix()
   endif

   return
end subroutine FindTheThirdLatticeVector

!> a function to reduce h, k, l by their GCD (Greatest common divisor)
subroutine gcd_reduce(h, k, l)
   use para, only : dp
   implicit none

   integer, intent(inout) :: h, k, l
   integer :: i
   real(dp) :: rmiller(3), r3(3), sumr
   rmiller(1)= dble(h)
   rmiller(2)= dble(k)
   rmiller(3)= dble(l)

   do i=1, 3
      if (rmiller(i)<1E-3) cycle
      r3= rmiller/rmiller(i)
      sumr= abs(mod(sum(r3), 1d0))
      sumr=min(abs(sumr-1d0), abs(sumr))
      if (sumr<1E-3) then
         h= int(r3(1)) 
         k= int(r3(2)) 
         l= int(r3(3)) 
         return
      endif
   enddo

   return
end subroutine gcd_reduce


 !> move the atoms into the home unitcell [0, 1)*[0, 1)*[0, 1)
subroutine transformtohomecell(pos)
   ! Transform the k points to the 1st BZ
   !
   ! By QuanSheng Wu
   !
   ! wuquansheng@gmail.com
   !
   ! Nov 9 2016  at ETHZ

   use para, only : dp

   integer :: i
   real(dp), intent(inout) :: pos(3)

   do i=1, 3
      do while (.true.)
         if (pos(i)>= 0.0d0 .and. pos(i)<1.0d0) then
            exit
         else if (pos(i)< 0.0d0) then
            pos(i)= pos(i)+ 1d0
         else if (pos(i)>=1.0d0) then
            pos(i)= pos(i)- 1d0
         endif
      enddo
   enddo

   return
end subroutine transformtohomecell

!====================================================================!
subroutine param_get_range_vector(keyword,inline,length,lcount,i_value)
!====================================================================!
!!   Read a range vector eg. 1,2,3,4-10  or 1 3 400:100           
!!   if(lcount) we return the number of states in length            
!!   From Wannier90, modified by QSWU
!====================================================================!

!> usage
!> first count howmany values
!> call param_get_range_vector('TOPATOMS','1,2,3,4-10',NTopAtoms,lcount=.true., TOPATOMS)
!> then get values
!> call param_get_range_vector('TOPATOMS','1,2,3,4-10',NTopAtoms,lcount=.false., TopAtoms)

    implicit none

    character(len=*), intent(in) :: keyword
    character(len=*), intent(inout) :: inline
    integer,           intent(inout) :: length
    !! Number of states
    logical,           intent(in)    :: lcount
    !! If T only count states
    integer, intent(out)   :: i_value(length)
    !! States specified in range vector

    integer   :: loop,num1,num2,i_punc
    integer   :: counter,i_digit,loop_r,range_size
    character(len=256) :: dummy
    character(len=10), parameter :: c_digit="0123456789"
    character(len=2) , parameter :: c_range="-:"
    character(len=3) , parameter :: c_sep=" ,;"
    character(len=5) , parameter :: c_punc=" ,;-:"
    character(len=2) , parameter :: comment_punc="#!%"
    character(len=5)  :: c_num1,c_num2

    !> remove the comment part
    i_punc= scan(inline,comment_punc)
    if (i_punc>0) then
       dummy= inline(1:i_punc-1)
    else
       dummy=adjustl(inline)
    endif

    dummy=adjustl(dummy)
    
    counter=0
    do 
       i_punc=scan(dummy,c_punc)
       if(i_punc==0) call printerrormsg('Error parsing keyword '//trim(keyword)) 
       c_num1=dummy(1:i_punc-1)
       read(c_num1,*,err=1201,end=1201) num1
       dummy=adjustl(dummy(i_punc:))
       !look for range
       if(scan(dummy,c_range)==1) then
          i_digit=scan(dummy,c_digit)
          dummy=adjustl(dummy(i_digit:))
          i_punc=scan(dummy,c_punc)
          c_num2=dummy(1:i_punc-1)
          read(c_num2,*,err=1201,end=1201) num2
          dummy=adjustl(dummy(i_punc:))
          range_size=abs(num2-num1)+1
          do loop_r=1,range_size
             counter=counter+1
             if(.not. lcount) i_value(counter)=min(num1,num2)+loop_r-1
          end do
       else
          counter=counter+1 
          if(.not. lcount) i_value(counter)=num1
       end if

       if(scan(dummy,c_sep)==1) dummy=adjustl(dummy(2:))
       if(scan(dummy,c_range)==1) call printerrormsg('Error parsing keyword '//trim(keyword)//' incorrect range') 
       if(index(dummy,' ')==1) exit
    end do

    if(lcount) length=counter
    if(.not.lcount) then
       do loop=1,counter-1
          do loop_r=loop+1,counter 
             if(i_value(loop)==i_value(loop_r)) &
                call printerrormsg('Error parsing keyword '//trim(keyword)//' duplicate values')
          end do
        end do
    end if

    return

1201 call printerrormsg('Error parsing keyword '//trim(keyword))


end  subroutine param_get_range_vector


!> Write out the POSCAR for a given cell
subroutine writeout_poscar(cell, poscarname)
   use para
   implicit none

   integer :: ia
   type(cell_type) :: cell
   character(*) :: poscarname

   !> print out the new basis
   outfileindex= outfileindex+ 1
   if (cpuid.eq.0) then
      open(outfileindex, file=poscarname)
      write(outfileindex, '(a)')"POSCAR for generated by WannierTools"
      write(outfileindex, '(a)')"1.0"
      write(outfileindex, '(3f12.6)') cell%Rua/Angstrom2atomic
      write(outfileindex, '(3f12.6)') cell%Rub/Angstrom2atomic
      write(outfileindex, '(3f12.6)') cell%Ruc/Angstrom2atomic
      write(outfileindex, '(30A6)') cell%Name_of_atomtype
      write(outfileindex, '(30i6)') cell%Num_atoms_eachtype
      write(outfileindex, '(a)')"Direct"
      do ia=1, cell%Num_atoms
         if(cpuid==0)write(outfileindex, '(3f12.6, a9)')cell%Atom_position_direct(:, ia), trim(adjustl(cell%Atom_name(ia)))
      enddo
      close(outfileindex)
   endif
   return
end subroutine writeout_poscar


!> generate the POSCAR for slab system
!> necessary input:
!> Nslab
!> SURFACE
!> Vacuum_thickness_in_Angstrom
subroutine generate_slab_poscar(cell)
   use para
   implicit none

   type(cell_type) :: cell

   integer :: i, it, ia
   real(dp) :: angle_t, ratio
   real(dp) :: R1(3), R2(3), R3(3), R3_slab(3), R12_cross(3)
   integer, allocatable :: Num_atoms_eachtype(:)
   real(dp), allocatable :: pos_cart(:, :)
   character(10), allocatable :: atom_name(:)
   integer :: Num_atoms_slab, num_atoms_primitive_cell
   real(dp), external :: norm, angle

   R1=cell%Rua
   R2=cell%Rub
   R3=cell%Ruc

   !> R12_cross=R1xR2
   call cross_product(R1, R2, R12_cross)

   !> angle of R12_cross and R3
   angle_t= angle (R12_cross, R3)
   angle_t= angle_t*pi/180d0

   ratio= Vacuum_thickness_in_Angstrom/cos(angle_t)/norm(R3)

   R3_slab= (Nslab+ ratio)*R3


   num_atoms_primitive_cell= cell%Num_atoms
   Num_atoms_slab= cell%Num_atoms*Nslab
   allocate(atom_name(Num_atoms_slab))
   allocate(Num_atoms_eachtype(cell%Num_atom_type))
   Num_atoms_eachtype= cell%Num_atoms_eachtype*Nslab

   allocate(pos_cart(3, Num_atoms_slab))
   pos_cart=0d0

   it= 0
   do ia=1, num_atoms_primitive_cell
      do i=1, Nslab
         it=it+1
         pos_cart(:, it)= cell%Atom_position_cart(:, ia)+ R3*(i-1d0+ratio/2d0)
         atom_name(it)= cell%atom_name(ia)
      enddo
   enddo

   !> print out the new basis
   outfileindex= outfileindex+ 1
   if (cpuid.eq.0) then
      open(outfileindex, file="POSCAR-slab")
      write(outfileindex, '(a)')"POSCAR for slab defined by SURFACE card and Nslab and Vacuum_thickness_in_Angstrom in wt.in by WannierTools"
      write(outfileindex, '(a)')"1.0"
      write(outfileindex, '(3f12.6)') R1/Angstrom2atomic
      write(outfileindex, '(3f12.6)') R2/Angstrom2atomic
      write(outfileindex, '(3f12.6)') R3_slab/Angstrom2atomic
      write(outfileindex, '(30A6)') cell%Name_of_atomtype
      write(outfileindex, '(30i6)') Num_atoms_eachtype
      write(outfileindex, '(a)')"Cartesian"
      do ia=1, Num_atoms_slab
         if(cpuid==0)write(outfileindex, '(3f12.6, a9)')pos_cart(:, ia)/Angstrom2atomic, trim(adjustl(atom_name(ia)))
      enddo
      close(outfileindex)
   endif
   return
end subroutine generate_slab_poscar

!> generate a POSCAR for supercell defined by Nslab1, Nslab2, Nslab3
!> necessary input:
!> Nslab1, Nslab2, Nslab3
subroutine generate_supercell_poscar()
   use para
   implicit none

   integer :: i, it, ia
   real(dp) :: angle_t, ratio
   real(dp) :: R1(3), R2(3), R3(3), R3_slab(3), R12_cross(3)
   integer, allocatable :: Num_atoms_eachtype(:)
   real(dp), allocatable :: pos_cart(:, :)
   character(10), allocatable :: atom_name(:)
   integer :: Num_atoms_slab, num_atoms_primitive_cell
   real(dp), external :: norm, angle

   R1=Cell_defined_by_surface%Rua
   R2=Cell_defined_by_surface%Rub
   R3=Cell_defined_by_surface%Ruc

   !> R12_cross=R1xR2
   call cross_product(R1, R2, R12_cross)

   !> angle of R12_cross and R3
   angle_t= angle (R12_cross, R3)
   angle_t= angle_t*pi/180d0

   ratio= Vacuum_thickness_in_Angstrom/cos(angle_t)/norm(R3)

   R3_slab= (Nslab+ ratio)*R3


   num_atoms_primitive_cell= Cell_defined_by_surface%Num_atoms
   Num_atoms_slab= Cell_defined_by_surface%Num_atoms*Nslab
   allocate(atom_name(Num_atoms_slab))
   allocate(Num_atoms_eachtype(Cell_defined_by_surface%Num_atom_type))
   Num_atoms_eachtype= Cell_defined_by_surface%Num_atoms_eachtype*Nslab

   allocate(pos_cart(3, Num_atoms_slab))
   pos_cart=0d0

   it= 0
   do ia=1, num_atoms_primitive_cell
      do i=1, Nslab
         it=it+1
         pos_cart(:, it)= Cell_defined_by_surface%Atom_position_cart(:, ia)+ R3*(i-1d0+ratio/2d0)
         atom_name(it)= Cell_defined_by_surface%atom_name(ia)
      enddo
   enddo

   !> print out the new basis
   outfileindex= outfileindex+ 1
   if (cpuid.eq.0) then
      open(outfileindex, file="POSCAR-slab")
      write(outfileindex, '(a)')"POSCAR for slab defined by SURFACE card and Nslab and Vacuum_thickness_in_Angstrom in wt.in by WannierTools"
      write(outfileindex, '(a)')"1.0"
      write(outfileindex, '(3f12.6)') R1/Angstrom2atomic
      write(outfileindex, '(3f12.6)') R2/Angstrom2atomic
      write(outfileindex, '(3f12.6)') R3_slab/Angstrom2atomic
      write(outfileindex, '(30A6)') Cell_defined_by_surface%Name_of_atomtype
      write(outfileindex, '(30i6)') Num_atoms_eachtype
      write(outfileindex, '(a)')"Cartesian"
      do ia=1, Num_atoms_slab
         if(cpuid==0)write(outfileindex, '(3f12.6, a9)')pos_cart(:, ia)/Angstrom2atomic, trim(adjustl(atom_name(ia)))
      enddo
      close(outfileindex)
   endif
   return
end subroutine generate_supercell_poscar

subroutine get_reciprocal_lattice(R1, R2, R3, K1, K2, K3)
   !> Get reciprocal lattice vectors with given direct lattice vectors
   !> volume= R1.(R2xR3)
   !> K1=2*pi* R2xR3/volume
   !> K2=2*pi* R3xR1/volume
   !> K3=2*pi* R1xR2/volume
   use para, only : dp, pi
   implicit none
   real(dp) :: volume
   real(dp), intent(in)  :: R1(3), R2(3), R3(3)
   real(dp), intent(out) :: K1(3), K2(3), K3(3)
   
   call cross_product(R2, R3, K1)
   call cross_product(R3, R1, K2)
   call cross_product(R1, R2, K3)
   volume=dot_product(K1, R1)

   K1=2d0*pi*K1/volume
   K2=2d0*pi*K2/volume
   K3=2d0*pi*K3/volume
   return
end subroutine get_reciprocal_lattice

subroutine get_volume(R1, R2, R3, volume)
   !> Get volume with three given vectors
   !> volume= R1.(R2xR3)
   use para, only : dp
   implicit none
   real(dp) :: R0(3)
   real(dp), intent(in)  :: R1(3), R2(3), R3(3)
   real(dp), intent(out) :: volume
   
   call cross_product(R2, R3, R0)
   volume=dot_product(R0, R1)
   return
end subroutine get_volume

