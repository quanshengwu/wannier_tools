!> some global parameters 
!> Copyright (c) 2010 QuanSheng Wu. All rights reserved.
!> add namelist for convenience  June 5th 2016 by QuanSheng Wu

  module wmpi
#if defined (MPI)
     include 'mpif.h'
#endif
  end module wmpi


  module para

     use wmpi
     implicit none

     integer,parameter :: stdout= 8

     !> define the file index to void the same index in different subroutines
     integer, public, save :: outfileindex= 10

     character(80) :: Hrfile
     character(80) :: Particle
     character(80) :: Package
     namelist / TB_FILE / Hrfile, Particle, Package

     
     !> control parameters
     logical :: BulkBand_calc    ! Flag for bulk energy band calculation
     logical :: BulkBand_points_calc    ! Flag for bulk energy band calculation
     logical :: BulkBand_plane_calc    ! Flag for bulk energy band calculation
     logical :: BulkFS_calc      ! Flag for bulk 3D fermi surface in 3D BZ calculation
     logical :: BulkFS_plane_calc ! Flag for bulk fermi surface for a fix k plane calculation
     logical :: BulkGap_cube_calc  ! Flag for Gap_cube calculation
     logical :: BulkGap_plane_calc ! Flag for Gap_plane calculation
     logical :: SlabBand_calc  ! Flag for 2D slab energy band calculation
     logical :: AHC_calc  ! Flag for Boltzmann tranport under magnetic field
     logical :: WireBand_calc  ! Flag for 1D wire energy band calculation
     logical :: SlabSS_calc    ! Flag for surface state ARPES spectrum calculation
     logical :: Dos_calc       ! Flag for density of state calculation
     logical :: JDos_calc      ! Flag for joint density of state calculation
     logical :: SlabArc_calc   ! Flag for surface state fermi-arc calculation
     logical :: SlabQPI_calc   ! Flag for surface state QPI spectrum calculation
     logical :: SlabSpintexture_calc ! Flag for surface state spin-texture calculation
     logical :: WannierCenter_calc  ! Flag for Wilson loop calculation
     logical :: Z2_3D_calc  ! Flag for Z2 number calculations of 6 planes
     logical :: WeylChirality_calc  ! Flag for Z2 number calculations of 6 planes
     logical :: Chern_3D_calc  ! Flag for Chern number calculations of 6 planes
     logical :: BerryPhase_calc   ! Flag for Berry phase calculation
     logical :: BerryCurvature_calc ! Flag for Berry curvature calculation
     logical :: EffectiveMass_calc  ! Flag for effective mass calculation
     logical :: FindNodes_calc  ! Flag for effective mass calculation
     

     namelist / Control / BulkBand_calc,BulkFS_calc, BulkGap_plane_calc, BulkFS_plane_calc, &
                          BulkGap_cube_calc, SlabBand_calc, WireBand_calc, &
                          SlabSS_calc, SlabArc_calc, SlabSpintexture_calc, &
                          WannierCenter_calc,BerryPhase_calc,BerryCurvature_calc, &
                          Z2_3D_calc, Chern_3D_calc, WeylChirality_calc, &
                          Dos_calc, JDos_calc, EffectiveMass_calc, SlabQPI_calc, FindNodes_calc, &
                          BulkBand_plane_calc, AHC_Calc, BulkBand_points_calc
 
     ! double precision  
     integer,parameter :: Dp=kind(1.0d0)

     ! number of slabs of Bi2Se3 
     ! slab=1 means there is a quintuple layer of Bi2Se3 system
     integer :: Nslab
     integer :: Nslab1
     integer :: Nslab2

     !> number of princple layers for surface green's function
     integer :: Np

     integer, public, save :: ijmax=6

     !> leading dimension of surface green's function 
     integer :: Ndim

     !> number of occupied bands for bulk unit cell
     integer :: Numoccupied

     !> number of electrons
     integer :: Ntotch

     integer :: Num_wann

     ! number of R points
     integer :: Nrpts

     ! number of k points used in ek_slab
     integer :: Nk  
     integer :: Nk1
     integer :: Nk2
     integer :: Nk3
     integer, parameter :: Nk2_max=10000  ! maximum number of k points 

     integer, public, save :: Nr1=5
     integer, public, save :: Nr2=5
     integer, public, save :: Nr3=2

     ! number of k points used for spintexture
     integer,parameter :: kmesh(2)=(/200 , 200/)
     integer,parameter :: knv=kmesh(1)*kmesh(2)


     ! a parameter to control soc
     ! Soc=0 means no spin-orbit coupling
     ! Soc!=0 means no spin-orbit coupling
     integer :: Soc

     ! used to calculate dos epsilon+i eta
     real(Dp) :: eta 
     real(Dp) :: Eta_Arc

     ! the number of omega
     integer :: OmegaNum 

     ! omega interval 
     real(dp) :: OmegaMin, OmegaMax

     ! Fermi energy for arc calculation
     real(Dp) :: E_arc

     ! threshold value for output the gap data for Gap3D
     real(Dp) :: Gap_threshold

     !> The largest distance between two WFs for which the Hamiltonian matrix element is retained and used in the band interpolation
     real(dp) :: bondlength_cutoff

     !> tolerance for wilson loop calculation in the integration over k
     real(dp) :: wcc_calc_tol


     !> namelist parameters
     namelist /PARAMETERS/ Eta_Arc, OmegaNum, OmegaMin, OmegaMax, &
        E_arc, Nk1, Nk2, Nk3, NP, Gap_threshold, wcc_calc_tol

     ! Fermi energy
     real(Dp) :: E_fermi

     !> surface onsite energy shift
     real(dp) :: surf_onsite

     !> magnetic field (Tesla)
     real(dp) :: Bx, By, Bz

     !> system parameters namelist
     namelist / SYSTEM / Soc, E_fermi, Bx, By, Bz, surf_onsite, &
        Nslab, Nslab1, Nslab2, Numoccupied, Ntotch, bondlength_cutoff

     !> e/2/h*a*a   a=1d-10m, h is the planck constant
     !> then the flux equals alpha*B*s
     real(dp),parameter :: alpha= 1.20736d0*1D-6

     !> some parameters related to atomic units
     real(dp),parameter :: bohr2atomic=0.529177211d0
     real(dp),parameter :: eV2Hartree= 1d0/27.211385d0

     ! circumference ratio pi  
     real(dp),parameter :: Pi= 3.14159265359d0
     real(dp),parameter :: half= 0.5d0
     real(dp),parameter :: zero= 0.0d0
     real(dp),parameter :: one = 1.0d0
     real(dp),parameter :: eps3= 1e-3
     real(dp),parameter :: eps6= 1e-6
     real(dp),parameter :: eps9= 1e-9

     real(Dp),parameter :: Ka(2)=(/1.0d0,0.0d0/)
     real(Dp),parameter :: Kb(2)=(/0.0d0,1.0d0/)

     real(Dp),public, save :: Ra2(2)
     real(Dp),public, save :: Rb2(2)

     real(Dp),public, save :: Ka2(2)
     real(Dp),public, save :: Kb2(2)

     ! three  primitive vectors in Cartsien coordinatec
     real(dp),public, save :: Rua(3)
     real(dp),public, save :: Rub(3)
     real(dp),public, save :: Ruc(3)

     !> three primitive vectors in new coordinate system, see slab part
     real(dp),public, save :: Rua_new(3)
     real(dp),public, save :: Rub_new(3)
     real(dp),public, save :: Ruc_new(3)

     ! three reciprocal primitive vectors  
     real(dp),public, save :: Kua(3)
     real(dp),public, save :: Kub(3)
     real(dp),public, save :: Kuc(3)

     real(dp),public, save :: Urot(3, 3)

     ! k list for 3D case band
     integer :: nk3lines
     integer :: nk3_band
     character(4), allocatable :: k3line_name(:)
     real(dp),allocatable :: k3line_stop(:)
     real(dp),allocatable :: k3line_start(:, :)
     real(dp),allocatable :: k3line_end(:, :)
     real(dp),allocatable :: K3list_band(:, :)
     real(dp),allocatable :: K3len(:)
     real(dp),allocatable :: K3points(:, :)

     !> k points in the point mode 
     integer :: Nk3_point_mode
     real(dp), allocatable :: k3points_pointmode_cart(:, :)
     real(dp), allocatable :: k3points_pointmode_direct(:, :)

     ! k path for berry phase, read from the input.dat
     ! in the KPATH_BERRY card
     integer :: NK_Berry
     character(10) :: DirectOrCart_Berry ! Whether direct coordinates or Cartisen coordinates
     real(dp), allocatable :: k3points_Berry(:, :) ! only in direct coordinates

     !>> top surface atoms
     integer :: NtopAtoms, NtopOrbitals
     integer, allocatable :: TopAtoms(:)
     integer, allocatable :: TopOrbitals(:)

     !>> bottom surface atoms
     integer :: NBottomAtoms, NBottomOrbitals
     integer, allocatable :: BottomAtoms(:)
     integer, allocatable :: BottomOrbitals(:)

     !>> effective mass

     !> k step for effective mass calculation
     real(dp), public, save :: dk_mass
     integer , public, save :: iband_mass
     real(dp), public, save :: k_mass(3)

     !>  klist for 2D case include all 2D system
     integer :: nk2lines
     integer :: knv2
     real(dp) :: kp(2, 32)
     real(dp) :: ke(2, 32)
     real(dp) :: k2line_stop(32)
     character(4) :: k2line_name(32)
     real(dp),allocatable :: k2len(:)
     real(dp),allocatable :: k2_path(:, :)

     !> A kpoint for 3D system--> only one k point
     real(dp), public, save :: Kpoint_3D_direct(3) ! the k point for effective mass calculation
     real(dp), public, save :: Kpoint_3D_cart(3) ! the k point for effective mass calculation
     
     character(10) :: DirectOrCart_SINGLE ! Whether direct coordinates or Cartisen coordinates
     real(dp), public, save :: Single_KPOINT_3D_DIRECT(3) ! the k point for effective mass calculation
     real(dp), public, save :: Single_KPOINT_3D_CART(3) ! the k point for effective mass calculation


     !> kpoints plane for 2D system--> arcs  
     real(dp) :: K2D_start(2)
     real(dp) :: K2D_vec1(2)
     real(dp) :: K2D_vec2(2)

     !> kpoints plane for 3D system --> gapshape
     real(dp) :: K3D_start(3)
     real(dp) :: K3D_vec1(3)
     real(dp) :: K3D_vec2(3)
     real(dp) :: K3D_vec3(3)

     !> kpoints plane for 3D system --> gapshape3D
     real(dp) :: K3D_start_cube(3)
     real(dp) :: K3D_vec1_cube(3)
     real(dp) :: K3D_vec2_cube(3)
     real(dp) :: K3D_vec3_cube(3)

     ! R coordinates  
     integer, allocatable     :: irvec(:,:)
     real(dp), allocatable    :: crvec(:,:)   ! R coordinates in Cartesian coordinates in units of Angstrom

     ! Hamiltonian m,n are band indexes
     complex(dp), allocatable :: HmnR(:,:,:)

     ! degree of degeneracy of R point 
     integer, allocatable     :: ndegen(:)
 
     ! complex constant 0+1*i
     complex(dp),parameter    :: zi=(0.0d0, 1.0d0)
     complex(dp),parameter    :: pi2zi=(0.0d0, 6.283185307179586d0)

     integer :: cpuid
     integer :: num_cpu

# if defined (MPI)
     integer, parameter :: mpi_in= mpi_integer
     integer, parameter :: mpi_dp= mpi_double_precision
     integer, parameter :: mpi_dc= mpi_double_complex
     integer, parameter :: mpi_cmw= mpi_comm_world
# endif 

     !> a matrix change old primitive cell to new primitive cell
     !> which can define new surface
     !> a 3*3 matrix
     real(dp), public, save :: Umatrix(3, 3)
     integer, public, save :: MillerIndices(3)    ! a matrix change old primitive cell to new primitive cell which can define new surface, it is a 3*3 matrix

     !> number of atoms in one primitive cell
     integer :: Num_atoms
     character(10) :: AngOrBohr
     character(10) :: DirectOrCart
     character(10), allocatable :: Atom_name(:)
     real(dp) :: CellVolume
     real(dp) :: kCubeVolume
     real(dp) :: PrimitiveCellVolume
     real(dp), allocatable :: Atom_position(:, :)
     real(dp), allocatable :: Atom_position_direct(:, :)
     real(dp), allocatable :: wannier_centers_cart(:, :)
     real(dp), allocatable :: wannier_centers_direct(:, :)

     integer :: max_projs
     integer, allocatable :: nprojs(:)
     character(10), allocatable :: proj_name(:, :)

     !> the start index for each atoms, only consider the spinless component
     integer, allocatable :: orbitals_start(:)

     !> symmetry operator apply on function basis
     complex(dp), allocatable :: inversion(:, :)
     complex(dp), allocatable :: mirror_x(:, :)
     complex(dp), allocatable :: mirror_z(:, :)
     complex(dp), allocatable :: glide(:, :)
     
     !> symmetry operator apply on coordinate system
     real(dp), allocatable :: inv_op(:, :)
     real(dp), allocatable :: mirror_z_op(:, :)
     real(dp), allocatable :: mirror_x_op(:, :)
     real(dp), allocatable :: mirror_y_op(:, :)
     real(dp), allocatable :: glide_y_op(:, :)

     !> weyl point information from the input.dat
     integer :: Num_Weyls
     character(10) :: DirectOrCart_Weyl ! Whether direct coordinates or Cartisen coordinates
     real(dp) :: kr0
     real(dp), allocatable :: weyl_position_cart(:, :)
     real(dp), allocatable :: weyl_position_direct(:, :)

     !> selected bands for magnetoresistance
     integer :: NumberofSelectedOrbitals
     integer, allocatable :: Selected_Orbitals(:)

     !> selected bands for magnetoresistance
     integer :: NumberofSelectedBands
     integer, allocatable :: Selected_band_index(:)

     !> time 
     character(8)  :: date_now
     character(10) :: time_now
     character(5)  :: zone_now

 end module para


 module wcc_module
    ! module for Wannier charge center (Wilson loop) calculations.
    use para
    implicit none

    
    type kline_wcc_type
       ! define a type of all properties at the k point to get the wcc
       real(dp) :: k(3)  ! coordinate
       real(dp) :: delta ! apart from the start point
       real(dp), allocatable :: wcc(:)  ! 1:Numoccupied, wannier charge center
       real(dp), allocatable :: gap(:)  ! 1:Numoccupied, the distance between  wcc neighbours
       real(dp) :: largestgap_pos_i     ! largest gap position
       real(dp) :: largestgap_pos_val   ! largest gap position value
       real(dp) :: largestgap_val       ! largest gap value
       logical  :: converged            ! converged or not 
       logical  :: calculated
    end type kline_wcc_type

    type kline_integrate_type
       ! define a type of all properties at the k point for integration
       real(dp) :: k(3)  ! coordinate
       real(dp) :: delta ! apart from the start point
       real(dp) :: b(3)  ! dis
       logical  :: calculated
       complex(dp), allocatable :: eig_vec(:, :)  !dim= (num_wann, num_wann)
    end type kline_integrate_type

 end module wcc_module


