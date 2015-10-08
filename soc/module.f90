! some global parameters 

  module para

     use mpi
     implicit none

     ! output file name

     character*40 :: outfilename
     character*80 :: infilename

     logical :: BulkBand_calc
     logical :: SlabBand_calc
     logical :: WireBand_calc
     logical :: SlabSS_calc
     logical :: SlabArc_calc
     logical :: SlabSpintexture_calc
     logical :: wanniercenter_calc
     logical :: berry_calc

     integer,parameter :: stdout= 6

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
     real(Dp) :: eta_arc

     ! the number of omega
     integer :: omeganum 

     ! omega interval 
     real(dp) :: omegamin, omegamax

     ! Fermi energy for arc calculation
     real(Dp) :: E_arc

     ! Fermi energy
     real(Dp) :: E_fermi

     !> surface onsite energy shift
     real(dp) :: surf_onsite

     !> magnetic field (Tesla)
     real(dp) :: Bx, By, Bz

     !> e/2/h*a*a   a=1d-10m, h is the planck constant
     !> then the flux equals alpha*B*s
     real(dp),parameter :: alpha= 1.20736d0*1D-6

     ! circumference ratio pi  
     real(dp),parameter :: Pi= 3.14159265359d0
     real(dp),parameter :: half= 0.5d0
     real(dp),parameter :: zero= 0.0d0
     real(dp),parameter :: one = 1.0d0
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

     ! k list for band
     integer :: nk3lines
     integer :: nk3_band
     character(4), allocatable :: k3line_name(:)
     real(dp),allocatable :: k3line_stop(:)
     real(dp),allocatable :: k3line_start(:, :)
     real(dp),allocatable :: k3line_end(:, :)
     real(dp),allocatable :: K3list_band(:, :)
     real(dp),allocatable :: K3len(:)
     real(dp),allocatable :: K3points(:, :)

     ! R coordinates  
     integer, allocatable     :: irvec(:,:)

     ! Hamiltonian m,n are band indexes
     complex(dp), allocatable :: HmnR(:,:,:)

     ! degree of degeneracy of R point 
     integer, allocatable     :: ndegen(:)
 
     ! complex constant 0+1*i
     complex(dp),parameter    :: zi=(0.0d0, 1.0d0)
     complex(dp),parameter    :: pi2zi=(0.0d0, 6.283185307179586d0)

     integer :: cpuid
     integer :: num_cpu
     integer, parameter :: mpi_in= mpi_integer
     integer, parameter :: mpi_dp= mpi_double_precision
     integer, parameter :: mpi_dc= mpi_double_complex
     integer, parameter :: mpi_cmw= mpi_comm_world

     !> a matrix change old primitive cell to new primitive cell
     !> which can define new surface
     !> a 3*3 matrix
     real(dp), public, save :: Umatrix(3, 3)

     !> number of atoms in one primitive cell
     integer :: Num_atoms
     character(10), allocatable :: atom_name(:)
     real(dp), allocatable :: Atom_position(:, :)

     integer :: max_projs
     integer, allocatable :: nprojs(:)
     character(10), allocatable :: proj_name(:, :)


 end module para
