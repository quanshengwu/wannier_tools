! some global parameters 

  module para

     use mpi
     implicit none

     ! output file name

     character*40 :: filename
     character*80 :: infilename

     logical :: BulkBand_calc
     logical :: SlabBand_calc
     logical :: WireBand_calc
     logical :: SlabSS_calc
     logical :: SlabArc_calc
     logical :: SlabSpintexture_calc
 
     ! double precision  
     integer,parameter :: Dp=kind(1.0d0)

     ! number of slabs of Bi2Se3 
     ! slab=1 means there is a quintuple layer of Bi2Se3 system
     integer :: Nslab

     !> number of princple layers for surface green's function
     integer :: Np

     integer, public, save :: ijmax=6

     !> leading dimension of surface green's function 
     integer :: Ndim

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

     ! the number of omega
     integer :: omeganum 

     ! omega interval 
     real(dp) :: omegamin, omegamax

     ! Fermi energy for arc calculation
     real(Dp) :: E_arc

     ! Fermi energy
     real(Dp) :: E_fermi


     ! circumference ratio pi  
     real(dp),parameter :: Pi= 3.14159265359d0

     real(Dp),parameter :: Ka(2)=(/1.0d0,0.0d0/)
     real(Dp),parameter :: Kb(2)=(/0.0d0,1.0d0/)

     real(Dp),parameter :: Ra2(2)=(/1d0,0.0d0/)
     real(Dp),parameter :: Rb2(2)=(/0.0d0,1d0/)

     real(Dp),parameter :: Ka2(2)=(/1d0,0.0d0/)
     real(Dp),parameter :: Kb2(2)=(/0.0d0,1d0/)

     ! three  primitive vectors  
     real(dp),public, save :: Rua(3)
     real(dp),public, save :: Rub(3)
     real(dp),public, save :: Ruc(3)

     ! three reciprocal primitive vectors  
     real(dp),public, save :: Kua(3)
     real(dp),public, save :: Kub(3)
     real(dp),public, save :: Kuc(3)

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


 end module para
