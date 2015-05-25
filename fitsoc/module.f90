! some global parameters 

  module para

     implicit none

     ! output file name
     character*40 :: filename
     character*80 :: infilename(2)

     ! double precision  
     integer,parameter :: Dp=kind(1.0d0)

     complex(dp), parameter :: zi= (0d0, 1d0)

	  integer :: Num_wann_nsoc
	  integer :: Num_wann_soc
	  integer :: num_bands_DFT

     ! number of R points
     integer :: Nrpts_nsoc
     integer :: Nrpts_soc

     !> soc= 1 : without soc in hr file 
     !> soc= 2 : consider soc in hr file
     integer :: soc

     ! R coordinates  
     integer, allocatable     :: irvec_soc(:,:)
     integer, allocatable     :: irvec_nsoc(:,:)

     ! Hamiltonian m,n are band indexes
     complex(dp), allocatable :: HmnR_soc(:,:,:)
     complex(dp), allocatable :: HmnR_nsoc(:,:,:)
     complex(dp), allocatable :: HmnR_nsoc_origin(:,:,:)

     ! degree of degeneracy of R point 
     integer, allocatable     :: ndegen_soc(:)
     integer, allocatable     :: ndegen_nsoc(:)

     integer :: Nk
     integer :: knv3

     !> k points coordinates
     real(dp), allocatable :: kpoints(:, :)

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

     !> weight for fitting of each kpoints and each band
     real(dp), allocatable :: weight(:, :)


     !> eigenvalue (Nwann, Nk)
     real(dp), allocatable :: eigval_soc(:, :)
     real(dp), allocatable :: eigval_soc_DFT(:, :)
     real(dp), allocatable :: eigval_nsoc(:, :)

     !> fermi level
     real(dp) :: E_fermi

     !> atoms information
     integer :: Num_atoms
     character(4), allocatable :: atom_name(:)
     real(dp), allocatable :: Atom_position(:, :)
     
     !> howmany projectors for each atom, with out spin degeneracy
     integer, allocatable :: nprojs(:)
     integer :: max_projs

     !> projectors name
     character(4), allocatable :: proj_name(:, :)

     !> spin orbital coupling strength
     real(dp), allocatable :: lambda_p(:)
     real(dp), allocatable :: lambda_d(:)

     real(dp), parameter :: pi=3.1415926536d0

 end module para
