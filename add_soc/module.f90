! some global parameters 

  module para

     implicit none

     ! output file name
     character*40 :: filename
     character*80 :: infilename

     ! double precision  
     integer,parameter :: Dp=kind(1.0d0)

     complex(dp), parameter :: zi= (0d0, 1d0)

	  integer :: Num_wann

     ! number of R points
     integer :: Nrpts

     !> soc= 1 : without soc in hr file 
     !> soc= 2 : consider soc in hr file
     integer :: soc

     ! R coordinates  
     integer, allocatable     :: irvec(:,:)
     integer, allocatable     :: irvec_2d(:,:)

     ! Hamiltonian m,n are band indexes
     complex(dp), allocatable :: HmnR(:,:,:)

     ! degree of degeneracy of R point 
     integer, allocatable     :: ndegen(:)

     !> fermi level
     real(dp) :: E_fermi

     !> atoms information
     integer :: Num_atoms
     character(4), allocatable :: atom_name(:)
     real(dp), allocatable :: Atom_position(:, :)
     
     integer :: Num_atom_type
     integer, allocatable :: atom_type(:)

     !> howmany projectors for each atom, with out spin degeneracy
     integer, allocatable :: nprojs(:)
     integer :: max_projs

     !> projectors name
     character(8), allocatable :: proj_name(:, :)

     !> spin orbital coupling strength
     !> for each atom's type
     real(dp), allocatable :: lambda_p(:)
     real(dp), allocatable :: lambda_d(:)

     !> for every atom's 
     real(dp), allocatable :: magmomx_p(:)
     real(dp), allocatable :: magmomx_d(:)
     real(dp), allocatable :: magmomy_p(:)
     real(dp), allocatable :: magmomy_d(:)
     real(dp), allocatable :: magmomz_p(:)
     real(dp), allocatable :: magmomz_d(:)

 end module para
