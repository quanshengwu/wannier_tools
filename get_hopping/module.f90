! some global parameters 

  module para

     implicit none

     ! output file name
     character*40 :: filename
     character*80 :: infilename

     ! double precision  
     integer,parameter :: Dp=kind(1.0d0)

	  integer :: Num_wann

     ! number of R points
     integer :: Nrpts

     !> soc= 1 : without soc in hr file 
     !> soc= 2 : consider soc in hr file
     integer :: soc

     ! R coordinates  
     integer, allocatable     :: irvec(:,:)

     ! Hamiltonian m,n are band indexes
     complex(dp), allocatable :: HmnR(:,:,:)

     ! degree of degeneracy of R point 
     integer, allocatable     :: ndegen(:)

     !> fermi level
     real(dp) :: E_fermi

     !> lattice vector
     real(dp),public, save :: Rua(3)
     real(dp),public, save :: Rub(3)
     real(dp),public, save :: Ruc(3)
     
     !> atoms information
     integer :: Num_atoms
     character(4), allocatable :: atom_name(:)
     real(dp), allocatable :: Atom_position(:, :)
     
     !> howmany projectors for each atom, with out spin degeneracy
     integer, allocatable :: nprojs(:)
     integer :: max_projs

     !> projectors name
     character(4), allocatable :: proj_name(:, :)

 end module para
