  module prec
     !>> A module controls the precision. 
     !> when the nnzmax is larger than 2,147,483,647 then li=8,
     !> otherwise, li=4. 
     !> warning: li=4 was tested, li=8 is not tested yet.
     integer,parameter :: li=4 ! long integer
     integer,parameter :: Dp=kind(1.0d0) ! double precision  
  end module prec


  module para
     !> Some global parameters 

     use prec
     implicit none

     character(80) :: version

     integer,parameter :: stdout= 8

     real(dp), parameter :: pi=atan(1.0)*4d0
     real(dp), parameter :: eps8=1D-8
     real(dp), parameter :: eps6=1D-6
     real(dp), parameter :: eps4=1D-4
     real(dp), parameter :: eps3=1D-3
     real(dp), parameter :: eps2=1D-2
     real(dp), parameter :: eps1=1D-1
     real(dp), parameter :: Rcut = 10d0 ! 10 Angstrom cutoff for the hopping calculation


     integer :: outfileindex

     !> number of graphene layers
     integer :: number_layers

     !> twisted angle index
     !> (m, m+1) see reference 
     !> Fatemeh Haddadi, QuanSheng Wu,* Alex J. Kruchkov,* and Oleg V. Yazyev*
     !> Moiré Flat Bands in Twisted Double Bilayer Graphene
     !> doi:10.1021/acs.nanolett.9b05117
     !> theta = acos((3*m^2 + 3*m + 0.5)/(3*m^2 + 3*m + 1));
     integer :: twisted_index_m
     real(dp) :: twisted_angle_degree

     !> twisted angle arrangement
     !> this array is used to define the twisted angle. 
     !> the unit is theta which is defined by twisted_index_m
     !> in this program, we only focus on the system with only one twisted angle.
     integer :: twisted_angle_array_input(100)
     integer, allocatable :: twisted_angle_array(:)

     !> stacking sequences
     !> A-B-C
     character(1) :: stacking_sequences_input(100)
     character(1), allocatable :: stacking_sequences(:)

     !> use POSCAR or not
     logical :: use_poscar

     !> generate hmnr or not
     logical :: hr_generate

     !> generate sparse hamiltonian or not
     logical :: gen_sparse_hr

     !> cutoff for hopping integrals
     real(dp) :: hr_cutoff
     real(dp) :: vpppi

     !> in-plane lattice constant of Graphene
     real(dp) :: lattice_constant_a

     !> out of plane lattice constant of Graphene
     real(dp) :: lattice_constant_c


     !>  iR_cut
     integer :: iR_cut

     namelist / PARAMETERS / number_layers, twisted_index_m, &
        twisted_angle_array_input, stacking_sequences_input , &
        use_poscar, hr_generate, gen_sparse_hr, hr_cutoff, vpppi, &
        iR_cut, lattice_constant_a, lattice_constant_c

 end module para


 subroutine now(time_now)
    !> Convert the current wall-time into a real number
    !
    !> The unit is second.
    use para, only : dp

    implicit none
    integer   :: time_new(8)
    real(dp)  :: time_now

    call Date_and_time(values=time_new)
    time_now= time_new(3)*24*3600+time_new(5)*3600+&
              time_new(6)*60+time_new(7)+time_new(8)/1000d0
 
    return
 end subroutine now



