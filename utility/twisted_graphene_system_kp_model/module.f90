  module prec
     !>> A module controls the precision. 
     !> when the nnzmax is larger than 2,147,483,647 then li=8,
     !> otherwise, li=4. 
     !> warning: li=4 was tested, li=8 is not tested yet.
     !> Author: Q.S Wu (wuquansheng@gmail.com)
     integer,parameter :: li=4 ! long integer
     integer,parameter :: Dp=kind(1.0d0) ! double precision  
  end module prec

  module wmpi
     use prec

#if defined (MPI)
     include 'mpif.h'
#endif

     integer :: cpuid  ! CPU id for mpi
     integer :: num_cpu  ! Number of processors for mpi

#if defined (MPI)
     integer, parameter :: mpi_in= mpi_integer
     integer, parameter :: mpi_dp= mpi_double_precision
     integer, parameter :: mpi_dc= mpi_double_complex
     integer, parameter :: mpi_cmw= mpi_comm_world
#endif 
  
     integer               :: BasisStart
     integer               :: BasisEnd
  
  end module wmpi

  module para
     !> Some global parameters 

     use wmpi
     use prec
     implicit none

     character(80) :: version

     integer,parameter :: stdout= 8

     real(dp), parameter :: pi=atan(1.0)*4d0
     real(dp), parameter :: sqrt3=dsqrt(3d0)
     real(dp), parameter :: eps8=1D-8
     real(dp), parameter :: eps6=1D-6
     real(dp), parameter :: eps4=1D-4
     real(dp), parameter :: eps3=1D-3
     real(dp), parameter :: eps2=1D-2
     real(dp), parameter :: eps1=1D-1
     real(dp), parameter :: Rcut = 10d0 ! 10 Angstrom cutoff for the hopping calculation

     complex(dp), parameter :: z0=(0d0, 0d0)
     complex(dp), parameter :: z1=(1d0, 0d0)
     complex(dp), parameter :: zi=(0d0, 1d0)
     complex(dp), parameter :: s0(2, 2)= reshape((/z1, z0, z0, z1/), (/2, 2/))

     !> three Pauli matrices
     complex(dp), parameter :: sx(2, 2)= reshape((/z0, z1, z1, z0/), (/2, 2/))
     complex(dp), parameter :: sy(2, 2)= reshape((/z0, zi,-zi, z0/), (/2, 2/))
     complex(dp), parameter :: sz(2, 2)= reshape((/z1, z0, z0, -z1/), (/2, 2/))


     integer :: outfileindex

     !> number of graphene layers
     integer :: number_layers

     !> twisted angle 
     integer :: twisted_index_m
     real(dp) :: twisted_angle_degree

     !> twisted angle arrangement
     !> this array is used to define the twisted angle. 
     !> the unit is theta which is defined by twisted_index_m
     !> in this program, we only focus on the system with only one twisted angle.
     integer :: twisted_angle_array_input(100)
     integer, allocatable :: twisted_angle_array(:)
     real(dp) :: interlayercoupling_ratio_array_input(100)
     real(dp), allocatable :: interlayercoupling_ratio_array(:)

     !> stacking sequences
     !> A-B-C
     character(1) :: stacking_sequences_input(100)
     character(1), allocatable :: stacking_sequences(:)

     !> coupling in the AA and AB range
     !> u_AA=0 gives the chiral limit E_n(k)=-E_{-n}(k)
     !> by default u_AA=79.7meV, u_AB=97.5meV from Phys. Rev. X 8, 031087
     !> vf=2.1354 eV* a
     real(dp) :: u_AA  ! hopping in AA region of TBG layer
     real(dp) :: u_AB  ! hopping in AB region of TBG layer
     !> vf=sqrt(3)/2*2.46*gamma_0
     real(dp) :: gamma_0, vf    ! Fermi velocity
     real(dp) :: vppsigma  ! nearest neighbour interlayer hopping
     !> v3=gamma_3*a*sqrt3/2
     real(dp) :: gamma_3, v3  ! next nearest neighbour interlayer hopping between two sublattices 
     !> v4=gamma_4*a*sqrt3/2
     real(dp) :: gamma_4, v4   ! next nearest neighbour interlayer hopping between the same sublattices
     real(dp) :: Electric_field  ! eV/Angstrom
     real(dp) :: lattice_constant_c  

     !> wave vector cutoff in the Hmnk construction
     !> see references
     integer :: Qcutoff

     !> number of k points per line
     integer :: Nk

     !> number of bands to be calculated
     integer :: Num_bands

     namelist / PARAMETERS / number_layers, twisted_angle_degree, &
        twisted_angle_array_input, stacking_sequences_input , &
        u_AA, u_AB, Qcutoff, vppsigma, Nk, gamma_0, vf, twisted_index_m, Num_bands, &
        gamma_3, gamma_4, &
        Electric_field, lattice_constant_c, interlayercoupling_ratio_array_input

     !> Graphene real lattice
     real(dp), parameter :: lattice_constant_graphene= 2.46d0  ! Angstrom
     real(dp), parameter :: a1(2)= (/0.5d0, sqrt3/2d0/)* lattice_constant_graphene
     real(dp), parameter :: a2(2)= (/-0.5d0, sqrt3/2d0/)* lattice_constant_graphene

     real(dp), parameter :: K_valley(2) = 4d0*pi/3d0/lattice_constant_graphene* (/1d0, 0d0/)
     real(dp), parameter :: Kp_valley(2) = 4d0*pi/3d0/lattice_constant_graphene* (/-1d0, 0d0/)

     !> define two reciprocal lattice vectors of morie cell, and qb, qtl, qtr
     real(dp) :: b1m(2), b2m(2), qb(2), qtl(2), qtr(2)

     !> define morie K valley, one is from top layer, the other one from bottom layer
     real(dp) :: K_valley1(2), K_valley2(2)

     !> rotation matrix defined by twisted angle
     real(dp) :: rot1(2, 2), rot2(2, 2)

     !> number of plane waves defined by Qcutoff
     integer :: Num_Qvectors

     !> dimension of the matrix
     !> Ndim= 2*num_Q*number_layers
     integer :: Ndim

     !> Q vectors 
     integer, allocatable :: Qvectors(:, :)

     ! k list for 3D case band
     integer :: nklines  ! Howmany k lines for bulk band calculation
     integer :: nk_band  ! Howmany k points for each k line
     character(4), allocatable :: kline_name(:) ! Name of the K points
     real(dp),allocatable :: kline_stop(:)  ! Connet points
     real(dp),allocatable :: kline_start(:, :) ! Start point for each k line
     real(dp),allocatable :: kline_end(:, :) ! End point for each k line
     real(dp),allocatable :: klist_band(:, :) ! coordinate of k points for bulk band calculation in kpath mode
     real(dp),allocatable :: klen(:)  ! put all k points in a line in order to plot the bands 
     real(dp),allocatable :: kpoints(:, :) ! coordinate of k points for bulk band calculation in cube mode

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



