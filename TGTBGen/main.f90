program main
   ! This program is used to generate Hamiltonian for Twisted Graphene system.
   ! constructed by Q.S.Wu On June 24, 2020
   ! wuquansheng@gmail.com
   ! Copyright (c) 2020 QuanSheng Wu. All rights reserved.

   use para
   implicit none

   !> file existence
   logical  :: exists
   integer :: ierr
   character(8) :: cht

   !> time measure
   real(dp) :: time_start, time_end, time_init

   ierr = 0
   !> initial the environment of mpi

   open(unit=stdout, file='hrgen.out')

   !call now(time_init)
   !call header

   !> print information for mpi
   !> readin the necessary parameters specified by user
   !> including the twisted angle, stacking pattern
   call readinput

   !> Generate the primitive unit cell 
   if (.not.use_poscar) call generate_crystal_structure()

   !> generate tight-binding Hamiltonian
   if (hr_generate) call generate_hr()
   
   !call now(time_end)

   !call print_time_cost(time_init, time_end, 'whole program')
   write(stdout,*)'Congratulations! you finished the calculation.'
   write(stdout,*)'See you next time :)'

end  !<< end of main program
