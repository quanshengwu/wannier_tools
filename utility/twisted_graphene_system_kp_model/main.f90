program main
   ! This program is used to study Twisted Graphene system using continuum model
   ! Constructed by Q.S.Wu On Oct. 28, 2020
   ! wuquansheng@gmail.com
   ! Copyright (c) 2020 QuanSheng Wu. All rights reserved.

   !> This code is used to generate results for paper
   !> ShengNan Zhang, Bo Xie, QuanSheng Wu, Jianpeng Liu, Oleg V. Yazyev, 
   !> Chiral Decomposition of Twisted Graphene Multilayers with Arbitrary Stacking
   !> arXiv:2012.11964 (2020)

   !> References:
   !> J. M. B. Lopes dos Santos, N. M. R. Peres, and A. H. Castro Neto, Phys. Rev. Lett. 99, 256802 (2007)
   !> Rafi Bistritzer and Allan H. MacDonald, PNAS July 26, 108 (30) 12233-12237 (2011)
   !> J. M. B. Lopes dos Santos, N. M. R. Peres, A. H. Castro Neto PHYSICAL REVIEW B 86, 155449 (2012)
   !> updates:
   !> 2021/Aug/30, This is included in WannierTools as an auxiliary tool.

   use wmpi
   use para
   implicit none

   !> file existence
   logical  :: exists
   integer :: ierr
   character(8) :: cht

   !> time measure
   real(dp) :: time_start, time_end, time_init

   ierr = 0
   cpuid= 0
   num_cpu= 1
   !> initial the environment of mpi
#if defined (MPI)
   call mpi_init(ierr)
   call mpi_comm_rank(mpi_cmw,cpuid,ierr)
   call mpi_comm_size(mpi_cmw,num_cpu,ierr)
#endif

   if (cpuid==0) open(unit=stdout, file='output.txt')

   !> if mpi initial wrong, alarm
   if (cpuid==0.and.ierr.ne.0)then
      write(stdout,*)'mpi initialize wrong'
      stop
   endif

   !call now(time_init)
   !call header

   !> print information for mpi
   if (cpuid==0) then
      write(stdout, '(1x, a, i5, a)')'You are using ', num_cpu, ' CPU cores'
      write(stdout, *)' '
   endif

   !> readin the necessary parameters specified by user
   !> including the twisted angle, stacking pattern
   call readinput

   call ek_bulk_line

   if (cpuid.eq.0)write(stdout,*)'Congratulations! you finished the calculation.'
   if (cpuid.eq.0)write(stdout,*)'See you next time :)'

#if defined (MPI)
   call mpi_finalize(ierr)
#endif

end  !<< end of main program
