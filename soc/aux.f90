!>> some auxilary subroutines for time or others

!>> Convert the current wall-time into a real number
!>  The unit is second.

  subroutine now(time_now)

     use wmpi
     use para, only : dp
     implicit none
     integer   :: time_new(8)
     real(dp)  :: time_now
     call Date_and_time(values=time_new)
     time_now= time_new(3)*24*3600+time_new(5)*3600+&
               time_new(6)*60+time_new(7)+time_new(8)/1000d0
  
     return
  end subroutine now

  !>  The unit is second.
  subroutine print_time_cost(time_start, time_end, subname)

     use wmpi
     use para, only : dp, stdout, cpuid
     implicit none
     real(dp), intent(in)   :: time_start
     real(dp), intent(in)   :: time_end
     character(*) :: subname

     if (cpuid==0)write(stdout, 101)'Time cost for ', subname, ' is about ', &
        time_end- time_start, ' s' 

     101 format(1x, 3a, f16.3, a)
     return
  end subroutine print_time_cost

  !> print header
  subroutine header
     use wmpi
     use para, only : stdout, cpuid
     implicit none
     if (cpuid==0) then
        write(stdout, '(a)') " -----------------------------------------------------------------------"
        write(stdout, '(a)') "  W            W            W              TTTTTTTTTTTTTTTTTTTT         "
        write(stdout, '(a)') "   W          W W          W                        TT                  "
        write(stdout, '(a)') "    W        W   W        W                         TT                  "
        write(stdout, '(a)') "     W      W     W      W                          TT                  "
        write(stdout, '(a)') "      W    W       W    W                           TT                  "
        write(stdout, '(a)') "       W  W         W  W                            TT                  "
        write(stdout, '(a)') "        WW           WW                             TT                  "
        write(stdout, '(a)') "        W             W                             TT                  "
        write(stdout, '(a)') "                                                                        "
        write(stdout, '(a)') "                        Welcome to Wannier_tools.                       "
        write(stdout, '(a)') "                             Version 2.2.1                              "
        write(stdout, '(a)') "                   Tools for topological novel materials.               "
        write(stdout, '(a)') "                          Enjoy it and good luck.                       "
        write(stdout, '(a)') "                           Author : QuanSheng Wu                        "
        write(stdout, '(a)') "                        Email : wuquansheng@gmail.com                   "
        write(stdout, '(a)') "               Find more information on www.wanniertools.com            "
        write(stdout, '(a)') " ======================================================================="
        write(stdout, '(a)') "                                                                        "
     endif
  end subroutine header

