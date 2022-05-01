!>> some auxilary subroutines for time or others

!>> Convert the current wall-time into a real number
!>  The unit is second.

  subroutine now(time_now)

     use wmpi
     use para, only : dp
     implicit none

     integer   :: time_new(8)
     real(dp), intent(inout)  :: time_now

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

  subroutine printallocationinfo(variablename, ierr)
     use para, only : stdout

     implicit none

     character(len=*), intent(in) :: variablename
     integer, intent(in) :: ierr

     if (ierr/=0) then
        write(*, *)"ERROR: no enough memory for ", variablename, '  STAT=',ierr
#if defined (MPI)
        call mpi_finalize(ierr)
#endif
        stop
     endif

     return
  end subroutine printallocationinfo

  subroutine printerrormsg(errormsg)
     use para, only : stdout
     use wmpi, only : cpuid
     implicit none
     character(*), intent(in) :: errormsg
     integer :: ierr

     if (cpuid==0) then
        write(stdout, *) trim(errormsg)
     endif

#if defined (MPI)
        call mpi_finalize(ierr)
#endif
     stop
     return
  end subroutine printerrormsg




  !> print header
  subroutine header
     use wmpi
     use para, only : stdout, cpuid, version
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
        write(stdout, '(a)') "                        Welcome to WannierTools.                       "
        write(stdout, '(a,a10)') "                             Version ", version                             
        write(stdout, '(a)') "                   Tools for novel topological materials.               "
        write(stdout, '(a)') "                          Enjoy it and good luck.                       "
        write(stdout, '(a)') "                           Author : QuanSheng Wu                        "
        write(stdout, '(a)') "                        Email : wuquansheng@gmail.com                   "
        write(stdout, '(a)') "               Find more information on www.wanniertools.com            "
        write(stdout, '(a)') " ======================================================================="
        write(stdout, '(a)') "                                                                        "
     endif
  end subroutine header

  !> print footer
  subroutine footer
     use wmpi
     use para
     implicit none
     if (cpuid==0) then
        write(stdout, '(2x,a)') ''
        write(stdout, '(2x,a)') ''
        write(stdout, '(2x,a)') "======================================================================="
        write(stdout, '(2x,a)') 'Congratulations! you finished the calculation.'
        write(stdout, '(2x,a)') "I hope you could find something useful from this calculation."
        write(stdout, '(2x,a)') "If so, I would like to ask you to cite our this paper:"
        write(stdout, '(2x,a)') ''
        write(stdout, '(2x,a)') "WannierTools : An open-source software package for novel topological materials"
        write(stdout, '(2x,a)') "QuanSheng Wu and ShengNan Zhang and Hai-Feng Song and Matthias Troyer and Alexey A. Soluyanov"
        write(stdout, '(2x,a)') "Computer Physics Communications 224, 405 (2018)"
        write(stdout, '(2x,a)') "https://doi.org/10.1016/j.cpc.2017.09.033"
        if (Boltz_OHE_calc .or. Boltz_k_calc .or. Boltz_evolve_k) then
           write(stdout, '(2x,a)') "Please also cite this paper for magnetoresistance calculation"
           write(stdout, '(2x,a)') "Magnetoresistance from Fermi surface topology" 
           write(stdout, '(2x,a)') "ShengNan Zhang, QuanSheng Wu, Yi Liu, and Oleg V. Yazyev"
           write(stdout, '(2x,a)') "Phys. Rev. B 99, 035142 (2019) DOI:10.1103/PhysRevB.99.035142"
        endif
        write(stdout, '(2x,a)') ''
        write(stdout, '(2x,a)') "For bugs, please report to wuquansheng@gmail.com                   "
        write(stdout, '(2x,a)') "or wanniertools@groups.google.com.                 "
        write(stdout, '(2x,a)') "More information could find on www.wanniertools.com            "
        write(stdout, '(2x,a)') 'See you next time :)'
        write(stdout, '(2x,a)') "======================================================================="
     endif
  end subroutine footer


  !> Sorting arr in ascending order
   subroutine sortheap(n, arr)
      use para, only : dp
      implicit none
      integer, intent(in) :: n
      real(dp), intent(inout) :: arr(n)

      !* local variables
      integer :: i

      do i=n/2, 1, -1
         call sift_down(i, n)
      enddo

      do i=n, 2, -1
         call swap(arr(1), arr(i))
         call sift_down(1, i-1)
      enddo
      contains
      subroutine sift_down(l, r)
         use para, only : dp
         integer, intent(in) :: l, r
         integer :: j, jold
         real(dp) :: a
         a= arr(l)
         jold= l
         j= l+ l

         do while (j<=r)
            if (j<r) then
               if (arr(j)<arr(j+1))j=j+1
            endif
            if (a>= arr(j))exit
            arr(jold)= arr(j)
            jold= j
            j= j+ j
         enddo
         arr(jold)= a
         return
      end subroutine sift_down
   end subroutine sortheap

   !>> swap two real numbers
   subroutine swap(a, b)
      use para, only : dp
      real(dp), intent(inout) :: a
      real(dp), intent(inout) :: b
      real(dp) :: c
      c=a
      a=b
      b=c
      return
   end subroutine swap


