  program main

     implicit none

     call readinput
     call readHmnR

     !> for InAs  14 spin-orbitals
      call addsoc_all

     !> for WTe2 88 spin-orbitals
    !call addsoc_pd_zjw

  end 
