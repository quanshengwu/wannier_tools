  program main

     implicit none

     call readinput

     !> for InAs  14 spin-orbitals
      call readHmnR
     !call readHmnR_updn
     !call addsoc_magmom_all_3d
      call addsoc_all

    !call readHmnR_2d
    !call addsoc_all_2d
    !call addsoc_magmom_all_2d
     !> for WTe2 88 spin-orbitals
    !call addsoc_pd_zjw

  end 
