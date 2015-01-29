!> get the hopping parameter for two-atoms wannier90_hr.dat
  subroutine parse
     use para
     implicit none

     integer :: i
     integer :: j
     integer :: ir

     integer :: natoms
     integer :: max_projs

     integer, allocatable :: R0(:)
     integer, allocatable :: R1(:)
     complex(dp), allocatable :: Hmn(:, :)

     allocate(R0(3), R1(3))
     allocate(Hmn(Num_wann, Num_wann))

     do ir=1, nrpts
        if (sum(abs(irvec(:, ir)))<0.1) then
           do i=1, Num_wann
              write(11, '(100f6.2)') real(HmnR(i, :, ir))
           enddo
        endif
     enddo


     open (unit=10, file='hopping.dat')
     write(10, '(a)')'# export hopping parameters from hr file'

     !> R= (0, 0, 0)
     R0= (/0, 0, 0/)
     R1= (/0, 0, 0/)
     do ir=1, nrpts
        if (irvec(1, ir)==0 .and. irvec(2, ir)==0 .and.irvec(3, ir)==0) then
           Hmn= HmnR(:, :, ir)
        endif
     enddo

     call parse_H(Hmn, R0, R1)



     close(10)


     !>--------------------------------------------------------------
     !>> onsite hopping atom1-atom2
     !>--------------------------------------------------------------
     return

  end subroutine parse

