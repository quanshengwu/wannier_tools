!
! this subroutine is used to read some paramters from
! input.dat
! constructed on 4/22/2010 by QS.Wu 


  subroutine readinput

     use para
     implicit none

     character*20 :: fname='gethopping.in'
     character*25 :: char_temp 

     logical :: exists
     integer :: i

     real(dp) :: pos(3)

     inquire(file=fname,exist=exists)
     if (exists)then
        write(*,*) 'read some paramters from input.dat'
        open(unit=1001,file=fname,status='old')
     else
        write(*,*)'file ' ,fname, ' dosnot exist'
        stop
     endif
 

     read(1001,*)infilename
     write(*,'(a,a25)')' input file:',infilename
     read(1001,*)Soc
     write(*,*)'soc',Soc
     read(1001,*)E_fermi
     write(*,*)'E_fermi',E_fermi

     !> lattice information
     read(1001, *)Rua
     read(1001, *)Rub
     read(1001, *)Ruc

     write(*, '(a)') '>> lattice information'
     write(*, '(3f10.6)')Rua
     write(*, '(3f10.6)')Rub
     write(*, '(3f10.6)')Ruc

     !> read atom position
     read(1001, *)Num_atoms
     allocate(atom_name(Num_atoms))
     allocate(Atom_position(3, Num_atoms))
     do i=1, Num_atoms
        read(1001, *) atom_name(i), Atom_position(:, i)
        write(*, '(a4,3f6.3)')atom_name(i), Atom_position(:, i)
        pos= Atom_position(:, i)
        Atom_position(:, i)= pos(1)*Rua+ pos(2)*Rub+ pos(3)*Ruc
     enddo


     !> read projectors
     allocate(nprojs(Num_atoms))
     nprojs= 0
     read(1001, *)nprojs

     max_projs= maxval(nprojs)
     allocate(proj_name(max_projs, Num_atoms))
     proj_name= ' '
     do i=1, Num_atoms
        read(1001, *)char_temp, proj_name(1:nprojs(i), i)
     enddo
     
     write(*,*)'read input.dat file successfully'

     return
  end subroutine



