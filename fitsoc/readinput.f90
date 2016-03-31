!
! this subroutine is used to read some paramters from
! input.dat
! constructed on 4/22/2010 by QS.Wu 


  subroutine readinput

     use para
     implicit none

     character*20 :: fname='fitsoc.in'
     character*25 :: char_temp 

     logical :: exists
     integer :: i
     integer :: j
     integer :: NN

     real(dp) :: t1

     real(dp) :: temp
     real(dp) :: cell_volume
     real(dp) :: pos(3)

     real(dp) :: k1(3), k2(3)
     real(dp) :: kstart(3), kend(3)

     inquire(file=fname,exist=exists)
     if (exists)then
        write(*,*) 'read some paramters from input.dat'
        open(unit=1001,file=fname,status='old')
     else
        write(*,*)'file ' ,fname, ' dosnot exist'
        stop
     endif
 

     read(1001,*)infilename(:)
     write(*,'(a,a25)')' input file:',infilename
     read(1001,*)E_fermi
     write(*,*)'E_fermi',E_fermi
     read(1001,*)Nk
     write(*,*)'Nk', Nk

     !> read atom position
     read(1001, *)Num_atoms
     allocate(atom_name(Num_atoms))
     allocate(Atom_position(3, Num_atoms))
     do i=1, Num_atoms
        read(1001, *) atom_name(i), Atom_position(:, i)
        write(*, '(a4,3f6.3)')atom_name(i), Atom_position(:, i)
        pos= Atom_position(:, i)
     enddo

     Num_atom_type= 1
     do i=1, Num_atoms- 1
        if (atom_name(i).ne.atom_name(i+1)) Num_atom_type= Num_atom_type+ 1
     enddo
     print *, Num_atom_type
     allocate(atom_type(Num_atoms))

     Num_atom_type= 1
     atom_type(1)= 1
     do i=1, Num_atoms- 1
        if (atom_name(i).ne.atom_name(i+1)) Num_atom_type= Num_atom_type+ 1
        atom_type(i+1)= Num_atom_type
     enddo

     !> lattice information
     read(1001, *)Rua
     read(1001, *)Rub
     read(1001, *)Ruc

     !> transform lattice from direct space to reciprocal space

     Kua= 0d0
     Kub= 0d0
     Kuc= 0d0
     cell_volume= Rua(1)*(Rub(2)*Ruc(3)- Rub(3)*Ruc(2)) &
                 +Rua(2)*(Rub(3)*Ruc(1)- Rub(1)*Ruc(3)) &
                 +Rua(3)*(Rub(1)*Ruc(2)- Rub(2)*Ruc(1)) 
     cell_volume= 2d0*3.1415926535d0/cell_volume
     Kua(1)= cell_volume*(Rub(2)*Ruc(3)- Rub(3)*Ruc(2))
     Kua(2)= cell_volume*(Rub(3)*Ruc(1)- Rub(1)*Ruc(3))
     Kua(3)= cell_volume*(Rub(1)*Ruc(2)- Rub(2)*Ruc(1))

     Kub(1)= cell_volume*(Ruc(2)*Rua(3)- Ruc(3)*Rua(2))
     Kub(2)= cell_volume*(Ruc(3)*Rua(1)- Ruc(1)*Rua(3))
     Kub(3)= cell_volume*(Ruc(1)*Rua(2)- Ruc(2)*Rua(1))

     Kuc(1)= cell_volume*(Rua(2)*Rub(3)- Rua(3)*Rub(2))
     Kuc(2)= cell_volume*(Rua(3)*Rub(1)- Rua(1)*Rub(3))
     Kuc(3)= cell_volume*(Rua(1)*Rub(2)- Rua(2)*Rub(1))

     write(*, '(a)') '>> lattice information'
     write(*, '(3f10.6)')Rua
     write(*, '(3f10.6)')Rub
     write(*, '(3f10.6)')Ruc

     write(*, '(a)') '>> Reciprocal lattice information'
     write(*, '(3f10.6)')Kua
     write(*, '(3f10.6)')Kub
     write(*, '(3f10.6)')Kuc

     !> kline for 3d band structure
     !> high symmetry k points
     read(1001, *) nk3lines
     allocate(k3line_start(3, nk3lines))
     allocate(k3line_end(3, nk3lines))
     allocate(k3line_name(nk3lines+1))
     allocate(k3line_stop(nk3lines+1))
     k3line_stop= 0d0
     k3line_start= 0d0
     k3line_end= 0d0
     k3line_name= ' '
     do i=1, nk3lines
        read(1001, *) k3line_name(i), k3line_start(:, i), &
                      char_temp, k3line_end(:, i)
     enddo
     k3line_name(nk3lines+1)= char_temp

     NN= Nk
     nk3_band= NN*nk3lines
     allocate(k3len(nk3_band))
     allocate(k3points(3, nk3_band))
     k3len=0d0
     k3points= 0d0
     t1= 0d0
     do j=1, nk3lines
        do i=1, NN
           kstart= k3line_start(:, j)
           kend  = k3line_end(:, j)
           k1= kstart(1)*Kua+ kstart(2)*Kub+ kstart(3)*Kuc
           k2= kend(1)*Kua+ kend(2)*Kub+ kend(3)*Kuc

           k3points(:, i+ (j-1)*NN)= kstart+ (kend- kstart)*dble(i-1)/dble(NN-1)
           
           temp= dsqrt((k2(1)- k1(1))**2 &
                 +(k2(2)- k1(2))**2  &
                 +(k2(3)- k1(3))**2)/dble(NN-1) 

           if (i.gt.1) then
              t1=t1+temp
           endif
           k3len(i+(j-1)*NN)= t1
        enddo
        k3line_stop(j+1)= t1
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
        write(*, '(100a6)') char_temp, proj_name(1:nprojs(i), i)
     enddo

     !> spin orbital coupling strength
     allocate(lambda_p(Num_atom_type))
     allocate(lambda_d(Num_atom_type))
     lambda_p= 0d0
     lambda_d= 0d0
     read(1001, *) lambda_p(:)
     write(*, '(a,100f7.4)')'lambda_p input', lambda_p
     read(1001, *) lambda_d(:)
     write(*, '(a,100f7.4)')'lambda_d input', lambda_d

     write(*,*)'read input.dat file successfully'

     return
  end subroutine



