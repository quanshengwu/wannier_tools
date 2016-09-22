!> get the hopping parameter for two-atoms wannier90_hr.dat
  subroutine parse
     use para
     implicit none

     integer :: i
     integer :: j
     integer :: n
     integer :: ia
     integer :: ib
     integer :: ic
     integer :: it
     integer :: ir
     integer :: iR0

     !> maximun order of nearest neighbour
     integer :: max_NN
     integer :: max_neighbours

     integer, allocatable :: R0(:)
     integer, allocatable :: R1(:)
     complex(dp), allocatable :: Hmn(:, :)

     type nearest_atoms
        !> the i'th nearest neighbour
        integer :: i

        !> howmany neighbours
        integer :: numbers

        integer, allocatable :: neighbour(:)
        character(4), allocatable :: neighbour_name(:)

        !> the distance between two atoms
        real(dp) :: distance

        !> get the lattice vector for each neighbours
        integer, allocatable :: iR(:)
     end type nearest_atoms

     !> the first dimension is the order ot the neighbour
     !> the second dimension is the atoms number
     type(nearest_atoms), allocatable :: neighbours(:,:)

     !> the distance between two atoms for each species of atom
     real(dp), allocatable :: distance(:, :)

     real(dp) :: dis
     real(dp) :: pos1(3)
     real(dp) :: pos2(3)
     max_NN= 10
     max_neighbours= 40

     allocate(neighbours(max_NN, Num_atoms))
     allocate(distance(nrpts*Num_atoms, Num_atoms))
     distance= 0d0

     do n=1, max_NN
        do i=1, Num_atoms
           allocate(neighbours(n, i)%iR(max_neighbours))
           allocate(neighbours(n, i)%neighbour(max_neighbours))
           allocate(neighbours(n, i)%neighbour_name(max_neighbours))
           neighbours(n, i)%neighbour_name= ' '
           neighbours(n, i)%neighbour= 0
           neighbours(n, i)%iR= 0
           neighbours(n, i)%distance= 0d0
           neighbours(n, i)%numbers = 0
        enddo
     enddo

     !> first step is to calculate all the distance from different atoms
     !> sweep over atoms in unit cell
     do i= 1, Num_atoms
        !> the first atom's position 
        pos1= Atom_position(:, i)
        it= 0
        do ir=1, nrpts
           ia= irvec(1, ir)
           ib= irvec(2, ir)
           ic= irvec(3, ir)
           do j=1, Num_atoms
              it= it+ 1
              !> the second atom's position
              pos2= ia*Rua+ ib*Rub+ ic*Ruc+ Atom_position(:, j)
              dis= (pos2(1)-pos1(1))**2+(pos2(2)-pos1(2))**2 + (pos2(3)-pos1(3))**2
              distance(it, i)= sqrt(dis)
           enddo ! j
        enddo ! ir
     enddo ! ia

     do i=1, Num_atoms
        call sortheap(nrpts*Num_atoms, distance(:, i))
     enddo

     !> get the distance for each kind of neighbour
     !> nearest and next nearest and next-next nearest neighbour
     do i=1, Num_atoms
        it= 0
        dis= 0
        do j=1, nrpts*Num_atoms
           if (abs(distance(j, i)-dis)>0.01d0)then
              it= it+ 1
              dis= distance(j, i)
              if (it> max_NN) exit
              neighbours(it, i)%i= it  ! the it'th nearest neighbour
              neighbours(it, i)%distance= distance(j, i)
           endif
        enddo ! j
     enddo ! i

 
     !> Calculate all the distance from different atoms
     !> sweep over atoms in unit cell
     !> to get all neighbours information
     do n=1, max_NN
     do i= 1, Num_atoms
        !> the first atom's position 
        pos1= Atom_position(:, i)
        it= 0
        do ir=1, nrpts
           ia= irvec(1, ir)
           ib= irvec(2, ir)
           ic= irvec(3, ir)
           do j=1, Num_atoms
              !> the second atom's position
              pos2= ia*Rua+ ib*Rub+ ic*Ruc+ Atom_position(:, j)
              dis= (pos2(1)-pos1(1))**2+(pos2(2)-pos1(2))**2 + (pos2(3)-pos1(3))**2
              dis= sqrt(dis)
              if (abs(dis- neighbours(n, i)%distance)<0.01)then
                 it= it+ 1
                 neighbours(n, i)%numbers= neighbours(n, i)%numbers+ 1
                 neighbours(n, i)%neighbour(it) = j
                 neighbours(n, i)%neighbour_name(it) = atom_name(j)
                 neighbours(n, i)%iR(it) = ir
              endif
           enddo ! j
        enddo ! ir
     enddo ! i    
     enddo ! n    


     do i=1, Num_atoms
     do n=1, max_NN
        write(*, *)' '
        write(*, '(a, i2, 2a)')"The", n,"'th nearest neighbour for ", atom_name(i)
        write(*, '(a, i3, a)')"There are", neighbours(n, 1)%numbers, " neighbours"
        write(*, *)(neighbours(n, i)%neighbour_name(j), j=1, neighbours(n, i)%numbers)
     enddo !n 
     enddo !i

     allocate(R0(3), R1(3))
     allocate(Hmn(Num_wann, Num_wann))

     open (unit=10, file='hopping.dat')
     write(10, '(a)')'# export hopping parameters from hr file'

     do ir=1, nrpts
        if (irvec(1, ir)==0.and. irvec(2, ir)==0 .and.irvec(3, ir)==0)then
           iR0= ir
        endif
     enddo


     !> sweep over atoms in R=0
     do i=1, Num_atoms 

        call parse2(iR0, i, i)

        !> for each atoms, sweep its neighbour order
        write(10, '(a)') '#-------------------------------------------------------- '
        write(10, '(a, a, a)')'# hopping bewteen ', trim(Atom_name(i)), &
           ' and its neighbours'
        do n=1, max_NN
           write(10, *) ' '
           write(10, '(a, i3, a)')"#the ", n, "'th order"
           write(10, '(a, i3, 80(6a))')"#There are", neighbours(n, i)%numbers, &
              " neighbours: ", &
           (neighbours(n, i)%neighbour_name(j), j=1, neighbours(n, i)%numbers)
           !> for each order, sweep its neighbours
           do j=1, neighbours(n, i)%numbers
              write(10, *)' '
              write(10, '(a, i2, a, i3)')'# ', j, ' in', neighbours(n, i)%numbers 
              write(10, '(5a, f18.7, a)')'# hopping bewteen ', trim(Atom_name(i)),&
                   '-', trim(neighbours(n, i)%neighbour_name(j)), '   distance :', &
                   neighbours(n, i)%distance, ' Ang'
              iR= neighbours(n, i)%iR(j)
              call parse2(iR, i, neighbours(n, i)%neighbour(j))
           enddo ! j
        enddo ! n
     enddo ! i

     close(10)

     !>--------------------------------------------------------------
     !>> onsite hopping atom1-atom2
     !>--------------------------------------------------------------
     return

  end subroutine parse

!> get the hopping parameter for two-atoms wannier90_hr.dat
  subroutine parse2(iR1, ia1, ia2)
     use para
     implicit none

     integer :: i
     integer :: j
     integer :: ir

     integer :: natoms

     integer :: iatom1
     integer :: iatom2
     integer :: col_start
     integer :: col_end
     integer :: row_start
     integer :: row_end
     integer :: row_diff
     integer :: col_diff
     integer :: row_offset
     integer :: col_offset
     integer :: R0(3)
     integer :: R1(3)

     !> atom species
     integer, intent(in) :: ia1
     integer, intent(in) :: ia2
     integer, intent(in) :: iR1
     complex(dp), allocatable :: Hmn(:, :)

     complex(dp), allocatable :: Hsub(:, :)

     allocate(Hmn(Num_wann, Num_wann))
     Hmn= 0d0

     natoms= Num_atoms
     Hmn= HmnR(:, :, iR1) !/ndegen(iR1)
     R1= irvec(:, iR1)
     write(10, '(a, 8i5)')'Hmn ', R1
     do i=1, Num_wann/2
        write(10, '(10f8.4)')(real(Hmn(i, j)), j=1, Num_wann/2)
     enddo

     allocate(Hsub(max_projs*soc, max_projs*soc))
     Hsub= 0d0

     R0= 0

     !> for Zincblende InSb  InAs GaSb AlSb, the projectors are
     !> In(Ga, Al) s, px, py, pz  Sb(As) px, py, pz

     100 format('# hopping between ', a, '-', a, ' from', &
        ' (', 3i3, ') to (', 3i3, ')')
     101 format(2x, 'Real', 2x, 100a12)
     1010 format(2x, 'Imag', 2x, 100a12)
     102 format(a8, 100f12.6 )


     !>--------------------------------------------------------------
     !>> onsite hopping atom1-atom2
     !>--------------------------------------------------------------
     iatom1= ia1  !> row
     iatom2= ia2  !> col
     row_offset= 0
     do i=1, ia1-1
        row_offset= row_offset+ nprojs(i)
     enddo

     col_offset= 0
     do i=1, ia2-1
        col_offset= col_offset+ nprojs(i)
     enddo

     !> for Hmn
     row_start= 1+ row_offset
     row_end= nprojs(iatom1)+ row_offset
     col_start= 1+ col_offset
     col_end= nprojs(iatom2)+ col_offset

     !> for Hsub
     row_diff= row_end- row_start+ 1
     col_diff= col_end- col_start+ 1
     do ir=1, nrpts
        if (sum(abs(irvec(:, ir)))<0.1) then
           Hsub(1: row_diff, 1:col_diff) &
              = Hmn(row_start:row_end, col_start:col_end)
           if (soc==2) then
           Hsub(row_diff+1:2*row_diff, col_diff+1:2*col_diff) &
              = Hmn(row_start+ Num_wann/2: row_end+ Num_wann/2, &
                     col_start+ Num_wann/2: col_end+ Num_wann/2)
           Hsub(1:row_diff, col_diff+1:2*col_diff) &
              = Hmn(row_start: row_end, &
                     col_start+ Num_wann/2: col_end+ Num_wann/2)
           Hsub(row_diff+1:2*row_diff, 1:col_diff) &
              = Hmn(row_start+ Num_wann/2: row_end+ Num_wann/2, &
                     col_start: col_end)
           endif ! soc
        endif ! R= (0, 0, 0)
     enddo ! ir

     !> Get Slater-Koster parameters
     call get_SK(Hmn, iR1, iatom1, iatom2)
     !> export onsite hopping for atom2
     write(10, 100)trim(atom_name(iatom1)), trim(atom_name(iatom2)), R0, R1
     !> no soc
     if (soc==1) then
        write(10, 101)proj_name(:, iatom2)
        do i=1, row_end- row_start+ 1
           write(10, 102)proj_name(i, iatom1), &
              real(Hsub(i, 1: col_diff))
        enddo
   
        write(10, 1010)proj_name(:, iatom2)
        do i=1, row_end- row_start+ 1
           write(10, 102)proj_name(i, iatom1), &
              aimag(Hsub(i, 1: col_diff))
        enddo
     else
        write(10, 101)proj_name(1:nprojs(iatom2), iatom2), &
           proj_name(1:nprojs(iatom2), iatom2)
        do i=1, row_diff
           write(10, 102)proj_name(i, iatom1), &
              real(Hsub(i, 1: 2*col_diff))
        enddo
        do i=1, row_diff
           write(10, 102)proj_name(i, iatom1), &
              real(Hsub(i+row_diff, 1: 2*col_diff))
        enddo
 
        write(10, 1010)proj_name(1:nprojs(iatom2), iatom2), &
           proj_name(1:nprojs(iatom2), iatom2)
        do i=1, row_diff
           write(10, 102)proj_name(i, iatom1), &
              aimag(Hsub(i, 1: 2*col_diff))
        enddo
        do i=1, row_diff
           write(10, 102)proj_name(i, iatom1), &
              aimag(Hsub(i+row_diff, 1: 2*col_diff))
        enddo
     
     endif


     !>--------------------------------------------------------------
     !>> onsite hopping atom1-atom2
     !>--------------------------------------------------------------
     return

  end subroutine parse2


!> get the hopping parameter for two-atoms wannier90_hr.dat
  subroutine parse_H(Hmn, R0, R1)
     use para
     implicit none

     integer :: i
     integer :: ir

     integer :: natoms

     integer :: iatom1
     integer :: iatom2
     integer :: col_start
     integer :: col_end
     integer :: row_start
     integer :: row_end
     integer :: row_diff
     integer :: col_diff
     integer :: row_offset
     integer :: col_offset

     integer, intent(in) :: R0(3)
     integer, intent(in) :: R1(3)
     complex(dp), intent(in) :: Hmn(Num_wann, Num_wann)

     complex(dp), allocatable :: Hsub(:, :)


     natoms= Num_atoms

     allocate(Hsub(max_projs*soc, max_projs*soc))

     !> for Zincblende InSb  InAs GaSb AlSb, the projectors are
     !> In(Ga, Al) s, px, py, pz  Sb(As) px, py, pz

     100 format('# hopping between ', a, '-', a, ' from', &
        ' (', 3i3, ') to (', 3i3, ')')
     101 format(2x, 'Real', 2x, 100a8)
     1010 format(2x, 'Imag', 2x, 100a8)
     102 format(a8, 100f18.7 )


     !>--------------------------------------------------------------
     !>> onsite hopping atom1-atom1
     !>--------------------------------------------------------------
     iatom1= 1
     iatom2= 1
     col_start= 1
     row_start= 1
     col_end= nprojs(iatom1)
     row_end= nprojs(iatom2)
     row_diff= row_end- row_start+ 1
     col_diff= col_end- col_start+ 1
     do ir=1, nrpts
        if (sum(abs(irvec(:, ir)))<0.1d0) then
           Hsub(1: row_diff, 1:col_diff) &
              = Hmn(row_start:row_end, col_start:col_end)
           if (soc==2) then
           Hsub(row_diff+1:2*row_diff, col_diff+1:2*col_diff) &
              = Hmn(row_start+ Num_wann/2: row_end+ Num_wann/2, &
                     col_start+ Num_wann/2: col_end+ Num_wann/2)
           Hsub(1:row_diff, col_diff+1:2*col_diff) &
              = Hmn(row_start: row_end, &
                     col_start+ Num_wann/2: col_end+ Num_wann/2)
           Hsub(row_diff+1:2*row_diff, 1:col_diff) &
              = Hmn(row_start+ Num_wann/2: row_end+ Num_wann/2, &
                     col_start: col_end)
           endif
        endif
     enddo

     !> export onsite hopping for atom2
     write(10, 100)trim(atom_name(iatom1)), trim(atom_name(iatom2)), R0, R1
     if (soc==1) then
        write(10, 101)proj_name(1:nprojs(iatom2), iatom2)
        do i=1, row_end- row_start+ 1
           write(10, 102)proj_name(i, iatom1), &
              real(Hsub(i, 1: col_diff))
        enddo
   
        write(10, 1010)proj_name(1:nprojs(iatom2), iatom2)
        do i=1, row_end- row_start+ 1
           write(10, 102)proj_name(i, iatom1), &
              aimag(Hsub(i, 1: col_diff))
        enddo
     else
        write(10, 101)proj_name(1:nprojs(iatom2), iatom2), &
           proj_name(1:nprojs(iatom2), iatom2)
        do i=1, row_diff
           write(10, 102)proj_name(i, iatom1), &
              real(Hsub(i, 1: 2*col_diff))
        enddo
        do i=1, row_diff
           write(10, 102)proj_name(i, iatom1), &
              real(Hsub(i+row_diff, 1: 2*col_diff))
        enddo
 
        write(10, 1010)proj_name(1:nprojs(iatom2), iatom2), &
           proj_name(1:nprojs(iatom2), iatom2)
        do i=1, row_diff
           write(10, 102)proj_name(i, iatom1), &
              aimag(Hsub(i, 1: 2*col_diff))
        enddo
        do i=1, row_diff
           write(10, 102)proj_name(i, iatom1), &
              aimag(Hsub(i+row_diff, 1: 2*col_diff))
        enddo
     
     endif

     !>--------------------------------------------------------------
     !>> onsite hopping atom2-atom2
     !>--------------------------------------------------------------
     iatom1= 2
     iatom2= 2
     row_offset= nprojs(1)
     col_offset= nprojs(1)

     !> for Hmn
     col_start= 1+ col_offset
     row_start= 1+ row_offset
     col_end= nprojs(iatom1)+ col_offset
     row_end= nprojs(iatom2)+ row_offset

     !> for Hsub
     row_diff= row_end- row_start+ 1
     col_diff= col_end- col_start+ 1
     Hsub= 0d0
     do ir=1, nrpts
        if (sum(abs(irvec(:, ir)))<0.1) then
           Hsub(1: row_diff, 1:col_diff) &
              = Hmn(row_start:row_end, col_start:col_end)
           if (soc==2) then
           Hsub(row_diff+1:2*row_diff, col_diff+1:2*col_diff) &
              = Hmn(row_start+ Num_wann/2: row_end+ Num_wann/2, &
                     col_start+ Num_wann/2: col_end+ Num_wann/2)
           Hsub(1:row_diff, col_diff+1:2*col_diff) &
              = Hmn(row_start: row_end, &
                     col_start+ Num_wann/2: col_end+ Num_wann/2)
           Hsub(row_diff+1:2*row_diff, 1:col_diff) &
              = Hmn(row_start+ Num_wann/2: row_end+ Num_wann/2, &
                     col_start: col_end)
           endif ! soc
        endif ! R= (0, 0, 0)
     enddo ! ir

     !> export onsite hopping for atom2
     write(10, *) ' '
     write(10, 100)trim(atom_name(iatom1)), trim(atom_name(iatom2)), R0, R1
     !> no soc
     if (soc==1) then
        write(10, 101)proj_name(:, iatom1)
        do i=1, row_end- row_start+ 1
           write(10, 102)proj_name(i, iatom2), &
              real(Hsub(i, 1: col_diff))
        enddo
   
        write(10, 1010)proj_name(:, iatom1)
        do i=1, row_end- row_start+ 1
           write(10, 102)proj_name(i, iatom2), &
              aimag(Hsub(i, 1: col_diff))
        enddo
     else
        write(10, 101)proj_name(1:nprojs(iatom1), iatom1), &
           proj_name(1:nprojs(iatom1), iatom1)
        do i=1, row_diff
           write(10, 102)proj_name(i, iatom2), &
              real(Hsub(i, 1: 2*col_diff))
        enddo
        do i=1, row_diff
           write(10, 102)proj_name(i, iatom2), &
              real(Hsub(i+row_diff, 1: 2*col_diff))
        enddo
 
        write(10, 1010)proj_name(1:nprojs(iatom1), iatom1), &
           proj_name(1:nprojs(iatom1), iatom1)
        do i=1, row_diff
           write(10, 102)proj_name(i, iatom2), &
              aimag(Hsub(i, 1: 2*col_diff))
        enddo
        do i=1, row_diff
           write(10, 102)proj_name(i, iatom2), &
              aimag(Hsub(i+row_diff, 1: 2*col_diff))
        enddo
     
     endif




     !>--------------------------------------------------------------
     !>> onsite hopping atom1-atom2
     !>--------------------------------------------------------------
     iatom1= 1
     iatom2= 2
     row_offset= 0
     col_offset= nprojs(1)

     !> for Hmn
     col_start= 1+ col_offset
     col_end= nprojs(iatom2)+ col_offset
     row_start= 1+ row_offset
     row_end= nprojs(iatom1)+ row_offset
     print *, row_start, row_end
     print *, col_start, col_end

     !> for Hsub
     row_diff= row_end- row_start+ 1
     col_diff= col_end- col_start+ 1
     Hsub= 0d0
     do ir=1, nrpts
        if (sum(abs(irvec(:, ir)))<0.1) then
           Hsub(1: row_diff, 1:col_diff) &
              = Hmn(row_start:row_end, col_start:col_end)
           if (soc==2) then
           Hsub(row_diff+1:2*row_diff, col_diff+1:2*col_diff) &
              = Hmn(row_start+ Num_wann/2: row_end+ Num_wann/2, &
                     col_start+ Num_wann/2: col_end+ Num_wann/2)
           Hsub(1:row_diff, col_diff+1:2*col_diff) &
              = Hmn(row_start: row_end, &
                     col_start+ Num_wann/2: col_end+ Num_wann/2)
           Hsub(row_diff+1:2*row_diff, 1:col_diff) &
              = Hmn(row_start+ Num_wann/2: row_end+ Num_wann/2, &
                     col_start: col_end)
           endif ! soc
        endif ! R= (0, 0, 0)
     enddo ! ir

     !> export onsite hopping for atom2
     write(10, *) ' '
     write(10, 100)trim(atom_name(iatom1)), trim(atom_name(iatom2)), R0, R1
     !> no soc
     if (soc==1) then
        write(10, 101)proj_name(:, iatom2)
        do i=1, row_end- row_start+ 1
           write(10, 102)proj_name(i, iatom1), &
              real(Hsub(i, 1: col_diff))
        enddo
   
        write(10, 1010)proj_name(:, iatom2)
        do i=1, row_end- row_start+ 1
           write(10, 102)proj_name(i, iatom1), &
              aimag(Hsub(i, 1: col_diff))
        enddo
     else
        write(10, 101)proj_name(1:nprojs(iatom2), iatom2), &
           proj_name(1:nprojs(iatom2), iatom2)
        do i=1, row_diff
           write(10, 102)proj_name(i, iatom1), &
              real(Hsub(i, 1: 2*col_diff))
        enddo
        do i=1, row_diff
           write(10, 102)proj_name(i, iatom1), &
              real(Hsub(i+row_diff, 1: 2*col_diff))
        enddo
 
        write(10, 1010)proj_name(1:nprojs(iatom2), iatom2), &
           proj_name(1:nprojs(iatom2), iatom2)
        do i=1, row_diff
           write(10, 102)proj_name(i, iatom1), &
              aimag(Hsub(i, 1: 2*col_diff))
        enddo
        do i=1, row_diff
           write(10, 102)proj_name(i, iatom1), &
              aimag(Hsub(i+row_diff, 1: 2*col_diff))
        enddo
     
     endif



     !>--------------------------------------------------------------
     !>> onsite hopping atom2-atom1
     !>--------------------------------------------------------------
     iatom1= 2
     iatom2= 1
     row_offset= nprojs(1)
     col_offset= 0

     !> for Hmn
     col_start= 1+ col_offset
     col_end= nprojs(iatom2)+ col_offset
     row_start= 1+ row_offset
     row_end= nprojs(iatom1)+ row_offset
     print *, row_start, row_end
     print *, col_start, col_end

     !> for Hsub
     row_diff= row_end- row_start+ 1
     col_diff= col_end- col_start+ 1
     Hsub= 0d0
     do ir=1, nrpts
        if (sum(abs(irvec(:, ir)))<0.1) then
           Hsub(1: row_diff, 1:col_diff) &
              = Hmn(row_start:row_end, col_start:col_end)
           if (soc==2) then
           Hsub(row_diff+1:2*row_diff, col_diff+1:2*col_diff) &
              = Hmn(row_start+ Num_wann/2: row_end+ Num_wann/2, &
                     col_start+ Num_wann/2: col_end+ Num_wann/2)
           Hsub(1:row_diff, col_diff+1:2*col_diff) &
              = Hmn(row_start: row_end, &
                     col_start+ Num_wann/2: col_end+ Num_wann/2)
           Hsub(row_diff+1:2*row_diff, 1:col_diff) &
              = Hmn(row_start+ Num_wann/2: row_end+ Num_wann/2, &
                     col_start: col_end)
           endif ! soc
        endif ! R= (0, 0, 0)
     enddo ! ir

     !> export onsite hopping for atom2
     write(10, *) ' '
     write(10, 100)trim(atom_name(iatom1)), trim(atom_name(iatom2)), R0, R1
     !> no soc
     if (soc==1) then
        write(10, 101)proj_name(:, iatom2)
        do i=1, row_end- row_start+ 1
           write(10, 102)proj_name(i, iatom1), &
              real(Hsub(i, 1: col_diff))
        enddo
   
        write(10, 1010)proj_name(:, iatom2)
        do i=1, row_end- row_start+ 1
           write(10, 102)proj_name(i, iatom1), &
              aimag(Hsub(i, 1: col_diff))
        enddo
     else
        write(10, 101)proj_name(1:nprojs(iatom2), iatom2), &
           proj_name(1:nprojs(iatom2), iatom2)
        do i=1, row_diff
           write(10, 102)proj_name(i, iatom1), &
              real(Hsub(i, 1: 2*col_diff))
        enddo
        do i=1, row_diff
           write(10, 102)proj_name(i, iatom1), &
              real(Hsub(i+row_diff, 1: 2*col_diff))
        enddo
 
        write(10, 1010)proj_name(1:nprojs(iatom2), iatom2), &
           proj_name(1:nprojs(iatom2), iatom2)
        do i=1, row_diff
           write(10, 102)proj_name(i, iatom1), &
              aimag(Hsub(i, 1: 2*col_diff))
        enddo
        do i=1, row_diff
           write(10, 102)proj_name(i, iatom1), &
              aimag(Hsub(i+row_diff, 1: 2*col_diff))
        enddo
     
     endif



     close(10)


     !>--------------------------------------------------------------
     !>> onsite hopping atom1-atom2
     !>--------------------------------------------------------------
     return

  end subroutine parse_H


  !> get Slater-Koster hopping parameters
  !> ref : doi:10.1103/PhysRevB.92.205108
  !> This version only works for InAs, GaSb, InSb, AlSb
  !> the orbital order is Cation s, px, py, pz, Anion px, py, pz
  subroutine get_SK(Hmn, iR, ia1, ia2)
     use para
     implicit none

     complex(dp), intent(in) :: Hmn(Num_wann, Num_wann)
     integer, intent(in) :: iR, ia1, ia2

     integer :: ia, ib, ic, i, j, it
     real(dp) :: pos1(3), pos2(3), vec(3)

     !> direction cosines
     !> l=x/r; m=y/r; n=z/r
     real(dp) :: l, m, n
     real(dp) :: cosines(3)
     !> norm of pos2-pos1
     real(dp) :: r

     logical :: ifs(Num_atoms)
     logical :: ifp(Num_atoms)

     !> the orbital index 
     integer :: sorbital(Num_atoms)
     integer :: porbital(Num_atoms)

     !> slater-koster parameters, only for s and p orbitals
     real(dp) :: Vsssigma
     real(dp) :: Vspsigma
     real(dp) :: Vppsigma, Vpppi

     real(dp) :: dis
     real(dp) :: delta(3,3)

     real(dp) :: diag_sum_pp
     real(dp) :: offdiag_sum_pp

     Vsssigma= 0d0
     Vspsigma= 0d0
     Vppsigma= 0d0
     Vpppi   = 0d0

     sorbital(1)=1   ! cation has  s orbital
     sorbital(2)=0   ! Anion has no s orbital

     porbital(1)=2
     porbital(2)=5

     ia= irvec(1, iR)
     ib= irvec(2, iR)
     ic= irvec(3, iR)
     pos1= Atom_position(:, ia1) 
     pos2= ia*Rua+ ib*Rub+ ic*Ruc+ Atom_position(:, ia2)
     vec= pos2- pos1
     dis= (pos2(1)-pos1(1))**2+(pos2(2)-pos1(2))**2 + (pos2(3)-pos1(3))**2
     dis= sqrt(dis)

     !> direction cosines
     r= sqrt(vec(1)*vec(1)+ vec(2)*vec(2)+ vec(3)*vec(3))
     if (r>0.001) then
        l= vec(1)/r
        m= vec(2)/r
        n= vec(3)/r
     else
        l=1d0
        m=1d0
        n=1d0
     endif
     cosines(1)=l
     cosines(2)=m
     cosines(3)=n

     !> check if there are s, p orbital in projector
     ifs=.FALSE.
     ifp=.FALSE.
     do ia=1, Num_atoms
        do j=1, nprojs(ia)
           if (index(proj_name(j, ia), 's')/=0) ifs(ia)=.TRUE.
           if (index(proj_name(j, ia), 'p')/=0) ifp(ia)=.TRUE.
        enddo !j
     enddo !ia

     !> ss sigma
     if (ifs(ia1).and. ifs(ia2)) Vsssigma= Hmn(sorbital(ia1), sorbital(ia2))

     !> sp sigma
     if (ifs(ia1).and. ifp(ia2)) then
        Vspsigma= 0d0
        j=0
        do i=0, 2
           if (abs(cosines(i+1))>0.001d0) then
              Vspsigma= Vspsigma+ Hmn(sorbital(ia1), porbital(ia2)+i )/cosines(i+1)  
              j= j+ 1
           endif
        enddo
        Vspsigma= Vspsigma/dble(j)
     endif

     !> sp sigma
     if (ifs(ia2).and. ifp(ia1)) then
        Vspsigma= 0d0
        j=0
        do i=0, 2
           if (abs(cosines(i+1))>0.001d0) then
              Vspsigma= Vspsigma+ Hmn(sorbital(ia2)+i, porbital(ia1) )/cosines(i+1)
              j= j+ 1
           endif
        enddo
        Vspsigma= Vspsigma/dble(j)
     endif

     !> pp hopping
     diag_sum_pp= 0d0
     if (ifp(ia1).and. ifp(ia2)) then
        do i=0, 2
           diag_sum_pp= diag_sum_pp+ Hmn(porbital(ia1)+i, porbital(ia2)+i)
        enddo
     endif

     offdiag_sum_pp= 0d0
     if (ifp(ia1).and. ifp(ia2)) then
        it= 0 
        do i=0, 2
           do j=0, 2
              if (i==j) cycle
              if (abs(cosines(i+1))<0.001 .or. abs(cosines(j+1))<0.001) cycle
              offdiag_sum_pp= offdiag_sum_pp+ &
                 Hmn(porbital(ia1)+i, porbital(ia2)+j)/cosines(i+1)/cosines(j+1)
              it= it+ 1
           enddo
        enddo
        if (it>0) offdiag_sum_pp= offdiag_sum_pp/dble(it)
     endif


     !> pp sigma
     Vppsigma= 0d0
     if (ifp(ia1).and. ifp(ia2)) then
       !do i=0, 2
       !   do j=0, 2
       !      Vppsigma= Vppsigma+ Hmn(porbital(ia1)+i, porbital(ia2)+j)*cosines(i+1)*cosines(j+1)
       !   enddo
       !enddo
        Vppsigma= (diag_sum_pp+ 2d0*offdiag_sum_pp)/3d0
     endif

     !> pp pi
     Vpppi= 0d0
     if (ifp(ia1).and. ifp(ia2)) then
       !do i=0, 2
       !   do j=0, 2
       !      Vpppi= Vpppi-Hmn(porbital(ia1)+i, porbital(ia2)+j)*cosines(i+1)*cosines(j+1)*0.5d0
       !   enddo
       !   Vpppi= Vpppi+ Hmn(porbital(ia1)+i, porbital(ia2)+i)*0.5d0
       !enddo

        Vpppi= (diag_sum_pp- offdiag_sum_pp)/3d0
     endif

     delta=0d0
     do i=1, 3
        delta(i, i)= 1d0
     enddo

     if (ifp(ia1).and. ifp(ia2)) then
        write(10, '(a)') '# P-P hoppings from slater-koster'
        do i=1, 3
           write(10, '(10f10.6)') ((Vppsigma-Vpppi)*cosines(i)*cosines(j)+ Vpppi*delta(i, j), j=1, 3)
        enddo
     endif

     !> write out the Slater-Koster parameters
     write(10, '(a, 3f10.6)')'# direction cosines: ', cosines
     write(10, '(a)')'# Slater koster parameters Vsssigma, Vspsigma, Vppsigma, Vpppi'
     write(10, '(5a, 10f10.6)')'#', trim(Atom_name(ia1)), '-', trim(Atom_name(ia2)), &
        '-SK : ', dis, Vsssigma, Vspsigma, Vppsigma, Vpppi

     return
  end subroutine get_SK

   !* heap sort algorithm see Numerical Reciples
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
         integer, intent(in) :: l, r
         integer :: j, jold
         real(8) :: a
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
      subroutine swap(a, b)
         real(8) :: a, b

         real(8) :: c
         c=a
         a=b
         b=c
      end subroutine swap
   end subroutine sortheap


