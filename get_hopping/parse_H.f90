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
     102 format(a8, 100f8.3 )


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

