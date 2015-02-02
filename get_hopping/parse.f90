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
     max_NN= 3
     max_neighbours= 20

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
           enddo
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

     !> R= (0, 0, 0)
     R0= (/0, 0, 0/)
     R1= (/0, 0, 0/)
     do ir=1, nrpts
        if (irvec(1, ir)==0 .and. irvec(2, ir)==0 .and.irvec(3, ir)==0) then
           Hmn= HmnR(:, :, ir)
        endif
     enddo

     call parse_H(Hmn, R0, R1)


     !> R= (0, 0, 1)


     close(10)


     !>--------------------------------------------------------------
     !>> onsite hopping atom1-atom2
     !>--------------------------------------------------------------
     return

  end subroutine parse

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

