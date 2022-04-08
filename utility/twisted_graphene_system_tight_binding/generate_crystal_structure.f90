   subroutine generate_crystal_structure()

      ! reference
      ! Graphene Bilayer with a Twist: Electronic Structure
      ! PRL 99, 256802 (2007)
      use para
      implicit none

      real(dp) :: lattice_constant, interlayer_distance, lc_z
      real(dp), parameter :: sqrt_3=dsqrt(3.0d0)
      real(dp) :: a1(3), a2(3), t1(3), t2(3), t3(3), A_site(3), B_site(3)
      real(dp) :: R1(3), R2(3), R3(3), Urot(3, 3)
      integer :: max_period, iter, i, j, m, ilayer, Nleft
      integer :: num_atom_per_layer, max_atoms_perlayer, num_atom_total
      real(dp), allocatable :: atoms_morie_cell_direct(:, :)

      real(dp), allocatable :: atoms_supercell_temp0(:, :)
      real(dp), allocatable :: atoms_supercell_temp(:, :)

      real(dp) :: theta, angle, tau(3), pos1(3), pos2(3), pos_cart(3), pos_direct(3)
      real(dp) :: xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz
      character(40) :: xyzfilename

      real(dp), external :: norm, get_angle
      real(dp) :: angle_t(160)
      integer :: m_index(160), number_of_atoms(160)


      !> primitive lattice vectors of Graphene 
      lattice_constant= lattice_constant_a
      interlayer_distance= lattice_constant_c
      a1(1)=lattice_constant*0.5d0*sqrt_3; a1(2)=lattice_constant*0.5d0; a1(3)= 0d0
      a2(1)=lattice_constant*0.5d0*sqrt_3; a2(2)=-lattice_constant*0.5d0; a2(3)= 0d0

      m=twisted_index_m
      theta= acos((3d0*m*m + 3d0*m + 0.5d0)/(3d0*m*m + 3d0*m + 1d0));

      do i=1,  160
         m_index(i)=i
         angle_t(i)= acos((3d0*i*i + 3d0*i + 0.5d0)/(3d0*i*i + 3d0*i + 1d0))*180d0/pi
         number_of_atoms(i)= 4*(3*i*i+3*i+1)
      enddo

      write(stdout, '(a)') ' '
      write(stdout, '(a)') '<< Index-angle table >>'
      do i=1, 16
         write(stdout, '(a10, 10i8)') ' m:', m_index((i-1)*10+1:i*10)
         write(stdout, '(a10, 10f8.3)')'theta:', angle_t((i-1)*10+1:i*10)
         write(stdout, '(a10, 10i8)')'No. atoms:', number_of_atoms((i-1)*10+1:i*10)
      write(stdout, '(a8)') ' '
      enddo
      write(stdout, '(a)') '>>Index-angle table<<'
      write(stdout, '(a)') ' '

      write(stdout, '(a, i4)'), 'm= ', m
      write(stdout, '(a, f9.2)'), 'twisted angle', theta*180/pi

      !> 20 angstrom is the thickness of vacuum.
      lc_z= interlayer_distance*(number_layers-1d0)+20d0
      !> primitive lattice vectors of Moire supercell
     !t1= m*a1 + (m+1)*a2;
      t1=-(m+1)*a1 + (2*m+1)*a2;
      t2= (2*m+1)*a1 - m*a2
      t3= (/0d0, 0d0, lc_z/)

      write(stdout, '(a, f8.3, a)')'length of lattice vector:', dsqrt(sum(t1**2)), ' Angstrom'

      !> calculated using the area ratio between moire supercell and Graphene primitive cell
      num_atom_per_layer= 2*(3*m*m+3*m+1)
      num_atom_total= num_atom_per_layer* number_layers
      write(stdout, '(a, i10)') 'number of atoms per cell', num_atom_total

      allocate(atoms_morie_cell_direct(3, num_atom_total))
      atoms_morie_cell_direct= 0d0

      !> build a supercell contains enough atoms
      if (m<2) then
         max_period= 9*m+4
      else
         max_period=3*m+4
      endif
      max_atoms_perlayer=max_period**2*2;
      allocate(atoms_supercell_temp0(3, max_atoms_perlayer))
      allocate(atoms_supercell_temp(3, max_atoms_perlayer*number_layers))
      atoms_supercell_temp0= 0d0
      atoms_supercell_temp= 0d0
      A_site= 0d0
      B_site= 0d0

      !> build a temporary supercell
      iter=0
      do i=1, max_period
         do j=1, max_period
            !> A site
            A_site= (i-max_period/2)*a1+ (j-max_period/2)*a2
           !A_site= (i-1)*a1+ (j-1)*a2
            iter=iter+1
            atoms_supercell_temp0(:, iter)=A_site
            B_site(2)=A_site(2)
            B_site(1)=A_site(1)+lattice_constant/sqrt_3
            iter=iter+1
            atoms_supercell_temp0(:, iter)=B_site
         enddo
      enddo

      !> write out 
     !call write_xyz(max_atoms_perlayer, atoms_supercell_temp, 'layer1.xyz')

      !> Let's find all the atoms inside the Moire supercell
      do ilayer=1, number_layers
         angle=theta*twisted_angle_array(ilayer)
         select case (stacking_sequences(ilayer))
            case ('A')
               tau= 0d0
            case ('B')
               tau(1)=-lattice_constant/sqrt_3; tau(2)=0d0; tau(3) = 0d0
            case ('C')
               tau(1)= lattice_constant/sqrt_3; tau(2)=0d0; tau(3) = 0d0
            case default
               tau= 0d0
         end select

         tau(3)= (ilayer-1d0)*interlayer_distance+lc_z/2d0-interlayer_distance*(number_layers-1d0)/2d0
         do i=1, max_atoms_perlayer
            pos1= atoms_supercell_temp0(:, i)+tau
            call rotate_z(angle, pos1, pos_cart)
            call cart_direct_real(t1, t2, t3, pos_cart, pos_direct)
           !pos_direct= mod(pos_direct, 1d0)
            call transformtohomecell(pos_direct)
            atoms_supercell_temp(:, i+(ilayer-1)*max_atoms_perlayer)=pos_direct
         enddo
        !write(xyzfilename, '(a,i1,a)')'layer',ilayer,'.xyz'
        !call write_xyz(max_atoms_perlayer, atoms_supercell_temp, xyzfilename)
      enddo
      call eliminate_duplicates(max_atoms_perlayer*number_layers, atoms_supercell_temp, Nleft)
      do i=1, Nleft
      !  write(stdout, '(3f12.6)')atoms_supercell_temp(1:3, i)
      enddo
      write(stdout, *) Nleft, num_atom_per_layer*number_layers
      call writeout_poscar(t1, t2, t3, Nleft, atoms_supercell_temp(:,1:Nleft))

      if (Nleft.ne.num_atom_per_layer*number_layers) then
         write(stdout, *)"Error: There is something wrong when generating the Morie cell"
         write(stdout, *) "Nleft, num_atom_per_layer*number_layers, max_atoms_perlayer*number_layers"
         write(stdout, *) Nleft, num_atom_per_layer*number_layers, max_atoms_perlayer*number_layers
         stop
      endif
      atoms_morie_cell_direct(:, :)= atoms_supercell_temp(:, 1:num_atom_total)


      !> write out crystal structure in LAMMPS format
      !> The three lattice vectors should be a=(xhi-xlo,0,0); b=(xy,yhi-ylo,0); c=(xz,yz,zhi-zlo)
      !> set new axis so that the first vector is along x direction
      !> define a rotation matrix
      !> ex'
      Urot(1, :)=t1/norm(t1)

      !> ez'
      call cross_product(t1, t2, Urot(3,:))
      Urot(3, :)=Urot(3, :)/norm(Urot(3, :))

      !> ey'= ez'\cross ex'
      call cross_product(Urot(3, :), Urot(1, :), Urot(2,:))
      Urot(2, :)=Urot(2, :)/norm(Urot(2, :))

      R1= Urot(:, 1)*t1(1)+ Urot(:, 2)*t1(2)+ Urot(:, 3)*t1(3)
      R2= Urot(:, 1)*t2(1)+ Urot(:, 2)*t2(2)+ Urot(:, 3)*t2(3)
      R3= Urot(:, 1)*t3(1)+ Urot(:, 2)*t3(2)+ Urot(:, 3)*t3(3)


      xlo= 0d0
      ylo= 0d0
      zlo= -40d0
      xhi=R1(1)
      xy =R2(1); yhi=R2(2)
      xz =R3(1); yz =R3(2); zhi=60d0

      !open file
      open(unit=1012, file='rigid_structure.lammps')
      !> a comment line
      write(1012, '(a, i3, a, a, f10.6, a)')'Twisted', number_layers, ' layers', &
         'theta= ', twisted_angle_degree, ' degree.'
      write(1012, '(i12, a)')num_atom_total, ' atoms'
      write(1012, '(a)')' '
      write(1012, '(i4, a)')number_layers, '  atom types'  
      write(1012, '(a)')' '
      write(1012, '(2f14.8, 2a4)')xlo, xhi, 'xlo', 'xhi'
      write(1012, '(2f14.8, 2a4)')ylo, yhi, 'ylo', 'yhi'
      write(1012, '(2f14.8, 2a4)')zlo, zhi, 'zlo', 'zhi'
      write(1012, '(3f14.8, a)')xy, xz, yz, ' xy xz yz'
      write(1012, '(a)')' '
      write(1012, '(a)')' Atoms'
      write(1012, '(a)')' '
      iter=0
      do i=1, number_layers
         do j=1, num_atom_per_layer
            iter=iter+1
            pos_direct= atoms_morie_cell_direct(:, iter)
            call direct_cart_real(t1, t2, t3, pos_direct, pos_cart)
            pos1= Urot(:, 1)*pos_cart(1)+ Urot(:, 2)*pos_cart(2)+ Urot(:, 3)*pos_cart(3)
            write(1012, '(2I5, 3f20.14)')iter, i, pos1
         enddo
      enddo
      close(1012)

      outfileindex= outfileindex+1
      open(outfileindex, file='POSCAR-rigid')
      write(outfileindex, '(a, i3, a, a, f10.6, a)')'! Twisted', number_layers, ' layers', &
      'theta= ', twisted_angle_degree, ' degree.'
      write(outfileindex, '(a)')"1.0"
      write(outfileindex, '(3f12.6)') R1
      write(outfileindex, '(3f12.6)') R2
      write(outfileindex, '(3f12.6)') R3
      write(outfileindex, '(30A6)') 'C'
      write(outfileindex, '(30i6)') num_atom_total
      write(outfileindex, '(a)')"Cartesian"
      do i=1, num_atom_total
         pos_direct= atoms_morie_cell_direct(:, i)
         call direct_cart_real(t1, t2, t3, pos_direct, pos_cart)
         pos1= Urot(:, 1)*pos_cart(1)+ Urot(:, 2)*pos_cart(2)+ Urot(:, 3)*pos_cart(3)
         write(outfileindex, '(3f12.6, a9)')pos1, 'C'
      enddo
      close(outfileindex)


      return
   end subroutine generate_crystal_structure

   subroutine write_xyz(num_atoms, atoms_position_cart, xyz_file_name)
      !> write out the positions of atoms into a xyz file
      !> <number of atoms>
      !> comment line
      !> <element> <X> <Y> <Z>
      !> <element> <X> <Y> <Z>
      !> ...

      use para
      implicit none
      integer, intent(in) :: num_atoms
      real(dp), intent(in) :: atoms_position_cart(3, num_atoms)
      character(*), intent(in) :: xyz_file_name

      integer :: ia

      open(unit=2423, file=xyz_file_name)

      write(2423, '(i10)') num_atoms
      write(2423, '(a)') 'All atoms are Carbon atoms'
      do ia=1, num_atoms
         write(2423, '(a3, 3f12.6)') 'C', atoms_position_cart(:, ia)
      enddo
      close(2423)

      return
   end subroutine write_xyz

   subroutine rotate_z(theta, R1, R2)
      !> Theta should be in SI unit. radian
      !> Rotate R1 vector by theta with respect to the z axis
      use para, only : dp
      implicit none

      real(dp), intent(in) :: theta 
      real(dp), intent(in) :: R1(3)
      real(dp), intent(out) :: R2(3)

      real(dp) :: rotation_matrix(3,3)

      rotation_matrix= 0d0
      rotation_matrix(1, :)= (/cos(theta), -sin(theta), 0d0/)
      rotation_matrix(2, :)= (/sin(theta),  cos(theta), 0d0/)
      rotation_matrix(3, :)= (/0d0       ,  0d0       , 1d0/)

      R2(1)=R1(1)*rotation_matrix(1, 1)+R1(2)*rotation_matrix(2, 1)
      R2(2)=R1(1)*rotation_matrix(1, 2)+R1(2)*rotation_matrix(2, 2)
      R2(3)=R1(3)

      return
   end subroutine rotate_z

   subroutine cart_direct_real(t1, t2, t3, pos_cart, pos_direct)
      !> transform from Cartesian coordinates to lattice vector basis
      !> in this subroutine, we only need to get the first two coordinates and discard the z axis.
      use para
      implicit none
      real(dp), intent(in) :: pos_cart(3)
      real(dp), intent(in) :: t1(3)
      real(dp), intent(in) :: t2(3)
      real(dp), intent(in) :: t3(3)
      real(dp), intent(inout) :: pos_direct(3)
   
      real(dp), allocatable :: mat(:, :)
      real(dp) :: denominator
   
      allocate(mat(2, 2))
      mat= 0d0
      pos_direct= 0d0
  
      denominator= 1d0/(t1(1)*t2(2)-t1(2)*t2(1))
      mat(1, 1)= t2(2)*denominator
      mat(1, 2)=-t1(2)*denominator
      mat(2, 1)=-t2(1)*denominator
      mat(2, 2)= t1(1)*denominator
   
      pos_direct(1:2)= pos_cart(1)*mat(1, :)+ pos_cart(2)*mat(2, :) 
      pos_direct(3)= pos_cart(3)/t3(3)
   
      deallocate(mat)
   
      return
   end subroutine cart_direct_real

   subroutine direct_cart_real(t1, t2, t3, pos_direct, pos_cart)
      !> transform from Cartesian coordinates to lattice vector basis
      !> in this subroutine, we only need to get the first two coordinates and discard the z axis.
      use para
      implicit none
      real(dp), intent(in) :: t1(3)
      real(dp), intent(in) :: t2(3)
      real(dp), intent(in) :: t3(3)
      real(dp), intent(inout) :: pos_cart(3)
      real(dp), intent(in) :: pos_direct(3)
   
      pos_cart= pos_direct(1)*t1+ pos_direct(2)*t2+pos_direct(3)*t3
   
      return
   end subroutine direct_cart_real


   subroutine eliminate_duplicates(isize, array, Nleft)
      ! Eliminate the duplicated entries.
      ! On output, array will be replaced.
      !
      ! By QuanSheng Wu 
      !
      ! wuquansheng@gmail.com
      !
      ! June 28 2020, modified from WannierTools.

      use para, only : dp, eps6
      integer, intent(in) :: isize
      integer, intent(out) :: Nleft
      real(dp), intent(inout) :: array(3, isize)
      real(dp), allocatable :: array_left(:, :)

      integer :: it, ik, ik1
      logical :: Logical_duplicate

      allocate(array_left(3, isize))
      array_left=0
      array_left(:, 1)= array(:, 1)

      Nleft = 1
      it= 1
      do ik=2, isize
         Logical_duplicate= .False.
         do ik1=1, Nleft
            if (sum(abs(array(:, ik)-array_left(:, ik1)))<eps6) then
               Logical_duplicate= .True.
               exit
            endif
         enddo
         if (.not.Logical_duplicate)then
            Nleft= Nleft+ 1
            array_left(:, Nleft)= array(:, ik)
         endif
      enddo

      array(:, 1:Nleft)= array_left(:, 1:Nleft)

      deallocate(array_left)

      return
   end subroutine eliminate_duplicates


   !> Write out the POSCAR for a given cell
   subroutine writeout_poscar(t1, t2, t3, num_atoms, atom_position_direct)
      use para, only : dp, outfileindex
      implicit none
   
      integer :: ia
      real(dp), intent(in) :: t1(3)
      real(dp), intent(in) :: t2(3)
      real(dp), intent(in) :: t3(3)
      integer, intent(in) :: num_atoms
      real(dp), intent(in) :: atom_position_direct(3, num_atoms)
      character(6) :: poscarname
  
      poscarname= 'POSCAR'
      !> print out the new basis
      outfileindex= outfileindex+ 1
      open(outfileindex, file=poscarname)
      write(outfileindex, '(a)')"POSCAR for twisted bilayer graphene"
      write(outfileindex, '(a)')"1.0"
      write(outfileindex, '(3f12.6)') t1
      write(outfileindex, '(3f12.6)') t2
      write(outfileindex, '(3f12.6)') t3
      write(outfileindex, '(30A6)') 'C'
      write(outfileindex, '(30i6)') num_atoms
      write(outfileindex, '(a)')"Direct"
      do ia=1, Num_atoms
         write(outfileindex, '(3f12.6, a9)')Atom_position_direct(:, ia)
      enddo
      close(outfileindex)
      return
   end subroutine writeout_poscar



    !> move the atoms into the home unitcell [0, 1)*[0, 1)*[0, 1)
   subroutine transformtohomecell(pos)
      ! Transform the k points to the 1st BZ
      !
      ! By QuanSheng Wu
      !
      ! wuquansheng@gmail.com
      !
      ! Nov 9 2016  at ETHZ
   
      use para, only : dp, pi
   
      integer :: i
      real(dp), intent(inout) :: pos(3)
      real(dp) :: irrational_shift(3)
      irrational_shift= (/ pi/1000d0, pi/1000d0, 0d0 /)
      pos= pos+ irrational_shift
   
      do i=1, 3
         do while (.true.)
            if (pos(i)>= 0.0d0 .and. pos(i)<1.0d0) then
               exit
            else if (pos(i)< 0.0d0) then
               pos(i)= pos(i)+ 1d0
            else if (pos(i)>=1.0d0) then
               pos(i)= pos(i)- 1d0
            endif
         enddo
      enddo
      pos= pos- irrational_shift
   
      return
   end subroutine transformtohomecell



!> cross product of two 3-value vectors
!> R3=R1 x R2
subroutine cross_product(R1, R2, R3)
   use para, only : dp
   implicit none

   real(dp), intent(in) :: R1(3), R2(3)
   real(dp), intent(out) :: R3(3)

   R3(1)= R1(2)*R2(3)- R1(3)*R2(2)
   R3(2)= R1(3)*R2(1)- R1(1)*R2(3)
   R3(3)= R1(1)*R2(2)- R1(2)*R2(1)

   return
end subroutine cross_product


!> Get the length of a given 3-value vector
function norm(R1)
   use para, only : dp

   implicit none
   real(dp), intent(in) :: R1(3)
   real(dp) :: norm1
   real(dp) :: norm

   norm1= R1(1)*R1(1)+ R1(2)*R1(2)+ R1(3)*R1(3)
   norm= sqrt(norm1)

   return
end function norm


!> Get the angle between two given vectors which are in cartesian coordinates
!> return in degree
function get_angle(R1, R2)
   use para, only : dp, pi

   implicit none
   real(dp), intent(in) :: R1(3), R2(3)
   real(dp) :: get_angle, adotb
   real(dp), external :: norm
 
   adotb=dot_product(R1, R2)
   get_angle= acos(adotb/norm(R1)/norm(R2))*180/pi

   return
end function get_angle


