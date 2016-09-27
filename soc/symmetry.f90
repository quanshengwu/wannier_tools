!> set symmetry 
  subroutine symmetry
     use para
     implicit none

     integer :: nwan
     integer :: ia, i, n

     !> get the atom afterr the mirror_x operation
     integer, allocatable :: iatom_mirror_x(:)
 
     !> set up symmetry operators
     !> here we assume that the wannier functions have the symmetry 
     !> as atomic orbitals
    
     nwan= Num_wann/2

     allocate(iatom_mirror_x(Num_atoms))
     allocate(inversion(Num_wann, Num_wann))
     allocate(mirror_z(Num_wann, Num_wann))
     allocate(mirror_x(Num_wann, Num_wann))
     inversion= 0d0
     mirror_z = 0d0
     mirror_x = 0d0
     iatom_mirror_x= 0

     !> inversion symmetry
     !> s-> s; p-> -p; d-> -d
     n= 0
     do ia=1, Num_atoms
        do i=1, nprojs(ia)
           n= n+ 1
           select case (proj_name(i, ia))
           case ('s', 'S')
              inversion(n, n)= 1
              inversion(n+ nwan, n+ nwan)= 1
           case ('px', 'Px', 'PX')
              inversion(n, n)=-1
              inversion(n+ nwan, n+ nwan)=-1
           case ('py', 'Py', 'PY')
              inversion(n, n)=-1
              inversion(n+ nwan, n+ nwan)=-1
           case ('pz', 'Pz', 'PZ')
              inversion(n, n)=-1
              inversion(n+ nwan, n+ nwan)=-1
           case ('dxy', 'Dxy', 'DXY')
              inversion(n, n)= 1
              inversion(n+ nwan, n+ nwan)= 1
           case ('dyz', 'Dyz', 'DYZ')
              inversion(n, n)= 1
              inversion(n+ nwan, n+ nwan)= 1
           case ('dxz', 'Dxz', 'DXZ', 'dzx', 'Dzx', 'DZX')
              inversion(n, n)= 1
              inversion(n+ nwan, n+ nwan)= 1
           case ('dx2-y2', 'Dx2-y2', 'DX2-Y2', 'dx2', 'DX2')
              inversion(n, n)= 1
              inversion(n+ nwan, n+ nwan)= 1
           case ('dz2', 'Dz2', 'DZ2')
              inversion(n, n)= 1
              inversion(n+ nwan, n+ nwan)= 1
           case default
              write(*, *) "ERROR: only support s px py pz dxy dyz dxz dx2-y2 dz2 orbitals"
              stop
           end select
        enddo ! i
     enddo ! ia

     !> mirror_x symmetry
     !> s-> s; px->-px, py->py, pz->  pz
     !> dxy-> -dxy, dyz->  dyz, dxz-> -dxz, dx2-> dx2 dz2->dz2 
     !> up-> -i dn  dn-> -i up
     !> here we throw away the phase -i, it is just a constant, leading M*M= -1

     n= 0
     mirror_x= 0d0
     do ia=1, Num_atoms
        do i=1, nprojs(ia)
           n= n+ 1
           select case (proj_name(i, ia))
           case ('s', 'S')
              mirror_x(n, n+ nwan)= 1d0
              mirror_x(n+ nwan, n)= 1d0
           case ('px', 'Px', 'PX')
              mirror_x(n, n+ nwan)=-1d0
              mirror_x(n+ nwan, n)=-1d0
           case ('py', 'Py', 'PY')
              mirror_x(n, n+ nwan)= 1d0
              mirror_x(n+ nwan, n)= 1d0
           case ('pz', 'Pz', 'PZ')
              mirror_x(n, n+ nwan)= 1d0
              mirror_x(n+ nwan, n)= 1d0
           case ('dxy', 'Dxy', 'DXY')
              mirror_x(n, n+ nwan)=-1d0
              mirror_x(n+ nwan, n)=-1d0
           case ('dyz', 'Dyz', 'DYZ')
              mirror_x(n, n+ nwan)= 1d0
              mirror_x(n+ nwan, n)= 1d0
           case ('dxz', 'Dxz', 'DXZ', 'dzx', 'Dzx', 'DZX')
              mirror_x(n, n+ nwan)=-1d0
              mirror_x(n+ nwan, n)=-1d0
           case ('dx2-', 'dx2-y2', 'Dx2-y2', 'DX2-Y2', 'dx2', 'DX2')
              mirror_x(n, n+ nwan)= 1d0
              mirror_x(n+ nwan, n)= 1d0
           case ('dz2', 'Dz2', 'DZ2')
              mirror_x(n, n+ nwan)= 1d0
              mirror_x(n+ nwan, n)= 1d0
           case default
              write(*, *) "ERROR: only support s px py pz dxy dyz dxz dx2-y2 dz2 orbitals"
              stop
           end select
        enddo ! i
     enddo ! ia

     do i=1, Num_wann
       !write(*, '(1000i2)')int(real(mirror_x(:, i)))
     enddo


     !> mirror_z symmetry
     !> s-> s; px->px, py->py, pz-> -pz
     !> dxy-> dxy, dyz-> -dyz, dxz-> -dxz, dx2-> dx2 dz2->dz2 
     !> up-> up  dn-> -dn
     n= 0
     do ia=1, Num_atoms
        do i=1, nprojs(ia)
           n= n+ 1
           select case (proj_name(i, ia))
           case ('s', 'S')
              mirror_z(n, n)= 1
              mirror_z(n+ nwan, n+ nwan)=-1
           case ('px', 'Px', 'PX')
              mirror_z(n, n)= 1
              mirror_z(n+ nwan, n+ nwan)=-1
           case ('py', 'Py', 'PY')
              mirror_z(n, n)= 1
              mirror_z(n+ nwan, n+ nwan)=-1
           case ('pz', 'Pz', 'PZ')
              mirror_z(n, n)=-1
              mirror_z(n+ nwan, n+ nwan)= 1
           case ('dxy', 'Dxy', 'DXY')
              mirror_z(n, n)= 1
              mirror_z(n+ nwan, n+ nwan)=-1
           case ('dyz', 'Dyz', 'DYZ')

              mirror_z(n, n)=-1
              mirror_z(n+ nwan, n+ nwan)= 1
           case ('dxz', 'Dxz', 'DXZ', 'dzx', 'Dzx', 'DZX')
              mirror_z(n, n)=-1
              mirror_z(n+ nwan, n+ nwan)= 1
           case ('dx2-', 'dx2-y2', 'Dx2-y2', 'DX2-Y2', 'dx2', 'DX2')
              mirror_z(n, n)= 1
              mirror_z(n+ nwan, n+ nwan)=-1
           case ('dz2', 'Dz2', 'DZ2')
              mirror_z(n, n)= 1
              mirror_z(n+ nwan, n+ nwan)=-1
           case default
              write(*, *) "ERROR: only support s px py pz dxy dyz dxz dx2-y2 dz2 orbitals"
              stop
           end select
        enddo ! i
     enddo ! ia


     !> set up symmetry operators
     allocate(mirror_x_op(3,3))
     allocate(mirror_y_op(3,3))
     allocate(mirror_z_op(3,3))

     !> for glide symmetry, (1:3, 1:3) shows the mirror operation, (1:3, 4) 
     !> gives the shift
     allocate(glide_y_op(3,4))

     mirror_x_op= 0d0
     mirror_y_op= 0d0
     mirror_z_op= 0d0
     glide_y_op= 0d0

     mirror_x_op(1, 1)=-1d0
     mirror_x_op(2, 2)= 1d0
     mirror_x_op(3, 3)= 1d0

     mirror_y_op(1, 1)= 1d0
     mirror_y_op(2, 2)=-1d0
     mirror_y_op(3, 3)= 1d0

     mirror_z_op(1, 1)= 1d0
     mirror_z_op(2, 2)= 1d0
     mirror_z_op(3, 3)=-1d0

     glide_y_op(1, 1) = 1d0
     glide_y_op(2, 2) =-1d0
     glide_y_op(3, 3) = 1d0
     glide_y_op(:, 4) = (/0.5d0, 0.0d0, 0.5d0/)

     return
  end subroutine symmetry
