!> set symmetry 
!> This subroutine is a temporary solution. 
!> It can be used in the mirror chern number calculation
!> In the future, we will added the automatically detection of the groups
!> and automatically generate the matrix representations for each operators 
!> in the basis of atomic like Wannier functions. 
!> If you want to study the mirror Chern number, you have to define the
!> operator matrix yourself. For example, mirror_x, mirror_z.
!> In the mirror Chern number calculation subroutine, we now only support 
!> very simple mirror symmetries mirror_x and mirror_z. So you have to rotate
!> your coordinates so that the mirror plane you concern is mirror_z or mirror_x
!> By default we use mirror_z. So if you want to use mirror_x, please modify 
!> wanniercenter.f90
!> If you have problems with that, please contact with wuquansheng@gmail.com
!> Sep/01/2018 EPFL Switzerland
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
   
     if (SOC<=0) return
     nwan= Num_wann/2

     allocate(iatom_mirror_x(Origin_cell%Num_atoms))
     allocate(inversion(Num_wann, Num_wann))
     allocate(mirror_z(Num_wann, Num_wann))
     allocate(mirror_y(Num_wann, Num_wann))
     allocate(mirror_x(Num_wann, Num_wann))
     allocate(C2yT(Num_wann, Num_wann))
     inversion= 0d0
     mirror_z = 0d0
     mirror_y = 0d0
     mirror_x = 0d0
     C2yT = 0d0
     iatom_mirror_x= 0

     !> symmetry operators for VASP
     !> for VASP, the orbital order is up up  up up  dn dn dn dn 
     if (index( Package, 'VASP')/=0.or. index( Package, 'Wien2k')/=0 &
        .or. index( Package, 'Abinit')/=0.or. index( Package, 'openmx')/=0) then
        !> inversion symmetry
        !> s-> s; p-> -p; d-> -d
        n= 0
        do ia=1, Origin_cell%Num_atoms
           do i=1, Origin_cell%nprojs(ia)
              n= n+ 1
              select case (Origin_cell%proj_name(i, ia))
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
        do ia=1, Origin_cell%Num_atoms
           do i=1, Origin_cell%nprojs(ia)
              n= n+ 1
              select case (Origin_cell%proj_name(i, ia))
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
 
        !> C2yT symmetry =  C_2y * i*sigma_y 
        !> s-> s; px->-px, py-> py, pz-> -pz
        !> dxy-> -dxy, dyz-> -dyz, dxz-> dxz, dx2-> dx2 dz2->dz2 
        !> up-> up  dn-> dn
        !> Here we keep the phase i of T, we also drop the Conjugation of T
      
        n= 0
        C2yT= 0d0
        do ia=1, Origin_cell%Num_atoms
           do i=1, Origin_cell%nprojs(ia)
              n= n+ 1
              select case (Origin_cell%proj_name(i, ia))
              case ('s', 'S')
                 C2yT(n, n)= 1d0
                 C2yT(n+ nwan, n+ nwan)= 1d0
              case ('px', 'Px', 'PX')
                 C2yT(n, n)=-1 
                 C2yT(n+ nwan, n+ nwan)=-1 
              case ('py', 'Py', 'PY')
                 C2yT(n, n)= 1 
                 C2yT(n+ nwan, n+ nwan)= 1 
              case ('pz', 'Pz', 'PZ')
                 C2yT(n, n)=-1 
                 C2yT(n+ nwan, n+ nwan)=-1 
              case ('dxy', 'Dxy', 'DXY')
                 C2yT(n, n)=-1 
                 C2yT(n+ nwan, n+ nwan)=-1 
              case ('dyz', 'Dyz', 'DYZ')
                 C2yT(n, n)=-1 
                 C2yT(n+ nwan, n+ nwan)=-1 
              case ('dxz', 'Dxz', 'DXZ', 'dzx', 'Dzx', 'DZX')
                 C2yT(n, n)= 1 
                 C2yT(n+ nwan, n+ nwan)= 1 
              case ('dx2-', 'dx2-y2', 'Dx2-y2', 'DX2-Y2', 'dx2', 'DX2')
                 C2yT(n, n)= 1 
                 C2yT(n+ nwan, n+ nwan)= 1 
              case ('dz2', 'Dz2', 'DZ2')
                 C2yT(n, n)= 1 
                 C2yT(n+ nwan, n+ nwan)= 1 
              case default
                 write(*, *) "ERROR: only support s px py pz dxy dyz dxz dx2-y2 dz2 orbitals"
                 stop
              end select
           enddo ! i
        enddo ! ia
      
     

        !> mirror_y symmetry = i*sigma_y*R_y
        !> s-> s; px->px, py->-py, pz->  pz
        !> dxy-> -dxy, dyz-> -dyz, dxz->  dxz, dx2-> dx2 dz2->dz2 
        !> up-> -zi*dn  dn-> zi*up
        !> Drop off phase i
      
        n= 0
        mirror_y= 0d0
        do ia=1, Origin_cell%Num_atoms
           do i=1, Origin_cell%nprojs(ia)
              n= n+ 1
              select case (Origin_cell%proj_name(i, ia))
              case ('s', 'S')
                 mirror_y(n, n+ nwan)=-zi
                 mirror_y(n+ nwan, n)= zi
              case ('px', 'Px', 'PX')
                 mirror_y(n, n+ nwan)=-zi
                 mirror_y(n+ nwan, n)= zi
              case ('py', 'Py', 'PY')
                 mirror_y(n, n+ nwan)= zi
                 mirror_y(n+ nwan, n)=-zi
              case ('pz', 'Pz', 'PZ')
                 mirror_y(n, n+ nwan)=-zi
                 mirror_y(n+ nwan, n)= zi
              case ('dxy', 'Dxy', 'DXY')
                 mirror_y(n, n+ nwan)= zi
                 mirror_y(n+ nwan, n)=-zi
              case ('dyz', 'Dyz', 'DYZ')
                 mirror_y(n, n+ nwan)= zi
                 mirror_y(n+ nwan, n)=-zi
              case ('dxz', 'Dxz', 'DXZ', 'dzx', 'Dzx', 'DZX')
                 mirror_y(n, n+ nwan)=-zi
                 mirror_y(n+ nwan, n)= zi
              case ('dx2-', 'dx2-y2', 'Dx2-y2', 'DX2-Y2', 'dx2', 'DX2')
                 mirror_y(n, n+ nwan)=-zi
                 mirror_y(n+ nwan, n)= zi
              case ('dz2', 'Dz2', 'DZ2')
                 mirror_y(n, n+ nwan)=-zi
                 mirror_y(n+ nwan, n)= zi
              case default
                 write(*, *) "ERROR: only support s px py pz dxy dyz dxz dx2-y2 dz2 orbitals"
                 stop
              end select
           enddo ! i
        enddo ! ia
      
        do i=1, Num_wann
          !write(*, '(1000i2)')int(real(mirror_x(:, i)))
        enddo
      
      
        !> mirror_z symmetry  i*sigma_z*R_z
        !> s-> s; px->px, py->py, pz-> -pz
        !> dxy-> dxy, dyz-> -dyz, dxz-> -dxz, dx2-> dx2 dz2->dz2 
        !> up-> up  dn-> -dn  Drop off phase i
        n= 0
        do ia=1, Origin_cell%Num_atoms
           do i=1, Origin_cell%nprojs(ia)
              n= n+ 1
              select case (Origin_cell%proj_name(i, ia))
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
 
     !> for QE, the orbital order is up dn  up dn  up dn up dn 
     elseif (index( Package, 'QE')/=0.or.index( Package, 'quantumespresso')/=0 &
         .or.index( Package, 'quantum-espresso')/=0.or.index( Package, 'pwscf')/=0) then
        !> inversion symmetry
        !> s-> s; p-> -p; d-> -d
        n= 0
        do ia=1, Origin_cell%Num_atoms
           do i=1, Origin_cell%nprojs(ia)
              n= n+ 1
              select case (Origin_cell%proj_name(i, ia))
              case ('s', 'S')
                 inversion(2*n-1, 2*n-1)= 1
                 inversion(2*n, 2*n)= 1
              case ('px', 'Px', 'PX')
                 inversion(2*n-1, 2*n-1)=-1
                 inversion(2*n, 2*n)=-1
              case ('py', 'Py', 'PY')
                 inversion(2*n-1, 2*n-1)=-1
                 inversion(2*n, 2*n)=-1
              case ('pz', 'Pz', 'PZ')
                 inversion(2*n-1, 2*n-1)=-1
                 inversion(2*n, 2*n)=-1
              case ('dxy', 'Dxy', 'DXY')
                 inversion(2*n-1, 2*n-1)= 1
                 inversion(2*n, 2*n)= 1
              case ('dyz', 'Dyz', 'DYZ')
                 inversion(2*n-1, 2*n-1)= 1
                 inversion(2*n, 2*n)= 1
              case ('dxz', 'Dxz', 'DXZ', 'dzx', 'Dzx', 'DZX')
                 inversion(2*n-1, 2*n-1)= 1
                 inversion(2*n, 2*n)= 1
              case ('dx2-y2', 'Dx2-y2', 'DX2-Y2', 'dx2', 'DX2')
                 inversion(2*n-1, 2*n-1)= 1
                 inversion(2*n, 2*n)= 1
              case ('dz2', 'Dz2', 'DZ2')
                 inversion(2*n-1, 2*n-1)= 1
                 inversion(2*n, 2*n)= 1
              case default
                 write(*, *) "ERROR: only support s px py pz dxy dyz dxz dx2-y2 dz2 orbitals"
                 stop
              end select
           enddo ! i
        enddo ! ia
       
        !> mirror_x symmetry i*sigma_x*R_x
        !> s-> s; px->-px, py->py, pz->  pz
        !> dxy-> -dxy, dyz->  dyz, dxz-> -dxz, dx2-> dx2 dz2->dz2 
        !> up-> i dn  dn-> i up
        !> here we throw away the phase i, it is just a constant, leading M*M= -1
      
        n= 0
        mirror_x= 0d0
        do ia=1, Origin_cell%Num_atoms
           do i=1, Origin_cell%nprojs(ia)
              n= n+ 1
              select case (Origin_cell%proj_name(i, ia))
              case ('s', 'S')
                 mirror_x(2*n-1, 2*n)= 1d0
                 mirror_x(2*n, 2*n-1)= 1d0
              case ('px', 'Px', 'PX')
                 mirror_x(2*n-1, 2*n)=-1d0
                 mirror_x(2*n, 2*n-1)=-1d0
              case ('py', 'Py', 'PY')
                 mirror_x(2*n-1, 2*n)= 1d0
                 mirror_x(2*n, 2*n-1)= 1d0
              case ('pz', 'Pz', 'PZ')
                 mirror_x(2*n-1, 2*n)= 1d0
                 mirror_x(2*n, 2*n-1)= 1d0
              case ('dxy', 'Dxy', 'DXY')
                 mirror_x(2*n-1, 2*n)=-1d0
                 mirror_x(2*n, 2*n-1)=-1d0
              case ('dyz', 'Dyz', 'DYZ')
                 mirror_x(2*n-1, 2*n)= 1d0
                 mirror_x(2*n, 2*n-1)= 1d0
              case ('dxz', 'Dxz', 'DXZ', 'dzx', 'Dzx', 'DZX')
                 mirror_x(2*n-1, 2*n)=-1d0
                 mirror_x(2*n, 2*n-1)=-1d0
              case ('dx2-', 'dx2-y2', 'Dx2-y2', 'DX2-Y2', 'dx2', 'DX2')
                 mirror_x(2*n-1, 2*n)= 1d0
                 mirror_x(2*n, 2*n-1)= 1d0
              case ('dz2', 'Dz2', 'DZ2')
                 mirror_x(2*n-1, 2*n)= 1d0
                 mirror_x(2*n, 2*n-1)= 1d0
              case default
                 write(*, *) "ERROR: only support s px py pz dxy dyz dxz dx2-y2 dz2 orbitals"
                 stop
              end select
           enddo ! i
        enddo ! ia
           
        !> mirror_y symmetry = i*sigma_y*R_z
        !> s-> s; px->px, py->-py, pz->  pz
        !> dxy-> -dxy, dyz-> -dyz, dxz-> dxz, dx2-> dx2 dz2->dz2 
        !> up-> -zi*dn  dn-> zi*up
        !> Drop off phase i
      
        n= 0
        mirror_y= 0d0
        do ia=1, Origin_cell%Num_atoms
           do i=1, Origin_cell%nprojs(ia)
              n= n+ 1
              select case (Origin_cell%proj_name(i, ia))
              case ('s', 'S')
                 mirror_y(2*n-1, 2*n)=-zi
                 mirror_y(2*n, 2*n-1)= zi
              case ('px', 'Px', 'PX')
                 mirror_y(2*n-1, 2*n)=-zi
                 mirror_y(2*n, 2*n-1)= zi
              case ('py', 'Py', 'PY')
                 mirror_y(2*n-1, 2*n)= zi
                 mirror_y(2*n, 2*n-1)=-zi
              case ('pz', 'Pz', 'PZ')
                 mirror_y(2*n-1, 2*n)=-zi
                 mirror_y(2*n, 2*n-1)= zi
              case ('dxy', 'Dxy', 'DXY')
                 mirror_y(2*n-1, 2*n)= zi
                 mirror_y(2*n, 2*n-1)=-zi
              case ('dyz', 'Dyz', 'DYZ')
                 mirror_y(2*n-1, 2*n)= zi
                 mirror_y(2*n, 2*n-1)=-zi
              case ('dxz', 'Dxz', 'DXZ', 'dzx', 'Dzx', 'DZX')
                 mirror_y(2*n-1, 2*n)=-zi
                 mirror_y(2*n, 2*n-1)= zi
              case ('dx2-', 'dx2-y2', 'Dx2-y2', 'DX2-Y2', 'dx2', 'DX2')
                 mirror_y(2*n-1, 2*n)=-zi
                 mirror_y(2*n, 2*n-1)= zi
              case ('dz2', 'Dz2', 'DZ2')
                 mirror_y(2*n-1, 2*n)=-zi
                 mirror_y(2*n, 2*n-1)= zi
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
        do ia=1, Origin_cell%Num_atoms
           do i=1, Origin_cell%nprojs(ia)
              n= n+ 1
              select case (Origin_cell%proj_name(i, ia))
              case ('s', 'S')
                 mirror_z(2*n-1, 2*n-1)= 1d0
                 mirror_z(2*n, 2*n)= -1d0
              case ('px', 'Px', 'PX')
                 mirror_z(2*n-1, 2*n-1)= 1d0
                 mirror_z(2*n, 2*n)= -1d0
              case ('py', 'Py', 'PY')
                 mirror_z(2*n-1, 2*n-1)= 1d0
                 mirror_z(2*n, 2*n)= -1d0
              case ('pz', 'Pz', 'PZ')
                 mirror_z(2*n-1, 2*n-1)=-1d0
                 mirror_z(2*n, 2*n)=  1d0
              case ('dxy', 'Dxy', 'DXY')
                 mirror_z(2*n-1, 2*n-1)= 1d0
                 mirror_z(2*n, 2*n)= -1d0
              case ('dyz', 'Dyz', 'DYZ')
                 mirror_z(2*n-1, 2*n-1)=-1d0
                 mirror_z(2*n, 2*n)=  1d0
              case ('dxz', 'Dxz', 'DXZ', 'dzx', 'Dzx', 'DZX')
                 mirror_z(2*n-1, 2*n-1)=-1d0
                 mirror_z(2*n, 2*n)=  1d0
              case ('dx2-', 'dx2-y2', 'Dx2-y2', 'DX2-Y2', 'dx2', 'DX2')
                 mirror_z(2*n-1, 2*n-1)= 1d0
                 mirror_z(2*n, 2*n)= -1d0
              case ('dz2', 'Dz2', 'DZ2')
                 mirror_z(2*n-1, 2*n-1)= 1d0
                 mirror_z(2*n, 2*n)= -1d0
              case default
                 write(*, *) "ERROR: only support s px py pz dxy dyz dxz dx2-y2 dz2 orbitals"
                 stop
              end select
           enddo ! i
        enddo ! ia
 
     endif


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
