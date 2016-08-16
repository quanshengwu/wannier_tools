
!> add spin-orbital term for all kinds of situation
!> we can get the orbital order from the input file fit_soc.in
!> for p orbitals, the name of the orbital should be px, py, pz
!> for d orbitals, the name of the orbital should be dxy, dyz, dxz, dx2-y2, dz2
!> the order of these orbitals can be specified by users
   subroutine  addsoc_all

      use para
      implicit none


      integer :: i, j, ir, ia
      integer :: i1, i2, j1, j2, it
      integer :: Num_wann_new

      !> includind s  px py pz dxy  dyz dxz dx2-y2 dz2
      integer :: LMMAX= 9

      real(dp), parameter :: sqrt3=1.732051d0

      real(dp) :: lambda0

      complex(dp), allocatable :: h_soc(:, :)

      complex(dp), allocatable :: soc0_p(:, :)
      complex(dp), allocatable :: soc0_d(:, :)
      complex(dp), allocatable :: soc_spd(:, :)

      integer, allocatable :: start(:)
      integer, allocatable :: orbital_index(:, :)
      real(dp), allocatable :: lambda(:, :)

     
      Num_wann_new= 2* Num_wann

      allocate(soc0_p(6, 6))
      allocate(soc0_d(10, 10))
      allocate(soc_spd(LMMAX*2, LMMAX*2))
      allocate(H_soc(Num_wann_new, Num_wann_new))
      allocate(start(Num_atoms))
      allocate(lambda(max_projs, Num_atoms))
      allocate(orbital_index(max_projs, Num_atoms))

      orbital_index= 0
      soc0_p= 0d0
      soc0_d= 0d0
      soc_spd= 0d0
      H_soc= 0d0 

      !> soc term for p orbitals
      !> 1:3 up  4:6 dn
      !> px py pz
      soc0_p(1,2)=-zi;    soc0_p(1,6)=1.0d0
      soc0_p(2,1)= zi;    soc0_p(2,6)=-zi
      soc0_p(3,4)=-1.0d0; soc0_p(3,5)= zi
      soc0_p(4,3)=-1.0d0; soc0_p(4,5)= zi
      soc0_p(5,3)=-zi;    soc0_p(5,4)=-zi
      soc0_p(6,1)= 1.0d0; soc0_p(6,2)= zi

      !> soc term for d orbitals
      !> 1:5 up 6:10 dn
      !> d orbitals should range like this
      !> dxy, dyz, dxz, dx2-y2 dz2
      soc0_d(1, 4)= 2d0*zi; soc0_d(1, 7)= 1d0; soc0_d(1, 8)= -zi
      soc0_d(2, 3)= zi; soc0_d(2, 6)=-1d0; soc0_d(2, 9)=-zi; soc0_d(2,10)= -sqrt3*zi
      soc0_d(3, 6)= zi; soc0_d(3, 9)= -1; soc0_d(3, 10)= sqrt3
      soc0_d(4, 7)= zi; soc0_d(4, 8)=  1
      soc0_d(5, 7)= sqrt3*zi; soc0_d(5, 8)= -sqrt3
      soc0_d(6, 9)= -2d0*zi
      soc0_d(7, 8)= -zi
      do i=1, 10
         do j=i+1, 10
            soc0_d(j, i)= conjg(soc0_d(i, j))
         enddo
      enddo

      !> get soc 
      soc_spd(2:4, 2:4)= soc0_p(1:3, 1:3)
      soc_spd(2:4, 11:13)= soc0_p(1:3, 4:6)
      soc_spd(11:13, 2:4)= soc0_p(4:6, 1:3)
      soc_spd(11:13, 11:13)= soc0_p(4:6, 4:6)

      !> get soc_d 
      soc_spd(5:9, 5:9)= soc0_d(1:5, 1:5)
      soc_spd(5:9, 14:18)= soc0_d(1:5, 6:10)
      soc_spd(14:18, 5:9)= soc0_d(6:10, 1:5)
      soc_spd(14:18, 14:18)= soc0_d(6:10, 6:10)

      lambda = 0d0
      orbital_index= 1
      !> the orbital order for default is
      !> s  px py pz dxy  dyz dxz dx2-y2 dz2
      !> 1  2  3  4  5    6   7   8      9

      !> get the orbital order from user's input 
      do ia=1, Num_atoms
         it= atom_type(ia)
         do i=1, nprojs(ia)
            select case (proj_name(i, ia))
            case ('s', 'S')
               orbital_index(i, ia)= 1
               lambda(i, ia)= 0d0
            case ('px', 'Px', 'PX')
               orbital_index(i, ia)= 2
               lambda(i, ia)= lambda_p(it)/2d0
            case ('py', 'Py', 'PY')
               orbital_index(i, ia)= 3
               lambda(i, ia)= lambda_p(it)/2d0
            case ('pz', 'Pz', 'PZ')
               orbital_index(i, ia)= 4
               lambda(i, ia)= lambda_p(it)/2d0
            case ('dxy', 'Dxy', 'DXY')
               orbital_index(i, ia)= 5
               lambda(i, ia)= lambda_d(it)/2d0
            case ('dyz', 'Dyz', 'DYZ')
               orbital_index(i, ia)= 6
               lambda(i, ia)= lambda_d(it)/2d0
            case ('dxz', 'Dxz', 'DXZ')
               orbital_index(i, ia)= 7
               lambda(i, ia)= lambda_d(it)/2d0
            case ('dx2-y2', 'Dx2-y2', 'DX2-Y2', 'dx2', 'DX2')
               orbital_index(i, ia)= 8
               lambda(i, ia)= lambda_d(it)/2d0
            case ('dz2', 'Dz2', 'DZ2')
               orbital_index(i, ia)= 9
               lambda(i, ia)= lambda_d(it)/2d0
            case default
               write(*, *) "ERROR: only support s px py pz dxy dyz dxz dx2-y2 dz2 &
                  orbitals"
               write(*, *)proj_name(i, ia)
               stop
            end select
         enddo
      enddo

      start(1)= 1
      do ia=2, Num_atoms
         start(ia)= start(ia-1)+ nprojs(ia-1)
      enddo

      do ia=1, Num_atoms
         do i=1, nprojs(ia)
            i1= start(ia)+ i- 1
            i2= orbital_index(i, ia) 
            lambda0= lambda(i, ia) ! soc strength for each orbital of atoms
            do j=1, nprojs(ia)
               j1= start(ia)+ j- 1
               j2= orbital_index(j, ia) 
               H_soc(i1, j1)= soc_spd(i2, j2)* lambda0
               H_soc(i1, j1+ Num_wann)= soc_spd(i2, j2+ LMMAX)* lambda0
               H_soc(i1+ Num_wann, j1)= soc_spd(i2+ LMMAX, j2)* lambda0
               H_soc(i1+ Num_wann, j1+ Num_wann)= &
                  soc_spd(i2+ LMMAX, j2+ LMMAX)* lambda0
            enddo ! j
         enddo ! i
      enddo ! ia

      do ir=1, nrpts
         if (irvec(1, ir).ne.0 .or. irvec(2, ir).ne.0 .or. irvec(3, ir).ne.0) cycle
         HmnR(:, :, ir)= HmnR(:, :, ir)+ H_soc* ndegen(ir)
      enddo

      !> write to new_hr.dat
      open(unit=105, file='new_hr.dat')
      write(105, *)'new HmnR from hand made soc'
      write(105, *)Num_wann_new
      write(105, *)nrpts
      write(105, '(15I5)')(ndegen(i), i=1, nrpts)
      do ir=1, nrpts
         do i=1, Num_wann_new
            do j=1, Num_wann_new
               write(105, '(5I5, 2f16.8)')irvec(:, ir), i, j, HmnR(i, j, ir)
            enddo
         enddo
      enddo
      close(105)


      return
   end subroutine  addsoc_all


!> add spin-orbital term for all kinds of situation
!> we can get the orbital order from the input file fit_soc.in
!> for p orbitals, the name of the orbital should be px, py, pz
!> for d orbitals, the name of the orbital should be dxy, dyz, dxz, dx2-y2, dz2
!> the order of these orbitals can be specified by users
   subroutine  addsoc_all_2D

      use para
      implicit none


      integer :: i, j, ir, ia
      integer :: i1, i2, j1, j2, it
      integer :: Num_wann_new

      !> includind s  px py pz dxy  dyz dxz dx2-y2 dz2
      integer :: LMMAX= 9

      real(dp), parameter :: sqrt3=1.732051d0

      real(dp) :: lambda0

      complex(dp), allocatable :: h_soc(:, :)

      complex(dp), allocatable :: soc0_p(:, :)
      complex(dp), allocatable :: soc0_d(:, :)
      complex(dp), allocatable :: soc_spd(:, :)

      integer, allocatable :: start(:)
      integer, allocatable :: orbital_index(:, :)
      real(dp), allocatable :: lambda(:, :)

     
      Num_wann_new= 2* Num_wann

      allocate(soc0_p(6, 6))
      allocate(soc0_d(10, 10))
      allocate(soc_spd(LMMAX*2, LMMAX*2))
      allocate(H_soc(Num_wann_new, Num_wann_new))
      allocate(start(Num_atoms))
      allocate(lambda(max_projs, Num_atoms))
      allocate(orbital_index(max_projs, Num_atoms))

      orbital_index= 0
      soc0_p= 0d0
      soc0_d= 0d0
      soc_spd= 0d0
      H_soc= 0d0 

      !> soc term for p orbitals
      !> 1:3 up  4:6 dn
      !> px py pz
      soc0_p(1,2)=-zi;    soc0_p(1,6)=1.0d0
      soc0_p(2,1)= zi;    soc0_p(2,6)=-zi
      soc0_p(3,4)=-1.0d0; soc0_p(3,5)= zi
      soc0_p(4,3)=-1.0d0; soc0_p(4,5)= zi
      soc0_p(5,3)=-zi;    soc0_p(5,4)=-zi
      soc0_p(6,1)= 1.0d0; soc0_p(6,2)= zi

      !> soc term for d orbitals
      !> 1:5 up 6:10 dn
      !> d orbitals should range like this
      !> dxy, dyz, dxz, dx2-y2 dz2
      soc0_d(1, 4)= 2d0*zi; soc0_d(1, 7)= 1d0; soc0_d(1, 8)= -zi
      soc0_d(2, 3)= zi; soc0_d(2, 6)=-1d0; soc0_d(2, 9)=-zi; soc0_d(2,10)= -sqrt3*zi
      soc0_d(3, 6)= zi; soc0_d(3, 9)= -1; soc0_d(3, 10)= sqrt3
      soc0_d(4, 7)= zi; soc0_d(4, 8)=  1
      soc0_d(5, 7)= sqrt3*zi; soc0_d(5, 8)= -sqrt3
      soc0_d(6, 9)= -2d0*zi
      soc0_d(7, 8)= -zi
      do i=1, 10
         do j=i+1, 10
            soc0_d(j, i)= conjg(soc0_d(i, j))
         enddo
      enddo

      !> get soc 
      soc_spd(2:4, 2:4)= soc0_p(1:3, 1:3)
      soc_spd(2:4, 11:13)= soc0_p(1:3, 4:6)
      soc_spd(11:13, 2:4)= soc0_p(4:6, 1:3)
      soc_spd(11:13, 11:13)= soc0_p(4:6, 4:6)

      !> get soc_d 
      soc_spd(5:9, 5:9)= soc0_d(1:5, 1:5)
      soc_spd(5:9, 14:18)= soc0_d(1:5, 6:10)
      soc_spd(14:18, 5:9)= soc0_d(6:10, 1:5)
      soc_spd(14:18, 14:18)= soc0_d(6:10, 6:10)

      lambda = 0d0
      orbital_index= 1
      !> the orbital order for default is
      !> s  px py pz dxy  dyz dxz dx2-y2 dz2
      !> 1  2  3  4  5    6   7   8      9

      !> get the orbital order from user's input 
      do ia=1, Num_atoms
         it= atom_type(ia)
         do i=1, nprojs(ia)
            select case (proj_name(i, ia))
            case ('s', 'S')
               orbital_index(i, ia)= 1
               lambda(i, ia)= 0d0
            case ('px', 'Px', 'PX')
               orbital_index(i, ia)= 2
               lambda(i, ia)= lambda_p(it)/2d0
            case ('py', 'Py', 'PY')
               orbital_index(i, ia)= 3
               lambda(i, ia)= lambda_p(it)/2d0
            case ('pz', 'Pz', 'PZ')
               orbital_index(i, ia)= 4
               lambda(i, ia)= lambda_p(it)/2d0
            case ('dxy', 'Dxy', 'DXY')
               orbital_index(i, ia)= 5
               lambda(i, ia)= lambda_d(it)/2d0
            case ('dyz', 'Dyz', 'DYZ')
               orbital_index(i, ia)= 6
               lambda(i, ia)= lambda_d(it)/2d0
            case ('dxz', 'Dxz', 'DXZ')
               orbital_index(i, ia)= 7
               lambda(i, ia)= lambda_d(it)/2d0
            case ('dx2-y2', 'Dx2-y2', 'DX2-Y2', 'dx2', 'DX2')
               orbital_index(i, ia)= 8
               lambda(i, ia)= lambda_d(it)/2d0
            case ('dz2', 'Dz2', 'DZ2')
               orbital_index(i, ia)= 9
               lambda(i, ia)= lambda_d(it)/2d0
            case default
               write(*, *) "ERROR: only support s px py pz dxy dyz dxz dx2-y2 dz2 &
                  orbitals"
               write(*, *)proj_name(i, ia)
               stop
            end select
         enddo
      enddo

      start(1)= 1
      do ia=2, Num_atoms
         start(ia)= start(ia-1)+ nprojs(ia-1)
      enddo

      do ia=1, Num_atoms
         do i=1, nprojs(ia)
            i1= start(ia)+ i- 1
            i2= orbital_index(i, ia) 
            lambda0= lambda(i, ia) ! soc strength for each orbital of atoms
            do j=1, nprojs(ia)
               j1= start(ia)+ j- 1
               j2= orbital_index(j, ia) 
               H_soc(i1, j1)= soc_spd(i2, j2)* lambda0
               H_soc(i1, j1+ Num_wann)= soc_spd(i2, j2+ LMMAX)* lambda0
               H_soc(i1+ Num_wann, j1)= soc_spd(i2+ LMMAX, j2)* lambda0
               H_soc(i1+ Num_wann, j1+ Num_wann)= &
                  soc_spd(i2+ LMMAX, j2+ LMMAX)* lambda0
            enddo ! j
         enddo ! i
      enddo ! ia

      do ir=1, nrpts
         if (irvec_2d(1, ir).ne.0 .or. irvec_2d(2, ir).ne.0) cycle
         HmnR(:, :, ir)= HmnR(:, :, ir)+ H_soc* ndegen(ir)
      enddo

      !> write to new_hr.dat
      open(unit=105, file='new_hr.dat')
      write(105, *)'new HmnR from hand made soc'
      write(105, *)Num_wann_new
      write(105, *)nrpts
      write(105, '(15I5)')(ndegen(i), i=1, nrpts)
      do ir=1, nrpts
         do i=1, Num_wann_new
            do j=1, Num_wann_new
               write(105, '(4I5, 2f16.8)')irvec_2d(:, ir), i, j, HmnR(i, j, ir)
            enddo
         enddo
      enddo
      close(105)


      return
   end subroutine  addsoc_all_2D

