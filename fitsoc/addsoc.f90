
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

     
      Num_wann_new= 2* Num_wann_nsoc

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
               H_soc(i1, j1+ Num_wann_nsoc)= soc_spd(i2, j2+ LMMAX)* lambda0
               H_soc(i1+ Num_wann_nsoc, j1)= soc_spd(i2+ LMMAX, j2)* lambda0
               H_soc(i1+ Num_wann_nsoc, j1+ Num_wann_nsoc)= &
                  soc_spd(i2+ LMMAX, j2+ LMMAX)* lambda0
            enddo ! j
         enddo ! i
      enddo ! ia

      do ir=1, nrpts_nsoc
         if (irvec_nsoc(1, ir).ne.0 .or. irvec_nsoc(2, ir).ne.0 .or. irvec_nsoc(3, ir).ne.0) cycle
         HmnR_nsoc(:, :, ir)= HmnR_nsoc_origin(:, :, ir)+ H_soc* ndegen_nsoc(ir)
      enddo

      return
   end subroutine  addsoc_all

   subroutine  addsoc_pd_zjw

      use para
      implicit none


      integer :: i, j, ir, ia
      integer :: Num_wann_new
      real(dp), parameter :: sqrt3=1.732051d0

      real(dp) :: lp
      real(dp) :: ld

      complex(dp), allocatable :: h_soc(:, :)

      complex(dp), allocatable :: soc_p(:, :)
      complex(dp), allocatable :: soc_d(:, :)

     
      Num_wann_new= 2* Num_wann_nsoc

      allocate(soc_p(6, 6))
      allocate(soc_d(10, 10))
      allocate (H_soc(Num_wann_new, Num_wann_new))

      soc_p= 0d0
      soc_d= 0d0
      H_soc= 0d0 

      !> soc term for p orbitals
      !> 1:3 up  4:6 dn
      !> px py pz
      soc_p(1,2)=-zi;    soc_p(1,6)=1.0d0
      soc_p(2,1)= zi;    soc_p(2,6)=-zi
      soc_p(3,4)=-1.0d0; soc_p(3,5)= zi
      soc_p(4,3)=-1.0d0; soc_p(4,5)= zi
      soc_p(5,3)=-zi;    soc_p(5,4)=-zi
      soc_p(6,1)= 1.0d0; soc_p(6,2)= zi

      !> 1:3 up  4:6 dn
      !> pz px py
      soc_p=0d0
      soc_p(2,3)=-zi;    soc_p(2,4)=1.0d0
      soc_p(3,2)= zi;    soc_p(3,4)=-zi
      soc_p(1,5)=-1.0d0; soc_p(1,6)= zi
      soc_p(5,1)=-1.0d0; soc_p(5,6)= zi
      soc_p(6,1)=-zi;    soc_p(6,5)=-zi
      soc_p(4,2)= 1.0d0; soc_p(4,3)= zi



      !> soc term for d orbitals
      !> 1:5 up 6:10 dn
      !> d orbitals should range like this
      !> dxy, dyz, dzx, dx2-y2 d3z2-r2
      soc_d(1, 4)= 2d0*zi; soc_d(1, 7)= 1d0; soc_d(1, 8)= -zi
      soc_d(2, 3)= zi; soc_d(2, 6)=-1d0; soc_d(2, 9)=-zi; soc_d(2,10)= -sqrt3*zi
      soc_d(3, 6)= zi; soc_d(3, 9)= -1; soc_d(3, 10)= sqrt3
      soc_d(4, 7)= zi; soc_d(4, 8)=  1
      soc_d(5, 7)= sqrt3*zi; soc_d(5, 8)= -sqrt3
      soc_d(6, 9)= -2d0*zi
      soc_d(7, 8)= -zi
      do i=1, 10
         do j=i+1, 10
            soc_d(j, i)= conjg(soc_d(i, j))
         enddo
      enddo

      !> d orbitals should range like this
      !> k3z2-r2, dxz, dyz, dx2-y2 dxy
      soc_d= 0d0
      soc_d(5, 4)= 2d0*zi; soc_d(5, 8)= 1d0; soc_d(5, 7)= -zi
      soc_d(3, 2)= zi; soc_d(3, 10)=-1d0; soc_d(3, 9)=-zi; soc_d(3,6)= -sqrt3*zi
      soc_d(2, 10)= zi; soc_d(2, 9)= -1; soc_d(2, 6)= sqrt3
      soc_d(4, 8)= zi; soc_d(4, 7)=  1
      soc_d(1, 8)= sqrt3*zi; soc_d(1, 7)= -sqrt3
      soc_d(10, 9)= -2d0*zi
      soc_d(8, 7)= -zi
 
      do i=1, 10
         do j=i+1, 10
            soc_d(j, i)= conjg(soc_d(i, j))
         enddo
      enddo

      ld= lambda_d(1)
      lp= lambda_p(5)

      !> 1:20 45:64 is the part of W orbitals
      do ia=1, 4
         H_soc((ia-1)*5+1:(ia-1)*5+5, (ia-1)*5+1:(ia-1)*5+5)= ld*soc_d(1:5, 1:5)/2d0
         H_soc((ia-1)*5+45:(ia-1)*5+49, (ia-1)*5+1:(ia-1)*5+5)= ld*soc_d(6:10, 1:5)/2d0
         H_soc((ia-1)*5+1:(ia-1)*5+5, (ia-1)*5+45:(ia-1)*5+49)= ld*soc_d(1:5, 6:10)/2d0
         H_soc((ia-1)*5+45:(ia-1)*5+49, (ia-1)*5+45:(ia-1)*5+49)= ld*soc_d(6:10, 6:10)/2d0
      enddo

      !> 21:44 65:88 is the part of Te orbitals
      do ia=1, 8
         H_soc((ia-1)*3+21:(ia-1)*3+23, (ia-1)*3+21:(ia-1)*3+23)= lp*soc_p(1:3, 1:3)/2d0
         H_soc((ia-1)*3+65:(ia-1)*3+67, (ia-1)*3+21:(ia-1)*3+23)= lp*soc_p(4:6, 1:3)/2d0
         H_soc((ia-1)*3+21:(ia-1)*3+23, (ia-1)*3+65:(ia-1)*3+67)= lp*soc_p(1:3, 4:6)/2d0
         H_soc((ia-1)*3+65:(ia-1)*3+67, (ia-1)*3+65:(ia-1)*3+67)= lp*soc_p(4:6, 4:6)/2d0
      enddo

      if (Num_wann_new.ne.88) stop' Num_wann should equal 44'
      do ir=1, nrpts_nsoc
         if (irvec_nsoc(1, ir).ne.0 .or. irvec_nsoc(2, ir).ne.0 .or. irvec_nsoc(3, ir).ne.0) cycle
         HmnR_nsoc(:, :, ir)= HmnR_nsoc_origin(:, :, ir)+ H_soc* ndegen_nsoc(ir)
      enddo

      return
   end subroutine  addsoc_pd_zjw




   subroutine  addsoc_pd

      use para
      implicit none


      integer :: i, j, ir, ia
      integer :: Num_wann_new
      real(dp), parameter :: sqrt3=1.732051d0

      real(dp) :: lp
      real(dp) :: ld

      complex(dp), allocatable :: h_soc(:, :)

      complex(dp), allocatable :: soc_p(:, :)
      complex(dp), allocatable :: soc_d(:, :)

     
      Num_wann_new= 2* Num_wann_nsoc

      allocate(soc_p(6, 6))
      allocate(soc_d(10, 10))
      allocate (H_soc(Num_wann_new, Num_wann_new))

      soc_p= 0d0
      soc_d= 0d0
      H_soc= 0d0 

      !> soc term for p orbitals
      !> 1:3 up  4:6 dn
      soc_p(1,2)=-zi;    soc_p(1,6)=1.0d0
      soc_p(2,1)= zi;    soc_p(2,6)=-zi
      soc_p(3,4)=-1.0d0; soc_p(3,5)= zi
      soc_p(4,3)=-1.0d0; soc_p(4,5)= zi
      soc_p(5,3)=-zi;    soc_p(5,4)=-zi
      soc_p(6,1)= 1.0d0; soc_p(6,2)= zi


      !> soc term for d orbitals
      !> 1:5 up 6:10 dn
      !> d orbitals should range like this
      !> dxy, dyz, dzx, dx2-y2 d3z2-r2
      soc_d(1, 4)= 2d0*zi; soc_d(1, 7)= 1d0; soc_d(1, 8)= -zi
      soc_d(2, 3)= zi; soc_d(2, 6)=-1d0; soc_d(2, 9)=-zi; soc_d(2,10)= -sqrt3*zi
      soc_d(3, 6)= zi; soc_d(3, 9)= -1; soc_d(3, 10)= sqrt3
      soc_d(4, 7)= zi; soc_d(4, 8)=  1
      soc_d(5, 7)= sqrt3*zi; soc_d(5, 8)= -sqrt3
      soc_d(6, 9)= -2d0*zi
      soc_d(7, 8)= -zi
      do i=1, 10
         do j=i+1, 10
            soc_d(j, i)= conjg(soc_d(i, j))
         enddo
      enddo

      !> d orbitals should range like this
      !>d3z2-r2, dx2-y2, dyz, dxy, dzx,  
      soc_d= 0d0
      soc_d(4, 2)= 2d0*zi; soc_d(4, 8)= 1d0; soc_d(4, 10)= -zi
      soc_d(3, 5)= zi; soc_d(3, 9)=-1d0; soc_d(3, 7)=-zi; soc_d(3,6)= -sqrt3*zi
      soc_d(5, 9)= zi; soc_d(5, 7)= -1; soc_d(5, 6)= sqrt3
      soc_d(2, 8)= zi; soc_d(2, 10)=  1
      soc_d(1, 8)= sqrt3*zi; soc_d(1, 10)= -sqrt3
      soc_d(9, 7)= -2d0*zi
      soc_d(8, 10)= -zi
      do i=1, 10
         do j=i+1, 10
            soc_d(j, i)= conjg(soc_d(i, j))
         enddo
      enddo

      ld= lambda_d(1)
      lp= lambda_p(5)

      !> 1:20 45:64 is the part of W orbitals
      do ia=1, 4
         H_soc((ia-1)*5+1:(ia-1)*5+5, (ia-1)*5+1:(ia-1)*5+5)= ld*soc_d(1:5, 1:5)/2d0
         H_soc((ia-1)*5+45:(ia-1)*5+49, (ia-1)*5+1:(ia-1)*5+5)= ld*soc_d(6:10, 1:5)/2d0
         H_soc((ia-1)*5+1:(ia-1)*5+5, (ia-1)*5+45:(ia-1)*5+49)= ld*soc_d(1:5, 6:10)/2d0
         H_soc((ia-1)*5+45:(ia-1)*5+49, (ia-1)*5+45:(ia-1)*5+49)= ld*soc_d(6:10, 6:10)/2d0
      enddo

      !> 21:44 65:88 is the part of Te orbitals
      do ia=1, 8
         H_soc((ia-1)*3+21:(ia-1)*3+23, (ia-1)*3+21:(ia-1)*3+23)= lp*soc_p(1:3, 1:3)/2d0
         H_soc((ia-1)*3+65:(ia-1)*3+67, (ia-1)*3+21:(ia-1)*3+23)= lp*soc_p(4:6, 1:3)/2d0
         H_soc((ia-1)*3+21:(ia-1)*3+23, (ia-1)*3+65:(ia-1)*3+67)= lp*soc_p(1:3, 4:6)/2d0
         H_soc((ia-1)*3+65:(ia-1)*3+67, (ia-1)*3+65:(ia-1)*3+67)= lp*soc_p(4:6, 4:6)/2d0
      enddo

      if (Num_wann_new.ne.88) stop' Num_wann should equal 44'
      do ir=1, nrpts_nsoc
         if (irvec_nsoc(1, ir).ne.0 .or. irvec_nsoc(2, ir).ne.0 .or. irvec_nsoc(3, ir).ne.0) cycle
         HmnR_nsoc(:, :, ir)= HmnR_nsoc_origin(:, :, ir)+ H_soc* ndegen_nsoc(ir)
      enddo

      return
   end subroutine  addsoc_pd


   subroutine  addsoc_p

      use para
      implicit none


      integer :: i, j, ir
      integer :: Num_wann_new

      complex(dp), allocatable :: h_soc(:, :)

      complex(dp), allocatable :: soc_p(:, :)

     
      Num_wann_new= 2* Num_wann_nsoc

      allocate(soc_p(6, 6))
      allocate (H_soc(Num_wann_new, Num_wann_new))

      soc_p= 0d0
      H_soc= 0d0 

      !> soc term for p orbitals
      !> 1:3 up  4:6 dn
      soc_p(1,2)=-zi;    soc_p(1,6)=1.0d0
      soc_p(2,1)= zi;    soc_p(2,6)=-zi
      soc_p(3,4)=-1.0d0; soc_p(3,5)= zi
      soc_p(4,3)=-1.0d0; soc_p(4,5)= zi
      soc_p(5,3)=-zi;    soc_p(5,4)=-zi
      soc_p(6,1)= 1.0d0; soc_p(6,2)= zi

      !> the orbital order like this
      !> s px py pz px py pz
      H_soc(2:4, 2:4)  = lambda_p(1)*soc_p(1:3, 1:3)/2d0
      H_soc(9:11, 9:11)= lambda_p(1)*soc_p(4:6, 4:6)/2d0
      H_soc(2:4, 9:11) = lambda_p(1)*soc_p(1:3, 4:6)/2d0
      H_soc(9:11, 2:4) = lambda_p(1)*soc_p(4:6, 1:3)/2d0

      H_soc(5:7, 5:7)    = lambda_p(2)*soc_p(1:3, 1:3)/2d0
      H_soc(12:14, 12:14)= lambda_p(2)*soc_p(4:6, 4:6)/2d0
      H_soc(5:7, 12:14)  = lambda_p(2)*soc_p(1:3, 4:6)/2d0
      H_soc(12:14, 5:7)  = lambda_p(2)*soc_p(4:6, 1:3)/2d0

      if (Num_wann_new.ne.14) stop' Num_wann should equal 7'
      do ir=1, nrpts_nsoc
         if (irvec_nsoc(1, ir).ne.0 .or. irvec_nsoc(2, ir).ne.0 .or. irvec_nsoc(3, ir).ne.0) cycle
         HmnR_nsoc(:, :, ir)= HmnR_nsoc_origin(:, :, ir)+ H_soc* ndegen_nsoc(ir)
      enddo

      do i=1, num_wann_soc
         write(10102, '(2000f6.2)')real(H_soc(:, i))
      enddo
      return
   end subroutine  addsoc_p

   subroutine  addsoc

      use para
      implicit none


      integer :: i, j, ir
      integer :: Num_wann_new

      complex(dp), allocatable :: h_soc(:, :)

      complex(dp), allocatable :: soc_p(:, :)

     
      Num_wann_new= 2* Num_wann_nsoc

      allocate(soc_p(6, 6))
      allocate (H_soc(Num_wann_new, Num_wann_new))

      soc_p= 0d0
      H_soc= 0d0 

      !> soc term for p orbitals
      !> 1:3 up  4:6 dn
      soc_p(1,2)=-zi;    soc_p(1,6)=1.0d0
      soc_p(2,1)= zi;    soc_p(2,6)=-zi
      soc_p(3,4)=-1.0d0; soc_p(3,5)= zi
      soc_p(4,3)=-1.0d0; soc_p(4,5)= zi
      soc_p(5,3)=-zi;    soc_p(5,4)=-zi
      soc_p(6,1)= 1.0d0; soc_p(6,2)= zi

      !> the orbital order like this
      !> s px py pz px py pz
      H_soc(2:4, 2:4)  = lambda_p(1)*soc_p(1:3, 1:3)/2d0
      H_soc(9:11, 9:11)= lambda_p(1)*soc_p(4:6, 4:6)/2d0
      H_soc(2:4, 9:11) = lambda_p(1)*soc_p(1:3, 4:6)/2d0
      H_soc(9:11, 2:4) = lambda_p(1)*soc_p(4:6, 1:3)/2d0

      H_soc(5:7, 5:7)    = lambda_p(2)*soc_p(1:3, 1:3)/2d0
      H_soc(12:14, 12:14)= lambda_p(2)*soc_p(4:6, 4:6)/2d0
      H_soc(5:7, 12:14)  = lambda_p(2)*soc_p(1:3, 4:6)/2d0
      H_soc(12:14, 5:7)  = lambda_p(2)*soc_p(4:6, 1:3)/2d0

      if (Num_wann_new.ne.14) stop' Num_wann should equal 7'
      do ir=1, nrpts_nsoc
         if (irvec_nsoc(1, ir).ne.0 .or. &
            irvec_nsoc(2, ir).ne.0 .or. &
            irvec_nsoc(3, ir).ne.0) cycle
         HmnR_nsoc(:, :, ir)= HmnR_nsoc_origin(:, :, ir)+ H_soc* ndegen_nsoc(ir)
      enddo

      deallocate(soc_p, H_soc)
      return
   end subroutine  addsoc
