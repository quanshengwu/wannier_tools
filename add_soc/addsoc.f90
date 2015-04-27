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

     
      Num_wann_new= 2* Num_wann
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

     
      Num_wann_new= 2* Num_wann

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
   end subroutine  addsoc_pd


   subroutine  addsoc_p

      use para
      implicit none


      integer :: i, j, ir
      integer :: Num_wann_new

      complex(dp), allocatable :: h_soc(:, :)

      complex(dp), allocatable :: soc_p(:, :)

     
      Num_wann_new= 2* Num_wann

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
               write(105, '(5I5, 2f16.8)')irvec(:, ir), HmnR(i, j, ir)
            enddo
         enddo
      enddo
      close(105)


      return
   end subroutine  addsoc_p

   subroutine  addsoc

      use para
      implicit none


      integer :: i, j, ir
      integer :: Num_wann_new

      complex(dp), allocatable :: h_soc(:, :)

      complex(dp), allocatable :: soc_p(:, :)

     
      Num_wann_new= 2* Num_wann

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
      do ir=1, nrpts
         if (irvec(1, ir).ne.0 .or. &
            irvec(2, ir).ne.0 .or. &
            irvec(3, ir).ne.0) cycle
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
               write(105, '(5I5, 2f16.8)')irvec(:, ir), HmnR(i, j, ir)
            enddo
         enddo
      enddo
      close(105)


      deallocate(soc_p, H_soc)
      return
   end subroutine  addsoc
