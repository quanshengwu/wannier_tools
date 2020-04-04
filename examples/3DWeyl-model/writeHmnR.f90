! 2-band 3D WSM model
! usage:
! compile and run
! gfortran writeHmnR.f90 -o writehmnr
! ./writehmnr
! > H=A(kx*s_x+ky*s_y)+(M0-M1(kx*kx+ky*ky+kz*kz))*s_z
   program writeHmnR

      implicit none

      integer, parameter :: dp=kind(1d0)
      complex(dp), parameter :: zi= (0d0, 1d0)
      complex(dp), parameter :: zzero= (0d0, 0d0)

      integer :: i, j
      integer :: ir
      integer :: nwann

      !> arrays for hamiltonian storage
      integer :: nrpts
      integer, allocatable :: ndegen(:)
      integer, allocatable :: irvec(:, :)
      complex(dp), allocatable :: hmnr(:, :, :)

      !> three lattice constants
      real(dp) :: A, M0, M1


      A=1d0
      M0=1d0
      M1=1d0

      nwann= 1
      nrpts=17
      allocate(irvec(3, nrpts))
      allocate(ndegen(nrpts))
      allocate(hmnr(nwann*2, nwann*2, nrpts))
      irvec=0
      ndegen=1
      hmnr= zzero


      ! 0 0 0
      ir= 1
      irvec(1, ir)= 0
      irvec(2, ir)= 0
      irvec(3, ir)= 0
      hmnr(1, 1, ir)= M0- 6d0*M1
      hmnr(2, 2, ir)=-M0+ 6d0*M1

      !1 0 0
      ir= ir+ 1
      irvec(1, ir)= 1
      irvec(2, ir)= 0
      irvec(3, ir)= 0
      hmnr(1, 1, ir)=M1
      hmnr(1, 2, ir)=zi*A/2d0
      hmnr(2, 1, ir)=zi*A/2d0
      hmnr(2, 2, ir)=-M1

      !-1 0 0
      ir= ir+ 1
      irvec(1, ir)=-1
      irvec(2, ir)= 0
      irvec(3, ir)= 0
      hmnr(1, 1, ir)=M1
      hmnr(1, 2, ir)=-zi*A/2d0
      hmnr(2, 1, ir)=-zi*A/2d0
      hmnr(2, 2, ir)=-M1

      ! 0 1 0
      ir= ir+ 1
      irvec(1, ir)= 0 
      irvec(2, ir)= 1
      irvec(3, ir)= 0
      hmnr(1, 1, ir)=M1
      hmnr(1, 2, ir)=A/2d0
      hmnr(2, 1, ir)=-A/2d0
      hmnr(2, 2, ir)=-M1

      !0 -1  0
      ir= ir+ 1
      irvec(1, ir)= 0 
      irvec(2, ir)=-1
      irvec(3, ir)= 0
      hmnr(1, 1, ir)=M1
      hmnr(1, 2, ir)=-A/2d0
      hmnr(2, 1, ir)=A/2d0
      hmnr(2, 2, ir)=-M1
 
      ! 0  0  1
      ir= ir+ 1
      irvec(1, ir)= 0 
      irvec(2, ir)= 0
      irvec(3, ir)= 1
      hmnr(1, 1, ir)= M1
      hmnr(2, 2, ir)=-M1

      ! 0  0 -1
      ir= ir+ 1
      irvec(1, ir)= 0 
      irvec(2, ir)= 0
      irvec(3, ir)=-1
      hmnr(1, 1, ir)= M1
      hmnr(2, 2, ir)=-M1

      nrpts = ir

      !> write to new_hr.dat
      open(unit=105, file='Weyl3D_hr.dat')
      write(105, *)'2-band 3D WSM toy model'
      write(105, *)nwann*2
      write(105, *)nrpts
      write(105, '(15I5)')(ndegen(i), i=1, nrpts)
      do ir=1, nrpts
         do i=1, nwann*2
            do j=1, nwann*2
               write(105, '(5I5, 2f16.8)')irvec(:, ir), i, j, HmnR(i, j, ir)
            enddo
         enddo
      enddo
      close(105)

   end ! end of program 
