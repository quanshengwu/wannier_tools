   program writeHmnR

      !> two band model
      !> H11(k)= -1 +2(mx+my+mz)- 2mx*cos(kx)- 2*my*cos(ky)- 2*mz*cos(kz)
      !> H22(k)=  1 -2(mx+my+mz)- 2mx*cos(kx)+ 2*my*cos(ky)+ 2*mz*cos(kz)
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
      real(dp) :: mx, my, mz


      mx=1d0
      my=1d0
      mz=1d0

      nwann= 2
      nrpts=17
      allocate(irvec(3, nrpts))
      allocate(ndegen(nrpts))
      allocate(hmnr(nwann, nwann, nrpts))
      irvec=0
      ndegen=1
      hmnr= zzero


      ! 0 0 0
      ir= 1
      irvec(1, ir)= 0
      irvec(2, ir)= 0
      irvec(3, ir)= 0
      hmnr(1, 1, ir)= -1.d0+ 2d0*(mx+my+mz)
      hmnr(2, 2, ir)=  1.d0- 2d0*(mx+my+mz)

      !1 0 0
      ir= ir+ 1
      irvec(1, ir)= 1
      irvec(2, ir)= 0
      irvec(3, ir)= 0
      hmnr(1, 1, ir)=-mx
      hmnr(2, 2, ir)=-mx

      !-1 0 0
      ir= ir+ 1
      irvec(1, ir)=-1
      irvec(2, ir)= 0
      irvec(3, ir)= 0
      hmnr(1, 1, ir)=-mx
      hmnr(2, 2, ir)=-mx

      ! 0 1 0
      ir= ir+ 1
      irvec(1, ir)= 0 
      irvec(2, ir)= 1
      irvec(3, ir)= 0
      hmnr(1, 1, ir)=-my
      hmnr(2, 2, ir)= my

      !0 -1  0
      ir= ir+ 1
      irvec(1, ir)= 0 
      irvec(2, ir)=-1
      irvec(3, ir)= 0
      hmnr(1, 1, ir)=-my
      hmnr(2, 2, ir)= my
 
      ! 0  0  1
      ir= ir+ 1
      irvec(1, ir)= 0 
      irvec(2, ir)= 0
      irvec(3, ir)= 1
      hmnr(1, 1, ir)=-mz
      hmnr(2, 2, ir)= mz

      ! 0  0 -1
      ir= ir+ 1
      irvec(1, ir)= 0 
      irvec(2, ir)= 0
      irvec(3, ir)=-1
      hmnr(1, 1, ir)=-mz
      hmnr(2, 2, ir)= mz

      nrpts = ir

      !> write to new_hr.dat
      open(unit=105, file='Free_compensate_model_hr.dat')
      write(105, *)'2-band electron-hole compensated model'
      write(105, *)nwann
      write(105, *)nrpts
      write(105, '(15I5)')(ndegen(i), i=1, nrpts)
      do ir=1, nrpts
         do i=1, nwann
            do j=1, nwann
               write(105, '(5I5, 2f16.8)')irvec(:, ir), i, j, HmnR(i, j, ir)
            enddo
         enddo
      enddo
      close(105)

   end ! end of program 
