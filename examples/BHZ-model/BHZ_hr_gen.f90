! > BHZ model
!> gfortran BHZ_hr_gen.f90 -o BHZ_hr_gen
!> ./BHZ_hr_gen
program main
      implicit none

      integer, parameter :: dp=kind(1d0)
      real(dp), parameter :: pi=atan(1d0)*4d0
      complex(dp), parameter :: zi=(0d0, 1d0)
      integer :: ir, nrpts, i, j, num_wann

      ! lattice constant
      real(dp), parameter :: A= 1d0
      real(dp), parameter :: B= 1d0
      real(dp), parameter :: Delta=  1.0d0
      real(dp), parameter :: alpha=  0.0d0
             
      ! R coordinates  
      integer, allocatable     :: Irvec(:,:)
      
      ! Hamiltonian m,n are band indexes
      complex(dp), allocatable :: HmnR(:,:,:)
      
      ! no of degeneracy of R point 
      integer, allocatable     :: ndegen(:)

      num_wann=4
      nrpts=7
      allocate(Irvec(3, nrpts))
      allocate(ndegen(nrpts))
      allocate(HmnR(num_wann, num_wann, nrpts))
 
      irvec=0
      ndegen=1
      hmnr= 0d0

      ! 0 0 0
      ir= 1
      irvec(1, ir)= 0
      irvec(2, ir)= 0
      hmnr(1, 1, ir)=  Delta- 4d0*B
      hmnr(2, 2, ir)= -Delta+ 4d0*B
      hmnr(3, 3, ir)=  Delta- 4d0*B
      hmnr(4, 4, ir)= -Delta+ 4d0*B

      ! 1 0 
      ir= 2
      irvec(1, ir)= 1
      irvec(2, ir)= 0
      hmnr(1, 1, ir)=   B
      hmnr(2, 2, ir)= - B
      hmnr(3, 3, ir)=   B
      hmnr(4, 4, ir)= - B
      hmnr(1, 3, ir)=   alpha/2d0
      hmnr(3, 1, ir)= - alpha/2d0
      hmnr(4, 2, ir)= - alpha/2d0
      hmnr(2, 4, ir)=   alpha/2d0
      hmnr(1, 2, ir)= 0.5d0*( alpha- zi*A)
      hmnr(2, 1, ir)= 0.5d0*(-alpha- zi*A)
      hmnr(3, 4, ir)= 0.5d0*( alpha+ zi*A)
      hmnr(4, 3, ir)= 0.5d0*(-alpha+ zi*A)

      ! 0 1 
      ir= 3
      irvec(1, ir)= 0
      irvec(2, ir)= 1
      hmnr(1, 1, ir)=   B
      hmnr(2, 2, ir)= - B
      hmnr(3, 3, ir)=   B
      hmnr(4, 4, ir)= - B
      hmnr(1, 3, ir)= - zi*alpha
      hmnr(3, 1, ir)= - zi*alpha
      hmnr(4, 2, ir)=   zi*alpha
      hmnr(2, 4, ir)=   zi*alpha
      hmnr(1, 2, ir)=  -A/2d0
      hmnr(2, 1, ir)=   A/2d0
      hmnr(3, 4, ir)=  -A/2d0
      hmnr(4, 3, ir)=   A/2d0

      !-1 0 
      ir= 4
      irvec(1, ir)=-1
      irvec(2, ir)= 0
      hmnr(1, 1, ir)=   B
      hmnr(2, 2, ir)= - B
      hmnr(3, 3, ir)=   B
      hmnr(4, 4, ir)= - B
      hmnr(1, 3, ir)= - alpha/2d0
      hmnr(3, 1, ir)=   alpha/2d0
      hmnr(4, 2, ir)=   alpha/2d0
      hmnr(2, 4, ir)= - alpha/2d0
      hmnr(1, 2, ir)= 0.5d0*(-alpha+ zi*A)
      hmnr(2, 1, ir)= 0.5d0*( alpha+ zi*A)
      hmnr(3, 4, ir)= 0.5d0*(-alpha- zi*A)
      hmnr(4, 3, ir)= 0.5d0*( alpha- zi*A)

      ! 0-1 
      ir= 5
      irvec(1, ir)= 0
      irvec(2, ir)=-1
      hmnr(1, 1, ir)=   B
      hmnr(2, 2, ir)= - B
      hmnr(3, 3, ir)=   B
      hmnr(4, 4, ir)= - B
      hmnr(1, 3, ir)=   zi*alpha
      hmnr(3, 1, ir)=   zi*alpha
      hmnr(4, 2, ir)= - zi*alpha
      hmnr(2, 4, ir)= - zi*alpha
      hmnr(1, 2, ir)=   A/2d0
      hmnr(2, 1, ir)=  -A/2d0
      hmnr(3, 4, ir)=   A/2d0
      hmnr(4, 3, ir)=  -A/2d0

      !> write to new_hr.dat
      open(unit=105, file='BHZ_hr.dat')
      write(105, *)'4-band BHZ model'
      write(105, *)'4'
      write(105, *)nrpts
      write(105, '(15I5)')(ndegen(i), i=1, nrpts)
      do ir=1, nrpts
         do i=1, 4
            do j=1, 4
               write(105, '(5I5, 2f16.8)')irvec(:, ir), i, j, HmnR(i, j, ir)
            enddo
         enddo
      enddo
      close(105)


end
