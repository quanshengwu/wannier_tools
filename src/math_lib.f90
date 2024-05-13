!> Here we put a set of basic mathematic operations
function det3(a)
   !> calculate the determinant of a real valued rank-3 matrix
   use para, only : dp
   implicit none
   real(dp) :: det3
   real(dp), intent(in) :: a(3, 3)

   det3= a(1,1)*(a(2,2)*a(3,3)-a(3,2)*a(2,3)) &
         +a(2,1)*(a(3,2)*a(1,3)-a(1,2)*a(3,3)) &
         +a(3,1)*(a(1,2)*a(2,3)-a(2,2)*a(1,3))

   return
end function det3

!> Get the angle between two given vectors which are in cartesian coordinates
!> return in degree
function angle(R1, R2)
   use para, only : dp, pi

   implicit none
   real(dp), intent(in) :: R1(3), R2(3)
   real(dp) :: angle, adotb
   real(dp), external :: norm
 
   adotb=dot_product(R1, R2)
   angle= acos(adotb/norm(R1)/norm(R2))*180d0/pi

   return
end function angle

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

!> shift the atom's position to the home unit cell
!> shift pos_direct_pc to the home unit cell (-0.5, 0.5]
subroutine in_home_cell_regularization(pos)
   use para, only : dp, eps3, pi

   implicit none

   real(dp), intent(inout) :: pos(3)

   integer :: i
   real(dp) :: irrational_shift(3)
   irrational_shift= (/ pi/1000d0, pi/1000d0, 0d0 /)
   pos= pos+ irrational_shift

   pos= pos-floor(pos)

   do i=1, 3
      if (abs(pos(i)-1)<0.03d0) pos(i)= 1d0
      if (abs(pos(i)-0.5)<0.03d0) pos(i)= 0.5d0
      if (pos(i)>0.5000000d0) pos(i)= pos(i)-1d0
   enddo
   pos= pos- irrational_shift

   return
end subroutine in_home_cell_regularization



!> the shortest difference betwerrn two vectors with respect to the lattice vectors
subroutine periodic_diff_1D(R2, R1, diff)
   !> diff= mod(R2-R1, 1)
   use para, only : dp

   implicit none

   real(dp), intent(in) :: R1, R2
   real(dp), intent(out) :: diff

   integer :: i

   diff= R2-R1
   diff= diff-floor(diff)

   if (diff>0.5000000d0) diff= diff-1d0

   return
end subroutine periodic_diff_1D

!> shift the atom's position to the home unit cell
!> shift pos_direct_pc to the home unit cell [-0.5, 0.5)
subroutine in_home_cell(R0)
   use para, only : dp

   implicit none

   real(dp), intent(inout) :: R0(3)

   integer :: i

   R0= R0-int8(R0)

   do i=1, 3
      if (R0(i)>=0.5000000d0) R0(i)= R0(i)-1d0
   enddo

   R0= R0+ 0.5d0
   R0= mod(R0, 1d0)

   return
end subroutine in_home_cell



!> the shortest difference betwerrn two vectors with respect to the lattice vectors
subroutine periodic_diff(R2, R1, diff)
   !> diff= mod(R2-R1, 1)
   use para, only : dp

   implicit none

   real(dp), intent(in) :: R1(3), R2(3)
   real(dp), intent(out) :: diff(3)

   integer :: i

   diff= R2-R1
   diff= diff-floor(diff)

   do i=1, 3
      if (diff(i)>0.5000000d0) diff(i)= diff(i)-1d0
   enddo

   return
end subroutine periodic_diff

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

Module Kronecker

contains

! Takes in Matrices A(i,j),B(k,l), assumed 2D, returns Kronecker Product
! C(i*k,j*l)
function KronProd(A,B) result(C)

   use para, only : Dp
   IMPLICIT NONE

   complex(Dp), dimension (:,:), intent(in)  :: A, B
   complex(Dp), dimension (:,:), allocatable :: C
   !real, dimension (:,:), intent(in)  :: A, B
   !real, dimension (:,:), allocatable :: C
   integer :: i = 0, j = 0, k = 0, l = 0
   integer :: m = 0, n = 0, p = 0, q = 0


   allocate(C(size(A,1)*size(B,1),size(A,2)*size(B,2)))
   C = 0

   do i = 1,size(A,1)
      do j = 1,size(A,2)
         n=(i-1)*size(B,1) + 1
         m=n+size(B,1) - 1
         p=(j-1)*size(B,2) + 1
         q=p+size(B,2) - 1
         C(n:m,p:q) = A(i,j)*B
      enddo
   enddo

end function KronProd
end module Kronecker