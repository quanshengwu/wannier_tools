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
