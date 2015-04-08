! Purpose :
!    Multidimensional minimization of the function funk(x) where !x(1:ndim)
!    is a vector in ndim dimensions, by the downhill simplex method of 
!    Nelder and Mead. The  matrix p(1:ndim+1,1:ndim) is input. Its ndim+1 rows 
!    are ndim-dimensional vectors which are the vertices of the starting 
!    simplex. Also input is the vector y(1:ndim+1), whose components must be 
!    pre-initialized to the values of funk evaluated at the ndim+1 vertices 
!    (rows) of p; and ftol the fractional convergence tolerance to be 
!    achieved in the  function value (n.b.!). On output, p and y will 
!    have been reset to ndim+1 new points all within ftol of a minimum
!    function value, and iter gives the number of function evaluations taken.
! Ref : 
!    Nelder, J.A., and Mead, R. 1965, Computer Journal, vol. 7, pp. 308–313. [1]
!    Yarbro, L.A., and Deming, S.N. 1974, Analytica Chimica Acta, vol. 73, pp.
!    391–398.
!    Jacoby, S.L.S, Kowalik, J.S., and Pizzo, J.T. 1972, Iterative Methods for
!    Nonlinear Optimization
!    Problems (Englewood Cliffs, NJ: Prentice-Hall).
!   
!    

   subroutine amoeba(ndim, p,y,ftol,func,iter)
      implicit none

      ! precision control 
      integer, parameter :: dp= kind(1d0)
      integer, parameter :: i4b= 4

      !>> inout variables
      ! dimension of our variable space
      integer(i4b), intent(in) :: ndim
      integer(i4b), intent(out) :: iter

      ! the fractional convergence tolerance to be achieved in the function
      ! value
      real(dp), intent(in) :: ftol

      real(dp), intent(inout) :: y(ndim+1)
      real(dp), intent(inout) :: p(ndim+1, ndim)

      interface
         function func(n, x)
         implicit none
         integer, parameter :: dp= kind(1d0)
         integer, intent(in) :: n
         real(dp), dimension(n), intent(in) :: x
         real(dp) :: func
         end function func
      end interface

      !>> local variables
      integer(i4b), parameter :: itmax=5000
      real(dp), parameter :: tiny=1.0e-10
      integer(i4b) :: ihi
      real(dp), dimension(size(p,2)) :: psum

      call amoeba_private
      contains

      subroutine amoeba_private
      implicit none
      integer(i4b) :: i,ilo,inhi
      real(dp) :: rtol,ysave,ytry,ytmp
      iter=0

      ! Enter here when starting or have just overall contracted.
      ! Recompute psum.
      psum(:)=sum(p(:,:),dim=1)

      do while (.true.)

         ! Enter here when have just changed a single point.
         ! Determine which point is the highest (worst) ihi, next-highest inhi,
         ! and lowest (best)  ilo
         ilo=iminloc(y(:))
         ihi=imaxloc(y(:))
         ytmp=y(ihi)
         y(ihi)=y(ilo)
         inhi=imaxloc(y(:))
         y(ihi)=ytmp

         ! Compute the fractional range from highest to lowest and return if
         ! satisfactory.
         ! If returning, put best point and value in slot 1.
         rtol=2.0_dp*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo))+tiny)
         if (rtol < ftol) then
            call swap1(y(1),y(ilo))
            call swap2(p(1,:),p(ilo,:))
            call printpy(p, y)
            print *, 'here'
            return
         end if
         if (iter >= itmax) then
            print*, 'itmax exceeded in amoeba'
            return
         endif
         
         ! begin a new iteration. first extrapolate by a factor ?1 through the
         ! face of the simplex
         ! across from the high point, i.e., reflect the simplex from the high
         ! point.
         ytry=amotry(-1.0_dp)
         iter=iter+1
         call printpy(p, y)
         if (ytry <= y(ilo)) then
            ! gives a result better than the best point, so
            ! try an additional extrapolation by a factor of 2
            ytry=amotry(2.0_dp)
            iter=iter+1
            call printpy(p, y)
         else if (ytry >= y(inhi)) then
            ! the reflected point is worse than the second
            ! highest, so look for an intermediate lower point, i.e., do a
            ! one-dimensional contraction.
            ysave=y(ihi)
            ytry=amotry(0.5_dp)
            iter=iter+1
            call printpy(p, y)
            if (ytry >= ysave) then
               ! can't seem to get rid of that high point. better contract
               ! around the lowest (best) point.
               p(:,:)=0.5_dp*(p(:,:)+spread(p(ilo,:),1,size(p,1)))
               do i=1,ndim+1
                 !if (i /= ilo) call func(ndim, p(i,:), y(i))
                  if (i /= ilo) y(i)=func(ndim, p(i,:))
               end do
               !keep track of function evaluations.
               iter=iter+ndim 
               psum(:)=sum(p(:,:),dim=1)
            end if
         end if

      end do !go back for the test of doneness and the next
   end subroutine amoeba_private
!bl

   !>>> Extrapolates by a factor fac through the face of the simplex across from
   ! the high point,
   ! tries it, and replaces the high point if the new point is better.
   function amotry(fac)
      implicit none
      real(dp), intent(in) :: fac
      real(dp) :: amotry
      real(dp) :: fac1,fac2,ytry
      real(dp), dimension(size(p,2)) :: ptry
      fac1=(1.0_dp-fac)/ndim
      fac2=fac1-fac
      ptry(:)=psum(:)*fac1-p(ihi,:)*fac2
      ytry=func(ndim, ptry)
     !call func(ndim, ptry, ytry)
      if (ytry < y(ihi)) then
         y(ihi)=ytry
         psum(:)=psum(:)-p(ihi,:)+ptry(:)
         p(ihi,:)=ptry(:)
      end if
      amotry=ytry
   end function amotry

   !>>> function that find the location of the minimum value of y
   function iminloc(y)
      implicit none
      real(dp), intent(in) :: y(ndim+1)
      integer(i4b) :: iminloc

      integer(i4b) :: i

      iminloc=1
      do i=1, ndim+1
         if (y(i)<y(iminloc)) iminloc=i
      enddo
   end function iminloc

   !>>> function that find the location of the maximum value of y
   function imaxloc(y)
      implicit none
      real(dp), intent(in) :: y(ndim+1)
      integer(i4b) :: imaxloc

      integer(i4b) :: i

      imaxloc=1
      do i=1, ndim+1
         if (y(i)>y(imaxloc)) imaxloc=i
      enddo
   end function imaxloc

   ! swap two values
   subroutine swap1(x, y)
      implicit none
      real(dp), intent(inout) :: x
      real(dp), intent(inout) :: y
      real(dp) :: z
      z= x
      x= y
      y= z
      return
   end subroutine swap1

   ! swap two vectors
   subroutine swap2(x, y)
      implicit none
      real(dp), intent(inout) :: x(ndim)
      real(dp), intent(inout) :: y(ndim)
      real(dp) :: z(ndim)
      z= x
      x= y
      y= z

      return
   end subroutine swap2

   subroutine printpy(p, y)

      implicit none

      real(dp), intent(in) :: p(ndim+1, ndim)
      real(dp), intent(in) :: y(ndim+1      )

      integer :: i

      ! output some informations
      write(6, '(a, i8)') 'Amoeba iter= ', iter
      write(6, '("x", 4f16.6)') (p(1, i), i=1, ndim)
      write(6, '("y", 4f16.6)') (y(   i), i=1, ndim+1)

      return
   end subroutine printpy
   end subroutine amoeba
