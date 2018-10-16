   subroutine Proteus(ndim, k,gap,ptol,func_gap,iter)
      ! This subroutine gives the local minimal of the func_gap in multidimensional
      ! parameters space
      ! Ref:  Nelder, J.A., and Mead, R. 1965, Computer Journal, vol. 7, pp. 308–313. [1]

      use para, only : dp, cpuid, stdout, eps9
      implicit none

      integer, intent(in) :: ndim  ! number of parameters
      integer, intent(out) :: iter ! number of iterations

      real(dp), intent(in) :: ptol  ! tolerance to be achieved in func_gap

      real(dp), intent(inout) :: gap(ndim+1)
      real(dp), intent(inout) :: k(ndim+1, ndim)

      ! function to calculate the energy gap for a given k point x
      interface
         function func_gap(n, x)
         implicit none
         integer, parameter :: dp= kind(1d0)
         integer, intent(in) :: n
         real(dp), dimension(n), intent(in) :: x
         real(dp) :: func_gap
         end function func_gap
      end interface

      integer, parameter :: itmax=500  ! maximum iterations
      real(dp), parameter :: tiny=1.0e-10
      integer :: iworst
      real(dp), dimension(size(k,2)) :: ksum

      call Proteus_0
      contains

      subroutine Proteus_0
      implicit none
      integer :: i,ibest,inworst
      real(dp) :: rtol,ysave,ytry,ytmp
      iter=0

      ! Enter here when starting or have just overall contracted.
      ! Recompute ksum.
      ksum(:)=sum(k(:,:),dim=1)

      do while (.true.)

         ! Enter here when have just changed a single point.
         ! Determine which point is the highest (worst) iworst, next-highest inworst,
         ! and lowest (best)  ibest
         ibest=iminloc(gap(:))
         iworst=imaxloc(gap(:))
         ytmp=gap(iworst)
         gap(iworst)=gap(ibest)
         inworst=imaxloc(gap(:))
         gap(iworst)=ytmp

         ! Compute the fractional range from highest to lowest and return if
         ! satisfactory.
         ! If returning, put best point and value in slot 1.
         rtol=2.0_dp*abs(gap(iworst)-gap(ibest))/(abs(gap(iworst))+abs(gap(ibest))+tiny)
         if (rtol < ptol.or.gap(ibest)<eps9) then
            call swap1(gap(1),gap(ibest))
            call swap2(k(1,:),k(ibest,:))
            call printkgap(k, gap)
            return
         end if
         if (iter >= itmax) then
            if (cpuid.eq.0) write(stdout, *)'itmax exceeded in Proteus'
            return
         endif
         
         ! start a new iteration. first extrapolate by a factor ?1 through the
         ! face of the simplex
         ! across from the high point, i.e., reflect the simplex from the high
         ! point.
         ytry=protry(-1.0_dp)
         iter=iter+1
         call printkgap(k, gap)
         if (ytry <= gap(ibest)) then
            ! gives a result better than the best point, so
            ! try an additional extrapolation by a factor of 2
            ytry=protry(2.0_dp)
            iter=iter+1
            call printkgap(k, gap)
         else if (ytry >= gap(inworst)) then
            ! the reflected point is worse than the second
            ! highest, so look for an intermediate lower point, i.e., do a
            ! one-dimensional contraction.
            ysave=gap(iworst)
            ytry=protry(0.5_dp)
            iter=iter+1
            call printkgap(k, gap)
            if (ytry >= ysave) then
               ! can't seem to get rid of that high point. better contract
               ! around the lowest (best) point.
               k(:,:)=0.5_dp*(k(:,:)+spread(k(ibest,:),1,size(k,1)))
               do i=1,ndim+1
                  if (i /= ibest) gap(i)=func_gap(ndim, k(i,:))
               end do
               !keep track of function evaluations.
               iter=iter+ndim 
               ksum(:)=sum(k(:,:),dim=1)
            end if
         end if

      end do !go back for the test of doneness and the next
   end subroutine Proteus_0

   function protry(fac)
      ! Extrapolates by a factor fac through the face of the simplex across from
      ! the high point,
      ! tries it, and replaces the high point if the new point is better.
      implicit none
      real(dp), intent(in) :: fac
      real(dp) :: protry
      real(dp) :: fac1,fac2,ytry
      real(dp), dimension(size(k,2)) :: ptry
      fac1=(1.0_dp-fac)/ndim
      fac2=fac1-fac
      ptry(:)=ksum(:)*fac1-k(iworst,:)*fac2
      ytry=func_gap(ndim, ptry)
      if (ytry < gap(iworst)) then
         gap(iworst)=ytry
         ksum(:)=ksum(:)-k(iworst,:)+ptry(:)
         k(iworst,:)=ptry(:)
      end if
      protry=ytry
   end function protry

   function iminloc(gap0)
   ! Find the location of the minimum value of y
      implicit none
      real(dp), intent(in) :: gap0(ndim+1)
      integer :: iminloc

      integer :: i

      iminloc=1
      do i=1, ndim+1
         if (gap0(i)<gap0(iminloc)) iminloc=i
      enddo
   end function iminloc

   function imaxloc(gap0)
   ! Find the location of the maximum value of gap
      implicit none
      real(dp), intent(in) :: gap0(ndim+1)
      integer :: imaxloc

      integer :: i

      imaxloc=1
      do i=1, ndim+1
         if (gap0(i)>gap0(imaxloc)) imaxloc=i
      enddo
   end function imaxloc

   subroutine swap1(x, y)
      ! swap x and y, real valued
      implicit none
      real(dp), intent(inout) :: x
      real(dp), intent(inout) :: y
      real(dp) :: z
      z= x
      x= y
      y= z
      return
   end subroutine swap1

   subroutine swap2(x, y)
      ! swap x and y, real vector valued
      implicit none
      real(dp), intent(inout) :: x(ndim)
      real(dp), intent(inout) :: y(ndim)
      real(dp) :: z(ndim)
      z= x
      x= y
      y= z

      return
   end subroutine swap2

   subroutine printkgap(k0, gap0)

      implicit none

      real(dp), intent(in) :: k0(ndim+1, ndim)
      real(dp), intent(in) :: gap0(ndim+1      )

      integer :: i
      
      if (mod(iter, 50).ne.1) return

      ! output some informations
      if (cpuid.eq.0) write(stdout, '(a, i8)') 'Proteus iter= ', iter
      if (cpuid.eq.0) write(stdout, '("k  ", 4f16.6)') (k0(1, i), i=1, ndim)
      if (cpuid.eq.0) write(stdout, '("gap", 4f16.6)') (gap0(   i), i=1, ndim+1)

      return
   end subroutine printkgap
   end subroutine Proteus
