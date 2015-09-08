! The arguments for this subroutine have the same significance as for VA05A/AD, so that F, X, are one dimensional arrays. This subroutine must be provided by the user and must calculate the values of the functions f (x ,x ,...x ), j 1 2 N j=1,2,...,m for the given variables x = (x ,x ,...,x ) . These function values must be planted in the array F 

   function func_old(N,X)
      use para
      implicit none

      ! no. of parameters
      integer,intent(in) :: N

      ! parameters
      real(Dp),intent(inout):: X(N)

      real(dp) :: func_old

      integer :: i
      ! eigenvalue for each k point
      integer :: ik
      real(dp) :: k(3)

      real(dp) :: efermi_nsoc
      real(dp), allocatable :: W(:)
      complex(dp), allocatable :: Hamk_nsoc(:, :)

      lambda_p(5:12)= X(1)
      lambda_d(1:4)= X(2)
      efermi_nsoc= X(3)

      !> add spin-orbital term onto the nsoc HmnR
      !> for a given spin orbitall coupling strength lambda_p
      call addsoc

      allocate(W(num_wann_soc))
      allocate(Hamk_nsoc(num_wann_soc, num_wann_soc))

      func_old= 0d0
      eigval_nsoc= 0d0
      do ik=1, knv3
         k= kpoints(:, ik)
         W=0d0
         Hamk_nsoc= 0d0
         call Hamk_bulk_nsoc(k, Hamk_nsoc)
         call eigensystem_c('N', 'U', num_wann_soc, Hamk_nsoc, W)
         eigval_nsoc(:, ik)= W

         do i=1, num_wann_soc
            func_old= func_old+ abs((eigval_soc(i, ik)- W(i)+ efermi_nsoc))**2
         enddo
      enddo
      func_old= func_old/knv3

      deallocate(W, Hamk_nsoc)

      return
   end function func_old

! The arguments for this subroutine have the same significance as for VA05A/AD, so that F, X, are one dimensional arrays. This subroutine must be provided by the user and must calculate the values of the functions f (x ,x ,...x ), j 1 2 N j=1,2,...,m for the given variables x = (x ,x ,...,x ) . These function values must be planted in the array F 

   function func(N,X)
      use para
      implicit none

      ! no. of parameters
      integer,intent(in) :: N

      ! parameters
      real(Dp),intent(inout):: X(N)
      real(dp) :: func


      integer :: i
      ! eigenvalue for each k point
      integer :: ik
      real(dp) :: k(3)

      real(dp) :: efermi_nsoc
      real(dp), allocatable :: W(:)
      complex(dp), allocatable :: Hamk_nsoc(:, :)

     !X(1)= 0
      lambda_p(1:Num_atom_type)= X(1:Num_atom_type)
      lambda_d(1:Num_atom_type)= X(1+Num_atom_type:2*Num_atom_type)
      efermi_nsoc= X(Num_atom_type*2+ 1)

      !> add spin-orbital term onto the nsoc HmnR
      !> for a given spin orbitall coupling strength lambda_p
      call addsoc_all
      !call addsoc_pd_zjw

      allocate(W(num_wann_soc))
      allocate(Hamk_nsoc(num_wann_soc, num_wann_soc))

      func= 0d0
      eigval_nsoc= 0d0
      do ik=1, knv3
     !do ik=20, 40
         k= kpoints(:, ik)
         W=0d0
         Hamk_nsoc= 0d0
         call Hamk_bulk_nsoc(k, Hamk_nsoc)
         call eigensystem_c('N', 'U', num_wann_soc, Hamk_nsoc, W)
         eigval_nsoc(:, ik)= W

         do i=1, num_wann_soc-4
        !do i= 51, 60
        !do i= 55, 58
        !do i= 1, 60
            func= func+ abs((eigval_soc(i, ik)- W(i)+ efermi_nsoc))**2 &
               * weight(i, ik)
         enddo
      enddo
 
      func= func/knv3

      deallocate(W, Hamk_nsoc)

      return
   end function func

