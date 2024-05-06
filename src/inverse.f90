  subroutine inv(ndim,Amat)


     implicit none

     integer,parameter :: dp=8

     integer           :: i
     integer           :: info

     integer,intent(in):: ndim

!    IPIV   : INTEGER. Array, DIMENSION at least max(1, n). The pivot indices that define
!    the permutation matrix P; row i of the matrix was interchanged with
!    row ipiv(i). Corresponds to the single precision factorization (if info=
!    0 and iter ≥ 0) or the double precision factorization (if info= 0 and
!    iter < 0).
     integer,allocatable   :: ipiv(:)


     complex(dp),parameter :: zone=(1.0d0,0.0d0)



!    Amat  :
!    Overwritten by the factors L and U from the factorization of A = P*L*U;
!    the unit diagonal elements of L are not stored.
!    If iterative refinement has been successfully used (info= 0 and iter≥
!    0), then A is unchanged.
!    If double precision factorization has been used (info= 0 and iter <
!    0), then the array A contains the factors L and U from the factorization
!    A = P*L*U; the unit diagonal elements of L are not stored.
     complex(dp),intent(inout):: Amat(ndim,ndim)

!    Bmat  :
!    Overwritten by the solution matrix X for dgesv, sgesv,zgesv,zgesv.
     complex(dp),allocatable :: Bmat(:,:)


     allocate(ipiv(ndim))
     allocate(Bmat(ndim,ndim))

     ipiv=0

     ! unit matrix
     Bmat= (0d0, 0d0)
     do i=1,ndim
        Bmat(i,i)= zone
     enddo

     call zgesv(ndim,ndim,Amat,ndim,ipiv,Bmat,ndim,info)

     if(info.ne.0)print *,'something wrong with zgesv'

     Amat=Bmat
     
     return
  end subroutine inv 

  !> get the inverse of a real matrix
  subroutine inv_r(ndim,Amat)

     implicit none

     integer,parameter :: dp=8

     integer           :: i
     integer           :: info

     integer,intent(in):: ndim

!    IPIV   : INTEGER. Array, DIMENSION at least max(1, n). The pivot indices that define
!    the permutation matrix P; row i of the matrix was interchanged with
!    row ipiv(i). Corresponds to the single precision factorization (if info=
!    0 and iter ≥ 0) or the double precision factorization (if info= 0 and
!    iter < 0).
     integer,allocatable   :: ipiv(:)


!    Amat  :
!    Overwritten by the factors L and U from the factorization of A = P*L*U;
!    the unit diagonal elements of L are not stored.
!    If iterative refinement has been successfully used (info= 0 and iter≥
!    0), then A is unchanged.
!    If double precision factorization has been used (info= 0 and iter <
!    0), then the array A contains the factors L and U from the factorization
!    A = P*L*U; the unit diagonal elements of L are not stored.
     real(dp),intent(inout):: Amat(ndim,ndim)

!    Bmat  :
!    Overwritten by the solution matrix X for dgesv, sgesv,zgesv,zgesv.
     real(dp),allocatable :: Bmat(:,:)


     allocate(ipiv(ndim))
     allocate(Bmat(ndim,ndim))

     ipiv=0

     ! unit matrix
     Bmat= 0d0 
     do i=1,ndim
        Bmat(i,i)= 1d0
     enddo

     call dgesv(ndim,ndim,Amat,ndim,ipiv,Bmat,ndim,info)

     if(info.ne.0)stop 'something wrong with dgesv in inverse.f90'

     Amat=Bmat
     
     return
  end subroutine inv_r


  subroutine sparse_solver(ndims, nnz, acsr, icsr, jcsr, b, x)
     !> this subroutine  solves A*x= b, where A is a sparse matrix stored in compress row format
     !> A is a Hermitian complex square matrix

     use para, only : dp
     implicit none
    
     !> the dimension of a square matrix
     integer, intent(in) :: ndims

     !> number of non-zeros elements
     integer, intent(in) :: nnz

     !> sparse matrix 
     ! acsr(K) = value of entry,
     ! icsr(K) = row of entry,
     ! jcsr(K) = column of entry.
     complex(dp), intent(in) :: acsr(nnz)
     integer, intent(in) :: icsr(ndims+1)
     integer, intent(in) :: jcsr(nnz)

     complex(dp), intent(out) :: x(ndims)
     complex(dp), intent(in)  :: b(ndims)
     real(dp), allocatable :: acsr_real(:)
     real(dp), allocatable :: acsr_imag(:)
     real(dp), allocatable :: x_real(:)
     real(dp), allocatable :: x_imag(:)
     real(dp), allocatable :: b_real(:)
     real(dp), allocatable :: b_imag(:)

     complex(dp), allocatable ::  acsc(:)
     integer, allocatable ::  icsc(:)
     integer, allocatable ::  jcsc(:)

     ! ITR_MAX, the maximum number of (outer) iterations to take.
     integer :: itr_max
 
     !MR, the maximum number of (inner) iterations to take.  0 < MR <= N.
     integer :: mr

     !TOL_ABS, an absolute tolerance applied to the current residual.
     !TOL_REL, a relative tolerance comparing the current residual to the initial residual.
     real(dp) :: tol_abs, tol_rel

     integer :: iopt, info, nrhs, i, ipos, ljob
     integer*8 :: factors


     !>--------------------------------------------------
     !>> call mgmres_st
     !> mgmres_st() applies restarted GMRES to a sparse triplet matrix.

     !allocate(acsr_real(nnz), acsr_imag(nnz))
     !allocate(x_real(ndims), x_imag(ndims))
     !allocate(b_real(ndims), b_imag(ndims))

     !acsr_real= real(acsr)
     !acsr_imag= aimag(acsr)
     !b_real= real(b)
     !b_imag= aimag(b)


     !itr_max=1000
     !tol_abs=1E-7
     !tol_rel=1E-7
     !mr= min(100, ndims)
    !call mgmres_csr ( ndims, nnz, icsr, jcsr, acsr_real, x_real, b_real, itr_max, mr, tol_abs, tol_rel )
    !call mgmres_csr ( ndims, nnz, icsr, jcsr, acsr_imag, x_imag, b_imag, itr_max, mr, tol_abs, tol_rel )
    
    !do i=1, ndims
    !   x(i) = cmplx(x_real(i), x_imag(i))
    !enddo
     !>--------------------------------------------------

     !>--------------------------------------------------
     !>> call superLU

#ifdef SUPERLU
     allocate( acsc(nnz), icsc(ndims+1), jcsc(nnz))
     !> convert csr to csc
     ljob = 1 ! fill c
     ipos = 1
     call csrcsc ( ndims, ljob, ipos, acsr, jcsr, icsr, acsc, jcsc, icsc )
     !> 
     x= b
! First, factorize the matrix. The factors are stored in *factors* handle.
      iopt = 1
      nrhs = 1
      call c_fortran_zgssv( iopt, ndims, nnz, nrhs, acsc, jcsc, icsc, &
           x, ndims, factors, info )
      
      if (info .eq. 0) then
      !  write (*,*) 'Factorization succeeded'
      else
      !  write(*,*) 'INFO from factorization = ', info
      endif
      
! Second, solve the system using the existing factors.
      iopt = 2
      call c_fortran_zgssv( iopt, ndims, nnz, nrhs, acsc, jcsc, icsc, &
                           x, ndims, factors, info )

      if (info .eq. 0) then
      !  write (*,*) 'Solve succeeded'
      else
      !  write(*,*) 'INFO from triangular solve = ', info
      endif

! Last, free the storage allocated inside SuperLU
      iopt = 3
      call c_fortran_zgssv( iopt, ndims, nnz, nrhs, acsc, jcsc, icsc, &
                           x, ndims, factors, info )
     !>--------------------------------------------------
#elif defined (INTELMKL)
!     !>--------------------------------------------------
     !>> call mkldss zgesv
      call zmat_mkldss_zgesv(ndims, nnz, acsr, jcsr, icsr, b, x)
     !>--------------------------------------------------
#else
     STOP "ERROR: Please add -DINTELMKL or -DSUPERLU in the Makefile to run the sparse_solver"
#endif


     return
  end subroutine sparse_solver

#if defined (INTELMKL)
      subroutine zmat_mkldss_zgesv(ndims, nnzmax, acsr, jcsr, icsr, xvec, yvec)
         use mkl_dss
         use para, only : dp, stdout, time_cost_t1, time_cost_t2, time_cost_t3

         implicit none

         ! external arguments
         ! dimension of matrix A
         integer, intent(in) :: ndims

         ! non-zero elements of matrix A
         integer, intent(in) :: nnzmax

         ! compressed sparse row storage of matrix A
         complex(dp), intent(in) :: acsr(nnzmax)
         integer, intent(in) :: jcsr(nnzmax)
         integer, intent(in) :: icsr(ndims+1)

         ! xvec vector with dimension "ndims"
         complex(dp), intent(in) :: xvec(ndims)

         ! yvec vector with dimension "ndims"
         complex(dp), intent(out) :: yvec(ndims)

         ! local variables
         type(mkl_dss_handle) :: handle
         integer :: nrhs
         integer :: ierr
         integer :: perm(1)
         real(dp) :: time_tic, time_toc

         !>>> step 0: setup the problem to be solved
         nrhs = 1
         perm(1) = 0

         !>>> step 1: initialize the solver.
         ierr = dss_create( handle, mkl_dss_defaults )
         if (ierr .ne. mkl_dss_success) then
            write(*,*) "mkl_dss_create returned error code", ierr; stop
         endif 

         !>>> step 2: define the non-zero structure of the matrix.
         ierr = dss_define_structure( handle, mkl_dss_symmetric_structure_complex, icsr, ndims, ndims, jcsr, nnzmax )
         !ierr = dss_define_structure( handle, MKL_DSS_NON_SYMMETRIC_COMPLEX, icsr, ndims, ndims, jcsr, nnzmax )
         if (ierr .ne. mkl_dss_success) then
            write(*,*) "mkl_dss_define_structure returned error code", ierr; stop
         endif ! back if (ierr .ne. mkl_dss_success} block

         call now(time_tic)
         !>>> step 3: reorder the matrix.
         ierr = dss_reorder( handle, mkl_dss_defaults, perm )
         if (ierr .ne. mkl_dss_success) then
            write(*,*) "mkl_dss_reorder returned error code", ierr; stop
         endif ! back if (ierr .ne. mkl_dss_success} block
         call now(time_toc)
         time_cost_t1= time_cost_t1+ time_toc- time_tic

         call now(time_tic)
         !>>> step 4: factor the matrix.
         ierr = dss_factor_complex( handle, mkl_dss_defaults, acsr )
         if (ierr .ne. mkl_dss_success) then
            write(*,*) "mkl_dss_factor_complex returned error code", ierr; stop
         endif ! back if (ierr .ne. mkl_dss_success} block
         call now(time_toc)
         time_cost_t2= time_cost_t2+ time_toc- time_tic

         call now(time_tic)
         ! allocate the solution vector and solve the problem.
         ierr = dss_solve_complex( handle, mkl_dss_defaults, xvec, nrhs, yvec )
         call now(time_toc)
         time_cost_t3= time_cost_t3+ time_toc- time_tic

         if (ierr .ne. mkl_dss_success) then
            write(*,*) "mkl_dss_solve_complex returned error code", ierr; stop
         endif ! back if (ierr .ne. mkl_dss_success} block

         ! deallocate solver storage and various local arrays.
         ierr = dss_delete( handle, mkl_dss_defaults )

         if (ierr .ne. mkl_dss_success) then
            write(*,*) "mkl_dss_solver returned error code", ierr; stop
         endif ! back if (ierr .ne. mkl_dss_success} block

         return
      end subroutine zmat_mkldss_zgesv
#endif

subroutine csrcsc ( n, job, ipos, a, ja, ia, ao, jao, iao )
 
!*****************************************************************************80
!
!! CSRCSC converts Compressed Sparse Row to Compressed Sparse Column.
!
!  Discussion:
!
!    This is essentially a transposition operation.  
!
!    It is NOT an in-place algorithm.
!
!    This routine transposes a matrix stored in a, ja, ia format.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) JOB, indicates whether or not to fill the values of the
!    matrix AO or only the pattern (IA, and JA).  Enter 1 for yes.
!
! ipos  = starting position in ao, jao of the transposed matrix.
!         the iao array takes this into account (thus iao(1) is set to ipos.)
!         Note: this may be useful if one needs to append the data structure
!         of the transpose to that of A. In this case use
!                call csrcsc (n,1,n+2,a,ja,ia,a,ja,ia(n+2))
!        for any other normal usage, enter ipos=1.
!
!    Input, real A(*), integer ( kind = 4 ) JA(*), IA(N+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Output, real AO(*), JAO(*), IAO(N+1), the matrix in CSC
!    Compressed Sparse Column format.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(*)
  complex ( kind = 8 ) ao(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(n+1)
  integer ( kind = 4 ) iao(n+1)
  integer ( kind = 4 ) ipos
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ja(*)
  integer ( kind = 4 ) jao(*)
  integer ( kind = 4 ) job
  integer ( kind = 4 ) k
  integer ( kind = 4 ) next
!
!  Compute lengths of rows of A'.
!
  iao(1:n+1) = 0

  do i = 1, n
    do k = ia(i), ia(i+1)-1
      j = ja(k) + 1
      iao(j) = iao(j) + 1
    end do
  end do
!
!  Compute pointers from lengths.
!
  iao(1) = ipos
  do i = 1, n
    iao(i+1) = iao(i) + iao(i+1)
  end do
!
!  Do the actual copying.
!
  do i = 1, n
    do k = ia(i), ia(i+1)-1
      j = ja(k)
      next = iao(j)
      if ( job == 1 ) then
        ao(next) = a(k)
      end if
      jao(next) = i
      iao(j) = next + 1
    end do
  end do
!
!  Reshift IAO and leave.
!
  do i = n, 1, -1
    iao(i+1) = iao(i)
  end do
  iao(1) = ipos

  return
end


