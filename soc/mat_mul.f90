! performs matrix-matrix multiply
! C=A*B
  subroutine mat_mul(nmatdim,A,B,C)
     
     use para, only : Dp
      implicit none


     integer,intent(in) :: nmatdim    

     complex(Dp) :: ALPHA
     complex(Dp) :: BETA 
 

     complex(Dp), intent(in)  :: A(nmatdim ,nmatdim)
     complex(Dp), intent(in)  :: B(nmatdim ,nmatdim)
     !complex(Dp) :: mat_mul(nmatdim,nmatdim)
     complex(Dp), intent(out) :: C(nmatdim,nmatdim)

     ALPHA=1.0d0 
     BETA=0.0D0

     C(:,:)=(0.0d0,0.0d0)

     call ZGEMM('N','N',nmatdim,nmatdim,nmatdim,ALPHA, &
               &  A,nmatdim,B,nmatdim,BETA,C,nmatdim)

     return
  end subroutine mat_mul

  !> ZGESVD computes the singular value decomposition (SVD) for GE matrices
  !> In this pack, we assume the matrix A is a square matrix, the dimension 
  !> of row and column are the same
  !> A = U * SIGMA * conjugate-transpose(V)
  !> VT= conjugate-transpose(V)
  subroutine zgesvd_pack(M, A, U, S, VT)

     use para, only : Dp
     implicit none

     integer, intent(in) :: M
     complex(dp), intent(inout) :: A(M, M)
     complex(dp), intent(out) :: U(M, M)
     real(dp)   , intent(out) :: S(M, M)
     complex(dp), intent(out) :: VT(M, M)

     character :: JOBU
     character :: JOBVT
     integer :: N
     integer :: LDA
     integer :: LDU
     integer :: LDVT
     integer :: LWORK
     complex(dp), allocatable :: WORK(:)
     real(dp), allocatable :: RWORK(:)
     integer :: INFO

     N= M
     LDA= M
     LDU= M
     LDVT= M
     allocate(RWORK(5*M))

     allocate(work(5*M))

     JOBU= 'A'
     JOBVT= 'A'

     LWORK = -1
     call zgesvd (JOBU, JOBVT, M, N, A, LDA, S, U, LDU, &
        VT, LDVT, WORK, LWORK, RWORK, INFO)
     if (INFO==0 .and. real(WORK(1))>0 )then
        LWORK= WORK(1)
        deallocate(work)
        allocate(WORK(LWORK))
     else
        write(*, *)'something wrong with zgesvd'
     endif


     call zgesvd (JOBU, JOBVT, M, N, A, LDA, S, U, LDU, &
        VT, LDVT, WORK, LWORK, RWORK, INFO)
     if (INFO /= 0) write(*, *)'something wrong with zgesvd'

     return
  end subroutine zgesvd_pack
