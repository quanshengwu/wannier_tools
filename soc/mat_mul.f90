! performs matrix-matrix multiply
! C=A*B
! function mat_mul(nmatdim,A,B)
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
!  end function mat_mul
