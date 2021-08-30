  subroutine eigensystem_c(JOBZ,UPLO,N,A,W)
     ! A pack of Lapack subroutine zheev, which is 
     ! a subroutine to calculate eigenvector and eigenvalue for a 
     ! complex hermite matrix
     !> copied from WannierTools
     !> Author: Q.S Wu (wuquansheng@gmail.com)

     use para, only : Dp, stdout
     implicit none

!  JOBZ    (input) CHARACTER*1
!          = 'N':  Compute eigenvalues only;
!          = 'V':  Compute eigenvalues and eigenvectors.
!
     character*1, intent(in) :: JOBZ

!  UPLO    (input) CHARACTER*1
!          = 'U':  Upper triangle of A is stored;
!          = 'L':  Lower triangle of A is stored.
     character*1, intent(in) :: UPLO

!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
     integer,   intent(in) :: N

!  A       (input/output) COMPLEX*16 array, dimension (LDA, N)
!          On entry, the Hermitian matrix A.  If UPLO = 'U', the
!          leading N-by-N upper triangular part of A contains the
!          upper triangular part of the matrix A.  If UPLO = 'L',
!          the leading N-by-N lower triangular part of A contains
!          the lower triangular part of the matrix A.
!          On exit, if JOBZ = 'V', then if INFO = 0, A contains the
!          orthonormal eigenvectors of the matrix A.
!          If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')
!          or the upper triangle (if UPLO='U') of A, including the
!          diagonal, is destroyed.

     complex(Dp),intent(inout) :: A(N,N)

!  W       (output) DOUBLE PRECISION array, dimension (N)
!          If INFO = 0, the eigenvalues in ascending order.

     real(Dp), intent(inout) :: W(N)
    
     integer :: info
     integer :: lwork

     real(Dp),allocatable ::  rwork(:)

     complex(Dp),allocatable :: work(:)

     lwork=16*N
     allocate(rwork(lwork))
     allocate( work(lwork))
     rwork= 0d0
     work= (0d0, 0d0)

     info=0
     W=0.0d0

     !> if N==1, you don't have to do the diagonalization
     if (N==1) then 
        W=A(1, 1)
        A(1, 1)= 1d0
        return
     endif

     call zheev( JOBZ, UPLO, N, A, N,  &
              W, work, lwork, rwork, info )

     if (info.ne.0) then
        write(stdout, *) 'ERROR : something wrong with zheev'
        stop
     endif

     deallocate(rwork, work)
     return
  end subroutine eigensystem_c

  subroutine eigensystem_r (JOBZ,UPLO,N,A,W)
     ! A pack of Lapack subroutine dsyev, which is 
     ! a subroutine to calculate eigenvector and eigenvalue for a 
     ! real symmetric matrix

     use para, only : Dp, stdout
     implicit none

!  JOBZ    (input) CHARACTER*1
!          = 'N':  Compute eigenvalues only;
!          = 'V':  Compute eigenvalues and eigenvectors.
!
     character*1, intent(in) :: JOBZ

!  UPLO    (input) CHARACTER*1
!          = 'U':  Upper triangle of A is stored;
!          = 'L':  Lower triangle of A is stored.
     character*1, intent(in) :: UPLO

!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
     integer,   intent(in) :: N

!  A       (input/output) COMPLEX*16 array, dimension (LDA, N)
!          On entry, the Hermitian matrix A.  If UPLO = 'U', the
!          leading N-by-N upper triangular part of A contains the
!          upper triangular part of the matrix A.  If UPLO = 'L',
!          the leading N-by-N lower triangular part of A contains
!          the lower triangular part of the matrix A.
!          On exit, if JOBZ = 'V', then if INFO = 0, A contains the
!          orthonormal eigenvectors of the matrix A.
!          If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')
!          or the upper triangle (if UPLO='U') of A, including the
!          diagonal, is destroyed.

     real(Dp),intent(inout) :: A(N,N)

!  W       (output) DOUBLE PRECISION array, dimension (N)
!          If INFO = 0, the eigenvalues in ascending order.

     real(Dp), intent(inout) :: W(N)
    
     integer :: info

     real(Dp),allocatable :: work(:)

     allocate( work(16*N))
     work= 0d0

     if (N==1) then 
        W=A(1, 1)
        A(1, 1)= 1d0
        return
     endif

     info=0
     W=0.0d0
     call dsyev( JOBZ, UPLO, N, A, N,  &
              W, work, 16*N, info )

     if (info.ne.0) then
        write(stdout, *) 'ERROR : something wrong with dsyev'
        stop
     endif

     deallocate(work)
     return
  end subroutine eigensystem_r

  subroutine zgeev_sys(N, A, W,JOBVL,VL,JOBVR,VR )
     ! a pack of Lapack subroutine zgeev, which is a subroutine to
     ! calculate eigenvector and eigenvalue of a complex hermite matrix

     use para, only : Dp
     implicit none

!  JOBVL   (input) CHARACTER*1
!  JOBVR   (input) CHARACTER*1
!          = 'N':  Compute eigenvalues only;
!          = 'V':  Compute eigenvalues and eigenvectors.
!
     character*1 :: JOBVL
     character*1 :: JOBVR

!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
     integer,   intent(in) :: N

!~         A is COMPLEX*16 array, dimension (LDA,N)
!~           On entry, the N-by-N matrix A.
!~           On exit, A has been overwritten.

     complex(Dp),intent(inout) :: A(N, N)

!  W       (output) DOUBLE PRECISION array, dimension (N)
!          If INFO = 0, the eigenvalues in ascending order.

!    eigenvalues
     complex(Dp), intent(out) :: W(N)
    
!    left eigenvectors
     complex(dp) :: VL(N,N)

!    right eigenvectors
     complex(dp) :: VR(N,N)

     integer :: info

     integer :: lda
     integer :: ldvl
     integer :: ldvr

     integer :: lwork

     real(dp),allocatable ::  rwork(:)

     complex(dp),allocatable :: work(:)

     !> only calculate eigenvalues
!~      JOBVL= 'N'
!~      JOBVR= 'N'

     lda=N
     ldvl=N
     ldvr=N
     lwork= 16*N

     allocate(rwork(16*N))
     allocate( work(lwork))
     VL= (0d0, 0d0)
     VR= (0d0, 0d0)
     rwork= 0d0
     work= (0d0, 0d0)

     info=0
     W=(0d0, 0d0)

     if (N==1) then 
        W=A(1, 1)
        VL(1, 1)= 1d0
        VR(1, 1)= 1d0
        return
     endif

     call  ZGEEV( JOBVL, JOBVR, N, A, LDA, W, VL, LDVL, VR, LDVR, &
                        WORK, LWORK, RWORK, INFO )
     if (info /= 0) then
        stop ">>> Error : something wrong happens in zgeev_pack" 
     endif


     deallocate(rwork)
     deallocate( work)

     return
  end subroutine 

   subroutine zgeev_pack(N,A,W)
   !pass the vars into zgeev_sys, without calculating jobvl
    use para, only : Dp
     implicit none
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
     integer,   intent(in) :: N

!  A       (input/output) complex, 
!~       in-> the matrix complex(16)
!~       out-> the eigenvectors(if jobvl=V

     complex(Dp),intent(inout) :: A(N, N)
     complex(dp),allocatable :: vl(:,:),vr(:,:)

!  W       (output) DOUBLE PRECISION array, dimension (N)
!          If INFO = 0, the eigenvalues in ascending order.

!    eigenvalues
     complex(Dp), intent(out) :: W(N)  
     allocate(vl(N,N),vr(N,N))
     call zgeev_sys(N,A,W,'N',vl,'N',vr)
     deallocate(vr,vl)
     return
   
   end subroutine zgeev_pack

  !> the IL-th through IU-th eigenvalues will be found.
  subroutine zheevx_pack(JOBZ, UPLO, N, il, iu, A, eigval, eigvec)
     use para, only : dp
     implicit none
    
     !> inout variables
     character(1), intent(in) :: JOBZ ! 'V' compute eigenvectors and eigenvalues; 'N' only eigenvalues
     character(1), intent(in) :: UPLO
     integer, intent(in) :: N  !> dimension of matrix A
     integer, intent(in) :: il !> lowest band to be calculated
     integer, intent(in) :: iu !> highest band to be calculated
     real(dp), intent(inout) :: eigval(iu-il+1) !> eigenvalues
     complex(dp), intent(inout) :: A(N, N) !> the input hermite complex matrix
     complex(dp), intent(inout) :: eigvec(N, iu-il+1) !> eigenvectors

     !> local variables
     real(dp), allocatable :: eigenvalues(:)
     integer , allocatable :: iwork(:)
     integer , allocatable :: ifail(:)
     real(dp), allocatable :: rwork(:)
     complex(dp), allocatable :: work(:)

     integer :: mdim
     integer :: lwork
     integer :: info
     real(dp) :: vl  !> not referenced in this subroutine
     real(dp) :: vu  !> not referenced in this subroutine

     real(dp), parameter :: abstol= 1D-10

     !real(dp), external :: DLAMCH
     !abstol= 2*DLAMCH('S')

     mdim= iu-il+1
     lwork= 64*N
     vl=0d0; vu=0d0 

     allocate(ifail(N))
     allocate(iwork(5*N))
     allocate(rwork(7*N))
     allocate(work(lwork))
     allocate(eigenvalues(N))
     ifail = 0
     iwork = 0
     rwork = 0d0
     work = 0d0

     eigenvalues= 0d0
     eigvec= 0d0

     call zheevx(JOBZ,'I',UPLO,N,A,N,vl,vu,il,iu,abstol,&
         mdim,eigenvalues,eigvec,N,work,-1,rwork,iwork,ifail,info)
     lwork=int(work(1))

     call zheevx(JOBZ,'I',UPLO,N,A,N,vl,vu,il,iu,abstol,&
         mdim,eigenvalues,eigvec,N,work,lwork,rwork,iwork,ifail,info)


     if (info/=0) then
        print *, ' Error info in zheevx: ', info
        stop ' Error happens in zheev_pack'
     endif

     eigval(1:mdim)= eigenvalues(1:mdim)

     deallocate(eigenvalues, work, iwork, ifail, rwork)

     return
  end subroutine  zheevx_pack


