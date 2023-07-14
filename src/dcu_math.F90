# if defined (DCU)
!> To generate dcu_math.o
!> hipfc -lhipblas -lhipsolver -lhipsparse -g -c dcu_math.F90
!> written by QSWu on Nov 18 2022 at IOP CAS, Beijing
!> Email: quansheng.wu@iphy.ac.cn or wuquansheng@gmail.com
!> tested with compiler/rocm/dtk-22.04.2 
!> hipfort was compiled with gcc and gfortran

!> A, B and C are double complex matrices, whose dimensionality is (ndim, ndim).
subroutine mat_mul_dcu_in(ndim, A, B, C)
   use dcu_matmul
   implicit none

   integer, intent(in) :: ndim
   complex(dp), intent(in), target  :: A(ndim, ndim)
   complex(dp), intent(in), target  :: B(ndim, ndim)
   complex(dp), intent(out), target :: C(ndim, ndim)

   integer :: istatus
   integer*8 :: isize
   complex(c_double_complex) :: alpha, beta

   alpha= 1d0
   beta = 0d0

   !> copy data from host to device
   !> d_A= A; d_B= B;
   istatus= hipMemcpy(d_A_mm, c_loc(A), sizeof(A), hipMemcpyHostToDevice)
   istatus= hipMemcpy(d_B_mm, c_loc(B), sizeof(A), hipMemcpyHostToDevice)

   istatus= hipblasZgemm(handle_mm, HIPBLAS_OP_N, HIPBLAS_OP_N, ndim, ndim, ndim, &
      alpha, d_A_mm, ndim, d_B_mm, ndim, beta, d_C_mm, ndim)

   !> copy data from device to host
   !> C= d_C
   call hipCheck(hipMemcpy(c_loc(C), d_C_mm, sizeof(A), hipMemcpyDeviceToHost))

   return
end subroutine mat_mul_dcu_in

!> A, B and C are double complex matrices, whose dimensionality is (ndim, ndim).
subroutine mat_mul_dcu_z(ndim, A, B, C)
   use iso_c_binding
   use hipfort
   use hipfort_check
   use hipfort_hipblas

   implicit none

   integer, parameter :: dp= 8
   integer, intent(in) :: ndim
   complex(dp), intent(in), dimension(ndim,ndim)  :: A
   complex(dp), intent(in), dimension(ndim,ndim)  :: B
   complex(dp), intent(out), dimension(ndim,ndim) :: C

   integer :: istatus
   integer*8 :: isize
   integer(kind(HIPBLAS_OP_N)), parameter :: transa = HIPBLAS_OP_N, transb =HIPBLAS_OP_N;
   complex(c_double_complex) :: alpha, beta

   !> variables on GPU device
   type(c_ptr) :: handle= c_null_ptr
   complex(dp), pointer, dimension(:, :) :: d_A
   complex(dp), pointer, dimension(:, :) :: d_B
   complex(dp), pointer, dimension(:, :) :: d_C

   alpha= 1d0
   beta = 0d0

   !> copy data from host to device
  !allocate(d_A(ndim, ndim), d_B(ndim, ndim), d_C(ndim, ndim))
   call hipCheck(hipMalloc(d_A, source=A))
   call hipCheck(hipMalloc(d_B, source=B))
   call hipCheck(hipMalloc(d_C, source=C))
  !istatus= hipMalloc(d_C, sizeof(C))

   !> d_A= A; d_B= B;
  !istatus= hipMemcpy(d_A, (A), hipMemcpyHostToDevice)
  !istatus= hipMemcpy(d_B, (B), hipMemcpyHostToDevice)

   istatus= hipblasCreate(handle)
   istatus= hipblasZgemm(handle, HIPBLAS_OP_N, HIPBLAS_OP_N, ndim, ndim, ndim, &
      alpha, d_A, ndim, d_B, ndim, beta, d_C, ndim)

!call hipblasCheck(hipblasCgemm(handle,transa,transb,m,n,k,alpha,da,size(da,1),db,size(db,1),beta,dc,size(dc,1)))

   !> copy data from device to host
   !> C= d_C
   call hipCheck(hipMemcpy(C, d_C, hipMemcpyDeviceToHost))

   !> Free memory on device
   istatus= hipblasDestroy(handle)
   !deallocate(d_A, d_B, d_C)
   istatus= hipFree(d_A)
   istatus= hipFree(d_B)
   istatus= hipFree(d_C)

   return
end subroutine mat_mul_dcu_z

!> A subroutine computes eigenvalues and eigenvectors of a symmetric or
!> Hermitian ndim x ndim matrix Amat(ndim, ndim)
!> input/output : Amat(ndim, ndim), Amat will be overwritten as the eigenvectors on output
!> use hipsolver on GPU device
subroutine hipsolver_zheev(ndim, Amat, eigval)

   use hipfort
   use hipfort_check
   use hipfort_hipblas
   use hipfort_hipsolver
   implicit none

   integer, parameter :: dp=kind(1d0)
   integer, intent(in) :: ndim
   real(dp), intent(out) :: eigval(ndim)
   complex(dp), intent(inout) :: Amat(ndim, ndim)

   integer :: istatus, info

   !> variables related to GPU
   type(c_ptr) :: handle
   integer(c_int), pointer :: devinfo
   complex(dp), pointer, dimension(:, :) :: d_A
   real(dp), pointer, dimension(:) :: d_eigval
   type(c_ptr) :: d_workspace
   integer :: lwork
   integer(c_size_t) :: lwork2

   integer(kind(HIPSOLVER_EIG_MODE_VECTOR)) :: jobz
   integer(kind(HIPBLAS_FILL_MODE_UPPER))   :: uplo

   jobz= HIPSOLVER_EIG_MODE_VECTOR  ! compute eigenvalues and eigenvectors.
   uplo= HIPBLAS_FILL_MODE_UPPER

   !> step 1: create hipsolver handle 
   istatus= hipsolverCreate(handle)
   
  
   !> step 2: allocate memory on device and copy data from host to device
   !allocate(d_A(ndim, ndim), d_eigval(ndim))
   istatus= hipMalloc(devinfo)
   istatus= hipMalloc(d_A, source=Amat)
   eigval= 0d0
   call hipCheck(hipMalloc(d_eigval, source=eigval))

   !> d_A= Amat
  !istatus= hipMemcpy(d_A, Amat, hipMemcpyHostToDevice)

   !> step 3: query working space of getrf 
   istatus= hipsolverZheevd_buffersize(handle, jobz, uplo, ndim, d_A, ndim, d_eigval, lwork)
   lwork2= lwork
   istatus= hipMalloc(d_workspace, lwork2)

   call hipCheck(hipsolverZheevd(handle, &
          jobz, uplo, ndim, d_A, ndim, d_eigval, d_workspace, lwork, devinfo))

   if (istatus/=0) print *, 'ERROR: something wrong with hipsolverZheevd on GPU'

   !> step 4: copy data from device back to host
   !> eigval= d_eigval
   !> Amat= d_A
   istatus= hipMemcpy(Amat, d_A, hipMemcpyDeviceToHost)
   istatus= hipMemcpy(eigval, d_eigval, hipMemcpyDeviceToHost)

   !> step 5: free memory on GPU
   !deallocate(d_A, d_eigval, d_workspace)
   istatus= hipsolverDestroy(handle)
   istatus= hipFree(d_A)
   istatus= hipFree(d_eigval)
   istatus= hipFree(d_workspace)

   return
end subroutine hipsolver_zheev

!> A subroutine performs y=A*x
!> A is  CSR sparse format stored, and the dimension is (ndim,ndim)
subroutine hipsparse_zcsrmv(ndim, nnz, x, iA, jA, matA, y)
    use hipfort
    use hipfort_hipsparse
    use iso_c_binding

    implicit none

    integer, parameter :: dp=kind(1d0)

    !> leading dimension of matrix matA
    integer, intent(in) :: ndim

    !> number of non-zeros entries of matrix matA
    integer, intent(in) :: nnz

    !> row index, column index, and values of matA
    !> ia and ja are one based.
    integer, target :: ia(ndim+1)
    integer, target :: ja(nnz)
    complex(dp), target :: matA(nnz)

    !> a dense vector
    complex(dp), target :: x(ndim)
    complex(dp), target :: y(ndim)

    integer :: istatus
    complex(dp), target :: alpha=1.0, beta=0.0

    integer*8, target :: buffersize

    !> define device variables
    type(c_ptr) :: buffer
    type(c_ptr) :: d_rowIdx, d_colIdx, d_matA

    type(c_ptr) :: handle
    type(c_ptr) :: matrixdescr
    type(c_ptr) :: VecXdescr
    type(c_ptr) :: VecYdescr

    type(c_ptr) :: d_x, d_y

    integer(c_int64_t) :: rows_li
    integer(c_int64_t) :: cols_li
    integer(c_int64_t) :: nnz_li
    
    rows_li= ndim
    cols_li= ndim
    nnz_li= nnz

    !> allocate memory on device, this works if compiled with HIP Fortran
    !> compiler such as nvfortran
    !allocate(d_x(ndim), d_y(ndim), d_csrRowPtr(ndim+1))
    !allocate(d_rowIdx(nnz), d_colIdx(nnz), d_matA(nnz))
    istatus= hipMalloc(d_x, sizeof(x))
    istatus= hipMalloc(d_y, sizeof(y))
    istatus= hipMalloc(d_rowIdx, sizeof(ia))
    istatus= hipMalloc(d_colIdx, sizeof(ja))
    istatus= hipMalloc(d_matA, sizeof(matA))

    y=0d0
    !> copy the data from host (CPU) to device (GPU)
    !> d_x= x; d_rowIdx= ia; d_colIdx= ja; d_matA= matA; d_y= y
    istatus= hipMemcpy(d_x, c_loc(x), sizeof(x), hipMemcpyHostToDevice)
    istatus= hipMemcpy(d_rowIdx, c_loc(ia), sizeof(ia), hipMemcpyHostToDevice)
    istatus= hipMemcpy(d_colIdx, c_loc(ja), sizeof(ja), hipMemcpyHostToDevice)
    istatus= hipMemcpy(d_matA, c_loc(matA), sizeof(matA), hipMemcpyHostToDevice)

    istatus=hipsparseCreate(handle)

    !> HIP_C_64F  =  5, /* complex as a pair of double numbers */
    !> HIPSPARSE_INDEX_32I = 2, ///< 32-bit signed integer for matrix/vector
    !> indices>
    istatus=hipsparseCreateCsr(matrixdescr, rows_li, cols_li, nnz_li, d_rowIdx, d_colIdx, d_matA,&
       HIPSPARSE_INDEX_32I, HIPSPARSE_INDEX_32I, HIPSPARSE_INDEX_BASE_ONE, HIP_C_64F)

    istatus=hipsparseCreateDnVec(vecXdescr, rows_li, d_x, HIP_C_64F)
    istatus=hipsparseCreateDnVec(vecYdescr, rows_li, d_y, HIP_C_64F)

    buffersize= 0
    istatus=hipsparseSpMV_bufferSize(handle, HIPSPARSE_OPERATION_NON_TRANSPOSE, &
       c_loc(alpha), matrixdescr, VecXdescr, c_loc(beta), vecYdescr, HIP_C_64F, HIPSPARSE_SPMV_ALG_DEFAULT, c_loc(buffersize))
    istatus=hipMalloc(buffer,buffersize)

    !> perform y=A*x on device
    istatus=hipsparseSpMV(handle, HIPSPARSE_OPERATION_NON_TRANSPOSE, c_loc(alpha), &
       matrixdescr, VecXdescr, c_loc(beta), vecYdescr, HIP_C_64F, HIPSPARSE_SPMV_ALG_DEFAULT, buffer)

    !> copy the data from device back to host
    !> y=d_y
    istatus= hipMemcpy(c_loc(y), d_y, sizeof(y), hipMemcpyDeviceToHost)

    istatus=hipsparseDestroy(handle)
    istatus=hipsparseDestroyDnVec(VecXdescr)
    istatus=hipsparseDestroyDnVec(VecYdescr)
    istatus=hipsparseDestroySpMat(matrixdescr)

    istatus=hipFree(buffer)
    istatus=hipFree(d_colIdx)
    istatus=hipFree(d_matA)
    istatus=hipFree(d_rowIdx)
    istatus=hipFree(d_x)
    istatus=hipFree(d_y)

    return
end subroutine hipsparse_zcsrmv



!> A subroutine performs y=A*x
!> A is  COO sparse format stored, and the dimension is (ndim,ndim)
subroutine hipsparse_zcoomv(ndim, nnz, x, iA, jA, matA, y)
    use hipfort
    use hipfort_hipsparse
    use iso_c_binding

    implicit none

    integer, parameter :: dp=kind(1d0)

    !> leading dimension of matrix matA
    integer, intent(in) :: ndim

    !> number of non-zeros entries of matrix matA
    integer, intent(in) :: nnz

    !> row index, column index, and values of matA
    !> ia and ja are one based.
    integer, target :: ia(nnz)
    integer, target :: ja(nnz)
    complex(dp), target :: matA(nnz)

    !> a dense vector
    complex(dp), target :: x(ndim)
    complex(dp), target :: y(ndim)

    integer :: istatus
    complex(dp), target :: alpha=1.0, beta=0.0

    integer*8, target :: buffersize
    integer, allocatable :: csrRow(:)

    !> define device variables
    type(c_ptr) :: buffer
    type(c_ptr) :: d_rowIdx, d_colIdx, d_csrRowPtr, d_matA

    type(c_ptr) :: handle
    type(c_ptr) :: matrixdescr
    type(c_ptr) :: VecXdescr
    type(c_ptr) :: VecYdescr

    type(c_ptr) :: d_x, d_y

    integer(c_int64_t) :: rows_li
    integer(c_int64_t) :: cols_li
    integer(c_int64_t) :: nnz_li
    
    rows_li= ndim
    cols_li= ndim
    nnz_li= nnz

    allocate(csrRow(ndim+1))
    csrRow= 1
 
    !> allocate memory on device, this works if compiled with HIP Fortran
    !> compiler such as nvfortran
    !allocate(d_x(ndim), d_y(ndim), d_csrRowPtr(ndim+1))
    !allocate(d_rowIdx(nnz), d_colIdx(nnz), d_matA(nnz))
    istatus= hipMalloc(d_x, sizeof(x))
    istatus= hipMalloc(d_y, sizeof(y))
    istatus= hipMalloc(d_csrRowPtr, sizeof(csrRow))
    istatus= hipMalloc(d_colIdx, sizeof(ia))
    istatus= hipMalloc(d_rowIdx, sizeof(ja))
    istatus= hipMalloc(d_matA, sizeof(matA))

    y=0d0
    !> copy the data from host (CPU) to device (GPU)
    !> d_x= x; d_rowIdx= ia; d_colIdx= ja; d_matA= matA; d_y= y
    istatus= hipMemcpy(d_x, c_loc(x), sizeof(x), hipMemcpyHostToDevice)
    istatus= hipMemcpy(d_rowIdx, c_loc(ia), sizeof(ia), hipMemcpyHostToDevice)
    istatus= hipMemcpy(d_colIdx, c_loc(ja), sizeof(ia), hipMemcpyHostToDevice)
    istatus= hipMemcpy(d_matA, c_loc(matA), sizeof(matA), hipMemcpyHostToDevice)

    istatus=hipsparseCreate(handle)

    !> convert coo to csr format on device
    istatus=hipsparseXcoo2csr(handle, d_rowIdx, nnz, ndim, d_csrRowPtr, HIPSPARSE_INDEX_BASE_ONE)

    !> HIP_C_64F  =  5, /* complex as a pair of double numbers */
    !> HIPSPARSE_INDEX_32I = 2, ///< 32-bit signed integer for matrix/vector
    !> indices>
    istatus=hipsparseCreateCsr(matrixdescr, rows_li, cols_li, nnz_li, d_csrRowPtr, d_colIdx, d_matA,&
       HIPSPARSE_INDEX_32I, HIPSPARSE_INDEX_32I, HIPSPARSE_INDEX_BASE_ONE, HIP_C_64F)

    istatus=hipsparseCreateDnVec(VecXdescr, rows_li, d_x, HIP_C_64F)
    istatus=hipsparseCreateDnVec(vecYdescr, rows_li, d_y, HIP_C_64F)

    buffersize= 0
    istatus=hipsparseSpMV_bufferSize(handle, HIPSPARSE_OPERATION_NON_TRANSPOSE, &
       c_loc(alpha), matrixdescr, VecXdescr, c_loc(beta), vecYdescr, HIP_C_64F, HIPSPARSE_SPMV_ALG_DEFAULT, c_loc(buffersize))
    istatus=hipMalloc(buffer,buffersize)

    !> perform y=A*x on device
    istatus=hipsparseSpMV(handle, HIPSPARSE_OPERATION_NON_TRANSPOSE, c_loc(alpha), &
       matrixdescr, VecXdescr, c_loc(beta), vecYdescr, HIP_C_64F, HIPSPARSE_SPMV_ALG_DEFAULT, buffer)

    !> copy the data from device back to host
    !> y=d_y
    istatus= hipMemcpy(c_loc(y), d_y, sizeof(y), hipMemcpyDeviceToHost)

    istatus=hipsparseDestroy(handle)
    istatus=hipsparseDestroyDnVec(VecXdescr)
    istatus=hipsparseDestroyDnVec(VecYdescr)
    istatus=hipsparseDestroySpMat(matrixdescr)

    istatus=hipFree(buffer)
    istatus=hipFree(d_colIdx)
    istatus=hipFree(d_csrRowPtr)
    istatus=hipFree(d_matA)
    istatus=hipFree(d_rowIdx)
    istatus=hipFree(d_x)
    istatus=hipFree(d_y)

    return
end subroutine hipsparse_zcoomv

!> subroutine to perform the inverse of a matrix Amat(ndim, ndim)
!> input/output : Amat(ndim, ndim), Amat will be overwritten on output
!> use hipsolver on GPU device
subroutine hipsolver_zgesv(ndim, Amat)

   use hipfort
   use hipfort_hipblas
   use hipfort_hipsolver
   implicit none

   integer, parameter :: dp=kind(1d0)
   integer, intent(in) :: ndim
   complex(dp), intent(inout), target :: Amat(ndim, ndim)
   complex(dp) :: A_t(ndim, ndim)

   integer :: trans, istatus, lwork, i, info
   integer, allocatable, target :: ipiv(:)
   complex(dp), allocatable, target :: Bmat(:, :)

   !> variables related to GPU
   type(c_ptr) :: handle
   integer(c_int), pointer :: devinfo
   type(c_ptr) :: d_devipiv
   type(c_ptr) :: d_A
   type(c_ptr) :: d_B
   type(c_ptr) :: d_workspace
   integer(c_size_t) :: lwork2

   allocate(ipiv(ndim))
   ipiv= 0

   !> step 1: create hipsolver handle 
   istatus= hipsolverCreate(handle)
  
   !> define the Bmat to be a unitary matrix as an input of hipsolverZgetrs
   allocate(Bmat(ndim, ndim))
   Bmat= 0d0
   do i=1, ndim
      Bmat(i, i)= 1d0
   enddo

   !> step 2: allocate memory on device and copy data from host to device
   !allocate(d_A(ndim, ndim), d_devipiv(ndim), d_B(ndim, ndim))
   istatus= hipMalloc(devinfo)
   istatus= hipMalloc(d_A, sizeof(Amat))
   istatus= hipMalloc(d_B, sizeof(Amat))
   istatus= hipMalloc(d_devipiv, sizeof(ipiv))

   !> d_A= Amat; d_devipiv= ipiv; d_B= Bmat
   istatus= hipMemcpy(d_A, c_loc(Amat), sizeof(Amat), hipMemcpyHostToDevice)
   istatus= hipMemcpy(d_B, c_loc(Bmat), sizeof(Bmat), hipMemcpyHostToDevice)
   istatus= hipMemcpy(d_devipiv, c_loc(ipiv), sizeof(ipiv), hipMemcpyHostToDevice)

   !> step 3: query working space of getrf 
   istatus= hipsolverZgetrf_bufferSize(handle, ndim, ndim, d_A, ndim, lwork)
   lwork2= lwork
   istatus= hipMalloc(d_workspace, lwork2)

   !> hipsolverZgetrf computes the LU factorization of a general mxn matrix
   istatus= hipsolverZgetrf(handle, &
         ndim, ndim, d_A, ndim, d_workspace, lwork, d_devipiv, devinfo)

   !> hipsolverZgetrs solves the system of linear equations A*x= b resulting from the LU
   ! factorization of a matrix using hipsolverZgetrf
   istatus= hipsolverZgetrs(handle, &
          HIPBLAS_OP_N, ndim, ndim, d_A, ndim, d_devipiv, d_B, ndim, &
          d_workspace, lwork, devinfo)

   if (istatus/=0) print *, 'ERROR: something wrong with zgesv on GPU'

   !> step 4: copy data from device back to host
   !> Amat= d_B
   istatus= hipMemcpy(c_loc(Amat), d_B, sizeof(Amat), hipMemcpyDeviceToHost)

   !> free memory on GPU
   istatus=hipsolverDestroy(handle)
   istatus= hipFree(d_A)
   istatus= hipFree(d_B)
   istatus= hipFree(d_devipiv)
   istatus= hipFree(d_workspace)
   deallocate(Bmat, ipiv)

   return
end subroutine hipsolver_zgesv

!> subroutine to perform the inverse of a matrix Amat(ndim, ndim)
!> input/output : Amat(ndim, ndim), Amat will be overwritten on output
!> use magmaf on GPU device
!> A*X= B
subroutine magmaf_zheevd_gpu_wt(JOBZ,UPLO,Ndim,A,W)

   use magma
   implicit none
   integer, parameter :: dp=kind(1d0)

   character*1, intent(in) :: JOBZ
   character*1, intent(in) :: UPLO
   integer,   intent(in) :: Ndim
   real(dp), intent(inout) :: W(Ndim)
   complex(dp),intent(inout) :: A(Ndim,Ndim)

   integer :: info, lwork, liwork, lrwork

   integer, allocatable :: iwork(:)

   real(dp),allocatable ::  rwork(:)

   complex(dp), allocatable :: work(:)

   !> variables related to GPU
   magma_devptr_t :: d_A, handle

   !> allocate GPU memory
   info = magmaf_zmalloc( d_A, ndim*ndim )
   
   !> copy a matrix from host to gpu device 
   call magmaf_queue_create( 0, handle)
   call magmaf_zsetmatrix( ndim, ndim, A, ndim, d_A, ndim, handle)
   call magmaf_queue_destroy( handle )

   if (Ndim==1) then
      W=A(1, 1)
      A(1, 1)= 1d0
      return
   endif

   liwork = 5*Ndim + 3
   lrwork = 2*Ndim*Ndim + 5*Ndim + 1
   lwork = Ndim*(Ndim+2)
   info = 0

   !> d_A= Amat
   allocate (work(lwork), rwork(lrwork), iwork(liwork))
   !> only A is on the gpu device
   call magmaf_zheevd_gpu( jobz, uplo, ndim, d_A, ndim, W, A, ndim, work, lwork, rwork,  &
        lrwork, iwork, liwork, info )

   !> step : copy data from device back to host
   call magmaf_queue_create( 0, handle)
   call magmaf_zgetmatrix(ndim, ndim, d_A, ndim, A, ndim, handle)
   call magmaf_queue_destroy( handle )

   if (info/=0) print *, 'ERROR: something wrong with magmaf_zheevd_gpu on GPU'

   !> free memory on GPU
   info = magmaf_free( d_A )
   deallocate(work, rwork, iwork)

   return
end subroutine magmaf_zheevd_gpu_wt



!> subroutine to perform the inverse of a matrix Amat(ndim, ndim)
!> input/output : Amat(ndim, ndim), Amat will be overwritten on output
!> use magmaf on GPU device
!> A*X= B
subroutine magmaf_zgesv_gpu_wt(ndim, Amat)

   use magma
   implicit none

   integer, parameter :: dp=kind(1d0)
   integer, intent(in) :: ndim
   complex(dp), intent(inout), target :: Amat(ndim, ndim)

   integer :: i, info
   integer, allocatable, target :: ipiv(:)
   complex(dp), allocatable, target :: Bmat(:, :)

   !> variables related to GPU
   magma_devptr_t :: d_A, d_B, handle

   allocate(ipiv(ndim))
   ipiv = 0
   info = 0

   !> allocate GPU memory
   info = magmaf_zmalloc( d_A, ndim*ndim )
   info = magmaf_zmalloc( d_B, ndim*ndim )
 
   !> copy a matrix from host to gpu device 
   call magmaf_queue_create( 0, handle)
   call magmaf_zsetmatrix( ndim, ndim, Amat, ndim, d_A, ndim, handle)
   call magmaf_queue_destroy( handle )

   !> define the Bmat to be a unitary matrix as an input of hipsolverZgetrs
   allocate(Bmat(ndim, ndim))
   Bmat= 0d0
   do i=1, ndim
      Bmat(i, i)= 1d0
   enddo

   !> B mat
   call magmaf_queue_create(0, handle)
   call magmaf_zsetmatrix( ndim, ndim, Bmat, ndim, d_B, ndim, handle)
   call magmaf_queue_destroy( handle )

   !> d_A= Amat; d_devipiv= ipiv; d_B= Bmat
   call magmaf_zgesv_gpu( ndim, ndim, d_A, ndim, ipiv, d_B, ndim, info )

   if (info/=0) print *, 'ERROR: something wrong with zgesv on GPU'

   !> step 4: copy data from device back to host
   call magmaf_queue_create( 0, handle)
   call magmaf_zgetmatrix(ndim, ndim, d_B, ndim, Amat, ndim, handle)
   call magmaf_queue_destroy( handle )

   !> free memory on GPU
   info = magmaf_free( d_A )
   info = magmaf_free( d_B )
   deallocate(Bmat, ipiv)

   return
end subroutine magmaf_zgesv_gpu_wt


! ! compile it with
! !> nvfortran -dcu -cudalib=hipsparse,hipsolver,hipblas -Mpreprocess -C -traceback cuda_math.cuf
! ! run it with 
! !  ./a.out
! !> modified by QSWu on Nov 17 2022
! !> tested with HIP Toolkit v11.8  HIP 22.5
! program product
!     use hipfort
!     use hipfort_hipblas
!     use hipfort_hipsparse
!     use iso_c_binding
! 
!     implicit none
! 
!     complex(8),allocatable :: x(:), matA(:), y(:), x2(:)
!     integer*4,allocatable :: rowIdx(:), colIdx(:),csrRowPtr(:)
! 
!     integer*4 :: nx, nz, status, i
!     complex(8), allocatable :: matA_dense(:, :)
!     complex(8), allocatable :: matB_dense(:, :)
!     complex(8), allocatable :: matC_dense(:, :)

!     real(8), allocatable :: eigval(:)
! 
!     integer*4 :: start, end, countrate
!     complex(8) :: t

!   interface hipsparse_zcsrmv
!     subroutine hipsparse_zcsrmv(ndim, nnz, x, iA, jA, matA, y)
!   use hipfort
!   use hipfort_hipblas
!   use hipfort_hipsparse
!   use iso_c_binding

!   implicit none

!   integer, parameter :: dp=kind(1d0)

!   !> leading dimension of matrix matA
!   integer, intent(in) :: ndim

!   !> number of non-zeros entries of matrix matA
!   integer, intent(in) :: nnz

!   !> row index, column index, and values of matA
!   !> ia and ja are one based.
!   integer, target :: ia(nnz)
!   integer, target :: ja(nnz)
!   complex(dp), target :: matA(nnz)

!   !> a dense vector
!   complex(dp), target :: x(ndim)
!   complex(dp), target :: y(ndim)
!     end subroutine

!     end interface


!     interface hipsparse_zcoomv
!     subroutine hipsparse_zcoomv(ndim, nnz, x, iA, jA, matA, y)
!   use hipfort
!   use hipfort_hipsparse
!   use iso_c_binding

!   implicit none

!   integer, parameter :: dp=kind(1d0)

!   !> leading dimension of matrix matA
!   integer, intent(in) :: ndim

!   !> number of non-zeros entries of matrix matA
!   integer, intent(in) :: nnz

!   !> row index, column index, and values of matA
!   !> ia and ja are one based.
!   integer, target :: ia(nnz)
!   integer, target :: ja(nnz)
!   complex(dp), target :: matA(nnz)

!   !> a dense vector
!   complex(dp), target :: x(ndim)
!   complex(dp), target :: y(ndim)
!     end subroutine

!     end interface

!    interface hipsolver_zgesv2
!    subroutine hipsolver_zgesv2(ndim, Amat)

!  use hipfort
!  use hipfort_hipsolver
!  implicit none
!  
!  integer, parameter :: dp=kind(1d0)
!  integer, intent(in) :: ndim
!  complex(dp), intent(inout), target :: Amat(ndim, ndim)
!
!    end subroutine hipsolver_zgesv2
!    end interface hipsolver_zgesv2
! 
!    interface hipsolver_zgesv
!    subroutine hipsolver_zgesv(ndim, Amat)

!  use hipfort
!  use hipfort_hipsolver
!  implicit none
!  
!  integer, parameter :: dp=kind(1d0)
!  integer, intent(in) :: ndim
!  complex(dp), intent(inout), target :: Amat(ndim, ndim)
!
!    end subroutine hipsolver_zgesv
!    end interface hipsolver_zgesv
! 
!     nx=4
!     nz=9
! 
!     allocate(x(nx), y(nx), csrRowPtr(nx+1), x2(nx))
!     allocate(rowIdx(nz), colIdx(nz), matA(nz))
!     allocate(matA_dense(nx, nx))
!     allocate(matB_dense(nx, nx))
!     allocate(matC_dense(nx, nx))
!     allocate(eigval(nx))
  
!     !> a standard matrix used in the example of hipsparse library
!     !> https://github.com/NVIDIA/HIPLibrarySamples/tree/master/cuSPARSE/spmv_csr
!     rowIdx=(/0, 0, 0, 1, 2, 2, 2, 3, 3/)
!     csrRowPtr=(/0, 3, 4, 7, 9/)
!     colIdx=(/0, 2, 3, 1, 0, 2, 3, 1, 3/)
!     matA=(/1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0/)
!     x=(/1.0, 2.0, 3.0, 4.0/)
!     !> the results of y=A*x=(/19, 8, 51, 52/)
! 
!     rowIdx=rowIdx+1
!     colIdx=colIdx+1
!     csrRowPtr= csrRowPtr+1
! 
!     call hipsparse_zcoomv(nx, nz, x, rowIdx, colIdx, matA, y)
!     print *, 'x:', x
!     print *, 'y:', y
! 
!     call hipsparse_zcsrmv(nx, nz, x, csrRowPtr, colIdx, matA, y)
!     print *, 'x:', x
!     print *, 'y:', y
! 
!     !>> Matrix multiplication of A*A
!     matA_dense(1, :)=(/0, 1, 0, 0/)
!     matA_dense(2, :)=(/1, 0, 1, 0/)
!     matA_dense(3, :)=(/0, 1, 0, 1/)
!     matA_dense(4, :)=(/0, 0, 1, 0/)
! 
!     !> gest for mat_mul_dcu_z
!     matB_dense= matA_dense
!     call mat_mul_dcu_z(nx, matA_dense, matB_dense, matC_dense)
!     print *, 'matC_dense=A*A'
!     do i=1, nx
!        write(*, '(100f8.2)')real(matC_dense(i, :))
!     enddo

!     !>> Inverse of a matrix 
!     matA_dense(1, :)=(/0, 1, 0, 0/)
!     matA_dense(2, :)=(/1, 0, 1, 0/)
!     matA_dense(3, :)=(/0, 1, 0, 1/)
!     matA_dense(4, :)=(/0, 0, 1, 0/)
! 
!     print *, 'matA_dense'
!     do i=1, nx
!        write(*, '(100f8.2)')real(matA_dense(i, :))
!     enddo
!     call hipsolver_zgesv (nx, matA_dense)
!     print *, 'inverse of matA_dense'
!     do i=1, nx
!        write(*, '(100f8.2)')real(matA_dense(i, :))
!     enddo
! 

!     !> get eigenvalue of A
!     matA_dense(1, :)=(/0, 1, 0, 0/)
!     matA_dense(2, :)=(/1, 0, 1, 0/)
!     matA_dense(3, :)=(/0, 1, 0, 1/)
!     matA_dense(4, :)=(/0, 0, 1, 0/)
!
!     call hipsolver_zheev(nx, matA_dense, eigval)
!     print *, 'eigval:', eigval
!     print *, ' '
!     print *, '>> Congratulations : test passed'
! 
! end
# endif 
