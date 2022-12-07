!> To generate cuda_blas.o
!> nvfortran -cuda -cudalib=cusparse,cusolver,cublas -Mpreprocess -C -traceback cuda_math.cuf
!> written by QSWu on Nov 18 2022 at IOP CAS, Beijing
!> Email: quansheng.wu@iphy.ac.cn or wuquansheng@gmail.com
!> tested with CUDA Toolkit v11.8  CUDA 22.5

!> A subroutine performs y=A*x
!> A is  COO sparse format stored, and the dimension is (ndim,ndim)
module cucsrmv_module
    use cudafor
    use cusparse
    use iso_c_binding

    implicit none

    integer, parameter :: dp_8=kind(1d0)
    integer :: istatus
    !> define device variables
    integer, allocatable, device :: buffer(:)
    integer, allocatable, device :: d_rowIdx(:), d_colIdx(:)
    complex(dp_8), allocatable, device :: d_x(:), d_matA(:), d_y(:)

    type(cusparseHandle) :: handle
    type(cusparseSpMatDescr) :: matrix
    type(cusparseDnVecDescr) :: vecX, vecY

contains
    subroutine cusparse_zcsrmv_allocate(ndim, nnz, A, ia, ja)
       implicit none
       integer, intent(in) :: ndim
       integer, intent(in) :: nnz
       integer, intent(in) :: ia(ndim+1)
       integer, intent(in) :: ja(nnz)
       complex(dp_8), intent(in) :: A(nnz)
       istatus= cusparseCreate(handle)
       allocate(d_x(ndim), d_y(ndim))
       allocate(d_rowIdx(ndim+1), d_colIdx(nnz), d_matA(nnz))
       d_matA= A; d_rowIdx= ia; d_colIdx= ja
       return
    end subroutine cusparse_zcsrmv_allocate
 
    subroutine cusparse_zcsrmv_deallocate
       implicit none
       istatus= cusparseDestroy(handle)
       istatus= cusparseDestroyDnVec(vecX)
       istatus= cusparseDestroyDnVec(vecY)
       istatus= cusparseDestroySpMat(matrix)
   
       deallocate(d_x, d_y, d_matA, d_rowIdx, d_colIdx)
       return
    end subroutine cusparse_zcsrmv_deallocate
end module cucsrmv_module
