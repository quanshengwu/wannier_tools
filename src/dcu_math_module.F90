#if defined (DCU)
module dcu_matmul
   use iso_c_binding
   use para, only : dp
   use hipfort
   use hipfort_check
   use hipfort_hipblas
   implicit none

   !> variables on GPU device
   type(c_ptr) :: handle_mm= c_null_ptr
   type(c_ptr) :: d_A_mm
   type(c_ptr) :: d_B_mm
   type(c_ptr) :: d_C_mm

   contains
   subroutine dcu_mm_allocate(ndim)
      integer, intent(in) :: ndim 
      call hipCheck(hipblasCreate(handle_mm))
      call hipCheck(hipMalloc(d_A_mm, int(ndim*16*ndim, c_size_t)))
      call hipCheck(hipMalloc(d_B_mm, int(ndim*16*ndim, c_size_t)))
      call hipCheck(hipMalloc(d_C_mm, int(ndim*16*ndim, c_size_t)))
   end subroutine dcu_mm_allocate

   subroutine dcu_mm_deallocate()
      call hipCheck(hipblasDestroy(handle_mm))
      call hipCheck(hipFree(d_A_mm))
      call hipCheck(hipFree(d_B_mm))
      call hipCheck(hipFree(d_C_mm))
   end subroutine dcu_mm_deallocate
end module dcu_matmul

module dcu_zheev
   use para, only : dp
   use hipfort
   use hipfort_check
   use hipfort_hipblas
   use hipfort_hipsolver
   implicit none

   !> variables related to GPU
   type(c_ptr) :: handle_zheev
   integer(c_int), pointer :: devinfo_zheev
   complex(dp), pointer, dimension(:, :) :: d_A_zheev
   real(dp), pointer, dimension(:) :: d_eigval_zheev

   contains
   subroutine dcu_zheev_allocate(ndim)
      integer, intent(in) :: ndim 

      !> step 1: create hipsolver handle 
      call hipCheck( hipsolverCreate(handle_zheev))
      
      !> step 2: allocate memory on device and copy data from host to device
      !allocate(d_A(ndim, ndim), d_eigval(ndim))
      call hipCheck(hipMalloc(devinfo_zheev))
      call hipCheck(hipMalloc(d_eigval_zheev, ndim))
      call hipCheck(hipMalloc(d_A_zheev, ndim, ndim))
   
   end subroutine dcu_zheev_allocate

   subroutine dcu_zheev_deallocate()
       call hipCheck(hipsolverDestroy(handle_zheev))
       call hipCheck(hipFree(d_A_zheev))
       call hipCheck(hipFree(devinfo_zheev))
       call hipCheck(hipFree(d_eigval_zheev))
  end subroutine dcu_zheev_deallocate
end module dcu_zheev

#endif
