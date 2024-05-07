!>  Define some CSR and ParCSR data structure
!>* modified by QuanSheng Wu on Apr. 28. 2015
!>* unstable version

module sparse

   use prec
   use WMPI
   use para
   implicit none


   !> define CSR(compressed sparse row) format matrix
   type WTCSR
      !> number of local rows, for different cpu, we get different rows
      integer(li) :: NumRows

      !> number of columns, for row-partioning method, this value is equals to
      !> the global columns of matrix, in our physics problem, NumCols=
      !> GlobalNumCols= GlobalNumRows
      integer(li) :: NumCols

      !> number of non-zeros entities
      integer(li) :: NumNonZeros

      !> array of length NumRows, record the row index
      !> because the index for ia is from 1 to NumRows, for this matrix, we need
      !> to know its row index
      integer(li), pointer :: RowIndex(:)

      !> integer array of length NumRows+1 containing index of first nonzero
      !> element of row
      integer(li), pointer :: ia(:)

      !> integer array of length NumNonZeros containing the column positions
      !> of the corresponding elements in a.
      integer(li), pointer :: ja(:)

      !> complex data for non-zeros entity
      complex(dp), pointer :: a(:)
   end type WTCSR

   !> Define parallel version of CSR format matrix ParCSR
   type WTParCSR

      !> mpi communicator
      integer :: comm

      !> global number of rows
      integer(li) :: GlobalNumRows

      !> global number of columns
      integer(li) :: GlobalNumCols

      !> number of non-zeros
      integer(li) :: NumNonZeros

      !> when doing matrix-vector multiplication, the diagonal block needn't
      !> the vector data from other procs, while the offdiagonal block needn't
      !> the local cpu's vector data. So, we can realize non-block communication
      !> when using ParCSR.
      !>* For realization, we
      !> store the diagonal and offdiagonal block as CSR matrix
      type(WTCSR), pointer :: Diag

      !> the column index of offd is not the original global column index,
      !> so we need another array ColMapOffd
      type(WTCSR), pointer :: Offd

      !> map column index of offd to global column index
      integer(li), pointer :: ColMapOffd(:)

      !> array of length NumCPUs, RowStarts(i) contains the
      !> global number of the first row on proc i,
      !> first_row_index = row_starts[my_id]
      !> row_starts[num_procs] = global_num_rows
      integer(li), pointer :: RowStarts(:)

      !> array of length num_procs+1, col_starts[i] contains the
      !> global number of the first column of diag on proc i,
      !> first_col_diag = col_starts[my_id],
      !> col_starts[num_procs] = global_num_cols */
      integer(li), pointer :: ColStarts(:)

      !> an flag that whether the matrix has generate sendrecv
      logical :: IsComm
      !> an object to store the send and recieve messages
      type(WTParCSRComm), pointer :: SendRecv

   end type WTParCSR

   !> define a type for sequential vector
   type WTSeqVec

      !> length of vector
      integer(li) :: length
      complex(dp), pointer :: data(:)

   end type WTSeqVec

   !> define a type for parallel vector
   type WTParVec

      !> mpi communicator
      integer :: comm

      !> length of vector
      integer(li) :: GlobalLength
      integer(li) :: FirstIndex
      integer(li) :: LastIndex
      type(WTSeqVec), pointer :: SeqVec

   end type WTParVec

   !> define some interfaces
   !> redefine some operators
   interface assignment (=)
      !* copy one parvec to another parvec
      module procedure WTParVecCopy
      module procedure WTParVecSetConstantValue_z
      module procedure WTParVecSetConstantValue_r
   end interface

   interface operator (*)
      !* perform two parvec multiply
      module procedure WTParVecInnerProd
      module procedure WTParVecScale_r1
      module procedure WTParVecScale_r2
      module procedure WTParVecScale_z1
      module procedure WTParVecScale_z2
   end interface

   interface operator (.DOT.)
      !* perform two parvec multiply
      module procedure WTParVecInnerProdWithNoReduce
   end interface

   interface WTParCSRMatVec
      module procedure WTParCSRMatVec_c
   end interface

   interface WTParVecSetConstantValue
      module procedure WTParVecSetConstantValue_z
      module procedure WTParVecSetConstantValue_r
   end interface

   interface WTParVecSetSingleValue
      module procedure WTParVecSetSingleValue_r
      module procedure WTParVecSetSingleValue_z
   end interface


   !> Maximum number of non-zero elements in the Hamiltonian matrix
   integer(li), public :: MaxNumNonZeros

   private
   public :: operator(*), assignment(=), operator(.dot.)
   public :: WTCSR
   public :: WTParCSR
   public :: WTParVec
   public :: ConvertCooToCsr
   public :: WTCSRCreate
   public :: WTCSRDestroy
   public :: WTCSRInitialize
   public :: WTSeqVecCreate
   public :: WTSeqVecDestroy
   public :: WTSeqVecInitialize
   public :: WTSeqVecPrint
   public :: WTSeqVecInnerProd

   public :: WTCSRToParCSR
   public :: WTParCSRCreate
   public :: WTParCSRInitialize
   public :: WTParCSRDestroy
   public :: WTParCSRMatVec

   public :: WTGenSendRecv
   public :: WTSendRecvDestroy

   public :: WTParVecCreate
   public :: WTParVecDestroy
   public :: WTParVecInitialize
   public :: WTParVecPrint
   public :: WTParVecCopy
   public :: WTParVecNorm
   public :: WTParVecAxpy
   public :: WTParVecAxpBy
   public :: WTParVecInnerProd
   public :: WTParVecInnerProdWithNoReduce
   public :: WTParVecNormalize
   public :: WTParVecGetSingleValue
   public :: WTParVecSetSingleValue
   public :: WTParVecSetRandomValue
   public :: WTParVecSetConstantValue
   public :: WTParVecSetSeqValue
   public :: WTParCSRMatrixCreate
   public :: csr_sort_indices
   public :: csr_sum_duplicates
   public :: arpack_sparse_coo_eigs
   public :: arpack_sparse_coo_eigs_nonorth
   public :: csrmv_z
   public :: coomv_z

contains

   !> Create a WTCSR type matrix, when create, we only give the basic
   !> information of the matrix, we don't allocate memory
   subroutine WTCSRCreate(nrows, ncols, maxnnz, HCSR)
      implicit none
      integer(li), intent(in) :: nrows
      integer(li), intent(in) :: ncols
      integer(li), intent(in) :: maxnnz
      type(WTCSR), intent(inout) :: HCSR

      integer :: ierr

      !if (.not. associated(HCSR)) allocate(HCSR)
      HCSR%NumRows= nrows
      HCSR%NumCols= ncols
      HCSR%NumNonZeros= maxnnz

      HCSR%RowIndex=> Null()
      if (.not.associated(HCSR%RowIndex)) then
         allocate(HCSR%RowIndex(nrows), stat=ierr)
         HCSR%RowIndex= 0
      endif

      !*initialize pointer
      HCSR%ia=> Null()
      HCSR%ja=> Null()
      HCSR%a=> Null()

      return
   end subroutine WTCSRCreate

   !> Free memory for CSR matrix
   subroutine WTCSRDestroy(HCSR)
      implicit none
      type(WTCSR) :: HCSR

      if (associated(HCSR%RowIndex))deallocate(HCSR%RowIndex)
      if (associated(HCSR%ia))deallocate(HCSR%ia)
      if (associated(HCSR%ja))deallocate(HCSR%ja)
      if (associated(HCSR%a))deallocate(HCSR%a)
      nullify(HCSR%RowIndex, HCSR%ia, HCSR%ja, HCSR%a)

      return
   end subroutine WTCSRDestroy

   !> initialize hamiltoinan matrix
   !> allocate memory
   subroutine WTCSRInitialize(HCSR)
      implicit none

      type(WTCSR) :: HCSR
      integer(li) :: nrows
      integer(li) :: maxnnz

      integer :: ierr
      real(dp) :: memory

      nrows= HCSR%NumRows
      maxnnz= HCSR%NumNonZeros

      !>  calculate memory required by ia, ja, H
      memory= ( nrows*4.0      & ! ia
         + maxnnz*4.0     & ! ja
         + maxnnz*16.0    & ! H
         )/1024.0/1024.0    ! MB

      if (stdout>0) then
         write(stdout,100)'Memory required by Hamiltonian :', memory, ' MB'
      endif


      if (.not.associated(HCSR%ia)) then
         allocate(HCSR%ia(nrows+1), stat=ierr)
         HCSR%ia= 0
      endif

      if (.not.associated(HCSR%ja) .and. maxnnz/=0) then
         allocate(HCSR%ja(maxnnz), stat=ierr)
         HCSR%ja= 0
      endif

      if (.not.associated(HCSR%a) .and. maxnnz/=0) then
         allocate(HCSR%a(maxnnz), stat=ierr)
         HCSR%a = cmplx(0d0, 0d0, dp)
      endif


100   format(2x,a,f9.1,a)

      return
   end subroutine WTCSRInitialize

   !> create a vector
   subroutine WTSeqVecCreate(length, vec)
      implicit none

      integer(li), intent(in) :: length
      type(WTSeqVec) :: vec

      vec%length= length

      return
   end subroutine WTSeqVecCreate

   !> destroy a vector
   subroutine WTSeqVecDestroy(vec)
      implicit none
      type(WTSeqVec) :: vec

      if (associated(vec%data)) deallocate(vec%data)
      nullify(vec%data)

      return
   end subroutine WTSeqVecDestroy

   !> initialize a vector
   subroutine WTSeqVecInitialize(vec)
      implicit none

      type(WTSeqVec) :: vec
      integer :: ierr

      vec%data=> Null()
      if (.not. associated(vec%data)) &
         allocate(vec%data(vec%length), stat= ierr)
      vec%data= (0d0, 0d0)

      return
   end subroutine WTSeqVecInitialize

   !> create a parallel vector
   subroutine WTParVecCreate(comm, GlobalLength, ParVec)
      implicit none

      integer, intent(in) :: comm
      integer(li), intent(in) :: GlobalLength
      type(WTParVec) :: ParVec

      integer(li) :: first
      integer(li) :: last
      integer(li) :: LocalLength
      integer :: NCPUs
      integer :: cpu_id
      integer :: ierr

#if defined (MPI)
      ParVec%Comm= comm
      call mpi_comm_size(comm, NCPUS, ierr)
      call mpi_comm_rank(comm, cpu_id, ierr)
#endif


      call WTGenerateLocalPartition(GlobalLength, NCPUs, &
         cpu_id, first, last)

      ParVec%GlobalLength= GlobalLength
      ParVec%FirstIndex= first
      ParVec%LastIndex= last

      LocalLength= last- first+ 1

      ParVec%SeqVec=> Null()
      if (.not.associated(ParVec%SeqVec))allocate(ParVec%SeqVec)
      call WTSeqVecCreate(LocalLength, ParVec%SeqVec)

      return
   end subroutine WTParVecCreate

   !> initialize a parallel vector
   subroutine WTParVecDestroy(ParVec)
      implicit none
      type(WTParVec) :: ParVec

      call WTSeqVecDestroy(ParVec%SeqVec)
      nullify(ParVec%SeqVec)

      return

   end subroutine WTParVecDestroy

   !> initialize a parallel vector
   subroutine WTParVecInitialize(ParVec)
      implicit none
      type(WTParVec) :: ParVec

      call WTSeqVecInitialize(ParVec%SeqVec)
      return
   end subroutine WTParVecInitialize

   !> Create ParCSR matrix
   subroutine WTParCSRCreate(comm, RowStarts, ColStarts, NumNonZerosDiag, &
      NumNonZerosOffd, HParCSR)

      integer, intent(in) :: comm
      integer(li) :: RowStarts(:)
      integer(li) :: ColStarts(:)
      integer(li), intent(in) :: NumNonZerosDiag
      integer(li), intent(in) :: NumNonZerosOffd
      type(WTParCSR) :: HParCSR

      integer(li) :: i
      integer(li) :: localnumrows
      integer(li) :: localnumcols
      integer(li) :: GlobalNumCols

      integer :: NCPUS
      integer :: icpu
      integer :: ierr

      HParCSR%Comm= comm
      HParCSR%Diag=> Null()
      HParCSR%Offd=> Null()
      HParCSR%ColMapOffd=> Null()
      HParCSR%RowStarts=> Null()
      HParCSR%ColStarts=> Null()
      HParCSR%SendRecv=> Null()

      !> get number of cpus and current cpu id
#if defined (MPI)
      call mpi_comm_size(comm, NCPUS, ierr)
      call mpi_comm_rank(comm, icpu, ierr)
#endif

      allocate(HParCSR%RowStarts(NCPUs+1))
      allocate(HParCSR%ColStarts(NCPUs+1))
      HParCSR%RowStarts= RowStarts  !< maybe some problem
      HParCSR%ColStarts= ColStarts

      HParCSR%GlobalNumRows= RowStarts(NCPUs+1)
      HParCSR%GlobalNumCols= ColStarts(NCPUs+1)

      localnumrows= RowStarts(icpu+2)- RowStarts(icpu+1)
      localnumcols= ColStarts(icpu+2)- ColStarts(icpu+1)

      if (.not. associated(HParCSR%diag)) allocate(HParCSR%diag)
      call WTCSRCreate(localnumrows, localnumcols, NumNonZerosDiag, HParCSR%diag)

      !> when create parcsr matrix, we don't know the details of the matrix,
      !> numrows and numcols for offd will be changed, some arrays will be
      !> reallocated if the new length is lager than the old one.
      if (.not. associated(HParCSR%offd)) allocate(HParCSR%offd)
      call WTCSRCreate(localnumrows, localnumcols, NumNonZerosOffd, HParCSR%Offd)

      !> assign data for RowIndex
      do i=RowStarts(icpu+1), RowStarts(icpu+2)-1
         HParCSR%Diag%RowIndex(i-RowStarts(icpu+1)+1)= i
      enddo
      HParCSR%Offd%RowIndex= HParCSR%Diag%RowIndex

      !> SendRecv is a pointer, when create, we allocate memory for it
      if (.not. associated(HParCSR%SendRecv)) allocate(HParCSR%SendRecv)

      return
   end subroutine WTParCSRCreate

   !> initialize ParCSR matrix
   subroutine WTParCSRInitialize(HParCSR)
      implicit none

      type(WTParCSR) :: HParCSR
      call WTCSRInitialize(HParCSR%diag)
      call WTCSRInitialize(HParCSR%offd)

      !> this one will be reallocate if the real columns is greater than the
      !> current one
      if (.not.associated(HParCSR%ColMapOffd).and.HParCSR%Offd%NumCols/=0) then
         allocate(HParCSR%ColMapOffd(HParCSR%Offd%NumCols))
      endif

      !* commpkg
      HParCSR%IsComm= .false.

      return
   end subroutine WTParCSRInitialize

   !>> Free memory for ParCSR matrix
   subroutine WTParCSRDestroy(HParCSR)
      type(WTParCSR)  :: HParCSR

      call WTCSRDestroy(HParCSR%diag)
      call WTCSRDestroy(HParCSR%offd)
      if (associated(HParCSR%ColMapOffd))deallocate(HParCSR%ColMapOffd)
      if (associated(HParCSR%RowStarts))deallocate(HParCSR%RowStarts)
      if (associated(HParCSR%ColStarts))deallocate(HParCSR%ColStarts)
      call WTSendRecvDestroy(HParCSR%sendrecv)
      Nullify(HParCSR%diag, HParCSR%offd)
      Nullify(HParCSR%ColStarts, HParCSR%sendrecv)
      Nullify(HParCSR%ColMapOffd, HParCSR%RowStarts, HParCSR%ColStarts)

      return
   end subroutine WTParCSRDestroy

   !> create a hamiltonian object
   subroutine WTParCSRMatrixCreate(Mdim, HParCSR)

      use para
      implicit none

      integer, intent(in) :: Mdim
      type (WTParCSR), intent(out) :: HParCSR

      integer(li), pointer :: RowStarts(:)
      integer(li), pointer :: ColStarts(:)

      integer :: ierr
      integer :: NCPUS

      NCPUS= 1
#if defined (MPI)
      call mpi_comm_size(mpi_cmw, NCPUS, ierr)
#endif


      RowStarts=> Null()
      ColStarts=> Null()
      allocate(RowStarts(NCPUs+1))
      allocate(ColStarts(NCPUs+1))
      call WTGeneratePartition(Mdim, NCPUS, RowStarts)
      ColStarts= RowStarts


      if (stdout>0) then
         write(stdout, *)' '
         write(stdout, '(2x,a)')'>> RowStarts'
         write(stdout, '(5i12)')RowStarts
         write(stdout, *)' '
      endif

#if defined (MPI)
      call WTParCSRCreate(mpi_cmw, RowStarts, ColStarts, 0, 0, HParCSR)
#endif

      deallocate(RowStarts, ColStarts)

      return
   end subroutine WTParCSRMatrixCreate


   !>> Converts CSR matrix to ParCSR matrix
   subroutine WTCSRToParCSR(HCSR, HParCSR)

      type(WTCSR) :: HCSR
      type(WTParCSR) :: HParCSR

      !> here numrows is the local rows
      integer(li) :: NumRows

      !> here numcols is the global columns and rows
      integer(li) :: NumCols

      !> number of columns for diagonal block
      integer(li) :: NumColsDiag

      !> number of columns for off-diagonal block, this value is unknow, and
      !> should be set by checking how many columns that have nonzero entities
      integer(li) :: NumColsOffd

      !> mark the diagonal block
      integer(li) :: FirstColDiag
      integer(li) :: LastColDiag

      integer(li) :: NumNonZeros

      integer(li), pointer :: a_i(:)
      integer(li), pointer :: a_j(:)
      complex(dp), pointer :: a(:)

      integer(li), pointer :: ColMapOffd(:)

      !> mark 1 if the j'th column has nonzeros value, length=NumCols
      integer, allocatable :: marker(:)

      !> For each convert, we generate the partition once.
      integer(li), allocatable :: RowStarts(:)
      integer(li), allocatable :: ColStarts(:)

      !> some loop index
      integer :: i, j, jo, jd

      integer :: counter
      integer :: NCPUS
      integer :: icpu
      integer :: comm
      integer :: ierr

      !a_i=> Null()
      !a_j=> Null()
      !a  => Null()

#if defined (MPI)
      comm= HParCSR%Comm
      !> get number of cpus and current cpu id
      call mpi_comm_size(comm, NCPUS, ierr)
      call mpi_comm_rank(comm, icpu, ierr)
#endif

      NumRows= HCSR%NumRows
      NumCols= HCSR%NumCols
      NumNonZeros= HCSR%NumNonZeros

      allocate(marker(NumCols))
      allocate(RowStarts(NCPUs+1))
      allocate(ColStarts(NCPUs+1))
      marker= 0
      RowStarts= HParCSR%RowStarts
      ColStarts= HParCSR%ColStarts

      NumColsDiag= ColStarts(icpu+2)- ColStarts(icpu+1)+1
      NumColsOffd= 0
      FirstColDiag= ColStarts(icpu+1)
      LastColDiag = ColStarts(icpu+2)- 1

      a_i=> HCSR%ia
      a_j=> HCSR%ja
      a  => HCSR%a


      if (NumCols.gt.NumColsDiag) then !< there are off-diagonal block

         !* here we don't know the nnz for diag and offd, so we just
         !* initialize ia and RowIndex
         !* PS: when create HParCSR, we should set maxnnz= 0, in order that
         !* we don't allocate memory for ja, a
         call WTCSRInitialize(HParCSR%diag)
         call WTCSRInitialize(HParCSR%offd)

         jo=1
         jd=1
         do i=1, NumRows
            HParCSR%offd%ia(i)= jo
            HParCSR%diag%ia(i)= jd

            !* sweep non-zero entries in each row
            do j=a_i(i), a_i(i+1)-1
               if (a_j(j)<FirstColDiag .or. a_j(j)>LastColDiag) then !< offd
                  if (marker(a_j(j))==0) then
                     marker(a_j(j))= 1
                     NumColsOffd= NumColsOffd+ 1
                  endif
                  jo= jo+ 1 !< number of non-zeros in the non-diagonal block
               else
                  jd= jd+ 1 !< number of non-zeros in the diagonal block
               endif
            enddo
         enddo
         HParCSR%offd%ia(NumRows+1)= jo
         HParCSR%diag%ia(NumRows+1)= jd

         if (NumColsOffd>0) then
            allocate(HParCSR%ColMapOffd(NumColsOffd))
            HParCSR%ColMapOffd= 0
         endif

         !* count how many non-zeros columns and
         !* record its global column index
         counter= 1
         do i=1, NumCols
            if (marker(i)==1) then
               HParCSR%ColMapOffd(counter)= i
               marker(i)= counter !< new column order
               counter= counter+ 1
            endif
         enddo

         !* till now, we know the nnz of diag and offd, so we can allocate
         !* memory for ja, a of diag and offd
         HParCSR%diag%NumNonZeros= jd
         call WTCSRInitialize(HParCSR%diag)

         HParCSR%offd%NumNonZeros= jo
         call WTCSRInitialize(HParCSR%offd)

         !* set number of columns for off-diagonal block of ParCSR
         HParCSR%offd%NumCols= NumColsOffd

         !* get data for diag and offd
         jo=1
         jd=1
         do i=1, NumRows
            do j=a_i(i), a_i(i+1)-1
               if (a_j(j)<FirstColDiag .or. a_j(j)>LastColDiag) then !< offd
                  HParCSR%offd%ja(jo)= marker(a_j(j))
                  HParCSR%offd%a(jo)= a(j)
                  jo= jo+ 1 !< number of non-zeros in the non-diagonal block
               else
                  HParCSR%diag%ja(jd)= a_j(j)- FirstColDiag+ 1
                  HParCSR%diag%a(jd)= a(j)
                  jd= jd+ 1 !< number of non-zeros in the diagonal block
               endif
            enddo
         enddo
      else !< there is no off-diagonal block
         HParCSR%diag%NumNonZeros= NumNonZeros
         HParCSR%offd%NumNonZeros= 0
         HParCSR%offd%NumCols= 0
         call WTCSRInitialize(HParCSR%diag)

         call WTCSRInitialize(HParCSR%offd)
         do i=1, NumNonZeros
            HParCSR%diag%ja(i)= a_j(i)
            HParCSR%diag%a(i)= a(i)
         enddo
         do i=1, NumRows+1
            HParCSR%diag%ia(i)= a_i(i)
            HParCSR%offd%ia(i)= 0
         enddo
      endif

      !> set rowindex for diag and offd
      !diag%RowIndex= HCSR%RowIndex
      !$OMP PARALLEL DO
      do i=1, NumRows
         HParCSR%diag%RowIndex(i)= HCSR%RowIndex(i)
      enddo
      !$OMP END PARALLEL DO

      !offd%RowIndex= HCSR%RowIndex
      !$OMP PARALLEL DO
      do i=1, NumRows
         HParCSR%offd%RowIndex(i)= HCSR%RowIndex(i)
      enddo
      !$OMP END PARALLEL DO

      deallocate(marker)

      return
   end subroutine WTCSRToParCSR


      !>> converts COO to CSR in place.
      subroutine ConvertCooToCsr( n, nnz, a, ia, ja, iwk)

         !*****************************************************************************80
         !   from http://people.sc.fsu.edu/~jburkardt/f_src/sparsekit/sparsekit.f90
         !
         !
         !! COOCSR_INPLACE converts COO to CSR in place.
         !
         !  Discussion:
         !
         !    This routine converts a matrix stored in coordinate format into
         !    the CSR format.  The conversion is done in place in that the arrays
         !    a,ja,ia of the result are overwritten onto the original arrays.
         !
         !    The entries of the output matrix are not sorted (the column
         !    indices in each are not in increasing order) use COOCSR
         !    if you want them sorted.
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
         !    Input, integer N, the row dimension of the matrix.
         !
         !    Input, integer NNZ, the number of nonzero elements in A.
         !
         !    Input, integer JOB.  When JOB = 1, the real values in A are
         !    filled.  Otherwise A is not touched and the structure of the
         !    array only (i.e. JA, IA)  is obtained.
         !
         !    Input/output, real A(NNZ).  On input, the matrix numeric values,
         !    stored in the COO format.  On output, the numeric values, stored
         !    in CSR format.
         !
         ! ja      = integer array of length nnz containing the column positions
         !         of the corresponding elements in a.
         !
         ! ia      = integer array of length nnz containing the row positions
         !         of the corresponding elements in a.
         !
         ! iwk      = integer work array of length n.
         !
         ! on return:
         !
         !    Output, real A(*), integer JA(*), IA(NROW+1), the matrix in CSR
         !    Compressed Sparse Row format.
         !
         implicit none

         ! >> inout variables
         integer(li), intent(in) :: n
         integer(li), intent(in) :: nnz
         integer(li), intent(inout) :: ia(nnz)
         integer(li), intent(inout) :: ja(nnz)
         integer(li), intent(inout) :: iwk(n+1)
         complex(dp), intent(inout) :: a(nnz)

         ! >> local variables
         integer(li) :: i
         integer(li) :: inext
         integer(li) :: init
         integer(li) :: ipos
         integer(li) :: j
         integer(li) :: jnext
         integer :: job
         integer(li) :: k
         complex(dp) ::  t
         complex(dp) ::  tnext
         logical :: values

         if (nnz.eq.0)return

         job=1
         values = (job == 1)
         !
         !  Find offdter array for resulting matrix.
         !
         iwk(1:n+1) = 0

         do k = 1, nnz
            i = ia(k)
            iwk(i+1) = iwk(i+1) + 1
         end do

         iwk(1) = 1
         do i = 2, n
            iwk(i) = iwk(i-1) + iwk(i)
         end do
         !
         !  Loop for a cycle in chasing process.
         !
         init = 1
         k = 0

5        continue

         if ( values ) then
            t = a(init)
         end if

         i = ia(init)
         j = ja(init)
         ia(init) = -1

6        continue
         k = k + 1
         !
         !  Current row number is I.  Determine where to go.
         !
         ipos = iwk(i)
         !
         !  Save the chased element.
         !
         if ( values ) then
            tnext = a(ipos)
         end if

         inext = ia(ipos)
         jnext = ja(ipos)
         !
         !  Then occupy its location.
         !
         if ( values ) then
            a(ipos) = t
!~             write(*,*) t
         end if

         ja(ipos) = j
         !
         !  Update pointer information for next element to come in row I.
         !
         iwk(i) = ipos + 1
         !
         !  Determine the next element to be chased.
         !
         if ( ia(ipos) < 0 ) then
            !if ( ia(ipos) <= 0 ) then !changed by QS.Wu
            go to 65
         end if

         t = tnext
         i = inext
         j = jnext
         ia(ipos) = -1

         if ( k < nnz ) then
            go to 6
         end if

         go to 70

65       continue

         init = init + 1

         if ( nnz < init ) then
            go to 70
         end if

         !if ( ia(init) <= 0 ) then !changed by QS.Wu
         if ( ia(init) < 0 ) then
            go to 65
         end if
         !
         !  Restart chasing.
         !
         go to 5

70       continue

         ia(1) = 1
         ia(2:n+1) = iwk(1:n)

         return
      end subroutine ConvertCooToCsr

      !> generate commpkg for parcsr matrix A
      !> added on May-12-2015 by QS.Wu
      subroutine WTGenSendRecv(HParCSR)
         implicit none

         !> hamiltonian in parcsr matrix format
         type(WTParCSR), intent(inout) :: HParCSR


         !* local variables
         integer(li) :: FirstColDiag
         integer :: NumColsOffd
         integer :: NumColsDiag

         integer, pointer :: SendCPUs(:)
         integer, pointer :: SendMapStarts(:)
         integer(li), pointer :: SendMapElements(:)
         integer, pointer :: RecvCPUs(:)
         integer, pointer :: RecvVecStarts(:)
         integer(li), pointer :: ColMapOffd(:)
         integer(li), pointer :: ColStarts(:)

         integer, pointer :: CPUMark(:)
         integer, pointer :: CPUAdd(:)
         integer, pointer :: tmp(:)
         integer, pointer :: RecvBuf(:)
         integer, pointer :: Displs(:)
         integer, pointer :: info(:)
         integer, pointer :: mpirequest(:)
         integer, pointer :: mpistatus(:)

         integer :: i, j
         integer :: ierr
         integer :: icpu
         integer :: CPUNo
         integer :: NCPUs
         integer :: cpu_id
         integer :: comm
         integer :: VecStart
         integer :: VecLen
         integer :: NumRecvs
         integer :: NumSends
         integer :: NumElements
         integer :: NumRequest
         integer :: localinfo
         integer(li) :: OffdCol

         !* communicator
         type(WTParCSRComm), pointer :: SendRecv


         !* initialize null pointers
         SendCPUs=> Null()
         SendMapStarts=> Null()
         SendMapElements=> Null()
         RecvCPUs=> Null()
         RecvVecStarts=> Null()
         ColMapOffd=> Null()
         ColStarts=> Null()
         CPUMark=> Null()
         CPUAdd=> Null()
         tmp=> Null()
         RecvBuf=> Null()
         Displs=> Null()
         info=> Null()
         mpirequest=> Null()
         mpistatus=> Null()
         SendRecv=> Null()


#if defined (MPI)
         comm= HParCSR%Comm
         call mpi_comm_size(comm, NCPUS, ierr)
         call mpi_comm_rank(comm, cpu_id, ierr)
#endif

         SendRecv=> HParCSR%SendRecv
         SendRecv%Comm= comm
         SendRecv%NumSends= 0
         SendRecv%NumRecvs= 0
         SendRecv%SendCPUs=> Null()
         SendRecv%SendMapStarts=> Null()
         SendRecv%SendMapElements=> Null()
         SendRecv%RecvCPUs=> Null()
         SendRecv%RecvVecStarts=> Null()

         NumColsDiag= HParCSR%Diag%NumCols
         NumColsOffd= HParCSR%Offd%NumCols
         ColMapOffd=> HParCSR%ColMapOffd
         ColStarts=> HParCSR%ColStarts
         FirstColDiag= ColStarts(cpu_id+1)

         allocate(CPUMark(NCPUs), CPUAdd(NCPUs), info(NCPUs), stat=ierr)
         CPUAdd= 0
         CPUMark= -1
         info= -1
         CPUNo= 0


         !* determine which cpus to receive from (set CPUMark) and num_recvs,
         !* at the end of the loop CPUMark(i) contains the number of elements to
         !* be received from the i-1'th cpu

         if (NumColsOffd>0) OffdCol= ColMapOffd(1)
         NumRecvs=0
         i=1
         do while(i<=NumColsOffd)

            !* get the cpu # where offdcol'th column is in
            if (NumColsDiag>0) CPUNo= min(NCPUs-1, OffdCol/NumColsDiag)
            do while (ColStarts(CPUNo+1)> OffdCol)
               CPUNo= CPUNo- 1
            enddo
            do while (ColStarts(CPUNo+2)-1< OffdCol)
               CPUNo= CPUNo+ 1
            enddo

            NumRecvs= NumRecvs+ 1
            CPUMark(NumRecvs)= CPUNo
            j= i

            do while (ColStarts(CPUNo+2)>OffdCol)
               CPUAdd(NumRecvs)= CPUAdd(NumRecvs)+ 1  !< recieve how many data
               if (j<NumColsOffd) then
                  j= j+ 1
                  OffdCol= ColMapOffd(j) !< next non-zero column
               else !< j=NumColsOffd, the last non-zero column for offd
                  j= j+ 1
                  OffdCol= ColStarts(NCPUs+1)
               endif
            enddo
            if (j<NumColsOffd+1) then
               i= j
            else !< exceed NumColsOffd, the loop will be ended
               i=j+1
            endif
         enddo !{while(i<=NumColsOffd)}

         !* the length of tmp array, the 2 comes from two arrays, one is
         !* CPUMark, another is CPUAdd.
         localinfo= 2*NumRecvs
         SendRecv%NumRecvs= NumRecvs

#if defined (MPI)
         call mpi_allgather(localinfo, 1, mpi_in, info, 1, mpi_in, comm, ierr)
#endif

         !* generate information to be sent: tmp contains for each RecvCPUs:
         !* id of RecvCPUs, number of elements to be received for this
         !* processor,
         !* indices of elements (in this order)

         allocate(Displs(NCPUs+1), stat=ierr)
         Displs(1)=1
         do i=2, NCPUs+1
            Displs(i)=Displs(i-1)+ info(i-1)
         enddo
         allocate(RecvBuf(Displs(NCPUs+1)), stat=ierr)
         RecvBuf= -1

         tmp=> null()
         if (NumRecvs>0) then
            allocate(tmp(localinfo), stat=ierr)
            allocate(SendRecv%RecvCPUs(NumRecvs), stat=ierr)
            RecvCPUs=> SendRecv%RecvCPUs
            tmp= -1
            RecvCPUs= -1
         else
            allocate(tmp(1)) !* allocate this array is only for mpi_allgatherv
            tmp= -1
         endif
         allocate(SendRecv%RecvVecStarts(NumRecvs+1), stat=ierr)
         RecvVecStarts=> SendRecv%RecvVecStarts
         RecvVecStarts= -1

         !* put recvcpus, cpuadd to array tmp
         j=1
         if (NumRecvs>0) RecvVecStarts(1)=1
         do i=1, NumRecvs
            NumElements= CPUAdd(i)
            RecvCPUs(i)= CPUMark(i)
            RecvVecStarts(i+1)= RecvVecStarts(i)+ NumElements
            tmp(j)= RecvCPUs(i)
            tmp(j+1)= CPUAdd(i)
            j= j+ 2
         enddo

#if defined (MPI)
         !* every cpu knows what are needed for other cpus
         call mpi_allgatherv(tmp, localinfo, mpi_in, RecvBuf, &
            info, Displs, mpi_in, comm, ierr)
#endif

         if (stdout>0) then


         endif

         !* get send information from RecvBuf
         NumSends= 0
         NumElements= 0
         CPUAdd(1)= 1
         do i=1, NCPUs
            j= Displs(i)
            do while(j<Displs(i+1))
               j= j+ 1
               if ( RecvBuf(j)==cpu_id) then
                  NumSends= NumSends+ 1
                  CPUMark(NumSends)= i-1
                  CPUAdd(NumSends+1)= CPUAdd(NumSends)+ RecvBuf(j+1)
                  exit
               endif
               j= j+ 1
            enddo
         enddo
         SendRecv%NumSends= NumSends

         !* determine SendCPUs and actual elements to be send (in SendMapElements)
         !* and SendMapStarts whose i-th entry points to the beginning of the
         !* elements to be send to i'th cpu
         if (NumSends>0) then
            allocate(SendRecv%SendCPUs(NumSends), stat=ierr)
            allocate(SendRecv%SendMapElements(CPUAdd(NumSends+1)-1), stat=ierr)
            SendCPUs=> SendRecv%SendCPUs
            SendMapElements=> SendRecv%SendMapElements
            SendCPUs= -1
            SendMapElements= -1
         endif
         allocate(SendRecv%SendMapStarts(NumSends+1), stat=ierr)
         SendMapStarts=> SendRecv%SendMapStarts
         SendMapStarts= -1

         NumRequest= NumRecvs+ NumSends
         if (NumRequest>0) then
            allocate(mpirequest(NumRequest), stat=ierr)
#if defined (MPI)
            allocate(mpistatus(NumRequest*mpi_status_size), stat=ierr)
#endif
            mpirequest= -1
            mpistatus= -1
         endif

         if (NumSends>0) SendMapStarts(1)= 1
         do i=1, NumSends
            SendMapStarts(i+1)= CPUAdd(i+1)
            SendCPUs(i)= CPUMark(i)
         enddo


         j=1
         do i=1, NumSends
            VecStart= SendMapStarts(i)
            VecLen= SendMapStarts(i+1)- VecStart
            icpu= SendCPUs(i)
#if defined (MPI)
            call mpi_irecv(SendMapElements(VecStart), VecLen, mpi_in, &
               icpu, 0, comm, mpirequest(j), ierr)
#endif
            j= j+ 1
         enddo

         do i=1, NumRecvs
            VecStart= RecvVecStarts(i)
            VecLen= RecvVecStarts(i+1)- VecStart
            icpu= RecvCPUs(i)
#if defined (MPI)
            call mpi_isend(ColMapOffd(VecStart), VecLen, mpi_in, &
               icpu, 0, comm, mpirequest(j), ierr)
#endif
            j= j+ 1
         enddo

         if (NumRequest>0) then
            ierr=1
#if defined (MPI)
            call mpi_waitall(NumRequest, mpirequest, mpistatus, ierr)
#endif
            if (ierr.eq.0)deallocate(mpirequest, mpistatus)
         endif

         !* minus the offset
         if (NumSends>0) then
            SendMapElements= SendMapElements- FirstColDiag+ 1
         endif

         !* deallocate memories
         if (associated(tmp)) deallocate(tmp)
         if (associated(CPUAdd)) deallocate(CPUAdd)
         if (associated(CPUMark)) deallocate(CPUMark)
         if (associated(RecvBuf)) deallocate(RecvBuf)
         if (associated(Displs)) deallocate(Displs)
         if (associated(info)) deallocate(info)

         HParCSR%IsComm= .true.

         !* only for debug
         !if (cpuid.eq.3.and.NumRequest>0) then
         !   print *,'SendCPUs', SendCPUs
         !   print *,'SendMapStarts', SendMapStarts
         !   print *,'SendMapElements', SendMapElements
         !   print *,'RecvCPUs', RecvCPUs
         !   print *,'RecvVecStarts', RecvVecStarts
         !endif

         return

      end subroutine WTGenSendRecv

      !> free memory for sendrecv, this subroutine will be called when
      !> the Hamiltonian is destroyed
      subroutine WTSendRecvDestroy(sendrecv)

         implicit none
         type(WTParCSRComm) :: sendrecv

         if (associated(sendrecv%SendCPUs)) &
            deallocate(sendrecv%SendCPUs)
         if (associated(sendrecv%SendMapStarts)) &
            deallocate(sendrecv%SendMapStarts)
         if (associated(sendrecv%SendMapElements)) &
            deallocate(sendrecv%SendMapElements)
         if (associated(sendrecv%RecvCPUs)) &
            deallocate(sendrecv%RecvCPUs)
         if (associated(sendrecv%RecvVecStarts)) &
            deallocate(sendrecv%RecvVecStarts)
         nullify(sendrecv%SendCPUs)
         nullify(sendrecv%SendMapStarts)
         nullify(sendrecv%SendMapElements)
         nullify(sendrecv%RecvCPUs)
         nullify(sendrecv%RecvVecStarts)

         return
      end subroutine WTSendRecvDestroy

      !> perform vector-vector multiply in Parallel vector
      !> GlobalResult=<bra|ket>
      complex(dp) function WTParVecInnerProd(bra, ket) result(GlobalResult)
         implicit none
         type(WTParVec), intent(in) :: bra
         type(WTParVec), intent(in) :: ket
         complex(dp) :: LocalResult

         integer :: ierr
         type(WTSeqVec), pointer :: x
         type(WTSeqVec), pointer :: y

         x=> bra%SeqVec
         y=> ket%SeqVec

         LocalResult= WTSeqVecInnerProd(x, y)

#if defined (MPI)
         call mpi_barrier(bra%comm, ierr)
         call mpi_allreduce(LocalResult, GlobalResult, 1, &
            mpi_dc, mpi_sum, bra%comm , ierr)
#endif

         return

      end function WTParVecInnerProd

      !> perform vector-vector multiply in Parallel vector, however, with no
      !> mpi_allreduce operation
      !> GlobalResult=<bra|ket>
      complex(dp) function WTParVecInnerProdWithNoReduce(bra, ket) result(LocalResult)
         implicit none
         type(WTParVec), intent(in) :: bra
         type(WTParVec), intent(in) :: ket

         integer :: ierr
         type(WTSeqVec), pointer :: x
         type(WTSeqVec), pointer :: y

         x=> bra%SeqVec
         y=> ket%SeqVec

         LocalResult= WTSeqVecInnerProd(x, y)

         return

      end function WTParVecInnerProdWithNoReduce

      !> perform vector-vector multiply in sequential vector
      complex(dp) function WTSeqVecInnerProd(bra, ket) result(LocalResult)
         type(WTSeqVec), intent(in) :: bra
         type(WTSeqVec), intent(in) :: ket

         complex(dp), pointer :: x(:)
         complex(dp), pointer :: y(:)

         integer(li) :: i
         integer(li) :: length

         length= bra%length
         x=> bra%data
         y=> ket%data
         LocalResult= (0d0, 0d0)
         !$OMP PARALLEL DO REDUCTION (+:LocalResult)
         do i=1, length
            LocalResult= LocalResult+ conjg(x(i))* y(i)
         enddo
         !$OMP END PARALLEL DO

         return
      end function WTSeqVecInnerProd

      !> print parallel vector on screen
      subroutine WTSeqVecPrint(SeqVec)
         implicit none
         type(WTSeqVec), intent(in) :: SeqVec

         integer(li) :: i
         integer(li) :: length

         length= SeqVec%length
         do i=1, length
            write(*, '(i8, 2f12.6)')i, SeqVec%data(i)
         enddo
         return
      end subroutine WTSeqVecPrint

      !> print parallel vector on screen
      subroutine WTParVecPrint(ParVec)
         implicit none
         type(WTParVec), intent(in) :: ParVec

         integer(li) :: i
         integer(li) :: length

         length= ParVec%SeqVec%length
         do i=1, length
            write(*, '(i8, 2f12.6)')i, ParVec%SeqVec%data(i)
         enddo
         return
      end subroutine WTParVecPrint

      !> Get value from seqvec, only one element
      function WTSeqVecGetSingleValue(SeqVec, i) result(val)
         implicit none
         integer(li), intent(in) :: i
         complex(dp) :: val
         type(WTSeqVec), intent(in) :: SeqVec

         val= SeqVec%data(i)
         return
      end function WTSeqVecGetSingleValue

      !> get value from ParVec, only one element
      function WTParVecGetSingleValue(ParVec, i) result(val)
         implicit none
         integer(li), intent(in) :: i
         type(WTParVec), intent(in) :: ParVec
         complex(dp) :: val

         val= WTSeqVecGetSingleValue(ParVec%SeqVec, i)

         return
      end function WTParVecGetSingleValue

      !> set constant value for seqvec, only one element
      subroutine WTSeqVecSetSingleValue_z(SeqVec, i, val)
         implicit none
         integer(li), intent(in) :: i
         complex(dp), intent(in) :: val
         type(WTSeqVec) :: SeqVec

         SeqVec%data(i)= val

         return
      end subroutine WTSeqVecSetSingleValue_z

      !> set constant value for seqvec, only one element
      subroutine WTSeqVecSetSingleValue_r(SeqVec, i, val)
         implicit none
         integer(li), intent(in) :: i
         real(dp), intent(in) :: val
         type(WTSeqVec) :: SeqVec

         SeqVec%data(i)= val

         return
      end subroutine WTSeqVecSetSingleValue_r

      !> set constant value for Parvec, only one element
      subroutine WTParVecSetSingleValue_z(ParVec, i, val)
         implicit none
         integer(li), intent(in) :: i
         complex(dp), intent(in) :: val
         type(WTParVec) :: ParVec

         call WTSeqVecSetSingleValue_z(ParVec%SeqVec, i, val)

         return
      end subroutine WTParVecSetSingleValue_z

      !> set constant value for Parvec, only one element
      subroutine WTParVecSetSingleValue_r(ParVec, i, val)
         implicit none
         integer(li), intent(in) :: i
         real(dp), intent(in) :: val
         type(WTSeqVec) :: SeqVec
         type(WTParVec) :: ParVec

         call WTSeqVecSetSingleValue_r(ParVec%SeqVec, i, val)

         return
      end subroutine WTParVecSetSingleValue_r

      !> set constant value for seqvec
      subroutine WTSeqVecSetConstantValue_z(SeqVec, val)
         implicit none
         complex(dp), intent(in) :: val
         type(WTSeqVec) :: SeqVec

         integer(li) :: i
         integer(li) :: length
         complex(dp), pointer :: x(:)

         length= SeqVec%length
         x=> SeqVec%data

         !$OMP PARALLEL DO PRIVATE(i)
         do i=1, length
            x(i)= val
         enddo
         !$OMP END PARALLEL DO

         return
      end subroutine WTSeqVecSetConstantValue_z

      !> set constant value for seqvec
      subroutine WTSeqVecSetConstantValue_r(SeqVec, val)
         implicit none
         real(dp), intent(in) :: val
         type(WTSeqVec) :: SeqVec

         integer(li) :: i
         integer(li) :: length
         complex(dp), pointer :: x(:)

         length= SeqVec%length
         x=> SeqVec%data

         !$OMP PARALLEL DO PRIVATE(i)
         do i=1, length
            x(i)= val
         enddo
         !$OMP END PARALLEL DO

         return
      end subroutine WTSeqVecSetConstantValue_r

      !> set random values for seqvec
      subroutine WTSeqVecSetRandomValue(SeqVec)
         implicit none
         type(WTSeqVec) :: SeqVec

         integer(li) :: i
         integer(li) :: length
         real(dp)    :: rann
         complex(dp), pointer :: x(:)

         length= SeqVec%length
         x=> SeqVec%data

         !$OMP PARALLEL DO PRIVATE(i, rann)
         do i=1, length
            call random_number(harvest= rann)
            x(i)=  rann
         enddo
         !$OMP END PARALLEL DO

         return
      end subroutine WTSeqVecSetRandomValue

      !> set constant value for parallel vector
      subroutine WTParVecSetConstantValue_z(ParVec, val)
         implicit none
         type(WTParVec), intent(out) :: ParVec
         complex(dp), intent(in) :: val

         call WTSeqVecSetConstantValue_z(ParVec%SeqVec, val)

         return
      end subroutine WTParVecSetConstantValue_z

      !> set constant value for parallel vector
      subroutine WTParVecSetConstantValue_r(ParVec, val)
         implicit none
         type(WTParVec), intent(out) :: ParVec
         real(dp), intent(in) :: val

         call WTSeqVecSetConstantValue_r(ParVec%SeqVec, val)

         return
      end subroutine WTParVecSetConstantValue_r

      !> set random value for parallel vector
      subroutine WTParVecSetRandomValue(ParVec)
         implicit none
         type(WTParVec) :: ParVec

         call WTSeqVecSetRandomValue(ParVec%SeqVec)

         return
      end subroutine WTParVecSetRandomValue

      !> set sequential values for seqvec
      subroutine WTSeqVecSetSeqValue(first, last, SeqVec)
         implicit none
         integer(li), intent(in) :: first
         integer(li), intent(in) :: last
         type(WTSeqVec) :: SeqVec

         integer(li) :: i
         complex(dp), pointer :: x(:)

         x=> SeqVec%data

         !$OMP PARALLEL DO PRIVATE(i)
         do i=first, last
            x(i-first+1)=  dble(i)
         enddo
         !$OMP END PARALLEL DO

         return
      end subroutine WTSeqVecSetSeqValue

      !> set Sequential value for parallel vector
      subroutine WTParVecSetSeqValue(ParVec)
         implicit none
         type(WTParVec) :: ParVec

         integer(li) :: first
         integer(li) :: last

         first= ParVec%FirstIndex
         last= ParVec%LastIndex

         call WTSeqVecSetSeqValue(first, last, ParVec%SeqVec)

         return
      end subroutine WTParVecSetSeqValue

      !> set sequential vector copy data from x to y,  y= x
      subroutine WTSeqVecCopy(y, x)
         implicit none
         type(WTSeqVec), intent(out) :: y
         type(WTSeqVec), intent(in) :: x

         integer(li) :: i

         !$OMP PARALLEL DO PRIVATE(i)
         do i=1, x%length
            y%data(i)= x%data(i)
         enddo
         !$OMP END PARALLEL DO

         return
      end subroutine WTSeqVecCopy

      !> set parallel vector copy data from x to y,  y= x
      !> notice the length of y should larger then x
      subroutine WTParVecCopy(y, x)
         implicit none
         type(WTParVec), intent(out) :: y
         type(WTParVec), intent(in) :: x

         call WTSeqVecCopy(y%SeqVec, x%SeqVec)

         return
      end subroutine WTParVecCopy

      !> perform a*x, a is complex number, x is seqvec
      function WTSeqVecScale_z1(a, x) result(y)

         implicit none
         complex(dp), intent(in) :: a
         type(WTSeqVec), intent(in), pointer :: x
         type(WTSeqVec), pointer :: y

         integer(li) :: i

         !$OMP PARALLEL DO PRIVATE(i)
         do i=1, x%length
            x%data(i)= x%data(i)* a
         enddo
         !$OMP END PARALLEL DO
         y=> x

         return
      end function WTSeqVecScale_z1

      !> perform a*x, a is complex number, x is parvec
      function WTParVecScale_z1(a, x) result(y)

         implicit none
         complex(dp), intent(in) :: a
         type(WTParVec), intent(in), pointer :: x
         type(WTParVec), pointer :: y
         type(WTSeqVec), pointer :: ySeq

         ySeq=> WTSeqVecScale_z1(a, x%SeqVec)
         y=> x

         return
      end function WTParVecScale_z1

      !> perform a*x , a is real number
      function WTSeqVecScale_r1(a, x) result(y)

         implicit none
         real(dp), intent(in) :: a
         type(WTSeqVec), pointer :: x
         type(WTSeqVec), pointer :: y

         integer(li) :: i

         !$OMP PARALLEL DO PRIVATE(i)
         do i=1, x%length
            x%data(i)= x%data(i)* a
         enddo
         !$OMP END PARALLEL DO
         y=> x

         return
      end function WTSeqVecScale_r1

      !> perform a*x, a is real number, x is parvec
      function WTParVecScale_r1(a, x) result(y)

         implicit none
         real(dp), intent(in) :: a
         type(WTParVec), intent(in), pointer :: x
         type(WTParVec), pointer :: y
         type(WTSeqVec), pointer :: ySeq

         ySeq=> WTSeqVecScale_r1(a, x%SeqVec)
         y=> x

         return
      end function WTParVecScale_r1

      !> perform a*x, a is complex number, x is seqvec
      function WTSeqVecScale_z2(x, a) result(y)

         implicit none
         complex(dp), intent(in) :: a
         type(WTSeqVec), pointer :: x
         type(WTSeqVec), pointer :: y

         integer(li) :: i

         !$OMP PARALLEL DO PRIVATE(i)
         do i=1, x%length
            x%data(i)= x%data(i)* a
         enddo
         !$OMP END PARALLEL DO

         y=> x

         return
      end function WTSeqVecScale_z2

      !> perform a*x, a is complex number, x is parvec
      function WTParVecScale_z2(x, a) result(y)

         implicit none
         complex(dp), intent(in) :: a
         type(WTParVec), intent(in), pointer :: x
         type(WTParVec), pointer :: y
         type(WTSeqVec), pointer :: ySeq

         ySeq=> WTSeqVecScale_z2(x%SeqVec, a)

         y=> x
         return
      end function WTParVecScale_z2

      !> perform a*x , a is real number
      function WTSeqVecScale_r2(x, a) result(y)

         implicit none
         real(dp), intent(in) :: a
         type(WTSeqVec), pointer :: x
         type(WTSeqVec), pointer :: y

         integer(li) :: i

         !$OMP PARALLEL DO PRIVATE(i)
         do i=1, x%length
            x%data(i)= x%data(i)* a
         enddo
         !$OMP END PARALLEL DO
         y=> x

         return
      end function WTSeqVecScale_r2

      !> perform a*x, a is real number, x is parvec
      function WTParVecScale_r2(x, a) result(y)

         implicit none
         real(dp), intent(in) :: a
         type(WTParVec), intent(in), pointer :: x
         type(WTParVec), pointer :: y
         type(WTSeqVec), pointer :: ySeq

         ySeq=> WTSeqVecScale_r2(x%SeqVec, a)
         y=> x

         return
      end function WTParVecScale_r2

      !> calculate norm for a sequential vector
      function WTSeqVecNorm(x) result(norm)
         implicit none

         real(dp) :: norm
         complex(dp) :: norm_z
         complex(dp) :: norm_zmpi
         type(WTSeqVec), intent(in) :: x

         integer :: ierr
         integer(li) :: i

         norm_z= zero
         do i=1, x%length
            norm_z= norm_z+ conjg(x%data(i))*x%data(i)
         enddo
#if defined (MPI)
         call mpi_allreduce(Norm_z, norm_zmpi, 1, &
            mpi_dc, mpi_sum, mpi_cmw , ierr)
#endif
         norm  = dble(norm_zmpi)
         norm  = dsqrt(norm)

         return
      end function WTSeqVecNorm

      !> calculate norm for a parallel vector
      function WTParVecNorm(x) result(norm)
         implicit none

         real(dp) :: norm
         complex(dp) :: norm_z
         type(WTParVec), intent(in) :: x

         norm_z= x*x
         norm  = dble(norm_z)
         norm  = dsqrt(norm)

         return
      end function WTParVecNorm

      !> normlize a vector
      subroutine WTParVecNormalize(x)
         implicit none

         type(WTParVec), pointer :: x

         real(dp) :: norm

         norm= WTParVecNorm(x)

         norm= 1/norm
         !x=> norm*x
         x=> WTParVecScale_r1(norm, x)

         return
      end subroutine WTParVecNormalize

      !> define two parallel vectors' add
      subroutine WTParVecAxpy(a, x, y)

         implicit none
         complex(dp), intent(in) :: a
         type (WTParVec), intent(in) :: x
         type (WTParVec), intent(inout):: y

         call WTSeqVecAxpy(a, x%SeqVec, y%SeqVec)

         return
      end subroutine WTParVecAxpy

      !> define two sequential vectors' add
      subroutine WTSeqVecAxpy(a, x, y)

         implicit none
         complex(dp), intent(in) :: a
         type (WTSeqVec), intent(in) :: x
         type (WTSeqVec), intent(inout) :: y

         integer(li) :: i

         !$OMP PARALLEL DO PRIVATE(i)
         do i=1, x%length
            y%data(i)= y%data(i)+ a*x%data(i)
         enddo
         !$OMP END PARALLEL DO

         return
      end subroutine WTSeqVecAxpy

      !> define two parallel vectors' add
      subroutine WTParVecAxpBy(a, x, b, y)

         implicit none
         complex(dp), intent(in) :: a
         complex(dp), intent(in) :: b
         type (WTParVec), intent(in) :: x
         type (WTParVec), intent(inout):: y

         call WTSeqVecAxpBy(a, x%SeqVec, b, y%SeqVec)

         return
      end subroutine WTParVecAxpBy

      !> define two sequential vectors' add
      subroutine WTSeqVecAxpBy(a, x, b, y)

         implicit none
         complex(dp), intent(in) :: a
         complex(dp), intent(in) :: b
         type (WTSeqVec), intent(in) :: x
         type (WTSeqVec), intent(inout) :: y

         integer(li) :: i

         !$OMP PARALLEL DO PRIVATE(i)
         do i=1, x%length
            y%data(i)= b*y%data(i)+ a*x%data(i)
         enddo
         !$OMP END PARALLEL DO

         return
      end subroutine WTSeqVecAxpBy

      !> perform matrix-vector multiplication in csr format
      subroutine WTCSRMatVec_c(alpha, HCSR, vec_in, beta, vec_out)

         implicit none

         complex(dp), intent(in) :: alpha
         complex(dp), intent(in) :: beta
         type(WTCSR), intent(in) :: HCSR
         type(WTSeqVec), intent(in) :: vec_in
         type(WTSeqVec), intent(out) :: vec_out

         integer(li) :: i, j
         integer(li) :: NumRows
         complex(dp) :: tempy

         NumRows= HCSR%NumRows

         if (abs(beta)<eps12) then
            !$OMP PARALLEL DO PRIVATE(i, j, tempy)
            do i=1, NumRows
               tempy= (0d0, 0d0)
               do j= HCSR%ia(i), HCSR%ia(i+1)-1
                  tempy= tempy+ HCSR%a(j)*vec_in%data(HCSR%ja(j))
               enddo
               vec_out%data(i)= tempy*alpha
            enddo
            !$OMP END PARALLEL DO
         elseif (abs(beta-1d0)<eps12) then
            !$OMP PARALLEL DO PRIVATE(i, j, tempy)
            do i=1, NumRows
               tempy= (0d0, 0d0)
               do j= HCSR%ia(i), HCSR%ia(i+1)-1
                  tempy= tempy+ HCSR%a(j)*vec_in%data(HCSR%ja(j))
               enddo
               vec_out%data(i)=vec_out%data(i)+ tempy*alpha
            enddo
            !$OMP END PARALLEL DO
         elseif (abs(beta-1d0)<eps12) then
            !$OMP PARALLEL DO PRIVATE(i, j, tempy)
            do i=1, NumRows
               tempy= (0d0, 0d0)
               do j= HCSR%ia(i), HCSR%ia(i+1)-1
                  tempy= tempy+ HCSR%a(j)*vec_in%data(HCSR%ja(j))
               enddo
               vec_out%data(i)=vec_out%data(i)*beta+ tempy*alpha
            enddo
            !$OMP END PARALLEL DO
         endif

         return
      end subroutine WTCSRMatVec_c

      !> perform matrix-vector multiplication in Parcsr format
      subroutine WTParCSRMatVec_c(HParCSR, vec_in, vec_out)
         implicit none

         !* inout variables
         type(WTParCSR), intent(in) :: HParCSR
         type(WTParVec), intent(in) :: vec_in
         type(WTParVec), intent(out) :: vec_out

         !* local variables

         type(WTCSR), pointer :: diag
         type(WTCSR), pointer :: offd
         type(WTSeqVec), pointer :: OffdVec
         type(WTSeqVec), pointer :: x
         type(WTSeqVec), pointer :: y

         complex(dp) :: alpha
         complex(dp) :: beta

         !> the data that should be sent to other cpus
         complex(dp), pointer :: SendData(:)

         !> the data that recieved from other cpus
         complex(dp), pointer :: RecvData(:)

         !> local data
         complex(dp), pointer :: LocalData(:)

         !> sendrecv data
         type(WTCommHandle), pointer :: CommHandle
         type(WTParCSRComm), pointer :: SendRecv
         integer, pointer :: SendCPUs(:)
         integer, pointer :: SendMapStarts(:)
         integer(li), pointer :: SendMapElements(:)
         integer, pointer :: RecvCPUs(:)
         integer, pointer :: RecvVecStarts(:)
         integer(li) :: NumColsOffd

         integer :: i, j, k
         integer :: NumSends
         integer :: VecStart
         integer :: ierr

         real(dp) :: time1
         real(dp) :: time2

         diag=> Null()
         offd=> Null()
         OffdVec=> Null()
         x=> Null()
         y=> Null()
         RecvData=> Null()
         SendData=> Null()
         LocalData=> Null()
         SendCPUs=> Null()
         SendMapStarts=> Null()
         SendMapElements=> Null()
         RecvCPUs=> Null()
         RecvVecStarts=> Null()

         diag=> HParCSR%Diag
         offd=> HParCSR%Offd
         x=> vec_in%SeqVec
         y=> vec_out%SeqVec
         LocalData=> x%data
         NumColsOffd= offd%NumCols

         if (.not.HParCSR%IsComm) then
            write(*,*)'you should generate sendrecv befor this subroutine calls'
            stop
         endif


         SendRecv=> HParCSR%SendRecv
         SendCPUs=> SendRecv%SendCPUs
         SendMapStarts=> SendRecv%SendMapStarts
         SendMapElements=> SendRecv%SendMapElements
         RecvCPUs=> SendRecv%RecvCPUs
         RecvVecStarts=> SendRecv%RecvVecStarts


         !> get send data from the vector on the current cpu
         NumSends= SendRecv%NumSends
         allocate(SendData(SendMapStarts(NumSends+1)-1))

         !> create a sequential vector to store the recieved data
         if (.not.associated(OffdVec)) allocate(OffdVec)
         call WTSeqVecCreate(NumColsOffd, OffdVec)
         call WTSeqVecInitialize(OffdVec)
         RecvData=> OffdVec%data

         k=1
         do i=1, NumSends
            VecStart= SendMapStarts(i)
            do j=VecStart, SendMapStarts(i+1)-1
               SendData(k)= LocalData(SendMapElements(j))
               k= k+ 1
            enddo
         enddo

#if defined (MPI)
         call mpi_barrier(vec_in%comm, ierr)
#endif
         !> send and recieve data from and to other cpus
         !> this is a non-block communication,
         !> the MatVec for Diag will be executed no matter
         !> the sendrecv operation finish or not
         allocate(CommHandle)
         call WTCommHandleCreate(SendRecv, SendData, RecvData, CommHandle)

         !> only for test, if not, we should place the below three lines to
         !> after the next three line
         !> wait all the cpus finish their communications
         call WTCommHandleDestroy(CommHandle)
         Nullify(CommHandle)
#if defined (MPI)
         call mpi_barrier(vec_in%comm, ierr)
#endif

         alpha= (1d0, 0d0)
         beta= (0d0, 0d0)
         call WTCSRMatVec_c(alpha, HParCSR%Diag, x, beta,  y)

         alpha= (1d0, 0d0)
         beta= (1d0, 0d0)
         call WTCSRMatVec_c(alpha, HParCSR%Offd, OffdVec, beta, y)

         if (associated(SendData)) deallocate(SendData)
         call WTSeqVecDestroy(OffdVec)
         Nullify(SendData, OffdVec)

         return
      end subroutine WTParCSRMatVec_c

      subroutine csr_sort_indices(n, nnz, Ap, Aj, Ax)
         ! Sort CSR column indices inplace
         use para, only : dp, stdout
         implicit none
         integer, intent(in) :: n
         integer, intent(in) :: nnz
         integer, intent(inout) :: Ap(n+1), Aj(nnz)
         complex(dp), intent(inout) :: Ax(nnz)
         integer :: i, r1, r2, l
         integer, allocatable :: idx(:)

         allocate(idx(nnz))

         if (cpuid.eq.0) then
         !  write(stdout, '(a, f16.1, a)') " Memory cost in csr_sort_indices : ", nnz/1024d0/1024d0*4d0, 'MB'
         endif

         do i = 1, N
            r1 = Ap(i)
            r2 = Ap(i+1)-1
            l = r2-r1+1
            call iargsort(l, Aj(r1:r2), idx(1:l))
            Aj(r1:r2) = Aj(r1+idx(1:l)-1)
            Ax(r1:r2) = Ax(r1+idx(1:l)-1)
         end do

         deallocate(idx)

         return
      end subroutine csr_sort_indices

      subroutine iargsort(n, a, b)
         implicit none
         ! Returns the indices that would sort an array.
         !
         ! Arguments
         ! ---------
         !
         integer, intent(in) :: n
         integer, intent(in) :: a(n)    ! array of numbers
         integer, intent(inout) :: b(n)         ! indices into the array 'a' that sort it
         !
         ! Example
         ! -------
         !
         ! iargsort([10, 9, 8, 7, 6])   ! Returns [5, 4, 3, 2, 1]

         integer :: i,imin,relimin(1)           ! indices: i, i of smallest, relative imin
         integer :: temp                        ! temporary
         integer, allocatable :: a2(:)

         allocate(a2(n))

         a2 = a
         do i = 1, N
            b(i) = i
         end do

         do i = 1, N-1
            ! find ith smallest in 'a'
            relimin = minloc(a2(i:))
            imin = relimin(1) + i - 1
            ! swap to position i in 'a' and 'b', if not already there
            if (imin /= i) then
               temp = a2(i); a2(i) = a2(imin); a2(imin) = temp
               temp = b(i); b(i) = b(imin); b(imin) = temp
            end if
         end do

         deallocate(a2)

         return
      end subroutine iargsort

      subroutine csr_sum_duplicates(n, nnz, ia, ja, A_val)
         ! Sum together duplicate column entries in each row of CSR matrix A
         ! The column indicies within each row must be in sorted order.
         ! Explicit zeros are retained.
         ! ia, ja, and A_val will be modified *inplace*
         !> ia row index dimension(N+1)
         !> ja column index dimension(nnz)
         !> A_val value of A dimension(nnz)
         !> n  rank of matrix A
         use para, only : dp
         implicit none
         integer, intent(in) :: n
         integer, intent(inout) :: nnz
         integer, intent(inout) :: ia(n+1), ja(nnz)
         complex(dp), intent(inout) :: A_val(nnz)

         integer :: i1, i2, i, j, jj
         complex(dp) :: valu

         nnz = 1
         i2 = 1
         do i = 1, n
            i1 = i2
            i2 = ia(i+1)
            jj = i1
            do while (jj < i2)
               j = ja(jj)
               valu = A_val(jj)
               jj = jj + 1
               do while (jj < i2)
                  if (ja(jj) == j) then
                     valu = valu + A_val(jj)
                     jj = jj + 1
                  else
                     exit
                  end if
               end do
               ja(nnz) = j
               A_val(nnz) = valu
               nnz = nnz + 1
            end do
            ia(i+1) = nnz
         end do

         !> modified by qswu
         nnz = nnz- 1

         return
      end subroutine csr_sum_duplicates

      !> csrmv_z multiplies a CSR matrix A times a vector x; y=A*x
      !> inputs:
      !> ndim, integer, the row dimension of the matrix
      !> nnz, integer, number of non-zero entries, the dimension of A_csr,
      !ja_csr
      !> complex A_csr(nnz), integer ia_csr(ndim+1), ja_csr(nnz), the matrix in
      !CSR Compresed Sparse Row format
      !> complex x(ndim)
      !> outputs:
      !> complex y(ndim)
      subroutine csrmv_z(ndim, nnz, A_csr, ia_csr, ja_csr, x, y)

         use para, only : dp
         implicit none

         integer, intent(in) :: ndim
         integer, intent(in) :: nnz
         complex(dp), intent(in) :: A_csr(nnz)
         integer, intent(in) :: ia_csr(ndim+1)
         integer, intent(in) :: ja_csr(nnz)
         complex(dp), intent(in) :: x(ndim)
         complex(dp), intent(out) :: y(ndim)

         integer :: i, k
         complex(dp) :: t_z

!        real(dp) :: time1, time2

!        call now(time1)
#if defined (INTELMKL)
         call mkl_zcsrgemv('N', ndim, a_csr, ia_csr, ja_csr, x, y)
#else
      !> A naive implementation
         do i=1, ndim
            t_z= 0d0
            do k=ia_csr(i), ia_csr(i+1)-1
               t_z= t_z+ A_csr(k) * x(ja_csr(k))
               
            enddo
            y(i) = t_z
         enddo
#endif

!        call now(time2)
!        time_total_debug= time_total_debug+ time2- time1
         return

      end subroutine csrmv_z

      !> coomv_z multiplies a COO matrix A times a vector x; y=A*x
      !> inputs:
      !> ndim, integer, the row dimension of the matrix
      !> nnz, integer, number of non-zero entries, the dimension of A_COO,
      !ja_COO
      !> complex A_COO(nnz), integer ia_COO(ndim+1), ja_COO(nnz), the matrix in
      !COO Compresed Sparse Row format
      !> complex x(ndim)
      !> outputs:
      !> complex y(ndim)
      subroutine coomv_z(ndim, nnz, A_coo, ia_coo, ja_coo, x, y)

         use para, only : dp
         implicit none

         integer, intent(in) :: ndim
         integer, intent(in) :: nnz
         complex(dp), intent(in) :: A_coo(nnz)
         integer, intent(in) :: ia_coo(ndim+1)
         integer, intent(in) :: ja_coo(nnz)
         complex(dp), intent(in) :: x(ndim)
         complex(dp), intent(out) :: y(ndim)

         integer :: i, k
         complex(dp) :: t_z

!        real(dp) :: time1, time2

!        call now(time1)
#if defined (INTELMKL)
         call mkl_zcoogemv('N', ndim, a_coo, ia_coo, ja_coo, nnz, x, y)
#else
      !> A naive implementation
         do i=1, nnz
            y(ia_coo(i)) = y(ia_coo(i))+ a_coo(i)* x(ja_coo(i))
         enddo
#endif

!        call now(time2)
!        time_total_debug= time_total_debug+ time2- time1
         return

      end subroutine coomv_z




      subroutine arpack_sparse_coo_eigs(ndims,nnzmax,nnz,acoo,jcoo,icoo,neval,nvecs,deval,sigma,zeigv, ritzvec)
         use para, only : dp
         implicit none
         ! external arguments
         ! dimension of matrix A
         integer, intent(inout) :: ndims

         ! maximum number of non-zero elements in matrix A
         integer, intent(in) :: nnzmax

         ! Number of non-zero elements in matrix A in COO format
         integer, intent(inout) :: nnz

         ! coordinate format storage of matrix A
         complex(dp), intent(inout) :: acoo(nnzmax)
         integer, intent(inout) :: jcoo(nnzmax)
         integer, intent(inout) :: icoo(nnzmax)

         ! number of selected eigenvals
         integer, intent(in) :: neval

         ! number of Arnoldi vectors
         integer, intent(in) :: nvecs

         !> calculate eigenvector or not
         logical, intent(in) :: ritzvec

         ! eigenvalues for selected "which"
         real(dp), intent(out) :: deval(neval)

         ! eigenvector for selected "which"

         complex(dp), intent(in) :: sigma


         ! compressed sparse row storage of matrix A
         ! if we need to write out zeigv, then the dimension is (ndims, neval)
         ! otherwise, it will not be used
         complex(dp),intent(out) :: zeigv(ndims, nvecs)
        
         zeigv= 0d0

#if defined (INTELMKL)
         !> zndrv2 needs a sparse solver to solve (A-sigma*I)*x=B with given A, sigma, and B
         !> get eigenvalues of a sparse matrix by calling arpack subroutine
         !> acoo, jcoo, icoo would be converted in to A-sigma*I, then converted into CSR format
         !> usually zndrv1 is about 10 times faster then zndrv2
         if (arpack_solver=='zndrv2') then
            call zmat_arpack_zndrv2(ndims, nnzmax, nnz, acoo, jcoo, icoo, sigma, neval, nvecs, deval, zeigv, ritzvec)
         else
            call zmat_arpack_zndrv1(ndims, nnzmax, nnz,  acoo, jcoo, icoo, sigma, neval, nvecs, deval, zeigv, ritzvec)
         endif
#else
         !> here acoo, icoo, jcoo are stored in COO format
         !> use matrix vector multiplication
         !> zndrv1 needs a matrix vector multiplication operator A*x
         call zmat_arpack_zndrv1(ndims, nnzmax, nnz,  acoo, jcoo, icoo, sigma, neval, nvecs, deval, zeigv, ritzvec)
#endif


         return
      end subroutine arpack_sparse_coo_eigs


      subroutine arpack_sparse_coo_eigs_nonorth(ndims, nnzmax, nnz, acoo_k, jcoo_k, icoo_k, &
             snnzmax, snnz, sacoo_k, sjcoo_k, sicoo_k, neval,nvecs,deval,sigma,zeigv, ritzvec)
         !> we want to get the eigen value of matrix A
         !> by solving A*x = lambda*S*x
         !> where A is a hermitian matrix, S is a overlap matrix between the non-orthogonal basis
         use para, only : dp
         implicit none
         ! external arguments
         ! dimension of matrix A
         integer, intent(inout) :: ndims

         ! maximum number of non-zero elements in matrix A
         integer, intent(in) :: nnzmax
         integer, intent(in) :: snnzmax

         ! Number of non-zero elements in matrix A in COO format
         integer, intent(inout) :: nnz
         integer, intent(inout) :: snnz

         ! coordinate format storage of matrix A
         complex(dp), intent(inout) :: acoo_k(nnzmax)
         !> column indices
         integer, intent(inout) :: jcoo_k(nnzmax)
         !> row indices
         integer, intent(inout) :: icoo_k(nnzmax)

         ! coordinate format storage of overlap matrix S
         complex(dp), intent(inout) :: sacoo_k(snnzmax)
         !> column indices
         integer, intent(inout) :: sjcoo_k(snnzmax)
         !> row indices
         integer, intent(inout) :: sicoo_k(snnzmax)

         ! number of selected eigenvals
         integer, intent(in) :: neval

         ! number of Arnoldi vectors
         integer, intent(in) :: nvecs

         !> calculate eigenvector or not
         logical, intent(in) :: ritzvec

         ! eigenvalues for selected "which"
         real(dp), intent(out) :: deval(neval)

         ! eigenvector for selected "which"

         complex(dp), intent(in) :: sigma


         ! compressed sparse row storage of matrix A
         ! if we need to write out zeigv, then the dimension is (ndims, neval)
         ! otherwise, it will not be used
         complex(dp),intent(out) :: zeigv(ndims, nvecs)
        
         zeigv= 0d0
         !> get eigenvalues of a sparse matrix by calling arpack subroutine
         !> here acoo, icoo, jcoo are stored in COO format
!        call zmat_arpack_zndrv3(ndims, nnzmax, nnz,  acoo_k, jcoo_k, icoo_k, &
!           snnzmax, snnz, sacoo_k, sjcoo_k, sicoo_k, sigma, neval, nvecs, deval, zeigv)

         !> acoo, jcoo, icoo would be converted in to A-sigma*S, then converted into CSR format
         call zmat_arpack_zndrv4(ndims, nnzmax, nnz,  acoo_k, jcoo_k, icoo_k, &
            snnzmax, snnz, sacoo_k, sjcoo_k, sicoo_k, sigma, neval, nvecs, deval, zeigv)


         return
      end subroutine arpack_sparse_coo_eigs_nonorth

      subroutine zmat_arpack_zndrv4(ndims, nnzmax, nnz, acsr, jcsr, icsr, &
            snnzmax, snnz, sacsr, sjcsr, sicsr, sigma, neval, nvecs, deval, zeigv)
         !> We want to solve A*x = lambda*S*x in shift-invert mode 
         !> where A is a hermitian matrix
         !> S is the overlap matrix between non-orthogonal basis
         use para, only : dp, stdout, cpuid
         implicit none

! external arguments
! dimension of matrix A
         integer, intent(in) :: ndims

! maximum number of non-zero elements of matrix A
         integer, intent(in) :: nnzmax

! number of non-zero elements of matrix A
         integer, intent(inout) :: nnz

! sparse storage of matrix A
!> the input sparse storage format is coo 
         complex(dp), intent(inout) :: acsr(nnzmax)
         integer, intent(inout) :: jcsr(nnzmax)
         integer, intent(inout) :: icsr(nnzmax)

! optional for overlap matrix, when Orthogonal_Basis=T, this can be omited
!> the input sparse storage format is coo 
         integer, intent(in) :: snnzmax
         integer, intent(inout) :: snnz
         integer, intent(inout) :: sjcsr(snnzmax)
         integer, intent(inout) :: sicsr(snnzmax)
         complex(dp), intent(inout) :: sacsr(snnzmax)


! shift energy
         complex(dp), intent(in) :: sigma
! number of selected eigenvals
         integer, intent(in) :: neval

! number of Arnoldi vectors
         integer, intent(in) :: nvecs

! eigenvalues for selected "which"
         real(dp), intent(out) :: deval(neval)

! eigenvector for selected "which"
         complex(dp), intent(out) :: zeigv(ndims, nvecs)

! loop index over neval
         integer :: ival, jval, i

! axuiliary integer variables
         integer :: itmp, iter

         real(dp) :: time_cost_mv, time_cost_zgesv, time_cost_convert,  time_cost_znaupd
         real(dp) :: time_start, time_end

! auxiliary real(dp) variables
         real(dp) :: dtmp

!
!     %--------------%
!     | Local Arrays |
!     %--------------%
!
         integer :: iparam(11), ipntr(14)
         complex(dp), allocatable :: ax(:)
         complex(dp), allocatable :: mx(:)
         complex(dp), allocatable :: d(:)
         complex(dp), allocatable :: v(:, :)
         complex(dp), allocatable :: workd(:)
         complex(dp), allocatable :: workev(:)
         complex(dp), allocatable :: resid(:)
         complex(dp), allocatable :: workl(:)
         logical   , allocatable :: select(:)

         integer, allocatable :: iwk(:)

         real(dp), allocatable :: rwork(:)
         real(dp), allocatable :: rd(:, :)

!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
         character  :: bmat*1, which*2
         integer    :: ido, n, nev, ncv, lworkl, info, j, &
            ierr, nconv, maxitr, ishfts, mode, ldv
         integer    :: istat
         integer    :: lwrkl, lwrkv, lwrkd
         real(dp)    :: tol
         logical    :: rvec
!
!     %-----------------------------%
!     | BLAS & LAPACK routines used |
!     %-----------------------------%
!
         real(dp), external :: dznrm2, dlapy2
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
!     %----------------------------------------------------%
!     | The number N is the dimension of the matrix.  A    |
!     | generalized eigenvalue problem is solved (BMAT =   |
!     | 'G').  NEV is the number of eigenvalues to be      |
!     | approximated.  The user can modify NEV, NCV, WHICH |
!     | to solve problems of different sizes, and to get   |
!     | different parts of the spectrum.  However, The     |
!     | following conditions must be satisfied:            |
!     |                    N <= MAXN,                      |
!     |                  NEV <= MAXNEV,                    |
!     |              NEV + 2 <= NCV <= MAXNCV              |
!     %----------------------------------------------------%
!
         n = ndims
         ldv = ndims
         nev = neval
         ncv = nvecs
         lwrkd = 3 * ndims; lwrkv = 3 * nvecs
         lwrkl = 3 * nvecs * nvecs  + 5 * nvecs


         ! allocate memory for arpack
         allocate(select(nvecs), stat=istat)
         allocate(ax(ndims), stat=istat)
         allocate(mx(ndims), stat=istat)
         allocate( d(neval), stat=istat)
        !allocate( v(ndims, nvecs), stat=istat)
         allocate(workd( lwrkd), stat=istat)
         allocate(workl( lwrkl), stat=istat)
         allocate(resid( ndims), stat=istat)
         allocate(workev(lwrkv), stat=istat)
         allocate(rwork(nvecs), stat=istat)
         allocate(rd(nvecs, 3), stat=istat)
         allocate(iwk(ndims+1))
         iwk= 0

         bmat  = 'G'
         which = 'LM'


         time_cost_mv= 0d0; time_cost_zgesv= 0d0; time_cost_convert= 0d0
         time_cost_znaupd= 0d0

         !
         !     %----------------------------------------------------%
         !     | Construct C = A - SIGMA*S                          |
         !     %----------------------------------------------------%
         do i=1, splen_overlap_input
            j=i+nnz
            icsr(j)= sicsr(i)
            jcsr(j)= sjcsr(i)
            acsr(j)= -sigma*sacsr(i)
         enddo
         !> number of non-zeros nnz will be updated
         nnz= nnz+ splen_overlap_input
         if (nnz>nnzmax) stop "ERROR: in zndrv4, nnz is larger than nnzmax"


         call now(time_start)
         !> prepare hamiltonian
         !> transform coo format to csr format
         call ConvertCooToCsr(ndims, nnz, acsr, icsr, jcsr, iwk)
         call csr_sort_indices(ndims, nnz, icsr, jcsr, acsr)
         !> eleminate the same entries in the sparse matrix
         call csr_sum_duplicates(ndims, nnz, icsr, jcsr, acsr)
         call csr_sort_indices(ndims, nnz, icsr, jcsr, acsr)

         !> transform coo format to csr format
         call ConvertCooToCsr(ndims, snnz, sacsr, sicsr, sjcsr, iwk)
         call csr_sort_indices(ndims, snnz, sicsr, sjcsr, sacsr)
         !> eleminate the same entries in the sparse matrix
         call csr_sum_duplicates(ndims, snnz, sicsr, sjcsr, sacsr)
         call csr_sort_indices(ndims, snnz, sicsr, sjcsr, sacsr)
         call now(time_end)
         time_cost_convert= time_end- time_start

!
!     %---------------------------------------------------%
!     | The work array WORKL is used in ZNAUPD as         |
!     | workspace.  Its dimension LWORKL is set as        |
!     | illustrated below.  The parameter TOL determines  |
!     | the stopping criterion. If TOL<=0, machine        |
!     | precision is used.  The variable IDO is used for  |
!     | reverse communication, and is initially set to 0. |
!     | Setting INFO=0 indicates that a random vector is  |
!     | generated to start the ARNOLDI iteration.         |
!     %---------------------------------------------------%
!
         lworkl = 3 * ncv * ncv + 5 * ncv
         tol    = 1.0D-7
         ido    = 0
         info   = 0

!     %---------------------------------------------------%
!     | This program uses exact shifts with respect to    |
!     | the current Hessenberg matrix (IPARAM(1) = 1).    |
!     | IPARAM(3) specifies the maximum number of Arnoldi |
!     | iterations allowed. Mode 3 of ZNAUPD  is used      |
!     | (IPARAM(7) = 3).  All these options can be        |
!     | changed by the user. For details see the          |
!     | documentation in ZNAUPD .                          |
!     %---------------------------------------------------%
!
         ishfts = 1
         maxitr = 50000
         mode   = 3

         iparam(1) = ishfts
         iparam(3) = maxitr
         iparam(7) = mode
         nconv=-1
!
!     %-------------------------------------------%
!     | M A I N   L O O P (Reverse communication) |
!     %-------------------------------------------%
!
         iter = 0
10       continue
!
!        %---------------------------------------------%
!        | Repeatedly call the routine ZNAUPD and take |
!        | actions indicated by parameter IDO until    |
!        | either convergence is indicated or maxitr   |
!        | has been exceeded.                          |
!        %---------------------------------------------%
!
         call now(time_start)
#if defined (ARPACK)
         call znaupd ( ido, bmat, n, which, nev, tol, resid, ncv, &
            zeigv, ldv, iparam, ipntr, workd, workl, lworkl, rwork, info )
#else
         STOP "ERROR : Please install WannierTools with ARPACK since you are diagonalizing a large sparse matrix"
#endif
         call now(time_end)
         time_cost_znaupd=  time_cost_znaupd+ time_end- time_start

         iter = iter + 1
         if (mod(iter,100) .eq. 0 .and. cpuid==0) then
            write(stdout, '(A, I10)') '>>> Iteration with znaupd', iter
         endif

         if (ido .eq. -1) then
!
!           %-------------------------------------------%
!           | Perform  y <--- OP*x = inv[A-SIGMA*S]*S*x |
!           | to force starting vector into the range   |
!           | of OP.   The user should supply his/her   |
!           | own matrix vector multiplication routine  |
!           | and a linear system solver.  The matrix   |
!           | vector multiplication routine should take |
!           | workd(ipntr(1)) as the input. The final   |
!           | result should be returned to              |
!           | workd(ipntr(2)).                          |
!           |    x = workd(ipntr(1))
!           |    y = workd(ipntr(2))
!           %-------------------------------------------%
!
 
            !> first perform S*x
            call now(time_start)
            !call mkl_zcsrgemv('N', ndims, sacsr, sicsr, sjcsr, workd(ipntr(1)), workd(ipntr(2)))
            call csrmv_z(ndims, snnz, sacsr, sicsr, sjcsr, workd(ipntr(1)), workd(ipntr(2)))
            call now(time_end)
            time_cost_mv=  time_cost_mv+ time_end- time_start
            call zcopy (ndims, workd(ipntr(2)),  1, workd(ipntr(1)), 1)
            !> then perform inv[A-SIGMA*S]*S*x
            call now(time_start)
            !call zmat_mkldss_zgesv(ndims, nnz, acsr, jcsr, icsr, workd(ipntr(1)), workd(ipntr(2)))
            call sparse_solver(ndims, nnz, acsr, icsr, jcsr, workd(ipntr(1)), workd(ipntr(2)))
            call now(time_end)
            time_cost_zgesv= time_cost_zgesv+ time_cost_mv+ time_end- time_start

!           %-----------------------------------------%
!           | L O O P   B A C K to call ZNAUPD again. |
!           %-----------------------------------------%
!
            go to 10
         else if ( ido .eq. 1) then
!
!           %-----------------------------------------%
!           | Perform y <-- OP*x = inv[A-sigma*S]*S*x |
!           | S*x has been saved in workd(ipntr(3)).  |
!           | The user only need the linear system    |
!           | solver here that takes workd(ipntr(3))  |
!           | as input, and returns the result to     |
!           | workd(ipntr(2)).                        |
!           |    S*x=workd(ipntr(3))
!           |    y  =workd(ipntr(2))
!           %-----------------------------------------%
!
            call now(time_start)
            !call zmat_mkldss_zgesv(ndims, nnz, acsr, jcsr, icsr, workd(ipntr(3)), workd(ipntr(2)))
            call sparse_solver(ndims, nnz, acsr, icsr, jcsr, workd(ipntr(3)), workd(ipntr(2)))
            call now(time_end)
            time_cost_zgesv= time_cost_zgesv+ time_cost_mv+ time_end- time_start

            go to 10
!
         else if ( ido .eq. 2) then
!
!           %-------------------------------------%
!           |        Perform  y <--- S*x          |
!           | The matrix vector multiplication    |
!           | routine should take workd(ipntr(1)) |
!           | as input and return the result to   |
!           | workd(ipntr(2)).                    |
!           |    x = workd(ipntr(1))
!           |    y = workd(ipntr(2))
!           %-------------------------------------%
!
            call now(time_start)
            !call mkl_zcsrgemv('N', ndims, sacsr, sicsr, sjcsr, workd(ipntr(1)), workd(ipntr(2)))
            call csrmv_z(ndims, snnz, sacsr, sicsr, sjcsr, workd(ipntr(1)), workd(ipntr(2)))
            call now(time_end)
            time_cost_mv=  time_cost_mv+ time_end- time_start
!
!           %-----------------------------------------%
!           | L O O P   B A C K to call ZNAUPD  again. |
!           %-----------------------------------------%
!
            go to 10
!
         end if
!
!     %----------------------------------------%
!     | Either we have convergence or there is |
!     | an error.                              |
!     %----------------------------------------%
!
         if ( info .lt. 0 ) then
!
!        %--------------------------%
!        | Error message, check the |
!        | documentation in ZNAUPD  |
!        %--------------------------%

            if (cpuid==0) write(stdout, *) ' '
            if (cpuid==0) write(stdout, *) ' Error with _naupd, info = ', info
            if (cpuid==0) write(stdout, *) ' Check the documentation of _naupd'
            if (cpuid==0) write(stdout, *) ' '
            stop

         else
!
!        %-------------------------------------------%
!        | No fatal errors occurred.                 |
!        | Post-Process using ZNEUPD.                |
!        |                                           |
!        | Computed eigenvalues may be extracted.    |
!        |                                           |
!        | Eigenvectors may also be computed now if  |
!        | desired.  (indicated by rvec = .true.)    |
!        %-------------------------------------------%
!
            rvec = .false.
            if (LandauLevel_wavefunction_calc.or.SlabBand_calc) rvec = .true.

#if defined (ARPACK)
            call zneupd (rvec, 'A', select, d, zeigv, ldv, sigma, &
               workev, bmat, n, which, nev, tol, resid, ncv,&
               zeigv, ldv, iparam, ipntr, workd, workl, lworkl, &
               rwork, ierr)
#else
         STOP "ERROR : Please install WannierTools with ARPACK since you are diagonalizing a large sparse matrix"
#endif
!
!        %----------------------------------------------%
!        | Eigenvalues are returned in the one          |
!        | dimensional array D.  The corresponding      |
!        | eigenvectors are returned in the first NCONV |
!        | (=IPARAM(5)) columns of the two dimensional  |
!        | array V if requested.  Otherwise, an         |
!        | orthogonal basis for the invariant subspace  |
!        | corresponding to the eigenvalues in D is     |
!        | returned in V.                               |
!        %----------------------------------------------%
!
            if ( ierr .ne. 0) then
!
!           %------------------------------------%
!           | Error condition:                   |
!           | Check the documentation of ZNEUPD. |
!           %------------------------------------%
!
               if (cpuid==0) write(stdout, *) ' '
               if (cpuid==0) write(stdout, *) ' Error with _neupd, info = ', ierr
               if (cpuid==0) write(stdout, *) ' Check the documentation of _neupd. '
               if (cpuid==0) write(stdout, *) ' '

            else

               nconv = iparam(5)
               do 20 j=1, nconv
!
!               %---------------------------%
!               | Compute the residual norm |
!               |                           |
!               |   ||  A*x - lambda*S*x ||   |
!               |                           |
!               | for the NCONV accurately  |
!               | computed eigenvalues and  |
!               | eigenvectors.  (iparam(5) |
!               | indicates how many are    |
!               | accurate to the requested |
!               | tolerance)                |
!               %---------------------------%
!

                  !call mkl_zcsrgemv('N', ndims, sacsr, sicsr, sjcsr, zeigV(1,j), mx)
                  call csrmv_z(ndims, snnz, sacsr, sicsr, sjcsr, zeigV(1,j), mx)
                  !call mkl_zcsrgemv('N', ndims, acsr, icsr, jcsr, zeigV(1,j), ax)
                  call csrmv_z(ndims, nnz, acsr, icsr, jcsr, zeigV(1,j), ax)
                  call zaxpy(n, -d(j), mx, 1, ax, 1)
                  rd(j,1) = dble(d(j))
                  rd(j,2) = dimag(d(j))
                  rd(j,3) = dznrm2(n, ax, 1)
                  rd(j,3) = rd(j,3) / dlapy2(rd(j,1),rd(j,2))
20             continue
!
!            %-----------------------------%
!            | Display computed residuals. |
!            %-----------------------------%
!
!              call dmout(6, nconv, 3, rd, nvecs, -6, &
!                 'Ritz values (Real, Imag) and relative residuals')
            end if
!
!        %-------------------------------------------%
!        | Print additional convergence information. |
!        %-------------------------------------------%
!
            if ( info .eq. 1) then
               if (cpuid==0) write(stdout, *) ' '
               if (cpuid==0) write(stdout, *) ' Maximum number of iterations reached.'
               if (cpuid==0) write(stdout, *) ' '
            else if ( info .eq. 3) then
               if (cpuid==0) write(stdout, *) ' '
               if (cpuid==0) write(stdout, *) ' No shifts could be applied during implicit &
               &Arnoldi update, try increasing NCV.'
               if (cpuid==0) write(stdout, *) ' '
            end if

            !print *, ' '
            if (cpuid==0) write(stdout, *) 'ZNDRV4'
            if (cpuid==0) write(stdout, *) '====== '
            if (cpuid==0) write(stdout, *) ' '
            if (cpuid==0) write(stdout, *) ' Size of the matrix is ', n
            if (cpuid==0) write(stdout, *) ' Get eigenvalues around sigma= ', real(sigma)
            if (cpuid==0) write(stdout, *) ' The number of Ritz values requested is ', nev
            if (cpuid==0) write(stdout, *) ' The number of Arnoldi vectors generated', &
               ' (NCV) is ', ncv
            if (cpuid==0) write(stdout, *) ' What portion of the spectrum: ', which
            if (cpuid==0) write(stdout, *) ' The number of converged Ritz values is ', &
               nconv
            if (cpuid==0) write(stdout, *) ' The number of Implicit Arnoldi update', &
               ' iterations taken is ', iparam(3)
            if (cpuid==0) write(stdout, *) ' The number of OP*x is ', iparam(9)
            if (cpuid==0) write(stdout, *) ' The convergence criterion is ', tol
            if (cpuid==0) write(stdout, *) ' '
            if (cpuid==0) write(stdout, '(a,f8.4, " s")') ' time_cost_zgesv: ', time_cost_zgesv
            if (cpuid==0) write(stdout, '(a,f8.4, " s")') ' time_cost_mv: ', time_cost_mv
            if (cpuid==0) write(stdout, '(a,f8.4, " s")') ' time_cost_t1: ', time_cost_t1
            if (cpuid==0) write(stdout, '(a,f8.4, " s")') ' time_cost_t2: ', time_cost_t2
            if (cpuid==0) write(stdout, '(a,f8.4, " s")') ' time_cost_t3: ', time_cost_t3
            if (cpuid==0) write(stdout, '(a,f8.4, " s")') ' time_cost_convert', time_cost_convert
            if (cpuid==0) write(stdout, '(a,f8.4, " s")') ' time_cost_znaupd:', time_cost_znaupd

         end if
!
!     %---------------------------%
!     | Done with program zndrv4. |
!     %---------------------------%
!

         do ival=1,neval
            deval(ival) = real(d(ival))
           !zeigv(1:ndims, ival) = V(1:ndims, ival)
         enddo ! over ival={1,neval} loop

!---------------------------------------------------------------------!
! use selection sort to minimize swaps of eigenvectors, ref: dsteqr.f !
!---------------------------------------------------------------------!
         do ival=1,neval-1
            itmp = ival; dtmp = deval(itmp)

            do jval=ival+1,neval
               if ((deval(jval)+1.0D-12) .lt. dtmp) then
                  itmp = jval; dtmp = deval(itmp)
               endif
            enddo ! over jval={ival+1,neval} loop

            if (itmp .ne. ival) then
               deval(itmp) = deval(ival); deval(ival) = dtmp
               call zswap(ndims, zeigv(1, ival), 1, zeigv(1,itmp), 1)
            endif
         enddo ! over ival={1,neval-1} loop

         return
      end subroutine zmat_arpack_zndrv4


      subroutine zmat_arpack_zndrv3(ndims, nnzmax, nnz, acsr, jcsr, icsr, snnzmax, snnz, sacsr, sjcsr, sicsr, sigma, neval, nvecs, deval, zeigv)
         use para, only : dp, stdout, cpuid
         implicit none

! external arguments
! dimension of matrix A
         integer, intent(in) :: ndims

! maximum number of non-zero elements of matrix A
         integer, intent(in) :: nnzmax

! number of non-zero elements of matrix A
         integer, intent(inout) :: nnz

! sparse storage of matrix A
!> the input sparse storage format is coo 
         complex(dp), intent(inout) :: acsr(nnzmax)
         integer, intent(inout) :: jcsr(nnzmax)
         integer, intent(inout) :: icsr(nnzmax)

! optional for overlap matrix, when Orthogonal_Basis=T, this can be omited
!> the input sparse storage format is coo 
         integer, intent(in) :: snnzmax
         integer, intent(inout) :: snnz
         integer, intent(inout) :: sjcsr(snnzmax)
         integer, intent(inout) :: sicsr(snnzmax)
         complex(dp), intent(inout) :: sacsr(snnzmax)


! shift energy
         complex(dp), intent(in) :: sigma
! number of selected eigenvals
         integer, intent(in) :: neval

! number of Arnoldi vectors
         integer, intent(in) :: nvecs

! eigenvalues for selected "which"
         real(dp), intent(out) :: deval(neval)

! eigenvector for selected "which"
         complex(dp), intent(out) :: zeigv(ndims, nvecs)

! loop index over neval
         integer :: ival, jval

! axuiliary integer variables
         integer :: itmp, iter

! auxiliary real(dp) variables
         real(dp) :: dtmp

!
!     %--------------%
!     | Local Arrays |
!     %--------------%
!
         integer :: iparam(11), ipntr(14)
         complex(dp), allocatable :: ax(:)
         complex(dp), allocatable :: mx(:)
         complex(dp), allocatable :: d(:)
         complex(dp), allocatable :: v(:, :)
         complex(dp), allocatable :: workd(:)
         complex(dp), allocatable :: workev(:)
         complex(dp), allocatable :: resid(:)
         complex(dp), allocatable :: workl(:)
         logical   , allocatable :: select(:)

         integer, allocatable :: iwk(:)

         real(dp), allocatable :: rwork(:)
         real(dp), allocatable :: rd(:, :)

!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
         character  :: bmat*1, which*2
         integer    :: ido, n, nev, ncv, lworkl, info, j, &
            ierr, nconv, maxitr, ishfts, mode, ldv
         integer    :: istat
         integer    :: lwrkl, lwrkv, lwrkd
         real(dp)    :: tol
         logical    :: rvec
!
!     %-----------------------------%
!     | BLAS & LAPACK routines used |
!     %-----------------------------%
!
         real(dp), external :: dznrm2, dlapy2
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
!     %----------------------------------------------------%
!     | The number N is the dimension of the matrix.  A    |
!     | generalized eigenvalue problem is solved (BMAT =   |
!     | 'G').  NEV is the number of eigenvalues to be      |
!     | approximated.  The user can modify NEV, NCV, WHICH |
!     | to solve problems of different sizes, and to get   |
!     | different parts of the spectrum.  However, The     |
!     | following conditions must be satisfied:            |
!     |                    N <= MAXN,                      |
!     |                  NEV <= MAXNEV,                    |
!     |              NEV + 2 <= NCV <= MAXNCV              |
!     %----------------------------------------------------%
!
         n = ndims
         ldv = ndims
         nev = neval
         ncv = nvecs
         lwrkd = 3 * ndims; lwrkv = 3 * nvecs
         lwrkl = 3 * nvecs * nvecs  + 5 * nvecs


         ! allocate memory for arpack
         allocate(select(nvecs), stat=istat)
         allocate(ax(ndims), stat=istat)
         allocate(mx(ndims), stat=istat)
         allocate( d(neval), stat=istat)
        !allocate( v(ndims, nvecs), stat=istat)
         allocate(workd( lwrkd), stat=istat)
         allocate(workl( lwrkl), stat=istat)
         allocate(resid( ndims), stat=istat)
         allocate(workev(lwrkv), stat=istat)
         allocate(rwork(nvecs), stat=istat)
         allocate(rd(nvecs, 3), stat=istat)
         allocate(iwk(ndims+1))
         iwk= 0

         bmat  = 'G'
         which = 'SM'

         !> prepare hamiltonian
         !> transform coo format to csr format
         call ConvertCooToCsr(ndims, nnz, acsr, icsr, jcsr, iwk)
         call csr_sort_indices(ndims, nnz, icsr, jcsr, acsr)
         !> eleminate the same entries in the sparse matrix
         call csr_sum_duplicates(ndims, nnz, icsr, jcsr, acsr)
         call csr_sort_indices(ndims, nnz, icsr, jcsr, acsr)

         !> transform coo format to csr format
         call ConvertCooToCsr(ndims, snnz, sacsr, sicsr, sjcsr, iwk)
         call csr_sort_indices(ndims, snnz, sicsr, sjcsr, sacsr)
         !> eleminate the same entries in the sparse matrix
         call csr_sum_duplicates(ndims, snnz, sicsr, sjcsr, sacsr)
         call csr_sort_indices(ndims, snnz, sicsr, sjcsr, sacsr)


!
!     %---------------------------------------------------%
!     | The work array WORKL is used in ZNAUPD as         |
!     | workspace.  Its dimension LWORKL is set as        |
!     | illustrated below.  The parameter TOL determines  |
!     | the stopping criterion. If TOL<=0, machine        |
!     | precision is used.  The variable IDO is used for  |
!     | reverse communication, and is initially set to 0. |
!     | Setting INFO=0 indicates that a random vector is  |
!     | generated to start the ARNOLDI iteration.         |
!     %---------------------------------------------------%
!
         lworkl = 3 * ncv * ncv + 5 * ncv
         tol    = 1.0D-7
         ido    = 0
         info   = 0
!
!     %---------------------------------------------------%
!     | This program uses exact shift with respect to     |
!     | the current Hessenberg matrix (IPARAM(1) = 1).    |
!     | IPARAM(3) specifies the maximum number of Arnoldi |
!     | iterations allowed.  Mode 1 of ZNAUPD is used     |
!     | (IPARAM(7) = 1). All these options can be changed |
!     | by the user. For details see the documentation in |
!     | ZNAUPD.                                           |
!     %---------------------------------------------------%

!  Mode 2:  A*x = lambda*M*x, M hermitian positive definite
!           ===> OP = inv[M]*A  and  B = M.
!           ===> (If M can be factored see remark 3 below)!
!  ISHIFT = 1: exact shifts with respect to the current
!              Hessenberg matrix H.  This is equivalent to
!              restarting the iteration from the beginning
!              after updating the starting vector with a linear
!              combination of Ritz vectors associated with the
!              "wanted" eigenvalues.
         ishfts = 1
         maxitr = 50000
         mode   = 2

         iparam(1) = ishfts
         iparam(3) = maxitr
         iparam(7) = mode
         nconv=-1
!
!     %-------------------------------------------%
!     | M A I N   L O O P (Reverse communication) |
!     %-------------------------------------------%
!
         iter = 0
10       continue
!
!        %---------------------------------------------%
!        | Repeatedly call the routine ZNAUPD and take |
!        | actions indicated by parameter IDO until    |
!        | either convergence is indicated or maxitr   |
!        | has been exceeded.                          |
!        %---------------------------------------------%
!
#if defined (ARPACK)
         call znaupd ( ido, bmat, n, which, nev, tol, resid, ncv, &
            zeigv, ldv, iparam, ipntr, workd, workl, lworkl, rwork, info )
#else
         STOP "ERROR : Please install WannierTools with ARPACK since you are diagonalizing a large sparse matrix"
#endif

         iter = iter + 1
         if (mod(iter,100) .eq. 0 .and. cpuid==0) then
            write(stdout, '(A, I10)') '>>> Iteration with znaupd', iter
         endif

         if (ido .eq. -1 .or. ido .eq. 1) then
!
!           %----------------------------------------%
!           | Perform  y <--- OP*x = inv[S]*A*x      |
!           | The user should supply his/her own     |
!           | matrix vector routine and a linear     |
!           | system solver.  The matrix-vector      |
!           | subroutine should take workd(ipntr(1)) |
!           | as input, and the final result should  |
!           | be returned to workd(ipntr(2)).        |
!           %----------------------------------------%
!
 
            !> first perform A*x
            !call mkl_zcsrgemv('N', ndims, acsr, icsr, jcsr, workd(ipntr(1)), workd(ipntr(2)))
            call csrmv_z(ndims, nnz, acsr, icsr, jcsr, workd(ipntr(1)), workd(ipntr(2)))
            call zcopy (ndims, workd(ipntr(2)),  1, workd(ipntr(1)), 1)
            !> first perform inv[S]*A*x
            !call zmat_mkldss_zgesv(ndims, snnz, sacsr, sjcsr, sicsr, workd(ipntr(1)), workd(ipntr(2)))
            call sparse_solver(ndims, snnz, sacsr, sicsr, sjcsr, workd(ipntr(1)), workd(ipntr(2)))
!
!           %-----------------------------------------%
!           | L O O P   B A C K to call ZNAUPD again. |
!           %-----------------------------------------%
!
            go to 10
         else if ( ido .eq. 2) then
!
!           %-------------------------------------%
!           |        Perform  y <--- M*x          |
!           | The matrix vector multiplication    |
!           | routine should take workd(ipntr(1)) |
!           | as input and return the result to   |
!           | workd(ipntr(2)).                    |
!           %-------------------------------------%
!
            !call mkl_zcsrgemv('N', ndims, sacsr, sicsr, sjcsr, workd(ipntr(1)), workd(ipntr(2)))
            call csrmv_z(ndims, snnz, sacsr, sicsr, sjcsr, workd(ipntr(1)), workd(ipntr(2)))
!
!           %-----------------------------------------%
!           | L O O P   B A C K to call ZNAUPD  again. |
!           %-----------------------------------------%
!
            go to 10
!
         end if
!
!     %----------------------------------------%
!     | Either we have convergence or there is |
!     | an error.                              |
!     %----------------------------------------%
!
         if ( info .lt. 0 ) then
!
!        %--------------------------%
!        | Error message, check the |
!        | documentation in ZNAUPD  |
!        %--------------------------%

            if (cpuid==0) write(stdout, *) ' '
            if (cpuid==0) write(stdout, *) ' Error with _naupd, info = ', info
            if (cpuid==0) write(stdout, *) ' Check the documentation of _naupd'
            if (cpuid==0) write(stdout, *) ' '
            stop

         else
!
!        %-------------------------------------------%
!        | No fatal errors occurred.                 |
!        | Post-Process using ZNEUPD.                |
!        |                                           |
!        | Computed eigenvalues may be extracted.    |
!        |                                           |
!        | Eigenvectors may also be computed now if  |
!        | desired.  (indicated by rvec = .true.)    |
!        %-------------------------------------------%
!
            rvec = .false.
            if (LandauLevel_wavefunction_calc.or.SlabBand_calc) rvec = .true.

#if defined (ARPACK)
            call zneupd (rvec, 'A', select, d, zeigv, ldv, sigma, &
               workev, bmat, n, which, nev, tol, resid, ncv,&
               zeigv, ldv, iparam, ipntr, workd, workl, lworkl, &
               rwork, ierr)
#else
         STOP "ERROR : Please install WannierTools with ARPACK since you are diagonalizing a large sparse matrix"
#endif
!
!        %----------------------------------------------%
!        | Eigenvalues are returned in the one          |
!        | dimensional array D.  The corresponding      |
!        | eigenvectors are returned in the first NCONV |
!        | (=IPARAM(5)) columns of the two dimensional  |
!        | array V if requested.  Otherwise, an         |
!        | orthogonal basis for the invariant subspace  |
!        | corresponding to the eigenvalues in D is     |
!        | returned in V.                               |
!        %----------------------------------------------%
!
            if ( ierr .ne. 0) then
!
!           %------------------------------------%
!           | Error condition:                   |
!           | Check the documentation of ZNEUPD. |
!           %------------------------------------%
!
               if (cpuid==0) write(stdout, *) ' '
               if (cpuid==0) write(stdout, *) ' Error with _neupd, info = ', ierr
               if (cpuid==0) write(stdout, *) ' Check the documentation of _neupd. '
               if (cpuid==0) write(stdout, *) ' '

            else

               nconv = iparam(5)
               do 20 j=1, nconv
!
!               %---------------------------%
!               | Compute the residual norm |
!               |                           |
!               |   ||  A*x - lambda*S*x ||   |
!               |                           |
!               | for the NCONV accurately  |
!               | computed eigenvalues and  |
!               | eigenvectors.  (iparam(5) |
!               | indicates how many are    |
!               | accurate to the requested |
!               | tolerance)                |
!               %---------------------------%
!

                  !call mkl_zcsrgemv('N', ndims, sacsr, sicsr, sjcsr, zeigV(1,j), mx)
                  call csrmv_z(ndims, snnz, sacsr, sicsr, sjcsr, zeigV(1,j), mx)
                  !call mkl_zcsrgemv('N', ndims, acsr, icsr, jcsr, zeigV(1,j), ax)
                  call csrmv_z(ndims, nnz, acsr, icsr, jcsr, zeigV(1,j), ax)
                  call zaxpy(n, -d(j), mx, 1, ax, 1)
                  rd(j,1) = dble(d(j))
                  rd(j,2) = dimag(d(j))
                  rd(j,3) = dznrm2(n, ax, 1)
                  rd(j,3) = rd(j,3) / dlapy2(rd(j,1),rd(j,2))
20             continue
!
!            %-----------------------------%
!            | Display computed residuals. |
!            %-----------------------------%
!
!              call dmout(6, nconv, 3, rd, nvecs, -6, &
!                 'Ritz values (Real, Imag) and relative residuals')
            end if
!
!        %-------------------------------------------%
!        | Print additional convergence information. |
!        %-------------------------------------------%
!
            if ( info .eq. 1) then
               if (cpuid==0) write(stdout, *) ' '
               if (cpuid==0) write(stdout, *) ' Maximum number of iterations reached.'
               if (cpuid==0) write(stdout, *) ' '
            else if ( info .eq. 3) then
               if (cpuid==0) write(stdout, *) ' '
               if (cpuid==0) write(stdout, *) ' No shifts could be applied during implicit &
               &Arnoldi update, try increasing NCV.'
               if (cpuid==0) write(stdout, *) ' '
            end if

            !print *, ' '
            if (cpuid==0) write(stdout, *) '_NDRV3'
            if (cpuid==0) write(stdout, *) '====== '
            if (cpuid==0) write(stdout, *) ' '
            if (cpuid==0) write(stdout, *) ' Size of the matrix is ', n
            if (cpuid==0) write(stdout, *) ' The number of Ritz values requested is ', nev
            if (cpuid==0) write(stdout, *) ' The number of Arnoldi vectors generated', &
               ' (NCV) is ', ncv
            if (cpuid==0) write(stdout, *) ' What portion of the spectrum: ', which
            if (cpuid==0) write(stdout, *) ' The number of converged Ritz values is ', &
               nconv
            if (cpuid==0) write(stdout, *) ' The number of Implicit Arnoldi update', &
               ' iterations taken is ', iparam(3)
            if (cpuid==0) write(stdout, *) ' The number of OP*x is ', iparam(9)
            if (cpuid==0) write(stdout, *) ' The convergence criterion is ', tol
            if (cpuid==0) write(stdout, *) ' '

         end if
!
!     %---------------------------%
!     | Done with program zndrv3. |
!     %---------------------------%
!

         do ival=1,neval
            deval(ival) = real(d(ival))
           !zeigv(1:ndims, ival) = V(1:ndims, ival)
         enddo ! over ival={1,neval} loop

!---------------------------------------------------------------------!
! use selection sort to minimize swaps of eigenvectors, ref: dsteqr.f !
!---------------------------------------------------------------------!
         do ival=1,neval-1
            itmp = ival; dtmp = deval(itmp)

            do jval=ival+1,neval
               if ((deval(jval)+1.0D-12) .lt. dtmp) then
                  itmp = jval; dtmp = deval(itmp)
               endif
            enddo ! over jval={ival+1,neval} loop

            if (itmp .ne. ival) then
               deval(itmp) = deval(ival); deval(ival) = dtmp
               call zswap(ndims, zeigv(1, ival), 1, zeigv(1,itmp), 1)
            endif
         enddo ! over ival={1,neval-1} loop

         return
      end subroutine zmat_arpack_zndrv3

      subroutine zmat_arpack_zndrv1(ndims, nnzmax, nnz, acsr, jcsr, icsr, sigma, neval, nvecs, deval, zeigv, ritzvec)
         use para, only : dp, stdout, cpuid
         implicit none

! external arguments
! dimension of matrix A
         integer, intent(in) :: ndims

! maximum number of non-zero elements of matrix A
         integer, intent(in) :: nnzmax

! number of non-zero elements of matrix A
         integer, intent(inout) :: nnz

! sparse storage of matrix A
!> the input sparse storage format is coo 
         complex(dp), intent(inout) :: acsr(nnzmax)
         integer, intent(inout) :: jcsr(nnzmax)
         integer, intent(inout) :: icsr(nnzmax)

! shift energy
         complex(dp), intent(in) :: sigma
! number of selected eigenvals
         integer, intent(in) :: neval

! number of Arnoldi vectors
         integer, intent(in) :: nvecs

! eigenvalues for selected "which"
         real(dp), intent(out) :: deval(neval)

! eigenvector for selected "which"
         complex(dp), intent(out) :: zeigv(ndims, nvecs)

!> calculate eigenvector or not
         logical, intent(in) :: ritzvec

! loop index over neval
         integer :: ival, jval

! axuiliary integer variables
         integer :: itmp, iter

! auxiliary real(dp) variables
         real(dp) :: dtmp

!
!     %--------------%
!     | Local Arrays |
!     %--------------%
!
         integer :: iparam(11), ipntr(14)
         complex(dp), allocatable :: ax(:)
         complex(dp), allocatable :: d(:)
         complex(dp), allocatable :: v(:, :)
         complex(dp), allocatable :: workd(:)
         complex(dp), allocatable :: workev(:)
         complex(dp), allocatable :: resid(:)
         complex(dp), allocatable :: workl(:)
         logical   , allocatable :: select(:)

         integer, allocatable :: iwk(:)

         real(dp), allocatable :: rwork(:)
         real(dp), allocatable :: rd(:, :)

!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
         character  :: bmat*1, which*2
         integer    :: ido, n, nev, ncv, lworkl, info, i, j, &
            ierr, nconv, maxitr, ishfts, mode, ldv
         integer    :: istat 
         integer    :: lwrkl, lwrkv, lwrkd
         real(dp)    :: tol
!
!     %-----------------------------%
!     | BLAS & LAPACK routines used |
!     %-----------------------------%
!
         real(dp), external :: dznrm2, dlapy2
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
!     %--------------------------------------------------%
!     | The number NX is the number of interior points   |
!     | in the discretization of the 2-dimensional       |
!     | convection-diffusion operator on the unit        |
!     | square with zero Dirichlet boundary condition.   |
!     | The number N(=NX*NX) is the dimension of the     |
!     | matrix.  A standard eigenvalue problem is        |
!     | solved (BMAT = 'I').  NEV is the number of       |
!     | eigenvalues to be approximated.  The user can    |
!     | modify NX, NEV, NCV, WHICH to solve problems of  |
!     | different sizes, and to get different parts of   |
!     | the spectrum.  However, The following            |
!     | conditions must be satisfied:                    |
!     |                   N <= MAXN                      |
!     |                 NEV <= MAXNEV                    |
!     |           NEV + 2 <= NCV <= MAXNCV               |
!     %--------------------------------------------------%
!
         n = ndims
         ldv = ndims
         nev = neval
         ncv = nvecs
         lwrkd = 3 * ndims; lwrkv = 3 * nvecs
         lwrkl = 3 * nvecs * nvecs  + 5 * nvecs


         ! allocate memory for arpack
         allocate(select(nvecs), stat=istat)
         allocate(ax(ndims), stat=istat)
         allocate( d(neval), stat=istat)
        !allocate( v(ndims, nvecs), stat=istat)
         allocate(workd( lwrkd), stat=istat)
         allocate(workl( lwrkl), stat=istat)
         allocate(resid( ndims), stat=istat)
         allocate(workev(lwrkv), stat=istat)
         allocate(rwork(nvecs), stat=istat)
         allocate(rd(nvecs, 3), stat=istat)
         allocate(iwk(ndims+1))
         iwk= 0

         bmat  = 'I'
         which = 'SM'

         !> added in the diagonal part, shift the spectrum by sigma
         do i=1, ndims
            j=i+nnz
            icsr(j)= i
            jcsr(j)= i
            acsr(j)= -sigma
         enddo
         nnz= nnz+ ndims

         !> prepare hamiltonian
         !> transform coo format to csr format
         call ConvertCooToCsr(ndims, nnz, acsr, icsr, jcsr, iwk)

         call csr_sort_indices(ndims, nnz, icsr, jcsr, acsr)
         !> eleminate the same entries in the sparse matrix
         call csr_sum_duplicates(ndims, nnz, icsr, jcsr, acsr)
         call csr_sort_indices(ndims, nnz, icsr, jcsr, acsr)

!
!     %---------------------------------------------------%
!     | The work array WORKL is used in ZNAUPD as         |
!     | workspace.  Its dimension LWORKL is set as        |
!     | illustrated below.  The parameter TOL determines  |
!     | the stopping criterion. If TOL<=0, machine        |
!     | precision is used.  The variable IDO is used for  |
!     | reverse communication, and is initially set to 0. |
!     | Setting INFO=0 indicates that a random vector is  |
!     | generated to start the ARNOLDI iteration.         |
!     %---------------------------------------------------%
!
         lworkl = 3 * ncv * ncv + 5 * ncv
         tol    = 1.0D-7
         ido    = 0
         info   = 0
!
!     %---------------------------------------------------%
!     | This program uses exact shift with respect to     |
!     | the current Hessenberg matrix (IPARAM(1) = 1).    |
!     | IPARAM(3) specifies the maximum number of Arnoldi |
!     | iterations allowed.  Mode 1 of ZNAUPD is used     |
!     | (IPARAM(7) = 1). All these options can be changed |
!     | by the user. For details see the documentation in |
!     | ZNAUPD.                                           |
!     %---------------------------------------------------%
!
         ishfts = 1
         maxitr = 50000
         mode   = 1

         iparam(1) = ishfts
         iparam(3) = maxitr
         iparam(7) = mode
!
!     %-------------------------------------------%
!     | M A I N   L O O P (Reverse communication) |
!     %-------------------------------------------%
!
         iter = 0
10       continue
!
!        %---------------------------------------------%
!        | Repeatedly call the routine ZNAUPD and take |
!        | actions indicated by parameter IDO until    |
!        | either convergence is indicated or maxitr   |
!        | has been exceeded.                          |
!        %---------------------------------------------%
!
#if defined (ARPACK)
         call znaupd ( ido, bmat, n, which, nev, tol, resid, ncv, &
            zeigv, ldv, iparam, ipntr, workd, workl, lworkl, rwork, info )
#else
         STOP "ERROR : Please install WannierTools with ARPACK since you are diagonalizing a large sparse matrix"
#endif

         iter = iter + 1
         if (mod(iter,100) .eq. 0 .and. cpuid==0) then
            write(stdout, '(A, I10)') '>>> Iteration with znaupd', iter
         endif

         if (ido .eq. -1 .or. ido .eq. 1) then
!
!           %-------------------------------------------%
!           | Perform matrix vector multiplication      |
!           |   y <--- OP*x =  A*x                      |
!           | The user should supply his/her own        |
!           | matrix vector multiplication routine here |
!           | that takes workd(ipntr(1)) as the input   |
!           | vector, and return the matrix vector      |
!           | product to workd(ipntr(2)).               |
!           %-------------------------------------------%
!
            !call mkl_zcsrgemv('N', ndims, acsr, icsr, jcsr, workd(ipntr(1)), workd(ipntr(2)))
            call csrmv_z(ndims, nnz, acsr, icsr, jcsr, workd(ipntr(1)), workd(ipntr(2)))
!
!           %-----------------------------------------%
!           | L O O P   B A C K to call ZNAUPD again. |
!           %-----------------------------------------%
!
            go to 10

         end if
!
!     %----------------------------------------%
!     | Either we have convergence or there is |
!     | an error.                              |
!     %----------------------------------------%
!
         if ( info .lt. 0 ) then
!
!        %--------------------------%
!        | Error message, check the |
!        | documentation in ZNAUPD  |
!        %--------------------------%

            if (cpuid==0) write(stdout, *) ' '
            if (cpuid==0) write(stdout, *) ' Error with _naupd, info = ', info
            if (cpuid==0) write(stdout, *) ' Check the documentation of _naupd'
            if (cpuid==0) write(stdout, *) ' '
            stop

         else
!
!        %-------------------------------------------%
!        | No fatal errors occurred.                 |
!        | Post-Process using ZNEUPD.                |
!        |                                           |
!        | Computed eigenvalues may be extracted.    |
!        |                                           |
!        | Eigenvectors may also be computed now if  |
!        | desired.  (indicated by rvec = .true.)    |
!        %-------------------------------------------%
!
#if defined (ARPACK)
            call zneupd (ritzvec, 'A', select, d, zeigv, ldv, sigma, &
               workev, bmat, n, which, nev, tol, resid, ncv,&
               zeigv, ldv, iparam, ipntr, workd, workl, lworkl, &
               rwork, ierr)
#else
         STOP "ERROR : Please install WannierTools with ARPACK"
#endif
!
!        %----------------------------------------------%
!        | Eigenvalues are returned in the one          |
!        | dimensional array D.  The corresponding      |
!        | eigenvectors are returned in the first NCONV |
!        | (=IPARAM(5)) columns of the two dimensional  |
!        | array V if requested.  Otherwise, an         |
!        | orthogonal basis for the invariant subspace  |
!        | corresponding to the eigenvalues in D is     |
!        | returned in V.                               |
!        %----------------------------------------------%
!
            if ( ierr .ne. 0) then
!
!           %------------------------------------%
!           | Error condition:                   |
!           | Check the documentation of ZNEUPD. |
!           %------------------------------------%
!
               if (cpuid==0) write(stdout, *) ' '
               if (cpuid==0) write(stdout, *) ' Error with _neupd, info = ', ierr
               if (cpuid==0) write(stdout, *) ' Check the documentation of _neupd. '
               if (cpuid==0) write(stdout, *) ' '

            else

               nconv = iparam(5)
               do 20 j=1, nconv
!
!               %---------------------------%
!               | Compute the residual norm |
!               |                           |
!               |   ||  A*x - lambda*x ||   |
!               |                           |
!               | for the NCONV accurately  |
!               | computed eigenvalues and  |
!               | eigenvectors.  (iparam(5) |
!               | indicates how many are    |
!               | accurate to the requested |
!               | tolerance)                |
!               %---------------------------%
!
                  !call mkl_zcsrgemv('N', ndims, acsr, icsr, jcsr, zeigV(1,j), ax)
                  call csrmv_z(ndims, nnz, acsr, icsr, jcsr, zeigV(1,j), ax)
                  call zaxpy(n, -d(j), zeigv(1,j), 1, ax, 1)
                  rd(j,1) = dble(d(j))
                  rd(j,2) = dimag(d(j))
                  rd(j,3) = dznrm2(n, ax, 1)
                  rd(j,3) = rd(j,3) / dlapy2(rd(j,1),rd(j,2))
20             continue
!
!            %-----------------------------%
!            | Display computed residuals. |
!            %-----------------------------%
!
!              call dmout(6, nconv, 3, rd, nvecs, -6, &
!                 'Ritz values (Real, Imag) and relative residuals')
            end if
!
!        %-------------------------------------------%
!        | Print additional convergence information. |
!        %-------------------------------------------%
!
            if ( info .eq. 1) then
               if (cpuid==0) write(stdout, *) ' '
               if (cpuid==0) write(stdout, *) ' Maximum number of iterations reached. try increasing NCV'
               if (cpuid==0) write(stdout, *) ' '
            else if ( info .eq. 3) then
               if (cpuid==0) write(stdout, *) ' '
               if (cpuid==0) write(stdout, *) ' No shifts could be applied during implicit &
               &Arnoldi update, try increasing NCV.'
               if (cpuid==0) write(stdout, *) ' '
            end if

            !print *, ' '
            if (cpuid==0) write(stdout, *) '_NDRV1'
            if (cpuid==0) write(stdout, *) '====== '
            if (cpuid==0) write(stdout, *) ' '
            if (cpuid==0) write(stdout, *) ' Size of the matrix is ', n
            if (cpuid==0) write(stdout, *) ' The number of Ritz values requested is ', nev
            if (cpuid==0) write(stdout, *) ' The number of Arnoldi vectors generated', &
               ' (NCV) is ', ncv
            if (cpuid==0) write(stdout, *) ' What portion of the spectrum: ', which
            if (cpuid==0) write(stdout, *) ' The number of converged Ritz values is ', &
               nconv
            if (cpuid==0) write(stdout, *) ' The number of Implicit Arnoldi update', &
               ' iterations taken is ', iparam(3)
            if (cpuid==0) write(stdout, *) ' The number of OP*x is ', iparam(9)
            if (cpuid==0) write(stdout, *) ' The convergence criterion is ', tol
            if (cpuid==0) write(stdout, *) ' '

         end if
!
!     %---------------------------%
!     | Done with program zndrv1. |
!     %---------------------------%
!
         !shift back the spectrum by sigma
         do ival=1,neval
            deval(ival) = real(d(ival)+sigma)
           !zeigv(1:ndims, ival) = V(1:ndims, ival)
         enddo ! over ival={1,neval} loop

!---------------------------------------------------------------------!
! use selection sort to minimize swaps of eigenvectors, ref: dsteqr.f !
!---------------------------------------------------------------------!
         do ival=1,neval-1
            itmp = ival; dtmp = deval(itmp)

            do jval=ival+1,neval
               if ((deval(jval)+1.0D-12) .lt. dtmp) then
                  itmp = jval; dtmp = deval(itmp)
               endif
            enddo ! over jval={ival+1,neval} loop

            if (itmp .ne. ival) then
               deval(itmp) = deval(ival); deval(ival) = dtmp
               call zswap(ndims, zeigv(1, ival), 1, zeigv(1,itmp), 1)
            endif
         enddo ! over ival={1,neval-1} loop

         return
      end subroutine zmat_arpack_zndrv1


      subroutine zmat_arpack_zndrv2(ndims, nnzmax, nnz, acsr, jcsr, icsr, sigma,neval, nvecs, deval, zeigv, ritzvec)
         use para, only : dp, stdout, cpuid, LandauLevel_wavefunction_calc, SlabBand_calc
         implicit none

         ! external arguments
         ! dimension of matrix A
         integer, intent(in) :: ndims

         ! maximum number of non-zero elements of matrix A
         integer, intent(in) :: nnzmax

         ! number of non-zero elements of matrix A
         integer, intent(inout) :: nnz

         !> input: coo format storage of matrix A
         !> output: csr format storage of matrix A-sigma*I
         complex(dp), intent(inout) :: acsr(nnzmax)
         integer, intent(inout) :: jcsr(nnzmax)
         integer, intent(inout) :: icsr(nnzmax)

         !> calculate eigenvalue close to sigma
         complex(dp), intent(in) :: sigma

         ! number of selected eigenvals
         integer, intent(in) :: neval

         ! number of Arnoldi vectors
         integer, intent(in) :: nvecs

         !> calculate eigenvector or not
         logical, intent(in) :: ritzvec

         ! eigenvalues for selected "which"
         real(dp), intent(out) :: deval(neval)

         ! eigenvector for selected "which"
         ! dimes= (ndims, neval)
         complex(dp), intent(out) :: zeigv(ndims, nvecs)

         ! loop index over neval
         integer :: ival, jval

         ! axuiliary integer variables
         integer :: itmp, iter, i

         ! auxiliary real(dp) variables
         real(dp) :: dtmp
         real(dp) :: memory

         !
         !     %------------------------%
         !     | Local Arrays for ARPACK|
         !     %------------------------%
         !
         integer :: iparam(11), ipntr(14)
         complex(dp), allocatable :: ax(:)
         complex(dp), allocatable :: d(:)
         complex(dp), allocatable :: workd(:)
         complex(dp), allocatable :: workev(:)
         complex(dp), allocatable :: resid(:)
         complex(dp), allocatable :: workl(:)
         complex(dp), allocatable :: v(:,:)
         logical   , allocatable :: select(:)

         integer, allocatable :: iwk(:)

         real(dp), allocatable :: rwork(:)
         real(dp), allocatable :: rd(:, :)

         !
         !     %---------------%
         !     | Local Scalars |
         !     %---------------%
         !
         character         bmat*1, which*2
         integer           ido, n, nev, ncv, lworkl, info, j, ierr, &
            nconv, maxitr, ishfts, mode, ldv
         integer    :: istat
         integer    :: lwrkl, lwrkv, lwrkd
         real(dp)    :: tol

         !
         !     %-----------------------------%
         !     | BLAS & LAPACK routines used |
         !     %-----------------------------%
         !
         real(dp), external :: dznrm2, dlapy2
         !
         !     %-----------------------%
         !     | Executable statements |
         !     %-----------------------%
         !
         !     %--------------------------------------------------%
         !     | The number N is the dimension of the matrix.  A  |
         !     | standard eigenvalue problem is solved (BMAT =    |
         !     | 'I').  NEV is the number of eigenvalues (closest |
         !     | to the shift SIGMA) to be approximated.  Since   |
         !     | the shift-invert mode is used, WHICH is set to   |
         !     | 'LM'.  The user can modify NEV, NCV, SIGMA to    |
         !     | solve problems of different sizes, and to get    |
         !     | different parts of the spectrum.  However, The   |
         !     | following conditions must be satisfied:          |
         !     |                 N <= MAXN,                       |
         !     |               NEV <= MAXNEV,                     |
         !     |           NEV + 2 <= NCV <= MAXNCV               |
         !     %--------------------------------------------------%
         !
         n = ndims
         ldv = ndims
         nev = neval
         ncv = nvecs
         if(ncv>n-1) ncv=n
         lwrkd = 3 * ndims; lwrkv = 3 * nvecs
         lwrkl = 3 * nvecs * nvecs  + 5 * nvecs
         ! allocate memory for arpack
         allocate(select(nvecs), stat=istat)
         allocate(ax(ndims), stat=istat)
         allocate( d(nvecs), stat=istat)
         !allocate( v(ndims, nvecs), stat=istat)
         allocate(workd( lwrkd), stat=istat)
         allocate(workl( lwrkl), stat=istat)
         allocate(resid( ndims), stat=istat)
         allocate(workev(lwrkv), stat=istat)
         allocate(rwork(nvecs), stat=istat)
         allocate(rd(nvecs, 3), stat=istat)
         allocate(iwk(ndims+1))

         memory= ndims/1024d0/1024d0*16d0*3+(lwrkl+lwrkd+lwrkv)/1024d0/1024d0*16d0
         if (cpuid.eq.0) then
            write(stdout, '(a, f16.1, a)')' > memory cost in zmat_arpack_zndrv2 is : ', memory, 'MB'
            write(stdout, '(a, I10, a)')' > We are calculating ', neval, ' eigenvalues'
            write(stdout, '(a, I10, a)')' > with ', ncv, ' Arnoldi Lanczos vectors.'
            write(stdout, '(a, I16)')' > Dimension of the matrix is: ', ndims
         endif


         
         iwk= 0
         bmat  = 'I'
         which = 'LM'

         !
         !     %----------------------------------------------------%
         !     | Construct C = A - SIGMA*I, factor C in complex     |
         !     | arithmetic (using LAPACK subroutine zgttrf). The   |
         !     | matrix A is chosen to be the tridiagonal matrix    |
         !     | derived from standard central difference of the    |
         !     | 1-d convection diffusion operator - u" + rho*u' on |
         !     | the interval [0, 1] with zero Dirichlet boundary   |
         !     | condition.                                         |
         !     %----------------------------------------------------%
         !

         !> added in the diagonal part
         do i=1, ndims
            j=i+nnz
            icsr(j)= i
            jcsr(j)= i
            acsr(j)= -sigma
         enddo
         nnz= nnz+ ndims


         !> transform coo format to csr format
         call ConvertCooToCsr(ndims, nnz, acsr, icsr, jcsr, iwk)
         call csr_sort_indices(ndims, nnz, icsr, jcsr, acsr)

         !> sum up the same entries in the sparse matrix
         call csr_sum_duplicates(ndims, nnz, icsr, jcsr, acsr)
         call csr_sort_indices(ndims, nnz, icsr, jcsr, acsr)

         !
         !     %-----------------------------------------------------%
         !     | The work array WORKL is used in ZNAUPD as           |
         !     | workspace.  Its dimension LWORKL is set as          |
         !     | illustrated below.  The parameter TOL determines    |
         !     | the stopping criterion. If TOL<=0, machine          |
         !     | precision is used.  The variable IDO is used for    |
         !     | reverse communication, and is initially set to 0.   |
         !     | Setting INFO=0 indicates that a random vector is    |
         !     | generated in ZNAUPD to start the Arnoldi iteration. |
         !     %-----------------------------------------------------%
         !
         lworkl = 3 * ncv * ncv + 5 * ncv
         tol    = 1.0D-7
         ido    = 0
         info   = 0
         !
         !     %---------------------------------------------------%
         !     | This program uses exact shifts with respect to    |
         !     | the current Hessenberg matrix (IPARAM(1) = 1).    |
         !     | IPARAM(3) specifies the maximum number of Arnoldi |
         !     | iterations allowed. Mode 3 of ZNAUPD is used      |
         !     | (IPARAM(7) = 3).  All these options can be        |
         !     | changed by the user. For details see the          |
         !     | documentation in ZNAUPD.                          |
         !     %---------------------------------------------------%
         !
         ishfts = 1
         maxitr = 50000
         mode   = 3

         iparam(1) = ishfts
         iparam(3) = maxitr
         iparam(7) = mode


         !
         !     %-------------------------------------------%
         !     | M A I N   L O O P (Reverse communication) |
         !     %-------------------------------------------%
         !
         iter = 0
20       continue
         !
         !        %---------------------------------------------%
         !        | Repeatedly call the routine ZNAUPD and take |
         !        | actions indicated by parameter IDO until    |
         !        | either convergence is indicated or maxitr   |
         !        | has been exceeded.                          |
         !        %---------------------------------------------%
         !
         !    write(*,*) 'ncv,nev',ncv,nev,n

         
#if defined (ARPACK)
         call znaupd ( ido, bmat, n, which, nev, tol, resid, ncv, &
            zeigv, ldv, iparam, ipntr, workd, workl, lworkl, rwork,info )
#else
         STOP "ERROR : Please install WannierTools with ARPACK since you are diagonalizing a large sparse matrix"
#endif

         iter = iter + 1
         if (mod(iter,10) .eq. 0 .and. cpuid==0 ) then
            write(stdout, '(A, I8)') '>>> Iteration with znaupd', iter
         endif ! back if (mod(iter,1) .eq. 0) block

         if (ido .eq. -1 .or. ido .eq. 1 ) then
            !
            !           %-------------------------------------------%
            !           | Perform  y <--- OP*x = inv[A-SIGMA*I]*x   |
            !           | The user should supply his/her own linear |
            !           | system solver here that takes             |
            !           | workd(ipntr(1)) as the input, and returns |
            !           | the result to workd(ipntr(2)).            |
            !           %-------------------------------------------%
            !
            call sparse_solver(ndims, nnz, acsr, icsr, jcsr, workd(ipntr(1)), workd(ipntr(2)))
            !
            !           %-----------------------------------------%
            !           | L O O P   B A C K to call ZNAUPD again. |
            !           %-----------------------------------------%
            !
            go to 20

         end if
         !
         !     %-----------------------------------------%
         !     | Either we have convergence, or there is |
         !     | an error.                               |
         !     %-----------------------------------------%
         !
         if ( info .lt. 0 ) then
            !
            !        %--------------------------%
            !        | Error message, check the |
            !        | documentation in ZNAUPD  |
            !        %--------------------------%
            !
            if (cpuid==0) write(stdout, *) ' '
            if (cpuid==0) write(stdout, *) ' Error with _naupd, info = ',info
            if (cpuid==0) write(stdout, *) ' Check the documentation in _naupd.'
            if (cpuid==0) write(stdout, *) ' '

         else
            !
            !        %-------------------------------------------%
            !        | No fatal errors occurred.                 |
            !        | Post-Process using ZNEUPD.                |
            !        |                                           |
            !        | Computed eigenvalues may be extracted.    |
            !        |                                           |
            !        | Eigenvectors may also be computed now if  |
            !        | desired.  (indicated by ritzvec = .true.)    |
            !        %-------------------------------------------%
            !

#if defined (ARPACK)
            call zneupd (ritzvec, 'A', select, d, zeigv, ldv, sigma, &
               workev, bmat, n, which, nev, tol, resid, ncv,&
               zeigv, ldv, iparam, ipntr, workd, workl, lworkl, &
               rwork, ierr)
#else
         STOP "ERROR : Please install WannierTools with ARPACK since you are diagonalizing a large sparse matrix"
#endif
            !
            !        %----------------------------------------------%
            !        | Eigenvalues are returned in the one          |
            !        | dimensional array D.  The corresponding      |
            !        | eigenvectors are returned in the first NCONV |
            !        | (=IPARAM(5)) columns of the two dimensional  |
            !        | array zeigV if requested.  Otherwise, an         |
            !        | orthogonal basis for the invariant subspace  |
            !        | corresponding to the eigenvalues in D is     |
            !        | returned in zeigV.                               |
            !        %----------------------------------------------%
            !
            if ( ierr .ne. 0) then
               !
               !           %------------------------------------%
               !           | Error condition:                   |
               !           | Check the documentation of ZNEUPD. |
               !           %------------------------------------%
               !
               if (cpuid==0) write(stdout, *) ' '
               if (cpuid==0) write(stdout, *) ' Error with _neupd, info = ', ierr
               if (cpuid==0) write(stdout, *) ' Check the documentation of _neupd. '
               if (cpuid==0) write(stdout, *) ' '
            else

               nconv = iparam(5)
               do 60 j=1, nconv
                  !
                  !              %---------------------------%
                  !              | Compute the residual norm |
                  !              |                           |
                  !              |   ||  A*x - lambda*x ||   |
                  !              |                           |
                  !              | for the NCONV accurately  |
                  !              | computed eigenvalues and  |
                  !              | eigenvectors.  (iparam(5) |
                  !              | indicates how many are    |
                  !              | accurate to the requested |
                  !              | tolerance)                |
                  !              %---------------------------%
                  !
                  !$ call av(n, zeigv(1,j), ax)
                  !call mkl_zcsrgemv('N', ndims, acsr, icsr, jcsr, zeigv(1,j), ax)
                  call csrmv_z(ndims, nnz, acsr, icsr, jcsr, zeigV(1,j), ax)
                  call zaxpy(n, -d(j), zeigv(1,j), 1, ax, 1)
                  rd(j,1) = dble(d(j))
                  rd(j,2) = dimag(d(j))
                  rd(j,3) = dznrm2(n, ax, 1)
                  rd(j,3) = rd(j,3) / dlapy2(rd(j,1),rd(j,2))
60             continue
               !
               !           %-----------------------------%
               !           | Display computed residuals. |
               !           %-----------------------------%
               !
               !                call dmout(6, nconv, 3, rd, ncv, -6, &
               !                    'Ritz values (Real, Imag) and relative residuals')

            end if
            !
            !        %-------------------------------------------%
            !        | if (cpuid==0) write(stdout, *)ditional convergence information. |
            !        %-------------------------------------------%
            !
            if ( info .eq. 1) then
               if (cpuid==0) write(stdout, *) ' '
               if (cpuid==0) write(stdout, *) ' Maximum number of iterations reached.'
               if (cpuid==0) write(stdout, *) ' '
            else if ( info .eq. 3) then
               if (cpuid==0) write(stdout, *) ' '
               if (cpuid==0) write(stdout, *) ' No shifts could be applied during implicit &
               & Arnoldi update, try increasing NCV.'
               if (cpuid==0) write(stdout, *) ' '
            end if
            !
            !        print *, ' '
            !        print *, '_NDRV2 '
            !        print *, '====== '
            !        print *, ' '
            !        print *, ' Size of the matrix is ', n
            !        print *, ' The number of Ritz values requested is ', nev
            !        print *, ' The number of Arnoldi vectors generated', &
            !            ' (NCV) is ', ncv
            !        print *, ' What portion of the spectrum: ', which
            !        print *, ' The number of converged Ritz values is ', &
            !            nconv
            !        print *, ' The number of Implicit Arnoldi update',   &
            !            ' iterations taken is ', iparam(3)
            !        print *, ' The number of OP*x is ', iparam(9)
            !        print *, ' The convergence criterion is ', tol
            !        print *, ' '

         end if
         !
         !     %---------------------------%
         !     | Done with program zndrv2. |
         !     %---------------------------%
         !

         do ival=1,neval
            deval(ival) = real(d(ival))
           !zeigv(1:ndims, ival) = zeigV(1:ndims, ival)
         enddo ! over ival={1,neval} loop

         !---------------------------------------------------------------------!
         ! use selection sort to minimize swaps of eigenvectors, ref: dsteqr.f !
         !---------------------------------------------------------------------!
         do ival=1,neval-1
            itmp = ival; dtmp = deval(itmp)

            do jval=ival+1,neval
               if ((deval(jval)+1.0D-12) .lt. dtmp) then
                  itmp = jval; dtmp = deval(itmp)
               endif
            enddo ! over jval={ival+1,neval} loop

            if (itmp .ne. ival) then
               deval(itmp) = deval(ival); deval(ival) = dtmp
               call zswap(ndims, zeigv(1, ival), 1, zeigv(1,itmp), 1)
            endif
         enddo ! over ival={1,neval-1} loop
!    write(*,*) 'well done zndrv2'

         deallocate(select)
         deallocate(ax)
         deallocate( d)
        !deallocate( v)
         deallocate(workd)
         deallocate(workl)
         deallocate(resid)
         deallocate(workev)
         deallocate(rwork)
         deallocate(rd)
         deallocate(iwk)

         return
      end subroutine zmat_arpack_zndrv2

   end module sparse


