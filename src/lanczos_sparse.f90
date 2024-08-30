 !> initialize lanczos vector with a constant vector
subroutine InitLanczosVector(vector)
   use prec
   use sparse
   implicit none
   real(dp) :: norm
   type (WTParVec), pointer :: vector

   call WTParVecSetSeqValue(vector)
   norm= WTParVecNorm(vector)
   vector= 1d0/norm*vector

   return
endsubroutine InitLanczosVector

 !> Lanczos algorithm
 !> NumLczVectors : number of Lanczos vectors
 !> Mdim : The dimension of the Ham when it becomes a dense matrix
 !> input: initial vector; Matrix stored in PARCSR format
 !> output : Alpha(n), Beta(n)
subroutine lanczos_seqsparse_cpu_z(NumLczVectors, NumLczVectors_out, Mdim, nnz, icsr, jcsr, acsr, InitialVector, Alpha, Betan)
   use prec
   use wmpi
   use sparse
   use para, only : stdout, eps12, OmegaNum
   implicit none

   !* inout variables
   !> number of lanczos vectors
   integer, intent(in) :: NumLczVectors
   integer, intent(inout) :: NumLczVectors_out
   integer, intent(in) :: Mdim
   integer, intent(in) :: nnz
   integer, intent(in) :: icsr(Mdim+1)
   integer, intent(in) :: jcsr(nnz)
   complex(dp), intent(in) :: InitialVector(Mdim)
   complex(dp), intent(in) :: acsr(nnz)

   complex(dp), intent(out) :: Alpha (NumLczVectors)
   complex(dp), intent(out) :: Betan (NumLczVectors)

   !* local variables
   !* loop index
   integer :: i
   real(dp) :: time_start, time_end, time_start0

   !* Lanczos vector
   complex(dp), allocatable  :: vec_0(:), vec_1(:), vec_2(:)

   !* temp variables
   real(dp) :: error
   logical :: converged
   complex(dp) :: t_z1, t_z2, zdotc
   real(dp), allocatable :: dos_old(:)

   !if (cpuid==0) write(stdout, '(2x,a)')'>> Lanczos algorithm'

   allocate(dos_old(OmegaNum))
   allocate(vec_0(Mdim), vec_1(Mdim), vec_2(Mdim))
   dos_old= 0d0
   vec_0= 0d0; vec_1= 0d0; vec_2= 0d0

   alpha= 0d0
   betan= 0d0

   call now(time_start0)
   time_start= time_start0

   !** run Lanczos procedure to get a set of lanczos vectors <<!

   !* initialize the first lanczos vector
   vec_0= InitialVector

   !* normalize \phi_0
   !* <\phi_0|\phi_0>
   t_z1= zdotc(Mdim, vec_0, 1, InitialVector, 1)
   if (dble(t_z1)<eps12) t_z1= eps12 !< avoiding computational error

   !> the norm of the initial Lanczos vector should not be zero
   betan(1)= dsqrt(dble(t_z1))

   !* vec_0= vec_0/betan(1)
   vec_0= 1d0/betan(1)* vec_0

   !* perform |vec(1,:)>=Ham |vec(0,:)>
!#if defined (INTELMKL)
!   call mkl_zcsrgemv('N', Mdim, acsr, icsr, jcsr, vec_0, vec_1)
!#endif
   call csrmv_z(Mdim, nnz, acsr, icsr, jcsr, vec_0, vec_1)

   !* a_0= <phi_0|H|phi_0>
   alpha(1) = zdotc(Mdim, vec_0, 1, vec_1, 1)

   !* |phi_1>'=H|phi_0>-a_0|phi_0>
   !call WTSeqVecAxpy(-alpha(1), vec_0, vec_1)
   call zaxpy(mdim, -alpha(1), vec_0, 1, vec_1, 1)

   betan(2) = zdotc(Mdim, vec_1, 1, vec_1, 1)

   if (dble(betan(2))<eps12) betan(2)= eps12 !< avoiding computational error
   betan(2) = dsqrt(dble(betan(2)))

   !* rescale vec_1
   vec_1= 1d0/betan(2)*vec_1

   !* time measurement
   call now(time_end)

   do i= 2, NumLczVectors
      if (Mdim>500000.and.mod(i, 100)==0.and.cpuid.eq.0) &
         write(stdout, '(2a, i9, "  /", i10, a, f10.1, "s", a, f10.1, "s")') &
         '  In lanczos_seqsparse_cpu_z ', ' i/NumLczVectors', i, NumLczVectors, &
         ' time elapsed: ', time_end-time_start0, &
         ' time left: ', (NumLczVectors-i+1d0)*(time_end-time_start)
      call now(time_start)


!#if defined (INTELMKL)
!      call mkl_zcsrgemv('N', Mdim, acsr, icsr, jcsr, vec_1, vec_2)
!#endif
      call csrmv_z(Mdim, nnz, acsr, icsr, jcsr, vec_1, vec_2)

      !* calculate alpha, betan
      t_z1= zdotc(Mdim, vec_1, 1, vec_2, 1)
      alpha(i)= dble(t_z1)

      !* vec_2= vec_2- betan(i)*vec_0- alpha(i)*vec_1
      !call WTSeqVecAxpy(-alpha(i), vec_1, vec_2)
      call zaxpy(mdim, -alpha(i), vec_1, 1, vec_2, 1)
      !call WTSeqVecAxpy(-betan(i), vec_0, vec_2)
      call zaxpy(mdim, -betan(i), vec_0, 1, vec_2, 1)


      if (i.eq.NumLczVectors)cycle

      t_z2= zdotc(Mdim, vec_2, 1, vec_2, 1)
      if (dble(t_z2)<eps12) t_z2= eps12 !< avoiding computational error
      betan(i+1)= dsqrt(dble(t_z2))

      vec_0= vec_1

      vec_2= 1/betan(i+1)* vec_2

      vec_1= vec_2

      !> check convergence
      converged= .false.
      if (mod(i, 10).eq.0)then
         call lanczos_converged(alpha, betan, NumLczVectors, i, dos_old,&
                             error, .true., converged)
      endif
      if (converged)then
         NumLczVectors_out= i
         if (cpuid==0) then
            write(stdout, *)' >> Lanczos procedure is converged at iteration: ', i
            write(stdout, '(a, E16.5)')' >> Estimated error is : ', error
         endif
         exit
      endif

      call now(time_end)
   enddo  ! iteration

1000 continue

   if (.not.converged)then
       NumLczVectors_out= NumLczVectors
       if (cpuid==0) then
          write(stdout, '(a, i16)')' >> Lanczos procedure is not converged at iteration: ', NumLczVectors
          write(stdout, '(a, E16.5)')' >> Estimated error is : ', error
       endif
    endif


   !if (cpuid==0) then
   !   write(stdout, '(a)') '  ilcz       alpha           betan'
   !   do i = 1, NumLczVectors
   !      write(stdout, '(i5,40f16.8)')i,dble(alpha(i)), dble(betan(i))
   !   enddo
   !endif

   return
end subroutine lanczos_seqsparse_cpu_z


 !> Lanczos algorithm
 !> NumLczVectors : number of Lanczos vectors
 !> Mdim : The dimension of the Ham when it becomes a dense matrix
 !> input: initial vector; Matrix stored in PARCSR format
 !> output : Alpha(n), Beta(n)
subroutine lanczos_parsparse_cpu_z(NumLczVectors, Mdim, Ham, InitialVector, Alpha, Betan)
   use prec
   use sparse
   use wmpi
   use para, only : stdout, eps12
   implicit none

   !* inout variables
   !> number of lanczos vectors
   integer, intent(in) :: NumLczVectors, Mdim
   type (WTParCSR), intent(in) :: Ham
   type (WTParVec), intent(in) :: InitialVector
   complex(dp), intent(out) :: Alpha (NumLczVectors)
   complex(dp), intent(out) :: Betan (NumLczVectors)

   !* local variables
   !* loop index
   integer :: i

   !* Lanczos vector
   type (WTParVec), pointer :: vec_0
   type (WTParVec), pointer :: vec_1
   type (WTParVec), pointer :: vec_2

   !* temp variables
   complex(dp) :: t_z1
   complex(dp) :: t_z2


   if (cpuid==0) write(stdout, '(2x,a)')'>> Lanczos algorithm'
   !print *, '%=============================================================%'
   !print *, '|             Lanczos on cpu with sparse matrix               |'
   !print *, '%=============================================================%'


   vec_0=> Null()
   vec_1=> Null()
   vec_2=> Null()
   allocate(vec_0, vec_1, vec_2)
#if defined (MPI)
   call WTParVecCreate(mpi_cmw, Mdim, vec_0)
   call WTParVecCreate(mpi_cmw, Mdim, vec_1)
   call WTParVecCreate(mpi_cmw, Mdim, vec_2)
#endif
   call WTParVecInitialize(vec_0)
   call WTParVecInitialize(vec_1)
   call WTParVecInitialize(vec_2)

   alpha= 0d0
   betan= 0d0


   !** run Lanczos procedure to get a set of lanczos vectors <<!

   !* initialize the first lanczos vector
   vec_0= InitialVector

   !* normalize \phi_0
   !* <\phi_0|\phi_0>
   t_z1 = vec_0*InitialVector
   if (dble(t_z1)<eps12) t_z1= eps12 !< avoiding computational error

   !> the norm of the initial Lanczos vector should not be zero
   betan(1)= dsqrt(dble(t_z1))

   !* vec_0= vec_0/betan(1)
   vec_0= 1d0/betan(1)* vec_0

   !* perform |vec(1,:)>=Ham |vec(0,:)>
   call WTParCSRMatVec(Ham, vec_0, vec_1)

   !* a_0= <phi_0|H|phi_0>
   alpha(1)= vec_0*vec_1

   !* |phi_1>'=H|phi_0>-a_0|phi_0>
   call WTParVecAxpy(-alpha(1), vec_0, vec_1)

   betan(2) = vec_1*vec_1

   if (dble(betan(2))<eps12) betan(2)= eps12 !< avoiding computational error
   betan(2) = dsqrt(dble(betan(2)))

   !* rescale vec_1
   vec_1= 1d0/betan(2)*vec_1


   do i= 2, NumLczVectors

      call WTParCSRMatVec(Ham, vec_1, vec_2)

      !* calculate alpha, betan
      t_z1= vec_1*vec_2
      alpha(i)= dble(t_z1)

      !* vec_2= vec_2- betan(i)*vec_0- alpha(i)*vec_1
      call WTParVecAxpy(-alpha(i), vec_1, vec_2)
      call WTParVecAxpy(-betan(i), vec_0, vec_2)


      if (i.eq.NumLczVectors)cycle

      t_z2= vec_2*vec_2
      if (dble(t_z2)<eps12) t_z2= eps12 !< avoiding computational error
      betan(i+1)= dsqrt(dble(t_z2))

      vec_0= vec_1

      vec_2= 1/betan(i+1)* vec_2

      vec_1= vec_2

   enddo

1000 continue

   if (cpuid==0) then
      write(stdout, '(a)') '  ilcz       alpha           betan'
      do i = 1, NumLczVectors
         write(stdout, '(i5,40f16.8)')i,dble(alpha(i)), dble(betan(i))
      enddo
   endif


   call WTParVecDestroy(vec_0)
   call WTParVecDestroy(vec_1)
   call WTParVecDestroy(vec_2)
   nullify(vec_0, vec_1, vec_2)


   return
end subroutine lanczos_parsparse_cpu_z

subroutine LandauLevel_B_dos_Lanczos
   !> we calculate the the spectrum with given magnetic field strength indicated by Magp and Magq,
   !> a serials of kpoints kpoints(3, NK) and a serial of energies Omega(OmegaNum).
   !> the output is dos_B_omega(Nk, OmegaNum)
   !> Bx and By are global variables, you have to set it up in the current subroutine
   use prec
   use sparse
   use wmpi
   use mt19937_64
   use para, only : Magq, Num_Wann, Bx, By, zi, pi, Fermi_broadening, iso_energy, &
      OmegaNum, OmegaMin, OmegaMax,  Magp, stdout, Magp_min, Magp_max, nnzmax_input, &
      outfileindex, Single_KPOINT_3D_DIRECT,splen,Is_Sparse_Hr, eV2Hartree, &
      MagneticSuperProjectedArea,ijmax,NumLCZVecs, NumRandomConfs, Add_Zeeman_Field
   implicit none

   !> magnetic field strength, this number should compatiable with the magnetic supercell
   !> which means B*AreaOfMagneticSupercell=2*pi
   !> The size of magnetic supercell is controled by Nq.
   !> here Magp should be integer from 1 to Nq

   !> energy interval
   !> OmegaNum is defined in the module.f90 and read from the input.dat or wt.in
   real, allocatable :: omega(:)
   real(dp) :: time_start, time_end, time_start0

   !> spectrum calculated
   real(dp), allocatable :: dos_B_omega(:, :, :), dos_B_omega_mpi(:, :, :)
   real(dp), allocatable ::  n_int(:, :)

   integer :: i,  Mdim, Nq, ie, ib, ierr, it, ieta, iter
   integer :: nnzmax, nnz, NumLczVectors, Nmag, NumLczVectors_out
   real(dp) :: k3(3), thetaj, B0, energy, B0Tesla_quantumflux_magsupcell, eta

   !> magnetic field strength
   real(dp), allocatable :: mag_Tesla(:), flux(:)

   !> left vector and right vector
   complex(dp), allocatable :: InitialVector(:), InitialVectors(:, :)

   !> Lanczos hamiltonian
   complex(dp), allocatable :: Alpha(:), Betan(:)

   integer, allocatable :: icsr(:), jcsr(:), iwk(:)
   complex(dp), allocatable :: acsr(:)
   complex(dp) :: zdotc, norm, theta
   real(dp) :: continued_fraction
   logical :: term

   integer :: NumberofEta, ie_Earc
   real(dp), allocatable :: eta_array(:), n_Earc(:)

   real(dp) :: Bx_in_Tesla, By_in_Tesla, Bz_in_Tesla
   integer(8) :: iseed
   real(dp) :: rseed0

   NumberofEta=9
   allocate(eta_array(NumberofEta))
   allocate(n_Earc(NumberofEta))
   eta_array=(/0.1d0, 0.2d0, 0.4d0, 0.8d0, 1.0d0, 2d0, 4d0, 8d0, 10d0/)
   eta_array= eta_array*Fermi_broadening


   Nq= Magq
   Nmag= Magp_max- Magp_min +1
   Mdim= Num_Wann*Nq
   if (nnzmax_input<0)then
      nnzmax= Num_wann*(2*ijmax+1+2)*Mdim
      if(Is_Sparse_Hr) nnzmax=splen*Nq
   else
      nnzmax= nnzmax_input
   endif
   if (cpuid==0) then
      write(stdout, '(a,i8)')' Magnetic supercell is Nq= ', Nq
      write(stdout, '(a,i18)')' Hamiltonian matrix dimension Mdim= ', Mdim
      write(stdout, '(a,i18)')' nnzmax for lanczos is ', nnzmax
      write(stdout, '(a)')' Maximal integer for 32-digit is 2,147,483,647'
      write(stdout, '(a, f20.1, a)')' Memory for dos matrix is ',  &
         (Nmag)/1024d0/1024d0*9*OmegaNum*16d0*2d0, ' MB'
      write(stdout, '(a, f20.1, a)')' Memory for sparse Hamiltonian is ', &
         nnzmax/1024d0/1024d0*(4d0+4d0+16d0), ' MB'
      write(stdout, '(a, f20.1, a)')' Memory for Lanczos vectors is ', &
         Mdim/1024d0/1024d0*(4d0+NumRandomConfs)*16d0, ' MB'
      write(stdout, '(a, f20.1, a)')' Memory for DOS storage is ', &
         Nmag*OmegaNum*NumberofEta*2d0/1024d0/1024d0*8d0, ' MB'
   endif

   NumLczVectors = NumLCZVecs
   NumLczVectors_out = NumLCZVecs
   if (NumLczVectors>Mdim) then
      NumLczVectors = Mdim-1
   endif

   allocate(omega(OmegaNum), mag_Tesla(Nmag), flux(Nmag))
   allocate(n_int(0:Omeganum, NumberofEta))
   allocate(dos_B_omega(Nmag, OmegaNum, NumberofEta), dos_B_omega_mpi(Nmag, OmegaNum, NumberofEta))
   dos_B_omega= 0d0; dos_B_omega_mpi= 0d0

   !> energy
   do ie=1, OmegaNum
      omega(ie)= OmegaMin+ (OmegaMax-OmegaMin)* (ie-1d0)/dble(OmegaNum-1)
   enddo ! ie

   !> deal with the magnetic field
   !> first transform the Bx By into B*Cos\theta, B*Sin\theta
   if (abs(By)<1e-8) then
      theta= 0d0
   elseif (By>0) then
      theta = atan(Bx/By)
   else
      theta = atan(Bx/By)+ pi
   endif
   if (abs(theta)<1e-8 .and. Bx<0 ) theta= theta+ pi

   !> when B0=2*pi/Nq, the magnetic flux in the magnetic supercell is 2*pi
   !> the magnetic flux in the primitive cell is 2*pi/(Nq)
   B0= 2d0*pi/dble(Nq)
   B0= abs(B0)

   !> transform it into Tesla
   !> B=2*pi*\phi_0/S0, where \phi_0 is the quantum flux, S0 is the projected area
   !> \phi_0 = h/2e  h=6.62607004*1E-34, e= 1.6*1E-19
   !> B0Tesla is the magnetic field when the flux of the magnetic-supercell is 2*pi
   B0Tesla_quantumflux_magsupcell= 6.62607004d0*1E-34/2d0/1.6021766208d0*1E19/MagneticSuperProjectedArea*1E20

   if (cpuid==0) then
      write(stdout, '(a, 2E18.8)')' Magnetic field B in Tesla that makes the flux through '
      write(stdout, '(a, 2E18.8)')' the magnetic supercell to be one quantum flux : ', B0Tesla_quantumflux_magsupcell
   endif


   !> The flux in the unit cell changes from 0 to 2*pi*Nmag/Nq
   do ib= Magp_min, Magp_max
      flux(ib-Magp_min+1)= B0* ib
      mag_Tesla(ib-Magp_min+1)= B0Tesla_quantumflux_magsupcell* ib
   enddo

   k3=Single_KPOINT_3D_DIRECT

   allocate(InitialVector(Mdim))
   allocate(InitialVectors(Mdim, NumRandomConfs))
   allocate(Alpha(NumLczVectors), Betan(NumLczVectors))
   allocate(icsr(nnzmax), jcsr(nnzmax), acsr(nnzmax))
   allocate(iwk (Mdim+1))

   !* Get the Hamiltonian
   time_start= 0d0
   time_start0= 0d0
   call now(time_start0)
   time_start= time_start0
   time_end  = time_start0

   !> Fix the initial vectors for each magnetic flux
   do it=1, NumRandomConfs
      call now(rseed0)
      iseed=cpuid*NumRandomConfs*3132+ it+123431+int(rseed0)
      call init_genrand64(iseed)
      do i=1, Mdim
        !call random_number(harvest= thetaj)
         thetaj= genrand64_real2()
         InitialVectors(i, it)= exp(zi*2d0*pi*thetaj)
     enddo
   enddo

   do iter=1+cpuid, Nmag*NumRandomConfs, num_cpu
  !do ib=1+cpuid, Nmag, num_cpu
  !   do it= 1, NumRandomConfs
      ib= (iter-1)/NumRandomConfs+ 1
      it= (iter-1 - (ib-1)*NumRandomConfs)+ 1
         if (cpuid.eq.0) &
            write(stdout, '(2a, i9, "  /", i10, a, f10.1, "s", a, f10.1, "s")') &
            'In LandauLevel_B_dos_Lanczos ', ' iter/Nmag*NumRandomConfs ', iter, Nmag*NumRandomConfs, &
            ' time elapsed: ', time_end-time_start0, &
            ' time left: ', ((Nmag*NumRandomConfs-iter-1d0)/dble(num_cpu))*(time_end-time_start)

         call now(time_start)

         !> here we set the magnetic field along the first vector of the SURFACE card
         Bx= -flux(ib)
         By= 0d0

         !> TODO: should be revised in future, here we only put magnetic field along z direction
         Bx_in_Tesla=0d0; By_in_Tesla= 0d0; Bz_in_Tesla=  mag_Tesla(ib)

         Alpha= 0d0; Betan= 0d0
         icsr=0; jcsr=0;acsr=0d0;iwk=0d0
         InitialVector(:)= InitialVectors(:, it)
         norm= zdotc(Mdim, InitialVector, 1, InitialVector, 1)

         InitialVector= InitialVector/dsqrt(dble(norm))

         nnz= nnzmax
         if(Is_Sparse_Hr) then

            !> this subroutine will change the tight binding model according to the magnetic field strength
            if (Add_Zeeman_Field) call add_zeeman_sparse_hr(Bx_in_Tesla, By_in_Tesla, Bz_in_Tesla)

            !> get the Hamiltonian for the magnetic supercell whose size is Nq.
            call ham_3Dlandau_sparseHR(nnz,Mdim,Nq,k3,acsr,jcsr,icsr)
         else
            call ham_3Dlandau_sparse1(nnz, Mdim, Nq, k3, acsr,jcsr,icsr)
         end if

         !> transform coo format to csr format
         call ConvertCooToCsr(Mdim, nnz, acsr, icsr, jcsr, iwk)
         call csr_sort_indices(Mdim, nnz, icsr, jcsr, acsr)

         !> eleminate the same entries in the sparse matrix
         call csr_sum_duplicates(Mdim, nnz, icsr, jcsr, acsr)
         call csr_sort_indices(Mdim, nnz, icsr, jcsr, acsr)

         call now(time_end)
         if (cpuid.eq.0) write(stdout, '(a, f10.1, "s")') &
            '  Hamiltonian construction time cost :', time_end-time_start   

         !* doing lanczos procedure in order to get alpha, beta
         call lanczos_seqsparse_cpu_z(NumLczVectors, NumLczVectors_out, Mdim, nnz, icsr, jcsr, acsr, InitialVector, Alpha, Betan)

         !> Betan(1) is meaningless.
         term = .True.
         do ieta=1, NumberofEta
            eta= eta_array(ieta)
            do ie=1, OmegaNum
               energy= omega(ie)
               dos_B_omega(ib, ie, ieta)= dos_B_omega(ib, ie, ieta)+ &
                  continued_fraction(Alpha, Betan, energy, eta, NumLczVectors_out, term)
            enddo
         enddo
         call now(time_end)
     !enddo ! it= 1, NumRandomConfs
  !enddo ! ib magnetic field
   enddo ! combing ib and it

#if defined (MPI)
   call mpi_allreduce(dos_B_omega, dos_B_omega_mpi, size(dos_B_omega), &
      mpi_dp, mpi_sum, mpi_cmw, ierr)
#else
   dos_B_omega_mpi= dos_B_omega
#endif
   dos_B_omega_mpi= dos_B_omega_mpi/NumRandomConfs

   outfileindex= outfileindex+ 1
   if (cpuid.eq.0)then
      open (unit=outfileindex, file='LandauLevel_B_dos.dat')
      open (unit=outfileindex+1,file='wannierdiagram.dat')
      write(outfileindex, '("#", a)')'Hofstadter butterfly (Landau level in (E,B) axis) '
      write(outfileindex, '("#", a, i6)')'Magnetic supercell size : ', Magq
      write(outfileindex+1, '("#", a)')'Wannier diagram '
      write(outfileindex+1, '("#", a, i6)')'Magnetic supercell size : ', Magq
      write(outfileindex, '("#", a, 3f12.6, a)')'k points in fractional coordinates (', k3, ')'
      write(outfileindex, '("# Column", I5, 100I16)')(i, i=1, 12)
      write(outfileindex, '("#", a15, 5a16)', advance='NO')'Flux', 'B(Tesla)', 'E(eV)'
      write(outfileindex+1, '("# Column", I5, 100I16)')(i, i=1, 20)
      write(outfileindex+1, '("#", a15, 2a16)', advance='NO')'Flux', 'B(Tesla)'
      do ieta=1, NumberofEta-1
         write(outfileindex, '(a16)', advance='NO') 'A(B, E)'
         write(outfileindex+1, '(2a16)', advance='NO') 'n   ', 'A(n, E)'
      enddo
      write(outfileindex, '(a16)') 'A(B, E)'
      write(outfileindex+1, '(2a16)') 'n   ', 'A(n, E)'

      write(outfileindex, '("#", a, 20X, 300f16.2)')'Broadening \eta (meV): ', Eta_array(:)*1000d0/eV2Hartree
      write(outfileindex+1, '("#", a, 5X, f26.2,300f32.2)')'Broadening \eta (meV): ', Eta_array(:)*1000d0/eV2Hartree

      !> find n_int at EF
      do ie=1, omeganum
         if (omega(ie)>iso_energy) then
            ie_Earc= ie- 1
            exit
         endif
      enddo

      do ib=1, Nmag
         n_int=0
         do ie=1, omeganum
            write(outfileindex, '(300f16.6)')flux(ib)/2d0/pi, mag_Tesla(ib), omega(ie)/eV2Hartree, dos_B_omega_mpi(ib, ie, :)
         enddo

         !> get number of electrons between the lowest energy level and iso_energy
         n_Earc= 0d0
         do ie=1, ie_Earc
            n_Earc(:)= n_Earc(:)+ dos_B_omega_mpi(ib, ie, :)
         enddo

         !> get number of electrons between the lowest energy level and omega(ie)
         do ie=1, omeganum
            n_int(ie, :)=n_int(ie-1, :)+ dos_B_omega_mpi(ib, ie, :)
         enddo

         !> set n(E)= n_int- n_Earc= \int_Earc^E \rho(\epsilon)d\epsilon
         do ie=1, omeganum
            n_int(ie, :)= n_int(ie, :)-n_Earc(:)
         enddo

        !do ieta=1, NumberofEta
        !   n_int(:, ieta)=n_int(:, ieta)/n_int(omeganum, ieta)
        !enddo

         do ie=1, omeganum
            write(outfileindex+1,'(300f16.6)')  flux(ib)/2d0/pi, mag_Tesla(ib), &
               (n_int(ie, ieta), dos_B_omega_mpi(ib, ie, ieta), ieta=1, NumberofEta)
         enddo

         write(outfileindex+1, *) ' '
         write(outfileindex, *) ' '
      enddo
      close(outfileindex)
      write(stdout,*)'calculate Landau level spectrum in B-E mode successfully'
   endif

   outfileindex= outfileindex+ 2
   if (cpuid.eq.0)then
      open (unit=outfileindex, file='LandauLevel_B_dos.gnu')
      write(outfileindex, '(a)')'#set terminal  postscript enhanced color font ",24"'
      write(outfileindex, '(a)')'set terminal pngcairo enhanced color font ",60" size 1920, 1680'
      write(outfileindex, '(a)')"set output 'LandauLevel_B_dos.png'"
      write(outfileindex, '(a)')'set pm3d'
      write(outfileindex, '(a)')'set palette rgb 21,22,23'
      write(outfileindex, '(a)')'#set isosamples 50,50'
      write(outfileindex, '(a)')'set size 0.9,1'
      write(outfileindex, '(a)')'set origin 0.05,0'
      write(outfileindex, '(a)')'set view map'
      write(outfileindex, '(a)')'unset ztics'
      write(outfileindex, '(a)')'unset surface'
      write(outfileindex, '(a)')'unset key'
      write(outfileindex, '(a)')'#set xtics font ",24"'
      write(outfileindex, '(a)')'#set ytics font ",24"'
      write(outfileindex, '(a)')'#set xlabel font ",24"'
      write(outfileindex, '(a)')'set ylabel "Energy (eV)"'
      write(outfileindex, '(a, i6,a)')'set title "Hofstadter butterfly with Nq=', Nq, '" font ",40"'
      write(outfileindex, '(a, f12.6, a, f12.6, a)')'set yrange [',OmegaMin/eV2Hartree,':',OmegaMax/eV2Hartree,']'
      write(outfileindex, '(a)')'set xlabel "Phi per unit cell"'
      write(outfileindex,*) '#set xlabel "{/Symbol F}/{/Symbol F}_0"'
      write(outfileindex, '(a, f12.6, a)')     'set xrange [ 0.0000 :', maxval(flux/2d0/pi) ,']'
      write(outfileindex, '(a)')     "splot 'LandauLevel_B_dos.dat' u 1:3:(log($4)) w pm3d"
      write(outfileindex, '(a)')'set xlabel "B (Tesla)"'
      write(outfileindex, '(a, f12.6, a)')     '#set xrange [ 0.0000 :', maxval(mag_Tesla) ,']'
      write(outfileindex, '(a)')     "#splot 'LandauLevel_B_dos.dat' u 2:3:(log($4)) w pm3d"
   endif

   outfileindex=outfileindex+1
   if(cpuid == 0) then
      open(unit=outfileindex, file='wannierdiagram.gnu')
      write(outfileindex,*) 'set terminal pngcairo enhanced color font ",60" size 1920, 1680'
      write(outfileindex,*) "set output 'wannierdiagram.png'"
      write(outfileindex,*) 'set pm3d'
      write(outfileindex, '(a)')'set size 0.9,1'
      write(outfileindex, '(a)')'set origin 0.05,0'
      write(outfileindex, '(a)')'set palette rgb 21,22,23'
      write(outfileindex,*) '#set isosamples 50,50'
      write(outfileindex,*) 'set view map'
      write(outfileindex,*) 'unset ztics'
      write(outfileindex,*) 'unset surface'
      write(outfileindex,*) 'unset key'

      write(outfileindex,*) 'set ylabel "n"'
      write(outfileindex, '(a, i6, a)') 'set title "Wannier diagram with Nq=', Nq, '"'
      write(outfileindex,*) '#set yrange [   ] noextend'
      write(outfileindex,*) '#set xlabel "{/Symbol F}/{/Symbol F}_0"'
      write(outfileindex,*) 'set xlabel "Phi/Phi_0 per unit cell"'
      write(outfileindex,*) '#set xlabel "{/Symbol F}/{/Symbol F}_0"'
      write(outfileindex,*) '#set xrange [ ] noextend'

      write(outfileindex,*) "splot 'wannierdiagram.dat' u 1:3:(log($4)) w pm3d #lc palette"
   end if

#if defined (MPI)
   call mpi_barrier(mpi_cmw, ierr)
#endif

   deallocate(InitialVector, Alpha, Betan, InitialVectors)
   deallocate(icsr, jcsr, acsr, iwk)

   return
end subroutine LandauLevel_B_dos_Lanczos


subroutine LandauLevel_k_dos_Lanczos
   !> we calculate the the spectrum with given magnetic field strength indicated by Magp and Magq,
   !> a serials of kpoints kpoints(3, NK) and a serial of energies Omega(OmegaNum).
   !> the output is dos_k_omega(Nk, OmegaNum)
   use prec
   use sparse
   use wmpi
   use para, only : Magq, Num_Wann, Bx, By, zi, pi, Fermi_broadening, Angstrom2atomic, &
      OmegaNum, OmegaMin, OmegaMax, nk3_band, Magp, stdout, kpath_3d, eV2Hartree, &
      outfileindex, K3len_mag,splen,Is_Sparse_Hr,ijmax,NumLCZVecs, MagneticSuperProjectedArea, &
      Nk3lines, k3line_mag_stop, k3line_name, NumRandomConfs, nnzmax_input
   implicit none

   !> magnetic field strength, this number should compatiable with the magnetic supercell
   !> which means B*AreaOfMagneticSupercell=2*pi
   !> The size of magnetic supercell is controled by Nq.
   !> here Magp should be integer from 1 to Nq

   !> energy interval
   !> OmegaNum is defined in the module.f90 and read from the input.dat or wt.in
   real, allocatable :: omega(:)
   real(dp) :: time_start, time_end, time_start0

   !> spectrum calculated
   real(dp), allocatable :: dos_k_omega(:, :), dos_k_omega_mpi(:, :)

   integer :: i,  Mdim, Nq, ie, ik, ierr, it
   integer :: nnzmax, nnz, NumLczVectors, NumLczVectors_out
   real(dp) :: k3(3), thetaj, B0, energy, B0Tesla, B0Tesla_quantumflux_magsupcell

   !> left vector and right vector
   complex(dp), allocatable :: InitialVector(:)

   !> Lanczos hamiltonian
   complex(dp), allocatable :: Alpha(:), Betan(:)

   integer, allocatable :: icsr(:), jcsr(:), iwk(:)
   complex(dp), allocatable :: acsr(:)
   complex(dp) :: zdotc, norm, theta
   real(dp) :: continued_fraction
   logical :: term

   Nq= Magq
   Mdim= Num_Wann*Magq
   !> need to be checked
   !if(Is_Sparse_Hr) nnzmax=Nq*splen
   if (nnzmax_input<0)then
      nnzmax= Num_wann*(2*ijmax+1+2)*Mdim
      if(Is_Sparse_Hr) nnzmax=splen*Nq
   else
      nnzmax= nnzmax_input
   endif
 
   NumLczVectors = NumLCZVecs
   NumLczVectors_out = NumLCZVecs
   if (NumLczVectors>Mdim) then
      NumLczVectors = Mdim-1
   endif

   allocate(omega(OmegaNum))
   allocate(dos_k_omega(nk3_band, OmegaNum), dos_k_omega_mpi(nk3_band, OmegaNum))
   dos_k_omega= 0d0; dos_k_omega_mpi= 0d0

   !> energy
   do ie=1, OmegaNum
      omega(ie)= OmegaMin+ (OmegaMax-OmegaMin)* (ie-1d0)/dble(OmegaNum-1)
   enddo ! ie


   !> deal with the magnetic field
   !> first transform the Bx By into B*Cos\theta, B*Sin\theta
   if (abs(By)<1e-8) then
      theta= 0d0
   elseif (By>0) then
      theta = atan(Bx/By)
   else
      theta = atan(Bx/By)+ pi
   endif
   if (abs(theta)<1e-8 .and. Bx<0 ) theta= theta+ pi

   !> the magnetic flux in the magnetic supercell is 2*pi*Magp
   B0= 2d0*pi/dble(Nq)* Magp
   B0= abs(B0)
   Bx= B0* Cos(theta)
   By= B0* Sin(theta)

   !> transform it into Tesla
   !> B=2*pi*\phi_0/S0, where \phi_0 is the quantum flux, S0 is the projected area
   !> \phi_0 = h/2e  h=6.62607004*1E-34, e= 1.6*1E-19
   !> B0Tesla_quantumflux_magsupcell is the magnetic field that makes the flux through the magnetic
   !> supercell to be one quantum flux
   B0Tesla_quantumflux_magsupcell= 6.62607004*1E-34/2d0/1.6021766208*1E19/MagneticSuperProjectedArea*1E20
   B0Tesla= 6.62607004*1E-34/2d0/1.6021766208*1E19/MagneticSuperProjectedArea*1E20*Magp

   if (cpuid==0) then
      write(stdout, '(a, 2E18.8)')' Magnetic field B in Tesla that makes the flux through '
      write(stdout, '(a, 2E18.8)')' the magnetic supercell to be one quantum flux : ', B0Tesla_quantumflux_magsupcell
      write(stdout, '(a, 2E18.8)')' Magnetic field B in Tesla that fits to Magp : ', B0Tesla
      write(stdout, '(a, 2f18.8)')' Magnetic field Bx, By= ', Bx, By
   endif


   allocate(InitialVector(Mdim))
   allocate(Alpha(NumLczVectors), Betan(NumLczVectors))
   allocate(icsr(nnzmax), jcsr(nnzmax), acsr(nnzmax))
   allocate(iwk (Mdim+1))

   !* Get the Hamiltonian
   time_start= 0d0
   time_start0= 0d0
   call now(time_start0)
   time_start= time_start0
   time_end  = time_start0
   do ik=1+cpuid, nk3_band, num_cpu
      do it=1, NumRandomConfs
         if (cpuid.eq.0) &
            write(stdout, '(2a, i9, "  /", i10, a, f10.1, "s", a, f10.1, "s")') &
            'In LandauLevel_k_dos_Lanczos ', ' ik/NK ', ik, nk3_band, &
            ' time elapsed: ', time_end-time_start0, &
            ' time left: ', ((nk3_band-ik)/dble(num_cpu))*(time_end-time_start)*NumRandomConfs

         call now(time_start)
         Alpha= 0d0; Betan= 0d0
         icsr=0; jcsr=0;acsr=0d0;iwk=0d0
         do i=1, Mdim
            call random_number(harvest= thetaj)
            InitialVector(i)= exp(zi*2d0*pi*thetaj)
            !InitialVector(i)= 1d0
         enddo
         norm= zdotc(Mdim, InitialVector, 1, InitialVector, 1)

         InitialVector= InitialVector/dsqrt(dble(norm))

         if (cpuid==0) write(stdout, '(a, 2i10)') 'LandauLevel_k_DOS', ik, nk3_band
         k3= kpath_3d(:, ik)
         nnz= nnzmax
         if(Is_Sparse_Hr) then
            call ham_3Dlandau_sparseHR(nnz,Mdim,NQ,k3,acsr,jcsr,icsr)
         else
            call ham_3Dlandau_sparse1(nnz, Mdim, Nq, k3, acsr,jcsr,icsr)
         end if

         !> transform coo format to csr format
         call ConvertCooToCsr(Mdim, nnz, acsr, icsr, jcsr, iwk)
         call csr_sort_indices(Mdim, nnz, icsr, jcsr, acsr)

         !> eleminate the same entries in the sparse matrix
         call csr_sum_duplicates(Mdim, nnz, icsr, jcsr, acsr)
         call csr_sort_indices(Mdim, nnz, icsr, jcsr, acsr)

         !* doing lanczos procedure in order to get alpha, beta
         call lanczos_seqsparse_cpu_z(NumLczVectors, NumLczVectors_out, Mdim, nnz, icsr, jcsr, acsr, InitialVector, Alpha, Betan)

         !> Betan(1) is meaningless.
         term = .False.
         do ie=1, OmegaNum
            energy= omega(ie)
            dos_k_omega(ik, ie)=dos_k_omega(ik, ie)+  &
               continued_fraction(Alpha,Betan,energy,Fermi_broadening,NumLczVectors_out, term)
         enddo
         call now(time_end)
      enddo ! it= 1, NumRandomConfs
   enddo ! ik

#if defined (MPI)
   call mpi_allreduce(dos_k_omega, dos_k_omega_mpi, size(dos_k_omega), &
      mpi_dp, mpi_sum, mpi_cmw, ierr)
#else
   dos_k_omega_mpi= dos_k_omega
#endif
   dos_k_omega_mpi= dos_k_omega_mpi/NumRandomConfs

   outfileindex= outfileindex+ 1
   if (cpuid.eq.0)then
      open (unit=outfileindex, file='LandauLevel_k_dos.dat')
      write(outfileindex, '("#", a14, f16.8)')'Magnetic field in Tesla: ', B0Tesla
      write(outfileindex, '("#", a14, 2a15, a)')'k', 'E(eV)', 'LDOS'
      do ik=1, nk3_band
         do ie=1, omeganum
            write(outfileindex, '(30f16.8)')K3len_mag(ik)*Angstrom2atomic, omega(ie)/eV2Hartree, dos_k_omega_mpi(ik, ie)
         enddo
         write(outfileindex, *) ' '
      enddo
      close(outfileindex)
      write(stdout,*)'calculate Landau level spectrum in k-E mode successfully'
   endif

   outfileindex= outfileindex+ 1
   if (cpuid.eq.0)then
      open (unit=outfileindex, file='LandauLevel_k_dos.gnu')
      write(outfileindex, '(a)')'#set terminal  postscript enhanced color font ",24"'
      write(outfileindex, '(a)')'set terminal pngcairo enhanced color font ",70" size 1920, 1680'
      write(outfileindex, '(a)')"set output 'LandauLevel_k_dos.png'"
      write(outfileindex, '(a)')'set pm3d'
      write(outfileindex, '(a)')'#set isosamples 50,50'
      write(outfileindex, '(a)')'set size 0.85, 1'
      write(outfileindex, '(a)')'set origin 0.1, 0'
      write(outfileindex, '(a)')'set view map'
      write(outfileindex, '(a)')'unset ztics'
      write(outfileindex, '(a)')'unset surface'
      write(outfileindex, '(a)')'unset key'
      write(outfileindex, '(a)')'#set xtics font ",24"'
      write(outfileindex, '(a)')'#set ytics font ",24"'
      write(outfileindex, '(a)')'#set xlabel font ",24"'
      write(outfileindex, '(a)')'set ylabel "Energy (eV)"'
      write(outfileindex, '(a, i6,a)')'set title "Landau level with Nq=', Nq, '"'
      write(outfileindex, '(a, f18.5, a)')'set xrange [0: ', maxval(K3len_mag*Angstrom2atomic), ']'
      write(outfileindex, '(a, f18.5, a, f18.5, a)')'set yrange [', OmegaMin/eV2Hartree, ':', OmegaMax/eV2Hartree, ']'
      write(outfileindex, 202, advance="no") (k3line_name(i), k3line_mag_stop(i), i=1, Nk3lines)
      write(outfileindex, 203)k3line_name(nk3lines+1), k3line_mag_stop(Nk3lines+1)

      do i=1, nk3lines-1
         write(outfileindex, 204)k3line_mag_stop(i+1), OmegaMin, k3line_mag_stop(i+1), OmegaMax
      enddo
      write(outfileindex, '(a)')'set pm3d interpolate 2,2'
      write(outfileindex, '(2a)')"splot 'LandauLevel_k_dos.dat' u 1:2:(log($3)) w pm3d"
      close(outfileindex)
   endif

202 format('set xtics (',:20('"',A1,'" ',F8.5,','))
203 format(A1,'" ',F8.5,')')
204 format('set arrow from ',F8.5,',',F10.5,' to ',F8.5,',',F10.5, ' nohead front lw 3')

#if defined (MPI)
   call mpi_barrier(mpi_cmw, ierr)
#endif
   deallocate(InitialVector, Alpha, Betan)
   deallocate(icsr, jcsr, acsr, iwk)

   return
end subroutine LandauLevel_k_dos_Lanczos


subroutine bulkbandk_dos_lanczos
   !> we calculate the the spectrum with given magnetic field strength indicated by Magp and Magq,
   !> a serials of kpoints kpoints(3, NK) and a serial of energies Omega(OmegaNum).
   !> the output is dos_k_omega(Nk, OmegaNum)
   use prec
   use sparse
   use wmpi
   use para, only : Magq, Num_Wann, Bx, By, zi, pi, Fermi_broadening, Angstrom2atomic, &
      OmegaNum, OmegaMin, OmegaMax, nk3_band, Magp, stdout, kpath_3d, nnzmax_input, &
      outfileindex, K3len,splen,Is_Sparse_Hr,ijmax,NumLCZVecs, eV2Hartree
   implicit none

   !> magnetic field strength, this number should compatiable with the magnetic supercell
   !> which means B*AreaOfMagneticSupercell=2*pi
   !> The size of magnetic supercell is controled by Nq.
   !> here Magp should be integer from 1 to Nq

   !> energy interval
   !> OmegaNum is defined in the module.f90 and read from the input.dat or wt.in
   real, allocatable :: omega(:)

   !> spectrum calculated
   real(dp), allocatable :: dos_k_omega(:, :), dos_k_omega_mpi(:, :)

   integer :: i, Mdim, Nq, ie, ik, ierr
   integer :: nnzmax, nnz, NumLczVectors, NumLczVectors_out
   real(dp) :: k3(3), thetaj, B0, energy

   !> left vector and right vector
   complex(dp), allocatable :: InitialVector(:)

   !> Lanczos hamiltonian
   complex(dp), allocatable :: Alpha(:), Betan(:)

   integer, allocatable :: icsr(:), jcsr(:), iwk(:)
   complex(dp), allocatable :: acsr(:)
   complex(dp) :: zdotc, norm, theta
   real(dp) :: continued_fraction
   logical :: term

   Nq= 1
   Mdim= Num_Wann*Nq
   if (nnzmax_input<0)then
      nnzmax= Num_wann*(2*ijmax+1+2)*Mdim
      if(Is_Sparse_Hr) nnzmax=splen*Nq
   else
      nnzmax= nnzmax_input
   endif
 
   !> need to be checked
   !if(Is_Sparse_Hr) nnzmax=Nq*splen

   NumLczVectors = NumLCZVecs
   if (NumLczVectors>Mdim) then
      NumLczVectors = Mdim-1
   endif

   allocate(omega(OmegaNum))
   allocate(dos_k_omega(nk3_band, OmegaNum), dos_k_omega_mpi(nk3_band, OmegaNum))
   dos_k_omega= 0d0; dos_k_omega_mpi= 0d0

   !> energy
   do ie=1, OmegaNum
      omega(ie)= OmegaMin+ (OmegaMax-OmegaMin)* (ie-1d0)/dble(OmegaNum-1)
   enddo ! ie


   !> deal with the magnetic field
   !> first transform the Bx By into B*Cos\theta, B*Sin\theta
   if (abs(By)<1e-8) then
      theta= 0d0
   elseif (By>0) then
      theta = atan(Bx/By)
   else
      theta = atan(Bx/By)+ pi
   endif
   if (abs(theta)<1e-8 .and. Bx<0 ) theta= theta+ pi

   !> the magnetic flux in the magnetic supercell is 2*pi*iq

   allocate(InitialVector(Mdim))
   allocate(Alpha(NumLczVectors), Betan(NumLczVectors))
   allocate(icsr(nnzmax), jcsr(nnzmax), acsr(nnzmax))
   allocate(iwk (Mdim+1))

   !* Get the Hamiltonian
   do ik=1+cpuid, nk3_band, num_cpu
      Alpha= 0d0; Betan= 0d0
      icsr=0; jcsr=0;acsr=0d0;iwk=0d0
      do i=1, Mdim
         call random_number(harvest= thetaj)
         InitialVector(i)= exp(zi*2d0*pi*thetaj)
      enddo
      norm= zdotc(Mdim, InitialVector, 1, InitialVector, 1)

      InitialVector= InitialVector/dsqrt(dble(norm))

      if (cpuid==0) write(stdout, '(a, 2i10)') 'LandauLevel_k_DOS', ik, nk3_band
      k3= kpath_3d(:, ik)
      nnz= nnzmax
      if(Is_Sparse_Hr) then
         call ham_bulk_coo_sparsehr(k3, acsr,jcsr,icsr)
      else
         stop " ERORR: we don't support this calculation"
      end if

      !> transform coo format to csr format
      call ConvertCooToCsr(Mdim, nnz, acsr, icsr, jcsr, iwk)
      call csr_sort_indices(Mdim, nnz, icsr, jcsr, acsr)

      !> eleminate the same entries in the sparse matrix
      call csr_sum_duplicates(Mdim, nnz, icsr, jcsr, acsr)
      call csr_sort_indices(Mdim, nnz, icsr, jcsr, acsr)

      !* doing lanczos procedure in order to get alpha, beta
      call lanczos_seqsparse_cpu_z(NumLczVectors, NumLczVectors_out, Mdim, nnz, icsr, jcsr, acsr, InitialVector, Alpha, Betan)

      !> Betan(1) is meaningless.
      term = .False.
      do ie=1, OmegaNum
         energy= omega(ie)
         dos_k_omega(ik, ie)=continued_fraction(Alpha,Betan,energy,Fermi_broadening,NumLczVectors_out, term)
      enddo
   enddo

#if defined (MPI)
   call mpi_allreduce(dos_k_omega, dos_k_omega_mpi, size(dos_k_omega), &
      mpi_dp, mpi_sum, mpi_cmw, ierr)
#else
   dos_k_omega_mpi= dos_k_omega
#endif

   outfileindex= outfileindex+ 1
   if (cpuid.eq.0)then
      open (unit=outfileindex, file='ekbulklcz.dat')
      do ik=1, nk3_band
         do ie=1, omeganum
            write(outfileindex, '(30f16.8)')K3len(ik)*Angstrom2atomic, omega(ie)/eV2Hartree, dos_k_omega_mpi(ik, ie)
         enddo
         write(outfileindex, *) ' '
      enddo
      close(outfileindex)
      write(stdout,*)'calculate Landau level spectrum in k-E mode successfully'
   endif

   outfileindex= outfileindex+ 1
   if (cpuid.eq.0)then
      open (unit=outfileindex, file='ekbulklcz.gnu')
      write(outfileindex,*)     '#set terminal  postscript enhanced color font ",24"'
      write(outfileindex,*)     'set terminal  pngcairo'
      write(outfileindex,*)     "set output 'LandauB_dos.png'"
      write(outfileindex,*)     'set pm3d'
      write(outfileindex,*)     '#set isosamples 50,50'
      write(outfileindex,*)     'set view map'
      write(outfileindex,*)     'unset ztics'
      write(outfileindex,*)     'unset surface'
      write(outfileindex,*)     'unset key'
      write(outfileindex, '(a)')'#set xtics font ",24"'
      write(outfileindex, '(a)')'#set ytics font ",24"'
      write(outfileindex, '(a)')'#set xlabel font ",24"'
      write(outfileindex, '(a)')'set xlabel "Phi per cell"'
      write(outfileindex,*) '#set xlabel "{/Symbol F}/{/Symbol F}_0"'
      write(outfileindex, '(a)')'set ylabel "Energy (eV)"'
      write(outfileindex, '(a, i6,a)')'set title "Hofstadter butterfly with Nq=', Nq, '" font ",40"'
      write(outfileindex,*)     'set yrange [',OmegaMin/eV2Hartree,':',OmegaMax/eV2Hartree,']'
      write(outfileindex,*)     'set xrange [ 0.0000 :', maxval(K3len*Angstrom2atomic) ,']'
      write(outfileindex,*)     "splot 'ekbulklcz.dat' u 1:2:(log($3)) w pm3d"

   endif

#if defined (MPI)
   call mpi_barrier(mpi_cmw, ierr)
#endif
   deallocate(InitialVector, Alpha, Betan)
   deallocate(icsr, jcsr, acsr, iwk)

   return
end subroutine bulkbandk_dos_lanczos


subroutine bulk_dos_lanczos
   !> calculate density of state using Lanczos method only for sparse_hr format
   !> The DOS is integrated in the KCUBE_BULK
   use prec
   use sparse
   use wmpi
   use para, only : Num_Wann, Bx, By, zi, pi, Fermi_broadening, &
      OmegaNum, OmegaMin, OmegaMax, nk3_band, Magp, stdout, nnzmax_input, &
      outfileindex,splen,Is_Sparse_Hr,ijmax,NumLCZVecs,Nk1,Nk2,Nk3,&
      K3D_start_cube,K3D_vec1_cube,K3D_vec2_cube,K3D_vec3_cube, &
      NumRandomConfs, Omega_array, eV2Hartree
   implicit none

   !> energy interval
   !> OmegaNum is defined in the module.f90 and read from the input.dat or wt.in
   real, allocatable :: omega(:)

   !> spectrum calculated
   real(dp), allocatable :: dos_k_omega(:, :), dos_k_omega_mpi(:, :),dos_lcz(:, :)

   integer :: NumberofEta, ieta
   integer :: i, Mdim, ie, ik, ierr, it
   integer :: nnzmax, nnz, NumLczVectors, NumLczVectors_out
   real(dp) :: k3(3), thetaj, energy, eta0
   real(dp), allocatable :: eta_array(:)

   !> left vector and right vector
   complex(dp), allocatable :: InitialVector(:)

   !> Lanczos hamiltonian
   complex(dp), allocatable :: Alpha(:), Betan(:)

   integer, allocatable :: icsr(:), jcsr(:), iwk(:)
   complex(dp), allocatable :: acsr(:)
   complex(dp) :: zdotc, norm
   real(dp) :: continued_fraction
   logical :: term

   integer :: knv3,ikx,iky,ikz
   real(dp) :: dk3, time_end, time_start

   NumberofEta=9
   Mdim= Num_Wann
   knv3= Nk1*Nk2*Nk3

   if (nnzmax_input<0)then
      nnzmax= Num_wann*(2*ijmax+1+2)*Mdim
      if (Is_Sparse_Hr) then
         nnzmax=splen
      else
         nnzmax= Num_wann* Num_wann
      endif
   else
      nnzmax= nnzmax_input
   endif
 
   !> need to be checked
   !if(Is_Sparse_Hr) nnzmax=Nq*splen

   NumLczVectors = NumLCZVecs
   NumLczVectors_out= NumLczVectors
   if (NumLczVectors>Mdim) then
      NumLczVectors = Mdim-1
   endif

   allocate(eta_array(NumberofEta))
   allocate(omega(OmegaNum),dos_lcz(OmegaNum, NumberofEta))
   allocate(dos_k_omega(OmegaNum, NumberofEta), dos_k_omega_mpi(OmegaNum, NumberofEta))
   dos_k_omega= 0d0; dos_k_omega_mpi= 0d0; dos_lcz= 0d0

   !> energy
   omega= Omega_array

   eta_array=(/0.1d0, 0.2d0, 0.4d0, 0.8d0, 1.0d0, 2d0, 4d0, 8d0, 10d0/)
   eta_array= eta_array*Fermi_broadening

   allocate(InitialVector(Mdim))
   allocate(Alpha(NumLczVectors), Betan(NumLczVectors))
   allocate(icsr(nnzmax), jcsr(nnzmax), acsr(nnzmax))
   allocate(iwk (Mdim+1))
   dk3= 1d0/dble(knv3)

   if (cpuid==0) write(stdout, '(a, 2i10)') '>>> LanczosDos_calc starting'
   if (cpuid==0) then
      write(stdout, '(2X, a, i10)')" Number of Lanczos vector we used is ", NumLczVectors
   endif

   time_start= 0d0
   time_end= 0d0
   do ik=1+cpuid, knv3, num_cpu
      do it= 1, NumRandomConfs
         Alpha= 0d0; Betan= 0d0
         icsr=0; jcsr=0;acsr=0d0;iwk=0d0
         do i=1, Mdim
            call random_number(harvest= thetaj)
            InitialVector(i)= exp(zi*2d0*pi*thetaj)
            !InitialVector(i)= 1d0
         enddo
         norm= zdotc(Mdim, InitialVector, 1, InitialVector, 1)

         InitialVector= InitialVector/dsqrt(dble(norm))


         if (cpuid.eq.0.and. mod(ik/num_cpu, 100).eq.0) &
            write(stdout, '(a, i18, "/", i18, a, f10.3, "s")') 'ik/knv3', &
            ik, knv3, ' time left', (knv3-ik)*(time_end-time_start)/num_cpu*NumRandomConfs
         call now(time_start)


         !> Get k coordinates
         ikx= (ik-1)/(nk2*nk3)+1
         iky= ((ik-1-(ikx-1)*Nk2*Nk3)/nk3)+1
         ikz= (ik-(iky-1)*Nk3- (ikx-1)*Nk2*Nk3)
         k3= K3D_start_cube+ K3D_vec1_cube*(ikx-1)/dble(nk1)  &
            + K3D_vec2_cube*(iky-1)/dble(nk2)  &
            + K3D_vec3_cube*(ikz-1)/dble(nk3)
         nnz= nnzmax

         !* Get the Hamiltonian
         if(Is_Sparse_Hr) then
            call ham_bulk_coo_sparsehr(k3, acsr,jcsr,icsr)
         else
            call ham_bulk_coo_densehr(k3,nnzmax, nnz,acsr,icsr,jcsr)
         end if

         !> transform coo format to csr format
         call ConvertCooToCsr(Mdim, nnz, acsr, icsr, jcsr, iwk)
         call csr_sort_indices(Mdim, nnz, icsr, jcsr, acsr)

         !> eleminate the same entries in the sparse matrix
         call csr_sum_duplicates(Mdim, nnz, icsr, jcsr, acsr)
         call csr_sort_indices(Mdim, nnz, icsr, jcsr, acsr)

         !* doing lanczos procedure in order to get alpha, beta
         call lanczos_seqsparse_cpu_z(NumLczVectors, NumLczVectors_out, Mdim, nnz, icsr, jcsr, acsr, InitialVector, Alpha, Betan)

         !> Betan(1) is meaningless.
         term = .False.
         do ie=1, OmegaNum
            energy= omega(ie)
            do ieta= 1, NumberofEta
               eta0= eta_array(ieta)
               dos_k_omega(ie, ieta)= dos_k_omega(ie, ieta)+ continued_fraction(Alpha,Betan,energy,eta0,NumLczVectors_out, term)
            enddo
         enddo
         call now(time_end)
      enddo ! sweep over different random initial Lanczos vector
   enddo ! sweep over k points

#if defined (MPI)
   call mpi_allreduce(dos_k_omega, dos_k_omega_mpi, size(dos_k_omega), &
      mpi_dp, mpi_sum, mpi_cmw, ierr)
#else
   dos_k_omega_mpi= dos_k_omega
#endif

   outfileindex= outfileindex+ 1
   if (cpuid.eq.0)then
      open (unit=outfileindex, file='bulkdoslcz.dat')
      write(outfileindex, *)'# Density of state of bulk system'
      write(outfileindex, '(a, 3i6)')'# Nk1*Nk2*Nk3 :', Nk1, Nk2, Nk3
      write(outfileindex, '(a16,a)')'# E(eV)', 'DOS(E) (states/eV/unit cell)'
      write(outfileindex, '(a16, 90f16.6)')'# E(eV)', eta_array
      dos_lcz= dos_k_omega_mpi/knv3/pi*Num_wann/NumRandomConfs
      do ie=1, omeganum
         write(outfileindex, '(20f16.8)') omega(ie)/eV2Hartree, dos_lcz(ie, :)
      enddo
      write(outfileindex, *) ' '

      close(outfileindex)
      write(stdout,*)'<<< calculate dos with Lanczos successfully'
   endif

   outfileindex= outfileindex+ 1
   if (cpuid.eq.0)then
      open (unit=outfileindex, file='bulkdoslcz.gnu')
      write(outfileindex,*)     '#set terminal  postscript enhanced color font ",24"'
      write(outfileindex,*)     '#set terminal  pngcairo enhanced color font ",40"'
      write(outfileindex,*)     'set terminal pdf enhanced color font ",20" size 4,4 '
      write(outfileindex,*)     "set output 'Lanczos_dos.pdf'"
      write(outfileindex,*)     'set style data linespoints'
      write(outfileindex,*)     '#unset ztics'
      write(outfileindex,*)     'set key box opaque samplen 0.5 textcolor variable font ",12"'
      write(outfileindex,*)     'set key title "Broadening" font ",12"'
      write(outfileindex, '(a)')'#set xtics font ",24"'
      write(outfileindex, '(a)')'#set ytics font ",24"'
      write(outfileindex, '(a)')'#set xlabel font ",24"'
      write(outfileindex, '(a)')'set xlabel "Energy (eV)"'
      write(outfileindex, '(a)')'set ylabel "Dos (states/eV/unit cell)"'
      write(outfileindex,*)     'set xrange [',OmegaMin/eV2Hartree,':',OmegaMax/eV2Hartree,']'
      write(outfileindex,*)     'set yrange [ 0.0000 :', maxval(dos_lcz)+0.2 ,']'
      write(outfileindex,'(a,f6.1,a)') "plot 'bulkdoslcz.dat' u 1:2 w lp ps 0.2 lw 1.0 title ' ",&
         eta_array(1)*1000," meV', \"
      do ieta=2, NumberofEta-1
         write(outfileindex,'(a,i1,a,f6.1,a)') "'' u 1:",ieta+1," w lp ps 0.2 lw 1.0 title ' ",&
            eta_array(ieta)*1000," meV', \"
      enddo
      write(outfileindex,'(a,f6.1,a)') "'' u 1:10 w lp ps 0.2 lw 1.0 title ' ",&
         eta_array(9)*1000," meV'"

   endif

#if defined (MPI)
   call mpi_barrier(mpi_cmw, ierr)
#endif
   deallocate(InitialVector, Alpha, Betan)
   deallocate(icsr, jcsr, acsr, iwk)

   return
end subroutine




subroutine SeqLanczosDOS
   use prec
   use sparse
   use wmpi
   use para, only : Magq, Num_Wann, Bx, By, zi, pi, Fermi_broadening, &
      eV2Hartree, OmegaNum, OmegaMin, OmegaMax, Omega_array
   implicit none

   integer :: i, Mdim, Nq, ie
   integer :: nnzmax, nnz, NumLczVectors, NumLczVectors_out
   real(dp) :: k3(3), thetaj, B0

   !> diagonal elements of Lanczos projective Hamiltonian
   real(dp), allocatable :: w(:), DOS(:)
   real(dp), allocatable :: omega(:)
   complex(dp), allocatable :: Hlcz(:,:)

   !> left vector and right vector
   complex(dp), allocatable :: InitialVector(:)

   !> Lanczos hamiltonian
   complex(dp), allocatable :: Alpha(:), Betan(:)

   integer, allocatable :: icsr(:), jcsr(:), iwk(:)
   complex(dp), allocatable :: acsr(:)
   complex(dp) :: zdotc, norm, theta
   real(dp) :: continued_fraction
   logical :: term

   Nq= Magq
   Mdim= Num_Wann*Magq
   nnzmax= Num_wann*13*Mdim
   NumLczVectors = NumLczVectors
   NumLczVectors_out= NumLczVectors


   allocate(omega(OmegaNum), DOS(OmegaNum))

   !> energy
   omega= Omega_array


   !> deal with the magnetic field
   !> first transform the Bx By into B*Cos\theta, B*Sin\theta
   if (abs(By)<1e-8) then
      theta= 0d0
   elseif (By>0) then
      theta = atan(Bx/By)
   else
      theta = atan(Bx/By)+ pi
   endif
   if (abs(theta)<1e-8 .and. Bx<0 ) theta= theta+ pi

   !> get the shortest bond in the home unit cell
   B0= 2d0*pi/ Nq
   B0= abs(B0)
   Bx= B0* Cos(theta)
   By= B0* Sin(theta)

   allocate(InitialVector(Mdim))
   allocate(Alpha(NumLczVectors), Betan(NumLczVectors))
   allocate(icsr(nnzmax), jcsr(nnzmax), acsr(nnzmax))
   allocate(iwk (Mdim+1))
   Alpha= 0d0
   Betan= 0d0

   do i=1, Mdim
      call random_number(harvest= thetaj)
      InitialVector(i)= exp(zi*2d0*pi*thetaj)
      !InitialVector(i)= 1d0
   enddo
   norm= zdotc(Mdim, InitialVector, 1, InitialVector, 1)

   InitialVector= InitialVector/dsqrt(dble(norm))

   !* Get the Hamiltonian
   k3= 0d0
   nnz= nnzmax
   call ham_3Dlandau_sparse1(nnz, Mdim, Nq, k3, acsr,jcsr,icsr)
   acsr=acsr/eV2Hartree

   !> transform coo format to csr format
   call ConvertCooToCsr(Mdim, nnz, acsr, icsr, jcsr, iwk)
   call csr_sort_indices(Mdim, nnz, icsr, jcsr, acsr)

   !> eleminate the same entries in the sparse matrix
   call csr_sum_duplicates(Mdim, nnz, icsr, jcsr, acsr)
   call csr_sort_indices(Mdim, nnz, icsr, jcsr, acsr)

   !* doing lanczos procedure in order to get alpha, beta
   call lanczos_seqsparse_cpu_z(NumLczVectors, NumLczVectors_out, Mdim, nnz, icsr, jcsr, acsr, InitialVector, Alpha, Betan)

   !> Betan(1) is meaningless.
   term = .True.
   do ie=1, OmegaNum
      dos(ie)=continued_fraction(Alpha,Betan,omega(ie),Fermi_broadening,NumLczVectors_out, term)
      write(10001, '(100f16.6)') omega(ie), dos(ie)
   enddo

   if (allocated(Hlcz)) deallocate(Hlcz)
   if (allocated(w)) deallocate(w)

   return

end subroutine SeqLanczosDOS

function continued_fraction(alphan,betan,omega,eta,NumLczVectors, term)
   !> Foa Torres L.E.F., Roche S., Charlier J.-C. - Introduction to Graphene-Based Nanomaterials_ From Electronic Structure to
   !Quantum Transport (2014, CUP).pdf
   use prec
   use para, only : stdout
   implicit none

   ! computes the continued fraction.
   integer, intent(in)  :: NumLczVectors
   complex(dp), intent(in) :: alphan(NumLczVectors)
   complex(dp), intent(in) :: betan(NumLczVectors)
   real(dp), intent(in) :: eta
   real(dp), intent(in) :: omega
   logical, intent(in) :: term

   real(dp) :: continued_fraction

   integer :: i, p,q
   complex(dp) :: res ,lastterm ! intermediate variable
   real(dp) :: aa, bb

   q=NumLczVectors/2
   IF (term) THEN
      aa=0.0
      bb=0.0
      DO p=1, q
         aa=aa+alphan(NumLczVectors-p)
         bb=bb+betan(NumLczVectors-p+1)
      ENDDO
      aa=aa/q
      bb=bb/q
      res=lastterm(aa-omega,bb*bb,eta)
   ELSE
      res = alphan(NumLczVectors)+ cmplx(-omega, eta, kind=dp)
   ENDIF
   DO i = 1, NumLczVectors -1
      res = alphan(NumLczVectors-i)+ cmplx(-omega, -eta, kind=dp) &
         -betan(NumLczVectors-i+1)*betan(NumLczVectors-i+1)/res
   ENDDO
   continued_fraction = AIMAG(1/res)
   return
end function continued_fraction

function lastterm(a,b,g)
   use prec, only : dp
   implicit none
   real(dp)    :: a, b, g, y1, y2, z1, z2, r
   complex(dp) :: lastterm

   y1 = a*a - g*g - 4*b
   y2 = -2*a*g
   r  = 0.5*sqrt(y1*y1 + y2*y2)

   if (g<0) then
      z1 =  a/2 + 0.5*sign(sqrt(y1/2 + r),y2)
      z2 = -g/2 + 0.5*sqrt(-y1/2 + r)
   else
      z1 =  a/2 - 0.5*sign(sqrt(y1/2 + r),y2)
      z2 = -g/2 - 0.5*sqrt(-y1/2+r)
   endif

   lastterm = cmplx(z1,z2,kind=dp)

end function lastterm

!> check whether the lanczos procedure is converged
subroutine lanczos_converged(alphan,betan, NumLczVectors, current_lcz_vectors, dos_old, &
                          error, term, converged)
  use para, only : dp, OmegaNum, Omega_array, Fermi_broadening, eps6
  implicit none

  !>> Lanczos elements
  integer, intent(in) :: current_lcz_vectors, NumLczVectors
  complex(dp), intent(in) :: alphan(NumLczVectors),betan(NumLczVectors)

  !> current number of Lanczos vectors
  logical, intent(in) :: term
  logical, intent(out) :: converged
  real(dp), intent(out) :: error
  real(dp), intent(inout) :: dos_old(OmegaNum)

  real(dp) :: tmp,area,e,diff
  real(dp), external :: continued_fraction
  integer  :: ie
  
  diff = 0d0
  area = 0d0
  do ie=1, OmegaNum
     e=Omega_array(ie)
     tmp = continued_fraction(alphan, betan, e, Fermi_broadening, current_lcz_vectors, term)
     diff = diff + abs(dos_old(ie)-tmp)
     area = area + abs(tmp)
     dos_old(ie) = tmp
  enddo
  diff = diff / area
  error= diff
  if(diff < eps6) then
     converged = .true.
  else
     converged = .false.
  endif

  return
end subroutine lanczos_converged
 


subroutine ParLanczosDOS
   use prec
   use sparse
   use wmpi
   use para
   implicit none

   integer :: i, Mdim
   integer :: NumLczVectors
   type (WTParCSR) :: Ham

   !> diagonal elements of Lanczos projective Hamiltonian
   complex(dp), allocatable :: Hlcz(:,:)
   real(dp), allocatable :: w(:)

   !> Lanczos hamiltonian
   complex(dp), allocatable :: Alphan(:), Betan(:)

   !> left vector and right vector
   type(WTParVec), pointer :: InitialVector

   Mdim= Num_Wann*Nslab1*Nslab2
   InitialVector=> Null()
   allocate(InitialVector)
   allocate(Alphan(NumLczVectors), Betan(NumLczVectors))
#if defined (MPI)
   call WTParVecCreate(mpi_cmw, Mdim, InitialVector)
#endif
   call WTParVecInitialize(InitialVector)

   !call InitLanczosVector(InitialVector)

   NumLczVectors = 100

   !* Get the Hamiltonian

   !* doing lanczos procedure in order to get alpha, beta
   call lanczos_parsparse_cpu_z(NumLczVectors, Mdim, Ham, InitialVector, Alphan, Betan)

   allocate(Hlcz(NumLczVectors, NumLczVectors))
   allocate(W(NumLczVectors))

   Hlcz= 0d0

   do i=1, NumLczVectors-1
      Hlcz(i,i)= dble(Alphan(i))
      Hlcz(i,i+1)= dble(betan(i+1))
      Hlcz(i+1,i)= dble(betan(i+1))
   enddo
   i=NumLczVectors
   Hlcz(i,i)= dble(Alphan(i))


   !* diagonalize Lanczos projective hamiltonian
   call eigensystem_c('N', 'U', NumLczVectors, Hlcz, W)

   !* free Hdiag, Hsubdiag
   if (allocated(Hlcz)) deallocate(Hlcz)
   if (allocated(w)) deallocate(w)

   return

end subroutine ParLanczosDOS

 !> Construct Hamiltonian for subspace of NumElectrons' Electrons
 !> provide a tempalette to construct a parallel CSR matrix
subroutine WTHamConstruct()

   use prec
   use sparse
   use wmpi
   use para, only : cpuid, num_cpu, Num_Wann, Nslab2, Nslab1, stdout
   implicit none

   integer(li) :: Mdim, NNZPerRow
   integer(li) :: RowStart, RowEnd, NumRows, NumCols, maxnnz, nnz

   !> these four arrays are used for construct hamiltonian in coo format
   !> and transform the coo format to csr format
   integer(li), allocatable, target :: iwk (:), icoo(:), jcoo(:)
   complex(dp), allocatable, target :: acoo(:)

   type(WTCSR), pointer :: HCSR
   type(WTParCSR), pointer :: HParCSR

   integer(li) :: i, j, ierr

   real(dp) :: time0, time1, time2, timegetham, memory

   call now(time0)

   !> for slab system
   Mdim= Num_Wann*Nslab1*Nslab2

   HCSR=> Null()

   if (.not.associated(HParCSR)) allocate(HParCSR)
   call WTParCSRMatrixCreate(Mdim, HParCSR)

   call WTGenerateLocalPartition(Mdim, num_cpu, cpuid, RowStart, RowEnd)

   NumRows= RowEnd- RowStart+ 1
   NumCols= Mdim
   NNZPerRow= 13
   if (NumCols<=num_cpu) then
      maxnnz= numcols
   else
      maxnnz= numcols* NNZPerRow/num_cpu
   endif

   !>  calculate memory required by ia, ja, H
   memory= ( NumRows*8.0      & ! iwk
      + maxnnz*8.0       & ! icoo
      + maxnnz*8.0       & ! jcoo
      + maxnnz*16.0      & ! acoo
      )/1024.0/1024.0    ! MB

   if (cpuid==0) then
      write(stdout,102)'Memory required by constructing Hamiltonian :', memory, ' MB'
   endif

102 format(2x,a,f9.1,a)

   allocate(iwk (NumRows+1), stat= ierr)
   allocate(icoo(maxnnz), stat= ierr)
   allocate(jcoo(maxnnz), stat= ierr)
   allocate(acoo(maxnnz), stat= ierr)

   call now(time1)
   !> The row index of the Hamiltonian is ranging from RowStart to RowEnd
   !> The coloum index of the Hamiltonian is ranging from 1 to Mdim
   !> Stored in COO sparse matrix format
   call WTGetHam(Mdim, RowStart, RowEnd, maxnnz, nnz, acoo, icoo, jcoo)
   call now(time2)
   timegetham= time2- time1
#if defined (MPI)
   call mpi_reduce(timegetham, time1, 1, mpi_dp, mpi_max, 0, mpi_cmw, ierr)
#endif

   call print_time_cost(time1, 0d0, 'WTGetHam')

   !$OMP PARALLEL DO PRIVATE(i)
   do i=1, maxnnz
      icoo(i)= icoo(i)- RowStart+ 1 !< make a shift
   enddo
   !$OMP END PARALLEL DO

   !> convert coo storage format to csr storage format
   call now(time1)
   call ConvertCooToCsr(NumRows, nnz, acoo, icoo, jcoo, iwk)
   call csr_sort_indices(Mdim, nnz, icoo, jcoo, acoo)
   call csr_sum_duplicates(Mdim, nnz, icoo, jcoo, acoo)
   call csr_sort_indices(Mdim, nnz, icoo, jcoo, acoo)

   call now(time2)
   timegetham= time2- time1
#if defined (MPI)
   call mpi_reduce(timegetham, time1, 1, mpi_dp, mpi_max, 0, mpi_cmw, ierr)
#endif

   call print_time_cost(time1, 0d0, 'WTGetHam')

   if (.not. associated(HCSR)) allocate(HCSR)
   call WTCSRCreate(NumRows, NumCols, nnz, HCSR)
   HCSR%ia=> icoo(1:NumRows+1)
   HCSR%ja=> jcoo(1:nnz)
   HCSR%a => acoo(1:nnz)

   !> this subroutine is tested to be right on May 13 2014 by QS.Wu
   call now(time1)
   call WTCSRToParCSR(HCSR, HParCSR)
   call now(time2)
   timegetham= time2- time1
#if defined (MPI)
   call mpi_reduce(timegetham, time1, 1, mpi_dp, mpi_max, 0, mpi_cmw, ierr)
#endif

   call print_time_cost(time1, 0d0, 'WTCSRToParCSR')

   call now(time1)
   call WTGenSendRecv(HParCSR)
   deallocate(icoo, jcoo, acoo, iwk)
   call now(time2)
   timegetham= time2- time1
#if defined (MPI)
   call mpi_reduce(timegetham, time1, 1, mpi_dp, mpi_max, 0, mpi_cmw, ierr)
#endif

   call print_time_cost(time1, 0d0, 'WTGenSendRecv')

   !* because pointer ia, ja, a are point to icoo, jcoo, acoo respectively,
   !* icoo, jcoo, acoo are deallocated, so we need to nullify ia, ja, a,
   !* then destroy CSR matrix
   nullify(HCSR%ia, HCSR%ja, HCSR%a)
   call WTCSRDestroy(HCSR)
   nullify(HCSR)

   call now(time2)
   call print_time_cost(time2, time0, 'WTHamConstruct')

   return
end subroutine WTHamConstruct


subroutine WTGetHam(Mdim, RowStart, RowEnd, maxnnz, nnz, acoo, icoo, jcoo)
   use para, only : dp
   implicit none

   integer, intent(in) :: Mdim
   integer, intent(in) :: RowStart
   integer, intent(in) :: RowEnd
   integer, intent(in) :: maxnnz
   integer, intent(in) :: nnz
   integer, intent(in) :: icoo(maxnnz)
   integer, intent(in) :: jcoo(maxnnz)
   complex(dp), intent(in) :: acoo(maxnnz)


   return
end subroutine WTGetHam


