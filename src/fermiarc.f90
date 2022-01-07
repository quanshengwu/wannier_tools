  subroutine SurfaceDOSkk
     ! This subroutine calculates surface states using
     ! iterative Green's function method  ---  J.Phys.F.Met.Phys.15(1985)851-858
     !> 1. Calculate surface local density of states in k-k mode
     !> 2. Calculate surface spin-texture in k-k mode
     !> 3. Calculation joint density of state which is also called the quasi-partical interference
     !> refs :http://science.sciencemag.org/content/sci/suppl/2016/03/09/351.6278.1184.DC1/Inoue.SM.pdf
     !> in this version, we also calculate the spin density jdos
     !  by Quan Sheng Wu on 4/20/2010
     !  mpi version      4/21/2010
     !  Add spintexture by QS.Wu on July/20/2010  10:26:13
     !
     !  Modified by QS.Wu on June/4/2016 in Jouvence Montreal Canada 
     !  Modified by Changming Yue on 12/4/2017 in IOP Beijing
     !  Merged fermiarc, spintexture and QPI on 10/09/2018
     ! License: GPL v3

     use wmpi
     use para
     implicit none

     integer :: arclfile, arcrfile, arcbulkfile
     integer :: arcljfile, arcrjfile
     integer :: arcljsfile, arcrjsfile
     integer :: spindosrfile, spindoslfile

     ! general loop index
     integer :: i, j, io, ierr, nkx, nky, Nwann

     ! kpoint loop index
     integer :: ikp, ik1, ik2, iq, Nk1_half, Nk2_half
     integer :: imin1, imax1, imin2, imax2, iq1, iq2

     real(dp) :: dos_l_max, dos_r_max
     real(dp) :: time_start, time_end
     real(dp) :: omega, k1min_shape, k1max_shape, k2min_shape, k2max_shape

     real(dp) :: k(2)
     real(dp) :: s0(3), s1(3)
     real(dp) :: sx_bulk, sy_bulk, sz_bulk

     integer , allocatable :: ik12(:,:)
     real(dp), allocatable :: k12(:,:), k12_shape(:,:)

     real(dp), allocatable :: dos_l(:,:), dos_l_mpi(:,:)
     real(dp), allocatable :: dos_r(:,:), dos_r_mpi(:,:)
     real(dp), allocatable :: dos_bulk(:,:), dos_bulk_mpi(:,:)
     real(dp), allocatable :: jdos_l(:), jdos_l_mpi(:)
     real(dp), allocatable :: jdos_r(:), jdos_r_mpi(:)
     real(dp), allocatable :: jsdos_l(:), jsdos_l_mpi(:)
     real(dp), allocatable :: jsdos_r(:), jsdos_r_mpi(:)
     real(dp), allocatable :: dos_l_only(:, :), dos_r_only(:, :)
     real(dp), allocatable :: sx_l(:, :), sy_l(:, :), sz_l(:, :)
     real(dp), allocatable :: sx_l_mpi(:, :), sy_l_mpi(:, :), sz_l_mpi(:, :)
     real(dp), allocatable :: sx_r(:, :), sy_r(:, :), sz_r(:, :)
     real(dp), allocatable :: sx_r_mpi(:, :), sy_r_mpi(:, :), sz_r_mpi(:, :)

     ! spin operator matrix spin_sigma_x,spin_sigma_y in spin_sigma_z representation
     complex(Dp),allocatable :: spin_sigma_x(:,:), spin_sigma_y(:,:), spin_sigma_z(:,:)
     complex(Dp),allocatable :: ctemp(:,:), ones(:,:)
     complex(dp), allocatable :: GLL(:,:), GRR(:,:), GB(:,:)
     complex(dp), allocatable :: H00(:,:), H01(:,:)

     !> Nk1 and Nk2 should be odd number so that the center of the kslice is (0,0)
     !> if you want to calculate the QPI
     nkx=Nk1; nky=Nk2
     if (mod(Nk1, 2)==0) nkx= Nk1+1
     if (mod(Nk2, 2)==0) nky= Nk2+1

     Nk1_half= (nkx-1)/2
     Nk2_half= (nky-1)/2
     Nwann= Num_wann/2

     allocate( ik12(2, nkx*nky), k12(2, nkx*nky), k12_shape(2, nkx*nky))
     allocate( dos_l(nkx,nky), dos_l_mpi(nkx,nky))
     allocate( dos_r(nkx,nky), dos_r_mpi(nkx,nky))
     allocate( dos_r_only(nkx,nky), dos_l_only(nkx,nky))
     allocate( dos_bulk(nkx,nky), dos_bulk_mpi(nkx,nky))
     allocate( jdos_l(nkx*nky), jdos_l_mpi(nkx*nky))
     allocate( jdos_r(nkx*nky), jdos_r_mpi(nkx*nky))
     allocate( GLL(ndim,ndim), GRR(ndim,ndim), GB(ndim, ndim))
     allocate( ctemp(ndim,ndim))
     if (SOC>0 .and. (SlabSpintexture_calc.or.SlabQPI_kplane_calc)) then
        allocate( jsdos_l(nkx*nky), jsdos_l_mpi(nkx*nky))
        allocate( jsdos_r(nkx*nky), jsdos_r_mpi(nkx*nky))
        allocate( sx_l(nkx,nky), sx_l_mpi(nkx,nky))
        allocate( sy_l(nkx,nky), sy_l_mpi(nkx,nky))
        allocate( sz_l(nkx,nky), sz_l_mpi(nkx,nky))
        allocate( sx_r(nkx,nky), sx_r_mpi(nkx,nky))
        allocate( sy_r(nkx,nky), sy_r_mpi(nkx,nky))
        allocate( sz_r(nkx,nky), sz_r_mpi(nkx,nky))
        allocate(spin_sigma_x(ndim,ndim), spin_sigma_y(ndim,ndim), spin_sigma_z(ndim,ndim))
     endif

     ik12=0
     k12=0d0; k12_shape=0d0
     dos_l=0d0; dos_l_mpi=0d0
     dos_r=0d0; dos_r_mpi=0d0
     jdos_l=0d0; jdos_l_mpi=1d-12
     jdos_r=0d0; jdos_r_mpi=1d-12
     if (SOC>0 .and. (SlabSpintexture_calc.or.SlabQPI_kplane_calc)) then
        spin_sigma_x=0.0d0; spin_sigma_y=0.0d0; spin_sigma_z=0.0d0
        jsdos_l=0d0; jsdos_l_mpi=1d-12
        jsdos_r=0d0; jsdos_r_mpi=1d-12
        sx_l=0d0; sy_l=0d0; sz_l=0d0
        sx_l_mpi=0d0; sy_l_mpi=0d0; sz_l_mpi=0d0
        sx_r=0d0; sy_r=0d0; sz_r=0d0
        sx_r_mpi=0d0; sy_r_mpi=0d0; sz_r_mpi=0d0
     endif

     ikp=0
     do i= 1, nkx
        do j= 1, nky
           ikp=ikp+1
           ik12(1, ikp)= i
           ik12(2, ikp)= j
           k12(:, ikp)=K2D_start+ (i-1)*K2D_vec1/dble(nkx-1) &
                      + (j-1)*K2D_vec2/dble(nky-1)
           k12_shape(:, ikp)= k12(1, ikp)* Ka2+ k12(2, ikp)* Kb2
        enddo
     enddo

     k1min_shape= minval(k12_shape(1, :))
     k2min_shape= minval(k12_shape(2, :))
     k1max_shape= maxval(k12_shape(1, :))
     k2max_shape= maxval(k12_shape(2, :))

     allocate(H00(Ndim, Ndim))
     allocate(H01(Ndim, Ndim))
     allocate(ones(Ndim, Ndim))
     GLL= 0d0; GRR= 0d0
     H00= 0d0; H01= 0d0; ones= 0d0

     do i=1,Ndim
        ones(i,i)=1.0d0
     enddo

     if (SOC>0 .and. (SlabSpintexture_calc.or.SlabQPI_kplane_calc)) then
        if (index(Particle,'phonon')/=0) then
           stop "ERROR: we don't support spintexture calculation for phonon system"
        endif
        Nwann= Num_wann/2
        !> spin operator matrix
        !> this part is package dependent. 
        if (index( Package, 'VASP')/=0.or. index( Package, 'Wien2k')/=0 &
           .or. index( Package, 'Abinit')/=0.or. index( Package, 'openmx')/=0) then
           do i=1, Np
              do j=1, Nwann
                 spin_sigma_x(Num_wann*(i-1)+j, Num_wann*(i-1)+Nwann+j)=1.0d0
                 spin_sigma_x(Num_wann*(i-1)+j+Nwann, Num_wann*(i-1)+j)=1.0d0
                 spin_sigma_y(Num_wann*(i-1)+j, Num_wann*(i-1)+Nwann+j)=-zi
                 spin_sigma_y(Num_wann*(i-1)+j+Nwann, Num_wann*(i-1)+j)=zi
                 spin_sigma_z(Num_wann*(i-1)+j, Num_wann*(i-1)+j)= 1d0
                 spin_sigma_z(Num_wann*(i-1)+j+Nwann, Num_wann*(i-1)+j+Nwann)=-1d0
              enddo
           enddo
        elseif (index( Package, 'QE')/=0.or.index( Package, 'quantumespresso')/=0 &
             .or.index( Package, 'quantum-espresso')/=0.or.index( Package, 'pwscf')/=0) then
           do i=1, Np
              do j=1, Nwann
                 spin_sigma_x(Num_wann*(i-1)+(2*j-1), Num_wann*(i-1)+2*j)=1.0d0
                 spin_sigma_x(Num_wann*(i-1)+2*j, Num_wann*(i-1)+(2*j-1))=1.0d0
                 spin_sigma_y(Num_wann*(i-1)+(2*j-1), Num_wann*(i-1)+2*j)=-zi
                 spin_sigma_y(Num_wann*(i-1)+2*j, Num_wann*(i-1)+(2*j-1))=zi
                 spin_sigma_z(Num_wann*(i-1)+(2*j-1), Num_wann*(i-1)+(2*j-1))=1.0d0
                 spin_sigma_z(Num_wann*(i-1)+2*j, Num_wann*(i-1)+2*j)=-1.0d0
              enddo
           enddo
        else
           if (cpuid.eq.0) write(stdout, *)'Error: please report your software generating tight binding and wannier90.wout to me'
           if (cpuid.eq.0) write(stdout, *)'wuquansheng@gmail.com'
           stop 'Error: please report your software and wannier90.wout to wuquansheng@gmail.com'
        endif
     endif

     omega = E_arc
     eta= eta_arc

     !> deal with phonon system
     !> for phonon system, omega should be changed to omega^2
     if (index(Particle,'phonon')/=0) then
        omega= omega*omega
     endif


     time_start= 0d0
     time_end= 0d0
     do ikp= 1+cpuid, nkx*nky, num_cpu
        if (cpuid==0.and. mod(ikp/num_cpu, 100)==0) &
           write(stdout, *) 'Arc, ik ', ikp, 'Nk',nkx*nky, 'time left', &
           (nkx*nky-ikp)*(time_end- time_start)/num_cpu, ' s'
        call now(time_start)
        k(1)= k12(1, ikp)
        k(2)= k12(2, ikp)

        if (index(Particle,'phonon')/=0.and.LOTO_correction) then
           call ham_qlayer2qlayer_LOTO(k,H00,H01)
        else
           call ham_qlayer2qlayer(k,H00,H01)
        endif


        !> calculate surface green function
        ! there are two method to calculate surface green's function
        ! the method in 1985 is better, you can find the ref in the
        ! subroutine
        call surfgreen_1985(omega,GLL,GRR,GB,H00,H01,ones)
        ! call surfgreen_1984(omega,GLL,GRR,H00,H01,ones)

        ik1= ik12(1, ikp)
        ik2= ik12(2, ikp)

        ! calculate spectral function
        do i= 1, NtopOrbitals
           io= TopOrbitals(i)
           dos_l(ik1, ik2)=dos_l(ik1, ik2)- aimag(GLL(io,io))/pi
        enddo ! i
        do i= 1, NBottomOrbitals
           io= Ndim- Num_wann+ BottomOrbitals(i)
           dos_r(ik1, ik2)=dos_r(ik1, ik2)- aimag(GRR(io,io))/pi
        enddo ! i
        do i= 1, Ndim
           dos_bulk(ik1, ik2)=dos_bulk(ik1, ik2)- AIMAG(GB(i,i))/pi
        enddo ! i

        if (SOC>0 .and. (SlabSpintexture_calc.or.SlabQPI_kplane_calc)) then
           !>> calculate spin-resolved bulk spectrum
           sx_bulk= 0d0
           call mat_mul(ndim,GB ,spin_sigma_x,ctemp)
           do i= 1, Ndim
              sx_bulk=sx_bulk- aimag(ctemp(i,i))/pi
           enddo ! i
   
           sy_bulk= 0d0
           call mat_mul(ndim,GB ,spin_sigma_y,ctemp)
           do i= 1, Ndim
              sy_bulk=sy_bulk- aimag(ctemp(i,i))/pi
           enddo ! i
   
           sz_bulk= 0d0
           call mat_mul(ndim,GB ,spin_sigma_z,ctemp)
           do i= 1, Ndim
              sz_bulk=sz_bulk- aimag(ctemp(i,i))/pi
           enddo ! i
   
           !>> calculate spin-resolved surface spectrum
   
           !ctemp=matmul(surfgreen,spin_sigma_x)
           call mat_mul(ndim,GLL,spin_sigma_x,ctemp)
           do i= 1, NtopOrbitals
              io= TopOrbitals(i)
              sx_l(ik1, ik2)=sx_l(ik1, ik2)- aimag(ctemp(io,io))/pi
           enddo ! i
          !sx_l(ik1, ik2)= sx_l(ik1, ik2)- sx_bulk
          !if (sx_l(ik1, ik2)<0) sx_l(ik1, ik2)= eps9
   
           !ctemp=matmul(surfgreen,spin_sigma_y)
           call mat_mul(ndim,GLL,spin_sigma_y,ctemp)
           do i= 1, NtopOrbitals
              io= TopOrbitals(i)
              sy_l(ik1, ik2)=sy_l(ik1, ik2)- aimag(ctemp(io,io))/pi
           enddo ! i
          !sy_l(ik1, ik2)= sy_l(ik1, ik2)- sy_bulk
          !if (sy_l(ik1, ik2)<0) sy_l(ik1, ik2)= eps9
   
   
           !ctemp=matmul(surfgreen,spin_sigma_z)
           call mat_mul(ndim,GLL,spin_sigma_z,ctemp)
           do i= 1, NtopOrbitals
              io= TopOrbitals(i)
              sz_l(ik1, ik2)=sz_l(ik1, ik2)- aimag(ctemp(io,io))/pi
           enddo ! i
          !sz_l(ik1, ik2)= sz_l(ik1, ik2)- sz_bulk
          !if (sz_l(ik1, ik2)<0) sz_l(ik1, ik2)= eps9
   
   
           !ctemp=matmul(surfgreen,spin_sigma_x)
           call mat_mul(ndim,GRR,spin_sigma_x,ctemp)
           do i= 1, NBottomOrbitals
              io= Ndim- Num_wann+ BottomOrbitals(i)
              sx_r(ik1, ik2)=sx_r(ik1, ik2)- aimag(ctemp(io,io))/pi
           enddo ! i
          !sx_r(ik1, ik2)= sx_r(ik1, ik2)- sz_bulk
          !if (sx_r(ik1, ik2)<0) sx_r(ik1, ik2)= eps9
   
   
           !ctemp=matmul(surfgreen,spin_sigma_y)
           call mat_mul(ndim,GRR,spin_sigma_y,ctemp)
           do i= 1, NBottomOrbitals
              io= Ndim- Num_wann+ BottomOrbitals(i)
              sy_r(ik1, ik2)=sy_r(ik1, ik2)- aimag(ctemp(io,io))/pi
           enddo ! i
          !sy_r(ik1, ik2)= sy_r(ik1, ik2)- sz_bulk
          !if (sy_r(ik1, ik2)<0) sy_r(ik1, ik2)= eps9
   
   
           !ctemp=matmul(surfgreen,spin_sigma_z)
           call mat_mul(ndim,GRR,spin_sigma_z,ctemp)
           do i= 1, NBottomOrbitals
              io= Ndim- Num_wann+ BottomOrbitals(i)
              sz_r(ik1, ik2)=sz_r(ik1, ik2)- aimag(ctemp(io,io))/pi
           enddo ! i
          !sz_r(ik1, ik2)= sz_r(ik1, ik2)- sz_bulk
          !if (sz_r(ik1, ik2)<0) sz_r(ik1, ik2)= eps9
   
        endif  !> SOC>0
        call now(time_end)

     enddo

#if defined (MPI)
     !> we don't have to do allreduce operation
     call mpi_allreduce(dos_l, dos_l_mpi, size(dos_l),mpi_double_precision,&
                     mpi_sum, mpi_comm_world, ierr)
     call mpi_allreduce(dos_r, dos_r_mpi, size(dos_r),mpi_double_precision,&
                     mpi_sum, mpi_comm_world, ierr)
     call mpi_allreduce(dos_bulk, dos_bulk_mpi, size(dos_bulk),mpi_double_precision,&
                     mpi_sum, mpi_comm_world, ierr)
     if (SOC>0 .and. (SlabSpintexture_calc.or.SlabQPI_kplane_calc)) then
        call mpi_allreduce(sx_l, sx_l_mpi, size(sx_l),mpi_double_precision,&
                        mpi_sum, mpi_comm_world, ierr)
        call mpi_allreduce(sy_l, sy_l_mpi, size(sy_l),mpi_double_precision,&
                        mpi_sum, mpi_comm_world, ierr)
        call mpi_allreduce(sz_l, sz_l_mpi, size(sz_l),mpi_double_precision,&
                        mpi_sum, mpi_comm_world, ierr)
        call mpi_allreduce(sx_r, sx_r_mpi, size(sx_r),mpi_double_precision,&
                        mpi_sum, mpi_comm_world, ierr)
        call mpi_allreduce(sy_r, sy_r_mpi, size(sy_r),mpi_double_precision,&
                        mpi_sum, mpi_comm_world, ierr)
        call mpi_allreduce(sz_r, sz_r_mpi, size(sz_r),mpi_double_precision,&
                        mpi_sum, mpi_comm_world, ierr)
     endif
#else
     dos_l_mpi= dos_l
     dos_r_mpi= dos_r
     dos_bulk_mpi= dos_bulk
     if (SOC>0 .and. (SlabSpintexture_calc.or.SlabQPI_kplane_calc)) then
        sx_l_mpi= sx_l
        sy_l_mpi= sy_l
        sz_l_mpi= sz_l
        sx_r_mpi= sx_r
        sy_r_mpi= sy_r
        sz_r_mpi= sz_r
     endif
#endif

     dos_l_max= maxval(dos_l_mpi)
     dos_r_max= maxval(dos_r_mpi)

     do ikp=1, nkx*nky
        ik1= ik12(1, ikp)
        ik2= ik12(2, ikp)

        dos_l_only(ik1, ik2)= dos_l_mpi(ik1, ik2)
        if (dos_l_only(ik1, ik2)<dos_l_max/10d0) then
           dos_l_only(ik1, ik2)=eps9
           if (SOC>0 .and. (SlabSpintexture_calc.or.SlabQPI_kplane_calc)) then
              sx_l_mpi(ik1, ik2)=eps9
              sy_l_mpi(ik1, ik2)=eps9
              sz_l_mpi(ik1, ik2)=eps9
           endif
        endif
        dos_r_only(ik1, ik2)= dos_r_mpi(ik1, ik2)
        if (dos_r_only(ik1, ik2)<dos_r_max/10d0) then
           dos_r_only(ik1, ik2)=eps9
           if (SOC>0 .and. (SlabSpintexture_calc.or.SlabQPI_kplane_calc)) then
              sx_r_mpi(ik1, ik2)=eps9
              sy_r_mpi(ik1, ik2)=eps9
              sz_r_mpi(ik1, ik2)=eps9
           endif
        endif
     enddo

     outfileindex= outfileindex+ 100
     arclfile= outfileindex
     outfileindex= outfileindex+ 1
     arcrfile= outfileindex
     outfileindex= outfileindex+ 1
     spindoslfile= outfileindex
     outfileindex= outfileindex+ 1
     spindosrfile= outfileindex
     outfileindex= outfileindex+ 1
     arcbulkfile= outfileindex

     outfileindex= outfileindex+ 1
     arcljfile= outfileindex
     outfileindex= outfileindex+ 1
     arcrjfile= outfileindex
     outfileindex= outfileindex+ 3
     arcljsfile= outfileindex
     outfileindex= outfileindex+ 3
     arcrjsfile= outfileindex
 
     if (cpuid.eq.0)then
        open (unit=arclfile, file='arc.dat_l')
        open (unit=arcrfile, file='arc.dat_r')
        open (unit=arcbulkfile, file='arc.dat_bulk')
        write(arclfile,'(a)')'# Surface density of states of bottom surface'
        write(arclfile,'(a)')'# The coordinate of k is redefined according to the SURFACE card'
        write(arclfile,'(a)')"# x axis is parallel to R1'"
        write(arclfile,'(a)')"# z axis is parallel to R1'xR2'"
        write(arclfile,'(a)')"# y axis is parallel to z x x"
        write(arclfile,'(30a16)')'#kx', 'ky', 'log(dos)'
        write(arcrfile,'(a)')'# Surface density of states of top surface'
        write(arcbulkfile,'(a)')'#Density of states of bulk system after projection onto the surface BZ'
        write(arcrfile,'(a)')'# The coordinate of k is redefined according to the SURFACE card'
        write(arcrfile,'(a)')"# x axis is parallel to R1'"
        write(arcrfile,'(a)')"# z axis is parallel to R1'xR2'"
        write(arcrfile,'(a)')"# y axis is parallel to z x x"
        write(arcrfile,'(30a16)')'#kx', 'ky', 'log(dos)'
        write(arcbulkfile,'(a)')'#Density of states of bulk system after projection onto the surface BZ'
        write(arcbulkfile,'(a)')'# The coordinate of k is redefined according to the SURFACE card'
        write(arcbulkfile,'(a)')"# x axis is parallel to R1'"
        write(arcbulkfile,'(a)')"# z axis is parallel to R1'xR2'"
        write(arcbulkfile,'(a)')"# y axis is parallel to z x x"
        write(arcbulkfile,'(30a16)')'#kx', 'ky', 'log(dos)'
        do ikp=1, nkx*nky
           write(arclfile, '(30f16.8)')k12_shape(:, ikp), log(dos_l_mpi(ik12(1, ikp), ik12(2, ikp))) 
           if (mod(ikp, nky)==0) write(arclfile, *)' '
           write(arcrfile, '(30f16.8)')k12_shape(:, ikp), log(dos_r_mpi(ik12(1, ikp), ik12(2, ikp))) 
           if (mod(ikp, nky)==0) write(arcrfile, *)' '
           write(arcbulkfile, '(30f16.8)')k12_shape(:, ikp), log(abs(dos_bulk_mpi(ik12(1, ikp), ik12(2, ikp))))
           if (mod(ikp, nky)==0) write(arcbulkfile, *)' '
        enddo
        close(arclfile)
        close(arcrfile)
        close(arcbulkfile)
     endif
           
     if (cpuid.eq.0.and.SOC>0 .and. (SlabSpintexture_calc.or.SlabQPI_kplane_calc))then
        open(spindoslfile,file='spindos.dat_l')
        open(spindosrfile,file='spindos.dat_r')
        write(spindoslfile,'(a)')'# The coordinates of sx,sy,sz are redefined according to the SURFACE card'
        write(spindoslfile,'(a)')"# x axis is parallel to R1'"
        write(spindoslfile,'(a)')"# z axis is parallel to R1'xR2'"
        write(spindoslfile,'(a)')"# y axis is parallel to z x x"
        write(spindosrfile,'(a)')'# The coordinates of sx,sy,sz are redefined according to the SURFACE card'
        write(spindosrfile,'(a)')"# x axis is parallel to R1'"
        write(spindosrfile,'(a)')"# z axis is parallel to R1'xR2'"
        write(spindosrfile,'(a)')"# y axis is parallel to z x x"
        write(spindoslfile,'(30a16)')'#kx', 'ky',  &
                            'sx', 'sy', 'sz'
        write(spindosrfile,'(30a16)')'#kx', 'ky',  &
                            'sx', 'sy', 'sz'
        do ikp=1, nkx*nky
           if (dos_l_only(ik12(1, ikp), ik12(2, ikp))>eps6)then
              s0(1)= (sx_l_mpi(ik12(1, ikp), ik12(2, ikp)))/dos_l_mpi(ik12(1, ikp), ik12(2, ikp))
              s0(2)= (sy_l_mpi(ik12(1, ikp), ik12(2, ikp)))/dos_l_mpi(ik12(1, ikp), ik12(2, ikp))
              s0(3)= (sz_l_mpi(ik12(1, ikp), ik12(2, ikp)))/dos_l_mpi(ik12(1, ikp), ik12(2, ikp))
              call rotate(s0, s1)
              write(spindoslfile, '(30f16.8)')k12_shape(:, ikp), s1
           endif
           if (dos_r_only(ik12(1, ikp), ik12(2, ikp))>eps6)then
              s0(1)= (sx_r_mpi(ik12(1, ikp), ik12(2, ikp)))/dos_r_mpi(ik12(1, ikp), ik12(2, ikp))
              s0(2)= (sy_r_mpi(ik12(1, ikp), ik12(2, ikp)))/dos_r_mpi(ik12(1, ikp), ik12(2, ikp))
              s0(3)= (sz_r_mpi(ik12(1, ikp), ik12(2, ikp)))/dos_r_mpi(ik12(1, ikp), ik12(2, ikp))
              call rotate(s0, s1)
              write(spindosrfile, '(30f16.8)')k12_shape(:, ikp), s1
           endif
        enddo
        close(spindoslfile)
        close(spindosrfile)
     endif ! SlabSpintexture_calc and SlabQPI_kplane_calc=T

     if (cpuid.eq.0)then
        write(stdout,*)'Ndim: ',ndim
        write(stdout,*)'Nk1,Nk2,eta: ',nkx, nky, eta
        write(stdout,*)'Calculated surface density of state successfully'
     endif

     IF (SlabQPI_kplane_calc) then

     !> calculate QPI (jdos)
     do iq= 1+ cpuid, nkx*nky, num_cpu
        iq1= ik12(1, iq)- Nk1_half
        iq2= ik12(2, iq)- Nk2_half
        if (cpuid==0.and. mod(iq/num_cpu, 100)==0) &
           write(stdout, *) 'JDOS, iq ', iq, 'Nq',nkx*nky, 'time left', &
           (nkx*nky-iq)*(time_end- time_start)/num_cpu, ' s'
        call now(time_start)
        imin1= max(-Nk1_half-iq1, -Nk1_half)+ Nk1_half+ 1
        imax1= min(Nk1_half-iq1, Nk1_half)+ Nk1_half+ 1
        imin2= max(-Nk2_half-iq2, -Nk2_half)+ Nk2_half+ 1
        imax2= min(Nk2_half-iq2, Nk2_half)+ Nk2_half+ 1
        do ik2= imin2, imax2
           do ik1= imin1, imax1
              jdos_l(iq)= jdos_l(iq)+ dos_l_only(ik1, ik2)* dos_l_only(ik1+iq1, ik2+iq2)
              jdos_r(iq)= jdos_r(iq)+ dos_r_only(ik1, ik2)* dos_r_only(ik1+iq1, ik2+iq2)
           enddo !ik1
        enddo !ik2

        !> Only works if spin orbital coupling is considered.
        if (SOC>0) then
           do ik2= imin2, imax2
              do ik1= imin1, imax1
                 jsdos_l(iq)= jsdos_l(iq)+ dos_l_only(ik1, ik2)* dos_l_only(ik1+iq1, ik2+iq2) &
                                         + sx_l_mpi(ik1, ik2)* sx_l_mpi(ik1+iq1, ik2+iq2) &
                                         + sy_l_mpi(ik1, ik2)* sy_l_mpi(ik1+iq1, ik2+iq2) &
                                         + sz_l_mpi(ik1, ik2)* sz_l_mpi(ik1+iq1, ik2+iq2)
                 jsdos_r(iq)= jsdos_r(iq)+ dos_r_only(ik1, ik2)* dos_r_only(ik1+iq1, ik2+iq2) &
                                         + sx_r_mpi(ik1, ik2)* sx_r_mpi(ik1+iq1, ik2+iq2) &
                                         + sy_r_mpi(ik1, ik2)* sy_r_mpi(ik1+iq1, ik2+iq2) &
                                         + sz_r_mpi(ik1, ik2)* sz_r_mpi(ik1+iq1, ik2+iq2)
              enddo !ik1
           enddo !ik2
        endif
        call now(time_end)
     enddo !iq

     jdos_l_mpi=1d-12
     jdos_r_mpi=1d-12
     jsdos_l_mpi=1d-12
     jsdos_r_mpi=1d-12

#if defined (MPI)
     call mpi_reduce(jdos_l, jdos_l_mpi, size(jdos_l),mpi_double_precision,&
                     mpi_sum, 0, mpi_comm_world, ierr)
     call mpi_reduce(jdos_r, jdos_r_mpi, size(jdos_r),mpi_double_precision,&
                     mpi_sum, 0, mpi_comm_world, ierr)
     if (SOC>0) then
        call mpi_reduce(jsdos_l, jsdos_l_mpi, size(jsdos_l),mpi_double_precision,&
                        mpi_sum, 0, mpi_comm_world, ierr)
        call mpi_reduce(jsdos_r, jsdos_r_mpi, size(jsdos_r),mpi_double_precision,&
                        mpi_sum, 0, mpi_comm_world, ierr)
     endif
#else
     jdos_l_mpi= jdos_l
     jdos_r_mpi= jdos_r
     if (SOC>0) then
        jsdos_l_mpi= jsdos_l
        jsdos_r_mpi= jsdos_r
     endif
#endif
     jdos_l_mpi= jdos_l_mpi/nkx/nky
     jdos_r_mpi= jdos_r_mpi/nkx/nky

     ENDIF ! SlabQPI_kplane_calc= T

     if (cpuid.eq.0.and.SlabQPI_kplane_calc)then
        write(stdout,*)'The calculation of joint density of state was done.'
        write(stdout,*)'Now it is ready to write out.'
        open (unit=arcljfile, file='arc.jdat_l')
        write(arcljfile,'(a)')'# Surface joint density of states of bottom surface'
        write(arcljfile,'(a)')'# The coordinate of k is redefined according to the SURFACE card'
        write(arcljfile,'(a)')"# x axis is parallel to R1'"
        write(arcljfile,'(a)')"# z axis is parallel to R1'xR2'"
        write(arcljfile,'(a)')"# y axis is parallel to z x x"
        write(arcljfile,'(30a16)')'#kx', 'ky', 'log(dos)'
        do ikp=1, nkx*nky
           write(arcljfile, '(30f16.8)')k12_shape(:, ikp),  log(abs(jdos_l_mpi(ikp)))
           if (mod(ikp, nky)==0) write(arcljfile, *)' '
        enddo
        close(arcljfile)

        open (unit=arcrjfile, file='arc.jdat_r')
        write(arcrjfile,'(a)')'# Surface joint density of states of top surface'
        write(arcrjfile,'(a)')'# The coordinate of k is redefined according to the SURFACE card'
        write(arcrjfile,'(a)')"# x axis is parallel to R1'"
        write(arcrjfile,'(a)')"# z axis is parallel to R1'xR2'"
        write(arcrjfile,'(a)')"# y axis is parallel to z x x"
        write(arcrjfile,'(30a16)')'#kx', 'ky', 'log(dos)'
        do ikp=1, nkx*nky
           write(arcrjfile, '(30f16.8)')k12_shape(:, ikp),  log(abs(jdos_r_mpi(ikp)))
           if (mod(ikp, nky)==0) write(arcrjfile, *)' '
        enddo
        close(arcrjfile)
     endif ! cpuid==0

     if (cpuid.eq.0.and.SOC>0.and.SlabQPI_kplane_calc)then
        open (unit=arcljsfile, file='arc.jsdat_l')
        write(arcljsfile,'(a)')'# Surface joint density of states(in consideration of spin contribution) of bottom surface'
        write(arcljsfile,'(a)')'# The coordinate of k is redefined according to the SURFACE card'
        write(arcljsfile,'(a)')"# x axis is parallel to R1'"
        write(arcljsfile,'(a)')"# z axis is parallel to R1'xR2'"
        write(arcljsfile,'(a)')"# y axis is parallel to z x x"
        write(arcljsfile,'(30a16)')'#kx', 'ky', 'log(dos)'
        do ikp=1, nkx*nky
           write(arcljsfile, '(30f16.8)')k12_shape(:, ikp), log(abs(jsdos_l_mpi(ikp)))
           if (mod(ikp, nky)==0) write(arcljsfile, *)' '
        enddo
        close(arcljsfile)
   
        open (unit=arcrjsfile, file='arc.jsdat_r')
        write(arcrjsfile,'(a)')'# Surface joint density of states(in consideration of spin contribution) of top surface'
        write(arcrjsfile,'(a)')'# The coordinate of k is redefined according to the SURFACE card'
        write(arcrjsfile,'(a)')"# x axis is parallel to R1'"
        write(arcrjsfile,'(a)')"# z axis is parallel to R1'xR2'"
        write(arcrjsfile,'(a)')"# y axis is parallel to z x x"
        write(arcrjsfile,'(30a16)')'#kx', 'ky', 'log(dos)'
        do ikp=1, nkx*nky
           write(arcrjsfile, '(30f16.8)')k12_shape(:, ikp), log(abs(jsdos_r_mpi(ikp)))
           if (mod(ikp, nky)==0) write(arcrjsfile, *)' '
        enddo
        close(arcrjsfile)
        write(stdout,*)'calculate joint density of state successfully'
     endif !> SOC>0

     outfileindex= outfileindex+ 1
     if (cpuid.eq.0.and.SOC>0.and.SlabQPI_kplane_calc)then
        open(unit=outfileindex, file='arc_l_jsdos.gnu')
        write(outfileindex, '(a)')"set encoding iso_8859_1"
        write(outfileindex, '(a)')'#set terminal  postscript enhanced color'
        write(outfileindex, '(a)')"#set output 'arc_l_jsdos.eps'"
        write(outfileindex, '(3a)')'#set terminal pngcairo truecolor enhanced', &
           '  font ",50" size 1920, 1680'
        write(outfileindex, '(3a)')'set terminal png truecolor enhanced', &
           ' font ",50" size 1920, 1680'
        write(outfileindex, '(a)')"set output 'arc_l_jsdos.png'"
        write(outfileindex,'(2a)') 'set palette defined (0  "white", ', &
           '6 "red", 20 "black" )'
        write(outfileindex, '(a)')'#set palette rgbformulae 33,13,10'
        write(outfileindex, '(a)')'unset ztics'
        write(outfileindex, '(a)')'unset key'
        write(outfileindex, '(a)')'set pm3d'
        write(outfileindex, '(a)')'set border lw 6'
        write(outfileindex, '(a)')'set size ratio -1'
        write(outfileindex, '(a)')'set view map'
        write(outfileindex, '(a)')'set xtics'
        write(outfileindex, '(a)')'set ytics'
        write(outfileindex, '(a)')'set xlabel "K_1 (1/{\305})"'
        write(outfileindex, '(a)')'set ylabel "K_2 (1/{\305})"'
        write(outfileindex, '(a)')'set ylabel offset 1, 0'
        write(outfileindex, '(a)')'set colorbox'
        write(outfileindex, '(a)')'unset cbtics'
        write(outfileindex, '(a, f8.5, a, f8.5, a)')'set xrange [', k1min_shape, ':', k1max_shape, ']'
        write(outfileindex, '(a, f8.5, a, f8.5, a)')'set yrange [', k2min_shape, ':', k2max_shape, ']'
        write(outfileindex, '(a)')'set pm3d interpolate 2,2'
        write(outfileindex, '(2a)')"splot 'arc.jsdat_l' u 1:2:(exp($3)) w pm3d"
        close(outfileindex)
     endif !> SOC>0
  
       !> write script for gnuplot
     outfileindex= outfileindex+ 1
     if (cpuid.eq.0.and.SOC>0.and.SlabQPI_kplane_calc)then
        open(unit=outfileindex, file='arc_r_jsdos.gnu')
        write(outfileindex, '(a)')"set encoding iso_8859_1"
        write(outfileindex, '(a)')'#set terminal  postscript enhanced color'
        write(outfileindex, '(a)')"#set output 'arc_r_jsdos.eps'"
        write(outfileindex, '(3a)')'#set terminal pngcairo truecolor enhanced', &
           '  font ",50" size 1920, 1680'
        write(outfileindex, '(3a)')'set terminal png truecolor enhanced', &
           ' font ",50" size 1920, 1680'
        write(outfileindex, '(a)')"set output 'arc_r_jsdos.png'"
        write(outfileindex,'(2a)') 'set palette defined (0  "white", ', &
           '6 "red", 20 "black" )'
        write(outfileindex, '(a)')'#set palette rgbformulae 33,13,10'
        write(outfileindex, '(a)')'unset ztics'
        write(outfileindex, '(a)')'unset key'
        write(outfileindex, '(a)')'set pm3d'
        write(outfileindex, '(a)')'set border lw 6'
        write(outfileindex, '(a)')'set size ratio -1'
        write(outfileindex, '(a)')'set view map'
        write(outfileindex, '(a)')'set xtics'
        write(outfileindex, '(a)')'set ytics'
        write(outfileindex, '(a)')'set xlabel "K_1 (1/{\305})"'
        write(outfileindex, '(a)')'set ylabel "K_2 (1/{\305})"'
        write(outfileindex, '(a)')'set ylabel offset 1, 0'
        write(outfileindex, '(a)')'set colorbox'
        write(outfileindex, '(a)')'unset cbtics'
        write(outfileindex, '(a, f8.5, a, f8.5, a)')'set xrange [', k1min_shape, ':', k1max_shape, ']'
        write(outfileindex, '(a, f8.5, a, f8.5, a)')'set yrange [', k2min_shape, ':', k2max_shape, ']'
        write(outfileindex, '(a)')'set pm3d interpolate 2,2'
        write(outfileindex, '(2a)')"splot 'arc.jsdat_r' u 1:2:(exp($3)) w pm3d"
        close(outfileindex)
     endif !> SOC>0



     !> write script for gnuplot
     outfileindex= outfileindex+ 1
     if (cpuid.eq.0.and.SOC>0.and.SlabQPI_kplane_calc)then
        open(unit=outfileindex, file='arc_l_jdos.gnu')
        write(outfileindex, '(a)')"set encoding iso_8859_1"
        write(outfileindex, '(a)')'#set terminal  postscript enhanced color'
        write(outfileindex, '(a)')"#set output 'arc_l_jdos.eps'"
        write(outfileindex, '(3a)')'#set terminal pngcairo truecolor enhanced', &
           '  font ",50" size 1920, 1680'
        write(outfileindex, '(3a)')'set terminal png truecolor enhanced', &
           ' font ",50" size 1920, 1680'
        write(outfileindex, '(a)')"set output 'arc_l_jdos.png'"
        write(outfileindex,'(2a)') 'set palette defined (0  "white", ', &
           '6 "red", 20 "black" )'
        write(outfileindex, '(a)')'#set palette rgbformulae 33,13,10'
        write(outfileindex, '(a)')'unset ztics'
        write(outfileindex, '(a)')'unset key'
        write(outfileindex, '(a)')'set pm3d'
        write(outfileindex, '(a)')'set border lw 6'
        write(outfileindex, '(a)')'set size ratio -1'
        write(outfileindex, '(a)')'set view map'
        write(outfileindex, '(a)')'set xtics'
        write(outfileindex, '(a)')'set ytics'
        write(outfileindex, '(a)')'set xlabel "K_1 (1/{\305})"'
        write(outfileindex, '(a)')'set ylabel "K_2 (1/{\305})"'
        write(outfileindex, '(a)')'set ylabel offset 1, 0'
        write(outfileindex, '(a)')'set colorbox'
        write(outfileindex, '(a)')'unset cbtics'
        write(outfileindex, '(a, f8.5, a, f8.5, a)')'set xrange [', k1min_shape, ':', k1max_shape, ']'
        write(outfileindex, '(a, f8.5, a, f8.5, a)')'set yrange [', k2min_shape, ':', k2max_shape, ']'
        write(outfileindex, '(a)')'set pm3d interpolate 2,2'
        write(outfileindex, '(2a)')"splot 'arc.jdat_l' u 1:2:(exp($3)) w pm3d"
        close(outfileindex)
     endif

     !> write script for gnuplot
     outfileindex= outfileindex+ 1
     if (cpuid.eq.0.and.SOC>0.and.SlabQPI_kplane_calc)then
        open(unit=outfileindex, file='arc_r_jdos.gnu')
        write(outfileindex, '(a)')"set encoding iso_8859_1"
        write(outfileindex, '(a)')'#set terminal  postscript enhanced color'
        write(outfileindex, '(a)')"#set output 'arc_r_jdos.eps'"
        write(outfileindex, '(3a)')'#set terminal pngcairo truecolor enhanced', &
           '  font ",50" size 1920, 1680'
        write(outfileindex, '(3a)')'set terminal png truecolor enhanced', &
           ' font ",50" size 1920, 1680'
        write(outfileindex, '(a)')"set output 'arc_r_jdos.png'"
        write(outfileindex,'(2a)') 'set palette defined (0  "white", ', &
           '6 "red", 20 "black" )'
        write(outfileindex, '(a)')'#set palette rgbformulae 33,13,10'
        write(outfileindex, '(a)')'unset ztics'
        write(outfileindex, '(a)')'unset key'
        write(outfileindex, '(a)')'set pm3d'
        write(outfileindex, '(a)')'set border lw 6'
        write(outfileindex, '(a)')'set size ratio -1'
        write(outfileindex, '(a)')'set view map'
        write(outfileindex, '(a)')'set xtics'
        write(outfileindex, '(a)')'set ytics'
        write(outfileindex, '(a)')'set xlabel "K_1 (1/{\305})"'
        write(outfileindex, '(a)')'set ylabel "K_2 (1/{\305})"'
        write(outfileindex, '(a)')'set ylabel offset 1, 0'
        write(outfileindex, '(a)')'set colorbox'
        write(outfileindex, '(a)')'unset cbtics'
        write(outfileindex, '(a, f8.5, a, f8.5, a)')'set xrange [', k1min_shape, ':', k1max_shape, ']'
        write(outfileindex, '(a, f8.5, a, f8.5, a)')'set yrange [', k2min_shape, ':', k2max_shape, ']'
        write(outfileindex, '(a)')'set pm3d interpolate 2,2'
        write(outfileindex, '(2a)')"splot 'arc.jdat_r' u 1:2:(exp($3)) w pm3d"
        close(outfileindex)
     endif

     !> write script for gnuplot
     outfileindex= outfileindex+ 1
     if (cpuid==0) then
        open(unit=outfileindex, file='arc_l.gnu')
        write(outfileindex, '(a)')"set encoding iso_8859_1"
        write(outfileindex, '(a)')'#set terminal  postscript enhanced color'
        write(outfileindex, '(a)')"#set output 'arc_l.eps'"
        write(outfileindex, '(3a)')'#set terminal pngcairo truecolor enhanced', &
           '  font ",50" size 1920, 1680'
        write(outfileindex, '(3a)')'set terminal png truecolor enhanced', &
           ' font ",50" size 1920, 1680'
        write(outfileindex, '(a)')"set output 'arc_l.png'"
        write(outfileindex,'(2a)') 'set palette defined ( -10 "#194eff", ', &
           '0 "white", 10 "red" )'
        write(outfileindex, '(a)')'#set palette rgbformulae 33,13,10'
        write(outfileindex, '(a)')'unset ztics'
        write(outfileindex, '(a)')'unset key'
        write(outfileindex, '(a)')'set pm3d'
        write(outfileindex, '(a)')'set border lw 6'
        write(outfileindex, '(a)')'set size ratio -1'
        write(outfileindex, '(a)')'set view map'
        write(outfileindex, '(a)')'set xtics'
        write(outfileindex, '(a)')'set ytics'
        write(outfileindex, '(a)')'set xlabel "K_1 (1/{\305})"'
        write(outfileindex, '(a)')'set ylabel "K_2 (1/{\305})"'
        write(outfileindex, '(a)')'set ylabel offset 1, 0'
        write(outfileindex, '(a)')'set colorbox'
        write(outfileindex, '(a, f8.5, a, f8.5, a)')'set xrange [', k1min_shape, ':', k1max_shape, ']'
        write(outfileindex, '(a, f8.5, a, f8.5, a)')'set yrange [', k2min_shape, ':', k2max_shape, ']'
        write(outfileindex, '(a)')'set pm3d interpolate 2,2'
        write(outfileindex, '(2a)')"splot 'arc.dat_l' u 1:2:3 w pm3d"

        close(outfileindex)
     endif

     !> write script for gnuplot
     outfileindex= outfileindex+ 1
     if (cpuid==0) then
        open(unit=outfileindex, file='arc_r.gnu')
        write(outfileindex, '(a)')"set encoding iso_8859_1"
        write(outfileindex, '(a)')'#set terminal  postscript enhanced color'
        write(outfileindex, '(a)')"#set output 'arc_r.eps'"
        write(outfileindex, '(3a)')'#set terminal pngcairo truecolor enhanced', &
           '  font ",50" size 1920, 1680'
        write(outfileindex, '(3a)')'set terminal png truecolor enhanced', &
           ' font ",50" size 1920, 1680'
        write(outfileindex, '(a)')"set output 'arc_r.png'"
        write(outfileindex,'(2a)') 'set palette defined ( -10 "#194eff", ', &
           '0 "white", 10 "red" )'
        write(outfileindex, '(a)')'#set palette rgbformulae 33,13,10'
        write(outfileindex, '(a)')'unset ztics'
        write(outfileindex, '(a)')'unset key'
        write(outfileindex, '(a)')'set pm3d'
        write(outfileindex, '(a)')'set border lw 6'
        write(outfileindex, '(a)')'set size ratio -1'
        write(outfileindex, '(a)')'set view map'
        write(outfileindex, '(a)')'set xtics'
        write(outfileindex, '(a)')'set ytics'
        write(outfileindex, '(a)')'set xlabel "K_1 (1/{\305})"'
        write(outfileindex, '(a)')'set ylabel "K_2 (1/{\305})"'
        write(outfileindex, '(a)')'set ylabel offset 1, 0'
        write(outfileindex, '(a)')'set colorbox'
        write(outfileindex, '(a, f8.5, a, f8.5, a)')'set xrange [', k1min_shape, ':', k1max_shape, ']'
        write(outfileindex, '(a, f8.5, a, f8.5, a)')'set yrange [', k2min_shape, ':', k2max_shape, ']'
        write(outfileindex, '(a)')'set pm3d interpolate 2,2'
        write(outfileindex, '(2a)')"splot 'arc.dat_r' u 1:2:3 w pm3d"

        close(outfileindex)
     endif

     !> write script for gnuplot for bulk green's function
     outfileindex= outfileindex+ 1
     if (cpuid==0) then
        open(unit=outfileindex, file='arc_bulk.gnu')
        write(outfileindex, '(a)')"set encoding iso_8859_1"
        write(outfileindex, '(a)')'#set terminal  postscript enhanced color'
        write(outfileindex, '(a)')"#set output 'arc_l.eps'"
        write(outfileindex, '(3a)')'#set terminal pngcairo truecolor enhanced', &
           '  font ",50" size 1920, 1680'
        write(outfileindex, '(3a)')'set terminal png truecolor enhanced', &
           ' font ",50" size 1920, 1680'
        write(outfileindex, '(a)')"set output 'arc_bulk.png'"
        write(outfileindex,'(2a)') 'set palette defined ( -10 "#194eff", ', &
           '0 "white", 10 "red" )'
        write(outfileindex, '(a)')'#set palette rgbformulae 33,13,10'
        write(outfileindex, '(a)')'unset ztics'
        write(outfileindex, '(a)')'unset key'
        write(outfileindex, '(a)')'set pm3d'
        write(outfileindex, '(a)')'set border lw 6'
        write(outfileindex, '(a)')'set size ratio -1'
        write(outfileindex, '(a)')'set view map'
        write(outfileindex, '(a)')'set xtics'
        write(outfileindex, '(a)')'set ytics'
        write(outfileindex, '(a)')'set xlabel "K_1 (1/{\305})"'
        write(outfileindex, '(a)')'set ylabel "K_2 (1/{\305})"'
        write(outfileindex, '(a)')'set ylabel offset 1, 0'
        write(outfileindex, '(a)')'set colorbox'
        write(outfileindex, '(a, f8.5, a, f8.5, a)')'set xrange [', k1min_shape, ':', k1max_shape, ']'
        write(outfileindex, '(a, f8.5, a, f8.5, a)')'set yrange [', k2min_shape, ':', k2max_shape, ']'
        write(outfileindex, '(a)')'set pm3d interpolate 2,2'
        write(outfileindex, '(2a)')"splot 'arc.dat_bulk' u 1:2:3 w pm3d"
        close(outfileindex)
     endif

     outfileindex= outfileindex+ 1
     if (cpuid.eq.0.and.SOC>0 .and. (SlabSpintexture_calc.or.SlabQPI_kplane_calc))then
     !> generate gnuplot scripts for plotting the spin texture
     if (cpuid.eq.0) then
        open(outfileindex,file='spintext_r.gnu')
        write(outfileindex, '(a)')"set encoding iso_8859_1"
        write(outfileindex, '(a)')'#set terminal pngcairo truecolor enhanced font ",80" size 3680, 3360'
        write(outfileindex, '(a)')'set terminal png truecolor enhanced font ",80" size 3680, 3360'
        write(outfileindex, '(a)')"set output 'spintext_r.png'"
        write(outfileindex, '(a)')'set palette defined ( -6 "white", 0 "white", 10 "black" )'
        write(outfileindex, '(a)')"set multiplot layout 1,1 "
  
        write(outfileindex, '(a)')"set origin 0.06, 0.0"
        write(outfileindex, '(a)')"set size 0.9, 0.9"
        write(outfileindex, '(a)')"set xlabel 'K_1'"
        write(outfileindex, '(a)')"set ylabel 'K_2'"
        write(outfileindex, '(a)')"unset key"
        write(outfileindex, '(a)')"set pm3d"
        write(outfileindex, '(a)')"set xtics nomirror scale 0.5"
        write(outfileindex, '(a)')"set ytics nomirror scale 0.5"
        write(outfileindex, '(a)')"set border lw 6"
        write(outfileindex, '(a)')"set size ratio -1"
        write(outfileindex, '(a)')"set view map"
        write(outfileindex, '(a)')"unset colorbox"
        write(outfileindex, '(a, f10.5, a, f10.5, a)')"set xrange [", minval(k12_shape(1, :)), ":", &
           maxval(k12_shape(1, :)), "]"
        write(outfileindex, '(a, f10.5, a, f10.5, a)')"set yrange [", minval(k12_shape(2, :)), ":", &
           maxval(k12_shape(2, :)), "]"
        write(outfileindex, '(a)')"set pm3d interpolate 2,2"
        write(outfileindex, '(a)')"set label 1 'Spin texture' at graph 0.25, 1.10 front"
        write(outfileindex, '(a)')"splot 'arc.dat_r' u 1:2:3 w pm3d, \"
        write(outfileindex, '(a)')"     'spindos.dat_r' u 1:2:(0):($3/5.00):($4/5.00):(0)  w vec  head lw 5 lc rgb 'orange' front"
  
        close(outfileindex)
     endif
   
     !> generate gnuplot scripts for plotting the spin texture
     outfileindex= outfileindex+ 1
     if (cpuid.eq.0) then
        open(outfileindex,file='spintext_l.gnu')
        write(outfileindex, '(a)')"set encoding iso_8859_1"
        write(outfileindex, '(a)')'#set terminal pngcairo truecolor enhanced font ",80" size 3680, 3360'
        write(outfileindex, '(a)')'set terminal png truecolor enhanced font ",80" size 3680, 3360'
        write(outfileindex, '(a)')"set output 'spintext_l.png'"
        write(outfileindex, '(a)')'set palette defined ( -6 "white", 0 "white", 10 "black" )'
        write(outfileindex, '(a)')"set multiplot layout 1,1 "
  
        write(outfileindex, '(a)')"set origin 0.06, 0.0"
        write(outfileindex, '(a)')"set size 0.9, 0.9"
        write(outfileindex, '(a)')"set xlabel 'K_1'"
        write(outfileindex, '(a)')"set ylabel 'K_2'"
        write(outfileindex, '(a)')"unset key"
        write(outfileindex, '(a)')"set pm3d"
        write(outfileindex, '(a)')"set xtics nomirror scale 0.5"
        write(outfileindex, '(a)')"set ytics nomirror scale 0.5"
        write(outfileindex, '(a)')"set border lw 6"
        write(outfileindex, '(a)')"set size ratio -1"
        write(outfileindex, '(a)')"set view map"
        write(outfileindex, '(a)')"unset colorbox"
        write(outfileindex, '(a, f10.5, a, f10.5, a)')"set xrange [", minval(k12_shape(1, :)), ":", &
           maxval(k12_shape(1, :)), "]"
        write(outfileindex, '(a, f10.5, a, f10.5, a)')"set yrange [", minval(k12_shape(2, :)), ":", &
           maxval(k12_shape(2, :)), "]"
        write(outfileindex, '(a)')"set pm3d interpolate 2,2"
        write(outfileindex, '(a)')"set label 1 'Spin texture' at graph 0.25, 1.10 front"
        write(outfileindex, '(a)')"splot 'arc.dat_l' u 1:2:3 w pm3d, \"
        write(outfileindex, '(a)')"     'spindos.dat_l' u 1:2:(0):($3/5.00):($4/5.00):(0)  w vec  head lw 5 lc rgb 'orange' front"
  
        close(outfileindex)
     endif
     endif !> SlabSpintexture_calc or SlabQPI_kplane_calc
   
     if (cpuid.eq.0)write(stdout,*)'calculate spintexture successfully' 
 


#if defined (MPI)
     call mpi_barrier(mpi_cmw, ierr)
#endif

     deallocate( ik12,  k12,  k12_shape)
     deallocate( dos_l,  dos_l_mpi)
     deallocate( dos_r,  dos_r_mpi)
     deallocate( dos_r_only,  dos_l_only)
     deallocate( dos_bulk,  dos_bulk_mpi)
     deallocate( jdos_l,  jdos_l_mpi)
     deallocate( jdos_r,  jdos_r_mpi)
     if (SOC>0 .and. (SlabSpintexture_calc.or.SlabQPI_kplane_calc)) then
        deallocate( jsdos_l,  jsdos_l_mpi)
        deallocate( jsdos_r,  jsdos_r_mpi)
        deallocate( sx_l,  sx_l_mpi)
        deallocate( sy_l,  sy_l_mpi)
        deallocate( sz_l,  sz_l_mpi)
        deallocate( sx_r,  sx_r_mpi)
        deallocate( sy_r,  sy_r_mpi)
        deallocate( sz_r,  sz_r_mpi)
        deallocate(spin_sigma_x, spin_sigma_y, spin_sigma_z)
     endif
     deallocate( GLL, GRR, GB)
     deallocate(ctemp, H00, H01, ones)

  return
  end subroutine SurfaceDOSkk

