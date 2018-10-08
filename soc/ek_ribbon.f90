! This subroutine is used to caculate energy dispersion for
! slab Bi2Se3
! Copyright (c) 2010 QuanSheng Wu. All rights reserved.

  subroutine ek_ribbon

     use wmpi
     use para
     implicit none

     integer :: mdim
     integer :: ndim1

     ! loop index
     integer :: i
     integer :: j
     integer :: m
     integer :: l


     ! wave vector
     real(Dp) :: k

     real(dp) :: emin, emax
     real(Dp) :: kmax=0.5d0

     integer :: lwork

     ! the indices of the smallest and largest
     ! eigenvalues to be returned
     integer :: il,iu

     !the lower and upper bounds of the interval
     !to be searched for eigenvalues
     real(Dp) :: vl,vu

     !The absolute error tolerance for the eigenvalues
     real(Dp) :: abstol

     ! parameters for zheev
     integer :: ierr


     ! energy dispersion
     real(Dp),allocatable :: ekribbon(:,:)
     real(Dp),allocatable :: ekribbon_mpi(:,:)
     real(Dp),allocatable ::  rwork(:)
     integer,allocatable ::  iwork(:)
     complex(Dp),allocatable :: work(:)

     !> color for plot, surface state weight
     real(dp), allocatable :: surf_weight(:, :)
     real(dp), allocatable :: surf_weight_mpi(:, :)

     ! eigenvalue
     real(Dp),allocatable :: eigenvalue(:)

     ! hamiltonian slab
     complex(Dp),allocatable ::z(:,:)
     complex(Dp),allocatable ::CHamk(:,:)


     Ndim1=Num_wann*nslab1*nslab2
     lwork=64*Ndim1

     allocate(ekribbon(Ndim1,Nk1))
     allocate(ekribbon_mpi(Ndim1,Nk1))
     allocate(z(Ndim1,Ndim1))
     allocate(CHamk(Ndim1,Ndim1))
     allocate(rwork(7*Ndim1))
     allocate(work(lwork))
     allocate(iwork(5*Ndim1))
     allocate(eigenvalue(Ndim1))

     allocate( surf_weight (Ndim1, Nk1))
     allocate( surf_weight_mpi (Ndim1, Nk1))

     ierr = 0

     ! sweep k
     ekribbon=0.0d0
     ekribbon_mpi=0.0d0
     surf_weight= 0d0
     surf_weight_mpi= 0d0

     kmax=0.5d0

     abstol=1e-10
     vl= omegamin
     vu= omegamax
     il= (NumOccupied-2)*Nslab1*Nslab2
     iu= (NumOccupied+2)*Nslab1*Nslab2
     mdim=iu-il+1

     if (cpuid==0) write(stdout, *)'number of bands calculating: ',mdim

     do i=1+cpuid, Nk1, num_cpu
        if (cpuid==0) write(stdout, *) "Ribbonek the i'th kpoint", i, Nk1
        k=kmax*real(i-1)/(Nk1-1)
        chamk=0.0d0
        call ham_ribbon(k,Chamk)
        eigenvalue=0.0d0

        ! diagonal Chamk
        call eigensystem_c('V', 'U', Ndim1, CHamk, eigenvalue)

        ! only eigenvalues are computed
        ! the eigenvalues with indices il through iu will be found
        !call zheevx('N','I','U',Ndim1,Chamk,Ndim1,vl,vu,il,iu,abstol,&
        !mdim,eigenvalue,z,Ndim1,work,lwork,rwork,iwork,ifail,info)

        ekribbon(:,i)=eigenvalue

        do j=1, Nslab1* Nslab2* Num_wann  !< bands
           do m=1, Nslab2
              do l=1, Num_wann  !< sum over orbitals
                 surf_weight(j, i)= surf_weight(j, i) &
                 + abs(CHamk((m-1)*Num_wann+ l , j))**2 & ! first slab
                 + abs(CHamk(Num_wann*Nslab1*Nslab2 -(m-1)*Num_wann- l+ 1, j))**2 !& ! last slab
              enddo ! m
           enddo ! l

           do m=1, Nslab1
              do l=1, Num_wann  !< sum over orbitals
                 surf_weight(j, i)= surf_weight(j, i) &
                 + abs(CHamk((m-1)*Num_wann*Nslab2+ l , j))**2 & ! first slab
                 + abs(CHamk(Num_wann*Nslab1*Nslab2 -(m-1)*Num_wann*Nslab2- l+ 1, j))**2 !& ! last slab
              enddo ! m
           enddo ! l

           surf_weight(j, i)= sqrt(surf_weight(j, i))
        enddo ! j

        if (cpuid.eq.0) write(stdout,'(a2,i4,f12.5,f10.2,a2)')'k',i,ekribbon(1,i)
     enddo
#if defined (MPI)
     call mpi_allreduce(ekribbon,ekribbon_mpi,size(ekribbon),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
     call mpi_allreduce(surf_weight, surf_weight_mpi,size(surf_weight),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
     ekribbon_mpi= ekribbon
     surf_weight_mpi= surf_weight
#endif
     surf_weight= surf_weight_mpi/ maxval(surf_weight_mpi)

     if (cpuid.eq.0) then
        open(unit=100, file='ribbonek.dat',status='unknown')
        do j=1, Ndim1
           do i=1,Nk1
              k=-kmax*real(Nk1-i)/(Nk1-1)
              write(100,'(2f15.7, i8)')k,ekribbon_mpi(j,Nk1-i+1), &
                 int(255-surf_weight(j, Nk1-i+1)*255d0)
           enddo
           do i=1,Nk1
              k=kmax*real(i-1)/(Nk1-1)
              write(100,'(2f15.7, i8)')k,ekribbon_mpi(j,i), &
                 int(255-surf_weight(j, i)*255d0)
           enddo
           write(100, *)' '
        enddo

        close(100)
        write(stdout,*) 'calculate energy band  done'
     endif

     emin= minval(ekribbon_mpi)-0.5d0
     emax= maxval(ekribbon_mpi)+0.5d0
     !> write script for gnuplot
     if (cpuid==0) then
        open(unit=118, file='ribbonek.gnu')
        write(118, '(a)')'#set terminal  postscript enhanced color'
        write(118, '(a)')"#set output 'ribbonek.eps'"
        write(118, '(3a)')'#set terminal  pngcairo truecolor enhanced', &
           '# font ",36" size 1920, 1680'
        write(118, '(3a)')'set terminal  png truecolor enhanced', &
           '  font ",36" size 1920, 1680'
        write(118, '(a)')"set output 'ribbonek.png'"
        write(118,'(2a)') 'set palette defined ( 0  "green", ', &
           '5 "yellow", 10 "red" )'
        write(118, '(a)')'set style data linespoints'
        write(118, '(a)')'unset ztics'
        write(118, '(a)')'unset key'
        write(118, '(a)')'set pointsize 0.8'
        write(118, '(a)')'set border lw 3 '
        write(118, '(a)')'set view 0,0'
        write(118, '(a)')'#set xtics offset 0, -1'
        write(118, '(a)')'set ylabel offset -1, 0 '
        write(118, '(a)')'set ylabel "Energy (eV)"'
        write(118, '(a)')'set xrange [-0.5:0.5]'
        write(118, '(a, f10.5, a, f10.5, a)')'set yrange [', emin, ':', emax, ']'
        write(118, '(a)')'rgb(r,g,b) = int(r)*65536 + int(g)*256 + int(b)'
        write(118, '(2a)')"plot 'ribbonek.dat' u 1:2:(rgb(255,$3, 3)) ",  &
            "w lp lw 2 pt 7  ps 1 lc rgb variable"

        close(118)
     endif

     return
  end subroutine ek_ribbon


  SUBROUTINE wirestat_jdos

! QPI calculates for 1D system using                                   !
! green's function method  ---  J.Phys.F.Met.Phys.15(1985)851-858      !
!
! History:
!
! Original version by Quan Sheng Wu on 4/20/2010                                !
!
!            mpi version      4/21/2010
!
!            Quansheng Wu on Jan 30 2015 at ETH Zurich
!
! Modified for 1D QPI by Jianzhou Zhao on 10/08/2018 @ ETH Zurich
!
! Copyright (c) 2010 QuanSheng Wu. All rights reserved.

    USE wmpi
    USE para, ONLY: omeganum, omegamin, omegamax, ndim, knv2, k2_path, outfileindex, &
                    BottomOrbitals, TopOrbitals, NBottomOrbitals, NtopOrbitals, stdout, &
                    k2len, Num_wann, eps9, zi, Np
    IMPLICIT NONE

    ! MPI error code
    INTEGER  :: ierr

    ! file id
    INTEGER  :: jdoslfile, jdosrfile, dosbulkfile
    INTEGER  :: doslfile, dosrfile, spindoslfile, spindosrfile

    ! general loop index
    INTEGER  :: i, j, io, iq, iq1, ik1, ikp
    INTEGER  :: Nk_half, imin1, imax1, nw_half
    REAL(DP) :: ktmp(2), eta
    ! string for integer
    CHARACTER(LEN=140) :: ichar, jchar, kchar, fmt

    ! start & stop time
    REAL(DP) :: time_start, time_end, time_togo

    ! Omega
    REAL(DP),    ALLOCATABLE  :: omega(:)
    ! DOS for surface and bulk
    REAL(DP),    ALLOCATABLE  :: dos_l(:,:), dos_r(:,:), dos_bulk(:,:)
    ! DOS for surface ONLY
    REAL(DP),    ALLOCATABLE  :: dos_l_only(:,:), dos_r_only(:,:)
    ! JDOS for both side
    REAL(DP),    ALLOCATABLE  :: jdos_l(:,:), jdos_r(:,:), jdos_l_only(:,:), jdos_r_only(:,:)
    ! Array for MPI
    REAL(DP),    ALLOCATABLE  :: dos_l_mpi(:,:), dos_r_mpi(:,:), dos_bulk_mpi(:,:)
    REAL(DP),    ALLOCATABLE  :: jdos_l_mpi(:,:), jdos_r_mpi(:,:), jdos_l_only_mpi(:,:), jdos_r_only_mpi(:,:)
    ! Green's function
    COMPLEX(DP), ALLOCATABLE  :: GLL(:,:), GRR(:,:), GB (:,:)
    COMPLEX(DP), ALLOCATABLE  :: H00(:,:), H01(:,:)
    ! Unit array
    COMPLEX(DP), ALLOCATABLE  :: ones(:,:)
    ! Spin resolved component
    REAL(DP),    ALLOCATABLE  :: sx_l(:, :), sy_l(:, :), sz_l(:, :)
    REAL(DP),    ALLOCATABLE  :: sx_r(:, :), sy_r(:, :), sz_r(:, :)
    REAL(DP),    ALLOCATABLE  :: sx_l_mpi(:, :), sy_l_mpi(:, :), sz_l_mpi(:, :)
    REAL(DP),    ALLOCATABLE  :: sx_r_mpi(:, :), sy_r_mpi(:, :), sz_r_mpi(:, :)
    COMPLEX(DP), ALLOCATABLE  :: sigma_x(:,:), sigma_y(:,:), sigma_z(:,:)
    COMPLEX(DP), ALLOCATABLE  :: ctemp(:,:)

    ALLOCATE( omega(omeganum) )
    ALLOCATE( dos_l(knv2, omeganum), dos_r(knv2, omeganum), dos_bulk(knv2, omeganum) )
    ALLOCATE( dos_l_only(knv2, omeganum), dos_r_only(knv2, omeganum) )
    ALLOCATE( dos_l_mpi(knv2, omeganum), dos_r_mpi(knv2, omeganum) )
    ALLOCATE( jdos_l(knv2, omeganum), jdos_r(knv2, omeganum) )
    ALLOCATE( jdos_l_only(knv2, omeganum), jdos_r_only(knv2, omeganum) )
    ALLOCATE( dos_bulk_mpi(knv2, omeganum))
    ALLOCATE( jdos_l_mpi(knv2, omeganum), jdos_r_mpi(knv2, omeganum) )
    ALLOCATE( jdos_l_only_mpi(knv2, omeganum), jdos_r_only_mpi(knv2, omeganum) )
    ALLOCATE( GLL(Ndim, Ndim),GRR(Ndim, Ndim),GB (Ndim, Ndim) )
    ALLOCATE( H00(Ndim, Ndim),H01(Ndim, Ndim) )
    ALLOCATE( ones(Ndim, Ndim) )
    ALLOCATE( sx_l(knv2, omeganum), sy_l(knv2, omeganum), sz_l(knv2, omeganum))
    ALLOCATE( sx_l_mpi(knv2, omeganum), sy_l_mpi(knv2, omeganum), sz_l_mpi(knv2, omeganum))
    ALLOCATE( sx_r(knv2, omeganum), sy_r(knv2, omeganum), sz_r(knv2, omeganum))
    ALLOCATE( sx_r_mpi(knv2, omeganum), sy_r_mpi(knv2, omeganum), sz_r_mpi(knv2, omeganum))
    ALLOCATE( sigma_x(ndim,ndim), sigma_y(ndim,ndim), sigma_z(ndim,ndim), ctemp(ndim,ndim))

    omega        = 0d0;      ones         = 0d0
    dos_l        = 0d0;      dos_r        = 0d0
    dos_l_only   = 0d0;      dos_r_only   = 0d0
    dos_l_mpi    = 0d0;      dos_r_mpi    = 0d0
    dos_bulk     = 0d0;      dos_bulk_mpi = 0d0
    GLL          = 0d0;      GRR          = 0d0;      GB           = 0d0
    H00          = 0d0;      H01          = 0d0
    sigma_x      = 0d0;      sigma_y      = 0d0;      sigma_z      = 0d0
    sx_l         = 0d0;      sy_l         = 0d0;      sz_l         = 0d0
    sx_r         = 0d0;      sy_r         = 0d0;      sz_r         = 0d0
    sx_l_mpi     = 0d0;      sy_l_mpi     = 0d0;      sz_l_mpi     = 0d0
    sx_r_mpi     = 0d0;      sy_r_mpi     = 0d0;      sz_r_mpi     = 0d0

    ! Broaden coeffient
    eta=(omegamax- omegamin)/DBLE(omeganum)
    ! omega list
    DO i = 1, omeganum
        omega(i) = omegamin+(i-1)*eta
    ENDDO
    eta = eta * 1.5d0

    DO i=1,Ndim
        ones(i,i) = 1.0d0
    ENDDO

    !> spin operator matrix
    !> Note, the basis here should be |↑↑↓↓>
    nw_half = Num_wann/2
    DO i=1, Np
        DO j=1, nw_half
            sigma_x( Num_wann*(i-1)+j        , Num_wann*(i-1)+j+nw_half ) =  1.0d0
            sigma_x( Num_wann*(i-1)+j+nw_half, Num_wann*(i-1)+j         ) =  1.0d0
            sigma_y( Num_wann*(i-1)+j        , Num_wann*(i-1)+j+nw_half ) = -zi
            sigma_y( Num_wann*(i-1)+j+nw_half, Num_wann*(i-1)+j         ) =  zi
            sigma_z( Num_wann*(i-1)+j        , Num_wann*(i-1)+j         ) =  1.0d0
            sigma_z( Num_wann*(i-1)+j+nw_half, Num_wann*(i-1)+j+nw_half ) = -1.0d0
        ENDDO
    ENDDO

    time_start = 0d0
    time_end   = 0d0
    DO ikp = 1+cpuid, knv2, num_cpu

        CALL now(time_start)

        ktmp = k2_path(ikp,:)

        CALL ham_qlayer2qlayer(ktmp,h00,h01)

        DO j = 1, omeganum
            !> calculate surface green function
            ! there are two method to calculate surface green's function
            ! the method in 1985 is better, you can find the ref in the
            ! subroutine
            CALL surfgreen_1985(omega(j),GLL,GRR,GB,H00,H01,ones)
            ! call surfgreen_1984(w,GLL,GRR,H00,H01,ones)
            ! calculate spectral function
            DO i = 1, NtopOrbitals
                io = TopOrbitals(i)
                dos_l_mpi(ikp, j) = dos_l_mpi(ikp,j) - AIMAG(GLL(io,io))
            ENDDO ! i
            DO i = 1, NBottomOrbitals
                io = Ndim- Num_wann+ BottomOrbitals(i)
                dos_r_mpi(ikp, j) = dos_r_mpi(ikp,j) - AIMAG(GRR(io,io))
            ENDDO ! i
            DO i = 1, Ndim
                dos_bulk_mpi(ikp, j) = dos_bulk_mpi(ikp,j) - AIMAG(GB(i,i))
            ENDDO ! i

            ! Spin resolved sprectrum
            CALL mat_mul(ndim,gll,sigma_x,ctemp)
            DO i = 1, NtopOrbitals
                io = TopOrbitals(i)
                sx_l_mpi(ikp, j) = sx_l_mpi(ikp, j)- AIMAG(ctemp(io,io))
            ENDDO ! i
            CALL mat_mul(ndim,gll,sigma_y,ctemp)
            DO i = 1, NtopOrbitals
                io = TopOrbitals(i)
                sy_l_mpi(ikp, j) = sy_l_mpi(ikp, j)- AIMAG(ctemp(io,io))
            ENDDO !
            CALL mat_mul(ndim,gll,sigma_z,ctemp)
            DO i = 1, NtopOrbitals
                io = TopOrbitals(i)
                sz_l_mpi(ikp, j) = sz_l_mpi(ikp, j)- AIMAG(ctemp(io,io))
            ENDDO ! i
            CALL mat_mul(ndim,grr,sigma_x,ctemp)
            DO i = 1, NtopOrbitals
                io = TopOrbitals(i)
                sx_r_mpi(ikp, j) = sx_r_mpi(ikp, j)- AIMAG(ctemp(io,io))
            ENDDO ! i
            CALL mat_mul(ndim,grr,sigma_y,ctemp)
            DO i = 1, NtopOrbitals
                io = TopOrbitals(i)
                sy_r_mpi(ikp, j) = sy_r_mpi(ikp, j)- AIMAG(ctemp(io,io))
            ENDDO !
            CALL mat_mul(ndim,grr,sigma_z,ctemp)
            DO i = 1, NtopOrbitals
                io = TopOrbitals(i)
                sz_r_mpi(ikp, j) = sz_r_mpi(ikp, j)- AIMAG(ctemp(io,io))
            ENDDO ! i

        ENDDO ! j
        CALL now(time_end)

        time_togo = (knv2-ikp)*(time_end-time_start)/num_cpu
        IF (ikp == 1) THEN
            WRITE(ichar, *) knv2
            WRITE(jchar, '("I", I0)') LEN(TRIM(ADJUSTL(ichar)))
            WRITE(ichar, *) int(time_togo)
            WRITE(kchar, '("F", I0, ".2")') LEN(TRIM(ADJUSTL(ichar)))+4
            fmt = "(' SSs ik/nk : ',"//TRIM(jchar)//",'/',"//TRIM(jchar)//",' , &
                  & Estimated time remaining : ', "//TRIM(kchar)//", ' s')"
        ENDIF

        IF ( cpuid==0 .AND. MOD(ikp/num_cpu, 10)==0 .AND. ikp/=1 ) &
        WRITE( stdout, fmt ) ikp, knv2, time_togo

    ENDDO ! ikp

!> we do have to do allreduce operation
#ifdef MPI
    CALL mpi_allreduce(dos_l_mpi, dos_l, SIZE(dos_l), mpi_double_precision,&
    mpi_sum, mpi_comm_world, ierr)
    CALL mpi_allreduce(dos_r_mpi, dos_r, SIZE(dos_r), mpi_double_precision,&
    mpi_sum, mpi_comm_world, ierr)
    CALL mpi_allreduce(dos_bulk_mpi, dos_bulk, SIZE(dos_bulk), mpi_double_precision,&
    mpi_sum, mpi_comm_world, ierr)
    CALL mpi_allreduce(sx_l_mpi, sx_l, SIZE(sx_l), mpi_double_precision,&
    mpi_sum, mpi_comm_world, ierr)
    CALL mpi_allreduce(sy_l_mpi, sy_l, SIZE(sy_l), mpi_double_precision,&
    mpi_sum, mpi_comm_world, ierr)
    CALL mpi_allreduce(sz_l_mpi, sz_l, SIZE(sz_l), mpi_double_precision,&
    mpi_sum, mpi_comm_world, ierr)
    CALL mpi_allreduce(sx_r_mpi, sx_r, SIZE(sx_r), mpi_double_precision,&
    mpi_sum, mpi_comm_world, ierr)
    CALL mpi_allreduce(sy_r_mpi, sy_r, SIZE(sy_r), mpi_double_precision,&
    mpi_sum, mpi_comm_world, ierr)
    CALL mpi_allreduce(sz_r_mpi, sz_r, SIZE(sz_r), mpi_double_precision,&
    mpi_sum, mpi_comm_world, ierr)
#else
    dos_l    = dos_l_mpi;     dos_r    = dos_r_mpi;     dos_bulk = dos_bulk_mpi
    sx_l     = sx_l_mpi;      sy_l     = sy_l_mpi;      sz_l     = sz_l_mpi
    sx_r     = sx_r_mpi;      sy_r     = sy_r_mpi;      sz_r     = sz_r_mpi
#endif

    DEALLOCATE( dos_l_mpi, dos_r_mpi, dos_bulk_mpi )
    DEALLOCATE( sx_l_mpi, sy_l_mpi, sz_l_mpi )
    DEALLOCATE( sx_r_mpi, sy_r_mpi, sz_r_mpi )
    DEALLOCATE( sigma_x, sigma_y, sigma_z, ctemp )

    DO ikp=1, knv2
        DO j=1, omeganum
            dos_l_only(ikp, j) = dos_l(ikp, j)- dos_bulk(ikp, j)
            IF (dos_l_only(ikp, j)<0) dos_l_only(ikp, j) = eps9
            dos_r_only(ikp, j) = dos_r(ikp, j)- dos_bulk(ikp, j)
            IF (dos_r_only(ikp, j)<0) dos_r_only(ikp, j) = eps9
        ENDDO
    ENDDO

    outfileindex = outfileindex+ 1
    doslfile     = outfileindex
    outfileindex = outfileindex+ 1
    dosrfile     = outfileindex
    outfileindex = outfileindex+ 1
    dosbulkfile  = outfileindex
    outfileindex = outfileindex+ 1
    spindoslfile = outfileindex
    outfileindex = outfileindex+ 1
    spindosrfile = outfileindex

    ! Write surface state to files
    IF (cpuid .eq. 0) THEN
        OPEN(unit=doslfile    , file='dos.dat_l')
        OPEN(unit=dosrfile    , file='dos.dat_r')
        OPEN(unit=dosbulkfile , file='dos.dat_bulk')
        OPEN(unit=spindoslfile, file='spindos.dat_l')
        OPEN(unit=spindosrfile, file='spindos.dat_r')
        DO ikp = 1, knv2
            DO j = 1, omeganum
                WRITE(doslfile,    2002) k2len(ikp), omega(j), dos_l(ikp, j), dos_l_only(ikp, j)
                WRITE(dosrfile,    2002) k2len(ikp), omega(j), dos_r(ikp, j), dos_r_only(ikp, j)
                WRITE(dosbulkfile, 2003) k2len(ikp), omega(j), dos_bulk(ikp, j)
                WRITE(spindoslfile,2001) k2len(ikp), omega(j), sx_l(ikp, j), sy_l(ikp, j), sz_l(ikp,j)
                WRITE(spindosrfile,2001) k2len(ikp), omega(j), sx_r(ikp, j), sy_r(ikp, j), sz_r(ikp,j)
            ENDDO
            WRITE(doslfile    , *)
            WRITE(dosrfile    , *)
            WRITE(dosbulkfile , *)
            WRITE(spindoslfile, *)
            WRITE(spindosrfile, *)
        ENDDO
        CLOSE(doslfile)
        CLOSE(dosrfile)
        CLOSE(dosbulkfile)
    ENDIF

    ! Begin to calculate JDOS from DOS before
    Nk_half         = (knv2-1)/2
    jdos_l_mpi      = 0.0d0
    jdos_r_mpi      = 0.0d0
    jdos_l_only_mpi = 0.0d0
    jdos_r_only_mpi = 0.0d0
    DO j = 1+cpuid, omeganum,num_cpu
        DO iq= 1, knv2
            iq1 = iq - Nk_half
            imin1 = MAX(-Nk_half-iq1, -Nk_half)+ Nk_half+ 1
            imax1 = MIN( Nk_half-iq1,  Nk_half)+ Nk_half+ 1
            DO ik1 = imin1, imax1
                jdos_l_mpi(iq,j) = jdos_l_mpi(iq,j)+ dos_l(ik1, j)* dos_l(ik1+iq1, j)
                jdos_r_mpi(iq,j) = jdos_r_mpi(iq,j)+ dos_r(ik1, j)* dos_r(ik1+iq1, j)
                jdos_l_only_mpi(iq,j) = jdos_l_only_mpi(iq,j)+ dos_l_only(ik1, j)* dos_l_only(ik1+iq1, j)
                jdos_r_only_mpi(iq,j) = jdos_r_only_mpi(iq,j)+ dos_r_only(ik1, j)* dos_r_only(ik1+iq1, j)
            ENDDO !ik1
        ENDDO !iq
    ENDDO

#ifdef MPI
    CALL mpi_reduce(jdos_l_mpi, jdos_l, SIZE(jdos_l), mpi_double_precision,&
    mpi_sum, 0, mpi_comm_world, ierr)
    CALL mpi_reduce(jdos_r_mpi, jdos_r, SIZE(jdos_r), mpi_double_precision,&
    mpi_sum, 0, mpi_comm_world, ierr)
    CALL mpi_reduce(jdos_l_only_mpi, jdos_l_only, SIZE(jdos_l_only), mpi_double_precision,&
    mpi_sum, 0, mpi_comm_world, ierr)
    CALL mpi_reduce(jdos_r_only_mpi, jdos_r_only, SIZE(jdos_r_only), mpi_double_precision,&
    mpi_sum, 0, mpi_comm_world, ierr)
#else
    jdos_l      = jdos_l_mpi;         jdos_r      = jdos_r_mpi
    jdos_l_only = jdos_l_only_mpi;    jdos_r_only = jdos_r_only_mpi
#endif

    outfileindex = outfileindex+ 1
    jdoslfile    = outfileindex
    outfileindex = outfileindex+ 1
    jdosrfile    = outfileindex

    IF(cpuid.eq.0) THEN
        OPEN(unit=jdoslfile, file='dos.jdat_l')
        OPEN(unit=jdosrfile, file='dos.jdat_r')
        DO ikp = 1, knv2
            DO j = 1, omeganum
                WRITE(jdoslfile, 2002) k2len(ikp), omega(j), jdos_l(ikp, j), jdos_l_only(ikp, j)
                WRITE(jdosrfile, 2002) k2len(ikp), omega(j), jdos_r(ikp, j), jdos_r_only(ikp, j)
            ENDDO
            WRITE(jdoslfile, *)
            WRITE(jdosrfile, *)
        ENDDO
        CLOSE(jdoslfile)
        CLOSE(jdosrfile)
        WRITE(stdout,*)
        WRITE(stdout,*)'calculate joint density of state successfully'
    ENDIF

    DEALLOCATE( ones, omega )
    DEALLOCATE( dos_l, dos_r, dos_l_only, dos_r_only, dos_bulk )
    DEALLOCATE( GLL, GRR, GB )
    DEALLOCATE( H00, H01 )
    DEALLOCATE( jdos_l, jdos_r, jdos_l_only, jdos_r_only )
    DEALLOCATE( jdos_l_mpi, jdos_r_mpi, jdos_l_only_mpi, jdos_r_only_mpi )
    DEALLOCATE( sx_l, sy_l, sz_l )
    DEALLOCATE( sx_r, sy_r, sz_r )

2001 FORMAT(5(1X,F16.8))
2002 FORMAT(4(1X,F16.8))
2003 FORMAT(3(1X,F16.8))

END SUBROUTINE wirestat_jdos

