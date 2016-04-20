
! this subroutine is used to calculate spin texture of
! the surface state of Bi2Se3 system
! 
! constructed by QS.Wu on July/20/2010  10:26:13


  subroutine spintext
     
    use para
    use mpi

    implicit none

    integer :: i, j, ikp

	 integer ::  nkz, nkx
    integer :: ierr

    integer :: Nband

    real(Dp)   :: k(2)
   
    real(Dp), allocatable   :: kxz(:,:) 
    
	 real(dp) :: kxmin, kxmax, kzmin, kzmax, omega

    real(Dp), allocatable   :: sx(:)
    real(Dp), allocatable   :: sy(:)
    real(Dp), allocatable   :: sz(:)
    real(Dp), allocatable   :: sx_mpi(:)
    real(Dp), allocatable   :: sy_mpi(:)
    real(Dp), allocatable   :: sz_mpi(:)
    real(Dp), allocatable   :: dos(:)
    real(Dp), allocatable   :: dos_mpi(:)
! spin operator matrix sigma_x,sigma_y in sigma_z representation
    complex(Dp),allocatable :: sigma_x(:,:) 
    complex(Dp),allocatable :: sigma_y(:,:) 
    complex(Dp),allocatable :: sigma_z(:,:) 

! surface green function

    complex(Dp),allocatable :: surfgreen(:,:)
    complex(Dp),allocatable :: GRR(:,:)
    complex(Dp),allocatable :: H00(:,:)
    complex(Dp),allocatable :: H01(:,:)
    complex(Dp),allocatable :: ones(:,:)
    complex(Dp),allocatable :: ctemp(:,:)

    allocate(sigma_x(ndim,ndim))
    allocate(sigma_y(ndim,ndim))
    allocate(sigma_z(ndim,ndim))
    allocate(surfgreen(ndim,ndim))
    allocate(ctemp(ndim,ndim))

    sigma_x   =0.0d0
    sigma_y   =0.0d0
    sigma_z   =0.0d0
    surfgreen =0.0d0
    ctemp     =0.0d0
    
    nkx = Nk 
    nkz = Nk

    allocate( kxz(2, nkx*nkz))
    allocate( dos(nkx*nkz))
    allocate( dos_mpi(nkx*nkz))
    allocate( GRR(ndim,ndim))
	 allocate(sx(knv2), sy(knv2), sz(knv2))
	 allocate(sx_mpi(knv2), sy_mpi(knv2), sz_mpi(knv2))
    allocate(H00(Ndim, Ndim))
    allocate(H01(Ndim, Ndim))
    allocate(ones(Ndim, Ndim))
    kxz=0d0
    dos=0d0
	 sx=0d0
	 sy=0d0
	 sz=0d0
	 sx_mpi=0d0
	 sy_mpi=0d0
	 sz_mpi=0d0
    GRR= 0d0
    H00= 0d0
    H01= 0d0
    ones= 0d0

    do i=1,Ndim
       ones(i,i)=1.0d0
    enddo
    
    ikp=0
    do i= 1, nkx
       do j= 1, nkz
          ikp=ikp+1
          kxz(:, ikp)=K2D_start1+ (i-1)*(K2D_end1- K2D_start1)/dble(nkx-1) &
                     +K2D_start2+ (j-1)*(K2D_end2- K2D_start2)/dble(nkz-1)
       enddo
    enddo


    Nband= Num_wann/2

    !> spin operator matrix
    do i=1, Np
       do j=1, Nband
          sigma_x(Num_wann*(i-1)+j, Num_wann*(i-1)+Nband+j)=1.0d0
          sigma_x(Num_wann*(i-1)+j+Nband, Num_wann*(i-1)+j)=1.0d0
          sigma_y(Num_wann*(i-1)+j, Num_wann*(i-1)+Nband+j)=-zi
          sigma_y(Num_wann*(i-1)+j+Nband, Num_wann*(i-1)+j)=zi
          sigma_z(Num_wann*(i-1)+j, Num_wann*(i-1)+j)= 1d0
          sigma_z(Num_wann*(i-1)+j+Nband, Num_wann*(i-1)+j+Nband)=-1d0
       enddo
    enddo

    if (cpuid.eq.0)write(stdout,*)'sigma_x'
	 do i=1, ndim
       if (cpuid.eq.0)write(stdout,'(240f3.0)')real(sigma_x(i,:)) 
	 enddo

    if (cpuid.eq.0)write(stdout,*)'sigma_y'
	 do i=1, ndim
       if (cpuid.eq.0)write(stdout,'(240f3.0)')aimag(sigma_y(i,:)) 
	 enddo

    if (cpuid.eq.0)write(stdout,*)'sigma_z'
	 do i=1, ndim
       if (cpuid.eq.0)write(stdout,'(240f3.0)')real(sigma_z(i,:)) 
	 enddo

    sx=0.0d0
    sy=0.0d0
    sz=0.0d0
 
    omega = E_arc
    eta= eta_arc

    do i=1+cpuid,knv2,num_cpu 
       if (cpuid==0) print *,i, knv2

       k=kxz(:,i) 

       !> get the hopping matrix between two principle layers
       call ham_qlayer2qlayer(k,H00,H01) 
      
       ! calculate surface green function
       call surfgreen_1985(omega,surfgreen,GRR,H00,H01,ones)

       !ctemp=matmul(surfgreen,sigma_x)       
       call mat_mul(ndim,surfgreen,sigma_x,ctemp)       
       do j=1,ndim 
          sx(i)=sx(i)-Aimag(ctemp(j,j))
       enddo

       !ctemp=matmul(surfgreen,sigma_y)       
       call mat_mul(ndim,surfgreen,sigma_y,ctemp)       
       do j=1,ndim 
          sy(i)=sy(i)-Aimag(ctemp(j,j))
       enddo

       !ctemp=matmul(surfgreen,sigma_z)       
       call mat_mul(ndim,surfgreen,sigma_z,ctemp)       
       do j=1,ndim 
          sz(i)=sz(i)-Aimag(ctemp(j,j))
       enddo

       do j=1,ndim
          dos(i)=dos(i)-Aimag(surfgreen(j,j))
       enddo


    enddo  ! i  kpoint
    
    call mpi_allreduce(sx,sx_mpi,size(sx),mpi_double_precision, &
                  mpi_sum , mpi_comm_world, ierr)
    call mpi_allreduce(sy,sy_mpi,size(sy),mpi_double_precision, &
                  mpi_sum , mpi_comm_world, ierr)
    call mpi_allreduce(sz,sz_mpi,size(sz),mpi_double_precision, &
                  mpi_sum , mpi_comm_world, ierr)
    call mpi_allreduce(dos,dos_mpi,size(dos),mpi_double_precision, &
                  mpi_sum , mpi_comm_world, ierr)

    sx = sx_mpi/pi
    sy = sy_mpi/pi
    sz = sz_mpi/pi
    dos= dos_mpi

    if (cpuid.eq.0)then
       open(100,file='spindos.dat')

       ikp= 0
       do i= 1, nkx
          do j= 1, nkz
             ikp=ikp+1
             k=kxz(:, ikp)
	   		 if (dos(ikp).gt.dos(ikp+1).and.dos(ikp).gt.dos(ikp-1).and. &
                real(log(dos(ikp)))>5.5) then
                write(100,'(6f16.8)')k,real(log(dos(ikp))), &
                   real(sx(ikp)),real(sy(ikp)),real(sz(ikp))                 
             else
                write(100,'(6f16.8)')k,real(log(dos(ikp))), &
                   real(sx(ikp)),real(sy(ikp)),real(sz(ikp))                 
	   	    endif
          enddo
          write(100, *)' '
       enddo

       close(100)
    endif

  
    if (cpuid.eq.0)write(stdout,*)'calculate spintexture successfully' 
 
    return
  end subroutine spintext
