! this subroutine is used to calculate spin texture of
! the surface state 
! 
! constructed by QS.Wu on July/20/2010  10:26:13
! constructed by QS.Wu on July/20/2010  10:26:13
!> modified by QS.Wu on June/4/2016 in Jouvence Montreal Canada 

  subroutine spintext
     
    use para
    use mpi

    implicit none

    integer :: i, j, ikp

    integer ::  nk2, nk1
    integer :: ierr

    integer :: Nband

    real(Dp)   :: k(2)
   
    real(Dp), allocatable   :: k12(:,:) 
    real(Dp), allocatable   :: k12_shape(:,:) 
    
	 real(dp) :: omega

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
    
    nk1 = Nk 
    nk2 = Nk

    allocate( k12(2, nk1*nk2))
    allocate( k12_shape(2, nk1*nk2))
    allocate( dos(nk1*nk2))
    allocate( dos_mpi(nk1*nk2))
    allocate( GRR(ndim,ndim))
    allocate(sx(nk1*nk2), sy(nk1*nk2), sz(nk1*nk2))
    allocate(sx_mpi(nk1*nk2), sy_mpi(nk1*nk2), sz_mpi(nk1*nk2))
    allocate(H00(Ndim, Ndim))
    allocate(H01(Ndim, Ndim))
    allocate(ones(Ndim, Ndim))
    k12=0d0
    k12_shape=0d0
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
    do i= 1, nk1
       do j= 1, nk2
          ikp=ikp+1
          k12(:, ikp)=K2D_start+ (i-1)*K2D_vec1/dble(nk1-1) &
                     + (j-1)*K2D_vec2/dble(nk2-1)
          k12_shape(:, ikp)= k12(1, ikp)* Ka2+ k12(2, ikp)* Kb2
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

    do ikp=1+cpuid,nk1*nk2,num_cpu 
       if (cpuid==0) write(stdout, *)'spintexture',ikp, nk1*nk2
       if (cpuid==0) write(*, *)'spintexture',ikp, nk1*nk2

       k=k12(:,ikp) 

       !> get the hopping matrix between two principle layers
       call ham_qlayer2qlayer(k,H00,H01) 
      
       ! calculate surface green function
       call surfgreen_1985(omega,surfgreen,GRR,H00,H01,ones)

       !> surface state spectrum for GLL
       do j=1,ndim
          dos(ikp)=dos(ikp)-Aimag(surfgreen(j,j))
       enddo


       !ctemp=matmul(surfgreen,sigma_x)       
       call mat_mul(ndim,surfgreen,sigma_x,ctemp)       
       do j=1,ndim 
          sx(ikp)=sx(ikp)-Aimag(ctemp(j,j))
       enddo

       !ctemp=matmul(surfgreen,sigma_y)       
       call mat_mul(ndim,surfgreen,sigma_y,ctemp)       
       do j=1,ndim 
          sy(ikp)=sy(ikp)-Aimag(ctemp(j,j))
       enddo

       !ctemp=matmul(surfgreen,sigma_z)       
       call mat_mul(ndim,surfgreen,sigma_z,ctemp)       
       do j=1,ndim 
          sz(ikp)=sz(ikp)-Aimag(ctemp(j,j))
       enddo
       sx(ikp)= sx(ikp)/dos(ikp)
       sy(ikp)= sy(ikp)/dos(ikp)
       sz(ikp)= sz(ikp)/dos(ikp)

    enddo  ! ikp  kpoint
    
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
       open(300,file='spindos.dat')
       open(301,file='spintext.dat')

       ikp= 0
       do i= 1, nk1
          do j= 1, nk2
             ikp=ikp+1
             k=k12_shape(:, ikp)
             write(300,'(60f16.8)')k, (log(dos(ikp))), &
                real(sx(ikp)),real(sy(ikp)),real(sz(ikp))                 
             if (dos(ikp).gt.dos(ikp+1).and.dos(ikp).gt.dos(ikp-1).and. &
                real(log(dos(ikp)))>0.5) then
                write(301,'(60f16.8)')k, (log(dos(ikp))), &
                   real(sx(ikp)),real(sy(ikp)),real(sz(ikp))                 
             endif
          enddo
          write(300, *)' '
       enddo

       close(300)
       close(301)
    endif


    !> generate gnuplot scripts for plotting the spin texture
    if (cpuid.eq.0) then
       open(302,file='spintext.gnu')
       write(302, '(a)')'set terminal pngcairo truecolor enhanced font ",80" size 3680, 3360'
       write(302, '(a)')"set output 'spintext.png'"
       write(302, '(a)')'set palette defined ( -6 "white", 0 "gray", 10 "blue" )'
       write(302, '(a)')"set multiplot layout 1,1 "

       write(302, '(a)')"set origin 0.06, 0.0"
       write(302, '(a)')"set size 0.9, 0.9"
       write(302, '(a)')"set xlabel 'K_1'"
       write(302, '(a)')"set ylabel 'K_2'"
       write(302, '(a)')"unset key"
       write(302, '(a)')"set pm3d"
       write(302, '(a)')"set xtics 0.1"
       write(302, '(a)')"set ytics 0.2"
       write(302, '(a)')"set border lw 6"
       write(302, '(a)')"set size ratio -1"
       write(302, '(a)')"set view map"
       write(302, '(a)')"unset colorbox"
       write(302, '(a, f10.5, a, f10.5, a)')"set xrange [", minval(k12_shape(1, :)), ":", &
          maxval(k12_shape(1, :)), "]"
       write(302, '(a, f10.5, a, f10.5, a)')"set yrange [", minval(k12_shape(2, :)), ":", &
          maxval(k12_shape(2, :)), "]"
       write(302, '(a)')"set pm3d interpolate 2,2"
       write(302, '(a)')"splot 'spindos.dat' u 1:2:3 w pm3d"
       write(302, '(a)')"unset xtics"
       write(302, '(a)')"unset ytics"
       write(302, '(a)')"unset xlabel"
       write(302, '(a)')"unset ylabel"
       write(302, '(a)')"set label 1 'Spin texture' at graph 0.25, 1.10 front"
       write(302, '(a)')"# the sencond plot"
       write(302, '(a)')"set origin 0.14, 0.115"
       write(302, '(a)')"set size 0.75, 0.70"
       write(302, '(a)')"plot 'spintext.dat' u 1:2:($4/5.00):($5/5.00)  w vec  head lw 4 lc rgb 'red' front"

       close(302)
    endif
  
    if (cpuid.eq.0)write(stdout,*)'calculate spintexture successfully' 
 
    return
  end subroutine spintext
