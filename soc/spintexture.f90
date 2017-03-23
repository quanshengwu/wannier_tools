! this subroutine is used to calculate spin texture of
! the surface state 
! 
! constructed by QS.Wu on July/20/2010  10:26:13
! constructed by QS.Wu on July/20/2010  10:26:13
!> modified by QS.Wu on June/4/2016 in Jouvence Montreal Canada 

  subroutine spintext
     
    use para
    use wmpi

    implicit none

    integer :: i, j, ikp

    integer :: ierr
    integer :: spindosfile, spintexturefile

    integer :: Nband

    real(Dp)   :: k(2)
   
    real(Dp), allocatable   :: k12(:,:) 
    real(Dp), allocatable   :: k12_shape(:,:) 
    
	 real(dp) :: omega
    real(dp) :: dos_l_max, dos_r_max

    real(Dp), allocatable   :: sx_l(:)
    real(Dp), allocatable   :: sy_l(:)
    real(Dp), allocatable   :: sz_l(:)
    real(Dp), allocatable   :: sx_l_mpi(:)
    real(Dp), allocatable   :: sy_l_mpi(:)
    real(Dp), allocatable   :: sz_l_mpi(:)
    real(Dp), allocatable   :: dos_l(:)
    real(Dp), allocatable   :: dos_bulk(:)
    real(Dp), allocatable   :: dos_l_mpi(:)
    real(dp), allocatable   :: dos_l_only(:)
    real(dp), allocatable   :: dos_r_only(:)
    real(Dp), allocatable   :: sx_r(:)
    real(Dp), allocatable   :: sy_r(:)
    real(Dp), allocatable   :: sz_r(:)
    real(Dp), allocatable   :: sx_r_mpi(:)
    real(Dp), allocatable   :: sy_r_mpi(:)
    real(Dp), allocatable   :: sz_r_mpi(:)
    real(Dp), allocatable   :: dos_r(:)
    real(Dp), allocatable   :: dos_r_mpi(:)

    ! spin operator matrix sigma_x,sigma_y in sigma_z representation
    complex(Dp),allocatable :: sigma_x(:,:) 
    complex(Dp),allocatable :: sigma_y(:,:) 
    complex(Dp),allocatable :: sigma_z(:,:) 

    ! surface green function
    complex(Dp),allocatable :: GLL(:,:)
    complex(Dp),allocatable :: GRR(:,:)
    complex(Dp),allocatable :: GB (:,:)
    complex(Dp),allocatable :: H00(:,:)
    complex(Dp),allocatable :: H01(:,:)
    complex(Dp),allocatable :: ones(:,:)
    complex(Dp),allocatable :: ctemp(:,:)

    allocate(sigma_x(ndim,ndim))
    allocate(sigma_y(ndim,ndim))
    allocate(sigma_z(ndim,ndim))
    allocate(GLL(ndim,ndim))
    allocate(ctemp(ndim,ndim))

    sigma_x   =0.0d0
    sigma_y   =0.0d0
    sigma_z   =0.0d0
    GLL =0.0d0
    ctemp     =0.0d0
    

    allocate( k12(2, nk1*nk2))
    allocate( k12_shape(2, nk1*nk2))
    allocate( dos_l(nk1*nk2))
    allocate( dos_l_mpi(nk1*nk2))
    allocate( dos_r(nk1*nk2))
    allocate( dos_bulk(nk1*nk2))
    allocate( dos_r_mpi(nk1*nk2))
    allocate( dos_l_only(nk1*nk2))
    allocate( dos_r_only(nk1*nk2))
    allocate( GRR(ndim,ndim))
    allocate( GB (ndim,ndim))
    allocate(sx_l(nk1*nk2), sy_l(nk1*nk2), sz_l(nk1*nk2))
    allocate(sx_l_mpi(nk1*nk2), sy_l_mpi(nk1*nk2), sz_l_mpi(nk1*nk2))
    allocate(sx_r(nk1*nk2), sy_r(nk1*nk2), sz_r(nk1*nk2))
    allocate(sx_r_mpi(nk1*nk2), sy_r_mpi(nk1*nk2), sz_r_mpi(nk1*nk2))
    allocate(H00(Ndim, Ndim))
    allocate(H01(Ndim, Ndim))
    allocate(ones(Ndim, Ndim))
    k12=0d0
    k12_shape=0d0
    dos_l=0d0
    dos_r=0d0
    dos_bulk=0d0
    sx_l=0d0
    sy_l=0d0
    sz_l=0d0
    sx_l_mpi=0d0
    sy_l_mpi=0d0
    sz_l_mpi=0d0
    sx_r=0d0
    sy_r=0d0
    sz_r=0d0
    sx_r_mpi=0d0
    sy_r_mpi=0d0
    sz_r_mpi=0d0
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
    !> this part is package dependent. 
    if (index( Package, 'VASP')/=0.or. index( Package, 'Wien2k')/=0 &
       .or. index( Package, 'Abinit')/=0.or. index( Package, 'openmx')/=0) then
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
    elseif (index( Package, 'QE')/=0.or.index( Package, 'quantumespresso')/=0 &
         .or.index( Package, 'quantum-espresso')/=0.or.index( Package, 'pwscf')/=0) then
       do i=1, Np
          do j=1, Nband
             sigma_x(Num_wann*(i-1)+(2*j-1), Num_wann*(i-1)+2*j)=1.0d0
             sigma_x(Num_wann*(i-1)+2*j, Num_wann*(i-1)+(2*j-1))=1.0d0
             sigma_y(Num_wann*(i-1)+(2*j-1), Num_wann*(i-1)+2*j)=-zi
             sigma_y(Num_wann*(i-1)+2*j, Num_wann*(i-1)+(2*j-1))=zi
             sigma_z(Num_wann*(i-1)+(2*j-1), Num_wann*(i-1)+(2*j-1))=1.0d0
             sigma_z(Num_wann*(i-1)+2*j, Num_wann*(i-1)+2*j)=-1.0d0
          enddo
       enddo
    else
       if (cpuid.eq.0) write(stdout, *)'Error: please report your software and wannier90.wout to me'
       if (cpuid.eq.0) write(stdout, *)'wuquansheng@gmail.com'
       stop 'Error: please report your software and wannier90.wout to wuquansheng@gmail.com'
    endif


   !if (cpuid.eq.0)write(stdout,*)'sigma_x'
    do i=1, ndim
      !if (cpuid.eq.0)write(stdout,'(240f3.0)')real(sigma_x(i,:)) 
    enddo

   !if (cpuid.eq.0)write(stdout,*)'sigma_y'
    do i=1, ndim
      !if (cpuid.eq.0)write(stdout,'(240f3.0)')aimag(sigma_y(i,:)) 
    enddo

   !if (cpuid.eq.0)write(stdout,*)'sigma_z'
    do i=1, ndim
      !if (cpuid.eq.0)write(stdout,'(240f3.0)')real(sigma_z(i,:)) 
    enddo

    omega = E_arc
    eta= eta_arc

    do ikp=1+cpuid,nk1*nk2,num_cpu 
       if (cpuid==0) write(stdout, *)'spintexture, ik, Nk1*Nk2 ',ikp, nk1*nk2

       k=k12(:,ikp) 

       !> get the hopping matrix between two principle layers
       call ham_qlayer2qlayer(k,H00,H01) 
      
       ! calculate surface green function
       call surfgreen_1985(omega,GLL,GRR,GB,H00,H01,ones)

       !> surface state spectrum for GLL
       do j=1,ndim
          dos_l(ikp)=dos_l(ikp)-Aimag(GLL(j,j))
          dos_r(ikp)=dos_r(ikp)-Aimag(GRR(j,j))
          dos_bulk(ikp)=dos_bulk(ikp)-Aimag(GB(j,j))
       enddo

       !ctemp=matmul(GLL,sigma_x)       
       call mat_mul(ndim,GLL,sigma_x,ctemp)       
       do j=1,ndim 
          sx_l(ikp)=sx_l(ikp)-Aimag(ctemp(j,j))
       enddo

       !ctemp=matmul(GLL,sigma_y)       
       call mat_mul(ndim,GLL,sigma_y,ctemp)       
       do j=1,ndim 
          sy_l(ikp)=sy_l(ikp)-Aimag(ctemp(j,j))
       enddo

       !ctemp=matmul(GLL,sigma_z)       
       call mat_mul(ndim,GLL,sigma_z,ctemp)       
       do j=1,ndim 
          sz_l(ikp)=sz_l(ikp)-Aimag(ctemp(j,j))
       enddo
       sx_l(ikp)= sx_l(ikp)/dos_l(ikp)
       sy_l(ikp)= sy_l(ikp)/dos_l(ikp)
       sz_l(ikp)= sz_l(ikp)/dos_l(ikp)


       !ctemp=matmul(GRR,sigma_x)       
       call mat_mul(ndim,GRR,sigma_x,ctemp)       
       do j=1,ndim 
          sx_r(ikp)=sx_r(ikp)-Aimag(ctemp(j,j))
       enddo

       !ctemp=matmul(GRR,sigma_y)       
       call mat_mul(ndim,GRR,sigma_y,ctemp)       
       do j=1,ndim 
          sy_r(ikp)=sy_r(ikp)-Aimag(ctemp(j,j))
       enddo

       !ctemp=matmul(GRR,sigma_z)       
       call mat_mul(ndim,GRR,sigma_z,ctemp)       
       do j=1,ndim 
          sz_r(ikp)=sz_r(ikp)-Aimag(ctemp(j,j))
       enddo
       sx_r(ikp)= sx_r(ikp)/dos_r(ikp)
       sy_r(ikp)= sy_r(ikp)/dos_r(ikp)
       sz_r(ikp)= sz_r(ikp)/dos_r(ikp)


    enddo  ! ikp  kpoint
    
#if defined (MPI)
    call mpi_allreduce(sx_l,sx_l_mpi,size(sx_l),mpi_double_precision, &
                  mpi_sum , mpi_comm_world, ierr)
    call mpi_allreduce(sy_l,sy_l_mpi,size(sy_l),mpi_double_precision, &
                  mpi_sum , mpi_comm_world, ierr)
    call mpi_allreduce(sz_l,sz_l_mpi,size(sz_l),mpi_double_precision, &
                  mpi_sum , mpi_comm_world, ierr)
    call mpi_allreduce(dos_l,dos_l_mpi,size(dos_l),mpi_double_precision, &
                  mpi_sum , mpi_comm_world, ierr)
    call mpi_allreduce(sx_r,sx_r_mpi,size(sx_r),mpi_double_precision, &
                  mpi_sum , mpi_comm_world, ierr)
    call mpi_allreduce(sy_r,sy_r_mpi,size(sy_r),mpi_double_precision, &
                  mpi_sum , mpi_comm_world, ierr)
    call mpi_allreduce(sz_r,sz_r_mpi,size(sz_r),mpi_double_precision, &
                  mpi_sum , mpi_comm_world, ierr)
    call mpi_allreduce(dos_r,dos_r_mpi,size(dos_r),mpi_double_precision, &
                  mpi_sum , mpi_comm_world, ierr)
#else
     sx_l_mpi= sx_l
     sy_l_mpi= sy_l
     sz_l_mpi= sz_l
     dos_l_mpi= dos_l
     sx_r_mpi= sx_r
     sy_r_mpi= sy_r
     sz_r_mpi= sz_r
     dos_r_mpi= dos_r
#endif

    sx_l = sx_l_mpi/pi
    sy_l = sy_l_mpi/pi
    sz_l = sz_l_mpi/pi
    dos_l= dos_l_mpi
    sx_r = sx_r_mpi/pi
    sy_r = sy_r_mpi/pi
    sz_r = sz_r_mpi/pi
    dos_r= dos_r_mpi

    do ikp= 1, Nk1*Nk2
       dos_l_only(ikp)= dos_l(ikp)- dos_bulk(ikp)
       if (dos_l_only(ikp)<0) dos_l_only(ikp)=eps9
       dos_r_only(ikp)= dos_r(ikp)- dos_bulk(ikp)
       if (dos_r_only(ikp)<0) dos_r_only(ikp)=eps9
    enddo


    outfileindex= outfileindex+ 1
    spindosfile= outfileindex
    outfileindex= outfileindex+ 1
    spintexturefile= outfileindex
    if (cpuid.eq.0)then
       open(spindosfile,file='spindos_l.dat')
       open(spintexturefile,file='spintext_l.dat')
       write(spindosfile,'(60a16)')'#kx', 'ky', 'log(dos)', &
                           'sx', 'sy', 'sz'
       write(spintexturefile,'(60a16)')'#kx', 'ky', 'log(dos)', &
                           'sx', 'sy', 'sz'

      !ikp= 0
      !do i= 1, nk1
      !   do j= 1, nk2
      !      ikp=ikp+1
      !      k=k12_shape(:, ikp)
      !      write(spindosfile,'(60f16.8)')k, (log(dos_l(ikp))), &
      !         real(sx_l(ikp)), real(sy_l(ikp)), real(sz_l(ikp))                 
      !      if (ikp>1 .and. ikp< (NK1*Nk2-1)) then
      !         if(dos_l(ikp)>dos_l(ikp+1).and.dos_l(ikp)>dos_l(ikp-1).and. &
      !         real(log(dos_l(ikp)))>0.5d0) & 
      !         write(spintexturefile,'(60f16.8)')k, (log(dos_l(ikp))), &
      !            real(sx_l(ikp)),real(sy_l(ikp)),real(sz_l(ikp))                 
      !      endif
      !   enddo
      !   write(spindosfile, *)' '
      !enddo

       dos_l_max= maxval(dos_l_only)/4d0
       ikp= 0
       do i= 1, nk1
          do j= 1, nk2
             ikp=ikp+1
             k=k12_shape(:, ikp)
             write(spindosfile,'(60f16.8)')k, (log(dos_l(ikp))), &
                real(sx_l(ikp)), real(sy_l(ikp)), real(sz_l(ikp))                 
             if (dos_l_only(ikp)>dos_l_max)then
                write(spintexturefile,'(60f16.8)')k, (log(dos_l(ikp))), &
                   real(sx_l(ikp)),real(sy_l(ikp)),real(sz_l(ikp))                 
             endif
          enddo
          write(spindosfile, *)' '
       enddo

       close(spindosfile)
       close(spintexturefile)
    endif

    outfileindex= outfileindex+ 1
    spindosfile= outfileindex
    outfileindex= outfileindex+ 1
    spintexturefile= outfileindex
    if (cpuid.eq.0)then
       open(spindosfile,file='spindos_r.dat')
       open(spintexturefile,file='spintext_r.dat')
       write(spindosfile,'(60a16)')'#kx', 'ky', 'log(dos)', &
                           'sx', 'sy', 'sz'
       write(spintexturefile,'(60a16)')'#kx', 'ky', 'log(dos)', &
                           'sx', 'sy', 'sz'

      !ikp= 0
      !do i= 1, nk1
      !   do j= 1, nk2
      !      ikp=ikp+1
      !      k=k12_shape(:, ikp)
      !      write(spindosfile,'(60f16.8)')k, (log(dos_r(ikp))), &
      !         real(sx_r(ikp)), real(sy_r(ikp)), real(sz_r(ikp))                 
      !      if (ikp>1 .and. ikp< (NK1*Nk2-1)) then
      !         if(dos_r(ikp)>dos_r(ikp+1).and.dos_r(ikp)>dos_r(ikp-1).and. &
      !         real(log(dos_r(ikp)))>0.5d0) & 
      !         write(spintexturefile,'(60f16.8)')k, (log(dos_r(ikp))), &
      !            real(sx_r(ikp)),real(sy_r(ikp)),real(sz_r(ikp))                 
      !      endif
      !   enddo
      !   write(spindosfile, *)' '
      !enddo

       dos_r_max= maxval(dos_r_only)/4d0
       ikp= 0
       do i= 1, nk1
          do j= 1, nk2
             ikp=ikp+1
             k=k12_shape(:, ikp)
             write(spindosfile,'(60f16.8)')k, (log(dos_r(ikp))), &
                real(sx_r(ikp)), real(sy_r(ikp)), real(sz_r(ikp))                 
             if (dos_r_only(ikp)>dos_r_max)then
                write(spintexturefile,'(60f16.8)')k, (log(dos_r(ikp))), &
                   real(sx_r(ikp)),real(sy_r(ikp)),real(sz_r(ikp))                 
             endif
          enddo
          write(spindosfile, *)' '
       enddo


       close(spindosfile)
       close(spintexturefile)
    endif


    !> generate gnuplot scripts for plotting the spin texture
    outfileindex= outfileindex+ 1
    if (cpuid.eq.0) then
       open(outfileindex,file='spintext_r.gnu')
       write(outfileindex, '(a)')"set encoding iso_8859_1"
       write(outfileindex, '(a)')'#set terminal pngcairo truecolor enhanced font ",80" size 3680, 3360'
       write(outfileindex, '(a)')'set terminal png truecolor enhanced font ",80" size 3680, 3360'
       write(outfileindex, '(a)')"set output 'spintext_r.png'"
       write(outfileindex, '(a)')'set palette defined ( -6 "white", 0 "gray", 10 "blue" )'
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
       write(outfileindex, '(a)')"splot 'spindos_r.dat' u 1:2:3 w pm3d"
       write(outfileindex, '(a)')"unset xtics"
       write(outfileindex, '(a)')"unset ytics"
       write(outfileindex, '(a)')"unset xlabel"
       write(outfileindex, '(a)')"unset ylabel"
       write(outfileindex, '(a)')"set label 1 'Spin texture' at graph 0.25, 1.10 front"
       write(outfileindex, '(a)')"# the sencond plot"
       write(outfileindex, '(a)')"set origin 0.14, 0.115"
       write(outfileindex, '(a)')"set size 0.75, 0.70"
       write(outfileindex, '(a)')"plot 'spintext_r.dat' u 1:2:($4/5.00):($5/5.00)  w vec  head lw 4 lc rgb 'red' front"

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
       write(outfileindex, '(a)')'set palette defined ( -6 "white", 0 "gray", 10 "blue" )'
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
       write(outfileindex, '(a)')"splot 'spindos_l.dat' u 1:2:3 w pm3d"
       write(outfileindex, '(a)')"unset xtics"
       write(outfileindex, '(a)')"unset ytics"
       write(outfileindex, '(a)')"unset xlabel"
       write(outfileindex, '(a)')"unset ylabel"
       write(outfileindex, '(a)')"set label 1 'Spin texture' at graph 0.25, 1.10 front"
       write(outfileindex, '(a)')"# the sencond plot"
       write(outfileindex, '(a)')"set origin 0.14, 0.115"
       write(outfileindex, '(a)')"set size 0.75, 0.70"
       write(outfileindex, '(a)')"plot 'spintext_l.dat' u 1:2:($4/5.00):($5/5.00)  w vec  head lw 4 lc rgb 'red' front"

       close(outfileindex)
    endif
  
    if (cpuid.eq.0)write(stdout,*)'calculate spintexture successfully' 
 
    return
  end subroutine spintext
