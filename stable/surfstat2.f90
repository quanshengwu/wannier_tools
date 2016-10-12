!+---------+---------+---------+---------+---------+---------+--------+!
! this subroutine is used to calculate surface state using             !
! green's function method  ---  J.Phys.F.Met.Phys.15(1985)851-858      !
! 
! History:
!         by Quan Sheng Wu on 4/20/2010                                !
!            mpi version      4/21/2010
!            change Kb to K=(Ka+Kb)/3 direction 4/22/2010
!            Quansheng Wu on Jan 30 2015 at ETH Zurich
! Copyright (c) 2010 QuanSheng Wu. All rights reserved.
!+---------+---------+---------+---------+---------+---------+--------+!
  subroutine surfstat2

     use wmpi
     use para
     implicit none
     
     integer :: ierr

     ! general loop index
     integer :: i,j 

     integer :: knv2

     ! kpoint loop index
     integer :: ikp

     integer :: NN, nlines

     integer :: Nwann

     real(dp) :: emin
     real(dp) :: emax
     real(dp) :: t1, temp
     real(dp) :: rx0, ry0, x0, y0
     real(dp) :: k(2), w, s
     real(dp) :: k1(2)
     real(dp) :: k2(2)

     real(dp) :: kp(16, 2)
     real(dp) :: ke(16, 2)
     real(dp) :: kpath_stop(16)
     character(4) :: kpath_name(17)

     real(dp) :: kstart(2), kend(2)

     real(dp), allocatable :: kpoint(:,:)
     real(dp), allocatable :: omega(:)
     
     real(dp), allocatable :: dos_l(:, :,:)
     real(dp), allocatable :: dos_r(:, :,:)
     real(dp), allocatable :: dos_l_mpi(:, :,:)
     real(dp), allocatable :: dos_r_mpi(:, :,:)

     complex(dp), allocatable :: GLL(:,:)
     complex(dp), allocatable :: GRR(:,:)
     complex(dp), allocatable :: H00(:,:)
     complex(dp), allocatable :: H01(:,:)
     complex(dp), allocatable :: ones(:,:)

     real(dp), allocatable :: k_len(:)

     kpath_name= ' '
     kp(1,:)=(/0.5d0, 0.00d0/)  ; kpath_name(1)= 'X'
     ke(1,:)=(/0.0d0, 0.00d0/)  
     kp(2,:)=(/0.0d0, 0.00d0/)  ; kpath_name(2)= 'G'
     ke(2,:)=(/0.0d0, 0.50d0/)  ! K
     kp(3,:)=(/0.0d0, 0.50d0/) ; kpath_name(3)= 'Y'     
     ke(3,:)=(/0.5d0, 0.5d0/)  ! K
     kp(4,:)=(/0.5d0, 0.5d0/)  ; kpath_name(4)= 'M'     
     ke(4,:)=(/0.0d0, 0.0d0/)  ; kpath_name(5)= 'G'  

    !kp(1,:)=(/0.5d0, 0.00d0/)  ; kpath_name(1)= 'X'
    !ke(1,:)=(/0.0d0, 0.00d0/)  
    !kp(1,:)=(/0.266666d0, 0.2666660d0/)  ; kpath_name(1)= 'G'
    !ke(1,:)=(/0.311111d0, 0.3111111d0/)  ! K
    !kp(2,:)=(/0.333333d0, 0.3333330d0/) ; kpath_name(2)= 'K'     
    !ke(2,:)=(/0.5d0, 0.0d0/)  ! K
    !kp(3,:)=(/0.5d0, 0.0d0/)  ; kpath_name(3)= 'Y'     
    !ke(3,:)=(/0.333333d0, 0.3333330d0/)  ! K
    !kp(5,:)=(/0.333333d0, 0.3333330d0/)  ! K
    !ke(5,:)=(/0.0d0, 0.0d0/)  ; kpath_name(5)= 'G'  


    !kp(1,:)=(/0.2d0, 0.2d0/)  ! Gamma
    !ke(1,:)=(/0.5d0, 0.5d0/)  ! Z

    !kp(1,:)=(/-0.5d0, 0.00d0/)  ; kpath_name(1)= 'A'
    !ke(1,:)=(/0.0d0, 0.00d0/)  
    !kp(2,:)=(/0.0d0, 0.00d0/)  ; kpath_name(2)= 'G'
    !ke(2,:)=(/0.5d0, 0.00d0/)  ; kpath_name(3)= 'A'
   
    !kp(1,:)=(/0.3d0, 0.00d0/)  ; kpath_name(1)= 'X'
    !ke(1,:)=(/0.0d0, 0.00d0/)  
    !kp(1, 1)= surf_onsite
    !ke(1, 2)= surf_onsite
    !kp(2, :)=(/0.0d0, 0.00d0/)  ; kpath_name(2)= 'G'

    !kp(1,:)=(/-0.5d0, 0.00d0/)  ; kpath_name(1)= 'T'
    !ke(1,:)=(/0.0d0, 0.20d0/)  
    !kp(2,:)=(/0.0d0, 0.20d0/)  ; kpath_name(2)= 'P'
    !ke(2,:)=(/0.5d0, 0.50d0/)  ; kpath_name(3)= 'M'

    !> 0.19 0.195  0.2  0.21  0.23  0.25  0.26   0.29
    ! 0.105 0.1076 0.11 0.116 0.127 0.138 0.1435 0.16
    !kp(1,:)=(/ 0.1380d0,-0.18d0/)  ; kpath_name(1)= 'T'
    !ke(1,:)=(/ 0.1380d0, 0.00d0/)  
    !kp(2,:)=(/ 0.1380d0, 0.00d0/)  ; kpath_name(2)= 'G'
    !ke(2,:)=(/ 0.1380d0, 0.18d0/)  

     kp(1,:)=(/ 0.0000d0, 0.00d0/)  ; kpath_name(1)= 'X'
     ke(1,:)=(/ 0.3000d0, 0.00d0/)  

     nlines=1
     NN= Nk
     knv2=NN*nlines
     allocate( kpoint(knv2, 2))
     allocate( k_len (knv2))
     kpoint= 0d0

     t1=0d0
     k_len=0d0
     kpath_stop= 0d0
     do j=1, nlines 
        do i=1, NN
           kstart= kp(j,:)
           kend  = ke(j,:)
           k1= kstart(1)*Ka2+ kstart(2)*Kb2
           k2= kend(1)*Ka2+ kend(2)*Kb2
           kpoint(i+(j-1)*NN,:)= kstart+ (kend-kstart)*(i-1)/dble(NN-1)
           
           temp= dsqrt((k2(1)- k1(1))**2 &
                 +(k2(2)- k1(2))**2)/dble(NN-1) 

           if (i.gt.1) then
              t1=t1+temp
           endif
           k_len(i+(j-1)*NN)= t1
        enddo
        kpath_stop(j+1)= t1
     enddo

     !> lines across two weyl points
    !k1= -0.5d0*ka2+0.5d0*kb2
    !t1= 0d0
    !do i=1, NN
    !   s= (i-1d0)/(NN-1d0)-0.5d0
    !   k2= k1
    !   kpoint(i, 1)= s
    !   kpoint(i, 2)= -2.818d0*s**3-0.2955*s
    !   k1= kpoint(i, 1)*ka2+ kpoint(i, 2)*kb2
    !   temp= dsqrt((k2(1)- k1(1))**2 &
    !         +(k2(2)- k1(2))**2)

    !   if (i.gt.1) then
    !      t1=t1+temp
    !   endif
    !   k_len(i)= t1
    !enddo

     !> k lines around one weyl points
     !> for WTe2, W1 (0.121480, 0.0447), W2= 0.121831, 0.0388)
     !> W0 (0.1216, 0.042)

     !> for MoTe2 W1 (0.10296,  0.05388) E=0.06eV
    !x0= 0.10296d0
    !y0= 0.05388d0

     !>           W2 (0.10516,  0.01378) E=0.005eV
    !x0= 0.10516d0
    !y0= 0.01378d0

    !rx0= 0.001d0
    !ry0= 0.005d0
    !s=0
    !k2(1)= x0+ rx0*cos(s+pi)
    !k2(2)= y0 + ry0*sin(s+pi)
    !k1= k2(1)*ka2+ k2(2)*kb2
    !t1=0d0
    !do i=1, NN
    !   k2= k1
    !   s= (i-1d0)/(NN-1d0)*2d0*pi+pi
    !   kpoint(i, 1)= x0 + rx0*cos(s)
    !   kpoint(i, 2)= y0 + ry0*sin(s)
    !   k1= kpoint(i, 1)*ka2+ kpoint(i, 2)*kb2
    !   temp= dsqrt((k2(1)- k1(1))**2 &
    !       + (k2(2)- k1(2))**2)

    !   if (i.gt.1) then
    !      t1=t1+temp
    !   endif
    !   k_len(i)= t1
    !enddo


     Nwann= Num_wann/2 
     if (SOC ==0 ) Nwann= Num_wann
     allocate( omega(omeganum))
     allocate( dos_l(Nwann, knv2, omeganum ))
     allocate( dos_r(Nwann, knv2, omeganum ))
     allocate( dos_l_mpi(Nwann, knv2, omeganum ))
     allocate( dos_r_mpi(Nwann, knv2, omeganum ))
     omega=0d0
     dos_l=0d0
     dos_r=0d0
     dos_l_mpi=0d0
     dos_r_mpi=0d0

     eta=(omegamax- omegamin)/dble(omeganum)*3.0d0

     do i= 1, omeganum
        omega(i)=omegamin+(i-1)*(omegamax-omegamin)/dble(omeganum)
     enddo

     allocate(GLL(Ndim, Ndim))
     allocate(GRR(Ndim, Ndim))
     allocate(H00(Ndim, Ndim))
     allocate(H01(Ndim, Ndim))
     allocate(ones(Ndim, Ndim))
     GLL= 0d0
     GRR= 0d0
     H00= 0d0
     H01= 0d0
     ones= 0d0

     do i=1,Ndim
        ones(i,i)=1.0d0
     enddo

     do ikp= 1+cpuid, knv2, num_cpu
        if (cpuid==0) write(*, *) ikp, 'ik', knv2
        if (cpuid==0) write(stdout, *) ikp, 'ik', knv2
        k= kpoint(ikp,:)

        !> get the hopping matrix between two principle layers
        call ham_qlayer2qlayer(k,H00,H01) 

        do j = 1, omeganum
           w=omega(j)

           !> calculate surface green function
           ! there are two method to calculate surface green's function 
           ! the method in 1985 is better, you can find the ref in the
           ! subroutine
            call surfgreen_1985(w,GLL,GRR,H00,H01,ones)
           !call surfgreen_1984(w,GLL,GRR,H00,H01,ones)

           ! calculate spectral function
           do i= 1, Nwann
              dos_l(i, ikp, j)=dos_l(i, ikp, j)- aimag(GLL(i,i)) &
                                        - aimag(GLL(i+Nwann,i+Nwann))
           enddo ! i
           do i=Ndim-Num_wann+1, Ndim- Nwann
              dos_r(i, ikp, j)=dos_r(i, ikp, j)- aimag(GRR(i,i)) &
                                        - aimag(GRR(i+Nwann,i+Nwann))
           enddo ! i
        enddo ! j
     enddo ! ikp

     !> we don't have to do allreduce operation
     call mpi_reduce(dos_l, dos_l_mpi, size(dos_l), mpi_double_precision,&
                     mpi_sum, 0, mpi_comm_world, ierr)
     call mpi_reduce(dos_r, dos_r_mpi, size(dos_r), mpi_double_precision,&
                     mpi_sum, 0, mpi_comm_world, ierr)
    !dos_l=log(dos_l_mpi)
    !dos_r=log(dos_r_mpi)
     dos_l= dos_l_mpi/maxval(dos_l_mpi)*255d0
     dos_r= dos_r_mpi/maxval(dos_r_mpi)*255d0

     if (cpuid.eq.0)then
        open (unit=16, file='dos.dat_l')
        open (unit=13, file='dos.dat_r')
        do ikp=1, knv2
           do j=1, omeganum 
              write(16, '(2f16.8, 30000i5)')k_len(ikp), omega(j), int8(dos_l(:, ikp, j))
              write(13, '(2f16.8, 30000i5)')k_len(ikp), omega(j), int8(dos_r(:, ikp, j))
           enddo
           write(16, *) ' '
           write(13, *) ' '
        enddo
        close(16)
        close(13)
        write(stdout,*)'ndim',ndim
        write(stdout,*) 'knv2,omeganum,eta',knv2, omeganum, eta
        write(stdout,*)'calculate density of state successfully'    
     endif

     emin= minval(omega)
     emax= maxval(omega)
     !> write script for gnuplot
     if (cpuid==0) then
        open(unit=104, file='surfdos_l.gnu')
        write(104, '(a)')'#set terminal  postscript enhanced color'
        write(104, '(a)')"#set output 'surfdos_l.eps'"
        write(104, '(3a)')'set terminal  pngcairo truecolor enhanced', &
           ' font ", 36" size 1920, 1680'
        write(104, '(a)')"set output 'surfdos_l.png'"
        write(104,'(2a)') 'set palette defined (-20 "#194eff", ', &
           '0 "white", 10 "red" )'
        write(104, '(a)')'#set palette rgbformulae 33,13,10'
        write(104, '(a)')'set style data linespoints'
        write(104, '(a)')'#set size ratio -1'
        write(104, '(a)')'unset ztics'
        write(104, '(a)')'unset key'
        write(104, '(a)')'set pointsize 0.8'
        write(104, '(a)')'set pm3d'
        write(104, '(a)')'#set view equal xyz'
        write(104, '(a)')'set view map'
        write(104, '(a)')'set border lw 3'
        write(104, '(a)')'set cbtics font ",48"'
        write(104, '(a)')'set xtics font ",48"'
        write(104, '(a)')'set ytics font ",48"'
        write(104, '(a)')'set ylabel font ",48"'
        write(104, '(a)')'set ylabel "Energy (eV)"'
        write(104, '(a)')'#set xtics offset 0, -1'
        write(104, '(a)')'#set ylabel offset -6, 0 '
        write(104, '(a, f8.5, a)')'set xrange [0: ', maxval(k_len), ']'
        write(104, '(a, f8.5, a, f8.5, a)')'set yrange [', emin, ':', emax, ']'
        write(104, 202, advance="no") (kpath_name(i), kpath_stop(i), i=1, nlines)
        write(104, 203)kpath_name(nlines+1), kpath_stop(nlines+1)

        do i=1, nlines-1
           write(104, 204)kpath_stop(i+1), emin, kpath_stop(i+1), emax
        enddo
        write(104, '(a)')'set pm3d interpolate 2,2'
        write(104, '(2a)')"splot 'dos.dat_l' u 1:2:3 w pm3d"

     endif

     !> write script for gnuplot
     if (cpuid==0) then
        open(unit=105, file='surfdos_r.gnu')
        write(105, '(a)')'#set terminal  postscript enhanced color'
        write(105, '(a)')"#set output 'surfdos_r.eps'"
        write(105, '(3a)')'set terminal  pngcairo truecolor enhanced', &
           ' font ", 36" size 1920, 1680'
        write(105, '(a)')"set output 'surfdos_r.png'"
        write(105,'(2a)') 'set palette defined (-20 "#194eff", ', &
           '0 "white", 10 "red" )'
        write(105, '(a)')'#set palette rgbformulae 33,13,10'
        write(105, '(a)')'set style data linespoints'
        write(105, '(a)')'unset ztics'
        write(105, '(a)')'unset key'
        write(105, '(a)')'set pointsize 0.8'
        write(105, '(a)')'set pm3d'
        write(105, '(a)')'set border lw 3'
        write(105, '(a)')'#set size ratio -1'
        write(105, '(a)')'#set view equal xyz'
        write(105, '(a)')'set view map'
        write(105, '(a)')'set cbtics font ",48"'
        write(105, '(a)')'set xtics font ",48"'
        write(105, '(a)')'set ytics font ",48"'
        write(105, '(a)')'set ylabel font ",48"'
        write(105, '(a)')'set ylabel "Energy (eV)"'
        write(105, '(a)')'#set xtics offset 0, -1'
        write(105, '(a)')'#set ylabel offset -6, 0 '
        write(105, '(a, f8.5, a)')'set xrange [0: ', maxval(k_len), ']'
        write(105, '(a, f8.5, a, f8.5, a)')'set yrange [', emin, ':', emax, ']'
        write(105, 202, advance="no") (kpath_name(i), kpath_stop(i), i=1, nlines)
        write(105, 203)kpath_name(nlines+1), kpath_stop(nlines+1)

        do i=1, nlines-1
           write(105, 204)kpath_stop(i+1), emin, kpath_stop(i+1), emax
        enddo
        write(105, '(a)')'set pm3d interpolate 2,2'
        write(105, '(2a)')"splot 'dos.dat_r' u 1:2:3 w pm3d"

     endif


     202 format('set xtics (',:20('"',A1,'" ',F8.5,','))
     203 format(A1,'" ',F8.5,')')
     204 format('set arrow from ',F8.5,',',F10.5,' to ',F8.5,',',F10.5, ' nohead front lw 3')
 

  return   
  end subroutine
