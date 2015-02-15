!+---------+---------+---------+---------+---------+---------+--------+!
! this subroutine is used to calculate surface state using             !
! green's function method  ---  J.Phys.F.Met.Phys.15(1985)851-858      !
! 
! History:
!         by Quan Sheng Wu on 4/20/2010                                !
!            mpi version      4/21/2010
!            change Kb to K=(Ka+Kb)/3 direction 4/22/2010
!            Quansheng Wu on Jan 30 2015 at ETH Zurich
!+---------+---------+---------+---------+---------+---------+--------+!
  subroutine surfstat

     use mpi
     use para
     implicit none
     
     integer :: ierr

     ! general loop index
     integer :: i,j 

     integer :: knv2

     ! kpoint loop index
     integer :: ikp

     integer :: NN, nlines

     real(dp) :: emin
     real(dp) :: emax
     real(dp) :: t1, temp
     real(dp) :: k(2), w
     real(dp) :: k1(2)
     real(dp) :: k2(2)

     real(dp) :: kp(16, 2)
     real(dp) :: ke(16, 2)
     real(dp) :: kpath_stop(16)
     character(4) :: kpath_name(17)

     real(dp) :: kstart(2), kend(2)

     real(dp), allocatable :: kpoint(:,:)
     real(dp), allocatable :: omega(:)
     
     real(dp), allocatable :: dos_l(:,:)
     real(dp), allocatable :: dos_r(:,:)
     real(dp), allocatable :: dos_l_mpi(:,:)
     real(dp), allocatable :: dos_r_mpi(:,:)

     complex(dp), allocatable :: GLL(:,:)
     complex(dp), allocatable :: GRR(:,:)
     complex(dp), allocatable :: H00(:,:)
     complex(dp), allocatable :: H01(:,:)
     complex(dp), allocatable :: ones(:,:)

     real(dp), allocatable :: k_len(:)

     kpath_name= ' '
     kp(1,:)=(/0.0d0, 0.5d0/)  ; kpath_name(1)= 'Y'
     ke(1,:)=(/0.0d0, 0.0d0/)  
     kp(2,:)=(/0.0d0, 0.0d0/)  ; kpath_name(2)= 'G'
     ke(2,:)=(/0.5d0, 0.00d0/)  ! K
     kp(3,:)=(/0.5d0, 0.00d0/) ; kpath_name(3)= 'X'     
     ke(3,:)=(/0.5d0, 0.5d0/)  ! K
     kp(4,:)=(/0.5d0, 0.5d0/)  ; kpath_name(4)= 'M'     
     ke(4,:)=(/0.0d0, 0.0d0/)  ; kpath_name(5)= 'G'  

     kp(5,:)=(/0.0d0, 0.0d0/)  ! Gamma
     kp(6,:)=(/5.0d0, 0.0d0/)  ! Z
    
     nlines=4
     NN=80
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


     allocate( omega(omeganum))
     allocate( dos_l(knv2, omeganum))
     allocate( dos_r(knv2, omeganum))
     allocate( dos_l_mpi(knv2, omeganum))
     allocate( dos_r_mpi(knv2, omeganum))
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
        if (cpuid==0) write(stdout, *) ikp, 'in', knv2
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
           ! call surfgreen_1984(w,GLL,GRR,H00,H01,ones)

           ! calculate spectral function
           do i= 1, ndim
              dos_l(ikp, j)=dos_l(ikp,j)- aimag(GLL(i,i))
              dos_r(ikp, j)=dos_r(ikp,j)- aimag(GRR(i,i))
           enddo ! i
        enddo ! j
     enddo ! ikp

     !> we don't have to do allreduce operation
     call mpi_reduce(dos_l, dos_l_mpi, size(dos_l), mpi_double_precision,&
                     mpi_sum, 0, mpi_comm_world, ierr)
     call mpi_reduce(dos_r, dos_r_mpi, size(dos_r), mpi_double_precision,&
                     mpi_sum, 0, mpi_comm_world, ierr)
     dos_l=log(dos_l_mpi)
     dos_r=log(dos_r_mpi)

     if (cpuid.eq.0)then
        open (unit=12, file='dos.dat_l')
        open (unit=13, file='dos.dat_r')
        do ikp=1, knv2
           do j=1, omeganum 
              write(12, '(3f16.8)')k_len(ikp), omega(j), dos_l(ikp, j)
              write(13, '(3f16.8)')k_len(ikp), omega(j), dos_r(ikp, j)
           enddo
           write(12, *) ' '
           write(13, *) ' '
        enddo
        close(12)
        close(13)
        write(stdout,*)'ndim',ndim
        write(stdout,*) 'knv2,omeganum,eta',knv2, omeganum, eta
        write(stdout,*)'calculate density of state successfully'    
     endif

     emin= minval(omega)
     emax= maxval(omega)
     !> write script for gnuplot
     if (cpuid==0) then
        open(unit=101, file='surfdos_l.gnu')
        write(101, '(a)')'#set terminal  postscript enhanced color'
        write(101, '(a)')"#set output 'surfdos_l.eps'"
        write(101, '(a)')'set terminal  png truecolor enhanced transarpent giant'
        write(101, '(a)')"set output 'surfdos_l.png'"
        write(101,'(2a)') '#set palette defined (-10 "green", ', &
           '0 "yellow", 10 "red" )'
        write(101, '(a)')'set palette rgbformulae 33,13,10'
        write(101, '(a)')'set style data linespoints'
        write(101, '(a)')'unset ztics'
        write(101, '(a)')'unset key'
        write(101, '(a)')'set pointsize 0.8'
        write(101, '(a)')'set pm3d'
        write(101, '(a)')'#set view equal xyz'
        write(101, '(a)')'set view map'
        write(101, '(a)')'set xtics font ",24"'
        write(101, '(a)')'set ytics font ",24"'
        write(101, '(a)')'set ylabel font ",24"'
        write(101, '(a)')'set ylabel "Energy (eV)"'
        write(101, '(a, f8.5, a)')'set xrange [0: ', maxval(k_len), ']'
        write(101, '(a, f8.5, a, f8.5, a)')'set yrange [', emin, ':', emax, ']'
        write(101, 202, advance="no") (kpath_name(i), kpath_stop(i), i=1, nlines)
        write(101, 203)kpath_name(nlines+1), kpath_stop(nlines+1)

        do i=1, nlines-1
           write(101, 204)kpath_stop(i+1), emin, kpath_stop(i+1), emax
        enddo
        write(101, '(a)')'set pm3d interpolate 2,2'
        write(101, '(2a)')"splot 'dos.dat_l' u 1:2:3 w pm3d"

     endif

     !> write script for gnuplot
     if (cpuid==0) then
        open(unit=101, file='surfdos_r.gnu')
        write(101, '(a)')'#set terminal  postscript enhanced color'
        write(101, '(a)')"#set output 'surfdos_r.eps'"
        write(101, '(a)')'set terminal  png truecolor enhanced transarpent giant'
        write(101, '(a)')"set output 'surfdos_r.png'"
        write(101,'(2a)') '#set palette defined (-10 "green", ', &
           '0 "yellow", 10 "red" )'
        write(101, '(a)')'set palette rgbformulae 33,13,10'
        write(101, '(a)')'set style data linespoints'
        write(101, '(a)')'unset ztics'
        write(101, '(a)')'unset key'
        write(101, '(a)')'set pointsize 0.8'
        write(101, '(a)')'set pm3d'
        write(101, '(a)')'#set view equal xyz'
        write(101, '(a)')'set view map'
        write(101, '(a)')'set xtics font ",24"'
        write(101, '(a)')'set ytics font ",24"'
        write(101, '(a)')'set ylabel font ",24"'
        write(101, '(a)')'set ylabel "Energy (eV)"'
        write(101, '(a, f8.5, a)')'set xrange [0: ', maxval(k_len), ']'
        write(101, '(a, f8.5, a, f8.5, a)')'set yrange [', emin, ':', emax, ']'
        write(101, 202, advance="no") (kpath_name(i), kpath_stop(i), i=1, nlines)
        write(101, 203)kpath_name(nlines+1), kpath_stop(nlines+1)

        do i=1, nlines-1
           write(101, 204)kpath_stop(i+1), emin, kpath_stop(i+1), emax
        enddo
        write(101, '(a)')'set pm3d interpolate 2,2'
        write(101, '(2a)')"splot 'dos.dat_r' u 1:2:3 w pm3d"

     endif


     202 format('set xtics (',:20('"',A3,'" ',F8.5,','))
     203 format(A3,'" ',F8.5,')')
     204 format('set arrow from ',F8.5,',',F10.5,' to ',F8.5,',',F10.5, ' nohead')
 

  return   
  end subroutine
