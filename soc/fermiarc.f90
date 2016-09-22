!+---------+---------+---------+---------+---------+---------+--------+!
! this subroutine is used to calculate surface state using             !
! green's function method  ---  J.Phys.F.Met.Phys.15(1985)851-858      !
! 
! History:
!         by Quan Sheng Wu on 4/20/2010                                !
!            mpi version      4/21/2010
!         Ca3PbO version     11/12/2011
!+---------+---------+---------+---------+---------+---------+--------+!
  subroutine fermiarc


     use wmpi
     use para
     implicit none
     


     integer :: ierr

     ! general loop index
     integer :: i,j 

     ! kpoint loop index
     integer :: ikp


     real(dp) :: k(2)

     real(dp) :: omega
     real(dp) :: k1min_shape, k1max_shape, k2min_shape, k2max_shape

     real(dp), allocatable :: k12(:,:)
     real(dp), allocatable :: k12_shape(:,:)
     
     real(dp), allocatable :: dos_l(:)
     real(dp), allocatable :: dos_l_mpi(:)
     real(dp), allocatable :: dos_r(:)
     real(dp), allocatable :: dos_r_mpi(:)

     complex(dp), allocatable :: ones(:,:)
     complex(dp), allocatable :: GLL(:,:), GRR(:,:)
     complex(dp), allocatable :: H00(:,:), H01(:,:)


     allocate( k12(2, nk1*nk2))
     allocate( k12_shape(2, nk1*nk2))
     allocate( dos_l(nk1*nk2))
     allocate( dos_l_mpi(nk1*nk2))
     allocate( dos_r(nk1*nk2))
     allocate( dos_r_mpi(nk1*nk2))
     allocate( GLL(ndim,ndim), GRR(ndim,ndim))
     k12=0d0
     k12_shape=0d0
     dos_l=0d0
     dos_l_mpi=0d0
     dos_r=0d0
     dos_r_mpi=0d0

     ikp=0
     do i= 1, nk1
        do j= 1, nk2
           ikp=ikp+1
           k12(:, ikp)=K2D_start+ (i-1)*K2D_vec1/dble(nk1-1) &
                      + (j-1)*K2D_vec2/dble(nk2-1)
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
     GLL= 0d0
     GRR= 0d0
     H00= 0d0
     H01= 0d0
     ones= 0d0

     do i=1,Ndim
        ones(i,i)=1.0d0
     enddo


     omega = E_arc
     eta= eta_arc

     do ikp= 1+cpuid, nk1*nk2, num_cpu
        if (cpuid==0) write(stdout, *) 'Arc, ik, knv2', ikp, nk1*nk2
        k(1)= k12(1, ikp)
        k(2)= k12(2, ikp)

        call ham_qlayer2qlayer(k,H00,H01) 

        !> calculate surface green function
        ! there are two method to calculate surface green's function 
        ! the method in 1985 is better, you can find the ref in the
        ! subroutine
        call surfgreen_1985(omega,GLL,GRR,H00,H01,ones)
        ! call surfgreen_1984(omega,GLL,GRR,H00,H01,ones)


        ! calculate spectral function
        do i= 1, ndim
           dos_l(ikp)=dos_l(ikp)- aimag(GLL(i,i))
           dos_r(ikp)=dos_r(ikp)- aimag(GRR(i,i))
        enddo
     enddo

     !> we don't have to do allreduce operation
     call mpi_reduce(dos_l, dos_l_mpi, size(dos_l),mpi_double_precision,&
                     mpi_sum, 0, mpi_comm_world, ierr)
     call mpi_reduce(dos_r, dos_r_mpi, size(dos_r),mpi_double_precision,&
                     mpi_sum, 0, mpi_comm_world, ierr)

     if (cpuid.eq.0)then
        open (unit=12, file='arc.dat_l')
        open (unit=13, file='arc.dat_r')
        do ikp=1, nk1*nk2
           write(12, '(3f16.8)')k12_shape(:, ikp), log(dos_l_mpi(ikp))
           if (mod(ikp, nk2)==0) write(12, *)' '
           write(13, '(3f16.8)')k12_shape(:, ikp), log(dos_r_mpi(ikp))
           if (mod(ikp, nk2)==0) write(13, *)' '
        enddo
        close(12)
        close(13)
        write(stdout,*)'ndim',ndim
        write(stdout,*) 'Nk1,Nk2,eta',Nk1, Nk2, eta
        write(stdout,*)'calculate density of state successfully'    
     endif

     !> write script for gnuplot
     if (cpuid==0) then
        open(unit=114, file='arc_l.gnu')
        write(114, '(a)')"set encoding iso_8859_1"
        write(114, '(a)')'#set terminal  postscript enhanced color'
        write(114, '(a)')"#set output 'arc_l.eps'"
        write(114, '(3a)')'#set terminal pngcairo truecolor enhanced', &
           '  font ",60" size 1920, 1680'
        write(114, '(3a)')'set terminal png truecolor enhanced', &
           ' font ",60" size 1920, 1680'
        write(114, '(a)')"set output 'arc_l.png'"
        write(114,'(2a)') 'set palette defined ( -10 "#194eff", ', &
           '0 "white", 10 "red" )'
        write(114, '(a)')'#set palette rgbformulae 33,13,10'
        write(114, '(a)')'unset ztics'
        write(114, '(a)')'unset key'
        write(114, '(a)')'set pm3d'
        write(114, '(a)')'set border lw 6'
        write(114, '(a)')'set size ratio -1'
        write(114, '(a)')'set view map'
        write(114, '(a)')'set xtics'
        write(114, '(a)')'set ytics'
        write(114, '(a)')'set xlabel "K_1 (1/{\305})"'
        write(114, '(a)')'set ylabel "K_2 (1/{\305})"'
        write(114, '(a)')'set colorbox'
        write(114, '(a, f8.5, a, f8.5, a)')'set xrange [', k1min_shape, ':', k1max_shape, ']'
        write(114, '(a, f8.5, a, f8.5, a)')'set yrange [', k2min_shape, ':', k2max_shape, ']'
        write(114, '(a)')'set pm3d interpolate 2,2'
        write(114, '(2a)')"splot 'arc.dat_l' u 1:2:3 w pm3d"

        close(114)
     endif

     !> write script for gnuplot
     if (cpuid==0) then
        open(unit=115, file='arc_r.gnu')
        write(115, '(a)')"set encoding iso_8859_1"
        write(115, '(a)')'#set terminal  postscript enhanced color'
        write(115, '(a)')"#set output 'arc_r.eps'"
        write(115, '(3a)')'#set terminal pngcairo truecolor enhanced', &
           '  font ",60" size 1920, 1680'
        write(115, '(3a)')'set terminal png truecolor enhanced', &
           ' font ",60" size 1920, 1680'
        write(115, '(a)')"set output 'arc_r.png'"
        write(115,'(2a)') 'set palette defined ( -10 "#194eff", ', &
           '0 "white", 10 "red" )'
        write(115, '(a)')'#set palette rgbformulae 33,13,10'
        write(115, '(a)')'unset ztics'
        write(115, '(a)')'unset key'
        write(115, '(a)')'set pm3d'
        write(115, '(a)')'set border lw 6'
        write(115, '(a)')'set size ratio -1'
        write(115, '(a)')'set view map'
        write(115, '(a)')'set xtics'
        write(115, '(a)')'set ytics'
        write(115, '(a)')'set xlabel "K_1 (1/{\305})"'
        write(115, '(a)')'set ylabel "K_2 (1/{\305})"'
        write(115, '(a)')'set colorbox'
        write(115, '(a, f8.5, a, f8.5, a)')'set xrange [', k1min_shape, ':', k1max_shape, ']'
        write(115, '(a, f8.5, a, f8.5, a)')'set yrange [', k2min_shape, ':', k2max_shape, ']'
        write(115, '(a)')'set pm3d interpolate 2,2'
        write(115, '(2a)')"splot 'arc.dat_r' u 1:2:3 w pm3d"

        close(115)
     endif


  return   
  end subroutine fermiarc


  !> calculation joint density of state 
  !> refs :http://science.sciencemag.org/content/sci/suppl/2016/03/09/351.6278.1184.DC1/Inoue.SM.pdf

  !> in this version, we also calculate the spin density jdos
  subroutine fermiarc_jdos


     use wmpi
     use para
     implicit none
     


     integer :: ierr

     ! general loop index
     integer :: i,j 

     integer :: Nwann
     ! kpoint loop index
     integer :: ikp, ik1, ik2, iq, Nk1_half, Nk2_half
     integer :: imin1, imax1, imin2, imax2, iq1, iq2

     real(dp) :: k(2), dis

     real(dp) :: omega
     real(dp) :: k1min_shape, k1max_shape, k2min_shape, k2max_shape

     integer , allocatable :: ik12(:,:)
     real(dp), allocatable :: k12(:,:)
     real(dp), allocatable :: k12_shape(:,:)
     
     real(dp), allocatable :: dos_l(:,:)
     real(dp), allocatable :: dos_l_mpi(:,:)
     real(dp), allocatable :: dos_r(:,:)
     real(dp), allocatable :: dos_r_mpi(:,:)
     real(dp), allocatable :: jdos_l(:)
     real(dp), allocatable :: jdos_l_mpi(:)
     real(dp), allocatable :: jdos_r(:)
     real(dp), allocatable :: jdos_r_mpi(:)
     real(dp), allocatable :: jsdos_l(:)
     real(dp), allocatable :: jsdos_l_mpi(:)
     real(dp), allocatable :: jsdos_r(:)
     real(dp), allocatable :: jsdos_r_mpi(:)
     real(dp), allocatable :: sx_l(:, :)
     real(dp), allocatable :: sy_l(:, :)
     real(dp), allocatable :: sz_l(:, :)
     real(dp), allocatable :: sx_l_mpi(:, :)
     real(dp), allocatable :: sy_l_mpi(:, :)
     real(dp), allocatable :: sz_l_mpi(:, :)
     real(dp), allocatable :: sx_r(:, :)
     real(dp), allocatable :: sy_r(:, :)
     real(dp), allocatable :: sz_r(:, :)
     real(dp), allocatable :: sx_r_mpi(:, :)
     real(dp), allocatable :: sy_r_mpi(:, :)
     real(dp), allocatable :: sz_r_mpi(:, :)

     ! spin operator matrix sigma_x,sigma_y in sigma_z representation
     complex(Dp),allocatable :: sigma_x(:,:) 
     complex(Dp),allocatable :: sigma_y(:,:) 
     complex(Dp),allocatable :: sigma_z(:,:) 
     complex(Dp),allocatable :: ctemp(:,:)
     complex(dp), allocatable :: ones(:,:)
     complex(dp), allocatable :: GLL(:,:), GRR(:,:)
     complex(dp), allocatable :: H00(:,:), H01(:,:)

     Nk1_half= (Nk1-1)/2
     Nk2_half= (Nk2-1)/2
     Nwann= Num_wann/2

     allocate( ik12(2, nk1*nk2))
     allocate( k12(2, nk1*nk2))
     allocate( k12_shape(2, nk1*nk2))
     allocate( dos_l(nk1,nk2))
     allocate( dos_l_mpi(nk1,nk2))
     allocate( dos_r(nk1,nk2))
     allocate( dos_r_mpi(nk1,nk2))
     allocate( jdos_l(nk1*nk2))
     allocate( jdos_l_mpi(nk1*nk2))
     allocate( jdos_r(nk1*nk2))
     allocate( jdos_r_mpi(nk1*nk2))
     allocate( jsdos_l(nk1*nk2))
     allocate( jsdos_l_mpi(nk1*nk2))
     allocate( jsdos_r(nk1*nk2))
     allocate( jsdos_r_mpi(nk1*nk2))
     allocate( sx_l(nk1,nk2))
     allocate( sx_l_mpi(nk1,nk2))
     allocate( sy_l(nk1,nk2))
     allocate( sy_l_mpi(nk1,nk2))
     allocate( sz_l(nk1,nk2))
     allocate( sz_l_mpi(nk1,nk2))
     allocate( sx_r(nk1,nk2))
     allocate( sx_r_mpi(nk1,nk2))
     allocate( sy_r(nk1,nk2))
     allocate( sy_r_mpi(nk1,nk2))
     allocate( sz_r(nk1,nk2))
     allocate( sz_r_mpi(nk1,nk2))
     allocate( GLL(ndim,ndim), GRR(ndim,ndim))
     allocate(sigma_x(ndim,ndim))
     allocate(sigma_y(ndim,ndim))
     allocate(sigma_z(ndim,ndim))
     allocate(ctemp(ndim,ndim))

     sigma_x   =0.0d0
     sigma_y   =0.0d0
     sigma_z   =0.0d0
     ik12=0
     k12=0d0
     k12_shape=0d0
     dos_l=0d0
     dos_l_mpi=0d0
     dos_r=0d0
     dos_r_mpi=0d0
     jdos_l=0d0
     jdos_l_mpi=1d-12
     jdos_r=0d0
     jdos_r_mpi=1d-12
     jsdos_l=0d0
     jsdos_l_mpi=1d-12
     jsdos_r=0d0
     jsdos_r_mpi=1d-12
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

     ikp=0
     do i= 1, nk1
        do j= 1, nk2
           ikp=ikp+1
           ik12(1, ikp)= i
           ik12(2, ikp)= j
           k12(:, ikp)=K2D_start+ (i-1)*K2D_vec1/dble(nk1-1) &
                      + (j-1)*K2D_vec2/dble(nk2-1)
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
     GLL= 0d0
     GRR= 0d0
     H00= 0d0
     H01= 0d0
     ones= 0d0

     do i=1,Ndim
        ones(i,i)=1.0d0
     enddo

     Nwann= Num_wann/2
 
     !> spin operator matrix
     do i=1, Np
        do j=1, Nwann
           sigma_x(Num_wann*(i-1)+j, Num_wann*(i-1)+Nwann+j)=1.0d0
           sigma_x(Num_wann*(i-1)+j+Nwann, Num_wann*(i-1)+j)=1.0d0
           sigma_y(Num_wann*(i-1)+j, Num_wann*(i-1)+Nwann+j)=-zi
           sigma_y(Num_wann*(i-1)+j+Nwann, Num_wann*(i-1)+j)=zi
           sigma_z(Num_wann*(i-1)+j, Num_wann*(i-1)+j)= 1d0
           sigma_z(Num_wann*(i-1)+j+Nwann, Num_wann*(i-1)+j+Nwann)=-1d0
        enddo
     enddo



     omega = E_arc
     eta= eta_arc

     do ikp= 1+cpuid, nk1*nk2, num_cpu
        if (cpuid==0) write(stdout, *) 'Arc, ik, knv2', ikp, nk1*nk2
        k(1)= k12(1, ikp)
        k(2)= k12(2, ikp)

        call ham_qlayer2qlayer(k,H00,H01) 

        !> calculate surface green function
        ! there are two method to calculate surface green's function 
        ! the method in 1985 is better, you can find the ref in the
        ! subroutine
        call surfgreen_1985(omega,GLL,GRR,H00,H01,ones)
        ! call surfgreen_1984(omega,GLL,GRR,H00,H01,ones)

        ik1= ik12(1, ikp)
        ik2= ik12(2, ikp)

        ! calculate spectral function
        do i= 1, ndim
           dos_l(ik1, ik2)=dos_l(ik1, ik2)- aimag(GLL(i,i))
           dos_r(ik1, ik2)=dos_r(ik1, ik2)- aimag(GRR(i,i))
        enddo
        !! calculate spectral function
        !do i= 1, Nwann
        !   dos_l(ik1, ik2)=dos_l(ik1, ik2)- aimag(GLL(i,i)) &
        !   - aimag(GLL(i+Nwann,i+Nwann))
        !enddo ! i
        !do i=Ndim-Num_wann+1, Ndim- Nwann
        !   dos_r(ik1, ik2)=dos_r(ik1, ik2)- aimag(GRR(i,i)) &
        !   - aimag(GRR(i+Nwann,i+Nwann))
        !enddo

        !ctemp=matmul(surfgreen,sigma_x)       
        call mat_mul(ndim,GLL,sigma_x,ctemp)       
        do j=1,ndim 
           sx_l(ik1, ik2)=sx_l(ik1, ik2)-Aimag(ctemp(j,j))
        enddo
  
        !ctemp=matmul(surfgreen,sigma_y)       
        call mat_mul(ndim,GLL,sigma_y,ctemp)       
        do j=1,ndim 
           sy_l(ik1, ik2)=sy_l(ik1, ik2)-Aimag(ctemp(j,j))
        enddo
  
        !ctemp=matmul(surfgreen,sigma_z)       
        call mat_mul(ndim,GLL,sigma_z,ctemp)       
        do j=1,ndim 
           sz_l(ik1, ik2)=sz_l(ik1, ik2)-Aimag(ctemp(j,j))
        enddo


        !ctemp=matmul(surfgreen,sigma_x)       
        call mat_mul(ndim,GRR,sigma_x,ctemp)       
        do j=1,ndim 
           sx_r(ik1, ik2)=sx_r(ik1, ik2)-Aimag(ctemp(j,j))
        enddo
  
        !ctemp=matmul(surfgreen,sigma_y)       
        call mat_mul(ndim,GRR,sigma_y,ctemp)       
        do j=1,ndim 
           sy_r(ik1, ik2)=sy_r(ik1, ik2)-Aimag(ctemp(j,j))
        enddo
  
        !ctemp=matmul(surfgreen,sigma_z)       
        call mat_mul(ndim,GRR,sigma_z,ctemp)       
        do j=1,ndim 
           sz_r(ik1, ik2)=sz_r(ik1, ik2)-Aimag(ctemp(j,j))
        enddo
     enddo

     !> we don't have to do allreduce operation
     call mpi_allreduce(dos_l, dos_l_mpi, size(dos_l),mpi_double_precision,&
                     mpi_sum, mpi_comm_world, ierr)
     call mpi_allreduce(dos_r, dos_r_mpi, size(dos_r),mpi_double_precision,&
                     mpi_sum, mpi_comm_world, ierr)
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

     if (cpuid.eq.0)then
        open (unit=12, file='arc.dat_l')
        open (unit=13, file='arc.dat_r')
        open(30014,file='spindos_l.dat')
        open(30015,file='spindos_r.dat')
        do ikp=1, nk1*nk2
           write(12, '(3f16.8)')k12_shape(:, ikp), log(dos_l_mpi(ik12(1, ikp), ik12(2, ikp)))
           if (mod(ikp, nk2)==0) write(12, *)' '
           write(13, '(3f16.8)')k12_shape(:, ikp), log(dos_r_mpi(ik12(1, ikp), ik12(2, ikp)))
           if (mod(ikp, nk2)==0) write(13, *)' '
           write(30014, '(30f16.8)')k12_shape(:, ikp), &
                                     log(dos_l_mpi(ik12(1, ikp), ik12(2, ikp))), &
                                     (sx_l_mpi(ik12(1, ikp), ik12(2, ikp)))/dos_l_mpi(ik12(1, ikp), ik12(2, ikp)), &
                                     (sy_l_mpi(ik12(1, ikp), ik12(2, ikp)))/dos_l_mpi(ik12(1, ikp), ik12(2, ikp)), &
                                     (sz_l_mpi(ik12(1, ikp), ik12(2, ikp)))/dos_l_mpi(ik12(1, ikp), ik12(2, ikp))
           write(30015, '(30f16.8)')k12_shape(:, ikp), &
                                     log(dos_r_mpi(ik12(1, ikp), ik12(2, ikp))), &
                                     (sx_r_mpi(ik12(1, ikp), ik12(2, ikp)))/dos_r_mpi(ik12(1, ikp), ik12(2, ikp)), &
                                     (sy_r_mpi(ik12(1, ikp), ik12(2, ikp)))/dos_r_mpi(ik12(1, ikp), ik12(2, ikp)), &
                                     (sz_r_mpi(ik12(1, ikp), ik12(2, ikp)))/dos_r_mpi(ik12(1, ikp), ik12(2, ikp))
           if (mod(ikp, nk2)==0) write(30014, *)' '
           if (mod(ikp, nk2)==0) write(30015, *)' '
        enddo
        close(12)
        close(13)
        close(30014)
        close(30015)
        write(stdout,*)'ndim',ndim
        write(stdout,*) 'Nk1,Nk2,eta',Nk1, Nk2, eta
        write(stdout,*)'calculate density of state successfully'    
     endif



     !> calculate jdos
     do iq= 1+ cpuid, Nk1*Nk2, num_cpu
        iq1= ik12(1, iq)- Nk1_half
        iq2= ik12(2, iq)- Nk2_half
        if (cpuid==0) write(stdout, *) 'Jdos, iq, knv2', iq, nk1*nk2
        imin1= max(-Nk1_half-iq1, -Nk1_half)+ Nk1_half+ 1
        imax1= min(Nk1_half-iq1, Nk1_half)+ Nk1_half+ 1
        imin2= max(-Nk2_half-iq2, -Nk2_half)+ Nk2_half+ 1
        imax2= min(Nk2_half-iq2, Nk2_half)+ Nk2_half+ 1
        do ik1= imin1, imax1
           do ik2= imin2, imax2
              jdos_l(iq)= jdos_l(iq)+ dos_l_mpi(ik1, ik2)* dos_l_mpi(ik1+iq1, ik2+iq2)
              jdos_r(iq)= jdos_r(iq)+ dos_r_mpi(ik1, ik2)* dos_r_mpi(ik1+iq1, ik2+iq2)
              jsdos_l(iq)= jsdos_l(iq)+ dos_l_mpi(ik1, ik2)* dos_l_mpi(ik1+iq1, ik2+iq2) &
                                      + sx_l_mpi(ik1, ik2)* sx_l_mpi(ik1+iq1, ik2+iq2) &
                                      + sy_l_mpi(ik1, ik2)* sy_l_mpi(ik1+iq1, ik2+iq2) &
                                      + sz_l_mpi(ik1, ik2)* sz_l_mpi(ik1+iq1, ik2+iq2)
              jsdos_r(iq)= jsdos_r(iq)+ dos_r_mpi(ik1, ik2)* dos_r_mpi(ik1+iq1, ik2+iq2) &
                                      + sx_r_mpi(ik1, ik2)* sx_r_mpi(ik1+iq1, ik2+iq2) &
                                      + sy_r_mpi(ik1, ik2)* sy_r_mpi(ik1+iq1, ik2+iq2) &
                                      + sz_r_mpi(ik1, ik2)* sz_r_mpi(ik1+iq1, ik2+iq2)
           enddo !ik
        enddo !ik
     enddo !ik

     jdos_l_mpi=1d-12
     jdos_r_mpi=1d-12
     jsdos_l_mpi=1d-12
     jsdos_r_mpi=1d-12
     call mpi_reduce(jdos_l, jdos_l_mpi, size(jdos_l),mpi_double_precision,&
                     mpi_sum, 0, mpi_comm_world, ierr)
     call mpi_reduce(jdos_r, jdos_r_mpi, size(jdos_r),mpi_double_precision,&
                     mpi_sum, 0, mpi_comm_world, ierr)
     call mpi_reduce(jsdos_l, jsdos_l_mpi, size(jsdos_l),mpi_double_precision,&
                     mpi_sum, 0, mpi_comm_world, ierr)
     call mpi_reduce(jsdos_r, jsdos_r_mpi, size(jsdos_r),mpi_double_precision,&
                     mpi_sum, 0, mpi_comm_world, ierr)

     if (cpuid.eq.0)then
        open (unit=1212, file='arc.jdat_l')
        open (unit=1313, file='arc.jdat_r')
        open (unit=1214, file='arc.jsdat_l')
        open (unit=1315, file='arc.jsdat_r')
        do ikp=1, nk1*nk2
           write(1212, '(30f16.8)')k12_shape(:, ikp), log(jdos_l_mpi(ikp))
           if (mod(ikp, nk2)==0) write(1212, *)' '
           write(1313, '(30f16.8)')k12_shape(:, ikp), log(jdos_r_mpi(ikp))
           if (mod(ikp, nk2)==0) write(1313, *)' '
           write(1214, '(30f16.8)')k12_shape(:, ikp), log(jsdos_l_mpi(ikp))
           if (mod(ikp, nk2)==0) write(1214, *)' '
           write(1315, '(30f16.8)')k12_shape(:, ikp), log(jsdos_r_mpi(ikp))
           if (mod(ikp, nk2)==0) write(1315, *)' '
        enddo
        close(1212)
        close(1313)
        close(1214)
        close(1315)
        write(stdout,*)'calculate joint density of state successfully'    
     endif

     !> write script for gnuplot
     if (cpuid==0) then
        open(unit=114, file='arc_l.gnu')
        write(114, '(a)')"set encoding iso_8859_1"
        write(114, '(a)')'#set terminal  postscript enhanced color'
        write(114, '(a)')"#set output 'arc_l.eps'"
        write(114, '(3a)')'#set terminal pngcairo truecolor enhanced', &
           '  font ",60" size 1920, 1680'
        write(114, '(3a)')'set terminal png truecolor enhanced', &
           ' font ",60" size 1920, 1680'
        write(114, '(a)')"set output 'arc_l.png'"
        write(114,'(2a)') 'set palette defined ( -10 "#194eff", ', &
           '0 "white", 10 "red" )'
        write(114, '(a)')'#set palette rgbformulae 33,13,10'
        write(114, '(a)')'unset ztics'
        write(114, '(a)')'unset key'
        write(114, '(a)')'set pm3d'
        write(114, '(a)')'set border lw 6'
        write(114, '(a)')'set size ratio -1'
        write(114, '(a)')'set view map'
        write(114, '(a)')'set xtics'
        write(114, '(a)')'set ytics'
        write(114, '(a)')'set xlabel "K_1 (1/{\305})"'
        write(114, '(a)')'set ylabel "K_2 (1/{\305})"'
        write(114, '(a)')'set colorbox'
        write(114, '(a, f8.5, a, f8.5, a)')'set xrange [', k1min_shape, ':', k1max_shape, ']'
        write(114, '(a, f8.5, a, f8.5, a)')'set yrange [', k2min_shape, ':', k2max_shape, ']'
        write(114, '(a)')'set pm3d interpolate 2,2'
        write(114, '(2a)')"splot 'arc.dat_l' u 1:2:3 w pm3d"

        close(114)
     endif

     !> write script for gnuplot
     if (cpuid==0) then
        open(unit=115, file='arc_r.gnu')
        write(115, '(a)')"set encoding iso_8859_1"
        write(115, '(a)')'#set terminal  postscript enhanced color'
        write(115, '(a)')"#set output 'arc_r.eps'"
        write(115, '(3a)')'#set terminal pngcairo truecolor enhanced', &
           '  font ",60" size 1920, 1680'
        write(115, '(3a)')'set terminal png truecolor enhanced', &
           ' font ",60" size 1920, 1680'
        write(115, '(a)')"set output 'arc_r.png'"
        write(115,'(2a)') 'set palette defined ( -10 "#194eff", ', &
           '0 "white", 10 "red" )'
        write(115, '(a)')'#set palette rgbformulae 33,13,10'
        write(115, '(a)')'unset ztics'
        write(115, '(a)')'unset key'
        write(115, '(a)')'set pm3d'
        write(115, '(a)')'set border lw 6'
        write(115, '(a)')'set size ratio -1'
        write(115, '(a)')'set view map'
        write(115, '(a)')'set xtics'
        write(115, '(a)')'set ytics'
        write(115, '(a)')'set xlabel "K_1 (1/{\305})"'
        write(115, '(a)')'set ylabel "K_2 (1/{\305})"'
        write(115, '(a)')'set colorbox'
        write(115, '(a, f8.5, a, f8.5, a)')'set xrange [', k1min_shape, ':', k1max_shape, ']'
        write(115, '(a, f8.5, a, f8.5, a)')'set yrange [', k2min_shape, ':', k2max_shape, ']'
        write(115, '(a)')'set pm3d interpolate 2,2'
        write(115, '(2a)')"splot 'arc.dat_r' u 1:2:3 w pm3d"

        close(115)
     endif


  return   
  end subroutine fermiarc_jdos



  !> calculation joint density of state 
  !> refs :http://science.sciencemag.org/content/sci/suppl/2016/03/09/351.6278.1184.DC1/Inoue.SM.pdf
  !> read dos_l.dat dos_r.dat from the file
  subroutine fermiarc_jdos2


     use wmpi
     use para
     implicit none
     


     integer :: ierr

     ! general loop index
     integer :: i,j 

     integer :: Nwann
     ! kpoint loop index
     integer :: ikp, ik1, ik2, iq, Nk1_half, Nk2_half
     integer :: imin1, imax1, imin2, imax2, iq1, iq2

     real(dp) :: k(2), dis
     real(dp) :: r1, r2

     real(dp) :: omega
     real(dp) :: k1min_shape, k1max_shape, k2min_shape, k2max_shape

     integer , allocatable :: ik12(:,:)
     real(dp), allocatable :: k12(:,:)
     real(dp), allocatable :: k12_shape(:,:)
     
     real(dp), allocatable :: dos_l(:,:)
     real(dp), allocatable :: dos_l_mpi(:,:)
     real(dp), allocatable :: dos_r(:,:)
     real(dp), allocatable :: dos_r_mpi(:,:)
     real(dp), allocatable :: jdos_l(:)
     real(dp), allocatable :: jdos_l_mpi(:)
     real(dp), allocatable :: jdos_r(:)
     real(dp), allocatable :: jdos_r_mpi(:)

     complex(dp), allocatable :: ones(:,:)
     complex(dp), allocatable :: GLL(:,:), GRR(:,:)
     complex(dp), allocatable :: H00(:,:), H01(:,:)

     Nk1_half= (Nk1-1)/2
     Nk2_half= (Nk2-1)/2
     Nwann= Num_wann/2

     allocate( ik12(2, nk1*nk2))
     allocate( k12(2, nk1*nk2))
     allocate( k12_shape(2, nk1*nk2))
     allocate( dos_l(nk1,nk2))
     allocate( dos_l_mpi(nk1,nk2))
     allocate( dos_r(nk1,nk2))
     allocate( dos_r_mpi(nk1,nk2))
     allocate( jdos_l(nk1*nk2))
     allocate( jdos_l_mpi(nk1*nk2))
     allocate( jdos_r(nk1*nk2))
     allocate( jdos_r_mpi(nk1*nk2))
     allocate( GLL(ndim,ndim), GRR(ndim,ndim))
     ik12=0
     k12=0d0
     k12_shape=0d0
     dos_l=0d0
     dos_l_mpi=0d0
     dos_r=0d0
     dos_r_mpi=0d0
     jdos_l=0d0
     jdos_l_mpi=1d-12
     jdos_r=0d0
     jdos_r_mpi=1d-12

     ikp=0
     do i= 1, nk1
        do j= 1, nk2
           ikp=ikp+1
           ik12(1, ikp)= i
           ik12(2, ikp)= j
           k12(:, ikp)=K2D_start+ (i-1)*K2D_vec1/dble(nk1-1) &
                      + (j-1)*K2D_vec2/dble(nk2-1)
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
     GLL= 0d0
     GRR= 0d0
     H00= 0d0
     H01= 0d0
     ones= 0d0

     do i=1,Ndim
        ones(i,i)=1.0d0
     enddo


     omega = E_arc
     eta= eta_arc

     if (cpuid==0) then
     !> read arc.dat_l and arc.dat_r 
     open(unit=13523, file='arc.dat_l')
     open(unit=13524, file='arc.dat_r')
     do ikp= 1, nk1*nk2
        if (mod(ikp, nk2)==0) write(stdout, *) 'Arc, ik, knv2', ikp, nk1*nk2
        k(1)= k12(1, ikp)
        k(2)= k12(2, ikp)

        read(13523, *)r1, r2, dos_l(ik12(1, ikp), ik12(2, ikp))
        if (mod(ikp, nk2)==0)read(13523, *)
        read(13524, *)r1, r2, dos_r(ik12(1, ikp), ik12(2, ikp))
        if (mod(ikp, nk2)==0)read(13524, *)
        dos_l= dexp(dos_l)
        dos_r= dexp(dos_r)

     enddo
     close(13523)
     close(13524)

     endif
     call MPI_bcast(dos_l,size(dos_l),mpi_dp,0,mpi_cmw,ierr)
     call MPI_bcast(dos_r,size(dos_r),mpi_dp,0,mpi_cmw,ierr)

     !> calculate jdos
     do iq= 1+ cpuid, Nk1*Nk2, num_cpu
        iq1= ik12(1, iq)- Nk1_half
        iq2= ik12(2, iq)- Nk2_half
        if (cpuid==0) write(stdout, *) 'Jdos, iq, knv2', iq, nk1*nk2
        imin1= max(-Nk1_half-iq1, -Nk1_half)+ Nk1_half+ 1
        imax1= min(Nk1_half-iq1, Nk1_half)+ Nk1_half+ 1
        imin2= max(-Nk2_half-iq2, -Nk2_half)+ Nk2_half+ 1
        imax2= min(Nk2_half-iq2, Nk2_half)+ Nk2_half+ 1
        do ik1= imin1, imax1
           do ik2= imin2, imax2
              jdos_l(iq)= jdos_l(iq)+ dos_l(ik1, ik2)* dos_l(ik1+iq1, ik2+iq2)
              jdos_r(iq)= jdos_r(iq)+ dos_r(ik1, ik2)* dos_r(ik1+iq1, ik2+iq2)
           enddo !ik
        enddo !ik
     enddo !ik

     jdos_l_mpi=0d0  
     jdos_r_mpi=0d0  
     call mpi_reduce(jdos_l, jdos_l_mpi, size(jdos_l),mpi_double_precision,&
                     mpi_sum, 0, mpi_comm_world, ierr)
     call mpi_reduce(jdos_r, jdos_r_mpi, size(jdos_r),mpi_double_precision,&
                     mpi_sum, 0, mpi_comm_world, ierr)

     if (cpuid.eq.0)then
        open (unit=1212, file='arc.jdat_l')
        open (unit=1313, file='arc.jdat_r')
        do ikp=1, nk1*nk2
           write(1212, '(2f10.5, 30f36.5)')k12_shape(:, ikp), (jdos_l_mpi(ikp))
           if (mod(ikp, nk2)==0) write(1212, *)' '
           write(1313, '(2f10.5, 30f36.5)')k12_shape(:, ikp), (jdos_r_mpi(ikp))
           if (mod(ikp, nk2)==0) write(1313, *)' '
        enddo
        close(1212)
        close(1313)
        write(stdout,*)'calculate joint density of state successfully'    
     endif

     !> write script for gnuplot
     if (cpuid==0) then
        open(unit=114, file='arc_l.gnu')
        write(114, '(a)')"set encoding iso_8859_1"
        write(114, '(a)')'#set terminal  postscript enhanced color'
        write(114, '(a)')"#set output 'arc_l.eps'"
        write(114, '(3a)')'#set terminal pngcairo truecolor enhanced', &
           '  font ",60" size 1920, 1680'
        write(114, '(3a)')'set terminal png truecolor enhanced', &
           ' font ",60" size 1920, 1680'
        write(114, '(a)')"set output 'arc_l.png'"
        write(114,'(2a)') 'set palette defined ( -10 "#194eff", ', &
           '0 "white", 10 "red" )'
        write(114, '(a)')'#set palette rgbformulae 33,13,10'
        write(114, '(a)')'unset ztics'
        write(114, '(a)')'unset key'
        write(114, '(a)')'set pm3d'
        write(114, '(a)')'set border lw 6'
        write(114, '(a)')'set size ratio -1'
        write(114, '(a)')'set view map'
        write(114, '(a)')'set xtics'
        write(114, '(a)')'set ytics'
        write(114, '(a)')'set xlabel "K_1 (1/{\305})"'
        write(114, '(a)')'set ylabel "K_2 (1/{\305})"'
        write(114, '(a)')'set colorbox'
        write(114, '(a, f8.5, a, f8.5, a)')'set xrange [', k1min_shape, ':', k1max_shape, ']'
        write(114, '(a, f8.5, a, f8.5, a)')'set yrange [', k2min_shape, ':', k2max_shape, ']'
        write(114, '(a)')'set pm3d interpolate 2,2'
        write(114, '(2a)')"splot 'arc.dat_l' u 1:2:3 w pm3d"

        close(114)
     endif

     !> write script for gnuplot
     if (cpuid==0) then
        open(unit=115, file='arc_r.gnu')
        write(115, '(a)')"set encoding iso_8859_1"
        write(115, '(a)')'#set terminal  postscript enhanced color'
        write(115, '(a)')"#set output 'arc_r.eps'"
        write(115, '(3a)')'#set terminal pngcairo truecolor enhanced', &
           '  font ",60" size 1920, 1680'
        write(115, '(3a)')'set terminal png truecolor enhanced', &
           ' font ",60" size 1920, 1680'
        write(115, '(a)')"set output 'arc_r.png'"
        write(115,'(2a)') 'set palette defined ( -10 "#194eff", ', &
           '0 "white", 10 "red" )'
        write(115, '(a)')'#set palette rgbformulae 33,13,10'
        write(115, '(a)')'unset ztics'
        write(115, '(a)')'unset key'
        write(115, '(a)')'set pm3d'
        write(115, '(a)')'set border lw 6'
        write(115, '(a)')'set size ratio -1'
        write(115, '(a)')'set view map'
        write(115, '(a)')'set xtics'
        write(115, '(a)')'set ytics'
        write(115, '(a)')'set xlabel "K_1 (1/{\305})"'
        write(115, '(a)')'set ylabel "K_2 (1/{\305})"'
        write(115, '(a)')'set colorbox'
        write(115, '(a, f8.5, a, f8.5, a)')'set xrange [', k1min_shape, ':', k1max_shape, ']'
        write(115, '(a, f8.5, a, f8.5, a)')'set yrange [', k2min_shape, ':', k2max_shape, ']'
        write(115, '(a)')'set pm3d interpolate 2,2'
        write(115, '(2a)')"splot 'arc.dat_r' u 1:2:3 w pm3d"

        close(115)
     endif


  return   
  end subroutine fermiarc_jdos2

