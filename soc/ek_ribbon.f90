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
     integer :: info,ierr,ifail
      
   
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
