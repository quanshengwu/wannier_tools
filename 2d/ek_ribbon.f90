!> This subroutine is used for calculating energy 
!> dispersion with wannier functions
! Copyright (c) 2010 QuanSheng Wu. All rights reserved.
  subroutine ek_ribbon
    
     use wmpi
     use para
     implicit none 

! loop index
     integer :: i     
     integer :: j
     integer :: l

     integer :: lwork

! wave vector 
     real(Dp) :: k

     real(dp) :: emin, emax
 
! parameters for zheev
     integer :: ierr
     real(Dp), allocatable ::  rwork(:)
     complex(Dp), allocatable :: work(:)
      
! eigenvalue 
     real(Dp), allocatable :: eigenvalue(:)
   
! energy dispersion
     real(Dp),allocatable :: ekribbon(:,:)
     real(Dp),allocatable :: ekribbon_mpi(:,:)


     !> color for plot, surface state weight
     real(dp), allocatable :: surf_weight(:, :)
     real(dp), allocatable :: surf_weight_mpi(:, :)

! hamiltonian ribbon
     complex(Dp),allocatable ::CHamk(:,:)

     lwork= 16*Nslab*Num_wann
     ierr = 0


     allocate(eigenvalue(nslab*Num_wann))
     allocate( surf_weight (Nslab* Num_wann, knv1))
     allocate( surf_weight_mpi (Nslab* Num_wann, knv1))
     allocate(ekribbon(Nslab*Num_wann,knv1))
     allocate(ekribbon_mpi(Nslab*Num_wann,knv1))
     allocate(CHamk(nslab*Num_wann,nslab*Num_wann))
     allocate(work(lwork))
     allocate(rwork(lwork))
 
     surf_weight= 0d0
     surf_weight_mpi= 0d0

     ! sweep k
     ekribbon=0.0d0
     ekribbon_mpi=0.0d0
     do i=1+cpuid,knv1,num_cpu
        if (cpuid==0) write(stdout, *), 'ribbonEk, ik ',  i, knv1
        k= k1_path(i)
        chamk=0.0d0 

        !> no magnetgic field
        if (abs(Bx)<eps9.and. abs(By)<eps9.and. abs(Bz)<eps9)then
           call ham_ribbon(k,Chamk)
        !> in-plane magnetic field
        elseif (abs(Bx)>eps9 .or. abs(By)>eps9)then
           write(stdout, *)'>> magnetic field is larger than zero'
           call ham_ribbon(k,Chamk)
          !call ham_slab_parallel_B(k,Chamk)
        !> vertical magnetic field
        else
           print *, 'Error: we only support in-plane magnetic field at present'
           stop 'please set Bz= 0'
        endif


        eigenvalue=0.0d0

        ! diagonal Chamk
        call eigensystem_c('V', 'U', Num_wann*Nslab, CHamk, eigenvalue)
       
        ekribbon(:,i)=eigenvalue

        do j=1, Nslab*Num_wann
           do l=1, Num_wann
              surf_weight(j, i)= surf_weight(j, i) &
                 + abs(CHamk(l, j))**2 & ! first slab
                 + abs(CHamk(Num_wann*Nslab- l+ 1, j))**2 !& ! last slab
                !+ abs(CHamk(Num_wann+ l, j))**2 & ! the second slab
                !+ abs(CHamk(Num_wann*(Nslab-1)- l, j))**2 ! last second slab
           enddo ! l
           surf_weight(j, i)= sqrt(surf_weight(j, i))
        enddo ! j 
     enddo ! i

     call mpi_allreduce(ekribbon,ekribbon_mpi,size(ekribbon),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
     call mpi_allreduce(surf_weight, surf_weight_mpi,size(surf_weight),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
 

     ekribbon=ekribbon_mpi
     if (maxval(surf_weight_mpi)<0.00001d0)surf_weight_mpi=1d0
     surf_weight= surf_weight_mpi/ maxval(surf_weight_mpi)
     

     if(cpuid==0)then
        open(unit=100, file='ribbonek.dat')
        do j=1, Num_wann*Nslab
           do i=1, knv2
             !write(100,'(3f15.7, i8)')k1len(i), ekslab(j,i), &
             !   (surf_weight(j, i))
              write(100,'(2f15.7, i8)')k1len(i), ekribbon(j,i), &
                 int(255-surf_weight(j, i)*255d0)
           enddo
           write(100 , *)''
        enddo
        close(100)
        write(stdout,*) 'calculate energy band  done'
     endif

     emin= minval(ekribbon)-0.5d0
     emax= maxval(ekribbon)+0.5d0
     !> write script for gnuplot
     if (cpuid==0) then
        open(unit=113, file='ribbonek.gnu')
        write(113, '(a)')"set encoding iso_8859_1"
        write(113, '(a)')'#set terminal  postscript enhanced color'
        write(113, '(a)')"#set output 'ribbonek.eps'"
        write(113, '(3a)')'#set terminal  pngcairo truecolor enhanced', &
           '  font ",60" size 1920, 1680'
        write(113, '(3a)')'set terminal  png truecolor enhanced', &
           '  font ",60" size 1920, 1680'
        write(113, '(a)')"set output 'ribbonek.png'"
        write(113,'(2a)') 'set palette defined ( 0  "green", ', &
           '5 "yellow", 10 "red" )'
        write(113, '(a)')'set style data linespoints'
        write(113, '(a)')'unset ztics'
        write(113, '(a)')'unset key'
        write(113, '(a)')'set pointsize 0.8'
        write(113, '(a)')'set border lw 3 '
        write(113, '(a)')'set view 0,0'
        write(113, '(a)')'#set xtics font ",36"'
        write(113, '(a)')'#set ytics font ",36"'
        write(113, '(a)')'#set ylabel font ",36"'
        write(113, '(a)')'#set xtics offset 0, -1'
        write(113, '(a)')'set ylabel offset -1, 0 '
        write(113, '(a)')'set ylabel "Energy (eV)"'
        write(113, '(a, f10.5, a)')'set xrange [0: ', maxval(k1len), ']'
        write(113, '(a, f10.5, a, f10.5, a)')'set yrange [', emin, ':', emax, ']'
        write(113, 202, advance="no") (trim(k1line_name(i)), k1line_stop(i), i=1, nk1lines)
        write(113, 203)trim(k1line_name(nk1lines+1)), k1line_stop(nk1lines+1)

        do i=1, nk1lines-1
           write(113, 204)k1line_stop(i+1), emin, k1line_stop(i+1), emax
        enddo
        write(113, '(a)')'rgb(r,g,b) = int(r)*65536 + int(g)*256 + int(b)'
        write(113, '(2a)')"plot 'ribbonek.dat' u 1:2:(rgb(255,$3, 3)) ",  &
            "w lp lw 2 pt 7  ps 1 lc rgb variable"

        !write(113, '(2a)')"splot 'ribbonek.dat' u 1:2:3 ",  &
        !   "w lp lw 2 pt 13 palette"
        close(113)
     endif

     202 format('set xtics (',:20('"',A3,'" ',F10.5,','))
     203 format(A3,'" ',F10.5,')')
     204 format('set arrow from ',F10.5,',',F10.5, &
        ' to ',F10.5,',',F10.5, ' nohead')
     
  return
  end subroutine ek_ribbon

