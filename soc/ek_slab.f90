! This subroutine is used to caculate energy dispersion with wannier functions

  subroutine ek_slab
    
     use mpi
     use para
     implicit none 

! loop index
     integer :: i     
     integer :: j
     integer :: l
     integer :: NN, nlines, knv3

     integer :: lwork

! wave vector 
     real(dp) :: k1(2)
     real(dp) :: k2(2)
     real(Dp) :: k(2), kstart(2), kend(2)
     real(dp) :: kp(16,2), ke(16,2), kpath_stop(16)
     character(4) :: kpath_name(17)

     real(dp) :: cell_volume
     real(dp) :: t1, temp
     real(dp) :: emin, emax
 
! parameters for zheev
     integer :: info,ierr
     real(Dp), allocatable ::  rwork(:)
     complex(Dp), allocatable :: work(:)
      
! eigenvalue 
     real(Dp) :: eigenvalue(nslab*Num_wann)
   
! energy dispersion
     real(Dp),allocatable :: ekslab(:,:)
     real(Dp),allocatable :: ekslab_mpi(:,:)

     real(dp), allocatable :: kpoint(:,:)
     real(dp), allocatable :: k_len(:)

     !> color for plot, surface state weight
     real(dp), allocatable :: surf_weight(:, :)
     real(dp), allocatable :: surf_weight_mpi(:, :)

! hamiltonian slab
     complex(Dp),allocatable ::CHamk(:,:)

     lwork= 16*Nslab*Num_wann
     ierr = 0

     kpath_name= ' '
     kp(1,:)=(/0.0d0, 0.5d0/)  ; kpath_name(1)= 'K'
     ke(1,:)=(/0.0d0, 0.0d0/)  
     kp(2,:)=(/0.0d0, 0.0d0/)  ; kpath_name(2)= 'G'
     ke(2,:)=(/0.5d0, 0.00d0/)  ! K
     kp(3,:)=(/0.5d0, 0.00d0/) ; kpath_name(3)= 'K'     
     ke(3,:)=(/0.5d0, 0.5d0/)  ! K
     kp(4,:)=(/0.5d0, 0.5d0/)  ; kpath_name(4)= 'M'     
     ke(4,:)=(/0.0d0, 0.0d0/)  ; kpath_name(5)= 'G'  

     kp(5,:)=(/0.0d0, 0.0d0/)  ! K
     ke(5,:)=(/0.5d0, 0.5d0/)  ! K
     kp(6,:)=(/0.0d0, 0.0d0/)  ! K
     ke(6,:)=(/1.0d0, 1.0d0/)  ! K
  


     nlines= 4
     NN=20
     knv3=NN*nlines
     allocate( kpoint(knv3, 3))
     allocate( k_len (knv3))
     allocate( surf_weight (Nslab* Num_wann, knv3))
     allocate( surf_weight_mpi (Nslab* Num_wann, knv3))
     allocate(ekslab(Nslab*Num_wann,knv3))
     allocate(ekslab_mpi(Nslab*Num_wann,knv3))
     allocate(CHamk(nslab*Num_wann,nslab*Num_wann))
     allocate(work(lwork))
     allocate(rwork(lwork))
 
     kpoint= 0d0
     surf_weight= 0d0
     surf_weight_mpi= 0d0

     t1=0d0
     k_len=0d0
     kpath_stop= 0d0
     do j=1, nlines 
        do i=1, NN
           k = kp(j,:)
           kstart=k
           k1= kstart(1)*Ka2+ kstart(2)*Kb2

           k = ke(j,:)
           kend=k
           k2= kend(1)*Ka2+ kend(2)*Kb2

           kpoint(i+(j-1)*NN,:)= kstart+ (kend-kstart)*dble(i-1)/dble(NN-1)
           
           temp= dsqrt((k2(1)- k1(1))**2 &
                 +(k2(2)- k1(2))**2)/dble(NN-1) 

           if (i.gt.1) then
              t1=t1+temp
           endif
           k_len(i+(j-1)*NN)= t1
        enddo
        kpath_stop(j+1)= t1
     enddo

     ! sweep k
     ekslab=0.0d0
     do i=1+cpuid,knv3,num_cpu
        k= kpoint(i, :)
        chamk=0.0d0 
        call ham_slab(k,Chamk)
        eigenvalue=0.0d0

        ! diagonal Chamk
        call eigensystem_c('V', 'U', Num_wann*Nslab, CHamk, eigenvalue)
       
        ekslab(:,i)=eigenvalue

        do j=1, Nslab* Num_wann
           do l=1, Num_wann
              surf_weight(j, i)= surf_weight(j, i) &
                 + abs(CHamk(l, j))**2 & ! first slab
                 + abs(CHamk(Num_wann*Nslab- l, j))**2 !& ! last slab
                !+ abs(CHamk(Num_wann+ l, j))**2 & ! the second slab
                !+ abs(CHamk(Num_wann*(Nslab-1)- l, j))**2 ! last second slab
           enddo ! l
           surf_weight(j, i)= sqrt(surf_weight(j, i))
        enddo ! j 
        if (cpuid==0) print *,'SlabEk,k', i, knv3
     enddo ! i
     call mpi_allreduce(ekslab,ekslab_mpi,size(ekslab),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
     call mpi_allreduce(surf_weight, surf_weight_mpi,size(surf_weight),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
 

     ekslab=ekslab_mpi
     surf_weight= surf_weight_mpi/ maxval(surf_weight_mpi)
        
     if(cpuid==0)then
        open(unit=100, file='slabek.dat')
        do j=1, Num_wann*Nslab
           do i=1, knv3
             !write(100,'(3f15.7, i8)')k_len(i), ekslab(j,i), &
             !   (surf_weight(j, i))
              write(100,'(2f15.7, i8)')k_len(i), ekslab(j,i), &
                 int(255-surf_weight(j, i)*255d0)
           enddo
           write(100 , *)''
        enddo
        close(100)
        write(*,*) 'calculate energy band  done'
     endif

     emin= minval(ekslab)-0.5d0
     emax= maxval(ekslab)+0.5d0
     !> write script for gnuplot
     if (cpuid==0) then
        open(unit=101, file='slabek.gnu')
        write(101, '(a)') 'set terminal  postscript enhanced color'
        write(101,'(2a)') 'set palette defined ( 0  "green", ', &
           '5 "yellow", 10 "red" )'
        write(101, '(a)')"set output 'slabek.eps'"
        write(101, '(a)')'set style data linespoints'
        write(101, '(a)')'unset ztics'
        write(101, '(a)')'unset key'
        write(101, '(a)')'set pointsize 0.8'
        write(101, '(a)')'set view 0,0'
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
        write(101, '(a)')'rgb(r,g,b) = int(r)*65536 + int(g)*256 + int(b)'
        write(101, '(2a)')"plot 'slabek.dat' u 1:2:(rgb(255,$3, 3)) ",  &
            "w lp lw 2 pt 7  ps 1 lc rgb variable"

        !write(101, '(2a)')"splot 'slabek.dat' u 1:2:3 ",  &
        !   "w lp lw 2 pt 13 palette"
     endif

     202 format('set xtics (',:20('"',A3,'" ',F8.5,','))
     203 format(A3,'" ',F8.5,')')
     204 format('set arrow from ',F8.5,',',F10.5,' to ',F8.5,',',F10.5, ' nohead')
     
  return
  end
