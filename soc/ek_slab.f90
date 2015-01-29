! This subroutine is used to caculate energy dispersion with wannier functions

  subroutine ek_slab
    
     use mpi
     use para
     implicit none 

! loop index
     integer :: i     
     integer :: j
     integer :: NN, nlines, knv3

     integer :: lwork

! wave vector 
     real(Dp) :: k(2), kstart(2), kend(2)
     real(dp) :: kp(16,2), ke(16,2)

     real(dp) :: cell_volume
     real(dp) :: t1, temp
 
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

! hamiltonian slab
     complex(Dp),allocatable ::CHamk(:,:)

     lwork= 16*Nslab*Num_wann
     ierr = 0

     kp(1,:)=(/0.0d0, 0.5d0/)  
     ke(1,:)=(/0.0d0, 0.0d0/)  
     kp(2,:)=(/0.0d0, 0.0d0/)  
     ke(2,:)=(/0.5d0, 0.00d0/)  ! K
     kp(3,:)=(/0.5d0, 0.00d0/)  ! K
     ke(3,:)=(/0.5d0, 0.5d0/)  ! K
     kp(4,:)=(/0.5d0, 0.5d0/)  ! K
     ke(4,:)=(/0.0d0, 0.0d0/)  ! K
     kp(5,:)=(/0.0d0, 0.0d0/)  ! K
     ke(5,:)=(/0.5d0, 0.5d0/)  ! K
     kp(6,:)=(/0.0d0, 0.0d0/)  ! K
     ke(6,:)=(/1.0d0, 1.0d0/)  ! K
    
     nlines= 4
     NN=Nk
     knv3=NN*nlines
     allocate( kpoint(knv3, 3))
     allocate( k_len (knv3))
     allocate(ekslab(Nslab*Num_wann,knv3))
     allocate(ekslab_mpi(Nslab*Num_wann,knv3))
     allocate(CHamk(nslab*Num_wann,nslab*Num_wann))
     allocate(work(lwork))
     allocate(rwork(lwork))
 
     kpoint= 0d0

     t1=0d0
     k_len=0d0
     do j=1, nlines 
        do i=1, NN
           k = kp(j,:)
           kstart=k

           k = ke(j,:)
           kend=k

           kpoint(i+(j-1)*NN,:)= kstart+ (kend-kstart)*dble(i-1)/dble(NN-1)
           
           temp= dsqrt((kstart(1)- kend(1))**2 &
                 +(kstart(2)- kend(2))**2)/dble(NN-1) 

           if (i.gt.1) then
              t1=t1+temp
           endif
           k_len(i+(j-1)*NN)= t1
        enddo
     enddo

     ! sweep k
     ekslab=0.0d0
     do i=1+cpuid,knv3,num_cpu
        k= kpoint(i, :)
        chamk=0.0d0 
        call ham_slab(k,Chamk)
        eigenvalue=0.0d0

        ! diagonal Chamk
        call eigensystem_c('N', 'U', Num_wann*nslab, CHamk, eigenvalue)
       
        ekslab(:,i)=eigenvalue
        print *,'cpuid,k',cpuid,i,ekslab(1,i)
     enddo
     call mpi_allreduce(ekslab,ekslab_mpi,size(ekslab),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
  
     ekslab=ekslab_mpi
        
     if(cpuid==0)then
        open(unit=100, file='slabek.dat')
        do j=1, Num_wann*Nslab
           do i=1, knv3
              write(100,'(10000f15.7)')k_len(i), ekslab(j,i)
           enddo
           write(100 , *)''
        enddo
        close(100)
        write(*,*) 'calculate energy band  done'
     endif
     
  return
  end
