!+---------+---------+---------+---------+---------+---------+--------+!
! this subroutine is used to calculate surface state using             !
! green's function method  ---  J.Phys.F.Met.Phys.15(1985)851-858      !
! 
! History:
!         by Quan Sheng Wu on 4/20/2010                                !
!            mpi version      4/21/2010
!            change Kb to K=(Ka+Kb)/3 direction 4/22/2010
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

     real(dp) :: t1, temp
     real(dp) :: k(2), w

     real(dp) :: kp(6, 2)
     real(dp) :: ke(6, 2)

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

     kp(1,:)=(/-0.5d0, 0.000d0/)  ! K
     ke(1,:)=(/0.0d0,  0.000d0/)  ! Gamma
     kp(2,:)=(/0.0d0,  0.000d0/)  ! X
     ke(2,:)=(/0.5d0,  0.000d0/)  ! M
     kp(3,:)=(/-0.0d0, 0.500d0/)  ! K
     ke(3,:)=(/0.0d0,  0.000d0/)  ! Gamma
     kp(4,:)=(/0.0d0,  0.000d0/)  ! X
     ke(4,:)=(/0.0d0,  0.500d0/)  ! M
     kp(5,:)=(/0.0d0, 0.0d0/)  ! Gamma
     kp(6,:)=(/5.0d0, 0.0d0/)  ! Z
    
     nlines=2
     NN=Nk
     knv2=NN*nlines
     allocate( kpoint(knv2, 2))
     allocate( k_len (knv2))
     kpoint= 0d0

     t1=0d0
     k_len=0d0
     do j=1, nlines 
        do i=1, NN
           kstart= kp(j,:)
           kend  = ke(j,:)
           kpoint(i+(j-1)*NN,:)= kstart+ (kend-kstart)*(i-1)/dble(NN-1)
           
           temp= dsqrt((ke(j,1)- kp(j,1))**2 &
                 +(ke(j,2)- kp(j,2))**2)/dble(NN-1) 

           if (i.gt.1) then
              t1=t1+temp
           endif
           k_len(i+(j-1)*NN)= t1
        enddo
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

     eta=(omegamax- omegamin)/dble(omeganum)*1.0d0

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
        if (cpuid==0) print *, ikp, 'in', knv2
        k= kpoint(ikp,:)

        !> get the hopping matrix between two principle layers
        call ham_qlayer2qlayer(k,H00,H01) 

        do j = 1, omeganum
           w=omega(j)

           !> calculate surface green function
           ! there are two method to calculate surface green's function 
           ! the method in 1985 is better, you can find the ref in the
           ! subroutine
           call surfgreen_1985(omega,GLL,GRR,H00,H01,ones)
           ! call surfgreen_1984(omega,GLL,GRR,H00,H01,ones)

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
     dos_l=dos_l_mpi
     dos_r=dos_r_mpi

     if (cpuid.eq.0)then
        open (unit=12, file='dos.dat_l')
        open (unit=13, file='dos.dat_r')
        do ikp=1, knv2
           do j=1, omeganum 
              write(12, '(3f16.8)')k_len(ikp), omega(j)*27.2114d0, dos_l(ikp, j)
              write(13, '(3f16.8)')k_len(ikp), omega(j)*27.2114d0, dos_r(ikp, j)
           enddo
        enddo
        close(12)
        close(13)
        write(*,*)'ndim',ndim
        write(*,*) 'knv2,omeganum,eta',knv2, omeganum, eta
        write(*,*)'calculate density of state successfully'    
     endif

  return   
  end subroutine
