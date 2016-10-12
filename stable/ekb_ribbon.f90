! This subroutine is used to caculate energy dispersion for 
! slab Bi2Se3
! Copyright (c) 2010 QuanSheng Wu. All rights reserved.

  subroutine ekb_ribbon
    
     use para
     implicit none 

     integer :: mdim
     integer :: ndim1

     ! loop index
     integer :: i     

     integer :: lwork

     ! the indices of the smallest and largest
     ! eigenvalues to be returned
     integer :: il,iu


     ! wave vector 
     real(Dp) :: k

     real(Dp) :: kmax=0.5d0
 
     !the lower and upper bounds of the interval
     !to be searched for eigenvalues
     real(Dp) :: vl,vu

     !The absolute error tolerance for the eigenvalues
     real(Dp) :: abstol

     ! parameters for zheev
     integer :: info,ierr,ifail

     ! number of magnetic fields
     integer :: Nb
      
     ! magnetic field
     real(dp) :: Bmag
     real(dp), allocatable :: magnetic(:)
   
     ! energy dispersion
     real(Dp),allocatable :: ekbribbon(:,:)
     real(Dp),allocatable ::  rwork(:)
     integer,allocatable ::  iwork(:)
     complex(Dp),allocatable :: work(:)

     ! eigenvalue 
     real(Dp),allocatable :: eigenvalue(:)

     ! hamiltonian slab
     complex(Dp),allocatable ::z(:,:)
     complex(Dp),allocatable ::CHamk(:,:)


     Nb= 400
     Ndim1=Num_wann*nslab1*nslab2
     lwork=64*Ndim1

     allocate(ekbribbon(Ndim1,Nb ))
     allocate(z(Ndim1,Ndim1))
     allocate(CHamk(Ndim1,Ndim1))
     allocate(rwork(7*Ndim1))
     allocate(work(lwork))
     allocate(iwork(5*Ndim1))
     allocate(eigenvalue(Ndim1))
     allocate(magnetic(Nb))
 

     ierr = 0
     
     ! sweep k
     ekbribbon=0.0d0

     kmax=0.5d0

     abstol=1e-10
     vl=-0.6/27.2114d0
     vu=0.5/27.2114d0
     il=17*Nslab1*Nslab2
     iu=20*Nslab1*Nslab2
     mdim=iu-il+1

     do i=1, Nb
        magnetic(i)= (i-1)*1.0d0/real(Nb)
     enddo

     do i=1,Nb
        Bmag= magnetic(i)
        k=0d0
        chamk=0.0d0 
        call ham_ribbon_b(bmag, k,Chamk)
        eigenvalue=0.0d0

        ! diagonal Chamk
        call eigensystem_c('N', 'U', Ndim1, CHamk, eigenvalue) 

        ! only eigenvalues are computed
        ! the eigenvalues with indices il through iu will be found

       !call zheevx('N','I','U',Ndim1,Chamk,Ndim1,vl,vu,il,iu,abstol,&
       !mdim,eigenvalue,z,Ndim1,work,lwork,rwork,iwork,ifail,info)


        ekbribbon(:,i)=eigenvalue

        write(stdout,'(a2,i4,f12.5,f10.2,a2)')'B, Nb', i, Nb, ekbribbon(1,i)     
     enddo
     
    !open(unit=100, file='ribbonek.dat',status='unknown')
    !do i=1,Nk1
    !   k=-kmax*real(Nk1-i)/(Nk1-1)
    !   write(100,'(60000f15.7)')k,ekbribbon(:,Nk1-i+1)
    !enddo 
  
    !do i=1,Nk1
    !   k=kmax*real(i-1)/(Nk1-1)
    !   write(100,'(60000f15.7)')k,ekbribbon(:,i)
    !enddo
    !close(100)

    open(unit=100, file='ribbonekb.dat',status='unknown')
    do i=1, Nb
       write(100, '(50000f15.7)')magnetic(i), ekbribbon(:, i)
    enddo
    !write(stdout,*) 'calculate energy band  done'

  return
  end
