! This subroutine is used to caculate energy dispersion for 
! slab Bi2Se3

  subroutine ek_ribbon

     use mpi
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


     integer:: Nk1

     ! wave vector 
     real(Dp) :: k

     real(Dp) :: kmax=0.5d0
 
     real(Dp) :: time1,time2,time3

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

     ! eigenvalue 
     real(Dp),allocatable :: eigenvalue(:)

     ! hamiltonian slab
     complex(Dp),allocatable ::z(:,:)
     complex(Dp),allocatable ::CHamk(:,:)


     Ndim1=Num_wann*nslab1*nslab2
     lwork=64*Ndim1

     Nk1= Nk

     allocate(ekribbon(Ndim1,Nk1))
     allocate(ekribbon_mpi(Ndim1,Nk1))
     allocate(z(Ndim1,Ndim1))
     allocate(CHamk(Ndim1,Ndim1))
     allocate(rwork(7*Ndim1))
     allocate(work(lwork))
     allocate(iwork(5*Ndim1))
     allocate(eigenvalue(Ndim1))
 

     ierr = 0
     
     ! sweep k
     ekribbon=0.0d0
     ekribbon_mpi=0.0d0

     kmax=0.5d0

     abstol=1e-10
     vl=-0.6/27.2114d0
     vu=0.5/27.2114d0
     il=17*Nslab1*Nslab2
     iu=20*Nslab1*Nslab2
     mdim=iu-il+1

     print *,'number of bands calculating: ',mdim

     do i=1+cpuid, Nk1, num_cpu
        if (cpuid.eq.0) print *, "Ribbonek the i'th kpoint", i, Nk1
        k=kmax*real(i-1)/(Nk1-1)
        chamk=0.0d0 
        call ham_ribbon(k,Chamk)
        eigenvalue=0.0d0

        ! diagonal Chamk
        call eigensystem_c('N', 'U', Ndim1, CHamk, eigenvalue) 
       
        ! only eigenvalues are computed
        ! the eigenvalues with indices il through iu will be found
        !call zheevx('N','I','U',Ndim1,Chamk,Ndim1,vl,vu,il,iu,abstol,&
        !mdim,eigenvalue,z,Ndim1,work,lwork,rwork,iwork,ifail,info)

        ekribbon(:,i)=eigenvalue

        if (cpuid.eq.0) write(stdout,'(a2,i4,f12.5,f10.2,a2)')'k',i,ekribbon(1,i),time2-time1,' s'
     enddo
     call mpi_allreduce(ekribbon,ekribbon_mpi,size(ekribbon),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
  
     if (cpuid.eq.0) then
        open(unit=100, file='ribbonek.dat',status='unknown')
        do i=1,Nk1
           k=-kmax*real(Nk1-i)/(Nk1-1)
           write(100,'(60000f15.7)')k,ekribbon_mpi(:,Nk1-i+1)
        enddo 
     
        do i=1,Nk1
           k=kmax*real(i-1)/(Nk1-1)
           write(100,'(60000f15.7)')k,ekribbon_mpi(:,i)
        enddo
        close(100)
        write(stdout,*) 'calculate energy band  done'
     endif

     return
  end subroutine ek_ribbon
