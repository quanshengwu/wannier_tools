! this subroutine is used to calculate surface state using             !
! green's function method  ---  J.Phys.F.Met.Phys.15(1985)851-858      !
! Using eq.(11) and eq.(13) (16), (17)
! History:
!         by Quan Sheng Wu on Oct/17/2012                                !
!+---------+---------+---------+---------+---------+---------+--------+!
  subroutine surfgreen_1985(omega,GLL,GRR,GB,H00,H01,ones, eta_broadening)
     use para
     implicit none

     ! inout variables     
     ! the factor 2 is induced by spin
     ! energy hbar omega
     real(Dp),intent(in) :: omega  
     real(Dp),intent(in) :: eta_broadening

     ! H00 Hamiltonian between nearest neighbour-quintuple-layers
     complex(Dp),intent(in) :: H00(Ndim,Ndim)

     ! H01 Hamiltonian between next-nearest neighbour-quintuple-layers
     complex(Dp),intent(in) :: H01(Ndim,Ndim)

     ! temp hamiltonian

     complex(Dp),intent(in)   :: ones(Ndim,Ndim)

     ! surface green function
     complex(Dp),intent(inout)  :: GLL(Ndim,Ndim)
     complex(Dp),intent(inout)  :: GRR(Ndim,Ndim)

     !> bulk green's function
     complex(Dp),intent(inout)  :: GB(Ndim,Ndim)

     ! >> local variables
     ! iteration number
     integer :: iter

     ! maximun iteration 
     integer ,parameter:: itermax=100

     ! accuracy control
     real(Dp) :: accuracy=1e-16

     ! a real type temp variable
     real(Dp) :: real_temp

     ! omegac=omega(i)+I * Fermi_broadening
     complex(Dp) :: omegac 


     ! some variables in Eq.(11)
     complex(Dp), allocatable :: alphai(:, :) 
     complex(Dp), allocatable :: betai(:, :) 
     complex(Dp), allocatable :: epsiloni(:, :) 
     complex(Dp), allocatable :: epsilons(:, :) 
     complex(Dp), allocatable :: epsilons_t(:, :) 

     complex(Dp), allocatable :: mat1 (:, :) 
     complex(Dp), allocatable :: mat2 (:, :) 

     ! g0= inv(w-e_i)
     complex(Dp), allocatable :: g0 (:, :) 

     ! allocate some variables
     allocate(alphai(Ndim, Ndim)) 
     allocate(betai (Ndim, Ndim)) 
     allocate(epsiloni (Ndim, Ndim)) 
     allocate(epsilons (Ndim, Ndim)) 
     allocate(epsilons_t(Ndim, Ndim)) 
     allocate(mat1(Ndim, Ndim)) 
     allocate(mat2(Ndim, Ndim)) 
     allocate(g0(Ndim, Ndim)) 

     epsiloni= H00
     epsilons= H00
     epsilons_t= H00
     alphai  = H01
     betai   = conjg(transpose(H01))
    !print *, sqrt(sum(abs(H00)**2)), 'H00'

     ! w+i*0^+
     omegac= dcmplx(omega, eta_broadening)
    !print *, omegac

     ! begin iteration
     do iter=1, itermax

        g0= omegac*ones- epsiloni
        call inv(Ndim, g0)

        ! a_i-1*(w-e_i-1)^-1
        call mat_mul(Ndim, alphai, g0, mat1 )
        
        ! b_i-1*(w-e_i-1)^-1
        call mat_mul(Ndim, betai, g0, mat2 )

        ! a_i-1*(w-e_i-1)^-1*b_i-1
        call mat_mul(Ndim, mat1, betai, g0)
        epsiloni= epsiloni+ g0
       !print *, sqrt(sum(abs(epsiloni)**2)), 'ei'
        ! es_i= es_i-1 + a_i-1*(w-e_i-1)^-1*b_i-1
        epsilons= epsilons+ g0
       !print *, sqrt(sum(abs(epsilons)**2)), 'es'
       !pause

        ! b_i-1*(w-e_i-1)^-1*a_i-1
        call mat_mul(Ndim, mat2, alphai, g0)
        epsiloni= epsiloni+ g0
        ! es_i= es_i-1 + a_i-1*(w-e_i-1)^-1*b_i-1
        epsilons_t= epsilons_t+ g0

        ! a_i= a_i-1*(w-e_i-1)^-1*a_i-1 
        call mat_mul(Ndim, mat1, alphai, g0)
        alphai= g0
        ! b_i= b_i-1*(w-e_i-1)^-1*b_i-1 
        call mat_mul(Ndim, mat2, betai, g0)
        betai= g0

       !real_temp=maxval(abs(alphai))   
        real_temp=sum(abs(alphai))   
       !if (cpuid.eq.0) print *, iter, real_temp
        if (real_temp.le.accuracy) exit

     enddo ! end of iteration

     ! calculate surface green's function
     GLL= omegac*ones- epsilons
     call inv(Ndim, GLL)

     GRR= omegac*ones- epsilons_t
     call inv(Ndim, GRR)

     GB = omegac*ones- epsiloni
     call inv(Ndim, GB)

     return
  end subroutine surfgreen_1985


!+---------+---------+---------+---------+---------+---------+--------+!
! this subroutine is used to calculate surface state using             !
! green's function method  ---  J.Phys.F.Met.Phys.14(1984)1205-1215    !
! Quick iterative scheme for the calculation of transfer matrices:
! application to Mo (100) 
! History:
!         by Quan Sheng Wu on 4/20/2010                                !
!            mpi version      4/21/2010
!            little change to cpu version in Zurich Swiss Jan 25 2015
!+---------+---------+---------+---------+---------+---------+--------+!
  subroutine surfgreen_1984(omega,GLL,GRR,H00,H01,ones, eta_broadening)

     use wmpi
     use para
     implicit none
     

     ! general loop index
     integer :: i,j 

     ! iteration loop index
     integer :: it

     ! iteration number
     integer :: iter

     ! maximun iteration 
     integer ,parameter :: itermax=100

     ! accuracy control
     real(Dp),parameter :: accuracy=1e-6

     ! a real type temp variable
     real(Dp) :: real_temp

     ! frequency 
     real(Dp),intent(in) :: omega

     ! energy energy=omega(i)+I * Fermi_broadening
     complex(Dp) :: energy

     ! surface green fuction  an output variable
     complex(Dp),intent(out)  :: GLL(Ndim,Ndim)
     complex(Dp),intent(out)  :: GRR(Ndim,Ndim)

     ! H00 Hamiltonian between nearest neighbour-quintuple-layers
     ! the factor 2 is induced by spin
     complex(Dp),intent(in)  :: H00(Ndim,Ndim)

     ! H01 Hamiltonian between next-nearest neighbour-quintuple-layers
     complex(Dp),intent(in)  :: H01(Ndim,Ndim)

     ! unit matrix ones
     complex(Dp),intent(in)  :: ones(Ndim,Ndim)

     !> infinite small value broadening
     real(dp), intent(in) :: eta_broadening

     complex(Dp),allocatable  :: H01dag(:,:)

     !> surface hamiltonian
     complex(Dp),allocatable  :: Hs(:,:)

     ! temp hamiltonian
     complex(Dp),allocatable  :: t0(:,:)
     complex(Dp),allocatable  :: Tmatrix(:,:)
     complex(Dp),allocatable  :: Tmatrixt(:,:)
     complex(Dp),allocatable  :: t0tilde(:,:)
     complex(Dp),allocatable  :: tnew(:,:)
     complex(Dp),allocatable  :: tnewtilde(:,:)
     complex(Dp) ,allocatable :: told(:,:)
     complex(Dp),allocatable  :: toldtilde(:,:)
     complex(Dp),allocatable  :: temp(:,:)
     complex(Dp),allocatable  :: Tmat_temp(:,:)
     complex(Dp),allocatable  :: Tmat_tempt(:,:)

     real(Dp),allocatable     :: abs_told(:,:)

     ! allocate variables
     allocate(Hs(ndim,ndim))
     allocate(t0(ndim,ndim))
     allocate(Tmatrix(ndim,ndim))
     allocate(Tmatrixt(ndim,ndim))
     allocate(t0tilde(ndim,ndim))
     allocate(tnew(ndim,ndim))
     allocate(tnewtilde(ndim,ndim))
     allocate(told(ndim,ndim))
     allocate(toldtilde(ndim,ndim))
     allocate(temp(ndim,ndim))
     allocate(Tmat_temp(ndim,ndim))
     allocate(Tmat_tempt(ndim,ndim))
     allocate(abs_told(ndim,ndim))
    
     allocate(H01dag(ndim,ndim))
 
     Hs=0.0d0
     t0=0.0d0
     t0tilde=0.0d0
     told=0.0d0
     toldtilde=0.0d0
     tnew=0.0d0
     tnewtilde=0.0d0
     Tmatrix=0.0d0
     Tmatrixt=0.0d0
     temp=0.0d0
     abs_told=0.0d0 
     Tmat_temp=0.0d0
	  GLL=0d0
	  GRR=0d0
     Tmat_tempt=0.0d0
     H01dag=0.0d0 
   
     ! H01dag=H01^dag
     do i=1,ndim
        do j=1,ndim
           H01dag(i,j)=conjg(H01(j,i))
        enddo
     enddo

     energy=omega+zi*eta_broadening
     temp=energy*ones-H00


     call inv(ndim,temp)
     t0= matmul(temp,H01dag)  
     t0tilde= matmul(temp,H01)
     told=t0
     toldtilde=t0tilde

     ! begin iteration
     iter=0    
     Tmatrix=0.0d0 
     Tmat_temp=ones
     Tmatrixt=0.0d0 
     Tmat_tempt=ones
     ITER1 : do it=1,itermax
        iter=iter+1

        temp=ones-matmul(told,toldtilde)-matmul(toldtilde,told)
        call inv(ndim,temp)

        tnew=matmul(temp,told)
        tnew=matmul(tnew,told)
        tnewtilde=matmul(temp,toldtilde)
        tnewtilde=matmul(tnewtilde,toldtilde)

        Tmat_temp=matmul(Tmat_temp,toldtilde)
        Tmatrix=Tmatrix+matmul(Tmat_temp,tnew)
        
		  Tmat_tempt=matmul(Tmat_tempt,told)
        Tmatrixt=Tmatrixt+matmul(Tmat_tempt,tnewtilde)
     
        told=tnew
        toldtilde=tnewtilde
      
        do i=1,ndim 
           do j=1,ndim
              abs_told(i,j)=abs(told(i,j))+abs(toldtilde(i,j)) 
           enddo 
        enddo 
        real_temp=maxval(abs_told)   
        if (real_temp.le.accuracy) exit ITER1 
     
     end do ITER1

     !print *,'iter,acc',iter,real_temp

     !> set up surface hamiltonian
     !> usually Hs= H00
     !> but you can add static potential on the surface
     Hs= H00
     do i=1, num_wann
        Hs(i, i)=Hs(i, i)+ surf_onsite
     enddo
  
     Tmatrix=t0+Tmatrix
     temp=energy*ones-Hs-matmul(H01,Tmatrix)    
     call inv(ndim,temp)
     
     ! g_00=(epsilon-Hs -H01*T)^-1
     GLL(1:ndim,1:ndim)=temp

     !> usually Hs= H00
     !> but you can add static potential on the surface
     Hs= H00
     do i=ndim-num_wann+1, ndim
        Hs(i, i)=Hs(i, i)+ surf_onsite
     enddo
 
     Tmatrixt=t0tilde+Tmatrixt
     temp=energy*ones-Hs -matmul(conjg(transpose(H01)),Tmatrixt)  
     call inv(ndim,temp)

     ! g_00=(epsilon-Hs -H01*T)^-1
     GRR(1:ndim,1:ndim)=temp

     return   
  end subroutine surfgreen_1984
