! This subroutine is used to caculate energy dispersion for 
! slab Bi2Se3

  subroutine psik()
    
     use para,only : Dp,Num_wann,Nslab
     implicit none 

! loop index
     integer     :: i,j     
     integer     :: info

! wave vector 
     real(Dp)    :: k(2)
 
! parameters for zheev
     real(Dp)    :: rwork(5*Num_wann*nslab)
     complex(Dp) :: work (5*Num_wann*nslab)
      
! eigenvalue 
     real(Dp)    :: eigenvalue (nslab*Num_wann)
     
   
! energy dispersion
     real(Dp)    :: ekslab(Nslab*Num_wann)
     real(Dp)    :: psi2  (Nslab)
     complex(Dp) :: psi   (Nslab*Num_wann)

! hamiltonian slab
     complex(Dp),allocatable :: CHamk(:,:)
     complex(Dp),allocatable :: eigenvector(:,:)


     allocate(CHamk(nslab*Num_wann,nslab*Num_wann))
     allocate(eigenvector(nslab*Num_wann,nslab*Num_wann))

! sweep k
     ekslab=0.0d0

! Ka direction
     k=0.0d0
     chamk=0.0d0 

! calculate Hamiltonian
     call ham_slab(k,Chamk)
     eigenvalue=0.0d0
     eigenvector=Chamk

! diagonal Chamk
     print *, 'here'
     call eigensystem_c('V', 'U', Num_wann*nslab, eigenvector, eigenvalue)
    
     ekslab=eigenvalue
     
     info=18*Nslab+10

     psi(:)=eigenvector(:, info)
     write(stdout,*) 'eigenvalue',info,ekslab(info)

     j=0
     psi2=0.0d0
     do i=1,Nslab
        do j=1,Num_wann/2
           psi2(i)=psi2(i)+abs(psi((i-1)*Num_wann+j))**2
        enddo
     enddo

     open(unit=100, file='squarepsi.dat')
     do i=1,Nslab
        write(100,'(2f16.9)')real(i),psi2(i)
     enddo
     close(100)
     
     write(stdout,*) 'calculate psi  done'
     
  end
