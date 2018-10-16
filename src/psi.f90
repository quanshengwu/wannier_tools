! This subroutine is used to caculate energy dispersion for 
! slab Bi2Se3

  subroutine psik()

     use para,only : Dp,Num_wann,Nslab, stdout, cpuid 
     
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
     call eigensystem_c('V', 'U', Num_wann*nslab, eigenvector, eigenvalue)
    
     ekslab=eigenvalue
    
     !> with given band index
     info=18*Nslab+2 

     psi(:)=eigenvector(:, info)
     if (cpuid.eq.0) write(stdout,*) 'eigenvalue',info,ekslab(info)

     j=0
     psi2=0.0d0
     do i=1,Nslab
        do j=1,Num_wann
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

! This subroutine is used to caculate energy dispersion for 
! bulk system

  subroutine psik_bulk()
    
     use para,only : Dp,Num_wann, stdout, Numoccupied, cpuid
     implicit none 

! loop index
     integer     :: i,ik, ib
     integer     :: info

! wave vector 
     real(Dp)    :: k(3)
     real(Dp)    :: kpoints(3, 8)
      
! eigenvalue 
     real(Dp)    :: eigenvalue (Num_wann)
     
   
! energy dispersion
     complex(Dp) :: psi   (Num_wann)

! hamiltonian slab
     complex(Dp),allocatable :: CHamk(:,:)
     complex(Dp),allocatable :: eigenvector(:,:)


     allocate(CHamk(Num_wann,Num_wann))
     allocate(eigenvector(Num_wann,Num_wann))


     kpoints(:, 1)= (/0.0d0, 0.0d0, 0.0d0/)
     kpoints(:, 2)= (/0.5d0, 0.0d0, 0.0d0/)
     kpoints(:, 3)= (/0.0d0, 0.5d0, 0.0d0/)
     kpoints(:, 4)= (/0.0d0, 0.0d0, 0.5d0/)
     kpoints(:, 5)= (/0.5d0, 0.5d0, 0.5d0/)
     kpoints(:, 6)= (/0.5d0, 0.5d0, 0.0d0/)
     kpoints(:, 7)= (/0.0d0, 0.5d0, 0.5d0/)
     kpoints(:, 8)= (/0.5d0, 0.0d0, 0.5d0/)

     ib= Numoccupied
     if (cpuid==0)open(unit=100, file='wavefunction.dat')
     do ik=1, 8

        k=kpoints(:, ik)
        ! calculate Hamiltonian
        call ham_bulk(k,Chamk)
        eigenvalue=0.0d0
        eigenvector=Chamk
        
        ! diagonal Chamk
        call eigensystem_c('V', 'U', Num_wann, eigenvector, eigenvalue)
       
        psi(:)=eigenvector(:, ib)
        if (cpuid==0)write(stdout,*) 'eigenvalue',info,eigenvalue(ib)
   
        if (cpuid==0)write(100, '(a,3f8.4)')'K point ', k
        do i=1, Num_wann
           if (cpuid==0)write(100,'(i5, 30f16.9)')i, eigenvector(i, ib-1), eigenvector(i, ib)
        enddo
        if (cpuid==0)write(100,*)' '
    
     enddo ! ik
        
     if (cpuid==0)close(100)
  end subroutine psik_bulk
