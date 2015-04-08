! This subroutine is used to caculate Hamiltonian for 
! bulk system . 

  subroutine hamk_bulk_nsoc(k,Hamk_bulk)
  
     use para
     implicit none

     ! loop index  
     integer :: i1,i2,ia,ib,ic,iR

     real(dp) :: kdotr

     ! wave vector in 2d
     real(Dp) :: k(3)      

     ! coordinates of R vector 
     real(Dp) :: R(3)      

     ! Hamiltonian of bulk system
     complex(Dp),intent(out) ::Hamk_bulk(num_wann_soc, num_wann_soc) 

     Hamk_bulk=0d0
     do iR=1, Nrpts_nsoc
        ia=irvec_nsoc(1,iR)  
        ib=irvec_nsoc(2,iR) 
        ic=irvec_nsoc(3,iR)

        R(1)=dble(ia)
        R(2)=dble(ib)
        R(3)=dble(ic)
        
        kdotr=k(1)*R(1) + k(2)*R(2) + k(3)*R(3)
      
        Hamk_bulk(:,:)=Hamk_bulk(:,:)&
        +HmnR_nsoc(:,:,iR)*Exp(2d0*pi*zi*kdotr)/ndegen_nsoc(iR)
     enddo


     !> check hermitcity
     do i1=1, Num_wann_soc
     do i2=1, Num_wann_soc
        if(abs(Hamk_bulk(i1,i2)-conjg(Hamk_bulk(i2,i1))).ge.1e-6)then
          write(*,*)'there are something wrong with Hamk_bulk'
          stop
        endif 
     enddo
     enddo

     return
  end subroutine Hamk_bulk_nsoc



  subroutine hamk_bulk_soc(k,Hamk_bulk)
  
     use para
     implicit none

     ! loop index  
     integer :: i1,i2,ia,ib,ic,iR

     real(dp) :: kdotr

     ! wave vector in 2d
     real(Dp) :: k(3)      

     ! coordinates of R vector 
     real(Dp) :: R(3)      

     ! Hamiltonian of bulk system
     complex(Dp),intent(out) ::Hamk_bulk(num_wann_soc, num_wann_soc) 

     Hamk_bulk=0d0
     do iR=1, Nrpts_soc
        ia=irvec_soc(1,iR)  
        ib=irvec_soc(2,iR) 
        ic=irvec_soc(3,iR)

        R(1)=dble(ia)
        R(2)=dble(ib)
        R(3)=dble(ic)
        
        kdotr=k(1)*R(1) + k(2)*R(2) + k(3)*R(3)
      
        Hamk_bulk(:,:)=Hamk_bulk(:,:)&
        +HmnR_soc(:,:,iR)*Exp(2d0*pi*zi*kdotr)/ndegen_soc(iR)
     enddo


     !> check hermitcity
     do i1=1, Num_wann_soc
     do i2=1, Num_wann_soc
        if(abs(Hamk_bulk(i1,i2)-conjg(Hamk_bulk(i2,i1))).ge.1e-6)then
          write(*,*)'there are something wrong with Hamk_bulk'
          stop
        endif 
     enddo
     enddo

     return
  end subroutine Hamk_bulk_soc
