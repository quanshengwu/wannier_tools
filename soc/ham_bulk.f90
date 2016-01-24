! This subroutine is used to caculate Hamiltonian for 
! bulk system . 

! History  
!        May/29/2011 by Quansheng Wu
  subroutine ham_bulk(k,Hamk_bulk)
  
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
     complex(Dp),intent(out) ::Hamk_bulk(Num_wann, Num_wann) 
     complex(dp), allocatable :: mat1(:, :)
     complex(dp), allocatable :: mat2(:, :)
     allocate(mat1(Num_wann, Num_wann))
     allocate(mat2(Num_wann, Num_wann))

     Hamk_bulk=0d0
     do iR=1, Nrpts
        ia=irvec(1,iR)  
        ib=irvec(2,iR) 
        ic=irvec(3,iR)

        R(1)=dble(ia)
        R(2)=dble(ib)
        R(3)=dble(ic)
        
        kdotr=k(1)*R(1) + k(2)*R(2) + k(3)*R(3)
      
        Hamk_bulk(:,:)=Hamk_bulk(:,:)&
        +HmnR(:,:,iR)*Exp(2d0*pi*zi*kdotr)/ndegen(iR)
     enddo

     call mat_mul(Num_wann, mirror_z, Hamk_bulk, mat1)
     call mat_mul(Num_wann, mat1, mirror_z, mat2)
     Hamk_bulk= (Hamk_bulk+ mat2)/2d0

     ! check hermitcity
     do i1=1, Num_wann
     do i2=1, Num_wann
        if(abs(Hamk_bulk(i1,i2)-conjg(Hamk_bulk(i2,i1))).ge.1e-6)then
          write(stdout,*)'there are something wrong with Hamk_bulk'
          stop
        endif 
     enddo
     enddo

  return
  end
