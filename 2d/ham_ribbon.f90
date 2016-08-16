! This subroutine is used to caculate Hamiltonian for 
! ribbon system . 

! History  
!        4/18/2010 by Quansheng Wu
! Copyright (c) 2010 QuanSheng Wu. All rights reserved.
  subroutine ham_ribbon(k,Hamk_ribbon)
  
     use para
     implicit none

     ! loop index  
     integer :: i1, i2, i, j
     integer :: ir

     ! wave vector in 2d
     real(Dp), intent(inout) :: k(2)      

     ! Hamiltonian of ribbon system
     complex(Dp),intent(out) ::Hamk_ribbon(Num_wann*nslab,Num_wann*nslab) 

     ! the factor 2 is induced by spin
     complex(Dp), allocatable :: Hij(:, :, :)

     allocate( Hij(-ijmax:ijmax,Num_wann,Num_wann))
     call ham_qlayer2qlayer2(k,Hij)

     Hamk_ribbon=0.0d0 
     ! i1 column index
     do i1=1, nslab
        ! i2 row index
        do i2=1, nslab
          if (abs(i2-i1).le.ijmax)then
            Hamk_ribbon((i2-1)*Num_wann+1:(i2-1)*Num_wann+Num_wann,&
                      (i1-1)*Num_wann+1:(i1-1)*Num_wann+Num_wann )&
            = Hij(i2-i1,1:Num_wann,1:Num_wann)
          endif 
        enddo ! i2
     enddo ! i1

     ! check hermitcity

     do i1=1,nslab*Num_wann
     do i2=1,nslab*Num_wann
        if(abs(Hamk_ribbon(i1,i2)-conjg(Hamk_ribbon(i2,i1))).ge.1e-6)then
          write(stdout,*)'there are something wrong with Hamk_ribbon'
          stop
        endif 
     enddo
     enddo

  return
  end subroutine ham_ribbon

