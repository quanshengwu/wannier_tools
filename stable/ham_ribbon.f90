! This subroutine is used to caculate Hamiltonian for 
! ribbon system . 

! History  
!        4/15/2010 by Quansheng Wu
! Copyright (c) 2010 QuanSheng Wu. All rights reserved.
  subroutine ham_ribbon(k,Hamk_ribbon)
  
     use para
     implicit none

! loop index  
     integer :: i1,j1,i2,j2

! wave vector in 2d
     real(Dp) :: k      

! Hamiltonian of slab system
     complex(Dp),intent(out) :: Hamk_ribbon(Num_wann*nslab1*nslab2, &
        Num_wann*nslab1*nslab2) 

! the factor 2 is induced by spin
     complex(Dp), allocatable :: Hij(:, :, :, :)

     allocate(Hij(-ijmax:ijmax,-ijmax:ijmax,Num_wann,Num_wann))
     call ham_qlayer2qlayerribbon(k,Hij)

     Hamk_ribbon=0.0d0 

! i1,j1 row index
     do i1=1,nslab1
     do j1=1,nslab2
! i2,j2 column index
     do i2=1,nslab1
     do j2=1,nslab2  
       if (abs(i2-i1).le.ijmax.and.abs(j2-j1).le.ijmax)then
         Hamk_ribbon((i1-1)*nslab2*Num_wann+(j1-1)*Num_wann+1: &
                     (i1-1)*nslab2*Num_wann+(j1-1)*Num_wann+Num_wann,&
                     (i2-1)*nslab2*Num_wann+(j2-1)*Num_wann+1: &
                     (i2-1)*nslab2*Num_wann+(j2-1)*Num_wann+Num_wann )&
         = Hij(i2-i1,j2-j1,1:Num_wann,1:Num_wann)

       endif 
     enddo
     enddo
     enddo
     enddo

  return
  end
