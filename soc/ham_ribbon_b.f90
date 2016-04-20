! This subroutine is used to caculate Hamiltonian for 
! ribbon system .  
! with magnetic field along x direction

! History  
!        4/15/2010 by Quansheng Wu
! Copyright (c) 2010 QuanSheng Wu. All rights reserved.
  subroutine ham_ribbon_b(bmag,k,Hamk_ribbon)
  
     use para
     implicit none

! loop index  
     integer :: i1,j1,i2,j2

! wave vector in 2d
     real(Dp), intent(in) :: k      
     real(dp), intent(in) :: bmag

! Hamiltonian of slab system
     complex(Dp),intent(out) ::Hamk_ribbon(Num_wann*nslab1*nslab2,Num_wann*nslab1*nslab2) 

     complex(Dp) :: Hij(Num_wann,Num_wann)

     integer :: ir

     integer :: ia, ib, ic
     integer :: new_ia, new_ib, new_ic

     real(dp) :: kdotr
     real(dp) :: phase
     complex(dp) :: ratio
     complex(dp) :: fac

     real(dp) :: xyz1(3)
     real(dp) :: xyz2(3)



     Hamk_ribbon=0.0d0 
 
     ! i1,j1 row index
     do i1=1,nslab1
     do j1=1,nslab2
 
     ! i2,j2 column index
     do i2=1,nslab1
     do j2=1,nslab2  

        Hij= zzero

        ! summation over R 
        do ir=1, nrpts
           ia= irvec(1, ir)
           ib= irvec(2, ir)
           ic= irvec(3, ir)
           new_ia= ia
           new_ib= ib
           new_ic= ic

           xyz1(1)= 0
           xyz1(2)= dble(i1)
           xyz1(3)= dble(j1)
          
           xyz2(1)= dble(ia)
           xyz2(2)= dble(i2)
           xyz2(3)= dble(j2)
           
           ! B along x direction 
          !phase= Bmag*(xyz2(3)+xyz1(3))*(xyz2(2)-xyz1(2))/2d0
          !fac= cos(phase)+ zi*sin(phase)
          
           ! B along z direction 
           phase=-Bmag*(xyz2(3)+xyz1(3))*(xyz2(1)-xyz1(1))/2d0
           fac= cos(phase)+ zi*sin(phase)


           kdotr= k*real(new_ia)
           ratio= cos(2d0*pi*kdotr)+ zi*sin(2d0*pi*kdotr)
           if (new_ib.eq.i2-i1.and.new_ic.eq.j2-j1)then
              Hij= Hij+ Hmnr(:, :, ir)/real(ndegen(ir))*fac*ratio
           endif

        enddo


       if (abs(i2-i1).le.ijmax.and.abs(j2-j1).le.ijmax)then
         Hamk_ribbon((i1-1)*nslab2*Num_wann+(j1-1)*Num_wann+1: &
                     (i1-1)*nslab2*Num_wann+(j1-1)*Num_wann+Num_wann,&
                     (i2-1)*nslab2*Num_wann+(j2-1)*Num_wann+1: &
                     (i2-1)*nslab2*Num_wann+(j2-1)*Num_wann+Num_wann )&
         = Hij

       endif 


     enddo
     enddo
     enddo
     enddo

  return
  end
