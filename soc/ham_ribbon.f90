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




! This subroutine is used to caculate Hamiltonian for 
! open system in all 3 directions.

! History  
!       11/15/2019 by Quansheng Wu
! Copyright (c) 2010 QuanSheng Wu. All rights reserved.
! subroutine ham_cube(k,Hamk_ribbon)
! 
!    use para
!    implicit none

!    integer :: i1,j1,i2,j2

!    ! Hamiltonian of cube system
!    complex(Dp),intent(out) :: Hamk_cube(Num_wann*nslab1*nslab2*nslab3, &
!       Num_wann*nslab1*nslab2*nslab3) 

!    Hamk_ribbon=0.0d0 

!    ! i1,j1 row index
!    do i1=1,nslab1
!    do j1=1,nslab2
!    do l1=1,nslab3

!    ! i2,j2 column index
!    do i2=1,nslab1
!    do j2=1,nslab2  
!    do l2=1,nslab3  
!      if (abs(i2-i1).le.ijmax.and.abs(j2-j1).le.ijmax)then
!        call get_ir(i2-i1, j2-j1, l2-l1, ir)
!        if (ir==0) cycle
!        Hamk_ribbon((i1-1)*nslab2*Num_wann+(j1-1)*Num_wann+1: &
!                    (i1-1)*nslab2*Num_wann+(j1-1)*Num_wann+Num_wann,&
!                    (i2-1)*nslab2*Num_wann+(j2-1)*Num_wann+1: &
!                    (i2-1)*nslab2*Num_wann+(j2-1)*Num_wann+Num_wann )&
!        = HmnR(1:Num_wann,1:Num_wann, ir)/ndegen(ir)
!      endif 
!    enddo
!    enddo
!    enddo
!    enddo
!    enddo
!    enddo

! return
! end ham_cube

  subroutine get_ir(ia, ib, ic, ir)
     use para, only : dp, irvec, Nrpts
     implicit none

     integer :: ia, ib, ic, ir
     integer :: i

     ir = 0
     do i=1, Nrpts
        if (irvec(1, ir)==ia.and.irvec(2, ir)==ib.and.irvec(3, ir)==ic) then
           ir=i
        endif
     enddo

     return

  end subroutine get_ir
