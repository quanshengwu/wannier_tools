! This subroutine is used to caculate Hamiltonian for 
! slab system . 

! History  
!        4/18/2010 by Quansheng Wu
  subroutine ham_slab(k,Hamk_slab)
  
     use para
     implicit none

! loop index  
     integer :: i1,i2, i
     integer :: ir

! wave vector in 2d
     real(Dp), intent(inout) :: k(2)      

! Hamiltonian of slab system
     complex(Dp),intent(out) ::Hamk_slab(Num_wann*nslab,Num_wann*nslab) 

! the factor 2 is induced by spin
     complex(Dp) :: Hij(-ijmax:ijmax,Num_wann,Num_wann)

     call ham_qlayer2qlayer2(k,Hij)


     Hamk_slab=0.0d0 
! i1,j1 column index
     do i1=1,nslab
! i2,j2 row index
     do i2=1,nslab
       if (abs(i2-i1).le.ijmax)then
         Hamk_slab((i2-1)*Num_wann+1:(i2-1)*Num_wann+Num_wann,&
                   (i1-1)*Num_wann+1:(i1-1)*Num_wann+Num_wann )&
         = Hij(i2-i1,1:Num_wann,1:Num_wann)

       endif 
       if ((i1.eq.1.and.i2.eq.1).or.(i1.eq.nslab.and.i2.eq.nslab))then
          do i=1, Num_wann
             Hamk_slab((i2-1)*Num_wann+i,(i1-1)*Num_wann+i )&
            =Hamk_slab((i2-1)*Num_wann+i,(i1-1)*Num_wann+i )
          enddo
       endif


     enddo
     enddo

 
! check hermitcity

     do i1=1,nslab*Num_wann
     do i2=1,nslab*Num_wann
        if(abs(Hamk_slab(i1,i2)-conjg(Hamk_slab(i2,i1))).ge.1e-6)then
          write(stdout,*)'there are something wrong with Hamk_slab'
          stop
        endif 
     enddo
     enddo

  return
  end subroutine ham_slab

!> This subroutine is used to caculate Hamiltonian for 
!> slab system under the in-plane magnetic field
!> the magnetic phase are chosen like this
!> phi_ij= alpha*[By*(xj-xi)*(zi+zj)-Bx*(yj-yi)*(zi+zj)] 
!> x, z are in unit of Angstrom, Bx, By are in unit of Tesla
!> History :
!        9/17/2015 by Quansheng Wu @ETH Zurich
  subroutine ham_slab_parallel_B(k,Hamk_slab)
  
     use para
     implicit none

! loop index  
     integer :: i1,i2, i
     integer :: ir

     integer :: ir

     integer :: ia, ib, ic
     integer :: new_ia, new_ib, new_ic

     real(dp) :: kdotr
     real(dp) :: phase
     complex(dp) :: ratio
     complex(dp) :: fac

     real(dp) :: xyz1(3)
     real(dp) :: xyz2(3)

! wave vector in 2d
     real(Dp), intent(inout) :: k(2)      

! Hamiltonian of slab system
     complex(Dp),intent(out) ::Hamk_slab(Num_wann*nslab,Num_wann*nslab) 

! the factor 2 is induced by spin
     complex(Dp), allocatable :: Hij(:, :, :)

     allocate(Hij(-ijmax:ijmax,Num_wann,Num_wann))
     Hij= 0d0

     !> construct hamiltonian for slab system under 
     Hamk_slab=0.0d0 
     ! i1 column index
     do i1=1,nslab
     ! i2 row index
     do i2=1,nslab
       if (abs(i2-i1).le.ijmax)then
         Hamk_slab((i2-1)*Num_wann+1:(i2-1)*Num_wann+Num_wann,&
                   (i1-1)*Num_wann+1:(i1-1)*Num_wann+Num_wann )&
         = Hij(i2-i1,1:Num_wann,1:Num_wann)

       endif 
     enddo
     enddo

     ! check hermitcity
     do i1=1,nslab*Num_wann
     do i2=1,nslab*Num_wann
        if(abs(Hamk_slab(i1,i2)-conjg(Hamk_slab(i2,i1))).ge.1e-6)then
          write(stdout,*)'there are something wrong with Hamk_slab'
          stop
        endif 
     enddo
     enddo

  return
  end subroutine ham_slab_parallel_B
