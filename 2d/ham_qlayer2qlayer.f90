! this subroutine is used to caculate Hamiltonian between
! slabs  
! 4/23/2010 by QS Wu
! Copyright (c) 2010 QuanSheng Wu. All rights reserved.

  subroutine ham_qlayer2qlayer(k,H00new,H01new)

     use para

     implicit none

! loop index
     integer :: i,j,iR

! index used to sign irvec     
     integer :: ia,ib

! 

! new index used to sign irvec     
     real(dp) :: new_ia,new_ib
     integer :: inew_ib

! wave vector k times lattice vector R  
     real(Dp) :: kdotr  

! input wave vector k's cooridinates
     real(Dp),intent(in) :: k

! H00 Hamiltonian between nearest neighbour-quintuple-layers
! the factor 2 is induced by spin

     complex(dp) :: ratio

     complex(Dp), allocatable :: Hij(:, :, :)


! H00 Hamiltonian between nearest neighbour-quintuple-layers
! the factor 2 is induced by spin


!     complex(Dp),allocatable,intent(out) :: H00new(:,:)
     complex(Dp),intent(out) :: H00new(Ndim,Ndim)

! H01 Hamiltonian between next-nearest neighbour-quintuple-layers
!     complex(Dp),allocatable,intent(out) :: H01new(:,:)
     complex(Dp),intent(out) :: H01new(Ndim,Ndim)

     allocate(Hij(-ijmax:ijmax,Num_wann,Num_wann))

     Hij=0.0d0
     do iR=1,Nrpts
        ia=irvec(1,iR)
        ib=irvec(2,iR)

        call latticetransform(ia, ib, new_ia, new_ib)
       
        inew_ib= int(new_ib)
        if (abs(new_ib).le.ijmax)then
           kdotr=k*new_ia
           ratio=cos(2d0*pi*kdotr)+zi*sin(2d0*pi*kdotr)

           Hij(inew_ib, 1:Num_wann, 1:Num_wann )&
           =Hij(inew_ib, 1:Num_wann, 1:Num_wann )&
           +HmnR(:,:,iR)*ratio/ndegen(iR)
        endif

     enddo

     H00new=0.0d0
     H01new=0.0d0

! nslab's principle layer 
! H00new
     do i=1,Np
     do j=1,Np
        if (abs(i-j).le.(ijmax)) then
          H00new(Num_wann*(i-1)+1:Num_wann*i,Num_wann*(j-1)+1:Num_wann*j)&
                =Hij(j-i,:,:)
        endif
     enddo
     enddo

! H01new
     do i=1,Np
     do j=Np+1,Np*2
        if (j-i.le.ijmax) then
           H01new(Num_wann*(i-1)+1:Num_wann*i,&
               Num_wann*(j-1-Np)+1:Num_wann*(j-Np))=Hij(j-i,:,:)
        endif
     enddo
     enddo

     do i=1,Ndim
     do j=1,Ndim
        if(abs(H00new(i,j)-conjg(H00new(j,i))).ge.1e-4)then
          write(stdout,*)'there are something wrong with ham_qlayer2qlayer'
        stop
        endif

     enddo
     enddo

  return
  end subroutine ham_qlayer2qlayer

  ! for slab
  subroutine ham_qlayer2qlayer2(k,Hij)

     use para

     implicit none

! loop index
     integer :: iR

! index used to sign irvec     
     integer :: ia,ib

! 

! new index used to sign irvec     
     real(dp) :: new_ia,new_ib
     integer :: inew_ib

! wave vector k times lattice vector R  
     real(Dp) :: kdotr  

! input wave vector k's cooridinates
     real(Dp),intent(in) :: k

! H00 Hamiltonian between nearest neighbour-quintuple-layers
! the factor 2 is induced by spin

     complex(dp) :: ratio

     complex(Dp), intent(out) :: Hij(-ijmax:ijmax,Num_wann,Num_wann)

     Hij=0.0d0
     do iR=1,Nrpts
        ia=irvec(1,iR)
        ib=irvec(2,iR)

        !> new lattice
        call latticetransform(ia, ib, new_ia, new_ib)

        inew_ib= int(new_ib)
        if (abs(new_ib).le.ijmax)then
           kdotr=k*new_ia
           ratio=cos(2d0*pi*kdotr)+zi*sin(2d0*pi*kdotr)

           Hij(inew_ib, 1:Num_wann, 1:Num_wann )&
           =Hij(inew_ib, 1:Num_wann, 1:Num_wann )&
           +HmnR(:,:,iR)*ratio/ndegen(iR)
        endif

     enddo

  return
  end subroutine ham_qlayer2qlayer2

  !> use Umatrix to get the new representation of a vector in new basis
  !> R= a*R1+b*R2+c*R3= x*R1'+y*R2'+z*R3'
  subroutine latticetransform(a, b, x, y)
     use para
     implicit none

     integer, intent(in)  :: a, b
     real(dp), intent(out) :: x, y

     real(dp) :: Uinv(2, 2)

     Uinv= Umatrix

     call inv_r(2, Uinv)

     x= a*Uinv(1, 1)+ b*Uinv(2, 1)
     y= a*Uinv(1, 2)+ b*Uinv(2, 2)

     return
  end subroutine latticetransform
