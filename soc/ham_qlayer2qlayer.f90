! this subroutine is used to caculate Hamiltonian between
! slabs  
! 4/23/2010 by QS Wu

  subroutine ham_qlayer2qlayer(k,H00new,H01new)

     use para

     implicit none

! loop index
     integer :: i,j,iR

! index used to sign irvec     
     integer :: ia,ib,ic

! 

! new index used to sign irvec     
     real(dp) :: new_ia,new_ib,new_ic
     integer :: inew_ia,inew_ib,inew_ic

! wave vector k times lattice vector R  
     real(Dp) :: kdotr  

! input wave vector k's cooridinates
     real(Dp),intent(in) :: k(2)

! H00 Hamiltonian between nearest neighbour-quintuple-layers
! the factor 2 is induced by spin

     complex(dp) :: ratio

     complex(Dp) :: Hij(-ijmax:ijmax,Num_wann,Num_wann)


! H00 Hamiltonian between nearest neighbour-quintuple-layers
! the factor 2 is induced by spin


!     complex(Dp),allocatable,intent(out) :: H00new(:,:)
     complex(Dp),intent(out) :: H00new(Ndim,Ndim)

! H01 Hamiltonian between next-nearest neighbour-quintuple-layers
!     complex(Dp),allocatable,intent(out) :: H01new(:,:)
     complex(Dp),intent(out) :: H01new(Ndim,Ndim)


     Hij=0.0d0
     do iR=1,Nrpts
        ia=irvec(1,iR)
        ib=irvec(2,iR)
        ic=irvec(3,iR)

        call latticetransform(ia, ib, ic, new_ia, new_ib, new_ic)
       
        inew_ic= int(new_ic)
        if (abs(new_ic).le.ijmax)then
           kdotr=k(1)*new_ia+ k(2)*new_ib
           ratio=cos(2d0*pi*kdotr)+zi*sin(2d0*pi*kdotr)


           Hij(inew_ic, 1:Num_wann, 1:Num_wann )&
           =Hij(inew_ic, 1:Num_wann, 1:Num_wann )&
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
     integer :: ia,ib,ic

! 

! new index used to sign irvec     
     real(dp) :: new_ia,new_ib,new_ic
     integer :: inew_ia,inew_ib,inew_ic

! wave vector k times lattice vector R  
     real(Dp) :: kdotr  

! input wave vector k's cooridinates
     real(Dp),intent(in) :: k(2)

! H00 Hamiltonian between nearest neighbour-quintuple-layers
! the factor 2 is induced by spin

     complex(dp) :: ratio

     complex(Dp), intent(out) :: Hij(-ijmax:ijmax,Num_wann,Num_wann)

     Hij=0.0d0
     do iR=1,Nrpts
        ia=irvec(1,iR)
        ib=irvec(2,iR)
        ic=irvec(3,iR)

        !> new lattice
        call latticetransform(ia, ib, ic, new_ia, new_ib, new_ic)

        inew_ic= int(new_ic)
        if (abs(new_ic).le.ijmax)then
           kdotr=k(1)*new_ia+k(2)*new_ib
           ratio=cos(2d0*pi*kdotr)+zi*sin(2d0*pi*kdotr)

           Hij(inew_ic, 1:Num_wann, 1:Num_wann )&
           =Hij(inew_ic, 1:Num_wann, 1:Num_wann )&
           +HmnR(:,:,iR)*ratio/ndegen(iR)
        endif

     enddo

  return
  end subroutine ham_qlayer2qlayer2

  !> for slab with in-plane magnetic field
  !> the magnetic phase are chosen like this
  !> phi_ij= alpha*[By*(xj-xi)*(zi+zj)-Bx*(yj-yi)*(zi+zj)] 
  !> x, z are in unit of Angstrom, Bx, By are in unit of Tesla
  !> History :
  !        9/17/2015 by Quansheng Wu @ETH Zurich
  subroutine ham_qlayer2qlayer_parallel_B(k,Hij)

     use para

     implicit none

! loop index
     integer :: iR

! index used to sign irvec     
     integer :: ia,ib,ic

! new index used to sign irvec     
     real(dp) :: new_ia,new_ib,new_ic
     integer :: inew_ia,inew_ib,inew_ic

! wave vector k times lattice vector R  
     real(Dp) :: kdotr  

! input wave vector k's cooridinates
     real(Dp),intent(in) :: k(2)
     complex(dp) :: ratio
     
     real(dp) :: phase
     complex(dp) :: fac

     real(dp) :: xyz1(3)
     real(dp) :: xyz2(3)


     complex(Dp), intent(out) :: Hij(-ijmax:ijmax,Num_wann,Num_wann)

     Hij=0.0d0
     do iR=1,Nrpts
        ia=irvec(1,iR)
        ib=irvec(2,iR)
        ic=irvec(3,iR)

        !> new lattice
        call latticetransform(ia, ib, ic, new_ia, new_ib, new_ic)

        inew_ic= int(new_ic)
        if (abs(new_ic).gt.ijmax) cycle

        do ia1=1, Num_atoms
        do ia2=1, Num_atoms
           xyz1(1)= 0
           xyz1(2)= dble(i1)
           xyz1(3)= dble(j1)
           
           xyz2(1)= dble(ia)
           xyz2(2)= dble(i2)
           xyz2(3)= dble(j2)
 

           phase=-Bmag*(xyz2(3)+xyz1(3))*(xyz2(1)-xyz1(1))/2d0
           fac= cos(phase)+ zi*sin(phase)

           kdotr=k(1)*new_ia+k(2)*new_ib
           ratio=cos(2d0*pi*kdotr)+zi*sin(2d0*pi*kdotr)

           Hij(inew_ic, 1:Num_wann, 1:Num_wann )&
           =Hij(inew_ic, 1:Num_wann, 1:Num_wann )&
           +HmnR(:,:,iR)*ratio/ndegen(iR)* fac
        enddo ! ia2
        enddo ! ia1
     enddo ! iR

     return
  end subroutine ham_qlayer2qlayer_parallel_B




  subroutine ham_qlayer2qlayerribbon(k,Hij)

     use para

     implicit none

! loop index
     integer :: iR

! index used to sign irvec     
     integer :: ia,ib,ic

! new index used to sign irvec     
     real(dp) :: new_ia,new_ib,new_ic
     integer :: inew_ia,inew_ib,inew_ic

! wave vector k times lattice vector R  
     real(Dp) :: kdotr

! input wave vector k's cooridinates
     real(Dp),intent(in) :: k

     complex(dp) :: ratio

     complex(Dp), intent(out) :: Hij(-ijmax:ijmax, &
        -ijmax:ijmax,Num_wann,Num_wann)

     Hij=0.0d0
     do iR=1,Nrpts
        ia=irvec(1,iR)
        ib=irvec(2,iR)
        ic=irvec(3,iR)

        !> new lattice
        call latticetransform(ia, ib, ic, new_ia, new_ib, new_ic)


        inew_ia= int(new_ia)
        inew_ib= int(new_ib)
        if (abs(new_ia).le.ijmax)then
        if (abs(new_ib).le.ijmax)then
           kdotr=k*new_ic
           ratio=cos(2d0*pi*kdotr)+zi*sin(2d0*pi*kdotr)

           Hij(inew_ia, inew_ib, 1:Num_wann, 1:Num_wann )&
           =Hij(inew_ia, inew_ib, 1:Num_wann, 1:Num_wann )&
           +HmnR(:,:,iR)*ratio/ndegen(iR)
        endif
        endif

     enddo

     return
  end subroutine ham_qlayer2qlayerribbon


  !> use Umatrix to get the new representation of a vector in new basis
  !> R= a*R1+b*R2+c*R3= x*R1'+y*R2'+z*R3'
  subroutine latticetransform(a, b, c, x, y, z)
     use para
     implicit none

     integer, intent(in)  :: a, b, c
     real(dp), intent(out) :: x, y, z

     real(dp) :: Uinv(3, 3)

     Uinv= Umatrix

     call inv_r(3, Uinv)

     x= a*Uinv(1, 1)+ b*Uinv(2, 1)+ c*Uinv(3, 1)
     y= a*Uinv(1, 2)+ b*Uinv(2, 2)+ c*Uinv(3, 2)
     z= a*Uinv(1, 3)+ b*Uinv(2, 3)+ c*Uinv(3, 3)

     return
  end subroutine latticetransform
