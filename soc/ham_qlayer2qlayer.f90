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
     integer :: new_ia,new_ib,new_ic

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
     complex(Dp),intent(out) :: H00new(Num_wann*Nslab,Num_wann*Nslab)

! H01 Hamiltonian between next-nearest neighbour-quintuple-layers
!     complex(Dp),allocatable,intent(out) :: H01new(:,:)
     complex(Dp),intent(out) :: H01new(Num_wann*Nslab,Num_wann*Nslab)



     Hij=0.0d0
     do iR=1,Nrpts
        ia=irvec(1,iR)
        ib=irvec(2,iR)
        ic=irvec(3,iR)

        ! 100 for fcc conventional cell
        new_ia=ia
        new_ib=ib
        new_ic=ib+ ic
        
        if (abs(new_ic).le.ijmax)then
           kdotr=k(1)*real(new_ia,Dp)+k(2)*real(new_ib,Dp)
           ratio=cos(2d0*pi*kdotr)+zi*sin(2d0*pi*kdotr)

           Hij(new_ic, 1:Num_wann, 1:Num_wann )&
           =Hij(new_ic, 1:Num_wann, 1:Num_wann )&
           +HmnR(:,:,iR)*ratio/ndegen(iR)
        endif

     enddo

     H00new=0.0d0
     H01new=0.0d0

! nslab's principle layer 
! H00new
     do i=1,Nslab
     do j=1,Nslab
        if (abs(i-j).le.(ijmax)) then
          H00new(Num_wann*(i-1)+1:Num_wann*i,Num_wann*(j-1)+1:Num_wann*j)&
                =Hij(j-i,:,:)
        endif
     enddo
     enddo

! H01new
     do i=1,Nslab
     do j=Nslab+1,Nslab*2
        if (j-i.le.ijmax) then
           H01new(Num_wann*(i-1)+1:Num_wann*i,&
               Num_wann*(j-1-Nslab)+1:Num_wann*(j-Nslab))=Hij(j-i,:,:)
        endif
     enddo
     enddo

     do i=1,Num_wann*Nslab
     do j=1,Num_wann*Nslab
        if(abs(H00new(i,j)-conjg(H00new(j,i))).ge.1e-4)then
          write(*,*)'there are something wrong with ham_qlayer2qlayer'
        stop
        endif

     enddo
     enddo

  return
  end

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
     integer :: new_ia,new_ib,new_ic

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

        !> 100 for fcc conventional cell
        new_ia=ia
        new_ib=ib
        new_ic=ib+ ic
        
       !new_ia=ia
       !new_ib=ib
       !new_ic=ic    
        if (abs(new_ic).le.ijmax)then
           kdotr=k(1)*real(new_ia,Dp)+k(2)*real(new_ib,Dp)
           ratio=cos(2d0*pi*kdotr)+zi*sin(2d0*pi*kdotr)

           Hij(new_ic, 1:Num_wann, 1:Num_wann )&
           =Hij(new_ic, 1:Num_wann, 1:Num_wann )&
           +HmnR(:,:,iR)*ratio/ndegen(iR)
        endif

     enddo

  return
  end subroutine ham_qlayer2qlayer2

  subroutine ham_qlayer2qlayerribbon(k,Hij)

     use para

     implicit none

! loop index
     integer :: iR

! index used to sign irvec     
     integer :: ia,ib,ic

! new index used to sign irvec     
     integer :: new_ia,new_ib,new_ic

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

        !> 100 for fcc conventional cell
        new_ia=ia
        new_ib=ib
        new_ic=ib+ ic
        

        if (abs(new_ib).le.ijmax)then
        if (abs(new_ic).le.ijmax)then
           kdotr=k*real(new_ia,Dp)
           ratio=cos(2d0*pi*kdotr)+zi*sin(2d0*pi*kdotr)

           Hij(new_ib, new_ic, 1:Num_wann, 1:Num_wann )&
           =Hij(new_ib, new_ic, 1:Num_wann, 1:Num_wann )&
           +HmnR(:,:,iR)*ratio/ndegen(iR)
        endif
        endif

     enddo

     return
  end subroutine ham_qlayer2qlayerribbon


