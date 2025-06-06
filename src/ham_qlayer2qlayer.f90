  subroutine ham_qlayer2qlayer(k,H00new,H01new)
     ! This subroutine caculates Hamiltonian between
     ! slabs For surface state calculation
     ! 4/23/2010 by QS Wu
     ! Copyright (c) 2010 QuanSheng Wu. All rights reserved.

     use para

     implicit none

     ! loop index
     integer :: i,j,iR

     ! index used to sign irvec     
     real(dp) :: ia,ib,ic

     ! new index used to sign irvec     
     real(dp) :: new_ia,new_ib,new_ic
     integer :: inew_ic
     integer :: int_ic

     ! wave vector k times lattice vector R  
     real(Dp) :: kdotr  

     ! input wave vector k's cooridinates
     real(Dp),intent(in) :: k(2)

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

     allocate(Hij(Num_wann,Num_wann,-ijmax:ijmax))

     Hij=0.0d0
     do iR=1,nrpts_surfacecell
        ia=irvec_surfacecell(1,iR)
        ib=irvec_surfacecell(2,iR)
        ic=irvec_surfacecell(3,iR)

        int_ic= int(ic)
        if (abs(ic).le.ijmax)then
           kdotr=k(1)*ia+ k(2)*ib
           ratio=cos(2d0*pi*kdotr)+zi*sin(2d0*pi*kdotr)

           Hij(1:Num_wann, 1:Num_wann, int_ic )&
           =Hij(1:Num_wann, 1:Num_wann, int_ic)&
           +HmnR_surfacecell(:,:,iR)*ratio/ndegen_surfacecell(iR)
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
                =Hij(:,:,j-i)
        endif
     enddo
     enddo

     ! H01new
     do i=1,Np
     do j=Np+1,Np*2
        if (j-i.le.ijmax) then
           H01new(Num_wann*(i-1)+1:Num_wann*i,&
               Num_wann*(j-1-Np)+1:Num_wann*(j-Np))=Hij(:,:,j-i)
        endif
     enddo
     enddo

     do i=1,Ndim
     do j=1,Ndim
        if(abs(H00new(i,j)-conjg(H00new(j,i))).ge.1e-4)then
       !  write(stdout,*)'there are something wrong with ham_qlayer2qlayer'
       !stop
        endif

     enddo
     enddo

     deallocate(Hij)

  return
  end subroutine ham_qlayer2qlayer

  subroutine s_qlayer2qlayer(k,H00new,H01new)
   ! This subroutine caculates Hamiltonian between
   ! slabs For surface state calculation
   ! 4/23/2010 by QS Wu
   ! Copyright (c) 2010 QuanSheng Wu. All rights reserved.

   use para

   implicit none

   ! loop index
   integer :: i,j,iR

   ! index used to sign irvec     
   real(dp) :: ia,ib,ic

   ! new index used to sign irvec     
   real(dp) :: new_ia,new_ib,new_ic
   integer :: inew_ic
   integer :: int_ic

   ! wave vector k times lattice vector R  
   real(Dp) :: kdotr  

   ! input wave vector k's cooridinates
   real(Dp),intent(in) :: k(2)

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

   allocate(Hij(Num_wann,Num_wann,-ijmax:ijmax))

   Hij=0.0d0
   do iR=1,nrpts_surfacecell
      ia=irvec_surfacecell(1,iR)
      ib=irvec_surfacecell(2,iR)
      ic=irvec_surfacecell(3,iR)

      int_ic= int(ic)
      if (abs(ic).le.ijmax)then
         kdotr=k(1)*ia+ k(2)*ib
         ratio=cos(2d0*pi*kdotr)+zi*sin(2d0*pi*kdotr)

         Hij(1:Num_wann, 1:Num_wann, int_ic )&
         =Hij(1:Num_wann, 1:Num_wann, int_ic)&
         +SmnR_surfacecell(:,:,iR)*ratio/ndegen_surfacecell(iR)
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
              =Hij(:,:,j-i)
      endif
   enddo
   enddo

   ! H01new
   do i=1,Np
   do j=Np+1,Np*2
      if (j-i.le.ijmax) then
         H01new(Num_wann*(i-1)+1:Num_wann*i,&
             Num_wann*(j-1-Np)+1:Num_wann*(j-Np))=Hij(:,:,j-i)
      endif
   enddo
   enddo

   do i=1,Ndim
   do j=1,Ndim
      if(abs(H00new(i,j)-conjg(H00new(j,i))).ge.1e-4)then
     !  write(stdout,*)'there are something wrong with ham_qlayer2qlayer'
     !stop
      endif

   enddo
   enddo

   deallocate(Hij)

 return
 end subroutine s_qlayer2qlayer


  subroutine ham_qlayer2qlayer_bak(k,H00new,H01new)
     ! This subroutine caculates Hamiltonian between
     ! slabs For surface state calculation
     ! 4/23/2010 by QS Wu
     ! Copyright (c) 2010 QuanSheng Wu. All rights reserved.

     use para

     implicit none

     ! loop index
     integer :: i,j,iR

     ! index used to sign irvec     
     real(dp) :: ia,ib,ic

     ! new index used to sign irvec     
     real(dp) :: new_ia,new_ib,new_ic
     integer :: inew_ic

     ! wave vector k times lattice vector R  
     real(Dp) :: kdotr  

     ! input wave vector k's cooridinates
     real(Dp),intent(in) :: k(2)

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
       !  write(stdout,*)'there are something wrong with ham_qlayer2qlayer'
       !stop
        endif

     enddo
     enddo

     deallocate(Hij)

  return
  end subroutine ham_qlayer2qlayer_bak

  subroutine ham_qlayer2qlayer_LOTO(k,H00new,H01new)
     ! This subroutine caculates Hamiltonian between
     ! slabs For surface state calculation
     ! Aug/01/2017 by QS Wu
     ! Copyright (c) 2017 QuanSheng Wu. All rights reserved.

     use para

     implicit none

     ! loop index
     integer :: i,j,iR

     ! index used to sign irvec     
     real(dp) :: ia,ib,ic
     integer  :: ii,jj,pp,qq

     ! new index used to sign irvec     
     real(dp) :: new_ia,new_ib,new_ic
     integer :: inew_ic

     ! wave vector k times lattice vector R  
     real(Dp) :: kdotr  

     ! input wave vector k's cooridinates
     real(Dp),intent(in) :: k(2)

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

     !> see eqn. (3) in J. Phys.: Condens. Matter 22 (2010) 202201
     complex(Dp),allocatable :: nac_correction(:, :, :) 

     real(dp) :: temp1(2), temp2, knew(2)
     real(dp) :: constant_t
     complex(dp), allocatable :: mat1(:, :)
     complex(dp), allocatable :: mat2(:, :)
 
     !> k times Born charge
     real(dp), allocatable :: kBorn(:, :)

     real(dp) :: nac_q

    
     allocate(kBorn(Origin_cell%Num_atoms, 3))
     allocate(mat1(Num_wann, Num_wann))
     allocate(mat2(Num_wann, Num_wann))
     allocate(nac_correction(Num_wann, Num_wann, Nrpts))
     allocate(Hij(-ijmax:ijmax,Num_wann,Num_wann))
     mat1 = 0d0
     mat2 = 0d0
     nac_correction= 0d0

     knew= k 

     !>  add loto splitting term
     temp1(1:2)= (/0.0,0.0/)
     temp2= 0.0
     if (abs((k(1)**2+k(2)**2)).le.0.0001) then  !> skip k=0
        knew(1)= 0.0001d0
        knew(2)= 0.0001d0
     endif

     !> see eqn. (3) in J. Phys.: Condens. Matter 22 (2010) 202201
     do qq= 1, 2
        temp1(qq)= knew(1)*Diele_Tensor(qq, 1)+knew(2)*Diele_Tensor(qq, 2) 
     enddo
     temp2= knew(1)*temp1(1)+ knew(2)*temp1(2)
     constant_t= 4d0*3.1415926d0/(temp2*Origin_cell%CellVolume)*VASPToTHZ/0.529177

     do ii=1, Origin_cell%Num_atoms
        do pp=1, 3
           kBorn(ii, pp)=  knew(1)*Born_Charge(ii,1,pp)+knew(2)*Born_Charge(ii,2,pp)
        enddo
     enddo
     
     nac_correction= 0d0
 
     Hij=0.0d0
     do iR=1,Nrpts
        ia=irvec(1,iR)
        ib=irvec(2,iR)
        ic=irvec(3,iR)

        !> new lattice
        call latticetransform(ia, ib, ic, new_ia, new_ib, new_ic)

        inew_ic= int(new_ic)
        if (abs(new_ic).le.ijmax)then
           kdotr=knew(1)*new_ia+knew(2)*new_ib
           ratio=cos(2d0*pi*kdotr)+zi*sin(2d0*pi*kdotr)

           do ii= 1,Origin_cell%Num_atoms
              do pp= 1, 3
                 do jj= 1, Origin_cell%Num_atoms
                    do qq= 1, 3
                       nac_q= kBorn(jj, qq)*kBorn(ii, pp)*constant_t/sqrt(Atom_Mass(ii)*Atom_Mass(jj))
                       nac_correction((ii-1)*3+pp,(jj-1)*3+qq,iR) = HmnR((ii-1)*3+pp,(jj-1)*3+qq,iR) + dcmplx(nac_q,0)
                    enddo  ! qq
                 enddo  ! jj
              enddo ! pp
           enddo  ! i
           Hij(inew_ic, 1:Num_wann, 1:Num_wann )&
           = Hij(inew_ic, 1:Num_wann, 1:Num_wann )&
           + nac_correction(:, :, iR)*ratio/ndegen(iR)
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
       !  write(stdout,*)'there are something wrong with ham_qlayer2qlayer'
       !stop
        endif

     enddo
     enddo

     deallocate(kBorn)
     deallocate(mat1)
     deallocate(mat2)
     deallocate(nac_correction)
     deallocate(Hij)

  return
  end subroutine ham_qlayer2qlayer_LOTO


  subroutine ham_qlayer2qlayer_velocity(k,Vij_x,Vij_y)
     ! This subroutine caculates velocity matrix between
     ! slabs  
     ! 06/Aug/2018 by QS Wu
     ! Copyright (c) 2018 QuanSheng Wu. All rights reserved.

     use para
     implicit none

     ! loop index
     real(dp) :: ia, ib, ic
     integer :: iR,  inew_ic

     ! new index used to sign irvec     
     real(dp) :: new_ia,new_ib,new_ic

     ! wave vector k times lattice vector R  
     real(Dp) :: kdotr, r0(3), r1(3)

     ! input wave vector k's cooridinates
     real(Dp),intent(in) :: k(2)

     complex(dp) :: ratio

     complex(Dp), intent(out) :: Vij_x(-ijmax:ijmax,Num_wann,Num_wann)
     complex(Dp), intent(out) :: Vij_y(-ijmax:ijmax,Num_wann,Num_wann)

     Vij_x=0.0d0; Vij_y= 0d0;
     do iR=1,Nrpts
        ia=irvec(1,iR)
        ib=irvec(2,iR)
        ic=irvec(3,iR)

        !> new lattice
        call latticetransform(ia, ib, ic, new_ia, new_ib, new_ic)
        r0= new_ia*Rua_newcell+ new_ib* Rub_newcell

        !> rotate the vector from the original coordinate to the new coordinates
        !> which defined like this: x is along R1', z is along R1'xR2', y is along z x y
        call rotate(r0, r1)

        inew_ic= int(new_ic)
        if (abs(new_ic).le.ijmax)then
           kdotr= k(1)*new_ia+ k(2)*new_ib
           ratio= cos(2d0*pi*kdotr)+ zi*sin(2d0*pi*kdotr)

           Vij_x(inew_ic, 1:Num_wann, 1:Num_wann ) &
           = Vij_x(inew_ic, 1:Num_wann, 1:Num_wann ) &
           + zi*r1(1)*HmnR(:,:,iR)*ratio/ndegen(iR)
           Vij_y(inew_ic, 1:Num_wann, 1:Num_wann ) &
           = Vij_y(inew_ic, 1:Num_wann, 1:Num_wann ) &
           + zi*r1(2)*HmnR(:,:,iR)*ratio/ndegen(iR)
        endif
     enddo

     return
  end subroutine ham_qlayer2qlayer_velocity


  subroutine ham_qlayer2qlayer2(k,Hij)
     ! This subroutine caculates Hamiltonian between
     ! slabs  
     ! 4/23/2010 by QS Wu
     ! Copyright (c) 2010 QuanSheng Wu. All rights reserved.

     use para

     implicit none

     ! loop index
     integer :: iR, inew_ic, int_ic
     real(dp) :: ia, ib, ic

     ! new index used to sign irvec     
     real(dp) :: new_ia,new_ib,new_ic

     ! wave vector k times lattice vector R  
     real(Dp) :: kdotr  

     ! input wave vector k's cooridinates
     real(Dp),intent(in) :: k(2)

     complex(dp) :: ratio

     complex(Dp), intent(out) :: Hij(-ijmax:ijmax,Num_wann,Num_wann)

     Hij=0.0d0

     do iR=1,nrpts_surfacecell
        ia=irvec_surfacecell(1,iR)
        ib=irvec_surfacecell(2,iR)
        ic=irvec_surfacecell(3,iR)

        int_ic= int(ic)
        if (abs(ic).le.ijmax)then
           kdotr=k(1)*ia+k(2)*ib
           ratio=cos(2d0*pi*kdotr)+zi*sin(2d0*pi*kdotr)

           Hij(int_ic, 1:Num_wann, 1:Num_wann )&
           =Hij(int_ic, 1:Num_wann, 1:Num_wann )&
           +HmnR_surfacecell(:,:,iR)*ratio/ndegen_surfacecell(iR)

        endif

     enddo

  return
  end subroutine ham_qlayer2qlayer2

  subroutine s_qlayer2qlayer2(k,Hij)
   ! This subroutine caculates Hamiltonian between
   ! slabs  
   ! 4/23/2010 by QS Wu
   ! Copyright (c) 2010 QuanSheng Wu. All rights reserved.

   use para

   implicit none

   ! loop index
   integer :: iR, inew_ic, int_ic
   real(dp) :: ia, ib, ic

   ! new index used to sign irvec     
   real(dp) :: new_ia,new_ib,new_ic

   ! wave vector k times lattice vector R  
   real(Dp) :: kdotr  

   ! input wave vector k's cooridinates
   real(Dp),intent(in) :: k(2)

   complex(dp) :: ratio

   complex(Dp), intent(out) :: Hij(-ijmax:ijmax,Num_wann,Num_wann)

   Hij=0.0d0

   do iR=1,nrpts_surfacecell
      ia=irvec_surfacecell(1,iR)
      ib=irvec_surfacecell(2,iR)
      ic=irvec_surfacecell(3,iR)

      int_ic= int(ic)
      if (abs(ic).le.ijmax)then
         kdotr=k(1)*ia+k(2)*ib
         ratio=cos(2d0*pi*kdotr)+zi*sin(2d0*pi*kdotr)

         Hij(int_ic, 1:Num_wann, 1:Num_wann )&
         =Hij(int_ic, 1:Num_wann, 1:Num_wann )&
         +SmnR_surfacecell(:,:,iR)*ratio/ndegen_surfacecell(iR)

      endif

   enddo

return
end subroutine s_qlayer2qlayer2

  subroutine ham_qlayer2qlayer2_LOTO(k,Hij)
     ! This subroutine caculates Hamiltonian between
     ! slabs for phonon system with LOTO coorection
     ! 4/23/2010 by QS Wu
     ! Copyright (c) 2010 QuanSheng Wu. All rights reserved.

     use para

     implicit none

     ! loop index
     integer :: iR

     ! index used to sign irvec     
     real(dp) :: ia,ib,ic
     integer  :: ii,jj,pp,qq

     ! 

     ! new index used to sign irvec     
     integer  :: inew_ic
     real(dp) :: new_ia,new_ib,new_ic

     ! wave vector k times lattice vector R  
     real(dp) :: kdotr  

     ! input wave vector k's cooridinates
     real(dp),intent(in) :: k(2)
    
     complex(dp) :: ratio

     ! H00 Hamiltonian between nearest neighbour-quintuple-layers
     ! the factor 2 is induced by spin
     complex(dp), intent(out) :: Hij(-ijmax:ijmax,Num_wann,Num_wann)

     !> see eqn. (3) in J. Phys.: Condens. Matter 22 (2010) 202201
     complex(Dp),allocatable :: nac_correction(:, :, :) 

     real(dp) :: temp1(2), temp2, knew(2)
     real(dp) ::  constant_t
     complex(dp), allocatable :: mat1(:, :)
     complex(dp), allocatable :: mat2(:, :)
 
     !> k times Born charge
     real(dp), allocatable :: kBorn(:, :)

     real(dp) :: nac_q

     allocate(kBorn(Origin_cell%Num_atoms, 3))
     allocate(mat1(Num_wann, Num_wann))
     allocate(mat2(Num_wann, Num_wann))
     allocate(nac_correction(Num_wann, Num_wann, Nrpts))
     mat1 = 0d0
     mat2 = 0d0
     nac_correction= 0d0

     knew= k 

     !>  add loto splitting term
     temp1(1:2)= (/0.0,0.0/)
     temp2= 0.0
     if (abs((k(1)**2+k(2)**2)).le.0.0001) then  !> skip k=0
        knew(1)= 0.0001d0
        knew(2)= 0.0001d0
     endif

     !> see eqn. (3) in J. Phys.: Condens. Matter 22 (2010) 202201
     do qq= 1, 2
        temp1(qq)= knew(1)*Diele_Tensor(qq, 1)+knew(2)*Diele_Tensor(qq, 2) 
     enddo
     temp2= knew(1)*temp1(1)+ knew(2)*temp1(2)
     constant_t= 4d0*3.1415926d0/(temp2*Origin_cell%CellVolume)*VASPToTHZ/0.529177

     do ii=1, Origin_cell%Num_atoms
        do pp=1, 3
           kBorn(ii, pp)=  knew(1)*Born_Charge(ii,1,pp)+knew(2)*Born_Charge(ii,2,pp)
        enddo
     enddo
     
     nac_correction= 0d0
 
     Hij=0.0d0
     do iR=1,Nrpts
        ia=irvec(1,iR)
        ib=irvec(2,iR)
        ic=irvec(3,iR)

        !> new lattice
        call latticetransform(ia, ib, ic, new_ia, new_ib, new_ic)

        inew_ic= int(new_ic)
        if (abs(new_ic).le.ijmax)then
           kdotr=knew(1)*new_ia+knew(2)*new_ib
           ratio=cos(2d0*pi*kdotr)+zi*sin(2d0*pi*kdotr)

           do ii= 1,Origin_cell%Num_atoms
              do pp= 1, 3
                 do jj= 1, Origin_cell%Num_atoms
                    do qq= 1, 3
                       nac_q= kBorn(jj, qq)*kBorn(ii, pp)*constant_t/sqrt(Atom_Mass(ii)*Atom_Mass(jj))
                       nac_correction((ii-1)*3+pp,(jj-1)*3+qq,iR) = HmnR((ii-1)*3+pp,(jj-1)*3+qq,iR) + dcmplx(nac_q,0)
                    enddo  ! qq
                 enddo  ! jj
              enddo ! pp
           enddo  ! i
           Hij(inew_ic, 1:Num_wann, 1:Num_wann )&
           = Hij(inew_ic, 1:Num_wann, 1:Num_wann )&
           + nac_correction(:, :, iR)*ratio/ndegen(iR)
        endif

     enddo

     deallocate(kBorn)
     deallocate(mat1)
     deallocate(mat2)
     deallocate(nac_correction)
  return
  end subroutine ham_qlayer2qlayer2_LOTO

  subroutine ham_qlayer2qlayerribbon(k,Hij)
     ! This subroutine caculates Hamiltonian between
     ! ribbon
     ! 4/23/2010 by QS Wu
     ! Copyright (c) 2010 QuanSheng Wu. All rights reserved.

     use para

     implicit none

     ! loop index
     integer :: iR

     ! index used to sign irvec     
     real(dp) :: ia,ib,ic

     ! new index used to sign irvec     
     real(dp) :: new_ia,new_ib,new_ic
     integer :: inew_ia,inew_ib

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


  subroutine latticetransform(a, b, c, x, y, z)
     !> use Umatrix to get the new representation of a vector in new basis
     !> R= a*R1+b*R2+c*R3= x*R1'+y*R2'+z*R3'
     use para
     implicit none

     real(dp), intent(in)  :: a, b, c
     real(dp), intent(out) :: x, y, z

     real(dp), allocatable :: Uinv(:, :)

     allocate(Uinv(3, 3))
     Uinv= Umatrix

     call inv_r(3, Uinv)

     x= a*Uinv(1, 1)+ b*Uinv(2, 1)+ c*Uinv(3, 1)
     y= a*Uinv(1, 2)+ b*Uinv(2, 2)+ c*Uinv(3, 2)
     z= a*Uinv(1, 3)+ b*Uinv(2, 3)+ c*Uinv(3, 3)

     return
  end subroutine latticetransform
  
    subroutine transform_r(v0, x, y, z)
     !> use Umatrix to get the new representation of a vector in new basis
     !> R= a*R1+b*R2+c*R3= x*R1'+y*R2'+z*R3'
     use para
     implicit none

     real(dp), intent(in)  :: v0(3)
     real(dp), intent(out) :: x, y, z

     real(dp), allocatable :: Uinv(:, :)

     allocate(Uinv(3, 3))
     Uinv= Umatrix

     call inv_r(3, Uinv)

     x= v0(1)*Uinv(1, 1)+ v0(2)*Uinv(2, 1)+ v0(3)*Uinv(3, 1)
     y= v0(1)*Uinv(1, 2)+ v0(2)*Uinv(2, 2)+ v0(3)*Uinv(3, 2)
     z= v0(1)*Uinv(1, 3)+ v0(2)*Uinv(2, 3)+ v0(3)*Uinv(3, 3)

     return
  end subroutine 
