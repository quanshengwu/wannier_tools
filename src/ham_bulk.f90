! This subroutine is used to caculate Hamiltonian for
! bulk system .

! History
!        May/29/2011 by Quansheng Wu
!        Aug/31/2018 by Quansheng Wu
subroutine ham_bulk(k,Hamk_bulk)

   use para
   implicit none

   real(dp), intent(in) :: k(3)

   ! Hamiltonian of bulk system
   complex(Dp),intent(out) ::Hamk_bulk(Num_wann, Num_wann)

   ! temporary integers
   integer :: i1,i2,ia,ib,ic,iR, ia1, ia2, nwann

   real(dp) :: kdotr, R(3)

   !> row start; row end; column start; column end
   integer :: rs, re, cs, ce

   real(dp) :: pos(3), pos1(3), pos2(3), pos_cart(3), pos_direct(3)

   !> distance between two atoms
   real(dp) :: dis

   real(dp), external :: norm
   complex(dp) :: ratio

   nwann= Num_wann/2
   Hamk_bulk=0d0
   if (soc>0 ) then
      !> the first atom in home unit cell
      do ia1=1, Num_atoms
         pos1= Atom_position_direct(:, ia1)
         rs= index_start(ia1)
         re= index_end(ia1)

         !> the second atom in unit cell R
         do ia2=1, Num_atoms
            pos2= Atom_position_direct(:, ia2)
            cs= index_start(ia2)
            ce= index_end(ia2)
            do iR=1,Nrpts

               pos_direct= irvec(:, iR)
               pos_direct= pos_direct+ pos2- pos1

               call direct_cart_real(pos_direct, pos_cart)

               dis= norm(pos_cart)
               if (dis> Rcut) cycle

               kdotr=k(1)*pos_direct(1) + k(2)*pos_direct(2) + k(3)*pos_direct(3)
               ratio= exp(pi2zi*kdotr)/ndegen(iR)

               Hamk_bulk(rs:re, cs:ce)=Hamk_bulk(rs:re, cs:ce)&
                  +HmnR(rs:re, cs:ce,iR)*ratio

               Hamk_bulk(rs:re, cs+nwann:ce+nwann)= &
                  Hamk_bulk(rs:re, cs+nwann:ce+nwann)+ &
                  HmnR(rs:re, cs+nwann:ce+nwann,iR)* &
                  ratio

               Hamk_bulk(rs+nwann:re+nwann, cs:ce)= &
                  Hamk_bulk(rs+nwann:re+nwann, cs:ce)+ &
                  HmnR(rs+nwann:re+nwann, cs:ce,iR)* &
                  ratio

               Hamk_bulk(rs+nwann:re+nwann, cs+nwann:ce+nwann)= &
                  Hamk_bulk(rs+nwann:re+nwann, cs+nwann:ce+nwann)+ &
                  HmnR(rs+nwann:re+nwann, cs+nwann:ce+nwann,iR)* &
                  ratio

            enddo ! iR
         enddo ! ia2
      enddo ! ia1
   else
      !> the first atom in home unit cell
      do ia1=1, Num_atoms
         pos1= Atom_position_direct(:, ia1)
         rs= index_start(ia1)
         re= index_end(ia1)

         !> the second atom in unit cell R
         do ia2=1, Num_atoms
            pos2= Atom_position_direct(:, ia2)
            cs= index_start(ia2)
            ce= index_end(ia2)
            do iR=1,Nrpts

               pos_direct= irvec(:, iR)
               pos_direct= pos_direct+ pos2- pos1

               call direct_cart_real(pos_direct, pos_cart)

               dis= norm(pos_cart)
               if (dis> Rcut) cycle

               kdotr=k(1)*pos_direct(1) + k(2)*pos_direct(2) + k(3)*pos_direct(3)
               ratio= exp(pi2zi*kdotr)/ndegen(iR)

               Hamk_bulk(rs:re, cs:ce)=Hamk_bulk(rs:re, cs:ce)&
                 +HmnR(rs:re, cs:ce,iR)*ratio

            enddo ! iR
         enddo ! ia2
      enddo ! ia1

   endif ! soc

! check hermitcity

   do i1=1, Num_wann
      do i2=1, Num_wann
         if(abs(Hamk_bulk(i1,i2)-conjg(Hamk_bulk(i2,i1))).ge.1e-6)then
            write(stdout,*)'there is something wrong with Hamk_bulk'
            write(stdout,*)'i1, i2', i1, i2
            write(stdout,*)'value at (i1, i2)', Hamk_bulk(i1, i2)
            write(stdout,*)'value at (i2, i1)', Hamk_bulk(i2, i1)
            !stop
         endif
      enddo
   enddo

   return
end subroutine ham_bulk


subroutine dHdk_atomicgauge(k, vx, vy, vz)
   use para, only : Nrpts, irvec, crvec, Atom_position, HmnR, ndegen, &
      Atom_position_direct, Num_wann, dp, Rcut, pi2zi, index_start, index_end, &
      zi, soc, Num_atoms
   implicit none

   !> momentum in 3D BZ
   real(dp), intent(in) :: k(3)

   !> velocity operator using lattice guage
   !> which don't take into account the atom's position
   !> this is a nature choice, while maybe not consistent with the symmetry
   complex(dp), intent(out) :: vx(Num_wann, Num_wann)
   complex(dp), intent(out) :: vy(Num_wann, Num_wann)
   complex(dp), intent(out) :: vz(Num_wann, Num_wann)

   integer :: iR, ia1, ia2, rs, re, cs, ce, nwann

   real(dp) :: pos(3), pos1(3), pos2(3), pos_cart(3), pos_direct(3)
   real(dp) :: kdotr, dis
   complex(dp) :: ratio
   real(dp), external :: norm

   nwann= Num_wann/2
   vx=0d0; vy=0d0; vz=0d0
   if (soc>0 ) then
      !> the first atom in home unit cell
      do ia1=1, Num_atoms
         pos1= Atom_position_direct(:, ia1)
         rs= index_start(ia1)
         re= index_end(ia1)

         !> the second atom in unit cell R
         do ia2=1, Num_atoms
            pos2= Atom_position_direct(:, ia2)
            cs= index_start(ia2)
            ce= index_end(ia2)
            do iR=1,Nrpts

               pos_direct= irvec(:, iR)
               pos_direct= pos_direct+ pos2- pos1

               call direct_cart_real(pos_direct, pos_cart)

               dis= norm(pos_cart)
               if (dis> Rcut) cycle

               kdotr=k(1)*pos_direct(1) + k(2)*pos_direct(2) + k(3)*pos_direct(3)
               ratio= exp(pi2zi*kdotr)/ndegen(iR)

               vx(rs:re, cs:ce)=vx(rs:re, cs:ce)+ &
                  zi*pos_cart(1)*HmnR(rs:re, cs:ce,iR)*ratio

               vx(rs:re, cs+nwann:ce+nwann)= &
                  vx(rs:re, cs+nwann:ce+nwann)+ &
                  zi*pos_cart(1)*HmnR(rs:re, cs+nwann:ce+nwann,iR)* ratio

               vx(rs+nwann:re+nwann, cs:ce)= &
                  vx(rs+nwann:re+nwann, cs:ce)+ &
                  zi*pos_cart(1)*HmnR(rs+nwann:re+nwann, cs:ce,iR)* ratio

               vx(rs+nwann:re+nwann, cs+nwann:ce+nwann)= &
                  vx(rs+nwann:re+nwann, cs+nwann:ce+nwann)+ &
                  zi*pos_cart(1)*HmnR(rs+nwann:re+nwann, cs+nwann:ce+nwann,iR)* ratio


               vy(rs:re, cs:ce)=vy(rs:re, cs:ce) + &
                  zi*pos_cart(2)*HmnR(rs:re, cs:ce,iR)*ratio

               vy(rs:re, cs+nwann:ce+nwann)= &
                  vy(rs:re, cs+nwann:ce+nwann)+ &
                  zi*pos_cart(2)*HmnR(rs:re, cs+nwann:ce+nwann,iR)* ratio

               vy(rs+nwann:re+nwann, cs:ce)= &
                  vy(rs+nwann:re+nwann, cs:ce)+ &
                  zi*pos_cart(2)*HmnR(rs+nwann:re+nwann, cs:ce,iR)* ratio

               vy(rs+nwann:re+nwann, cs+nwann:ce+nwann)= &
                  vy(rs+nwann:re+nwann, cs+nwann:ce+nwann)+ &
                  zi*pos_cart(2)*HmnR(rs+nwann:re+nwann, cs+nwann:ce+nwann,iR)* ratio


               vz(rs:re, cs:ce)=vz(rs:re, cs:ce)+ &
                  zi*pos_cart(3)*HmnR(rs:re, cs:ce,iR)*ratio

               vz(rs:re, cs+nwann:ce+nwann)= &
                  vz(rs:re, cs+nwann:ce+nwann)+ &
                  zi*pos_cart(3)*HmnR(rs:re, cs+nwann:ce+nwann,iR)* ratio

               vz(rs+nwann:re+nwann, cs:ce)= &
                  vz(rs+nwann:re+nwann, cs:ce)+ &
                  zi*pos_cart(3)*HmnR(rs+nwann:re+nwann, cs:ce,iR)* ratio

               vz(rs+nwann:re+nwann, cs+nwann:ce+nwann)= &
                  vz(rs+nwann:re+nwann, cs+nwann:ce+nwann)+ &
                  zi*pos_cart(3)*HmnR(rs+nwann:re+nwann, cs+nwann:ce+nwann,iR)* ratio

            enddo ! iR
         enddo ! ia2
      enddo ! ia1
   else
      !> the first atom in home unit cell
      do ia1=1, Num_atoms
         pos1= Atom_position_direct(:, ia1)
         rs= index_start(ia1)
         re= index_end(ia1)

         !> the second atom in unit cell R
         do ia2=1, Num_atoms
            pos2= Atom_position_direct(:, ia2)
            cs= index_start(ia2)
            ce= index_end(ia2)
            do iR=1,Nrpts

               pos_direct= irvec(:, iR)
               pos_direct= pos_direct+ pos2- pos1

               call direct_cart_real(pos_direct, pos_cart)

               dis= norm(pos_cart)
               if (dis> Rcut) cycle

               kdotr=k(1)*pos_direct(1) + k(2)*pos_direct(2) + k(3)*pos_direct(3)
               ratio= exp(pi2zi*kdotr)/ndegen(iR)

               vx(rs:re, cs:ce)=vx(rs:re, cs:ce)+ &
                  zi*pos_cart(1)*HmnR(rs:re, cs:ce,iR)*ratio

               vy(rs:re, cs:ce)=vy(rs:re, cs:ce)+ &
                  zi*pos_cart(2)*HmnR(rs:re, cs:ce,iR)*ratio

               vz(rs:re, cs:ce)=vz(rs:re, cs:ce)+ &
                  zi*pos_cart(3)*HmnR(rs:re, cs:ce,iR)*ratio

            enddo ! iR
         enddo ! ia2
      enddo ! ia1

   endif ! soc


   return
end subroutine dHdk_atomicgauge


subroutine dHdk_latticegauge(k, vx, vy, vz)
   use para, only : Nrpts, irvec, crvec, Atom_position, &
      HmnR, ndegen, Num_wann, zi, pi2zi, dp
   implicit none

   !> momentum in 3D BZ
   real(dp), intent(in) :: k(3)

   !> velocity operator using lattice guage
   !> which don't take into account the atom's position
   !> this is a nature choice, while maybe not consistent with the symmetry
   complex(dp), intent(out) :: vx(Num_wann, Num_wann)
   complex(dp), intent(out) :: vy(Num_wann, Num_wann)
   complex(dp), intent(out) :: vz(Num_wann, Num_wann)

   integer :: iR

   real(dp) :: kdotr
   complex(dp) :: ratio

   do iR= 1, Nrpts
      kdotr= k(1)*irvec(1,iR) + k(2)*irvec(2,iR) + k(3)*irvec(3,iR)
      ratio= Exp(pi2zi*kdotr)
      vx= vx+ zi*crvec(1, iR)*HmnR(:,:,iR)*ratio/ndegen(iR)
      vy= vy+ zi*crvec(2, iR)*HmnR(:,:,iR)*ratio/ndegen(iR)
      vz= vz+ zi*crvec(3, iR)*HmnR(:,:,iR)*ratio/ndegen(iR)
   enddo ! iR

   return
end subroutine dHdk_latticegauge


  subroutine ham_bulk_latticegauge(k,Hamk_bulk)
     ! This subroutine caculates Hamiltonian for
     ! bulk system without the consideration of the atom's position
     !
     ! History
     !
     !        May/29/2011 by Quansheng Wu

     use para
     implicit none

! loop index
     integer :: i1,i2,ia,ib,ic,iR
     integer :: nwann

     real(dp) :: kdotr

! wave vector in 2d
     real(Dp) :: k(3)

     ! coordinates of R vector
     real(Dp) :: R(3), R1(3), R2(3)

     complex(dp) :: factor

! Hamiltonian of bulk system
     complex(Dp),intent(out) ::Hamk_bulk(Num_wann, Num_wann)
!    complex(dp), allocatable :: mat1(:, :)
!    complex(dp), allocatable :: mat2(:, :)
!    allocate(mat1(Num_wann, Num_wann))
!    allocate(mat2(Num_wann, Num_wann))

     Hamk_bulk=0d0
     do iR=1, Nrpts
        ia=irvec(1,iR)
        ib=irvec(2,iR)
        ic=irvec(3,iR)

        R(1)=dble(ia)
        R(2)=dble(ib)
        R(3)=dble(ic)
        kdotr=k(1)*R (1) + k(2)*R (2) + k(3)*R (3)
        factor= exp(pi2zi*kdotr)

        Hamk_bulk(:, :)=&
        Hamk_bulk(:, :) &
        + HmnR(:, :, iR)*factor/ndegen(iR)
     enddo ! iR

    !call mat_mul(Num_wann, mirror_z, Hamk_bulk, mat1)
    !call mat_mul(Num_wann, mat1, mirror_z, mat2)
    !Hamk_bulk= (Hamk_bulk+ mat2)/2d0

     ! check hermitcity
     do i1=1, Num_wann
     do i2=1, Num_wann
        if(abs(Hamk_bulk(i1,i2)-conjg(Hamk_bulk(i2,i1))).ge.1e-6)then
          write(stdout,*)'there are something wrong with Hamk_bulk'
          write(stdout,*)'i1, i2', i1, i2
          write(stdout,*)'value at (i1, i2)', Hamk_bulk(i1, i2)
          write(stdout,*)'value at (i2, i1)', Hamk_bulk(i2, i1)
         !stop
        endif
     enddo
     enddo

  return
  end subroutine ham_bulk_latticegauge


  subroutine ham_bulk_LOTO(k,Hamk_bulk)
     ! This subroutine caculates Hamiltonian for
     ! bulk system without the consideration of the atom's position
     ! with the LOTO correction for phonon system
     !
     ! History
     !
     !        July/15/2017 by TianTian Zhang

     use para
     implicit none

     ! loop index
     integer :: i1,i2,ia,ib,ic,iR
     integer  :: ii,jj,mm,nn,pp,qq,xx,yy,zz

     real(dp) :: kdotr

     ! wave vector in 2d
     real(Dp) :: k(3)

     ! coordinates of R vector
     real(Dp) :: R(3), R1(3), R2(3)

     ! Hamiltonian of bulk system
     complex(Dp),intent(out) ::Hamk_bulk(Num_wann, Num_wann)

     !> see eqn. (3) in J. Phys.: Condens. Matter 22 (2010) 202201
     complex(Dp),allocatable :: nac_correction(:, :, :)

     complex(dp), allocatable :: mat1(:, :)
     complex(dp), allocatable :: mat2(:, :)
     real(dp) :: temp1(3)=(/0.0,0.0,0.0/)
     real(dp) :: temp2=0.0
     real(dp) :: temp3(30),constant_t
     real(dp) ::A_ii(3)=(/0.0,0.0,0.0/)
     real(dp) ::A_jj(3)=(/0.0,0.0,0.0/)

     !> k times Born charge
     real(dp), allocatable :: kBorn(:, :)

     real(dp) :: nac_q

     allocate(kBorn(Num_atoms, 3))
     allocate(mat1(Num_wann, Num_wann))
     allocate(mat2(Num_wann, Num_wann))
     allocate(nac_correction(Num_wann, Num_wann, Nrpts))
     mat1 = 0d0
     mat2 = 0d0
     nac_correction= 0d0

     Hamk_bulk=0d0

     !>  add loto splitting term
     temp1(1:3)= (/0.0,0.0,0.0/)
     temp2= 0.0
     if (abs((k(1)**2+k(2)**2+k(3)**2)).le.0.00001)then  !> skip k=0
        k(1)=0.00001d0
        k(2)=0.00001d0
        k(3)=0.00001d0
     endif

     !> see eqn. (3) in J. Phys.: Condens. Matter 22 (2010) 202201
     do qq= 1, 3
        temp1(qq)= k(1)*Diele_Tensor(qq, 1)+k(2)*Diele_Tensor(qq, 2)+k(3)*Diele_Tensor(qq, 3)
     enddo
     temp2= k(1)*temp1(1)+ k(2)*temp1(2)+ k(3)*temp1(3)
     constant_t= 4d0*3.1415926d0/(temp2*CellVolume)*VASPToTHZ

     do ii=1, Num_atoms
        do pp=1, 3
           kBorn(ii, pp)=  k(1)*Born_Charge(ii,1,pp)+k(2)*Born_Charge(ii,2,pp)+k(3)*Born_Charge(ii,3,pp)
        enddo
     enddo

     nac_correction= 0d0
     do iR=1, Nrpts
        R(1)=dble(irvec(1,iR))
        R(2)=dble(irvec(2,iR))
        R(3)=dble(irvec(3,iR))
        kdotr=k(1)*R(1) + k(2)*R(2) + k(3)*R(3)

        do ii= 1,Num_atoms
           do pp= 1, 3
              do jj= 1, Num_atoms
                 do qq= 1,3
                    nac_q= kBorn(jj, qq)*kBorn(ii, pp)*constant_t/sqrt(Atom_Mass(ii)*Atom_Mass(jj))
                    nac_correction((ii-1)*3+pp,(jj-1)*3+qq,iR) = HmnR((ii-1)*3+pp,(jj-1)*3+qq,iR) + dcmplx(nac_q,0)
                 enddo  ! qq
              enddo  ! jj
           enddo ! pp
        enddo  ! ii

        Hamk_bulk(:, :)= Hamk_bulk(:, :) &
        + nac_correction(:, :, iR)*Exp(2d0*pi*zi*kdotr)/ndegen(iR)
     enddo ! iR

     ! check hermitcity
     do i1=1, Num_wann
     do i2=1, Num_wann
        if(abs(Hamk_bulk(i1,i2)-conjg(Hamk_bulk(i2,i1))).ge.1e-6)then
          write(stdout,*)'there are something wrong with Hamk_bulk'
          write(stdout,*)'i1, i2', i1, i2
          write(stdout,*)'value at (i1, i2)', Hamk_bulk(i1, i2)
          write(stdout,*)'value at (i2, i1)', Hamk_bulk(i2, i1)
         !stop
        endif
     enddo
     enddo

     deallocate(kBorn)
     deallocate(mat1)
     deallocate(mat2)
     deallocate(nac_correction)
  return
  end subroutine ham_bulk_LOTO

  subroutine ham_bulk_kp(k,Hamk_bulk)
     ! This subroutine caculates Hamiltonian for
     ! bulk system without the consideration of the atom's position
     !
     ! History
     !
     !        May/29/2011 by Quansheng Wu

     use para
     implicit none

     integer :: i1,i2,ia,ib,ic,iR, nwann

     ! coordinates of R vector
     real(Dp) :: R(3), R1(3), R2(3), kdotr, kx, ky, kz, m, shift

     complex(dp) :: factor

     real(dp), intent(in) :: k(3)

     ! Hamiltonian of bulk system
     complex(Dp),intent(out) ::Hamk_bulk(Num_wann, Num_wann)
     if (Num_wann/=3) then
        print *, "Error : in this kp model, num_wann should be 3"
        print *, 'Num_wann', Num_wann

        stop
     endif

     kx=k(1)
     ky=k(2)
     kz=k(3)
     m=-0.60d0

     shift= 0.625d0
     !> kp
     Hamk_bulk(1, 1)= kx**3+ shift
     Hamk_bulk(1, 2)= ky
     Hamk_bulk(1, 3)= kz
     Hamk_bulk(2, 1)= ky
     Hamk_bulk(2, 2)= -kx+ kx**2+ m
     Hamk_bulk(2, 3)= -1d0/4d0*ky*kz
     Hamk_bulk(3, 1)= kz
     Hamk_bulk(3, 2)= -1d0/4d0*ky*kz
     Hamk_bulk(3, 3)= -kx- kx**2- m

    !Hamk_bulk(1, 1)= sin(kx)**3
    !Hamk_bulk(1, 2)= sin(ky)
    !Hamk_bulk(1, 3)= sin(kz)
    !Hamk_bulk(2, 1)= sin(ky)
    !Hamk_bulk(2, 2)= -sin(kx)+ (1-cos(kx))*2d0+ m
    !Hamk_bulk(2, 3)= -1d0/4d0*sin(ky)*sin(kz)
    !Hamk_bulk(3, 1)= sin(kz)
    !Hamk_bulk(3, 2)= -1d0/4d0*sin(ky)*sin(kz)
    !Hamk_bulk(3, 3)= -sin(kx)- (1-cos(kx))*2d0- m

     ! check hermitcity
     do i1=1, Num_wann
     do i2=1, Num_wann
        if(abs(Hamk_bulk(i1,i2)-conjg(Hamk_bulk(i2,i1))).ge.1e-6)then
          write(stdout,*)'there are something wrong with Hamk_bulk'
          write(stdout,*)'i1, i2', i1, i2
          write(stdout,*)'value at (i1, i2)', Hamk_bulk(i1, i2)
          write(stdout,*)'value at (i2, i1)', Hamk_bulk(i2, i1)
         !stop
        endif
     enddo
     enddo

     return
  end subroutine ham_bulk_kp


