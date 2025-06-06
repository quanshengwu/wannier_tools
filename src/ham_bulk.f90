! This subroutine is used to caculate Hamiltonian for
! bulk system .

! History
!        May/29/2011 by Quansheng Wu

subroutine ham_bulk_atomicgauge(k,Hamk_bulk)
   ! This subroutine caculates Hamiltonian for
   ! bulk system with the consideration of the atom's position
   !
   ! History
   !
   !        May/29/2011 by Quansheng Wu
   !  Atomic gauge Guan Yifei 2019
   !  Lattice gauge Hl
   !  Atomic gauge Ha= U* Hl U 
   !  where U = e^ik.wc(i) on diagonal

   use para
   implicit none

   integer :: i1,i2,iR

   ! wave vector in 3d
   real(Dp) :: k(3), kdotr, pos0(3), dis

   complex(dp) :: factor

   real(dp) :: pos(3), pos1(3), pos2(3), pos_cart(3), pos_direct(3)
   ! Hamiltonian of bulk system
   complex(Dp),intent(out) ::Hamk_bulk(Num_wann, Num_wann)
   complex(dp), allocatable :: mat1(:, :)
   real(dp), external :: norm

   allocate(mat1(Num_wann, Num_wann))

   Hamk_bulk=0d0
  !do iR=1, Nrpts
  !   kdotr= k(1)*irvec(1,iR) + k(2)*irvec(2,iR) + k(3)*irvec(3,iR)
  !   factor= exp(pi2zi*kdotr)

  !   Hamk_bulk(:, :)= Hamk_bulk(:, :) &
  !      + HmnR(:, :, iR)*factor/ndegen(iR)
  !enddo ! iR
 
  !mat1=0d0
  !do i1=1,Num_wann
  !   pos0=Origin_cell%wannier_centers_direct(:, i1)
  !   kdotr= k(1)*pos0(1)+ k(2)*pos0(2)+ k(3)*pos0(3)
  !   mat1(i1,i1)= exp(pi2zi*kdotr)
  !enddo
  !Hamk_bulk=matmul(conjg(mat1),matmul(Hamk_bulk,mat1))


   !> the first atom in home unit cell
   do iR=1, Nrpts
      do i2=1, Num_wann
         pos2= Origin_cell%wannier_centers_direct(:, i2)
         !> the second atom in unit cell R
         do i1=1, Num_wann
            pos1= Origin_cell%wannier_centers_direct(:, i1)
            pos_direct= irvec(:, iR)
            pos_direct= pos_direct+ pos2- pos1

            call direct_cart_real(pos_direct, pos_cart, Origin_cell%lattice)

            dis= norm(pos_cart)
            if (dis> Rcut) cycle

            kdotr=k(1)*pos_direct(1) + k(2)*pos_direct(2) + k(3)*pos_direct(3)
            factor= (cos(twopi*kdotr)+zi*sin(twopi*kdotr))

            Hamk_bulk(i1, i2)= Hamk_bulk(i1, i2) &
               + HmnR(i1, i2, iR)*factor/ndegen(iR)
         enddo ! i1
      enddo ! i2
   enddo ! iR

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
end subroutine ham_bulk_atomicgauge


subroutine valley_k_atomicgauge(k,valley_k)
   ! This subroutine performs the Fourier transform of avalley operator
   ! History
   !        Nov/5/2023 by Quansheng Wu

   use para
   implicit none

   integer :: i1,i2,iR

   ! wave vector in 3d
   real(Dp) :: k(3), kdotr, pos0(3), dis

   complex(dp) :: factor

   real(dp) :: pos(3), pos1(3), pos2(3), pos_cart(3), pos_direct(3)
   ! Hamiltonian of bulk system
   complex(Dp),intent(out) :: valley_k(Num_wann, Num_wann)
   real(dp), external :: norm
   

   valley_k= 0d0
   !> the first atom in home unit cell
   do iR=1, Nrpts_valley
      do i2=1, Num_wann
         pos2= Origin_cell%wannier_centers_direct(:, i2)
         !> the second atom in unit cell R
         do i1=1, Num_wann
            pos1= Origin_cell%wannier_centers_direct(:, i1)
            pos_direct= irvec_valley(:, iR)
            pos_direct= pos_direct+ pos2- pos1

            call direct_cart_real(pos_direct, pos_cart, Origin_cell%lattice)

            dis= norm(pos_cart)
            if (dis> Rcut) cycle

            kdotr=k(1)*pos_direct(1) + k(2)*pos_direct(2) + k(3)*pos_direct(3)
            factor= (cos(twopi*kdotr)+zi*sin(twopi*kdotr))

            valley_k(i1, i2)= valley_k(i1, i2) &
               + valley_operator_R(i1, i2, iR)*factor
         enddo ! i1
      enddo ! i2
   enddo ! iR

   ! check hermitcity
   do i1=1, Num_wann
      do i2=1, Num_wann
         if(abs(valley_k(i1,i2)-conjg(valley_k(i2,i1))).ge.1e-6)then
            write(stdout,*)'there is something wrong with Hamk_bulk'
            write(stdout,*)'i1, i2', i1, i2
            write(stdout,*)'value at (i1, i2)', valley_k(i1, i2)
            write(stdout,*)'value at (i2, i1)', valley_k(i2, i1)
            !stop
         endif
      enddo
   enddo


   return
end subroutine valley_k_atomicgauge

subroutine d2Hdk2_atomicgauge(k, DHDk2_wann)
   !> second derivatve of H(k)
   use para, only : Nrpts, irvec, crvec, Origin_cell, HmnR, ndegen, &
       Num_wann, dp, Rcut, pi2zi, zi, twopi, pi
   implicit none

   !> momentum in 3D BZ
   real(dp), intent(in) :: k(3)

   !> second derivate of H(k)
   complex(dp), intent(out) :: DHDk2_wann(Num_wann, Num_wann, 3, 3)

   integer :: iR, i1, i2, i, j

   real(dp) :: pos(3), pos1(3), pos2(3), pos_cart(3), pos_direct(3)
   real(dp) :: kdotr, dis
   complex(dp) :: ratio
   real(dp), external :: norm

   DHDk2_wann= 0d0
   !> the first atom in home unit cell
   do iR=1, Nrpts
      do i2=1, Num_wann
         pos2= Origin_cell%wannier_centers_direct(:, i2)
         !> the second atom in unit cell R
         do i1=1, Num_wann
            pos1= Origin_cell%wannier_centers_direct(:, i1)
            pos_direct= irvec(:, iR)
            pos_direct= pos_direct+ pos2- pos1

            call direct_cart_real(pos_direct, pos_cart, Origin_cell%lattice)

            dis= norm(pos_cart)
            if (dis> Rcut) cycle

            kdotr=k(1)*pos_direct(1) + k(2)*pos_direct(2) + k(3)*pos_direct(3)
            ratio= (cos(twopi*kdotr)+zi*sin(twopi*kdotr))/ndegen(iR)

            do j=1, 3
               do i=1, 3
                  DHDk2_wann(i1, i2, i, j)=DHDk2_wann(i1, i2, i, j) &
                     -pos_cart(i)*pos_cart(j)*HmnR(i1, i2, iR)*ratio
               enddo ! j 
            enddo ! i
         enddo ! i1
      enddo ! i2
   enddo ! iR

   return
end subroutine d2Hdk2_atomicgauge

subroutine d2Hdk2_atomicgauge_wann(k, D2HDk2_wann)
   !> second derivatve of H(k)
   use para, only : Nrpts, irvec, crvec, Origin_cell, HmnR, ndegen, &
       Num_wann, dp, Rcut, pi2zi, zi, twopi, pi
   implicit none

   !> momentum in 3D BZ
   real(dp), intent(in) :: k(3)

   !> second derivate of H(k)
   complex(dp), intent(out) :: D2HDk2_wann(Num_wann, Num_wann, 3, 3)

   integer :: iR, i1, i2, i, j

   real(dp) :: pos(3), pos1(3), pos2(3), pos_cart(3), pos_direct(3)
   real(dp) :: kdotr, dis
   complex(dp) :: ratio
   real(dp), external :: norm

   D2HDk2_wann= 0d0
   !> the first atom in home unit cell
   do iR=1, Nrpts
      do i2=1, Num_wann
         pos2= Origin_cell%wannier_centers_direct(:, i2)
         !> the second atom in unit cell R
         do i1=1, Num_wann
            pos1= Origin_cell%wannier_centers_direct(:, i1)
            pos_direct= irvec(:, iR)
            pos_direct= pos_direct+ pos2- pos1

            call direct_cart_real(pos_direct, pos_cart, Origin_cell%lattice)

            dis= norm(pos_cart)
            if (dis> Rcut) cycle

            kdotr=k(1)*pos_direct(1) + k(2)*pos_direct(2) + k(3)*pos_direct(3)
            ratio= (cos(twopi*kdotr)+zi*sin(twopi*kdotr))/ndegen(iR)

            do j=1, 3
               do i=1, 3
                  D2HDk2_wann(i1, i2, i, j)=D2HDk2_wann(i1, i2, i, j) &
                     -pos_cart(i)*pos_cart(j)*HmnR(i1, i2, iR)*ratio
               enddo ! j 
            enddo ! i
         enddo ! i1
      enddo ! i2
   enddo ! iR

   return
end subroutine d2Hdk2_atomicgauge_wann

subroutine d2Hdk2_atomicgauge_Ham(UU, dHdkdk, D2HDk2_Ham)
   !> Transform d2H/dk^2 from Wannier gauge to Hamiltonian gauge
   use para, only : Num_wann, dp
   implicit none

   ! Inputs
   complex(dp), intent(in)  :: UU(Num_wann, Num_wann)
   complex(dp), intent(in)  :: dHdkdk(Num_wann, Num_wann, 3, 3)

   ! Output
   complex(dp), intent(out) :: D2HDk2_Ham(Num_wann, Num_wann, 3, 3)

   ! Locals
   integer :: i, j
   complex(dp) :: mat_tmp(Num_wann, Num_wann)

   do j = 1, 3
      do i = 1, 3
         call rotation_to_Ham_basis(UU, dHdkdk(:, :, i, j), mat_tmp)
         D2HDk2_Ham(:, :, i, j) = mat_tmp
      enddo ! i
   enddo ! j

   return
end subroutine d2Hdk2_atomicgauge_Ham

subroutine dHdk_atomicgauge(k, velocity_Wannier)
   !> Velocity operator in Wannier basis using atomic gauge
   !> First derivate of H(k); dH(k)/dk
   use para, only : Nrpts, irvec, Origin_cell, HmnR, ndegen, &
       Num_wann, dp, Rcut, pi2zi,  &
       zi, soc, zzero, twopi
   implicit none

   !> momentum in 3D BZ
   real(dp), intent(in) :: k(3)

   !> velocity operator in Wannier basis using atomic gauge 
   complex(dp), intent(out) :: velocity_Wannier(Num_wann, Num_wann, 3)

   integer :: iR, i1, i2, i

   real(dp) :: pos1(3), pos2(3), pos_cart(3), pos_direct(3)
   real(dp) :: kdotr, dis
   complex(dp) :: ratio
   real(dp), external :: norm

   velocity_Wannier= zzero
   do iR=1, Nrpts
      do i2=1, Num_wann
         pos2= Origin_cell%wannier_centers_direct(:, i2)
         !> the second atom in unit cell R
         do i1=1, Num_wann
            !> the first atom in home unit cell
            pos1= Origin_cell%wannier_centers_direct(:, i1)
            pos_direct= irvec(:, iR)+ pos2- pos1

            call direct_cart_real(pos_direct, pos_cart, Origin_cell%lattice)

            dis= norm(pos_cart)
            if (dis> Rcut) cycle

            kdotr=k(1)*pos_direct(1) + k(2)*pos_direct(2) + k(3)*pos_direct(3)
            ratio= (cos(twopi*kdotr)+zi*sin(twopi*kdotr))/ndegen(iR)

            do i=1, 3
               velocity_Wannier(i1, i2, i)=velocity_Wannier(i1, i2, i)+ &
                  zi*pos_cart(i)*HmnR(i1, i2, iR)*ratio
            enddo ! i

         enddo ! i2
      enddo ! i1
   enddo ! iR

   return
end subroutine dHdk_atomicgauge

subroutine dHdk_atomicgauge_Ham(k, eigvec, Vmn_Ham)
   !> Velocity operator in Hamiltonian basis using atomic gauge
   !> see https://www.wanniertools.org/theory/tight-binding-model/
   use para, only : Num_wann, dp
   implicit none

   !> momentum in 3D BZ
   real(dp), intent(in) :: k(3)

   !> eigenvectors of H, H*eigvec(:, n)= E(n)*eigvec(:, n)
   complex(dp), intent(in) :: eigvec(Num_wann, Num_wann)

   !> velocity operator in the diagonalized Hamiltonian basis using lattice gauge
   complex(dp), intent(out) :: Vmn_Ham(Num_wann, Num_wann, 3)

   !> velocity operator in Wannier basis using lattice gauge
   complex(dp), allocatable :: Vmn_wann(:, :, :)

   integer :: i

   allocate(Vmn_wann(Num_wann, Num_wann, 3))
   Vmn_Ham= 0d0
   call dHdk_atomicgauge(k, Vmn_wann)
   do i=1, 3
      call rotation_to_Ham_basis(eigvec, Vmn_wann(:, :, i), Vmn_Ham(:, :, i))
   enddo
   deallocate(Vmn_wann)

   return
end subroutine dHdk_atomicgauge_Ham

subroutine ham_bulk_latticegauge(k,Hamk_bulk)
   ! This subroutine caculates Hamiltonian for
   ! bulk system without the consideration of the atom's position
   !
   ! History
   !
   !        May/29/2011 by Quansheng Wu

   use para, only : dp, pi2zi, HmnR, ndegen, nrpts, irvec, Num_wann, stdout, twopi, zi
   implicit none

   ! loop index
   integer :: i1,i2,iR

   real(dp) :: kdotr

   complex(dp) :: factor

   real(dp), intent(in) :: k(3)

   ! Hamiltonian of bulk system
   complex(Dp),intent(out) ::Hamk_bulk(Num_wann, Num_wann)

   Hamk_bulk=0d0
   do iR=1, Nrpts
      kdotr= k(1)*irvec(1,iR) + k(2)*irvec(2,iR) + k(3)*irvec(3,iR)
      factor= (cos(twopi*kdotr)+zi*sin(twopi*kdotr))

      Hamk_bulk(:, :)= Hamk_bulk(:, :)+ HmnR(:, :, iR)*factor/ndegen(iR)
   enddo ! iR
   
   !call mat_mul(Num_wann, mirror_z, Hamk_bulk, mat1)
   !call mat_mul(Num_wann, mat1, mirror_z, mat2)
   !Hamk_bulk= (Hamk_bulk+ mat2)/2d0

   ! check hermitcity
  !打印出k点
   ! write(*,*)'check Hamk_bulk hermitcity'
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
end subroutine ham_bulk_latticegauge


subroutine S_bulk_latticegauge(k,Sk_bulk)
   ! This subroutine caculates Hamiltonian for
   ! bulk system without the consideration of the atom's position
   !
   ! History
   !
   !        May/29/2011 by Quansheng Wu

   use para, only : dp, pi2zi, SmnR, ndegen, nrpts, irvec, Num_wann, stdout, twopi, zi
   implicit none

   ! loop index
   integer :: i1,i2,iR

   real(dp) :: kdotr

   complex(dp) :: factor

   real(dp), intent(in) :: k(3)

   ! Hamiltonian of bulk system
   complex(Dp),intent(out) ::Sk_bulk(Num_wann, Num_wann)

   Sk_bulk=0d0
   do iR=1, Nrpts
      kdotr= k(1)*irvec(1,iR) + k(2)*irvec(2,iR) + k(3)*irvec(3,iR)
      factor= (cos(twopi*kdotr)+zi*sin(twopi*kdotr))

      Sk_bulk(:, :)= Sk_bulk(:, :)+ SmnR(:, :, iR)*factor
   enddo ! iR
   
   !call mat_mul(Num_wann, mirror_z, Hamk_bulk, mat1)
   !call mat_mul(Num_wann, mat1, mirror_z, mat2)
   !Hamk_bulk= (Hamk_bulk+ mat2)/2d0

   ! check hermitcity
   ! write(*,*)'check Sk_bulk'
  do i1=1, Num_wann
    do i2=1, Num_wann
       if(abs(Sk_bulk(i1,i2)-conjg(Sk_bulk(i2,i1))).ge.1e-6)then
          write(stdout,*)'there is something wrong with Sk_bulk'
          write(stdout,*)'i1, i2', i1, i2
          write(stdout,*)'value at (i1, i2)', Sk_bulk(i1, i2)
          write(stdout,*)'value at (i2, i1)', Sk_bulk(i2, i1)
          !stop
       endif
    enddo
  enddo

   return
end subroutine S_bulk_latticegauge

!---------------------------------------------------------------
! Surroutines:change the non-orthogonal basis H and S to orthogonal basis
!   H_new = U * H_old * U
!   U = S_vec * M_inv * S_vec^H
!---------------------------------------------------------------
subroutine orthogonalize_hamiltonian(Hamk_bulk, Sk_bulk , N)
   use para, only: Dp, stdout,E_fermi,eV2Hartree
   implicit none
 
   integer, intent(in)       :: N
   complex(Dp), intent(in )    :: Sk_bulk(N,N)
   complex(Dp), intent(inout)  :: Hamk_bulk(N,N)
 
   !-- 局部 --
   integer           :: i, info, j
   real(Dp), allocatable    :: S_vals(:)
   complex(Dp), allocatable :: S_vec(:,:), M_inv(:,:), U(:,:)
 

   allocate(S_vals(N))
   allocate(S_vec (N,N))
   allocate(M_inv (N,N))
   allocate(U     (N,N))
   U = (0.0_dp,0.0_dp)
   ! 1) Sk_bulk =  S_vec * diag(S_vals) * S_vec^H 
   S_vec = Sk_bulk
   call eigensystem_c('V','U', N, S_vec, S_vals)
 
   ! 2)  M_inv = diag(1/√S_vals)
   M_inv = (0.0_dp,0.0_dp)
   ! do i = 1, N
   !   M_inv(i,i) = 1.0_dp / sqrt(S_vals(i))
   ! end do

   do i = 1, N
      if (S_vals(i) <= 0.0_dp) then
         write(stdout,'(a,i5,1pe13.4)') ' WARNING: λ <= 0, i =', i, S_vals(i)
      else
         M_inv(i,i) = 1.0_dp / sqrt(S_vals(i))
      endif
   enddo
 
   ! 3) U = S_vec * M_inv * S_vec^H
   ! U = matmul( S_vec, matmul(M_inv, conjg(transpose(S_vec))) )
   U = matmul( S_vec, M_inv )
   U = matmul(U, conjg(transpose(S_vec)))

   ! call mat_mul(N, S_vec, M_inv, U)
   ! call mat_mul(N, U, conjg(transpose(S_vec)), U)

 
   ! 4)Hamiltonian = U * Hamk_bulk * U
   ! Hamk_bulk = matmul( U, matmul(Hamk_bulk, U) )
   Hamk_bulk = matmul(U, Hamk_bulk)
   Hamk_bulk = matmul(Hamk_bulk, U)
   ! call mat_mul(N, U, Hamk_bulk, Hamk_bulk)
   ! call mat_mul(N, Hamk_bulk, U, Hamk_bulk)

   do i = 1, N
      Hamk_bulk(i,i) = Hamk_bulk(i,i) - E_fermi*eV2Hartree
   enddo

   deallocate(S_vals, S_vec, M_inv, U)
end subroutine orthogonalize_hamiltonian
 


subroutine dHdk_latticegauge_wann(k, velocity_Wannier)
   use para, only : Nrpts, irvec, crvec, Origin_cell, &
      HmnR, ndegen, Num_wann, zi, pi2zi, dp, twopi
   implicit none

   !> momentum in 3D BZ
   real(dp), intent(in) :: k(3)

   !> velocity operator in Wannier basis using lattice gauge
   complex(dp), intent(out) :: velocity_Wannier(Num_wann, Num_wann, 3)

   integer :: iR, i

   real(dp) :: kdotr
   complex(dp) :: ratio

   do i=1, 3
      do iR= 1, Nrpts
         kdotr= k(1)*irvec(1,iR) + k(2)*irvec(2,iR) + k(3)*irvec(3,iR)
         ratio= (cos(twopi*kdotr)+zi*sin(twopi*kdotr))
         velocity_Wannier(:, :, i)= velocity_Wannier(:, :, i)+ &
            zi*crvec(i, iR)*HmnR(:,:,iR)*ratio/ndegen(iR)
      enddo ! iR
   enddo

   return
end subroutine dHdk_latticegauge_wann

subroutine dHdk_latticegauge_Ham(k, eigval, eigvec, Vmn_Ham)
   use para, only : Nrpts, irvec, crvec, Origin_cell, &
      HmnR, ndegen, Num_wann, zi, pi2zi, dp, zzero, twopi
   implicit none

   !> momentum in 3D BZ
   real(dp), intent(in) :: k(3)

   real(dp), intent(in) :: eigval(Num_wann)
   complex(dp), intent(in) :: eigvec(Num_wann, Num_wann)

   !> velocity operator in the diagonalized Hamiltonian basis using lattice gauge
   complex(dp), intent(out) :: Vmn_Ham(Num_wann, Num_wann, 3)

   !> velocity operator in Wannier basis using lattice gauge
   complex(dp), allocatable :: Vmn_wann(:, :, :), Vmn_Ham0(:, :, :)

   !> wnm=En-Em
   complex(dp), allocatable :: wnm(:, :), temp(:, :)

   !> it's a diagonal matrix for each axis
   complex(dp), allocatable :: itau(:, :, :)

   integer :: i, m, n, l

   allocate(Vmn_wann(Num_wann, Num_wann, 3))
   allocate(Vmn_Ham0(Num_wann, Num_wann, 3))
   Vmn_wann= 0d0; Vmn_Ham0= 0d0

   call dHdk_latticegauge_wann(k, Vmn_wann)
   do i=1, 3
      call rotation_to_Ham_basis(eigvec, Vmn_wann(:, :, i), Vmn_Ham0(:, :, i))
   enddo
 
   allocate(wnm(Num_wann, Num_wann))
   allocate(temp(Num_wann, Num_wann))
   allocate(itau(Num_wann, Num_wann, 3))
   wnm=0d0; temp=0d0; itau= zzero

   !\  -i\tau
   do i=1, 3
      do m=1, Num_wann
         itau(m, m, i)= -zi*Origin_cell%wannier_centers_cart(i, m)
      enddo
   enddo

   do m=1, Num_wann
      do n=1, Num_wann
         wnm(m, n)= eigval(n)-eigval(m)
      enddo
   enddo
  
   !> \sum_l (-i*\tau_{l})*conjg(vec(l, m))*vec(l, n)
   do i=1, 3
      call mat_mul(Num_wann, itau(:, :, i), eigvec, temp)
      call mat_mul(Num_wann, conjg(transpose(eigvec)), temp, Vmn_Ham(:, :, i))
   enddo
 
   do i=1, 3
      do n=1, Num_wann
         do m=1, Num_wann
            temp(m, n) = wnm(m, n)*Vmn_Ham(m, n, i)
         enddo
      enddo
      Vmn_Ham(:, :, i)= temp(:, :)
   enddo


   Vmn_Ham= Vmn_Ham+ Vmn_Ham0
  !Vmn_Ham = Vmn_Ham0

!  print *, Vmn_Ham(1, 1, 1)

  !Vmn_Ham= 0d0
  !do m=1, Num_wann
  !   do n=1, Num_wann
  !      do i=1, 3
  !         do l=1, Num_wann
  !         Vmn_Ham(m, n, i)= Vmn_Ham(m, n, i)- &
  !            (eigval(n)- eigval(m))*zi*Origin_cell%wannier_centers_cart(i, l)*conjg(eigvec(l, m))*eigvec(l, n)
  !         enddo ! summation over l
  !      enddo
  !   enddo
  !enddo
  !Vmn_Ham= Vmn_Ham+ Vmn_Ham0

!  print *, Vmn_Ham(1, 1, 1)
!  stop

   return
end subroutine dHdk_latticegauge_Ham


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
   complex(dp) :: ratio

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

   allocate(kBorn(Origin_cell%Num_atoms, 3))
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
   if (abs((k(1)**2+k(2)**2+k(3)**2)).le.0.0001)then  !> skip k=0
      k(1)=0.0001d0
      k(2)=0.0001d0
      k(3)=0.0001d0
   endif

   !> see eqn. (3) in J. Phys.: Condens. Matter 22 (2010) 202201
   do qq= 1, 3
      temp1(qq)= k(1)*Diele_Tensor(qq, 1)+k(2)*Diele_Tensor(qq, 2)+k(3)*Diele_Tensor(qq, 3)
   enddo
   temp2= k(1)*temp1(1)+ k(2)*temp1(2)+ k(3)*temp1(3)
   constant_t= 4d0*3.1415926d0/(temp2*Origin_cell%CellVolume)*VASPToTHZ

   do ii=1, Origin_cell%Num_atoms
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

      do ii= 1,Origin_cell%Num_atoms
         do pp= 1, 3
            do jj= 1, Origin_cell%Num_atoms
               do qq= 1,3
                  nac_q= kBorn(jj, qq)*kBorn(ii, pp)*constant_t/sqrt(Atom_Mass(ii)*Atom_Mass(jj))
                  nac_correction((ii-1)*3+pp,(jj-1)*3+qq,iR) = HmnR((ii-1)*3+pp,(jj-1)*3+qq,iR) + dcmplx(nac_q,0)
               enddo  ! qq
            enddo  ! jj
         enddo ! pp
      enddo  ! ii

      ratio= (cos(twopi*kdotr)+zi*sin(twopi*kdotr))
      Hamk_bulk(:, :)= Hamk_bulk(:, :) &
         + nac_correction(:, :, iR)*ratio/ndegen(iR)
   enddo ! iR

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

   deallocate(kBorn)
   deallocate(mat1)
   deallocate(mat2)
   deallocate(nac_correction)
   return
end subroutine ham_bulk_LOTO

subroutine ham_bulk_kp_abcb_graphene(k,Hamk_bulk)
   ! > construct the kp model at K point of ABCB tetralayer graphene

   use para, only : Num_wann, dp, stdout, zi, zzero, One_complex, &
      Symmetrical_Electric_field_in_eVpA, eV2Hartree, Angstrom2atomic, &
      Electric_field_in_eVpA, Origin_cell, cpuid
   implicit none

   integer :: i1,i2,ia,ib,ic,iR, nwann, io, i
   integer :: add_electric_field
   real(dp) :: pos(Origin_cell%Num_atoms)
   real(dp) :: static_potential

   ! coordinates of R vector
   real(Dp) :: R(3), R1(3), R2(3), kdotr, kx, ky, kz
   real(dp) :: v,v3,v4,v6,v7,d1,d2,d3,d4,d5,d6,d7,d8,g1,g2,g3,g4,g5,g6

   complex(dp) :: factor, km, kp

   real(dp), intent(in) :: k(3)

   !> the k point that the kp model is constructed
   !> in fractional coordinate
   real(dp) :: k_center_direct(3), l(19)

   ! Hamiltonian of bulk system
   complex(Dp),intent(out) ::Hamk_bulk(Num_wann, Num_wann)

   if (Num_wann/=8) then
      print *, "Error : in this kp model, num_wann should be 8"
      print *, 'Num_wann', Num_wann
      stop
   endif
  !l=(/ 1.26416117e+01, -5.57664351e+00, -5.87940524e+00,  1.60839693e+01, &
  ! 8.13551651e+00,  2.82753440e-01,  9.16947940e-03, -2.61922210e-03, &
  !-2.77215942e-03,  7.01879554e-03, -3.38271168e-01, -6.64427749e-02, &
  !-3.08356540e-02, -8.82150889e-03, -2.57865658e-03,  1.83192652e-01, &
  ! 4.35695252e-03, -2.28829893e-01, -1.82821963e-02/)

   l=(/1.57937270e+01, -2.49047067e+00, -1.82376120e+00,  1.21854647e-02, &
    -5.75783391e+00, -2.09787893e-01,  7.48348077e-03,  2.03490675e-01, &
    -1.82294202e-04,  3.11506594e-02, -1.32325756e-01,  2.72568864e-01, &
     2.27263042e-02, -2.57806046e-01, -2.01920757e-01, -2.46168084e-01, &
    -1.45799538e-01, -1.59225389e-01, -2.20335989e-01/)

   v=l(1); v3=l(2); v4=l(3); v6=l(4); v7=l(5);
   d1=l(6); d2= l(7); d3= l(8); d4=l(9); d5=l(10); d6= l(11); d7=l(12); d8= l(19)
   g1= l(13); g2= l(14); g3= l(15); g4= l(16); g5= l(17); g6= l(18);

   k_center_direct= (/1d0/3d0, 1d0/3d0, 0d0/)
   kx=k(1) -k_center_direct(1)
   ky=k(2) -k_center_direct(2)
   kz=k(3) -k_center_direct(3)

   km=-kx+ zi*ky
   kp=-kx- zi*ky
   !> write out the parameters
   if (cpuid.eq.0) then
   !  write(stdout, '(a, 8f8.4)') '> Onsite energies: ', d1, d2, d3, d4, d5, d6, d7, d8
   endif

   Hamk_bulk= 0d0
   !> kp
   Hamk_bulk(:, 1)= (/d1*One_complex, v*kp, -v4*kp, v3*km, v6*km, g6*One_complex, zzero, zzero/)
   Hamk_bulk(:, 2)= (/v*km, d2*One_complex, g1*One_complex, -v4*kp, v7*kp, v6*km, zzero, zzero/)
   Hamk_bulk(:, 3)= (/-v4*km, g1*One_complex, d3*One_complex, v*kp, -v4*kp, v3*km, g4*One_complex, zzero/)
   Hamk_bulk(:, 4)= (/v3*kp, -v4*km, v*km, d4*One_complex, g2*One_complex, -v4*kp, zzero, g5*One_complex/)
   Hamk_bulk(:, 5)= (/v6*kp, v7*km, -v4*km, g2*One_complex, d5*One_complex, v*kp, -v4*km, g3*One_complex/)
   Hamk_bulk(:, 6)= (/g6*One_complex, v6*kp, v3*kp, -v4*km, v*km, d6*One_complex, v3*kp, -v4*km/)
   Hamk_bulk(:, 7)= (/zzero, zzero, g4*One_complex, zzero, -v4*kp, v3*km, d7*One_complex, v*kp/)
   Hamk_bulk(:, 8)= (/zzero, zzero, zzero, g5*One_complex, g3*One_complex, -v4*kp, v*km, d8*One_complex/)

   !> add electrical field

   call get_stacking_direction_and_pos(add_electric_field, pos)
 
   if (add_electric_field>0) then
      io=0
      do ia=1, Origin_cell%Num_atoms
         !static_potential= pos(ia)*Origin_cell%cell_parameters(add_electric_field)*Electric_field_in_eVpA
         !if (Inner_symmetrical_Electric_Field) then
         !   static_potential= abs(pos(ia)*Origin_cell%cell_parameters(add_electric_field))*Electric_field_in_eVpA
         !endif
         static_potential= abs(pos(ia)*Origin_cell%cell_parameters(add_electric_field))&
            *Symmetrical_Electric_field_in_eVpA*eV2Hartree/Angstrom2atomic+ &
            (pos(ia)*Origin_cell%cell_parameters(add_electric_field))*&
            Electric_field_in_eVpA*eV2Hartree/Angstrom2atomic

         do i=1, Origin_cell%nprojs(ia)
            io=io+1
            Hamk_bulk(io, io)= Hamk_bulk(io, io)+ static_potential
         enddo ! nproj
      enddo ! ia
   endif  ! add electric field or not

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
   Hamk_bulk= Hamk_bulk*eV2Hartree

   return
end subroutine ham_bulk_kp_abcb_graphene


subroutine ham_bulk_kp(k,Hamk_bulk)
   ! > construct the kp model at K point of WC system
   ! > space group is 187. The little group is C3h
   ! Sep/10/2018 by Quansheng Wu

   use para, only : Num_wann, dp, stdout, zi
   implicit none

   integer :: i1,i2,ia,ib,ic,iR, nwann

   ! coordinates of R vector
   real(Dp) :: R(3), R1(3), R2(3), kdotr, kx, ky, kz
   real(dp) :: A1, A2, B1, B2, C1, C2, D1, D2
   real(dp) :: m1x, m2x, m3x, m4x, m1z, m2z, m3z, m4z
   real(dp) :: E1, E2, E3, E4

   complex(dp) :: factor, kminus, kplus

   real(dp), intent(in) :: k(3)

   ! Hamiltonian of bulk system
   complex(Dp),intent(out) ::Hamk_bulk(Num_wann, Num_wann)

   if (Num_wann/=4) then
      print *, "Error : in this kp model, num_wann should be 4"
      print *, 'Num_wann', Num_wann
      stop
   endif

   kx=k(1)
   ky=k(2)
   kz=k(3)
   E1= 1.25d0
   E2= 0.85d0
   E3= 0.25d0
   E4=-0.05d0

   A1= 0.10d0
   A2= 0.30d0
   B1= 0.05d0
   B2= 0.1d0

   C1= -1.211d0
   C2= 1.5d0
   D1=-0.6d0
   D2= 0.6d0

   m1x= -1.8d0
   m2x= -1.6d0
   m3x=  1.2d0
   m4x=  1.4d0
   m1z= 2d0
   m2z= 3d0
   m3z=-1d0
   m4z=-1d0

   kminus= kx- zi*ky
   kplus= kx+ zi*ky

   Hamk_bulk= 0d0
   !> kp
   Hamk_bulk(1, 1)= E1+ m1x*(kx*kx+ky*ky)+ m1z*kz*kz
   Hamk_bulk(2, 2)= E2+ m2x*(kx*kx+ky*ky)+ m2z*kz*kz
   Hamk_bulk(3, 3)= E3+ m3x*(kx*kx+ky*ky)+ m3z*kz*kz
   Hamk_bulk(4, 4)= E4+ m4x*(kx*kx+ky*ky)+ m4z*kz*kz

   Hamk_bulk(1, 2)=-zi*D1*kplus*kz
   Hamk_bulk(2, 1)= zi*D1*kminus*kz
   Hamk_bulk(3, 4)= zi*D2*kminus*kz
   Hamk_bulk(4, 3)=-zi*D2*kplus*kz

   Hamk_bulk(1, 4)= zi*C1*kz
   Hamk_bulk(2, 3)= zi*C2*kz
   Hamk_bulk(3, 2)=-zi*C2*kz
   Hamk_bulk(4, 1)=-zi*C1*kz

   Hamk_bulk(1, 3)=  A1*kplus+ B1*kminus*kminus !+ D*kplus*kplus*kplus
   Hamk_bulk(2, 4)=  A2*kminus+ B2*kplus*kplus !+ D*kminus*kminus*kminus
   Hamk_bulk(3, 1)= conjg(Hamk_bulk(1, 3))
   Hamk_bulk(4, 2)= conjg(Hamk_bulk(2, 4))

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
end subroutine ham_bulk_kp

subroutine ham_bulk_coo_densehr(k,nnzmax, nnz, acoo,icoo,jcoo)
   !> This subroutine use sparse hr format
   ! History
   !        Dec/10/2018 by Quansheng Wu
   use para
   implicit none

   real(dp), intent(in) :: k(3)
   integer, intent(in) :: nnzmax
   integer, intent(out) :: nnz
   integer,intent(inout):: icoo(nnzmax),jcoo(nnzmax)
   complex(dp),intent(inout) :: acoo(nnzmax)

   ! loop index
   integer :: i1,i2,ia,ib,ic,iR, iz
   integer :: nwann

   real(dp) :: kdotr

   ! wave vector in 3D BZ
   ! coordinates of R vector
   real(Dp) :: R(3), R1(3), R2(3)

   complex(dp) :: factor

   ! Hamiltonian of bulk system
   complex(Dp), allocatable :: Hamk_bulk(:, :)

   allocate( Hamk_bulk(Num_wann, Num_wann))
   Hamk_bulk= 0d0

   do iR=1, Nrpts
      ia=irvec(1,iR)
      ib=irvec(2,iR)
      ic=irvec(3,iR)

      R(1)=dble(ia)
      R(2)=dble(ib)
      R(3)=dble(ic)
      kdotr=k(1)*R (1) + k(2)*R (2) + k(3)*R (3)
      factor= (cos(twopi*kdotr)+zi*sin(twopi*kdotr))

      Hamk_bulk(:, :)=&
         Hamk_bulk(:, :) &
         + HmnR(:, :, iR)*factor/ndegen(iR)
   enddo ! iR

   !> transform Hamk_bulk into sparse COO format

   nnz= 0
   do i1=1, Num_wann
      do i2=1, Num_wann
         if(abs(Hamk_bulk(i1,i2)).ge.1e-6)then
            nnz= nnz+ 1
            if (nnz>nnzmax) stop 'nnz is larger than nnzmax in ham_bulk_coo_densehr'
            icoo(nnz) = i1
            jcoo(nnz) = i2
            acoo(nnz) = Hamk_bulk(i1, i2)
         endif
      enddo
   enddo

   ! check hermitcity
   do i1=1, Num_wann
      do i2=1, Num_wann
         if(abs(Hamk_bulk(i1,i2)-conjg(Hamk_bulk(i2,i1))).ge.1e-6)then
            write(stdout,*)'there is something wrong with Hamk_bulk'
            write(stdout,*)'i1, i2', i1, i2
            write(stdout,*)'value at (i1, i2)', Hamk_bulk(i1, i2)
            write(stdout,*)'value at (i2, i1)', Hamk_bulk(i2, i1)
         endif
      enddo
   enddo

   return
end subroutine ham_bulk_coo_densehr


subroutine ham_bulk_coo_sparsehr_latticegauge(k,acoo,icoo,jcoo)
   !> This subroutine use sparse hr format
   !> Here we use atomic gauge which means the atomic position is taken into account
   !> in the Fourier transformation
   use para
   implicit none

   real(dp) :: k(3), posij(3)
   real(dp) :: kdotr
   integer,intent(inout):: icoo(splen),jcoo(splen)!,splen
   complex(dp),intent(inout) :: acoo(splen)
   complex(dp) ::  ratio
   integer :: i,j,ir

   do i=1,splen
      icoo(i)=hicoo(i)
      jcoo(i)=hjcoo(i)
      posij= hirv(1:3, i)
      kdotr=posij(1)*k(1)+posij(2)*k(2)+posij(3)*k(3)
      ratio= (cos(twopi*kdotr)+zi*sin(twopi*kdotr))
      acoo(i)=ratio*hacoo(i)
   end do

   return
end subroutine ham_bulk_coo_sparsehr_latticegauge

subroutine ham_bulk_coo_sparsehr(k,acoo,icoo,jcoo)
   !> This subroutine use sparse hr format
   !> Here we use atomic gauge which means the atomic position is taken into account
   !> in the Fourier transformation
   use para
   implicit none

   real(dp) :: k(3), posij(3)
   real(dp) :: kdotr
   integer,intent(inout):: icoo(splen),jcoo(splen)!,splen
   complex(dp),intent(inout) :: acoo(splen)
   complex(dp) ::  ratio
   integer :: i,j,ir

   do i=1, splen
      icoo(i)=hicoo(i)
      jcoo(i)=hjcoo(i)
      posij= hirv(1:3, i)+ Origin_cell%wannier_centers_direct(:, jcoo(i))- Origin_cell%wannier_centers_direct(:, icoo(i))
      kdotr=posij(1)*k(1)+posij(2)*k(2)+posij(3)*k(3)
      ratio= (cos(twopi*kdotr)+zi*sin(twopi*kdotr))
      acoo(i)=ratio*hacoo(i)
   end do

   return
end subroutine ham_bulk_coo_sparsehr


subroutine overlap_bulk_coo_sparse(k, acoo, icoo, jcoo)
   !> This subroutine use sparse hr format
   !> Here we use atomic gauge which means the atomic position is taken into account
   !> in the Fourier transformation
   use para
   implicit none

   real(dp) :: k(3), posij(3)
   real(dp) :: kdotr
   integer,intent(inout) :: icoo(splen_overlap_input),jcoo(splen_overlap_input)
   complex(dp),intent(inout) :: acoo(splen_overlap_input)
   complex(dp) ::  ratio
   integer :: i,j,ir

   do i=1, splen_overlap_input
      icoo(i)=sicoo(i)
      jcoo(i)=sjcoo(i)
      posij= sirv(1:3, i)+ Origin_cell%wannier_centers_direct(:, sjcoo(i))- Origin_cell%wannier_centers_direct(:, sicoo(i))
      kdotr=posij(1)*k(1)+posij(2)*k(2)+posij(3)*k(3)
      ratio= (cos(twopi*kdotr)+zi*sin(twopi*kdotr))
      acoo(i)=ratio*sacoo(i)
   end do

   return
end subroutine overlap_bulk_coo_sparse

subroutine ham_bulk_tpaw(mdim, k, Hamk_moire)
   ! This subroutine caculates Hamiltonian for
   ! a twist system
   !
   !  Lattice gauge Hl
   !  Atomic gauge Ha= U* Hl U 
   !  where U = e^ik.wc(i) on diagonal
   ! contact wuquansheng@gmail.com
   ! 2024 Jan 19

   use para
   implicit none

   !> leading dimension of Ham_moire
   integer, intent(in) :: mdim
   real(dp), intent(in) :: k(3)
   complex(dp), intent(inout) :: hamk_moire(mdim, mdim)

   ! wave vector in 3d
   real(Dp) :: kdotr, pos0(3), dis

   complex(dp) :: factor

   real(dp) :: pos(3), pos1(3), pos2(3), pos_cart(3), pos_direct(3)

   integer :: Num_gvectors

   !> dimension 3*Num_gvectors 
   real(dp), allocatable :: gvectors(:, :)

   integer :: i1, i2

   !> for a test, we assume a twist homo-bilayer system
   ! mdim= Num_gvectors* Folded_cell%NumberOfspinorbitals*2

   hamk_moire=0d0

   !> first, set G vectors

   !> construct T matrix
   !> the dimension of the T matrix is num_wann
   

   ! check hermitcity
   do i1=1, mdim
      do i2=1, mdim
         if(abs(hamk_moire(i1,i2)-conjg(hamk_moire(i2,i1))).ge.1e-6)then
            write(stdout,*)'there is something wrong with hamk_moire'
            write(stdout,*)'i1, i2', i1, i2
            write(stdout,*)'value at (i1, i2)', hamk_moire(i1, i2)
            write(stdout,*)'value at (i2, i1)', hamk_moire(i2, i1)
            !stop
         endif
      enddo
   enddo

   return
end subroutine ham_bulk_tpaw


subroutine valley_k_coo_sparsehr(nnz, k,acoo,icoo,jcoo)
   !> This subroutine use sparse hr format
   !> Here we use atomic gauge which means the atomic position is taken into account
   !> in the Fourier transformation
   use para
   implicit none

   real(dp) :: k(3), posij(3)
   real(dp) :: kdotr
   integer, intent(in) :: nnz
   integer,intent(inout) :: icoo(nnz),jcoo(nnz)
   complex(dp),intent(inout) :: acoo(nnz)
   complex(dp) ::  ratio
   integer :: i,j,ir

   do i=1,nnz
      ir= valley_operator_irv(i)
      icoo(i)=valley_operator_icoo(i)
      jcoo(i)=valley_operator_jcoo(i)
      posij=irvec_valley(:, ir)+ Origin_cell%wannier_centers_direct(:, jcoo(i))- Origin_cell%wannier_centers_direct(:, icoo(i))
      kdotr=posij(1)*k(1)+posij(2)*k(2)+posij(3)*k(3)
      ratio= (cos(twopi*kdotr)+zi*sin(twopi*kdotr))
      acoo(i)=ratio*valley_operator_acoo(i)
   enddo

   return
end subroutine valley_k_coo_sparsehr


subroutine rotation_to_Ham_basis(UU, mat_wann, mat_ham)
   !> this subroutine rotate the matrix from Wannier basis to Hamiltonian basis
   !> UU are the eigenvectors from the diagonalization of Hamiltonian
   !> mat_ham=UU_dag*mat_wann*UU
   use para, only : dp, Num_wann
   implicit none
   complex(dp), intent(in) :: UU(Num_wann, Num_wann)
   complex(dp), intent(in) :: mat_wann(Num_wann, Num_wann)
   complex(dp), intent(out) :: mat_ham(Num_wann, Num_wann)
   complex(dp), allocatable :: UU_dag(:, :), mat_temp(:, :)

   allocate(UU_dag(Num_wann, Num_wann), mat_temp(Num_wann, Num_wann))
   UU_dag= conjg(transpose(UU))

   call mat_mul(Num_wann, mat_wann, UU, mat_temp) 
   call mat_mul(Num_wann, UU_dag, mat_temp, mat_ham) 

   return
end subroutine rotation_to_Ham_basis



