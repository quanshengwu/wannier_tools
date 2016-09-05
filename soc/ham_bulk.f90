! This subroutine is used to caculate Hamiltonian for 
! bulk system . 

! History  
!        May/29/2011 by Quansheng Wu
! Copyright (c) 2010 QuanSheng Wu. All rights reserved.
  subroutine ham_bulk(k,Hamk_bulk)
  
     use para
     implicit none

! loop index  
     integer :: i1,i2,ia,ib,ic,iR, ia1, ia2
     integer :: istart1, istart2, iend1, iend2
     integer :: nwann

     real(dp) :: kdotr

! wave vector in 2d
     real(Dp) :: k(3)      

     ! coordinates of R vector 
     real(Dp) :: R(3), R1(3), R2(3), R3(3)

! Hamiltonian of bulk system
     complex(Dp),intent(out) ::Hamk_bulk(Num_wann, Num_wann) 
     complex(dp), allocatable :: mat1(:, :)
     complex(dp), allocatable :: mat2(:, :)
     allocate(mat1(Num_wann, Num_wann))
     allocate(mat2(Num_wann, Num_wann))

     nwann= Num_wann/2
     Hamk_bulk=0d0
     do iR=1, Nrpts
        ia=irvec(1,iR)  
        ib=irvec(2,iR) 
        ic=irvec(3,iR)

        R(1)=dble(ia)
        R(2)=dble(ib)
        R(3)=dble(ic)


        do ia1= 1, Num_atoms
           istart1= sum(nprojs(1:ia1-1))+1
           iend1= sum(nprojs(1:ia1))
           do ia2= 1, Num_atoms
              istart2= sum(nprojs(1:ia2-1))+1
              iend2= sum(nprojs(1:ia2))
              R1= Atom_position_direct(:, ia1)
              R2= Atom_position_direct(:, ia2)
              R3= R+ R1- R2
              kdotr=k(1)*R3(1) + k(2)*R3(2) + k(3)*R3(3)
             !kdotr=k(1)*R (1) + k(2)*R (2) + k(3)*R (3)

              Hamk_bulk(istart1:iend1,istart2:iend2)=&
              Hamk_bulk(istart1:iend1,istart2:iend2) &
              +HmnR(istart1:iend1,istart2:iend2,iR)*Exp(2d0*pi*zi*kdotr)/ndegen(iR)

              if (soc>0) then
                  Hamk_bulk(istart1+nwann:iend1+nwann,istart2:iend2)=&
                  Hamk_bulk(istart1+nwann:iend1+nwann,istart2:iend2) &
                  +HmnR(istart1+nwann:iend1+nwann,istart2:iend2,iR)*Exp(2d0*pi*zi*kdotr)/ndegen(iR)

                  Hamk_bulk(istart1+nwann:iend1+nwann,istart2+nwann:iend2+nwann)=&
                  Hamk_bulk(istart1+nwann:iend1+nwann,istart2+nwann:iend2+nwann) &
                  +HmnR(istart1+nwann:iend1+nwann,istart2+nwann:iend2+nwann,iR)*Exp(2d0*pi*zi*kdotr)/ndegen(iR)

                  Hamk_bulk(istart1:iend1,istart2+nwann:iend2+nwann)=&
                  Hamk_bulk(istart1:iend1,istart2+nwann:iend2+nwann) &
                  +HmnR(istart1:iend1,istart2+nwann:iend2+nwann,iR)*Exp(2d0*pi*zi*kdotr)/ndegen(iR)

              endif

           enddo ! ia2
        enddo ! ia1
     enddo ! iR

    !call mat_mul(Num_wann, mirror_z, Hamk_bulk, mat1)
    !call mat_mul(Num_wann, mat1, mirror_z, mat2)
    !Hamk_bulk= (Hamk_bulk+ mat2)/2d0

     ! check hermitcity
     do i1=1, Num_wann
     do i2=1, Num_wann
        if(abs(Hamk_bulk(i1,i2)-conjg(Hamk_bulk(i2,i1))).ge.1e-6)then
          write(stdout,*)'there are something wrong with Hamk_bulk'
          stop
        endif 
     enddo
     enddo

  return
  end  subroutine ham_bulk



! This subroutine is used to caculate Hamiltonian for 
! bulk system . 

! History  
!        May/29/2011 by Quansheng Wu
  subroutine ham_bulk_old(k,Hamk_bulk)
  
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

! Hamiltonian of bulk system
     complex(Dp),intent(out) ::Hamk_bulk(Num_wann, Num_wann) 
     complex(dp), allocatable :: mat1(:, :)
     complex(dp), allocatable :: mat2(:, :)
     allocate(mat1(Num_wann, Num_wann))
     allocate(mat2(Num_wann, Num_wann))

     Hamk_bulk=0d0
     do iR=1, Nrpts
        ia=irvec(1,iR)  
        ib=irvec(2,iR) 
        ic=irvec(3,iR)

        R(1)=dble(ia)
        R(2)=dble(ib)
        R(3)=dble(ic)
        kdotr=k(1)*R (1) + k(2)*R (2) + k(3)*R (3)

        Hamk_bulk(:, :)=&
        Hamk_bulk(:, :) &
        + HmnR(:, :, iR)*Exp(2d0*pi*zi*kdotr)/ndegen(iR)
     enddo ! iR

    !call mat_mul(Num_wann, mirror_z, Hamk_bulk, mat1)
    !call mat_mul(Num_wann, mat1, mirror_z, mat2)
    !Hamk_bulk= (Hamk_bulk+ mat2)/2d0

     ! check hermitcity
     do i1=1, Num_wann
     do i2=1, Num_wann
        if(abs(Hamk_bulk(i1,i2)-conjg(Hamk_bulk(i2,i1))).ge.1e-6)then
          write(stdout,*)'there are something wrong with Hamk_bulk'
          stop
        endif 
     enddo
     enddo

  return
  end subroutine ham_bulk_old

