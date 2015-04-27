! This subroutine is used to caculate Hamiltonian for 
! bulk system . 

! History  
!        May/29/2011 by Quansheng Wu
  subroutine ham_bulk(k,Hamk_bulk)
  
     use para
     implicit none

! loop index  
     integer :: i1,i2,ia,ib,ic,iR

     integer :: ia1, ia2

     integer :: nwann

     real(dp) :: kdotr

! wave vector in 2d
     real(Dp) :: k(3)      

     ! coordinates of R vector 
     real(Dp) :: R(3)      

     !> row start; row end; column start; column end
     real(dp) :: rs, re
     real(dp) :: cs, ce

     real(dp) :: pos(3)
     real(dp) :: pos1(3)
     real(dp) :: pos2(3)

     !> distance between two atoms
     real(dp) :: dis

     ! Hamiltonian of bulk system
     complex(Dp),intent(out) ::Hamk_bulk(Num_wann, Num_wann) 

     real(dp), external :: norm

     nwann= Num_wann/2
     Hamk_bulk=0d0
     if (soc>0) then
        !> the first atom in home unit cell
        do ia1=1, Num_atoms
           pos1= Atom_position(:, ia1)
           rs= index_start(ia1)
           re= index_end(ia1)
           !> the second atom in unit cell R
           do ia2=1, Num_atoms
              pos2= Atom_position(:, ia2)
              cs= index_start(ia2)
              ce= index_end(ia2)
              do iR=1,Nrpts
                 ia=irvec(1,iR)
                 ib=irvec(2,iR)
                 ic=irvec(3,iR)
         
                 pos= ia*Rua+ ib*Rub+ ic* Ruc
                 pos= pos+ pos2- pos1
                 dis= norm(pos)
         
                 if (dis> Rcut) cycle
        
                 kdotr=k(1)*ia + k(2)*ib + k(3)*ic
               
                 Hamk_bulk(rs:re, cs:ce)=Hamk_bulk(rs:re, cs:ce)&
                 +HmnR(rs:re, cs:ce,iR)*Exp(2d0*pi*zi*kdotr)/ndegen(iR)

                 Hamk_bulk(rs:re, cs+nwann:ce+nwann)= &
                 Hamk_bulk(rs:re, cs+nwann:ce+nwann)+ &
                 HmnR(rs:re, cs+nwann:ce+nwann,iR)* &
                 Exp(2d0*pi*zi*kdotr)/ndegen(iR)

                 Hamk_bulk(rs+nwann:re+nwann, cs:ce)= &
                 Hamk_bulk(rs+nwann:re+nwann, cs:ce)+ &
                 HmnR(rs+nwann:re+nwann, cs:ce,iR)* &
                 Exp(2d0*pi*zi*kdotr)/ndegen(iR)

                 Hamk_bulk(rs+nwann:re+nwann, cs+nwann:ce+nwann)= &
                 Hamk_bulk(rs+nwann:re+nwann, cs+nwann:ce+nwann)+ &
                 HmnR(rs+nwann:re+nwann, cs+nwann:ce+nwann,iR)* &
                 Exp(2d0*pi*zi*kdotr)/ndegen(iR)

              enddo ! iR
           enddo ! ia2
        enddo ! ia1
     else
        !> the first atom in home unit cell
        do ia1=1, Num_atoms
           pos1= Atom_position(:, ia1)
           rs= index_start(ia1)
           re= index_end(ia1)
           !> the second atom in unit cell R
           do ia2=1, Num_atoms
              pos2= Atom_position(:, ia2)
              cs= index_start(ia2)
              ce= index_end(ia2)
              do iR=1,Nrpts
                 ia=irvec(1,iR)
                 ib=irvec(2,iR)
                 ic=irvec(3,iR)
         
                 pos= ia*Rua+ ib*Rub+ ic* Ruc
                 pos= pos+ pos2- pos1
                 dis= norm(pos)
         
                 if (dis> Rcut) cycle
        
                 kdotr=k(1)*R(1) + k(2)*R(2) + k(3)*R(3)
               
                 Hamk_bulk(rs:re, cs:ce)=Hamk_bulk(rs:re, cs:ce)&
                 +HmnR(rs:re, cs:ce,iR)*Exp(2d0*pi*zi*kdotr)/ndegen(iR)

              enddo ! iR
           enddo ! ia2
        enddo ! ia1

     endif

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
  end
