  subroutine ham_slab(k,Hamk_slab)
     ! This subroutine is used to caculate Hamiltonian for 
     ! slab system . 
     ! 
     ! History  
     !        4/18/2010 by Quansheng Wu
     ! Copyright (c) 2010 QuanSheng Wu. All rights reserved.
  
     use para
     implicit none

     ! loop index  
     integer :: i1, i2, i, j
     integer :: ir

     ! wave vector in 2d
     real(Dp), intent(inout) :: k(2)      

     ! Hamiltonian of slab system
     complex(Dp),intent(out) ::Hamk_slab(Num_wann*nslab,Num_wann*nslab) 

     ! the factor 2 is induced by spin
     complex(Dp), allocatable :: Hij(:, :, :)

     allocate( Hij(-ijmax:ijmax,Num_wann,Num_wann))


     !> deal with phonon system
     if (index(Particle,'phonon')/=0.and.LOTO_correction) then
        call ham_qlayer2qlayer2_LOTO(k,Hij)
     else
        call ham_qlayer2qlayer2(k,Hij)
     endif


     Hamk_slab=0.0d0 
     ! i1 column index
     do i1=1, nslab
        ! i2 row index
        do i2=1, nslab
          if (abs(i2-i1).le.ijmax)then
            Hamk_slab((i2-1)*Num_wann+1:(i2-1)*Num_wann+Num_wann,&
                      (i1-1)*Num_wann+1:(i1-1)*Num_wann+Num_wann )&
            = Hij(i2-i1,1:Num_wann,1:Num_wann)
          endif 
        enddo ! i2
     enddo ! i1

     ! check hermitcity

     do i1=1,nslab*Num_wann
     do i2=1,nslab*Num_wann
        if(abs(Hamk_slab(i1,i2)-conjg(Hamk_slab(i2,i1))).ge.1e-6)then
         !write(stdout,*)'there are something wrong with Hamk_slab'
         !stop
        endif 
     enddo
     enddo

     deallocate( Hij)
  return
  end subroutine ham_slab


  subroutine ham_slab_parallel_B(k,Hamk_slab)
     ! This subroutine is used to caculate Hamiltonian for 
     ! slab system . 
     !> for slab with in-plane magnetic field
     !> the magnetic phase are chosen like this
     !> phi_ij= alpha*[By*(xj-xi)*(zi+zj)-Bx*(yj-yi)*(zi+zj)] 
     !> x, z are in unit of Angstrom, Bx, By are in unit of Tesla
     !> History :
     !        9/21/2015 by Quansheng Wu @ETH Zurich
  
     use para
     implicit none

     ! loop index  
     integer :: i1, i2, i, j

     ! wave vector in 2d
     real(Dp), intent(inout) :: k(2)      

     ! loop index
     integer :: iR

     ! index used to sign irvec     
     integer :: ia,ib,ic
     integer :: ia1, ia2

     integer :: istart1, istart2
     integer :: iend1, iend2

     integer :: inew_ic

     !> nwann= Num_wann/2
     integer :: nwann
     
     integer, allocatable :: orbital_start(:)

     ! new index used to sign irvec     
     real(dp) :: new_ia,new_ib,new_ic

     ! wave vector k times lattice vector R  
     real(Dp) :: kdotr  
     real(dp) :: phase
     complex(dp) :: ratio
     complex(dp) :: fac

     real(dp) :: R(3)
     real(dp) :: Rp1(3)
     real(dp) :: Rp2(3)
     real(dp) :: R1(3)
     real(dp) :: R2(3)
     real(dp) :: Ri(3)
     real(dp) :: Rj(3)
     real(dp) :: tau1(3)
     real(dp) :: tau2(3)


     ! Hamiltonian of slab system
     complex(Dp),intent(out) ::Hamk_slab(Num_wann*nslab,Num_wann*nslab) 

     nwann= Num_wann/2
     allocate( orbital_start(Num_atoms+ 1))
     orbital_start= 0
     orbital_start(1)= 1
     do ia=1, Num_atoms
        orbital_start(ia+1)= orbital_start(ia)+ nprojs(ia)
     enddo

     Hamk_slab=0.0d0 
     ! i1 column index
     do i1=1, Nslab
        ! i2 row index
        do i2=1, Nslab
           !> sum over R points to get H(k1, k2)
           do iR=1, Nrpts
              ia=irvec(1,iR)
              ib=irvec(2,iR)
              ic=irvec(3,iR)
      
              !> new lattice
              call latticetransform(ia, ib, ic, new_ia, new_ib, new_ic)
      
              inew_ic= int(new_ic)
              if (inew_ic /= (i2-i1)) cycle
      
              !> exp(i k.R)
              kdotr= k(1)*new_ia+ k(2)*new_ib
              ratio= cos(2d0*pi*kdotr)+ zi*sin(2d0*pi*kdotr)
      
              R1= (i1-1)*Ruc_new
              R2= new_ia*Rua_new+ new_ib*Rub_new+ (i2-1)*Ruc_new
      
              do ia1=1, Num_atoms
              do ia2=1, Num_atoms
                 R1= Atom_position(:, ia1)
                 R2= Atom_position(:, ia2)
                 call rotate(R1, tau1)
                 call rotate(R2, tau2)
      
                
                 Ri= Rp1+ tau1
                 Rj= Rp2+ tau2
      
                 phase= alpha*By*(Rj(3)+Ri(3))*(Rj(1)-Ri(1))  &
                      - alpha*Bx*(Rj(3)+Ri(3))*(Rj(2)-Ri(2))
                 fac= cos(phase)+ zi*sin(phase)
      
                !write(*, '(a, 4i5   )') 'iR, ia ib ic', ir, ia, ib, ic
                !write(*, '(a, 4f10.5)') '            ', new_ia, new_ib, new_ic
                !write(*, '(a, 3f10.5)') 'Ri', Ri
                !write(*, '(a, 3f10.5)') 'Rj', Rj
                !write(*, '(a, 3f10.5)') 'phase', phase
      
                 istart1= (i2-1)*Num_wann+ orbital_start(ia1)
                 iend1= (i2-1)*Num_wann+ orbital_start(ia1+1)- 1 
                 istart2= (i1-1)*Num_wann+ orbital_start(ia2)
                 iend2= (i1-1)*Num_wann+ orbital_start(ia2+1)- 1
                 
                 Hamk_slab( istart1:iend1, istart2:iend2) &
                 = Hamk_slab( istart1:iend1, istart2:iend2) &
                 + HmnR( istart1- (i2-1)*Num_wann:iend1- (i2-1)*Num_wann, &
                 istart2- (i1-1)*Num_wann:iend2- (i1-1)*Num_wann, iR)*ratio/ndegen(iR)* fac
      
                 !> there is soc term in the hr file
                 if (soc>0) then
                    istart1= (i2-1)*Num_wann+ orbital_start(ia1) + Nwann
                    iend1= (i2-1)*Num_wann+ orbital_start(ia1+1)- 1 + Nwann 
                    istart2= (i1-1)*Num_wann+ orbital_start(ia2)
                    iend2= (i1-1)*Num_wann+ orbital_start(ia2+1)- 1
                    
                    Hamk_slab( istart1:iend1, istart2:iend2) &
                    = Hamk_slab( istart1:iend1, istart2:iend2) &
                    + HmnR( istart1- (i2-1)*Num_wann:iend1- (i2-1)*Num_wann, &
                    istart2- (i1-1)*Num_wann:iend2- (i1-1)*Num_wann, iR)*ratio/ndegen(iR)* fac
      
                    istart1= (i2-1)*Num_wann+ orbital_start(ia1)
                    iend1= (i2-1)*Num_wann+ orbital_start(ia1+1)- 1 
                    istart2= (i1-1)*Num_wann+ orbital_start(ia2) + Nwann
                    iend2= (i1-1)*Num_wann+ orbital_start(ia2+1)- 1 + Nwann
                    
                    Hamk_slab( istart1:iend1, istart2:iend2) &
                    = Hamk_slab( istart1:iend1, istart2:iend2) &
                    + HmnR( istart1- (i2-1)*Num_wann:iend1- (i2-1)*Num_wann, &
                    istart2- (i1-1)*Num_wann:iend2- (i1-1)*Num_wann, iR)*ratio/ndegen(iR)* fac
      
                    istart1= (i2-1)*Num_wann+ orbital_start(ia1) + Nwann
                    iend1= (i2-1)*Num_wann+ orbital_start(ia1+1)- 1 + Nwann 
                    istart2= (i1-1)*Num_wann+ orbital_start(ia2) + Nwann
                    iend2= (i1-1)*Num_wann+ orbital_start(ia2+1)- 1 + Nwann
                    
                    Hamk_slab( istart1:iend1, istart2:iend2) &
                    = Hamk_slab( istart1:iend1, istart2:iend2) &
                    + HmnR( istart1- (i2-1)*Num_wann:iend1- (i2-1)*Num_wann, &
                    istart2- (i1-1)*Num_wann:iend2- (i1-1)*Num_wann, iR)*ratio/ndegen(iR)* fac
                 endif ! soc
              enddo ! ia2
              enddo ! ia1
           enddo ! iR
        enddo ! i2
     enddo ! i1

 

     !> check hermitcity
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

