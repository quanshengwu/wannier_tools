  !> calculate landau levels in 3D system in special k line
  !> fix B field
  !> the magnetic field is in the a1 a2 plane
  !> B= (B cos\theta, B sin\theta, 0)
  !> construct on Dec 8 2015
  !> By QuanSheng Wu at ETH Zurich Honggerberg
! Copyright (c) 2010 QuanSheng Wu. All rights reserved.
  subroutine landau_level_k
     use para
     implicit none


     !> magnetic supercell size, perpendicular to the magnetic field
     integer :: Nq
     integer :: ik
     integer :: ia1
     integer :: ia2
     integer :: i, j

     !> Ndimq= Nq* Num_wann
     integer :: Ndimq

     integer :: nlines

     integer :: ierr

     real(dp) :: B0
     real(dp) :: theta


     real(dp) :: t1, temp, dis, dis1
     ! wave vector 
     real(dp) :: k3(3)

     !> dim= Ndimq, knv3
     real(dp), allocatable :: W(:)
     real(dp), allocatable :: eigv(:, :)
     real(dp), allocatable :: eigv_mpi(:, :)

     !> dim= Ndimq*Ndimq
     complex(dp), allocatable :: ham_landau(:, :)


     Nq= nslab
     Ndimq= Num_wann* Nq
     allocate( ham_landau(Ndimq, Ndimq))
     allocate( W( Ndimq))
     allocate( eigv( Ndimq, knv2))
     allocate( eigv_mpi( Ndimq, knv2))
     eigv_mpi= 0d0
     eigv    = 0d0
     ham_landau= 0d0


     !> deal with the magnetic field
     !> first transform the Bx By into B*Cos\theta, B*Sin\theta
     if (abs(By)<1e-8) then
        theta= 0d0
     elseif (By>0) then
        theta = atan(Bx/By)
     else
        theta = atan(Bx/By)+ pi
     endif
     if (abs(theta)<1e-8 .and. Bx<0 ) theta= theta+ pi 

     !> get the shortest bond in the home unit cell
     dis1= sqrt(sum((Rua_new)**2))
     do ia1=1, Num_atoms
        do ia2=ia1+1, Num_atoms
           dis = (Rua_new(2))*Cos(theta)- (Rua_new(1))*Sin(theta) 
           if (abs(dis)< abs(dis1) .and. abs(dis)>1e-5) dis1= abs(dis)
           dis = (Rub_new(2))*Cos(theta)- (Rub_new(1))*Sin(theta) 
           if (abs(dis)< abs(dis1) .and. abs(dis)>1e-5) dis1= abs(dis)
        enddo
     enddo

     !> The flux in the supercell should be 2*pi
     if (dis1< 1e-9) stop 'something wrong with the atom position'
     B0= 2d0*pi/ (dis1* Ruc_new(3)*Nq)
     B0= abs(B0)
     Bx= B0* Cos(theta)
     By= B0* Sin(theta)

     if (cpuid==0) then
        write(stdout, '(a, 2f18.8)')' Bx, By= ', Bx, By
     endif

     !> calculate the landau levels along special k line
     k3= 0
     do ik=1+ cpuid, knv2, num_cpu

        if (cpuid==0) print *, ik, knv2
        k3(1:2) = k2_path(ik, :)
        call ham_3Dlandau2(Ndimq, Nq, k3, ham_landau)
        
        !> diagonalization by call zheev in lapack
        W= 0d0
        call eigensystem_c( 'N', 'U', Ndimq ,ham_landau, W)

        eigv(:, ik)= W
     enddo !ik

     call mpi_allreduce(eigv,eigv_mpi,size(eigv),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)

     if (cpuid.eq.0) then
        open(unit=100, file='landaulevel_k.dat')
        do j=1, Ndimq
           do i=1, knv2
              write(100,'(2f15.7, i8)')k2len(i), eigv_mpi(j, i)
           enddo
           write(100 , *)''
        enddo
        close(100)
        write(stdout,*) 'calculate landau level done'
     endif

     return
  end subroutine landau_level_k

  !> calculate landau levels in 3D system for different B
  !> fix k point, usually Gamma point
  subroutine landau_level_B
     use para
     implicit none


     !> magnetic supercell size, perpendicular to the magnetic field
     integer :: Nq
     integer :: ia1
     integer :: ia2
     integer :: ib
     integer :: i, j

     !> Ndimq= Nq* Num_wann
     integer :: Ndimq

     integer :: Nmag

     integer :: ierr

     real(dp) :: B0
     real(dp) :: theta

     real(dp) :: dis, dis1
     ! wave vector 
     real(dp) :: k3(3)

     !> dim= Ndimq, knv3
     real(dp), allocatable :: mag(:)
     real(dp), allocatable :: W(:)
     real(dp), allocatable :: eigv(:, :)
     real(dp), allocatable :: eigv_mpi(:, :)

     !> dim= Ndimq*Ndimq
     complex(dp), allocatable :: ham_landau(:, :)

     Nq= nslab
     Nmag= Nq
     Ndimq= Num_wann* Nq
     allocate( ham_landau(Ndimq, Ndimq))
     allocate( W( Ndimq))
     allocate( eigv( Ndimq, Nmag))
     allocate( eigv_mpi( Ndimq, Nmag))
     allocate( mag(Nq))
     mag= 0d0
     eigv_mpi= 0d0
     eigv    = 0d0
     ham_landau= 0d0

     !> deal with the magnetic field
     !> first transform the Bx By into B*Cos\theta, B*Sin\theta
     if (abs(By)<1e-8) then
        theta= 0d0
     elseif (By>0) then
        theta = atan(Bx/By)
     else
        theta = atan(Bx/By)+ pi
     endif
     if (abs(theta)<1e-8 .and. Bx<0 ) theta= theta+ pi 

     !> get the shortest bond in the home unit cell
     dis1= sqrt(sum((Rua_new)**2))
     do ia1=1, Num_atoms
        do ia2=ia1+1, Num_atoms
           dis= (Rua_new(2))*Cos(theta)- (Rua_new(1))*Sin(theta) 
           if (abs(dis)< abs(dis1) .and. abs(dis)>1e-5) dis1= abs(dis)
           dis= (Rub_new(2))*Cos(theta)- (Rub_new(1))*Sin(theta) 
           if (abs(dis)< abs(dis1) .and. abs(dis)>1e-5) dis1= abs(dis)
        enddo
     enddo


     !> The flux in the supercell should be 2*pi
     if (dis1< 1e-9) stop 'something wrong with the atom position'
     B0= 2d0*pi/ (dis1* Ruc_new(3)*Nq)
     B0= abs(B0)

     do ib=1, Nmag
        mag(ib)= B0* ib
     enddo

     k3= 0d0
     !> calculate the landau levels along special k line
     do ib=1+ cpuid, Nmag, num_cpu

        Bx= mag(ib)* Cos(theta)
        By= mag(ib)* Sin(theta)
        
        if (cpuid==0) print *, ib, Nmag
        call ham_3Dlandau2(Ndimq, Nq, k3, ham_landau)
        
        !> diagonalization by call zheev in lapack
        W= 0d0
        call eigensystem_c( 'N', 'U', Ndimq ,ham_landau, W)

        eigv(:, ib)= W
     enddo !ik

     call mpi_allreduce(eigv,eigv_mpi,size(eigv),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)

     if (cpuid.eq.0) then
        open(unit=101, file='landaulevel_B.dat')
        do j=1, Ndimq
           do i=1, Nmag
              write(101,'(2f15.7, i8)')mag(i), eigv_mpi(j, i)
           enddo
           write(101 , *)''
        enddo
        close(101)
        write(stdout,*) 'calculate landau level done'
     endif

     return

  end subroutine landau_level_B

  !> calculate hamiltonian for landau levels
  !> consider the internal atom's position
  subroutine ham_3Dlandau(Ndimq, Nq, k, ham_landau)
     use para
     implicit none

     integer, intent(in) :: Ndimq
     integer, intent(in) :: Nq
     real(dp), intent(in) :: k(3)
     complex(dp), intent(out) :: ham_landau(Ndimq, Ndimq)

     !> inta-hopping for the supercell
     complex(dp), allocatable :: H00(:, :)

     !> inter-hopping for the supercell
     complex(dp), allocatable :: H01(:, :)

     ! loop index  
     integer :: i1, i2

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

     real(dp) :: Rp1(3)
     real(dp) :: Rp2(3)
     real(dp) :: R1(3)
     real(dp) :: R2(3)
     real(dp) :: Ri(3)
     real(dp) :: Rj(3)
     real(dp) :: tau1(3)
     real(dp) :: tau2(3)

     allocate( H00( Ndimq, Ndimq))
     allocate( H01( Ndimq, Ndimq))

     nwann= Num_wann/2
     allocate( orbital_start(Num_atoms+ 1))
     orbital_start= 0
     orbital_start(1)= 1
     do ia=1, Num_atoms
        orbital_start(ia+1)= orbital_start(ia)+ nprojs(ia)
     enddo

     !> calculate intra-hopping
     H00=0.0d0 
     ! i1 column index
     do i1=1, Nq
        ! i2 row index
        do i2=1, Nq
           if (abs(i2-i1)> ijmax) cycle
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
              Rp1= R1
              Rp2= R2
      
              do ia1=1, Num_atoms
              do ia2=1, Num_atoms
                 R1= Atom_position(:, ia1)
                 R2= Atom_position(:, ia2)
                 call rotate(R1, tau1)
                 call rotate(R2, tau2)
      
                
                 Ri= Rp1+ tau1
                 Rj= Rp2+ tau2

                !write(*, '(3i3, 3(2x, 3f7.3))')ia, ib, ic, Ri, &
                !Rj-Ri, Ri+Rj
                !write(*, '(3i3, 3(2x, 3f7.3))')ia, ib, ic, Rp2, tau2, Rj
      
                 phase= By*(Rj(3)+Ri(3))*(Rj(1)-Ri(1))  &
                      - Bx*(Rj(3)+Ri(3))*(Rj(2)-Ri(2))
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
                 
                 H00( istart1:iend1, istart2:iend2) &
                 = H00( istart1:iend1, istart2:iend2) &
                 + HmnR( istart1- (i2-1)*Num_wann:iend1- (i2-1)*Num_wann, &
                 istart2- (i1-1)*Num_wann:iend2- (i1-1)*Num_wann, iR)*ratio/ndegen(iR)* fac
      
                 !> there is soc term in the hr file
                 if (soc>0) then
                    istart1= (i2-1)*Num_wann+ orbital_start(ia1) + Nwann
                    iend1= (i2-1)*Num_wann+ orbital_start(ia1+1)- 1 + Nwann 
                    istart2= (i1-1)*Num_wann+ orbital_start(ia2)
                    iend2= (i1-1)*Num_wann+ orbital_start(ia2+1)- 1
                    
                    H00( istart1:iend1, istart2:iend2) &
                    = H00( istart1:iend1, istart2:iend2) &
                    + HmnR( istart1- (i2-1)*Num_wann:iend1- (i2-1)*Num_wann, &
                    istart2- (i1-1)*Num_wann:iend2- (i1-1)*Num_wann, iR)*ratio/ndegen(iR)* fac
      
                    istart1= (i2-1)*Num_wann+ orbital_start(ia1)
                    iend1= (i2-1)*Num_wann+ orbital_start(ia1+1)- 1 
                    istart2= (i1-1)*Num_wann+ orbital_start(ia2) + Nwann
                    iend2= (i1-1)*Num_wann+ orbital_start(ia2+1)- 1 + Nwann
                    
                    H00( istart1:iend1, istart2:iend2) &
                    = H00( istart1:iend1, istart2:iend2) &
                    + HmnR( istart1- (i2-1)*Num_wann:iend1- (i2-1)*Num_wann, &
                    istart2- (i1-1)*Num_wann:iend2- (i1-1)*Num_wann, iR)*ratio/ndegen(iR)* fac
      
                    istart1= (i2-1)*Num_wann+ orbital_start(ia1) + Nwann
                    iend1= (i2-1)*Num_wann+ orbital_start(ia1+1)- 1 + Nwann 
                    istart2= (i1-1)*Num_wann+ orbital_start(ia2) + Nwann
                    iend2= (i1-1)*Num_wann+ orbital_start(ia2+1)- 1 + Nwann
                    
                    H00( istart1:iend1, istart2:iend2) &
                    = H00( istart1:iend1, istart2:iend2) &
                    + HmnR( istart1- (i2-1)*Num_wann:iend1- (i2-1)*Num_wann, &
                    istart2- (i1-1)*Num_wann:iend2- (i1-1)*Num_wann, iR)*ratio/ndegen(iR)* fac
                 endif ! soc
              enddo ! ia2
              enddo ! ia1
           enddo ! iR
        enddo ! i2
     enddo ! i1

 

     !> check hermitcity
     do i1=1,Nq*Num_wann
     do i2=1,Nq*Num_wann
       !if(abs(H00(i1,i2)-conjg(H00(i2,i1))).ge.1e-6)then
       !  write(stdout,*)'there are something wrong with H00'
       !  print *, i1, i2, H00(i1, i2)
       !  stop
       !endif 
     enddo
     enddo

     !>> calculate inter-hopping
     H01=0.0d0 
     ! i1 column index
     do i1=1, Nq
        ! i2 row index
        do i2=1, Nq
           !> sum over R points to get H(k1, k2)
           do iR=1, Nrpts
              ia=irvec(1,iR)
              ib=irvec(2,iR)
              ic=irvec(3,iR)
      
              !> new lattice
              call latticetransform(ia, ib, ic, new_ia, new_ib, new_ic)
      
              inew_ic= int(new_ic)
              if (inew_ic /= (i2+ Nq -i1)) cycle
      
              !> exp(i k.R)
              kdotr= k(1)*new_ia+ k(2)*new_ib
              ratio= cos(2d0*pi*kdotr)+ zi*sin(2d0*pi*kdotr)
      
              R1= (i1-1)*Ruc_new
              R2= new_ia*Rua_new+ new_ib*Rub_new+ (i2-1+ Nq)*Ruc_new
              Rp1= R1
              Rp2= R2
             !call rotate(R1, Rp1)
             !call rotate(R2, Rp2)
      
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
                 
                 H01( istart1:iend1, istart2:iend2) &
                 = H01( istart1:iend1, istart2:iend2) &
                 + HmnR( istart1- (i2-1)*Num_wann:iend1- (i2-1)*Num_wann, &
                 istart2- (i1-1)*Num_wann:iend2- (i1-1)*Num_wann, iR)*ratio/ndegen(iR)* fac
      
                 !> there is soc term in the hr file
                 if (soc>0) then
                    istart1= (i2-1)*Num_wann+ orbital_start(ia1) + Nwann
                    iend1= (i2-1)*Num_wann+ orbital_start(ia1+1)- 1 + Nwann 
                    istart2= (i1-1)*Num_wann+ orbital_start(ia2)
                    iend2= (i1-1)*Num_wann+ orbital_start(ia2+1)- 1
                    
                    H01( istart1:iend1, istart2:iend2) &
                    = H01( istart1:iend1, istart2:iend2) &
                    + HmnR( istart1- (i2-1)*Num_wann:iend1- (i2-1)*Num_wann, &
                    istart2- (i1-1)*Num_wann:iend2- (i1-1)*Num_wann, iR)*ratio/ndegen(iR)* fac
      
                    istart1= (i2-1)*Num_wann+ orbital_start(ia1)
                    iend1= (i2-1)*Num_wann+ orbital_start(ia1+1)- 1 
                    istart2= (i1-1)*Num_wann+ orbital_start(ia2) + Nwann
                    iend2= (i1-1)*Num_wann+ orbital_start(ia2+1)- 1 + Nwann
                    
                    H01( istart1:iend1, istart2:iend2) &
                    = H01( istart1:iend1, istart2:iend2) &
                    + HmnR( istart1- (i2-1)*Num_wann:iend1- (i2-1)*Num_wann, &
                    istart2- (i1-1)*Num_wann:iend2- (i1-1)*Num_wann, iR)*ratio/ndegen(iR)* fac
      
                    istart1= (i2-1)*Num_wann+ orbital_start(ia1) + Nwann
                    iend1= (i2-1)*Num_wann+ orbital_start(ia1+1)- 1 + Nwann 
                    istart2= (i1-1)*Num_wann+ orbital_start(ia2) + Nwann
                    iend2= (i1-1)*Num_wann+ orbital_start(ia2+1)- 1 + Nwann
                    
                    H01( istart1:iend1, istart2:iend2) &
                    = H01( istart1:iend1, istart2:iend2) &
                    + HmnR( istart1- (i2-1)*Num_wann:iend1- (i2-1)*Num_wann, &
                    istart2- (i1-1)*Num_wann:iend2- (i1-1)*Num_wann, iR)*ratio/ndegen(iR)* fac
                 endif ! soc
              enddo ! ia2
              enddo ! ia1
           enddo ! iR
        enddo ! i2
     enddo ! i1


     !> periodic boundary
     ham_landau= H00+ H01* exp(zi*k(3)) + conjg(transpose(H01))* exp(-zi*k(3))

     !> check hermitcity
     do i1=1,Nq*Num_wann
     do i2=1,Nq*Num_wann
       !if(abs(H00(i1,i2)-conjg(H00(i2,i1))).ge.1e-6)then
       !  write(stdout,*)'there are something wrong with ham_landau'
       !  stop
       !endif 
     enddo
     enddo

     return
  end subroutine ham_3Dlandau


  !> calculate hamiltonian for landau levels
  !> don't consider the internal atom's position
  subroutine ham_3Dlandau2(Ndimq, Nq, k, ham_landau)
     use para
     implicit none

     integer, intent(in) :: Ndimq
     integer, intent(in) :: Nq
     real(dp), intent(in) :: k(3)
     complex(dp), intent(out) :: ham_landau(Ndimq, Ndimq)

     !> inta-hopping for the supercell
     complex(dp), allocatable :: H00(:, :)

     !> inter-hopping for the supercell
     complex(dp), allocatable :: H01(:, :)

     ! loop index  
     integer :: i1, i2

     ! loop index
     integer :: iR

     ! index used to sign irvec     
     integer :: ia,ib,ic

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

     real(dp) :: Ri(3)
     real(dp) :: Rj(3)

     allocate( H00( Ndimq, Ndimq))
     allocate( H01( Ndimq, Ndimq))

     nwann= Num_wann/2
     allocate( orbital_start(Num_atoms+ 1))
     orbital_start= 0
     orbital_start(1)= 1
     do ia=1, Num_atoms
        orbital_start(ia+1)= orbital_start(ia)+ nprojs(ia)
     enddo

     !> calculate intra-hopping
     H00=0.0d0 
     ! i1 column index
     do i1=1, Nq
        ! i2 row index
        do i2=1, Nq
           if (abs(i2-i1)> ijmax) cycle
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
      
              Ri= (i1-1)*Ruc_new
              Rj= new_ia*Rua_new+ new_ib*Rub_new+ (i2-1)*Ruc_new

              phase= By*(Rj(3)+Ri(3))*(Rj(1)-Ri(1))  &
                   - Bx*(Rj(3)+Ri(3))*(Rj(2)-Ri(2))
              fac= cos(phase)+ zi*sin(phase)
      
              istart1= (i2-1)*Num_wann+ 1
              iend1= (i2-1)*Num_wann+ Num_wann
              istart2= (i1-1)*Num_wann+ 1
              iend2= (i1-1)*Num_wann+ Num_wann

              H00(istart1:iend1, istart2:iend2)= &
              H00(istart1:iend1, istart2:iend2)+ &
              HmnR(:, :, iR)/dble(ndegen(iR))* ratio*fac

           enddo ! iR
          !pause
        enddo ! i2
     enddo ! i1

 

     !> check hermitcity
     do i1=1,Nq*Num_wann
     do i2=1,Nq*Num_wann
        if(abs(H00(i1,i2)-conjg(H00(i2,i1))).ge.1e-6)then
          write(stdout,*)'there are something wrong with H00'
          write(*,*)'there are something wrong with H00'
          print *, i1, i2, H00(i1, i2)
          stop
        endif 
     enddo
     enddo

     !> calculate inter-hopping
     H01=0.0d0 
     ! i1 column index
     do i1=1, Nq
        ! i2 row index
        do i2=1, Nq
           if (abs(i2+Nq-i1)> ijmax) cycle
           !> sum over R points to get H(k1, k2)
           do iR=1, Nrpts
              ia=irvec(1,iR)
              ib=irvec(2,iR)
              ic=irvec(3,iR)
      
              !> new lattice
              call latticetransform(ia, ib, ic, new_ia, new_ib, new_ic)
      
              inew_ic= int(new_ic)
              if (inew_ic /= (i2+Nq-i1)) cycle
      
              !> exp(i k.R)
              kdotr= k(1)*new_ia+ k(2)*new_ib
              ratio= cos(2d0*pi*kdotr)+ zi*sin(2d0*pi*kdotr)
      
              Ri= (i1-1)*Ruc_new
              Rj= new_ia*Rua_new+ new_ib*Rub_new+ (i2+Nq-1)*Ruc_new

              phase= By*(Rj(3)+Ri(3))*(Rj(1)-Ri(1))  &
                   - Bx*(Rj(3)+Ri(3))*(Rj(2)-Ri(2))
              fac= cos(phase)+ zi*sin(phase)
      
              istart1= (i2-1)*Num_wann+ 1
              iend1= (i2-1)*Num_wann+ Num_wann
              istart2= (i1-1)*Num_wann+ 1
              iend2= (i1-1)*Num_wann+ Num_wann

              H01(istart1:iend1, istart2:iend2)= &
              H01(istart1:iend1, istart2:iend2)+ &
              HmnR(:, :, iR)/dble(ndegen(ir))* ratio*fac

           enddo ! iR
        enddo ! i2
     enddo ! i1


     !> periodic boundary
     ham_landau= H00 + H01* exp(zi*k(3)) + conjg(transpose(H01))* exp(-zi*k(3))

     !> check hermitcity
     do i1=1,Nq*Num_wann
     do i2=1,Nq*Num_wann
       !if(abs(H00(i1,i2)-conjg(H00(i2,i1))).ge.1e-6)then
       !  write(stdout,*)'there are something wrong with ham_landau'
       !  stop
       !endif 
     enddo
     enddo

     return
  end subroutine ham_3Dlandau2

