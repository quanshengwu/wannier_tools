!
! this subroutine is used to read some paramters from
! input.dat
! constructed on 4/22/2010 by QS.Wu 


  subroutine readinput

     use para
     implicit none

     character*12 :: fname='input.dat'
     character*25 :: char_temp 
     logical ::  exists
     real(dp) :: cell_volume
     real(dp) :: cell_volume2

     integer  :: i
     integer  :: j
     integer  :: NN
     real(dp) :: t1, temp
     real(dp) :: pos(3)
     real(dp) :: k1(3), k2(3)
     real(dp) :: kstart(3), kend(3)
     real(dp) :: R1(3), R2(3), R3(3) 
     real(dp) :: Urot(3, 3)
     real(dp), external :: norm
    
     inquire(file=fname,exist=exists)
     if (exists)then
        if(cpuid==0)write(stdout,*) 'read some paramters from input.dat'
        open(unit=1001,file=fname,status='old')
     else
        if(cpuid==0)write(stdout,*)'file' ,fname, 'dosnot exist'
        stop
     endif
 

     read(1001,*)infilename
     if(cpuid==0)write(stdout,'(a,a25)')' input file:',infilename
     read(1001,*)outfilename
     if(cpuid==0)write(stdout,'(2a)')' output file:',outfilename
     if(cpuid==0)open(unit=stdout, file=outfilename)
     read(1001,*) BulkBand_calc
     read(1001,*) SlabBand_calc
     read(1001,*) WireBand_calc
     read(1001,*) SlabSS_calc
     read(1001,*) SlabArc_calc
     read(1001,*) SlabSpintexture_calc
     read(1001,*) wanniercenter_calc
     read(1001,*)Nk
     if(cpuid==0)write(stdout,*)'Nk',Nk
     read(1001,*)omeganum
     if(cpuid==0)write(stdout,*)'omeganum',omeganum
     read(1001,*)omegamin, omegamax
     if(cpuid==0)write(stdout,*)'omegamin, omegamax', omegamin, omegamax
     read(1001,*)E_arc
     if(cpuid==0)write(stdout,*)'E_arc', E_arc
     read(1001,*)nslab
     if(cpuid==0)write(stdout,*)'nslab',nslab
     read(1001,*)Np
     if(cpuid==0)write(stdout,*)'Np',Np
     read(1001,*)Numoccupied
     if(cpuid==0)write(stdout,*)'Numoccupied', Numoccupied
     read(1001,*)Soc
     if(cpuid==0)write(stdout,*)'soc',Soc
     read(1001,*)eta_arc
     if(cpuid==0)write(stdout,*)'eta_arc',eta_arc
     read(1001,*)E_fermi
     if(cpuid==0)write(stdout,*)'E_fermi',E_fermi

     Nslab1= Nslab
     Nslab2= Np

     !> lattice information
     read(1001, *)Rua
     read(1001, *)Rub
     read(1001, *)Ruc

     !> transform lattice from direct space to reciprocal space

     Kua= 0d0
     Kub= 0d0
     Kuc= 0d0
     cell_volume= Rua(1)*(Rub(2)*Ruc(3)- Rub(3)*Ruc(2)) &
                 +Rua(2)*(Rub(3)*Ruc(1)- Rub(1)*Ruc(3)) &
                 +Rua(3)*(Rub(1)*Ruc(2)- Rub(2)*Ruc(1)) 
     cell_volume= 2d0*3.1415926535d0/cell_volume
     Kua(1)= cell_volume*(Rub(2)*Ruc(3)- Rub(3)*Ruc(2))
     Kua(2)= cell_volume*(Rub(3)*Ruc(1)- Rub(1)*Ruc(3))
     Kua(3)= cell_volume*(Rub(1)*Ruc(2)- Rub(2)*Ruc(1))

     Kub(1)= cell_volume*(Ruc(2)*Rua(3)- Ruc(3)*Rua(2))
     Kub(2)= cell_volume*(Ruc(3)*Rua(1)- Ruc(1)*Rua(3))
     Kub(3)= cell_volume*(Ruc(1)*Rua(2)- Ruc(2)*Rua(1))

     Kuc(1)= cell_volume*(Rua(2)*Rub(3)- Rua(3)*Rub(2))
     Kuc(2)= cell_volume*(Rua(3)*Rub(1)- Rua(1)*Rub(3))
     Kuc(3)= cell_volume*(Rua(1)*Rub(2)- Rua(2)*Rub(1))

     if(cpuid==0)write(stdout, '(a)') '>> lattice information'
     if(cpuid==0)write(stdout, '(3f10.6)')Rua
     if(cpuid==0)write(stdout, '(3f10.6)')Rub
     if(cpuid==0)write(stdout, '(3f10.6)')Ruc

     if(cpuid==0)write(stdout, '(a)') '>> Reciprocal lattice information'
     if(cpuid==0)write(stdout, '(3f10.6)')Kua
     if(cpuid==0)write(stdout, '(3f10.6)')Kub
     if(cpuid==0)write(stdout, '(3f10.6)')Kuc

     !> read atom position
     read(1001, *)Num_atoms
     allocate(atom_name(Num_atoms))
     allocate(Atom_position(3, Num_atoms))
     do i=1, Num_atoms
        read(1001, *) atom_name(i), Atom_position(:, i)
        if(cpuid==0)write(stdout, '(a4,3f6.3)')atom_name(i), Atom_position(:, i)
        pos= Atom_position(:, i)
        Atom_position(:, i)= pos(1)*Rua+ pos(2)*Rub+ pos(3)*Ruc
     enddo
     if(cpuid==0)write(stdout,'(a)')'Atom position in cartisen coordinate'
     do i=1, Num_atoms
        if(cpuid==0)write(stdout, '(a4,3f6.3)')atom_name(i), Atom_position(:, i)
     enddo


     !> read projectors
     allocate(nprojs(Num_atoms))
     nprojs= 0
     read(1001, *)nprojs

     max_projs= maxval(nprojs)
     allocate(proj_name(max_projs, Num_atoms))
     proj_name= ' '
     do i=1, Num_atoms
        read(1001, *)char_temp, proj_name(1:nprojs(i), i)
     enddo




     !> kline for 3d band structure
     !> high symmetry k points
     read(1001, *) nk3lines
     allocate(k3line_start(3, nk3lines))
     allocate(k3line_end(3, nk3lines))
     allocate(k3line_name(nk3lines+1))
     allocate(k3line_stop(nk3lines+1))
     k3line_stop= 0d0
     k3line_start= 0d0
     k3line_end= 0d0
     k3line_name= ' '
     do i=1, nk3lines
        read(1001, *) k3line_name(i), k3line_start(:, i), &
                      char_temp, k3line_end(:, i)
     enddo
     k3line_name(nk3lines+1)= char_temp

     NN= Nk
     nk3_band= NN*nk3lines
     allocate(k3len(nk3_band))
     allocate(k3points(3, nk3_band))
     k3len=0d0
     k3points= 0d0
     t1= 0d0
     do j=1, nk3lines
        do i=1, NN
           kstart= k3line_start(:, j)
           kend  = k3line_end(:, j)
           k1= kstart(1)*Kua+ kstart(2)*Kub+ kstart(3)*Kuc
           k2= kend(1)*Kua+ kend(2)*Kub+ kend(3)*Kuc

           k3points(:, i+ (j-1)*NN)= kstart+ (kend- kstart)*dble(i-1)/dble(NN-1)
           
           temp= dsqrt((k2(1)- k1(1))**2 &
                 +(k2(2)- k1(2))**2  &
                 +(k2(3)- k1(3))**2)/dble(NN-1) 

           if (i.gt.1) then
              t1=t1+temp
           endif
           k3len(i+(j-1)*NN)= t1
        enddo
        k3line_stop(j+1)= t1
     enddo


     !> read information for new lattice 
     !> in order to get different surface state
     !> R1'=U11*R1+U12*R2+U13*R3
     !> R2'=U21*R1+U22*R2+U23*R3
     !> R3'=U31*R1+U32*R2+U33*R3
     read(1001, *)Umatrix(1, :)
     read(1001, *)Umatrix(2, :)
     read(1001, *)Umatrix(3, :)

     !> check whether Umatrix is right
     !> the volume of the new cell should be the same as the old ones
     R1= Umatrix(1, 1)*Rua+ Umatrix(1, 2)*Rub+ Umatrix(1, 3)*Ruc
     R2= Umatrix(2, 1)*Rua+ Umatrix(2, 2)*Rub+ Umatrix(2, 3)*Ruc
     R3= Umatrix(3, 1)*Rua+ Umatrix(3, 2)*Rub+ Umatrix(3, 3)*Ruc


     cell_volume2= R1(1)*(R2(2)*R3(3)- R2(3)*R3(2)) &
                 +R1(2)*(R2(3)*R3(1)- R2(1)*R3(3)) &
                 +R1(3)*(R2(1)*R3(2)- R2(2)*R3(1)) 
     cell_volume2= 2d0*3.1415926535d0/cell_volume2

     if (abs(abs(cell_volume2)-abs(cell_volume))> 0.001d0.and.cpuid==0) then
        write(stdout, *)' ERROR The Umatrix is wrong, the new cell', &
           'volume should be the same as the old ones'
        stop
     endif

     !> get the surface vector, we should set the new coordinate system
     !> set R1 to the new x direction ex'
     !> set R1\cross R2 to the new z direction ez'
     !> set ey'= ez'\cross ex'
     !> then e_i'= \sum_j U_ij e_j
     Urot= 0d0
     !> e_x'
     Urot(1, :)= R1/norm(R1)

     !> e_z'
     Urot(3, 1)= (R1(2)*R2(3)- R1(3)*R2(2))
     Urot(3, 2)= (R1(3)*R2(1)- R1(1)*R2(3))
     Urot(3, 3)= (R1(1)*R2(2)- R1(2)*R2(1))
     Urot(3, :)= Urot(3, :)/norm(Urot(3, :))

     !> e_y'= e_z'\cross e_x'
     Urot(2, 1)= (Urot(3, 2)*Urot(1, 3)- Urot(3, 3)*Urot(1, 2))
     Urot(2, 2)= (Urot(3, 3)*Urot(1, 1)- Urot(3, 1)*Urot(1, 3))
     Urot(2, 3)= (Urot(3, 1)*Urot(1, 2)- Urot(3, 2)*Urot(1, 1))
     Urot(2, :)= Urot(2, :)/norm(Urot(2, :))

     !> then transform R1, R2 to the new coordinates
     !> R1'_j= \sum_i U_ij R_i
     !> because the z direction is perpendicular to R1, R2, 
     !> so the z coordinates for R1, R2 in the new axis are zero
     Ra2(1)= Urot(1, 1)*R1(1)+ Urot(1, 2)*R1(2)+ Urot(1, 3)*R1(3)
     Ra2(2)= Urot(2, 1)*R1(1)+ Urot(2, 2)*R1(2)+ Urot(2, 3)*R1(3)
     Rb2(1)= Urot(1, 1)*R2(1)+ Urot(1, 2)*R2(2)+ Urot(1, 3)*R2(3)
     Rb2(2)= Urot(2, 1)*R2(1)+ Urot(2, 2)*R2(2)+ Urot(2, 3)*R2(3)

     !> get the surface reciprocal vector
     cell_volume=Ra2(1)*Rb2(2)- Rb2(1)*Ra2(2)
     cell_volume= abs(cell_volume)
     Ka2(1)= 2d0*pi/cell_volume*Rb2(2)
     Ka2(2)=-2d0*pi/cell_volume*Rb2(1)
     Kb2(1)=-2d0*pi/cell_volume*Ra2(2)
     Kb2(2)= 2d0*pi/cell_volume*Ra2(1)



     close(1001)

     eta=(omegamax- omegamin)/omeganum*eta

     if(cpuid==0)write(stdout,*)'read input.dat file successfully'

     return
  end subroutine readinput

  function norm(R1)
     use para, only : dp

     implicit none
     real(dp), intent(in) :: R1(3)
     real(dp) :: norm

     norm= sqrt(R1(1)*R1(1)+ R1(2)*R1(2)+ R1(3)*R1(3))

     return
  end function norm


