!
! this subroutine is used to read some paramters from
! input.dat
! constructed on 4/22/2010 by QS.Wu 


  subroutine readinput

     use wmpi
     use para
     implicit none

     character*12 :: fname='input.dat'
     character*25 :: char_temp 
     character*80 :: inline
     logical ::  exists
     logical ::  lfound
     real(dp) :: cell_volume
     real(dp) :: cell_volume2

     integer  :: stat
     integer  :: i, ia, n
     integer  :: j, io
     integer  :: NN
     integer :: nwann
     real(dp) :: t1, temp
     real(dp) :: pos(3)
     real(dp) :: k1(3), k2(3)
     real(dp) :: kstart(3), kend(3)
     real(dp) :: R1(3), R2(3), R3(3) 
     real(dp), external :: norm

    
     inquire(file=fname,exist=exists)
     if (exists)then
        if(cpuid==0)write(stdout,*) '  '
        if(cpuid==0)write(stdout,*) '>>>Read some paramters from input.dat'
        open(unit=1001,file=fname,status='old')
     else
        if(cpuid==0)write(stdout,*)'file' ,fname, 'dosnot exist'
        stop
     endif

    !read(1001,*)Hrfile
     read(1001, TB_FILE, iostat= stat)
     if (stat/=0) then
        Hrfile='wannier90_hr.dat'
        inquire(file='wannier90_hr.dat',exist=exists)
        if (.not.exists) stop "TB_FIlE namelist should be given or wannier90_hr.dat should exist"
     endif
     if(cpuid==0)write(stdout,'(1x, a, a25)')"Tight binding Hamiltonian file: ",Hrfile
    

     BulkBand_calc         = .FALSE.
     BulkFS_calc           = .FALSE.
     BulkGap_cube_calc     = .FALSE.
     BulkGap_plane_calc    = .FALSE.
     SlabBand_calc         = .FALSE.
     WireBand_calc         = .FALSE.
     SlabSS_calc           = .FALSE.
     SlabArc_calc          = .FALSE.
     SlabQPI_calc          = .FALSE.
     SlabSpintexture_calc  = .FALSE.
     wanniercenter_calc    = .FALSE.
     BerryPhase_calc       = .FALSE.
     BerryCurvature_calc   = .FALSE.
     Dos_calc              = .FALSE.
     JDos_calc             = .FALSE.
     EffectiveMass_calc    = .FALSE.
     
     read(1001, CONTROL, iostat=stat)

     if (stat/=0) then
        write(*, *)"ERROR: namelist CONTROL should be set"
        write(*, *)"You should set one of these functions to be T"
        write(*, *)"BulkBand_calc,BulkFS_calc,BulkGap_cube_calc,BulkGap_plane_calc"
        write(*, *)"SlabBand_calc,WireBand_calc,SlabSS_calc,SlabArc_calc "
        write(*, *)"SlabSpintexture,wanniercenter_calc"
        write(*, *)"BerryPhase_calc,BerryCurvature_calc"
        write(*, *)"The default Vaule is F"
        stop
     endif

     !> control parameters
     if (cpuid==0) then
        write(stdout, *) "  "
        write(stdout, *) ">>>Control parameters: " 
        write(stdout, *) "BulkBand_calc       : ",  BulkBand_calc
        write(stdout, *) "BulkFS_calc         : ",  BulkFS_calc
        write(stdout, *) "BulkGap_cube_calc   : ",  BulkGap_cube_calc
        write(stdout, *) "BulkGap_plane_calc  : ",  BulkGap_plane_calc
        write(stdout, *) "SlabBand_calc       : ",  SlabBand_calc
        write(stdout, *) "SlabSS_calc         : ",  SlabSS_calc
        write(stdout, *) "SlabArc_calc        : ",  SlabArc_calc
        write(stdout, *) "SlabSpintexture_calc: ",  SlabSpintexture_calc
        write(stdout, *) "wanniercenter_calc  : ", wanniercenter_calc
        write(stdout, *) "BerryPhase_calc     : ", BerryPhase_calc
        write(stdout, *) "BerryCurvature_calc : ", BerryCurvature_calc
     endif

     !> set system parameters by default
     Nslab= 10
     Nslab1= 1 
     Nslab2= 1 
     Numoccupied = 0
     Ntotch = 0
     SOC = -1
     E_FERMI = 0
     Bx = 0
     By = 0
     Bz = 0
     surf_onsite = 0
     
     !> read system parameters from file
     read(1001, SYSTEM, iostat=stat)
     if (stat/=0) then
        write(*, *)"ERROR: namelist SYSTEM is wrong and should be set correctly"
       !stop
     endif

     if (SOC == -1) then
        write(*, *)"ERROR: you should set SOC in namelist SYSTEM correctly"
        stop
     endif

     if (Numoccupied == 0) then
        write(*, *)"ERROR: you should set Numoccupied in namelist SYSTEM correctly"
        stop
     endif

     if (Ntotch == 0) then
        Ntotch = Numoccupied
     endif

     if (cpuid==0) then
        write(stdout, *) "  "
        write(stdout, *) ">>>System parameters: " 
        write(stdout, '(1x, a, i6 )')"NumSlabs :", Nslab
        write(stdout, '(1x, a, i6)')"Nslab1 for ribbon  :", Nslab1
        write(stdout, '(1x, a, i6)')"Nslab2 for ribbon  :", Nslab2
        write(stdout, '(1x, a, i6)')"Number of Occupied bands:", NumOccupied
        write(stdout, '(1x, a, i6)')"Number of total electrons:", Ntotch
        write(stdout, '(1x, a, i6)')"With SOC or not in Hrfile:", SOC
        write(stdout, '(1x, a, 3f16.6)')"Fermi energy :", E_FERMI
        write(stdout, '(1x, a, 3f16.6)')"Bx, By, Bz :", Bx, By, Bz
        write(stdout, '(1x, a, 3f16.6)')"surf_onsite :", surf_onsite
     endif

     !> set up parameters for calculation
     E_arc = 0.0d0
     Eta_Arc= 0.001d0
     OmegaNum = 100
     OmegaMin = -1d0
     OmegaMax =  1d0
     Nk1 = 50
     Nk2 = 50
     Nk3 = 50
     NP = 2
     Gap_threshold= 0.01d0
    
     read(1001, PARAMETERS, iostat= stat)

     if (cpuid==0) then
        write(stdout, *) "  "
        write(stdout, *) ">>>calculation parameters : " 
        write(stdout, '(1x, a, f16.5)')'E_arc : ', E_arc
        write(stdout, '(1x, a, f16.5)')'Eta_arc : ', Eta_arc
        write(stdout, '(1x, a, f16.5)')'Gap_threshold', Gap_threshold
        write(stdout, '(1x, a, f16.5)')'OmegaMin : ', OmegaMin
        write(stdout, '(1x, a, f16.5)')'OmegaMax : ', OmegaMax
        write(stdout, '(1x, a, i6   )')'OmegaNum : ', OmegaNum
        write(stdout, '(1x, a, i6   )')'Nk1 : ', Nk1
        write(stdout, '(1x, a, i6   )')'Nk2 : ', Nk2
        write(stdout, '(1x, a, i6   )')'Nk3 : ', Nk3
        write(stdout, '(1x, a, i6   )')'NP number of principle layers  : ', Np
     endif
   
     NK = Nk1

     !> read lattice information
     rewind(1001)
     lfound = .false.
     do while (.true.)
        read(1001, *, end= 100)inline
        if (trim(adjustl(inline))=='LATTICE') then
           lfound= .true.
           if (cpuid==0) write(stdout, *)' '
           if (cpuid==0) write(stdout, *)'We found LATTICE card'
           exit
        endif
     enddo
     100 continue

     if (lfound) then
        read(1001, *)inline   ! The unit of lattice vector
        AngOrBohr=trim(adjustl(inline))
        read(1001, *)Rua
        read(1001, *)Rub
        read(1001, *)Ruc
     else 
        stop 'ERROR: please set lattice information'
     endif

     if (index(AngOrBohr, 'Bohr')>0) then
        Rua= Rua*0.529177d0
        Rub= Rub*0.529177d0
        Ruc= Ruc*0.529177d0
     endif

     !> transform lattice from direct space to reciprocal space

     Kua= 0d0
     Kub= 0d0
     Kuc= 0d0
     cell_volume= Rua(1)*(Rub(2)*Ruc(3)- Rub(3)*Ruc(2)) &
                 +Rua(2)*(Rub(3)*Ruc(1)- Rub(1)*Ruc(3)) &
                 +Rua(3)*(Rub(1)*Ruc(2)- Rub(2)*Ruc(1)) 
     CellVolume= cell_volume
     cell_volume= 2d0*3.1415926535d0/cell_volume

     PrimitiveCellVolume= cell_volume
     Kua(1)= cell_volume*(Rub(2)*Ruc(3)- Rub(3)*Ruc(2))
     Kua(2)= cell_volume*(Rub(3)*Ruc(1)- Rub(1)*Ruc(3))
     Kua(3)= cell_volume*(Rub(1)*Ruc(2)- Rub(2)*Ruc(1))

     Kub(1)= cell_volume*(Ruc(2)*Rua(3)- Ruc(3)*Rua(2))
     Kub(2)= cell_volume*(Ruc(3)*Rua(1)- Ruc(1)*Rua(3))
     Kub(3)= cell_volume*(Ruc(1)*Rua(2)- Ruc(2)*Rua(1))

     Kuc(1)= cell_volume*(Rua(2)*Rub(3)- Rua(3)*Rub(2))
     Kuc(2)= cell_volume*(Rua(3)*Rub(1)- Rua(1)*Rub(3))
     Kuc(3)= cell_volume*(Rua(1)*Rub(2)- Rua(2)*Rub(1))

     if(cpuid==0)write(stdout, '(a)') '>> lattice information (Angstrom)'
     if(cpuid==0)write(stdout, '(3f12.6)')Rua
     if(cpuid==0)write(stdout, '(3f12.6)')Rub
     if(cpuid==0)write(stdout, '(3f12.6)')Ruc

     if(cpuid==0)write(stdout, '(a)') '>> Reciprocal lattice information (1/Angstrom)'
     if(cpuid==0)write(stdout, '(3f12.6)')Kua
     if(cpuid==0)write(stdout, '(3f12.6)')Kub
     if(cpuid==0)write(stdout, '(3f12.6)')Kuc

     !> Read atom positions information
     rewind(1001)
     lfound = .false.
     do while (.true.)
        read(1001, *, end= 101)inline
        if (trim(adjustl(inline))=='ATOM_POSITIONS') then
           lfound= .true.
           if (cpuid==0) write(stdout, *)' '
           if (cpuid==0) write(stdout, *)'We found ATOM_POSITIONS card'
           exit
        endif
     enddo
     101 continue

     if (lfound) then
        read(1001, *)Num_atoms   ! The unit of lattice vector
        if(cpuid==0)write(stdout, '(a, i)')'Num_atoms', Num_atoms
        allocate(atom_name(Num_atoms))
        allocate(Atom_position(3, Num_atoms))
        allocate(Atom_position_direct(3, Num_atoms))
        read(1001, *)inline   ! The unit of lattice vector
        DirectOrCart= trim(adjustl(inline))
        do i=1, Num_atoms
           read(1001, *) atom_name(i), Atom_position(:, i)
           if(cpuid==0)write(stdout, '(a4,3f12.6)')atom_name(i), Atom_position(:, i)
           if (index(DirectOrCart, "D")>0)then
              pos= Atom_position(:, i)
              Atom_position(:, i)= pos(1)*Rua+ pos(2)*Rub+ pos(3)*Ruc
           endif
        enddo
        if(cpuid==0)write(stdout,'(a)')'Atom position in cartisen coordinate'
        do i=1, Num_atoms
           if(cpuid==0)write(stdout, '(a4,3f12.6)')atom_name(i), Atom_position(:, i)
        enddo
     
        if(cpuid==0)write(stdout,'(a)')'Atom position in direct coordinate'
        do ia=1, Num_atoms
           call cart_direct_real(Atom_position(:, ia), Atom_position_direct(:, ia))
           if(cpuid==0)write(stdout, '(a4,3f12.6)')atom_name(ia), Atom_position_direct(:, ia)
        enddo
     else
        stop "ERROR: please set atom's positions information"
     endif

     !> Read projectors information
     rewind(1001)
     lfound = .false.
     do while (.true.)
        read(1001, *, end= 102)inline
        if (trim(adjustl(inline))=='PROJECTORS'.or.&
            trim(adjustl(inline))=='PROJECTOR') then
           lfound= .true.
           if (cpuid==0) write(stdout, *)' '
           if (cpuid==0) write(stdout, *)'We found PROJECTORS card'
           exit
        endif
     enddo
     102 continue

     if (lfound) then
        allocate(nprojs(Num_atoms))
        nprojs= 0
        read(1001, *)nprojs
        if(cpuid==0)write(stdout, '(a, 40i5)')'nprojs', nprojs

        max_projs= maxval(nprojs)
        allocate(proj_name(max_projs, Num_atoms))
        proj_name= ' '
        do i=1, Num_atoms
           read(1001, *)char_temp, proj_name(1:nprojs(i), i)
           if(cpuid==0)write(stdout, '(40a8)') &
              char_temp, proj_name(1:nprojs(i), i)
        enddo
     else
        stop "ERROR: please set projectors for Wannier functions information"
     endif

     !> set up orbitals_start
     allocate(orbitals_start(Num_atoms))
     orbitals_start= 1
     do i=1, Num_atoms-1
        orbitals_start(i+1)= orbitals_start(i)+ nprojs(i)
     enddo


     !> read Wannier centres
     nwann= sum(nprojs)
     if (SOC>0) nwann= 2*nwann
     allocate(wannier_centers_cart(3, Nwann))
     allocate(wannier_centers_direct(3, Nwann))
     wannier_centers_direct= 0d0
     wannier_centers_cart= 0d0
     !> default wannier centers 
     i= 0   
     do ia= 1, Num_atoms
        do j= 1, nprojs(ia)
           i= i+ 1
           wannier_centers_cart(:, i)=  Atom_position(:, ia)
           call cart_direct_real(wannier_centers_cart(:, i),  &
              wannier_centers_direct(:, i))
           if (SOC>0) then
              wannier_centers_cart(:, i+Nwann/2)=  Atom_position(:, ia)
              call cart_direct_real(wannier_centers_cart(:, i+Nwann/2),  &
                 wannier_centers_direct(:, i+Nwann/2))
           endif
        enddo ! j
     enddo ! ia

     rewind(1001)
     lfound = .false.
     do while (.true.)
        read(1001, *, end= 110)inline
        if (trim(adjustl(inline))=='WANNIER_CENTERS' &
           .or. trim(adjustl(inline))=='WANNIER_CENTRES') then
           lfound= .true.
           if (cpuid==0) write(stdout, *)' '
           if (cpuid==0) write(stdout, *)'We found wannier_centers card'
           exit
        endif
     enddo

     if (lfound) then
        read(1001, *)inline   ! The unit of lattice vector
        DirectOrCart= trim(adjustl(inline))

        if (index(DirectOrCart, "D")>0)then
           do i=1, Nwann
              read(1001, *) wannier_centers_direct(:, i)
              call direct_cart_real(wannier_centers_direct(:, i), &
                 wannier_centers_cart(:, i))
           enddo

        else
           do i=1, Nwann
              read(1001, *) wannier_centers_cart(:, i)
              call cart_direct_real(wannier_centers_cart(:, i), &
                 wannier_centers_direct(:, i))
           enddo
        endif
     endif ! found wannier_centers card
 
     110 continue

     if (lfound) then
        if (cpuid==0) then
           write(stdout, *)" "
           write(stdout, *)">> Wannier centers from input.dat, in unit of reciprocal lattice vector"
           write(stdout, '(a6, 4a10)')'iwann', 'R1', 'R2', 'R3'
           do i=1, Nwann
              write(stdout, '(i6, 3f10.6)')i, wannier_centers_direct(:, i)
             !write(stdout, '(i6, 3f10.6)')i, wannier_centers_cart(:, i)
           enddo
        endif
     else
        if (cpuid==0) then
           write(stdout, *)" "
           write(stdout, *)">> Wannier centers by default, in unit of reciprocal lattice vector"
           write(stdout, '(a6, 4a10)')'iwann', 'R1', 'R2', 'R3'
           do i=1, Nwann
              write(stdout, '(i6, 3f10.6)')i, wannier_centers_direct(:, i)
             !write(stdout, '(i6, 3f10.6)')i, wannier_centers_cart(:, i)
           enddo
        endif
     endif


     !> read surface information
     rewind(1001)
     lfound = .false.
     do while (.true.)
        read(1001, *, end= 103)inline
        if (trim(adjustl(inline))=='SURFACE') then
           lfound= .true.
           if (cpuid==0) write(stdout, *)' '
           if (cpuid==0) write(stdout, *)'We found SURFACE card'
           exit
        endif
     enddo
     103 continue

     if (.not.lfound) then
        print *, inline
        stop 'ERROR: please set surface information'
     endif

     !> read information for new lattice 
     !> in order to get different surface state
     !> R1'=U11*R1+U12*R2+U13*R3
     !> R2'=U21*R1+U22*R2+U23*R3
     !> R3'=U31*R1+U32*R2+U33*R3
     read(1001, *)Umatrix(1, :)
     read(1001, *)Umatrix(2, :)
     read(1001, *)Umatrix(3, :)

     if (cpuid==0) then
        write(stdout, '(a)')'>> new vectors to define the surface (in unit of lattice vector)' 
        write(stdout, '(a, 3f12.6)')' The 1st vector on surface  :', Umatrix(1, :)
        write(stdout, '(a, 3f12.6)')' The 2nd vector on surface  :', Umatrix(2, :)
        write(stdout, '(a, 3f12.6)')' The 3rd vector out surface :', Umatrix(3, :)
     endif

     !> check whether Umatrix is right
     !> the volume of the new cell should be the same as the old ones
     R1= Umatrix(1, 1)*Rua+ Umatrix(1, 2)*Rub+ Umatrix(1, 3)*Ruc
     R2= Umatrix(2, 1)*Rua+ Umatrix(2, 2)*Rub+ Umatrix(2, 3)*Ruc
     R3= Umatrix(3, 1)*Rua+ Umatrix(3, 2)*Rub+ Umatrix(3, 3)*Ruc


     cell_volume2= R1(1)*(R2(2)*R3(3)- R2(3)*R3(2)) &
                 +R1(2)*(R2(3)*R3(1)- R2(1)*R3(3)) &
                 +R1(3)*(R2(1)*R3(2)- R2(2)*R3(1)) 
     cell_volume2= 2d0*3.1415926535d0/cell_volume2

     if (cell_volume2<0) then
        R3=-R3
        Umatrix(3, :)= -Umatrix(3, :)
     endif

     if (abs(abs(cell_volume2)-abs(cell_volume))> 0.001d0.and.cpuid==0) then
        write(stdout, *)' ERROR The Umatrix is wrong, the new cell', &
           'volume should be the same as the old ones'
        write(stdout, '(a,2f10.4)')'cell_volume vs cell_volume-new', cell_volume, cell_volume2
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

     call rotate(R1, Rua_new)
     call rotate(R2, Rub_new)
     call rotate(R3, Ruc_new)

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

     if (abs(cell_volume)<1e-6) stop 'cell_volume equal zero'

     Ka2(1)= 2d0*pi/cell_volume*Rb2(2)
     Ka2(2)=-2d0*pi/cell_volume*Rb2(1)
     Kb2(1)=-2d0*pi/cell_volume*Ra2(2)
     Kb2(2)= 2d0*pi/cell_volume*Ra2(1)

     if (cpuid==0) then
        write(stdout, *)'Cell_Volume: ', Cell_Volume
        write(stdout, *)'Ra2, Rb2'
        write(stdout, '(3f10.4)')Ra2
        write(stdout, '(3f10.4)')Rb2
        write(stdout, *)'Ka2, Kb2'
        write(stdout, '(3f10.4)')ka2
        write(stdout, '(3f10.4)')kb2
     endif


     !> read kpath_bulk information
     rewind(1001)
     lfound = .false.
     do while (.true.)
        read(1001, *, end= 104)inline
        if (trim(adjustl(inline))=='KPATH_BULK') then
           lfound= .true.
           if (cpuid==0) write(stdout, *)' '
           if (cpuid==0) write(stdout, *)'We found KPATH_BULK card'
           exit
        endif
     enddo

     !> kline for 3d band structure
     !> high symmetry k points
     read(1001, *) nk3lines
     if(cpuid==0)write(stdout, '(a, 40i5)')'Number of K lines : ', nk3lines
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
        if(cpuid==0)write(stdout, '(a5, 3f9.4, 2x, a5, 3f9.4)')&
          k3line_name(i), k3line_start(:, i), &
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
          !k1= kstart
          !k2= kend

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

     104 continue
     if (.not.lfound .and. BulkBand_calc == .TRUE.) then
        stop 'ERROR: please set KPATH_BULK for bulk band structure calculation'
     endif

     !> read kpath_slab information
     rewind(1001)
     lfound = .false.
     do while (.true.)
        read(1001, *, end= 105)inline
        if (trim(adjustl(inline))=='KPATH_SLAB') then
           lfound= .true.
           if (cpuid==0) write(stdout, *)' '
           if (cpuid==0) write(stdout, *)'We found KPATH_SLAB card'
           exit
        endif
     enddo

     !> read in k lines for 2D system
     k2line_name= ' '
     if (cpuid==0) write(stdout, *)'k lines for 2D system'
     read(1001, *)nk2lines
     if (cpuid==0) write(stdout, *)'No. of k lines for 2D system : ', nk2lines
     do i=1, nk2lines
        read(1001, *) k2line_name(i), kp(:, i), &
                      char_temp, ke(:, i)
        if(cpuid==0) write(stdout, '(a6, 2f9.5, 4x, a6, 2f9.5)')&
                      k2line_name(i), kp(:, i), &
                      char_temp, ke(:, i)
     enddo
     k2line_name(nk2lines+1) = char_temp
 
     NN= Nk
     knv2= NN*nk2lines
     allocate( k2_path(knv2, 2))
     allocate( k2len (knv2))
     k2_path= 0d0
     k2len= 0d0

     t1=0d0
     k2len=0d0
     k2line_stop= 0d0
     do j=1, nk2lines 
        do i=1, NN
           kstart(1:2)= kp(:, j)
           kend(1:2)  = ke(:, j)
           k1(1:2)= kstart(1)*Ka2+ kstart(2)*Kb2
           k2(1:2)= kend(1)*Ka2+ kend(2)*Kb2
           k2_path(i+(j-1)*NN,:)= kstart(1:2)+ (kend(1:2)-kstart(1:2))*(i-1)/dble(NN-1)
           
           temp= dsqrt((k2(1)- k1(1))**2 &
                 +(k2(2)- k1(2))**2)/dble(NN-1) 

           if (i.gt.1) then
              t1=t1+temp
           endif
           k2len(i+(j-1)*NN)= t1
        enddo
        k2line_stop(j+1)= t1

     enddo


     105 continue
     if (.not.lfound .and.(SlabBand_calc == .TRUE. .or. SlabSS_calc==.TRUE.)) then
        stop 'ERROR: please set KPATH_SLAB for slab band structure calculation'
     endif


     !> read kplane_slab information
     !> default value for KPLANE_SLAB
     K2D_start= (/-0.5, -0.5/)
     K2D_vec1 = (/ 1.0,  0.0/)
     K2D_vec2 = (/ 0.0,  1.0/)

     rewind(1001)
     lfound = .false.
     do while (.true.)
        read(1001, *, end= 106)inline
        if (trim(adjustl(inline))=='KPLANE_SLAB') then
           lfound= .true.
           if (cpuid==0) write(stdout, *)' '
           if (cpuid==0) write(stdout, *)'We found KPLANE_SLAB card'
           exit
        endif
     enddo

     !> kpoints plane for 2D system--> arcs  
     read(1001, *)K2D_start
     read(1001, *)K2D_vec1
     read(1001, *)K2D_vec2
     106 continue

     if (cpuid==0) write(stdout, *)'>> Kpoints plane for 2D system--> arcs  '
     if (cpuid==0) write(stdout, '((a, 2f8.4))')'K2D_start:', K2D_start
     if (cpuid==0) write(stdout, '((a, 2f8.4))')'The first vector: ', K2D_vec1
     if (cpuid==0) write(stdout, '((a, 2f8.4))')'The second vector: ', K2D_vec2
     if (.not.lfound .and.(SlabArc_calc == .TRUE. .or. SlabSpintexture_calc==.TRUE.)) then
        stop 'ERROR: please set KPLANE_SLAB for arc or spintexture calculations'
     endif



     !> read kplane_bulk information
     !> default value for KPLANE_BULK
     K3D_start= (/ 0.0,  0.0,   0.0/)
     K3D_vec1 = (/ 1.0,  0.0,   0.0/)
     K3D_vec2 = (/ 0.0,  0.5,   0.0/)

     rewind(1001)
     lfound = .false.
     do while (.true.)
        read(1001, *, end= 107)inline
        if (trim(adjustl(inline))=='KPLANE_BULK') then
           lfound= .true.
           if (cpuid==0) write(stdout, *)' '
           if (cpuid==0) write(stdout, *)'We found KPLANE_BULK card'
           exit
        endif
     enddo

     !> kpoints plane for 3D system--> gapshape
     read(1001, *)K3D_start
     read(1001, *)K3D_vec1
     read(1001, *)K3D_vec2
     107 continue

     if (cpuid==0) write(stdout, *)'>> Kpoints plane for 3D system--> gapshape  '
     if (cpuid==0) write(stdout, '((a, 3f8.4))')'k3D_start : ', K3D_start
     if (cpuid==0) write(stdout, '((a, 3f8.4))')'The 1st vector: ', K3D_vec1
     if (cpuid==0) write(stdout, '((a, 3f8.4))')'The 2nd vector: ', K3D_vec2
     if (.not.lfound .and.(BulkGap_plane_calc == .TRUE. .or. wanniercenter_calc==.TRUE.)) then
        stop 'ERROR: please set KPLANE_bulk for gap or WCC calculations'
     endif



     !> read kcube_bulk information
     !> default value for KCUBE_BULK
     K3D_start_cube= (/-0.5, -0.5,  -0.5/)
     K3D_vec1_cube = (/ 1.0,  0.0,   0.0/)
     K3D_vec2_cube = (/ 0.0,  1.0,   0.0/)
     K3D_vec3_cube = (/ 0.0,  0.0,   1.0/)

     rewind(1001)
     lfound = .false.
     do while (.true.)
        read(1001, *, end= 108)inline
        if (trim(adjustl(inline))=='KCUBE_BULK') then
           lfound= .true.
           if (cpuid==0) write(stdout, *)' '
           if (cpuid==0) write(stdout, *)'We found KCUBE_BULK card'
           exit
        endif
     enddo

     !> kpoints plane for 3D system--> gapshape
     read(1001, *)K3D_start_cube
     read(1001, *)K3D_vec1_cube
     read(1001, *)K3D_vec2_cube
     read(1001, *)K3D_vec3_cube

     kCubeVolume= K3D_vec1_cube(1)*(K3D_vec2_cube(2)*K3D_vec3_cube(3) &
                  - K3D_vec2_cube(3)*K3D_vec3_cube(2)) &
                 + K3D_vec1_cube(2)*(K3D_vec2_cube(3)*K3D_vec3_cube(1) &
                  - K3D_vec2_cube(1)*K3D_vec3_cube(3)) &
                 + K3D_vec1_cube(3)*(K3D_vec2_cube(1)*K3D_vec3_cube(2) &
                  - K3D_vec2_cube(2)*K3D_vec3_cube(1)) 

     kCubeVolume= kCubeVolume*PrimitiveCellVolume 

     108 continue
     if (cpuid==0) write(stdout, *)'>> Kpoints plane for 3D system--> gapshape3D  '
     if (cpuid==0) write(stdout, '((a, 3f8.4))')'k3D_start :', K3D_start_cube
     if (cpuid==0) write(stdout, '((a, 3f8.4))')'The 1st vector: ', K3D_vec1_cube
     if (cpuid==0) write(stdout, '((a, 3f8.4))')'The 2nd vector: ', K3D_vec2_cube
     if (cpuid==0) write(stdout, '((a, 3f8.4))')'The 3rd vector: ', K3D_vec3_cube
     if (.not.lfound .and.(BulkGap_cube_calc == .TRUE.)) then
        stop 'ERROR: please set KCUBE_BULK for gap3D calculations'
     endif

     !> default parameters for effective mass calculation
     dk_mass= 0.01  ! in unit of 1/Ang
     iband_mass= NumOccupied
     k_mass= 0 
     rewind(1001)
     lfound = .false.
     do while (.true.)
        read(1001, *, end= 109)inline
        if (trim(adjustl(inline))=='EFFECTIVE_MASS') then
           lfound= .true.
           if (cpuid==0) write(stdout, *)' '
           if (cpuid==0) write(stdout, *)'We found EFFECTIVE_MASS card'
           exit
        endif
     enddo

     read(1001, *)iband_mass
     read(1001, *)dk_mass
     read(1001, *)k_mass

     109 continue
     if (cpuid==0) write(stdout, *)' '
     if (.not.lfound.and.cpuid==0)write(stdout, *)'>> Using default parameters for effective mass calculation'
     if (cpuid==0) write(stdout, *)'>> Effective mass calculation parameters  '
     if (cpuid==0) write(stdout, '(a, i5, a)')'>> The ', iband_mass, "'th band"
     if (cpuid==0) write(stdout, '(a, f7.4, a)')'>> k step ', dk_mass, " in unit of 1/Angstrom"
     if (cpuid==0) write(stdout, '(a, 3f7.4, a)')'>> k points ', k_mass, " in unit of reciprocal primitive cell"
     k1=k_mass
     call direct_cart_rec(k1, k_mass)
     if (cpuid==0) write(stdout, '(a, 3f7.4, a)')'>> k points ', k_mass, " in unit of 1/Angstrom"
    
     !> setup the atoms on top and bottom surface that used for output the 
     !> surface-state spectrum
     !> by default we output all the atoms' weight
     NtopAtoms   = Num_atoms 
     NbottomAtoms= Num_atoms 
     allocate(TopAtoms(NtopAtoms))
     allocate(BottomAtoms(NbottomAtoms))
     do i=1, NTopAtoms
        TopAtoms(i)= i
     enddo
     do i=1, NBottomAtoms
        BottomAtoms(i)= i
     enddo

     rewind(1001)
     lfound = .false.
     do while (.true.)
        read(1001, *, end= 112)inline
        if (trim(adjustl(inline))=='TOPBOTTOMATOMS') then
           lfound= .true.
           if (cpuid==0) write(stdout, *)' '
           if (cpuid==0) write(stdout, *)'We found TOPBOTTOMATOMS card'
           exit
        endif
     enddo
     read(1001, *)NtopAtoms
     deallocate(TopAtoms)
     allocate(TopAtoms(NtopAtoms))
     read(1001, *)TopAtoms
     read(1001, *)NbottomAtoms
     deallocate(BottomAtoms)
     allocate(BottomAtoms(NbottomAtoms))
     read(1001, *)BottomAtoms

     112 continue
     if (.not.lfound.and.cpuid==0)write(stdout, *)'>> Output all atoms weight for surface state spectrum'
     if (cpuid==0) write(stdout, *)'> NtopAtoms ', NtopAtoms
     if (cpuid==0) write(stdout, '(a,100i3)')'> TopAtoms ', TopAtoms
     if (cpuid==0) write(stdout, '(a,100i3)')'> NbottomAtoms ', NbottomAtoms
     if (cpuid==0) write(stdout, '(a,100i3)')'> BottomAtoms ',BottomAtoms 

     NtopOrbitals=0
     do i=1, NTopAtoms
        NtopOrbitals= NtopOrbitals+ nprojs(TopAtoms(i))
     enddo
     if (SOC>0) NtopOrbitals= NtopOrbitals*2
     allocate(TopOrbitals(NtopOrbitals))
     TopOrbitals= 1

     !> set up top surface orbitals for output the surface spectrum
     io=0
     do i=1, NTopAtoms
        do j=1, nprojs(TopAtoms(i))
           io =io+ 1
           TopOrbitals(io)= orbitals_start(i)+ j- 1
           if (SOC>0)TopOrbitals(io+ NtopOrbitals/2 )= orbitals_start(i)+ j- 1+ Nwann/2
        enddo ! j
     enddo ! i

     NBottomOrbitals=0
     do i=1, NBottomAtoms
        NBottomOrbitals= NBottomOrbitals+ nprojs(BottomAtoms(i))
     enddo
     if (SOC>0) NBottomOrbitals= NBottomOrbitals*2
     allocate(BottomOrbitals(NBottomOrbitals))
     BottomOrbitals= 1

     !> set up Bottom surface orbitals for output the surface spectrum
     io=0
     do i=1, NBottomAtoms
        do j=1, nprojs(BottomAtoms(i))
           io =io+ 1
           BottomOrbitals(io)= orbitals_start(i)+ j- 1
           if (SOC>0)BottomOrbitals(io+ NBottomOrbitals/2)= orbitals_start(i)+ j- 1+ Nwann/2
        enddo ! j
     enddo ! i

     if (cpuid==0) write(stdout, *)'> NtopOrbitals ', NtopOrbitals
     if (cpuid==0) write(stdout, '(a,999i4)')'> TopOrbitals ', TopOrbitals
     if (cpuid==0) write(stdout, '(a,999i4)')'> NBottomOrbitals ', NBottomOrbitals
     if (cpuid==0) write(stdout, '(a,999i4)')'> BottomOrbitals ',BottomOrbitals


     !> close input.dat
     close(1001)

     eta=(omegamax- omegamin)/omeganum*2d0

     if(cpuid==0)write(stdout,*)'<<<Read input.dat file successfully'


     return
  end subroutine readinput

  function norm(R1)
     use para, only : dp

     implicit none
     real(dp), intent(in) :: R1(3)
     real(dp) :: norm1
     real(dp) :: norm

     norm1= R1(1)*R1(1)+ R1(2)*R1(2)+ R1(3)*R1(3)
     norm= sqrt(norm1)

     return
  end function norm


  !> rotate a vector to the new coordinate system
  subroutine rotate(R1, R2)
     use para, only : dp, Urot
     implicit none
     real(dp), intent(in) :: R1(3)
     real(dp), intent(inout) :: R2(3)

     R2(1)= Urot(1, 1)*R1(1)+ Urot(1, 2)*R1(2)+ Urot(1, 3)*R1(3)
     R2(2)= Urot(2, 1)*R1(1)+ Urot(2, 2)*R1(2)+ Urot(2, 3)*R1(3)
     R2(3)= Urot(3, 1)*R1(1)+ Urot(3, 2)*R1(2)+ Urot(3, 3)*R1(3)

     return
  end subroutine rotate


  !> transform from Cartesian coordinates to direct lattice vector basis
   subroutine cart_direct_real(R1, R2)
      use para
      implicit none
      real(dp), intent(in) :: R1(3)
      real(dp), intent(inout) :: R2(3)
      real(dp) :: mata(3, 3)

      mata(1, :)= Rua
      mata(2, :)= Rub
      mata(3, :)= Ruc

      call inv_r(3, mata)
      R2= R1(1)*mata(1, :)+ R1(2)*mata(2, :)+ R1(3)*mata(3, :)

      return
   end subroutine cart_direct_real

  !> transform from direct lattice vector basis to Cartesian coordinates
   subroutine direct_cart_real(R1, R2)
      use para
      implicit none
      real(dp), intent(in) :: R1(3)
      real(dp), intent(inout) :: R2(3)

      R2= R1(1)*Rua+ R1(2)*Rub+ R1(3)*Ruc

      return
   end subroutine direct_cart_real


  !> transform from Cartesian coordinates to reciprocal lattice vector basis
   subroutine cart_direct_rec(k1, k2)
      use para
      implicit none
      real(dp), intent(in) :: k1(3)
      real(dp), intent(inout) :: k2(3)
      real(dp) :: mata(3, 3)

      mata(1, :)= Kua
      mata(2, :)= Kub
      mata(3, :)= Kuc

      call inv_r(3, mata)
      K2= k1(1)*mata(1, :)+ k1(2)*mata(2, :)+ k1(3)*mata(3, :)

      return
   end subroutine cart_direct_rec

   subroutine direct_cart_rec(k1, k2)
      use para
      implicit none
      real(dp), intent(in) :: k1(3)
      real(dp), intent(inout) :: k2(3)

      K2= k1(1)*Kua+ k1(2)*Kub+ k1(3)*Kuc

      return
   end subroutine direct_cart_rec
