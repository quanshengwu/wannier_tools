!
! this subroutine is used to read some paramters from
! input.dat
! constructed on 4/22/2010 by QS.Wu 


  subroutine readinput

     use wmpi
     use para
     implicit none

     character*12 :: fname='wt.in'
     character*25 :: char_temp 
     character*80 :: inline
     logical ::  exists
     logical ::  lfound
     real(dp) :: cell_volume
     real(dp) :: cell_volume2

     integer  :: stat
     integer  :: i, ia, n, ik
     integer  :: j, io, it
     integer  :: NN
     integer :: nwann
     real(dp) :: t1, temp
     real(dp) :: pos(3)
     real(dp) :: k1(3), k2(3), k(3)
     real(dp) :: kstart(3), kend(3)
     real(dp) :: R1(3), R2(3), R3(3), Rt(3)
     real(dp), external :: norm

     real(dp), allocatable :: mass_temp(:)
     real(dp), allocatable :: Born_Charge_temp(:, :, :)


    
     inquire(file=fname,exist=exists)
     if (exists)then
        if(cpuid==0)write(stdout,*) '  '
        if(cpuid==0)write(stdout,*) '>>>Read some paramters from wt.in'
        open(unit=1001,file=fname,status='old')
     else
        if(cpuid==0)write(stdout,*)'file  ' ,fname, 'dosnot exist'
        stop
     endif

     Particle='electron'
     Package= 'VASP'
     KPorTB = 'TB'
     read(1001, TB_FILE, iostat= stat)
     if (stat/=0) then
        Hrfile='wannier90_hr.dat'
        Particle='electron'
        inquire(file='wannier90_hr.dat',exist=exists)
        if (.not.exists) stop "TB_FIlE namelist should be given or wannier90_hr.dat should exist"
     endif
     if(cpuid==0)write(stdout,'(1x, a, a6, a)')"You are using : ", KPorTB, " model"
     if(cpuid==0)write(stdout,'(1x, a, a25)')"Tight binding Hamiltonian file: ",Hrfile
     if(cpuid==0)write(stdout,'(1x, a, a25)')"System of particle: ", Particle
     if(cpuid==0)write(stdout,'(1x, a, a25)')"Tight binding Hamiltonian obtained from : ",Package

     if (index(Particle, 'electron')==0 .and. index(Particle, 'phonon')==0 &
        .and. index(Particle, 'photon')==0) then
        write(stdout, *)' ERROR: Particle shoule equal either "electron", &
           "phonon", or "photon"' 
        stop
     endif
    

     BulkBand_calc         = .FALSE.
     BulkBand_points_calc  = .FALSE.
     BulkBand_plane_calc   = .FALSE.
     BulkFS_calc           = .FALSE.
     BulkGap_cube_calc     = .FALSE.
     BulkGap_plane_calc    = .FALSE.
     BulkFS_Plane_calc     = .FALSE.
     SlabBand_calc         = .FALSE.
     WireBand_calc         = .FALSE.
     SlabSS_calc           = .FALSE.
     SlabArc_calc          = .FALSE.
     SlabQPI_calc          = .FALSE.
     SlabSpintexture_calc  = .FALSE.
     WannierCenter_calc    = .FALSE.
     WeylChirality_calc    = .FALSE.
     Z2_3D_calc            = .FALSE.
     Chern_3D_calc         = .FALSE.
     BerryPhase_calc       = .FALSE.
     BerryCurvature_calc   = .FALSE.
     Dos_calc              = .FALSE.
     JDos_calc             = .FALSE.
     EffectiveMass_calc    = .FALSE.
     FindNodes_calc        = .FALSE.
     LOTO_correction       = .FALSE.
     AHC_Calc              = .FALSE.
     Translate_to_WS_calc  = .FALSE.
     
     read(1001, CONTROL, iostat=stat)


     if (stat/=0) then
        backspace(1001)
        read(1001,fmt='(A)') inline
        write(*,'(A)') &
           'Invalid line in namelist: '//trim(inline)
        write(*, *)"ERROR: namelist CONTROL should be set"
        write(*, *)"You should set one of these functions to be T"
        write(*, *)"BulkBand_points_calc, BulkBand_calc,BulkGap_plane_calc"
        write(*, *)"BulkFS_calc,BulkGap_cube_calc"
        write(*, *)"SlabBand_calc,WireBand_calc,SlabSS_calc,SlabArc_calc "
        write(*, *)"SlabSpintexture,wanniercenter_calc"
        write(*, *)"BerryPhase_calc,BerryCurvature_calc, Z2_3D_calc"
        write(*, *)"Dos_calc, JDos_calc, FindNodes_calc, LOTO_correction "
        write(*, *)"BulkFS_plane_calc, BulkBand_plane_calc"
        write(*, *)"Z2_3D_calc"
        write(*, *)"Chern_3D_calc"
        write(*, *)"WeylChirality_calc"
        write(*, *)"Translate_to_WS_calc"
        write(*, *)"The default Vaule is F"
        write(*, *)"Or there should be some unknonw tags set with = "
        stop
     endif

     !> control parameters
     if (cpuid==0) then
        write(stdout, *) "  "
        write(stdout, *) ">>>Control parameters: " 
        write(stdout, *) "BulkBand_calc       : ",  BulkBand_calc
        write(stdout, *) "BulkBand_points_calc : ",  BulkBand_points_calc
        write(stdout, *) "BulkBand_plane_calc : ",  BulkBand_plane_calc
        write(stdout, *) "BulkFS_calc         : ",  BulkFS_calc
        write(stdout, *) "BulkFS_Plane_calc   : ",  BulkFS_Plane_calc
        write(stdout, *) "BulkGap_cube_calc   : ",  BulkGap_cube_calc
        write(stdout, *) "BulkGap_plane_calc  : ",  BulkGap_plane_calc
        write(stdout, *) "SlabBand_calc       : ",  SlabBand_calc
        write(stdout, *) "SlabSS_calc         : ",  SlabSS_calc
        write(stdout, *) "SlabArc_calc        : ",  SlabArc_calc
        write(stdout, *) "SlabSpintexture_calc: ",  SlabSpintexture_calc
        write(stdout, *) "wanniercenter_calc  : ",  wanniercenter_calc
        write(stdout, *) "Z2_3D_calc          : ",  Z2_3D_calc
        write(stdout, *) "Chern_3D_calc       : ",  Chern_3D_calc
        write(stdout, *) "WeylChirality_calc  : ",  WeylChirality_calc
        write(stdout, *) "BerryPhase_calc     : ",  BerryPhase_calc
        write(stdout, *) "Dos_calc            : ",  DOS_calc
        write(stdout, *) "AHC_calc            : ",  AHC_calc
        write(stdout, *) "JDos_calc           : ",  JDOS_calc
        write(stdout, *) "FindNodes_calc      : ",  FindNodes_calc
        write(stdout, *) "LOTO_correction     : ", LOTO_correction
        write(stdout, *) "BerryCurvature_calc : ",  BerryCurvature_calc
        write(stdout, *) "EffectiveMass_calc  : ",  EffectiveMass_calc 
        write(stdout, *) "Translate_to_WS_calc: ",  Translate_to_WS_calc
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
     bondlength_cutoff= 1000d0
     
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
        write(stdout, '(1x, a, i6)')"Nslab1 for nonawire  :", Nslab1
        write(stdout, '(1x, a, i6)')"Nslab2 for nanowire  :", Nslab2
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
        if(cpuid==0)write(stdout, '(a, i5)')'Num_atoms', Num_atoms
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

     !> orbital index order 
     allocate(index_start(Num_atoms))
     allocate(index_end  (Num_atoms))
     index_start= 0
     index_end= 0
     index_start(1)= 1
     index_end(1)= nprojs(1)
     do i=2, Num_atoms
        index_start(i)= index_start(i-1)+ nprojs(i-1)
        index_end(i)= index_end(i-1)+ nprojs(i)
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

        it= 0
        if (index(DirectOrCart, "D")>0)then
           do i=1, Nwann
              read(1001, *, end=207, err=207) wannier_centers_direct(:, i)
              it= it+ 1
              call direct_cart_real(wannier_centers_direct(:, i), &
                 wannier_centers_cart(:, i))
           enddo

        else
           do i=1, Nwann
              read(1001, *, end=207, err=207) wannier_centers_cart(:, i)
              it= it+ 1
              call cart_direct_real(wannier_centers_cart(:, i), &
                 wannier_centers_direct(:, i))
           enddo
        endif
     endif ! found wannier_centers card
     207 continue
     if (it< Nwann.and.cpuid==0) then
        write(stdout, *)' '
        write(stdout, *)' >>>> Error happens in Wannier_centres card'
        write(stdout, *)' Error: the number of wannier_centers lines should ' 
        write(stdout, *)' equal to the number wannier functions (include spin)'
        write(stdout, *)' Num_wann', Nwann, ' the centres lines you given ', it
        write(stdout, *)' Otherwise, if you do not know the meaning of this,' 
        write(stdout, *)' please delete this card'
        stop
     endif


 
     110 continue

     if (lfound) then
        if (cpuid==0) then
           write(stdout, *)" "
           write(stdout, *)">> Wannier centers from wt.in, in unit of reciprocal lattice vector"
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


     !> read surface information by Miller indices
     rewind(1001)
     lfound = .false.
     do while (.true.)
        read(1001, *, end= 224)inline
        if (trim(adjustl(inline))=='MILLER_INDICES'  &
           .or. trim(adjustl(inline))=='MILLER_INDEX') then
           lfound= .true.
           if (cpuid==0) write(stdout, *)' '
           if (cpuid==0) write(stdout, *)'We found MILLER_INDEX card for slab calculations'
           exit
        endif
     enddo
     224 continue

     MillerIndices= 0
     if (lfound) then
        read(1001, *, err=225, iostat=stat) MillerIndices(:)
        225 continue
        if (stat/=0) stop "Something wrong with setting of MillerIndices, they should be like this 1 0 0"
        if (cpuid.eq.0) then
           write(stdout, '(a, 3i6)')'  Miller indices are :', MillerIndices
        endif
        call MillerIndicestoumatrix()
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

     if (.not.lfound.and.sum(abs(MillerIndices))==0) then
        print *, inline
        write(stdout, *) 'ERROR: please set surface information by setting SURFACE card or', &
             'MILLER_INDEX card'
        write(*, *) 'ERROR: please set surface information by setting SURFACE card or', &
             'MILLER_INDEX card'
        stop
     endif

     if (lfound) then
        !> read information for new lattice 
        !> in order to get different surface state
        !> R1'=U11*R1+U12*R2+U13*R3
        !> R2'=U21*R1+U22*R2+U23*R3
        !> R3'=U31*R1+U32*R2+U33*R3
        read(1001, *)Umatrix(1, :)
        read(1001, *)Umatrix(2, :)

        Umatrix(3,:)=(/0.0,0.0,1.0/)
        read(1001, *, err=260, iostat=stat)Umatrix(3, :)
        260 continue
   
        if (cpuid==0) then
           write(stdout, '(a)')' '
           write(stdout, '(a)')'>> new vectors to define the surface (in unit of lattice vector)' 
           write(stdout, '(a, 3f12.6)')' The 1st vector on surface     :', Umatrix(1, :)
           write(stdout, '(a, 3f12.6)')' The 2nd vector on surface     :', Umatrix(2, :)
           write(stdout, '(a, 3f12.6)')' The 3rd vector out of surface :', Umatrix(3, :)
        endif

        if (sum(abs(MillerIndices))>0 .and. cpuid==0) then
           write(stdout, '(a, a, a)')' Attention: you specified surface information twice.' , &
           ' However, we only take the information from SURFACE card, ', &
           ' and omit the settings by MILLER_INDEX card'
        endif
     endif

     !> check whether Umatrix is right
     !> the volume of the new cell should be the same as the old ones
     !> Here R1, R2, R3 are vectors defined by SURFACE CARD in original cartesian coordinates
     R1= Umatrix(1, 1)*Rua+ Umatrix(1, 2)*Rub+ Umatrix(1, 3)*Ruc
     R2= Umatrix(2, 1)*Rua+ Umatrix(2, 2)*Rub+ Umatrix(2, 3)*Ruc
     R3= Umatrix(3, 1)*Rua+ Umatrix(3, 2)*Rub+ Umatrix(3, 3)*Ruc

     cell_volume2= R1(1)*(R2(2)*R3(3)- R2(3)*R3(2)) &
                 +R1(2)*(R2(3)*R3(1)- R2(1)*R3(3)) &
                 +R1(3)*(R2(1)*R3(2)- R2(2)*R3(1)) 

     if (abs(cell_volume2)>0.001d0) cell_volume2= 2d0*3.1415926535d0/cell_volume2

     if (cell_volume2<0) then
        R3=-R3
        Umatrix(3, :)= -Umatrix(3, :)
     endif

     if (abs(abs(cell_volume2)-abs(cell_volume))> 0.001d0.and.cpuid==0) then
        write(stdout, '(a)')' '
        write(stdout, '(2a)')' Warnning: The Umatrix is wrongly set, the new cell', &
           'volume should be the same as the old ones. '
        write(stdout, '(a,2f10.4)')' cell_volume vs cell_volume-new', cell_volume, cell_volume2
        write(stdout, '(a)')" However, don't worry, WannierTools will help you to find a suitable rotation matrix."
        write(stdout, '(a)')" I am looking for new unit cell atuomatically: "
     endif     !> print out the new basis
     if (abs(abs(cell_volume2)-abs(cell_volume))> 0.001d0) then
        call FindTheThirdLatticeVector()
        R1= Umatrix(1, 1)*Rua+ Umatrix(1, 2)*Rub+ Umatrix(1, 3)*Ruc
        R2= Umatrix(2, 1)*Rua+ Umatrix(2, 2)*Rub+ Umatrix(2, 3)*Ruc
        R3= Umatrix(3, 1)*Rua+ Umatrix(3, 2)*Rub+ Umatrix(3, 3)*Ruc
     endif     !> print out the new basis

     if (cpuid.eq.0) then
        write(stdout, *)" "
        write(stdout, *)"The rotated new unit cell : "
        write(stdout, '(3f12.6)') R1
        write(stdout, '(3f12.6)') R2
        write(stdout, '(3f12.6)') R3
        write(stdout, *)" "
     endif

     if (cpuid.eq.0) then
        write(stdout, *)"Fractional coordinates of atoms in units of new lattice vectors : "
        do ia=1, Num_atoms
           call rotate_newlattice(Atom_position_direct(:, ia), Rt)
           if(cpuid==0)write(stdout, '(a4,3f12.6)')atom_name(ia), Rt
        enddo
        write(stdout, *)" "
     endif

     !> print out the new basis
     outfileindex= outfileindex+ 1
     if (cpuid.eq.0) then
        open(outfileindex, file="POSCAR-rotated")
        write(outfileindex, '(a)')"Rotated POSCAR by SURFACE card in wt.in by WannierTools"
        write(outfileindex, '(a)')"1.0"
        write(outfileindex, '(3f12.6)') R1
        write(outfileindex, '(3f12.6)') R2
        write(outfileindex, '(3f12.6)') R3
        write(outfileindex, '(1000i5)') Num_atoms
        write(outfileindex, '(a)')"Direct"
        do ia=1, Num_atoms
           call rotate_newlattice(Atom_position_direct(:, ia), Rt)
           if(cpuid==0)write(outfileindex, '(3f12.6, a9)')Rt, trim(adjustl(atom_name(ia)))  
        enddo
        close(outfileindex)
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
     it=0
     do i=1, nk3lines
        read(1001, *, err=201) k3line_name(i), k3line_start(:, i), &
                      char_temp, k3line_end(:, i)
        it= it+ 1
        if(cpuid==0)write(stdout, '(a5, 3f9.4, 2x, a5, 3f9.4)')&
          k3line_name(i), k3line_start(:, i), &
          char_temp, k3line_end(:, i)

     enddo
     201 continue
     if (it< nk3lines.and.cpuid==0) then
        write(stdout, *)' '
        write(stdout, *)' >>>> Error happens in KPATH_BULK card'
        write(stdout, *)' Error: the number of kpath lines should consistent'
        write(stdout, *)' Nk3lines ', nk3lines, ' the k lines you given ', it
        write(stdout, *)' Please set nk3lines to be  ', it
        stop
     endif

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
     if (.not.lfound .and. BulkBand_calc) then
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
     it = 0
     do i=1, nk2lines
        read(1001, *, err=202) k2line_name(i), kp(:, i), &
                      char_temp, ke(:, i)
        it= it+ 1
        if (cpuid==0) write(stdout, '(a6, 2f9.5, 4x, a6, 2f9.5)')&
                      k2line_name(i), kp(:, i), &
                      char_temp, ke(:, i)
     enddo
     202 continue

     if (it< nk2lines.and.cpuid==0) then
        write(stdout, *)' '
        write(stdout, *)' >>>> Error happens in KPATH_SLAB card'
        write(stdout, *)' Error: the number of kpath lines should consistent'
        write(stdout, *)' Nk2lines ', nk2lines, ' the k lines you given ', it
        write(stdout, *)' Please set nk2lines to be  ', it
        stop
     endif


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
     if (.not.lfound .and.(SlabBand_calc .or. SlabSS_calc)) then
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
     it = 0
     read(1001, *, err=203)K2D_start
     it= it+ 1
     read(1001, *, err=203)K2D_vec1
     it= it+ 1
     read(1001, *, err=203)K2D_vec2
     it= it+ 1
     203 continue
     if (it< 3.and.cpuid==0) then
        write(stdout, *)' '
        write(stdout, *)' >>>> Error happens in KPLANE_SLAB card'
        write(stdout, *)' Error: There are three lines in this card to specify the start point'
        write(stdout, *)' of vectors, and two lines to assign two vectors.'
        write(stdout, *)" If you don't know the meaning of this card, please delete this card"
        stop
     endif


     106 continue

     if (cpuid==0) write(stdout, *)'>> Kpoints plane for 2D system--> arcs  '
     if (cpuid==0) write(stdout, '((a, 2f8.4))')'K2D_start:', K2D_start
     if (cpuid==0) write(stdout, '((a, 2f8.4))')'The first vector: ', K2D_vec1
     if (cpuid==0) write(stdout, '((a, 2f8.4))')'The second vector: ', K2D_vec2
     if (.not.lfound .and.(SlabArc_calc  .or. SlabSpintexture_calc)) then
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
     it= 0
     read(1001, *, err=206)K3D_start
     it= it+ 1
     read(1001, *, err=206)K3D_vec1
     it= it+ 1
     read(1001, *, err=206)K3D_vec2
     it= it+ 1
     206 continue
     if (it< 3.and.cpuid==0) then
        write(stdout, *)' '
        write(stdout, *)' >>>> Error happens in KPLANE_BULK card'
        write(stdout, *)' Error: There are three lines in this card to specify the start point'
        write(stdout, *)' of vectors, and two lines to assign two vectors.'
        write(stdout, *)" If you don't know the meaning of this card, please delete this card"
        stop
     endif


     107 continue

     if (cpuid==0) write(stdout, *)'>> Kpoints plane for 3D system--> gapshape  '
     if (cpuid==0) write(stdout, '((a, 3f8.4))')'k3D_start : ', K3D_start
     if (cpuid==0) write(stdout, '((a, 3f8.4))')'The 1st vector: ', K3D_vec1
     if (cpuid==0) write(stdout, '((a, 3f8.4))')'The 2nd vector: ', K3D_vec2
     if (.not.lfound .and.(BulkGap_plane_calc  .or. wanniercenter_calc)) then
        stop 'ERROR: please set KPLANE_bulk for gap or WCC calculations'
     endif



     !> read kcube_bulk information
     !> default value for KCUBE_BULK
     K3D_start_cube= (/ 0.0,  0.0,   0.0/)
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
     it= 0
     read(1001, *, err=204)K3D_start_cube
     it= it+ 1
     read(1001, *, err=204)K3D_vec1_cube
     it= it+ 1
     read(1001, *, err=204)K3D_vec2_cube
     it= it+ 1
     read(1001, *, err=204)K3D_vec3_cube
     it= it+ 1
     204 continue
     if (it< 3.and.cpuid==0) then
        write(stdout, *)' '
        write(stdout, *)' >>>> Error happens in KCUBE_BULK card'
        write(stdout, *)' Error: There are four lines in this card to specify the start point'
        write(stdout, *)' of vectors, and three lines to assign two vectors.'
        write(stdout, *)" If you don't know the meaning of this card, please delete this card"
        stop
     endif


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
     if (.not.lfound .and.(BulkGap_cube_calc)) then
        stop 'ERROR: please set KCUBE_BULK for gap3D calculations'
     endif

     !> set default parameters for Berry phase calculation
     NK_Berry= 2
     allocate(k3points_Berry(3, NK_Berry))
     DirectOrCart_Berry='Direct'
     k3points_Berry= 0d0
     rewind(1001)
     lfound = .false.
     do while (.true.)
        read(1001, *, end= 113)inline
        if (trim(adjustl(inline))=='KPATH_BERRY') then
           lfound= .true.
           if (cpuid==0) write(stdout, *)' '
           if (cpuid==0) write(stdout, *)'We found KPATH_BERRY card'
           exit
        endif
     enddo

     read(1001, *, end=208, err=208)NK_Berry
     if (cpuid==0) write(stdout, '(a, i10)')'NK_Berry', NK_Berry
     read(1001, *, end=208, err=208)inline   ! The unit of lattice vector
     DirectOrCart_Berry= trim(adjustl(inline))

     deallocate(k3points_Berry)
     allocate(k3points_Berry(3, NK_Berry))
     k3points_Berry= 0d0

     it= 0
     if (index(DirectOrCart_Berry, "D")>0)then
        do ik=1, NK_Berry
           read(1001, *, end=208, err=208)k3points_Berry(:, ik)   ! The unit of lattice vector
           it = it+ 1
        enddo
     else
        do ik=1, NK_Berry
           read(1001, *, end=208, err=208)k    ! The unit of lattice vector
           call cart_direct_rec(k, k3points_Berry(:, ik))
           it = it+ 1
        enddo
     endif  ! Direct or Cart coordinates
     208 continue

     if (it< NK_Berry.and. cpuid==0) then
        write(stdout, *)"ERROR: something wrong in the KPATH_BERRY card"
        write(stdout, *)"No. of kpoints for berry is not consistent with No. of lines"
        write(stdout, *)"I found ", it, " lines"
        write(stdout, *)"while you set NK_Berry to be ", NK_Berry
        stop
     endif
      
     113 continue

     if (cpuid==0) write(stdout, *)' '
     if (.not.lfound.and.cpuid==0.and.BerryPhase_calc)then
        write(stdout, *)'Error : you have to set KPATH_BERRY card with a list of k points'
        stop
     endif

     !< end of Berry phase setting


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

     it= 0
     read(1001, *, end=205, err=205)iband_mass
     it= it+ 1
     read(1001, *, end=205, err=205)dk_mass
     it= it+ 1
     read(1001, *, end=205, err=205)k_mass
     it= it+ 1
     205 continue
     if (it< 3.and.cpuid==0) then
        write(stdout, *)' '
        write(stdout, *)' >>>> Error happens in EFFECTIVE_MASS card'
        write(stdout, *)' Error: There are three lines in this card to specify the iband_mass'
        write(stdout, *)" , dk_mass and k_mass, like this: "
        write(stdout, *)"EFFECTIVE_MASS"
        write(stdout, *)" 6       ! the 6'th band"
        write(stdout, *)" 0.01    ! in unit of 1/Bohr"
        write(stdout, *)" 0 0 0   ! k point"
        stop
     endif


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

 

     !>> setting up a series of k points in 3D BZ
     rewind(1001)
     lfound = .false.
     do while (.true.)
        read(1001, *, end= 311)inline
        if (trim(adjustl(inline))=='KPOINTS_3D') then
           lfound= .true.
           if (cpuid==0) write(stdout, *)' '
           if (cpuid==0) write(stdout, *)'We found KPOINTS_3D card'
           exit
        endif
     enddo

     read(1001, *, end=319, err=319, iostat=stat)Nk3_point_mode   ! The unit of lattice vector
     allocate(k3points_pointmode_cart(3, Nk3_point_mode))
     allocate(k3points_pointmode_direct(3, Nk3_point_mode))
     k3points_pointmode_cart= 0d0
     k3points_pointmode_direct= 0d0

     read(1001, *, end=319, err=319, iostat=stat)inline   ! The unit of lattice vector
     if (index(trim(adjustl(inline)), "D")>0)then
        do ik= 1, Nk3_point_mode
           read(1001, *, end=319, err=319, iostat=stat)k3points_pointmode_direct(:, ik)
           call direct_cart_rec(k3points_pointmode_direct(:, ik), &
              k3points_pointmode_cart(:, ik))
        enddo
     else
        do ik= 1, Nk3_point_mode
           read(1001, *, end=319, err=319, iostat=stat)k3points_pointmode_cart(:, ik)
           call cart_direct_rec(k3points_pointmode_cart(:, ik), &
              k3points_pointmode_direct(:, ik))
        enddo
     endif
     
     319 continue
     if (stat/=0 .and. cpuid==0) then
        write(stdout, '(8f10.5)') "ERROR: there are something wrong in KPOINTS_3D card"
        write(stdout, '(8f10.5)') "It should be like this:"
        write(stdout, '(8f10.5)') "The number of lines below 'Direct' "
        write(stdout, '(8f10.5)') "should be the same as the number of k points"
        write(stdout, '(8f10.5)') "KPOINTS_3D"
        write(stdout, '(8f10.5)') "4 ! number of k points"
        write(stdout, '(8f10.5)') "Direct"
        write(stdout, '(8f10.5)') "0.0  0.0  0.0"
        write(stdout, '(8f10.5)') "0.5  0.0  0.0"
        write(stdout, '(8f10.5)') "0.0  0.5  0.0"
        write(stdout, '(8f10.5)') "0.0  0.0  0.5"
     endif

     !> print out the single kpoint positions
     if (cpuid==0) then
        write(stdout, '(a)')" "
        write(stdout, '(a)')"KPOINTS_3D positions"
        do ik=1, Nk3_point_mode
           write(stdout, '(8a10)')"kx", 'ky', 'kz', 'k1', 'k2', 'k3'
           write(stdout, '(8f10.5)')k3points_pointmode_cart(:, ik), k3points_pointmode_direct(:, ik)
        enddo
        write(stdout, '(a)')" "
     endif

     321 continue
     if (cpuid==0) write(stdout, *)' '
     if (.not. lfound) then
        Nk3_point_mode = 1
        allocate(k3points_pointmode_cart(3, Nk3_point_mode))
        allocate(k3points_pointmode_direct(3, Nk3_point_mode))
        k3points_pointmode_cart= 0d0
        k3points_pointmode_direct= 0d0
     endif
     if (.not.lfound.and.cpuid==0)write(stdout, *)'>> We use the default values for k3points_pointmode_direct=[0,0,0]' 
     if (.not.lfound.and.cpuid==0)write(stdout, *)'>> and Nk3_point_mode = 1'





     !>> setting up a single k points in 3D BZ
     Single_KPOINT_3D_CART= [0.d0, 0d0, 0d0]
     Single_KPOINT_3D_DIRECT= [0.d0, 0d0, 0d0]
     rewind(1001)
     lfound = .false.
     do while (.true.)
        read(1001, *, end= 311)inline
        if (trim(adjustl(inline))=='SINGLEKPOINT_3D') then
           lfound= .true.
           if (cpuid==0) write(stdout, *)' '
           if (cpuid==0) write(stdout, *)'We found SINGLEKPOINT_3D card'
           exit
        endif
     enddo

     read(1001, *, end=309, err=309, iostat=stat)inline   ! The unit of lattice vector
     DirectOrCart_SINGLE= trim(adjustl(inline))
     if (index(DirectOrCart_SINGLE, "D")>0)then
        read(1001, *, end=309, err=309, iostat=stat)Single_KPOINT_3D_DIRECT(1:3)
        call direct_cart_rec(Single_KPOINT_3D_DIRECT, Single_KPOINT_3D_CART)
     else
        read(1001, *, end=309, err=309, iostat=stat)Single_KPOINT_3D_CART(1:3)
        call cart_direct_rec(Single_KPOINT_3D_CART, Single_KPOINT_3D_DIRECT)
     endif
     
     309 continue
     if (stat/=0 .and. cpuid==0) then
        write(stdout, '(8f10.5)') "ERROR: there are something wrong in SINGLEKPOINT_3D card"
        write(stdout, '(8f10.5)') "It should be like this:"
        write(stdout, '(8f10.5)') "SINGLEKPOINT_3D"
        write(stdout, '(8f10.5)') "Direct"
        write(stdout, '(8f10.5)') "0.0  0.0  0.0"
     endif

     !> print out the single kpoint positions
     if (cpuid==0) then
        write(stdout, '(a)')" "
        write(stdout, '(a)')"Single_KPOINT_3D positions"
        write(stdout, '(8a10)')"kx", 'ky', 'kz', 'k1', 'k2', 'k3'
        write(stdout, '(8f10.5)') Single_KPOINT_3D_CART, Single_KPOINT_3D_DIRECT
        write(stdout, '(a)')" "
     endif

     311 continue
     if (cpuid==0) write(stdout, *)' '
     if (.not.lfound.and.cpuid==0)write(stdout, *)'>> We use the default values for Single_KPOINT_3D_DIRECT=[0,0,0]' 


   
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


     !> setup for Weyl points chirality calculation
     !> default
     Num_Weyls= 0 ! in unit of 1/Ang
     kr0=0
     rewind(1001)
     lfound = .false.
     do while (.true.)
        read(1001, *, end= 111)inline
        if (trim(adjustl(inline))=='WEYL_CHIRALITY') then
           lfound= .true.
           if (cpuid==0) write(stdout, *)' '
           if (cpuid==0) write(stdout, *)'We found WEYL_CHIRALITY card'
           exit
        endif
     enddo

     read(1001, *, end=209, err=209)Num_Weyls
     allocate(weyl_position_direct(3, Num_Weyls))
     allocate(weyl_position_cart(3, Num_Weyls))
     weyl_position_direct= 0d0
     weyl_position_cart= 0d0
     read(1001, *, end=209, err=209)inline   ! The unit of lattice vector
     DirectOrCart_Weyl= trim(adjustl(inline))
     read(1001, *, end=209, err=209)kr0
     it= 0
     do i=1, Num_Weyls
        if (index(DirectOrCart_Weyl, "D")>0)then
           read(1001, *, end=209, err=209)weyl_position_direct(:, i)
           call direct_cart_rec(weyl_position_direct(:, i), weyl_position_cart(:, i))
           it = it+ 1
        else
           read(1001, *, end=209, err=209)weyl_position_cart(:, i)
           call cart_direct_rec(weyl_position_cart(:, i), weyl_position_direct(:, i))
           it = it+ 1
        endif
     enddo
     
     209 continue

     !> print out the Weyl positions
     if (cpuid==0) then
        write(stdout, '(a)')" "
        write(stdout, '(a)')"Weyl point positions"
        write(stdout, '(8a10)')"kx", 'ky', 'kz', 'k1', 'k2', 'k3'
        do i=1, Num_Weyls
           write(stdout, '(8f10.5)')weyl_position_cart(:, i), weyl_position_direct(:, i)
        enddo

        write(stdout, '(a)')" "
     endif

 
     if (it< Num_Weyls.and.cpuid==0) then
        write(stdout, *)' Error: Num_Weyls should the same as ', &
           ' the number of weyl position lines'
        write(stdout, *)' Num_Weyls = ', Num_Weyls
        write(stdout, *)' Num of pos lines = ', it
        stop
     endif

     111 continue
     if (cpuid==0) write(stdout, *)' '
     if (.not.lfound.and.cpuid==0)write(stdout, *)'>> We do not calculate chirality for weyl points'
     if (.not.lfound.and.WeylChirality_calc.and.cpuid==0) then
        write(stdout, *) 'ERROR: you should specify the WEYL_CHIRALITY card, see documentation'
     endif

     !> parameters for selectedorbitals
     rewind(1001)
     lfound = .false.
     stat=0
     do while (.true.)
        read(1001, *, end= 331)inline
        if (trim(adjustl(inline))=='SELECTEDORBITALS') then
           lfound= .true.
           if (cpuid==0) write(stdout, *)' '
           if (cpuid==0) write(stdout, *)'We found SELECTEDORBITALS card'
           exit
        endif
     enddo

     it= 0
     stat= 0
     NumberofSelectedOrbitals= 0
     read(1001, *, err=331, iostat=stat)NumberofSelectedOrbitals

     if (NumberofSelectedOrbitals>0) then
        allocate(Selected_Orbitals(NumberofSelectedOrbitals))
        Selected_Orbitals= 0
        read(1001, *, err=331, iostat=stat) (Selected_Orbitals(i), i=1, NumberofSelectedOrbitals)
     else
        stop 'NumberofSelectedOrbitals should be an integer and larger than zero'
     endif
   
     331 continue
     if (stat/=0) then
        if (cpuid==0) write(stdout, *)' '
        if (cpuid==0) write(stdout, *)' '
        if (cpuid==0) write(stdout, '(a)')'Error: Please set the right number of orbitals and orbitals in wt.in like this:'
        if (cpuid==0) write(stdout, *)' '
        if (cpuid==0) write(stdout, '(a)')'SelectedOrbitals'
        if (cpuid==0) write(stdout, '(a)')'4  ! number of selected orbitals'
        if (cpuid==0) write(stdout, '(a)')'4 5 6 7 ! orbitals indices '
        stop 'Errors happen in the WT.in, please check informations in the WT.out'
     endif


     !> parameters for selectedbands
     rewind(1001)
     lfound = .false.
     stat=0
     do while (.true.)
        read(1001, *, end= 231)inline
        if (trim(adjustl(inline))=='SELECTEDBANDS') then
           lfound= .true.
           if (cpuid==0) write(stdout, *)' '
           if (cpuid==0) write(stdout, *)'We found SELECTEDBANDS card'
           exit
        endif
     enddo

     it= 0
     stat= 0
     NumberofSelectedBands= 0
     read(1001, *, err=231, iostat=stat)NumberofSelectedBands

     if (NumberofSelectedBands>0) then
        allocate(Selected_band_index(NumberofSelectedBands))
        Selected_band_index= 0
        read(1001, *, err=231, iostat=stat) (Selected_band_index(i), i=1, NumberofSelectedBands)
     else
        stop 'NumberofSelectedBands should be an integer and larger than zero'
     endif

     231 continue
     if (stat/=0) then
        if (cpuid==0) write(stdout, *)' '
        if (cpuid==0) write(stdout, *)' '
        if (cpuid==0) write(stdout, '(a)')'Error: Please set the right number of bands and band indices in wt.in like this:'
        if (cpuid==0) write(stdout, *)' '
        if (cpuid==0) write(stdout, '(a)')'SELECTEDBANDS'
        if (cpuid==0) write(stdout, '(a)')'4  ! number of selected bands'
        if (cpuid==0) write(stdout, '(a)')'4 5 6 7 ! band indices '
        stop 'Errors happen in the WT.in, please check informations in the WT.out'
     endif


     !> for phonon system,  LO-TO correction, by T.T Zhang
     !> Atomic MASS in unit of g/mol
     rewind(1001)
     lfound = .false.
     do while (.true.)
        read(1001, *, end= 221)inline
        if (trim(adjustl(inline))=='ATOM_MASS') then
           lfound= .true.
           if (cpuid==0) write(stdout, *)' '
           if (cpuid==0) write(stdout, *)'We found ATOM_MASS card'
           exit
        endif
     enddo
     221 continue

     if (lfound) then
        read(1001,*)Num_atom_type
        if(cpuid==0)write(stdout,'(a,i10,a)')'There are', Num_atom_type, 'kind of atoms'
        if (.not. allocated(Num_atoms_eachtype)) allocate(Num_atoms_eachtype(Num_atom_type))
        if (.not. allocated(mass_temp)) allocate(mass_temp(Num_atom_type))
        if (.not. allocated(ATOM_MASS)) allocate(ATOM_MASS(Num_atoms))
        read(1001,*)Num_atoms_eachtype(1:Num_atom_type)
        read(1001,*)mass_temp(1:Num_atom_type)
        do i= 1, Num_atom_type
           write(stdout,'(a,i10,a)')'Each type have', Num_atoms_eachtype(i), ' atoms'
           write(stdout,'(a,f12.6)')'And their mass is', mass_temp(i)
        enddo
        it=0
        do i=1, Num_atom_type   
           do j=1, Num_atoms_eachtype(i)
              it=it+1
              ATOM_MASS(it)=mass_temp(i)
           enddo
        enddo
     else 
        if (LOTO_correction)stop "ERROR: please set ATOM_MASS card for LOTO correction of phonon spectrum"
     endif

     !>> setup Dielectric tensor for a given material
     rewind(1001)
     lfound = .false.
     do while (.true.)
        read(1001, *, end= 222)inline
        if (trim(adjustl(inline))=='LOTO_DT') then
           lfound= .true.
           if (cpuid==0) write(stdout, *)' '
           if (cpuid==0) write(stdout, *)'We found LOTO_DT card for LOTO correction'
           exit
        endif
     enddo
     222 continue

     if (lfound) then
           read(1001, *)Diele_Tensor(1,:)   ! Diele_Tensor is a 3*3 tensor for a material
           write(stdout,'(a,3f12.5)')'Diele_tensor(1,:)',Diele_Tensor(1,:)
           read(1001, *)Diele_Tensor(2,:)
           write(stdout,'(a,3f12.5)')'Diele_tensor(2,:)',Diele_Tensor(2,:)
           read(1001, *)Diele_Tensor(3,:)
           write(stdout,'(a,3f12.5)')'Diele_tensor(3,:)',Diele_Tensor(3,:)
     else 
        if (LOTO_correction) then
           if (cpuid==0) then
              write(stdout, *)"ERROR: please set LOTO_DT card for LOTO correction of phonon spectrum"
              write(stdout, *)"ERROR: please set Dielectronic Tensor information"
              stop "ERROR: Check error messages in WT.out"
           endif
        endif
     endif

     !>> setup Born charge for a given material
     rewind(1001)
     lfound = .false.
     do while (.true.)
        read(1001, *, end= 223)inline
        if (trim(adjustl(inline))=='LOTO_BC') then
           lfound= .true.
           if (cpuid==0) write(stdout, *)' '
           if (cpuid==0) write(stdout, *)'We found LOTO_BC card for LOTO correction'
           exit
        endif
     enddo
     223 continue

     if (lfound) then
        it=0
        allocate(Born_Charge(Num_atoms,3,3))
        allocate(Born_Charge_temp(Num_atom_type,3,3))
        do i=1,Num_atom_type  
           read(1001, *)Born_Charge_temp(i,1,:)
           read(1001, *)Born_Charge_temp(i,2,:)
           read(1001, *)Born_Charge_temp(i,3,:)
           do j=1,Num_atoms_eachtype(i)
              it=it+1
              Born_Charge(it,:,:)=Born_Charge_temp(i,:,:)
              write(stdout,'(a,i3,2X,a6)')'Born_Charge for atom ', it, atom_name(it)
              write(stdout,'(3f12.5)')Born_Charge(it,1,:)
              write(stdout,'(3f12.5)')Born_Charge(it,2,:)
              write(stdout,'(3f12.5)')Born_Charge(it,3,:)
           enddo
        enddo
     else 
        if (LOTO_correction) then
           if (cpuid==0) then
              write(stdout, *)"ERROR: please set LOTO_BC card for LOTO correction of phonon spectrum"
              write(stdout, *)"ERROR: please set Born charge information"
              stop "ERROR: Check error messages in WT.out"
           endif
        endif
     endif




     !> close wt.in 
     close(1001)

     eta=(omegamax- omegamin)/omeganum*2d0

     if(cpuid==0)write(stdout,*)'<<<Read wt.in file successfully'


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

  !> rotate a vector in unit of  the original lattice vector into the new lattice
  !> vector defined by Umatrix
  subroutine rotate_newlattice(R1, R2)
     use para, only : dp, Umatrix
     implicit none
     real(dp), intent(in) :: R1(3)
     real(dp), intent(inout) :: R2(3)
     real(dp) :: Umatrix_inv(3, 3)
     Umatrix_inv= Umatrix

     call inv_r(3, Umatrix_inv)

     R2(1)= Umatrix_inv(1, 1)*R1(1)+ Umatrix_inv(2, 1)*R1(2)+ Umatrix_inv(3, 1)*R1(3)
     R2(2)= Umatrix_inv(1, 2)*R1(1)+ Umatrix_inv(2, 2)*R1(2)+ Umatrix_inv(3, 2)*R1(3)
     R2(3)= Umatrix_inv(1, 3)*R1(1)+ Umatrix_inv(2, 3)*R1(2)+ Umatrix_inv(3, 3)*R1(3)

     return
  end subroutine rotate_newlattice




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



   !> define a new unit cell with the given MillerIndices [hkl]
   subroutine MillerIndicestoumatrix()
      use para
      implicit none
      integer :: i1, i2, i3, j1, j2, j3, h, k, l, it
      real(dp) :: R(3), R1(3), R2(3), R3(3), Rhkl(3), dot

      integer, allocatable :: vector_on_hkl_surface(:, :)
      integer :: Nvectors_on_hkl_surface
      real(dp) :: smallest_area, area
      real(dp) :: largestangle, angle
      real(dp) :: smallest_volume, cell_volume
      real(dp) :: smallest_length
      real(dp) :: norm_1, norm_2, norm_3

      integer :: iRmax
      iRmax= 10

      allocate(vector_on_hkl_surface(3, (2*iRmax+1)**3))
      vector_on_hkl_surface= 0

      h= MillerIndices(1)
      k= MillerIndices(2)
      l= MillerIndices(3)
      Rhkl= h*Rua+ k*Rub+ l*Ruc

      !> Firstly, find all vectors that are orthorgonal to hkl
      it= 0
      do i1=-iRmax, iRmax
         do i2=-iRmax, iRmax
            do i3=-iRmax, iRmax
               if (i1==0 .and. i2==0 .and. i3==0) cycle
              !R= i1*Rua+i2*Rub+i3*Ruc
              !dot=abs(R(1)*Rhkl(1)+ R(2)*Rhkl(2)+ R(3)*Rhkl(3))
               dot=abs(i1*h+i2*k+i3*l)
               if (dot<eps9) then
                  it= it+1
                  vector_on_hkl_surface(1, it)= i1
                  vector_on_hkl_surface(2, it)= i2
                  vector_on_hkl_surface(3, it)= i3
               endif
            enddo
         enddo
      enddo
      Nvectors_on_hkl_surface= it

      !> secondly, find the smallest area and the largest
      !> angle of two vectors
      smallest_area= 99999999d0
      largestangle= 0
      do i1=1, Nvectors_on_hkl_surface
         do i2=i1+1, Nvectors_on_hkl_surface
            R1= vector_on_hkl_surface(1, i1)*Rua+ &
                vector_on_hkl_surface(2, i1)*Rub+ &
                vector_on_hkl_surface(3, i1)*Ruc  
            R2= vector_on_hkl_surface(1, i2)*Rua+ &
                vector_on_hkl_surface(2, i2)*Rub+ &
                vector_on_hkl_surface(3, i2)*Ruc  
            dot= R1(1)*R2(1)+ R1(2)*R2(2)+ R1(3)*R2(3)
            norm_1= sqrt(R1(1)*R1(1)+ R1(2)*R1(2)+ R1(3)*R1(3))
            norm_2= sqrt(R2(1)*R2(1)+ R2(2)*R2(2)+ R2(3)*R2(3))
            if (norm_1*norm_2<eps9) stop 'norm of vector should larger than zero'
            R3(1)= R1(2)*R2(3)- R1(3)*R2(2)
            R3(2)= R1(3)*R2(1)- R1(1)*R2(3)
            R3(3)= R1(1)*R2(2)- R1(2)*R2(1)
            area= sqrt(R3(1)*R3(1)+ R3(2)*R3(2)+ R3(3)*R3(3))
            if (area<eps6) cycle
            if (area< smallest_area)smallest_area= area
         enddo
      enddo

      !> thirdly, find the largest
      !> angle in those two vectors which have smallest area 
      largestangle= 0
      do i1=1, Nvectors_on_hkl_surface
         do i2=i1+1, Nvectors_on_hkl_surface
            R1= vector_on_hkl_surface(1, i1)*Rua+ &
                vector_on_hkl_surface(2, i1)*Rub+ &
                vector_on_hkl_surface(3, i1)*Ruc  
            R2= vector_on_hkl_surface(1, i2)*Rua+ &
                vector_on_hkl_surface(2, i2)*Rub+ &
                vector_on_hkl_surface(3, i2)*Ruc  

            R3(1)= R1(2)*R2(3)- R1(3)*R2(2)
            R3(2)= R1(3)*R2(1)- R1(1)*R2(3)
            R3(3)= R1(1)*R2(2)- R1(2)*R2(1)
            area= dsqrt(R3(1)*R3(1)+ R3(2)*R3(2)+ R3(3)*R3(3))

            dot= R1(1)*R2(1)+ R1(2)*R2(2)+ R1(3)*R2(3)
            norm_1= dsqrt(R1(1)*R1(1)+ R1(2)*R1(2)+ R1(3)*R1(3))
            norm_2= dsqrt(R2(1)*R2(1)+ R2(2)*R2(2)+ R2(3)*R2(3))
            angle= dacos(dot/norm_1/norm_2)
            if (angle>pi/2d0) angle= abs(angle-pi)

            if (dabs(area- smallest_area)<eps9)then 
               if (angle> largestangle)largestangle= angle
            endif
         enddo
      enddo 

      !> thirdly, find the two vectors which have smallest area and largest
      !> angle
  l1: do i1=1, Nvectors_on_hkl_surface
         do i2=i1+1, Nvectors_on_hkl_surface
            R1= vector_on_hkl_surface(1, i1)*Rua+ &
                vector_on_hkl_surface(2, i1)*Rub+ &
                vector_on_hkl_surface(3, i1)*Ruc  
            R2= vector_on_hkl_surface(1, i2)*Rua+ &
                vector_on_hkl_surface(2, i2)*Rub+ &
                vector_on_hkl_surface(3, i2)*Ruc  
            dot= R1(1)*R2(1)+ R1(2)*R2(2)+ R1(3)*R2(3)
            norm_1= dsqrt(R1(1)*R1(1)+ R1(2)*R1(2)+ R1(3)*R1(3))
            norm_2= dsqrt(R2(1)*R2(1)+ R2(2)*R2(2)+ R2(3)*R2(3))
            angle= dacos(dot/norm_1/norm_2)
            R3(1)= R1(2)*R2(3)- R1(3)*R2(2)
            R3(2)= R1(3)*R2(1)- R1(1)*R2(3)
            R3(3)= R1(1)*R2(2)- R1(2)*R2(1)
            area= dsqrt(R3(1)*R3(1)+ R3(2)*R3(2)+ R3(3)*R3(3))
            if (angle>pi/2d0) angle= abs(angle-pi)
            if (dabs(area- smallest_area)<eps9 .and. dabs(angle-largestangle)<eps9)then 
               Umatrix(1, :)= vector_on_hkl_surface(:, i1)
               Umatrix(2, :)= vector_on_hkl_surface(:, i2)
               exit l1
            endif
         enddo 
      enddo l1
 
      !> The last step, find the third vector that makes the new unit cell has
      !> the same volume as the old unit cell
      smallest_volume= 9999999d0
      R1= Umatrix(1, 1)*Rua+  Umatrix(1, 2)*Rub+ Umatrix(1, 3)*Ruc  
      R2= Umatrix(2, 1)*Rua+  Umatrix(2, 2)*Rub+ Umatrix(2, 3)*Ruc  
      do i1=-iRmax, iRmax
         do i2=-iRmax, iRmax
            do i3=-iRmax, iRmax
               if (i1==0 .and. i2==0 .and. i3==0) cycle
               R3= i1*Rua+i2*Rub+i3*Ruc
               cell_volume= R1(1)*(R2(2)*R3(3)- R2(3)*R3(2)) &
                 +R1(2)*(R2(3)*R3(1)- R2(1)*R3(3)) &
                 +R1(3)*(R2(1)*R3(2)- R2(2)*R3(1)) 
               cell_volume= dabs(cell_volume)
               if (cell_volume< eps9) cycle
               if (cell_volume< smallest_volume) smallest_volume= cell_volume
            enddo
         enddo
      enddo

      smallest_length= 9999999d0
      do i1=-iRmax, iRmax
         do i2=-iRmax, iRmax
            do i3=-iRmax, iRmax
               if (i1==0 .and. i2==0 .and. i3==0) cycle
               R3= i1*Rua+i2*Rub+i3*Ruc
               cell_volume= R1(1)*(R2(2)*R3(3)- R2(3)*R3(2)) &
                 +R1(2)*(R2(3)*R3(1)- R2(1)*R3(3)) &
                 +R1(3)*(R2(1)*R3(2)- R2(2)*R3(1)) 
               cell_volume= dabs(cell_volume)
               if (dabs(cell_volume- smallest_volume)<eps9) then
                  norm_3= dsqrt(R3(1)*R3(1)+ R3(2)*R3(2)+ R3(3)*R3(3))
                  if (norm_3< smallest_length) smallest_length= norm_3
               endif
            enddo
         enddo
      enddo 

      l2: do i1=-iRmax, iRmax
         do i2=-iRmax, iRmax
            do i3=-iRmax, iRmax
               if (i1==0 .and. i2==0 .and. i3==0) cycle
               R3= i1*Rua+i2*Rub+i3*Ruc
               cell_volume= R1(1)*(R2(2)*R3(3)- R2(3)*R3(2)) &
                 +R1(2)*(R2(3)*R3(1)- R2(1)*R3(3)) &
                 +R1(3)*(R2(1)*R3(2)- R2(2)*R3(1)) 
               cell_volume= dabs(cell_volume)
               norm_3= dsqrt(R3(1)*R3(1)+ R3(2)*R3(2)+ R3(3)*R3(3))
               if (dabs(cell_volume- smallest_volume)<eps9 .and. &
                   dabs(norm_3- smallest_length)<eps9) then
                  Umatrix(3, 1)= i1
                  Umatrix(3, 2)= i2
                  Umatrix(3, 3)= i3
                  exit l2
               endif
            enddo
         enddo
      enddo l2


      R1= Umatrix(1, 1)*Rua+  Umatrix(1, 2)*Rub+ Umatrix(1, 3)*Ruc  
      R2= Umatrix(2, 1)*Rua+  Umatrix(2, 2)*Rub+ Umatrix(2, 3)*Ruc  
      R3= Umatrix(3, 1)*Rua+  Umatrix(3, 2)*Rub+ Umatrix(3, 3)*Ruc  
      cell_volume= R1(1)*(R2(2)*R3(3)- R2(3)*R3(2)) &
                 + R1(2)*(R2(3)*R3(1)- R2(1)*R3(3)) &
                 + R1(3)*(R2(1)*R3(2)- R2(2)*R3(1)) 
      if (cell_volume<0) Umatrix(3, :)= -Umatrix(3, :)

      R1= Umatrix(1, 1)*Rua+  Umatrix(1, 2)*Rub+ Umatrix(1, 3)*Ruc  
      R2= Umatrix(2, 1)*Rua+  Umatrix(2, 2)*Rub+ Umatrix(2, 3)*Ruc  
      R3= Umatrix(3, 1)*Rua+  Umatrix(3, 2)*Rub+ Umatrix(3, 3)*Ruc  
      cell_volume= R1(1)*(R2(2)*R3(3)- R2(3)*R3(2)) &
                 + R1(2)*(R2(3)*R3(1)- R2(1)*R3(3)) &
                 + R1(3)*(R2(1)*R3(2)- R2(2)*R3(1)) 

      if (abs(cell_volume- CellVolume)< eps9 .and. cpuid==0) then
         write(stdout, *)'  Congratulations, you got a unit cell that has ', &
         ' the same volume as the original unit cell '
         write(stdout, *)' The unitary rotation matrix is : '
         write(stdout, '(3f10.3)')Umatrix(1,:)
         write(stdout, '(3f10.3)')Umatrix(2,:)
         write(stdout, '(3f10.3)')Umatrix(3,:)

         write(stdout, *)' '
         write(stdout, *)'The lattice vectors for new cell are : '
         write(stdout, '(a,3f10.3)')' R1=', R1
         write(stdout, '(a,3f10.3)')' R2=', R2
         write(stdout, '(a,3f10.3)')' R3=', R3
         write(stdout, *)' Where R1 and R2 are in (hkl) plane'
      endif

      return
   end subroutine MillerIndicestoumatrix


   !> define a new unit cell with the given MillerIndices [hkl]
   subroutine FindTheThirdLatticeVector()
      use para
      implicit none
      integer :: i1, i2, i3, j1, j2, j3, h, k, l, it
      real(dp) :: R(3), R1(3), R2(3), R3(3), cross(3), dot

      integer, allocatable :: vectors_parallel_umatrix1(:, :)
      integer, allocatable :: vectors_parallel_umatrix2(:, :)
      integer :: Nvectors_parallel_umatrix1, Nvectors_parallel_umatrix2
      real(dp) :: smallest_volume, cell_volume
      real(dp) :: smallest_length
      real(dp) :: norm_1, norm_2, norm_3

      integer :: iRmax
      iRmax= 6 

      allocate(vectors_parallel_umatrix1(3, (2*iRmax+1)**3))
      allocate(vectors_parallel_umatrix2(3, (2*iRmax+1)**3))
      vectors_parallel_umatrix1= 0
      vectors_parallel_umatrix2= 0

      !> Firstly, find all vectors that are parallel to Umatrix(:,1)
      it= 0
      do i1=-iRmax, iRmax
         do i2=-iRmax, iRmax
            do i3=-iRmax, iRmax
               if (i1==0 .and. i2==0 .and. i3==0) cycle
              !R= i1*Rua+i2*Rub+i3*Ruc
              !dot=abs(R(1)*Rhkl(1)+ R(2)*Rhkl(2)+ R(3)*Rhkl(3))
               cross(1)= Umatrix(2, 1)*i3- i2*Umatrix(3, 1)
               cross(2)= Umatrix(3, 1)*i1- i3*Umatrix(1, 1)
               cross(3)= Umatrix(1, 1)*i2- i1*Umatrix(2, 1)
               dot = i1*Umatrix(1, 1)+ i2*Umatrix(2, 1)+ i3*Umatrix(3, 1)
               if ((abs(cross(1))+abs(cross(2))+abs(cross(3)))<eps9.and.dot>0) then
                  it= it+1
                  vectors_parallel_umatrix1(1, it)= i1
                  vectors_parallel_umatrix1(2, it)= i2
                  vectors_parallel_umatrix1(3, it)= i3
               endif
            enddo
         enddo
      enddo
      Nvectors_parallel_umatrix1= it

      !> and select the shortest vectors_parallel_umatrix1
      smallest_length= 9999999d0
      do it= 1, Nvectors_parallel_umatrix1
         R1= vectors_parallel_umatrix1(1, it)*Rua+ &
             vectors_parallel_umatrix1(2, it)*Rub+ &
             vectors_parallel_umatrix1(3, it)*Ruc  
         norm_3= dsqrt(R1(1)*R1(1)+ R1(2)*R1(2)+ R1(3)*R1(3))
         if (norm_3< smallest_length) then
            smallest_length= norm_3
            Umatrix(:, 1) = vectors_parallel_umatrix1(:, it)
         endif
      enddo 


      !> secondly, find all vectors that are parallel to Umatrix(:,2)
      it= 0
      do i1=-iRmax, iRmax
         do i2=-iRmax, iRmax
            do i3=-iRmax, iRmax
               if (i1==0 .and. i2==0 .and. i3==0) cycle
              !R= i1*Rua+i2*Rub+i3*Ruc
              !dot=abs(R(1)*Rhkl(1)+ R(2)*Rhkl(2)+ R(3)*Rhkl(3))
               cross(1)= Umatrix(2, 2)*i3- i2*Umatrix(3, 2)
               cross(2)= Umatrix(3, 2)*i1- i3*Umatrix(1, 2)
               cross(3)= Umatrix(1, 2)*i2- i1*Umatrix(2, 2)
               dot = i1*Umatrix(1, 2)+ i2*Umatrix(2, 2)+ i3*Umatrix(3, 2)
               if ((abs(cross(1))+abs(cross(2))+abs(cross(3)))<eps9.and.dot>0) then
                  it= it+1
                  vectors_parallel_umatrix2(1, it)= i1
                  vectors_parallel_umatrix2(2, it)= i2
                  vectors_parallel_umatrix2(3, it)= i3
               endif
            enddo
         enddo
      enddo
      Nvectors_parallel_umatrix2= it

      !> and select the shortest vectors_parallel_umatrix1
      smallest_length= 9999999d0
      do it= 1, Nvectors_parallel_umatrix2
         R2= vectors_parallel_umatrix2(1, it)*Rua+ &
             vectors_parallel_umatrix2(2, it)*Rub+ &
             vectors_parallel_umatrix2(3, it)*Ruc  
         norm_3= dsqrt(R2(1)*R2(1)+ R2(2)*R2(2)+ R2(3)*R2(3))
         if (norm_3< smallest_length) then
            smallest_length= norm_3
            Umatrix(:, 2) = vectors_parallel_umatrix2(:, it)
         endif
      enddo 

      !> The last step, find the third vector that makes the new unit cell has
      !> the same volume as the old unit cell
      smallest_volume= 9999999d0
      R1= Umatrix(1, 1)*Rua+  Umatrix(1, 2)*Rub+  Umatrix(1, 3)*Ruc
      R2= Umatrix(2, 1)*Rua+  Umatrix(2, 2)*Rub+  Umatrix(2, 3)*Ruc
      do i1=-iRmax, iRmax
         do i2=-iRmax, iRmax
            do i3=-iRmax, iRmax
               if (i1==0 .and. i2==0 .and. i3==0) cycle
               R3= i1*Rua+i2*Rub+i3*Ruc
               cell_volume= R1(1)*(R2(2)*R3(3)- R2(3)*R3(2)) &
                 +R1(2)*(R2(3)*R3(1)- R2(1)*R3(3)) &
                 +R1(3)*(R2(1)*R3(2)- R2(2)*R3(1)) 
               cell_volume= dabs(cell_volume)
               if (cell_volume< eps9) cycle
               if (cell_volume< smallest_volume) smallest_volume= cell_volume
            enddo
         enddo
      enddo

      !> find the third vector with the shortest length
      smallest_length= 9999999d0
      do i1=-iRmax, iRmax
         do i2=-iRmax, iRmax
            do i3=-iRmax, iRmax
               if (i1==0 .and. i2==0 .and. i3==0) cycle
               R3= i1*Rua+i2*Rub+i3*Ruc
               cell_volume= R1(1)*(R2(2)*R3(3)- R2(3)*R3(2)) &
                 +R1(2)*(R2(3)*R3(1)- R2(1)*R3(3)) &
                 +R1(3)*(R2(1)*R3(2)- R2(2)*R3(1)) 
               cell_volume= dabs(cell_volume)
               if (dabs(cell_volume- smallest_volume)<eps9) then
                  norm_3= dsqrt(R3(1)*R3(1)+ R3(2)*R3(2)+ R3(3)*R3(3))
                  if (norm_3< smallest_length) smallest_length= norm_3
               endif
            enddo
         enddo
      enddo 

      l2: do i1=-iRmax, iRmax
         do i2=-iRmax, iRmax
            do i3=-iRmax, iRmax
               if (i1==0 .and. i2==0 .and. i3==0) cycle
               R3= i1*Rua+i2*Rub+i3*Ruc
               cell_volume= R1(1)*(R2(2)*R3(3)- R2(3)*R3(2)) &
                 +R1(2)*(R2(3)*R3(1)- R2(1)*R3(3)) &
                 +R1(3)*(R2(1)*R3(2)- R2(2)*R3(1)) 
               cell_volume= dabs(cell_volume)
               norm_3= dsqrt(R3(1)*R3(1)+ R3(2)*R3(2)+ R3(3)*R3(3))
               if (dabs(cell_volume- smallest_volume)<eps9 .and. &
                   dabs(norm_3- smallest_length)<eps9) then
                  Umatrix(3, 1)= i1
                  Umatrix(3, 2)= i2
                  Umatrix(3, 3)= i3
                  exit l2
               endif
            enddo
         enddo
      enddo l2

      R1= Umatrix(1, 1)*Rua+  Umatrix(1, 2)*Rub+ Umatrix(1, 3)*Ruc  
      R2= Umatrix(2, 1)*Rua+  Umatrix(2, 2)*Rub+ Umatrix(2, 3)*Ruc  
      R3= Umatrix(3, 1)*Rua+  Umatrix(3, 2)*Rub+ Umatrix(3, 3)*Ruc  
      cell_volume= R1(1)*(R2(2)*R3(3)- R2(3)*R3(2)) &
                 + R1(2)*(R2(3)*R3(1)- R2(1)*R3(3)) &
                 + R1(3)*(R2(1)*R3(2)- R2(2)*R3(1)) 
      if (cell_volume<0) Umatrix(3, :)= -Umatrix(3, :)

      R1= Umatrix(1, 1)*Rua+  Umatrix(1, 2)*Rub+ Umatrix(1, 3)*Ruc  
      R2= Umatrix(2, 1)*Rua+  Umatrix(2, 2)*Rub+ Umatrix(2, 3)*Ruc  
      R3= Umatrix(3, 1)*Rua+  Umatrix(3, 2)*Rub+ Umatrix(3, 3)*Ruc  
      cell_volume= R1(1)*(R2(2)*R3(3)- R2(3)*R3(2)) &
                 + R1(2)*(R2(3)*R3(1)- R2(1)*R3(3)) &
                 + R1(3)*(R2(1)*R3(2)- R2(2)*R3(1)) 

      if (abs(cell_volume- CellVolume)< eps9 .and. cpuid==0) then
         write(stdout, *)'  Congratulations, you got a unit cell that has ', &
         ' the same volume as the original unit cell '
         write(stdout, *)' The unitary rotation matrix is : '
         write(stdout, '(3f10.3)')Umatrix(1,:)
         write(stdout, '(3f10.3)')Umatrix(2,:)
         write(stdout, '(3f10.3)')Umatrix(3,:)

         write(stdout, *)' '
         write(stdout, *)'The lattice vectors for new cell are : '
         write(stdout, '(a,3f10.3)')' R1=', R1
         write(stdout, '(a,3f10.3)')' R2=', R2
         write(stdout, '(a,3f10.3)')' R3=', R3
         write(stdout, *)' Where R1, R2, R3 are in cartesian coordinates'
      endif

      return
   end subroutine FindTheThirdLatticeVector

