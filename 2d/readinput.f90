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
     integer  :: j
     integer  :: NN
     integer :: nwann
     real(dp) :: t1, temp
     real(dp) :: pos(2)
     real(dp) :: k1(2), k2(2)
     real(dp) :: kstart(2), kend(2), kstart1, kend1
     real(dp) :: R1(2), R2(2) 
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

     !> inout file



    !read(1001,*)Hrfile
     read(1001, TB_FILE, iostat= stat)
     if (stat/=0) then
        Hrfile='wannier90_hr.dat'
        inquire(file='wannier90_hr.dat',exist=exists)
        if (.not.exists) stop "TB_FIlE namelist should be given or wannier90_hr.dat should exist"
     endif
     if(cpuid==0)write(stdout,'(1x, a, a25)')"Tight binding Hamiltonian file: ",Hrfile
    

     BulkBand_calc         = .FALSE.
     RibbonBand_calc         = .FALSE.
     SlabSS_calc           = .FALSE.
     Dos_calc           = .FALSE.
     wanniercenter_calc    = .FALSE.
     
     read(1001, CONTROL, iostat=stat)

     if (stat/=0) then
        write(*, *)"ERROR: namelist CONTROL should be set"
        write(*, *)"You should set one of these functions to be T"
        write(*, *)"BulkBand_calc,BulkFS_calc,BulkGap_cube_calc,BulkGap_plane_calc"
        write(*, *)"RibbonBand_calc,WireBand_calc,SlabSS_calc,SlabArc_calc "
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
        write(stdout, *) "RibbonBand_calc       : ",  RibbonBand_calc
        write(stdout, *) "SlabSS_calc         : ",  SlabSS_calc
        write(stdout, *) "wanniercenter_calc  : ", wanniercenter_calc
     endif

     !> set system parameters by default
     Nslab= 10
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
     else 
        stop 'ERROR: please set lattice information'
     endif

     if (index(AngOrBohr, 'Bohr')>0) then
        Rua= Rua*0.529177d0
        Rub= Rub*0.529177d0
     endif

     !> transform lattice from direct space to reciprocal space

     Kua= 0d0
     Kub= 0d0
     cell_volume=Rua(1)*Rub(2)- Rub(1)*Rua(2)
     cell_volume= abs(cell_volume)
     cell_volume= 2d0*3.1415926535d0/cell_volume

     Kua(1)= cell_volume*Rub(2)
     Kua(2)=-cell_volume*Rub(1)

     Kub(1)=-cell_volume*Rua(2)                
     Kub(2)= cell_volume*Rua(1)                

     if(cpuid==0)write(stdout, '(a)') '>> lattice information (Angstrom)'
     if(cpuid==0)write(stdout, '(2f12.6)')Rua
     if(cpuid==0)write(stdout, '(2f12.6)')Rub

     if(cpuid==0)write(stdout, '(a)') '>> Reciprocal lattice information (1/Angstrom)'
     if(cpuid==0)write(stdout, '(2f12.6)')Kua
     if(cpuid==0)write(stdout, '(2f12.6)')Kub

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
        allocate(Atom_position(2, Num_atoms))
        allocate(Atom_position_direct(2, Num_atoms))
        read(1001, *)inline   ! The unit of lattice vector
        DirectOrCart= trim(adjustl(inline))
        do i=1, Num_atoms
           read(1001, *) atom_name(i), Atom_position(:, i)
           if(cpuid==0)write(stdout, '(a4,2f12.6)')atom_name(i), Atom_position(:, i)
           if (index(DirectOrCart, "D")>0)then
              pos= Atom_position(:, i)
              Atom_position(:, i)= pos(1)*Rua+ pos(2)*Rub            
           endif
        enddo
        if(cpuid==0)write(stdout,'(a)')'Atom position in cartisen coordinate'
        do i=1, Num_atoms
           if(cpuid==0)write(stdout, '(a4,2f12.6)')atom_name(i), Atom_position(:, i)
        enddo
     
        if(cpuid==0)write(stdout,'(a)')'Atom position in direct coordinate'
        do ia=1, Num_atoms
           call cart_direct_real(Atom_position(:, ia), Atom_position_direct(:, ia))
           if(cpuid==0)write(stdout, '(a4,2f12.6)')atom_name(ia), Atom_position_direct(:, ia)
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


     !> read edge information
     rewind(1001)
     lfound = .false.
     do while (.true.)
        read(1001, *, end= 103)inline
        if (trim(adjustl(inline))=='EDGE') then
           lfound= .true.
           if (cpuid==0) write(stdout, *)' '
           if (cpuid==0) write(stdout, *)'We found EDGE card'
           exit
        endif
     enddo
     103 continue

     if (.not.lfound) then
        print *, inline
        stop 'ERROR: please set edge information'
     endif

     !> read information for new lattice 
     !> in order to get different edge state
     !> R1'=U11*R1+U12*R2
     !> R2'=U21*R1+U22*R2
     read(1001, *)Umatrix(1, :)
     read(1001, *)Umatrix(2, :)

     if (cpuid==0) then
        write(stdout, '(a)')'>> new vectors to define the edge (in unit of lattice vector)' 
        write(stdout, '(a, 2f12.6)')' The 1st vector on edge  :', Umatrix(1, :)
        write(stdout, '(a, 2f12.6)')' The 2nd vector on edge  :', Umatrix(2, :)
     endif

     !> check whether Umatrix is right
     !> the volume of the new cell should be the same as the old ones
     R1= Umatrix(1, 1)*Rua+ Umatrix(1, 2)*Rub
     R2= Umatrix(2, 1)*Rua+ Umatrix(2, 2)*Rub


     cell_volume2= R1(1)*R2(2)- R1(2)*R2(1)
     cell_volume2= 2d0*3.1415926535d0/cell_volume2

     if (cell_volume2<0) then
        R2=-R2
        Umatrix(2, :)= -Umatrix(2, :)
     endif

     if (abs(abs(cell_volume2)-abs(cell_volume))> 0.001d0.and.cpuid==0) then
        write(stdout, *)' ERROR The Umatrix is wrong, the new cell', &
           'volume should be the same as the old ones'
        write(stdout, '(a,2f10.4)')'cell_volume vs cell_volume-new', cell_volume, cell_volume2
        stop
     endif

     !> read kpath_slab information
     rewind(1001)
     lfound = .false.
     do while (.true.)
        read(1001, *, end= 105)inline
        if (trim(adjustl(inline))=='KPATH_BULK') then
           lfound= .true.
           if (cpuid==0) write(stdout, *)' '
           if (cpuid==0) write(stdout, *)'We found KPATH_BULK card'
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
           k1(1:2)= kstart(1)*Kua+ kstart(2)*Kub
           k2(1:2)= kend(1)*Kua+ kend(2)*Kub
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

     !> read kpath_ribbon information
     rewind(1001)
     lfound = .false.
     do while (.true.)
        read(1001, *, end= 116)inline
        if (trim(adjustl(inline))=='KPATH_RIBBON') then
           lfound= .true.
           if (cpuid==0) write(stdout, *)' '
           if (cpuid==0) write(stdout, *)'We found KPATH_RIBBON card'
           exit
        endif
     enddo

     !> read in k lines for 1D system
     k1line_name= ' '
     if (cpuid==0) write(stdout, *)'k lines for 1D system'
     read(1001, *)nk1lines
     if (cpuid==0) write(stdout, *)'No. of k lines for 1D system : ', nk1lines
     do i=1, nk1lines
        read(1001, *) k1line_name(i), kp1(i), &
                      char_temp, ke1(i)
        if(cpuid==0) write(stdout, '(a6, 1f9.5, 4x, a6, 1f9.5)')&
                      k1line_name(i), kp1(i), &
                      char_temp, ke1(i)
     enddo
     k1line_name(nk1lines+1) = char_temp
 
     NN= Nk
     knv1= NN*nk1lines
     allocate( k1_path(knv2))
     allocate( k1len (knv2))
     k1_path= 0d0
     k1len= 0d0

     t1=0d0
     k1len=0d0
     k1line_stop= 0d0
     do j=1, nk1lines 
        do i=1, NN
           kstart1= kp1(j)
           kend1  = ke1(j)
           k1_path(i+(j-1)*NN)= kstart1+ (kend1-kstart1)*(i-1)/dble(NN-1)
           
           temp= abs((kend1- kstart1))/dble(NN-1) 

           if (i.gt.1) then
              t1=t1+temp
           endif
           k1len(i+(j-1)*NN)= t1
        enddo
        k1line_stop(j+1)= t1

     enddo


     116 continue
     if (.not.lfound .and.RibbonBand_calc) then
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
        if (trim(adjustl(inline))=='KPLANE_BULK') then
           lfound= .true.
           if (cpuid==0) write(stdout, *)' '
           if (cpuid==0) write(stdout, *)'We found KPLANE_BULK card'
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
     real(dp) :: norm

     norm= sqrt(R1(1)*R1(1)+ R1(2)*R1(2)+ R1(3)*R1(3))

     return
  end function norm



   subroutine cart_direct_real(R1, R2)
      use para
      implicit none
      real(dp), intent(in) :: R1(2)
      real(dp), intent(inout) :: R2(2)
      real(dp) :: mata(2, 2)

      mata(1, :)= Rua
      mata(2, :)= Rub

      call inv_r(2, mata)
      R2= R1(1)*mata(1, :)+ R1(2)*mata(2, :)

      return
   end subroutine cart_direct_real

   subroutine direct_cart_real(R1, R2)
      use para
      implicit none
      real(dp), intent(in) :: R1(2)
      real(dp), intent(inout) :: R2(2)

      R2= R1(1)*Rua+ R1(2)*Rub

      return
   end subroutine direct_cart_real
