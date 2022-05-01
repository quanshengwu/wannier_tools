!--------+--------+--------+--------+--------+--------+--------+------!
! Main program of WannierTools based on tight binding model formated
! as wannier90_hr.dat defined in Wannier90 software package.
!
! Ref:
! WannierTools : An open-source software package for novel topological materials
! QuanSheng Wu and ShengNan Zhang and Hai-Feng Song and Matthias Troyer and Alexey A. Soluyanov
! Computer Physics Communications 224, 405 (2018)
! Magnetoresistance from Fermi surface topology, 
! ShengNan Zhang, QuanSheng Wu, Yi Liu, and Oleg V. Yazyev, 
! Phys. Rev. B 99, 035142 (2019)
!
! constructed by Q.S.Wu on 4/9/2010
! change      by Q.S.Wu on 4/22/2010
! changed     by Q.S.wu on July/15/2010
! Jan 25 2015 by Q.S.Wu at ETH Zurich
! version     2.2.1  At EPFL, Switzerland, Sep. 14. 2017
! version     2.4.0  At EPFL, Switzerland, Aug. 31. 2018
! version     2.4.1  At EPFL, Switzerland, Oct. 15. 2018
! version     2.4.2  At EPFL, Switzerland, July. 9. 2019
! version     2.5.0  At EPFL, Switzerland, Dec. 9. 2019, magnetoresistance, band unfolding
! version     2.5.1  At EPFL, Switzerland, Mar. 6. 2020, For WannierTools tutorial 2020
! version     2.6.0  At EPFL, Switzerland, Feb.15. 2021, Landau level, sparse Hamiltonian, TBG
! version     2.6.1  At Beijing, China, April 10. 2022 clean for the Wannier90 tutorial 2022
!
! Corresponding to : wuquansheng@gmail.com
!
! License: GPL V3
!--------+--------+--------+--------+--------+--------+--------+------!

  program main

     use wmpi
     use para
     implicit none

     !> file existence
     logical  :: exists
     integer :: ierr
     character(8) :: cht

     !> time measure
     real(Dp) :: time_start, time_end, time_init


     !> version of WannierTools
     version='2.6.1'

     ierr = 0
     cpuid= 0
     num_cpu= 1
     !> initial the environment of mpi
#if defined (MPI)
     call mpi_init(ierr)
     call mpi_comm_rank(mpi_cmw,cpuid,ierr)
     call mpi_comm_size(mpi_cmw,num_cpu,ierr)
#endif

     if (cpuid==0) open(unit=stdout, file='WT.out')

     !> if mpi initial wrong, alarm
     if (cpuid==0.and.ierr.ne.0)then
        write(stdout,*)'mpi initialize wrong'
        stop
     endif

     call now(time_init)
     call header

     !> print information for mpi
     if (cpuid==0) then
        write(stdout, '(1x, a, i5, a)')'You are using ', num_cpu, ' CPU cores'
        write(stdout, *)' '
     endif

   !> readin the control parameters for this program
   call readinput

   !> determine whether you are using kp model or TB model
   IF (index(KPorTB, 'KP')/=0) then
      Num_wann= sum(Origin_cell%nprojs)
      if (SOC>0) num_wann= 2*num_wann
   ELSE  ! you are using TB model
      !> open file Hmn_R.data to get Num_wann and Nrpts
      if (cpuid.eq.0)then
         write(stdout,*)''
         inquire (file =Hrfile, EXIST = exists)
         if (exists)then
            !> the standard format defined by wannier90
            if (index(Hrfile, 'HWR')==0) then
               write(stdout,'(2x,a,a,a)')'File ',trim(Hrfile), &
                  ' exist, We are using HmnR from wannier90'
               open(unit=1001,file=Hrfile,status='old')
               read(1001,*)
               !> for sparse hr format, this line is the number of non-zero lines
               if(Is_Sparse_Hr) read(1001,*)
               read(1001,*) Num_wann
               read(1001,*) Nrpts
               write(stdout,*)'>> Num_wann', Num_wann
               write(stdout,*)'>> NRPTS', NRPTS
               close(1001)
            !> the format defined by OpenMX
            else
               write(stdout,'(2x, 3a)')'File ',trim(Hrfile), &
                  ' exist, We are using HmnR from HWR'
               open(unit=1001,file=Hrfile,status='old')
               read(1001,*)
               read(1001,'(a26, i3)')cht, Num_wann
               read(1001,'(a32,i10)')cht,Nrpts
               write(stdout,*)'>> Num_wann', Num_wann
               write(stdout,*)'>> NRPTS', NRPTS
               close(1001)
            endif ! hwr or not
         else
            write(stdout,'(2x,a25)')'>>> Error : no HmnR input'
            stop
         endif ! exists or not

         !> We need to extend the spinless hamiltonian to spinfull hamiltonian if
         !> we want to add Zeeman field when spin-orbit coupling is not included in 
         !> the hr file.
         if (Add_Zeeman_Field.and.SOC==0)then
            Num_wann= Num_wann*2
            if (cpuid==0) then
               write(stdout,*)'>> Num_wann is doubled due to the consideration of Zeeman effect'
               write(stdout,*)">> Num_wann : ", Num_wann
            endif
         endif

         if ((soc==0 .and. sum(Origin_cell%nprojs)/=Num_wann .and. .not.Add_Zeeman_Field) .or. &
            (soc>0 .and. sum(Origin_cell%nprojs)/=Num_wann/2))then
            print *, 'sum(Origin_cell%nprojs), num_wann, num_wann/2'
            print *, sum(Origin_cell%nprojs), num_wann, num_wann/2
            print *, "ERROR: Maybe the SOC tags in the SYSTEM is wrongly set"
            stop "ERROR: the summation of all projectors times spin degeneracy is not equal to num_wann"
         endif
      endif ! cpuid

      !> broadcast and Nrpts to every cpu
#if defined (MPI)
     call MPI_bcast(Num_wann,1,mpi_in,0,mpi_cmw,ierr)
     call MPI_bcast(Nrpts,1,mpi_in,0,mpi_cmw,ierr)
#endif
     !> dimension for surface green's function
     Ndim= Num_wann* Np

      if (cpuid==0)then
         write(stdout,*) ' >> Begin to read Hmn_R.data'
      endif

      allocate(irvec(3,nrpts))
      allocate(ndegen(nrpts))
      irvec= 0
      ndegen=1

      if(Is_HrFile) then
         !> allocate necessary arrays for tight binding hamiltonians
         !> normal hmnr file
         if(.not. Is_Sparse_Hr) then
            !> for the dense hr file, we allocate HmnR
            allocate(HmnR(num_wann,num_wann,nrpts))
            HmnR= 0d0
            call readNormalHmnR()
         !> sparse hmnr input
         else
            call readSparseHmnR()
         end if
      else
         stop "We only support Is_HrFile=.true. for this version"
      end if

      if (cpuid==0)then
         write(stdout,*) ' << Read Hmn_R.data successfully'
      endif

     !> broadcast data to every cpu
#if defined (MPI)
     !call MPI_bcast(irvec,size(irvec),mpi_in,0,mpi_cmw,ierr)
     !call MPI_bcast(HmnR,size(HmnR),mpi_dc,0,mpi_cmw,ierr)
     !call MPI_bcast(ndegen,size(ndegen),mpi_in,0,mpi_cmw,ierr)
#endif

     ENDIF  ! end if the choice of kp model or TB model

   !> import symmetry
   call symmetry


   !> unfold bulk band line mode
   if (BulkBand_unfold_line_calc) then
      if(cpuid.eq.0)write(stdout, *)' '
      if(cpuid.eq.0)write(stdout, *)'>> Start of unfolding bulk band'
      call now(time_start)
      call unfolding_kpath
      call now(time_end)
      call print_time_cost(time_start, time_end, 'BulkBand_unfold_line_calc')
      if(cpuid.eq.0)write(stdout, *)'<< End of unfolding bulk band'
   endif


   !> unfold bulk band kplane mode
   if (BulkBand_unfold_plane_calc.or.QPI_unfold_plane_calc) then
      if(cpuid.eq.0)write(stdout, *)' '
      if(cpuid.eq.0)write(stdout, *)'>> Start of unfolding bulk band in plane mode'
      call now(time_start)
      call unfolding_kplane
      call now(time_end)
      call print_time_cost(time_start, time_end, 'BulkBand_unfold_plane_calc')
      if(cpuid.eq.0)write(stdout, *)'<< End of unfolding bulk band in plane mode'
   endif

   !> bulk band
   if (BulkBand_calc.or.BulkBand_line_calc) then
      if(cpuid.eq.0)write(stdout, *)' '
      if(cpuid.eq.0)write(stdout, *)'>> Start of calculating bulk band'
      call now(time_start)
      if(Is_Sparse_Hr) then
#if defined (INTELMKL)
            call sparse_ekbulk
#endif
      else
         call ek_bulk_line
      end if
      call now(time_end)
      call print_time_cost(time_start, time_end, 'BulkBand')
      if(cpuid.eq.0)write(stdout, *)'<< End of calculating bulk band'
   endif

     !> bulk band of a series k points.
     if (BulkBand_points_calc) then
        if(cpuid.eq.0)write(stdout, *)' '
        if(cpuid.eq.0)write(stdout, *)'>> Start of calculating the bulk band in points mode'
        call now(time_start)
        call ek_bulk_point_mode
        call now(time_end)
        call print_time_cost(time_start, time_end, 'BulkBand_points')
        if(cpuid.eq.0)write(stdout, *)'<< End of calculating the bulk band in points mode'
     endif


     !> bulk band in a plane. For Dirac or Weyl cone
   if (BulkBand_plane_calc) then
      if(cpuid.eq.0)write(stdout, *)' '
      if(cpuid.eq.0)write(stdout, *)'>> Start of calculating the bulk band in plane'
      call now(time_start)

      if (Is_Sparse_Hr) then
#if defined (INTELMKL)
         call sparse_ekbulk_plane
#endif
      else
         call ek_bulk_plane
      endif

     !call ek_bulk_plane_C2yT
      call now(time_end)
      call print_time_cost(time_start, time_end, 'BulkBand_plane')
      if(cpuid.eq.0)write(stdout, *)'<< End of calculating the bulk band in plane'
   endif


   if (LandauLevel_B_dos_calc) then
      if(cpuid.eq.0)write(stdout, *)' '
      if(cpuid.eq.0)write(stdout, *)'>> Start of calculating Landau level spectrum'
      call now(time_start)
#if defined (INTELMKL)
      call LandauLevel_B_dos_Lanczos
#endif
      call now(time_end)
      call print_time_cost(time_start, time_end, 'LandauLevel_B_dos_calc')
      if(cpuid.eq.0)write(stdout, *)'<< End of calculating Landau level spectrum'
   endif


   if (LandauLevel_k_dos_calc) then
      if(cpuid.eq.0)write(stdout, *)' '
      if(cpuid.eq.0)write(stdout, *)'>> Start of calculating Landau level spectrum'
      call now(time_start)
#if defined (INTELMKL)
      call LandauLevel_k_dos_Lanczos
#endif
      call now(time_end)
      call print_time_cost(time_start, time_end, 'LandauLevel_k_dos_calc')
      if(cpuid.eq.0)write(stdout, *)'<< End of calculating Landau level spectrum'
   endif
   
   if (Hof_Butt_calc.or.LandauLevel_B_calc) then
      if(cpuid.eq.0)write(stdout, *)' '
      if(cpuid.eq.0)write(stdout, *)'>> start of calculating the Hofstader butterfly '
      call now(time_start)
      if (Is_HrFile) then
         if(Is_Sparse_Hr) then
#if defined (INTELMKL)
            call sparse_landau_level_B
#endif
         else
            call landau_level_B
         end if
      endif
      call now(time_end)
      call print_time_cost(time_start, time_end, 'Hof_Butt_calc')
      if(cpuid.eq.0)write(stdout, *)'<< End of calculating the Hofstader butterfly'
   endif


   !> calculate LandauLevel along kpath with fixed B set by Nslab
   if (LandauLevel_kplane_calc)then
      if(cpuid.eq.0)write(stdout, *)' '
      if(cpuid.eq.0)write(stdout, *)'>> Start to calculate LandauLevel_kplane_calc'
      call now(time_start)
      if (Is_HrFile) then
         call landau_level_kplane
      endif
      call now(time_end)
      call print_time_cost(time_start, time_end, 'LandauLevel_kplane_calc')
      if(cpuid.eq.0)write(stdout, *)'End of LandauLevel_kplane_calc calculation'
   endif


   !> calculate LandauLevel along kpath with fixed B set by Magq
   if (LandauLevel_k_calc)then
      if(cpuid.eq.0)write(stdout, *)' '
      if(cpuid.eq.0)write(stdout, *)'>> Start to calculate LandauLevel_k_calc'
      call now(time_start)
      if (Is_HrFile) then
         if(Is_Sparse_Hr.or.Num_wann*Magq>2000) then
#if defined (INTELMKL)
            call sparse_landau_level_k
#endif
         else
            call landau_level_k
         endif
      else
      endif
      call now(time_end)
      call print_time_cost(time_start, time_end, 'LandauLevel_k_calc')
      if(cpuid.eq.0)write(stdout, *)'End of LandauLevel_k_calc calculation'
   endif


   if (BulkBand_cube_calc) then
      if(cpuid.eq.0)write(stdout, *)' '
      if(cpuid.eq.0)write(stdout, *)'>> Start of calculating the bulk band in a cube of BZ'
      call now(time_start)
      call ek_bulk_cube
      call now(time_end)
      call print_time_cost(time_start, time_end, 'BulkBand_cube')
      if(cpuid.eq.0)write(stdout, *)'<< End of calculating the bulk band in a cube of BZ'
   endif


     !> Find nodes in BZ
     if (FindNodes_calc) then
        if(cpuid.eq.0)write(stdout, *)' '
        if(cpuid.eq.0)write(stdout, *)'>> Start of nodes searching'
        call now(time_start)
        call FindNodes
        call now(time_end)
        call print_time_cost(time_start, time_end, 'FindNodes')
        if(cpuid.eq.0)write(stdout, *)'<< End of nodes searching'
     endif

     !> calculate Fermi surface on a k plane
     if (BulkFS_Plane_calc) then
        if(cpuid.eq.0)write(stdout, *)' '
        if(cpuid.eq.0)write(stdout, *)'>> Start of calculating the bulk FS in a k plane'
        call now(time_start)
        call fermisurface_kplane
        call now(time_end)
        call print_time_cost(time_start, time_end, 'BulkFS_Plane')
        if(cpuid.eq.0)write(stdout, *)'<< End of calculating the bulk FS in a k plane'
     endif

     if (BulkFS_Plane_stack_calc) then
        if(cpuid.eq.0)write(stdout, *)' '
        if(cpuid.eq.0)write(stdout, *)'>> Start of calculating the bulk FS in a k plane stacking'
        call now(time_start)
        call fermisurface_stack
        call now(time_end)
        call print_time_cost(time_start, time_end, 'BulkFS_Plane_stack')
        if(cpuid.eq.0)write(stdout, *)'<< End of calculating the bulk FS in a k plane stacking'
     endif


     !> calculate 3D Fermi surface
     if (BulkFS_calc) then
        if(cpuid.eq.0)write(stdout, *)' '
        if(cpuid.eq.0)write(stdout, *)'>> Start of calculating the bulk FS'
        call now(time_start)
        call fermisurface3D
        call now(time_end)
        call print_time_cost(time_start, time_end, 'BulkFS')
        if(cpuid.eq.0)write(stdout, *)'<< End of calculating the bulk FS'
     endif

     !> calculate density of state and joint density of state
     if (JDos_calc.and.Dos_calc) then
        if(cpuid.eq.0)write(stdout, *)' '
        if(cpuid.eq.0)write(stdout, *)'>> Start of calculating DOS and Jdos for bulk system'
        call now(time_start)
        call dos_joint_dos
        call now(time_end)
        call print_time_cost(time_start, time_end, 'Dos_calc and Jdos_calc')
        if(cpuid.eq.0)write(stdout, *)'<< End of calculating the DOS and Jdos for bulk system'
     else
        if (Dos_calc) then
           if(cpuid.eq.0)write(stdout, *)' '
           if(cpuid.eq.0)write(stdout, *)'>> Start of calculating DOS for bulk system'
           call now(time_start)
           if(.not. Is_Sparse_Hr) then
              call dos_sub
           else
#if defined (INTELMKL)
              call dos_sparse
#endif     
           end if
           call now(time_end)
           call print_time_cost(time_start, time_end, 'Dos_calc')
           if(cpuid.eq.0)write(stdout, *)'<< End of calculating the DOS for bulk system'
        endif

        if (JDos_calc) then
           if(cpuid.eq.0)write(stdout, *)' '
           if(cpuid.eq.0)write(stdout, *)'>> Start of calculating JDOS for bulk system'
           call now(time_start)
           call Joint_dos
           call now(time_end)
           call print_time_cost(time_start, time_end, 'JDos_calc')
           if(cpuid.eq.0)write(stdout, *)'<< End of calculating the JDOS for bulk system'
        endif
     endif

     !> effective mass
     if (EffectiveMass_calc) then
        if(cpuid.eq.0)write(stdout, *)' '
        if(cpuid.eq.0)write(stdout, *)'>> Start of calculating the effective mass'
        call now(time_start)
        call effective_mass_calc
        call now(time_end)
        call print_time_cost(time_start, time_end, 'EffectiveMass_calc')
        if(cpuid.eq.0)write(stdout, *)'<< End of calculating the effective mass'
     endif


     if (BulkGap_plane_calc) then
        if(cpuid.eq.0)write(stdout, *)' '
        if(cpuid.eq.0)write(stdout, *)'>> Start of calculating the bulk gap in plane'
        call now(time_start)
       !call psik_bulk
       !call ek_bulk_polar
       !call ek_bulk_fortomas
       !call ek_bulk_spin
       !call ek_bulk2D
       !call ek_bulk2D_spin
        call gapshape
        call now(time_end)
        call print_time_cost(time_start, time_end, 'BulkGap_plane')
        if(cpuid.eq.0)write(stdout, *)'<< End of calculating the bulk gap in plane'
     endif

     if (BulkGap_Cube_calc) then
        if(cpuid.eq.0)write(stdout, *)' '
        if(cpuid.eq.0)write(stdout, *)'>> start of calculating the bulk gap in Cube'
        call now(time_start)
        call gapshape3D
        call now(time_end)
        call print_time_cost(time_start, time_end, 'BulkGap_Cube')
        if(cpuid.eq.0)write(stdout, *)'<< End of calculating the bulk gap in Cube'
     endif

     !> slab band kplane mode
     if (SlabBand_plane_calc)then
        if(cpuid.eq.0)write(stdout, *)' '
        if(cpuid.eq.0)write(stdout, *)'>> Start of calculating the slab band structure in k plane mode'
        call now(time_start)
        call ek_slab_kplane
        call now(time_end)
        call print_time_cost(time_start, time_end, 'SlabBand_plane')
        if(cpuid.eq.0)write(stdout, *)'<< End of calculating the slab band structure in k plane mode'
     endif

     !> Wave function of Slab 
     if (SlabBandWaveFunc_calc)then
        if(cpuid.eq.0)write(stdout, *)' '
        if(cpuid.eq.0)write(stdout, *)'>> Start of calculating the slab band wave function'
        call now(time_start)
        call psik_slab
        call now(time_end)
        call print_time_cost(time_start, time_end, 'SlabBandWaveFunc_calc')
        if(cpuid.eq.0)write(stdout, *)'<< End of calculating the slab band wave function'
     endif

     !> slab band kpath mode
     if (SlabBand_calc)then
        if(cpuid.eq.0)write(stdout, *)' '
        if(cpuid.eq.0)write(stdout, *)'>> Start of calculating the slab band structure'
        call now(time_start)
        if (Is_Sparse_Hr) then
#if defined (INTELMKL)
          call ek_slab_sparseHR
#endif
        else
           call ek_slab
        endif
        call now(time_end)
        call print_time_cost(time_start, time_end, 'SlabBand_calc')
        if(cpuid.eq.0)write(stdout, *)'<< End of calculating the slab band structure'
     endif

     if (BerryCurvature_slab_calc)then
        if(cpuid.eq.0)write(stdout, *)' '
        if(cpuid.eq.0)write(stdout, *)'>> Start of calculating the Berry curvature for a slab system'
        call now(time_start)
        call berry_curvarture_slab
        call now(time_end)
        call print_time_cost(time_start, time_end, 'BerryCurvature_slab')
        if(cpuid.eq.0)write(stdout, *)'End of calculating the Berry curvature for a slab system'
     endif

     if (Berrycurvature_kpath_EF_calc.or.BerryCurvature_kpath_Occupied_calc)then
        if(cpuid.eq.0)write(stdout, *)' '
        if(cpuid.eq.0)write(stdout, *)'>> Start of calculating the Berry curvature'
        call now(time_start)
        call berry_curvarture_line
        call now(time_end)
        call print_time_cost(time_start, time_end, 'BerryCurvature')
        if(cpuid.eq.0)write(stdout, *)'End of calculating the Berry curvature'
     endif

     if (BerryCurvature_calc.or.Berrycurvature_EF_calc)then
        if(cpuid.eq.0)write(stdout, *)' '
        if(cpuid.eq.0)write(stdout, *)'>> Start of calculating the Berry curvature'
        call now(time_start)
        call berry_curvarture_plane_full
        call now(time_end)
        call print_time_cost(time_start, time_end, 'BerryCurvature')
        if(cpuid.eq.0)write(stdout, *)'End of calculating the Berry curvature'
     endif

     if (BerryCurvature_Cube_calc)then
        if(cpuid.eq.0)write(stdout, *)' '
        if(cpuid.eq.0)write(stdout, *)'>> Start of calculating the Berry curvature in a k-cube'
        call now(time_start)
        call berry_curvarture_cube
        call now(time_end)
        call print_time_cost(time_start, time_end, 'BerryCurvature_Cube')
        if(cpuid.eq.0)write(stdout, *)'End of calculating the Berry curvature in a cube'
     endif


     if (WireBand_calc) then
        if(cpuid.eq.0)write(stdout, *)' '
        if(cpuid.eq.0)write(stdout, *)'>> Start of calculating the wire band'
        call now(time_start)
        call ek_ribbon
        call now(time_end)
        call print_time_cost(time_start, time_end, 'WireBand')
        if(cpuid.eq.0)write(stdout, *)'End of calculating the wire band'
     endif

     !> Chirality of Weyl points calculation
     if (WeylChirality_calc)then
        if(cpuid.eq.0)write(stdout, *)' '
        if(cpuid.eq.0)write(stdout, *)'>> Start of chirality of Weyl points calculating'
        call now(time_start)
        call wannier_center3D_weyl
        call now(time_end)
        call print_time_cost(time_start, time_end, 'WeylChirality_calc')
        if(cpuid.eq.0)write(stdout, *)'<< End of chirality of Weyl points calculating'
     endif


     !> wannier center calculate
     if (wanniercenter_calc)then
        if(cpuid.eq.0)write(stdout, *)' '
        if(cpuid.eq.0)write(stdout, *)'>> Start of calculating the Wilson loop'
        call now(time_start)
        call wannier_center3D_plane_adaptive
        call now(time_end)
        call print_time_cost(time_start, time_end, 'WannierCenter')
        if(cpuid.eq.0)write(stdout, *)'<< End of calculating the Wilson loop'
     endif

     !> mirror chern number calculation
     if (MirrorChern_calc)then
        if(cpuid.eq.0)write(stdout, *)' '
        if(cpuid.eq.0)write(stdout, *)'>> Start of calculating the mirror chern number'
        call now(time_start)
        call wannier_center3D_plane_mirror
        call now(time_end)
        call print_time_cost(time_start, time_end, 'MirrorChern_calc')
        if(cpuid.eq.0)write(stdout, *)'<< End of calculating the mirror chern number'
     endif

     !> wannier center calculattion for the whole BZ, 6 planes
     if (Z2_3D_calc)then
        if(cpuid.eq.0)write(stdout, *)' '
        if(cpuid.eq.0)write(stdout, *)'>> Start of calculating Z2 number for the bulk'
        call now(time_start)
        call Z2_3D_adaptive
        call now(time_end)
        call print_time_cost(time_start, time_end, 'Z2_calc')
        if(cpuid.eq.0)write(stdout, *)'<< End of calculating Z2 number for the bulk'
     endif

     !> wannier center calculattion for the whole BZ, 6 planes
     if (Chern_3D_calc)then
        if(cpuid.eq.0)write(stdout, *)' '
        if(cpuid.eq.0)write(stdout, *)'>> Start of calculating Chern number for the bulk'
        call now(time_start)
        call Chern_3D
        call now(time_end)
        call print_time_cost(time_start, time_end, 'Chern_3D_calc')
        if(cpuid.eq.0)write(stdout, *)'<< End of calculating Chern number for the bulk'
     endif


     if (BerryPhase_calc) then
        if(cpuid.eq.0)write(stdout, *)' '
        if(cpuid.eq.0)write(stdout, *)'>> Start of calculating the Berry phase'
        call now(time_start)
        call berryphase
        call now(time_end)
        call print_time_cost(time_start, time_end, 'BerryPhase')
        if(cpuid.eq.0)write(stdout, *)'End of calculating the Berry phase'
     endif

     !> calculate ordinary hall effect with Boltzmann transport
     if (Boltz_evolve_k)then
        if(cpuid.eq.0)write(stdout, *)' '
        if(cpuid.eq.0)write(stdout, *)'>> Start to calculate ordinary hall effects for different k'
        call now(time_start)
        call evolve_k_ohe
        call now(time_end)
        call print_time_cost(time_start, time_end, 'Boltz_evolve_k')
        if(cpuid.eq.0)write(stdout, *)'End of OHE calculation'
     endif
  
     !> calculate  with Boltzmann transport
     if (Boltz_k_calc)then
        if(cpuid.eq.0)write(stdout, *)' '
        if(cpuid.eq.0)write(stdout, *)'>> Start to calculate ordinary hall effects for different k'
        call now(time_start)
        call sigma_k_ohe
        call now(time_end)
        call print_time_cost(time_start, time_end, 'Boltz_k_calc')
        if(cpuid.eq.0)write(stdout, *)'End of OHE calculation'
     endif
  
     !> calculate ordinary hall effect with Boltzmann transport
     if (Boltz_OHE_calc)then
        if(cpuid.eq.0)write(stdout, *)' '
        if(cpuid.eq.0)write(stdout, *)'>> Start to calculate ordinary hall effects'
        call now(time_start)
        call sigma_resistivity
        call now(time_end)
        call print_time_cost(time_start, time_end, 'Boltz_OHE_calc')
        if(cpuid.eq.0)write(stdout, *)'End of OHE calculation'
     endif


     !> calculate anomalouls hall conductivity
     if (AHC_calc)then
        if(cpuid.eq.0)write(stdout, *)' '
        if(cpuid.eq.0)write(stdout, *)'>> Start to calculate anomalouls hall conductivity'
        call now(time_start)
        call sigma_AHC
        call now(time_end)
        call print_time_cost(time_start, time_end, 'AHC_calc')
        if(cpuid.eq.0)write(stdout, *)'End of AHC calculation'
     endif


     !> surface state
     if (SlabSS_calc) then
        if(cpuid.eq.0)write(stdout, *)' '
        if(cpuid.eq.0)write(stdout, *)'>> Start of calculating the surface state'
        call now(time_start)
        call surfstat
        call now(time_end)
        call print_time_cost(time_start, time_end, 'SlabSS_calc')
        if(cpuid.eq.0)write(stdout, *)'End of calculating the surface state'
     endif

     !> fermi arc
     if (SlabArc_calc) then
        if(cpuid.eq.0)write(stdout, *)' '
        if(cpuid.eq.0)write(stdout, *)'>> Start of calculating the surface arc'
        call now(time_start)
        call SurfaceDOSkk
        call now(time_end)
        call print_time_cost(time_start, time_end, 'SlabArc')
        if(cpuid.eq.0)write(stdout, *)'End of calculating the surface arc'
     endif

     !> fermi arc QPI in kpath mode
     if (SlabQPI_kpath_calc) then
        if(cpuid.eq.0)write(stdout, *)' '
        if(cpuid.eq.0)write(stdout, *)'>> Start of calculating the surface QPI in kpath mode'
        call now(time_start)
        call surfstat_jdos
        call now(time_end)
        call print_time_cost(time_start, time_end, 'SlabQPI')
        if(cpuid.eq.0)write(stdout, *)'End of calculating the surface QPI in kpath mode'
     endif

     !> Surface State QPI in kplane mode
     if (SlabQPI_calc.or.SlabQPI_kplane_calc) then
        if(cpuid.eq.0)write(stdout, *)' '
        if(cpuid.eq.0)write(stdout, *)'>> Start of calculating the surface QPI'
        call now(time_start)
        call SurfaceDOSkk
        call now(time_end)
        call print_time_cost(time_start, time_end, 'SlabQPI')
        if(cpuid.eq.0)write(stdout, *)'End of calculating the surface QPI'
     endif

     !> calculate spin-texture
     if (SlabSpintexture_calc)then
        if(cpuid.eq.0)write(stdout, *)' '
        if(cpuid.eq.0)write(stdout, *)'>> Start of calculating the spin texture for surface'
        call now(time_start)
       !call spintext
        call SurfaceDOSkk
        call now(time_end)
        call print_time_cost(time_start, time_end, 'SlabSpintexture')
        if(cpuid.eq.0)write(stdout, *)'End of calculating the spin texture for surface'
     endif

     !> calculate spin-texture for bulk bands
     if (BulkSpintexture_calc)then
        if(cpuid.eq.0)write(stdout, *)' '
        if(cpuid.eq.0)write(stdout, *)'>> Start of calculating the spin texture for surface'
        call now(time_start)
        call fermisurface_kplane
        call now(time_end)
        call print_time_cost(time_start, time_end, 'SlabSpintexture')
        if(cpuid.eq.0)write(stdout, *)'End of calculating the spin texture for surface'
     endif

     call now(time_end)

     if(cpuid.eq.0)write(stdout, *)' '
     call print_time_cost(time_init, time_end, 'whole program')
     call footer


#if defined (MPI)
     call mpi_finalize(ierr)
#endif

  end  !<< end of main program
