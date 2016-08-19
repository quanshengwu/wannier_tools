!--------+--------+--------+--------+--------+--------+--------+------!
! main program of a set of tools based on Wannier90 TB
! constructed by Q.S.Wu on 4/9/2010
! change      by Q.S.Wu on 4/22/2010
! changed     by Q.S.wu on July/15/2010
! version     HmnR.data  contains soc
! mpi-version is not test yet, take you own risk.
! mpi-version is  tested , please report bugs to QSWU
! Jan 25 2015 by Q.S.Wu at ETH Zurich 
! wuquansheng@gmail.com
! Copyright (c) 2010 QuanSheng Wu. All rights reserved.
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
     real(dp) :: time_start, time_end, time_init

     !> initial the environment of mpi 
     call mpi_init(ierr)
     call mpi_comm_rank(mpi_cmw,cpuid,ierr)
     call mpi_comm_size(mpi_cmw,num_cpu,ierr)

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

     !> open file Hmn_R.data to get Num_wann and Nrpts
     if (cpuid.eq.0)then
        write(stdout,*)''
        inquire (file =Hrfile, EXIST = exists)
        if (exists)then
           if (index(Hrfile, 'HWR')==0) then
              write(stdout,'(2x,a,a,a)')'File ',trim(Hrfile), &
                 ' exist, We are using HmnR from wannier90'
              open(unit=1001,file=Hrfile,status='old')
              read(1001,*)
              read(1001,'(i)') Num_wann
              read(1001,'(i)') Nrpts
              write(stdout,*)'>> Num_wann', Num_wann 
              write(stdout,*)'>> NRPTS', NRPTS
              close(1001)
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
        if ((soc==0 .and. sum(nprojs)/=Num_wann) .or. &
           (soc>0 .and. sum(nprojs)/=Num_wann/2))then
           print *, sum(nprojs), num_wann, num_wann/2
           stop 'projectors are wrong'
        endif
     endif ! cpuid


     !> broadcast and Nrpts to every cpu
     call MPI_bcast(Num_wann,1,mpi_in,0,mpi_cmw,ierr)
     call MPI_bcast(Nrpts,1,mpi_in,0,mpi_cmw,ierr)
     
     !> dimension for surface green's function
     Ndim= Num_wann* Np


     !> allocate necessary arrays for tight binding hamiltonians
     allocate(irvec(3,nrpts))
     allocate(ndegen(nrpts))
     allocate(HmnR(num_wann,num_wann,nrpts))

     if (cpuid==0)then
        write(stdout,*) ' >> Begin to read Hmn_R.data'
     endif
     call readHmnR() 
     if (cpuid==0)then
        write(stdout,*) ' << Read Hmn_R.data successfully'
     endif


     !> broadcast data to every cpu
     call MPI_bcast(irvec,size(irvec),mpi_in,0,mpi_cmw,ierr)
     call MPI_bcast(HmnR,size(HmnR),mpi_dc,0,mpi_cmw,ierr)
     call MPI_bcast(ndegen,size(ndegen),mpi_in,0,mpi_cmw,ierr)

     !> import symmetry 
    !call symmetry

     !> bulk band
     if (BulkBand_calc) then
        if(cpuid.eq.0)write(stdout, *)' '
        if(cpuid.eq.0)write(stdout, *)'>> Start of calculating bulk band'
        call now(time_start)
        call ek_bulk
       !call ek_bulk_mirror_x
       !call ek_bulk_mirror_z
        call now(time_end)
        call print_time_cost(time_start, time_end, 'BulkBand')
        if(cpuid.eq.0)write(stdout, *)'<< End of calculating bulk band'
     endif

     if (BulkFS_calc) then
        if(cpuid.eq.0)write(stdout, *)' '
        if(cpuid.eq.0)write(stdout, *)'>> Start of calculating the bulk FS'
        call now(time_start)
       !call fermisurface3D
        call fermisurface
        call now(time_end)
        call print_time_cost(time_start, time_end, 'BulkFS')
        if(cpuid.eq.0)write(stdout, *)'<< End of calculating the bulk FS'
     endif

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
           call dos_sub
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



     if (BulkGap_plane_calc) then
        if(cpuid.eq.0)write(stdout, *)' '
        if(cpuid.eq.0)write(stdout, *)'>> Start of calculating the bulk gap in plane'
        call now(time_start)
       !call psik_bulk
       !call ek_bulk_polar
       !call ek_bulk_fortomas
       !call ek_bulk_spin
       !call dos_calc
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
       !call landau_level_k
       !call landau_level_B
        call now(time_end)
        call print_time_cost(time_start, time_end, 'BulkGap_Cube')
        if(cpuid.eq.0)write(stdout, *)'<< End of calculating the bulk gap in Cube'
     endif

     !> slab band
     if (SlabBand_calc)then
        if(cpuid.eq.0)write(stdout, *)' '
        if(cpuid.eq.0)write(stdout, *)'>> Start of calculating the slab band structure'
        call now(time_start)
        call ek_slab
       !call psik     
        call now(time_end)
        call print_time_cost(time_start, time_end, 'SlabBand')
        if(cpuid.eq.0)write(stdout, *)'<< End of calculating the slab band structure'
     endif
 
     if (BerryCurvature_calc)then
        if(cpuid.eq.0)write(stdout, *)' '
        if(cpuid.eq.0)write(stdout, *)'>> Start of calculating the Berry curvature'
        call now(time_start)
        call berry_curvarture 
        call now(time_end)
        call print_time_cost(time_start, time_end, 'BerryCurvature')
        if(cpuid.eq.0)write(stdout, *)'End of calculating the Berry curvature'
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
 
     !> wannier center calculate
     if (wanniercenter_calc)then
        if(cpuid.eq.0)write(stdout, *)' '
        if(cpuid.eq.0)write(stdout, *)'>> Start of calculating the Wilson loop'
        call now(time_start)
        call wannier_center3D_plane
       !call wannier_center2D
       !call wannier_center2D_alt
       !call wannier_center3D
       !call wannier_center3D_plane_mirror_plus
       !call wannier_center3D_plane_mirror_minus
        call now(time_end)
        call print_time_cost(time_start, time_end, 'WannierCenter')
        if(cpuid.eq.0)write(stdout, *)'<< End of calculating the Wilson loop'
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

     !> surface state
     if (SlabSS_calc) then
        if(cpuid.eq.0)write(stdout, *)' '
        if(cpuid.eq.0)write(stdout, *)'>> Start of calculating the surface state'
        call now(time_start)
        call surfstat
        call now(time_end)
        call print_time_cost(time_start, time_end, 'BerryCurvature')
        if(cpuid.eq.0)write(stdout, *)'End of calculating the surface state'
     endif
   
     !> fermi arc
     if (SlabArc_calc) then
        if(cpuid.eq.0)write(stdout, *)' '
        if(cpuid.eq.0)write(stdout, *)'>> Start of calculating the surface arc'
        call now(time_start)
        call fermiarc
        call now(time_end)
        call print_time_cost(time_start, time_end, 'SlabArc')
        if(cpuid.eq.0)write(stdout, *)'End of calculating the surface arc'
     endif
     
     !> calculate spin-texture     
     if (SlabSpintexture_calc)then
        if(cpuid.eq.0)write(stdout, *)' '
        if(cpuid.eq.0)write(stdout, *)'>> Start of calculating the spin texture for surface'
        call now(time_start)
        call spintext
        call now(time_end)
        call print_time_cost(time_start, time_end, 'SlabSpintexture')
        if(cpuid.eq.0)write(stdout, *)'End of calculating the spin texture for surface'
     endif

     call now(time_end)

     if(cpuid.eq.0)write(stdout, *)' '
     call print_time_cost(time_init, time_end, 'whole program')
     if (cpuid.eq.0)write(stdout,*)'Congratulations! you finished the calculation.'
     if (cpuid.eq.0)write(stdout,*)'See you next time :)'

     
     call mpi_finalize(ierr)

  end  !<< end of main program
