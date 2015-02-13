!--------+--------+--------+--------+--------+--------+--------+------!
! main program of TB on BiSe infinte surface 
! constructed by Q.S.Wu on 4/9/2010
! change      by Q.S.Wu on 4/22/2010
! changed     by Q.S.wu on July/15/2010
! version     HmnR.data  contains soc
! mpi-version is not test yet, take you own risk.
! Jan 25 2015 by Q.S.Wu at ETH Zurich 
! 
!--------+--------+--------+--------+--------+--------+--------+------!

  program main

    use para
    use mpi
    implicit none

    ! file existence
    logical  :: exists

	  character*4 :: c_temp
	  integer :: i_temp
	  real(dp) :: r_temp

     integer :: namelen
     integer :: ierr

     ! initial the environment of mpi 
     call mpi_init(ierr)
     call mpi_comm_rank(mpi_cmw,cpuid,ierr)
     call mpi_comm_size(mpi_cmw,num_cpu,ierr)

     ! if mpi initial wrong, alarm 
     if (cpuid==0.and.ierr.ne.0)then
        write(*,*)'mpi initialize wrong'
        stop 
     endif

     call readinput

     ! open file Hmn_R.data to get Num_wann and Nrpts
     if (cpuid.eq.0)then
        write(*,*)''
        inquire (file =infilename, EXIST = exists)
        if (exists)then
           if (.not.index(infilename, 'HWR')) then
              write(*,'(2x,a,a,a)')'File ',trim(infilename), &
                 ' exist, We are using HmnR from wannier90'
              open(unit=1001,file=infilename,status='old')
              read(1001,*)
              read(1001,'(i)') Num_wann
              read(1001,'(i)') Nrpts
              write(*,*)'>> Num_wann', Num_wann 
              write(*,*)'>> NRPTS', NRPTS
              close(1001)
           else 
              write(*,'(2x, 3a)')'File ',trim(infilename), &
                 ' exist, We are using HmnR from HWR'
              open(unit=1001,file=infilename,status='old')
              read(1001,*)
              read(1001,'(a26, i3)') c_temp, Num_wann
              read(1001,'(a32,i10)')c_temp,Nrpts
              write(*,*)'>> Num_wann', Num_wann 
              write(*,*)'>> NRPTS', NRPTS
              close(1001)
           endif ! hwr or not
        else
           write(*,'(2x,a25)')'>>> Error : no HmnR input'
           stop
        endif ! exists or not
     endif ! cpuid


	  ! broadcast and Nrpts to every cpu
     call MPI_bcast(Num_wann,1,mpi_in,0,mpi_cmw,ierr)
     call MPI_bcast(Nrpts,1,mpi_in,0,mpi_cmw,ierr)
	  
     !> dimension for surface green's function
     Ndim= Num_wann* Np


	  allocate(irvec(3,nrpts))
	  allocate(ndegen(nrpts))
     allocate(HmnR(num_wann,num_wann,nrpts))

     if(cpuid==0)then
       write(*,*) 'begin reading Hmn_R.data'
       call readHmnR() 
       write(*,*) 'read Hmn_R.data successfully'
     endif


	  ! broadcast data to every cpu
     call MPI_bcast(irvec,size(irvec),mpi_in,0,mpi_cmw,ierr)
     call MPI_bcast(HmnR,size(HmnR),mpi_dc,0,mpi_cmw,ierr)
     call MPI_bcast(ndegen,size(ndegen),mpi_in,0,mpi_cmw,ierr)

     !> bulk band
	  if(cpuid.eq.0)print *,'begin to calculate bulk band'
     if (BulkBand_calc)call ek_bulk
     if(cpuid.eq.0)print *,'end calculate bulk band'

     !> slab band
     if (SlabBand_calc)call ek_slab
     

     !> surface state
	  if(cpuid.eq.0)print *,'begin to calculate surface state'
     if (SlabSS_calc)call surfstat
     if(cpuid.eq.0)print *,'end calculate surface state'
    
     !> fermi arc
	  if(cpuid.eq.0)print *,'begin to calculate fermi arc'
     if (SlabArc_calc)call fermiarc
     if(cpuid.eq.0)print *,'end calculate fermi arc'
     
     !> calculate spin-texture     
     if(cpuid.eq.0)print *,'begin to calculate spin texture'
     if (SlabSpintexture_calc)call spintext
     if(cpuid.eq.0)print *,'end calculate spin texture'
   
     if (cpuid.eq.0)write(*,*)'Congratulations! you finished the calculation.'
     
     call mpi_finalize(ierr)

  end
