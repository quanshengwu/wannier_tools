!--------+--------+--------+--------+--------+--------+--------+------!
! main program of TB on BiSe infinte surface 
! constructed by Q.S.Wu on 4/9/2010
! change      by Q.S.Wu on 4/22/2010
! changed     by Q.S.wu on July/15/2010
! version     HmnR.data  contains soc
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
              write(*,'(2x,a,a,a)')'File ',infilename, &
                 ' exist, We are using HmnR from wannier90'
              open(unit=1001,file=infilename,status='old')
              read(1001,*)
              read(1001,'(i)') Num_wann
              read(1001,'(i)') Nrpts
              write(*,*)'>> Num_wann', Num_wann 
              write(*,*)'>> NRPTS', NRPTS
              close(1001)
           else 
              write(*,'(2x,a,a,a)')'File ',infilename, &
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
     call MPI_bcast(Nrpts,1,mpi_in,0,mpi_cmw,ierr)
	  
	  ! broadcast ndim,Nk,omeganum,Maxomega,nslab,soc,eta to every cpu
     call MPI_bcast(ndim,1,mpi_in,0,mpi_cmw,ierr)
     call MPI_bcast(Nk,1,mpi_in,0,mpi_cmw,ierr)
     call MPI_bcast(nslab,1,mpi_in,0,mpi_cmw,ierr)
     call MPI_bcast(omeganum,1,mpi_in,0,mpi_cmw,ierr)
     call MPI_bcast(soc,1,mpi_in,0,mpi_cmw,ierr)
     call MPI_bcast(eta,1,mpi_dp,0,mpi_cmw,ierr)
     call MPI_bcast(omegamin,1,mpi_dp,0,mpi_cmw,ierr)
     call MPI_bcast(omegamax,1,mpi_dp,0,mpi_cmw,ierr)
     call MPI_bcast(Rua,3,mpi_dp,0,mpi_cmw,ierr)
     call MPI_bcast(Rub,3,mpi_dp,0,mpi_cmw,ierr)
     call MPI_bcast(Ruc,3,mpi_dp,0,mpi_cmw,ierr)
     call MPI_bcast(Kua,3,mpi_dp,0,mpi_cmw,ierr)
     call MPI_bcast(Kub,3,mpi_dp,0,mpi_cmw,ierr)
     call MPI_bcast(Kuc,3,mpi_dp,0,mpi_cmw,ierr)


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

	  call ek_bulk
     call ek_slab
     

	  If (nslab.ge.1)then
	     ! call ek_slab 
	     if(cpuid.eq.0)print *,'begin to calculate surface state'
       !call surfstat
        if(cpuid.eq.0)print *,'end calculate surface state'
        
	     if(cpuid.eq.0)print *,'begin to calculate fermi arc'
       !call fermiarc
        if(cpuid.eq.0)print *,'end calculate fermi arc'
        
       ! calculate spin-texture     
        if(cpuid.eq.0)print *,'begin to calculate spin texture'
        !call spintext
        if(cpuid.eq.0)print *,'end calculate spin texture'
    
     Endif

     if (cpuid.eq.0)write(*,*)'Congratulations! you finished the calculation.'
     
     call mpi_finalize(ierr)

  end
