
! read data from HmnR.data   constructed by quansheng wu 4/2/2010

  subroutine readHmnR()

     use para

     implicit none

! file existenct
     logical :: exists

     integer :: i, j, ir
     integer :: i1, i2, i3
     integer :: i4, i5    
     integer :: nwan
     real(dp) :: r1, r2

     !> read in HmnR without spin orbital coupling
     inquire (file =infilename(1), EXIST = exists)
     if (exists)then
        open(12, file=infilename(1), status='old')
        read(12, *)
        read(12, *)nwan
        read(12, *)Nrpts_nsoc
        Num_wann_nsoc= nwan
        allocate(irvec_nsoc(3,nrpts_nsoc))
        allocate(ndegen_nsoc(nrpts_nsoc))
        allocate(HmnR_nsoc(num_wann_nsoc*2,num_wann_nsoc*2,nrpts_nsoc))
        allocate(HmnR_nsoc_origin(num_wann_nsoc*2,num_wann_nsoc*2,nrpts_nsoc))
        HmnR_nsoc_origin= 0d0
        HmnR_nsoc= 0d0
        ndegen_nsoc= 1
        irvec_nsoc= 0
        read(12,*)(ndegen_nsoc(i), i=1, nrpts_nsoc)
        do ir=1, nrpts_nsoc
           do i=1, nwan
              do j=1, nwan
                 read(12,*)i1, i2, i3, i4, i5, r1, r2
                 HmnR_nsoc_origin(i4,i5,ir)= dcmplx(r1, r2) ! in eV
                 !write(*,'(5i5,2f10.5)')i1, i2, i3, i4, i5, r1, r2
              enddo
           enddo
           irvec_nsoc(1, ir)=i1
           irvec_nsoc(2, ir)=i2
           irvec_nsoc(3, ir)=i3
        enddo
        close(12)

        do ir=1, nrpts_nsoc
           HmnR_nsoc_origin(1+Num_wann_nsoc:2*Num_wann_nsoc, &
                            1+Num_wann_nsoc:2*Num_wann_nsoc, ir)= &
           HmnR_nsoc_origin(1:Num_wann_nsoc, 1:Num_wann_nsoc, ir)
        enddo

        HmnR_nsoc= HmnR_nsoc_origin
     endif
   
     !> read in HmnR with spin orbital coupling
     inquire (file =infilename(2), EXIST = exists)
     if (exists)then
        open(12, file=infilename(2), status='old')
        read(12, *)
        read(12, *)nwan
        read(12, *)Nrpts_soc
        Num_wann_soc= nwan
        allocate(irvec_soc(3,nrpts_soc))
        allocate(ndegen_soc(nrpts_soc))
        allocate(HmnR_soc(num_wann_soc,num_wann_soc,nrpts_soc))
        HmnR_soc= 0d0
        ndegen_soc= 1
        irvec_soc= 0
        read(12,*)(ndegen_soc(i), i=1, nrpts_soc)
        do ir=1, nrpts_soc
           do i=1, nwan
              do j=1, nwan
                 read(12,*)i1, i2, i3, i4, i5, r1, r2
                 HmnR_soc(i4,i5,ir)= dcmplx(r1, r2) ! in eV
                 !write(*,'(5i5,2f10.5)')i1, i2, i3, i4, i5, r1, r2
              enddo
           enddo
           irvec_soc(1, ir)=i1
           irvec_soc(2, ir)=i2
           irvec_soc(3, ir)=i3
        enddo
        close(12)
        do ir=1, nrpts_soc
           if (sum(abs(irvec_soc(:, ir)))<0.1) then
              do i=1, nwan
                 HmnR_soc(i, i, ir)= HmnR_soc(i, i, ir)- E_fermi
              enddo
           endif
        enddo
     endif
           
     write(*, *)' read hmnr file sucessfully'
     return
  end subroutine readHmnR
