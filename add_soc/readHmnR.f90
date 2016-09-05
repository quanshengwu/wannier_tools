
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

     !> for InSb soc-hyb 666

     inquire (file =infilename, EXIST = exists)
     if (exists)then
        open(12, file=infilename, status='old')
        read(12, *)
        read(12, *)nwan
        read(12, *)nrpts
        num_wann= nwan
        allocate(irvec(3,nrpts))
        allocate(ndegen(nrpts))
        allocate(HmnR(num_wann*2,num_wann*2,nrpts))
        read(12,*)(ndegen(i), i=1, nrpts)
        do ir=1, nrpts
           do i=1, nwan
              do j=1, nwan
                 read(12,*)i1, i2, i3, i4, i5, r1, r2
                 HmnR(i4,i5,ir)= dcmplx(r1, r2) ! in eV
                 HmnR(i4+nwan,i5+nwan,ir)= dcmplx(r1,-r2) ! in eV
                 !write(*,'(5i5,2f10.5)')i1, i2, i3, i4, i5, r1, r2
              enddo
           enddo
           irvec(1, ir)=i1
           irvec(2, ir)=i2
           irvec(3, ir)=i3
        enddo
        close(12)
        do ir=1, nrpts
           if (sum(abs(irvec(:, ir)))<0.1) then
              do i=1, nwan*2
                 HmnR(i, i, ir)= HmnR(i, i, ir)- E_fermi
              enddo
           endif
        enddo
     endif
         
     return
  end subroutine readHmnR


! read data from HmnR.data   constructed by quansheng wu 4/2/2010

  subroutine readHmnR_updn()

     use para

     implicit none

! file existenct
     logical :: exists

     integer :: i, j, ir
     integer :: i1, i2, i3
     integer :: i4, i5    
     integer :: nwan
     real(dp) :: r1, r2

     !> for InSb soc-hyb 666

     inquire (file ='wannier90_hr.up.dat', EXIST = exists)
     if (.not.exists) stop 'wannier90_hr.up.dat is not existed'
     inquire (file ='wannier90_hr.dn.dat', EXIST = exists)
     if (.not.exists) stop 'wannier90_hr.dn.dat is not existed'

     if (exists)then
        open(12, file='wannier90_hr.up.dat', status='old')
        open(13, file='wannier90_hr.dn.dat', status='old')
        read(12, *)
        read(12, *)nwan
        read(12, *)nrpts
        read(13, *)
        read(13, *)nwan
        read(13, *)nrpts
        num_wann= nwan
        allocate(irvec(3,nrpts))
        allocate(ndegen(nrpts))
        allocate(HmnR(num_wann*2,num_wann*2,nrpts))
        read(12,*)(ndegen(i), i=1, nrpts)
        read(13,*)(ndegen(i), i=1, nrpts)
        do ir=1, nrpts
           do i=1, nwan
              do j=1, nwan
                 read(12,*)i1, i2, i3, i4, i5, r1, r2
                 HmnR(i4,i5,ir)= dcmplx(r1, r2) ! in eV
                 read(13,*)i1, i2, i3, i4, i5, r1, r2
                 HmnR(i4+nwan,i5+nwan,ir)= dcmplx(r1, r2) ! in eV
                 !write(*,'(5i5,2f10.5)')i1, i2, i3, i4, i5, r1, r2
              enddo
           enddo
           irvec(1, ir)=i1
           irvec(2, ir)=i2
           irvec(3, ir)=i3
        enddo
        close(12)
        do ir=1, nrpts
           if (sum(abs(irvec(:, ir)))<0.1) then
              do i=1, nwan*2
                 HmnR(i, i, ir)= HmnR(i, i, ir)- E_fermi
              enddo
           endif
        enddo
     endif
         
     return
  end subroutine readHmnR_updn

! read data from HmnR.data   constructed by quansheng wu 4/2/2010

  subroutine readHmnR_2d()

     use para

     implicit none

! file existenct
     logical :: exists

     integer :: i, j, ir
     integer :: i1, i2, i3
     integer :: i4, i5    
     integer :: nwan
     real(dp) :: r1, r2

     !> for InSb soc-hyb 666

     inquire (file =infilename, EXIST = exists)
     if (exists)then
        open(12, file=infilename, status='old')
        read(12, *)
        read(12, *)nwan
        read(12, *)nrpts
        num_wann= nwan
        allocate(irvec_2d(2,nrpts))
        allocate(ndegen(nrpts))
        allocate(HmnR(num_wann*2,num_wann*2,nrpts))
        read(12,*)(ndegen(i), i=1, nrpts)
        do ir=1, nrpts
           do i=1, nwan
              do j=1, nwan
                 read(12,*)i1, i2, i4, i5, r1, r2
                 HmnR(i4,i5,ir)= dcmplx(r1, r2) ! in eV
                 HmnR(i4+nwan,i5+nwan,ir)= dcmplx(r1,-r2) ! in eV
                 !write(*,'(5i5,2f10.5)')i1, i2, i3, i4, i5, r1, r2
              enddo
           enddo
           irvec_2d(1, ir)=i1
           irvec_2d(2, ir)=i2
        enddo
        close(12)
        do ir=1, nrpts
           if (sum(abs(irvec_2d(:, ir)))<0.1) then
              do i=1, nwan*2
                 HmnR(i, i, ir)= HmnR(i, i, ir)- E_fermi
              enddo
           endif
        enddo
     endif
         
     return
  end subroutine readHmnR_2d
