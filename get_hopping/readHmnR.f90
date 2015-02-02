
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
        allocate(HmnR(num_wann,num_wann,nrpts))
        read(12,*)(ndegen(i), i=1, nrpts)
        do ir=1, nrpts
           do i=1, nwan
              do j=1, nwan
                 read(12,*)i1, i2, i3, i4, i5, r1, r2
                 HmnR(i4,i5,ir)= dcmplx(r1, r2) ! in eV
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
              do i=1, nwan
                 HmnR(i, i, ir)= HmnR(i, i, ir)- E_fermi
              enddo
           endif
        enddo
     endif
         
     return
  end subroutine readHmnR
