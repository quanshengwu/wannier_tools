! read data from HmnR.data   constructed by quansheng wu 4/2/2010
! Copyright (c) 2010 QuanSheng Wu. All rights reserved.

  subroutine readHmnR()

     use para

     implicit none

     character*4 :: c_temp

! file existenct
     logical :: exists

     integer :: i, j, ir, ia
     integer :: n, m
     integer :: i1, i2, i3
     integer :: i4, i5    
     integer :: nwan
     real(dp) :: r1, r2

! real and imag of HmnR
     real(dp) :: rh,ih
     
     if(cpuid.eq.0)write(stdout,*)''
     open(12, file=Hrfile)
     if (index(Hrfile, 'HWR')==0) then
        read(12, *)
        read(12, *)nwan
        read(12, *)nrpts
        read(12,*)(ndegen(i), i=1, nrpts)
        do ir=1, nrpts
           do i=1, nwan
              do j=1, nwan
                 read(12,*)i1, i2, i3, i4, i5, r1, r2
                !read(12 , '(5I5, 2f12.3)')i1, i2, i3, i4, i5, r1, r2
                !r1=(anint(r1*100000))/100000
                !r2=(anint(r2*100000))/100000
                 HmnR(i4,i5,ir)= dcmplx(r1, r2) ! in eV
                 !write(*,'(5i5,2f10.5)')i1, i2, i3, i4, i5, r1, r2
              enddo
           enddo
           irvec(1, ir)=i1
           irvec(2, ir)=i2
           irvec(3, ir)=i3
        enddo

     else

        !File *.HWR exist, We are using HmnR from WHM
        ! skip 8 lines
        do i=1,8
           read(12,*)
        enddo
        read(12,'(a11,f18.7)')c_temp,E_fermi
        do iR=1,Nrpts
           read(12,'(a3,3i5,a3,i4)')c_temp,irvec(1:3,iR),c_temp,ndegen(iR)
           do i=1, Num_wann*Num_wann
              read(12,*)n,m,rh,ih
              HmnR(n,m,iR)=rh+ zi*ih   ! in Hartree
           enddo
           if (sum(abs(irvec(:, ir)))==0) then
              do i=1, Num_wann
                 HmnR(i,i,iR)= HmnR(i,i,iR)-E_fermi  ! in Hartree
              enddo
           endif
        enddo
        HmnR= HmnR*27.2114d0

        nwan= Num_wann
        if (cpuid==0) then
           open(unit=105, file='wannier90_hr.dat')
           write(105, *)'hr file transformed from HWR'
           write(105, *)Nwan
           write(105, *)nrpts
           write(105, '(15I5)')(ndegen(i), i=1, nrpts)
           do ir=1, nrpts
              do i=1, Nwan
                 do j=1, Nwan
                    write(105, '(5I5, 2f16.8)')irvec(:, ir), i, j, HmnR(i, j, ir)
                 enddo
              enddo
           enddo
           close(105)
        endif


     endif ! HWR or not
     close(12)

    !call get_fermilevel

     do iR=1,Nrpts
        if (Irvec(1,iR).eq.0.and.Irvec(2,iR).eq.0.and.Irvec(3,iR).eq.0)then
           do i=1, Num_wann
              HmnR(i,i,iR)=HmnR(i,i,iR)-E_fermi
           enddo
        endif
     enddo
   
     return
  end subroutine readHmnR
