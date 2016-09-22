
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
     integer :: order(14)
     real(dp) :: r1, r2
     complex(dp), allocatable :: HmnR_new(:, :, :)

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
        allocate(HmnR_new(num_wann,num_wann,nrpts))
        read(12,*)(ndegen(i), i=1, nrpts)
        do ir=1, nrpts
           do i=1, nwan
              do j=1, nwan
                 read(12,*)i1, i2, i3, i4, i5, r1, r2
                 HmnR(i4,i5,ir)= dcmplx(r1, r2)/ndegen(ir) ! in eV
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

     !> reorder the orbital 
     !> in my convention, I use s px py pz. But in Georg's convention, 
     !> he use s py pz px. 
     !> for InAs, the symmetrized hamiltonian has the order like this
     !> s, pz, px, py, pz, px, py
     !> I need change it to 
     !> s, px, py, pz, px, py, pz
     order(1)=1
     order(2)=3
     order(3)=4
     order(4)=2
     order(5)=6
     order(6)=7
     order(7)=5
     order(8)=8
     order(9)=10
     order(10)=11
     order(11)=9 
     order(12)=13
     order(13)=14
     order(14)=12

     if (num_wann/=14) stop 'only works for InAs with 14 orbitals'
     do ir=1, Nrpts
        do i=1, num_wann
        do j=1, num_wann
    !      HmnR_new(i, j, ir)= HmnR(order(i), order(j), ir)
        enddo
        enddo
     enddo

    !HmnR=HmnR_new
         
     return
  end subroutine readHmnR
