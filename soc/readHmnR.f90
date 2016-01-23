! read data from HmnR.data   constructed by quansheng wu 4/2/2010

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
     open(12, file=infilename)
     if (.not.index(infilename, 'HWR')) then
        read(12, *)
        read(12, *)nwan
        read(12, *)nrpts
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

     else

        !File *.HWR exist, We are using HmnR from WHM
        ! skip 8 lines
        do i=1,8
           read(12,*)
        enddo
        read(12,'(a11,f18.7)')c_temp,E_fermi
        print*,E_fermi
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
   
     
     !> set up symmetry operators
     !> here we assume that the wannier functions have the symmetry 
     !> as atomic orbitals
    
     nwan= Num_wann/2
     allocate(inversion(Num_wann, Num_wann))
     allocate(mirror_z(Num_wann, Num_wann))
     inversion= 0d0
     mirror_z = 0d0

     !> inversion symmetry
     !> s-> s; p-> -p; d-> -d
     n= 0
     do ia=1, Num_atoms
        do i=1, nprojs(ia)
           n= n+ 1
           select case (proj_name(i, ia))
           case ('s', 'S')
              inversion(n, n)= 1
              inversion(n+ nwan, n+ nwan)= 1
           case ('px', 'Px', 'PX')
              inversion(n, n)=-1
              inversion(n+ nwan, n+ nwan)=-1
           case ('py', 'Py', 'PY')
              inversion(n, n)=-1
              inversion(n+ nwan, n+ nwan)=-1
           case ('pz', 'Pz', 'PZ')
              inversion(n, n)=-1
              inversion(n+ nwan, n+ nwan)=-1
           case ('dxy', 'Dxy', 'DXY')
              inversion(n, n)= 1
              inversion(n+ nwan, n+ nwan)= 1
           case ('dyz', 'Dyz', 'DYZ')
              inversion(n, n)= 1
              inversion(n+ nwan, n+ nwan)= 1
           case ('dxz', 'Dxz', 'DXZ', 'dzx', 'Dzx', 'DZX')
              inversion(n, n)= 1
              inversion(n+ nwan, n+ nwan)= 1
           case ('dx2-y2', 'Dx2-y2', 'DX2-Y2', 'dx2', 'DX2')
              inversion(n, n)= 1
              inversion(n+ nwan, n+ nwan)= 1
           case ('dz2', 'Dz2', 'DZ2')
              inversion(n, n)= 1
              inversion(n+ nwan, n+ nwan)= 1
           case default
              write(*, *) "ERROR: only support s px py pz dxy dyz dxz dx2-y2 dz2 &
                 orbitals"
              stop
           end select
        enddo ! i
     enddo ! ia


     !> mirror_z symmetry
     !> s-> s; px->px, py->py, pz-> -pz
     !> dxy-> dxy, dyz-> -dyz, dxz-> -dxz, dx2-> dx2 dz2->dz2 
     !> up-> up  dn-> -dn
     n= 0
     do ia=1, Num_atoms
        do i=1, nprojs(ia)
           n= n+ 1
           select case (proj_name(i, ia))
           case ('s', 'S')
              mirror_z(n, n)= 1
              mirror_z(n+ nwan, n+ nwan)=-1
           case ('px', 'Px', 'PX')
              mirror_z(n, n)= 1
              mirror_z(n+ nwan, n+ nwan)=-1
           case ('py', 'Py', 'PY')
              mirror_z(n, n)= 1
              mirror_z(n+ nwan, n+ nwan)=-1
           case ('pz', 'Pz', 'PZ')
              mirror_z(n, n)=-1
              mirror_z(n+ nwan, n+ nwan)= 1
           case ('dxy', 'Dxy', 'DXY')
              mirror_z(n, n)= 1
              mirror_z(n+ nwan, n+ nwan)=-1
           case ('dyz', 'Dyz', 'DYZ')

              mirror_z(n, n)=-1
              mirror_z(n+ nwan, n+ nwan)= 1
           case ('dxz', 'Dxz', 'DXZ', 'dzx', 'Dzx', 'DZX')
              mirror_z(n, n)=-1
              mirror_z(n+ nwan, n+ nwan)= 1
           case ('dx2-y2', 'Dx2-y2', 'DX2-Y2', 'dx2', 'DX2')
              mirror_z(n, n)= 1
              mirror_z(n+ nwan, n+ nwan)=-1
           case ('dz2', 'Dz2', 'DZ2')
              mirror_z(n, n)= 1
              mirror_z(n+ nwan, n+ nwan)=-1
           case default
              write(*, *) "ERROR: only support s px py pz dxy dyz dxz dx2-y2 dz2 &
                 orbitals"
              stop
           end select
        enddo ! i
     enddo ! ia




     return
  end subroutine readHmnR
