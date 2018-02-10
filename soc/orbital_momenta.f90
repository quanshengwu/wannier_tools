  subroutine orbital_momenta(k, eigvec)
     ! Calculate the orbital momenta for each band at Gamma point
     ! 
     ! 2016 by QuanSheng Wu @ ETHZ
     !
     ! wuquansheng@gmail.com
     use para
     implicit none

     real(dp), intent(in) :: k(2)
     complex(dp), intent(in) :: eigvec(Num_wann, Num_wann)

     complex(dp), allocatable :: Lx(:, :)
     complex(dp), allocatable :: Ly(:, :)
     complex(dp), allocatable :: Lz(:, :)
     complex(dp), allocatable :: sx(:, :)
     complex(dp), allocatable :: sy(:, :)
     complex(dp), allocatable :: sz(:, :)
     complex(dp), allocatable :: jx(:, :)
     complex(dp), allocatable :: jy(:, :)
     complex(dp), allocatable :: jz(:, :)
     complex(dp), allocatable :: j_tot(:, :)
     complex(dp), allocatable :: jtemp(:, :)
     complex(dp), allocatable :: eigvec_dag(:, :)
     complex(dp), allocatable :: sub(:, :)
     real(dp), allocatable :: eig2(:)

     integer :: i, j
     integer :: nwann

     return
     if (SOC==0) return

     nwann=Num_wann/2
     if (nwann.ne.7) return

     allocate(Lx(Num_wann, Num_wann))
     allocate(Ly(Num_wann, Num_wann))
     allocate(Lz(Num_wann, Num_wann))
     allocate(sx(Num_wann, Num_wann))
     allocate(sy(Num_wann, Num_wann))
     allocate(sz(Num_wann, Num_wann))
     allocate(jx(Num_wann, Num_wann))
     allocate(jy(Num_wann, Num_wann))
     allocate(jz(Num_wann, Num_wann))
     allocate(j_tot(Num_wann, Num_wann))
     allocate(jtemp(Num_wann, Num_wann))
     allocate(eigvec_dag(Num_wann, Num_wann))
     allocate(sub(2, 2))
     allocate(eig2(2))
     Lx= 0d0
     Ly= 0d0
     Lz= 0d0
     sx= 0d0
     sy= 0d0
     sz= 0d0
     jx= 0d0
     jtemp= 0d0
     jy= 0d0
     jz= 0d0
     j_tot= 0d0
     sub= 0d0
     eig2= 0d0

     if (Num_wann<7) return

     !> only for AlSb/InAs/GaSb/AlSb system
     !>> p orbital operators
     !> spin up
     Lz(   3,   4)= -zi
     Lz(   4,   3)=  zi
     Lz(   6,   7)= -zi
     Lz(   7,   6)=  zi

     Lx(   2,   4)=  zi
     Lx(   4,   2)= -zi
     Lx(   5,   7)=  zi
     Lx(   7,   5)= -zi

     Ly(   2,   3)= -zi
     Ly(   3,   2)=  zi
     Ly(   5,   6)= -zi
     Ly(   6,   5)=  zi

     !> spin down
     Lz(   3+ nwann,   4+ nwann)= -zi
     Lz(   4+ nwann,   3+ nwann)=  zi
     Lz(   6+ nwann,   7+ nwann)= -zi
     Lz(   7+ nwann,   6+ nwann)=  zi

     Lx(   2+ nwann,   4+ nwann)=  zi
     Lx(   4+ nwann,   2+ nwann)= -zi
     Lx(   5+ nwann,   7+ nwann)=  zi
     Lx(   7+ nwann,   5+ nwann)= -zi

     Ly(   2+ nwann,   3+ nwann)= -zi
     Ly(   3+ nwann,   2+ nwann)=  zi
     Ly(   5+ nwann,   6+ nwann)= -zi
     Ly(   6+ nwann,   5+ nwann)=  zi

    !Lx=  Lx
    !Ly=  Ly
    !Lz=  Lz

    !write(*, *)'Lx'
     do i=1, Num_wann
    !   write(*, '(200i3)')int8(real(Lx(i, :)))
     enddo

    !write(*, *)'Ly'
     do i=1, Num_wann
    !   write(*, '(200i3)')int8(real(Ly(i, :)))
     enddo

    !write(*, *)'Lz'
     do i=1, Num_wann
    !   write(*, '(200i3)')int8(real(Lz(i, :)))
     enddo

     !>> spin operators
     do i=1, nwann
        !> spin up
        sx(   i,          i+ nwann)= 1d0
        sx(   i+ nwann,   i       )= 1d0

        sy(   i,          i+ nwann)= -zi
        sy(   i+ nwann,   i       )=  zi

        sz(   i,          i       )= 1d0
        sz(   i+ nwann,   i+ nwann)= -1d0
     enddo

    !write(*, *)'Sx'
     do i=1, Num_wann
    !   write(*, '(200i3)')int8(real(Sx(i, :)))
     enddo

    !write(*, *)'Sy'
     do i=1, Num_wann
    !   write(*, '(200i3)')int8(real(Sy(i, :)))
     enddo

    !write(*, *)'Sz'
     do i=1, Num_wann
    !   write(*, '(200i3)')int8(real(Sz(i, :)))
     enddo


     !> j=l+s/2
     jx= 1*Lx+ 1*sx/2
     jy= 1*Ly+ 1*sy/2
     jz= 1*Lz+ 1*sz/2

     call mat_mul(Num_wann, jx, jx, jtemp)
     j_tot= jtemp
     call mat_mul(Num_wann, jy, jy, jtemp)
     j_tot= jtemp+ j_tot
     call mat_mul(Num_wann, jz, jz, jtemp)
     j_tot= jtemp+ j_tot

     jx= 1*Lx+ 1*sx/2
     jy= 1*Ly+ 1*sy/2
     jz= 1*Lz+ 1*sz/2

     eigvec_dag= conjg(transpose(eigvec))

     call mat_mul(Num_wann, jx, eigvec, jtemp)
     call mat_mul(Num_wann, eigvec_dag, jtemp, jx)

     !> due to the degeneracy of the bands, so we need to diagonalize the block matrix 
     !> for each degenerate bands
     do i=1, Num_wann/2
        sub(1, 1)= jx(2*i-1, 2*i-1)
        sub(1, 2)= jx(2*i-1, 2*i)
        sub(2, 1)= jx(2*i, 2*i-1)
        sub(2, 2)= jx(2*i, 2*i)
        call eigensystem_c('N', 'U', 2, sub, eig2)
        jx(2*i-1, 2*i-1)= eig2(1)
        jx(2*i, 2*i)= eig2(2)
     enddo


     call mat_mul(Num_wann, jy, eigvec, jtemp)
     call mat_mul(Num_wann, eigvec_dag, jtemp, jy)

     do i=1, Num_wann/2
        sub(1, 1)= jy(2*i-1, 2*i-1)
        sub(1, 2)= jy(2*i-1, 2*i)
        sub(2, 1)= jy(2*i, 2*i-1)
        sub(2, 2)= jy(2*i, 2*i)
        call eigensystem_c('N', 'U', 2, sub, eig2)
        jy(2*i-1, 2*i-1)= eig2(1)
        jy(2*i, 2*i)= eig2(2)
     enddo
  
     call mat_mul(Num_wann, jz, eigvec, jtemp)
     call mat_mul(Num_wann, eigvec_dag, jtemp, jz)

    !write(*, *)'jz'
     do i=1, Num_wann
    !   write(*, '(200f4.1)')(real(jz(i, :)))
     enddo



     do i=1, Num_wann/2
        sub(1, 1)= jz(2*i-1, 2*i-1)
        sub(1, 2)= jz(2*i-1, 2*i)
        sub(2, 1)= jz(2*i, 2*i-1)
        sub(2, 2)= jz(2*i, 2*i)
        call eigensystem_c('V', 'U', 2, sub, eig2)
        jz(2*i-1, 2*i-1)= eig2(1)
        jz(2*i, 2*i-1)= 0d0
        jz(2*i-1, 2*i)= 0d0
        jz(2*i, 2*i)= eig2(2)
      enddo

 
     call mat_mul(Num_wann, j_tot, eigvec, jtemp)
     call mat_mul(Num_wann, eigvec_dag, jtemp, j_tot)

    !write(*, *)'j_tot'
     do i=1, Num_wann
    !   write(*, '(200f4.1)')(real(j_tot(i, :)))
     enddo



     do i=1, Num_wann/2
        sub(1, 1)= j_tot(2*i-1, 2*i-1)
        sub(1, 2)= j_tot(2*i-1, 2*i)
        sub(2, 1)= j_tot(2*i, 2*i-1)
        sub(2, 2)= j_tot(2*i, 2*i)
        call eigensystem_c('V', 'U', 2, sub, eig2)
        j_tot(2*i-1, 2*i-1)= eig2(1)
        j_tot(2*i, 2*i)= eig2(2)
      enddo

    !write(*, *)'j_tot'
     do i=1, Num_wann
    !   write(*, '(200f4.1)')(real(j_tot(i, :)))
     enddo



     if (cpuid.eq.0) then
        open(unit=53243, file='orbital.txt')
        write(53243, *)"# top valance band : ", Numoccupied
        write(53243, '(a, 2f10.4)')"# k point :", k
        write(53243, '(a6, 10a8)')"# Band ", 'jx', 'jy', 'jz', 'j_tot'
        do i=1, Num_wann
           write(53243, '(i6, 3000f8.3)') i, real(jx(i, i)), real(jy(i, i)), real(jz(i, i)), &
          !(real(jx(i, i)))**2+ (real(jy(i, i)))**2+ (real(jz(i, i)))**2
           (real(j_tot(i, i)) )
        enddo
        close(53243)

     endif

     return
  end subroutine orbital_momenta
