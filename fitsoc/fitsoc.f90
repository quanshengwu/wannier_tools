!  > using minimization algorithm to search spin-orbital strength 
   subroutine fitsoc

      use para
      implicit none

      integer :: i
      integer :: j
      integer :: l
      integer :: n
      integer :: ik
      real(dp) :: k(3)
      real(dp), allocatable :: x(:)

      integer :: npara

      integer  :: maxfun

      ! output frequency
      integer  :: iprint

      integer  :: iexit

      integer  :: Nitermax
      integer  :: iter

      real(dp) :: ftol

      real(dp), allocatable :: W(:)
      complex(dp), allocatable :: Hamk_soc(:, :)

      real(Dp),allocatable :: y(:)
      real(Dp),allocatable :: p(:, :)

      interface
         function func(n, x)
         implicit none
         integer, parameter :: dp= kind(1d0)
         integer, intent(in) :: n
         real(dp), dimension(n), intent(in) :: x
         real(dp) :: func
         end function func
      end interface

      allocate(W(num_wann_soc))
      allocate(Hamk_soc(num_wann_soc, num_wann_soc))
      W=0d0
      Hamk_soc= 0d0

      n= 3


      knv3= Nk*Nk*Nk
      allocate(kpoints(3, knv3))
      kpoints= 0d0

      ik= 0
      do i=1, Nk
         do j=1, Nk
            do l=1, Nk
               ik= ik+ 1
               kpoints(1, ik)= dble(i)/Nk
               kpoints(2, ik)= dble(j)/Nk
               kpoints(3, ik)= dble(l)/Nk
            enddo !l 
         enddo !j 
      enddo !i 


      !> using kline mode
      deallocate(kpoints)
      allocate(kpoints(3, nk3_band))
      kpoints= k3points

      knv3= nk3_band
      num_bands_DFT= 82

      allocate(eigval_soc(num_wann_soc, knv3))
      allocate(eigval_soc_DFT(num_bands_DFT, knv3))
      allocate(eigval_nsoc(num_wann_soc, knv3))
      eigval_soc= 0d0
      eigval_soc_DFT= 0d0
     
      !> get eigenvalues for soc hamiltonian
     !do ik=1, knv3
     !   k= kpoints(:, ik)
     !   W=0d0
     !   Hamk_soc= 0d0
     !   call Hamk_bulk_soc(k, Hamk_soc)
     !   call eigensystem_c('N', 'U', num_wann_soc, Hamk_soc, W)
     !   eigval_soc(:, ik)= W
     !enddo

      !> get eigenvalues for soc hamiltonian from DFT
      call EigvalFromDFT
     !eigval_soc=eigval_soc_DFT(19:18+num_wann_soc, :)
      eigval_soc(1:num_bands_DFT-17+1, :)=eigval_soc_DFT(17:num_bands_DFT, :)

      npara= Num_atom_type*2+ 1
 
      allocate(x(npara))
      allocate(y(npara+1))
      allocate(p(npara+1, npara))
      x=0d0
      y=0d0
      p=0d0

      !> initial value for spin-orbital coupling
      x(1:Num_atom_type)= lambda_p(1:Num_atom_type)
      x(1+Num_atom_type:2*Num_atom_type)= lambda_d(1:Num_atom_type)
      x(Num_atom_type*2+ 1)= E_fermi

      Nitermax=10000
      do i=1, npara+1
         p(i, :)= x
      enddo
      do i=1, npara
         p(i+1, i)= p(i+1, i)+0.02d0*i
      enddo

      do i=1, npara+1
         y(i)= func(npara, p(i, :))
      enddo

      ftol=1d-6 
      call amoeba(npara, p, y, ftol, func, iter)

      x= p(1, :)
      write(*,*)'  Iprint',iprint, 'iter', iter
      write(*, '(3f10.4)')x
      write(*, *)' '

      open(unit=14, file='bulkek-soc.dat')
     !do i=1, num_wann_soc
      do i=1, num_bands_DFT-16
         do ik=1, knv3
            write(14, '(1000f19.9)')k3len(ik),eigval_soc(i, ik)
         enddo
         write(14, *)' '
      enddo
      close(14)
 
      open(unit=15, file='bulkek-nsoc.dat')
      do i=1, num_wann_soc
         do ik=1, knv3
            write(15, '(1000f19.9)')k3len(ik),eigval_nsoc(i, ik)- X(3)
         enddo
         write(15, *)' '
      enddo
      close(15)
 

      return
   end subroutine fitsoc

   subroutine fitsoc_p

      use para
      implicit none

      integer :: i
      integer :: j
      integer :: l
      integer :: n
      integer :: ik
      real(dp) :: k(3)
      real(dp), allocatable :: x(:)

      integer :: npara

      integer  :: maxfun

      ! output frequency
      integer  :: iprint

      integer  :: iexit

      integer  :: Nitermax
      integer  :: iter

      real(dp) :: ftol

      real(dp), allocatable :: W(:)
      complex(dp), allocatable :: Hamk_soc(:, :)

      real(Dp),allocatable :: y(:)
      real(Dp),allocatable :: p(:, :)

      interface
         function func(n, x)
         implicit none
         integer, parameter :: dp= kind(1d0)
         integer, intent(in) :: n
         real(dp), dimension(n), intent(in) :: x
         real(dp) :: func
         end function func
      end interface

      allocate(W(num_wann_soc))
      allocate(Hamk_soc(num_wann_soc, num_wann_soc))
      W=0d0
      Hamk_soc= 0d0

      knv3= Nk*Nk*Nk
      allocate(kpoints(3, knv3))
      kpoints= 0d0

      ik= 0
      do i=1, Nk
         do j=1, Nk
            do l=1, Nk
               ik= ik+ 1
               kpoints(1, ik)= dble(i)/Nk
               kpoints(2, ik)= dble(j)/Nk
               kpoints(3, ik)= dble(l)/Nk
            enddo !l 
         enddo !j 
      enddo !i 


      !> using kline mode
      deallocate(kpoints)
      allocate(kpoints(3, nk3_band))
      kpoints= k3points

      knv3= nk3_band
      num_bands_DFT= 82
      allocate(weight(num_wann_soc, knv3))
      weight=1d0

      allocate(eigval_soc(num_wann_soc, knv3))
      allocate(eigval_soc_DFT(num_bands_DFT, knv3))
      allocate(eigval_nsoc(num_wann_soc, knv3))
      eigval_soc= 0d0
      eigval_soc_DFT= 0d0
     
      !> get eigenvalues for soc hamiltonian
      do ik=1, knv3
         k= kpoints(:, ik)
         if (abs(k(1))<0.2d0 .and. abs(k(2))<0.2d0 .and. abs(k(3))<0.2d0)then
        !   weight(:, ik)=30d0
         endif
         if (abs(k(1))<0.1d0 .and. abs(k(2))<0.1d0 .and. abs(k(3))<0.1d0)then
        !   weight(:, ik)=300d0
         endif
         if (abs(k(1))<0.02d0 .and. abs(k(2))<0.02d0 .and. abs(k(3))<0.02d0)then
            weight(:, ik)=3000d0
         endif
         W=0d0
         Hamk_soc= 0d0
         call Hamk_bulk_soc(k, Hamk_soc)
         call eigensystem_c('N', 'U', num_wann_soc, Hamk_soc, W)
         eigval_soc(:, ik)= W
      enddo

      !> get eigenvalues for soc hamiltonian from DFT
     !call EigvalFromDFT
     !eigval_soc=eigval_soc_DFT(19:18+num_wann_soc, :)
     !eigval_soc(1:num_bands_DFT-17+1, :)=eigval_soc_DFT(17:num_bands_DFT, :)

      npara= 2*Num_atom_type+ 1
      allocate(x(npara))
      allocate(y(npara+1))
      allocate(p(npara+1, npara))
      x= 0d0
      y=0d0
      p=0d0

      !> initial value for spin-orbital coupling
      x(1:Num_atom_type)= lambda_p(1:Num_atom_type)
      x(1+Num_atom_type:2*Num_atom_type)= lambda_d(1:Num_atom_type)
      x(Num_atom_type*2+ 1)= E_fermi

      Nitermax=10000
      do i=1, npara+1
         p(i, :)= x
      enddo
      do i=1, npara
         p(i+1, i)= p(i+1, i)+0.02d0*i
      enddo

      do i=1, npara+1
         y(i)= func(npara, p(i, :))
      enddo

      ftol=1d-6 
      call amoeba(npara, p, y, ftol, func, iter)

      x= p(1, :)
      write(*,*)'iter', iter
      write(*, '(30f10.4)')x
      write(*, *)' '

      open(unit=14, file='bulkek-soc.dat')
      do i=1, num_wann_soc
     !do i=1, num_bands_DFT-16
         do ik=1, knv3
            write(14, '(1000f19.9)')k3len(ik),eigval_soc(i, ik)
         enddo
         write(14, *)' '
      enddo
      close(14)
 
      open(unit=15, file='bulkek-nsoc.dat')
      do i=1, num_wann_soc
         do ik=1, knv3
            write(15, '(1000f19.9)')k3len(ik),eigval_nsoc(i, ik)- X(npara)
         enddo
         write(15, *)' '
      enddo
      close(15)
 

      return
   end subroutine fitsoc_p


   !> take eigenvalues from DFT
   !> read ek_nxy.dat
   subroutine EigvalFromDFT

      use para
      implicit none

      integer :: ik
      real(dp) :: r1

      real(dp), allocatable :: r2(:, :)

      !> only for WTe2
      allocate(r2(num_bands_DFT,180))
      r2=0d0

      open(unit=140, file='ek_nxy.dat')

      do ik=1, 100
         read(140, *)r1, r2(:, ik)
      enddo

      do ik=1, nk3_band
         eigval_soc_DFT(:, ik)= r2(:, ik)
      enddo

      close(140)

      return
   end subroutine EigvalFromDFT

