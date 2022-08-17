   subroutine FindNodes
      ! We try to find all the local minimal of the energy gap in 3D BZ,
      ! basically, we generate a serials of starting k point for FindNode_k0
      !
      ! By QuanSheng Wu 
      !
      ! wuquansheng@gmail.com
      !
      ! Nov 9 2016  at ETHZ

      use para
      implicit none

      integer :: ik, ikx, iky, ikz, knv3, ierr, Nleft
      real(dp) :: k(3), k_cart(3), k_out(3), gap_out
      real(dp), allocatable :: kabc_minimal(:, :)
      real(dp), allocatable :: gap_minimal(:)
      real(dp), allocatable :: kabc_minimal_mpi(:, :)
      real(dp), allocatable :: gap_minimal_mpi(:)
  
      interface
         function func_energy(x)
         implicit none
         integer, parameter :: dp= kind(1d0)
         real(dp), dimension(3), intent(in) :: x
         real(dp) :: func_energy
         end function func_energy
      end interface

      logical :: In_KCUBE

      knv3= Nk1*Nk2*Nk3
      allocate(kabc_minimal(3, knv3))
      allocate(gap_minimal(knv3))
      allocate(kabc_minimal_mpi(3, knv3))
      allocate(gap_minimal_mpi(knv3))

      kabc_minimal= 0d0
      gap_minimal= 0d0
     

      do ik=1+cpuid, knv3, num_cpu
         ikx= (ik- 1)/(Nk2*Nk3)+ 1
         iky= (ik- (ikx-1)*Nk2*Nk3- 1)/Nk3+ 1
         ikz= ik- (ikx-1)*Nk2*Nk3- (iky-1)*Nk3 
         if (cpuid.eq.0) then
            write(stdout, *)'>>>Find nodes for ik', ik, 'in knv3', knv3
         endif
         
         k= K3D_start_cube+ K3D_vec1_cube*(ikx-1)/dble(nk1)  &
                   + K3D_vec2_cube*(iky-1)/dble(nk2)  &
                   + K3D_vec3_cube*(ikz-1)/dble(nk3)

         call FindNode_k0(k, k_out, gap_out)
         kabc_minimal(:, ik)= k_out
         gap_minimal(ik)= gap_out
 
      enddo 

#if defined (MPI)
      kabc_minimal_mpi= 0
      gap_minimal_mpi= 0
      call mpi_reduce(kabc_minimal, kabc_minimal_mpi, size(kabc_minimal), &
                      mpi_dp, mpi_sum, 0, mpi_cmw, ierr)
      call mpi_reduce(gap_minimal, gap_minimal_mpi, size(gap_minimal), &
                      mpi_dp, mpi_sum, 0, mpi_cmw, ierr)
#else
      kabc_minimal_mpi= kabc_minimal
      gap_minimal_mpi= gap_minimal
#endif


      outfileindex= outfileindex+ 1 
      if (cpuid==0)then
         open(unit=outfileindex, file='Nodes.dat')
         write (outfileindex, '(a)') '# local minimal position and the related energy gap'
         write(outfileindex, '("#", A10, 80A14)') 'kx', 'ky', 'kz', 'gap', 'E', 'k1', 'k2', 'k3'

         if (index(KPorTB, 'KP')==0)then
            ! transform k into [-0.5:0.5]*[-0.5:0,5]*[-0.5:0.5]
            do ik=1, knv3
               call transformto1BZ(kabc_minimal_mpi(:, ik))
            enddo
         endif

         ! eliminate the duplicated k point
         call eliminate_duplicates(knv3, kabc_minimal_mpi, gap_minimal_mpi, Nleft)

         if (index(KPorTB, 'KP')==0)then
            ! move the k points into Wigner-Seitz cell
            do ik=1, Nleft
               call moveinto_wigner_seitzcell(kabc_minimal_mpi(:, ik))
            enddo
         endif

         do ik=1, Nleft
           !if (gap_minimal_mpi(ik)< Gap_threshold.and.In_KCUBE(kabc_minimal_mpi(:, ik))) then
            if (gap_minimal_mpi(ik)< Gap_threshold) then
               call direct_cart_rec(kabc_minimal_mpi(:, ik), k_cart)
               write(outfileindex, '(80f14.8)') k_cart*Angstrom2atomic, gap_minimal_mpi(ik)/eV2Hartree, &
                  func_energy(kabc_minimal_mpi(:, ik))/eV2Hartree, kabc_minimal_mpi(:, ik)
            endif
         enddo
         close(outfileindex)
      endif
       
      !> write script for gnuplot
      outfileindex= outfileindex+ 1
      if (cpuid==0) then
         open(unit=outfileindex, file='Nodes.gnu')
         write(outfileindex, '(a)')"set encoding iso_8859_1"
         write(outfileindex, '(a)')'#set terminal  postscript enhanced color'
         write(outfileindex, '(a)')"#set output 'gap.eps'"
         write(outfileindex, '(3a)')'set terminal  png      truecolor enhanced', &
            ' size 1920, 1680 font ",36"'
         write(outfileindex, '(a)')"set output 'Nodes.png'"
         write(outfileindex, '(a)')'unset ztics'
         write(outfileindex, '(a)')'unset key'
         write(outfileindex, '(a)')'set ticslevel 0'
         write(outfileindex, '(a)')'#set view equal xyz'
         write(outfileindex, '(a)')'set border 4095 front lt black linewidth 2 dashtype solid'
         write(outfileindex, '(a)')'set xlabel "k_1"'
         write(outfileindex, '(a)')'set ylabel "k_2"'
         write(outfileindex, '(a)')'set zlabel "k_3"'
         write(outfileindex, '(a)')'set xtics nomirror scale 0.5'
         write(outfileindex, '(a)')'set ytics nomirror scale 0.5'
         write(outfileindex, '(a)')'set ztics nomirror scale 0.5'
         write(outfileindex, '(a)')'set xtics offset  0.0,-1.0 , 0'
         write(outfileindex, '(a)')'set ytics offset -2.5,   0 , 0'
         write(outfileindex, '(a)')'set size ratio -1'
         write(outfileindex, '(a)')'#set view 60, 140, 1, 1'
         write(outfileindex, '(a, f10.5, a, f10.5, a)')'set xrange [', -0.5, ':', 0.5, ']'
         write(outfileindex, '(a, f10.5, a, f10.5, a)')'set yrange [', -0.5, ':', 0.5, ']'
         write(outfileindex, '(a, f10.5, a, f10.5, a)')'set zrange [', -0.5, ':', 0.5, ']'
         write(outfileindex, '(2a)')"splot 'Nodes.dat' u 6:7:8 w p pt 7 ps 2"
         close(outfileindex)
     
      endif
      
#if defined (MPI)
      call mpi_barrier(mpi_cmw, ierr)
#endif

      deallocate(kabc_minimal)
      deallocate(gap_minimal)
      deallocate(kabc_minimal_mpi)
      deallocate(gap_minimal_mpi)

      return
   end subroutine FindNodes

   !> https://math.stackexchange.com/questions/1472049/check-if-a-point-is-inside-a-rectangular-shaped-area-3d?noredirect=1&lq=1
   !> Check if a k point is inside a rectangular shaped area (3D)?
   function In_KCUBE(k)
      use para
      implicit none

      real(dp), intent(in) :: k(3)
      logical :: In_KCUBE

      real(dp) :: a, b, c, p0(3), p1(3), p2(3), p3(3), u(3), v(3), w(3)
      p0= K3D_start_cube
      u= K3D_vec1_cube
      v= K3D_vec2_cube
      w= K3D_vec3_cube
      p1= p0+ u
      p2= p0+ v
      p3= p0+ w

      !> check whether the dot product u.x is between u.P0 and u.P1
      a=dot_product(u, k)
      b=dot_product(u, p0)
      c=dot_product(u, p1)
      In_KCUBE= (b<a).and.(a<c)

      !> check whether the dot product v.x is between v.P0 and v.P2
      a=dot_product(v, k)
      b=dot_product(v, p0)
      c=dot_product(v, p2)
      In_KCUBE= In_KCUBE.and.(b<a).and.(a<c)

      !> check whether the dot product w.x is between w.P0 and w.P3
      a=dot_product(w, k)
      b=dot_product(w, p0)
      c=dot_product(w, p3)
      In_KCUBE= In_KCUBE.and.(b<a).and.(a<c)

      return
   end function In_KCUBE


   subroutine FindNode_k0(k0, k_out, gap_out)
      ! We find the local minimal of the energy gap in 3D BZ, 
      ! the result depends on the initial point
      !
      ! By QuanSheng Wu 
      !
      ! wuquansheng@gmail.com
      !
      ! Nov 9 2016  at ETHZ

      use para
      implicit none

      real(dp), intent(in)  :: k0(3)
      real(dp), intent(out) :: k_out(3)
      real(dp)              :: gap_out

      integer :: i

      ! parameters for Proteus
      integer :: npara ! number of parameters to be fit
      integer  :: iter
      real(dp) :: ptol
      real(dp),allocatable :: x(:)
      real(dp),allocatable :: gap(:)
      real(dp),allocatable :: k(:, :)

      interface
         function func_gap(n, x)
         implicit none
         integer, parameter :: dp= kind(1d0)
         integer, intent(in) :: n
         real(dp), dimension(n), intent(in) :: x
         real(dp) :: func_gap
         end function func_gap
      end interface


      npara = 3

      allocate(x(npara))
      allocate(gap(npara+1))
      allocate(k(npara+1, npara))
      x=0d0
      gap=0d0
      k=0d0

      x= k0  ! initial k point

      do i=1, npara+1
         k(i, :)= x
      enddo
      do i=1, npara
         k(i+1, i)= k(i+1, i)+0.02d0*i
      enddo

      do i=1, npara+1
         gap(i)= func_gap(npara, k(i, :))
      enddo

      ptol=1d-5 
      call Proteus(npara, k, gap, ptol, func_gap, iter)

      x= k(1, :)
      gap_out= func_gap(npara, x)

      k_out= x
      if (cpuid.eq.0) then
         write(stdout,*)'iter', iter
         write(stdout, '(4A12)') 'k1', 'k2', 'k3', 'gap'
         write(stdout, '(4f12.6)')x, gap_out
         write(stdout, *)' '
      endif

      return
   end subroutine FindNode_k0

   function func_energy(X)
      ! this function calculates the energy gap at a given k point
      !
      ! By QuanSheng Wu 
      !
      ! wuquansheng@gmail.com
      !
      ! Nov 9 2016  at ETHZ

      use para
      implicit none

      real(Dp),intent(in):: X(3) ! k point coordinates

      real(dp) :: func_energy  ! energy at a given k point

      real(dp), allocatable :: W(:)
      complex(dp), allocatable :: Hamk_bulk(:, :)

      allocate(W(num_wann))
      allocate(Hamk_bulk(num_wann, num_wann))

      ! generate bulk Hamiltonian
      if (index(KPorTB, 'KP')/=0)then
         call ham_bulk_kp  (X, Hamk_bulk)
      else
         call ham_bulk_atomicgauge     (X, Hamk_bulk)
      endif

      ! diagonalization by call zheev in lapack
      W= 0d0
      call eigensystem_c( 'N', 'U', Num_wann ,Hamk_bulk, W)

      func_energy= W(Numoccupied)

      !> deal with phonon system
      !> sign(A, B) returns the value of A with the sign of B.
      if (index(Particle,'phonon')/=0) then
         func_energy= sqrt(abs(func_energy))*sign(1d0, func_energy)
      endif

      deallocate(W, Hamk_bulk)

      return
   end function func_energy

   function func_gap(N,X)
      ! this function calculates the energy gap at a given k point
      !
      ! By QuanSheng Wu 
      !
      ! wuquansheng@gmail.com
      !
      ! Nov 9 2016  at ETHZ

      use para
      implicit none

      integer,intent(in) :: N ! no. of parameters, N=3

      real(Dp),intent(in):: X(N) ! k point coordinates

      real(dp) :: func_gap  ! energy gap at a given k point

      real(dp), allocatable :: W(:)
      complex(dp), allocatable :: Hamk_bulk(:, :)

      allocate(W(num_wann))
      allocate(Hamk_bulk(num_wann, num_wann))
      ! generate bulk Hamiltonian
      if (index(KPorTB, 'KP')/=0)then
         call ham_bulk_kp(X, Hamk_bulk)
      else
         call ham_bulk_atomicgauge  (X, Hamk_bulk)
      endif
 
      ! diagonalization by call zheev in lapack
      W= 0d0
      call eigensystem_c( 'N', 'U', Num_wann ,Hamk_bulk, W)

      func_gap= W(Numoccupied+1)- W(Numoccupied)

      deallocate(W, Hamk_bulk)

      return
   end function func_gap


   subroutine transformto1BZ(kin)
      ! Transform the k points to the 1st BZ
      !
      ! By QuanSheng Wu 
      !
      ! wuquansheng@gmail.com
      !
      ! Nov 9 2016  at ETHZ

      use para, only : dp

      integer :: i
      real(dp), intent(inout) :: kin(3)

      do i=1, 3
         do while (.true.)
            if (kin(i)>=-0.5d0 .and. kin(i)<0.5d0) then
               exit
            else if (kin(i)<-0.5d0) then
               kin(i)= kin(i)+ 1d0
            else if (kin(i)>=0.5d0) then
               kin(i)= kin(i)- 1d0
            endif
         enddo
      enddo

      return
   end subroutine transformto1BZ


   subroutine eliminate_duplicates(knv3, kabc_minimal, gap_minimal, Nleft)
      ! Eliminate the duplicated k points
      !
      ! By QuanSheng Wu 
      !
      ! wuquansheng@gmail.com
      !
      ! Nov 9 2016  at ETHZ

      use para, only : dp, eps6
      integer, intent(in) :: knv3
      integer, intent(out) :: Nleft
      real(dp), intent(inout) :: kabc_minimal(3, knv3)
      real(dp), intent(inout) :: gap_minimal(knv3)
      real(dp), allocatable :: kabc_minimal_left(:, :)
      real(dp), allocatable :: gap_minimal_left(:)

      integer :: it, ik, ik1
      logical :: Logical_duplicate

      allocate(kabc_minimal_left(3, knv3))
      allocate(gap_minimal_left(knv3))
      kabc_minimal_left=0
      kabc_minimal_left(:, 1)= kabc_minimal(:, 1)
      gap_minimal_left(1)= gap_minimal(1)

      Nleft = 1
      it= 1
      do ik=2, knv3
         Logical_duplicate= .False.
         do ik1=1, Nleft
            if (sum(abs(kabc_minimal(:, ik)-kabc_minimal_left(:, ik1)))<eps6) then
               Logical_duplicate= .True.
               exit
            endif
         enddo
         if (.not.Logical_duplicate)then
            Nleft= Nleft+ 1
            kabc_minimal_left(:, Nleft)= kabc_minimal(:, ik)
            gap_minimal_left(Nleft)= gap_minimal(ik)
         endif
      enddo

      kabc_minimal(:, 1:Nleft)= kabc_minimal_left(:, 1:Nleft)
      gap_minimal(1:Nleft)= gap_minimal_left(1:Nleft)

      deallocate(gap_minimal_left, kabc_minimal_left)

      return
   end subroutine eliminate_duplicates

   !> shift the k points into the Wigner-Seitz cell centered at Gamma point.
   !> input : k(3) in fractional units
   !> output : k(3) in fractional units
   subroutine moveinto_wigner_seitzcell(k)
      use para, only : dp, Origin_cell
      implicit none

      integer :: ik, l, m, n, ik0
      real(dp) :: k0(3)
      real(dp) :: klen
      real(dp), intent(inout) :: k(3)

      real(dp) :: smallest_vec(3)
      real(dp) :: smallest_vec_len
      real(dp), allocatable :: shiftedk(:, :)

      allocate(shiftedk(3, 125))
      shiftedk= 0d0

      ik0= 0
      do l=-2, 2
         do m=-2, 2
            do n=-2, 2
               ik0= ik0+ 1
               shiftedk(1, ik0)= k(1)+ l
               shiftedk(2, ik0)= k(2)+ m
               shiftedk(3, ik0)= k(3)+ n
            enddo ! n
         enddo ! m
      enddo ! l

      smallest_vec_len= 999d0
      ! find the smallest k vector in 125 ks
      do ik0=1, 125
         k0= shiftedk(1, ik0)* Origin_cell%Kua+ shiftedk(2, ik0)*Origin_cell%Kub+ shiftedk(3, ik0)*Origin_cell%Kuc
         klen= dsqrt(k0(1)*k0(1)+ k0(2)*k0(2)+ k0(3)*k0(3))
         if (klen< smallest_vec_len) then
            smallest_vec_len= klen
            smallest_vec= shiftedk(:, ik0)
         endif
      enddo

      k = smallest_vec

      return
   end subroutine moveinto_wigner_seitzcell

   function if_in_the_WS(kin)
      ! check whether kin is in the Wigner-Seitz cell
      !
      ! By QuanSheng Wu 
      !
      ! wuquansheng@gmail.com
      !
      ! Jun 11 2018  at EPFL
      ! kin(3) must be in unit of reciprocal lattice

      use para, only : dp, eps6
      implicit none

      integer :: i
      logical :: if_in_the_WS
      real(dp), intent(in) :: kin(3)
      real(dp) :: k(3)

      k= kin
      call moveinto_wigner_seitzcell(k)

      !> check whethe k is the same as kin
      if_in_the_WS= .False.
      if (sum(abs(k-kin))<eps6) if_in_the_WS= .True.

      return
   end function if_in_the_WS


