   !> this suboutine is used for wannier center calculation for 3D system
   !> only for one plane

   subroutine  wannier_center2D_plane
      use para
      use wmpi
      implicit none

      integer :: i
      integer :: j
      integer :: l 
      integer :: m 
      integer :: ia
      integer :: nfill
      integer :: imax

      integer :: ik1
      integer :: ik2

      integer :: ierr

      !> k points in kx-ky plane
      real(dp), allocatable :: kpoints(:, :, :)

      !> hamiltonian for each k point
      !> and also the eigenvector of hamiltonian after eigensystem_c
      complex(dp), allocatable :: Hamk(:, :)
      complex(dp), allocatable :: Hamk_dag(:, :)

      !> eigenvector for each kx
      complex(dp), allocatable :: Eigenvector(:, :, :)

      !> Mmnkb=<u_n(k)|u_m(k+b)>
      !> |u_n(k)> is the periodic part of wave function
      complex(dp), allocatable :: Mmnkb(:, :)
      complex(dp), allocatable :: Mmnkb_com(:, :)
      complex(dp), allocatable :: Mmnkb_full(:, :)

      !> 
      complex(dp), allocatable :: Lambda_eig(:)
      complex(dp), allocatable :: Lambda(:, :)
      complex(dp), allocatable :: Lambda0(:, :)

      !> three matrix for SVD 
      !> M= U.Sigma.V^\dag
      !> VT= V^\dag
      complex(dp), allocatable :: U(:, :)
      real   (dp), allocatable :: Sigma(:, :)
      complex(dp), allocatable :: VT(:, :)
   
      !> wannier centers for each ky, bands
      real(dp), allocatable :: WannierCenterKy(:, :)
      real(dp), allocatable :: WannierCenterKy_mpi(:, :)

      !> larges gap of each two wannier centers for a given k point
      !> dim= Nky
      real(dp), allocatable :: largestgap(:)
      real(dp), allocatable :: largestgap_mpi(:)

      !> eigenvalue
      real(dp), allocatable :: eigenvalue(:)

      !> for each orbital, it correspond to an atom
      !> dim= Num_wann
      integer, allocatable :: AtomIndex_orbital(:)

      real(dp) :: Umatrix_t(2, 2)

      !> b.r
      real(dp) :: br

      !> exp(-i*b.R)
      complex(dp) :: ratio

      real(dp) :: k(2)
      real(dp) :: b(2)

      real(dp) :: maxgap
      real(dp) :: maxgap0

      !> Z2 calculation for time reversal invariant system
      integer :: Z2
      integer :: Delta

      real(dp) :: g
      real(dp) :: phi1
      real(dp) :: phi2
      real(dp) :: phi3
      real(dp) :: zm
      real(dp) :: zm1
      real(dp) :: xnm1
      real(dp) :: Deltam
      real(dp), allocatable :: xnm(:)
      real(dp) :: k0(2), k1(2), k2(2)



      nfill= Numoccupied

      allocate(kpoints(2, Nk1, Nk2))
      kpoints= 0d0

      allocate(Lambda_eig(nfill))
      allocate(Lambda(nfill, nfill))
      allocate(Lambda0(nfill, nfill))
      allocate(Mmnkb(nfill, nfill))
      allocate(Mmnkb_com(nfill, nfill))
      allocate(Mmnkb_full(Num_wann, Num_wann))
      allocate(hamk(Num_wann, Num_wann))
      allocate(hamk_dag(Num_wann, Num_wann))
      allocate(Eigenvector(Num_wann, Num_wann, Nk1))
      allocate(eigenvalue(Num_wann))
      allocate(U(nfill, nfill))
      allocate(Sigma(nfill, nfill))
      allocate(VT(nfill, nfill))
      allocate(WannierCenterKy(nfill, Nk2))
      allocate(WannierCenterKy_mpi(nfill, Nk2))
      allocate(AtomIndex_orbital(Num_wann))
      allocate(xnm(nfill))
      allocate(largestgap(Nk2))
      allocate(largestgap_mpi(Nk2))
      largestgap= 0d0
      largestgap_mpi= 0d0
      WannierCenterKy= 0d0
      WannierCenterKy_mpi= 0d0
      hamk=0d0
      eigenvalue=0d0
      Eigenvector=0d0
      Mmnkb_full=0d0
      Mmnkb=0d0
      Mmnkb_com=0d0
      Lambda =0d0
      Lambda0=0d0
      U= 0d0
      Sigma= 0d0
      VT= 0d0

      !> set k plane
      !> the first dimension should be in one primitive cell, [0, 2*pi]
      !> the first dimension is the integration direction
      !> the WCCs are calculated along the second k line
      k0= K2D_start ! 
      k1= K2D_vec1   !  
      k2= K2D_vec2   ! 

      do ik2=1, Nk2
         do ik1=1, Nk1
            kpoints(:, ik1, ik2)= k0+ k1*(ik1-1d0)/dble(Nk1)+ k2*(ik2-1d0)/dble(Nk2-1)
         enddo
      enddo
      b= k1/dble(Nk1)
      b= b(1)*kua+b(2)*kub

      !> set up atom index for each orbitals in the basis
      if (soc>0) then  !> with spin orbital coupling
         l= 0
         do ia=1, Num_atoms  !> spin up
            do j=1, nprojs(ia)
               l= l+ 1
               AtomIndex_orbital(l)= ia
            enddo ! l
         enddo ! ia
         do ia=1, Num_atoms  !> spin down
            do j=1, nprojs(ia)
               l= l+ 1
               AtomIndex_orbital(l)= ia
            enddo ! l
         enddo ! ia
      else  !> without spin orbital coupling
         l= 0
         do ia=1, Num_atoms  !> spin down
            do j=1, nprojs(ia)
               l= l+ 1
               AtomIndex_orbital(l)= ia
            enddo ! l
         enddo ! ia

      endif

      Umatrix_t= transpose(Umatrix)
      call inv_r(2, Umatrix_t)

      !>> Get wannier center for ky=0 plane
      !> for each ky, we can get wanniercenter
      do ik2=1+ cpuid, Nk2, num_cpu
         if (cpuid.eq.0) write(stdout, *)' Wilson loop ',  'ik, Nk', ik2, nk2
         Lambda0=0d0
         do i=1, nfill
            Lambda0(i, i)= 1d0
         enddo

         !> for each k1, we get the eigenvectors
         do ik1=1, Nk1
            !if (cpuid==0) print *, ik1, ik2, Nk1, Nk2
            k= kpoints(:, ik1, ik2)

            call ham_bulk_old(k,hamk)
           !call ham_bulk    (k,hamk)

            !> diagonal hamk
            call eigensystem_c('V', 'U', Num_wann, hamk, eigenvalue)

            Eigenvector(:, :, ik1)= hamk
         enddo

         !> sum over k1 to get wanniercenters
         do ik1=1, Nk1
            !> <u_k|u_k+1>
            Mmnkb= 0d0
            hamk_dag= Eigenvector(:, :, ik1)
            if (ik1==Nk1) then
               hamk= Eigenvector(:, :, 1)
            else
               hamk= Eigenvector(:, :, ik1+ 1)
            endif
            do m=1, Num_wann
               ia= AtomIndex_orbital(m)
               br= b(1)*Atom_position(1, ia)+ &
                   b(2)*Atom_position(2, ia)
               ratio= cos(br)- zi* sin(br)
              !ratio= 1d0
         
               do j=1, nfill
                  do i=1, nfill
                     Mmnkb(i, j)=  Mmnkb(i, j)+ &
                        conjg(hamk_dag(m, i))* hamk(m, j)* ratio
                  enddo ! i
               enddo ! j
            enddo ! m

            !> perform Singluar Value Decomposed of Mmnkb
            call zgesvd_pack(nfill, Mmnkb, U, Sigma, VT)

            !> after the calling of zgesvd_pack, Mmnkb becomes a temporal matrix
            U= conjg(transpose(U))
            VT= conjg(transpose(VT))
            call mat_mul(nfill, VT, U, Mmnkb)

            call mat_mul(nfill, Mmnkb, Lambda0, Lambda)
            Lambda0 = Lambda
         enddo  !< ik1

         !> diagonalize Lambda to get the eigenvalue 
         call zgeev_pack(nfill, Lambda, Lambda_eig)
         do i=1, nfill
            WannierCenterKy(i, ik2)= aimag(log(Lambda_eig(i)))/2d0/pi
            WannierCenterKy(i, ik2)= mod(WannierCenterKy(i, ik2)+10d0, 1d0) 
         enddo

         call sortheap(nfill, WannierCenterKy(:, ik2))

         maxgap0= -99999d0
         imax= nfill
         do i=1, nfill
            if (i/=nfill) then
               maxgap= WannierCenterKy(i+1, ik2)- WannierCenterKy(i, ik2)
            else
               maxgap=1d0+ WannierCenterKy(1, ik2)- WannierCenterKy(nfill, ik2)
            endif

            if (maxgap>maxgap0) then
               maxgap0= maxgap
               imax= i
            endif

         enddo

         if (imax==nfill) then
            largestgap(ik2)= (WannierCenterKy(1, ik2)+ &
               WannierCenterKy(nfill, ik2) +1d0)/2d0
            largestgap(ik2)= mod(largestgap(ik2), 1d0)
         else
            largestgap(ik2)= (WannierCenterKy(imax+1, ik2)+ &
               WannierCenterKy(imax, ik2))/2d0
         endif


      enddo !< ik2

      WannierCenterKy_mpi= 0d0
      largestgap_mpi= 0d0
      call mpi_allreduce(WannierCenterKy, WannierCenterKy_mpi, &
           size(WannierCenterKy), mpi_dp, mpi_sum, mpi_cmw, ierr)
      call mpi_allreduce(largestgap, largestgap_mpi, &
           size(largestgap), mpi_dp, mpi_sum, mpi_cmw, ierr)

      if (cpuid==0) then
         open(unit=110, file='wcc.dat')

         do ik2=1, Nk2
            write(110, '(10000f16.8)') dble(ik2-1)/dble(Nk2-1)/2d0, &
               largestgap_mpi(ik2), dmod(sum(WannierCenterKy_mpi(:, ik2)), 1d0), & 
               WannierCenterKy_mpi(:, ik2)
         enddo
         close(110)
      endif



      !> Z2 calculation Alexey Soluyanov arXiv:1102.5600

      Delta= 0
      !> for each iky, we get a Deltam
      do ik2=1, nk2-1
      
         !> largestgap position
         zm= largestgap_mpi(ik2)
        !if (ik2==nk2) then
        !   zm1= largestgap_mpi(1)
        !   xnm= WannierCenterKy_mpi(1:nfill, 1)
        !else
            zm1= largestgap_mpi(ik2+1)
            xnm= WannierCenterKy_mpi(1:nfill, ik2+1)
        !endif
         Deltam= 1
         do i=1, nfill
            xnm1= xnm(i)
            phi1= 2d0*pi*zm
            phi2= 2d0*pi*zm1
            phi3= 2d0*pi*xnm1
            
            g= sin(phi2-phi1)+ sin(phi3-phi2)+  sin(phi1-phi3) 
            Deltam= Deltam* sign(1d0, g)
         enddo !i 
         if (Deltam<0) then
            Delta= Delta+ 1
         endif
      enddo !ik2

      Z2= mod(Delta, 2)

      if (cpuid==0) write(stdout, *)'Z2 for the plane you choose: ', Z2


      !> generate gnu script for wannier charge center plots
      if (cpuid==0) then
         open(unit=301, file='wcc.gnu')
         write(301, '(a)')"set encoding iso_8859_1"
         write(301, '(a)')'set terminal  postscript enhanced color font ",30"'
         write(301, '(a)')"set output 'wcc.eps'"
         write(301, '(a)')'unset key '
         write(301, '(a)')'set border lw 3 '
         write(301, '(a)')'set xtics offset 0, 0.2'
         write(301, '(a)')'set xtics format "%4.1f" nomirror out '
         write(301, '(a)')'set xlabel "k" '
         write(301, '(a)')'set xlabel offset 0, 0.7 '
         write(301, '(a)')'set ytics 0.5 '
         write(301, '(a)')'set ytics format "%4.1f" nomirror out'
         write(301, '(a)')'set ylabel "WCC"'
         write(301, '(a)')'set ylabel offset 2, 0.0 '
         write(301, '(a)')'set xrange [0: 0.5]'
         write(301, '(a)')'set yrange [0:1]'
         write(301, '(a)')"plot 'wcc.dat' u 1:2 w l lw 2  lc 'blue', \"   
         write(301, '(a, i5, a)')" for [i=4: ", nfill+3, "] 'wcc.dat' u 1:i w p  pt 7  ps 1.1 lc 'red'"
         close(301)
      endif

      return
   end subroutine  wannier_center2D_plane


   subroutine sortheap(n, arr)
      implicit none
      integer, intent(in) :: n
      real(8), intent(inout) :: arr(n)

      !* local variables
      integer :: i

      do i=n/2, 1, -1
         call sift_down(i, n)
      enddo

      do i=n, 2, -1
         call swap(arr(1), arr(i))
         call sift_down(1, i-1)
      enddo
      contains
      subroutine sift_down(l, r)
         integer, intent(in) :: l, r
         integer :: j, jold
         real(8) :: a
         a= arr(l)
         jold= l
         j= l+ l

         do while (j<=r)
            if (j<r) then
               if (arr(j)<arr(j+1))j=j+1
            endif
            if (a>= arr(j))exit
            arr(jold)= arr(j)
            jold= j
            j= j+ j
         enddo
         arr(jold)= a
         return
      end subroutine sift_down
   end subroutine sortheap

   !>> swap two real numbers
   subroutine swap(a, b)
      real(8), intent(inout) :: a
      real(8), intent(inout) :: b
      real(8) :: c
      c=a
      a=b
      b=c
      return
   end subroutine swap


