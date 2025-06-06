
   subroutine  Z2_3D_adaptive
      ! this suboutine is used for wannier center calculation for 3D system
      use para
      use wmpi
      implicit none

      integer :: ik2, i, j 
      integer :: Nk2_adaptive

      real(dp) :: kstart(3)
      real(dp) :: kvec1(3)
      real(dp) :: kvec2(3)

      integer :: Z2
      integer :: Z2_all(6)

      !> wannier centers for each ky, bands
      real(dp), allocatable :: wcc(:, :)
      real(dp), allocatable :: kpath_wcc(:)
      real(dp), allocatable :: wcc_all(:, :, :)

      !> larges gap of each two wannier centers for a given k point
      !> dim= Nky
      real(dp), allocatable :: largestgap(:)
      real(dp), allocatable :: largestgap_all(:,:)

      allocate(wcc(NumberofSelectedOccupiedBands, Nk2_max))
      allocate(wcc_all(NumberofSelectedOccupiedBands, Nk2_max, 6))
      allocate(largestgap(Nk2_max))
      allocate(largestgap_all(Nk2_max,6))
      allocate(kpath_wcc(Nk2_max))
      largestgap= 0d0
      wcc= 0d0
      wcc_all= 0d0

      !> integration over kc (\bar{c}}, wcc along kb, fixed k1=0
      kstart= (/0.0d0, 0.0d0, 0.0d0/)
      kvec1 = (/0.0d0, 0.0d0, 1.0d0/)
      kvec2 = (/0.0d0, 0.5d0, 0.0d0/)

      call  wannier_center3D_plane_adaptive_func(Kstart, Kvec1, Kvec2, &
         largestgap, wcc, Z2, Nk2_adaptive, kpath_wcc)
      Z2_all(1)= Z2

      outfileindex= outfileindex+ 1
      if (cpuid==0) then
         open(unit=outfileindex, file='wanniercenter3D_Z2_1.dat')
         do i=1, NumberofSelectedOccupiedBands
            do ik2=1, Nk2_adaptive
               write(outfileindex, '(10000f16.8)') kpath_wcc(ik2), &
                  (dmod((wcc(i, ik2)), 1d0))
            enddo
            write(outfileindex, *)' '
         enddo
         close(outfileindex)
      endif


      !> integration over kc (\bar{c}}, wcc along kb, fixed k1=ka/2
      kstart= (/0.5d0, 0.0d0, 0.0d0/)
      kvec1 = (/0.0d0, 0.0d0, 1.0d0/)
      kvec2 = (/0.0d0, 0.5d0, 0.0d0/)

      call  wannier_center3D_plane_adaptive_func(Kstart, Kvec1, Kvec2, &
         largestgap, wcc, Z2, Nk2_adaptive, kpath_wcc)
      Z2_all(2)= Z2

      outfileindex= outfileindex+ 1
      if (cpuid==0) then
         open(unit=outfileindex, file='wanniercenter3D_Z2_2.dat')
         do i=1, NumberofSelectedOccupiedBands
            do ik2=1, Nk2_adaptive
               write(outfileindex, '(10000f16.8)') kpath_wcc(ik2), &
                  (dmod((wcc(i, ik2)), 1d0))
            enddo
            write(outfileindex, *)' '
         enddo
         close(outfileindex)
      endif



      !> integration over kc (\bar{c}}, wcc along ka, fixed k2=0
      kstart= (/0.0d0, 0.0d0, 0.0d0/)
      kvec1 = (/0.0d0, 0.0d0, 1.0d0/)
      kvec2 = (/0.5d0, 0.0d0, 0.0d0/)

      call  wannier_center3D_plane_adaptive_func(Kstart, Kvec1, Kvec2, &
         largestgap, wcc, Z2, Nk2_adaptive, kpath_wcc)
      Z2_all(3)= Z2

      outfileindex= outfileindex+ 1
      if (cpuid==0) then
         open(unit=outfileindex, file='wanniercenter3D_Z2_3.dat')
         do i=1, NumberofSelectedOccupiedBands
            do ik2=1, Nk2_adaptive
               write(outfileindex, '(10000f16.8)') kpath_wcc(ik2), &
                  (dmod((wcc(i, ik2)), 1d0))
            enddo
            write(outfileindex, *)' '
         enddo
         close(outfileindex)
      endif

      !> integration over kc (\bar{c}}, wcc along ka, fixed k2=kb/2
      kstart= (/0.0d0, 0.5d0, 0.0d0/)
      kvec1 = (/0.0d0, 0.0d0, 1.0d0/)
      kvec2 = (/0.5d0, 0.0d0, 0.0d0/)

      call  wannier_center3D_plane_adaptive_func(Kstart, Kvec1, Kvec2, &
         largestgap, wcc, Z2, Nk2_adaptive, kpath_wcc)
      Z2_all(4)= Z2

      outfileindex= outfileindex+ 1
      if (cpuid==0) then
         open(unit=outfileindex, file='wanniercenter3D_Z2_4.dat')
         do i=1, NumberofSelectedOccupiedBands
            do ik2=1, Nk2_adaptive
               write(outfileindex, '(10000f16.8)') kpath_wcc(ik2), &
                  (dmod((wcc(i, ik2)), 1d0))
            enddo
            write(outfileindex, *)' '
         enddo
         close(outfileindex)
      endif

      !> integration over ka (\bar{a}}, wcc along kb, fixed k3=0
      kstart= (/0.0d0, 0.0d0, 0.0d0/)
      kvec1 = (/1.0d0, 0.0d0, 0.0d0/)
      kvec2 = (/0.0d0, 0.5d0, 0.0d0/)

      call  wannier_center3D_plane_adaptive_func(Kstart, Kvec1, Kvec2, &
         largestgap, wcc, Z2, Nk2_adaptive, kpath_wcc)
      Z2_all(5)= Z2

      outfileindex= outfileindex+ 1
      if (cpuid==0) then
         open(unit=outfileindex, file='wanniercenter3D_Z2_5.dat')
         do i=1, NumberofSelectedOccupiedBands
            do ik2=1, Nk2_adaptive
               write(outfileindex, '(10000f16.8)') kpath_wcc(ik2), &
                  (dmod((wcc(i, ik2)), 1d0))
            enddo
            write(outfileindex, *)' '
         enddo
         close(outfileindex)
      endif

      !> integration over ka (\bar{a}}, wcc along kb, fixed k3=kc/2
      kstart= (/0.0d0, 0.0d0, 0.5d0/)
      kvec1 = (/1.0d0, 0.0d0, 0.0d0/)
      kvec2 = (/0.0d0, 0.5d0, 0.0d0/)

      call  wannier_center3D_plane_adaptive_func(Kstart, Kvec1, Kvec2, &
         largestgap, wcc, Z2, Nk2_adaptive, kpath_wcc)
      Z2_all(6)= Z2

      outfileindex= outfileindex+ 1
      if (cpuid==0) then
         open(unit=outfileindex, file='wanniercenter3D_Z2_6.dat')
         do i=1, NumberofSelectedOccupiedBands
            do ik2=1, Nk2_adaptive
               write(outfileindex, '(10000f16.8)') kpath_wcc(ik2), &
                  (dmod((wcc(i, ik2)), 1d0))
            enddo
            write(outfileindex, *)' '
         enddo
         close(outfileindex)
      endif

      outfileindex= outfileindex+ 1
      if (cpuid==0) then
         open(unit=outfileindex, file='wanniercenter3D_Z2.gnu')

         write(outfileindex,'(a)')' #gnuplot version>5.4 '
         write(outfileindex,'(a)')' set encoding iso_8859_1'
         write(outfileindex,'(a)')' set terminal pdf enhanced color font ",12" size 5,4'
         write(outfileindex,'(a)')' set output "wanniercenter3D_Z2.pdf"'
         write(outfileindex,'(a)')' set size 0.6,1.0'
         write(outfileindex,'(a)')' set multiplot '
         write(outfileindex,'(a)')' unset key '
         write(outfileindex,'(a)')' set border lw 1 '
         write(outfileindex,'(3a)')' NOXTICS = "set format x ', "''" ,&
            '; unset xtics; unset xlabel"'
         write(outfileindex,'(3a)')' XTICS = "set xtics format ', "'%4.1f'", &
            '; set xtics 0.5 nomirror in offset 0, 0.3; set mxtics 5;"'
         write(outfileindex,'(3a)')' NOYTICS = "set format y '," '';", 'unset ylabel"'
         write(outfileindex,'(3a)')' YTICS = "set ytics format '," '%1.0f' ",&
            '1.0  nomirror in offset 0.7,0; set mytics 2;"'
         write(outfileindex,'(a)')' TMARGIN = "set tmargin at screen 0.96; set bmargin at screen 0.72"'
         write(outfileindex,'(a)')' MMARGIN = "set tmargin at screen 0.63; set bmargin at screen 0.39"'
         write(outfileindex,'(a)')' BMARGIN = "set tmargin at screen 0.30; set bmargin at screen 0.06"'
         write(outfileindex,'(a)')' LMARGIN = "set lmargin at screen 0.20; set rmargin at screen 0.45"'
         write(outfileindex,'(a)')' RMARGIN = "set lmargin at screen 0.50; set rmargin at screen 0.75"'
         write(outfileindex,'(a)')' TITLE = "offset 0, -0.7"'
         write(outfileindex,'(3a)')' LCOLOR = "rgb '," '#696969'",'"'
         write(outfileindex,'(a)')' '
         write(outfileindex,'(3a)')' POS = "at graph -0.23,1.0 font', " ',12' ",'"'
         write(outfileindex,'(3a)')' POS2 = "at graph -0.15,1.0 font'," ',12' ",'"'
         write(outfileindex,'(a)')' '
         write(outfileindex,'(a)')' set xrange [0: 0.5]'
         write(outfileindex,'(a)')' set yrange [0:1]'
         write(outfileindex,'(a)')' @TMARGIN; @LMARGIN'
         write(outfileindex,'(a)')' @XTICS; @YTICS'
         write(outfileindex,'(a)')' #set title "k_1=0.0" @TITLE'
         write(outfileindex,'(a)')' set xlabel "k_2" offset 0,1.6'
         write(outfileindex,'(a)')' set ylabel "c" rotate by  0 offset 1.8,0'
         write(outfileindex,'(a)')' set label 1 "(a)"  @POS front '
         write(outfileindex,'(a)')' plot "wanniercenter3D_Z2_1.dat" u ($1/2):2 w p  pt 7  ps 0.2 lc @LCOLOR'
         write(outfileindex,'(a)')' '
         write(outfileindex,'(a)')' @TMARGIN; @RMARGIN'
         write(outfileindex,'(a)')' @NOYTICS; @XTICS '
         write(outfileindex,'(a)')' #set title "k_1=0.5" @TITLE'
         write(outfileindex,'(a)')' set label 1 "(b)" @POS2 front'
         write(outfileindex,'(a)')' set xlabel "k_2" offset 0,1.6'
         write(outfileindex,'(a)')' unset ylabel'
         write(outfileindex,'(a)')' plot "wanniercenter3D_Z2_2.dat" u ($1/2):2 w p  pt 7  ps 0.2 lc @LCOLOR'
         write(outfileindex,'(a)')' '
         write(outfileindex,'(a)')' @MMARGIN; @LMARGIN'
         write(outfileindex,'(a)')' @YTICS; @XTICS '
         write(outfileindex,'(a)')' #set title "k_2=0.0" @TITLE'
         write(outfileindex,'(a)')' set label 1 "(c)" @POS front'
         write(outfileindex,'(a)')' set xlabel "k_1" offset 0,1.6'
         write(outfileindex,'(a)')' set ylabel "c" rotate by  0 offset 1.8,0'
         write(outfileindex,'(a)')' plot "wanniercenter3D_Z2_3.dat" u ($1/2):2 w p  pt 7  ps 0.2 lc @LCOLOR'
         write(outfileindex,'(a)')' '
         write(outfileindex,'(a)')' @MMARGIN; @RMARGIN'
         write(outfileindex,'(a)')' @NOYTICS; @XTICS '
         write(outfileindex,'(a)')' set label 1 "(d)" @POS2 front'
         write(outfileindex,'(a)')' #set title "k_2=0.5" @TITLE'
         write(outfileindex,'(a)')' set xlabel "k_1" offset 0,1.6'
         write(outfileindex,'(a)')' plot "wanniercenter3D_Z2_4.dat" u ($1/2):2 w p  pt 7  ps 0.2 lc @LCOLOR'
         write(outfileindex,'(a)')''
         write(outfileindex,'(a)')' @BMARGIN; @LMARGIN'
         write(outfileindex,'(a)')' @YTICS; @XTICS '
         write(outfileindex,'(a)')' #set title "k_3=0.0" @TITLE'
         write(outfileindex,'(a)')' set label 1 "(e)" @POS front'
         write(outfileindex,'(a)')' set ylabel "a" rotate by  0 offset 1.8,0'
         write(outfileindex,'(a)')' set xlabel "k_2" offset 0,1.6'
         write(outfileindex,'(a)')' plot "wanniercenter3D_Z2_5.dat" u ($1/2):2 w p  pt 7  ps 0.2 lc @LCOLOR'
         write(outfileindex,'(a)')' '
         write(outfileindex,'(a)')' '
         write(outfileindex,'(a)')' @BMARGIN; @RMARGIN'
         write(outfileindex,'(a)')' @NOYTICS; @XTICS '
         write(outfileindex,'(a)')' #set title "k_3=0.5" @TITLE'
         write(outfileindex,'(a)')' set label 1 "(f)" @POS2 front'
         write(outfileindex,'(a)')' set xlabel "k_2" offset 0,1.6'
         write(outfileindex,'(a)')' plot "wanniercenter3D_Z2_6.dat" u ($1/2):2 w p  pt 7  ps 0.2 lc @LCOLOR'
 
         close(outfileindex)
      endif


      if (cpuid==0) then
         write(stdout, *)'# Notes: Please check the Wilson loop plots gnuplot wanniercenter3D_Z2.gnu '
         write(stdout, *)'# z2 number for 6 planes'
         write(stdout, *)'k1=0.0, k2-k3 plane: ', Z2_all(1)
         write(stdout, *)'k1=0.5, k2-k3 plane: ', Z2_all(2)
         write(stdout, *)'k2=0.0, k1-k3 plane: ', Z2_all(3)
         write(stdout, *)'k2=0.5, k1-k3 plane: ', Z2_all(4)
         write(stdout, *)'k3=0.0, k1-k2 plane: ', Z2_all(5)
         write(stdout, *)'k3=0.5, k1-k2 plane: ', Z2_all(6)
      endif

      return
   end subroutine  Z2_3D_adaptive


   subroutine  wannier_center3D_plane_adaptive
      !> this suboutine is used for wannier center calculation for 3D system
      !> only for one plane

      use para
      use wmpi
      implicit none
      integer :: ik2

      real(dp), allocatable :: wcc(:, :)

      !> larges gap of each two wannier centers for a given k point
      !> dim= Nky
      real(dp), allocatable :: largestgap(:)
      real(dp), allocatable :: kpath_wcc(:)

      !> Z2 calculation for time reversal invariant system
      integer :: Z2

      integer :: Nk2_adaptive

      allocate(wcc(NumberofSelectedOccupiedBands, Nk2_max))
      allocate(largestgap(Nk2_max))
      allocate(kpath_wcc(Nk2_max))
      kpath_wcc= 0d0
      largestgap= 0d0
      wcc= 0d0

      call  wannier_center3D_plane_adaptive_func(K3D_start, K3D_vec1, K3D_vec2, &
         largestgap, wcc, Z2, Nk2_adaptive, kpath_wcc)


      outfileindex= outfileindex+ 1
      if (cpuid==0) then
         open(unit=outfileindex, file='wcc.dat')
         write(outfileindex, '(4A16, A)')'#      k', 'largestgap', 'sum(wcc(:,ik))', &
                                          'wcc(i, ik)', '(i=1, NumberofSelectedOccupiedBands)'
         do ik2=1, Nk2_adaptive
            write(outfileindex, '(10000f16.8)')kpath_wcc(ik2), &
               largestgap(ik2), dmod(sum(wcc(:, ik2)), 1d0), & 
               wcc(:, ik2)
         enddo
         close(outfileindex)
      endif

      if (cpuid==0) write(stdout, *)'Z2 for the plane you choose: ', Z2
   
      !> generate gnu script for wannier charge center plots
      outfileindex= outfileindex+ 1
      if (cpuid==0) then
         open(unit=outfileindex, file='wcc.gnu')
         write(outfileindex, '(a)')"# gnuplot version > 5.4"
         write(outfileindex, '(a)')"set encoding iso_8859_1"
         write(outfileindex, '(a)')'set terminal pdf enhanced color font ",24" size 5, 4'
         write(outfileindex, '(a)')"set output 'wcc.pdf'"
         write(outfileindex, '(a)')'unset key '
         write(outfileindex, '(a)')'set border lw 3 '
         write(outfileindex, '(a)')'set xtics offset 0, 0.2'
         write(outfileindex, '(a)')'set xtics format "%4.1f" nomirror out '
         write(outfileindex, '(a)')'set xlabel "k in unit of line 3 of KPLANE BULK" font ",18"'
         write(outfileindex, '(a)')'set xlabel offset 0, 0.7 '
         write(outfileindex, '(a)')'set ytics 0.5 '
         write(outfileindex, '(a)')'set ytics format "%4.1f" nomirror out'
         write(outfileindex, '(a)')'set ylabel "WCC"'
         write(outfileindex, '(a)')'set ylabel offset 2, 0.0 '
         write(outfileindex, '(a)')'set xrange [0: 1.0]'
         write(outfileindex, '(a)')'set yrange [0:1]'
         write(outfileindex, '(a, i5, a)')"plot for [i=4: ", NumberofSelectedOccupiedBands+3, &
            "] 'wcc.dat' u 1:i w p  pt 7  ps 0.5 lc 'black'"
         close(outfileindex)
      endif

      return
   end subroutine  wannier_center3D_plane_adaptive


   subroutine  wannier_center3D_plane_adaptive_func(kstart, kvec1, kvec2, &
         largestgap, wcc, Z2, Nk2_adaptive, kpath_wcc)
      !> this suboutine is used for wannier center calculation for 3D system
      !> only for one plane

      use para
      use wcc_module
      use wmpi
      implicit none

      integer :: i, j, l , m , ia, imax
      integer :: ik, ik1, ik2
      integer :: ierr

      real(dp), intent(in) :: kstart(3)
      real(dp), intent(in) :: kvec1(3)  ! the integration direction
      real(dp), intent(in) :: kvec2(3)
      real(dp), intent(out) :: largestgap(Nk2_max)
      real(dp), intent(out) :: wcc(NumberofSelectedOccupiedBands, Nk2_max)
      real(dp), intent(out) :: kpath_wcc(Nk2_max)
      integer , intent(out) :: Nk2_adaptive  ! number of k points in end

      !> k points in kx-ky plane
      real(dp), allocatable :: kpoints(:, :)
  
      !> wannier centers for each ky, bands
      real(dp), allocatable :: wcc_one_k(:)
      real(dp), allocatable :: wcc_gap(:)
      real(dp), allocatable :: wcc_current(:)
      real(dp), allocatable :: wcc_next(:)

      !> some arrays for mpi
      real(dp), allocatable :: wcc_k(:, :)
      real(dp), allocatable :: wcc_k_mpi(:, :)
      real(dp), allocatable :: largestgap_val_arr(:)
      real(dp), allocatable :: largestgap_val_arr_mpi(:)
      real(dp), allocatable :: largestgap_pos_val_arr(:)
      real(dp), allocatable :: largestgap_pos_val_arr_mpi(:)
      real(dp), allocatable :: largestgap_pos_i_arr(:)
      real(dp), allocatable :: largestgap_pos_i_arr_mpi(:)

      !> maximal tolerance of the wcc between two neighbours
      real(dp) :: neighbour_tol, wcc_tol
      real(dp) :: k2(3)
      real(dp) :: largestgap_val, largestgap_pos_i, largestgap_pos_val

      !> define kline
      type(kline_wcc_type) :: kline_wcc(Nk2_max)

      !> Z2 calculation for time reversal invariant system
      integer :: Z2
      integer :: Delta

      integer :: iter

      logical :: converged, exceed

      real(dp) :: g, phi1, phi2, phi3, zm
      real(dp) :: zm1, xnm1, Deltam, dis
      real(dp), allocatable :: xnm(:)

      allocate(xnm(NumberofSelectedOccupiedBands))
      allocate(wcc_one_k(NumberofSelectedOccupiedBands))
      allocate(wcc_gap(NumberofSelectedOccupiedBands))
      allocate(wcc_current(NumberofSelectedOccupiedBands))
      allocate(wcc_next(NumberofSelectedOccupiedBands))
      allocate(wcc_k(NumberofSelectedOccupiedBands, Nk2_max))
      allocate(wcc_k_mpi(NumberofSelectedOccupiedBands, Nk2_max))
      allocate(largestgap_val_arr(Nk2_max))
      allocate(largestgap_val_arr_mpi(Nk2_max))
      allocate(largestgap_pos_i_arr(Nk2_max))
      allocate(largestgap_pos_i_arr_mpi(Nk2_max))
      allocate(largestgap_pos_val_arr(Nk2_max))
      allocate(largestgap_pos_val_arr_mpi(Nk2_max))
     
    
      !> first we set an uniform mesh along kvec2
      allocate(kpoints(3, Nk2))
      kpoints= 0d0

      !> initial k points
      do ik2=1, Nk2
         kpoints(:, ik2)= kstart+ kvec2*(ik2-1d0)/dble(Nk2-1)
         kpath_wcc(ik2)= (ik2-1d0)/(Nk2-1)
      enddo

      !> fill data for kline_wcc
      do ik2=1, Nk2
         kline_wcc(ik2)%k= kpoints(:, ik2)
         kline_wcc(ik2)%delta= kpath_wcc(ik2)
         kline_wcc(ik2)%largestgap_pos_i= 0
         kline_wcc(ik2)%largestgap_pos_val= 0d0
         kline_wcc(ik2)%largestgap_val= 0d0
         kline_wcc(ik2)%converged= .False.
         kline_wcc(ik2)%calculated= .False.
         allocate(kline_wcc(ik2)%wcc(NumberofSelectedOccupiedBands))
         allocate(kline_wcc(ik2)%gap(NumberofSelectedOccupiedBands))
         kline_wcc(ik2)%wcc= 0d0
         kline_wcc(ik2)%gap= 0d0
      enddo
      

      !>> Get wannier center for ky=0 plane, for each ky, we can get wanniercenter
      !> adaptively increase the k points until it converges
      Nk2_adaptive= Nk2
      
      iter= 0
      neighbour_tol= wcc_neighbour_tol
      converged=.False.
      Do while (.true.)
         if (converged .or. Nk2_adaptive>=Nk2_max) exit
         iter= iter+ 1
         largestgap_val_arr=0d0
         largestgap_val_arr_mpi=0d0
         largestgap_pos_i_arr=0d0
         largestgap_pos_i_arr_mpi=0d0
         largestgap_pos_val_arr=0d0
         largestgap_pos_val_arr_mpi=0d0
         wcc_k= 0d0
         wcc_k_mpi= 0d0

         do ik2=1+ cpuid, Nk2_adaptive, num_cpu
            k2= kline_wcc(ik2)%k
            if (kline_wcc(ik2)%calculated) then
               largestgap_val_arr(ik2)= kline_wcc(ik2)%largestgap_val
               largestgap_pos_i_arr(ik2)= kline_wcc(ik2)%largestgap_pos_i
               largestgap_pos_val_arr(ik2)= kline_wcc(ik2)%largestgap_pos_val
               wcc_k(:, ik2)= kline_wcc(ik2)%wcc
               cycle ! < skip the calculated point
            endif
            call Wcc_integrate_func(k2, kvec1, wcc_one_k, &
               wcc_gap, largestgap_val, largestgap_pos_i, largestgap_pos_val)
            largestgap_val_arr(ik2)= largestgap_val
            largestgap_pos_i_arr(ik2)= largestgap_pos_i
            largestgap_pos_val_arr(ik2)= largestgap_pos_val
            wcc_k(:, ik2)= wcc_one_k
         enddo !< ik2

#if defined (MPI)
         call mpi_allreduce(wcc_k, wcc_k_mpi, &
              size(wcc_k), mpi_dp, mpi_sum, mpi_cmw, ierr)
         call mpi_allreduce(largestgap_val_arr, largestgap_val_arr_mpi, &
              size(largestgap_val_arr), mpi_dp, mpi_sum, mpi_cmw, ierr)
         call mpi_allreduce(largestgap_pos_i_arr, largestgap_pos_i_arr_mpi, &
              size(largestgap_pos_i_arr), mpi_dp, mpi_sum, mpi_cmw, ierr)
         call mpi_allreduce(largestgap_pos_val_arr, largestgap_pos_val_arr_mpi, &
              size(largestgap_pos_val_arr), mpi_dp, mpi_sum, mpi_cmw, ierr)
#else
         wcc_k_mpi= wcc_k
         largestgap_val_arr_mpi= largestgap_val_arr
         largestgap_pos_i_arr_mpi= largestgap_pos_i_arr
         largestgap_pos_val_arr_mpi= largestgap_pos_val_arr
#endif

         do ik2=1, Nk2_adaptive
            kline_wcc(ik2)%wcc= wcc_k_mpi(:, ik2)
            kline_wcc(ik2)%largestgap_val= largestgap_val_arr_mpi(ik2)
            kline_wcc(ik2)%largestgap_pos_i= largestgap_pos_i_arr_mpi(ik2)
            kline_wcc(ik2)%largestgap_pos_val= largestgap_pos_val_arr_mpi(ik2)
            kline_wcc(ik2)%calculated= .TRUE.
         enddo !< ik2

         kline_wcc(Nk2_adaptive)%converged= .TRUE. !< the last point is converged 

         !> check the difference between the neighbours
         do ik2=1, Nk2_adaptive-1
            largestgap_val= max(kline_wcc(ik2)%largestgap_pos_val, &
               kline_wcc(ik2+1)%largestgap_pos_val)  
            wcc_current= dmod(kline_wcc(ik2)%wcc+10d0-largestgap_val, 1d0)
            wcc_next= dmod(kline_wcc(ik2+1)%wcc+10d0-largestgap_val, 1d0)
            call sortheap(NumberofSelectedOccupiedBands, wcc_current)
            call sortheap(NumberofSelectedOccupiedBands, wcc_next)

            !> get the largest gap
            wcc_tol=-1d0
            do i=1, NumberofSelectedOccupiedBands
               call dis_phase(wcc_next(i), wcc_current(i), dis)
               if ( dis>wcc_tol) wcc_tol= dis
            enddo ! i

            largestgap_val= min(kline_wcc(ik2)%largestgap_val, &
               kline_wcc(ik2+1)%largestgap_val)  
            if (wcc_tol< neighbour_tol*largestgap_val) then
               kline_wcc(ik2)%converged= .TRUE.
            else
               kline_wcc(ik2)%converged= .FALSE.
            endif
         enddo !> ik2

         !> add more k points to get a converged k line

         !> first check if the added k points will exceed the maximal No. k points
         ik= Nk2_adaptive
         exceed= .false.
         do ik2=1, Nk2_adaptive-1
            if (.not.kline_wcc(ik2)%converged)then ! < add one k point between ik2 and ik2+1
               ik= ik+ 1 !> add to the bottom of kline_wcc
               if (ik> Nk2_max) exceed= .true.
            endif
         enddo

         if (exceed) then
            write(stdout, *)'Wcc calculation is not converged, may be there are some nodal points'
            exit
         endif

         !> if not, then add k points
         ik= Nk2_adaptive
         do ik2=1, Nk2_adaptive-1
            if (.not.kline_wcc(ik2)%converged)then ! < add one k point between ik2 and ik2+1
               ik= ik+ 1 !> add to the bottom of kline_wcc
               kline_wcc(ik)%k= (kline_wcc(ik2)%k+ kline_wcc(ik2+1)%k)/2d0 
               kline_wcc(ik)%delta= (kline_wcc(ik2)%delta+ kline_wcc(ik2+1)%delta)/2d0 
               if (.not.allocated(kline_wcc(ik)%wcc))allocate(kline_wcc(ik)%wcc(NumberofSelectedOccupiedBands))
               if (.not.allocated(kline_wcc(ik)%gap))allocate(kline_wcc(ik)%gap(NumberofSelectedOccupiedBands))
               kline_wcc(ik)%calculated= .False.
               kline_wcc(ik)%converged= .False.
               kline_wcc(ik)%largestgap_val= 0d0
               kline_wcc(ik)%largestgap_pos_val= 0d0
               kline_wcc(ik)%largestgap_pos_i= 0
            endif
         enddo ! ik2

         Nk2_adaptive= ik

         converged= .true.
         do ik2=1, Nk2_adaptive
            if (.not.kline_wcc(ik2)%converged)converged=.false.
         enddo

         !> sorted the kline_wcc according to kline_wcc(ik)%delta
         call kline_wcc_sorted(kline_wcc, Nk2_adaptive)
         if (cpuid.eq.0) then
            write(stdout, *)'Wcc iter:', iter
            do ik2=1, Nk2_adaptive
               write(stdout, '(i5,3f10.5,2L2,2f10.5)')ik2, kline_wcc(ik2)%k, &
                  kline_wcc(ik2)%calculated, kline_wcc(ik2)%converged,&
                  kline_wcc(ik2)%largestgap_val*neighbour_tol
            enddo
         endif

      Enddo !< converge criterion

      do ik2=1, Nk2_adaptive
         wcc(:, ik2)= kline_wcc(ik2)%wcc
         kpath_wcc(ik2)= kline_wcc(ik2)%delta
         largestgap(ik2)= kline_wcc(ik2)%largestgap_pos_val
      enddo

      !> Z2 calculation Alexey Soluyanov arXiv:1102.5600

      Delta= 0
      !> for each iky, we get a Deltam
      do ik2=1, Nk2_adaptive-1
      
         !> largestgap position
         zm= kline_wcc(ik2)%largestgap_pos_val
         zm1= kline_wcc(ik2+1)%largestgap_pos_val         
         xnm= kline_wcc(ik2+1)%wcc
         Deltam= 1
         do i=1, NumberofSelectedOccupiedBands
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

      return
   end subroutine  wannier_center3D_plane_adaptive_func


   subroutine  Wcc_integrate_func(k2, kvec1, wcc, &
         wcc_gap, largestgap_val, largestgap_pos_i, largestgap_pos_val)
      !> this suboutine is used for wannier center calculation for 3D system
      !> only for one plane

      use sparse
      use para
      use wcc_module
      use wmpi
      implicit none

      real(dp), intent(in) :: kvec1(3)  ! the integration direction
      real(dp), intent(in) :: k2(3)
      real(dp), intent(out) :: wcc(NumberofSelectedOccupiedBands)
      real(dp), intent(out) :: wcc_gap(NumberofSelectedOccupiedBands)
      real(dp), intent(out) :: largestgap_val
      real(dp), intent(out) :: largestgap_pos_i
      real(dp), intent(out) :: largestgap_pos_val

      logical :: not_in
      integer :: i, j, l , m , ia, it, imax

      integer :: ik1, ik, Nk_start, Nk_Max, Nk_adaptive

      integer :: ierr, iter

      real(dp) :: dis,  max_diff,  wcc_tol

      !> k points in kx-ky plane
      real(dp), allocatable :: kpoints(:, :)

      !> hamiltonian for each k point
      !> and also the eigenvector of hamiltonian after eigensystem_c
      complex(dp), allocatable :: Hamk(:, :), Hamk_dag(:, :),Sk(:, :)

      !> for sparse hamiltonian
      !> dim= Num_wann*Num_wann
      integer :: nnzmax, nnz
      integer, allocatable :: jcoo(:), icoo(:)
      complex(dp), allocatable :: acoo(:)
      complex(dp), allocatable :: zeigv(:, :)

      !number of ARPACK eigenvalues
      integer :: neval
   
      ! number of Arnoldi vectors
      integer :: nvecs
   
      !> calculate eigenvector or not
      logical :: ritzvec
   
      !shift-invert shiftsigma
      complex(dp) :: shiftsigma=(0d0,0d0)

      !> Mmnkb=<u_n(k)|u_m(k+b)>
      !> |u_n(k)> is the periodic part of wave function
      complex(dp), allocatable :: Mmnkb(:, :), Mmnkb_com(:, :)

      !> 
      complex(dp), allocatable :: Lambda_eig(:), Lambda(:, :), Lambda0(:, :)

      !> three matrix for SVD 
      !> M= U.Sigma.V^\dag
      !> VT= V^\dag
      complex(dp), allocatable :: U(:, :)
      real   (dp), allocatable :: Sigma(:, :)
      complex(dp), allocatable :: VT(:, :)
   
      !> wannier centers for each ky, bands
      real(dp), allocatable :: wcc_old(:)
      real(dp), allocatable :: wcc_new(:)

      !> eigenvalue
      real(dp), allocatable :: eigenvalue(:)

      !> b.r
      real(dp) :: br

      !> exp(-i*b.R)
      complex(dp) :: ratio

      real(dp) :: k(3), b(3)

      real(dp) :: maxgap, maxgap0, maxgap_new, maxgap_old

      type(kline_integrate_type) :: kline_integrate(Nk2_max)

      wcc_tol= wcc_calc_tol

      allocate(Lambda_eig(NumberofSelectedOccupiedBands))
      allocate(Lambda(NumberofSelectedOccupiedBands, NumberofSelectedOccupiedBands))
      allocate(Lambda0(NumberofSelectedOccupiedBands, NumberofSelectedOccupiedBands))
      allocate(Mmnkb(NumberofSelectedOccupiedBands, NumberofSelectedOccupiedBands))
      allocate(Mmnkb_com(NumberofSelectedOccupiedBands, NumberofSelectedOccupiedBands))
      allocate(eigenvalue(Num_wann))
      allocate(U(NumberofSelectedOccupiedBands, NumberofSelectedOccupiedBands))
      allocate(Sigma(NumberofSelectedOccupiedBands, NumberofSelectedOccupiedBands))
      allocate(VT(NumberofSelectedOccupiedBands, NumberofSelectedOccupiedBands))
      allocate(wcc_new(NumberofSelectedOccupiedBands))
      allocate(wcc_old(NumberofSelectedOccupiedBands))
      if (Is_Sparse) then
         allocate(hamk(Num_wann, NumberofSelectedOccupiedBands))
         allocate(hamk_dag(Num_wann, NumberofSelectedOccupiedBands))
      else
         allocate(hamk(Num_wann, Num_wann))
         allocate(hamk_dag(Num_wann, Num_wann))
      endif
      if (.not. Orthogonal_Basis) then
         if (Is_Sparse) then
            stop 'wannier center calculation for northogonal basis is not implemented for sparse hamiltonian now'
         else
            !> Sk is the overlap matrix
            !> Sk= <u_n(k)|u_m(k+b)>
            !> |u_n(k)> is the periodic part of wave function
            allocate(Sk(Num_wann, Num_wann))
            Sk= 0d0
         endif
      endif
         

      wcc_old= 0d0
      hamk=0d0
      eigenvalue=0d0
      Mmnkb=0d0
      Mmnkb_com=0d0
      Lambda =0d0
      Lambda0=0d0
      U= 0d0
      Sigma= 0d0
      VT= 0d0

      if (Is_Sparse) then

         if (OmegaNum==0) OmegaNum= Num_wann
         if (NumSelectedEigenVals==0) NumSelectedEigenVals=OmegaNum
      
         neval= NumSelectedEigenVals
      
         if (neval>Num_wann-2) then
            neval= Num_wann- 2
            nvecs= Num_wann
         endif
      
         !> ncv= NumLCZVecs if specfied in wt.in
         if (NumLCZVecs.ne.0) then
            nvecs= NumLCZVecs
         else
            nvecs=int(2*neval)
            if (nvecs<50) nvecs= int(6*neval)
         endif
      
         if (neval+2>=nvecs) neval= nvecs-2
      
         if (nvecs>Num_wann) nvecs= Num_wann

         shiftsigma=(1d0,0d0)*iso_energy
         nnzmax=splen+Num_wann
         nnz=splen
         allocate( acoo(nnzmax))
         allocate( jcoo(nnzmax))
         allocate( icoo(nnzmax))
         allocate( zeigv(Num_wann,nvecs))
      endif

      !> the first dimension should be in one primitive cell, [0, 1] 
      !> the first dimension is the integration direction
      !> the WCCs are calculated along the second k line
      !> For chern number calculation, the second k line should be in one primitive
      !> vector. 
      !> For  Z2 number clculation, the second k line should 
      !> from one TRIM to another TRIM, usually, we study half of the
      !> reciprocal lattice vector

      !>> at the begining, we set the kmesh to be 4
      !>> initial the integration k line
      Nk_start= 32
      allocate(kpoints(3, Nk_start))
      do ik1=1, Nk_start
         kpoints(:, ik1)= k2+ kvec1*(ik1-1d0)/dble(Nk_start)
      enddo

      do ik1=1, Nk_start
         kline_integrate(ik1)%k= kpoints(:, ik1)
         kline_integrate(ik1)%delta= (ik1-1d0)/dble(Nk_start)
         if (ik1==Nk_start) then
            b= kpoints(:, 2)- kpoints(:, 1)
            b= b(1)*Origin_cell%Kua+b(2)*Origin_cell%Kub+b(3)*Origin_cell%Kuc
            kline_integrate(ik1)%b = b
         else
            b= kpoints(:, ik1+1)- kpoints(:, ik1)
            b= b(1)*Origin_cell%Kua+b(2)*Origin_cell%Kub+b(3)*Origin_cell%Kuc
            kline_integrate(ik1)%b= b
         endif
         if (.not.allocated(kline_integrate(ik1)%eig_vec)) &
            allocate(kline_integrate(ik1)%eig_vec(Num_wann, NumberofSelectedOccupiedBands))
         kline_integrate(ik1)%eig_vec= (0d0, 0d0)
         kline_integrate(ik1)%calculated= .False. ! not calculated
      enddo

      !>> increase the k mesh until the wcc converges
      Nk_max= 1025
      Nk_adaptive= Nk_start
      max_diff= 999
      iter= 0
      Do while (Nk_adaptive<Nk_max.and. max_diff> wcc_tol)
         iter= iter+ 1
         wcc_old = wcc

         Nk_adaptive= Nk_adaptive*2  ! double the k mesh at each iteration

         deallocate(kpoints)
         allocate(kpoints(3, Nk_adaptive))
         do ik1=1, Nk_adaptive
            kpoints(:, ik1)= k2+ kvec1*(ik1-1d0)/dble(Nk_adaptive)
         enddo


         !> check whether this k point is in the kline_integrate or not.
         it = Nk_adaptive/2
         do ik1=1, Nk_adaptive
            not_in= .TRUE.
            do ik=1, Nk_adaptive/2
               if ((abs((ik1-1)/dble(Nk_adaptive)-kline_integrate(ik)%delta))<eps6) not_in= .False.
            enddo
            !> if this k point is not in the quenue, add this k point to the end of kline_integrate
            if (not_in) then
               it= it+ 1
               kline_integrate(it)%k= kpoints(:, ik1)
               kline_integrate(it)%delta= (ik1-1d0)/Nk_adaptive
               if (.not. allocated(kline_integrate(it)%eig_vec))&
                   allocate(kline_integrate(it)%eig_vec(Num_wann, NumberofSelectedOccupiedBands))
               kline_integrate(it)%eig_vec= 0d0
               kline_integrate(it)%calculated= .False.
            endif
         enddo
         if (it/=Nk_adaptive) stop 'Error: Nk_adaptive should equal to it'
 
         !>> Secondly, sort the kline
         call kline_integrate_sorted(kline_integrate, Nk_adaptive)
 
         !>> Firstly, for each k1, we get the eigenvectors

         do ik1=1, Nk_adaptive
            if (kline_integrate(ik1)%calculated) then
               cycle
            endif
            k= kline_integrate(ik1)%k
            if (Is_Sparse) then
               call ham_bulk_coo_sparsehr_latticegauge(k,acoo,icoo,jcoo)
               nnz= splen
              
               ritzvec= .true.
               call arpack_sparse_coo_eigs(Num_wann,nnzmax,nnz,acoo,jcoo,icoo,neval,nvecs,eigenvalue,shiftsigma, zeigv, ritzvec)
               kline_integrate(ik1)%eig_vec(1:Num_wann, 1:NumberofSelectedOccupiedBands)= &
                  zeigv(1:Num_wann, Selected_Occupiedband_index(1):Selected_Occupiedband_index(NumberofSelectedOccupiedBands))

            else
               !> get the TB hamiltonian in k space
               ! call ham_bulk_latticegauge(k,hamk)
               if (Orthogonal_Basis) then
                  call ham_bulk_latticegauge(k, hamk)
               else
                  
                  Sk= 0d0
                  call S_bulk_latticegauge(k, Sk)
                  call ham_bulk_latticegauge(k, hamk)
                  call orthogonalize_hamiltonian(hamk, Sk, Num_wann)
                  do i = 1, Num_wann
                     hamk(i, i)= hamk(i, i) - E_fermi
                  enddo !
               endif
               !> diagonal hamk
               call eigensystem_c('V', 'U', Num_wann, hamk, eigenvalue)
               kline_integrate(ik1)%eig_vec(1:Num_wann, 1:NumberofSelectedOccupiedBands)= &
                  hamk(1:Num_wann, Selected_Occupiedband_index(1):Selected_Occupiedband_index(NumberofSelectedOccupiedBands))
 
            endif
   
         enddo  ! ik1

         do ik1=1, Nk_adaptive
            kline_integrate(ik1)%calculated= .TRUE.
         enddo

         Lambda0=0d0
         do i=1, NumberofSelectedOccupiedBands
            Lambda0(i, i)= 1d0
         enddo
  
         !>> Thirdly, sum over k1 to get wanniercenters
         do ik1=1, Nk_adaptive
            !> <u_k|u_k+1>
            Mmnkb= 0d0
            hamk_dag(:, 1:NumberofSelectedOccupiedBands)= kline_integrate(ik1)%eig_vec
            if (ik1==Nk_adaptive) then
               hamk(:, 1:NumberofSelectedOccupiedBands)= kline_integrate(1)%eig_vec     
            else
               hamk(:, 1:NumberofSelectedOccupiedBands)= kline_integrate(ik1+1)%eig_vec 
            endif

            b= kline_integrate(ik1)%b
            do m=1, Num_wann
               br= b(1)*Origin_cell%wannier_centers_cart(1, m)+ &
                   b(2)*Origin_cell%wannier_centers_cart(2, m)+ &
                   b(3)*Origin_cell%wannier_centers_cart(3, m)
               ratio= cos(br)- zi* sin(br)
         
               do j=1, NumberofSelectedOccupiedBands
                  do i=1, NumberofSelectedOccupiedBands
                     Mmnkb(i, j)=  Mmnkb(i, j)+ &
                        conjg(hamk_dag(m, i))* hamk(m, j)* ratio
                  enddo ! i
               enddo ! j
            enddo ! m
   
            !> perform Singluar Value Decomposed of Mmnkb
            call zgesvd_pack(NumberofSelectedOccupiedBands, Mmnkb, U, Sigma, VT)
   
            !> after the calling of zgesvd_pack, Mmnkb becomes a temporal matrix
            U = conjg(transpose(U))
            VT= conjg(transpose(VT))
            call mat_mul(NumberofSelectedOccupiedBands, VT, U, Mmnkb)
   
            call mat_mul(NumberofSelectedOccupiedBands, Mmnkb, Lambda0, Lambda)
            Lambda0 = Lambda
         enddo  !< ik1
   
         !> diagonalize Lambda to get the eigenvalue 
         call zgeev_pack(NumberofSelectedOccupiedBands, Lambda, Lambda_eig)
         do i=1, NumberofSelectedOccupiedBands
            wcc(i)= aimag(log(Lambda_eig(i)))/2d0/pi
            wcc(i)= dmod(wcc(i)+10d0, 1d0) 
         enddo
  
         call sortheap(NumberofSelectedOccupiedBands, wcc(:))

         wcc_new= wcc
         call largest_gap_find(wcc_new, maxgap_new)
         call largest_gap_find(wcc_old, maxgap_old)

         maxgap_new= max(maxgap_new, maxgap_old)
         wcc_new= dmod(wcc_new+ 10d0- maxgap_new, 1d0)
         wcc_old= dmod(wcc_old+ 10d0- maxgap_old, 1d0)
         call sortheap(NumberofSelectedOccupiedBands, wcc_new)
         call sortheap(NumberofSelectedOccupiedBands, wcc_old)

         !>> check the convergence
         max_diff= -1d0
         do i=1, NumberofSelectedOccupiedBands
            call dis_phase(wcc_new(i), wcc_old(i), dis)
            if (dis>max_diff) max_diff= dis
         enddo


      Enddo  !> increase the k meth for integration

      if (cpuid.eq.0) write(stdout, *)'Wcc integration max_diff, Nk_adaptive', max_diff, Nk_adaptive


      maxgap0= -99999d0
      imax= NumberofSelectedOccupiedBands
      do i=1, NumberofSelectedOccupiedBands
         if (i/=NumberofSelectedOccupiedBands) then
            maxgap= wcc(i+1)- wcc(i)
         else
            maxgap=1d0+ wcc(1)- wcc(NumberofSelectedOccupiedBands)
         endif
         wcc_gap(i)= maxgap

         if (maxgap>maxgap0) then
            maxgap0= maxgap
            imax= i
         endif
      enddo

      if (imax==NumberofSelectedOccupiedBands) then
         largestgap_pos_val= (wcc(1)+ &
            wcc(NumberofSelectedOccupiedBands) +1d0)/2d0
         largestgap_pos_val= mod(largestgap_pos_val, 1d0)
      else
         largestgap_pos_val= (wcc(imax+1)+ &
            wcc(imax))/2d0
      endif
      largestgap_pos_i= imax
      largestgap_val= maxgap0


      return
   end subroutine Wcc_integrate_func


   subroutine kline_integrate_sorted(kline_integrate, Nk_adaptive)
      ! sort the kline according to delta
      ! use bubble sorting algorithm
      use para
      use wcc_module
      use wmpi
      implicit none

      type(kline_integrate_type) :: kline_integrate(Nk_adaptive)
      integer, intent(in) :: Nk_adaptive

      integer :: i, j
      real(dp) :: k(3)  ! coordinate
      real(dp) :: delta ! apart from the start point
      real(dp) :: b(3)  ! dis
      logical  :: calculated
      complex(dp), allocatable :: eig_vec(:, :)  !dim= (num_wann, NumberofSelectedOccupiedBands)

      allocate(eig_vec(num_wann, NumberofSelectedOccupiedBands))
      eig_vec= 0d0

      do i= Nk_adaptive-1, 1, -1  ! > start n-1 sweeping
         do j=1, i
            if ( kline_integrate(j)%delta> kline_integrate(j+1)%delta) then
               !> first swap delta
               delta= kline_integrate(j)%delta
               kline_integrate(j)%delta= kline_integrate(j+1)%delta
               kline_integrate(j+1)%delta= delta

               !> then swap k, calculated, eig_vec
               k= kline_integrate(j)%k
               kline_integrate(j)%k= kline_integrate(j+1)%k
               kline_integrate(j+1)%k= k

               calculated= kline_integrate(j)%calculated
               kline_integrate(j)%calculated= kline_integrate(j+1)%calculated
               kline_integrate(j+1)%calculated= calculated

               eig_vec= kline_integrate(j)%eig_vec
               kline_integrate(j)%eig_vec= kline_integrate(j+1)%eig_vec
               kline_integrate(j+1)%eig_vec= eig_vec
            endif

         enddo ! j
      enddo ! i

      do i=1, Nk_adaptive
         if (i==Nk_adaptive) then
            b= kline_integrate(2)%k- kline_integrate(1)%k
            b= b(1)*Origin_cell%Kua+b(2)*Origin_cell%Kub+b(3)*Origin_cell%Kuc
            kline_integrate(i)%b = b
         else
            b= kline_integrate(i+ 1)%k- kline_integrate(i)%k
            b= b(1)*Origin_cell%Kua+b(2)*Origin_cell%Kub+b(3)*Origin_cell%Kuc
            kline_integrate(i)%b= b
         endif
      enddo

      deallocate(eig_vec)

      return
   end subroutine kline_integrate_sorted


   subroutine kline_wcc_sorted(kline_wcc, Nk_adaptive)
      ! sort the kline according to delta
      ! use bubble sorting algorithm
      use para
      use wcc_module
      use wmpi
      implicit none

      type(kline_wcc_type) :: kline_wcc(Nk_adaptive)
      integer, intent(in) :: Nk_adaptive

      integer :: i, j
      real(dp) :: k(3)  ! coordinate
      real(dp) :: delta ! apart from the start point
      logical  :: calculated
      logical  :: converged
      real(dp) :: largestgap_pos_i     ! largest gap position
      real(dp) :: largestgap_pos_val   ! largest gap position value
      real(dp) :: largestgap_val       ! largest gap value
      real(dp), allocatable :: wcc(:)  !dim= (NumberofSelectedOccupiedBands)
      real(dp), allocatable :: gap(:)  !dim= (NumberofSelectedOccupiedBands)

      allocate(wcc(NumberofSelectedOccupiedBands))
      allocate(gap(NumberofSelectedOccupiedBands))
      wcc= 0d0
      gap= 0d0

      do i= Nk_adaptive-1, 1, -1  ! > start n-1 sweeping
         do j=1, i
            if ( kline_wcc(j)%delta> kline_wcc(j+1)%delta) then
               !> first swap delta
               delta= kline_wcc(j)%delta
               kline_wcc(j)%delta= kline_wcc(j+1)%delta
               kline_wcc(j+1)%delta= delta

               !> then swap k, calculated, wcc, gap
               k= kline_wcc(j)%k
               kline_wcc(j)%k= kline_wcc(j+1)%k
               kline_wcc(j+1)%k= k

               calculated= kline_wcc(j)%calculated
               kline_wcc(j)%calculated= kline_wcc(j+1)%calculated
               kline_wcc(j+1)%calculated= calculated

               converged= kline_wcc(j)%converged
               kline_wcc(j)%converged= kline_wcc(j+1)%converged
               kline_wcc(j+1)%converged= converged

               wcc= kline_wcc(j)%wcc
               kline_wcc(j)%wcc= kline_wcc(j+1)%wcc
               kline_wcc(j+1)%wcc= wcc

               gap= kline_wcc(j)%gap
               kline_wcc(j)%gap= kline_wcc(j+1)%gap
               kline_wcc(j+1)%gap= gap

               !> then swap largestgap_val, largestgap_pos_val, largestgap_pos_i
               largestgap_val= kline_wcc(j)%largestgap_val
               kline_wcc(j)%largestgap_val= kline_wcc(j+1)%largestgap_val
               kline_wcc(j+1)%largestgap_val= largestgap_val

               largestgap_pos_i= kline_wcc(j)%largestgap_pos_i
               kline_wcc(j)%largestgap_pos_i= kline_wcc(j+1)%largestgap_pos_i
               kline_wcc(j+1)%largestgap_pos_i= largestgap_pos_i

               largestgap_pos_val= kline_wcc(j)%largestgap_pos_val
               kline_wcc(j)%largestgap_pos_val= kline_wcc(j+1)%largestgap_pos_val
               kline_wcc(j+1)%largestgap_pos_val= largestgap_pos_val

            endif

         enddo ! j
      enddo ! i

      deallocate(wcc, gap)

      return
   end subroutine kline_wcc_sorted

   subroutine largest_gap_find(wcc, maxgap0)

      use para, only : dp, NumberofSelectedOccupiedBands
      use wmpi
      implicit none

      real(dp), intent(in) :: wcc(NumberofSelectedOccupiedBands)
      real(dp) :: maxgap0, maxgap

      integer :: i

      maxgap0= -99999d0
      do i=1, NumberofSelectedOccupiedBands
         if (i/=NumberofSelectedOccupiedBands) then
            maxgap= wcc(i+1)- wcc(i)
         else
            maxgap=1d0+ wcc(1)- wcc(NumberofSelectedOccupiedBands)
         endif

         if (maxgap>maxgap0) then
            maxgap0= maxgap
         endif
      enddo

      return
   end subroutine largest_gap_find

   subroutine dis_phase(x, y, dis)
      use para, only : dp
      use wmpi
      real(dp), intent(in) :: x, y
      real(dp), intent(out) :: dis

      real(dp) :: a, b
      a=x
      b=y

      a= dmod(a, 1d0)
      b= dmod(b, 1d0)

      dis= min( dmod(abs(1d0+ a- b), 1d0), dmod( abs(1d0- a+ b), 1d0))

      return
   end subroutine dis_phase
