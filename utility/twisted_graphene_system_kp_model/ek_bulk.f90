subroutine ek_bulk_line
   ! Calculate bulk's energy bands for twisted graphene system using BM continuum model
   ! Line mode
   ! Copied from WannierTools  https://github.com/quanshengwu/wannier_tools
   !> Author: Q.S Wu (wuquansheng@gmail.com)
   ! Copyright (c) 2020 QuanSheng Wu. All rights reserved.

   use wmpi
   use para

   implicit none

   real(dp) :: km(2)

   real(Dp), allocatable :: W(:)
   integer :: ik, io, i, j, knv3, ierr, ilayer, valley, il, iu
   real(dp) :: emin, emax

   ! Hamiltonian of bulk system
   complex(Dp), allocatable :: hmnk(:, :), eigvec(:, :)

   ! eigen value of H
   real(dp), allocatable :: eigv(:,:,:), eigv_mpi(:,:,:)
   real(dp), allocatable :: weight(:,:,:,:), weight_mpi(:,:,:,:)
   integer, external :: stacking_chirality

   il= Ndim/2-Num_bands/2+1
   iu= Ndim/2+Num_bands/2

   knv3= nk_band
   allocate(W(Ndim))
   allocate(hmnk(Ndim, Ndim))
   allocate( eigvec(Ndim, Num_bands))
   allocate( eigv    (Num_bands, knv3, 2))
   allocate( eigv_mpi(Num_bands, knv3, 2))
   allocate( weight    (number_layers, Num_bands, knv3, 2))
   allocate( weight_mpi(number_layers, Num_bands, knv3, 2))
   eigvec=0d0
   eigv    = 0d0
   eigv_mpi= 0d0
   weight = 0d0
   weight_mpi = 0d0

   do ik= 1+cpuid, knv3, num_cpu
      if (cpuid.eq.0) write(stdout, '(a, 2i8)')' ik, knv3 :', ik, knv3

      km = kpoints(:, ik)
      
      hmnk= 0d0

      !> first calculate band structure of K valley
      valley=1
      call get_hk(km, valley, hmnk) 
      W= 0d0
     !call eigensystem_c('V', 'U', Ndim ,hmnk, W)
      call zheevx_pack('V', 'U', Ndim, il, iu, hmnk, W, eigvec)
      eigv(:, ik, 1)= W(1:Num_bands)
      do j=1, Num_bands  !> band
         do ilayer=1, number_layers
            do i=1, Num_Qvectors*2
               io= i+ (ilayer-1)*Num_Qvectors*2
               weight(ilayer, j, ik, 1)= weight(ilayer, j, ik, 1)+ abs(eigvec(io, j))**2 
            enddo
         enddo
      enddo 

      !> Then calculate band structure of K' valley
      hmnk= 0d0
      valley=-1
      call get_hk(km, valley, hmnk) 
      W= 0d0
     !call eigensystem_c('V', 'U', Ndim ,hmnk, W)
      call zheevx_pack('V', 'U', Ndim, il, iu, hmnk, W, eigvec)
      eigv(:, ik, 2)= W(1:Num_bands)
      do j=1, Num_bands  !> band
         do ilayer=1, number_layers
            do i=1, Num_Qvectors*2
               io= i+ (ilayer-1)*Num_Qvectors*2
               weight(ilayer, j, ik, 2)= weight(ilayer, j, ik, 2)+ abs(eigvec(io, j))**2 
            enddo
         enddo
      enddo 
   enddo ! ik

#if defined (MPI)
   call mpi_allreduce(eigv,eigv_mpi,size(eigv),&
      mpi_dp,mpi_sum,mpi_cmw,ierr)
   call mpi_allreduce(weight, weight_mpi,size(weight),&
      mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
   eigv_mpi= eigv
   weight_mpi= weight
#endif
   emin= minval(eigv_mpi)
   emax= maxval(eigv_mpi)

   outfileindex= outfileindex+ nklines+1
   if (cpuid==0)then
      open(unit=outfileindex, file='bulkek.dat-valley-K-matlab')

      write(outfileindex, "('%column', i5, 200i16)")(i, i=1, 1+Num_bands)
      do ik=1, knv3
         write(outfileindex, '(200f16.9)')klen(ik),eigv_mpi(:, ik, 1)
      enddo
      close(outfileindex)
   endif


   outfileindex= outfileindex+ nklines+1
   if (cpuid==0)then
      open(unit=outfileindex, file='bulkek.dat-valley-K')

      write(outfileindex, &
         "('#', a12, a14, 1X, '| projection', 100(3X,'|layer', i2, ': A '))")&
         'klen', 'E', (i, i=1, number_layers)
      write(outfileindex, "('#column', i5, 200i16)")(i, i=1, 2+number_layers)
      do i=1, Num_bands
         do ik=1, knv3
            write(outfileindex, '(200f16.9)')klen(ik),eigv_mpi(i, ik, 1), &
               weight_mpi(:, i, ik, 1)
         enddo
         write(outfileindex, *)' '
      enddo
      close(outfileindex)
   endif


   outfileindex= outfileindex+ nklines+1
   if (cpuid==0)then
      open(unit=outfileindex, file='bulkek.dat-valley-Kprime')

      write(outfileindex, &
         "('#', a12, a14, 1X, '| projection', 100(3X,'|layer', i2, ': A '))")&
         'klen', 'E', (i, i=1, number_layers)
      write(outfileindex, "('#column', i5, 200i16)")(i, i=1, 2+number_layers)
      do i=1, Num_bands
         do ik=1, knv3
            write(outfileindex, '(200f16.9)')klen(ik),eigv_mpi(i, ik, 2), &
               weight_mpi(:, i, ik, 2)
         enddo
         write(outfileindex, *)' '
      enddo
      close(outfileindex)
   endif

   !> write script for gnuplot
   outfileindex= outfileindex+ 1
   if (cpuid==0) then
      open(unit=outfileindex, file='bulkek.gnu')
      write(outfileindex, '(a)') 'set terminal pdf enhanced color font ",30" size 4, 5'
      write(outfileindex,'(2a)') 'set palette defined ( 0  "green", ', &
         '5 "yellow", 10 "red" )'
      write(outfileindex, '(a)')"set output 'bulkek.pdf'"
      write(outfileindex, '(a)')'set style data linespoints'
      write(outfileindex, '(a)')'unset ztics'
      write(outfileindex, '(a)')'unset key'
      write(outfileindex, '(a)')'set pointsize 0.8'
      write(outfileindex, '(a)')'set view 0,0'
      write(outfileindex, '(a)')'set xtics font ",24"'
      write(outfileindex, '(a)')'set ytics font ",24"'
      write(outfileindex, '(a)')'set ylabel font ",24"'
      write(outfileindex, '(a)')'set ylabel offset 1.5,0'
      write(outfileindex, '(a, f10.5, a)')'set xrange [0: ', maxval(klen), ']'
      write(outfileindex, '(a,f12.6)')'emin=', emin
      write(outfileindex, '(a,f12.6)')'emax=', emax
      write(outfileindex, '(a)')'set ylabel "Energy (eV)"'
      write(outfileindex, '(a)')'set yrange [ emin : emax ]'
      write(outfileindex, 202, advance="no") (kline_name(i), kline_stop(i), i=1, nklines)
      write(outfileindex, 203)kline_name(nklines+1), kline_stop(nklines+1)

      do i=1, nklines-1
         write(outfileindex, 204)kline_stop(i+1), 'emin', kline_stop(i+1), 'emax'
      enddo
      write(outfileindex, '(2a)')"# please comment the following lines to plot the fatband "
      write(outfileindex, '(2a)')"plot 'bulkek.dat-valley-K' u 1:2 ",  &
         " w lp lw 2 pt 7  ps 0.1 lc rgb 'blue', \"
      write(outfileindex, '(2a)')"     'bulkek.dat-valley-Kprime' u 1:2 ",  &
         " w lp lw 2 pt 7  ps 0.1 lc rgb 'red'"
      write(outfileindex, '(2a)')" " 
      write(outfileindex, '(2a)')"# uncomment the following lines to plot the fatband "
      write(outfileindex, '(a)')'#set cbrange [0:1]'
      write(outfileindex, '(a)')'#set cbtics 1'
      write(outfileindex, '(2a)')"#plot 'bulkek.dat-valley-K' u 1:2:3 ",  &
         " w lp lw 2 pt 7  ps 0.1 lc palette"
      close(outfileindex)
   endif

202 format('set xtics (',20('"',A3,'" ',F10.5,','))
203 format(A3,'" ',F10.5,')')
204 format('set arrow from ',F10.5,',',A5,' to ',F10.5,',',A5, ' nohead')

   deallocate(W)
   deallocate(hmnk)
   deallocate( eigv    )
   deallocate( eigv_mpi)
   deallocate( weight    )
   deallocate( weight_mpi)

   return
end subroutine ek_bulk_line


