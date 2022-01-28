  subroutine sigma_AHC
     !> Calculate anomalous hall conductivity AHC
     !
     !> ref : Physical Review B 74, 195118(2006)
     !
     !> eqn (35)
     !
     !> Dec. 05 2017 by Quansheng Wu @ EPFL
     !> modified on June. 26 2018 by Quansheng Wu @ Airplane from Beijing to Zurich
     !
     ! Copyright (c) 2018 QuanSheng Wu. All rights reserved.

     use wmpi
     use para
     implicit none
    
     integer :: iR, ik, ikx, iky, ikz
     integer :: m, n, i, j, ie
     integer :: ierr, knv3

     real(dp) :: kdotr, mu, Beta_fake
     real(dp) :: k(3)

     real(dp) :: time_start, time_end

     ! eigen value of H
     real(dp), allocatable :: W(:)
     complex(dp), allocatable :: Hamk_bulk(:, :)
     complex(dp), allocatable :: UU(:, :)

     complex(dp) :: ratio

     !> conductivity  dim= OmegaNum
     real(dp), allocatable :: energy(:)
     real(dp), allocatable :: sigma_tensor_ahc(:, :)
     real(dp), allocatable :: sigma_tensor_ahc_mpi(:, :)
     
     !> Fermi-Dirac distribution
     real(dp), external :: fermi

     !> D_mn^H=V_mn/(En-Em) for m!=n
     !> D_nn^H=0 
     complex(dp), allocatable :: Dmn_Ham(:, :, :)
     complex(dp), allocatable :: Vmn_Ham(:, :, :)

     !> Berry curvature vectors for all bands
     real(dp),allocatable :: Omega_BerryCurv(:, :)
     real(dp),allocatable :: Omega_BerryCurv_t(:, :)

     allocate(Dmn_Ham(Num_wann, Num_wann, 3))
     allocate(Vmn_Ham(Num_wann, Num_wann, 3))
     allocate(Omega_BerryCurv(Num_wann, 3))
     allocate(Omega_BerryCurv_t(Num_wann, 3))

     allocate( W (Num_wann))
     allocate( Hamk_bulk(Num_wann, Num_wann))
     allocate( UU(Num_wann, Num_wann))
     allocate( energy(OmegaNum))
     allocate( sigma_tensor_ahc    (3, OmegaNum))
     allocate( sigma_tensor_ahc_mpi(3, OmegaNum))
     sigma_tensor_ahc    = 0d0
     sigma_tensor_ahc_mpi= 0d0
     Hamk_bulk=0d0
     UU= 0d0
     
     !> energy
     do ie=1, OmegaNum
        if (OmegaNum>1) then
           energy(ie)= OmegaMin+ (OmegaMax-OmegaMin)* (ie-1d0)/dble(OmegaNum-1)
        else
           energy= OmegaMin
        endif
     enddo ! ie

     knv3= Nk1*Nk2*Nk3

     call now(time_start) 
     do ik= 1+ cpuid, knv3, num_cpu
        if (cpuid.eq.0.and. mod(ik/num_cpu, 100).eq.0) then
           call now(time_end) 
           write(stdout, '(a, i18, "/", i18, a, f10.2, "s")') 'ik/knv3', &
           ik, knv3, '  time left', (knv3-ik)*(time_end-time_start)/num_cpu/100d0
           time_start= time_end
        endif

        ikx= (ik-1)/(nk2*nk3)+1
        iky= ((ik-1-(ikx-1)*Nk2*Nk3)/nk3)+1
        ikz= (ik-(iky-1)*Nk3- (ikx-1)*Nk2*Nk3)
        k= K3D_start_cube+ K3D_vec1_cube*(ikx-1)/dble(nk1)  &
         + K3D_vec2_cube*(iky-1)/dble(nk2)  &
         + K3D_vec3_cube*(ikz-1)/dble(nk3)

        ! calculation bulk hamiltonian by a direct Fourier transformation of HmnR
        call ham_bulk_atomicgauge(k, Hamk_bulk)
       !call ham_bulk_latticegauge(k, Hamk_bulk)
   
        !> diagonalization by call zheev in lapack
        UU=Hamk_bulk
        call eigensystem_c( 'V', 'U', Num_wann, UU, W)
  
        !> get velocity operator in Hamiltonian basis
        call dHdk_atomicgauge_Ham(k, UU, Vmn_Ham)
       !call dHdk_latticegauge_Ham(k, W, UU, Vmn_Ham)

        call get_Dmn_Ham(W, Vmn_Ham, Dmn_Ham)

        !> calculate Berry curvature at a single k point for all bands
        !> \Omega_n^{\gamma}(k)=i\sum_{\alpha\beta}\epsilon_{\gamma\alpha\beta}(D^{\alpha\dag}D^{\beta})_{nn}
        call berry_curvarture_singlek_allbands(Dmn_Ham, Omega_BerryCurv)
 
        !> consider the Fermi-distribution according to the broadening Earc_eta
        Beta_fake= 1d0/Eta_Arc

        do ie=1, OmegaNum
           mu = energy(ie)
           do m= 1, Num_wann
              Omega_BerryCurv_t(m, :)= Omega_BerryCurv(m, :)*fermi(W(m)-mu, Beta_fake)
           enddo
           sigma_tensor_ahc_mpi(:, ie)= sigma_tensor_ahc_mpi(:, ie)+ &
              (sum(Omega_BerryCurv_t(:, :), dim=1))
        enddo ! ie
     enddo ! ik

#if defined (MPI)
     call mpi_allreduce(sigma_tensor_ahc_mpi,sigma_tensor_ahc,size(sigma_tensor_ahc),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
     sigma_tensor_ahc= sigma_tensor_ahc_mpi
#endif

     !> in the version older than 2.5.2, we use Angstrom as the unit as length
     !sigma_tensor_ahc= sigma_tensor_ahc/dble(knv3)/Origin_cell%CellVolume*24341d0*kCubeVolume/Origin_cell%ReciprocalCellVolume  ! in (Omega*cm)^-1
     
     !> in the latest version, we use the atomic unit
     !> in unit of (Ohm*m)^-1
     sigma_tensor_ahc= sigma_tensor_ahc/dble(knv3)/Origin_cell%CellVolume*Echarge**2/hbar/Bohr_radius*kCubeVolume/Origin_cell%ReciprocalCellVolume
     !> in unit of (Ohm*cm)^-1
     sigma_tensor_ahc= sigma_tensor_ahc/100d0


     outfileindex= outfileindex+ 1
     if (cpuid.eq.0) then
        open(unit=outfileindex, file='sigma_ahc.txt')
        write(outfileindex, '("#",a)')' Anomalous hall conductivity in unit of (Ohm*cm)^-1'
        write(outfileindex, '("#",a13, 20a16)')'Eenergy (eV)', '\sigma_xy', '\sigma_yz', '\sigma_zx'
        do ie=1, OmegaNum
           write(outfileindex, '(200E16.8)')energy(ie)/eV2Hartree, sigma_tensor_ahc(3, ie), &
                                                        sigma_tensor_ahc(1, ie), &
                                                        sigma_tensor_ahc(2, ie)

        enddo
        close(outfileindex)
     endif

     !> write script for gnuplot
     outfileindex= outfileindex+ 1
     if (cpuid==0) then
        open(unit=outfileindex, file='sigma_ahc.gnu')
        write(outfileindex, '(a)') 'set terminal pdf enhanced color font ",20"'
        write(outfileindex, '(a)')"set output 'sigma_ahc.pdf'"
        write(outfileindex, '(a)')'set key samplen 0.8'
        write(outfileindex, '(a)')'set ylabel offset 0.0,0'
        write(outfileindex, '(a, f10.5, a, f10.5, a)')'set xrange [', OmegaMin/eV2Hartree, ':', OmegaMax/eV2Hartree, ']'
        write(outfileindex, '(a)')'set xlabel "Energy (eV)"'
        write(outfileindex, '(a)')'set ylabel "\sigma (1/(Ohm*cm))"'
        write(outfileindex, '(2a)')"plot 'sigma_ahc.txt' u 1:2 w l title '\sigma_{xy}' lc rgb 'red' lw 4, \"
        write(outfileindex, '(2a)')"     'sigma_ahc.txt' u 1:3 w l title '\sigma_{yz}' lc rgb 'blue' lw 4, \"
        write(outfileindex, '(2a)')"     'sigma_ahc.txt' u 1:4 w l title '\sigma_{zx}' lc rgb 'orange' lw 4 "
        close(outfileindex)
     endif

202 format('set xtics (',20('"',A3,'" ',F10.5,','))
203 format(A3,'" ',F10.5,')')
204 format('set arrow from ',F10.5,',',A5,' to ',F10.5,',',A5, ' nohead')


     deallocate( W, Hamk_bulk, UU, energy)
     deallocate( sigma_tensor_ahc, sigma_tensor_ahc_mpi)
 
     return
  end subroutine sigma_AHC
