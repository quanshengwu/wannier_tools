!-------------------------------------------------------------------!
!> This file contains subroutines that calculate anomalous          !
!> transport properties, including anomalous Hall conductivity,     !
!> spin Hall conductivity and anomalous Nernst coefficient.         !
!>                                                                  !
!> References :                                                     !
!> [1] Physical Review B 74, 195118(2006)        (AHC1)             !
!> [2] Physcial Review Letter 95, 156601 (2005)  (SHC1)             !
!> [3] Physcial Review B 98, 214402 (2018)       (SHC2)             !
!> [4] Physical Review Letter 97, 026603 (2006)  (ANC1)             !
!> restructed by Hanqi Pi on May 2023
!-------------------------------------------------------------------!

   subroutine sigma_AHC
   !------------------------------------------------------------------!
   !> This subroutine is to calculate the anomalous Hall conductivity !
   !>                                                                 !
   !> It produces the following output files:                         !
   !>                                                                 !
   !> 1. sigma_ahc_eta***meV.txt: sigma tensor for several eta values !
   !> 2. sigma_ahc.gnu: gnuplot script to plot sigma tensor           !
   !>                                                                 !
   !> The output sigma is in the unit of S/cm.                        !
   !>                                                                 !
   !> References: eq(35) of AHC1                                      !
   !>                                                                 !
   !> Dec. 05 2017 by Quansheng Wu @ EPFL                             !
   !> modified on June. 26 2018 by Quansheng Wu @ Airplane from       !
   !> Beijing to Zurich                                               !
   !------------------------------------------------------------------!

     use wmpi
     use para
     implicit none
    
     integer :: iR, ik, ikx, iky, ikz
     integer :: i, m, ie, ieta, iT
     integer :: ierr, knv3
     integer :: NumberofEta

     !> conductivity  dim= OmegaNum
     real(dp), allocatable :: energy(:)
     real(dp), allocatable :: sigma_tensor_ahc(:, :, :)

     real(dp), allocatable :: eta_array(:)
     real(dp), allocatable :: T_list(:)

     character*40 :: ahcfilename, etaname

     NumberofEta = 9 

     allocate(eta_array(NumberofEta))
     allocate(T_list(NumT))
     allocate( energy(OmegaNum))
     allocate( sigma_tensor_ahc    (3, OmegaNum, NumberofEta))
     sigma_tensor_ahc    = 0d0
      
     eta_array=(/0.1d0, 0.2d0, 0.4d0, 0.8d0, 1.0d0, 2d0, 4d0, 8d0, 10d0/)
     eta_array= eta_array*Fermi_broadening


     !> energy
     do ie=1, OmegaNum
        if (OmegaNum>1) then
           energy(ie)= OmegaMin+ (OmegaMax-OmegaMin)* (ie-1d0)/dble(OmegaNum-1)
        else
           energy= OmegaMin
        endif
     enddo ! ie

     !> temperature
     do iT=1, NumT
         if (NumT>1) then
            T_list(iT)= Tmin+(Tmax-Tmin)*(iT-1d0)/dble(NumT-1)
         else
            T_list= Tmin
         endif
     enddo ! iT

     call  sigma_ahc_vary_ChemicalPotential(OmegaNum, energy, NumberofEta, eta_array, sigma_tensor_ahc)

     outfileindex= outfileindex+ 1
     if (cpuid.eq.0) then
        do ieta=1, NumberofEta
           write(etaname, '(f12.2)')eta_array(ieta)*1000d0/eV2Hartree
           write(ahcfilename, '(7a)')'sigma_ahc_eta', trim(adjustl(etaname)), 'meV.txt'
           open(unit=outfileindex, file=ahcfilename)
           write(outfileindex, '("#",10a)')' Anomalous hall conductivity in unit of S/cm,', 'Brodening eta= ',  trim(adjustl(etaname)), ' meV'
           write(outfileindex, "('#column', i5, 3000i16)")(i, i=1, 4)
           write(outfileindex, '("#",a13, 20a16)')'Eenergy (eV)', '\sigma_xy', '\sigma_yz', '\sigma_zx'
           do ie=1, OmegaNum
              write(outfileindex, '(200E16.8)')energy(ie)/eV2Hartree, sigma_tensor_ahc(3, ie, ieta), &
                                                                      sigma_tensor_ahc(1, ie, ieta), &
                                                                      sigma_tensor_ahc(2, ie, ieta)
           
           enddo ! ie
           close(outfileindex)
        enddo
     endif

     !> write script for gnuplot
     outfileindex= outfileindex+ 1
     if (cpuid==0) then
        write(etaname, '(f12.2)')Fermi_broadening*1000d0/eV2Hartree
        write(ahcfilename, '(7a)')'sigma_ahc_eta', trim(adjustl(etaname)), 'meV.txt'
        open(unit=outfileindex, file='sigma_ahc.gnu')
        write(outfileindex, '(a)') 'set terminal pdf enhanced color font ",20"'
        write(outfileindex, '(a)')"set output 'sigma_ahc.pdf'"
        write(outfileindex, '(a)')'set key samplen 0.8'
        write(outfileindex, '(a)')'set ylabel offset 0.0,0'
        write(outfileindex, '(a, f10.5, a, f10.5, a)')'set xrange [', OmegaMin/eV2Hartree, ':', OmegaMax/eV2Hartree, ']'
        write(outfileindex, '(a)')'set xlabel "Energy (eV)"'
        write(outfileindex, '(a)')'set ylabel "AHC (S/cm)"'
        write(outfileindex, '(5a)')"#plot '",  trim(adjustl(ahcfilename)),  "' u 1:2 w l title '\sigma_{xy}' lc rgb 'red' lw 4, \"
        write(outfileindex, '(5a)')"'",  trim(adjustl(ahcfilename)), "' u 1:3 w l title '\sigma_{yz}' lc rgb 'blue' lw 4, \"
        write(outfileindex, '(5a)')"'", trim(adjustl(ahcfilename)), "' u 1:4 w l title '\sigma_{zx}' lc rgb 'orange' lw 4 "
        write(outfileindex, '(5a)')"plot '",  trim(adjustl(ahcfilename)),  "' u 1:2 w l title '{/Symbol s}_{xy}' lc rgb 'red' lw 4, \"
        write(outfileindex, '(5a)')"'",  trim(adjustl(ahcfilename)), "' u 1:3 w l title '{/Symbol s}_{yz}' lc rgb 'blue' lw 4, \"
        write(outfileindex, '(5a)')"'", trim(adjustl(ahcfilename)), "' u 1:4 w l title '{/Symbol s}_{zx}' lc rgb 'orange' lw 4 "
        close(outfileindex)
     endif

     deallocate(energy, sigma_tensor_ahc)

     return
  end subroutine sigma_AHC

  subroutine alpha_ANE
   !------------------------------------------------------------------!
   !> This subroutine is to calculate the anomalous Nernst coefficient!
   !>                                                                 !
   !> It produces the following output files:                         !
   !>                                                                 !
   !> 1. alpha_ane_eta***meV.txt: ane tensor for several eta values   !
   !> 2. alpha_ane.gnu: gnuplot script to plot ane tensor             !
   !>                                                                 !
   !> The first type of files are output on a grid of (T, mu)         !
   !> T is the temperature and mu is the chemical potential.          !
   !> The output coefficient is in the unit of A/(mK).                !
   !>                                                                 !
   !> References: eq(8) of ANC1                                       !
   !>                                                                 !
   !> Dec 14 2022  by Hanqi Pi @ iopcas                               !
   !> Modified at Apr 30 2023 by Hanqi Pi @ iopcas                    !
   !------------------------------------------------------------------!
   use para 
   implicit none 
   
   ! energy refers to the variable to be integrated
   ! mu refers to the chemical potential
   real(dp) :: energy,  omega, T, mu, time_start, time_end

   real(dp), parameter :: ANE_int_interval = 1.0d0*eV2Hartree
   real(dp), parameter :: ANE_int_step     = 0.001d0*eV2Hartree
   
   real(dp), allocatable :: T_list(:)
   real(dp), allocatable :: mu_list(:) ! chemical potential alpha_ane
   real(dp), allocatable :: energy_list(:) ! chemical potential for sigma_ahc (num_step+1)
   real(dp), allocatable :: nernst(:,:,:,:) ! (coordinate,Temp, mu, NumberofEta)
   real(dp), allocatable :: ahc(:,:,:) ! (coordinate, energy, NumberofEta) 
   real(dp), allocatable :: minus_dfde(:,:,:) ! derivative of Fermi-Dirac (num_step+1, Temp, mu)
   integer :: i, ie, iT, imu,ieta, NumberofEta

   ! integral region is (OmegaMin-ANE_int_interval,OmegaMax+ANE_int_interval)
   ! num_step = (OmegaMax-OmegaMin+2*ANE_int_interval)/ANE_int_step
   ! Note: make sure integral region includes Efermi
   integer :: num_step 

   real(dp), allocatable :: eta_array(:)

   character*40 :: anefilename, etaname

   num_step = int((OmegaMax-OmegaMin+2*ANE_int_interval)/ANE_int_step)
   NumberofEta = 9 

   allocate(eta_array(NumberofEta))
   allocate( T_list(NumT) )
   allocate( mu_list(OmegaNum) )
   allocate( energy_list(num_step+1) )
   allocate( nernst(3, NumT, OmegaNum, NumberofEta)) 
   allocate( ahc(3, num_step+1, NumberofEta)) 
   allocate( minus_dfde(num_step+1, NumT, OmegaNum) )
   T_list = 0.0d0
   energy_list = 0.0d0
   nernst = 0.0d0
   ahc = 0.0d0
   minus_dfde = 0.0d0

   eta_array=(/0.1d0, 0.2d0, 0.4d0, 0.8d0, 1.0d0, 2d0, 4d0, 8d0, 10d0/)
   eta_array= eta_array*Fermi_broadening

   !> energy to be integrated
   do ie = 1,  num_step+1
       energy_list(ie) = (OmegaMin-ANE_int_interval)+ANE_int_step*(ie-1)
   enddo

   !> chemical potential
   do imu=1, OmegaNum
       if (OmegaNum>1) then
          mu_list(imu)= OmegaMin+ (OmegaMax-OmegaMin)* (imu-1d0)/dble(OmegaNum-1)
       else
          mu_list= OmegaMin
       endif
   enddo ! imu

   !> temperature
   do iT=1, NumT
       if (NumT>1) then
          T_list(iT)= Tmin+(Tmax-Tmin)*(iT-1d0)/dble(NumT-1)
       else
          T_list= Tmin
       endif
   enddo ! iT

   !> get the anomalous hall conductivity
   call now(time_start)
   call sigma_ahc_vary_ChemicalPotential(num_step+1, energy_list, NumberofEta, eta_array, ahc)
   call now(time_end)
   call print_time_cost(time_start, time_end, 'calculation of AHC')

   !> get the derivative of Fermi-Dirac
   call now(time_start)
   do imu = 1, OmegaNum
      mu = mu_list(imu)
      do iT = 1, NumT
         T = T_list(iT)
         do ie = 1, num_step+1
            energy = energy_list(ie)
            call minusdfde_calc_single(energy/eV2Hartree, T*8.617333262E-5, &
                                       mu/eV2Hartree, minus_dfde(ie, iT, imu))
         enddo ! imu
      enddo ! iT
   enddo ! ie
   call now(time_end)
   call print_time_cost(time_start, time_end, 'calculation of derivative of Fermi-Dirac')

   ! --------------------------------------------------------------------
   !  nernst = -int dE 1/e* df/dmu * ahc * (E - mu) / T
   !
   !  We use Hartree as the unit of energy and convert them to eV with eV2Hartree
   !  Therefore the unit of nernst is 1/e * eV * 1/eV * S/cm * eV * 1/K = A/cm/K
   !  To convert it to A/m/K, we need to multiply it by 100
   ! --------------------------------------------------------------------

   !> get the anomalous Nernst coefficient
   call now(time_start)
   do ieta = 1, NumberofEta
      do imu = 1, OmegaNum
         mu = mu_list(imu)
         do iT = 1, NumT
            T = T_list(iT)
            do ie = 1, num_step+1
               energy = energy_list(ie)
               omega = energy - mu
               nernst(:, iT, imu,ieta) = nernst(:, iT, imu,ieta)- minus_dfde(ie,iT,imu)*ahc(:,ie,ieta)*&
                                    (omega/eV2Hartree)/T*(ANE_int_step/eV2Hartree)*100
            enddo ! ie
         enddo ! iT
      enddo ! ie
   enddo ! ieta
   call now(time_end)
   call print_time_cost(time_start, time_end, 'integration of AHC')

   if(cpuid .eq. 0) then
      do ieta =1, NumberofEta
         outfileindex = outfileindex+1
         write(etaname,'(f12.2)') eta_array(ieta)*1000d0/eV2Hartree
         write(anefilename,'(7a)') 'alpha_ane_eta',trim(adjustl(etaname)),'meV.txt'
         open(unit=outfileindex, file=anefilename)
         write(outfileindex, '("#",a)')'Anomalous Nernst coefficient  (A(m-1 K-1))'
         write(outfileindex, "('#column', i7, i13,3i16)")(i, i=1, 5)
         write(outfileindex,'("#",2a13, 3a16)')'Mu(eV)', 'T(K)', '\alpha_xy', '\alpha_yz', '\alpha_zx'
         do imu = 1, OmegaNum
            do iT = 1, NumT
                  write(outfileindex,'(F14.6,F13.3,3E16.6)') mu_list(imu)/eV2Hartree,T_list(iT), nernst(3, iT, imu,ieta),&
                                                                                             nernst(1, iT, imu, ieta),&
                                                                                             nernst(2, iT, imu, ieta)
            enddo
         enddo
         close(outfileindex)
      enddo
   endif

   !> write script for gnuplot
   outfileindex = outfileindex+1
   if(cpuid .eq. 0) then
         write(etaname, '(f12.2)')Fermi_broadening*1000d0/eV2Hartree
         write(anefilename, '(7a)')'alpha_ane_eta', trim(adjustl(etaname)), 'meV.txt'
         open(unit=outfileindex, file='alpha_ane.gnu')
         write(outfileindex, '(a)') 'set terminal pdfcairo enhanced color font ",30" size 13, 6'
         write(outfileindex, '(a)') "set output 'ane.pdf'"
         write(outfileindex, '(a)') 'set key outside'
         write(outfileindex, '(a)') "set palette defined (0 'red', 1 'green')"
         write(outfileindex, '(a)') 'unset colorbox'
         write(outfileindex, '(a, f6.2)') 'Tmin = ',Tmin
         write(outfileindex, '(a, f6.2)') 'Tmax = ',Tmax
         write(outfileindex, '(a, I4)') 'NumT = ',NumT
         write(outfileindex, '(a, f6.2)') 'OmegaMin = ',OmegaMin/eV2Hartree
         write(outfileindex, '(a, f6.2)') 'OmegaMax = ',OmegaMax/eV2Hartree
         write(outfileindex, '(a, I4)') 'OmegaNum = ',OmegaNum
         write(outfileindex, '(a)') ''
         write(outfileindex, '(a)') '#plot the temperature-dependent alpha_yx for the first six chemical potentials'
         write(outfileindex, '(a)') 'set xlabel "T (K)"'
         write(outfileindex, '(a)') 'set ylabel "{/Symbol a}_{yx} (A/(mK))"'
         write(outfileindex, '(a)') 'set ylabel offset 0.0,0'
         write(outfileindex, '(5a)')& 
                                    "plot for [i=0:5] '",trim(adjustl(anefilename)),"' every ::i*NumT::i*NumT+(NumT-1)", &
                                    "u 2:(-$3) w l lt palette frac i/5. title sprintf('{/Symbol m}=%.3f eV', ",&
                                    "OmegaMin+(OmegaMax-OmegaMin)/(OmegaNum*1.0-1.0)*i )"
         write(outfileindex, '(a)') ''
         write(outfileindex, '(a)') '#plot the chemical potential dependent alpha_yx'
         write(outfileindex, '(a)') '#set xlabel "E-E_f (eV)"'
         write(outfileindex, '(a)') '#set ylabel "{/Symbol a}_{yx} (A/(mK))"'
         write(outfileindex, '(a)') '#set ylabel offset 0.0,0'
         write(outfileindex, '(5a)')& 
                                    "#plot for [i=0:NumT-1] '",trim(adjustl(anefilename)),"' every NumT::i", &
                                    " u 1:(-$3) w l lt palette frac i/(NumT*1.0-1.0) title sprintf('T=%.3f K',",&
                                    "Tmin+(Tmax-Tmin)/(NumT*1.0-1.0)*i)"
         close(outfileindex)
   endif

   deallocate(mu_list, energy_list, T_list, nernst, ahc)
   return

end subroutine alpha_ANE

subroutine sigma_ahc_vary_ChemicalPotential(NumOfmu, mulist, NumberofEta, eta_array, sigma_tensor_ahc)
   !------------------------------------------------------------------!
   !> This subroutine is to calculate the anomalous Hall conductivity !
   !> as a function of chemical potential                             !
   !>                                                                 !
   !> This subroutine is called by sigma_AHC and alpha_ANE            !
   !>                                                                 !
   !> Dec. 05 2017 by Quansheng Wu @ EPFL                             !
   !> modified on June. 26 2018 by Quansheng Wu @ Airplane from       !
   !> Beijing to Zurich                                               !
   !> modified on Apr. 29. 2023 by Hanqi Pi @ Beijing                 !
   !------------------------------------------------------------------!

     use wmpi
     use para
     implicit none

     integer, intent(in)    :: NumOfmu
     real(dp), intent(in)    :: mulist(NumOfmu)
     integer, intent(in)    :: NumberofEta
     real(dp), intent(in)    :: eta_array(NumberofEta)
     real(dp), intent(inout) :: sigma_tensor_ahc(3, NumOfmu, NumberofEta)
     
     integer :: iR, ik, ikx, iky, ikz
     integer :: i, m, ie, ieta
     integer :: ierr, knv3

     real(dp) :: mu, Beta_fake, eta_local
     real(dp) :: k(3)

     real(dp) :: time_start, time_end

     ! eigen value of H
     real(dp), allocatable :: W(:)
     complex(dp), allocatable :: Hamk_bulk(:, :)
     complex(dp), allocatable :: UU(:, :)

     !> conductivity  dim= OmegaNum
     real(dp), allocatable :: energy(:)
     real(dp), allocatable :: sigma_tensor_ahc_mpi(:, :, :)
     
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
     allocate( sigma_tensor_ahc_mpi(3, NumOfmu, NumberofEta))
     sigma_tensor_ahc_mpi= 0d0
     Hamk_bulk=0d0
     UU= 0d0
      
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
        call Berry_curvature_singlek_allbands(Dmn_Ham, Omega_BerryCurv)
 
        do ieta= 1, NumberofEta
           eta_local = eta_array(ieta)
           !> consider the Fermi-distribution according to the broadening Earc_eta
           Beta_fake= 1d0/eta_local
   
           do ie=1, NumOfmu
              mu = mulist(ie)
              do m= 1, Num_wann
                 Omega_BerryCurv_t(m, :)= Omega_BerryCurv(m, :)*fermi(W(m)-mu, Beta_fake)
              enddo
              sigma_tensor_ahc_mpi(:, ie, ieta)= sigma_tensor_ahc_mpi(:, ie, ieta)- &
                 (sum(Omega_BerryCurv_t(:, :), dim=1))
           enddo ! ie
        enddo ! ieta
     enddo ! ik

#if defined (MPI)
     call mpi_allreduce(sigma_tensor_ahc_mpi,sigma_tensor_ahc,size(sigma_tensor_ahc),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
     sigma_tensor_ahc= sigma_tensor_ahc_mpi
#endif
   ! --------------------------------------------------------------------
   ! At this point ahc contains
   !
   !           sum_k Omega =  -N*V_c*int dk/[8(pi)^3] Omega
   !
   ! (N is the number of kpoints, V_c is the cell volume). We want
   !
   !           ahc = -(e^2/hbar) int dk/[4(pi)^3] Omega
   ! 
   ! in the version older than 2.5.2, we use Angstrom as the unit as length
   ! Hence we need to multiply by  e^2/hbar/N/V_c = 0.00024341/N/V_c
   ! where 'V_c' in Angstroms^3, and 'e', 'hbar' in SI units.
   ! Then it has units of 
   !
   !           C^2 (J*s)^-1 Angstrom^-3  Angstrom^2 = S/Angstrom                         
   ! 
   ! in the latest version, we use the atomic unit
   ! Similary we need to multiply by  e^2/hbar/N/V_c/Bohr_radius 
   ! and obtain the sigma_tensor_ahc in unit of (Ohm*m)^-1
   ! --------------------------------------------------------------------

      !> in the version older than 2.5.2, we use Angstrom as the unit as length
      !sigma_tensor_ahc= sigma_tensor_ahc/dble(knv3)/Origin_cell%CellVolume*24341d0*kCubeVolume/Origin_cell%ReciprocalCellVolume  ! in (Omega*cm)^-1
      
      !> in the latest version, we use the atomic unit
      sigma_tensor_ahc= sigma_tensor_ahc/dble(knv3)/Origin_cell%CellVolume*Echarge**2/hbar/&
         Bohr_radius*kCubeVolume/Origin_cell%ReciprocalCellVolume
      !> in unit of (Ohm*cm)^-1
      sigma_tensor_ahc= sigma_tensor_ahc/100d0
      
      ! deallocate( W, Hamk_bulk, UU)
      ! deallocate( sigma_tensor_ahc_mpi)
  end subroutine sigma_ahc_vary_ChemicalPotential

  subroutine sigma_SHC
   !------------------------------------------------------------------!
   !> This subroutine is to calculate the spin Hall conductivity      !
   !>                                                                 !
   !> It produces the following output files:                         !
   !>                                                                 !
   !> 1. sigma_shc_eta***meV.txt: sigma tensor for several eta values !
   !> 2. sigma_shc.gnu: gnuplot script to plot sigma tensor           !
   !>                                                                 !
   !> The output sigma is in the unit of (hbar/e)S/cm.                !
   !>                                                                 !
   !> References: SHC1 and SHC2                                       !
   !>                                                                 !
   !> Dec. 05 2022 by Quansheng Wu @ Beijing                          !
   !------------------------------------------------------------------!

     use wmpi
     use para
     implicit none
    
     integer :: ik, ikx, iky, ikz, ieta
     integer :: m, n, i, j, ie, ialpha, ibeta, igamma
     integer :: ierr, knv3, nwann
     integer :: NumberofEta

     real(dp) :: mu, Beta_fake, deno_fac, eta_local
     real(dp) :: k(3)

     real(dp) :: time_start, time_end

     ! eigen value of H
     real(dp), allocatable :: W(:)
     complex(dp), allocatable :: Hamk_bulk(:, :)
     complex(dp), allocatable :: UU(:, :)

     !> energy  dim= OmegaNum
     real(dp), allocatable :: energy(:)

     !> sigma^gamma_{alpha, beta}, alpha, beta, gamma=1,2,3 for x, y, z
     !>  sigma_tensor_shc(ie, igamma, ialpha, ibeta, ieta)
     real(dp), allocatable :: sigma_tensor_shc(:, :, :, :, :)
     real(dp), allocatable :: sigma_tensor_shc_mpi(:, :, :, :, :)
     
     !> Fermi-Dirac distribution
     real(dp), external :: fermi

     complex(dp), allocatable :: Vmn_Ham(:, :, :)
     complex(dp), allocatable :: Vmn_wann(:, :, :)
     complex(dp), allocatable :: spin_sigma(:, :)
     complex(dp), allocatable :: j_spin_gamma_alpha(:, :)
     complex(dp), allocatable :: mat_t(:, :)

     !> Berry curvature vectors for all bands
     real(dp),allocatable :: Omega_spin(:)
     real(dp),allocatable :: Omega_spin_t(:)

     ! spin operator matrix spin_sigma_x,spin_sigma_y in spin_sigma_z representation
     complex(Dp),allocatable :: pauli_matrices(:, :, :) 
    
     real(dp), allocatable :: eta_array(:)
     character(80) :: shcfilename, etaname

     NumberofEta=9

     allocate(Vmn_Ham(Num_wann, Num_wann, 3))
     allocate(Vmn_wann(Num_wann, Num_wann, 3))
     allocate(spin_sigma(Num_wann, Num_wann))
     allocate(j_spin_gamma_alpha(Num_wann, Num_wann))
     allocate(Omega_spin(Num_wann))
     allocate(Omega_spin_t(Num_wann))

     allocate( W (Num_wann))
     allocate( Hamk_bulk(Num_wann, Num_wann))
     allocate( UU(Num_wann, Num_wann))
     allocate( mat_t(Num_wann, Num_wann))
     allocate( energy(OmegaNum))
     allocate(eta_array(NumberofEta))
     allocate( sigma_tensor_shc    (OmegaNum, 3, 3, 3, NumberofEta))
     allocate( sigma_tensor_shc_mpi(OmegaNum, 3, 3, 3, NumberofEta))

     allocate(pauli_matrices(Num_wann, Num_wann, 3))
     spin_sigma= 0d0
     sigma_tensor_shc    = 0d0
     sigma_tensor_shc_mpi= 0d0
     Hamk_bulk=0d0
     UU= 0d0
     Vmn_wann= 0d0
     pauli_matrices= 0d0
 
     eta_array=(/0.1d0, 0.2d0, 0.4d0, 0.8d0, 1.0d0, 2d0, 4d0, 8d0, 10d0/)
     eta_array= eta_array*Fermi_broadening

     nwann= Num_wann/2
     !> spin operator matrix
     !> this part is package dependent. 
    !if (index( Package, 'VASP')/=0.or. index( Package, 'Wien2k')/=0 &
    !   .or. index( Package, 'Abinit')/=0.or. index( Package, 'openmx')/=0) then
        do j=1, nwann
           pauli_matrices(j, nwann+j, 1)=1.0d0
           pauli_matrices(j+nwann, j, 1)=1.0d0
           pauli_matrices(j, nwann+j, 2)=-zi
           pauli_matrices(j+nwann, j, 2)=zi
           pauli_matrices(j, j, 3)= 1d0
           pauli_matrices(j+nwann, j+nwann, 3)=-1d0
        enddo
    !else
    !   if (cpuid.eq.0) write(stdout, *)'Error: please report your software and wannier90.wout to me'
    !   if (cpuid.eq.0) write(stdout, *)'wuquansheng@gmail.com'
    !   stop 'Error: please report your software and wannier90.wout to wuquansheng@gmail.com'
    !endif

   
     !> energy range (chemical potential range)
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
        UU= Hamk_bulk
        call eigensystem_c( 'V', 'U', Num_wann, UU, W)
  
        !> get velocity operator in Wannier basis
        !> \partial_k H_nm
        call dHdk_atomicgauge(k, Vmn_wann)

        do ialpha= 1, 3
           call rotation_to_Ham_basis(UU, Vmn_wann(:, :, ialpha), Vmn_Ham(:, :, ialpha))
        enddo

        !> spin axis igamma= x, y, z
        do igamma= 1, 3
           !> set Pauli matrix
           spin_sigma= pauli_matrices(:, :, igamma)
   
           !> calculate spin current operator j_spin_gamma_alpha^l_alpha= 1/2*{Sigma_gamma, v_alpha} 
           ! in order to calculate SHC^z_xy
           do ialpha= 1, 3
              !> in Wannier basis
              j_spin_gamma_alpha= 0d0
              call mat_mul(Num_wann, spin_sigma, Vmn_wann(:, :, ialpha), j_spin_gamma_alpha(:, :))
              call mat_mul(Num_wann, Vmn_wann(:, :, ialpha), spin_sigma, mat_t)
              j_spin_gamma_alpha(:, :)= j_spin_gamma_alpha(:, :)+ mat_t
              j_spin_gamma_alpha= j_spin_gamma_alpha/2d0
         
              !> rotate to Hamiltonian basis
              mat_t= j_spin_gamma_alpha(:, :)
              call rotation_to_Ham_basis(UU, mat_t, j_spin_gamma_alpha(:, :))

              do ieta= 1, NumberofEta
                 eta_local= eta_array(ieta)
                 !> \Omega_spin^l_n^{\gamma}(k)=-2\sum_{m}*aimag(Im({js(\gamma),v(\alpha)}/2)_nm*v_beta_mn))/((w(n)-w(m))^2+Fermi_broadening^2)
                 do ibeta= 1, 3
                    Omega_spin= 0d0
                    do n= 1, Num_wann
                       do m= 1, Num_wann
                          if (abs(W(n)-W(m))<eps9 ) cycle
                          deno_fac= -2d0/((W(n)-W(m))**2+ eta_local**2)
                          Omega_spin(n)= Omega_spin(n)+ &
                             aimag(j_spin_gamma_alpha(n, m)*Vmn_Ham(m, n, ibeta))*deno_fac
                       enddo
                    enddo
        
                    !> consider the Fermi-distribution according to the broadening Earc_eta
                    Beta_fake= 1d0/eta_local
            
                    do ie=1, OmegaNum
                       mu = energy(ie)
                       do n= 1, Num_wann
                          Omega_spin_t(n)= Omega_spin(n)*fermi(W(n)-mu, Beta_fake)
                       enddo
   
                       !> sum over all "spin" Berry curvature below chemical potential mu
                       sigma_tensor_shc_mpi(ie, igamma, ialpha, ibeta, ieta)= &
                          sigma_tensor_shc_mpi(ie, igamma, ialpha, ibeta, ieta)+ &
                          sum(Omega_spin_t(:))
                    enddo ! ie
                 enddo ! ibeta  v
              enddo ! ieta=1, NumberofEta
           enddo ! ialpha  j
        enddo ! igamma  spin
     enddo ! ik

#if defined (MPI)
     call mpi_allreduce(sigma_tensor_shc_mpi,sigma_tensor_shc,size(sigma_tensor_shc),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
     sigma_tensor_shc= sigma_tensor_shc_mpi
#endif

     !> in the latest version, we use the atomic unit
     !> in unit of ((hbar/e)(Ohm*m)^-1
     sigma_tensor_shc= sigma_tensor_shc/dble(knv3)/Origin_cell%CellVolume*&
        Echarge**2/hbar/Bohr_radius*kCubeVolume/Origin_cell%ReciprocalCellVolume/2d0

     !> in unit of ((hbar/e)(Ohm*cm)^-1
     sigma_tensor_shc= sigma_tensor_shc/100d0


     outfileindex= outfileindex+ 1
     if (cpuid.eq.0) then
        do ieta=1, NumberofEta
           write(etaname, '(f12.2)')eta_array(ieta)*1000d0/eV2Hartree
           write(shcfilename, '(7a)')'sigma_shc_eta', trim(adjustl(etaname)), 'meV.txt'
           open(unit=outfileindex, file=shcfilename)
           write(outfileindex, '("#",10a)')' Spin hall conductivity in unit of (hbar/e)S/cm,', 'Brodening eta= ',  trim(adjustl(etaname)), ' meV'
           write(outfileindex, "('#column', i5, 3000i16)")(i, i=1, 28)
           write(outfileindex, '("#",a13, 27a16)')'Eenergy (eV)', &
             'xx^x', 'xy^x', 'xz^x', 'yx^x', 'yy^x', 'yz^x', 'zx^x', 'zy^x', 'zz^x', &
             'xx^y', 'xy^y', 'xz^y', 'yx^y', 'yy^y', 'yz^y', 'zx^y', 'zy^y', 'zz^y', &
             'xx^z', 'xy^z', 'xz^z', 'yx^z', 'yy^z', 'yz^z', 'zx^z', 'zy^z', 'zz^z'
           do ie=1, OmegaNum
              write(outfileindex, '(E16.8)', advance='no')energy(ie)/eV2Hartree
              do igamma=1, 3
              do ialpha=1, 3
              do ibeta =1, 3
                 if (ialpha*ibeta*igamma/=27)then
                    write(outfileindex, '(200E16.8)', advance='no') sigma_tensor_shc(ie, igamma, ialpha, ibeta, ieta)
                 else
                    write(outfileindex, '(200E16.8)', advance='yes') sigma_tensor_shc(ie, igamma, ialpha, ibeta, ieta)
                 endif
              enddo
              enddo
              enddo
           enddo
        enddo ! ieta
        close(outfileindex)
     endif

     !> write script for gnuplot
     outfileindex= outfileindex+ 1
     if (cpuid==0) then
        write(etaname, '(f12.2)')Fermi_broadening*1000d0/eV2Hartree
        write(shcfilename, '(7a)')'sigma_shc_eta', trim(adjustl(etaname)), 'meV.txt'
        open(unit=outfileindex, file='sigma_shc.gnu')
        write(outfileindex, '(a)') 'set terminal pdf enhanced color font ",20"'
        write(outfileindex, '(a)')"set output 'sigma_shc.pdf'"
        write(outfileindex, '(a)')'set key samplen 0.8'
        write(outfileindex, '(a)')'set ylabel offset 0.0,0'
        write(outfileindex, '(a, f10.5, a, f10.5, a)')'set xrange [', OmegaMin/eV2Hartree, ':', OmegaMax/eV2Hartree, ']'
        write(outfileindex, '(a)')'set xlabel "Energy (eV)"'
        write(outfileindex, '(a)')'set ylabel "SHC (\hbar/e)S/cm"'
        write(outfileindex, '(5a)')"#plot '",  trim(adjustl(shcfilename)),  "' u 1:21 w l title '\sigma_{xy}^z' lc rgb 'red' lw 4, \"
        write(outfileindex, '(5a)')"'",  trim(adjustl(shcfilename)), "' u 1:17 w l title '\sigma_{zx}^y' lc rgb 'blue' lw 4, \"
        write(outfileindex, '(5a)')"'", trim(adjustl(shcfilename)), "' u 1:13 w l title '\sigma_{xz}^y' lc rgb 'orange' lw 4 "
        write(outfileindex, '(5a)')"plot '",  trim(adjustl(shcfilename)),  "' u 1:21 w l title '{/Symbol s}_{xy}^z' lc rgb 'red' lw 4, \"
        write(outfileindex, '(5a)')"'",  trim(adjustl(shcfilename)), "' u 1:17 w l title '{/Symbol s}_{zx}^y' lc rgb 'blue' lw 4, \"
        write(outfileindex, '(5a)')"'", trim(adjustl(shcfilename)), "' u 1:13 w l title '{/Symbol s}_{xz}^y' lc rgb 'orange' lw 4 "
        close(outfileindex)
     endif



     deallocate( W, Hamk_bulk, UU, energy)
     deallocate( sigma_tensor_shc, sigma_tensor_shc_mpi)
 
     return
  end subroutine sigma_SHC
