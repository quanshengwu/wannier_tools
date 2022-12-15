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
     integer :: i, m, ie, ieta
     integer :: ierr, knv3
     integer :: NumberofEta

     real(dp) :: mu, Beta_fake, eta_local
     real(dp) :: k(3)

     real(dp) :: time_start, time_end

     ! eigen value of H
     real(dp), allocatable :: W(:)
     complex(dp), allocatable :: Hamk_bulk(:, :)
     complex(dp), allocatable :: UU(:, :)

     !> conductivity  dim= OmegaNum
     real(dp), allocatable :: energy(:)
     real(dp), allocatable :: sigma_tensor_ahc(:, :, :)
     real(dp), allocatable :: sigma_tensor_ahc_mpi(:, :, :)
     
     !> Fermi-Dirac distribution
     real(dp), external :: fermi

     real(dp), allocatable :: eta_array(:)

     !> D_mn^H=V_mn/(En-Em) for m!=n
     !> D_nn^H=0 
     complex(dp), allocatable :: Dmn_Ham(:, :, :)
     complex(dp), allocatable :: Vmn_Ham(:, :, :)

     !> Berry curvature vectors for all bands
     real(dp),allocatable :: Omega_BerryCurv(:, :)
     real(dp),allocatable :: Omega_BerryCurv_t(:, :)

     character*40 :: ahefilename, etaname

     NumberofEta = 9 

     allocate(Dmn_Ham(Num_wann, Num_wann, 3))
     allocate(Vmn_Ham(Num_wann, Num_wann, 3))
     allocate(Omega_BerryCurv(Num_wann, 3))
     allocate(Omega_BerryCurv_t(Num_wann, 3))

     allocate(eta_array(NumberofEta))
     allocate( W (Num_wann))
     allocate( Hamk_bulk(Num_wann, Num_wann))
     allocate( UU(Num_wann, Num_wann))
     allocate( energy(OmegaNum))
     allocate( sigma_tensor_ahc    (3, OmegaNum, NumberofEta))
     allocate( sigma_tensor_ahc_mpi(3, OmegaNum, NumberofEta))
     sigma_tensor_ahc    = 0d0
     sigma_tensor_ahc_mpi= 0d0
     Hamk_bulk=0d0
     UU= 0d0
      
     eta_array=(/0.1d0, 0.2d0, 0.4d0, 0.8d0, 1.0d0, 2d0, 4d0, 8d0, 10d0/)
     eta_array= eta_array*Eta_Arc


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
 
        do ieta= 1, NumberofEta
           eta_local = eta_array(ieta)
           !> consider the Fermi-distribution according to the broadening Earc_eta
           Beta_fake= 1d0/eta_local
   
           do ie=1, OmegaNum
              mu = energy(ie)
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

     !> in the version older than 2.5.2, we use Angstrom as the unit as length
     !sigma_tensor_ahc= sigma_tensor_ahc/dble(knv3)/Origin_cell%CellVolume*24341d0*kCubeVolume/Origin_cell%ReciprocalCellVolume  ! in (Omega*cm)^-1
     
     !> in the latest version, we use the atomic unit
     !> e^2/h*1/a0: the factor of unit conversion from atomic to SI 
     !> 1/knv3/CellVolume : the factor of summation of k
     !> 1/knv3/CellVolume\sum_k= \int_dk^3 /(2\pi)^3
     !> in unit of (Ohm*m)^-1
     sigma_tensor_ahc= sigma_tensor_ahc/dble(knv3)/Origin_cell%CellVolume*Echarge**2/hbar/&
        Bohr_radius*kCubeVolume/Origin_cell%ReciprocalCellVolume
     !> in unit of (Ohm*cm)^-1
     sigma_tensor_ahc= sigma_tensor_ahc/100d0


     outfileindex= outfileindex+ 1
     if (cpuid.eq.0) then
        do ieta=1, NumberofEta
           write(etaname, '(f12.2)')eta_array(ieta)*1000d0/eV2Hartree
           write(ahefilename, '(7a)')'sigma_ahe_eta', trim(adjustl(etaname)), 'meV.txt'
           open(unit=outfileindex, file=ahefilename)
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
        write(etaname, '(f12.2)')Eta_Arc*1000d0/eV2Hartree
        write(ahefilename, '(7a)')'sigma_ahe_eta', trim(adjustl(etaname)), 'meV.txt'
        open(unit=outfileindex, file='sigma_ahc.gnu')
        write(outfileindex, '(a)') 'set terminal pdf enhanced color font ",20"'
        write(outfileindex, '(a)')"set output 'sigma_ahc.pdf'"
        write(outfileindex, '(a)')'set key samplen 0.8'
        write(outfileindex, '(a)')'set ylabel offset 0.0,0'
        write(outfileindex, '(a, f10.5, a, f10.5, a)')'set xrange [', OmegaMin/eV2Hartree, ':', OmegaMax/eV2Hartree, ']'
        write(outfileindex, '(a)')'set xlabel "Energy (eV)"'
        write(outfileindex, '(a)')'set ylabel "AHC (S/cm)"'
        write(outfileindex, '(5a)')"plot '",  trim(adjustl(ahefilename)),  "' u 1:2 w l title '\sigma_{xy}' lc rgb 'red' lw 4, \"
        write(outfileindex, '(5a)')"'",  trim(adjustl(ahefilename)), "' u 1:3 w l title '\sigma_{zx}' lc rgb 'blue' lw 4, \"
        write(outfileindex, '(5a)')"'", trim(adjustl(ahefilename)), "' u 1:4 w l title '\sigma_{xz}' lc rgb 'orange' lw 4 "
        close(outfileindex)
     endif


     deallocate( W, Hamk_bulk, UU, energy)
     deallocate( sigma_tensor_ahc, sigma_tensor_ahc_mpi)
 
     return
  end subroutine sigma_AHC


  subroutine sigma_SHC
     !> Calculate spin hall conductivity SHC
     !
     !> refs : 
     ! [1] Y. Yao and Z. Fang, Phys.Rev.Lett.95, 156601 (2005).
     ! [2] Junfeng Qiao et al., Plys.Rev.B 98, 214402 (2018)
     !
     !> Dec. 05 2022 by Quansheng Wu @ Beijing
     !
     ! Copyright (c) 2022 QuanSheng Wu. All rights reserved.

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
     eta_array= eta_array*Eta_Arc

     nwann= Num_wann/2
     !> spin operator matrix
     !> this part is package dependent. 
     if (index( Package, 'VASP')/=0.or. index( Package, 'Wien2k')/=0 &
        .or. index( Package, 'Abinit')/=0.or. index( Package, 'openmx')/=0) then
        do j=1, nwann
           pauli_matrices(j, nwann+j, 1)=1.0d0
           pauli_matrices(j+nwann, j, 1)=1.0d0
           pauli_matrices(j, nwann+j, 2)=-zi
           pauli_matrices(j+nwann, j, 2)=zi
           pauli_matrices(j, j, 3)= 1d0
           pauli_matrices(j+nwann, j+nwann, 3)=-1d0
        enddo
     elseif (index( Package, 'QE')/=0.or.index( Package, 'quantumespresso')/=0 &
        .or.index( Package, 'quantum-espresso')/=0.or.index( Package, 'pwscf')/=0) then
        do j=1, nwann
           pauli_matrices((2*j-1), 2*j, 1)=1.0d0
           pauli_matrices(2*j, (2*j-1), 1)=1.0d0
           pauli_matrices((2*j-1), 2*j, 2)=-zi
           pauli_matrices(2*j, (2*j-1), 2)=zi
           pauli_matrices((2*j-1), (2*j-1), 3)=1.0d0
           pauli_matrices(2*j, 2*j, 3)=-1.0d0
        enddo
     else
        if (cpuid.eq.0) write(stdout, *)'Error: please report your software and wannier90.wout to me'
        if (cpuid.eq.0) write(stdout, *)'wuquansheng@gmail.com'
        stop 'Error: please report your software and wannier90.wout to wuquansheng@gmail.com'
     endif

   
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
        UU=Hamk_bulk
        call eigensystem_c( 'V', 'U', Num_wann, UU, W)
  
        !> get velocity operator in Wannier basis
        !> \partial_k H_nm
        call dHdk_atomicgauge(k, Vmn_wann)

        !> spin axis igamma= x, y, z
        do igamma= 1, 3
           !> set Pauli matrix
           spin_sigma= pauli_matrices(:, :, igamma)
   
           !> calculate spin current operator j_spin_gamma_alpha^l_alpha= 1/2*{Sigma_l, v_alpha} 
           ! in order to calculate SHC^z_xy
           j_spin_gamma_alpha= 0d0
           do ialpha= 1, 3
              !> in Wannier basis
              call mat_mul(Num_wann, spin_sigma, Vmn_wann(:, :, ialpha), j_spin_gamma_alpha(:, :))
              call mat_mul(Num_wann, Vmn_wann(:, :, ialpha), spin_sigma, mat_t)
              j_spin_gamma_alpha(:, :)= j_spin_gamma_alpha(:, :)+ mat_t
              j_spin_gamma_alpha= j_spin_gamma_alpha/2d0
         
              !> rotate to Hamiltonian basis
              mat_t= j_spin_gamma_alpha(:, :)
              call rotation_to_Ham_basis(UU, mat_t, j_spin_gamma_alpha(:, :))
              call rotation_to_Ham_basis(UU, Vmn_wann(:, :, ialpha), Vmn_Ham(:, :, ialpha))

              do ieta= 1, NumberofEta
                 eta_local = eta_array(ieta)
                 !> \Omega_spin^l_n^{\gamma}(k)=-2\sum_{m}*aimag(Im({js(\gamma),v(\alpha)}/2)_nm*v_beta_mn))/((w(n)-w(m))^2+eta_arc^2)
                 do ibeta= 1, 3
                    Omega_spin= 0d0
                    do n= 1, Num_wann
                       do m= 1, Num_wann
                          if (m==n) cycle
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
             'xx^x', 'xy^x', 'xz^x', 'yx^x', 'yy^x', 'yz^x', 'zx^x', 'yy^x', 'zz^x', &
             'xx^y', 'xy^y', 'xz^y', 'yx^y', 'yy^y', 'yz^y', 'zx^y', 'yy^y', 'zz^y', &
             'xx^z', 'xy^z', 'xz^z', 'yx^z', 'yy^z', 'yz^z', 'zx^z', 'yy^z', 'zz^z'
           do ie=1, OmegaNum
              write(outfileindex, '(E16.8)', advance='no')energy(ie)/eV2Hartree
              do igamma=1, 3
              do ialpha=1, 3
              do ibeta=1, 3
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
        write(etaname, '(f12.2)')Eta_Arc*1000d0/eV2Hartree
        write(shcfilename, '(7a)')'sigma_shc_eta', trim(adjustl(etaname)), 'meV.txt'
        open(unit=outfileindex, file='sigma_shc.gnu')
        write(outfileindex, '(a)') 'set terminal pdf enhanced color font ",20"'
        write(outfileindex, '(a)')"set output 'sigma_shc.pdf'"
        write(outfileindex, '(a)')'set key samplen 0.8'
        write(outfileindex, '(a)')'set ylabel offset 0.0,0'
        write(outfileindex, '(a, f10.5, a, f10.5, a)')'set xrange [', OmegaMin/eV2Hartree, ':', OmegaMax/eV2Hartree, ']'
        write(outfileindex, '(a)')'set xlabel "Energy (eV)"'
        write(outfileindex, '(a)')'set ylabel "SHC (\hbar/e)S/cm"'
        write(outfileindex, '(5a)')"plot '",  trim(adjustl(shcfilename)),  "' u 1:21 w l title '\sigma_{xy}^z' lc rgb 'red' lw 4, \"
        write(outfileindex, '(5a)')"'",  trim(adjustl(shcfilename)), "' u 1:17 w l title '\sigma_{zx}^y' lc rgb 'blue' lw 4, \"
        write(outfileindex, '(5a)')"'", trim(adjustl(shcfilename)), "' u 1:13 w l title '\sigma_{xz}^y' lc rgb 'orange' lw 4 "
        close(outfileindex)
     endif



     deallocate( W, Hamk_bulk, UU, energy)
     deallocate( sigma_tensor_shc, sigma_tensor_shc_mpi)
 
     return
  end subroutine sigma_SHC
