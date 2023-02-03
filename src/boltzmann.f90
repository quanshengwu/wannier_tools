!>-------------------------------------------------------------------------------<!
!> this subroutine is used for anomalous Nernst calculation 
!> Hanqi Pi (hqpi1999@gmail.com)
!> Dec 14 2022 @ iopcas
!>-------------------------------------------------------------------------------<!

subroutine ANE_varyT
    ! PRL 97, 026603 (2006) eq(8)
    use para 
    implicit none 
    
    ! energy refers to the variable to be integrated
    ! mu refers to the chemical potential
    real(dp) :: energy,  omega, T
    
    real(dp), allocatable :: T_list(:)
    real(dp), allocatable :: energy_list(:)
    real(dp), allocatable :: nernst(:,:) 
    real(dp), allocatable :: ahc(:,:) 
    integer :: ie, iT

    ! integral region is (-ANE_int_interval, ANE_int_interval)
    ! num_step = 2*ANE_int_interval/ANE_int_step
    ! Note: make sure integrated region includes Efermi
    integer :: num_step 

    ! derivative of Fermi-Dirac function on energy
    real(dp), external :: partial_fermi

    num_step = (2*ANE_int_interval)/ANE_int_step

    allocate( T_list(NumT) )
    allocate( energy_list(num_step+1) )
    allocate( nernst(3, NumT)) 
    allocate( ahc(3, num_step+1)) 
    T_list = 0.0d0
    energy_list = 0.0d0
    nernst = 0.0d0
    ahc = 0.0d0

    !> energy to be integrated
    do ie = 1,  num_step+1
        energy_list(ie) = -ANE_int_interval+ANE_int_step*(ie-1)
    enddo

    !> temperature
    do iT=1, NumT
        if (NumT>1) then
           T_list(iT)= Tmin+(Tmax-Tmin)*(iT-1d0)/dble(NumT-1)
        else
           T_list= Tmin
        endif
    enddo ! iT

    !> get the anomalous hall conductivity
    call ahc_zerotmp(num_step, energy_list, ahc)

    do iT = 1, NumT
        T = T_list(iT)
        do ie = 1, num_step+1
            energy = energy_list(ie)
            nernst(:, iT) = nernst(:, iT) + partial_fermi(energy, T)*ahc(:,ie)&
                            *(energy/eV2Hartree)/T*(ANE_int_step/eV2Hartree)*100
        enddo
    enddo

    outfileindex = outfileindex+1
    if(cpuid .eq. 0) then
        open(unit=outfileindex, file='ane.txt')
        write(outfileindex, '("#",a)')'Anomalous Nernst coefficient (A(m-1 K-1))'
        write(outfileindex,'("#",a13, 3a16)')'Temperature', '\alpha_xy', '\alpha_yz', '\alpha_zx'
        do iT = 1, NumT
            write(outfileindex,'(F9.3,3E16.8)') T_list(iT), nernst(3, iT),&
                                                            nernst(1, iT),&
                                                            nernst(2, iT)
        
        enddo
        close(outfileindex)
    endif
    deallocate(T_list, energy_list, nernst, ahc)
    return


end subroutine ANe_varyT

subroutine ANE_varyEf
    ! PRL 97, 026603 (2006) eq(8)
    use para 
    implicit none 
    
    ! energy refers to the variable to be integrated
    ! mu refers to the chemical potential
    real(dp) :: energy,  omega, T, Ef
    
    ! real(dp), allocatable :: T_list(:)
    real(dp), allocatable :: Ef_list(:)
    real(dp), allocatable :: energy_list(:)
    real(dp), allocatable :: nernst(:,:) 
    real(dp), allocatable :: ahc(:,:) 
    integer :: ie, iT, ief

    ! integral region is (OmegaMin-ANE_int_interval,OmegaMax+ANE_int_interval)
    ! num_step = 2*ANE_int_interval/ANE_int_step
    ! Note: make sure integral region includes Efermi
    integer :: num_step 

    ! derivative of Fermi-Dirac function on energy
    real(dp), external :: partial_fermi

    num_step = (OmegaMax-OmegaMin+2*ANE_int_interval)/ANE_int_step

    ! allocate( T_list(NumT) )
    allocate( Ef_list(OmegaNum) )
    allocate( energy_list(num_step+1) )
    allocate( nernst(3, OmegaNum)) 
    allocate( ahc(3, num_step+1)) 
    ! T_list = 0.0d0
    energy_list = 0.0d0
    nernst = 0.0d0
    ahc = 0.0d0

    !> energy to be integrated
    do ie = 1,  num_step+1
        energy_list(ie) = (OmegaMin-ANE_int_interval)+ANE_int_step*(ie-1)
    enddo

     !> energy
    do ie=1, OmegaNum
        if (OmegaNum>1) then
           Ef_list(ie)= OmegaMin+ (OmegaMax-OmegaMin)* (ie-1d0)/dble(OmegaNum-1)
        else
           Ef_list= OmegaMin
        endif
    enddo ! ie

    !> get the anomalous hall conductivity
    call ahc_zerotmp(num_step, energy_list, ahc)

    do ief = 1, OmegaNum
        Ef = Ef_list(ief)
        do ie = 1, num_step+1
            energy = energy_list(ie)
            omega = energy - Ef
            nernst(:, ief) = nernst(:, ief)+ partial_fermi(omega, Tmin)*ahc(:,ie)&
                             *(omega/eV2Hartree)/Tmin*(ANE_int_step/eV2Hartree)*100
        enddo
    enddo

    outfileindex = outfileindex+1
    if(cpuid .eq. 0) then
        open(unit=outfileindex, file='ane.txt')
        write(outfileindex, '("#",a)')'Anomalous Nernst coefficient  (A(m-1 K-1))'
        write(outfileindex,'("#",a13, 3a16)')'Energy', '\alpha_xy', '\alpha_yz', '\alpha_zx'
        do ief = 1, OmegaNum
            write(outfileindex,'(F9.3,3E16.8)') Ef_list(ief)/eV2Hartree, nernst(3, ief),&
                                                            nernst(1, ief),&
                                                            nernst(2, ief)
        
        enddo
        close(outfileindex)
    endif
    deallocate(Ef_list, energy_list, nernst, ahc)
    return


end subroutine ANE_varyEf

subroutine ahc_zerotmp(num_step, energy_list, ahc)
    use para
    use wmpi 
    implicit none

    integer, intent(in) :: num_step
    real(dp), intent(in) :: energy_list(num_step+1)
    real(dp), intent(inout) :: ahc(3, num_step+1)

    real(dp), external :: fermi,partial_fermi
    
    real(dp) :: time_start, time_end
    real(dp) :: energy
    real(dp) :: ahc_mpi(3, num_step+1) 
    real(dp) :: k(3)
    real(dp) :: kdotr,Beta_fake
    ! eigen value/vector of H
    real(dp), allocatable :: W(:)
    complex(dp), allocatable :: Hamk_bulk(:,:)
    complex(dp), allocatable :: UU(:,:)
    complex(dp), allocatable :: UU_dag(:,:)
    complex(dp), allocatable :: Amat(:,:)

    ! velocities
    complex(dp), allocatable :: vx(:,:), vy(:,:), vz(:,:)

    ! berry curvature
    ! real(dp), allocatable :: omega_x(:), omega_y(:), omega_z(:)
    ! real(dp), allocatable :: Omega_x_t(:), Omega_y_t(:), Omega_z_t(:)
    complex(dp), allocatable :: omega_x(:), omega_y(:), omega_z(:)
    complex(dp), allocatable :: Omega_x_t(:), Omega_y_t(:), Omega_z_t(:)

    
    integer :: ik, knum, ikx, iky, ikz, iR, m, n, ierr, ie 
    
    ! initialize ahc to 0.0d0 before every iteration
    ahc = 0.0d0
    ahc_mpi = 0.0d0

    allocate(W(Num_wann))
    allocate(Hamk_bulk(Num_wann, Num_wann))
    allocate(UU(Num_wann, Num_wann))
    allocate(UU_dag(Num_wann, Num_wann))
    allocate(Amat(Num_wann, Num_wann))
    allocate(vx(Num_wann, num_wann))
    allocate(vy(Num_wann, num_wann))
    allocate(vz(Num_wann, num_wann))
    allocate(omega_x(Num_wann), omega_y(Num_wann), omega_z(Num_wann))
    allocate(omega_x_t(Num_wann), omega_y_t(Num_wann), omega_z_t(Num_wann))
    
    W = 0.0d0
    Hamk_bulk = 0.0d0
    UU = 0.0d0
    UU_dag = 0.0d0
    Amat = 0.0d0
    omega_x_t=0.0d0
    omega_y_t=0.0d0
    omega_z_t=0.0d0
    

    knum = Nk1*Nk2*Nk3

    call now(time_start)

    do ik = 1+cpuid, knum, num_cpu 
        if (cpuid.eq.0.and. mod(ik/num_cpu, 100).eq.0) then
            call now(time_end) 
            write(stdout, '(a, i18, "/", i18, a, f10.2, "s")') 'ik/knv3', &
            ik, knum, '  time left', (knum-ik)*(time_end-time_start)/num_cpu/100d0
            time_start= time_end
        endif

        ikx = (ik-1)/(Nk2*Nk3)+1
        iky = ((ik-1)-(ikx-1)*Nk2*Nk3)/Nk3+1        
        ikz = ik-(ikx-1)*Nk2*Nk3-(iky-1)*Nk3 
        k = K3D_start_cube+ K3D_vec1_cube*(ikx-1)/dble(Nk1) &
            + K3D_vec2_cube*(iky-1)/dble(Nk2) &
            + K3D_vec3_cube*(ikz-1)/dble(Nk3)
        
        
        ! get Hamilton
        call ham_bulk_latticegauge(k, Hamk_bulk)
        UU=Hamk_bulk
        call eigensystem_c('V','U', Num_wann, UU, W)
        UU_dag = conjg(transpose(UU))
        vx=0.0d0
        vy=0.0d0
        vz=0.0d0
        do iR = 1, Nrpts
            kdotr = dot_product(k,irvec(:,iR))
            vx = vx+zi*crvec(1,iR)*HmnR(:,:,iR)*exp(pi2zi*kdotr)/ndegen(iR)
            vy = vy+zi*crvec(2,iR)*HmnR(:,:,iR)*exp(pi2zi*kdotr)/ndegen(iR)
            vz = vz+zi*crvec(3,iR)*HmnR(:,:,iR)*exp(pi2zi*kdotr)/ndegen(iR)
        enddo
        !> unitility rotate velocity
        call mat_mul(Num_wann, vx, UU, Amat) 
        call mat_mul(Num_wann, UU_dag, Amat, vx) 
        call mat_mul(Num_wann, vy, UU, Amat) 
        call mat_mul(Num_wann, UU_dag, Amat, vy) 
        call mat_mul(Num_wann, vz, UU, Amat) 
        call mat_mul(Num_wann, UU_dag, Amat, vz)

        omega_x=0.0d0
        omega_y=0.0d0
        omega_z=0.0d0
        
        Beta_fake= 1d0/Eta_Arc

        !> work 
        do m= 1, Num_wann
            do n= 1, Num_wann
               if (abs(W(m)-W(n))<eps9) cycle
               Omega_x(m)= Omega_x(m)+ vy(n, m)*vz(m, n)/((W(m)-W(n))**2)
               Omega_y(m)= Omega_y(m)+ vz(n, m)*vx(m, n)/((W(m)-W(n))**2)
               Omega_z(m)= Omega_z(m)+ vx(n, m)*vy(m, n)/((W(m)-W(n))**2)
            enddo ! m
        enddo ! n
    
        Omega_x= -Omega_x*2d0*zi
        Omega_y= -Omega_y*2d0*zi
        Omega_z= -Omega_z*2d0*zi
        
        do ie=1, num_step+1
            energy = energy_list(ie)
            do m= 1, Num_wann
               Omega_x_t(m)= Omega_x(m)*fermi(W(m)-energy, Beta_fake)
               Omega_y_t(m)= Omega_y(m)*fermi(W(m)-energy, Beta_fake)
               Omega_z_t(m)= Omega_z(m)*fermi(W(m)-energy, Beta_fake)
            enddo
            ahc_mpi(1, ie)= ahc_mpi(1, ie)+ real(sum(Omega_x_t))
            ahc_mpi(2, ie)= ahc_mpi(2, ie)+ real(sum(Omega_y_t))
            ahc_mpi(3, ie)= ahc_mpi(3, ie)+ real(sum(Omega_z_t))
        enddo ! ie


    enddo

#if defined (MPI)
        call mpi_allreduce(ahc_mpi, ahc, size(ahc),&
                            mpi_dp, mpi_sum, mpi_cmw, ierr)
#else
        ahc = ahc_mpi
#endif

   ! --------------------------------------------------------------------
   ! At this point ahc contains
   !
   !           sum_k Omega =  N*V_c*int dk/[8(pi)^3] Omega
   !
   ! (N is the number of kpoints, V_c is the cell volume). We want
   !
   !           ahc = -(e^2/hbar) int dk/[4(pi)^3] Omega
   !
   ! Hence we need to multiply by  -e^2/hbar/N/V_c = 0.00024341/N/V_c
   ! where 'V_c' in Angstroms^3, and 'e', 'hbar' in SI units.
   ! Then it has units of 
   !
   !           C^2 (J*s)^-1 Angstrom^-3  Angstrom^2 = S/Angstrom                         
   ! 
   ! Thus, to have a result in S/cm, we need to multiply ahc with 24341/N/V_c
   ! --------------------------------------------------------------------


    !> in the version older than 2.5.2, we use Angstrom as the unit as length
    !ahc = -ahc/dble(knum)/Origin_cell%CellVolume*24341d0*kCubeVolume/Origin_cell%ReciprocalCellVolume ! in (Omega*cm)^-1
     
     !> in the latest version, we use the atomic unit
     !> in unit of (Ohm*m)^-1
    ahc= -ahc/dble(knum)/Origin_cell%CellVolume*Echarge**2/hbar/Bohr_radius*kCubeVolume/Origin_cell%ReciprocalCellVolume
     !> in unit of (Ohm*cm)^-1
    ahc= ahc/100d0

    ! > this part is for checking 
    outfileindex = outfileindex+1
    if(cpuid .eq. 0) then
        open(unit=outfileindex, file='ahc.txt')
        write(outfileindex, '("#",a)')'Anomalous Hall effect'
        write(outfileindex,'("#",a13, 3a16)')'energy', '\sigma_xy', '\sigma_yz', '\sigma_zx'
        do ie = 1, num_step+1
            energy = energy_list(ie)
            write(outfileindex,'(F9.3,3E16.8)') energy/eV2Hartree, ahc(3, ie),&
                                                        ahc(1, ie), &
                                                        ahc(2, ie)
        
        enddo
        close(outfileindex)
    endif
    deallocate(W, Hamk_bulk, UU, UU_dag, Amat, vx, vy, vz, omega_x, omega_y,omega_z)
    return
end subroutine ahc_zerotmp


function partial_fermi(omega, T) result(value)
    !this function sets the derivative of Fermi-Dirac function of energy
    !note that: omega = energy - chmical_potential

    use para
    implicit none

    real(dp), intent(in) :: omega, T 
    real(dp) :: value

    real(dp), parameter :: k_b = 8.617333262E-5!Boltzmann constant
    real(dp) :: beta_tmp

    beta_tmp = 1/(k_b*T)  

    if (abs(beta_tmp*omega/eV2Hartree) .gt. 37D0)  then
        value = 0.0d0
    else  
        value = -beta_tmp/((1.0d0+exp(omega/eV2Hartree*beta_tmp))*(1.0d0+exp(-omega/eV2Hartree*beta_tmp)))
    endif 

    return 
end function partial_fermi

