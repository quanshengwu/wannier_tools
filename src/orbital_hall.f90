subroutine Intra_Orbital_hall_conductivity
    use wmpi
    use para
    implicit none

    integer :: knv3, ik, ii, jj, mm, nn, ss, tt, ierr, imu, i, j
    integer :: ikx, iky, ikz
    real(dp) :: yita
    complex(dp) :: cmplx_i, cmplx_1, cmplx_0, omega
    real(dp) :: k(3)
    real(dp) :: time_start, time_end
    real(dp) :: Fermi_energy

    real(dp), allocatable :: W(:)
    complex(dp), allocatable :: Hamk_bulk(:, :)
    complex(dp), allocatable :: UU(:, :)

    complex(dp), allocatable :: D_Ham(:, :, :)
    complex(dp), allocatable :: V_Ham(:, :, :)
    
    complex(dp), allocatable :: L_wann_z(:, :)
    complex(dp), allocatable :: L_Ham_z(:, :)
    
    complex(dp), allocatable :: L_wann_x(:, :)
    complex(dp), allocatable :: L_Ham_x(:, :)

    complex(dp), allocatable :: L_wann_y(:, :)
    complex(dp), allocatable :: L_Ham_y(:, :)

    complex(dp), allocatable :: Orbital_Current_Lz(:, :, :)
    complex(dp), allocatable :: Orbital_Current_Lx(:, :, :)
    complex(dp), allocatable :: Orbital_Current_Ly(:, :, :)
    
    complex(dp), allocatable :: Orbital_HC_Lz(:, :, :)
    complex(dp), allocatable :: Orbital_HC_Lz_mpi(:, :, :)
    complex(dp), allocatable :: Orbital_HC_Lx(:, :, :)
    complex(dp), allocatable :: Orbital_HC_Lx_mpi(:, :, :)
    complex(dp), allocatable :: Orbital_HC_Ly(:, :, :)
    complex(dp), allocatable :: Orbital_HC_Ly_mpi(:, :, :)
    

    real(dp), allocatable :: mulist(:)
    real(dp), allocatable :: occ(:)
    real(dp) :: OHC_out(54)

    knv3= Nk1*Nk2*Nk3
    cmplx_i = (0.0d0, 1.0d0)
    cmplx_1 = (1.0d0, 0.0d0)
    cmplx_0 = (0.0d0, 0.0d0)


    allocate( W(Num_wann) )
    allocate(Hamk_bulk(Num_wann, Num_wann))
    allocate( UU(Num_wann, Num_wann))


    allocate(mulist(OmegaNum))
    allocate(D_Ham(Num_wann, Num_wann, 3))
    allocate(V_Ham(Num_wann, Num_wann, 3))
    allocate(L_Ham_z(Num_wann, Num_wann),L_wann_z(Num_wann, Num_wann))
    allocate(L_Ham_x(Num_wann, Num_wann),L_wann_x(Num_wann, Num_wann))
    allocate(L_Ham_y(Num_wann, Num_wann),L_wann_y(Num_wann, Num_wann))
    allocate(Orbital_HC_Lz(OmegaNum, 3, 3))
    allocate(Orbital_HC_Lz_mpi(OmegaNum, 3, 3))
    allocate(Orbital_HC_Lx(OmegaNum, 3, 3))
    allocate(Orbital_HC_Lx_mpi(OmegaNum, 3, 3))
    allocate(Orbital_HC_Ly(OmegaNum, 3, 3))
    allocate(Orbital_HC_Ly_mpi(OmegaNum, 3, 3))
    allocate(Orbital_Current_Lz(Num_wann, Num_wann, 3))
    allocate(Orbital_Current_Lx(Num_wann, Num_wann, 3))
    allocate(Orbital_Current_Ly(Num_wann, Num_wann, 3))
    allocate( occ(Num_wann))


    yita = 0.002
    Orbital_HC_Lz = cmplx_0
    Orbital_HC_Lz_mpi = cmplx_0
    
    Orbital_HC_Lx = cmplx_0
    Orbital_HC_Lx_mpi = cmplx_0
    Orbital_HC_Ly = cmplx_0
    Orbital_HC_Ly_mpi = cmplx_0
    L_wann_z = cmplx_0
    L_wann_x = cmplx_0
    L_wann_y = cmplx_0
    
    !LZ
    !L_wann_z(3,3) = 1
    !L_wann_z(4,4) = 1
    !L_wann_z(7,7) = 1
    !L_wann_z(8,8) = 1
    
    !L_wann_z(3,4) = cmplx_i
    !L_wann_z(4,3) = -cmplx_i
    !L_wann_z(6,7) = cmplx_i
    !L_wann_z(7,6) = -cmplx_i

    

    !Lx
    !L_wann_x(4,2) = cmplx_i
    !L_wann_x(2,4) = -cmplx_i
    !L_wann_x(7,5) = cmplx_i
    !L_wann_x(5,7) = -cmplx_i



    !Ly
    !L_wann_y(2,3) = cmplx_i
    !L_wann_y(3,2) = -cmplx_i
    !L_wann_y(5,6) = cmplx_i
    !L_wann_y(6,5) = -cmplx_i

    call read_OAM_operator(L_wann_x, L_wann_y, L_wann_z)


    

    


    




    

    do imu= 1, OmegaNum 
        mulist(imu) = OmegaMin + ( OmegaMax - OmegaMin) * imu / OmegaNum
    enddo

    write(stdout,*) 'kmesh check start!'
    do ik = 1, knv3
        ikx= (ik-1)/(nk2*nk3)+1
        iky= ((ik-1-(ikx-1)*Nk2*Nk3)/nk3)+1
        ikz= (ik-(iky-1)*Nk3- (ikx-1)*Nk2*Nk3)
        !write(stdout,*) ikx, iky, ikz
    enddo
    write(stdout,*) 'kmesh check end!'

    call now(time_start) 
    !do ik= 1, knv3
    do ik= 1+ cpuid, knv3, num_cpu
        if (cpuid.eq.0.and. mod(ik/num_cpu, 100).eq.0) then
            call now(time_end) 
            write(stdout, '(a, i18, "/", i18, a, f10.2, "s")') 'ik/knv3', &
            ik, knv3, '  time left', (knv3-ik)*(time_end-time_start)/num_cpu/100d0
            time_start= time_end
        endif

        Hamk_bulk = cmplx_0
        V_Ham = cmplx_0
        D_Ham = cmplx_0
        L_Ham_z = cmplx_0
        L_Ham_x = cmplx_0
        L_Ham_y = cmplx_0
        
        UU = cmplx_0
        W = 0.0_dp
        Orbital_Current_Lz = cmplx_0 
        Orbital_Current_Lx = cmplx_0 
        Orbital_Current_Ly = cmplx_0 

        


        ikx= (ik-1)/(nk2*nk3)+1
        iky= ((ik-1-(ikx-1)*Nk2*Nk3)/nk3)+1
        ikz= (ik-(iky-1)*Nk3- (ikx-1)*Nk2*Nk3)
        k = K3D_start_cube+ K3D_vec1_cube*(ikx-1)/dble(nk1)  &
            + K3D_vec2_cube*(iky-1)/dble(nk2)  &
            + K3D_vec3_cube*(ikz-1)/dble(nk3)
        
        call ham_bulk_atomicgauge(k, Hamk_bulk)

        UU=Hamk_bulk
        call eigensystem_c( 'V', 'U', Num_wann, UU, W)

        call dHdk_atomicgauge_Ham(k, UU, V_Ham)
        !call get_Dmn_Ham(W, V_Ham, D_Ham)

        !call get_Lmn(W, D_Ham, V_Ham, L_Ham_z)

        call rotation_to_Ham_basis(UU, L_wann_z, L_Ham_z)
        call rotation_to_Ham_basis(UU, L_wann_x, L_Ham_x)
        call rotation_to_Ham_basis(UU, L_wann_y, L_Ham_y)
        !L_Ham_z(:,:)=L_Wann_z(:,:)
        do ss=1,3
            do mm = 1,Num_wann
                do nn = 1,Num_wann
                    do tt = 1,Num_wann
                    
                        Orbital_Current_Lz(mm,nn,ss) = Orbital_Current_Lz(mm,nn,ss) + &
                                0.5*( V_Ham(mm,tt,ss) * L_Ham_z(tt,nn) + L_Ham_z(mm,tt) * V_Ham(tt,nn,ss) )
                        
                        Orbital_Current_Lx(mm,nn,ss) = Orbital_Current_Lx(mm,nn,ss) + &
                                0.5*( V_Ham(mm,tt,ss) * L_Ham_x(tt,nn) + L_Ham_x(mm,tt) * V_Ham(tt,nn,ss) )

                        Orbital_Current_Ly(mm,nn,ss) = Orbital_Current_Ly(mm,nn,ss) + &
                                0.5*( V_Ham(mm,tt,ss) * L_Ham_y(tt,nn) + L_Ham_y(mm,tt) * V_Ham(tt,nn,ss) )
                    enddo
                enddo
            enddo
        enddo

       !Orbital_Current_Lz(:,:,:) = V_Ham(:,:,:) 
        

        do imu = 1, OmegaNum    
            Fermi_energy = mulist(imu)
            !write(stdout,*) 'Fermi_energy', Fermi_energy 

            occ = 0.0_dp
            do i = 1, Num_wann
                if (W(i) < Fermi_energy) occ(i) = 1
                !
            enddo
            !write(stdout,*) 'occ', occ(:) 
            do ii = 1, Num_wann
                do jj = 1, Num_wann
                    if(ii==jj) cycle
                    !if(ii>jj) cycle
                    !if( abs(W(ii)-W(jj))<eps12*1) cycle
                    do mm = 1, 3
                        do nn = 1, 3
                            Orbital_HC_Lz_mpi(imu,mm,nn) = Orbital_HC_Lz_mpi(imu,mm,nn) + &
                                -cmplx_i*(occ(ii)-occ(jj))*V_Ham(ii,jj,mm)*Orbital_Current_Lz(jj,ii,nn)/(W(ii)-W(jj)+cmplx_i*yita)**2

                            !Orbital_HC_Lz_mpi(imu,1) = Orbital_HC_Lz_mpi(imu,1) + &
                            !    -cmplx_i*(occ(ii)-occ(jj))*V_Ham(ii,jj,2)*Orbital_Current_Lz(jj,ii,3)/(W(ii)-W(jj)+cmplx_i*yita)**2

                            !Orbital_HC_Lz_mpi(imu,2) = Orbital_HC_Lz_mpi(imu,2) + &
                            !    -cmplx_i*(occ(ii)-occ(jj))*V_Ham(ii,jj,3)*Orbital_Current_Lz(jj,ii,1)/(W(ii)-W(jj)+cmplx_i*yita)**2
                    

                    
                            Orbital_HC_Lx_mpi(imu,mm,nn) = Orbital_HC_Lx_mpi(imu,mm,nn) + &
                                -cmplx_i*(occ(ii)-occ(jj))*V_Ham(ii,jj,mm)*Orbital_Current_Lx(jj,ii,nn)/(W(ii)-W(jj)+cmplx_i*yita)**2

                            !Orbital_HC_Lx_mpi(imu,1) = Orbital_HC_Lx_mpi(imu,1) + &
                            !    -cmplx_i*(occ(ii)-occ(jj))*V_Ham(ii,jj,2)*Orbital_Current_Lx(jj,ii,3)/(W(ii)-W(jj)+cmplx_i*yita)**2

                            !Orbital_HC_Lx_mpi(imu,2) = Orbital_HC_Lx_mpi(imu,2) + &
                            !    -cmplx_i*(occ(ii)-occ(jj))*V_Ham(ii,jj,3)*Orbital_Current_Lx(jj,ii,1)/(W(ii)-W(jj)+cmplx_i*yita)**2

                            
                            
                            Orbital_HC_Ly_mpi(imu,mm,nn) = Orbital_HC_Ly_mpi(imu,mm,nn) + &
                                -cmplx_i*(occ(ii)-occ(jj))*V_Ham(ii,jj,mm)*Orbital_Current_Ly(jj,ii,nn)/(W(ii)-W(jj)+cmplx_i*yita)**2

                            !Orbital_HC_Ly_mpi(imu,1) = Orbital_HC_Ly_mpi(imu,1) + &
                            !    -cmplx_i*(occ(ii)-occ(jj))*V_Ham(ii,jj,2)*Orbital_Current_Ly(jj,ii,3)/(W(ii)-W(jj)+cmplx_i*yita)**2

                            !Orbital_HC_Ly_mpi(imu,2) = Orbital_HC_Ly_mpi(imu,2) + &
                            !    -cmplx_i*(occ(ii)-occ(jj))*V_Ham(ii,jj,3)*Orbital_Current_Ly(jj,ii,1)/(W(ii)-W(jj)+cmplx_i*yita)**2
                        enddo
                    enddo
                enddo
            enddo
        
        enddo

    enddo

#if defined (MPI)
    call mpi_allreduce(Orbital_HC_Lz_mpi, Orbital_HC_Lz, size(Orbital_HC_Lz), mpi_dc,mpi_sum,mpi_cmw,ierr)
    call mpi_allreduce(Orbital_HC_Lx_mpi, Orbital_HC_Lx, size(Orbital_HC_Lx), mpi_dc,mpi_sum,mpi_cmw,ierr)
    call mpi_allreduce(Orbital_HC_Ly_mpi, Orbital_HC_Ly, size(Orbital_HC_Ly), mpi_dc,mpi_sum,mpi_cmw,ierr)
#else
    Orbital_HC_Lz = Orbital_HC_Lz_mpi
    Orbital_HC_Lx = Orbital_HC_Lx_mpi
    Orbital_HC_Ly = Orbital_HC_Ly_mpi
#endif

    !Orbital_HC_Lz = Orbital_HC_Lz_mpi


    if (cpuid.eq.0) then
        outfileindex= outfileindex+ 1
        open(unit=outfileindex, file='Intra_Orbital_Hall_Conductivity.dat')
        !write(outfileindex, '("#",10a)')' the symmetric part of linear conductivity'
        write(outfileindex, "('#column', i5, 3000i16)")(i, i=1, 55)
        write(outfileindex, '("#",a13, 20a16)') 'Omega', 'Re[\sigma_Lx_xx]', 'Im[\sigma_Lx_xx]', 'Re[\sigma_Lx_xy]', &
                            'Im[\sigma_Lx_xy]', 'Re[\sigma_Lx_xz]', 'Im[\sigma_Lx_xz]', &
                            'Re[\sigma_Lx_yx]', 'Im[\sigma_Lx_yx]', 'Re[\sigma_Lx_yy]', &
                            'Im[\sigma_Lx_yy]', 'Re[\sigma_Lx_yz]', 'Im[\sigma_Lx_yz]', &
                            'Re[\sigma_Lx_zx]', 'Im[\sigma_Lx_zx]', 'Re[\sigma_Lx_zy]', &
                            'Im[\sigma_Lx_zy]', 'Re[\sigma_Lx_zz]', 'Im[\sigma_Lx_zz]', &
                            'Re[\sigma_Ly_xx]', 'Im[\sigma_Ly_xx]', 'Re[\sigma_Ly_xy]', &
                            'Im[\sigma_Ly_xy]', 'Re[\sigma_Ly_xz]', 'Im[\sigma_Ly_xz]', &
                            'Re[\sigma_Ly_yx]', 'Im[\sigma_Ly_yx]', 'Re[\sigma_Ly_yy]', &
                            'Im[\sigma_Ly_yy]', 'Re[\sigma_Ly_yz]', 'Im[\sigma_Ly_yz]', &
                            'Re[\sigma_Ly_zx]', 'Im[\sigma_Ly_zx]', 'Re[\sigma_Ly_zy]', &
                            'Im[\sigma_Ly_zy]', 'Re[\sigma_Ly_zz]', 'Im[\sigma_Ly_zz]', &
                            'Re[\sigma_Lz_xx]', 'Im[\sigma_Lz_xx]', 'Re[\sigma_Lz_xy]', &
                            'Im[\sigma_Lz_xy]', 'Re[\sigma_Lz_xz]', 'Im[\sigma_Lz_xz]', &
                            'Re[\sigma_Lz_yx]', 'Im[\sigma_Lz_yx]', 'Re[\sigma_Lz_yy]', &
                            'Im[\sigma_Lz_yy]', 'Re[\sigma_Lz_yz]', 'Im[\sigma_Lz_yz]', &
                            'Re[\sigma_Lz_zx]', 'Im[\sigma_Lz_zx]', 'Re[\sigma_Lz_zy]', &
                            'Im[\sigma_Lz_zy]', 'Re[\sigma_Lz_zz]', 'Im[\sigma_Lz_zz]'
        do imu=1,OmegaNum
            

            OHC_out(1) = real(Orbital_HC_Lx(imu,1,1)) /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Origin_cell%cell_parameters(1)/Angstrom2atomic*1E-8!*Echarge**2/hbar
            OHC_out(2) = aimag(Orbital_HC_Lx(imu,1,1)) /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Origin_cell%cell_parameters(1)/Angstrom2atomic*1E-8!*Echarge**2/hbar
            OHC_out(3) = real(Orbital_HC_Lx(imu,1,2)) /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Origin_cell%cell_parameters(1)/Angstrom2atomic*1E-8!*Echarge**2/hbar
            OHC_out(4) = aimag(Orbital_HC_Lx(imu,1,2)) /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Origin_cell%cell_parameters(1)/Angstrom2atomic*1E-8!*Echarge**2/hbar
            OHC_out(5) = real(Orbital_HC_Lx(imu,1,3)) /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Origin_cell%cell_parameters(1)/Angstrom2atomic*1E-8!*Echarge**2/hbar
            OHC_out(6) = aimag(Orbital_HC_Lx(imu,1,3))/dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Origin_cell%cell_parameters(1)/Angstrom2atomic*1E-8!*Echarge**2/hbar
            OHC_out(7) = real(Orbital_HC_Lx(imu,2,1)) /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Origin_cell%cell_parameters(1)/Angstrom2atomic*1E-8!*Echarge**2/hbar
            OHC_out(8) = aimag(Orbital_HC_Lx(imu,2,1)) /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Origin_cell%cell_parameters(1)/Angstrom2atomic*1E-8!*Echarge**2/hbar
            OHC_out(9) = real(Orbital_HC_Lx(imu,2,2)) /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Origin_cell%cell_parameters(1)/Angstrom2atomic*1E-8!*Echarge**2/hbar
            OHC_out(10) = aimag(Orbital_HC_Lx(imu,2,2)) /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Origin_cell%cell_parameters(1)/Angstrom2atomic*1E-8!*Echarge**2/hbar
            OHC_out(11) = real(Orbital_HC_Lx(imu,2,3)) /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Origin_cell%cell_parameters(1)/Angstrom2atomic*1E-8!*Echarge**2/hbar
            OHC_out(12) = aimag(Orbital_HC_Lx(imu,2,3))/dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Origin_cell%cell_parameters(1)/Angstrom2atomic*1E-8!*Echarge**2/hbar
            OHC_out(13) = real(Orbital_HC_Lx(imu,3,1)) /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Origin_cell%cell_parameters(1)/Angstrom2atomic*1E-8!*Echarge**2/hbar
            OHC_out(14) = aimag(Orbital_HC_Lx(imu,3,1)) /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Origin_cell%cell_parameters(1)/Angstrom2atomic*1E-8!*Echarge**2/hbar
            OHC_out(15) = real(Orbital_HC_Lx(imu,3,2)) /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Origin_cell%cell_parameters(1)/Angstrom2atomic*1E-8!*Echarge**2/hbar
            OHC_out(16) = aimag(Orbital_HC_Lx(imu,3,2)) /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Origin_cell%cell_parameters(1)/Angstrom2atomic*1E-8!*Echarge**2/hbar
            OHC_out(17) = real(Orbital_HC_Lx(imu,3,3)) /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Origin_cell%cell_parameters(1)/Angstrom2atomic*1E-8!*Echarge**2/hbar
            OHC_out(18) = aimag(Orbital_HC_Lx(imu,3,3))/dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Origin_cell%cell_parameters(1)/Angstrom2atomic*1E-8!*Echarge**2/hbar

            OHC_out(19) = real(Orbital_HC_Ly(imu,1,1)) /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Origin_cell%cell_parameters(2)/Angstrom2atomic*1E-8!*Echarge**2/hbar
            OHC_out(20) = aimag(Orbital_HC_Ly(imu,1,1)) /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Origin_cell%cell_parameters(2)/Angstrom2atomic*1E-8!*Echarge**2/hbar
            OHC_out(21) = real(Orbital_HC_Ly(imu,1,2)) /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Origin_cell%cell_parameters(2)/Angstrom2atomic*1E-8!*Echarge**2/hbar
            OHC_out(22) = aimag(Orbital_HC_Ly(imu,1,2)) /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Origin_cell%cell_parameters(2)/Angstrom2atomic*1E-8!*Echarge**2/hbar
            OHC_out(23) = real(Orbital_HC_Ly(imu,1,3)) /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Origin_cell%cell_parameters(2)/Angstrom2atomic*1E-8!*Echarge**2/hbar
            OHC_out(24) = aimag(Orbital_HC_Ly(imu,1,3))/dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Origin_cell%cell_parameters(2)/Angstrom2atomic*1E-8!*Echarge**2/hbar
            OHC_out(25) = real(Orbital_HC_Ly(imu,2,1)) /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Origin_cell%cell_parameters(2)/Angstrom2atomic*1E-8!*Echarge**2/hbar
            OHC_out(26) = aimag(Orbital_HC_Ly(imu,2,1)) /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Origin_cell%cell_parameters(2)/Angstrom2atomic*1E-8!*Echarge**2/hbar
            OHC_out(27) = real(Orbital_HC_Ly(imu,2,2)) /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Origin_cell%cell_parameters(2)/Angstrom2atomic*1E-8!*Echarge**2/hbar
            OHC_out(28) = aimag(Orbital_HC_Ly(imu,2,2)) /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Origin_cell%cell_parameters(2)/Angstrom2atomic*1E-8!*Echarge**2/hbar
            OHC_out(29) = real(Orbital_HC_Ly(imu,2,3)) /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Origin_cell%cell_parameters(2)/Angstrom2atomic*1E-8!*Echarge**2/hbar
            OHC_out(30) = aimag(Orbital_HC_Ly(imu,2,3))/dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Origin_cell%cell_parameters(2)/Angstrom2atomic*1E-8!*Echarge**2/hbar
            OHC_out(31) = real(Orbital_HC_Ly(imu,3,1)) /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Origin_cell%cell_parameters(2)/Angstrom2atomic*1E-8!*Echarge**2/hbar
            OHC_out(32) = aimag(Orbital_HC_Ly(imu,3,1)) /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Origin_cell%cell_parameters(2)/Angstrom2atomic*1E-8!*Echarge**2/hbar
            OHC_out(33) = real(Orbital_HC_Ly(imu,3,2)) /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Origin_cell%cell_parameters(2)/Angstrom2atomic*1E-8!*Echarge**2/hbar
            OHC_out(34) = aimag(Orbital_HC_Ly(imu,3,2)) /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Origin_cell%cell_parameters(2)/Angstrom2atomic*1E-8!*Echarge**2/hbar
            OHC_out(35) = real(Orbital_HC_Ly(imu,3,3)) /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Origin_cell%cell_parameters(2)/Angstrom2atomic*1E-8!*Echarge**2/hbar
            OHC_out(36) = aimag(Orbital_HC_Ly(imu,3,3))/dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Origin_cell%cell_parameters(2)/Angstrom2atomic*1E-8!*Echarge**2/hbar

            OHC_out(37) = real(Orbital_HC_Lz(imu,1,1)) /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Origin_cell%cell_parameters(3)/Angstrom2atomic*1E-8!*Echarge**2/hbar
            OHC_out(38) = aimag(Orbital_HC_Lz(imu,1,1)) /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Origin_cell%cell_parameters(3)/Angstrom2atomic*1E-8!*Echarge**2/hbar
            OHC_out(39) = real(Orbital_HC_Lz(imu,1,2)) /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Origin_cell%cell_parameters(3)/Angstrom2atomic*1E-8!*Echarge**2/hbar
            OHC_out(40) = aimag(Orbital_HC_Lz(imu,1,2)) /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Origin_cell%cell_parameters(3)/Angstrom2atomic*1E-8!*Echarge**2/hbar
            OHC_out(41) = real(Orbital_HC_Lz(imu,1,3)) /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Origin_cell%cell_parameters(3)/Angstrom2atomic*1E-8!*Echarge**2/hbar
            OHC_out(42) = aimag(Orbital_HC_Lz(imu,1,3))/dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Origin_cell%cell_parameters(3)/Angstrom2atomic*1E-8!*Echarge**2/hbar
            OHC_out(43) = real(Orbital_HC_Lz(imu,2,1)) /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Origin_cell%cell_parameters(3)/Angstrom2atomic*1E-8!*Echarge**2/hbar
            OHC_out(44) = aimag(Orbital_HC_Lz(imu,2,1)) /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Origin_cell%cell_parameters(3)/Angstrom2atomic*1E-8!*Echarge**2/hbar
            OHC_out(45) = real(Orbital_HC_Lz(imu,2,2)) /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Origin_cell%cell_parameters(3)/Angstrom2atomic*1E-8!*Echarge**2/hbar
            OHC_out(46) = aimag(Orbital_HC_Lz(imu,2,2)) /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Origin_cell%cell_parameters(3)/Angstrom2atomic*1E-8!*Echarge**2/hbar
            OHC_out(47) = real(Orbital_HC_Lz(imu,2,3)) /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Origin_cell%cell_parameters(3)/Angstrom2atomic*1E-8!*Echarge**2/hbar
            OHC_out(48) = aimag(Orbital_HC_Lz(imu,2,3))/dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100**Origin_cell%cell_parameters(3)/Angstrom2atomic*1E-8!*Echarge**2/hbar
            OHC_out(49) = real(Orbital_HC_Lz(imu,3,1)) /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Origin_cell%cell_parameters(3)/Angstrom2atomic*1E-8!*Echarge**2/hbar
            OHC_out(50) = aimag(Orbital_HC_Lz(imu,3,1)) /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Origin_cell%cell_parameters(3)/Angstrom2atomic*1E-8!*Echarge**2/hbar
            OHC_out(51) = real(Orbital_HC_Lz(imu,3,2)) /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Origin_cell%cell_parameters(3)/Angstrom2atomic*1E-8!*Echarge**2/hbar
            OHC_out(52) = aimag(Orbital_HC_Lz(imu,3,2)) /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Origin_cell%cell_parameters(3)/Angstrom2atomic*1E-8!*Echarge**2/hbar
            OHC_out(53) = real(Orbital_HC_Lz(imu,3,3)) /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Origin_cell%cell_parameters(3)/Angstrom2atomic*1E-8!*Echarge**2/hbar
            OHC_out(54) = aimag(Orbital_HC_Lz(imu,3,3))/dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Origin_cell%cell_parameters(3)/Angstrom2atomic*1E-8!*Echarge**2/hbar

            Fermi_energy = mulist(imu)/eV2Hartree
            write(outfileindex, '(200E16.8)')  Fermi_energy  , OHC_out
        enddo
        close(outfileindex)
    endif




end subroutine Intra_Orbital_hall_conductivity




















subroutine Intra_Orbital_hall_curvature_line
    use wmpi
    use para
    implicit none

    integer :: knv3, ik, ii, jj, mm, nn, ss, tt, ierr, imu, i, j
    integer :: ikx, iky, ikz
    real(dp) :: yita
    complex(dp) :: cmplx_i, cmplx_1, cmplx_0, omega
    real(dp) :: k(3)
    real(dp) :: time_start, time_end
    real(dp) :: Fermi_energy

    real(dp), allocatable :: W(:)
    complex(dp), allocatable :: Hamk_bulk(:, :)
    complex(dp), allocatable :: UU(:, :)

    complex(dp), allocatable :: D_Ham(:, :, :)
    complex(dp), allocatable :: V_Ham(:, :, :)
    
    complex(dp), allocatable :: L_wann_z(:, :)
    complex(dp), allocatable :: L_Ham_z(:, :)
    
    complex(dp), allocatable :: L_wann_x(:, :)
    complex(dp), allocatable :: L_Ham_x(:, :)

    complex(dp), allocatable :: L_wann_y(:, :)
    complex(dp), allocatable :: L_Ham_y(:, :)

    complex(dp), allocatable :: Orbital_Current_Lz(:, :, :)
    complex(dp), allocatable :: Orbital_Current_Lx(:, :, :)
    complex(dp), allocatable :: Orbital_Current_Ly(:, :, :)
    
    complex(dp), allocatable :: Orbital_HC_Lz(:, :)
    complex(dp), allocatable :: Orbital_HC_Lz_mpi(:, :)
    complex(dp), allocatable :: Orbital_HC_Lx(:, :)
    complex(dp), allocatable :: Orbital_HC_Lx_mpi(:, :)
    complex(dp), allocatable :: Orbital_HC_Ly(:, :)
    complex(dp), allocatable :: Orbital_HC_Ly_mpi(:, :)
    

    real(dp), allocatable :: mulist(:)
    real(dp), allocatable :: occ(:)
    real(dp) :: OHC_out(18)

    knv3= nk3_band
    cmplx_i = (0.0d0, 1.0d0)
    cmplx_1 = (1.0d0, 0.0d0)
    cmplx_0 = (0.0d0, 0.0d0)


    allocate( W(Num_wann) )
    allocate(Hamk_bulk(Num_wann, Num_wann))
    allocate( UU(Num_wann, Num_wann))


    allocate(mulist(OmegaNum))
    allocate(D_Ham(Num_wann, Num_wann, 3))
    allocate(V_Ham(Num_wann, Num_wann, 3))
    allocate(L_Ham_z(Num_wann, Num_wann),L_wann_z(Num_wann, Num_wann))
    allocate(L_Ham_x(Num_wann, Num_wann),L_wann_x(Num_wann, Num_wann))
    allocate(L_Ham_y(Num_wann, Num_wann),L_wann_y(Num_wann, Num_wann))

    allocate(Orbital_HC_Lz(knv3, 3))
    allocate(Orbital_HC_Lz_mpi(knv3, 3))
    allocate(Orbital_HC_Lx(knv3, 3))
    allocate(Orbital_HC_Lx_mpi(knv3, 3))
    allocate(Orbital_HC_Ly(knv3, 3))
    allocate(Orbital_HC_Ly_mpi(knv3, 3))
    
    allocate(Orbital_Current_Lz(Num_wann, Num_wann, 3))
    allocate(Orbital_Current_Lx(Num_wann, Num_wann, 3))
    allocate(Orbital_Current_Ly(Num_wann, Num_wann, 3))
    allocate( occ(Num_wann))


    yita = 0.001
    Orbital_HC_Lz = cmplx_0
    Orbital_HC_Lz_mpi = cmplx_0
    Orbital_HC_Lx = cmplx_0
    Orbital_HC_Lx_mpi = cmplx_0
    Orbital_HC_Ly = cmplx_0
    Orbital_HC_Ly_mpi = cmplx_0
    
    L_wann_z = cmplx_0
    L_wann_x = cmplx_0
    L_wann_y = cmplx_0
    
    L_wann_z(3,4) = cmplx_i
    L_wann_z(4,3) = -cmplx_i
    L_wann_z(6,7) = cmplx_i
    L_wann_z(7,6) = -cmplx_i

    

    !Lx
    L_wann_x(4,2) = cmplx_i
    L_wann_x(2,4) = -cmplx_i
    L_wann_x(7,5) = cmplx_i
    L_wann_x(5,7) = -cmplx_i



    !Ly
    L_wann_y(2,3) = cmplx_i
    L_wann_y(3,2) = -cmplx_i
    L_wann_y(5,6) = cmplx_i
    L_wann_y(6,5) = -cmplx_i

    

    


    




    

    do imu= 1, OmegaNum 
        mulist(imu) = OmegaMin + ( OmegaMax - OmegaMin) * imu / OmegaNum
    enddo



    call now(time_start) 
    !do ik= 1, knv3
    do ik= 1+ cpuid, knv3, num_cpu
        if (cpuid.eq.0.and. mod(ik/num_cpu, 100).eq.0) then
            call now(time_end) 
            write(stdout, '(a, i18, "/", i18, a, f10.2, "s")') 'ik/knv3', &
            ik, knv3, '  time left', (knv3-ik)*(time_end-time_start)/num_cpu/100d0
            time_start= time_end
        endif

        Hamk_bulk = cmplx_0
        V_Ham = cmplx_0
        D_Ham = cmplx_0
        L_Ham_z = cmplx_0
        L_Ham_x = cmplx_0
        L_Ham_y = cmplx_0
        
        UU = cmplx_0
        W = 0.0_dp
        Orbital_Current_Lz = cmplx_0 
        Orbital_Current_Lx = cmplx_0 
        Orbital_Current_Ly = cmplx_0 

        


        k = kpath_3d(:, ik)
        
        call ham_bulk_atomicgauge(k, Hamk_bulk)

        UU=Hamk_bulk
        call eigensystem_c( 'V', 'U', Num_wann, UU, W)

        call dHdk_atomicgauge_Ham(k, UU, V_Ham)
        !call get_Dmn_Ham(W, V_Ham, D_Ham)

        !call get_Lmn(W, D_Ham, V_Ham, L_Ham_z)

        call rotation_to_Ham_basis(UU, L_wann_z, L_Ham_z)
        call rotation_to_Ham_basis(UU, L_wann_x, L_Ham_x)
        call rotation_to_Ham_basis(UU, L_wann_y, L_Ham_y)
        !L_Ham_z(:,:)=L_Wann_z(:,:)
        do ss=1,3
            do mm = 1,Num_wann
                do nn = 1,Num_wann
                    do tt = 1,Num_wann
                    
                        Orbital_Current_Lz(mm,nn,ss) = Orbital_Current_Lz(mm,nn,ss) + &
                                0.5*( V_Ham(mm,tt,ss) * L_Ham_z(tt,nn) + L_Ham_z(mm,tt) * V_Ham(tt,nn,ss) )
                        
                        Orbital_Current_Lx(mm,nn,ss) = Orbital_Current_Lx(mm,nn,ss) + &
                                0.5*( V_Ham(mm,tt,ss) * L_Ham_x(tt,nn) + L_Ham_x(mm,tt) * V_Ham(tt,nn,ss) )

                        Orbital_Current_Ly(mm,nn,ss) = Orbital_Current_Ly(mm,nn,ss) + &
                                0.5*( V_Ham(mm,tt,ss) * L_Ham_y(tt,nn) + L_Ham_y(mm,tt) * V_Ham(tt,nn,ss) )
                    enddo
                enddo
            enddo
        enddo

       !Orbital_Current_Lz(:,:,:) = V_Ham(:,:,:) 
        

        
        Fermi_energy = 0
        !write(stdout,*) 'Fermi_energy', Fermi_energy 

        occ = 0.0_dp
        do i = 1, Num_wann
            if (W(i) < Fermi_energy) occ(i) = 1
            !
        enddo
        !write(stdout,*) 'occ', occ(:) 
        do ii = 1, Num_wann
            do jj = 1, Num_wann
                if(ii==jj) cycle
                !if(occ(ii) < occ(jj)) cycle
                !if(ii>jj) cycle
                !if( abs(W(ii)-W(jj))<eps12*1) cycle
                Orbital_HC_Lz_mpi(ik,3) = Orbital_HC_Lz_mpi(ik,3) + &
                    -cmplx_i*(occ(ii)-occ(jj))*V_Ham(ii,jj,1)*Orbital_Current_Lz(jj,ii,2)/(W(ii)-W(jj)+cmplx_i*yita)**2

                Orbital_HC_Lz_mpi(ik,1) = Orbital_HC_Lz_mpi(ik,1) + &
                    -cmplx_i*(occ(ii)-occ(jj))*V_Ham(ii,jj,2)*Orbital_Current_Lz(jj,ii,3)/(W(ii)-W(jj)+cmplx_i*yita)**2

                Orbital_HC_Lz_mpi(ik,2) = Orbital_HC_Lz_mpi(ik,2) + &
                    -cmplx_i*(occ(ii)-occ(jj))*V_Ham(ii,jj,3)*Orbital_Current_Lz(jj,ii,1)/(W(ii)-W(jj)+cmplx_i*yita)**2
                

                
                Orbital_HC_Lx_mpi(ik,3) = Orbital_HC_Lx_mpi(ik,3) + &
                    -cmplx_i*(occ(ii)-occ(jj))*V_Ham(ii,jj,1)*Orbital_Current_Lx(jj,ii,2)/(W(ii)-W(jj)+cmplx_i*yita)**2

                Orbital_HC_Lx_mpi(ik,1) = Orbital_HC_Lx_mpi(ik,1) + &
                    -cmplx_i*(occ(ii)-occ(jj))*V_Ham(ii,jj,2)*Orbital_Current_Lx(jj,ii,3)/(W(ii)-W(jj)+cmplx_i*yita)**2

                Orbital_HC_Lx_mpi(ik,2) = Orbital_HC_Lx_mpi(ik,2) + &
                    -cmplx_i*(occ(ii)-occ(jj))*V_Ham(ii,jj,3)*Orbital_Current_Lx(jj,ii,1)/(W(ii)-W(jj)+cmplx_i*yita)**2

                
                
                Orbital_HC_Ly_mpi(ik,3) = Orbital_HC_Ly_mpi(ik,3) + &
                    -cmplx_i*(occ(ii)-occ(jj))*V_Ham(ii,jj,1)*Orbital_Current_Ly(jj,ii,2)/(W(ii)-W(jj)+cmplx_i*yita)**2

                Orbital_HC_Ly_mpi(ik,1) = Orbital_HC_Ly_mpi(ik,1) + &
                    -cmplx_i*(occ(ii)-occ(jj))*V_Ham(ii,jj,2)*Orbital_Current_Ly(jj,ii,3)/(W(ii)-W(jj)+cmplx_i*yita)**2

                Orbital_HC_Ly_mpi(ik,2) = Orbital_HC_Ly_mpi(ik,2) + &
                    -cmplx_i*(occ(ii)-occ(jj))*V_Ham(ii,jj,3)*Orbital_Current_Ly(jj,ii,1)/(W(ii)-W(jj)+cmplx_i*yita)**2
            enddo
        enddo

    enddo

#if defined (MPI)
    call mpi_allreduce(Orbital_HC_Lz_mpi, Orbital_HC_Lz, size(Orbital_HC_Lz), mpi_dc,mpi_sum,mpi_cmw,ierr)
    call mpi_allreduce(Orbital_HC_Lx_mpi, Orbital_HC_Lx, size(Orbital_HC_Lx), mpi_dc,mpi_sum,mpi_cmw,ierr)
    call mpi_allreduce(Orbital_HC_Ly_mpi, Orbital_HC_Ly, size(Orbital_HC_Ly), mpi_dc,mpi_sum,mpi_cmw,ierr)
#else
    Orbital_HC_Lz = Orbital_HC_Lz_mpi
    Orbital_HC_Lx = Orbital_HC_Lx_mpi
    Orbital_HC_Ly = Orbital_HC_Ly_mpi
#endif

    !Orbital_HC_Lz = Orbital_HC_Lz_mpi


    if (cpuid.eq.0) then
        outfileindex= outfileindex+ 1
        open(unit=outfileindex, file='Intra_Orbital_Hall_Curvature_line.dat')
        !write(outfileindex, '("#",10a)')' the symmetric part of linear conductivity'
        write(outfileindex, "('#column', i5, 3000i16)")(i, i=1, 19)
        write(outfileindex, '("#",a13, 20a16)') 'Omega', 'Re[\sigma_Lz_xy]', 'Im[\sigma_Lz_xy]', 'Re[\sigma_Lz_yz]', &
                            'Im[\sigma_Lz_yz]', 'Re[\sigma_Lz_zx]', 'Im[\sigma_Lz_zx]', &
                            'Re[\sigma_Lx_xy]', 'Im[\sigma_Lx_xy]', 'Re[\sigma_Lx_yz]', &
                            'Im[\sigma_Lx_yz]', 'Re[\sigma_Lx_zx]', 'Im[\sigma_Lx_zx]', &
                            'Re[\sigma_Ly_xy]', 'Im[\sigma_Ly_xy]', 'Re[\sigma_Ly_yz]', &
                            'Im[\sigma_Ly_yz]', 'Re[\sigma_Ly_zx]', 'Im[\sigma_Ly_zx]'
        do ik=1,knv3
            OHC_out(1) = real(Orbital_HC_Lz(ik,3)) /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Echarge**2/hbar
            OHC_out(2) = aimag(Orbital_HC_Lz(ik,3)) /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Echarge**2/hbar
            OHC_out(3) = real(Orbital_HC_Lz(ik,1)) /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Echarge**2/hbar
            OHC_out(4) = aimag(Orbital_HC_Lz(ik,1)) /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Echarge**2/hbar
            OHC_out(5) = real(Orbital_HC_Lz(ik,2)) /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Echarge**2/hbar
            OHC_out(6) = aimag(Orbital_HC_Lz(ik,2))/dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Echarge**2/hbar

            OHC_out(7) = real(Orbital_HC_Lx(ik,3)) /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Echarge**2/hbar
            OHC_out(8) = aimag(Orbital_HC_Lx(ik,3)) /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Echarge**2/hbar
            OHC_out(9) = real(Orbital_HC_Lx(ik,1)) /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Echarge**2/hbar
            OHC_out(10) = aimag(Orbital_HC_Lx(ik,1)) /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Echarge**2/hbar
            OHC_out(11) = real(Orbital_HC_Lx(ik,2)) /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Echarge**2/hbar
            OHC_out(12) = aimag(Orbital_HC_Lx(ik,2))/dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Echarge**2/hbar

            OHC_out(13) = real(Orbital_HC_Ly(ik,3)) /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Echarge**2/hbar
            OHC_out(14) = aimag(Orbital_HC_Ly(ik,3)) /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Echarge**2/hbar
            OHC_out(15) = real(Orbital_HC_Ly(ik,1)) /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Echarge**2/hbar
            OHC_out(16) = aimag(Orbital_HC_Ly(ik,1)) /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Echarge**2/hbar
            OHC_out(17) = real(Orbital_HC_Ly(ik,2)) /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Echarge**2/hbar
            OHC_out(18) = aimag(Orbital_HC_Ly(ik,2))/dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Echarge**2/hbar

            
            write(outfileindex, '(200E16.8)')  ik  , OHC_out
        enddo
        close(outfileindex)
    endif




end subroutine Intra_Orbital_hall_curvature_line



















subroutine Intra_Orbital_texture
    use wmpi
    use para
    implicit none

    integer :: knv3, ik, ii, jj, mm, nn, ss, tt, ierr, imu, i, j
    integer :: ikx, iky, ikz
    real(dp) :: yita
    complex(dp) :: cmplx_i, cmplx_1, cmplx_0, omega
    real(dp) :: k(3)
    real(dp) :: time_start, time_end
    real(dp) :: Fermi_energy

    real(dp), allocatable :: W(:)
    complex(dp), allocatable :: Hamk_bulk(:, :)
    complex(dp), allocatable :: UU(:, :)

    complex(dp), allocatable :: D_Ham(:, :, :)
    complex(dp), allocatable :: V_Ham(:, :, :)
    
    complex(dp), allocatable :: L_wann_z(:, :)
    complex(dp), allocatable :: L_Ham_z(:, :)
    
    complex(dp), allocatable :: L_wann_x(:, :)
    complex(dp), allocatable :: L_Ham_x(:, :)

    complex(dp), allocatable :: L_wann_y(:, :)
    complex(dp), allocatable :: L_Ham_y(:, :)

    complex(dp), allocatable :: Orbital_Current_Lz(:, :, :)
    complex(dp), allocatable :: Orbital_Current_Lx(:, :, :)
    complex(dp), allocatable :: Orbital_Current_Ly(:, :, :)
    
    complex(dp), allocatable :: Orbital_Curvature_EF(:,:)


    complex(dp), allocatable :: L_texture_x(:, :)
    complex(dp), allocatable :: L_texture_y(:, :)
    complex(dp), allocatable :: L_texture_z(:, :)
    real(dp), allocatable :: Energy_list(:,:)

    real(dp), allocatable :: mulist(:)
    real(dp), allocatable :: occ(:)
    real(dp) :: OHC_out(7)

    knv3= nk3_band
    cmplx_i = (0.0d0, 1.0d0)
    cmplx_1 = (1.0d0, 0.0d0)
    cmplx_0 = (0.0d0, 0.0d0)


    allocate( W(Num_wann) )
    allocate(Hamk_bulk(Num_wann, Num_wann))
    allocate( UU(Num_wann, Num_wann))


    allocate(mulist(OmegaNum))
    allocate(D_Ham(Num_wann, Num_wann, 3))
    allocate(V_Ham(Num_wann, Num_wann, 3))
    allocate(L_Ham_z(Num_wann, Num_wann),L_wann_z(Num_wann, Num_wann))
    allocate(L_Ham_x(Num_wann, Num_wann),L_wann_x(Num_wann, Num_wann))
    allocate(L_Ham_y(Num_wann, Num_wann),L_wann_y(Num_wann, Num_wann))

    allocate(L_texture_x(knv3,Num_wann), L_texture_y(knv3,Num_wann), L_texture_z(knv3,Num_wann))
    allocate(Energy_list(knv3,Num_wann))
    allocate( occ(Num_wann))

    allocate(Orbital_Current_Lz(Num_wann, Num_wann, 3))
    allocate(Orbital_Current_Lx(Num_wann, Num_wann, 3))
    allocate(Orbital_Current_Ly(Num_wann, Num_wann, 3))
    allocate(Orbital_Curvature_EF(knv3, 3))


    yita = 0.002

    
    L_wann_z = cmplx_0
    L_wann_x = cmplx_0
    L_wann_y = cmplx_0

    L_texture_x = cmplx_0
    L_texture_y = cmplx_0
    L_texture_z = cmplx_0
    
    Orbital_Curvature_EF = cmplx_0

    !p
    !L_wann_z(3,3) = 1
    !L_wann_z(4,4) = 1
    !L_wann_z(2,2) = 1
    !L_wann_z(7,7) = 1
    !L_wann_z(8,8) = 1
    !L_wann_z(6,6) = 1

    L_wann_z(2,3) = cmplx_i
    L_wann_z(3,2) = -cmplx_i
    L_wann_z(6,7) = cmplx_i
    L_wann_z(7,6) = -cmplx_i

    

    !Lx
    L_wann_x(3,4) = cmplx_i
    L_wann_x(4,3) = -cmplx_i
    L_wann_x(7,8) = cmplx_i
    L_wann_x(8,7) = -cmplx_i



    !Ly
    L_wann_y(4,2) = cmplx_i
    L_wann_y(2,4) = -cmplx_i
    L_wann_y(8,6) = cmplx_i
    L_wann_y(6,8) = -cmplx_i
    !s
    !L_wann_y(1,1) = 1
    !L_wann_y(5,5) = 1

    do imu= 1, OmegaNum 
        mulist(imu) = OmegaMin + ( OmegaMax - OmegaMin) * imu / OmegaNum
    enddo



    call now(time_start) 
    !do ik= 1, knv3
    do ik= 1, knv3
        call now(time_end) 
        time_start= time_end
        

        Hamk_bulk = cmplx_0
        V_Ham = cmplx_0
        D_Ham = cmplx_0
        L_Ham_z = cmplx_0
        L_Ham_x = cmplx_0
        L_Ham_y = cmplx_0 
        
        Orbital_Current_Lz = cmplx_0 
        Orbital_Current_Lx = cmplx_0 
        Orbital_Current_Ly = cmplx_0 

        UU = cmplx_0
        W = 0.0_dp


        


        k = kpath_3d(:, ik)
        
        call ham_bulk_atomicgauge(k, Hamk_bulk)

        UU=Hamk_bulk
        call eigensystem_c( 'V', 'U', Num_wann, UU, W)

        call dHdk_atomicgauge_Ham(k, UU, V_Ham)
        !call get_Dmn_Ham(W, V_Ham, D_Ham)

        !call get_Lmn(W, D_Ham, V_Ham, L_Ham_z)

        call rotation_to_Ham_basis(UU, L_wann_z, L_Ham_z)
        call rotation_to_Ham_basis(UU, L_wann_x, L_Ham_x)
        call rotation_to_Ham_basis(UU, L_wann_y, L_Ham_y)
        

        do ss=1,3
            do mm = 1,Num_wann
                do nn = 1,Num_wann
                    do tt = 1,Num_wann
                    
                        Orbital_Current_Lz(mm,nn,ss) = Orbital_Current_Lz(mm,nn,ss) + &
                                0.5*( V_Ham(mm,tt,ss) * L_Ham_z(tt,nn) + L_Ham_z(mm,tt) * V_Ham(tt,nn,ss) )
                        
                        Orbital_Current_Lx(mm,nn,ss) = Orbital_Current_Lx(mm,nn,ss) + &
                                0.5*( V_Ham(mm,tt,ss) * L_Ham_x(tt,nn) + L_Ham_x(mm,tt) * V_Ham(tt,nn,ss) )

                        Orbital_Current_Ly(mm,nn,ss) = Orbital_Current_Ly(mm,nn,ss) + &
                                0.5*( V_Ham(mm,tt,ss) * L_Ham_y(tt,nn) + L_Ham_y(mm,tt) * V_Ham(tt,nn,ss) )
                    enddo
                enddo
            enddo
        enddo
       
        occ = 0.0_dp
        do i = 1, Num_wann
            !if (W(i) < Fermi_energy) occ(i) = 1
            
            do j = 1, Num_wann
                if(i==j) cycle
                L_texture_z(ik, i) = L_texture_z(ik, i) - 1*cmplx_i*V_Ham(i,j,1)*Orbital_Current_Lz(j,i,2)/(W(i)-W(j)+cmplx_i*yita)**2 !L_Ham_z(i,i)!*L_Ham_z(i+2,i)
                L_texture_x(ik, i) = L_texture_x(ik, i) - 1*cmplx_i*V_Ham(i,j,1)*Orbital_Current_Lx(j,i,2)/(W(i)-W(j)+cmplx_i*yita)**2!L_Ham_x(i,i)!*L_Ham_x(i+2,i)
                L_texture_y(ik, i) = L_texture_y(ik, i) - 1*cmplx_i*V_Ham(i,j,1)*Orbital_Current_Ly(j,i,2)/(W(i)-W(j)+cmplx_i*yita)**2!L_Ham_y(i,i)!*L_Ham_y(i+2,i)
                Energy_list(ik, i) = W(i)

            enddo
        enddo
        Fermi_energy = 0
        do ii = 1, Num_wann
            if(W(ii)<Fermi_energy) then
                Orbital_Curvature_EF(ik,1) =  Orbital_Curvature_EF(ik,1) + L_texture_x(ik,ii)
                Orbital_Curvature_EF(ik,2) =  Orbital_Curvature_EF(ik,2) + L_texture_y(ik,ii)
                Orbital_Curvature_EF(ik,3) =  Orbital_Curvature_EF(ik,3) + L_texture_z(ik,ii)
            endif
        enddo

        
        

    enddo

    

    !Orbital_HC_Lz = Orbital_HC_Lz_mpi


    if (cpuid.eq.0) then
        outfileindex= outfileindex+ 1
        open(unit=outfileindex, file='Intra_Orbital_Texture_line.dat')
        !write(outfileindex, '("#",10a)')' the symmetric part of linear conductivity'
        write(outfileindex, "('#column', i5, 3000i16)")(i, i=1, 8)
        write(outfileindex, '("#",a13, 20a16)') 'k_index','energy','Lx.real', 'Lx.imag', 'Ly.real', 'Ly.imag', 'Lz.real', 'Lz.imag'
        do ii=1,Num_wann
            write(outfileindex, *)  'band', ii
            do ik=1,knv3
                OHC_out(1) = Energy_list(ik,ii)/eV2Hartree
                OHC_out(2) = real(L_texture_x(ik,ii)) /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Echarge**2/hbar
                OHC_out(3) = aimag(L_texture_x(ik,ii)) /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Echarge**2/hbar
                OHC_out(4) = real(L_texture_y(ik,ii)) /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Echarge**2/hbar
                OHC_out(5) = aimag(L_texture_y(ik,ii)) /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Echarge**2/hbar
                OHC_out(6) = real(L_texture_z(ik,ii)) /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Echarge**2/hbar
                OHC_out(7) = aimag(L_texture_z(ik,ii))/dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Echarge**2/hbar

                

                
                write(outfileindex, '(200E16.8)')  ik  , OHC_out
            enddo
        enddo
        close(outfileindex)


        outfileindex= outfileindex+ 1
        open(unit=outfileindex, file='Intra_Orbital_curvature_bench.dat')
        !write(outfileindex, '("#",10a)')' the symmetric part of linear conductivity'
        write(outfileindex, "('#column', i5, 3000i16)")(i, i=1, 8)
        write(outfileindex, '("#",a13, 20a16)') 'k_index','Lx.real', 'Lx.imag', 'Ly.real', 'Ly.imag', 'Lz.real', 'Lz.imag', '000'
        
        do ik=1,knv3
            OHC_out(1) = real(Orbital_Curvature_EF(ik,1)) /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Echarge**2/hbar
            OHC_out(2) = aimag(Orbital_Curvature_EF(ik,1)) /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Echarge**2/hbar
            OHC_out(3) = real(Orbital_Curvature_EF(ik,2)) /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Echarge**2/hbar
            OHC_out(4) = aimag(Orbital_Curvature_EF(ik,2)) /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Echarge**2/hbar
            OHC_out(5) = real(Orbital_Curvature_EF(ik,3)) /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Echarge**2/hbar
            OHC_out(6) = aimag(Orbital_Curvature_EF(ik,3))/dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume/Bohr_radius/100*Echarge**2/hbar
            OHC_out(7) = 0
            

            
            write(outfileindex, '(200E16.8)')  ik  , OHC_out
        
        enddo
        close(outfileindex)
    endif




end subroutine Intra_Orbital_texture


