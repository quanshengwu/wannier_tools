  

  subroutine sigma_AHC_old
     !> Calculate anomalous hall conductivity AHC
     !
     !> ref : Physical Review B 74, 195118(2006)
     !
     !> eqn (34)
     !
     !> Dec. 05 2017 by Quansheng Wu @ EPFL
     !
     ! Copyright (c) 2017 QuanSheng Wu. All rights reserved.

     use wmpi
     use para
     implicit none
    
     integer :: iR, ik, ikx, iky, ikz
     integer :: m, n, i, j, ie
     integer :: ierr, knv3

     real(dp) :: kdotr, mu
     real(dp) :: k(3)

     real(dp) :: time_start, time_end

     ! eigen value of H
	  real(dp), allocatable :: W(:)
     complex(dp), allocatable :: Hamk_bulk(:, :)
     complex(dp), allocatable :: Amat(:, :)
     complex(dp), allocatable :: UU(:, :)
     complex(dp), allocatable :: UU_dag(:, :)

     !> velocities
     complex(dp), allocatable :: vx(:, :)
     complex(dp), allocatable :: vy(:, :)
     complex(dp), allocatable :: vz(:, :)
     complex(dp), allocatable :: DHDk(:, :, :)
     complex(dp), allocatable :: DHDkdag(:, :, :)
    
     !> Berry curvature
     complex(dp) :: Omega
     complex(dp) :: ratio

     !> conductivity  dim= OmegaNum
     real(dp), allocatable :: energy(:)
     real(dp), allocatable :: sigma_tensor_ahc(:, :, :)
     real(dp), allocatable :: sigma_tensor_ahc_mpi(:, :, :)

     allocate( W       (Num_wann))
     allocate( vx      (Num_wann, Num_wann))
     allocate( vy      (Num_wann, Num_wann))
     allocate( vz      (Num_wann, Num_wann))
     allocate( Hamk_bulk(Num_wann, Num_wann))
     allocate( Amat(Num_wann, Num_wann))
     allocate( UU(Num_wann, Num_wann))
     allocate( UU_dag(Num_wann, Num_wann))
     allocate( DHDk    (Num_wann, Num_wann, 3))
     allocate( DHDkdag (Num_wann, Num_wann, 3))
     allocate( energy(OmegaNum))
     allocate( sigma_tensor_ahc    (3, 3, OmegaNum))
     allocate( sigma_tensor_ahc_mpi(3, 3, OmegaNum))
     sigma_tensor_ahc    = 0d0
     sigma_tensor_ahc_mpi= 0d0
     omega= 0d0
     vx=0d0
     vy=0d0
     vz=0d0
     Hamk_bulk=0d0
     Amat= 0d0
     UU_dag=0d0
     UU= 0d0
     DHDk= 0d0
     DHDkdag= 0d0
     
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


        ! calculation bulk hamiltonian
        call ham_bulk_latticegauge(k, Hamk_bulk)

        !> diagonalization by call zheev in lapack
        W= 0d0
        UU=Hamk_bulk
        call eigensystem_c( 'V', 'U', Num_wann, UU, W)
       !call zhpevx_pack(hamk_bulk,Num_wann, W, UU)

        UU_dag= conjg(transpose(UU))

        vx= 0d0
        vy= 0d0
        vz= 0d0
        do iR= 1, Nrpts
           kdotr= k(1)*irvec(1,iR) + k(2)*irvec(2,iR) + k(3)*irvec(3,iR)
           ratio= zi*exp(pi2zi*kdotr)/ndegen(iR)
           vx= vx+ crvec(1, iR)*HmnR(:,:,iR)*ratio
           vy= vy+ crvec(2, iR)*HmnR(:,:,iR)*ratio
           vz= vz+ crvec(3, iR)*HmnR(:,:,iR)*ratio
        enddo ! iR

        !> unitility rotate velocity
        call mat_mul(Num_wann, vx, UU, Amat) 
        call mat_mul(Num_wann, UU_dag, Amat, vx) 
        call mat_mul(Num_wann, vy, UU, Amat) 
        call mat_mul(Num_wann, UU_dag, Amat, vy) 
        call mat_mul(Num_wann, vz, UU, Amat) 
        call mat_mul(Num_wann, UU_dag, Amat, vz) 

        !> calculate conductivity for each chemical potential
        do ie= 1, OmegaNum
           mu= energy(ie) 

           DHDk= 0d0
           DHDkdag= 0d0
           do m= 1, Num_wann
              do n= 1, Num_wann
                 if (W(n)> mu .and. W(m)< mu) then
                !if (n> Numoccupied .and. m<= Numoccupied) then  ! "=" a bug reported by Linlin Wang
                    DHDk(n, m, 1)= zi*vx(n, m)/(W(m)-W(n))
                    DHDk(n, m, 2)= zi*vy(n, m)/(W(m)-W(n))
                    DHDk(n, m, 3)= zi*vz(n, m)/(W(m)-W(n))
                 else
                    DHDk(n, m, 1)= 0d0
                    DHDk(n, m, 2)= 0d0
                    DHDk(n, m, 3)= 0d0
                 endif
              enddo ! m
           enddo ! n
   
           do m= 1, Num_wann
              do n= 1, Num_wann
                 if (W(m)> mu .and. W(n)< mu) then
                !if (m>Numoccupied .and. n<=Numoccupied) then  ! "=" a bug reported by Linlin Wang
                    DHDkdag(n, m, 1)= zi*vx(n, m)/(W(m)-W(n))
                    DHDkdag(n, m, 2)= zi*vy(n, m)/(W(m)-W(n))
                    DHDkdag(n, m, 3)= zi*vz(n, m)/(W(m)-W(n))
                 else              
                    DHDkdag(n, m, 1)= 0d0
                    DHDkdag(n, m, 2)= 0d0
                    DHDkdag(n, m, 3)= 0d0
                 endif
              enddo ! m
           enddo ! n
   
           !> rotate DHDk and DHDkdag to diagonal basis
           do i=1, 3
              call mat_mul(Num_wann, DHDk(:, :, i), UU_dag, Amat) 
              call mat_mul(Num_wann, UU, Amat, DHDk(:, :, i)) 
              call mat_mul(Num_wann, DHDkdag(:, :, i), UU_dag, Amat) 
              call mat_mul(Num_wann, UU, Amat, DHDkdag(:, :, i)) 
           enddo
   
           do i=1, 3
              do j=1, 3
                 call mat_mul(Num_wann, DHDk(:, :, i), DHDkdag(:, :, j), Amat)
                 call Im_trace(Num_wann, Amat, Omega)
                 sigma_tensor_ahc_mpi(i, j, ie)= sigma_tensor_ahc_mpi(i, j, ie)+ Omega
              enddo ! j
           enddo ! i
   
        enddo ! ie
     enddo ! ik

#if defined (MPI)
     call mpi_allreduce(sigma_tensor_ahc_mpi,sigma_tensor_ahc,size(sigma_tensor_ahc),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
     sigma_tensor_ahc= sigma_tensor_ahc_mpi
#endif
     sigma_tensor_ahc= sigma_tensor_ahc/dble(knv3)/Origin_cell%CellVolume*24341d0*2d0  ! in (Omega*cm)^-1

     outfileindex= outfileindex+ 1
     if (cpuid.eq.0) then
        open(unit=outfileindex, file='sigma_ahe.txt')
        write(outfileindex, '("#",20a16)')'E(eV)', 'xx', 'xy', 'xz', 'yx', 'yy', 'yz', 'zx', 'zy','zz'
        do ie=1, OmegaNum
           write(outfileindex, '(200f16.8)')energy(ie), ((sigma_tensor_ahc(i, j, ie), i=1, 3), j=1, 3)
        enddo
        close(outfileindex)
     endif

     deallocate( W       )
     deallocate( vx      )
     deallocate( vy      )
     deallocate( vz      )
     deallocate( Hamk_bulk)
     deallocate( Amat)
     deallocate( UU)
     deallocate( UU_dag)
     deallocate( DHDk    )
     deallocate( DHDkdag )
     deallocate( energy)
     deallocate( sigma_tensor_ahc    )
     deallocate( sigma_tensor_ahc_mpi)
 
     return
  end subroutine sigma_AHC_old

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
     complex(dp), allocatable :: Amat(:, :)
     complex(dp), allocatable :: UU(:, :)
     complex(dp), allocatable :: UU_dag(:, :)

     !> velocities
     complex(dp), allocatable :: vx(:, :), vy(:, :), vz(:, :)
    
     !> Berry curvature
     complex(dp), allocatable :: Omega_x(:), Omega_y(:), Omega_z(:)
     complex(dp), allocatable :: Omega_x_t(:), Omega_y_t(:), Omega_z_t(:)
     complex(dp) :: ratio

     !> conductivity  dim= OmegaNum
     real(dp), allocatable :: energy(:)
     real(dp), allocatable :: sigma_tensor_ahc(:, :)
     real(dp), allocatable :: sigma_tensor_ahc_mpi(:, :)
     
     !> Fermi-Dirac distribution
     real(dp), external :: fermi


     allocate( W (Num_wann))
     allocate( vx(Num_wann, Num_wann), vy(Num_wann, Num_wann), vz(Num_wann, Num_wann))
     allocate( Hamk_bulk(Num_wann, Num_wann))
     allocate( Amat(Num_wann, Num_wann))
     allocate( UU(Num_wann, Num_wann))
     allocate( UU_dag(Num_wann, Num_wann))
     allocate( energy(OmegaNum))
     allocate( sigma_tensor_ahc    (3, OmegaNum))
     allocate( sigma_tensor_ahc_mpi(3, OmegaNum))
     allocate(Omega_x(Num_wann), Omega_y(Num_wann), Omega_z(Num_wann))
     allocate(Omega_x_t(Num_wann), Omega_y_t(Num_wann), Omega_z_t(Num_wann))
     sigma_tensor_ahc    = 0d0
     sigma_tensor_ahc_mpi= 0d0
     vx=0d0
     vy=0d0
     vz=0d0
     Hamk_bulk=0d0
     Amat= 0d0
     UU_dag=0d0
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
        call ham_bulk_latticegauge(k, Hamk_bulk)
   
        !> diagonalization by call zheev in lapack
        UU=Hamk_bulk
        call eigensystem_c( 'V', 'U', Num_wann, UU, W)
  
        vx= 0d0; vy= 0d0; vz= 0d0
        UU_dag= conjg(transpose(UU))
        do iR= 1, Nrpts
           kdotr= k(1)*irvec(1,iR) + k(2)*irvec(2,iR) + k(3)*irvec(3,iR)
           vx= vx+ zi*crvec(1, iR)*HmnR(:,:,iR)*Exp(pi2zi*kdotr)/ndegen(iR)
           vy= vy+ zi*crvec(2, iR)*HmnR(:,:,iR)*Exp(pi2zi*kdotr)/ndegen(iR)
           vz= vz+ zi*crvec(3, iR)*HmnR(:,:,iR)*Exp(pi2zi*kdotr)/ndegen(iR)
        enddo ! iR
   
        !> unitility rotate velocity
        call mat_mul(Num_wann, vx, UU, Amat) 
        call mat_mul(Num_wann, UU_dag, Amat, vx) 
        call mat_mul(Num_wann, vy, UU, Amat) 
        call mat_mul(Num_wann, UU_dag, Amat, vy) 
        call mat_mul(Num_wann, vz, UU, Amat) 
        call mat_mul(Num_wann, UU_dag, Amat, vz) 
   
        Omega_x=0d0;Omega_y=0d0; Omega_z=0d0
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
   
        !> consider the Fermi-distribution according to the brodening Earc_eta
        Beta_fake= 1d0/Eta_Arc

        do ie=1, OmegaNum
           mu = energy(ie)
           do m= 1, Num_wann
              Omega_x_t(m)= Omega_x(m)*fermi(W(m)-mu, Beta_fake)
              Omega_y_t(m)= Omega_y(m)*fermi(W(m)-mu, Beta_fake)
              Omega_z_t(m)= Omega_z(m)*fermi(W(m)-mu, Beta_fake)
           enddo
           sigma_tensor_ahc_mpi(1, ie)= sigma_tensor_ahc_mpi(1, ie)+ real(sum(Omega_x_t))
           sigma_tensor_ahc_mpi(2, ie)= sigma_tensor_ahc_mpi(2, ie)+ real(sum(Omega_y_t))
           sigma_tensor_ahc_mpi(3, ie)= sigma_tensor_ahc_mpi(3, ie)+ real(sum(Omega_z_t))
        enddo ! ie
     enddo ! ik

#if defined (MPI)
     call mpi_allreduce(sigma_tensor_ahc_mpi,sigma_tensor_ahc,size(sigma_tensor_ahc),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
     sigma_tensor_ahc= sigma_tensor_ahc_mpi
#endif

     sigma_tensor_ahc= sigma_tensor_ahc/dble(knv3)/Origin_cell%CellVolume*24341d0*kCubeVolume/Origin_cell%ReciprocalCellVolume ! in (Omega*cm)^-1

     outfileindex= outfileindex+ 1
     if (cpuid.eq.0) then
        open(unit=outfileindex, file='sigma_ahc.txt')
        write(outfileindex, '("#",a)')' Anomalous hall conductivity in unit of (Ohm*cm)^-1 (Omega*cm)^-1, and e^2/h'
        write(outfileindex, '("#",a13, 20a16)')'Eenergy (eV)', '\sigma_xy', '\sigma_yz', '\sigma_zx'
        do ie=1, OmegaNum
           write(outfileindex, '(200E16.8)')energy(ie), sigma_tensor_ahc(3, ie), &
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
        write(outfileindex, '(a, f10.5, a, f10.5, a)')'set xrange [', OmegaMin, ':', OmegaMax, ']'
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


     deallocate( W , vx, vy, vz, Hamk_bulk, Amat, UU, UU_dag, energy)
     deallocate( sigma_tensor_ahc, sigma_tensor_ahc_mpi)
 
     return
  end subroutine sigma_AHC
