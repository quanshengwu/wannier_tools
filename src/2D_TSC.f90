subroutine ham_slab_surface_zeeman(k, Bz_BdG_phase, Hamk_slab_surf_zeeman)
    ! This subroutine is used to caculate Hamiltonian for
    ! slab system with surface Zeeman splitting.
    !
    ! History
    !        4/18/2010 by Quansheng Wu
    !        6/21/2022 by Aiyun Luo
    !        4/07/2024 by Jingnan Hu
    ! Copyright (c) 2010 QuanSheng Wu. All rights reserved.

    use para
    implicit none

    ! loop index
    integer :: i1, i2, i, ii, j

    ! wave vector in 2d
    real(Dp), intent(in) :: k(2)

    real(Dp), intent(in) :: Bz_BdG_phase

    ! A vector perpendicular to the surface

    real(Dp) :: ez(3)
    real(dp), external :: norm

    ! Hamiltonian of slab system with surface Zeeman splitting
    complex(Dp), intent(out) :: Hamk_slab_surf_zeeman(Num_wann*nslab,Num_wann*nslab)

    ! the factor 2 is induced by spin
    complex(Dp), allocatable :: Hij(:, :, :)

    allocate( Hij(-ijmax:ijmax,Num_wann,Num_wann))

    call ham_qlayer2qlayer2(k,Hij)

    Hamk_slab_surf_zeeman= 0.0d0
    ! i1 column index
    do i1=1, nslab
       ! i2 row index
       do i2=1, nslab
         if (abs(i2-i1).le.ijmax)then
            Hamk_slab_surf_zeeman((i2-1)*Num_wann+1:(i2-1)*Num_wann+Num_wann,&
                     (i1-1)*Num_wann+1:(i1-1)*Num_wann+Num_wann )&
           = Hij(i1-i2,1:Num_wann,1:Num_wann)
         endif
       enddo ! i2
    enddo ! i1

    !> There are several types of Zeeman splitting

    ez(1)= (Rua_newcell(2)*Rub_newcell(3)-Rua_newcell(3)*Rub_newcell(2))
    ez(2)= (Rua_newcell(3)*Rub_newcell(1)-Rua_newcell(1)*Rub_newcell(3))
    ez(3)= (Rua_newcell(1)*Rub_newcell(2)-Rua_newcell(2)*Rub_newcell(1))
    ez(:)= ez(:)/norm(ez(:))

    if (cpuid==0) then
        write(stdout, *)'A vector perpendicular to the slab: '
        write(stdout, "(3f12.6)")ez(:)
    endif


    !> 1. Add Zeeman splitting only in the bottom slab
    if(Add_surf_zeeman_field == 1) then
        do i=1, Num_wann/2
            ii= i
            !>Bz_surf
            Hamk_slab_surf_zeeman(ii, ii)= Hamk_slab_surf_zeeman(ii, ii)+ Bz_BdG_phase * ez(3) * eV2Hartree
            Hamk_slab_surf_zeeman(ii+Num_wann/2, ii+Num_wann/2)= Hamk_slab_surf_zeeman(ii+Num_wann/2, ii+Num_wann/2)- Bz_BdG_phase * ez(3) * eV2Hartree
            !>Bx_surf, By_surf
            Hamk_slab_surf_zeeman(ii, ii+Num_wann/2)= Hamk_slab_surf_zeeman(ii, ii+Num_wann/2)+ Bz_BdG_phase * ez(1) * eV2Hartree- zi*Bz_BdG_phase * ez(2) * eV2Hartree
            Hamk_slab_surf_zeeman(ii+Num_wann/2, ii)= Hamk_slab_surf_zeeman(ii+Num_wann/2, ii)+ Bz_BdG_phase * ez(1) * eV2Hartree+ zi*Bz_BdG_phase * ez(2) * eV2Hartree
       enddo
    endif

    !> 2. Add Zeeman splitting only in the top slab
    if(Add_surf_zeeman_field == 2) then
        do i=1, Num_wann/2
            ii= (nslab-1)*Num_wann+ i
            !>Bz_surf
            Hamk_slab_surf_zeeman(ii, ii)= Hamk_slab_surf_zeeman(ii, ii)+ Bz_BdG_phase * ez(3) * eV2Hartree
            Hamk_slab_surf_zeeman(ii+Num_wann/2, ii+Num_wann/2)= Hamk_slab_surf_zeeman(ii+Num_wann/2, ii+Num_wann/2)- Bz_BdG_phase * ez(3) * eV2Hartree
            !>Bx_surf, By_surf
            Hamk_slab_surf_zeeman(ii, ii+Num_wann/2)= Hamk_slab_surf_zeeman(ii, ii+Num_wann/2) + Bz_BdG_phase * ez(1) * eV2Hartree- zi*Bz_BdG_phase * ez(2) * eV2Hartree
            Hamk_slab_surf_zeeman(ii+Num_wann/2, ii)= Hamk_slab_surf_zeeman(ii+Num_wann/2, ii) + Bz_BdG_phase * ez(1) * eV2Hartree+ zi*Bz_BdG_phase * ez(2) * eV2Hartree
       enddo
    endif    

    !> 3. Add Zeeman splitting in whole slab
    if(Add_surf_zeeman_field == 3) then
        do j=1, nslab
           do i=1, Num_wann/2
              ii= i+(j-1)*Num_wann
              !>Bz_surf
              Hamk_slab_surf_zeeman(ii, ii)= Hamk_slab_surf_zeeman(ii, ii)+ Bz_BdG_phase * ez(3) * eV2Hartree
              Hamk_slab_surf_zeeman(ii+Num_wann/2, ii+Num_wann/2)= Hamk_slab_surf_zeeman(ii+Num_wann/2, ii+Num_wann/2)- Bz_BdG_phase * ez(3) * eV2Hartree
              !>Bx_surf, By_surf
              Hamk_slab_surf_zeeman(ii, ii+Num_wann/2)= Hamk_slab_surf_zeeman(ii, ii+Num_wann/2)+ Bz_BdG_phase * ez(1) * eV2Hartree- zi*Bz_BdG_phase * ez(2) * eV2Hartree
              Hamk_slab_surf_zeeman(ii+Num_wann/2, ii)= Hamk_slab_surf_zeeman(ii+Num_wann/2, ii)+ Bz_BdG_phase * ez(1) * eV2Hartree+ zi*Bz_BdG_phase * ez(2) * eV2Hartree
           enddo
        enddo
    endif   

    !Hamk_slab_surf_zeeman(1, 1)=Hamk_slab_surf_zeeman(1, 1)+5*eV2Hartree
    !Hamk_slab_surf_zeeman(2, 2)=Hamk_slab_surf_zeeman(2, 2)+5*eV2Hartree
    !Hamk_slab_surf_zeeman(8, 8)=Hamk_slab_surf_zeeman(8, 8)+5*eV2Hartree
    !Hamk_slab_surf_zeeman(9, 9)=Hamk_slab_surf_zeeman(9, 9)+5*eV2Hartree

    !Hamk_slab_surf_zeeman((nslab-1)*Num_wann+6, (nslab-1)*Num_wann+6)=Hamk_slab_surf_zeeman((nslab-1)*Num_wann+6, (nslab-1)*Num_wann+6)-5*eV2Hartree
    !Hamk_slab_surf_zeeman((nslab-1)*Num_wann+7, (nslab-1)*Num_wann+7)=Hamk_slab_surf_zeeman((nslab-1)*Num_wann+7, (nslab-1)*Num_wann+7)-5*eV2Hartree
    !Hamk_slab_surf_zeeman((nslab-1)*Num_wann+13, (nslab-1)*Num_wann+13)=Hamk_slab_surf_zeeman((nslab-1)*Num_wann+13, (nslab-1)*Num_wann+13)-5*eV2Hartree
    !Hamk_slab_surf_zeeman((nslab-1)*Num_wann+14, (nslab-1)*Num_wann+14)=Hamk_slab_surf_zeeman((nslab-1)*Num_wann+14, (nslab-1)*Num_wann+14)-5*eV2Hartree


    ! check hermitcity
    do i1=1,nslab*Num_wann
    do i2=1,nslab*Num_wann
       if(abs(Hamk_slab_surf_zeeman(i1,i2)-conjg(Hamk_slab_surf_zeeman(i2,i1))).ge.1e-6)then
        write(stdout,*)'there are something wrong with Hamk_slab_surf_zeeman'
        stop
       endif
    enddo
    enddo

    deallocate(Hij)
 return
 end subroutine ham_slab_surface_zeeman

 subroutine ham_slab_BdG(k,Bz_BdG_phase,Hamk_slab_BdG)
    ! This subroutine is used to caculate Hamiltonian for
    ! slab BdG system .
    !
    ! History
    !       03/30/2022 by Aiyun Luo
    !        4/07/2024 by Jingnan Hu

    use para
    use Kronecker, only : KronProd     ! Kroneker product

    implicit none

    ! loop index
    integer :: i1, i2, i

    ! wave vector in 2d
    real(Dp), intent(in) :: k(2)

    real(Dp), intent(in) :: Bz_BdG_phase

    ! Delta_decay coefficient
    real(Dp) :: kd

    ! Hamiltonian of slab BdG system
    complex(Dp), intent(out) :: Hamk_slab_BdG(Num_wann_BdG*nslab,Num_wann_BdG*nslab)

    ! Hamiltonian of slab system with Zeeman splitting in the top or bottom surface
    complex(Dp), allocatable :: Hamk_slab(:,:)
    complex(Dp), allocatable :: Hamk_slab_minus(:,:)

    ! superconducting pairing strength, here we only consider the onsite s-wave pairing
    complex(Dp), allocatable :: Hamk_Delta(:,:)

    ! Pauli matrices
    complex(Dp) :: sigmax(2,2), sigmay(2,2), sigmaz(2,2), sigma0(2,2)

    ! indetify diagnoal matrix for construct s-wave pairing
    complex(Dp), allocatable :: I_nslab(:,:)
    complex(Dp), allocatable :: I_norb(:,:)
    complex(Dp), allocatable :: I_mu(:,:)
    complex(Dp), allocatable :: Hamk_temp(:,:)

    allocate(Hamk_slab(Num_wann*nslab,Num_wann*nslab))
    allocate(Hamk_slab_minus(Num_wann*nslab,Num_wann*nslab))
    allocate(I_nslab(nslab, nslab))
    allocate(I_norb(Num_wann/2, Num_wann/2))
    allocate(I_mu(nslab*Num_wann, nslab*Num_wann))
    allocate(Hamk_temp(Num_wann, Num_wann))
    allocate(Hamk_Delta(Num_wann*Nslab, Num_wann*Nslab))
    
    Hamk_slab_BdG  = 0.0d0
    Hamk_slab      = 0.0d0;    Hamk_slab_minus = 0.0d0
    I_nslab        = 0.0d0;    I_norb          = 0.0d0
    I_mu           = 0.0d0
    Hamk_temp      = 0.0d0;    Hamk_Delta      = 0.0d0

    sigmax(1,1)= (0.0d0, 0.0d0);   sigmax(1,2)= (1.0d0, 0.0d0)
    sigmax(2,1)= (1.0d0, 0.0d0);   sigmax(2,2)= (0.0d0, 0.0d0)
    
    sigmay(1,1)= (0.0d0, 0.0d0);   sigmay(1,2)= (0.0d0,-1.0d0)
    sigmay(2,1)= (0.0d0, 1.0d0);   sigmay(2,2)= (0.0d0, 0.0d0)

    sigmaz(1,1)= (1.0d0, 0.0d0);   sigmaz(1,2)= ( 0.0d0, 0.0d0)
    sigmaz(2,1)= (0.0d0, 0.0d0);   sigmaz(2,2)= (-1.0d0, 0.0d0)

    sigma0(1,1)= (1.0d0, 0.0d0);   sigma0(1,2)= (0.0d0, 0.0d0)
    sigma0(2,1)= (0.0d0, 0.0d0);   sigma0(2,2)= (1.0d0, 0.0d0)
    
    call ham_slab_surface_zeeman( k, Bz_BdG_phase, Hamk_slab)
    call ham_slab_surface_zeeman(-k, Bz_BdG_phase, Hamk_slab_minus)

    !> s-wave superconducting pairing for basis: up up dn dn (wannier90_hr.dat from vasp)
    !> s-wave pairing: Kron(I_nslab, Kron(i*sigmay*Delta_BdG, I_norb))
    call eye_mat(Num_wann/2, I_norb)
    call eye_mat(nslab*Num_wann,I_mu)

    !> Delta_decay (Ref: Theory of the Superconducting Proximity Effect)
    kd=sqrt((t_BdG/xi_BdG)*((t_BdG/xi_BdG)+(3/l_BdG)))


    if(Add_Delta_BdG == 1) then
       do i=1 ,nslab
          I_nslab(i,i)=cosh(kd*(i-1-nslab)*Ruc_newcell(3)/(10*Angstrom2atomic))/cosh(kd*(nslab)*Ruc_newcell(3)/(10*Angstrom2atomic))
       enddo !i
    endif

    if(Add_Delta_BdG == 2) then
       do i=1 ,nslab
          I_nslab(i,i)=cosh(kd*(i)*Ruc_newcell(3)/(10*Angstrom2atomic))/cosh(kd*(nslab)*Ruc_newcell(3)/(10*Angstrom2atomic))
       enddo !i
    endif

    Hamk_temp= KronProd(zi*Delta_BdG*sigmay*eV2Hartree, I_norb)
    Hamk_Delta= KronProd(I_nslab, Hamk_temp)

    !> constructing the slab BdG Hamiltonian with top or bottom surface exchange field
    !> basis: C1^dag, C2^dag, C1, C2 
    Hamk_slab_BdG(1:nslab*Num_wann,1:nslab*Num_wann)= Hamk_slab-mu_BdG*eV2Hartree*I_mu
    Hamk_slab_BdG(nslab*Num_wann+1:nslab*Num_wann_BdG,nslab*Num_wann+1:nslab*Num_wann_BdG)= -1.0d0*conjg(Hamk_slab_minus)+mu_BdG*eV2Hartree*I_mu
    !> add onsite s-wave pairing into slab BdG Hamiltonian
    Hamk_slab_BdG(1:nslab*Num_wann, nslab*Num_wann+1: nslab*Num_wann_BdG)= Hamk_Delta
    Hamk_slab_BdG(nslab*Num_wann+1: nslab*Num_wann_BdG, 1:nslab*Num_wann)= transpose(conjg(Hamk_Delta))

    ! check hermitcity
    do i1=1,nslab*Num_wann_BdG
    do i2=1,nslab*Num_wann_BdG
       if(abs(Hamk_slab_BdG(i1,i2)-conjg(Hamk_slab_BdG(i2,i1))).ge.1e-6)then
        write(stdout,*)'there are something wrong with Hamk_slab_BdG'
        stop
       endif
    enddo
    enddo

    deallocate(Hamk_slab)
    deallocate(Hamk_slab_minus)
    deallocate(I_nslab)
    deallocate(I_norb)
    deallocate(Hamk_temp)
    deallocate(Hamk_Delta)
 return
end subroutine ham_slab_BdG

subroutine eye_mat(ndim, A)
   use para, only : Dp
   implicit none

   integer, intent(in) :: ndim
   complex(Dp), intent(out) :: A(ndim, ndim)

   !> loop index
   integer :: i, j

   do i=1, ndim
   do j=1, ndim
       if(i.eq.j) then
          A(i, j)= 1.0d0
       else
          A(i, j)= 0.0d0
       endif
   enddo ! i
   enddo ! j

   return
end subroutine eye_mat

subroutine ek_slab_BdG
   !> This subroutine is used for calculating BdG energy 
   !> dispersion with wannier functions for 2D slab system
   !> Added by Aiyun Luo at 2022/03
   !
   ! Copyright (c) 2010 QuanSheng Wu. All rights reserved.
  
   use wmpi
   use para
   implicit none 

   ! loop index
   integer :: i, j, l, lwork, ierr, io

   real(Dp) :: k(2), emin, emax, maxweight

   ! time measurement
   real(dp) :: time_start, time_end, time_start0

   ! Delta_decay coefficient
   real(Dp) :: kd

   ! Delta_decay
   real(Dp),allocatable :: Delta_decay(:)

   ! parameters for zheev
   real(Dp), allocatable ::  rwork(:)
   complex(Dp), allocatable :: work(:)

   ! eigenvalue 
   real(Dp), allocatable :: eigenvalue_BdG(:)

   ! energy dispersion
   real(Dp),allocatable :: ekslab_BdG(:,:), ekslab_BdG_mpi(:,:)

   !> color for plot, surface state weight
   real(dp), allocatable :: surf_l_weight_BdG(:, :), surf_l_weight_BdG_mpi(:, :)
   real(dp), allocatable :: surf_r_weight_BdG(:, :), surf_r_weight_BdG_mpi(:, :)

   ! hamiltonian slab
   complex(Dp),allocatable ::CHamk_BdG(:,:)

   lwork= 16*Nslab*Num_wann_BdG
   ierr = 0


   allocate(eigenvalue_BdG(nslab*Num_wann_BdG))
   allocate( surf_l_weight_BdG (Nslab* Num_wann_BdG, knv2))
   allocate( surf_l_weight_BdG_mpi (Nslab* Num_wann_BdG, knv2))
   allocate( surf_r_weight_BdG (Nslab* Num_wann_BdG, knv2))
   allocate( surf_r_weight_BdG_mpi (Nslab* Num_wann_BdG, knv2))
   allocate(ekslab_BdG(Nslab*Num_wann_BdG,knv2))
   allocate(ekslab_BdG_mpi(Nslab*Num_wann_BdG,knv2))
   allocate(CHamk_BdG(nslab*Num_wann_BdG,nslab*Num_wann_BdG))
   allocate(work(lwork))
   allocate(rwork(lwork))
   allocate(Delta_decay(nslab))

   surf_l_weight_BdG= 0d0
   surf_l_weight_BdG_mpi= 0d0
   surf_r_weight_BdG= 0d0
   surf_r_weight_BdG_mpi= 0d0

   ! sweep k
   ekslab_BdG=0.0d0
   ekslab_BdG_mpi=0.0d0
   time_start= 0d0
   time_start0= 0d0
   call now(time_start0)
   time_start= time_start0
   time_end  = time_start0
   do i= 1+cpuid, knv2, num_cpu
      if (cpuid==0.and. mod(i/num_cpu, 4)==0) &
         write(stdout, '(a, i9, "  /", i10, a, f10.1, "s", a, f10.1, "s")') &
         ' Slabek: ik', i, knv2, ' time left', &
         (knv2-i)*(time_end- time_start)/num_cpu, &
         ' time elapsed: ', time_end-time_start0 

      call now(time_start)

      k= k2_path(i, :)
      chamk_BdG=0.0d0 

      call ham_slab_BdG(k,Bz_surf,Chamk_BdG)


      eigenvalue_BdG=0.0d0

      ! diagonal Chamk
      call eigensystem_c('V', 'U', Num_wann_BdG*Nslab, CHamk_BdG, eigenvalue_BdG)
     
      ekslab_BdG(:,i)=eigenvalue_BdG

      ! H*chamk(:,n)=E(n)*chamk(:,n)
      !> Nslab*Num_wann
      !> rho(:)=abs(chamk(:,n))**2
      !> (a1 o1, o2 o3, a2, o1, o2, o3; a1 o1, o2 o3, a2, o1, o2, o3), (a1 o1, o2 o3, a2, o1, o2, o3; a1 o1, o2 o3, a2, o1, o2, o3), (a1 o1, o2 o3, a2, o1, o2, o3; a1 o1, o2 o3, a2, o1, o2, o3), 
      do j=1, Nslab* Num_wann_BdG
         !> left is the bottom surface
         do l= 1, NBottomOrbitals
            io= BottomOrbitals(l) 
            surf_l_weight_BdG(j, i)= surf_l_weight_BdG(j, i) &
               + abs(CHamk_BdG(io, j))**2 & ! first slab -- electron
               + abs(CHamk_BdG(Num_wann*Nslab+io, j))**2 ! first slab -- hole
         enddo ! l sweeps the selected orbitals

         !> right is the top surface
         do l= 1, NTopOrbitals
            io= Num_wann*(Nslab-1)+ TopOrbitals(l)   
            surf_r_weight_BdG(j, i)= surf_r_weight_BdG(j, i) &
               + abs(CHamk_BdG(io, j))**2 & ! first slab -- electron
               + abs(CHamk_BdG(Num_wann*Nslab+io, j))**2 ! first slab -- hole
         enddo ! l sweeps the selected orbitals

      enddo ! j 
      call now(time_end)
   enddo ! i

#if defined (MPI)
   call mpi_allreduce(ekslab_BdG,ekslab_BdG_mpi,size(ekslab_BdG),&
                     mpi_dp,mpi_sum,mpi_cmw,ierr)
   call mpi_allreduce(surf_l_weight_BdG, surf_l_weight_BdG_mpi,size(surf_l_weight_BdG),&
                     mpi_dp,mpi_sum,mpi_cmw,ierr)
   call mpi_allreduce(surf_r_weight_BdG, surf_r_weight_BdG_mpi,size(surf_r_weight_BdG),&
                     mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
   ekslab_BdG_mpi= ekslab_BdG
   surf_l_weight_BdG_mpi= surf_l_weight_BdG
   surf_r_weight_BdG_mpi= surf_r_weight_BdG
#endif

   ekslab_BdG_mpi= ekslab_BdG_mpi/eV2Hartree

   ekslab_BdG=ekslab_BdG_mpi

   maxweight=maxval(surf_r_weight_BdG_mpi+ surf_l_weight_BdG_mpi)
   surf_l_weight_BdG= surf_l_weight_BdG_mpi/ maxweight
   surf_r_weight_BdG= surf_r_weight_BdG_mpi/ maxweight
   
   outfileindex= outfileindex+ 1
   if(cpuid==0)then
      open(unit=outfileindex, file='slabek_BdG.dat')
      write(outfileindex, "('#', a10, a15, 5X, 2a16 )")'# k', ' E', 'BS weight', 'TS weight'
      do j=1, Num_wann_BdG*Nslab
         do i=1, knv2
           !write(outfileindex,'(3f15.7, i8)')k2len(i), ekslab(j,i), &
           !   (surf_weight(j, i))
            write(outfileindex,'(2f15.7, 2f16.7)')k2len(i)*Angstrom2atomic, ekslab_BdG(j,i), &
               (surf_l_weight_BdG(j, i)), &
               (surf_r_weight_BdG(j, i))
         enddo
         write(outfileindex , *)''
      enddo
      close(outfileindex)
      write(stdout,*) 'calculate energy band  done'
   endif

   emin= minval(ekslab_BdG)-0.5d0
   emax= maxval(ekslab_BdG)+0.5d0
   !> write script for gnuplot
   outfileindex= outfileindex+ 1
   if (cpuid==0) then
      open(unit=outfileindex, file='slabek_BdG.gnu')
      write(outfileindex, '(a)')"set encoding iso_8859_1"
      write(outfileindex, '(a)')'#set terminal  postscript enhanced color'
      write(outfileindex, '(a)')"#set output 'slabek_BdG.eps'"
      write(outfileindex, '(3a)')'#set terminal  pngcairo truecolor enhanced', &
         '  font ",60" size 1920, 1680'
      write(outfileindex, '(3a)')'set terminal  png truecolor enhanced', &
         '  font ",60" size 1920, 1680'
      write(outfileindex, '(a)')"set output 'slabek_BdG.png'"
      write(outfileindex,'(2a)') 'set palette defined ( 0  "green", ', &
         '5 "yellow", 10 "red" )'
      write(outfileindex, '(a)')'set style data linespoints'
      write(outfileindex, '(a)')'unset ztics'
      write(outfileindex, '(a)')'unset key'
      write(outfileindex, '(a)')'set pointsize 0.8'
      write(outfileindex, '(a)')'set border lw 3 '
      write(outfileindex, '(a)')'set view 0,0'
      write(outfileindex, '(a)')'#set xtics font ",36"'
      write(outfileindex, '(a)')'#set ytics font ",36"'
      write(outfileindex, '(a)')'#set ylabel font ",36"'
      write(outfileindex, '(a)')'#set xtics offset 0, -1'
      write(outfileindex, '(a)')'set ylabel offset -1, 0 '
      write(outfileindex, '(a, f10.5, a)')'set xrange [0: ', maxval(k2len)*Angstrom2atomic, ']'
      if (index(Particle,'phonon')/=0) then
         write(outfileindex, '(a, f10.5, a)')'set yrange [0:', emax, ']'
         write(outfileindex, '(a)')'set ylabel "Frequency (THz)"'
      else
         write(outfileindex, '(a)')'set ylabel "Energy (eV)"'
         write(outfileindex, '(a, f10.5, a, f10.5, a)')'set yrange [', emin, ':', emax, ']'
      endif
      write(outfileindex, 202, advance="no") (trim(k2line_name(i)), k2line_stop(i)*Angstrom2atomic, i=1, nk2lines)
      write(outfileindex, 203)trim(k2line_name(nk2lines+1)), k2line_stop(nk2lines+1)*Angstrom2atomic

      do i=1, nk2lines-1
         if (index(Particle,'phonon')/=0) then
            write(outfileindex, 204)k2line_stop(i+1)*Angstrom2atomic, 0.0, k2line_stop(i+1)*Angstrom2atomic, emax
         else
            write(outfileindex, 204)k2line_stop(i+1)*Angstrom2atomic, emin, k2line_stop(i+1)*Angstrom2atomic, emax
         endif
      enddo
      write(outfileindex, '(a)')'#rgb(r,g,b) = int(r)*65536 + int(g)*256 + int(b)'
      write(outfileindex, '(2a)')"#plot 'slabek_BdG.dat' u 1:2:(rgb(255,$3, 3)) ",  &
          "w lp lw 2 pt 7  ps 1 lc rgb variable"
      write(outfileindex, '(2a)')"# (a) "
      write(outfileindex, '(2a)')"# plot the top and bottom surface's weight together"
      write(outfileindex, '(2a)')"#plot 'slabek_BdG.dat' u 1:2:($3+$4) ",  &
          "w lp lw 2 pt 7  ps 1 lc palette"
      write(outfileindex, '(2a)')"# (b) "
      write(outfileindex, '(2a)') &
         "# plot top and bottom surface's weight with red and blue respectively"
      write(outfileindex,'(2a)') 'set palette defined ( -1  "blue", ', &
         '0 "grey", 1 "red" )'
      write(outfileindex, '(2a)')"plot 'slabek_BdG.dat' u 1:2:($4-$3) ",  &
          "w lp lw 2 pt 7  ps 1 lc palette"

      !write(outfileindex, '(2a)')"splot 'slabek.dat' u 1:2:3 ",  &
      !   "w lp lw 2 pt 13 palette"
      close(outfileindex)
   endif

   !> Delta_decay (Ref: Theory of the Superconducting Proximity Effect)
   kd=sqrt((t_BdG/xi_BdG)*((t_BdG/xi_BdG)+(3/l_BdG)))

   do i = 1, Nslab
      Delta_decay(i)=Delta_BdG*cosh(kd*(i-1-nslab)*Ruc_newcell(3)/(10*Angstrom2atomic))/cosh(kd*(nslab)*Ruc_newcell(3)/(10*Angstrom2atomic))
   enddo

   outfileindex= outfileindex+ 1
   if(cpuid==0)then
      open(unit=outfileindex, file='Delta_decay.dat')
      write(outfileindex, "('#', a15, a15 )")'d(layer)', ' Delta'
      do i=1, Nslab
         write(outfileindex,'(i10, f16.7)') i, Delta_decay(i)
      enddo
      close(outfileindex)
      write(stdout,*) 'Plot the Delta decay curve'
   endif

   outfileindex= outfileindex+ 1
   if (cpuid==0) then
      open(unit=outfileindex, file='Delta_decay.gnu')
      write(outfileindex, '(a)')'set terminal pdf enhanced color font ",24"'
      write(outfileindex, '(a)')"set output 'Delta_decay.pdf'"
      write(outfileindex, '(a)')"set style data linespoints"
      write(outfileindex, '(a)')"unset key"
      write(outfileindex, '(a)')"set pointsize 0.8"
      write(outfileindex, '(a)')"set border lw 2"
      write(outfileindex, '(a, i8, a)')'set xrange [1: ', Nslab, ']'
      write(outfileindex, '(a)')'set xlabel "Distance(layer)"'
      write(outfileindex, '(a)')'set ylabel "Delta(eV)"'
      write(outfileindex, '(a, f10.5, a)')'set yrange [0: ', Delta_BdG, ']'
      write(outfileindex, '(a)')"plot 'Delta_decay.dat' u 1:2  w lp lw 2 pt 7  ps 0.2 lc rgb 'red', 0 w l lw 2 dt 2"
      close(outfileindex)
   endif

   202 format('set xtics (',:20('"',A3,'" ',F10.5,','))
   203 format(A3,'" ',F10.5,')')
   204 format('set arrow from ',F10.5,',',F10.5, &
      ' to ',F10.5,',',F10.5, ' nohead')
 
   deallocate(eigenvalue_BdG)
   deallocate( surf_l_weight_BdG )
   deallocate( surf_l_weight_BdG_mpi )
   deallocate( surf_r_weight_BdG )
   deallocate( surf_r_weight_BdG_mpi )
   deallocate(ekslab_BdG)
   deallocate(ekslab_BdG_mpi)
   deallocate(CHamk_BdG)
   deallocate(work)
   deallocate(rwork)
 
return
end subroutine ek_slab_BdG

subroutine  wannier_center2D_BdG
   ! This suboutine is used for wannier center calculation for slab system
   ! 
   ! Copyright (c) 2010 QuanSheng Wu. All rights reserved.

    use para
    use wmpi
    implicit none

    integer :: Nkx
    integer :: Nky

    integer :: i, j, l, ia, ia1, m
    integer :: nfill

    integer :: ikx
    integer :: iky

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

    !> eigenvalue
    real(dp), allocatable :: eigenvalue(:)

    !> 2D surface BZ
    real(dp) :: kx
    real(dp) :: ky
    real(dp) :: k(2), b(2)
    real(dp) :: k0(2), k1(2), k2(2)

    !> b.R
    real(dp) :: br

    !> exp(-i*b.R)
    complex(dp) :: ratio

    real(dp) :: slab_Rua(2)
    real(dp) :: slab_Rub(2)
    real(dp) :: slab_Kua(2)
    real(dp) :: slab_Kub(2)
    real(dp) :: cell_slab

    !> for each orbital, it correspond to an atom
    !> dim= Num_wann_BdG
    integer, allocatable :: AtomIndex_orbital(:)

    !> atom position in the unit cell
    !> for slab BdG system, dim=Nslab*Origin_cell%Num_atoms
    real(dp), allocatable :: AtomsPosition_unitcell(:, :)
    real(dp), allocatable :: AtomsPosition_supercell(:,:)

    real(dp) :: Umatrix_t(3,3)


    Nkx= Nk1
    Nky= Nk2

    nfill= Num_wann*Nslab

    allocate(kpoints(2, Nkx, Nky))
    kpoints= 0d0

    allocate(Lambda_eig(nfill))
    allocate(Lambda(nfill, nfill))
    allocate(Lambda0(nfill, nfill))
    allocate(Mmnkb(nfill, nfill))
    allocate(Mmnkb_com(nfill, nfill))
    allocate(Mmnkb_full(Num_wann_BdG*Nslab, Num_wann_BdG*Nslab))
    allocate(hamk(Num_wann_BdG*Nslab, Num_wann_BdG*Nslab))
    allocate(hamk_dag(Num_wann_BdG*Nslab, Num_wann_BdG*Nslab))
    allocate(Eigenvector(Num_wann_BdG*Nslab, Num_wann_BdG*Nslab, Nkx))
    allocate(eigenvalue(Num_wann_BdG*Nslab))
    allocate(U(nfill, nfill))
    allocate(Sigma(nfill, nfill))
    allocate(VT(nfill, nfill))
    allocate(WannierCenterKy(nfill, Nky))
    allocate(WannierCenterKy_mpi(nfill, Nky))
    allocate(AtomIndex_orbital(Num_wann_BdG*Nslab))
    allocate(AtomsPosition_unitcell(3, Origin_cell%Num_atoms))
    allocate(AtomsPosition_supercell(3, Nslab*Origin_cell%Num_atoms))
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

    slab_Rua= 0.0d0
    slab_Rub= 0.0d0
    slab_Kua= 0.0d0
    slab_kub= 0.0d0
    cell_slab= 0.0d0

    slab_Rua= Rua_new(1:2)
    slab_Rub= Rub_new(1:2)
    cell_slab= slab_Rua(1)* slab_Rub(2)- slab_Rua(2)* slab_Rub(1)
    cell_slab= abs(cell_slab)

    if (abs(cell_slab)< 1e-6) stop "cell_volume equal to 0"

    slab_Kua(1)= 2d0*pi/cell_slab*slab_Rub(2)
    slab_Kua(2)=-2d0*pi/cell_slab*slab_Rub(1)
    slab_Kub(1)=-2d0*pi/cell_slab*slab_Rua(2)
    slab_Kub(2)= 2d0*pi/cell_slab*slab_Rua(1)

    if (cpuid==0) then
        write(stdout, *)'2D Primitive Cell_Volume: ', cell_slab/(Angstrom2atomic**2)
        write(stdout, *)'slab_Rua, slab_Rub'
        write(stdout, '(3f10.4)')slab_Rua/Angstrom2atomic
        write(stdout, '(3f10.4)')slab_Rub/Angstrom2atomic
        write(stdout, *)'slab_Kua, slab_Kub'
        write(stdout, '(3f10.4)')slab_Kua/Angstrom2atomic
        write(stdout, '(3f10.4)')slab_Kub/Angstrom2atomic
    endif

    b= 0.0d0
    b(1)= 1.d0/real(Nkx)
    b(2)= 0.d0
    b= b(1)*slab_Kua+ b(2)*slab_Kub

    k0= K2D_start
    k1= K2D_vec1
    k2= K2D_vec2
    do iky=1, Nky
       do ikx=1, Nkx
          kpoints(:, ikx, iky)= k0+ k1*dble(ikx-1.d0)/dble(Nkx)+ k2*(iky-1)/dble(Nky-1d0)
       enddo
    enddo

    !> set up atom index for each orbitals in the basis
    if (soc>0) then   !> with spin-orbit coupling
      l= 0
      do i=1, Nslab
          do ia=1, Origin_cell%Num_atoms    !> spin up
              do j=1, Origin_cell%nprojs(ia)
                  l= l+ 1
                  AtomIndex_orbital(l)= ia+ (i-1)*Origin_cell%Num_atoms  !> electron
                  AtomIndex_orbital(l+ Num_wann*Nslab)= ia+ (i-1)*Origin_cell%Num_atoms  !> hole
              enddo ! j
          enddo ! ia
          do ia=1, Origin_cell%Num_atoms    !> spin down
              do j=1, Origin_cell%nprojs(ia)
                  l= l+ 1
                  AtomIndex_orbital(l)= ia+ (i-1)*Origin_cell%Num_atoms  !> electron
                  AtomIndex_orbital(l+ Num_wann*Nslab)= ia+ (i-1)*Origin_cell%Num_atoms  !> hole
              enddo ! j
          enddo ! ia
      enddo ! i
    endif

    if (cpuid==0) then
        write(stdout, *)'AtomIndex_orbital: '
        write(stdout, *)AtomIndex_orbital
    endif

    Umatrix_t= transpose(Umatrix)
    call inv_r(3, Umatrix_t)

    !> set up atoms' position in the unit cell of the new basis
    !> only for 2D slab system
    AtomsPosition_unitcell=0.0d0

    !do ia=1, Origin_cell%Num_atoms
    !   do i=1, 3
    !      do j=1,3
    !          AtomsPosition_unitcell(i,ia)= AtomsPosition_unitcell(i, ia)+ &
    !              Umatrix_t(i,j)*Origin_cell%Atom_position_cart(j, ia)
    !      enddo ! j
    !   enddo ! i
    !enddo ! ia

    do ia=1, Origin_cell%Num_atoms
       do i=1, 3
          AtomsPosition_unitcell(i,ia)= Cell_defined_by_surface%Atom_position_cart(i,ia)
       enddo ! i
    enddo ! ia


    if (cpuid==0) then
        write(stdout, *)'AtomPosition_unitcell: '
        do ia= 1, Origin_cell%Num_atoms
           write(stdout, "(3f12.6)")AtomsPosition_unitcell(:,ia)
        enddo
    endif

    !> set up atoms' position in the supercell
    !> actually, we only need the first two corordinates: x, y
    AtomsPosition_supercell=0.0d0
    ia1= 0
    do i=1, Nslab
       do ia=1, Origin_cell%Num_atoms
          ia1= ia1+ 1
          AtomsPosition_supercell(:, ia1)= AtomsPosition_unitcell(:,ia) !+Ruc_new*(i-1d0)
       enddo ! ia
    enddo ! i

    if (cpuid==0) then
        write(stdout, *) 'AtomPosition_supercell: '
        do ia=1, Nslab*Origin_cell%Num_atoms
           write(stdout, "(3f12.6)")AtomsPosition_supercell(:,ia)
        enddo
    endif


    !> for each ky, we can get wanniercenter
    do iky=1+ cpuid, nky, num_cpu
       Lambda0=0d0
       do i=1, nfill
          Lambda0(i, i)= 1d0 ! lam0=I
       enddo

       if (cpuid==0) print *, iky, nky
       !> for each kx, we get the eigenvectors
       do ikx=1, nkx
          k(1)= kpoints(1, ikx, iky)
          k(2)= kpoints(2, ikx, iky)

          call ham_slab_BdG(k,Bz_surf,hamk)

          !> diagonal hamk
          call eigensystem_c('V', 'U', Num_wann_BdG*Nslab, hamk, eigenvalue)

          Eigenvector(:, :, ikx)= hamk
       enddo

       !> <u_k|u_k+1>
       !> sum over kx to get wanniercenters
       do ikx=1, nkx
          Mmnkb= 0d0
          hamk_dag= Eigenvector(:, :, ikx)
          if (ikx==nkx) then
             hamk= Eigenvector(:, :, 1)
          else
             hamk= Eigenvector(:, :, ikx+ 1)
          endif
          
          do l=1, Nslab*2
             do m=1, Num_wann
                ia= AtomIndex_orbital(m+ (l-1)*Num_wann)
                br= b(1)*AtomsPosition_supercell(1, ia)+ &
                    b(2)*AtomsPosition_supercell(2, ia)
                ratio= cos(br)- zi* sin(br)

                do i= 1, nfill
                   do j= 1, nfill
                      Mmnkb(i, j)= Mmnkb(i, j)+ &
                          conjg(hamk_dag((l-1)*Num_wann+m, i))* &
                          hamk((l-1)*Num_wann+m, j)* ratio
                   enddo ! j
                 enddo ! i
               enddo ! m
          enddo ! l

         !> <u_k|u_k+1>
         !call mat_mul(Num_wann_BdG*Nslab, hamk_dag, hamk, Mmnkb_full)
         !Mmnkb= Mmnkb_full(1:nfill, 1:nfill)

         !Mmnkb_com= 0d0
         !hamk_dag= Eigenvector(:, :, ikx)
         !hamk= Eigenvector(:, :, ikx+1)
         !do i=1, nfill
         !   do j=1, nfill
         !      do l= 1, Num_wann*Nslab
         !         Mmnkb_com(i, j)=  Mmnkb_com(i, j)+ conjg(hamk_dag(l, i))* hamk(l, j)
         !      enddo
         !   enddo
         !enddo

         !print *, maxval(real(Mmnkb-Mmnkb_com))
         !stop


          !> perform Singluar Value Decomposed of Mmnkb
          call zgesvd_pack(nfill, Mmnkb, U, Sigma, VT)

          !> after the calling of zgesvd_pack, Mmnkb becomes a temporal matrix
          U= conjg(transpose(U))
          VT= conjg(transpose(VT))
          call mat_mul(nfill, VT, U, Mmnkb)

         !> check hermicity
         !do i=1, nfill
         !   do j=i, nfill
         !      if (abs(Mmnkb(i, j)-conjg(Mmnkb(j, i)))>0.0001d0)then
         !         print *, 'Mmnkb is not Hermitian'
         !         print*, i, j, Mmnkb(i, j), Mmnkb(j, i)

         !      endif
         !   enddo
         !enddo

         !stop


          call mat_mul(nfill, Mmnkb, Lambda0, Lambda)
          Lambda0 = Lambda
       enddo  !< ikx

       !> diagonalize Lambda to get the eigenvalue 
       call zgeev_pack(nfill, Lambda, Lambda_eig)
       do i=1, nfill
          WannierCenterKy(i, iky)= aimag(log(Lambda_eig(i)))/2d0/pi
          WannierCenterKy(i, iky)= mod(WannierCenterKy(i, iky)+10d0, 1d0)
       enddo

    enddo !< iky

#if defined (MPI)
    call mpi_allreduce(WannierCenterKy, WannierCenterKy_mpi, &
         size(WannierCenterKy), mpi_dp, mpi_sum, mpi_cmw, ierr)
#else
    WannierCenterKy_mpi= WannierCenterKy
#endif


    outfileindex= outfileindex+ 1
    if (cpuid==0) then
       open(unit=outfileindex, file='wanniercenter_BdG.dat')
       do i=1, nfill
       do iky=1, Nky
          write(outfileindex, '(10000f16.8)') kpoints(2, 1, iky), &
             dmod(WannierCenterKy_mpi(i, iky), 1d0)
       enddo
       enddo
       close(outfileindex)
    endif


    outfileindex= outfileindex+ 1
    if (cpuid==0) then
       open(unit=outfileindex, file='wanniercenter_BdG_total.dat')
       do iky=1, Nky
          write(outfileindex, '(10000f16.8)') kpoints(2, 1, iky), &
             dmod(sum(WannierCenterKy_mpi(:, iky)), 1d0) 
       enddo
       close(outfileindex)
    endif

    outfileindex= outfileindex+ 1
    if (cpuid==0) then
       open(unit=outfileindex, file='wcc_slab_BdG.gnu')
 
       write(outfileindex,*) 'set encoding iso_8859_1'
       write(outfileindex,*) 'set terminal  postscript enhanced color font "Roman,36" '
       write(outfileindex,*) "set output 'wcc_slab_BdG.eps'"
       write(outfileindex,*) 'set size ratio -1 '
       write(outfileindex,*) 'set multiplot '
       write(outfileindex,*) 'unset key'
       write(outfileindex,*) 'set border lw 1 '
       write(outfileindex,*) 'set xtics 0.5 nomirror '
       write(outfileindex,*) 'set xtics ("k_y" 0, "-{/Symbol p}" -0.5, "{/Symbol p}" 0.5) '
       write(outfileindex,*) 'set ytics 0.5 nomirror '
       write(outfileindex,*) 'set xrange [-0.50: 0.5]'
       write(outfileindex,*) 'set yrange [-0.00: 1.0]'
       write(outfileindex,*) 'set ylabel "{/Symbol q}(2{/Symbol p})" rotate by  90 offset 2.8,0 '
       write(outfileindex,*) 'plot "wanniercenter_BdG.dat" u 1:2 w p  pt 7  ps 0.6 lc rgb "blue"'
       close(outfileindex)
    endif    

    outfileindex= outfileindex+ 1
    if (cpuid==0) then
       open(unit=outfileindex, file='wcc_slab_BdG_total.gnu')

       write(outfileindex,*) 'set encoding iso_8859_1'
       write(outfileindex,*) 'set terminal  postscript enhanced color font "Roman,36" '
       write(outfileindex,*) "set output 'wcc_slab_BdG_total.eps'"
       write(outfileindex,*) 'set size ratio -1 '
       write(outfileindex,*) 'set multiplot '
       write(outfileindex,*) 'unset key'
       write(outfileindex,*) 'set border lw 1 '
       write(outfileindex,*) 'set xtics 0.5 nomirror '
       write(outfileindex,*) 'set xtics ("k_y" 0, "-{/Symbol p}" -0.5, "{/Symbol p}" 0.5) '
       write(outfileindex,*) 'set ytics 0.5 nomirror '
       write(outfileindex,*) 'set xrange [-0.50: 0.5]'
       write(outfileindex,*) 'set yrange [-0.00: 1.0]'
       write(outfileindex,*) 'set ylabel "{/Symbol q}(2{/Symbol p})" rotate by 90 offset 2.8,0 '
       write(outfileindex,*) 'plot "wanniercenter_BdG_total.dat" u 1:2 w p  pt 7  ps 0.6 lc rgb "blue"'
       close(outfileindex)
    endif

    return
 end subroutine  wannier_center2D_BdG

 subroutine ek_slab_BdG_phase
  
   use wmpi
   use para
   implicit none 

   ! loop index
   integer :: i, j, l, lwork, ierr

   real(Dp) :: k(2), emin, emax

   real(Dp) :: Bz_BdG_phase(Nk1)
   real(Dp) :: Bz_BdG_i

   ! time measurement
   real(dp) :: time_start, time_end, time_start0

   ! parameters for zheev
   real(Dp), allocatable ::  rwork(:)
   complex(Dp), allocatable :: work(:)

   ! eigenvalue 
   real(Dp), allocatable :: eigenvalue_BdG(:)

   ! energy dispersion
   real(Dp), allocatable :: ekslab_BdG(:,:), ekslab_BdG_mpi(:,:)

   ! hamiltonian slab
   complex(Dp), allocatable :: CHamk_BdG(:,:)

   lwork= 16*Nslab*Num_wann_BdG
   ierr = 0

   allocate(eigenvalue_BdG(nslab* Num_wann_BdG))
   allocate(ekslab_BdG(nslab* Num_wann_BdG, Nk1))
   allocate(ekslab_BdG_mpi(nslab* Num_wann_BdG, Nk1))
   allocate(CHamk_BdG(nslab* Num_wann_BdG, nslab* Num_wann_BdG))

   allocate(work(lwork))
   allocate(rwork(lwork))

   ! sweep k
   ekslab_BdG=0.0d0
   ekslab_BdG_mpi=0.0d0

   time_start= 0d0
   time_start0= 0d0
   call now(time_start0)
   time_start= time_start0
   time_end  = time_start0

   Bz_BdG_phase=0d0
   do i=1, Nk1
      Bz_BdG_phase(i)= real(i-1d0)/real(Nk1-1d0)*Bz_surf
   enddo

   do i= 1+cpuid, Nk1, num_cpu
      if (cpuid==0.and. mod(i/num_cpu, 4)==0) &
         write(stdout, '(a, i9, "  /", i10, a, f10.1, "s", a, f10.1, "s")') &
         ' Slabek: ik', i, Nk1, ' time left', &
         (Nk1-i)*(time_end- time_start)/num_cpu, &
         ' time elapsed: ', time_end-time_start0 

      call now(time_start)

      k= Single_KPOINT_2D_DIRECT

      Bz_BdG_i= Bz_BdG_phase(i)
      
      CHamk_BdG=0.0d0
      call ham_slab_BdG(k, Bz_BdG_i, CHamk_BdG)

      eigenvalue_BdG=0.0d0

      ! diagonal Chamk_BdG
      call eigensystem_c('V', 'U', Num_wann_BdG*Nslab, CHamk_BdG, eigenvalue_BdG)
     
      ekslab_BdG(:,i)=eigenvalue_BdG

      call now(time_end)
   enddo ! i

#if defined (MPI)
   call mpi_allreduce(ekslab_BdG,ekslab_BdG_mpi,size(ekslab_BdG),&
                     mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
   ekslab_BdG_mpi= ekslab_BdG
#endif

   ekslab_BdG_mpi= ekslab_BdG_mpi/eV2Hartree

   ekslab_BdG=ekslab_BdG_mpi
   
   outfileindex= outfileindex+ 1
   if(cpuid==0)then
      open(unit=outfileindex, file='slabek_BdG_phase.dat')
      write(outfileindex, "('#', a10, a15, 5X)")'# k', ' E'
      do j=1, Num_wann_BdG*Nslab
         do i=1, Nk1
            write(outfileindex,'(2f15.7)')Bz_BdG_phase(i), ekslab_BdG(j,i)
         enddo
         write(outfileindex , *)''
      enddo
      close(outfileindex)
      write(stdout,*) 'calculate BdG_phase done'
   endif


   emin= -Delta_BdG*2
   emax= Delta_BdG*2
   !> write script for gnuplot
   outfileindex= outfileindex+ 1
   if (cpuid==0) then
      open(unit=outfileindex, file='slabek_BdG_phase.gnu')
      write(outfileindex, '(a)')'set terminal pdf enhanced color font ",24"'
      write(outfileindex, '(a)')"set output 'slabek_BdG_phase.pdf'"
      write(outfileindex, '(a)')"set style data linespoints"
      write(outfileindex, '(a)')"unset key"
      write(outfileindex, '(a)')"set pointsize 0.8"
      write(outfileindex, '(a)')"set border lw 2"
      write(outfileindex, '(a, f10.5, a)')'set xrange [0: ', Bz_surf, ']'
      write(outfileindex, '(a)')'set xlabel "Bz surf(eV)"'
      write(outfileindex, '(a)')'set ylabel "Energy(eV)"'
      write(outfileindex, '(a, f10.5, a, f10.5, a)')'set yrange [', emin, ':', emax, ']'
      write(outfileindex, '(a)')"plot 'slabek_BdG_phase.dat' u 1:2  w lp lw 2 pt 7  ps 0.2 lc rgb 'red', 0 w l lw 2 dt 2"
      close(outfileindex)
   endif


   202 format('set xtics (',:20('"',A3,'" ',F10.5,','))
   203 format(A3,'" ',F10.5,')')
   204 format('set arrow from ',F10.5,',',F10.5, &
      ' to ',F10.5,',',F10.5, ' nohead')
 

   deallocate(eigenvalue_BdG)
   deallocate(ekslab_BdG)
   deallocate(ekslab_BdG_mpi)
   deallocate(CHamk_BdG)

   deallocate(work)
   deallocate(rwork)
 
return
end subroutine ek_slab_BdG_phase
