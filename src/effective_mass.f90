   subroutine effective_mass_calc()
      !> effective_mass_calc calculates effective mass for a
      !> specified  band index in the framework of Wannier tight binding
      !
      !> QuanSheng Wu (wuquansheng@gmail.com)
      !
      !> Aug 30 2016 @ ETH Zurich
      !
      !> refs: 
      !
      !> The effective mass calculation with VASP has been done with Alexandr Fonari
      !
      !> see https://afonari.com/emc.html
      !
      !> Abramowitz, M.; Stegun, I. Handbook of Mathematical Functions: with Formulas, Graphs, and
      !> Mathematical Tables; Dover, 1972.
      !
      !> here we implement the TB version
      use wmpi
      use para
      implicit none

      !> k points for effective mass calculation
      integer :: Nk_mass
      integer :: ik1, ik2, ik3, i, j

      real(dp) :: k1(3)
      real(dp) :: k(3)

      !> energy level for kpoints at the given band index
      real(dp), allocatable :: W(:)
      real(dp), allocatable :: eigval(:, :, :)
      complex(dp), allocatable :: Hamk_bulk(:, :)

      real(dp), allocatable :: Emass(:, :)

      allocate(W(Num_wann))
      allocate(Emass(3, 3))
      allocate(eigval(-2:2, -2:2, -2:2))
      allocate(Hamk_bulk(Num_wann, Num_wann))
      W= 0d0
      eigval= 0d0
      Hamk_bulk= 0d0

   
      !> get the eigenvalues
      do ik1=-2, 2
      do ik2=-2, 2
      do ik3=-2, 2
         k1(1)= ik1*dk_mass+ k_mass(1)
         k1(2)= ik2*dk_mass+ k_mass(2)
         k1(3)= ik3*dk_mass+ k_mass(3)
         call cart_direct_rec(k1, k)
         call ham_bulk_latticegauge(k, Hamk_bulk)

         W= 0d0
         call eigensystem_c( 'V', 'U', Num_wann ,Hamk_bulk, W)
         eigval(ik1, ik2, ik3)= W(iband_mass)
      enddo !ik3
      enddo !ik2
      enddo !ik1

      !> use atomic unit
      eigval = eigval*eV2Hartree
     
      !> xx
      Emass(1, 1) = (-eigval(-2, 0, 0) + 16.0d0*eigval(-1, 0, 0) - 30.0d0*eigval(0,0,0)&
         &+16.0d0*eigval(1, 0, 0) - eigval(2, 0, 0))/(12.0d0*dk_mass*dk_mass)

      !> yy
      Emass(2, 2) = (-eigval(0, -2, 0) + 16.0d0*eigval(0, -1, 0) - 30.0d0*eigval(0,0,0)&
         &+16.0d0*eigval(0, 1, 0) - eigval(0, 2, 0))/(12.0d0*dk_mass*dk_mass)

      !> zz
      Emass(3, 3) = (-eigval(0, 0, -2) + 16.0d0*eigval(0, 0, -1) - 30.0d0*eigval(0,0,0)&
         &+16.0d0*eigval(0, 0, 1) - eigval(0, 0, 2))/(12.0d0*dk_mass*dk_mass)

      Emass(2, 1)=&
         &(-63.0d0*(eigval( 1,-2,0)+ eigval( 2,-1, 0)+ eigval(-2, 1,0)+ eigval(-1,2,0))&
         &+ 63.0d0*(eigval(-1,-2,0)+ eigval(-2,-1, 0)+ eigval( 1, 2,0)+ eigval( 2,1,0))&
         &+ 44.0d0*(eigval( 2,-2,0)+ eigval(-2, 2, 0)- eigval(-2,-2,0)- eigval( 2,2,0))&
         &+ 74.0d0*(eigval(-1,-1,0)+ eigval( 1, 1, 0)- eigval( 1,-1,0)- eigval(-1,1,0)))&
         &/(600.0d0*dk_mass*dk_mass)
 
      Emass(3, 1)=&
         &(-63.0d0*(eigval( 1,0,-2)+ eigval( 2, 0,-1)+ eigval(-2,0, 1)+ eigval(-1,0,2))&
         &+ 63.0d0*(eigval(-1,0,-2)+ eigval(-2, 0,-1)+ eigval( 1,0, 2)+ eigval( 2,0,1))&
         &+ 44.0d0*(eigval( 2,0,-2)+ eigval(-2, 0, 2)- eigval(-2,0,-2)- eigval( 2,0,2))&
         &+ 74.0d0*(eigval(-1,0,-1)+ eigval( 1, 0, 1)- eigval( 1,0,-1)- eigval(-1,0,1)))&
         &/(600.0d0*dk_mass*dk_mass)
      !yz
      Emass(3, 2)=&
         &(-63.0D0*(eigval(0, 1,-2)+ eigval(0, 2,-1)+ eigval(0,-2, 1)+ eigval(0,-1,2))&
         &+ 63.0D0*(eigval(0,-1,-2)+ eigval(0,-2,-1)+ eigval(0, 1, 2)+ eigval(0, 2,1))&
         &+ 44.0D0*(eigval(0, 2,-2)+ eigval(0,-2, 2)- eigval(0,-2,-2)- eigval(0, 2,2))&
         &+ 74.0D0*(eigval(0,-1,-1)+ eigval(0, 1, 1)- eigval(0, 1,-1)- eigval(0,-1,1)))&
         &/(600.0D0*dk_mass*dk_mass)

      Emass(1, 2)= Emass(2, 1)
      Emass(1, 3)= Emass(3, 1)
      Emass(2, 3)= Emass(3, 2)

      if (cpuid==0) write(stdout,*) ">> Inverse of effective mass tensor:"
      if (cpuid==0) write(stdout,"(3F15.8)") ((Emass(i,j), j=1,3),i=1,3)
      if (cpuid==0) write(stdout,*)

      call inv_r(3, Emass)
      Emass= Emass/Angstrom2atomic/Angstrom2atomic
      if (cpuid==0) write(stdout,'(a, i5, a)') ">> The effective mass tensor for ", iband_mass, " 'th band"
      if (cpuid==0) write(stdout,*) ">> Effective mass tensor (in unit of bare electron mass):"
      if (cpuid==0) write(stdout,"(3F15.8)") ((Emass(i,j), j=1,3),i=1,3)
      if (cpuid==0) write(stdout,*)

      deallocate(W)
      deallocate(Emass)
      deallocate(eigval)
      deallocate(Hamk_bulk)
      return
   end subroutine effective_mass_calc

