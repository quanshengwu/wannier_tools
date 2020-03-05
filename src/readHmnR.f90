subroutine readNormalHmnR()
   !>> Read in the tight-binding model from wannier90_hr.dat
   !>  The format is defined by the wannier90 software
   ! Constructed by quansheng wu 4/2/2010
   !
   ! Yifei Guan added the sparse hr file parsing June/2018
   ! License: GPL V3

   use para

   implicit none

   character*4 :: c_temp

   ! file existence
   logical :: exists

   integer :: i, j, ir, ia, ib, ic, irp, iap, ibp, icp, io
   integer :: i1, i2, i3, i4, i5
   integer :: n, m, ia1, ia2, ir0
   integer :: add_electric_field
   integer :: rs, re, cs, ce, nwann, nwann_nsoc

   real(dp) :: r1, r2, static_potential
   real(dp) :: tot, rh, ih, dis
   real(dp) :: new_ia, new_ib, new_ic

   !> for newcell
   real(dp) :: vrpt(3),apos1d(3),apos2d(3), pos(Origin_cell%Num_atoms)

   !>count newcell nrpts
   integer,parameter :: mxir=20
   integer, allocatable  :: nzir(:, :, :)
   integer :: nir1,nir2,nir3
   integer :: nrpts_new

   !> all new irs
   integer, allocatable :: allirs(:, :, :, :)
   integer :: irn1(3),irn2(3)

   allocate( nzir(-mxir:mxir,-mxir:mxir,-mxir:mxir))
   allocate( allirs(nrpts, num_wann, num_wann, 3))


   !> add a check for NumOccupied parameter
   if (NumOccupied<=0 .or. NumOccupied>Num_wann) then
      write(stdout, '(a, i6, a)')">>> ERROR: NumOccupied should be in [1, ", Num_wann, " ]"
      write(stdout, '(a)')">>> Usually, it is the number of occupied Wannier bands."
      stop
   endif

   if(cpuid.eq.0)write(stdout,*)' '
   open(12, file=Hrfile, status='OLD')

   if (index(Hrfile, 'HWR')==0) then
      !> for Normal HmnR obtained from Wannier90 or sparse HmnR

      !> skip a comment line
      read(12, *)
   
      !> number of Wannier orbitals in the hr file
      nwann=0
      read(12, *)nwann
      if (nwann==0) then
         stop "ERROR : num_wann is zero in hr file"
      endif
      nwann_nsoc=nwann
      if (SOC>0) nwann_nsoc= nwann/2
   
      !> number of lattice vectors taken into account
      read(12, *)nrpts
   
      !> The degeneracy of each R point, this is caused by the Fourier transform in the Wigner-Seitz cell
      read(12, *)(ndegen(i), i=1, nrpts)
      ir=0
      do ir=1,nrpts
         do n=1,nwann
            do m=1,nwann
               read(12,*,end=1001)i1, i2, i3, i4, i5, r1, r2
               irvec(1,ir)=i1
               irvec(2,ir)=i2
               irvec(3,ir)=i3
               HmnR(i4,i5,ir)=dcmplx(r1,r2)
            end do
         enddo
      enddo

   else

      !File *.HWR exist, We are using HmnR from WHM
      ! skip 8 lines
      do i=1,8
         read(12,*)
      enddo
      read(12,'(a11,f18.7)')c_temp,E_fermi
      do iR=1,Nrpts
         read(12,'(a3,3i5,a3,i4)')c_temp,irvec(1:3,iR),c_temp,ndegen(iR)
         do i=1, Num_wann*Num_wann
            read(12,*)n,m,rh,ih
            HmnR(n,m,iR)=rh+ zi*ih   ! in Hartree
         enddo
         if (sum(abs(irvec(:, ir)))==0) then
            do i=1, Num_wann
               HmnR(i,i,iR)= HmnR(i,i,iR)-E_fermi  ! in Hartree
            enddo
         endif
      enddo
      HmnR= HmnR*27.2114d0 ! from Hartree to eV

      if (cpuid==0) then
         open(unit=105, file='wannier90_hr.dat')
         write(105, *)'hr file transformed from HWR'
         write(105, *)Nwann
         write(105, *)nrpts
         write(105, '(15I5)')(ndegen(i), i=1, nrpts)
         do ir=1, nrpts
            do i=1, Nwann
               do j=1, Nwann
                  write(105, '(5I5, 2f16.8)')irvec(:, ir), i, j, HmnR(i, j, ir)
               enddo
            enddo
         enddo
         close(105)
      endif


   endif ! HWR or not
   1001 continue
   close(12)

   !call get_fermilevel
   !check sum rule
   tot= 0d0
   do ir=1, Nrpts
      tot= tot+ 1d0/ndegen(ir)
   enddo

   !> get the Cartesian coordinates for R points
   allocate( crvec(3, nrpts))

   !> get R coordinates
   do iR=1, Nrpts
      crvec(:, iR)= Origin_cell%Rua*irvec(1,iR) + Origin_cell%Rub*irvec(2,iR) + Origin_cell%Ruc*irvec(3,iR)
   enddo

   !>> hamiltonian for the new cell
   allocate(HmnR_newcell(Num_wann, Num_wann, Nrpts))
   allocate(irvec_newcell(3, Nrpts))
   allocate(ndegen_newcell(Nrpts))
   HmnR_newcell= HmnR
   ndegen_newcell= ndegen

   !> extract the fermi energy
   do iR=1,Nrpts
      if (Irvec(1,iR).eq.0.and.Irvec(2,iR).eq.0.and.Irvec(3,iR).eq.0)then
         do i=1, Num_wann
            HmnR(i,i,iR)=HmnR(i,i,iR)-E_fermi
         enddo
      endif
   enddo

   !> Adding zeeman field
   !> Bx=Bdirection(1)
   !> By=Bdirection(2)
   !> Bz=Bdirection(3)
   !> Hz= Zeeman_energy_in_eV*(Bx*sx+By*sy+Bz*sz)/2d0
   !> sx, sy, sz are Pauli matrices.
   if (Add_Zeeman_Field) then
      if (abs(Zeeman_energy_in_eV)<eps6)then
         Zeeman_energy_in_eV= Bmagnitude*Effective_gfactor*Bohr_magneton
         if (cpuid==0)then
            write(stdout, '(1x, a, 3f16.6)')"Zeeman_energy_in_eV: ",  Zeeman_energy_in_eV
         endif
      endif

      !> After considering the Zeeman field, we already extended the spin space to spin-full.
      SOC = 1
   endif ! Add_Zeeman_Field
 
   call get_stacking_direction_and_pos(add_electric_field, pos)
   if (add_electric_field>0) then
      ir0=0
      do ir=1, nrpts
         if (irvec(1, ir)==0.and.irvec(2, ir)==0.and.irvec(3, ir)==0) ir0=ir
      enddo
      io=0
      do ia=1, Origin_cell%Num_atoms
         static_potential= pos(ia)*Origin_cell%cell_parameters(add_electric_field)*Electric_field_in_eVpA
         if (Inner_symmetrical_Electric_Field) then
            static_potential= abs(pos(ia)*Origin_cell%cell_parameters(add_electric_field))*Electric_field_in_eVpA
         endif
         do i=1, Origin_cell%nprojs(ia)
            io=io+1
            HmnR(io, io, ir0)= HmnR(io, io, ir0)+ static_potential
            if (SOC>0) then
               HmnR(io+Num_wann/2, io+Num_wann/2, ir0)= HmnR(io+Num_wann/2, io+Num_wann/2, ir0)+ static_potential
            endif ! SOC
         enddo ! nproj
      enddo ! ia
   endif  ! add electric field or not


   return
end subroutine readNormalHmnR

!>-----------------------------------------------------------------------
! The sparse hamiltonian is stored as follows
! start with a comment line
! 10  ! splen: number of non-zero Hmn(R) lines 
! 20  ! nwan:  number of orbitals (including spin-degeneracy if SOC is included) 
! 9   ! nrpts: number of R points
! 1    0    0    1    1    1.0    0.0  ! i1, i2, i3, n, m, rh, ih  Hnm(i1, i2, i3)
! ...  ! in total there are splen lines
!> Here the unit of energy is eV
!>-----------------------------------------------------------------------
subroutine readSparseHmnR
   !> This subroutine not just read the sparse hr file, but also can read
   !> the standard hr file defined in the Wannier90.
   use para
   implicit none
   integer:: i,j,nwann,nwann_nsoc,i1,i2,i3,i4,i5,ir, ia, io
   integer :: q,n,m, ir0, irb
   real(8) :: r1,r2, pos(Origin_cell%Num_atoms), dxmax(3), center, static_potential

   !> the direction which adding electric field which is also the stacking direction
   integer :: add_electric_field
   logical :: lfound
   real(dp) :: Bx_in_Tesla, By_in_Tesla, Bz_in_Tesla
   complex(dp) :: h_value

   open(12, file=Hrfile)

   !> skip a comment line
   read(12, *)

   !> comparing with the standard hr file, we add another line to show howmany
   !> lines that Hmn(R) is not zero.
   if(Is_Sparse_Hr) then
      read(12,*) splen_input   !> number of non-zeros lines
   end if

   !> number of Wannier orbitals in the hr file
   nwann=0
   read(12, *)nwann
   if (nwann==0) then
      stop "ERROR : num_wann is zero in hr file"
   endif
   nwann_nsoc=nwann
   if (SOC>0) nwann_nsoc= nwann/2

   !> number of lattice vectors taken into account
   read(12, *)nrpts

   !> The degeneracy of each R point
   read(12, *)(ndegen(i), i=1, nrpts)
   ir=0
  
   !> whether we need to add the electric field potential
   add_electric_field=0
   call get_stacking_direction_and_pos(add_electric_field, pos)

   if (add_electric_field>0) then
      !> with Electric_field
      if (Add_Zeeman_Field) then
         if (SOC==0) then
            splen= splen_input*2+ 6*nwann_nsoc
         else
            splen= splen_input+ 6*nwann_nsoc
         endif
      else
         if (SOC==0) then
            splen= splen_input+ nwann_nsoc
         else
            splen= splen_input+ nwann_nsoc*2
         endif
      endif
   else
      !> without electric field
      if (Add_Zeeman_Field) then
         if (SOC==0) then
            splen= splen_input*2+ 4*nwann_nsoc
         else
            splen= splen_input+ 4*nwann_nsoc
         endif
      else
         splen= splen_input
      endif
   end if

   allocate(hacoo(splen),hicoo(splen),hjcoo(splen),hirv(splen))
   hacoo=(0d0, 0d0)
   hicoo=0
   hjcoo=0
   hirv=0  

   j=0
   ir0=0
   ir=0
   do q=1,nrpts
      do n=1,nwann
         do m=1,nwann
            read(12,*,end=1001)i1, i2, i3, i4, i5, r1, r2
            if(ir==0) then
               ir=ir+1
               irvec(1,ir)=i1
               irvec(2,ir)=i2
               irvec(3,ir)=i3
            else
               i=abs(i1-irvec(1,ir))+abs(i2-irvec(2,ir))+abs(i3-irvec(3,ir))
               if(i/=0) ir=ir+1
               2011 continue
               irvec(1,ir)=i1
               irvec(2,ir)=i2
               irvec(3,ir)=i3
            endif
            j=j+1
            hicoo(j)=i4
            hjcoo(j)=i5
            hirv (j)=ir
            hacoo(j)=dcmplx(r1,r2)
            if (i1==0.and.i2==0.and.i3==0.and.i4==i5) then
               hacoo(j)= hacoo(j)- E_FERMI
               ir0= ir
            endif
         end do
      enddo
   enddo


   1001 continue 

   !> Adding zeeman field
   !> Bx=Bdirection(1)
   !> By=Bdirection(2)
   !> Bz=Bdirection(3)
   !> Hz= Zeeman_energy_in_eV*(Bx*sx+By*sy+Bz*sz)/2d0
   !> sx, sy, sz are Pauli matrices.
   if (Add_Zeeman_Field) then
      if (abs(Zeeman_energy_in_eV)<eps6)then
         Zeeman_energy_in_eV= Bmagnitude*Effective_gfactor*Bohr_magneton
         if (cpuid==0)then
            write(stdout, '(1x, a, 3f16.6)')"Zeeman_energy_in_eV: ",  Zeeman_energy_in_eV
         endif
      endif

      !> for sparse hr format
      !> extend from spinless to spinfull
      if (SOC_in==0) then
         do i=1, splen_input
            j= i+ splen_input
            hicoo(j)=hicoo(i)+nwann_nsoc
            hjcoo(j)=hjcoo(i)+nwann_nsoc
            hirv (j)=hirv (i)
            hacoo(j)=hacoo(i)
         enddo
      endif
   
      !> zeeman energy
      if (SOC_in==0) then
         j= splen_input*2
      else
         j= splen_input
      endif

      Bx_in_Tesla= Bx
      By_in_Tesla= By
      Bz_in_Tesla= Bz
      call add_zeeman_sparse_hr(Bx_in_Tesla, By_in_Tesla, Bz_in_Tesla)
      j=j+ nwann_nsoc*4
      !> After adding zeeman term, SOC=1
      SOC = 1
   endif ! Add_Zeeman_Field

   !> Adding electric field
   if (add_electric_field>0) then
      if (Add_Zeeman_Field) then
         !> continue with the j used in the Add_Zeeman_Field
         io=0
         do ia=1, Origin_cell%Num_atoms
            static_potential= pos(ia)*Origin_cell%cell_parameters(add_electric_field)*Electric_field_in_eVpA
            if (Inner_symmetrical_Electric_Field) then
               static_potential= abs(pos(ia)*Origin_cell%cell_parameters(add_electric_field))*Electric_field_in_eVpA
            endif
            do i=1, Origin_cell%nprojs(ia)
               io=io+1
               !> spin up
               j=j+1
               hicoo(j)= io 
               hjcoo(j)= io 
               hirv (j)= ir0
               hacoo(j)= static_potential
               !> spin down
               j=j+1
               hicoo(j)= io + nwann_nsoc
               hjcoo(j)= io + nwann_nsoc
               hirv (j)= ir0
               hacoo(j)= static_potential
            enddo ! nproj
         enddo ! ia
      else
         j=splen_input
         io=0
         do ia=1, Origin_cell%Num_atoms
            static_potential= pos(ia)*Origin_cell%cell_parameters(add_electric_field)*Electric_field_in_eVpA
            if (Inner_symmetrical_Electric_Field) then
               static_potential= abs(pos(ia)*Origin_cell%cell_parameters(add_electric_field))*Electric_field_in_eVpA
            endif
            do i=1, Origin_cell%nprojs(ia)
               io=io+1
               j=j+1
               hicoo(j)= io 
               hjcoo(j)= io 
               hirv (j)= ir0
               hacoo(j)= static_potential
               if (SOC>0) then  ! with SOC
                  j=j+1
                  hicoo(j)= io + nwann_nsoc
                  hjcoo(j)= io + nwann_nsoc
                  hirv (j)= ir0
                  hacoo(j)= static_potential
               endif
            enddo ! nproj
         enddo ! ia
      endif ! Add_Zeeman_Field or not
   endif ! add_electric_field
 
   !> transform it into dense format only when number of wannier orbitals is less than 100

   outfileindex= outfileindex+ 1
   if (cpuid.eq.0.and. nwann<100) then
      open (unit=outfileindex, file='hr.dat-dense')
      write(outfileindex, *) ' ! HmnR file from sparse hr file'
      write(outfileindex, '(I10, a)') nwann, '  ! Num_wann: number of wannier orbitals'
      write(outfileindex, '(I10, a)') nrpts, '  ! NRPTS: number of R points'
      write(outfileindex, '(10I5)') ndegen(:) 
      do iR=1, NRPTS
         do n=1, nwann
            do m=1, nwann
               lfound=.false.
               h_value=0d0
               do i=1, splen
                  if (hicoo(i)==n.and.hjcoo(i)==m.and.hirv(i)==iR)then
                     h_value= h_value+ hacoo(i)
                  endif
               enddo !i
               write(outfileindex, '(5I5, 2f12.6)')irvec(:, iR), n, m, h_value
            enddo !m
         enddo !n
      enddo !iR
      close(outfileindex)
   endif




   return
end subroutine readSparseHmnR


!> return the maximum distance between two neighbour elements of array x
!> dx_max= mod(dx_max, 1)
!> x is in [0, 1)
subroutine get_max_nn_distance(n, x, dx_max)
   use para, only: dp
   implicit none

   !> dimension of arr
   integer, intent(in) :: n

   !> an array
   real(dp), intent(in) :: x(n)
   real(dp) :: x1(n)

   real(dp), intent(out) :: dx_max

   integer :: i
   real(dp) :: dx

   dx_max=0d0
   do i=1, n-1
      dx= x1(i+1)-x1(i)
      if (dx>dx_max) dx_max= dx
   enddo
   dx= 1d0+x(1)-x(n)
   if (dx>dx_max) dx_max= dx

   return
end subroutine get_max_nn_distance


!> output: add_electric_field is the direction which the E field is adding on.
!> if add_electric_field=0, we don't add E field.
!> pos are the direct coordinate of atoms along the direction where the Electric_field applied.
subroutine get_stacking_direction_and_pos(add_electric_field, pos)
   use para
   implicit none

   integer, intent(inout) :: add_electric_field
   real(dp), intent(out) :: pos(Origin_cell%Num_atoms)

   integer :: i, ia
   real(dp) :: dxmax(3), center

   !> Adding electric field 
   !> The electric field is only applied on the 2D system.
   !> First we determine whether we have an vacuum and which direction
   !> Here we only support adding electric field along the unit lattice vectors.
   add_electric_field=0
   pos=0d0
   if (abs(Electric_field_in_eVpA)>eps6) then
      do i=1, 3
         pos= Origin_cell%Atom_position_direct(i, :)
         !> shift all the values into [0, 1)
         do ia=1, Origin_cell%Num_atoms
            pos(ia)=mod(pos(ia), 1d0)
         enddo

         call sortheap(Origin_cell%Num_atoms, pos)
         call get_max_nn_distance(Origin_cell%Num_atoms, pos, dxmax(i))
         dxmax(i)= dxmax(i)*Origin_cell%cell_parameters(i)
         if (dxmax(i)>9d0) then
            add_electric_field=i
         endif
      enddo
   endif ! add electric field or not

   if (add_electric_field==0.and.cpuid==0) then
      write(stdout, '(a)') " "
      write(stdout, '(a)') "  Note: This is not a 2D system, so we can't add electric field."
      write(stdout, '(a)') " "
   endif

   !> shift all the positions together and centered at zero along the electric field
   if (add_electric_field>0) then
      pos=Origin_cell%Atom_position_direct(add_electric_field, :)
      pos= mod(pos, 1d0)-0.5
      center= (maxval(pos)+minval(pos))/2d0
      pos=pos-center
   endif  ! add electric field or not

   return
end subroutine get_stacking_direction_and_pos


subroutine add_zeeman_sparse_hr(Bx_in_Tesla, By_in_Tesla, Bz_in_Tesla)
   !> add Zeeman energy on the sparse hmnr based on the magnetic field strength
   !> magnetic_field_in_tesla is in Tesla
   !> Here the coordinate system is the same as user-specified in the LATTICE card.
   use para
   implicit none

   real(dp), intent(in) :: Bx_in_Tesla, By_in_Tesla, Bz_in_Tesla 

   integer :: nwann_nsoc, i, j, ir0, ir
   real(dp) :: Zeeman_energy_in_eV_factor

   nwann_nsoc= Num_wann/2

   if (.not.Add_Zeeman_Field) return

   ir0=0
   do ir=1, nrpts
      if (irvec(1,ir)==0.and.irvec(2,ir)==0.and.irvec(3,ir)==0) then
         ir0= ir
         exit
      endif
   enddo
   if (ir0==0) stop 'something wrong with irvec in subroutine add_zeeman_sparse_hr'

   Zeeman_energy_in_eV_factor= Effective_gfactor*Bohr_magneton
   !> zeeman energy
   if (SOC_in==0) then
      j= splen_input*2
   else
      j= splen_input
   endif

   do i=1, nwann_nsoc
      j= j+ 1
      hicoo(j)= i
      hjcoo(j)= i
      hirv (j)= ir0
      hacoo(j)= Zeeman_energy_in_eV_factor/2d0*Bz_in_Tesla
   enddo
   do i=1, nwann_nsoc
      j= j+ 1
      hicoo(j)= i+ nwann_nsoc
      hjcoo(j)= i+ nwann_nsoc
      hirv (j)= ir0
      hacoo(j)= -Zeeman_energy_in_eV_factor/2d0*Bz_in_Tesla
   enddo
   do i=1, nwann_nsoc
      j= j+ 1
      hicoo(j)= i
      hjcoo(j)= i+ nwann_nsoc
      hirv (j)= ir0
      hacoo(j)=  Zeeman_energy_in_eV_factor/2d0*(By_in_Tesla*zi+ Bx_in_Tesla)
   enddo
   do i=1, nwann_nsoc
      j= j+ 1
      hicoo(j)= i+ nwann_nsoc
      hjcoo(j)= i
      hirv (j)= ir0
      hacoo(j)=  Zeeman_energy_in_eV_factor/2d0*(-By_in_Tesla*zi+ Bx_in_Tesla)
   enddo
   if (j>splen) stop "ERROR happend in setting zeeman energy in readhmnr.f90"

   return

end subroutine add_zeeman_sparse_hr



!subroutine add_zeeman_normal_hr(Bx_in_Tesla, By_in_Tesla, Bz_in_Tesla)
!  !> add Zeeman energy on the sparse hmnr based on the magnetic field strength
!  !> magnetic_field_in_tesla is in Tesla
!  !> Here the coordinate system is the same as user-specified in the LATTICE card.
!  use para
!  implicit none

!  real(dp), intent(in) :: Bx_in_Tesla, By_in_Tesla, Bz_in_Tesla 

!  integer :: nwann_nsoc, i, j, ir0, ir
!  real(dp) :: Zeeman_energy_in_eV_factor

!  nwann_nsoc= Num_wann/2

!  if (.not.Add_Zeeman_Field) return

!  ir0=0
!  do ir=1, nrpts
!     if (irvec(1,ir)==0.and.irvec(2,ir)==0.and.irvec(3,ir)==0) then
!        ir0= ir
!     endif
!  enddo
!  if (ir0==0) stop 'something wrong with irvec in subroutine add_zeeman_sparse_hr'

!  Zeeman_energy_in_eV_factor= Effective_gfactor*Bohr_magneton

!  if (SOC==0) then
!     do ir=1, nrpts
!        do n=1, nwann_nsoc
!           do m=1, nwann_nsoc
!              if (irvec(1, ir)==0.and.irvec(2, ir)==0.and.irvec(3, ir)==0.and.n==m) then
!                 HmnR(n+nwann_nsoc, m+nwann_nsoc, ir)= HmnR(n, m, ir)- Zeeman_energy_in_eV_factor/2d0*Bz_in_Tesla
!                 HmnR(n, m, ir)            = HmnR(n, m, ir)+ Zeeman_energy_in_eV_factor/2d0*Bz_in_Tesla
!                 HmnR(n      , m+nwann_nsoc, ir)=                 Zeeman_energy_in_eV_factor/2d0*Bx_in_Tesla &
!                                                           - zi*Zeeman_energy_in_eV_factor/2d0*By_in_Tesla 
!                 HmnR(n+nwann_nsoc, m      , ir)=                 Zeeman_energy_in_eV_factor/2d0*Bx_in_Tesla &
!                                                           + zi*Zeeman_energy_in_eV_factor/2d0*By_in_Tesla  
!              else
!                 HmnR(n+nwann_nsoc, m+nwann_nsoc, ir)= HmnR(n, m, ir)
!              endif ! R=0
!           enddo ! m
!        enddo ! n
!     enddo ! ir
!  else
!     do ir=1, nrpts
!        if (irvec(1, ir)/=0.or.irvec(2, ir)/=0.or.irvec(3, ir)/=0) cycle
!        do n=1, nwann_nsoc
!           do m=1, nwann_nsoc
!              if (n==m) then
!                 HmnR(n+nwann_nsoc, m+nwann_nsoc, ir)= HmnR(n+nwann_nsoc, &
!                 m+nwann_nsoc, ir)-Zeeman_energy_in_eV_factor/2d0*Bz_in_Tesla
!                 HmnR(n      , m+nwann_nsoc, ir)= HmnR(n, m+nwann_nsoc, ir)+ &
!                   Zeeman_energy_in_eV_factor/2d0*Bx_in_Tesla &
!                                                           -zi*Zeeman_energy_in_eV_factor/2d0*By_in_Tesla  
!                 HmnR(n+nwann_nsoc, m      , ir)= HmnR(n+nwann_nsoc, m, ir)+  &
!                  Zeeman_energy_in_eV_factor/2d0*Bx_in_Tesla &
!                                                           +zi*Zeeman_energy_in_eV_factor/2d0*By_in_Tesla  
!                 HmnR(n, m, ir)= HmnR(n, m, ir)+Zeeman_energy_in_eV_factor/2d0*Bz_in_Tesla
!              endif ! R=0
!           enddo ! m
!        enddo ! n
!     enddo ! ir
!  endif ! SOC


!  return

!end subroutine add_zeeman_normal_hr
