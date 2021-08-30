

subroutine readNormalHmnR()
   !>> Read in the tight-binding model from wannier90_hr.dat
   !>  The format is defined by the wannier90 software
   ! Constructed by quansheng wu 4/2/2010
   !
   ! Yifei Guan added the sparse hr file parsing June/2018
   ! License: GPL V3

   use para
   !> in: N of wann
   !> out : nth atom

   implicit none

   character*4 :: c_temp

   integer :: i, j, ir, ia, io 
   integer :: i1, i2, i3, i4, i5
   integer :: n, m, ir0
   integer :: add_electric_field
   integer :: nwann, nwann_nsoc

   real(dp) :: static_potential
   real(dp) :: tot, rh, ih
   real(dp) :: pos(Origin_cell%Num_atoms)


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
      read(12, *)Nrpts

      !> The degeneracy of each R point, this is caused by the Fourier transform in the Wigner-Seitz cell
      read(12, *)(ndegen(i), i=1, Nrpts)
      ir=0
      do ir=1,Nrpts
         do n=1,nwann
            do m=1,nwann
               read(12,*,end=1001)i1, i2, i3, i4, i5, rh, ih
               irvec(1,ir)=i1
               irvec(2,ir)=i2
               irvec(3,ir)=i3
               HmnR(i4,i5,ir)=dcmplx(rh,ih)
            end do
         enddo
      enddo

      !> extract the fermi energy
      do iR=1,Nrpts
         if (Irvec(1,iR).eq.0.and.Irvec(2,iR).eq.0.and.Irvec(3,iR).eq.0)then
            do i=1, Num_wann
               HmnR(i,i,iR)=HmnR(i,i,iR)-E_fermi
            enddo
         endif
      enddo
      !> WannierTools codes use Hatree atomic units
      if (index(Particle,'phonon')/=0) then
         HmnR= HmnR*eV2Hartree*eV2Hartree ! from eV to Hartree
      else
         HmnR= HmnR*eV2Hartree ! from eV to Hartree
      endif

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

      if (cpuid==0) then
         open(unit=105, file='wannier90_hr.dat')
         write(105, *)'hr file transformed from HWR'
         write(105, *)Num_wann
         write(105, *)nrpts
         write(105, '(15I5)')(ndegen(i), i=1, nrpts)
         do ir=1, nrpts
            do i=1, Num_wann
               do j=1, Num_wann
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

   !> change from "up dn up dn" to "up up dn dn"
   if (index( Package, 'QE')/=0.or.index( Package, 'quantumespresso')/=0 &
         .or.index( Package, 'quantum-espresso')/=0.or.index( Package, 'pwscf')/=0) then
      call reorder_wannierbasis

      if (cpuid==0.and.export_newhr) then
         !> write to new_hr.dat
         outfileindex= outfileindex+ 1
         open(unit=outfileindex, file='wannier90_hr_standard.dat')
         write(outfileindex, '(a,1X,a,1X,a,1X, a, a)')  &
            'HmnR transformed from QE ', time_now, date_now, 'UTC', zone_now
         write(outfileindex, *)Num_wann
         write(outfileindex, *)nrpts
         write(outfileindex, '(15I5)')ndegen
         do ir=1, nrpts
            do j=1, Num_wann
               do i=1, Num_wann
                  write( outfileindex, '(5I5, 2f16.6)') &
                     irvec(:, ir), i, j, HmnR(i, j, ir)/eV2Hartree
               end do
            end do
         end do
         close(outfileindex)
      endif
   endif


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
            write(stdout, '(1x, a, 3f16.6)')"Zeeman_energy_in_eV: ",  Zeeman_energy_in_eV/eV2Hartree
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
         !static_potential= pos(ia)*Origin_cell%cell_parameters(add_electric_field)*Electric_field_in_eVpA
         !if (Inner_symmetrical_Electric_Field) then
         !   static_potential= abs(pos(ia)*Origin_cell%cell_parameters(add_electric_field))*Electric_field_in_eVpA
         !endif
         static_potential= abs(pos(ia)*Origin_cell%cell_parameters(add_electric_field))&
            *Symmetrical_Electric_field_in_eVpA*eV2Hartree/Angstrom2atomic+ &
            (pos(ia)*Origin_cell%cell_parameters(add_electric_field))*&
            Electric_field_in_eVpA*eV2Hartree/Angstrom2atomic

         do i=1, Origin_cell%nprojs(ia)
            io=io+1
            HmnR(io, io, ir0)= HmnR(io, io, ir0)+ static_potential
            if (SOC>0) then
               HmnR(io+Num_wann/2, io+Num_wann/2, ir0)= HmnR(io+Num_wann/2, io+Num_wann/2, ir0)+ static_potential
            endif ! SOC
         enddo ! nproj
      enddo ! ia
   endif  ! add electric field or not

   call get_hmnr_cell(Cell_defined_by_surface)

   return
end subroutine readNormalHmnR


subroutine get_hmnr_cell(cell)
   !> Get new hmnr for a new cell with the same size as the previous one
   use para
   implicit none

   type(cell_type) :: cell

   !type(dense_tb_hr) :: cell_hr

   integer :: ir, i, j, iter
   real(dp) :: shift_vec_direct(3)

   !> for newcell
   real(dp) :: apos1d(3),apos2d(3)
   !>count newcell nrpts
   integer :: max_ir
   integer :: nir1,nir2,nir3, ir_cell
   integer :: nrpts_new, nrpts_max
   real(dp) :: new_ia, new_ib, new_ic, max_val

   !> all new irs
   integer, allocatable  :: rpts_array(:, :, :), rpts_map(:, :, :)
   integer, allocatable :: allirs(:, :, :, :)
   integer :: irn1(3),irn2(3)

   max_ir=8
   nrpts_max=(2*max_ir+1)**3
   allocate( rpts_array(-max_ir:max_ir,-max_ir:max_ir,-max_ir:max_ir))
   allocate( rpts_map(-max_ir:max_ir,-max_ir:max_ir,-max_ir:max_ir))

   allocate(allirs(nrpts, num_wann, num_wann, 3))

   call cart_direct_real(shift_to_topsurface_cart, shift_vec_direct, cell%lattice)

   call date_and_time(DATE=date_now,ZONE=zone_now, TIME=time_now)
   !> Get new Hrs
   rpts_array=0
   nrpts_new=0
   rpts_map= 0

   !> get number of R points for the new cell first 
   do ir=1, nrpts
      do i=1, Num_wann
         do j=1, Num_wann
            apos1d= Origin_cell%wannier_centers_direct(:, i)- shift_vec_direct
            apos2d= Origin_cell%wannier_centers_direct(:, j)- shift_vec_direct+irvec(:, ir)
            call latticetransform(apos1d(1),apos1d(2),apos1d(3),new_ia,new_ib,new_ic)
            irn1=floor([new_ia,new_ib,new_ic])

            call latticetransform(apos2d(1),apos2d(2),apos2d(3),new_ia,new_ib,new_ic)
            irn2=floor([new_ia,new_ib,new_ic])

            nir1=irn2(1)-irn1(1)
            nir2=irn2(2)-irn1(2)
            nir3=irn2(3)-irn1(3)
            if (abs(nir1)>max_ir .or. abs(nir2)>max_ir .or. abs(nir3)>max_ir) cycle
            rpts_array(nir1,nir2,nir3)=1
         enddo
      enddo
   enddo

   !> find all irvec
   nrpts_new= sum(rpts_array)
   allocate(irvec_newcell(3, Nrpts_new))
   iter= 0
   do nir3=-max_ir,max_ir
      do nir2=-max_ir,max_ir
         do nir1=-max_ir,max_ir
            if (rpts_array(nir1, nir2, nir3)==1) then
               iter=iter+1
               irvec_newcell(:, iter)=[nir1, nir2, nir3]
               rpts_map(nir1, nir2, nir3)=iter
            endif
         enddo
      enddo
   enddo

   !>> hamiltonian for the new cell
   allocate(HmnR_newcell(Num_wann, Num_wann, Nrpts_new))
   allocate(ndegen_newcell(Nrpts_new))
   HmnR_newcell= 0d0
   ndegen_newcell= 1

   !> get number of R points for the new cell first 
   do ir=1, nrpts
      do j=1, Num_wann
         do i=1, Num_wann
           !ia1=Origin_cell%spinorbital_to_atom_index(i)
           !ia2=Origin_cell%spinorbital_to_atom_index(j)
           !apos1d= Origin_cell%Atom_position_direct(:, ia1)- shift_vec_direct
           !apos2d= Origin_cell%Atom_position_direct(:, ia2)- shift_vec_direct+irvec(:, ir)
            apos1d= Origin_cell%wannier_centers_direct(:, i)- shift_vec_direct
            apos2d= Origin_cell%wannier_centers_direct(:, j)- shift_vec_direct+irvec(:, ir)
            call latticetransform(apos1d(1),apos1d(2),apos1d(3),new_ia,new_ib,new_ic)
            irn1=floor([new_ia,new_ib,new_ic])

            call latticetransform(apos2d(1),apos2d(2),apos2d(3),new_ia,new_ib,new_ic)
            irn2=floor([new_ia,new_ib,new_ic])

            nir1=irn2(1)-irn1(1)
            nir2=irn2(2)-irn1(2)
            nir3=irn2(3)-irn1(3)
            if (abs(nir1)>max_ir .or. abs(nir2)>max_ir .or. abs(nir3)>max_ir) cycle
            ir_cell= rpts_map(nir1, nir2, nir3)
            HmnR_newcell(i,j,ir_cell)=HmnR(i,j,ir)/ndegen(ir)
         enddo
      enddo
   enddo

   !> do cut-off according to the hopping value
   iter= 0 
   do ir=1, nrpts_new
      max_val=maxval(abs(HmnR_newcell(:, :, ir)))/eV2Hartree
      if (max_val<eps4) then
         iter= iter+ 1
      endif
   enddo

   if (cpuid==0.and.export_newhr) then
      !> write to new_hr.dat
      nrpts_new=sum(rpts_array)-iter
      outfileindex= outfileindex+ 1
      open(unit=outfileindex, file='wannier90_hr_newcell.dat')
      write(outfileindex, '(a,1X,a,1X,a,1X, a, a)')  &
         'HmnR for new cell generated at ', time_now, date_now, 'UTC', zone_now
      write(outfileindex, *)Num_wann
      write(outfileindex, *)nrpts_new
      write(outfileindex, '(15I5)')(1 , i=1, nrpts_new)
      do ir=1, nrpts_new
         max_val=maxval(abs(HmnR_newcell(:, :, ir)))/eV2Hartree
        !if (max_val<eps4) cycle
         do j=1, Num_wann
            do i=1, Num_wann
               write( outfileindex, '(5I5, 2f16.6)') &
                  irvec_newcell(:, ir), i, j, HmnR_newcell(i, j, ir)/eV2Hartree
            end do
         end do
      end do
      close(outfileindex)
   endif

   return
end subroutine get_hmnr_cell

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
   integer :: q,n,m, ir0 
   real(8) :: r1,r2, pos(Origin_cell%Num_atoms), static_potential

   !> the direction which adding electric field which is also the stacking direction
   integer :: add_electric_field
   real(dp) :: Bx_in_au, By_in_au, Bz_in_au
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

   !> in order to include the Fermi level
   splen=splen+nwann

   allocate(hacoo(splen),hicoo(splen),hjcoo(splen),hirv(splen))
   hacoo=(0d0, 0d0)
   hicoo=0
   hjcoo=0
   hirv=0

   j=0
   ir0=0
   ir=1
   irvec=-9999
   read(12,*,end=1001)i1, i2, i3, i4, i5, r1, r2
   irvec(1,ir)=i1
   irvec(2,ir)=i2
   irvec(3,ir)=i3
   !> will reread the above line
   backspace(12)
   do q=1,nrpts
      do n=1,nwann
         do m=1,nwann
            read(12,*,end=1001)i1, i2, i3, i4, i5, r1, r2
            if (sum(abs(irvec(:,ir)-[i1,i2,i3]))/=0) ir=ir+1
            j=j+1
            hicoo(j)=i4
            hjcoo(j)=i5
            hirv (j)=ir
            hacoo(j)=dcmplx(r1,r2)*eV2Hartree
            if (i1==0.and.i2==0.and.i3==0.and.i4==i5) then
               ir0= ir
            endif
            irvec(1,ir)=i1
            irvec(2,ir)=i2
            irvec(3,ir)=i3
         end do
      enddo
   enddo
1001 continue
   !> correct nrpts
   Nrpts=ir

   if (cpuid.eq.0) write(stdout, '(a, i6)')' >> NRPTS is ', Nrpts

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
            write(stdout, '(1x, a, 3f16.6)')"Zeeman_energy_in_eV: ",  Zeeman_energy_in_eV/eV2Hartree
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

      Bx_in_au= Bx
      By_in_au= By
      Bz_in_au= Bz
      call add_zeeman_sparse_hr(Bx_in_au, By_in_au, Bz_in_au)
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
            !static_potential= pos(ia)*Origin_cell%cell_parameters(add_electric_field)*Electric_field_in_eVpA
            !if (Inner_symmetrical_Electric_Field) then
            !   static_potential= abs(pos(ia)*Origin_cell%cell_parameters(add_electric_field))*Electric_field_in_eVpA
            !endif
            static_potential= abs(pos(ia)*Origin_cell%cell_parameters(add_electric_field))*Symmetrical_Electric_field_in_eVpA+ &
               (pos(ia)*Origin_cell%cell_parameters(add_electric_field))*Electric_field_in_eVpA
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
            !static_potential= pos(ia)*Origin_cell%cell_parameters(add_electric_field)*Electric_field_in_eVpA
            !if (Inner_symmetrical_Electric_Field) then
            !   static_potential= abs(pos(ia)*Origin_cell%cell_parameters(add_electric_field))*Electric_field_in_eVpA
            !endif
            static_potential= abs(pos(ia)*Origin_cell%cell_parameters(add_electric_field))*Symmetrical_Electric_field_in_eVpA+ &
               (pos(ia)*Origin_cell%cell_parameters(add_electric_field))*Electric_field_in_eVpA
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

   !> add Fermi level
   do n=1,num_wann
      j=j+1
      hicoo(j)=n
      hjcoo(j)=n
      hirv (j)=ir0
      hacoo(j)= - E_FERMI*eV2Hartree
   enddo

   !> transform it into dense format only when number of wannier orbitals is less than 100

   outfileindex= outfileindex+ 1
   if (cpuid.eq.0.and. nwann<100) then
      open (unit=outfileindex, file='hr.dat-dense')
      write(outfileindex, *) ' ! HmnR file from sparse hr file'
      write(outfileindex, '(I10, a)') nwann, '  ! Num_wann: number of wannier orbitals'
      write(outfileindex, '(I10, a)') Nrpts, '  ! NRPTS: number of R points'
      write(outfileindex, '(10I5)') ndegen(:)
      do iR=1, NRPTS
         do n=1, nwann
            do m=1, nwann
               h_value=0d0
               do i=1, splen
                  if (hicoo(i)==n.and.hjcoo(i)==m.and.hirv(i)==iR)then
                     h_value= h_value+ hacoo(i)
                  endif
               enddo !i
               write(outfileindex, '(5I5, 2f12.6)')irvec(:, iR), n, m, h_value/eV2Hartree
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
   if (abs(Symmetrical_Electric_field_in_eVpA)>eps6.or.abs(Electric_field_in_eVpA)>eps6) then
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
      write(stdout, '(a)') "  Note: This is not a 2D system or E=0, so we can't add electric field."
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


subroutine add_zeeman_sparse_hr(Bx_in_au, By_in_au, Bz_in_au)
   !> add Zeeman energy on the sparse hmnr based on the magnetic field strength
   !> magnetic_field_in_au is in au
   !> Here the coordinate system is the same as user-specified in the LATTICE card.
   use para
   implicit none

   real(dp), intent(in) :: Bx_in_au, By_in_au, Bz_in_au

   integer :: nwann_nsoc, i, j, ir0, ir
   real(dp) :: Zeeman_energy_in_hartree_factor

   nwann_nsoc= Num_wann/2

   if (.not.Add_Zeeman_Field) return
   if (cpuid==0) write(stdout, *)'>> We are adding Zeeman term on sparse Hamiltonian'

   ir0=0
   do ir=1, Nrpts
      if (irvec(1,ir)==0.and.irvec(2,ir)==0.and.irvec(3,ir)==0) then
         ir0= ir
         exit
      endif
   enddo
   if (ir0==0) stop 'something wrong with irvec in subroutine add_zeeman_sparse_hr'

   Zeeman_energy_in_hartree_factor= Effective_gfactor*Bohr_magneton
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
      hacoo(j)= Zeeman_energy_in_hartree_factor/2d0*Bz_in_au
   enddo
   do i=1, nwann_nsoc
      j= j+ 1
      hicoo(j)= i+ nwann_nsoc
      hjcoo(j)= i+ nwann_nsoc
      hirv (j)= ir0
      hacoo(j)= -Zeeman_energy_in_hartree_factor/2d0*Bz_in_au
   enddo
   do i=1, nwann_nsoc
      j= j+ 1
      hicoo(j)= i
      hjcoo(j)= i+ nwann_nsoc
      hirv (j)= ir0
      hacoo(j)=  Zeeman_energy_in_hartree_factor/2d0*(By_in_au*zi+ Bx_in_au)
   enddo
   do i=1, nwann_nsoc
      j= j+ 1
      hicoo(j)= i+ nwann_nsoc
      hjcoo(j)= i
      hirv (j)= ir0
      hacoo(j)=  Zeeman_energy_in_hartree_factor/2d0*(-By_in_au*zi+ Bx_in_au)
   enddo
   if (j>splen) stop "ERROR happend in setting zeeman energy in readhmnr.f90"

   return

end subroutine add_zeeman_sparse_hr

subroutine reorder_wannierbasis
   !> change from "up dn up dn" to "up up dn dn"
   use para
   implicit none

   integer :: i, j, ir, nwann
   integer, allocatable :: orderi(:)
   complex(dp), allocatable :: hmn_temp(:, :)

   if (SOC==0) return

   allocate(hmn_temp(num_wann,num_wann))

   allocate(orderi(num_wann))
   hmn_temp=0d0
   orderi=1

   nwann=num_wann/2
   do i=1, nwann
      orderi(i)= 2*i-1
      orderi(i+nwann)= 2*i  
   enddo

   do ir=1, Nrpts
      do i=1, num_wann
         do j=1, num_wann
            hmn_temp(i, j)= HmnR(orderi(i), orderi(j), ir)
         enddo
      enddo
      HmnR(:, :, ir)= hmn_temp
   enddo

   return
end subroutine reorder_wannierbasis





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
!  do ir=1, Nrpts
!     if (irvec(1,ir)==0.and.irvec(2,ir)==0.and.irvec(3,ir)==0) then
!        ir0= ir
!     endif
!  enddo
!  if (ir0==0) stop 'something wrong with irvec in subroutine add_zeeman_sparse_hr'

!  Zeeman_energy_in_eV_factor= Effective_gfactor*Bohr_magneton

!  if (SOC==0) then
!     do ir=1, Nrpts
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
!     do ir=1, Nrpts
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
