

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
   integer :: stat, idx, nRused
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
      if (.not. Is_Sparse_Hr) then

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

         if ((soc==0 .and. sum(Origin_cell%nprojs)/=nwann .and. .not.Add_Zeeman_Field) .or. &
         (soc>0 .and. sum(Origin_cell%nprojs)/=nwann/2))then
         print *, 'sum(Origin_cell%nprojs), num_wann, num_wann/2'
         print *, sum(Origin_cell%nprojs), nwann, nwann/2
         print *, "ERROR: Maybe the SOC tags in the SYSTEM is wrongly set"
         stop "ERROR: the summation of all projectors times spin degeneracy is not equal to num_wann"
         endif
   
         !> number of lattice vectors taken into account
         read(12, *)Nrpts
         if (.not. allocated(HmnR)) allocate(HmnR(Num_wann,Num_wann,nrpts))
         allocate(irvec(3,nrpts))
         allocate(ndegen(nrpts))
         irvec= 0
         ndegen=1

         !> The degeneracy of each R point, this is caused by the Fourier transform in the Wigner-Seitz cell
         read(12, *)(ndegen(i), i=1, Nrpts)
      else
         !> for Sparse HmnR
         !> skip a comment line
         read(12, *)
         read(12,*) splen_input   !> number of non-zeros lines
         if (splen_input<=0) then
            stop "ERROR : splen_input is zero in hr file"
         endif

         !> 读取 Wannier 轨道数
         read(12, *) nwann
         if (nwann /= Num_wann) then
            write(*, *) '>>> ERROR: Num_wann mismatch in Hr file: ', nwann, ' vs ', Num_wann
            stop
         endif
         nwann_nsoc = nwann
         if (SOC > 0) nwann_nsoc = nwann/2

         !> 读取 R 点个数
         read(12, *) Nrpts
         if (.not. allocated(HmnR)) allocate(HmnR(Num_wann,Num_wann,nrpts))
         allocate(irvec(3, Nrpts))
         allocate(ndegen(Nrpts))
         ndegen = 1
      endif

      nRused = 0
      do 
         read(12,*,iostat=stat) i1, i2, i3, i4, i5, rh, ih
         if (stat<0) exit   ! EOF
         if (stat>0) stop "格式错误"
         !—— 1) 在已有 R 点中查找 ——
         idx = 0
         do j = 1, nRused
            if (irvec(1,j)==i1 .and. irvec(2,j)==i2 .and. irvec(3,j)==i3) then
               idx = j
               exit
            end if
         end do
      
         !—— 2) 如果没找到，就新增一个 R 点 ——
         if (idx == 0) then
            nRused = nRused + 1
            idx = nRused
            irvec(:,idx) = (/ i1, i2, i3 /)
            ! write(*, '(a, i5, a, i5, a, i5)') 'New R point found: i1 =', i1, ', i2=', i2, ', i3=', i3
         end if
      
         !—— 3) 填写 overlap 矩阵 ——
         HmnR(i4, i5, idx) = dcmplx(rh, ih)

      end do

      ! ir=0
      ! do ir=1,Nrpts
      !    do n=1,nwann
      !       do m=1,nwann
      !          read(12,*,end=1001)i1, i2, i3, i4, i5, rh, ih
      !          irvec(1,ir)=i1
      !          irvec(2,ir)=i2
      !          irvec(3,ir)=i3
      !          HmnR(i4,i5,ir)=dcmplx(rh,ih)
      !       end do
      !    enddo
      ! enddo

      !> extract the fermi energy
      if (Orthogonal_Basis) then
         do iR=1,Nrpts
            if (Irvec(1,iR).eq.0.and.Irvec(2,iR).eq.0.and.Irvec(3,iR).eq.0)then
               do i=1, Num_wann
                  HmnR(i,i,iR)=HmnR(i,i,iR)-E_fermi
                  ! write(*, '(a, 5i5, a, f18.7)') 'HmnR at R=0', iR, i, irvec(:,iR), ' is ', HmnR(i,i,iR)
               enddo
            endif
         enddo
      endif

      !直接结束所有程序
      ! stop "ERROR: Orthogonal_Basis is not supported in this version of Wannier90"
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
         .or.index( Package, 'quantum-espresso')/=0.or.index(Package, 'VASP6')/=0.or.index( Package, 'pwscf')/=0) then
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

   !> get the ir index that the R(ir)=(0, 0, 0)
   ir0=0
   do ir=1, nrpts
      if (irvec(1, ir)==0.and.irvec(2, ir)==0.and.irvec(3, ir)==0) ir0=ir
   enddo

   !> Adding zeeman field
   !> Bx=Bdirection(1)
   !> By=Bdirection(2)
   !> Bz=Bdirection(3)
   !> Hz= Zeeman_energy_in_eV*(Bx*sx+By*sy+Bz*sz)/2d0
   !> sx, sy, sz are Pauli matrices.
   if (Add_Zeeman_Field) then
      call add_zeeman_normal_hr()
      !> After considering the Zeeman field, we already extended the spin space to spin-full.
      SOC = 1
   endif ! Add_Zeeman_Field

   call get_stacking_direction_and_pos(add_electric_field, pos)
 
   if (add_electric_field>0) then
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

   !> write out Hmn(R=0)
   if (cpuid.eq.0 .and. Num_wann< 300)then
      write(stdout, '(a)')" "
      write(stdout, '(a)')" >> Hopping parameters in the home unit cell"
      write(stdout, '(a)')" >> H00= Hmn(R=0) real part"
      do i=1, Num_wann
         write(stdout, '(50000f7.3)') real(HmnR(i, :, ir0))/eV2Hartree
      enddo
      write(stdout, '(a)')" "
      write(stdout, '(a)')" >> H00= Hmn(R=0) imagary part"
      do i=1, Num_wann
         write(stdout, '(50000f7.3)') aimag(HmnR(i, :, ir0))/eV2Hartree
      enddo
      write(stdout, '(a)')" "
   endif

   ! call get_hmnr_cell(Cell_defined_by_surface)

   return
end subroutine readNormalHmnR


subroutine read_OAM_operator(Lx_wann, Ly_wann, Lz_wann)
   !>> Read in the valley operator from valley_operator.dat
   ! Constructed by quansheng wu 04 Nov. 2023
   ! License: GPL V3

   use para

   implicit none


   integer :: n, m, ir0
   integer :: add_electric_field
   integer :: nwann, nwann_nsoc
   logical :: exists



   integer :: i, j, ios
   real :: re, im
   character(len=2) :: current_matrix
   character(len=256) :: line

   complex(dp), intent(out) :: Lx_wann(num_wann, num_wann), Ly_wann(num_wann, num_wann), Lz_wann(num_wann, num_wann)

   !allocate(Lx_wann(num_wann, num_wann))
   !allocate(Ly_wann(num_wann, num_wann))
   !allocate(Lz_wann(num_wann, num_wann))

   if(cpuid.eq.0)write(stdout,*)'****************OAM******************** '
   inquire (file ="oam.dat", EXIST = exists)
   if (exists)then
      write(stdout,*)'****************OAM******************** '
      open(12, file="oam.dat", status='OLD', iostat = ios)
      !if (ios /= 0) then
      !  write(*,*) "Error opening file!"
      !  stop
      !end if
   else
      STOP ">> for OAM projection , you have to prepare a file oam.dat"
   endif


   
   !Lx_wann = 0
   !Ly_wann = 0
   !Lz_wann = 0
   read(12,*) nwann
   write(stdout,*) 'num_wann_OAM', num_wann
   do
        read(12, '(A)', iostat=ios) line
        if (ios /= 0) exit

        ! 检查是否是矩阵标识行
        if (line(1:2) == 'Lx' .or. line(1:2) == 'Ly' .or. line(1:2) == 'Lz') then
            current_matrix = line(1:2)
            cycle
        end if

        ! 跳过空行
        if (len_trim(line) == 0) cycle

        ! 读取数据行
        read(line, *, iostat=ios) i, j, re, im
        write(stdout,*) i,j,re,im
        if (ios /= 0) then
            write(*,*) "Error reading line:", trim(line)
            cycle
        end if


         if (current_matrix == 'Lx') then
                Lx_wann(i,j) = cmplx(re, im)
         else if (current_matrix == 'Ly') then
                Ly_wann(i,j) = cmplx(re, im)
         else if (current_matrix == 'Lz') then
                Lz_wann(i,j) = cmplx(re, im)
         else
            write(*,*) "Unknown matrix identifier:", current_matrix
        end if
    end do




   !1001 continue
   close(12)

   do i=1,num_wann
      do j=1,num_wann
      if (abs(Lx_wann(i,j))>0.0001) write(stdout,*) i,j, Lx_wann(i,j)
      enddo
   enddo

   if (cpuid.eq.0) write(stdout, *) ">>> finished reading of OAM operator"

end subroutine read_OAM_operator

subroutine readNormalSmnR()
   !>> Read in the overlap matrix from wannier90_sr.dat
   !>  The format is defined by the Wannier90 software
   !  Constructed by quansheng wu 4/2/2010
   !  Modified for overlap (SmnR)

   use para
   implicit none

   integer :: ir, n, m, i, j, k, idx, nRused
   integer :: i1, i2, i3, i4, i5
   integer :: nwann, nwann_nsoc
   real(dp) :: rh, ih
   ! complex(dp), allocatable :: SmnR(:,:,:)  ! overlap matrix
   logical :: exists

   if(cpuid.eq.0)write(stdout,*)' '
   if(cpuid.eq.0)write(stdout,'(2a)')' >> Now we are reading the overlap files: ', trim(adjustl(Overlapfile))
   inquire(file= Overlapfile, EXIST = exists)
   if (exists)then
      open(13, file= Overlapfile, status='OLD')
   else
      STOP ">> for non-orthogonal basis , you have to prepare a file like overlap.dat"
   endif

   if(.not. Is_Sparse_Hr) then
      write(*, *) '>>> ERROR: The overlap file is not sparse. Please use the sparse format.'
      stop
   else
      !> 跳过注释行
      read(13, *) 
      
      read(13,*) splen_overlap_input   !> number of non-zeros lines

      !> 读取 Wannier 轨道数
      read(13, *) nwann
      if (nwann /= Num_wann) then
         write(*, *) '>>> ERROR: Num_wann mismatch in Sr file: ', nwann, ' vs ', Num_wann
         stop
      endif
      nwann_nsoc = nwann
      if (SOC > 0) nwann_nsoc = nwann/2

      !> 读取 R 点个数
      read(13, *) Nrpts

      !> 分配数组
      allocate(SmnR(Num_wann, Num_wann, Nrpts))
      ! allocate(irvec(3,    Nrpts))
      ! allocate(ndegen(    Nrpts))

      nRused = Nrpts
      do k = 1, splen_overlap_input
         read(13,*) i1, i2, i3, i4, i5, rh, ih
      
         !—— 1) 在已有 R 点中查找 ——
         idx = 0
         do j = 1, nRused
            if (irvec(1,j)==i1 .and. irvec(2,j)==i2 .and. irvec(3,j)==i3) then
               idx = j
               exit
            end if
         end do

         !—— 2) 如果没找到，就新增一个 R 点 ——
         if (idx == 0) then
            ! nRused = nRused + 1
            ! idx = nRused
            ! irvec(:,idx) = (/ i1, i2, i3 /)
            ! write(*, '(a, i5, a, i5, a, i5)') 'New R point found: ', i1, i2, i3
            stop "ERROR: R point not found in SmnR reading"
         end if

         !—— 3) 填写 overlap 矩阵 ——

         SmnR(i4, i5, idx) = dcmplx(rh, ih)
         if (i4 == i5) then
            !打印i1,i2,i3,i4,i5以及值
            ! write(*, '(a, 5i5, a, f18.7)') 'SmnR at R=', i1, i2, i3, i4, i5, ' is ', SmnR(i4, i5, idx)
         endif
      end do

      do iR=1,Nrpts
         if (Irvec(1,iR).eq.0.and.Irvec(2,iR).eq.0.and.Irvec(3,iR).eq.0)then
            do i=1, Num_wann
               !打印ir, i, SmnR(i, i, iR)
               ! write(*, '(a, i5, a, i5, a, f18.7)') 'SmnR at R=0', iR, ' for i=', i, ' is ', SmnR(i, i, iR)
            enddo
         endif
      enddo

   endif



   close(13)
end subroutine readNormalSmnR




subroutine read_valley_operator()
   !>> Read in the valley operator from valley_operator.dat
   ! Constructed by quansheng wu 04 Nov. 2023
   ! License: GPL V3

   use para

   implicit none

   character*4 :: c_temp

   integer :: i, j, ir, ia, io 
   integer :: i1, i2, i3, i4, i5
   integer :: n, m, ir0
   integer :: add_electric_field
   integer :: nwann, nwann_nsoc
   logical :: exists

   real(dp) :: static_potential
   real(dp) :: tot, rh, ih
   real(dp) :: pos(Origin_cell%Num_atoms)


   if(cpuid.eq.0)write(stdout,*)' '
   inquire (file ="valley_operator.dat", EXIST = exists)
   if (exists)then
      open(12, file="valley_operator.dat", status='OLD')
   else
      STOP ">> for valley projection , you have to prepare a file valley_operator.dat"
   endif

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
   read(12, *)Nrpts_valley

   !> The degeneracy of each R point, this is caused by the Fourier transform in the Wigner-Seitz cell
   read(12, *)(j, i=1, Nrpts_valley)
   allocate(irvec_valley(3, Nrpts_valley))
   allocate(valley_operator_R(num_wann, num_wann, Nrpts_valley))
   ir=0
   do ir=1,Nrpts_valley
      do n=1,nwann
         do m=1,nwann
            read(12,*,end=1001)i1, i2, i3, i4, i5, rh, ih
            irvec_valley(1,ir)=i1
            irvec_valley(2,ir)=i2
            irvec_valley(3,ir)=i3
            valley_operator_R(i4,i5,ir)=dcmplx(rh,ih)
         end do
      enddo
   enddo

   1001 continue
   close(12)

   if (cpuid.eq.0) write(stdout, *) ">>> finished reading of valley operator"

end subroutine read_valley_operator



subroutine get_hmnr_cell(cell)
   !> Get new hmnr for a new cell with the same size as the previous one
   use para
   implicit none

   type(cell_type) :: cell

   !type(dense_tb_hr) :: cell_hr

   integer :: ir, i, j, iter
   real(dp) :: shift_vec_direct(3)

   !> for newcell
   real(dp) :: apos1d(3),apos2d(3),apos1dprime(3),apos2dprime(3)
   !>count newcell nrpts
   integer :: max_ir
   integer :: nir1,nir2,nir3, ir_cell
   integer :: nrpts_new, nrpts_max
   real(dp) :: new_ia, new_ib, new_ic, max_val

   !> all new irs
   integer, allocatable  :: rpts_array(:, :, :), rpts_map(:, :, :)
   integer, allocatable :: allirs(:, :, :, :)
   integer :: irn1(3),irn2(3)

   !> The Allocate Process
   integer :: ia1,ia2,ia1prime,ia2prime

   max_ir=8
   nrpts_max=(2*max_ir+1)**3
   allocate( rpts_array(-max_ir:max_ir,-max_ir:max_ir,-max_ir:max_ir))
   allocate( rpts_map(-max_ir:max_ir,-max_ir:max_ir,-max_ir:max_ir))

   allocate(allirs(nrpts, num_wann, num_wann, 3))

   allocate(Rmn_old(3))
   allocate(Rmn_new(3))
   allocate(irvec_new(3))
   allocate(irvec_new_int(3))

   ! call cart_direct_real(shift_to_topsurface_cart, shift_vec_direct, cell%lattice)

   call date_and_time(DATE=date_now,ZONE=zone_now, TIME=time_now)
   !> Get new Hrs
   rpts_array=0
   nrpts_surfacecell=0
   rpts_map= 0

   !> 1. Get number of R points for the new cell
   do ir=1, nrpts
      do j=1, Num_wann
         do i=1, Num_wann
            ia1 = Origin_cell%spinorbital_to_atom_index(i)
            ia2 = Origin_cell%spinorbital_to_atom_index(j)
            apos1d = Origin_cell%Atom_position_direct(:, ia1)
            apos2d = Origin_cell%Atom_position_direct(:, ia2)

            ia1prime = Cell_defined_by_surface%spinorbital_to_atom_index(i)
            ia2prime = Cell_defined_by_surface%spinorbital_to_atom_index(j)
            apos1dprime = Cell_defined_by_surface%Atom_position_direct(:,ia1prime)
            apos2dprime = Cell_defined_by_surface%Atom_position_direct(:,ia2prime)

            Rmn_old = irvec(:, ir) + apos2d - apos1d
            call latticetransform(Rmn_old(1),Rmn_old(2),Rmn_old(3),Rmn_new(1),Rmn_new(2),Rmn_new(3))
            irvec_new = Rmn_new - (apos2dprime - apos1dprime)

            !> Due to the accuracy of computing, need to rounding  (but always tiny)
            irvec_new_int(1) = ANINT(irvec_new(1))
            irvec_new_int(2) = ANINT(irvec_new(2))
            irvec_new_int(3) = ANINT(irvec_new(3))

            if (abs(irvec_new_int(1))>max_ir .or. abs(irvec_new_int(2))>max_ir .or. abs(irvec_new_int(3))>max_ir) cycle
            rpts_array(irvec_new_int(1),irvec_new_int(2),irvec_new_int(3))=1

         enddo
      enddo
   enddo

   !> The total number of lattice points searched above
   nrpts_surfacecell= sum(rpts_array)

   allocate(irvec_surfacecell(3, nrpts_surfacecell))
   irvec_surfacecell = 0

   !> 2. Create an order map to sign the R points generated in step1

   iter = 0
   do nir3=-max_ir,max_ir
      do nir2=-max_ir,max_ir
         do nir1=-max_ir,max_ir
            if (rpts_array(nir1, nir2, nir3)==1) then
               iter=iter+1
               irvec_surfacecell(:, iter)=[nir1, nir2, nir3]
               rpts_map(nir1, nir2, nir3)=iter
            endif
         enddo
      enddo
   enddo


   allocate(HmnR_surfacecell(Num_wann, Num_wann, nrpts_surfacecell))
   allocate(ndegen_surfacecell(nrpts_surfacecell))
   HmnR_surfacecell= 0d0
   ndegen_surfacecell= 1

   if (.not. Orthogonal_Basis) then
      allocate(SmnR_surfacecell(Num_wann, Num_wann, nrpts_surfacecell))
      SmnR_surfacecell = 0d0
   endif

   !> 3. Allocate the old HmnR to the new HmnR.
   !>  Note: The cell can't be treated as a mass point. We need to consider the relative coordinates of atoms.
   do ir=1, nrpts
      do j=1, Num_wann
         do i=1, Num_wann

            ia1 = Origin_cell%spinorbital_to_atom_index(i)
            ia2 = Origin_cell%spinorbital_to_atom_index(j)
            apos1d = Origin_cell%Atom_position_direct(:, ia1)
            apos2d = Origin_cell%Atom_position_direct(:, ia2) 

            ia1prime = Cell_defined_by_surface%spinorbital_to_atom_index(i)
            ia2prime = Cell_defined_by_surface%spinorbital_to_atom_index(j)
            apos1dprime = Cell_defined_by_surface%Atom_position_direct(:,ia1prime)
            apos2dprime = Cell_defined_by_surface%Atom_position_direct(:,ia2prime)

            !> R'_mn = R' + tau'_2 - tau'_1
            !> R_mn = R + tau_2 - tau_1
            !> R'_mn = Pinv * R_mn
            !> R' = (Pinv * R_mn) - (tau'_2 - tau'_1)

            Rmn_old = irvec(:, ir) + apos2d - apos1d
            call latticetransform(Rmn_old(1),Rmn_old(2),Rmn_old(3),Rmn_new(1),Rmn_new(2),Rmn_new(3))
            irvec_new = Rmn_new - (apos2dprime - apos1dprime)

            !> For safety, we perform a rounding operation due to the accuracy of computing
            irvec_new_int(1) = ANINT(irvec_new(1))
            irvec_new_int(2) = ANINT(irvec_new(2))
            irvec_new_int(3) = ANINT(irvec_new(3))

            if (abs(irvec_new_int(1))>max_ir .or. abs(irvec_new_int(2))>max_ir .or. abs(irvec_new_int(3))>max_ir) cycle
            ir_cell = rpts_map(irvec_new_int(1), irvec_new_int(2), irvec_new_int(3))
            HmnR_surfacecell(i,j,ir_cell) = HmnR(i,j,ir)/ndegen(ir)
            
            if (.not. Orthogonal_Basis) then
               SmnR_surfacecell(i,j,ir_cell) = SmnR(i,j,ir)/ndegen(ir)
            endif 

         enddo
      enddo
   enddo

   ! !> do cut-off according to the hopping value
   ! iter= 0 
   ! do ir=1, nrpts_new
   !    max_val=maxval(abs(HmnR_newcell(:, :, ir)))/eV2Hartree
   !    if (max_val<eps4) then
   !       iter= iter+ 1
   !    endif
   ! enddo

   if (cpuid==0.and.export_newhr) then
      !> Problem: Just cut the nrpts_surfacecell, but the irvec is not well ordered. The order (iter) can be anywhere.
      !> write to new_hr.dat
      ! nrpts_surfacecell=sum(rpts_array)-iter

      nrpts_surfacecell=sum(rpts_array)
      outfileindex= outfileindex+ 1
      open(unit=outfileindex, file='wannier90_hr_newcell.dat')
      write(outfileindex, '(a,1X,a,1X,a,1X, a, a)')  &
         'HmnR for new cell generated at ', time_now, date_now, 'UTC', zone_now
      write(outfileindex, *)Num_wann
      write(outfileindex, *)nrpts_surfacecell
      write(outfileindex, '(15I5)')(1 , i=1, nrpts_surfacecell)
      do ir=1, nrpts_surfacecell
         max_val=maxval(abs(HmnR_surfacecell(:, :, ir)))/eV2Hartree
        !if (max_val<eps4) cycle
         do j=1, Num_wann
            do i=1, Num_wann
               write( outfileindex, '(5I5, 2f16.6)') &
                  irvec_surfacecell(:, ir), i, j, HmnR_surfacecell(i, j, ir)/eV2Hartree
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
         !> no electrical field and Zeeman field
         splen= splen_input
      endif
   end if

   !> in order to include the Fermi level
   !> for non-orthogonal_Basis, we can't include the fermi level in to H'
   if (.not.Orthogonal_Basis) then
      splen=splen
   else
      splen=splen+nwann
   endif

   allocate(hacoo(splen),hicoo(splen),hjcoo(splen),hirv(3, splen))
   hacoo=(0d0, 0d0)
   hicoo=0
   hjcoo=0
   hirv=0

   !> will reread the above line
   do j=1, splen_input
      read(12,*,end=1001)i1, i2, i3, i4, i5, r1, r2
      hirv (1, j)=i1
      hirv (2, j)=i2
      hirv (3, j)=i3
      hicoo(j)=i4
      hjcoo(j)=i5

      if (trim(adjustl(Package))=="OPENMX") then
         hacoo(j)=dcmplx(r1,r2)
      else
         hacoo(j)=dcmplx(r1,r2)*eV2Hartree
      endif
   enddo
   j=j-1
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
            hirv (1, j)=hirv (1, i)
            hirv (2, j)=hirv (2, i)
            hirv (3, j)=hirv (3, i)
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
            static_potential= abs(pos(ia)*Origin_cell%cell_parameters(add_electric_field))&
            *Symmetrical_Electric_field_in_eVpA*eV2Hartree/Angstrom2atomic+ &
            (pos(ia)*Origin_cell%cell_parameters(add_electric_field))*&
            Electric_field_in_eVpA*eV2Hartree/Angstrom2atomic
            do i=1, Origin_cell%nprojs(ia)
               io=io+1
               !> spin up
               j=j+1
               hicoo(j)= io
               hjcoo(j)= io
               hirv (1, j)= 0
               hirv (2, j)= 0
               hirv (3, j)= 0
               hacoo(j)= static_potential
               !> spin down
               j=j+1
               hicoo(j)= io + nwann_nsoc
               hjcoo(j)= io + nwann_nsoc
               hirv (1, j)= 0
               hirv (2, j)= 0
               hirv (3, j)= 0
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
            static_potential= abs(pos(ia)*Origin_cell%cell_parameters(add_electric_field))&
            *Symmetrical_Electric_field_in_eVpA*eV2Hartree/Angstrom2atomic+ &
            (pos(ia)*Origin_cell%cell_parameters(add_electric_field))*&
            Electric_field_in_eVpA*eV2Hartree/Angstrom2atomic

            do i=1, Origin_cell%nprojs(ia)
               io=io+1
               j=j+1
               hicoo(j)= io
               hjcoo(j)= io
               hirv (1, j)= 0
               hirv (2, j)= 0
               hirv (3, j)= 0
               hacoo(j)= static_potential
               if (SOC>0) then  ! with SOC
                  j=j+1
                  hicoo(j)= io + nwann_nsoc
                  hjcoo(j)= io + nwann_nsoc
                  hirv (1, j)= 0
                  hirv (2, j)= 0
                  hirv (3, j)= 0
                  hacoo(j)= static_potential
               endif
            enddo ! nproj
         enddo ! ia
      endif ! Add_Zeeman_Field or not
   endif ! add_electric_field

   !> the Fermi level can be added to the Hamiltonian directly only with in the orthogonal basis
   if (Orthogonal_Basis) then
      !> add Fermi level
      do n=1,num_wann
         j=j+1
         hicoo(j)=n
         hjcoo(j)=n
         hirv (1, j)= 0
         hirv (2, j)= 0
         hirv (3, j)= 0
         hacoo(j)= - E_FERMI*eV2Hartree
      enddo
   endif

   !> transform it into dense format only when number of wannier orbitals is less than 100

   outfileindex= outfileindex+ 1
   if (cpuid.eq.0.and. nwann<500) then
!     open (unit=outfileindex, file='hr.dat-dense')
!     write(outfileindex, *) ' ! HmnR file from sparse hr file'
!     write(outfileindex, '(I10, a)') nwann, '  ! Num_wann: number of wannier orbitals'
!     write(outfileindex, '(I10, a)') Nrpts, '  ! NRPTS: number of R points'
!     write(outfileindex, '(10I5)') ndegen(:)
!     do iR=1, NRPTS
!        i1= irvec(1, iR)
!        i2= irvec(2, iR)
!        i3= irvec(3, iR)
!        do n=1, nwann
!           do m=1, nwann
!              h_value=0d0
!              do i=1, splen
!                 if (hicoo(i)==n.and.hjcoo(i)==m.and.hirv(1, i)==i1.and.hirv(2, i)==i2.and.hirv(3, i)==i3)then
!                    h_value= h_value+ hacoo(i)
!                 endif
!              enddo !i
!              write(outfileindex, '(5I5, 2f12.6)')irvec(:, iR), n, m, h_value/eV2Hartree
!           enddo !m
!        enddo !n
!     enddo !iR
!     close(outfileindex)
   endif

   return
end subroutine readSparseHmnR

subroutine readsparse_overlap
   !> This subroutine reads the overlap matrix S between the non-orthogonal basis
   !> 2023 Dec. 4th
   !> by QSWU (wuquansheng@gmail.com)
   use para
   implicit none
   integer:: i,j,nwann,nwann_nsoc,i1,i2,i3,i4,i5,ir,n,m 
   real(dp) :: r1, r2
   logical :: exists

   if(cpuid.eq.0)write(stdout,*)' '
   if(cpuid.eq.0)write(stdout,'(2a)')' >> Now we are reading the overlap files: ', trim(adjustl(Overlapfile))
   inquire(file= Overlapfile, EXIST = exists)
   if (exists)then
      open(13, file= Overlapfile, status='OLD')
   else
      STOP ">> for non-orthogonal basis , you have to prepare a file like overlap.dat"
   endif


   !> skip a comment line
   read(13, *)

   !> comparing with the standard overlap file, we add another line to show howmany
   !> lines that S is not zero.
!  if (Is_Sparse_Hr) then
   read(13,*) splen_overlap_input   !> number of non-zeros lines
!  end if

   !> number of Wannier orbitals in the hr file
   nwann=0
   read(13, *)nwann
   if (nwann==0.or.nwann.ne.Num_wann) then
      print *, 'in readsparse_overlap'
      print *, 'nwann, num_wann', nwann, num_wann
      stop "ERROR : num_wann is zero or num_wann is not consistent with hr.dat"
   endif
   nwann_nsoc=nwann
   if (SOC>0) nwann_nsoc= nwann/2

   !> overlap matrix in sparse format
   allocate(sacoo(splen_overlap_input))
   allocate(sicoo(splen_overlap_input))
   allocate(sjcoo(splen_overlap_input))
   allocate(sirv(3, splen_overlap_input))
   sacoo=(0d0, 0d0)
   sicoo=0
   sjcoo=0
   sirv=0

   !> will reread the above line
   do i=1, splen_overlap_input
      read(13,*,end=1003)i1, i2, i3, i4, i5, r1, r2
      sirv(1, i)= i1
      sirv(2, i)= i2
      sirv(3, i)= i3
      sicoo(i)=i4
      sjcoo(i)=i5
      sacoo(i)=dcmplx(r1,r2)
   enddo
   1003 continue
   close(13)

   if (cpuid.eq.0) write(stdout, '(a, i20)')' >> Number of non-zeros splen_overlap_input', splen_overlap_input
   if (cpuid.eq.0) write(stdout, '(a)')' >> Sparse overlap matrix reading finished '

   return
end subroutine readsparse_overlap



subroutine readsparse_valley_operator
   !> This subroutine not just read the sparse valley operator
   !> 2023 Nov. 4
   use para
   implicit none
   integer:: i,j,nwann,nwann_nsoc,i1,i2,i3,i4,i5,ir,n,m 
   real(dp) :: r1, r2
   logical :: exists

   !> the direction which adding electric field which is also the stacking direction
   integer :: add_electric_field
   real(dp) :: Bx_in_au, By_in_au, Bz_in_au
   complex(dp) :: h_value

   if(cpuid.eq.0)write(stdout,*)' '
   inquire (file ="valley_operator.dat", EXIST = exists)
   if (exists)then
      open(13, file="valley_operator.dat", status='OLD')
   else
      STOP ">> for valley projection , you have to prepare a file valley_operator.dat"
   endif


   !> skip a comment line
   read(13, *)

   !> comparing with the standard hr file, we add another line to show howmany
   !> lines that Hmn(R) is not zero.
   if(Is_Sparse_Hr) then
      read(13,*) splen_valley_input   !> number of non-zeros lines
   end if

   !> number of Wannier orbitals in the hr file
   nwann=0
   read(13, *)nwann
   if (nwann==0.or.nwann.ne.Num_wann) then
      print *, 'nwann, num_wann', nwann, num_wann
      stop "ERROR : num_wann is zero or num_wann is not consistent with hr.dat"
   endif
   nwann_nsoc=nwann
   if (SOC>0) nwann_nsoc= nwann/2

   !> number of lattice vectors taken into account
   read(13, *)Nrpts_valley
   allocate(irvec_valley(3, Nrpts_valley))

   !> The degeneracy of each R point
   read(13, *)(j, i=1, nrpts)
   ir=0

   allocate(valley_operator_acoo(splen_valley_input))
   allocate(valley_operator_icoo(splen_valley_input))
   allocate(valley_operator_jcoo(splen_valley_input))
   allocate(valley_operator_irv(splen_valley_input))
   valley_operator_acoo=(0d0, 0d0)
   valley_operator_icoo=0
   valley_operator_jcoo=0
   valley_operator_irv=0

   j=0
   ir=1
   irvec_valley=-9999
   read(13,*,end=1002)i1, i2, i3, i4, i5, r1, r2
   irvec_valley(1,ir)=i1
   irvec_valley(2,ir)=i2
   irvec_valley(3,ir)=i3
   !> will reread the above line
   backspace(13)
   do i=1,Nrpts_valley
      do n=1,nwann
         do m=1,nwann
            read(13,*,end=1002)i1, i2, i3, i4, i5, r1, r2
            if (sum(abs(irvec_valley(:,ir)-[i1,i2,i3]))/=0) ir=ir+1
            j=j+1
            valley_operator_icoo(j)=i4
            valley_operator_jcoo(j)=i5
            valley_operator_irv (j)=ir
            valley_operator_acoo(j)=dcmplx(r1,r2)
            irvec_valley(1, ir)= i1
            irvec_valley(2, ir)= i2
            irvec_valley(3, ir)= i3
         end do
      enddo
   enddo
1002 continue
   close(13)
   !> correct nrpts
   Nrpts_valley=ir

   if (cpuid.eq.0) write(stdout, '(a, i6)')' >> Nrpts_valley is ', Nrpts_valley
   if (cpuid.eq.0) write(stdout, '(a, 2i10)')' >> splen_valley_input', splen_valley_input, j
   if (cpuid.eq.0) write(stdout, '(a)')' >> Valley operator reading finished '

   return
end subroutine readsparse_valley_operator


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
   real(dp), allocatable :: pos_selected_for_center(:, :)

   integer :: i, ia, NumberofSelectedAtoms_center, iter, ia_g, ig
   real(dp) :: dxmax(3), center

   !> sum over all selected atoms
   NumberofSelectedAtoms_center= sum(NumberofSelectedAtoms(:))
   allocate(pos_selected_for_center(3, NumberofSelectedAtoms_center))
   pos_selected_for_center= 0d0

   iter= 0
   do ig=1,  NumberofSelectedAtoms_groups
      do ia_g= 1, NumberofSelectedAtoms(ig)
         iter=iter+1
         ia= Selected_Atoms(ig)%iarray(ia_g)
         pos_selected_for_center(:, iter)= Origin_cell%Atom_position_direct(:, ia)
      enddo
   enddo
   if (cpuid.eq.0)then
      write(stdout, '(a, 3f12.6)') ">>: Center position for zero electrical potential for adding electrical field", &
         sum(pos_selected_for_center(:, :), dim=2)/NumberofSelectedAtoms_center
   endif

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
      pos= mod(pos, 1d0)-0.5d0
      pos_selected_for_center= mod(pos_selected_for_center, 1d0)-0.5d0
      if (sum(center_atom_for_electric_field)>0) then
         center= (Origin_cell%Atom_position_direct(add_electric_field, center_atom_for_electric_field(1)) &
                + Origin_cell%Atom_position_direct(add_electric_field, center_atom_for_electric_field(2)))/2d0
         center= mod(center, 1d0)-0.5d0
      else
        !center= (maxval(pos)+minval(pos))/2d0
         center= sum(pos_selected_for_center(add_electric_field, :))/NumberofSelectedAtoms_center
      endif
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
      hirv (1, j)= 0
      hirv (2, j)= 0
      hirv (3, j)= 0
      hacoo(j)= Zeeman_energy_in_hartree_factor/2d0*Bz_in_au
   enddo
   do i=1, nwann_nsoc
      j= j+ 1
      hicoo(j)= i+ nwann_nsoc
      hjcoo(j)= i+ nwann_nsoc
      hirv (1, j)= 0
      hirv (2, j)= 0
      hirv (3, j)= 0
      hacoo(j)= -Zeeman_energy_in_hartree_factor/2d0*Bz_in_au
   enddo
   do i=1, nwann_nsoc
      j= j+ 1
      hicoo(j)= i
      hjcoo(j)= i+ nwann_nsoc
      hirv (1, j)= 0
      hirv (2, j)= 0
      hirv (3, j)= 0
      hacoo(j)=  Zeeman_energy_in_hartree_factor/2d0*(By_in_au*zi+ Bx_in_au)
   enddo
   do i=1, nwann_nsoc
      j= j+ 1
      hicoo(j)= i+ nwann_nsoc
      hjcoo(j)= i
      hirv (1, j)= 0
      hirv (2, j)= 0
      hirv (3, j)= 0
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





 subroutine add_zeeman_normal_hr()
   !> add Zeeman energy on the sparse hmnr based on Zeeman_energy_in_eV and the magnetic field direction
   !> Bx=Bdirection(1)
   !> By=Bdirection(2)
   !> Bz=Bdirection(3)
   !> magnetic_field_in_tesla is in Tesla
   !> Here the coordinate system is the same as user-specified in the LATTICE card.
   !> H_zeeman
   use para
   implicit none

   real(dp) :: Bdirection_x, Bdirection_y, Bdirection_z

   integer :: nwann_nsoc, n, m, ir0, ir
   real(dp) :: Zeeman_energy_in_au

   if (.not.Add_Zeeman_Field) return

   nwann_nsoc= Num_wann/2

   ir0=0
   do ir=1, Nrpts
      if (irvec(1,ir)==0.and.irvec(2,ir)==0.and.irvec(3,ir)==0) then
         ir0= ir
      endif
   enddo
   if (ir0==0) stop 'something wrong with irvec in subroutine add_zeeman_normal_hr'

   if (abs(Zeeman_energy_in_eV)<eps9)then
      !> in atomic unit
      Zeeman_energy_in_au= Effective_gfactor*Bohr_magneton*Bmagnitude
   else
      Zeeman_energy_in_au= Zeeman_energy_in_eV*eV2Hartree
   endif
   if (cpuid==0)then
      write(stdout, '(1x, a, 1f16.6, a)')" >>> Zeeman_energy_in_au used in WT: ",  Zeeman_energy_in_au, ' Hartree'
      write(stdout, '(1x, a, 1f16.6, a)')" >>> Zeeman_energy_in_eV used in WT: ",  Zeeman_energy_in_au/eV2Hartree, ' eV'
      write(stdout, '(1x, a, 3f16.3, a)')" >>> Bdirection: ",  Bdirection 
   endif

   if (SOC==0) then
      !> We need to extend the spin up to the spin down part
      do n=1, nwann_nsoc
         HmnR(n, n, ir0)            = HmnR(n, n, ir0)+ Zeeman_energy_in_au/2d0*Bdirection(3)
         HmnR(n+nwann_nsoc, n+nwann_nsoc, ir0)= HmnR(n, n, ir0)- Zeeman_energy_in_au/2d0*Bdirection(3)
         HmnR(n      , n+nwann_nsoc, ir0)= Zeeman_energy_in_au/2d0*Bdirection(1) &
                                         - zi*Zeeman_energy_in_au/2d0*Bdirection(2)
         HmnR(n+nwann_nsoc, n      , ir0)= Zeeman_energy_in_au/2d0*Bdirection(1) &
                                         + zi*Zeeman_energy_in_au/2d0*Bdirection(2)
      enddo
 
      do ir=1, Nrpts
         do m=1, nwann_nsoc
            do n=1, nwann_nsoc
               HmnR(n+nwann_nsoc, m+nwann_nsoc, ir)= HmnR(n, m, ir)
            enddo ! m
         enddo ! n
      enddo ! ir
   else
      do n=1, nwann_nsoc
         HmnR(n, n, ir0)            = HmnR(n, n, ir0)+ Zeeman_energy_in_au/2d0*Bdirection(3)
         HmnR(n+nwann_nsoc, n+nwann_nsoc, ir0)= HmnR(n+nwann_nsoc, n+nwann_nsoc, ir0)- Zeeman_energy_in_au/2d0*Bdirection(3)
         HmnR(n      , n+nwann_nsoc, ir0)= HmnR(n      , n+nwann_nsoc, ir0)+ Zeeman_energy_in_au/2d0*Bdirection(1) &
                                         - zi*Zeeman_energy_in_au/2d0*Bdirection(2)
         HmnR(n+nwann_nsoc, n      , ir0)= HmnR(n+nwann_nsoc, n      , ir0)+ Zeeman_energy_in_au/2d0*Bdirection(1) &
                                         + zi*Zeeman_energy_in_au/2d0*Bdirection(2)
      enddo
   endif ! SOC

   return

 end subroutine add_zeeman_normal_hr
