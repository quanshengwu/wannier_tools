!------------------------------------------------------------------------------
! In this file, we are going to find the symmetry operations according to the 
! unit cell defined by lattice vectors Origin_cell%Rua, Origin_cell%Rub, Origin_cell%Ruc and atomic positions. The
! steps are as following:
! 1. define the basic point group operations
! 2. try to find the generators 
! 3. use the generators to get all the operations
! The operations are composed by the point group operations pgops
! and the translations tau. 
! April 12 2018 in EPFL 
! QuanSheng Wu (wuquansheng@gmail.com)
!------------------------------------------------------------------------------

  subroutine DefineBasicOperations(BasicOperations_space, BasicOperations_spin, BasicOperations_inversion)
     !> We adopt the definations from Z.Fang, There are 18 basic operations.
     !> The cartesian coordinates are used.
     !> E, I, 
     !> BasicOperations_space is the group operation in real space
     !> BasicOperations_spin is the group operators in spin space 
     !> BasicOperations_inversion is a indicator to check whether this operation
     !> contains inversion or not
     use para, only : dp, zi, pi
     implicit none

     real(dp), intent(out) :: BasicOperations_space(3, 3, 0:48)
     complex(dp), intent(out) :: BasicOperations_spin(2, 2, 0:48)
     real(dp), intent(out) :: BasicOperations_inversion(0:48)

     real(dp)    :: dcos_60, dsin_60, dcos_120, dsin_120
     complex(dp) :: z0, ze, z30, z_30, z45, z_45, z60, z_60, z90, z_90, &
                    z120, z_120, z180, z_180, cos_30, sin_30, &
                    cos_45, sin_45, cos_60, sin_60, cos_120, sin_120

     z0 = dcmplx(0.0d0, 0.0d0)
     ze = dcmplx(1.0d0, 0.0d0)
     z30  = dcmplx(dcos(pi/6.d0), dsin(pi/6.d0))
     z_30 = dcmplx(dcos(pi/6.d0),-dsin(pi/6.d0))
     z45  = dcmplx(dcos(pi/4.d0), dsin(pi/4.d0))
     z_45 = dcmplx(dcos(pi/4.d0),-dsin(pi/4.d0))
     z60  = dcmplx(dcos(pi/3.d0), dsin(pi/3.d0))
     z_60 = dcmplx(dcos(pi/3.d0),-dsin(pi/3.d0))
     z90  = dcmplx(dcos(pi/2.d0), dsin(pi/2.d0))
     z_90 = dcmplx(dcos(pi/2.d0),-dsin(pi/2.d0))
     z120 = dcmplx(dcos(2.d0*pi/3.d0), dsin(2.d0*pi/3.d0))
     z_120= dcmplx(dcos(2.d0*pi/3.d0),-dsin(2.d0*pi/3.d0))
     z180 = dcmplx(dcos(pi/1.d0), dsin(pi/1.d0))
     z_180= dcmplx(dcos(pi/1.d0),-dsin(pi/1.d0))
     cos_30=dcmplx(dcos(pi/6.d0),0.d0)
     sin_30=dcmplx(dsin(pi/6.d0),0.d0)
     cos_45=dcmplx(dcos(pi/4.d0),0.d0)
     sin_45=dcmplx(dsin(pi/4.d0),0.d0)
     cos_60=dcmplx(dcos(pi/3.d0),0.d0)
     sin_60=dcmplx(dsin(pi/3.d0),0.d0)
     cos_120=dcmplx(dcos(2.d0*pi/3.d0),0.d0)
     sin_120=dcmplx(dsin(2.d0*pi/3.d0),0.d0)
  
     dcos_60=dcos(pi/3.d0)
     dsin_60=dsin(pi/3.d0)
     dcos_120=dcos(2.d0*pi/3.d0)
     dsin_120=dsin(2.d0*pi/3.d0)

     BasicOperations_space=0.d0

     !--    No. 0 = 0
     BasicOperations_space(1,:,0) = (/ 0.0d0 ,   0.0d0 ,   0.0d0  /)  
     BasicOperations_space(2,:,0) = (/ 0.0d0 ,   0.0d0 ,   0.0d0  /)
     BasicOperations_space(3,:,0) = (/ 0.0d0 ,   0.0d0 ,   0.0d0  /)
     BasicOperations_spin(1,:,0) = (/ z0 ,  z0 /)
     BasicOperations_spin(2,:,0) = (/ z0 ,  z0 /)
     BasicOperations_inversion(0)= 1

     !--    No. 1 = E
     BasicOperations_space(1,:,1) = (/ 1.0d0 ,   0.0d0 ,   0.0d0  /)  
     BasicOperations_space(2,:,1) = (/ 0.0d0 ,   1.0d0 ,   0.0d0  /)
     BasicOperations_space(3,:,1) = (/ 0.0d0 ,   0.0d0 ,   1.0d0  /)
     BasicOperations_spin(1,:,1) = (/ ze ,  z0 /)
     BasicOperations_spin(2,:,1) = (/ z0 ,  ze /)
     BasicOperations_inversion(1)= 1

     !--    No. 2 = I
     BasicOperations_space(1,:,2) = (/-1.0d0 ,   0.0d0 ,   0.0d0  /)  
     BasicOperations_space(2,:,2) = (/ 0.0d0 ,  -1.0d0 ,   0.0d0  /)
     BasicOperations_space(3,:,2) = (/ 0.0d0 ,   0.0d0 ,  -1.0d0  /)
     BasicOperations_spin(1,:,2) = (/ ze ,  z0 /)
     BasicOperations_spin(2,:,2) = (/ z0 ,  ze /)
     BasicOperations_inversion(2)=-1

     !--    No. 3 = R4X
     BasicOperations_space(1,:,3) = (/ 1.0d0 ,   0.0d0 ,   0.0d0  /) 
     BasicOperations_space(2,:,3) = (/ 0.0d0 ,   0.0d0 ,  -1.0d0  /)
     BasicOperations_space(3,:,3) = (/ 0.0d0 ,   1.0d0 ,   0.0d0  /)
     BasicOperations_spin(1,:,3) = (/  ze*cos_45,  zi*sin_45 /)
     BasicOperations_spin(2,:,3) = (/  zi*sin_45,  ze*cos_45 /)
     BasicOperations_inversion(3)= 1

     !--    No. 4 = R4Y
     BasicOperations_space(1,:,4) = (/ 0.0d0 ,   0.0d0 ,   1.0d0  /)  
     BasicOperations_space(2,:,4) = (/ 0.0d0 ,   1.0d0 ,   0.0d0  /)
     BasicOperations_space(3,:,4) = (/-1.0d0 ,   0.0d0 ,   0.0d0  /)
     BasicOperations_spin(1,:,4) = (/ ze*dcos(pi/4.d0), ze*dsin(pi/4.d0) /)
     BasicOperations_spin(2,:,4) = (/-ze*dsin(pi/4.d0), ze*dcos(pi/4.d0) /)
     BasicOperations_inversion(4)= 1

     !--    No. 5 = R4Z
     BasicOperations_space(1,:,5) = (/ 0.0d0 ,  -1.0d0 ,   0.0d0  /)  
     BasicOperations_space(2,:,5) = (/ 1.0d0 ,   0.0d0 ,   0.0d0  /)
     BasicOperations_space(3,:,5) = (/ 0.0d0 ,   0.0d0 ,   1.0d0  /)
     BasicOperations_spin(1,:,5) = (/  z45, z0 /)
     BasicOperations_spin(2,:,5) = (/ z0,  z_45 /)
     BasicOperations_inversion(5)= 1

     !--    No. 6 = R2X
     BasicOperations_space(1,:,6) = (/ 1.0d0 ,   0.0d0 ,   0.0d0  /)  
     BasicOperations_space(2,:,6) = (/ 0.0d0 ,  -1.0d0 ,   0.0d0  /)
     BasicOperations_space(3,:,6) = (/ 0.0d0 ,   0.0d0 ,  -1.0d0  /)
     BasicOperations_spin(1,:,6) = (/  z0, zi /)
     BasicOperations_spin(2,:,6) = (/  zi, z0 /)
     BasicOperations_inversion(6)= 1
 
     !--    No. 7 = R2Y
     BasicOperations_space(1,:,7) = (/-1.0d0 ,   0.0d0 ,   0.0d0  /)  
     BasicOperations_space(2,:,7) = (/ 0.0d0 ,   1.0d0 ,   0.0d0  /)
     BasicOperations_space(3,:,7) = (/ 0.0d0 ,   0.0d0 ,  -1.0d0  /)
     BasicOperations_spin(1,:,7) = (/  z0, ze /)
     BasicOperations_spin(2,:,7) = (/ -ze, z0 /)
     BasicOperations_inversion(7)= 1

     !--    No. 8 = R2Z
     BasicOperations_space(1,:,8) = (/-1.0d0 ,   0.0d0 ,   0.0d0  /)  
     BasicOperations_space(2,:,8) = (/ 0.0d0 ,  -1.0d0 ,   0.0d0  /)
     BasicOperations_space(3,:,8) = (/ 0.0d0 ,   0.0d0 ,   1.0d0  /)
     BasicOperations_spin(1,:,8) = (/ zi, z0 /)
     BasicOperations_spin(2,:,8) = (/ z0, -zi /)
     BasicOperations_inversion(8)= 1

     !--    No. 9 = R3XYZ
     BasicOperations_space(1,:,9) = (/ 0.0d0 ,   0.0d0 ,   1.0d0  /)
     BasicOperations_space(2,:,9) = (/ 1.0d0 ,   0.0d0 ,   0.0d0  /)
     BasicOperations_space(3,:,9) = (/ 0.0d0 ,   1.0d0 ,   0.0d0  /)
     BasicOperations_spin(1,:,9)=(/dcmplx(0.5d0,0.5d0),dcmplx(0.5d0,0.5d0)/)
     BasicOperations_spin(2,:,9)=(/-dcmplx(0.5d0,-0.5d0),dcmplx(0.5d0,-0.5d0)/)
     BasicOperations_inversion(9)= 1

     !--    No. 10 = mirror_x
     BasicOperations_space(1,:,10) = (/-1.0d0 ,   0.0d0 ,   0.0d0  /)  
     BasicOperations_space(2,:,10) = (/ 0.0d0 ,   1.0d0 ,   0.0d0  /)
     BasicOperations_space(3,:,10) = (/ 0.0d0 ,   0.0d0 ,   1.0d0  /)
     BasicOperations_spin(1,:,10) = (/ z0,zi /)
     BasicOperations_spin(2,:,10) = (/ zi, z0 /)
     BasicOperations_inversion(10)=-1

     !--    No. 11 = mirror_y
     BasicOperations_space(1,:,11) = (/ 1.0d0 ,   0.0d0 ,   0.0d0  /)  
     BasicOperations_space(2,:,11) = (/ 0.0d0 ,  -1.0d0 ,   0.0d0  /)
     BasicOperations_space(3,:,11) = (/ 0.0d0 ,   0.0d0 ,   1.0d0  /)
     BasicOperations_spin(1,:,11) = (/  z0,  ze /)
     BasicOperations_spin(2,:,11) = (/ -ze, z0 /)
     BasicOperations_inversion(11)=-1

     !--    No. 12 = mirror_z
     BasicOperations_space(1,:,12) = (/ 1.0d0 ,   0.0d0 ,   0.0d0  /)  
     BasicOperations_space(2,:,12) = (/ 0.0d0 ,   1.0d0 ,   0.0d0  /)
     BasicOperations_space(3,:,12) = (/ 0.0d0 ,   0.0d0 ,  -1.0d0  /)
     BasicOperations_spin(1,:,12) = (/ zi, z0 /)
     BasicOperations_spin(2,:,12) = (/ z0,-zi /)
     BasicOperations_inversion(12)=-1

     !--    No. 13 = mirror_xy
     BasicOperations_space(1,:,13) = (/ 0.0d0 ,   1.0d0 ,   0.0d0  /)  
     BasicOperations_space(2,:,13) = (/ 1.0d0 ,   0.0d0 ,   0.0d0  /)
     BasicOperations_space(3,:,13) = (/ 0.0d0 ,   0.0d0 ,   1.0d0  /)
     BasicOperations_spin(1,:,13) = (/ z0, z45 /)
     BasicOperations_spin(2,:,13) = (/-z_45, z0 /)
     BasicOperations_inversion(13)=-1

     !--    No. 14 = S4 (z)
     BasicOperations_space(1,:,14) = (/ 0.0d0 ,  -1.0d0 ,   0.0d0  /)  
     BasicOperations_space(2,:,14) = (/ 1.0d0 ,   0.0d0 ,   0.0d0  /)
     BasicOperations_space(3,:,14) = (/ 0.0d0 ,   0.0d0 ,  -1.0d0  /)
     BasicOperations_spin(1,:,14) = (/ z_45, z0 /)
     BasicOperations_spin(2,:,14) = (/ z0, z45 /)
     BasicOperations_inversion(14)=-1

     !--    No. 15 = R2 (xy)
     BasicOperations_space(1,:,15) = (/ 0.0d0 ,   1.0d0 ,   0.0d0  /)  
     BasicOperations_space(2,:,15) = (/ 1.0d0 ,   0.0d0 ,   0.0d0  /)
     BasicOperations_space(3,:,15) = (/ 0.0d0 ,   0.0d0 ,  -1.0d0  /)
     BasicOperations_spin(1,:,15) = (/ z0, z_45 /)
     BasicOperations_spin(2,:,15) = (/-z45, z0 /)
     BasicOperations_inversion(15)=+1

     !--    No. 16 = R2 (-xy)
     BasicOperations_space(1,:,16) = (/ 0.0d0 ,  -1.0d0 ,   0.0d0  /)  
     BasicOperations_space(2,:,16) = (/-1.0d0 ,   0.0d0 ,   0.0d0  /)
     BasicOperations_space(3,:,16) = (/ 0.0d0 ,   0.0d0 ,  -1.0d0  /)
     BasicOperations_spin(1,:,16) = (/ z0, z45 /)
     BasicOperations_spin(2,:,16) = (/-z_45, z0 /)
     BasicOperations_inversion(16)=+1

     !--    No. 17 = R3Z
     BasicOperations_space(1,:,17)  = (/ dcos_120,-dsin_120, 0.0d0  /)
     BasicOperations_space(2,:,17)  = (/ dsin_120, dcos_120,   0.0d0  /)
     BasicOperations_space(3,:,17)  = (/ 0.0d0 ,   0.0d0 ,   1.0d0  /)
     BasicOperations_spin(1,:,17) = (/  z60, z0 /)
     BasicOperations_spin(2,:,17) = (/ z0,  z_60 /)
     BasicOperations_inversion(17)= 1

     !--    No. 18 = R6Z
     BasicOperations_space(1,:,18) = (/ dcos_60, -dsin_60,   0.0d0  /)
     BasicOperations_space(2,:,18) = (/ dsin_60,  dcos_60,   0.0d0  /)
     BasicOperations_space(3,:,18) = (/ 0.0d0 ,   0.0d0 ,   1.0d0  /)
     BasicOperations_spin(1,:,18) = (/  z30, z0 /)
     BasicOperations_spin(2,:,18) = (/ z0,  z_30 /)
     BasicOperations_inversion(18)= 1

     return
  end subroutine DefineBasicOperations

  subroutine FindGenerators(BasicOperations_space, BasicOperations_inversion, num_gen, generators_find, tau_find )
     !> In this subroutine, we are going to find the group generators with the given crystal stOrigin_cell%Ructure
     !> Here we assume that we have time reversal symmetry.
     
#if defined (MPI)
     use para, only : dp, Origin_cell, mpi_cmw, mpi_dp, mpi_in, mpi_sum, &
#else
     use para, only : dp, Origin_cell, &
#endif
                      cpuid, num_cpu, stdout, &
                      spatial_inversion

     implicit none

     !> input : 18 basic operators, and the lattice vectors Origin_cell%Rua, Origin_cell%Rub, Origin_cell%Ruc, and atomic positions
     real(dp), intent(in)  :: BasicOperations_space(3,3,0:48)
     real(dp), intent(in)  :: BasicOperations_inversion(0:48)

     !> output : number of generators and the generators itself. The generators are
     !> stored as integers which is the indices of BasicOperations_space
     !> Only generators_find(1:num_gen) and tau_find(1:3, 1:num_gen) are meaningful
     integer , intent(out) :: num_gen
     integer , intent(out) :: generators_find(48)
     real(dp), intent(out) :: tau_find(3,48)

     integer  :: i,j,k,ncell,ia,ia2,ic,na,it,nt,nop,ierr
     integer  :: natom_typ(Origin_cell%Num_atom_type),natom_min
     integer  :: typ_min, ntaup(48), ntaup_mpi(48)
     real(dp) :: temp,del
     real(dp) :: pos1(3)
     real(dp) :: Tmat(3,3),Tmat_inv(3,3)
     real(dp), external :: det3
     
     integer, allocatable :: typ_all(:, :)
     integer, allocatable :: op_tau(:, :), op_tau_mpi(:, :)

     real(dp), allocatable :: pos_new(:, :), pos_new2(:, :)
     real(dp), allocatable :: op_in(:, :, :)
     real(dp), allocatable :: tau_p(:, :), taut(:, :)
     real(dp), allocatable :: pos_min(:, :), pos_typ(:, :, :)
     real(dp), allocatable :: pos_all(:, :, :), post(:, :)
     real(dp), allocatable :: tau_find_t(:, :, :), tau_find_t_mpi(:, :, :)

     allocate(typ_all(Origin_cell%Num_atoms, 2197))
     allocate(pos_new(8*Origin_cell%Num_atoms, 3), pos_new2(8*Origin_cell%Num_atoms, 3))
     allocate(op_tau(48, 50), op_tau_mpi(48, 50))
     allocate(tau_find_t(3, 48, 50), tau_find_t_mpi(3, 48, 50))
     allocate(op_in(3, 3, 0:48))
     allocate(tau_p(3, 2197*Origin_cell%Num_atoms), taut(3, 2197*Origin_cell%Num_atoms))
     allocate(pos_min(3, 2197*Origin_cell%Num_atoms), pos_typ(3, 2197*Origin_cell%Num_atoms, Origin_cell%Num_atom_type))
     allocate(pos_all(Origin_cell%Num_atoms, 3, 2197), post(Origin_cell%Num_atoms, 3))
     del=1.d-3
     op_tau=0
     op_tau_mpi=0
     ntaup=1
     ntaup_mpi=0
     generators_find=1
     tau_find=0.d0
     tau_find_t=0.d0
     tau_find_t_mpi=0.d0
     taut=0.d0
     tau_p=0.d0
     if(cpuid.eq.0)write(stdout,*)'   '
     if(cpuid.eq.0)write(stdout,*)'   <<  Symmetry Finding starts >> '

     !-----------------------------------------------------------------------------
     !>  To get all the atom position in the 3x3 super-cell 
     !>  post: Origin_cell%Atom_position_direct for all the Origin_cell%Num_atoms atoms, rather than only Origin_cell%Num_atoms atoms
     !> typ_all: the atom type for all atoms in the 3x3 cell 
     !> natom_typ(Origin_cell%Num_atom_type): the total number of atoms for each type in the 3x3 cell
     !> pos_typ: the Origin_cell%Atom_position_direct for each type of atoms
     !-----------------------------------------------------------------------------
     na=0 

     !> shift the atom's position into the home unit cell [0,1)
     !> position in unit of Origin_cell%Rua, Origin_cell%Rub, Origin_cell%Ruc
     do ia=1,Origin_cell%Num_atoms
        post(ia,1)=Origin_cell%Atom_position_direct(1,ia)+1.d0-int(Origin_cell%Atom_position_direct(1,ia)+1.000001d0)
        post(ia,2)=Origin_cell%Atom_position_direct(2,ia)+1.d0-int(Origin_cell%Atom_position_direct(2,ia)+1.000001d0)
        post(ia,3)=Origin_cell%Atom_position_direct(3,ia)+1.d0-int(Origin_cell%Atom_position_direct(3,ia)+1.000001d0)
     enddo

     ncell=0 
     natom_typ=0
     do i=-6,6
     do j=-6,6
     do k=-6,6
       ncell=ncell+1
       do ia=1,Origin_cell%Num_atoms
         pos_all(ia,1,ncell)=post(ia,1)+dble(i)
         pos_all(ia,2,ncell)=post(ia,2)+dble(j)
         pos_all(ia,3,ncell)=post(ia,3)+dble(k)
         it=Origin_cell%itype_atom(ia)
         typ_all(ia,ncell)=it
         natom_typ(it)=natom_typ(it)+1
         nt=natom_typ(it)
         pos_typ(:,nt,it)=pos_all(ia,:,ncell)
       enddo
     enddo
     enddo
     enddo

     if(cpuid.eq.0)then
          write(stdout,*)'       Type         num_atom    '
        do it=1,Origin_cell%Num_atom_type
          write(stdout,'(5x,2I10)')it,natom_typ(it)
        enddo
     end if

     !>  To find the atom type, which has minimum number of atoms --
     natom_min=2197*Origin_cell%Num_atoms
     do it=1,Origin_cell%Num_atom_type
        if (natom_typ(it).le.natom_min)then
           natom_min=natom_typ(it)
           typ_min=it
        end if
     enddo
     if(cpuid.eq.0)write(stdout,*)'      typ_min, natom_min=',typ_min,natom_min


     !-----------------------------------------------------------------------------
     !>  To get all the symmetry operations for abc basis (primitive cell) 
     !>  BasicOperations_space:  input operation under cardisian coordiants
     !>  op_in:   the operation under abc basis  
     !-----------------------------------------------------------------------------
     Tmat(:,1)= Origin_cell%Rua
     Tmat(:,2)= Origin_cell%Rub
     Tmat(:,3)= Origin_cell%Ruc
     Tmat_inv= Tmat
     call inv_r(3, Tmat_inv)
     do nop=1,48
       op_in(:,:,nop)=matmul(Tmat_inv,matmul(BasicOperations_space(:,:,nop),Tmat))
     end do

     !> Start the symmetry operation on atoms 
     DO 100 nop=cpuid+1,48,num_cpu
        if(dabs(det3(op_in(:,:,nop))).lt.del) goto 100
        !>  Rotation on Origin_cell%Atom_position_direct 
        ncell=0
        do i=0,1
        do j=0,1
        do k=0,1
        do ia=1,Origin_cell%Num_atoms
          pos1(1)=post(ia,1)+dble(i)
          pos1(2)=post(ia,2)+dble(j)
          pos1(3)=post(ia,3)+dble(k)
          pos_new(ia+ncell*Origin_cell%Num_atoms,:)=MATMUL(op_in(:,:,nop),pos1(:))
        end do
        ncell=ncell+1
        end do
        end do
        end do
        !> to find all the possible tau 
        ntaup(nop)=1
        DO ia2=1,Origin_cell%Num_atoms
           IF(Origin_cell%itype_atom(ia2).eq.typ_min)then
              pos1(:)=pos_new(ia2,:)
            
              do ia=1,natom_min
                tau_p(:,ia)=pos_typ(:,ia,typ_min)-pos1(:)
                tau_p(1,ia)=tau_p(1,ia)+10.d0 -int(tau_p(1,ia)+10.000001d0)
                tau_p(2,ia)=tau_p(2,ia)+10.d0 -int(tau_p(2,ia)+10.000001d0)
                tau_p(3,ia)=tau_p(3,ia)+10.d0 -int(tau_p(3,ia)+10.000001d0)
              end do
            
              do 101 ia=1,natom_min
                do i=1,ntaup(nop)
                  temp=dsqrt( (tau_p(1,ia)-taut(1,i))**2  &
                                  +(tau_p(2,ia)-taut(2,i))**2  &
                                  +(tau_p(3,ia)-taut(3,i))**2 )
                  if(temp.lt.del.or.ntaup(nop).ge.50) go to 101
                end do
                ntaup(nop)=ntaup(nop)+1
                taut(:,ntaup(nop))=tau_p(:,ia)
              101 continue
           END IF
        END DO
       
        !>  To check the symmetry for each taut 
        do 104 nt=1,min(50,ntaup(nop))
           do ia=1,8*Origin_cell%Num_atoms
             pos_new2(ia,1)=pos_new(ia,1)+taut(1,nt)
             pos_new2(ia,2)=pos_new(ia,2)+taut(2,nt)
             pos_new2(ia,3)=pos_new(ia,3)+taut(3,nt)
           end do
          
           do ia=1,8*Origin_cell%Num_atoms
              it=Origin_cell%itype_atom(mod(ia+Origin_cell%Num_atoms-1,Origin_cell%Num_atoms)+1)
              do ia2=1,natom_typ(it)
                 temp=dsqrt( (pos_new2(ia,1)-pos_typ(1,ia2,it))**2  &
                                 +(pos_new2(ia,2)-pos_typ(2,ia2,it))**2  &
                                 +(pos_new2(ia,3)-pos_typ(3,ia2,it))**2 )
                 if(temp.lt.del) go to 103
                 if(ia2.eq.natom_typ(it)) go to 104
              end do
              103 continue
              if (ia.eq.8*Origin_cell%Num_atoms) then
                 op_tau(nop,nt)=1
                 tau_find_t(:,nop,nt)=taut(:,nt)
              end if
           end do
        104 continue
     100 CONTINUE

#if defined (MPI)
     call mpi_barrier(mpi_cmw,ierr)
     call mpi_allreduce(tau_find_t,tau_find_t_mpi,size(tau_find_t),&
                       mpi_dp,mpi_sum,mpi_cmw,ierr)
     call mpi_allreduce(op_tau,op_tau_mpi,size(op_tau),&
                       mpi_in,mpi_sum,mpi_cmw,ierr)
     call mpi_allreduce(ntaup,ntaup_mpi,size(ntaup),&
                       mpi_in,mpi_sum,mpi_cmw,ierr)
     call mpi_barrier(mpi_cmw,ierr)
#else
     ntaup_mpi= ntaup
     op_tau_mpi= op_tau
     tau_find_t_mpi= tau_find_t
#endif

     !-----------------------------------------------------------------------------
     !>  Transfer tau: from abc basis (primitive cell)
     !>                to conventioinal cell
     !-----------------------------------------------------------------------------
     if (cpuid.eq.0) write(stdout, '(3X, a)')'>> Here are the generators we found:'
     num_gen=0
     DO 120 nop=1,18
     IF (SUM(op_tau_mpi(nop,:)).ge.1)then
        if (cpuid.eq.0)then
           write(stdout, '(a, f5.2)')'    Inversion       ', BasicOperations_inversion(nop)
           SELECT CASE(nop)
              CASE(1)
                 write(stdout,*)'   nop  = 1    E '
                 write(stdout,'(7x,3f12.6)')BasicOperations_space(1,:,1)
                 write(stdout,'(7x,3f12.6)')BasicOperations_space(2,:,1)
                 write(stdout,'(7x,3f12.6)')BasicOperations_space(3,:,1)
              CASE(2)
                 write(stdout,*)'   nop  = 2    I '
                 write(stdout,'(7x,3f12.6)')BasicOperations_space(1,:,2)
                 write(stdout,'(7x,3f12.6)')BasicOperations_space(2,:,2)
                 write(stdout,'(7x,3f12.6)')BasicOperations_space(3,:,2)
              CASE(3)
                 write(stdout,*)'   nop  = 3    R4X '
                 write(stdout,'(7x,3f12.6)')BasicOperations_space(1,:,3)
                 write(stdout,'(7x,3f12.6)')BasicOperations_space(2,:,3)
                 write(stdout,'(7x,3f12.6)')BasicOperations_space(3,:,3)
              CASE(4)
                 write(stdout,*)'   nop  = 4    R4Y '
                 write(stdout,'(7x,3f12.6)')BasicOperations_space(1,:,4)
                 write(stdout,'(7x,3f12.6)')BasicOperations_space(2,:,4)
                 write(stdout,'(7x,3f12.6)')BasicOperations_space(3,:,4)
              CASE(5)
                 write(stdout,*)'   nop  = 5    R4Z '
                 write(stdout,'(7x,3f12.6)')BasicOperations_space(1,:,5)
                 write(stdout,'(7x,3f12.6)')BasicOperations_space(2,:,5)
                 write(stdout,'(7x,3f12.6)')BasicOperations_space(3,:,5)
              CASE(6)
                 write(stdout,*)'   nop  = 6    R2X '
                 write(stdout,'(7x,3f12.6)')BasicOperations_space(1,:,6)
                 write(stdout,'(7x,3f12.6)')BasicOperations_space(2,:,6)
                 write(stdout,'(7x,3f12.6)')BasicOperations_space(3,:,6)
              CASE(7)
                 write(stdout,*)'   nop  = 7    R2Y '
                 write(stdout,'(7x,3f12.6)')BasicOperations_space(1,:,7)
                 write(stdout,'(7x,3f12.6)')BasicOperations_space(2,:,7)
                 write(stdout,'(7x,3f12.6)')BasicOperations_space(3,:,7)
              CASE(8)
                 write(stdout,*)'   nop  = 8    R2Z '
                 write(stdout,'(7x,3f12.6)')BasicOperations_space(1,:,8)
                 write(stdout,'(7x,3f12.6)')BasicOperations_space(2,:,8)
                 write(stdout,'(7x,3f12.6)')BasicOperations_space(3,:,8)
              CASE(9)
                 write(stdout,*)'   nop  = 9    R3XYZ '
                 write(stdout,'(7x,3f12.6)')BasicOperations_space(1,:,9)
                 write(stdout,'(7x,3f12.6)')BasicOperations_space(2,:,9)
                 write(stdout,'(7x,3f12.6)')BasicOperations_space(3,:,9)
              CASE(10)
                 write(stdout,*)'   nop  = 10    Mirror_x '
                 write(stdout,'(7x,3f12.6)')BasicOperations_space(1,:,10)
                 write(stdout,'(7x,3f12.6)')BasicOperations_space(2,:,10)
                 write(stdout,'(7x,3f12.6)')BasicOperations_space(3,:,10)
              CASE(11)
                 write(stdout,*)'   nop  = 11    Mirror_y '
                 write(stdout,'(7x,3f12.6)')BasicOperations_space(1,:,11)
                 write(stdout,'(7x,3f12.6)')BasicOperations_space(2,:,11)
                 write(stdout,'(7x,3f12.6)')BasicOperations_space(3,:,11)
              CASE(12)
                 write(stdout,*)'   nop  = 12    Mirror_z '
                 write(stdout,'(7x,3f12.6)')BasicOperations_space(1,:,12)
                 write(stdout,'(7x,3f12.6)')BasicOperations_space(2,:,12)
                 write(stdout,'(7x,3f12.6)')BasicOperations_space(3,:,12)
              CASE(13)
                 write(stdout,*)'   nop  = 13    Mirror_xy '
                 write(stdout,'(7x,3f12.6)')BasicOperations_space(1,:,13)
                 write(stdout,'(7x,3f12.6)')BasicOperations_space(2,:,13)
                 write(stdout,'(7x,3f12.6)')BasicOperations_space(3,:,13)
              CASE(14)
                 write(stdout,*)'   nop  = 14    S4Z '
                 write(stdout,'(7x,3f12.6)')BasicOperations_space(1,:,14)
                 write(stdout,'(7x,3f12.6)')BasicOperations_space(2,:,14)
                 write(stdout,'(7x,3f12.6)')BasicOperations_space(3,:,14)
              CASE(15)
                 write(stdout,*)'   nop  = 15    R2XY '
                 write(stdout,'(7x,3f12.6)')BasicOperations_space(1,:,15)
                 write(stdout,'(7x,3f12.6)')BasicOperations_space(2,:,15)
                 write(stdout,'(7x,3f12.6)')BasicOperations_space(3,:,15)
              CASE(16)
                 write(stdout,*)'   nop  = 16    R2-XY '
                 write(stdout,'(7x,3f12.6)')BasicOperations_space(1,:,16)
                 write(stdout,'(7x,3f12.6)')BasicOperations_space(2,:,16)
                 write(stdout,'(7x,3f12.6)')BasicOperations_space(3,:,16)
              CASE(17)
                 write(stdout,*)'   nop  = 17    R3Z '
                 write(stdout,'(7x,3f12.6)')BasicOperations_space(1,:,17)
                 write(stdout,'(7x,3f12.6)')BasicOperations_space(2,:,17)
                 write(stdout,'(7x,3f12.6)')BasicOperations_space(3,:,17)
              CASE(18)
                write(stdout,*)'   nop  = 18    R6Z '
                write(stdout,'(7x,3f12.6)')BasicOperations_space(1,:,18)
                write(stdout,'(7x,3f12.6)')BasicOperations_space(2,:,18)
                write(stdout,'(7x,3f12.6)')BasicOperations_space(3,:,18)
              CASE Default
                write(stdout,*)'   nop  = ',nop
                write(stdout,'(7x,3f12.6)')BasicOperations_space(1,:,nop)
                write(stdout,'(7x,3f12.6)')BasicOperations_space(2,:,nop)
                write(stdout,'(7x,3f12.6)')BasicOperations_space(3,:,nop)
           END SELECT
        end if
  
       num_gen=num_gen+1
       generators_find(num_gen)=nop
       spatial_inversion(num_gen) = BasicOperations_inversion(nop)
  
       do nt=1,min(50,ntaup_mpi(nop))
          if (op_tau_mpi(nop,nt).eq.1)then
             tau_find_t_mpi(:,nop,nt)=matmul(Tmat,tau_find_t_mpi(:,nop,nt))
             tau_find_t_mpi(:,nop,nt)=matmul(Tmat_inv,tau_find_t_mpi(:,nop,nt))
             tau_find(:,num_gen)=tau_find_t_mpi(:,nop,nt)
             if (cpuid.eq.0) then
                write(stdout,'(2x,"  Tau",3f12.6)')tau_find_t_mpi(:,nop,nt)
                write(stdout, *)' '
             end if
          end if
       end do

     END IF
     120 CONTINUE

#if defined (MPI)
     call mpi_barrier(mpi_cmw,ierr)
#endif

     return

  end subroutine FindGenerators


  subroutine FindSymmetryOperators(BasicOperations_space, BasicOperations_inversion)
     use para
     implicit none

     !> input : 18 basic operators, and the lattice vectors Origin_cell%Rua, Origin_cell%Rub, Origin_cell%Ruc, and atomic positions
     real(dp), intent(in)  :: BasicOperations_space(3,3,0:48)
     real(dp), intent(in)  :: BasicOperations_inversion(0:48)

     !> input : number of generators and the generators itself. The generators are
     !> stored as integers which is the indices of BasicOperations_space
     !> Only generators_find(1:num_gen) and tau_find(1:3, 1:num_gen) are meaningful


     integer :: i, j, k
     real(dp) :: tau_new(3), op_cart_new(3, 3)
     real(dp) :: altv(3, 3), TT(3, 3), T1(3, 3), T2(3, 3)

     altv(:, 1)= Origin_cell%Rua
     altv(:, 2)= Origin_cell%Rub
     altv(:, 3)= Origin_cell%Ruc
     TT= altv
     call inv_r(3, TT)


     do i=1, number_group_generators
        pggen_cart(:, :, i)= BasicOperations_space(:, :, generators_find(i))
        !> transform cartesian coordinates to direct/fractional coordinates
        pggen_direct(:, :, i)= matmul(TT, matmul(pggen_cart(:, :, i), altv))
        pgop_cart(:, :, i)= pggen_cart(:, :, i)
        pgop_direct(:, :, i)= pggen_direct(:, :, i)
        tau_direct(:, i)= tau_find(:, i)
     enddo

     !> use the generators to find all the space group operators
     number_group_operators= number_group_generators
     do i=1, number_group_generators
        sec_loop: do j=1, number_group_generators
           !> operators multiplication
           op_cart_new= matmul(pggen_cart(:, :, i), pggen_cart(:, :, j))
           tau_new= matmul(pggen_direct(:, :, i), tau_find(:, j))+ tau_find(:, i)

           !> transform to the home unit cell [0,1)
           tau_new(1)= tau_new(1)+1.0d0-int(tau_new(1)+1.0001d0)
           tau_new(2)= tau_new(2)+1.0d0-int(tau_new(2)+1.0001d0)
           tau_new(3)= tau_new(3)+1.0d0-int(tau_new(3)+1.0001d0)

           !> check whether the new opertors is already in the opertor list
           do k=1, number_group_operators
              if (sum(dabs(op_cart_new(:, :)- pgop_cart(:, :, k)))<eps6)then
                 cycle sec_loop
              endif
           enddo
           number_group_operators= number_group_operators+ 1
           pgop_cart(:, :, number_group_operators)= op_cart_new
           tau_direct(:, number_group_operators)= tau_new
           spatial_inversion(number_group_operators)= spatial_inversion(i)*spatial_inversion(j)
        enddo sec_loop
     enddo

     do i=1, number_group_operators 
        pgop_direct(:, :, i)= matmul(TT, matmul(pgop_cart(:, :, i), altv))
        tau_cart(:, i)= matmul(altv, tau_direct(:, i))
     enddo

     !> do it again in order to get all the operations
     number_group_generators= number_group_operators
     pggen_cart= pgop_cart
     pggen_direct= pgop_direct
     tau_find= tau_direct
     do i=1, number_group_generators
        sec_loop2: do j=1, number_group_generators
           !> operators multiplication
           op_cart_new= matmul(pggen_cart(:, :, i), pggen_cart(:, :, j))
           tau_new= matmul(pggen_direct(:, :, i), tau_find(:, j))+ tau_find(:, i)

           !> transform to the home unit cell [0,1)
           tau_new(1)= tau_new(1)+1.0d0-int(tau_new(1)+1.0001d0)
           tau_new(2)= tau_new(2)+1.0d0-int(tau_new(2)+1.0001d0)
           tau_new(3)= tau_new(3)+1.0d0-int(tau_new(3)+1.0001d0)

           !> check whether the new opertors is already in the opertor list
           do k=1, number_group_operators
              if (sum(dabs(op_cart_new(:, :)- pgop_cart(:, :, k)))<eps6)then
                 cycle sec_loop2
              endif
           enddo
           number_group_operators= number_group_operators+ 1
           pgop_cart(:, :, number_group_operators)= op_cart_new
           tau_direct(:, number_group_operators)= tau_new
           spatial_inversion(number_group_operators)= spatial_inversion(i)*spatial_inversion(j)
        enddo sec_loop2
     enddo

     !> till now, we found all the operators
     do i=1, number_group_operators 
        pgop_direct(:, :, i)= matmul(TT, matmul(pgop_cart(:, :, i), altv))
        tau_cart(:, i)= matmul(altv, tau_direct(:, i))
     enddo

     if (cpuid.eq.0) then
        write(stdout,'(a)')" "
        write(stdout, '(5x,a,i10)')">> Assuming there is time-reversal symmetry."
        write(stdout, '(5x,a,i10)')"We found total number of operators: ", number_group_operators 
        do i=1, number_group_operators 
           write(stdout, '(5x,a, i10)')"No. operators: ", i
           write(stdout, '(5x,a, f5.1)')"Inversion:  ", spatial_inversion(i)
           write(stdout, '(5x,a24, 13x, a24)')"Direct", "Cartesian"
           write(stdout,'(8x,3f12.6, 3x,3f12.6)')pgop_direct(1, :, i), pgop_cart(1, :, i)
           write(stdout,'(8x,3f12.6, 3x,3f12.6)')pgop_direct(2, :, i), pgop_cart(2, :, i)
           write(stdout,'(8x,3f12.6, 3x,3f12.6)')pgop_direct(3, :, i), pgop_cart(3, :, i)
           write(stdout,'(5x, "Tau",6f12.6)')tau_direct(:, i), tau_cart(:, i)
           write(stdout,'(a)')" "
        enddo
     endif

     return
  end subroutine FindSymmetryOperators

  subroutine Magnetic_adapted_Operators
     !> In this subroutine, we delete the operators which is not consist with the 
     !> magnetic moment and the magnetic field
     !> After the calling of this subroutine, the following information will be updated
     !> number_group_operators, pgop_cart, pgop_direct, spatial_inversion
     !> The principle is as follows:
     !> given a operator R and a pseudovector M at atom A, check whether RM is the same as the 
     !> pseudovector in another atom B which is the image of A related by the operator R
     use para

     implicit none

     integer :: iatomA, iatomB, iop, i, j, i1, i2, i3, ia, ifind, number_group_operators_left
     real(dp) :: op_rotate_cart(3, 3), op_rotate_direct(3, 3), op_inv, p1(3), MM_B(3), MM_OP_A(3), p2(3)
     real(dp) :: atom_position_A_cart(3), atom_position_B_cart(3)
     real(dp) :: atom_position_A_direct(3), atom_position_B_direct(3)

     !> magnetic_op_left= .False. if the operator is not consit with the Atom_magnetic_moment_field
     logical :: magnetic_op_left(48)

     !> magnetic moment after considering the magnetic field
     real(dp), allocatable :: Atom_magnetic_moment_field(:, :)

     !> space group operators left after consider the magnetic in cartesian and direct coordinate
     real(dp), allocatable :: pgop_cart_left(:, :, :)
     real(dp), allocatable :: pgop_direct_left(:, :, :)
     real(dp), allocatable :: tau_cart_left(:,:)
     real(dp), allocatable :: tau_direct_left(:,:)
     real(dp), allocatable :: spatial_inversion_left(:)
     integer, allocatable :: imap_sym_local(:, :)
     integer, allocatable :: imap_sym_local_left(:, :)
     allocate(pgop_cart_left(3, 3, 48), pgop_direct_left(3, 3, 48))
     allocate(tau_cart_left(3, 48), tau_direct_left(3, 48), spatial_inversion_left(48))
     allocate(Atom_magnetic_moment_field(3, Origin_cell%Num_atoms))
     allocate(imap_sym_local(Origin_cell%Num_atoms , number_group_operators))
     allocate(imap_sym_local_left(Origin_cell%Num_atoms, number_group_operators))
     imap_sym_local= -1

     !> magnetic moment after considering the magnetic field
     do ia=1, Origin_cell%Num_atoms
        Atom_magnetic_moment_field(1, ia)= Origin_cell%Atom_magnetic_moment(1, ia)+ Bx*Magneticfluxdensity_atomic
        Atom_magnetic_moment_field(2, ia)= Origin_cell%Atom_magnetic_moment(2, ia)+ By*Magneticfluxdensity_atomic
        Atom_magnetic_moment_field(3, ia)= Origin_cell%Atom_magnetic_moment(3, ia)+ Bz*Magneticfluxdensity_atomic
     enddo
     if (cpuid==0) write(stdout, '(a)', advance='no')'>>> Check the consistance of the symmetry'
     if (cpuid==0) write(stdout, '(a)')  'operators with the magnetic configuration.'

     magnetic_op_left(1:number_group_operators)= .True.
     do iop= 1, number_group_operators
        op_rotate_direct= pgop_direct(:, :, iop)  ! operator matrix in unit of lattice vector
        op_rotate_cart= pgop_cart(:, :, iop)  ! operator matrix in unit of lattice vector
        op_inv = spatial_inversion(iop) ! whether have inversion or not
        do iatomA= 1, Origin_cell%Num_atoms
           !> apply operator R to the atom A's position
           atom_position_A_direct= Origin_cell%Atom_position_direct(:, iatomA)
           atom_position_B_direct= matmul(op_rotate_direct, atom_position_A_direct)+ tau_direct(:, iop)
           !> now find the atomic index for atom B
           ifind= 0
           findB: do iatomB=1, Origin_cell%Num_atoms
              p1= Origin_cell%Atom_position_direct(:, iatomB)
              do i1=-3, 3
              do i2=-3, 3
              do i3=-3, 3
                 p2(1)= p1(1)+ i1
                 p2(2)= p1(2)+ i2
                 p2(3)= p1(3)+ i3
                 if (sum(abs(p2-atom_position_B_direct))<eps3) then
                    ifind= iatomB
                    exit findB
                 endif
              enddo
              enddo
              enddo
           enddo findB

           if (ifind==0) then
              write(*, *) 'ERROR: There are some bugs in the symmetry module!'
              write(*, *) "ERROR: I can't find the image of an atom after some symmetry operation."
              write(*, *) 'ERROR: Please turn off the symmetry functionlity!'
              stop
           endif

           !> after the symmtry operation, the iatomA is mapped to ifind
           imap_sym_local(iatomA, iop)= ifind

           !> magnetic moment of atom B which is the image after the symmetry operation R to A
           MM_B= Atom_magnetic_moment_field(:, ifind)

           !> Apply the operator onto the magnetic moment of atom A.
           MM_OP_A= matmul(op_rotate_cart*op_inv, Atom_magnetic_moment_field(:, iatomA))

           !> if MM_OP_A is not equal MM_B then the operator is not consist with the magnetic moment
           if (sum(abs(MM_OP_A- MM_B))> eps3) then
              magnetic_op_left(iop)= .False.
           endif

           if (cpuid==0) then
             !write(stdout, '(a, i8)')' No. of operators: ', iop 
             !write(stdout, '(a, 3f8.4)')' Operators in fractional unit with tau = ', tau_direct(:, iop)
             !write(stdout, '(a, f5.1)')' Inversion : ', spatial_inversion(iop)
             !write(stdout, '(2a24)')'Direct', 'Cartesian'
             !do i=1, 3
             !   write(stdout, '(10f8.4)') op_rotate_direct(i, :), op_rotate_cart(i, :)
             !enddo
             !write(stdout, '(a, i5, 5X, 3f8.4)')' Applying on atom A :', iatomA, atom_position_A_direct
             !write(stdout, '(a, i5, 5X, 3f8.4)')' Getting  a  atom B :', ifind, Origin_cell%Atom_position_direct(:, ifind)
             !write(stdout, '(a, 3f8.4)')'Magnetic moment on atom A :', Atom_magnetic_moment_field(:, iatomA)
             !write(stdout, '(a, 3f8.4)')'Magnetic moment on atom B :', Atom_magnetic_moment_field(:, ifind)
             !write(stdout, '(a, 3f8.4)')'Symmetry operation on magnetic moment of atom A ', MM_OP_A
             !write(stdout, *)' '
           endif

        enddo ! iatomA
     enddo  ! iop

     !> remove the non-compatiable operators
     number_group_operators_left=0
     do iop= 1, number_group_operators
        if (magnetic_op_left(iop))then
           number_group_operators_left= number_group_operators_left+ 1
           pgop_direct_left(:, :, number_group_operators_left)= pgop_direct(:, :, iop)
           pgop_cart_left(:, :, number_group_operators_left)= pgop_cart(:, :, iop)
           tau_direct_left(:, number_group_operators_left)= tau_direct(:, iop) 
           tau_cart_left(:, number_group_operators_left)= tau_cart(:, iop) 
           spatial_inversion_left(number_group_operators_left)= spatial_inversion(number_group_operators_left)
           imap_sym_local_left(:, number_group_operators_left)= imap_sym_local(:, iop) 
        endif
     enddo
     number_group_operators= number_group_operators_left
     if (allocated(imap_sym)) then
        deallocate(imap_sym)
        allocate(imap_sym(Origin_cell%Num_atoms, number_group_operators))
     endif

     !> only take the compatiable operators
     if (number_group_operators>0) then
        pgop_cart(:, :, 1:number_group_operators)= pgop_cart_left(:, :, 1:number_group_operators)
        pgop_direct(:, :, 1:number_group_operators)= pgop_direct_left(:, :, 1:number_group_operators)
        tau_cart(:, 1:number_group_operators)= tau_cart_left(:, 1:number_group_operators)
        tau_direct(:, 1:number_group_operators)= tau_direct_left(:, 1:number_group_operators)
        spatial_inversion(1:number_group_operators)= spatial_inversion_left(1:number_group_operators)
        imap_sym(:, 1:number_group_operators) = imap_sym_local_left(:, 1:number_group_operators)
     endif

     if (cpuid.eq.0) then
        write(stdout, '(a)')'>>> Space group operators'
        do iop= 1, number_group_operators
           write(stdout, '(a)')' '
           write(stdout, '(a, i3)')' No. of operators: OP=', iop 
           write(stdout, '(a, 3f8.4)')' Operators in fractional unit with tau = ', tau_direct(:, iop)
           write(stdout, '(a, f5.1)')' Inversion : ', spatial_inversion(iop)
           write(stdout, '(a)')' Rotation matrix:'
           write(stdout, '(a16, 5x, a24)')'Direct', 'Cartesian'
           do i=1, 3
              write(stdout, '(3f8.4, 4X, 3f8.4)') pgop_direct(i, :, iop), pgop_cart(i, :, iop)
           enddo
        enddo  ! iop
   
        write(stdout, '(a)')'>> A map between a atom in unit cell and the atom after space group operation'
        write(stdout, '(4x, a4, a8, 2x, 48("  OP=", i3))')' Name ', 'index ', (iop, iop=1, number_group_operators)
        do iatomA= 1, Origin_cell%Num_atoms
           write(stdout, '(3x, a4, 100i8)') trim(adjustl(Origin_cell%Atom_name(iatomA))), iatomA, imap_sym(iatomA, 1:number_group_operators)
        enddo ! iatomA
        write(stdout, '(a)')' '
     endif

     return
  end subroutine Magnetic_adapted_Operators

  subroutine symmetry
     !>> Set symmetry, not finished yet.
     !
     use para
     implicit none

     integer :: nwan, Nk_reduced, k_deg, Nkmesh
     integer :: ia, i, n, j, ik, ik1, ik2, ik3, i1, i2, i3, j1

     real(dp) :: k3_out_cart(3)
     real(dp) :: k3_in_cart(3), symm_matrix(3, 3), k3_out_direct(3)
     real(dp) :: altv(3, 3)

     !> get the atom afterr the mirror_x operation
     integer, allocatable :: iatom_mirror_x(:)
 
     !> set up symmetry operators
     !> here we assume that the wannier functions have the symmetry 
     !> as atomic orbitals
  
     
     real(dp), allocatable :: BasicOperations_space(:, :, :)
     complex(dp), allocatable :: BasicOperations_spin(:, :, :)
     real(dp), allocatable :: BasicOperations_inversion(:)

     !> generators_find, tau_find are defined in module.f90
     allocate(generators_find(48), tau_find(3, 48))

     allocate(BasicOperations_space(3, 3, 0:48))
     allocate(BasicOperations_spin(2, 2, 0:48))
     allocate(BasicOperations_inversion(0:48))
     allocate(pggen_cart(3, 3, 48), pggen_direct(3, 3, 48))
     allocate(pgop_cart(3, 3, 48), pgop_direct(3, 3, 48))
     allocate(tau_cart(3, 48), tau_direct(3, 48), spatial_inversion(48))
     allocate(imap_sym(Origin_cell%Num_atoms, 48))
     imap_sym= -1
     BasicOperations_space= 0d0
     BasicOperations_spin= 0d0
     BasicOperations_inversion= 0d0
     pgop_direct= 0d0; pgop_cart= 0d0
     pggen_cart= 0d0; pggen_direct= 0d0
     tau_cart= 0d0; tau_direct= 0d0 ; spatial_inversion= 0d0

     !> This is the unitary operation (unchange)
     do j=1, 3
        pgop_cart(j, j, 1)= 1d0
     enddo

     if (Symmetry_Import_calc) then

        !> define 18 basic operations
        call DefineBasicOperations(BasicOperations_space, BasicOperations_spin, BasicOperations_inversion)
   
        !> get generators for the space group
        call FindGenerators(BasicOperations_space, BasicOperations_inversion, number_group_generators, generators_find, tau_find)
   
        !> setup pggen_direct, pggen_cart, pgop_direct, pgop_cart, tau_cart, tau_direct
        !> get all the space group operators and the number of operators number_group_operators
        call FindSymmetryOperators(BasicOperations_space, BasicOperations_inversion)
   
        !> Reconsider the space group operators according to the atom's magnetic moment
        call Magnetic_adapted_Operators()
     endif

     if (.not.Symmetry_Import_calc) number_group_operators = 1
     
     Nk_reduced= 0
     if (Boltz_OHE_calc.or.DOS_calc) then
        !> setup KCube3D_symm
        KCube3D_symm%Nk_total= Nk1*Nk2*Nk3
        Nkmesh= (Nk1+1)*(Nk2+1)*(Nk3+1)
   
        !> The coordinate of the k point can be deduced from ik_array
        !> ikx= (ik-1)/(nk2*nk3)+1
        !> iky= ((ik-1-(ikx-1)*Nk2*Nk3)/nk3)+1
        !> ikz= (ik-(iky-1)*Nk3- (ikx-1)*Nk2*Nk3)
        !> k= K3D_start_cube+ K3D_vec1_cube*(ikx-1)/dble(nk1-1)  &
        !>  + K3D_vec2_cube*(iky-1)/dble(nk2-1)  &
        !>  + K3D_vec3_cube*(ikz-1)/dble(nk3-1)
        allocate(KCube3D_symm%ik_array_symm(Nkmesh))
        allocate(KCube3D_symm%ik_relate(Nkmesh))
        allocate(KCube3D_symm%weight_k(Nkmesh))
        KCube3D_symm%weight_k= 0d0
        do ik=1, Nk1*Nk2*Nk3
           KCube3D_symm%ik_relate(ik)=  ik
        enddo
   
        !> reduce k points in the first BZ with symmetries
        k_deg= 1
        do ik=1, Nk1*Nk2*Nk3
           !> ik= 1+(ik3-1)+(ik2-1)*Nk3+(ik1-1)*Nk3*Nk2
           ik1= (ik-1)/(Nk2*Nk3)+1
           ik2= ((ik-1-(ik1-1)*Nk2*Nk3)/Nk3)+1
           ik3= (ik-(ik2-1)*Nk3- (ik1-1)*Nk2*Nk3)
          !write(1001, '(4i5)')ik1, ik2, ik3, ik
           k3_in_cart = Origin_cell%Kua*(ik1-1d0)/dble(Nk1)  &
              + Origin_cell%Kub*(ik2-1d0)/dble(Nk2)  &
              + Origin_cell%Kuc*(ik3-1d0)/dble(Nk3)
           if (KCube3D_symm%ik_relate(ik).eq.ik) then
              Nk_reduced= Nk_reduced+ 1
              KCube3D_symm%ik_array_symm(Nk_reduced)= ik
              k_deg= 1
              do j=1, number_group_operators
                 symm_matrix= pgop_cart(:, :, j)
                 !> k3_out_cart= symm_matrix*k3_in_cart
                 call symm_operation(symm_matrix,k3_in_cart, k3_out_cart)
                 if (ik==1) goto 130
                 call cart_direct_rec(k3_out_cart, k3_out_direct)
   
                 i1= int8(k3_out_direct(1)*dble(Nk1)+ eps6)+1
                 i2= int8(k3_out_direct(2)*dble(Nk2)+ eps6)+1
                 i3= int8(k3_out_direct(3)*dble(Nk3)+ eps6)+1
                 j1= 1+(i3-1)+(i2-1)*Nk3+(i1-1)*Nk3*Nk2
   
                 if (j1/=ik .and. KCube3D_symm%ik_relate(j1)==j1) then
                    k_deg= k_deg+ 1
                 endif
                 KCube3D_symm%ik_relate(j1)= Nk_reduced
              enddo ! j=1, number_group_operators
              130 continue
              KCube3D_symm%weight_k(Nk_reduced) = dble(k_deg)/dble(KCube3D_symm%Nk_total)
           endif
        enddo  ! ik

        !> total number of reduced k points after symmetry operation
        KCube3D_symm%Nk_total_symm= Nk_reduced
        KCube3D_symm%ik_array_symm(Nk_reduced+1:KCube3D_symm%Nk_total)= 0
     endif

     if (cpuid.eq.0) then
        write(stdout, *)' '
        write(stdout, *)'>> We finished the symmetry reducing of the kpoints in the kcube'
        write(stdout, '(2x, a, i10)')"Number of kpoints mesh is ", Nk1*Nk2*Nk3
        write(stdout, '(2x, a, i10)')"Number of symmetry reduced kpoints is ", Nk_reduced
        write(stdout, *)' '
     endif

     if (SOC.eq.0) return
     nwan= Num_wann/2

     allocate(iatom_mirror_x(Origin_cell%Num_atoms))
     allocate(inversion(Num_wann, Num_wann))
     allocate(mirror_z(Num_wann, Num_wann))
     allocate(mirror_y(Num_wann, Num_wann))
     allocate(mirror_x(Num_wann, Num_wann))
     allocate(C2yT(Num_wann, Num_wann))
     inversion= 0d0
     mirror_z = 0d0
     mirror_y = 0d0
     mirror_x = 0d0
     C2yT = 0d0
     iatom_mirror_x= 0

     !> symmetry operators for VASP
     !> we always assume the orbital order is up up up up dn dn dn dn 
    !if (index( Package, 'VASP')/=0.or. index( Package, 'Wien2k')/=0 &
    !   .or. index( Package, 'Abinit')/=0.or. index( Package, 'openmx')/=0) then
        !> inversion symmetry
        !> s-> s; p-> -p; d-> d; f->-f
        n= 0
        do ia=1, Origin_cell%Num_atoms
           do i=1, Origin_cell%nprojs(ia)
              n= n+ 1
              select case (Origin_cell%proj_name(i, ia))
              case ('s', 'S')
                 inversion(n, n)= 1
                 inversion(n+ nwan, n+ nwan)= 1
              case ('px', 'Px', 'PX')
                 inversion(n, n)=-1
                 inversion(n+ nwan, n+ nwan)=-1
              case ('py', 'Py', 'PY')
                 inversion(n, n)=-1
                 inversion(n+ nwan, n+ nwan)=-1
              case ('pz', 'Pz', 'PZ')
                 inversion(n, n)=-1
                 inversion(n+ nwan, n+ nwan)=-1
              case ('dxy', 'Dxy', 'DXY')
                 inversion(n, n)= 1
                 inversion(n+ nwan, n+ nwan)= 1
              case ('dyz', 'Dyz', 'DYZ')
                 inversion(n, n)= 1
                 inversion(n+ nwan, n+ nwan)= 1
              case ('dxz', 'Dxz', 'DXZ', 'dzx', 'Dzx', 'DZX')
                 inversion(n, n)= 1
                 inversion(n+ nwan, n+ nwan)= 1
              case ('dx2-y2', 'Dx2-y2', 'DX2-Y2', 'dx2', 'DX2')
                 inversion(n, n)= 1
                 inversion(n+ nwan, n+ nwan)= 1
              case ('dz2', 'Dz2', 'DZ2')
                 inversion(n, n)= 1
                 inversion(n+ nwan, n+ nwan)= 1
              case ('FZ3', 'FZ2')
                 inversion(n, n)= -1
                 inversion(n+ nwan, n+ nwan)= -1
              case ('FXZ2')
                 inversion(n, n)= -1
                 inversion(n+ nwan, n+ nwan)= -1
              case ('FYZ2')
                 inversion(n, n)= -1
                 inversion(n+ nwan, n+ nwan)= -1
              case ('FXYZ')
                 inversion(n, n)= -1
                 inversion(n+ nwan, n+ nwan)= -1
              case ('FZ(X2-Y2)', 'FZX2')
                 inversion(n, n)= -1
                 inversion(n+ nwan, n+ nwan)= -1
              case ('FX(X2-3Y2)', 'FX3', "FX3Y2")
                 inversion(n, n)= -1
                 inversion(n+ nwan, n+ nwan)= -1
              case ('FY(3X2-Y2)', 'FY3X2', 'FY3')
                 inversion(n, n)= -1
                 inversion(n+ nwan, n+ nwan)= -1
              case default
                 write(*, *) "ERROR: only support s px py pz dxy dyz dxz dx2-y2 dz2, "
                 write(*, *) "fz3, fxz2, fyz2, fxyz, fzx2, fx3y2, fx3 orbitals"
                 stop
              end select
           enddo ! i
        enddo ! ia
 
        !> mirror_x symmetry
        !> s-> s; px->-px, py->py, pz->  pz
        !> dxy-> -dxy, dyz->  dyz, dxz-> -dxz, dx2-> dx2 dz2->dz2 
        !> up-> -i dn  dn-> -i up
        !> here we throw away the phase -i, it is just a constant, leading M*M= -1
      
        n= 0
        mirror_x= 0d0
        do ia=1, Origin_cell%Num_atoms
           do i=1, Origin_cell%nprojs(ia)
              n= n+ 1
              select case (Origin_cell%proj_name(i, ia))
              case ('s', 'S')
                 mirror_x(n, n+ nwan)= 1d0
                 mirror_x(n+ nwan, n)= 1d0
              case ('px', 'Px', 'PX')
                 mirror_x(n, n+ nwan)=-1d0
                 mirror_x(n+ nwan, n)=-1d0
              case ('py', 'Py', 'PY')
                 mirror_x(n, n+ nwan)= 1d0
                 mirror_x(n+ nwan, n)= 1d0
              case ('pz', 'Pz', 'PZ')
                 mirror_x(n, n+ nwan)= 1d0
                 mirror_x(n+ nwan, n)= 1d0
              case ('dxy', 'Dxy', 'DXY')
                 mirror_x(n, n+ nwan)=-1d0
                 mirror_x(n+ nwan, n)=-1d0
              case ('dyz', 'Dyz', 'DYZ')
                 mirror_x(n, n+ nwan)= 1d0
                 mirror_x(n+ nwan, n)= 1d0
              case ('dxz', 'Dxz', 'DXZ', 'dzx', 'Dzx', 'DZX')
                 mirror_x(n, n+ nwan)=-1d0
                 mirror_x(n+ nwan, n)=-1d0
              case ('dx2-', 'dx2-y2', 'Dx2-y2', 'DX2-Y2', 'dx2', 'DX2')
                 mirror_x(n, n+ nwan)= 1d0
                 mirror_x(n+ nwan, n)= 1d0
              case ('dz2', 'Dz2', 'DZ2')
                 mirror_x(n, n+ nwan)= 1d0
                 mirror_x(n+ nwan, n)= 1d0
               case ('FZ3', 'FZ2')
                 mirror_x(n, n+ nwan)= 1d0
                 mirror_x(n+ nwan, n)= 1d0
              case ('FXZ2')
                 mirror_x(n, n+ nwan)=-1d0
                 mirror_x(n+ nwan, n)=-1d0
              case ('FYZ2')
                 mirror_x(n, n+ nwan)= 1d0
                 mirror_x(n+ nwan, n)= 1d0
              case ('FXYZ')
                 mirror_x(n, n+ nwan)=-1d0
                 mirror_x(n+ nwan, n)=-1d0
              case ('FZ(X2-Y2)', 'FZX2')
                 mirror_x(n, n+ nwan)= 1d0
                 mirror_x(n+ nwan, n)= 1d0
              case ('FX(X2-3Y2)', 'FX3', "FX3Y2")
                 mirror_x(n, n+ nwan)=-1d0
                 mirror_x(n+ nwan, n)=-1d0
              case ('FY(3X2-Y2)', 'FY3X2', 'FY3')
                 mirror_x(n, n+ nwan)= 1d0
                 mirror_x(n+ nwan, n)= 1d0
 
              case default
                 write(*, *) "ERROR: only support s px py pz dxy dyz dxz dx2-y2 dz2, "
                 write(*, *) "fz3, fxz2, fyz2, fxyz, fzx2, fx3y2, fx3 orbitals"
                 stop
              end select
           enddo ! i
        enddo ! ia
 
        !> C2yT symmetry =  C_2y * i*sigma_y 
        !> s-> s; px->-px, py-> py, pz-> -pz
        !> dxy-> -dxy, dyz-> -dyz, dxz-> dxz, dx2-> dx2 dz2->dz2 
        !> up-> up  dn-> dn
        !> Here we keep the phase i of T, we also drop the Conjugation of T
      
        n= 0
        C2yT= 0d0
        do ia=1, Origin_cell%Num_atoms
           do i=1, Origin_cell%nprojs(ia)
              n= n+ 1
              select case (Origin_cell%proj_name(i, ia))
              case ('s', 'S')
                 C2yT(n, n)= 1d0
                 C2yT(n+ nwan, n+ nwan)= 1d0
              case ('px', 'Px', 'PX')
                 C2yT(n, n)=-1 
                 C2yT(n+ nwan, n+ nwan)=-1 
              case ('py', 'Py', 'PY')
                 C2yT(n, n)= 1 
                 C2yT(n+ nwan, n+ nwan)= 1 
              case ('pz', 'Pz', 'PZ')
                 C2yT(n, n)=-1 
                 C2yT(n+ nwan, n+ nwan)=-1 
              case ('dxy', 'Dxy', 'DXY')
                 C2yT(n, n)=-1 
                 C2yT(n+ nwan, n+ nwan)=-1 
              case ('dyz', 'Dyz', 'DYZ')
                 C2yT(n, n)=-1 
                 C2yT(n+ nwan, n+ nwan)=-1 
              case ('dxz', 'Dxz', 'DXZ', 'dzx', 'Dzx', 'DZX')
                 C2yT(n, n)= 1 
                 C2yT(n+ nwan, n+ nwan)= 1 
              case ('dx2-', 'dx2-y2', 'Dx2-y2', 'DX2-Y2', 'dx2', 'DX2')
                 C2yT(n, n)= 1 
                 C2yT(n+ nwan, n+ nwan)= 1 
              case ('dz2', 'Dz2', 'DZ2')
                 C2yT(n, n)= 1 
                 C2yT(n+ nwan, n+ nwan)= 1 
              case ('FZ3', 'FZ2')
                 C2yT(n, n)=-1 
                 C2yT(n+ nwan, n+ nwan)=-1 
              case ('FXZ2')
                 C2yT(n, n)=-1 
                 C2yT(n+ nwan, n+ nwan)=-1 
              case ('FYZ2')
                 C2yT(n, n)= 1 
                 C2yT(n+ nwan, n+ nwan)= 1 
              case ('FXYZ')
                 C2yT(n, n)= 1 
                 C2yT(n+ nwan, n+ nwan)= 1 
              case ('FZ(X2-Y2)', 'FZX2')
                 C2yT(n, n)=-1 
                 C2yT(n+ nwan, n+ nwan)=-1 
              case ('FX(X2-3Y2)', 'FX3', "FX3Y2")
                 C2yT(n, n)=-1 
                 C2yT(n+ nwan, n+ nwan)=-1 
              case ('FY(3X2-Y2)', 'FY3X2', 'FY3')
                 C2yT(n, n)= 1 
                 C2yT(n+ nwan, n+ nwan)= 1 
 
              case default
                 write(*, *) "ERROR: only support s px py pz dxy dyz dxz dx2-y2 dz2, "
                 write(*, *) "fz3, fxz2, fyz2, fxyz, fzx2, fx3y2, fx3 orbitals"
                 stop
              end select
           enddo ! i
        enddo ! ia
      
     

        !> mirror_y symmetry = i*sigma_y*R_y
        !> s-> s; px->px, py->-py, pz->  pz
        !> dxy-> -dxy, dyz-> -dyz, dxz->  dxz, dx2-> dx2 dz2->dz2 
        !> up-> -zi*dn  dn-> zi*up
        !> Drop off phase i
      
        n= 0
        mirror_y= 0d0
        do ia=1, Origin_cell%Num_atoms
           do i=1, Origin_cell%nprojs(ia)
              n= n+ 1
              select case (Origin_cell%proj_name(i, ia))
              case ('s', 'S')
                 mirror_y(n, n+ nwan)=-zi
                 mirror_y(n+ nwan, n)= zi
              case ('px', 'Px', 'PX')
                 mirror_y(n, n+ nwan)=-zi
                 mirror_y(n+ nwan, n)= zi
              case ('py', 'Py', 'PY')
                 mirror_y(n, n+ nwan)= zi
                 mirror_y(n+ nwan, n)=-zi
              case ('pz', 'Pz', 'PZ')
                 mirror_y(n, n+ nwan)=-zi
                 mirror_y(n+ nwan, n)= zi
              case ('dxy', 'Dxy', 'DXY')
                 mirror_y(n, n+ nwan)= zi
                 mirror_y(n+ nwan, n)=-zi
              case ('dyz', 'Dyz', 'DYZ')
                 mirror_y(n, n+ nwan)= zi
                 mirror_y(n+ nwan, n)=-zi
              case ('dxz', 'Dxz', 'DXZ', 'dzx', 'Dzx', 'DZX')
                 mirror_y(n, n+ nwan)=-zi
                 mirror_y(n+ nwan, n)= zi
              case ('dx2-', 'dx2-y2', 'Dx2-y2', 'DX2-Y2', 'dx2', 'DX2')
                 mirror_y(n, n+ nwan)=-zi
                 mirror_y(n+ nwan, n)= zi
              case ('dz2', 'Dz2', 'DZ2')
                 mirror_y(n, n+ nwan)=-zi
                 mirror_y(n+ nwan, n)= zi
              case ('FZ3', 'FZ2')
                 mirror_y(n, n+ nwan)=-zi
                 mirror_y(n+ nwan, n)= zi
              case ('FXZ2')
                 mirror_y(n, n+ nwan)=-zi
                 mirror_y(n+ nwan, n)= zi
              case ('FYZ2')
                 mirror_y(n, n+ nwan)= zi
                 mirror_y(n+ nwan, n)=-zi
              case ('FXYZ')
                 mirror_y(n, n+ nwan)= zi
                 mirror_y(n+ nwan, n)=-zi
              case ('FZ(X2-Y2)', 'FZX2')
                 mirror_y(n, n+ nwan)=-zi
                 mirror_y(n+ nwan, n)= zi
              case ('FX(X2-3Y2)', 'FX3', "FX3Y2")
                 mirror_y(n, n+ nwan)=-zi
                 mirror_y(n+ nwan, n)= zi
              case ('FY(3X2-Y2)', 'FY3X2', 'FY3')
                 mirror_y(n, n+ nwan)= zi
                 mirror_y(n+ nwan, n)=-zi
 
              case default
                 write(*, *) "ERROR: only support s px py pz dxy dyz dxz dx2-y2 dz2, "
                 write(*, *) "fz3, fxz2, fyz2, fxyz, fzx2, fx3y2, fx3 orbitals"
                 stop
              end select
           enddo ! i
        enddo ! ia
      
        do i=1, Num_wann
          !write(*, '(1000i2)')int(real(mirror_x(:, i)))
        enddo
      
      
        !> mirror_z symmetry  i*sigma_z*R_z, but here, we omit the i
        !> s-> s; px->px, py->py, pz-> -pz
        !> dxy-> dxy, dyz-> -dyz, dxz-> -dxz, dx2-> dx2 dz2->dz2 
        !> fz3-> -fz3, fxz2->  fxz2, fyz2-> fyz2, fxyz-> -fxyz, fzx2-> -fzx2, fx3y2->  fx3y2, fy3x2->  fy3x2
        !> up-> up  dn-> -dn  Drop off phase i
        n= 0
        do ia=1, Origin_cell%Num_atoms
           do i=1, Origin_cell%nprojs(ia)
              n= n+ 1
              select case (Origin_cell%proj_name(i, ia))
              case ('s', 'S')
                 mirror_z(n, n)= 1
                 mirror_z(n+ nwan, n+ nwan)=-1
              case ('px', 'Px', 'PX')
                 mirror_z(n, n)= 1
                 mirror_z(n+ nwan, n+ nwan)=-1
              case ('py', 'Py', 'PY')
                 mirror_z(n, n)= 1
                 mirror_z(n+ nwan, n+ nwan)=-1
              case ('pz', 'Pz', 'PZ')
                 mirror_z(n, n)=-1
                 mirror_z(n+ nwan, n+ nwan)= 1
              case ('dxy', 'Dxy', 'DXY')
                 mirror_z(n, n)= 1
                 mirror_z(n+ nwan, n+ nwan)=-1
              case ('dyz', 'Dyz', 'DYZ')
                 mirror_z(n, n)=-1
                 mirror_z(n+ nwan, n+ nwan)= 1
              case ('dxz', 'Dxz', 'DXZ', 'dzx', 'Dzx', 'DZX')
                 mirror_z(n, n)=-1
                 mirror_z(n+ nwan, n+ nwan)= 1
              case ('dx2-', 'dx2-y2', 'Dx2-y2', 'DX2-Y2', 'dx2', 'DX2')
                 mirror_z(n, n)= 1
                 mirror_z(n+ nwan, n+ nwan)=-1
              case ('dz2', 'Dz2', 'DZ2')
                 mirror_z(n, n)= 1
                 mirror_z(n+ nwan, n+ nwan)=-1
              case ('FZ3', 'FZ2')
                 mirror_z(n, n)= -1
                 mirror_z(n+ nwan, n+ nwan)=  1
              case ('FXZ2')
                 mirror_z(n, n)=  1
                 mirror_z(n+ nwan, n+ nwan)= -1
              case ('FYZ2')
                 mirror_z(n, n)=  1
                 mirror_z(n+ nwan, n+ nwan)= -1
              case ('FXYZ')
                 mirror_z(n, n)= -1
                 mirror_z(n+ nwan, n+ nwan)=  1
              case ('FZ(X2-Y2)', 'FZX2')
                 mirror_z(n, n)= -1
                 mirror_z(n+ nwan, n+ nwan)=  1
              case ('FX(X2-3Y2)', 'FX3', "FX3Y2")
                 mirror_z(n, n)=  1
                 mirror_z(n+ nwan, n+ nwan)= -1
              case ('FY(3X2-Y2)', 'FY3X2', 'FY3')
                 mirror_z(n, n)=  1
                 mirror_z(n+ nwan, n+ nwan)= -1
              case default
                 write(*, *) "ERROR: only support s px py pz dxy dyz dxz dx2-y2 dz2, "
                 write(*, *) "fz3, fxz2, fyz2, fxyz, fzx2, fx3y2, fx3 orbitals"
                 stop
              end select
           enddo ! i
        enddo ! ia
        do i=1, Num_wann
          !write(*, '(1000i2)')int(real(mirror_z(:, i)))
        enddo
 
    !endif


     !> set up symmetry operators
     allocate(mirror_x_op(3,3))
     allocate(mirror_y_op(3,3))
     allocate(mirror_z_op(3,3))

     !> for glide symmetry, (1:3, 1:3) shows the mirror operation, (1:3, 4) 
     !> gives the shift
     !> not finished yet
     allocate(glide_y_op(3,4))

     mirror_x_op= 0d0
     mirror_y_op= 0d0
     mirror_z_op= 0d0
     glide_y_op= 0d0

     mirror_x_op(1, 1)=-1d0
     mirror_x_op(2, 2)= 1d0
     mirror_x_op(3, 3)= 1d0

     mirror_y_op(1, 1)= 1d0
     mirror_y_op(2, 2)=-1d0
     mirror_y_op(3, 3)= 1d0

     mirror_z_op(1, 1)= 1d0
     mirror_z_op(2, 2)= 1d0
     mirror_z_op(3, 3)=-1d0

     glide_y_op(1, 1) = 1d0
     glide_y_op(2, 2) =-1d0
     glide_y_op(3, 3) = 1d0
     glide_y_op(:, 4) = (/0.5d0, 0.0d0, 0.5d0/)

     return
  end subroutine symmetry


  subroutine symm_operation(symm_matrix,k3_in, k3_out)
     !> apply symmetry operator on a k point 
     !> k3_out= symm_matrix*k3_in

     use para, only : dp, Origin_cell, stdout
     implicit none
     real(dp), intent(in ) :: k3_in(3)
     real(dp), intent(out) :: k3_out(3)
     real(dp), intent(in ) :: symm_matrix(3, 3)
     
     integer :: i, j
     real(dp) :: alpha,beta,gamma,total

     !> k3_in, k3_out, symm_matrix are in cartesian coordinate
     k3_out= 0.0d0
     do i= 1, 3
        do j= 1, 3
           k3_out(i)= k3_out(i)+ symm_matrix(i,j)*k3_in(j)
        end do
     end do

     !> get the direct coordinate (alpha, beta, gamma)
     total=Origin_cell%Kua(1)*Origin_cell%Kub(2)*Origin_cell%Kuc(3) &
          +Origin_cell%Kua(2)*Origin_cell%Kub(3)*Origin_cell%Kuc(1) &
          +Origin_cell%Kua(3)*Origin_cell%Kub(1)*Origin_cell%Kuc(2) &
          -Origin_cell%Kua(1)*Origin_cell%Kub(3)*Origin_cell%Kuc(2) &
          -Origin_cell%Kua(2)*Origin_cell%Kub(1)*Origin_cell%Kuc(3) &
          -Origin_cell%Kua(3)*Origin_cell%Kub(2)*Origin_cell%Kuc(1)
     alpha=(k3_out(1)*Origin_cell%Kub(2)*Origin_cell%Kuc(3) &
           +k3_out(2)*Origin_cell%Kub(3)*Origin_cell%Kuc(1) &
           +k3_out(3)*Origin_cell%Kub(1)*Origin_cell%Kuc(2) &
           -k3_out(1)*Origin_cell%Kub(3)*Origin_cell%Kuc(2) &
           -k3_out(2)*Origin_cell%Kub(1)*Origin_cell%Kuc(3) &
           -k3_out(3)*Origin_cell%Kub(2)*Origin_cell%Kuc(1))/total

     total=Origin_cell%Kub(1)*Origin_cell%Kua(2)*Origin_cell%Kuc(3) &
          +Origin_cell%Kub(2)*Origin_cell%Kua(3)*Origin_cell%Kuc(1) &
          +Origin_cell%Kub(3)*Origin_cell%Kua(1)*Origin_cell%Kuc(2) &
          -Origin_cell%Kub(1)*Origin_cell%Kua(3)*Origin_cell%Kuc(2) &
          -Origin_cell%Kub(2)*Origin_cell%Kua(1)*Origin_cell%Kuc(3) &
          -Origin_cell%Kub(3)*Origin_cell%Kua(2)*Origin_cell%Kuc(1)
     beta=(k3_out(1)*Origin_cell%Kua(2)*Origin_cell%Kuc(3) &
          +k3_out(2)*Origin_cell%Kua(3)*Origin_cell%Kuc(1) &
          +k3_out(3)*Origin_cell%Kua(1)*Origin_cell%Kuc(2) &
          -k3_out(1)*Origin_cell%Kua(3)*Origin_cell%Kuc(2) &
          -k3_out(2)*Origin_cell%Kua(1)*Origin_cell%Kuc(3) &
          -k3_out(3)*Origin_cell%Kua(2)*Origin_cell%Kuc(1))/total

     total=Origin_cell%Kuc(1)*Origin_cell%Kub(2)*Origin_cell%Kua(3) &
          +Origin_cell%Kuc(2)*Origin_cell%Kub(3)*Origin_cell%Kua(1) &
          +Origin_cell%Kuc(3)*Origin_cell%Kub(1)*Origin_cell%Kua(2) &
          -Origin_cell%Kuc(1)*Origin_cell%Kub(3)*Origin_cell%Kua(2) &
          -Origin_cell%Kuc(2)*Origin_cell%Kub(1)*Origin_cell%Kua(3) &
          -Origin_cell%Kuc(3)*Origin_cell%Kub(2)*Origin_cell%Kua(1)
     gamma=(k3_out(1)*Origin_cell%Kub(2)*Origin_cell%Kua(3) &
           +k3_out(2)*Origin_cell%Kub(3)*Origin_cell%Kua(1) &
           +k3_out(3)*Origin_cell%Kub(1)*Origin_cell%Kua(2) &
           -k3_out(1)*Origin_cell%Kub(3)*Origin_cell%Kua(2) &
           -k3_out(2)*Origin_cell%Kub(1)*Origin_cell%Kua(3) &
           -k3_out(3)*Origin_cell%Kub(2)*Origin_cell%Kua(1))/total

     !> shift it into primitive cell [0, 1)
     alpha=alpha+10.0d0
     alpha=alpha-int(alpha+1.0d-10)

     beta=beta+10.0d0
     beta=beta-int(beta+1.0d-10)

     gamma=gamma+10.0d0
     gamma=gamma-int(gamma+1.0d-10)

     !> finally in cartesian coordinate
     k3_out= alpha*Origin_cell%Kua+beta*Origin_cell%Kub+gamma*Origin_cell%Kuc

     return
  end subroutine symm_operation


