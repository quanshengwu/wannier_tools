!include "mkl.fi"
!~ include "mkl_dss.f90"


subroutine ham_3Dlandau_sparse1(nnz, Ndimq, Nq, k, acoo,jcoo,icoo)
!> We generate the TB hamiltonian under magnetic field with dense formated HmnR.
!> At the begining, we generate coo-formated hamlitonian. However, in the end,
!> the coo-formated hamlitonian will be transform into csr format
   use para
   implicit none

   integer, intent(inout) :: nnz
   integer, intent(in) :: Ndimq
   integer, intent(in) :: Nq
   real(dp), intent(in) :: k(3)
   !    complex(dp) :: ham_landau(Ndimq, Ndimq)
   complex(dp), intent(inout) :: acoo(nnz)
   integer, intent(inout) :: jcoo(nnz)
   integer, intent(inout) :: icoo(nnz)

   !> inta-hopping for the supercell
   complex(dp), allocatable :: H00(:, :)

   !> inter-hopping for the supercell
   complex(dp), allocatable :: H01(:, :)
   complex(dp), allocatable :: H02(:, :)
   ! loop index
   integer :: i1, i2,i, i_t, j_t
   integer :: coo1,coo2,ncoo
   ! loop index
   integer :: iR

   ! index used to sign irvec
   real(dp) :: ia,ib,ic
   integer :: ia1, ia2

   integer :: istart1, istart2
   integer :: iend1, iend2

   integer :: inew_ic

   !> nwann= Num_wann/2
   integer :: nwann

   integer, allocatable :: orbital_start(:)

   !> new index used to sign irvec
   real(dp) :: new_ia,new_ib,new_ic

   !> wave vector k times lattice vector R
   real(Dp) :: kdotr
   real(dp) :: phase
   complex(dp) :: ratio
   complex(dp) :: fac,tmp

   real(dp) :: Rp1(3)
   real(dp) :: Rp2(3)
   real(dp) :: R1(3)
   real(dp) :: R2(3)
   real(dp) :: Ri(3)
   real(dp) :: Rj(3)
   real(dp) :: tau1(3)
   real(dp) :: tau2(3)

   !> calculate phase induced by magnetic field
   real(dp),external :: phase1,phase2

   !> for leaking purpose, only will be allocated when nnz is not big enough
   !complex(dp), allocatable :: acoo_extend(:)
   !integer, allocatable :: icoo_extend(:)
   !integer, allocatable :: jcoo_extend(:)

   !> a work array for coocsr
   integer(li), allocatable, target :: iwk (:)

   allocate( H00( Num_wann, Num_wann))
   allocate( H01( Num_wann, Num_wann))
   allocate( H02( Num_wann, Num_wann))
   H00= zzero
   H01= zzero
   H02= zzero

   nwann= Num_wann/2
   allocate( orbital_start(Origin_cell%Num_atoms+ 1))
   orbital_start= 0
   orbital_start(1)= 1
   do ia1=1, Origin_cell%Num_atoms
      orbital_start(ia1+1)= orbital_start(ia1)+ Origin_cell%nprojs(ia1)
   enddo

   !    write(*,*) Origin_cell%Num_atoms
   !> calculate intra-hopping
   acoo=zzero
   ncoo=0

   ! i1 column index, sweep magnetic supercells
   do i1=1, Nq
      ! i2 row index, sweep magnetic supercells
      do i2=1, Nq
         H00=zzero
         if (abs(i2-i1)> ijmax) cycle

         !> sum over R points to get H(k1, k2)
         do iR=1, Nrpts ! iR'th in lattice vector
            ia=dble(irvec(1,iR))
            ib=dble(irvec(2,iR))
            ic=dble(irvec(3,iR))

            !> new lattice
            call latticetransform(ia, ib, ic, new_ia, new_ib, new_ic)

            inew_ic= int(new_ic)
            if (inew_ic /= (i2-i1)) cycle

            !> exp(i k.R)
            kdotr= k(1)*new_ia+ k(2)*new_ib
            ratio= cos(2d0*pi*kdotr)+ zi*sin(2d0*pi*kdotr)

            Rp1= (i1-1)*Ruc_new
            Rp2= new_ia*Rua_new+ new_ib*Rub_new+ (i2-1)*Ruc_new

            do ia1=1, Origin_cell%Num_atoms
               do ia2=1, Origin_cell%Num_atoms

                  R1= Origin_cell%Atom_position_cart(:, ia1)
                  R2= Origin_cell%Atom_position_cart(:, ia2)

                  call rotate(R1, tau1)
                  call rotate(R2, tau2)


                  Ri= Rp1+ tau1
                  Rj= Rp2+ tau2

                  phase=phase2(Ri,Rj)
                  fac= cos(phase)+ zi*sin(phase)

                  istart1= orbital_start(ia1)
                  iend1= orbital_start(ia1+1)- 1
                  istart2= orbital_start(ia2)
                  iend2= orbital_start(ia2+1)- 1

                  H00 ( istart1:iend1, istart2:iend2) &
                     = H00 ( istart1:iend1, istart2:iend2) &
                     + HmnR( istart1:iend1, istart2:iend2, iR)*ratio/ndegen(iR)* fac

                  !> there is soc term in the hr file
                  if (soc>0) then
                     istart1=  orbital_start(ia1) + Nwann
                     iend1=  orbital_start(ia1+1)- 1 + Nwann
                     istart2=  orbital_start(ia2)
                     iend2=  orbital_start(ia2+1)- 1

                     H00 ( istart1:iend1, istart2:iend2) &
                        = H00 ( istart1:iend1, istart2:iend2) &
                        + HmnR( istart1:iend1, istart2:iend2, iR)*ratio/ndegen(iR)* fac

                     istart1=  orbital_start(ia1)
                     iend1=  orbital_start(ia1+1)- 1
                     istart2=  orbital_start(ia2) + Nwann
                     iend2=  orbital_start(ia2+1)- 1 + Nwann

                     H00 ( istart1:iend1, istart2:iend2) &
                        = H00 ( istart1:iend1, istart2:iend2) &
                        + HmnR( istart1:iend1, istart2:iend2, iR)*ratio/ndegen(iR)* fac

                     istart1=  orbital_start(ia1) + Nwann
                     iend1=  orbital_start(ia1+1)- 1 + Nwann
                     istart2=  orbital_start(ia2) + Nwann
                     iend2=  orbital_start(ia2+1)- 1 + Nwann

                     H00 ( istart1:iend1, istart2:iend2) &
                        = H00 ( istart1:iend1, istart2:iend2) &
                        + HmnR( istart1:iend1, istart2:iend2, iR)*ratio/ndegen(iR)* fac

                  endif ! soc


               enddo ! ia2
            enddo ! ia1
         enddo ! iR

         !> get the sparse matrix elements
         do coo2=1,Num_wann
            do coo1=1,Num_wann
               tmp=H00(coo1,coo2)
               if(abs(tmp)/eV2Hartree >eps6) then
                  j_t=coo1+(i1-1)*Num_wann
                  i_t=coo2+(i2-1)*Num_wann
                  !if (i_t>=j_t) then
                  ncoo=ncoo+1
                  jcoo(ncoo)=j_t
                  icoo(ncoo)=i_t
                  acoo(ncoo)=tmp !0.01
                  !endif
               end if
            end do
         end do
         H00=0d0
      enddo ! i2

   enddo ! i1


   if (cpuid.eq.0) write(stdout,*) 'nnz, nnzmax after H00 : ', ncoo, nnz
   if (nnz< ncoo) then
      write(*, *)'Please increase nnz in the sparse.f90'
      write(*,*) 'nnz, nnz after H00 : ', ncoo, nnz
      stop
   endif

   !> if we want to calculate the Chern number, we use open boundary such
   !> that we can get the Chern number from the edge states.
   if(landau_chern_calc) goto 328

   !>> calculate inter-hopping between layers
   ! i1 column index
   do i1=1, Nq
      ! i2 row index
      do i2=1, Nq
         H01=zzero
         if (abs(i2+Nq -i1)> ijmax) cycle
         !> sum over R points to get H(k1, k2)
         do iR=1, Nrpts
            ia=dble(irvec(1,iR))
            ib=dble(irvec(2,iR))
            ic=dble(irvec(3,iR))

            !> new lattice
            call latticetransform(ia, ib, ic, new_ia, new_ib, new_ic)
            inew_ic= int(new_ic)
            if (inew_ic /= (i2+ Nq -i1)) cycle

            !> exp(i k.R)
            kdotr= k(1)*new_ia+ k(2)*new_ib
            ratio= cos(2d0*pi*kdotr)+ zi*sin(2d0*pi*kdotr)

            Rp1= (i1-1)*Ruc_new
            Rp2= new_ia*Rua_new+ new_ib*Rub_new+ (i2-1+ Nq)*Ruc_new

            do ia1=1, Origin_cell%Num_atoms
               do ia2=1, Origin_cell%Num_atoms

                  R1= Origin_cell%Atom_position_cart(:, ia1)
                  R2= Origin_cell%Atom_position_cart(:, ia2)
                  call rotate(R1, tau1)
                  call rotate(R2, tau2)

                  Ri= Rp1+ tau1
                  Rj= Rp2+ tau2
                  phase=phase2(Ri,Rj)

                  fac= cos(phase)+ zi*sin(phase)

                  istart1=orbital_start(ia1)
                  iend1= orbital_start(ia1+1)- 1
                  istart2= orbital_start(ia2)
                  iend2= orbital_start(ia2+1)- 1

                  H01 ( istart1:iend1, istart2:iend2) &
                     = H01 ( istart1:iend1, istart2:iend2) &
                     + HmnR( istart1:iend1, istart2:iend2, iR)*ratio/ndegen(iR)* fac

                  !> there is soc term in the hr file
                  if (soc>0) then
                     istart1=  orbital_start(ia1) + Nwann
                     iend1=  orbital_start(ia1+1)- 1 + Nwann
                     istart2=  orbital_start(ia2)
                     iend2=  orbital_start(ia2+1)- 1

                     H01 ( istart1:iend1, istart2:iend2) &
                        = H01 ( istart1:iend1, istart2:iend2) &
                        + HmnR( istart1:iend1, istart2:iend2, iR)*ratio/ndegen(iR)* fac

                     istart1=  orbital_start(ia1)
                     iend1=  orbital_start(ia1+1)- 1
                     istart2=  orbital_start(ia2) + Nwann
                     iend2=  orbital_start(ia2+1)- 1 + Nwann

                     H01 ( istart1:iend1, istart2:iend2) &
                        = H01 ( istart1:iend1, istart2:iend2) &
                        + HmnR( istart1:iend1, istart2:iend2, iR)*ratio/ndegen(iR)* fac

                     istart1=  orbital_start(ia1) + Nwann
                     iend1=  orbital_start(ia1+1)- 1 + Nwann
                     istart2=  orbital_start(ia2) + Nwann
                     iend2=  orbital_start(ia2+1)- 1 + Nwann

                     H01 ( istart1:iend1, istart2:iend2) &
                        = H01 ( istart1:iend1, istart2:iend2) &
                        + HmnR( istart1:iend1, istart2:iend2, iR)*ratio/ndegen(iR)* fac
                  endif ! soc
               enddo ! ia2
            enddo ! ia1
         enddo ! iR

         do coo2=1,Num_wann!*(i1-1),Num_wann*i1
            do coo1=1,Num_wann!*(i2-1),Num_wann*i2
               tmp=H01(coo1,coo2)
               if(abs(tmp)/eV2Hartree >eps6) then

                  !> only store the upper triangle matrix for mkl_dss purpose
                  j_t=coo1+(i1-1)*Num_wann
                  i_t=coo2+(i2-1)*Num_wann
                  !if (i_t>=j_t) then
                  ncoo=ncoo+1
                  jcoo(ncoo)=j_t
                  icoo(ncoo)=i_t
                  acoo(ncoo)=tmp*exp(zi*k(3))
                  !endif

                  !> add the transpose-conjugate part
                  i_t=coo1+(i1-1)*Num_wann
                  j_t=coo2+(i2-1)*Num_wann
                  !if (i_t>=j_t) then
                  ncoo=ncoo+1
                  jcoo(ncoo)=j_t
                  icoo(ncoo)=i_t
                  acoo(ncoo)=conjg(tmp)*exp(-zi*k(3)) !0.01
                  !endif
               end if
            end do
         end do
      enddo ! i2

   enddo ! i1

328 continue
   if (cpuid.eq.0) write(stdout,*) 'nnz, nnzmax after H01:', ncoo, nnz
   if (nnz< ncoo) then
      write(*, *)'Please increase nnzmax in the sparse.f90'
      write(*,*) 'nnz, nnz after H01 : ', ncoo, nnz
      stop
   endif

   nnz= ncoo

   return
end subroutine ham_3Dlandau_sparse1

subroutine ham_3Dlandau_sparseHR(nnz, Ndimq, Nq, k, acoo,jcoo,icoo)
!> We generate the TB hamiltonian under magnetic field with sparse formated HmnR.
!> At the begining, we generate coo-formated hamlitonian. However, in the end,
!> the coo-formated hamlitonian will be transform into csr format
   use para
   implicit none

   !> input: nnz is the maximum number of non-zeros entries
   !> output: nnz is the number of non-zeros entries of acoo
   integer, intent(inout) :: nnz
   integer, intent(in) :: Ndimq
   integer, intent(in) :: Nq
   real(dp), intent(in) :: k(3)

   !> output hamiltonian stored as COO sparse matrix format
   complex(dp), intent(inout) :: acoo(nnz)
   integer, intent(inout) :: jcoo(nnz)
   integer, intent(inout) :: icoo(nnz)
   integer,allocatable :: rxyz(:,:),nonzero_counter(:,:,:)
   

   ! loop index
   integer :: i1, i2,i, i_t, j_t,i3
   integer :: coo1,coo2,ncoo

   ! loop index
   integer :: iR,ims

   ! index used to sign irvec
   real(dp) :: ia,ib,ic
   integer :: ia1, ia2

   integer :: inew_ic,inew_ia,inew_ib

   !> new index used to sign irvec
   real(dp) :: new_ia,new_ib,new_ic

   !> wave vector k times lattice vector R
   real(Dp) :: kdotr
   real(dp) :: phase
   complex(dp) :: ratio
   complex(dp) :: fac,tmp

   real(dp) :: Rp1(3), Rp2(3), R1(3), R2(3)
   real(dp) :: Ri(3), Rj(3), tau1(3), tau2(3)

   !> calculate phase induced by magnetic field
   real(dp),external :: phase1,phase2

   !> a work array for coocsr
   integer(li), allocatable, target :: iwk (:)

   if(export_maghr) then
      allocate(rxyz(nnz,3))
      allocate(nonzero_counter(-ijmax:ijmax,-ijmax:ijmax,-ijmax:ijmax))
      nonzero_counter=0
   endif
   

   !> calculate intra-hopping
   acoo=zzero
   ncoo=0
   tmp=0d0
   ! i1 column index, sweep magnetic supercells
   do i1=1, Nq
      ! i2 row index, sweep magnetic supercells
      do i2=1, Nq
         if (abs(i2-i1)> ijmax) cycle

         !> sum over R points to get H(k1, k2)
         do ims=1,splen
            ia=dble(hirv(1,ims))
            ib=dble(hirv(2,ims))
            ic=dble(hirv(3,ims))

            !> new lattice
            call latticetransform(ia, ib, ic, new_ia, new_ib, new_ic)

            inew_ic= int(new_ic)
            inew_ib=int(new_ib);inew_ia=int(new_ia)
            if (inew_ic /= (i2-i1)) cycle

            !> exp(i k.R)
            kdotr= k(1)*new_ia+ k(2)*new_ib
            ratio= cos(2d0*pi*kdotr)+ zi*sin(2d0*pi*kdotr)

            Rp1= (i1-1)*Ruc_new
            Rp2= new_ia*Rua_new+ new_ib*Rub_new+ (i2-1)*Ruc_new

            ia1=hicoo(ims)
            ia2=hjcoo(ims)

            R1= Origin_cell%wannier_centers_cart(:, ia1)
            R2= Origin_cell%wannier_centers_cart(:, ia2)

            call rotate(R1, tau1)
            call rotate(R2, tau2)


            Ri= Rp1+ tau1
            Rj= Rp2+ tau2

            phase=phase2(Ri,Rj)
            fac= cos(phase)+ zi*sin(phase)

            tmp=hacoo(ims)*ratio* fac
            if(abs(tmp)/eV2Hartree > 1e-6) then
               ncoo=ncoo+1
               icoo(ncoo)=hicoo(ims)+(i1-1)*Num_wann
               jcoo(ncoo)=hjcoo(ims)+(i2-1)*Num_wann
               acoo(ncoo)=acoo(ncoo)+tmp
              !if (export_maghr) then
              !   rxyz(ncoo,:)=[inew_ia,inew_ib,inew_ic]
              !   nonzero_counter(inew_ia,inew_ib,inew_ic)=1
              !endif
            endif
         enddo ! iR

      enddo ! i2
   enddo ! i1

   !> if we want to calculate the Chern number, we use open boundary such
   !> that we can get the Chern number from the edge states.
   if(landau_chern_calc) goto 239

   do i1=1, Nq
      ! i2 row index, sweep magnetic supercells
      do i2=1, Nq
         if (abs(i2-i1+Nq)> ijmax) cycle

         !> sum over R points to get H(k1, k2)
         do ims=1,splen
            ia=dble(hirv(1,ims))
            ib=dble(hirv(2,ims))
            ic=dble(hirv(3,ims))

            !> new lattice
            call latticetransform(ia, ib, ic, new_ia, new_ib, new_ic)

            inew_ic= int(new_ic)
            inew_ib=int(new_ib);inew_ia=int(new_ia)
            if (inew_ic /= (i2-i1+Nq)) cycle

            !> exp(i k.R)
            kdotr= k(1)*new_ia+ k(2)*new_ib+ k(3)
            ratio= cos(2d0*pi*kdotr)+ zi*sin(2d0*pi*kdotr)

            Rp1= (i1-1)*Ruc_new
            Rp2= new_ia*Rua_new+ new_ib*Rub_new+ (i2-1+Nq)*Ruc_new

            ia1=hicoo(ims)
            ia2=hjcoo(ims)

            R1= Origin_cell%wannier_centers_cart(:, ia1)
            R2= Origin_cell%wannier_centers_cart(:, ia2)


            call rotate(R1, tau1)
            call rotate(R2, tau2)


            Ri= Rp1+ tau1
            Rj= Rp2+ tau2

            phase=phase2(Ri,Rj)
            fac= cos(phase)+ zi*sin(phase)

            tmp=hacoo(ims)*ratio* fac
            if(abs(tmp)/eV2Hartree  > 1e-6) then
               ncoo=ncoo+1
               icoo(ncoo)=hicoo(ims)+(i1-1)*Num_wann
               jcoo(ncoo)=hjcoo(ims)+(i2-1)*Num_wann
               acoo(ncoo)=acoo(ncoo)+tmp
               if(export_maghr) then
               rxyz(ncoo,:)=[inew_ia,inew_ib,inew_ic]
               nonzero_counter(inew_ia,inew_ib,inew_ic)=1
               endif
            endif
            tmp=conjg(tmp)
            if(abs(tmp)/eV2Hartree  > 1e-6) then
               ncoo=ncoo+1
               icoo(ncoo)=hjcoo(ims)+(i2-1)*Num_wann
               jcoo(ncoo)=hicoo(ims)+(i1-1)*Num_wann
               acoo(ncoo)=acoo(ncoo)+tmp
              !if (export_maghr) then
              !   rxyz(ncoo,:)=-[inew_ia,inew_ib,inew_ic]
              !   nonzero_counter(-inew_ia,-inew_ib,-inew_ic)=1
              !endif
            endif

         enddo ! iR

      enddo ! i2
   enddo ! i1

239 continue
   if (cpuid.eq.0) write(stdout,*) 'nnz, nnzmax after H01:', ncoo, nnz
   if (nnz< ncoo) then
      write(*, *)'>>> Error : Please increase nnz in the sparse.f90'
      write(*, *) 'nnz, nnzmax after H01 : ', ncoo, nnz
      stop
   endif

   nnz= ncoo

   !> usually, we don't export magnetic hr file 
   if(export_maghr) then
      outfileindex= outfileindex+ 1
      open(unit=outfileindex,file='hr_mag.dat')
      write(outfileindex,*) '!magnetic supercell hr automated'
      write(outfileindex,*) ncoo
      write(outfileindex,*) ndimq
      write(outfileindex,*) sum(nonzero_counter)
      write(outfileindex,*) ('1   ',i=1,sum(nonzero_counter))
      do i1=-ijmax,ijmax
         do i2=-ijmax,ijmax
            do i3=-ijmax,ijmax
               if(nonzero_counter(i1,i2,i3)==0) cycle
               do i=1,ncoo
                  if(sum(abs(rxyz(i,:)-[i1,i2,i3]))/=0) cycle
                  write(outfileindex,*) rxyz(i,:),icoo(i),jcoo(i),real(acoo(i)),imag(acoo(i))
               end do
            end do
         end do
      end do
   endif
   if(export_maghr) deallocate(rxyz)
   return
end subroutine ham_3Dlandau_sparseHR

!use this to make eigs from coo
!> use solver from ARPACK


subroutine sparse_landau_level_B
   use para
   use sparse
   implicit none

   !> magnetic supercell size, perpendicular to the magnetic field
   !> Ndimq= Nq* Num_wann
   integer :: Nq, Ndimq, Nmag, Nmag1

   integer :: ia1, ia2, ib, i, j, kk1, kk2, ie, iq, ierr, ig

   real(dp) :: B0, B0Tesla, B0Tesla_quantumflux_magsupcell, theta, dis, dis1

   ! wave vector
   real(dp) :: k3(3)

   !> dim= Ndimq, knv3
   real(dp), allocatable :: mag(:), mag_Tesla(:), W(:), eigv(:, :), eigv_mpi(:, :)

   real(dp) :: time_start, time_end, time_start0

   !> dim= Ndimq*Ndimq
   !> maximum number of non-zeros matrix elements
   integer :: nnzmax, nnz

   complex(dp), allocatable :: acoo(:)
   integer, allocatable :: icoo(:), jcoo(:)

   !number of ARPACK eigenvalues to be obtained
   integer :: neval

   ! number of Arnoldi vectors
   integer :: nvecs

   !shift-invert sigma
   complex(dp) :: sigma

   !> time measurement
   real(dp) :: time1, time2, time3, times, timee

   !> eigenvector of the sparse matrix acoo. Dim=(Ndimq, neval)
   complex(dp), allocatable :: psi(:)

   !> the dimension of zeigv is (Ndimq, nvecs)
   !> however, only (1:Ndimq, 1:neval) are useful
   complex(dp), allocatable :: zeigv(:, :)

   !> print the weight for the Selected_WannierOrbitals
   real(dp), allocatable :: dos_selected(:, :, :), dos_selected_mpi(:, :, :)


   Nq= Magq
   Ndimq= Num_wann* Nq
   if (NumSelectedEigenVals==0) NumSelectedEigenVals=OmegaNum

   neval=NumSelectedEigenVals
   if (neval>=Ndimq) neval= Ndimq- 2

   !> ncv
   nvecs=int(3*neval)

   if (nvecs<50) nvecs= 50
   if (nvecs>Ndimq) nvecs= Ndimq


   sigma=(1d0,0d0)*iso_energy

   Nmag= Magp-1
   Nmag1=Nmag/1

   if (nnzmax_input<0)then
      nnzmax= Num_wann*(2*ijmax+1)*Ndimq+Ndimq
      if(Is_Sparse_Hr) nnzmax=splen*nq+Ndimq
   else
      nnzmax= nnzmax_input
   endif
 
   allocate( acoo(nnzmax), stat=ierr)
   call printallocationinfo('acoo', ierr)
   allocate( jcoo(nnzmax), stat=ierr)
   call printallocationinfo('jcoo', ierr)
   allocate( icoo(nnzmax), stat=ierr)
   call printallocationinfo('icoo', ierr)
   allocate( W( neval))
   allocate( eigv( neval, Nmag1+1), stat=ierr)
   call printallocationinfo('eigv', ierr)
   allocate( eigv_mpi( neval, Nmag1+1), stat=ierr)
   call printallocationinfo('eigv_mpi', ierr)
   allocate( mag(Nmag+1), mag_Tesla(Nmag+ 1))
   allocate( psi(ndimq))

   allocate( zeigv(ndimq,nvecs), stat=ierr)
   call printallocationinfo('zeigv', ierr)
   allocate( dos_selected     (neval,   Nmag+ 1, NumberofSelectedOrbitals_groups), stat=ierr)
   call printallocationinfo('dos_selected', ierr)
   allocate( dos_selected_mpi (neval,   Nmag+ 1, NumberofSelectedOrbitals_groups), stat=ierr)
   call printallocationinfo('dos_selected_mpi', ierr)
   zeigv= 0d0
   dos_selected= 0d0
   dos_selected_mpi= 0d0

   mag= 0d0; mag_Tesla= 0d0
   eigv_mpi= 0d0
   eigv    = 0d0
   acoo=0d0
   jcoo=0
   icoo=0


   if (cpuid.eq.0) write(stdout,*) 'sigma=',sigma,Ndimq,NumSelectedEigenVals
   !> deal with the magnetic field
   !> first transform the Bx By into B*Cos\theta, B*Sin\theta
   if (abs(By)<1e-8) then
      if (Bx<0) then
         theta= pi
      else
         theta= 0d0
      endif
   elseif (By>0) then
      theta = atan(Bx/By)
   else
      theta = atan(Bx/By)+ pi
   endif

   !> the magnetic flux in the magnetic supercell is 2*pi*Magp
   B0= 2d0*pi/dble(Nq)
   B0= abs(B0)
   Bx= B0* Cos(theta)
   By= B0* Sin(theta)

   !> transform it into Tesla
   !> B=2*pi*\phi_0/S0, where \phi_0 is the quantum flux, S0 is the projected area
   !> \phi_0 = h/2e  h=6.62607004*1E-34, e= 1.6*1E-19
   !> B0Tesla_quantumflux_magsupcell is the magnetic field that makes the flux through the magnetic
   !> supercell to be one quantum flux
   B0Tesla_quantumflux_magsupcell= 6.62607004*1E-34/1.6021766208*1E19/MagneticSuperProjectedArea*1E20
   B0Tesla= 6.62607004*1E-34/1.6021766208*1E19/MagneticSuperProjectedArea*1E20*Magp

   if (cpuid==0) then
      write(stdout, '(a, 2E18.8)')' Magnetic field B in Tesla that makes the flux through '
      write(stdout, '(a, 2E18.8)')' the magnetic supercell to be one quantum flux : ', B0Tesla_quantumflux_magsupcell
      write(stdout, '(a, 2E18.8)')' Magnetic field B in Tesla that fits to Magp : ', B0Tesla
      write(stdout, '(a, 2f18.8)')' Magnetic field Bx, By= ', Bx, By
   endif


   k3=Single_KPOINT_3D_DIRECT

   !> calculate the landau levels along special k line

   !> The flux in the unit cell changes from 0 to 2*pi
   do ib=1, Nmag+ 1
      mag(ib)= B0* (ib-1)
      mag_Tesla(ib)= B0Tesla_quantumflux_magsupcell* (ib-1)
   enddo

   time_start= 0d0
   time_start0= 0d0
   call now(time_start0)
   time_start= time_start0
   time_end  = time_start0
   do ib=1+ cpuid, Nmag1+1, num_cpu
      call now(times)
      !Bx= mag(ib)* Cos(theta)
      !By= mag(ib)* Sin(theta)
      Bx= -mag(ib)
      By= 0d0
      nnz= nnzmax
      if (cpuid.eq.0) &
         write(stdout, '(2a, i9, "  /", i10, a, f10.1, "s", a, f10.1, "s")') &
         'In sparse_landau_level_B ', ' ib/Nmag ', ib, Nmag+1, &
         ' time elapsed: ', time_end-time_start0, &
         ' time left: ', ((Nmag-ib)/num_cpu)*(time_end-time_start)

      call now(time_start)

      call now(time1)
      if(Is_Sparse_Hr) then
         call ham_3Dlandau_sparseHR(nnz, Ndimq, Nq, k3, acoo,jcoo,icoo)
      else
         call ham_3Dlandau_sparse1(nnz, Ndimq, Nq, k3, acoo,jcoo,icoo)
      endif
      acoo=acoo/eV2Hartree
      call now(time2)

      !> diagonalization by call zheev in lapack
      W= 0d0
      call arpack_sparse_coo_eigs(Ndimq,nnzmax,nnz,acoo,jcoo,icoo,neval,nvecs,W,sigma, zeigv, LandauLevel_wavefunction_calc)
      call now(time3)
      eigv(:, ib)= W
      if (cpuid==0)write(stdout, '(a, f20.2, a)')'  >> Time cost for constructing H: ', time2-time1, ' s'
      if (cpuid==0)write(stdout, '(a, f20.2, a)')'  >> Time cost for diagonalize H: ', time3-time2, ' s'

      if (LandauLevel_wavefunction_calc) then
         do ie= 1, neval
            psi(:)= zeigv(:, ie)
            do ig=1, NumberofSelectedOrbitals_groups
               do iq=1, Nq
                  do i= 1, NumberofSelectedOrbitals(ig)
                     j= Num_wann*(iq-1)+ Selected_WannierOrbitals(ig)%iarray(i)
                     dos_selected(ie, ib, ig)= dos_selected(ie, ib, ig)+ abs(psi(j))**2
                  enddo ! sweep the selected orbitals
               enddo ! iq sweep the magnetic supercell
            enddo ! ig group of selected orbitals
         enddo ! ib sweep the eigenvalue
      endif


      call now(time_end)
   enddo !ib

#if defined (MPI)
   call mpi_allreduce(eigv,eigv_mpi,size(eigv),&
      mpi_dp,mpi_sum,mpi_cmw,ierr)
   call mpi_allreduce(dos_selected, dos_selected_mpi,size(dos_selected),&
      mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
   eigv_mpi= eigv
   dos_selected_mpi= dos_selected
#endif

   outfileindex= outfileindex+ 1
   if (cpuid.eq.0) then
      open(unit=outfileindex, file='landaulevel_B.dat')
      if (LandauLevel_wavefunction_calc) then
         write(outfileindex, '("#", a14, 2a15, a)')'Phi per cell', 'B (Tesla)', ' Eig of LL', ' Weight on the selected orbitals'
         do j=1, neval
            do i=1, Nmag1+1
               write(outfileindex,'(4f19.7)')mag(i)/2d0/pi, mag_Tesla(i), eigv_mpi(j, i),  &
                  (dos_selected_mpi(j, i, ig), ig=1, NumberofSelectedOrbitals_groups)
            enddo
            write(outfileindex , *)''
         enddo
      else
         write(outfileindex, '("#", a14, 2a15)')'Phi per cell', 'B (Tesla)', ' Eig of LL'
         do j=1, neval
            do i=1, Nmag1+1
               write(outfileindex,'(3f19.7)')mag(i)/2d0/pi, mag_Tesla(i), eigv_mpi(j, i)
            enddo
            write(outfileindex , *)''
         enddo
      endif

      close(outfileindex)
      write(stdout,*) 'calculate landau level done'
   endif

   outfileindex= outfileindex+ 1
   if (cpuid==0) then
      open(unit=outfileindex, file='landaulevel_B.gnu')
      write(outfileindex, '(a)')"set encoding iso_8859_1"
      write(outfileindex, '(a)')'#set terminal  postscript enhanced color'
      write(outfileindex, '(a)')"#set output 'landaulevel_B.eps'"
      write(outfileindex, '(3a)')'#set terminal  pngcairo truecolor enhanced', &
         '  font ",60" size 1920, 1680'
      write(outfileindex, '(3a)')'set terminal  png truecolor enhanced', &
         '  font ",60" size 1920, 1680'
      write(outfileindex, '(a)')"set output 'landaulevel_B.png'"
      write(outfileindex, '(a)')'unset key'
      write(outfileindex, '(a)')'set pointsize 0.8'
      write(outfileindex, '(a)')'set palette rgb 21,22,23'
      write(outfileindex, '(a)')'set border lw 3 '
      write(outfileindex, '(a)')'#set xtics font ",36"'
      write(outfileindex, '(a)')'#set ytics font ",36"'
      write(outfileindex, '(a)')'#set xlabel font ",36"'
      write(outfileindex, '(a)')'set xlabel "Phi per cell"'
      write(outfileindex,*) '#set xlabel "{/Symbol F}/{/Symbol F}_0"'
      write(outfileindex, '(a)')'set ylabel "Energy (eV)"'
      write(outfileindex, '(a, i6,a)')'set title "Landau level with Nq=', Nq, '"'
      write(outfileindex, '(a)')'#set xtics offset 0, -1'
      write(outfileindex, '(a)')'#set ylabel offset -1, 0 '
      write(outfileindex, '(a)')'rgb(r,g,b) = int(r)*65536 + int(g)*256 + int(b)'
      if (LandauLevel_wavefunction_calc) then
         write(outfileindex, '(2a)')"plot 'landaulevel_B.dat' u 1:3:(rgb(255,255-255*$3, 3)) ",  &
            "w p  pt 7  ps 1 lc rgb variable"
      else
         write(outfileindex, '(2a)')"plot 'landaulevel_B.dat' u 1:3",  &
            " w p  pt 7  ps 1"
      endif
      close(outfileindex)
   endif

#if defined (MPI)
   call mpi_barrier(mpi_cmw, ierr)
#endif

   deallocate( acoo)
   deallocate( jcoo)
   deallocate( icoo)
   deallocate( W)
   deallocate( eigv)
   deallocate( eigv_mpi)
   deallocate( mag)

   return

end subroutine sparse_landau_level_B


subroutine sparse_landau_level_k
   use para
   use sparse
   implicit none


   !> magnetic supercell size, perpendicular to the magnetic field
   !> Ndimq= Nq* Num_wann
   integer :: Nq, Ndimq
   integer :: Nmag,Nmag1

   !> some temporary integers
   integer :: ik, ia1, ia2, i, j, ierr, ib, iq, ig

   real(dp) :: B0, B0Tesla, B0Tesla_quantumflux_magsupcell
   real(dp) :: theta

   real(dp) :: t1, temp, dis, dis1

   ! wave vector
   real(dp) :: k3(3)

   !> dim= Ndimq, knv3
   real(dp), allocatable :: W(:)
   real(dp), allocatable :: eigv(:, :)
   real(dp), allocatable :: eigv_mpi(:, :)

   real(dp) :: emin, emax
   real(dp) :: time_start, time_end, time_start0

   !> dim= Ndimq*Ndimq
   integer :: nnzmax, nnz
   complex(dp), allocatable :: acoo(:)
   integer, allocatable :: jcoo(:)
   integer, allocatable :: icoo(:)

   !> eigenvector of the sparse matrix acoo. Dim=(Ndimq, neval)
   complex(dp), allocatable :: psi(:)
   complex(dp), allocatable :: zeigv(:, :)

   !> print the weight for the Selected_WannierOrbitals
   real(dp), allocatable :: dos_selected(:, :, :)
   real(dp), allocatable :: dos_selected_mpi(:, :, :)

   real(dp), allocatable :: dos_l_selected(:, :, :)
   real(dp), allocatable :: dos_l_selected_mpi(:, :, :)
   real(dp), allocatable :: dos_r_selected(:, :, :)
   real(dp), allocatable :: dos_r_selected_mpi(:, :, :)

   !number of ARPACK eigenvalues
   integer :: neval

   ! number of Arnoldi vectors
   integer :: nvecs

   !shift-invert sigma
   complex(dp) :: sigma

   !> time measurement
   real(dp) :: time1, time2, time3

   Nq= Magq
   Nmag= Nq
   Ndimq= Num_wann* Nq
   if (nnzmax_input<0)then
      nnzmax= Num_wann*(2*ijmax+1)*Ndimq+Ndimq
      if(Is_Sparse_Hr) nnzmax=splen*nq+Ndimq
   else
      nnzmax= nnzmax_input
   endif
 

   if (NumSelectedEigenVals==0) NumSelectedEigenVals=OmegaNum
   neval=NumSelectedEigenVals
   if (neval>=Ndimq) neval= Ndimq- 2

   !> ncv
   nvecs=int(3*neval)

   if (nvecs<50) nvecs= 50
   if (nvecs>Ndimq) nvecs= Ndimq

   sigma=(1d0,0d0)*iso_energy

   allocate( acoo(nnzmax), stat= ierr)
   call printallocationinfo('acoo', ierr)
   allocate( jcoo(nnzmax), stat= ierr)
   call printallocationinfo('jcoo', ierr)
   allocate( icoo(nnzmax), stat= ierr)
   call printallocationinfo('icoo', ierr)
   allocate( W( neval), stat= ierr)
   allocate( eigv( neval, nk3_band))
   call printallocationinfo('eigv', ierr)
   allocate( eigv_mpi( neval, nk3_band), stat= ierr)
   call printallocationinfo('eigv_mpi', ierr)
   allocate( psi(ndimq))
   allocate( zeigv(ndimq,nvecs), stat= ierr)
   call printallocationinfo('zeigv', ierr)
   allocate( dos_selected     (neval,   nk3_band, NumberofSelectedOrbitals_groups), stat= ierr)
   call printallocationinfo('dos_selected', ierr)
   allocate( dos_selected_mpi (neval,   nk3_band, NumberofSelectedOrbitals_groups), stat= ierr)
   call printallocationinfo('dos_selected_mpi', ierr)
   if (landau_chern_calc) then
      allocate( dos_l_selected     (neval,   nk3_band, NumberofSelectedOrbitals_groups), stat= ierr)
      call printallocationinfo('dos_l_selected', ierr)
      allocate( dos_l_selected_mpi (neval,   nk3_band, NumberofSelectedOrbitals_groups), stat= ierr)
      call printallocationinfo('dos_l_selected_mpi', ierr)
      allocate( dos_r_selected     (neval,   nk3_band, NumberofSelectedOrbitals_groups), stat= ierr)
      call printallocationinfo('dos_r_selected', ierr)
      allocate( dos_r_selected_mpi (neval,   nk3_band, NumberofSelectedOrbitals_groups), stat= ierr)
      call printallocationinfo('dos_r_selected_mpi', ierr)
      dos_l_selected= 0d0
      dos_l_selected_mpi= 0d0
      dos_r_selected= 0d0
      dos_r_selected_mpi= 0d0
   endif
   dos_selected= 0d0
   dos_selected_mpi= 0d0

   if (cpuid==0) then
      write(stdout, '(a, I16)')' nnzmax is ',  &
         nnzmax
      write(stdout, '(a, I16)')' Dimension of magnetic supercell Ndimq is ',  &
         ndimq
      write(stdout, '(a, I16)')' Number of ritz vectors nvecs is ',  &
         nvecs
      write(stdout, '(a, f20.1, a)')' Memory for acoo, icoo, jcoo is ',  &
         nnzmax*(16d0+4d0+4d0)/1024d0/1024d0, ' MB'
      write(stdout, '(a, f20.1, a)')' Memory for zeigv matrix is ',  &
         ndimq/1024d0/1024d0*nvecs*16d0, ' MB'
      write(stdout, '(a, f20.1, a)')' Memory for eigv matrix is ',  &
         neval*nk3_band*16d0/1024d0/1024d0, ' MB'
      if (landau_chern_calc) then
         write(stdout, '(a, f20.1, a)')' Memory for dos_selected is ',  &
            neval*nk3_band*6*NumberofSelectedOrbitals_groups*8d0/1024d0/1024d0, ' MB'
      else
         write(stdout, '(a, f20.1, a)')' Memory for dos_selected is ',  &
            neval*nk3_band*2*NumberofSelectedOrbitals_groups*8d0/1024d0/1024d0, ' MB'
      endif
      write(stdout, *)' '
   endif
   eigv_mpi= 0d0
   eigv    = 0d0
   acoo= 0d0


   !> deal with the magnetic field
   !> first transform the Bx By into B*Cos\theta, B*Sin\theta
   if (abs(By)<1e-8) then
      theta= 0d0
   elseif (By>0) then
      theta = atan(Bx/By)
   else
      theta = atan(Bx/By)+ pi
   endif
   if (abs(theta)<1e-8 .and. Bx<0 ) theta= theta+ pi

   !> get the shortest bond in the home unit cell
!    if (dis1< 1e-9) stop 'something wrong with the atom position'
   !> the flux in the unit cell not the magnetic supercell
   B0= 2d0*pi/dble(Nq)*Magp
   Bx= B0* dcos(theta)
   By= B0* dsin(theta)


   !> transform it into Tesla
   !> B=2*pi*\phi_0/S0, where \phi_0 is the quantum flux, S0 is the projected area
   !> \phi_0 = h/2e  h=6.62607004*1E-34, e= 1.6*1E-19
   !> B0Tesla_quantumflux_magsupcell is the magnetic field that makes the flux through the magnetic
   !> supercell to be one quantum flux
   B0Tesla_quantumflux_magsupcell= 6.62607004*1E-34/1.6021766208*1E19/MagneticSuperProjectedArea*1E20
   B0Tesla= 6.62607004*1E-34/1.6021766208*1E19/MagneticSuperProjectedArea*1E20*Magp

   if (cpuid==0) then
      write(stdout, '(a, 2E18.8)')' Magnetic field B in Tesla that makes the flux through '
      write(stdout, '(a, 2E18.8)')' the magnetic supercell to be one quantum flux : ', B0Tesla_quantumflux_magsupcell
      write(stdout, '(a, 2E18.8)')' Magnetic field B in Tesla that fits to Magp : ', B0Tesla
      write(stdout, '(a, 2f18.8)')' Magnetic field Bx, By= ', Bx, By
   endif


   !> calculate the landau levels along special k line
   k3= 0
   time_start= 0d0
   time_start0= 0d0
   call now(time_start0)
   time_start= time_start0
   time_end  = time_start0
   do ik=1+ cpuid, nk3_band, num_cpu
      if (cpuid.eq.0) &
         write(stdout, '(2a, i9, "  /", i10, a, f10.1, "s", a, f10.1, "s")') &
         'In sparse_landau_level_k ', ' ik/NK ', ik, nk3_band, &
         ' time elapsed: ', time_end-time_start0, &
         ' time left: ', ((nk3_band-ik)/num_cpu)*(time_end-time_start)

      call now(time_start)

      if (cpuid==0) write(stdout, '(a, 2i10)') 'LandauLevel_k_calc', ik,nk3_band
      k3 = kpath_3d(:, ik)
      nnz= nnzmax
      call now(time1)
      acoo= 0d0
      if(.not. Is_Sparse_Hr) then
         call ham_3Dlandau_sparse1(nnz, Ndimq, Nq, k3, acoo,jcoo,icoo)
      else
         call ham_3Dlandau_sparseHR(nnz, Ndimq, Nq, k3, acoo,jcoo,icoo)
      end if
      acoo=acoo/eV2Hartree
      call now(time2)

      !> diagonalization by call zheev in lapack
      W= 0d0
      call arpack_sparse_coo_eigs(Ndimq,nnzmax,nnz,acoo,jcoo,icoo,neval,nvecs,W,sigma, zeigv, LandauLevel_wavefunction_calc)

      call now(time3)
      eigv(1:neval, ik)= W(1:neval)

      !> calculate the weight on the selected orbitals
      do ib= 1, neval
         psi(:)= zeigv(:, ib)  !> the eigenvector of ib'th band
         do ig=1, NumberofSelectedOrbitals_groups
            do iq=1, Nq
               do i= 1, NumberofSelectedOrbitals(ig)
                  j= Num_wann*(iq-1)+ Selected_WannierOrbitals(ig)%iarray(i)
                  dos_selected(ib, ik, ig)= dos_selected(ib, ik, ig)+ abs(psi(j))**2

               enddo ! sweep the selected orbitals
            enddo ! iq sweep the magnetic supercell

            if (landau_chern_calc) then
               do iq=1, 2  ! edge states
                  if (iq>Nq) cycle
                  do i= 1, NumberofSelectedOrbitals(ig)
                     j= Num_wann*(iq-1)+ Selected_WannierOrbitals(ig)%iarray(i)
                     dos_l_selected(ib, ik, ig)= dos_l_selected(ib, ik, ig)+ abs(psi(j))**2
                  enddo ! sweep the selected orbitals
               enddo ! iq sweep the magnetic supercell
               do iq=Nq-1, Nq ! edge states
                  if (iq<1) cycle
                  do i= 1, NumberofSelectedOrbitals(ig)
                     j= Num_wann*(iq-1)+ Selected_WannierOrbitals(ig)%iarray(i)
                     dos_r_selected(ib, ik, ig)= dos_r_selected(ib, ik, ig)+ abs(psi(j))**2
                  enddo ! sweep the selected orbitals
               enddo ! iq sweep the magnetic supercell
            endif
         enddo ! ig

      enddo ! ib sweep the eigenvalue

      if (cpuid==0)write(stdout, '(a, f20.2, a)')'  >> Time cost for constructing H: ', time2-time1, ' s'
      if (cpuid==0)write(stdout, '(a, f20.2, a)')'  >> Time cost for diagonalize H: ', time3-time2, ' s'
      call now(time_end)
   enddo !ik

#if defined (MPI)
   call mpi_allreduce(eigv,eigv_mpi,size(eigv),&
      mpi_dp,mpi_sum,mpi_cmw,ierr)
   call mpi_allreduce(dos_selected, dos_selected_mpi,size(dos_selected),&
      mpi_dp,mpi_sum,mpi_cmw,ierr)

   if (landau_chern_calc) then
      call mpi_allreduce(dos_l_selected, dos_l_selected_mpi,size(dos_l_selected),&
         mpi_dp,mpi_sum,mpi_cmw,ierr)
      call mpi_allreduce(dos_r_selected, dos_r_selected_mpi,size(dos_r_selected),&
         mpi_dp,mpi_sum,mpi_cmw,ierr)
   endif
#else
   eigv_mpi= eigv
   dos_selected_mpi= dos_selected
   if (landau_chern_calc) then
      dos_l_selected_mpi= dos_l_selected
      dos_r_selected_mpi= dos_r_selected
   endif
#endif

   !> minimum and maximum value of energy bands
   emin= minval(eigv_mpi)-0.5d0
   emax= maxval(eigv_mpi)+0.5d0


   outfileindex= outfileindex+ 1
   if (cpuid.eq.0) then
      open(unit=outfileindex, file='landaulevel_k.dat')
      if (LandauLevel_wavefunction_calc) then
         write(outfileindex, '("#", a14, f16.8)')'Magnetic field in Tesla: ', B0Tesla
         write(outfileindex, '("#", a14, a15, a)')'k ', ' Eig of LL', ' Weight on the selected orbitals'
         do j=1, neval
            do i=1,nk3_band
               if (landau_chern_calc) then
                  write(outfileindex,'(20f16.8)')k3len_mag(i)*Angstrom2atomic, eigv_mpi(j, i), &
                     (dos_l_selected_mpi(j, i, ig), dos_r_selected_mpi(j, i, ig), &
                     ig=1, NumberofSelectedOrbitals_groups)
               else
                  write(outfileindex,'(20f16.8)')k3len_mag(i)*Angstrom2atomic, eigv_mpi(j, i), &
                     (dos_selected_mpi(j, i, ig), ig=1, NumberofSelectedOrbitals_groups)
               endif
            enddo
            write(outfileindex , *)''
         enddo
      else
         write(outfileindex, '("#", a14, f16.8)')'Magnetic field in Tesla: ', B0Tesla
         write(outfileindex, '("#", a14, a15)')'k ', ' Eig of LL'
         do j=1, neval
            do i=1,nk3_band
               write(outfileindex,'(20f16.8)')K3len_mag(i)*Angstrom2atomic, eigv_mpi(j, i)
            enddo
            write(outfileindex , *)''
         enddo
      endif
      close(outfileindex)
      write(stdout,*) 'calculate landau level done'
   endif

   outfileindex= outfileindex+ 1
   if (cpuid==0) then
      open(unit=outfileindex, file='landaulevel_k.gnu')
      write(outfileindex, '(a)')"set encoding iso_8859_1"
      write(outfileindex, '(a)')'#set terminal  postscript enhanced color'
      write(outfileindex, '(a)')"#set output 'landaulevel_k.eps'"
      write(outfileindex, '(3a)')'#set terminal  pngcairo truecolor enhanced', &
         '  font ",60" size 1920, 1680'
      write(outfileindex, '(3a)')'set terminal  png truecolor enhanced', &
         '  font ",60" size 1920, 1680'
      write(outfileindex, '(a)')"set output 'landaulevel_k.png'"
      write(outfileindex, '(a)')'unset key'
      write(outfileindex, '(a)')'set pointsize 0.8'
      write(outfileindex, '(a)')'set palette rgb 21,22,23'
      write(outfileindex, '(a)')'set border lw 3 '
      write(outfileindex, '(a)')'#set xtics font ",36"'
      write(outfileindex, '(a)')'#set ytics font ",36"'
      write(outfileindex, '(a)')'#set xlabel font ",36"'
      write(outfileindex, '(a)')'set ylabel "Energy (eV)"'
      write(outfileindex, '(a, f10.5, a)')'set xrange [0: ', maxval(k3len_mag*Angstrom2atomic), ']'
      write(outfileindex, '(a, i6, a, i6, a)')'set title "Landau level with p/q=', magp, "/", Nq,  '"'
      write(outfileindex, '(a)')'#set xtics offset 0, -1'
      write(outfileindex, '(a)')'#set ylabel offset -1, 0 '
      write(outfileindex, '(a)')'rgb(r,g,b) = int(r)*65536 + int(g)*256 + int(b)'
      write(outfileindex, 202, advance="no") (k3line_name(i), k3line_mag_stop(i)*Angstrom2atomic, i=1, nk3lines)
      write(outfileindex, 203)k3line_name(nk3lines+1), k3line_mag_stop(nk3lines+1)*Angstrom2atomic
      write(outfileindex, '(a, f10.5, a, f10.5, a)')'set yrange [', emin, ':', emax, ']'

      do i=1, nk3lines-1
         write(outfileindex, 204)k3line_mag_stop(i+1)*Angstrom2atomic, emin, k3line_mag_stop(i+1)*Angstrom2atomic, emax
      enddo

202   format('set xtics (',:20('"',A3,'" ',F10.5,','))
203   format(A3,'" ',F10.5,')')
204   format('set arrow from ',F10.5,',',F10.5,' to ',F10.5,',',F10.5, ' nohead')

      if (LandauLevel_wavefunction_calc) then
         if (landau_chern_calc) then
            write(outfileindex,'(2a)') 'set palette defined ( -1  "blue", ', &
               '0 "grey", 1 "red" )'
            write(outfileindex, '(2a)')"plot 'landaulevel_k.dat' u 1:2:($4-$3)",  &
               "w p pt 7  ps 2 lc palette"
         else
            write(outfileindex, '(2a)')"plot 'landaulevel_k.dat' u 1:2:(rgb(255,255-255*$3, 3)) ",  &
               " w p  pt 7  ps 2 lc rgb variable"
         endif
      else
         write(outfileindex, '(2a)')"plot 'landaulevel_k.dat' u 1:2",  &
            " w p  pt 7  ps 1"
      endif
      close(outfileindex)
   endif



#if defined (MPI)
   call mpi_barrier(mpi_cmw, ierr)
#endif

   deallocate( acoo)
   deallocate( jcoo)
   deallocate( icoo)
   deallocate( W)
   deallocate( eigv)
   deallocate( eigv_mpi)
   deallocate( zeigv)
   deallocate( dos_selected)
   deallocate( dos_selected_mpi)

   return
end subroutine sparse_landau_level_k


!> calculate density of states
subroutine sparse_landau_dos
   use para
   use sparse
   implicit none

   integer,parameter :: nkk=7
   !> magnetic supercell size, perpendicular to the magnetic field
   integer :: Nq
   integer :: ia1
   integer :: ia2
   integer :: ib
   integer :: i, j



   !> Ndimq= Nq* Num_wann
   integer :: Ndimq

   integer :: Nmag,Nmag1

   integer :: ierr

   real(dp) :: B0
   real(dp) :: theta
   real(dp) :: eta_broadening

   integer :: nnzmax, nnz
   real(dp) :: dis, dis1
   ! wave vector
   real(dp) :: k3(3),dk(nkk*2-1)
   real(dp) :: times,timee
   !> dim= Ndimq, knv3
   real(dp), allocatable :: mag(:)
   real(dp), allocatable :: W(:)
   real(dp), allocatable :: eigv(:, :)
   real(dp), allocatable :: eigv_mpi(:, :)
   real(dp), allocatable :: eigv_mpi2(:, :)
   complex(dp), allocatable :: zeigv(:, :)

   !> dim= Ndimq*Ndimq
   complex(dp), allocatable :: acoo(:)
   integer, allocatable :: jcoo(:)
   integer, allocatable :: icoo(:)

   !number of ARPACK eigenvalues
   integer :: neval
   ! number of Arnoldi vectors
   integer :: nvecs
   !shift-invert sigma
   complex(dp) :: sigma=(0d0,0d0)


   !k-path calculating

   integer :: knv3
   integer :: ikx,iky,ikz
   integer :: ik,ie,NE,iv
   real(dp) :: emin,x
   real(dp) :: emax
   real(dp),allocatable::omega(:),dos_mpi(:,:),dos(:,:)
   real(dp), external :: delta



   allocate(zeigv(ndimq,nvecs))
   knv3= Nk1*Nk2*Nk3
   NE= OmegaNum
   emin= OmegaMin
   emax= OmegaMax
   eta_broadening= (emax- emin)/ dble(NE)*1d0
   if (cpuid.eq.0) write(stdout,*) 'eta_broadening=',eta_broadening
   allocate(dos_mpi(knv3,NE))
   allocate(dos(knv3,NE))
   allocate(omega(NE))


   if (NumSelectedEigenVals==0) NumSelectedEigenVals=Ndimq
   neval=NumSelectedEigenVals
   nvecs=2*neval
   sigma=(1d0,0d0)*iso_energy
   Nq= Magq
   Nmag= Nq
   Nmag1=20
   Ndimq= Num_wann* Nq
   nnzmax= Num_wann*(2*ijmax+1)*Ndimq+Ndimq
   nnz=nnzmax
   allocate( acoo(nnzmax))
   allocate( jcoo(nnzmax))
   allocate( icoo(nnzmax))
   allocate( W( neval))
   allocate( eigv( neval, Nmag1))
   allocate( eigv_mpi( neval, Nmag1))
   allocate( eigv_mpi2( neval, Nmag1))
   allocate( mag(Nq))
   mag= 0d0
   eigv_mpi= 0d0
   eigv    = 0d0
   acoo=0d0
   jcoo=0
   icoo=0
   if (cpuid.eq.0) write(stdout,*) 'shift energy sigma=',real(sigma),Ndimq,neval
   do i=1,nkk*2-1
      dk(i)=(i-nkk)/(2d0*nkk-1d0)
   end do
   !> deal with the magnetic field
   !> first transform the Bx By into B*Cos\theta, B*Sin\theta
   if (abs(By)<1e-8) then
      theta= 0d0
   elseif (By>0) then
      theta = atan(Bx/By)
   else
      theta = atan(Bx/By)+ pi
   endif
   if (abs(theta)<1e-8 .and. Bx<0 ) theta= theta+ pi

   B0=2d0*pi/(1000)
   B0=abs(B0)


   !    do kk1=1,nkk*2-1
   !        do kk2=1,nkk*2-1
   !            k3=[0d0,dk(kk1),dk(kk2)]
   !    k3=[0d0,0.33d0,0.33d0]
   !> calculate the landau levels along special k line
   eigv_mpi2=0
   !        do kk1=5,0,-1

   do ib=1, Nmag
      mag(ib)= B0* (ib)!+2d0*kk1*pi
   enddo
   do ie=1, NE
      omega(ie)= emin+ (emax-emin)* (ie-1d0)/dble(NE-1)
   enddo ! ie

   !    do ib=1+ cpuid, Nmag1, num_cpu
   dos_mpi=0
   dos=0
   do ik=1+cpuid,knv3,num_cpu
      ib=1
      ikx= (ik-1)/(nk2*nk3)+1
      iky= ((ik-1-(ikx-1)*Nk2*Nk3)/nk3)+1
      ikz= (ik-(iky-1)*Nk3- (ikx-1)*Nk2*Nk3)
      k3= K3D_start_cube+ K3D_vec1_cube*(ikx-1)/dble(nk1)  &
         + K3D_vec2_cube*(iky-1)/dble(nk2)  &
         + K3D_vec3_cube*(ikz-1)/dble(nk3)
      if(cpuid==0) call now(times)
      Bx= mag(ib)* Cos(theta)
      By= mag(ib)* Sin(theta)
      if (cpuid==0) print *, ib, Nmag
      call ham_3Dlandau_sparse1(nnz, Ndimq, Nq, k3, acoo,jcoo,icoo)
      acoo=acoo/eV2Hartree

      !> diagonalization by call zheev in lapack
      W= 0d0
      !        call eigensystem_c( 'N', 'U', Ndimq ,ham_landau, W)
      call arpack_sparse_coo_eigs(Ndimq,nnzmax,nnz,acoo,jcoo,icoo,neval,nvecs,W,sigma, zeigv, LandauLevel_wavefunction_calc)
      do ie= 1, NE
         do iv= 1, neval
            x= omega(ie)- W(iv)
            dos_mpi(ik,ie) = dos_mpi(ik,ie)+ delta(eta_broadening, x)
         enddo ! iv
      enddo ! ie
      if(cpuid==0) then
         call now(timee)
!            write(*,*) timee-times
      endif

   enddo !ik


   do ie=1,NE
      !write(231,*) omega(ie),mag(ib),(dos_mpi(ik,ie),ik=1,knv3)
   end do


#if defined (MPI)
   call mpi_allreduce(dos_mpi,dos,size(dos),&
      mpi_dp,mpi_sum,mpi_cmw,ierr)
#else
   dos= dos_mpi
#endif
!    do ie=1,NE
!        write(231,*) omega(ie),mag(ib),dos(ie)
!    end do

   !    end do
!    write(*,*) 'all eigvalue'


#if defined (MPI)
   call mpi_barrier(mpi_cmw, ierr)
#endif
   go to 21

21 return

end subroutine



subroutine sparse_export_maghr
   use para
   use sparse
   implicit none


   !> magnetic supercell size, perpendicular to the magnetic field
   !> Ndimq= Nq* Num_wann
   integer :: Nq, Ndimq
   integer :: Nmag,Nmag1

   !> some temporary integers
   integer :: ik, ia1, ia2, i, j, ierr, ib, iq, ig

   real(dp) :: B0, B0Tesla, B0Tesla_quantumflux_magsupcell
   real(dp) :: theta

   real(dp) :: t1, temp, dis, dis1

   ! wave vector
   real(dp) :: k3(3)

   !> dim= Ndimq, knv3
   real(dp), allocatable :: W(:)
   real(dp), allocatable :: eigv(:, :)
   real(dp), allocatable :: eigv_mpi(:, :)

   real(dp) :: emin, emax
   real(dp) :: time_start, time_end, time_start0

   !> dim= Ndimq*Ndimq
   integer :: nnzmax, nnz
   complex(dp), allocatable :: acoo(:)
   integer, allocatable :: jcoo(:)
   integer, allocatable :: icoo(:)

   !> eigenvector of the sparse matrix acoo. Dim=(Ndimq, neval)
   complex(dp), allocatable :: psi(:)
   complex(dp), allocatable :: zeigv(:, :)

   !> print the weight for the Selected_WannierOrbitals
   real(dp), allocatable :: dos_selected(:, :, :)
   real(dp), allocatable :: dos_selected_mpi(:, :, :)

   real(dp), allocatable :: dos_l_selected(:, :, :)
   real(dp), allocatable :: dos_l_selected_mpi(:, :, :)
   real(dp), allocatable :: dos_r_selected(:, :, :)
   real(dp), allocatable :: dos_r_selected_mpi(:, :, :)

   !number of ARPACK eigenvalues
   integer :: neval

   ! number of Arnoldi vectors
   integer :: nvecs

   !shift-invert sigma
   complex(dp) :: sigma

   !> time measurement
   real(dp) :: time1, time2, time3

   Nq= Magq
   Nmag= Nq
   Ndimq= Num_wann* Nq
   nnzmax= Num_wann*(2*ijmax+1)*Ndimq+Ndimq
   if(Is_Sparse_Hr) nnzmax=splen*nq+Ndimq
   if (NumSelectedEigenVals==0) NumSelectedEigenVals=OmegaNum
   neval=NumSelectedEigenVals
   if (neval>=Ndimq) neval= Ndimq- 2

   !> ncv
   nvecs=int(3*neval)

   if (nvecs<50) nvecs= 50
   if (nvecs>Ndimq) nvecs= Ndimq

   sigma=(1d0,0d0)*iso_energy

   allocate( acoo(nnzmax), stat= ierr)
   call printallocationinfo('acoo', ierr)
   allocate( jcoo(nnzmax), stat= ierr)
   call printallocationinfo('jcoo', ierr)
   allocate( icoo(nnzmax), stat= ierr)
   call printallocationinfo('icoo', ierr)
   allocate( W( neval), stat= ierr)
   allocate( eigv( neval, nk3_band))
   call printallocationinfo('eigv', ierr)
   allocate( eigv_mpi( neval, nk3_band), stat= ierr)
   call printallocationinfo('eigv_mpi', ierr)
   allocate( psi(ndimq))
   allocate( zeigv(ndimq,nvecs), stat= ierr)
   call printallocationinfo('zeigv', ierr)
   allocate( dos_selected     (neval,   nk3_band, NumberofSelectedOrbitals_groups), stat= ierr)
   call printallocationinfo('dos_selected', ierr)
   allocate( dos_selected_mpi (neval,   nk3_band, NumberofSelectedOrbitals_groups), stat= ierr)
   call printallocationinfo('dos_selected_mpi', ierr)
   if (landau_chern_calc) then
      allocate( dos_l_selected     (neval,   nk3_band, NumberofSelectedOrbitals_groups), stat= ierr)
      call printallocationinfo('dos_l_selected', ierr)
      allocate( dos_l_selected_mpi (neval,   nk3_band, NumberofSelectedOrbitals_groups), stat= ierr)
      call printallocationinfo('dos_l_selected_mpi', ierr)
      allocate( dos_r_selected     (neval,   nk3_band, NumberofSelectedOrbitals_groups), stat= ierr)
      call printallocationinfo('dos_r_selected', ierr)
      allocate( dos_r_selected_mpi (neval,   nk3_band, NumberofSelectedOrbitals_groups), stat= ierr)
      call printallocationinfo('dos_r_selected_mpi', ierr)
      dos_l_selected= 0d0
      dos_l_selected_mpi= 0d0
      dos_r_selected= 0d0
      dos_r_selected_mpi= 0d0
   endif
   dos_selected= 0d0
   dos_selected_mpi= 0d0

   if (cpuid==0) then
      write(stdout, '(a, I16)')' nnzmax is ',  &
         nnzmax
      write(stdout, '(a, I16)')' Dimension of magnetic supercell Ndimq is ',  &
         ndimq
      write(stdout, '(a, I16)')' Number of ritz vectors nvecs is ',  &
         nvecs
      write(stdout, '(a, f20.1, a)')' Memory for acoo, icoo, jcoo is ',  &
         nnzmax*(16d0+4d0+4d0)/1024d0/1024d0, ' MB'
      write(stdout, '(a, f20.1, a)')' Memory for zeigv matrix is ',  &
         ndimq/1024d0/1024d0*nvecs*16d0, ' MB'
      write(stdout, '(a, f20.1, a)')' Memory for eigv matrix is ',  &
         neval*nk3_band*16d0/1024d0/1024d0, ' MB'
      if (landau_chern_calc) then
         write(stdout, '(a, f20.1, a)')' Memory for dos_selected is ',  &
            neval*nk3_band*6*NumberofSelectedOrbitals_groups*8d0/1024d0/1024d0, ' MB'
      else
         write(stdout, '(a, f20.1, a)')' Memory for dos_selected is ',  &
            neval*nk3_band*2*NumberofSelectedOrbitals_groups*8d0/1024d0/1024d0, ' MB'
      endif
      write(stdout, *)' '
   endif
   eigv_mpi= 0d0
   eigv    = 0d0
   acoo= 0d0


   !> deal with the magnetic field
   !> first transform the Bx By into B*Cos\theta, B*Sin\theta
   if (abs(By)<1e-8) then
      theta= 0d0
   elseif (By>0) then
      theta = atan(Bx/By)
   else
      theta = atan(Bx/By)+ pi
   endif
   if (abs(theta)<1e-8 .and. Bx<0 ) theta= theta+ pi

   !> get the shortest bond in the home unit cell
!    if (dis1< 1e-9) stop 'something wrong with the atom position'
   !> the flux in the unit cell not the magnetic supercell
   B0= 2d0*pi/dble(Nq)*Magp
   Bx= B0* dcos(theta)
   By= B0* dsin(theta)


   !> transform it into Tesla
   !> B=2*pi*\phi_0/S0, where \phi_0 is the quantum flux, S0 is the projected area
   !> \phi_0 = h/2e  h=6.62607004*1E-34, e= 1.6*1E-19
   !> B0Tesla_quantumflux_magsupcell is the magnetic field that makes the flux through the magnetic
   !> supercell to be one quantum flux
   B0Tesla_quantumflux_magsupcell= 6.62607004*1E-34/1.6021766208*1E19/MagneticSuperProjectedArea*1E20
   B0Tesla= 6.62607004*1E-34/1.6021766208*1E19/MagneticSuperProjectedArea*1E20*Magp

   if (cpuid==0) then
      write(stdout, '(a, 2E18.8)')' Magnetic field B in Tesla that makes the flux through '
      write(stdout, '(a, 2E18.8)')' the magnetic supercell to be one quantum flux : ', B0Tesla_quantumflux_magsupcell
      write(stdout, '(a, 2E18.8)')' Magnetic field B in Tesla that fits to Magp : ', B0Tesla
      write(stdout, '(a, 2f18.8)')' Magnetic field Bx, By= ', Bx, By
   endif


   !> calculate the landau levels along special k line
   k3= 0
   time_start= 0d0
   time_start0= 0d0
   call now(time_start0)
   time_start= time_start0
   time_end  = time_start0

   if (cpuid.eq.0) &
      write(stdout, '(2a, i9, "  /", i10, a, f10.1, "s", a, f10.1, "s")') &
      'In sparse_landau_level_k ', ' ik/NK ', ik, nk3_band, &
      ' time elapsed: ', time_end-time_start0, &
      ' time left: ', ((nk3_band-ik)/num_cpu)*(time_end-time_start)

   call now(time_start)


   if (cpuid==0) write(stdout, '(a, 2i10)') 'maghr', ik,nk3_band
   k3=Single_KPOINT_3D_DIRECT
   nnz= nnzmax
   call now(time1)
   if(.not. Is_Sparse_Hr) then
      call ham_3Dlandau_sparse1(nnz, Ndimq, Nq, k3, acoo,jcoo,icoo)
   else
      call ham_3Dlandau_sparseHR(nnz, Ndimq, Nq, k3, acoo,jcoo,icoo)
   end if
   acoo=acoo/eV2Hartree
   call now(time2)
   return




#if defined (MPI)
   call mpi_allreduce(eigv,eigv_mpi,size(eigv),&
      mpi_dp,mpi_sum,mpi_cmw,ierr)
   call mpi_allreduce(dos_selected, dos_selected_mpi,size(dos_selected),&
      mpi_dp,mpi_sum,mpi_cmw,ierr)

   if (landau_chern_calc) then
      call mpi_allreduce(dos_l_selected, dos_l_selected_mpi,size(dos_l_selected),&
         mpi_dp,mpi_sum,mpi_cmw,ierr)
      call mpi_allreduce(dos_r_selected, dos_r_selected_mpi,size(dos_r_selected),&
         mpi_dp,mpi_sum,mpi_cmw,ierr)
   endif
#else
   eigv_mpi= eigv
   dos_selected_mpi= dos_selected
   if (landau_chern_calc) then
      dos_l_selected_mpi= dos_l_selected
      dos_r_selected_mpi= dos_r_selected
   endif
#endif

   !> minimum and maximum value of energy bands
   emin= minval(eigv_mpi)-0.5d0
   emax= maxval(eigv_mpi)+0.5d0


   outfileindex= outfileindex+ 1
   if (cpuid.eq.0) then
      open(unit=outfileindex, file='landaulevel_k.dat')
      if (LandauLevel_wavefunction_calc) then
         write(outfileindex, '("#", a14, f16.8)')'Magnetic field in Tesla: ', B0Tesla
         write(outfileindex, '("#", a14, a15, a)')'k ', ' Eig of LL', ' Weight on the selected orbitals'
         do j=1, neval
            do i=1,nk3_band
               if (landau_chern_calc) then
                  write(outfileindex,'(20f16.8)')k3len_mag(i)*Angstrom2atomic, eigv_mpi(j, i), &
                     (dos_l_selected_mpi(j, i, ig), dos_r_selected_mpi(j, i, ig), &
                     ig=1, NumberofSelectedOrbitals_groups)
               else
                  write(outfileindex,'(20f16.8)')k3len_mag(i)*Angstrom2atomic, eigv_mpi(j, i), &
                     (dos_selected_mpi(j, i, ig), ig=1, NumberofSelectedOrbitals_groups)
               endif
            enddo
            write(outfileindex , *)''
         enddo
      else
         write(outfileindex, '("#", a14, f16.8)')'Magnetic field in Tesla: ', B0Tesla
         write(outfileindex, '("#", a14, a15)')'k ', ' Eig of LL'
         do j=1, neval
            do i=1,nk3_band
               write(outfileindex,'(20f16.8)')K3len_mag(i)*Angstrom2atomic, eigv_mpi(j, i)
            enddo
            write(outfileindex , *)''
         enddo
      endif
      close(outfileindex)
      write(stdout,*) 'calculate landau level done'
   endif

   outfileindex= outfileindex+ 1
   if (cpuid==0) then
      open(unit=outfileindex, file='landaulevel_k.gnu')
      write(outfileindex, '(a)')"set encoding iso_8859_1"
      write(outfileindex, '(a)')'#set terminal  postscript enhanced color'
      write(outfileindex, '(a)')"#set output 'landaulevel_k.eps'"
      write(outfileindex, '(3a)')'#set terminal  pngcairo truecolor enhanced', &
         '  font ",60" size 1920, 1680'
      write(outfileindex, '(3a)')'set terminal  png truecolor enhanced', &
         '  font ",60" size 1920, 1680'
      write(outfileindex, '(a)')"set output 'landaulevel_k.png'"
      write(outfileindex, '(a)')'unset key'
      write(outfileindex, '(a)')'set pointsize 0.8'
      write(outfileindex, '(a)')'set palette rgb 21,22,23'
      write(outfileindex, '(a)')'set border lw 3 '
      write(outfileindex, '(a)')'#set xtics font ",36"'
      write(outfileindex, '(a)')'#set ytics font ",36"'
      write(outfileindex, '(a)')'#set xlabel font ",36"'
      write(outfileindex, '(a)')'set ylabel "Energy (eV)"'
      write(outfileindex, '(a, f10.5, a)')'set xrange [0: ', maxval(k3len_mag*Angstrom2atomic), ']'
      write(outfileindex, '(a, i6, a, i6, a)')'set title "Landau level with p/q=', magp, "/", Nq,  '"'
      write(outfileindex, '(a)')'#set xtics offset 0, -1'
      write(outfileindex, '(a)')'#set ylabel offset -1, 0 '
      write(outfileindex, '(a)')'rgb(r,g,b) = int(r)*65536 + int(g)*256 + int(b)'
      write(outfileindex, 202, advance="no") (k3line_name(i), k3line_mag_stop(i)*Angstrom2atomic, i=1, nk3lines)
      write(outfileindex, 203)k3line_name(nk3lines+1), k3line_mag_stop(nk3lines+1)*Angstrom2atomic
      write(outfileindex, '(a, f10.5, a, f10.5, a)')'set yrange [', emin, ':', emax, ']'

      do i=1, nk3lines-1
         write(outfileindex, 204)k3line_mag_stop(i+1)*Angstrom2atomic, emin, k3line_mag_stop(i+1)*Angstrom2atomic, emax
      enddo

202   format('set xtics (',:20('"',A3,'" ',F10.5,','))
203   format(A3,'" ',F10.5,')')
204   format('set arrow from ',F10.5,',',F10.5,' to ',F10.5,',',F10.5, ' nohead')

      if (LandauLevel_wavefunction_calc) then
         if (landau_chern_calc) then
            write(outfileindex,'(2a)') 'set palette defined ( -1  "blue", ', &
               '0 "grey", 1 "red" )'
            write(outfileindex, '(2a)')"plot 'landaulevel_k.dat' u 1:2:($4-$3)",  &
               "w p pt 7  ps 2 lc palette"
         else
            write(outfileindex, '(2a)')"plot 'landaulevel_k.dat' u 1:2:(rgb(255,255-255*$3, 3)) ",  &
               " w p  pt 7  ps 2 lc rgb variable"
         endif
      else
         write(outfileindex, '(2a)')"plot 'landaulevel_k.dat' u 1:2",  &
            " w p  pt 7  ps 1"
      endif
      close(outfileindex)
   endif



#if defined (MPI)
   call mpi_barrier(mpi_cmw, ierr)
#endif

   deallocate( acoo)
   deallocate( jcoo)
   deallocate( icoo)
   deallocate( W)
   deallocate( eigv)
   deallocate( eigv_mpi)
   deallocate( zeigv)
   deallocate( dos_selected)
   deallocate( dos_selected_mpi)

   return
end subroutine sparse_export_maghr
