  subroutine psik_slab()
     ! Psik calculates the eigenvector for a given k point and selected bands
     ! for 2D slab system

     use para,only : Dp,Num_wann,Nslab, stdout, cpuid , outfileindex, Single_KPOINT_2D_DIRECT, &
        NumberofSelectedBands, Selected_band_index, Origin_cell
     
     implicit none 

     ! loop index
     integer     :: i,j, mdim, ib, iband

     ! wave vector 
     real(Dp)    :: k(2)
      
     ! eigenvalue 
     real(Dp), allocatable   :: W(:)
     
     !> norm of psi |\psi|^2
     real(Dp), allocatable   :: psi2  (:, :)
     !real(Dp), allocatable   :: psi_atom  (:, :)

     !> wave function for the given band and k point
     complex(Dp), allocatable:: psi(:, :)

     ! hamiltonian slab
     complex(Dp),allocatable :: hamk_slab(:,:), hamk_slab_t(:,:)

     mdim= Nslab*Num_wann

     allocate(psi2(Nslab, NumberofSelectedBands))
     ! allocate(psi_atom(Nslab*Origin_cell%Num_atoms, NumberofSelectedBands))
     allocate(psi (mdim, 1))
     allocate(W(mdim))
     allocate(hamk_slab(mdim,mdim), hamk_slab_t(mdim,mdim))
     psi2=0d0;  psi=0d0; W=0d0
     ! psi_atom = 0d0
     hamk_slab_t= 0d0; hamk_slab= 0d0


     !> Single_KPOINT_2D_DIRECT is provided in the wt.in or input.dat
     k= Single_KPOINT_2D_DIRECT
     hamk_slab=0.0d0 

     ! calculate Hamiltonian
     call ham_slab(k,hamk_slab)

     psi2=0.0d0
     do ib=1, NumberofSelectedBands
        iband = Selected_band_index(ib)
        if (iband > mdim.or. iband<0) then
           write(*, *)"ERROR: selected bands should be smaller than ", mdim
           stop
        endif

        !> diagonal hamk_slab, only get one eigenvalue
        W=0.0d0; psi= 0d0
        hamk_slab_t= hamk_slab
        call zheevx_pack('V', 'U', mdim, iband, iband, hamk_slab_t, W, psi)
       
        if (cpuid.eq.0) write(stdout,'(2X, a, i8, a, f16.6)') 'Eigenvalue for band ', iband, ' is', W(1)
   
        j=0
        do i=1,Nslab
           do j=1,Num_wann
              psi2(i, ib)=psi2(i, ib)+abs(psi((i-1)*Num_wann+j, 1))**2
           enddo
        enddo

       !it = 0
       !do i=1,Nslab
       !   do ia=1, Origin_cell%Num_atoms 
       !      it = it + 1
       !      do j=1, Origin_cell%nprojs(ia)
       !         psi_atom(it, ib)=psi_atom(it, ib) + abs(psi((i-1)*Num_wann+j+sum(Origin_cell%nprojs(1:ia))-Origin_cell%nprojs(1), 1))**2
       !      enddo !j
       !      if (SOC>0) then
       !         do j=1, Origin_cell%nprojs(ia)
       !            psi_atom(it, ib)=psi_atom(it, ib) + abs(psi((i-1)*Num_wann+j+sum(Origin_cell%nprojs(1:ia))-Origin_cell%nprojs(1),
       !            1)+Num_wann/2)**2
       !         enddo !j
       !      endif
       !   enddo !ia
       !enddo ! i

     enddo

     outfileindex= outfileindex+ 1
     if (cpuid==0) then
        open(unit=outfileindex, file='psi_abs.txt')
        write(outfileindex, '(a8, a5, 2000i16 )')'#islab', 'band', Selected_band_index(:)
        do i=1,Nslab
           write(outfileindex,'(i8, 5X, 2000f16.9)')i,psi2(i, :)
        enddo
        close(outfileindex)
        write(stdout,*) '<< Calculating psi done'
     endif
 
     deallocate(psi2)
     deallocate(psi )
     deallocate(hamk_slab)

     return
    
  end subroutine psik_slab

  subroutine psik_bulk()
     ! Psik calculates the eigenvector for a given k point and energe band
     ! for 3D bulk system
    
     use para,only : Dp,Num_wann, stdout, Numoccupied, cpuid, outfileindex
     implicit none 

! loop index
     integer     :: i,ik, ib
     integer     :: info

! wave vector 
     real(Dp)    :: k(3)
     real(Dp)    :: kpoints(3, 8)
      
! eigenvalue 
     real(Dp)    :: eigenvalue (Num_wann)
     
   
! energy dispersion
     complex(Dp) :: psi   (Num_wann)

! hamiltonian slab
     complex(Dp),allocatable :: CHamk(:,:)
     complex(Dp),allocatable :: eigenvector(:,:)


     allocate(CHamk(Num_wann,Num_wann))
     allocate(eigenvector(Num_wann,Num_wann))


     kpoints(:, 1)= (/0.0d0, 0.0d0, 0.0d0/)
     kpoints(:, 2)= (/0.5d0, 0.0d0, 0.0d0/)
     kpoints(:, 3)= (/0.0d0, 0.5d0, 0.0d0/)
     kpoints(:, 4)= (/0.0d0, 0.0d0, 0.5d0/)
     kpoints(:, 5)= (/0.5d0, 0.5d0, 0.5d0/)
     kpoints(:, 6)= (/0.5d0, 0.5d0, 0.0d0/)
     kpoints(:, 7)= (/0.0d0, 0.5d0, 0.5d0/)
     kpoints(:, 8)= (/0.5d0, 0.0d0, 0.5d0/)

     ib= Numoccupied
     outfileindex= outfileindex+ 1
     if (cpuid==0)open(unit=outfileindex, file='wavefunction.dat')
     do ik=1, 8

        k=kpoints(:, ik)
        ! calculate Hamiltonian
        call ham_bulk_latticegauge(k,Chamk)
        eigenvalue=0.0d0
        eigenvector=Chamk
        
        ! diagonal Chamk
        call eigensystem_c('V', 'U', Num_wann, eigenvector, eigenvalue)
       
        psi(:)=eigenvector(:, ib)
        if (cpuid==0)write(stdout,*) 'eigenvalue',info,eigenvalue(ib)
   
        if (cpuid==0)write(outfileindex, '(a,3f8.4)')'K point ', k
        do i=1, Num_wann
           if (cpuid==0)write(outfileindex,'(i5, 30f16.9)')i, eigenvector(i, ib-1), eigenvector(i, ib)
        enddo
        if (cpuid==0)write(outfileindex,*)' '
    
     enddo ! ik
        
     if (cpuid==0)close(outfileindex)
  end subroutine psik_bulk
