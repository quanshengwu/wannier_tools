
subroutine readinput
   ! Read in the control paramters from input.dat,
   ! and set default values if not specified in the input.dat
   !
   ! Constructed on 4/22/2020 by QS.Wu

   use wmpi
   use para
   implicit none

   integer :: stat
   character*12 :: fname='system.in'
   character*25 :: char_temp
   character*256 :: inline
   logical ::  exists, lfound

   integer :: i, j, m, it, iter, NN
   real(dp) :: theta, kD, t1, temp
   real(dp) :: k1(2), k2(2), kstart(2), kend(2)
   interlayercoupling_ratio_array_input= 1d0

   if (cpuid.eq.0) then
      write(stdout, *)' The system.in file should looks like:'
      write(stdout, *)'  &PARAMETERS '
      write(stdout, *)'  number_layers = 2     '
      write(stdout, *)' !twisted_index_m= 1      '
      write(stdout, *)'  twisted_angle_degree= 1.08   '
      write(stdout, *)'  twisted_angle_array_input = 0 1  '
      write(stdout, *)'  interlayercoupling_ratio_array= 1    '
      write(stdout, *)'  stacking_sequences_input = "A" "B"  '
      write(stdout, *)'  u_AA=0.0797  ! eV '
      write(stdout, *)'  u_AB=0.0975  ! eV '
      write(stdout, *)'  vppsigma=0.30d0  ! eV '
      write(stdout, *)'  Qcutoff = 4'
      write(stdout, *)'  Electric_field= 0.1 ! eV/Angstrom'
      write(stdout, *)'  lattice_constant_c = 3.36 ' 
      write(stdout, *)'  / '
   endif


   inquire(file=fname,exist=exists)
   if (exists)then
      if(cpuid==0)write(stdout,*) '  '
      if(cpuid==0)write(stdout,*) '>>>Read some paramters from system.in'
      open(unit=1001,file=fname,status='old')
   else
      if(cpuid==0)write(stdout,*)'file' ,fname, 'dosnot exist'
      stop
   endif

   !> initialization
   number_layers= 2
   twisted_index_m= -1
   twisted_angle_degree=1.08d0 ! degree
   twisted_angle_array_input(1:2)=(/0, 1/)
   stacking_sequences_input(1:2)=(/'A', 'B'/)
   u_AA= 0.0797d0
   u_AB= 0.0975d0
   vf= 2.1354*lattice_constant_graphene
   gamma_3=0d0
   gamma_4=0d0
   vppsigma= 0.3d0
   Nk=100
   Qcutoff = 4
   Num_bands= 0
   Electric_field= 0d0
   lattice_constant_c= 3.36d0 ! 3.36 Angstrom

   read(1001, PARAMETERS, iostat=stat)
   
   if (stat/=0) then
      backspace(1001)
      read(1001,fmt='(A)') inline
      write(*,'(A)') &
         '>>> ERROR : Invalid line in namelist PARAMETERS : '//trim(inline)
      stop
   endif


   !> twisted angle array
   allocate(twisted_angle_array(number_layers))
   twisted_angle_array=twisted_angle_array_input(1:number_layers)

   !> stacking sequences
   allocate(stacking_sequences(number_layers))
   stacking_sequences= stacking_sequences_input(1:number_layers)

   !> introduce one parameter to control the interlayer coupling
   allocate(interlayercoupling_ratio_array(number_layers-1))
   interlayercoupling_ratio_array= interlayercoupling_ratio_array_input(1:number_layers-1)


   if (twisted_index_m>0) then
      m= twisted_index_m
      twisted_angle_degree= acos((3d0*m*m + 3d0*m + 0.5d0)/(3d0*m*m + 3d0*m + 1d0));
      twisted_angle_degree= twisted_angle_degree*180d0/pi
   endif


   if (cpuid==0) then
      write(stdout, '(a, f6.2, a)')">> twisted angle: ", twisted_angle_degree, ' degree'
      write(stdout, '(a, 100i2)')">> twisted angle array: ", twisted_angle_array
      write(stdout, '(a, 100a3)')">> stacking_sequences: ", stacking_sequences
      write(stdout, '(a, f12.6, a)')">> vf= ", vf*1000d0, ' meV*Anstrom'
      write(stdout, '(a, f12.6, a)')">> u_AA= ", u_AA*1000d0, ' meV'
      write(stdout, '(a, f12.6, a)')">> u_AB= ", u_AB*1000d0, ' meV'
      write(stdout, '(a, f12.6, a)')">> gamma_3= ", gamma_3, ' eV'
      write(stdout, '(a, f12.6, a)')">> gamma_4= ", gamma_4, ' eV'
      write(stdout, '(a, f12.6, a)')">> vppsigma= ", vppsigma*1000d0, ' meV'
      write(stdout, '(a, i6)')">> Qcutoff= ", Qcutoff
      write(stdout, '(a, f12.6)')">> Electric_field= ", Electric_field
      write(stdout, '(a, f12.6)')">> lattice_constant_c= ", lattice_constant_c
      write(stdout, *) ' '
   endif

   v3=gamma_3*sqrt(3d0)/2d0*lattice_constant_graphene
   v4=gamma_4*sqrt(3d0)/2d0*lattice_constant_graphene

   !> setup moire unit cell parameters
   theta=twisted_angle_degree*pi/180d0
   rot1=reshape((/cos(theta/2d0),-sin(theta/2d0),  sin(theta/2d0), cos(theta/2d0)/), (/2, 2/))
   rot2=reshape((/cos(theta/2d0), sin(theta/2d0), -sin(theta/2d0), cos(theta/2d0)/), (/2, 2/))

   !> one edge of Moire BZ
   !     /\
   ! -> |  |
   !     \/
   qb= matmul(rot2, K_valley)- matmul(rot1, K_valley)
   kD= sqrt(sum(abs(qb)**2)) 

   if (cpuid==0) write( stdout, '(a, 2f10.6)')"qb= ", qb
   if (cpuid==0) write( stdout, '(a, 2f10.6)')"length of qb= ", kd
   if (cpuid==0) write( stdout, *)" "
  

   !> moire reciprocal lattice vectors
   b1m= (/-0.5d0, -sqrt3/2d0/)*sqrt3*kD
   b2m= (/1.0d0, 0.0d0/)*sqrt3*kD
   if (cpuid==0) write( stdout, '(a)')"Lattice vectors of moire supercell: "
   if (cpuid==0) write( stdout, '(2f12.6)')a1 
   if (cpuid==0) write( stdout, '(2f12.6)')a2 
   if (cpuid==0) write( stdout, '(a)')"Reciprocal lattice vectors of moire supercell: "
   if (cpuid==0) write( stdout, '(2f12.6)')b1m
   if (cpuid==0) write( stdout, '(2f12.6)')b2m

   !> moire valley of Graphene K valley. It sould be multiplied by -1 for graphene K' valley.
   K_valley1= (/sqrt3/2d0, -0.5d0/)*kD
   K_valley2= (/sqrt3/2d0,  0.5d0/)*kD
   if (cpuid==0) write( stdout, '(a)')" Morie valley for Graphene K valley"
   if (cpuid==0) write( stdout, '(2f12.6)')K_valley1
   if (cpuid==0) write( stdout, '(2f12.6)')K_valley2
   if (cpuid==0) write( stdout, '(a)')" Morie valley for Graphene K' valley"
   if (cpuid==0) write( stdout, '(2f12.6)')K_valley1*(-1d0)
   if (cpuid==0) write( stdout, '(2f12.6)')K_valley2*(-1d0)


   !> number of wave vectors
   !> Qvectors is in unit of b1m and b2m
   Num_Qvectors= (2*Qcutoff+1)**2
   allocate(Qvectors(2, Num_Qvectors))
   iter=0
   do i=1, 2*Qcutoff+1
      do j=1, 2*Qcutoff+1
         iter=iter+1
         Qvectors(1, iter)= i-Qcutoff-1
         Qvectors(2, iter)= j-Qcutoff-1
      enddo
   enddo


   !> dimension of Hamiltonian
   Ndim = 2*Num_Qvectors*number_layers
   if (Num_bands==0.or.Num_bands>Ndim) Num_bands=Ndim
   if (mod(Num_bands, 2).ne.0) Num_bands= Num_bands-1

   if (cpuid==0) write( stdout, '(a, i10)')"Number of Q vectors Num_Qvectors=: ", Num_Qvectors
   if (cpuid==0) write( stdout, '(a, i10)')"Dimension of the Hamiltonian matrix Ndim=: ", Ndim
   if (cpuid==0) write( stdout, '(a, i6, a)')"We are calculating ", Num_bands, ' bands close to E=0.'

   !> read kpath_bulk information
   rewind(1001)
   lfound = .false.
   do while (.true.)
      read(1001, *, end= 104)inline
      inline=upper(inline)
      if (trim(adjustl(inline))=='KPATH_BULK') then
         lfound= .true.
         if (cpuid==0) write(stdout, *)' '
         if (cpuid==0) write(stdout, *)'We found KPATH_BULK card'
         exit
      endif
   enddo

   !> kline for 3d band stOrigin_cell%Ructure
   !> high symmetry k points
   read(1001, *) nklines
   if(cpuid==0)write(stdout, '(a, 40i5)')'Number of K lines : ', nklines
   if(cpuid==0)write(stdout, '(a, 40i5)')'Number of Kpoints per line : ', Nk
   allocate(kline_start(2, nklines))
   allocate(kline_end(2, nklines))
   allocate(kline_name(nklines+1))
   allocate(kline_stop(nklines+1))
   kline_stop= 0d0
   kline_start= 0d0
   kline_end= 0d0
   kline_name= ' '
   it=0
   do i=1, nklines
      read(1001, *, err=201) kline_name(i), kline_start(:, i), &
         char_temp, kline_end(:, i)
      it= it+ 1

   enddo
201 continue
   kline_name(nklines+1)= char_temp

   if (it< nklines.and.cpuid==0) then
      write(stdout, *)' '
      write(stdout, *)' >>>> Error happens in KPATH_BULK card'
      write(stdout, *)' Error: the number of kpath lines should consistent'
      write(stdout, *)' Nklines ', nklines, ' the k lines you given ', it
      write(stdout, *)' Please set nklines to be  ', it
      stop
   endif


104 continue
   if (.not.lfound ) then
      if (cpuid==0) then
         write(stdout, *)' '
         write(stdout, *)">> You didn't set KPATH_BULK, we will use the default one"
         write(stdout, *)">> G-M-K-G-K'"
         write(stdout, *)' '
      endif
      nklines=4
      if(cpuid==0)write(stdout, '(a, 40i5)')'Number of K lines : ', nklines
      if(cpuid==0)write(stdout, '(a, 40i5)')'Number of Kpoints per line : ', Nk
      allocate(kline_start(2, nklines))
      allocate(kline_end(2, nklines))
      allocate(kline_name(nklines+1))
      allocate(kline_stop(nklines+1))
      kline_stop= 0d0
      kline_start= 0d0
      kline_end= 0d0
      kline_name= ' '
      kline_name(:)=(/'G ', 'M ', 'K ', 'G ', "K'"/)
      kline_start(:, 1)= (/0d0, 0d0/)
      kline_start(:, 2)= (/0d0, 0.5d0/)
      kline_start(:, 3)= (/0.33333d0, 0.33333d0/)
      kline_start(:, 4)= (/0d0, 0d0/)
      kline_end(:, 1)= (/0d0, 0.5d0/)
      kline_end(:, 2)= (/0.33333d0, 0.33333d0/)
      kline_end(:, 3)= (/0d0, 0d0/)
      kline_end(:, 4)= (/-0.33333d0, -0.33333d0/)
   endif

   if (cpuid==0) then
      write(stdout, *) ' '
      write(stdout, *) 'KPATH: '
      do i=1, nklines
         write(stdout, '(a5, 2f10.6, 3x, a5, 2f10.6)') kline_name(i), kline_start(:, i), kline_name(i+1), kline_end(:, i)
      enddo
      write(stdout, *) ' '
   endif
   NN= Nk
   nk_band= NN*nklines

   allocate(klen(nk_band))
   allocate(kpoints(2, nk_band))
   klen=0d0
   kpoints= 0d0
   t1= 0d0
   do j=1, nklines
      do i=1, NN
         kstart= kline_start(:, j)
         kend  = kline_end(:, j)
         k1= kstart(1)*b1m+ kstart(2)*b2m
         k2= kend(1)*b1m+ kend(2)*b2m

         kpoints(:, i+ (j-1)*NN)= kstart+ (kend- kstart)*dble(i-1)/dble(NN-1)

         temp= dsqrt((k2(1)- k1(1))**2 &
            +(k2(2)- k1(2))**2)/dble(NN-1)

         if (i.gt.1) then
            t1=t1+temp
         endif
         klen(i+(j-1)*NN)= t1
      enddo
      kline_stop(j+1)= t1
   enddo

   close(1001)
   outfileindex= 32432

   contains
function upper(s1) result (s2)
   character(*)       :: s1
   character(len(s1)) :: s2
   character          :: ch
   integer, parameter :: DUC = ICHAR('A') - ICHAR('a')
   integer            :: i

   do i = 1,LEN(s1)
      ch = s1(i:i)
      if (ch >= 'a'.AND.ch <= 'z') ch = CHAR(ICHAR(ch)+DUC)
      s2(i:i) = ch
   enddo
end function upper


end subroutine readinput



