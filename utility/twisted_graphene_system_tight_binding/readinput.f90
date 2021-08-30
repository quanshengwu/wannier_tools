
subroutine readinput
   ! Read in the control paramters from input.dat,
   ! and set default values if not specified in the input.dat
   !
   ! Constructed on 4/22/2010 by QS.Wu

   use para
   implicit none

   integer :: stat
   character*12 :: fname='system.in'
   character*25 :: char_temp
   character*256 :: inline
   logical ::  exists, lfound

   ! an example for Twisted A-AB system
   write(stdout, *) ' '
   write(stdout, *) '&PARAMETERS'
   write(stdout, *) 'number_layers = 3    ! Two layers for TBG'
   write(stdout, *) 'twisted_index_m= 1   ! twisted index m, theta= acos((3d0*m*m + 3d0*m + 0.5d0)/(3d0*m*m + 3d0*m + 1d0));'
   write(stdout, *) 'twisted_angle_array_input = 0 1 1  ! twisted angle array, unit is theta; Number_of_layers numbers'
   write(stdout, *) 'stacking_sequences_input = "A" "A"  "B" ! AB stacking sequences', &
      ' only three values "A", "B", "C"; Number_of_layers numbers'
   write(stdout, *) 'use_poscar = F    ! use POSCAR or not, =T we will generate POSCAR '
   write(stdout, *) 'hr_generate = F    ! use POSCAR or not, =T we will generate hr.dat'
   write(stdout, *) 'gen_sparse_hr = F   ! use POSCAR or not, =T we will generate hr.dat in sparse format'
   write(stdout, *) 'hr_cutoff=0.10000  ! set HmnR=0 if HmnR<hr_cutoff'
   write(stdout, *) 'vpppi=-2.81       ! pi bond of p orbital '
   write(stdout, *) 'iR_cut = 1   ! R is in [-iR_cut, -iR_cut+1, ..., iR_cut]'
   write(stdout, *) '/'
   write(stdout, *) ' '


   inquire(file=fname,exist=exists)
   if (exists)then
      write(stdout,*) '  '
      write(stdout,*) '>>>Read some paramters from system.in'
      open(unit=1001,file=fname,status='old')
   else
      write(stdout,*)'file' ,fname, 'dosnot exist'
      stop
   endif

   !> initialization
   number_layers= 2
   twisted_index_m= 10
   vpppi=-2.81d0
   hr_cutoff=eps6
   iR_cut=1

   read(1001, PARAMETERS, iostat=stat)
   
   if (stat/=0) then
      backspace(1001)
      read(1001,fmt='(A)') inline
      write(stdout,'(A)') &
         '>>> ERROR : Invalid line in namelist PARAMETERS : '//trim(inline)
      stop
   endif


   twisted_angle_degree= acos((3d0*twisted_index_m*twisted_index_m + 3d0*twisted_index_m +&
   0.5d0)/(3d0*twisted_index_m*twisted_index_m + 3d0*twisted_index_m + 1d0))*180/pi

   !> twisted angle array
   allocate(twisted_angle_array(number_layers))
   twisted_angle_array=twisted_angle_array_input(1:number_layers)
   write(stdout, '(a, 100i2)')">> twisted angle array: ", twisted_angle_array

   !> stacking sequences
   allocate(stacking_sequences(number_layers))
   stacking_sequences= stacking_sequences_input(1:number_layers)
   write(stdout, '(a, 100a3)')">> stacking_sequences: ", stacking_sequences


   !> a logical number: whether use the existed POSCAR or not
   if (use_poscar) then
      write(stdout, *)">> We are using the user-provided POSCAR"
   else
      write(stdout, *)">> We will generate POSCAR"
   endif

   !> 
   write(stdout, '(a, f12.6)')">> hr_cutoff= ", hr_cutoff
   write(stdout, '(a, f12.6)')">> vpppi= ", vpppi
   write(stdout, '(a, i6)')">> iR_cut= ", iR_cut


   if (use_poscar) then
      write(stdout, *)">> We are going to generate tight-binding model"
   else
      write(stdout, *)">> We don't generate tight-binding model"
   endif

   close(1001)

   outfileindex= 32432
   return
end subroutine readinput


