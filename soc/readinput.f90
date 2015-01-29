!
! this subroutine is used to read some paramters from
! input.dat
! constructed on 4/22/2010 by QS.Wu 


  subroutine readinput

     use para
     implicit none

     character*12 :: fname='input.dat'
     character*25 :: char_temp 
     logical ::  exists
     real(dp) :: cell_volume

     integer  :: i
     integer  :: j
     integer  :: NN
     real(dp) :: t1, temp
     real(dp) :: k1(3), k2(3)
     real(dp) :: kstart(3), kend(3)
    
     inquire(file=fname,exist=exists)
     if (exists)then
        write(*,*) 'read some paramters from input.dat'
        open(unit=1001,file=fname,status='old')
     else
        write(*,*)'file' ,fname, 'dosnot exist'
        stop
     endif
 

     read(1001,*)infilename
     write(*,'(a,a25)')' input file:',infilename
     read(1001,*)filename
     write(*,'(2a)')' output file:',filename
     read(1001,*)Nk
     write(*,*)'Nk',Nk
     read(1001,*)omeganum
     write(*,*)'omeganum',omeganum
     read(1001,*)omegamin, omegamax
     write(*,*)'omegamin, omegamax', omegamin, omegamax
     read(1001,*)nslab
     write(*,*)'nslab',nslab
     read(1001,*)Soc
     write(*,*)'soc',Soc
     read(1001,*)eta
     write(*,*)'eta',eta
     read(1001,*)E_fermi
     write(*,*)'E_fermi',E_fermi

     !> lattice information
     allocate(Rua(3), Rub(3), Ruc(3))
     allocate(Kua(3), Kub(3), Kuc(3))
     read(1001, *)Rua
     read(1001, *)Rub
     read(1001, *)Ruc

     !> transform lattice from direct space to reciprocal space

     Kua= 0d0
     Kub= 0d0
     Kuc= 0d0
     cell_volume= Rua(1)*(Rub(2)*Ruc(3)- Rub(3)*Ruc(2)) &
                 +Rua(2)*(Rub(3)*Ruc(1)- Rub(1)*Ruc(3)) &
                 +Rua(3)*(Rub(1)*Ruc(2)- Rub(2)*Ruc(1)) 
     cell_volume= 2d0*3.1415926535d0/cell_volume
     Kua(1)= cell_volume*(Rub(2)*Ruc(3)- Rub(3)*Ruc(2))
     Kua(2)= cell_volume*(Rub(3)*Ruc(1)- Rub(1)*Ruc(3))
     Kua(3)= cell_volume*(Rub(1)*Ruc(2)- Rub(2)*Ruc(1))

     Kub(1)= cell_volume*(Ruc(2)*Rua(3)- Ruc(3)*Rua(2))
     Kub(2)= cell_volume*(Ruc(3)*Rua(1)- Ruc(1)*Rua(3))
     Kub(3)= cell_volume*(Ruc(1)*Rua(2)- Ruc(2)*Rua(1))

     Kuc(1)= cell_volume*(Rua(2)*Rub(3)- Rua(3)*Rub(2))
     Kuc(2)= cell_volume*(Rua(3)*Rub(1)- Rua(1)*Rub(3))
     Kuc(3)= cell_volume*(Rua(1)*Rub(2)- Rua(2)*Rub(1))

     write(*, '(a)') '>> lattice information'
     write(*, '(3f10.6)')Rua
     write(*, '(3f10.6)')Rub
     write(*, '(3f10.6)')Ruc

     write(*, '(a)') '>> Reciprocal lattice information'
     write(*, '(3f10.6)')Kua
     write(*, '(3f10.6)')Kub
     write(*, '(3f10.6)')Kuc

     !> kline for 3d band structure
     !> high symmetry k points
     read(1001, *) nk3lines
     allocate(k3line_start(3, nk3lines))
     allocate(k3line_end(3, nk3lines))
     do i=1, nk3lines
        read(1001, *) char_temp, k3line_start(:, i), &
                      char_temp, k3line_end(:, i)
     enddo
     close(1001)

     NN= 20
     nk3_band= NN*nk3lines
     allocate(k3len(nk3_band))
     allocate(k3points(3, nk3_band))
     k3len=0d0
     k3points= 0d0
     t1= 0d0
     do j=1, nk3lines
        do i=1, NN
           kstart= k3line_start(:, j)
           kend  = k3line_end(:, j)
           k1= kstart(1)*Kua+ kstart(2)*Kub+ kstart(3)*Kuc
           k2= kend(1)*Kua+ kend(2)*Kub+ kend(3)*Kuc

           k3points(:, i+ (j-1)*NN)= kstart+ (kend- kstart)*dble(i-1)/dble(NN-1)
           
           temp= dsqrt((k2(1)- k1(1))**2 &
                 +(k2(2)- k1(2))**2  &
                 +(k2(3)- k1(3))**2)/dble(NN-1) 

           if (i.gt.1) then
              t1=t1+temp
           endif
           k3len(i+(j-1)*NN)= t1
        enddo
     enddo

	  omegamin=omegamin
	  omegamax=omegamax
 
     eta=(omegamax- omegamin)/omeganum*eta
	  print * ,'eta', eta
     close(1001)

     write(*,*)'read input.dat file successfully'

     return
  end subroutine



