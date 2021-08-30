
  subroutine get_hk(km, valley, hmnk)
     !> We generate hamiltonian for a given k in the Morie BZ
     !> km is in fractional unit of b1m and b2m
     ! Copied from WannierTools  https://github.com/quanshengwu/wannier_tools
     !> Author: Q.S Wu (wuquansheng@gmail.com)
     use para
     implicit none

     !> Graphene valley index, +1 for K, -1 for K'
     integer, intent(in) :: valley

     !> k points in the moire BZ
     real(dp), intent(in) :: km(2)

     !> hamiltonian
     complex(dp), intent(inout) :: Hmnk(Ndim, Ndim)

     logical :: is_tbg_layer(number_layers)
     logical :: is_bottom_layer(number_layers)

     integer :: iter, ilayer, i, iq, jq, q1(2), q2(2)
     integer :: idx1, idx2

     real(dp) :: qvec(2), phi, kvec(2)

     !> mono-layer graphene Hamiltonian
     !> H(k)=-vf*(k*valley*sx+k*sy)
     complex(dp) :: H_SLG(2, 2), kplus, kminus

     !> inter-layer hopping defined by the stacking_chirality without twist 
     complex(dp) :: H_interlayer(2, 2)

     !> three inter-layer hopping with twisting
     complex(dp) :: T1(2, 2), T2(2, 2), T3(2, 2)

     ! a function to determinate the stacking chirality based on their sequence
     integer, external :: stacking_chirality

     is_tbg_layer = .false.
     iter = 0
     is_bottom_layer = .true.
     do ilayer=1, number_layers-1
        if (twisted_angle_array(ilayer).ne. twisted_angle_array(ilayer+1)) then
           iter= iter+1
           is_tbg_layer(ilayer)=.true.
           is_tbg_layer(ilayer+1)=.true.
           is_bottom_layer(ilayer+1:number_layers)=.false.
        endif
     enddo

     if (iter>1) stop ">> ERROR: We only can deal with a system of only a single twist!"


     ! three tunnelling matrices under twisting
     phi=valley*2d0*pi/3d0
     T1=reshape((/u_AA, u_AB, u_AB, u_AA/), (/2, 2/))
     T2=reshape((/u_AA*z1, u_AB*exp( zi*phi), u_AB*exp(-zi*phi), u_AA*z1/), (/2, 2/))
     T3=reshape((/u_AA*z1, u_AB*exp(-zi*phi), u_AB*exp( zi*phi), u_AA*z1/), (/2, 2/))

     !> three kind of tunnelling matrices without twisting
     !> stacking chirality is 1: AB, BC, CA stacking
!    H_interlayer(:, :, 1)= reshape((/ 0d0, vppsigma, 0d0, 0d0/) , (/2, 2/))

     !> stacking chirality is 0: AA, BB, CC stacking
!    H_interlayer(:, :, 0)= reshape((/ vppsigma, 0d0, 0d0, vppsigma/) , (/2, 2/))

     !> stacking chirality is -1: BA, CB, AC stacking
!    H_interlayer(:, :,-1)= reshape((/ 0d0, 0d0, vppsigma, 0d0/) , (/2, 2/))

     !> diagonal term: intra-layer hopping
     do ilayer=1, number_layers
        do iq=1, Num_Qvectors
           
           qvec= Qvectors(:, iq)
           if (is_bottom_layer(ilayer)) then
              kvec= - K_valley1*valley+ (km(1)+qvec(1))*b1m+ (km(2)+qvec(2))*b2m
           else 
              kvec= - K_valley2*valley+ (km(1)+qvec(1))*b1m+ (km(2)+qvec(2))*b2m
           endif

           H_SLG=  vf*(kvec(1)*valley*sx+kvec(2)*sy)

           ! Electrical field 
           H_SLG= H_SLG- Electric_field*(ilayer-1d0-number_layers/2.0)*lattice_constant_c*s0
           

           idx1= (ilayer-1)*Num_Qvectors*2+ (iq-1)*2

           Hmnk(idx1+1:idx1+2, idx1+1:idx1+2)=H_SLG
        enddo
     enddo

     !> diagonal term: inter-layer hopping for non-TBG pair layers
     do ilayer=1, number_layers-1
        if (.not.(is_tbg_layer(ilayer).and.is_tbg_layer(ilayer+1))) then
           do iq=1, Num_Qvectors
              idx1= (ilayer-1)*Num_Qvectors*2+ (iq-1)*2
              idx2= ilayer*Num_Qvectors*2+ (iq-1)*2

              qvec= Qvectors(:, iq)
              if (is_bottom_layer(ilayer)) then
                 kvec= - K_valley1*valley+ (km(1)+qvec(1))*b1m+ (km(2)+qvec(2))*b2m
              else 
                 kvec= - K_valley2*valley+ (km(1)+qvec(1))*b1m+ (km(2)+qvec(2))*b2m
              endif
              kminus=-(kvec(1)*valley-zi*kvec(2))
              kplus=-(kvec(1)*valley+zi*kvec(2))
   
              i= stacking_chirality(stacking_sequences(ilayer), stacking_sequences(ilayer+1))

              !> AB, BC, CA
              if (i==1) then
                 H_interlayer(:, :)= reshape((/ v4*kminus, vppsigma*z1, v3*kplus, v4*kminus/) , (/2, 2/))
              !> BA, CB, AC
              elseif (i==-1) then
                 H_interlayer(:, :)= reshape((/ v4*kplus, v3*kminus, vppsigma*z1, v4*kplus/) , (/2, 2/))
              !> AA, BB, CC
              elseif (i==0) then
                 H_interlayer(:, :)= reshape((/ vppsigma, 0d0, 0d0, vppsigma/) , (/2, 2/))
              endif
              Hmnk(idx1+1:idx1+2, idx2+1:idx2+2)=H_interlayer(:, :)*interlayercoupling_ratio_array(ilayer)
              Hmnk(idx2+1:idx2+2, idx1+1:idx1+2)=transpose(conjg(H_interlayer(:, :)))*interlayercoupling_ratio_array(ilayer)
           enddo
        endif
     enddo

     !> non-diagonal term: inter-layer hopping for TBG layers
     do ilayer=1, number_layers-1
        if (is_tbg_layer(ilayer).and.is_tbg_layer(ilayer+1)) then
           do iq=1, Num_Qvectors
              q1=Qvectors(:, iq)
              do jq=1, Num_Qvectors
                 q2=Qvectors(:, jq)

                 idx1= (ilayer-1)*Num_Qvectors*2+ (iq-1)*2
                 idx2= ilayer*Num_Qvectors*2+ (jq-1)*2
                 if (iq==jq) then
                    Hmnk(idx1+1:idx1+2, idx2+1:idx2+2)=T1
                    Hmnk(idx2+1:idx2+2, idx1+1:idx1+2)=transpose(conjg(T1))
                 endif

                 if ((q2(1)-q1(1)==-valley).and.(q1(2)==q2(2))) then
                    Hmnk(idx1+1:idx1+2, idx2+1:idx2+2)=T2
                    Hmnk(idx2+1:idx2+2, idx1+1:idx1+2)=transpose(conjg(T2))
                 endif

                 if (q2(1)-q1(1)==-valley.and.q2(2)-q1(2)==-valley) then
                    Hmnk(idx1+1:idx1+2, idx2+1:idx2+2)=T3
                    Hmnk(idx2+1:idx2+2, idx1+1:idx1+2)=transpose(conjg(T3))
                 endif
              enddo ! jq
           enddo ! iq
        endif ! is_tbg_layer
     enddo ! ilayer

     return
  end subroutine get_hk


  function stacking_chirality(bottom, top)
     !> Determinate the stacking chirality with given stacking sequence.
     !> 1 for AB, BC, CA  
     !> 0 for AA, BB, CC  
     !> -1 for BA, CB, AC  
     !> Author: Q.S Wu (wuquansheng@gmail.com)

     implicit none

     integer :: stacking_chirality
     character(*), intent(inout) :: top, bottom
     character(2) :: stacking_sqeuence


     write(stacking_sqeuence, '(2A1)') bottom, top

     select case (stacking_sqeuence)
        case ('AB')
           stacking_chirality= 1
        case ('BC')
           stacking_chirality= 1
        case ('CA')
           stacking_chirality= 1

        case ('BA')
           stacking_chirality=-1
        case ('CB')
           stacking_chirality=-1
        case ('AC')
           stacking_chirality=-1

        case ('AA')
           stacking_chirality= 0
        case ('BB')
           stacking_chirality= 0
        case ('CC')
           stacking_chirality= 0
        case default
           stop 'ERROR in stacking_chirality: we only support AB BC CA BA CB AC AA BB CC stacking'
     end select

  end function stacking_chirality

