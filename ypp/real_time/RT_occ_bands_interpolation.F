!
!        Copyright (C) 2000-2016 the YAMBO team
!              http://www.yambo-code.org
!
! Authors (see AUTHORS file for details): DS
! 
! This file is distributed under the terms of the GNU 
! General Public License. You can redistribute it and/or 
! modify it under the terms of the GNU General Public 
! License as published by the Free Software Foundation; 
! either version 2, or (at your option) any later version.
!
! This program is distributed in the hope that it will 
! be useful, but WITHOUT ANY WARRANTY; without even the 
! implied warranty of MERCHANTABILITY or FITNESS FOR A 
! PARTICULAR PURPOSE.  See the GNU General Public License 
! for more details.
!
! You should have received a copy of the GNU General Public 
! License along with this program; if not, write to the Free 
! Software Foundation, Inc., 59 Temple Place - Suite 330,Boston, 
! MA 02111-1307, USA or visit http://www.gnu.org/copyleft/gpl.txt.
!
subroutine RT_occ_bands_interpolation(En,kpt,n_T_steps)
 !
 use pars,           ONLY:SP,schlen
 use units,          ONLY:HA2EV,FS2AUT
 use electrons,      ONLY:levels,n_sp_pol
 use R_lattice,      ONLY:bz_samp
 use com,            ONLY:msg,of_open_close,error
 use YPP,            ONLY:BANDS_steps,Nel_fac,coo_in,coo_out,k_transform,interp_grid,&
&                         l_RT_life_bands,l_RT_occ_bands,electrons_bands_interpolate
 use stderr,         ONLY:real2ch,intc
 use SC,             ONLY:SC_bands
 use parser_m,        ONLY:parser
 use rt_ctl,         ONLY:OCC_T_range,OCC_deltaT,RT_apply,RT_carriers_object
 use real_time,      ONLY:RT_carriers
 use interpolate,    ONLY:eval_interpolation_coeff,bz_interpolation,reset_interpolation
 !
 implicit none
 !
 integer,       intent(in) :: n_T_steps
 type(bz_samp), intent(in) :: kpt
 type(levels),  intent(in) :: En
 !
 ! Work space
 !
 type(bz_samp)         :: k_bands
 type(levels)          :: E_bands
 !
 type(bz_samp)         :: USER_K,CIRCUIT_K,INTERP_K
 real(SP), pointer     :: bands_interpolated(:,:,:) => null()
 !
 real(SP), allocatable :: fit_on_circuit(:,:,:,:),variable_to_fit(:,:,:)
 real(SP)              :: io_Time,Nel,Nhole,Real_Nhole,Real_Nel,Max_occ,values(3)
 character(3)          :: message
 character(schlen)     :: headings(3),file_name_rt
 logical               :: l_normalize_Nel,l_g_space_interp,l_nearest_neig_interp
 integer               :: ib,ik,i_sp_pol,i_RT,i_T,ID_interp
 !
 !
 call section('*','Interpolation tool')
 !=====================================
 !
 if(all(interp_grid>0).or.BANDS_steps<=0) &
&  call error(" Interpolation of occupations/lifetimes only on a band circuit")
 !
 if(n_sp_pol>1) &
&  call error('Occupations for spin polarized systems not implemented')
 !
 ! Map RT_carriers%kpt ==> k_bands%kpt
 !     RT_carriers%E   ==> E_bands%E 
 !======================================
 call RT_occ_bands_kpts_and_En_map(En,kpt,E_bands,k_bands)
 !
 ! Perform the energies interpolation
 !======================================
 call electrons_bands_interpolate(E_bands,K_bands,ID_interp,USER_k,INTERP_k,CIRCUIT_k,SC_bands,bands_interpolated)
 !
 ! Perform the occupations interpolation
 !======================================
 Max_occ=0._SP
 !
 ! Define the kind of interpolation
 call parser('NNInterp',l_nearest_neig_interp)
 l_g_space_interp=.not.l_nearest_neig_interp
 !
 ! Check if Nel has to be renormalized
 call parser('NormN',l_normalize_Nel)
 !
 allocate(variable_to_fit(SC_bands(1):SC_bands(2),1:n_sp_pol,RT_carriers%nk))
 !
 allocate(fit_on_circuit(SC_bands(1):SC_bands(2),1:n_sp_pol,CIRCUIT_k%nibz,n_T_steps))
 !
 do i_T=1,n_T_steps
   !
   io_Time=min(OCC_T_range(1)+(i_T-1)*OCC_deltaT,OCC_T_range(2))
   !
   call RT_apply(SC_bands,En,kpt,Time=io_Time,keep_RT_carriers=.TRUE.)
   !
   Real_Nhole=0._SP
   Real_Nel=0._SP
   do i_RT=1,RT_carriers%nk
     ib      =RT_carriers%table(i_RT,1)
     !ik_small_grid=RT%table(i_RT,2)
     ik      =RT_carriers%table(i_RT,3)
     i_sp_pol=RT_carriers%table(i_RT,4)
     !
     if(ib<=E_bands%nbf) Real_Nhole=Real_Nhole-RT_carriers%delta_f(i_RT)*k_bands%weights(ik)
     if(ib> E_bands%nbf) Real_Nel  =Real_Nel  +RT_carriers%delta_f(i_RT)*k_bands%weights(ik)
     !
   enddo
   Max_occ=maxval((/Max_occ,RT_carriers%delta_f/))
   !
   if(l_g_space_interp) then
     !
     ! STEP 1: define coefficients for occupations/lifetimes
     !
     do i_RT=1,RT_carriers%nk
       ib      =RT_carriers%table(i_RT,1)
       !ik_small_grid=RT%table(i_RT,2)
       ik      =RT_carriers%table(i_RT,3)
       i_sp_pol=RT_carriers%table(i_RT,4)
       !
       if(l_RT_occ_Bands ) variable_to_fit(ib,i_sp_pol,ik)=    RT_carriers%delta_f(i_RT)
       if(l_RT_life_Bands) variable_to_fit(ib,i_sp_pol,ik)=sum(RT_carriers%GAMMA_bare(i_RT,:))
       !
     enddo
     !
     call eval_interpolation_coeff(R2D=variable_to_fit,k=k_bands,Nk=k_bands%nibz,ID=ID_interp)
     !
     ! STEP 2
     ! 1: interpolate on the circuit      if (BANDS_steps>0)
     !
     if (BANDS_steps> 0) call bz_interpolation(USER_k=CIRCUIT_k,R2D=fit_on_circuit(:,:,:,i_T),ID=ID_interp)
     !
     call reset_interpolation(ID_interp)
     !
   elseif(l_nearest_neig_interp) then
     !
     message='RT '
     if(i_T==1) message='RT0'
     !
     if(l_RT_occ_Bands)  then
       call Nearest_kpt_interpolation(RT_carriers%nk,CIRCUIT_k%nibz,RT_carriers%nb,&
&                                       RT_carriers%nstates,1,1,RT_carriers%table,RT_carriers%kpt, &
&                                       CIRCUIT_k%pt,RT_carriers%k_weight,CIRCUIT_k%weights,&
&                                       RT_carriers%delta_f,fit_on_circuit(:,:,:,i_T),message,l_normalize_Nel)
     endif
     if(l_RT_life_Bands) then
       call Nearest_kpt_interpolation(RT_carriers%nk,CIRCUIT_k%nibz,&
&                                       RT_carriers%nb,RT_carriers%nstates,1,1,RT_carriers%table,RT_carriers%kpt, &
&                                       CIRCUIT_k%pt,RT_carriers%k_weight,CIRCUIT_k%weights,&
&                                       RT_carriers%GAMMA_bare,fit_on_circuit(:,:,:,i_T),message,.false.)
     endif
     !
   endif
   !
   ! Normilize the occupations
   !
   if(l_RT_occ_Bands) then
     Nhole=-sum(fit_on_circuit(SC_bands(1):E_bands%nbf,:,:,i_T))/real(CIRCUIT_k%nibz)
     Nel=   sum(fit_on_circuit(E_bands%nbf+1:SC_bands(2),:,:,i_T))/real(CIRCUIT_k%nibz)
     !
     if(.not.l_nearest_neig_interp) then
       if(Nhole>0.and.l_normalize_Nel) &
&          fit_on_circuit(SC_bands(1):E_bands%nbf,:,:,i_T)=&
&            fit_on_circuit(SC_bands(1):E_bands%nbf,:,:,i_T)/Nhole*Real_Nhole
       if(Nel>0  .and.l_normalize_Nel) &
&          fit_on_circuit(E_bands%nbf+1:SC_bands(2),:,:,i_T)=&
&            fit_on_circuit(E_bands%nbf+1:SC_bands(2),:,:,i_T)/Nel*Real_Nel
     endif
     !
   endif
   !
   if(l_g_space_interp) call reset_interpolation(ID_interp)
   !
 enddo 
 !
 call RT_carriers_object(RT_carriers,WHAT='free')
 !
 deallocate(variable_to_fit)
 !
 ! Output: bands on circuit
 !==========================
 !
 headings(1)="Step"
 if(l_RT_occ_Bands)  headings(3)=' occ'
 if(l_RT_life_Bands) headings(3)='LIFE'
 !
 coo_out=coo_in
 !
 do i_T=1,n_T_steps
   !
   io_Time=min(OCC_T_range(1)+(i_T-1)*OCC_deltaT,OCC_T_range(2))
   !
   if (l_RT_occ_Bands ) file_name_rt="YPP-RT_occ_bands_iT"//intc(i_T)
   if (l_RT_life_Bands) file_name_rt="YPP-RT_REF_LIFE_bands"
   !
   call of_open_close(file_name_rt,'ot')
   !
   call msg('o '//file_name_rt," "," ",INDENT=0)
   call msg('o '//file_name_rt,"#",' t= '//trim(real2ch( io_Time/FS2AUT )),INDENT=0)
   !
   do ib=SC_bands(1),SC_bands(2)
     !
     call msg('o '//file_name_rt,"#"," ",INDENT=0)
     !
     headings(2)=' b'//intc(ib) 
     !
     call msg('o '//file_name_rt,"#",headings,INDENT=0,USE_TABS=.true.)
     call msg('o '//file_name_rt,"#"," ",INDENT=0)
     !
     do ik=1,CIRCUIT_K%nibz
       !
       values(1)=real(ik,SP)
       values(2)=bands_interpolated(ib,1,ik)*HA2EV
       if(l_RT_occ_Bands)  values(3)=fit_on_circuit(ib,1,ik,i_T)/Max_occ*Nel_fac
       if(l_RT_life_Bands) values(3)=fit_on_circuit(ib,1,ik,i_T)*HA2EV*1000._SP
       call msg('o '//file_name_rt,' ',values,INDENT=0,USE_TABS=.true.)
       !
     enddo
     !
   enddo
   !
   call of_open_close(file_name_rt)
   !
 enddo
 !
 deallocate(fit_on_circuit)
 !
end subroutine RT_occ_bands_interpolation