!
! License-Identifier: GPL
!
! Copyright (C) 2015 The Yambo Team
!
! Authors (see AUTHORS file for details): DS
!
subroutine DIPOLE_x_real_space(Xen,Xk,Dip,idir_todo)
 !
 use drivers,             ONLY:l_sc_run
 use pars,                ONLY:SP,cI,cZERO
 use D_lattice,           ONLY:x_cc,molecule_position,a
 use wave_func,           ONLY:WF
 use interfaces,          ONLY:WF_load,WF_free
 use wrapper_omp,         ONLY:Vstar_dot_V_omp
 use FFT_m,               ONLY:fft_size
 use R_lattice,           ONLY:bz_samp
 use electrons,           ONLY:levels,n_sp_pol,n_spinor
 use DIPOLES,             ONLY:DIPOLE_t,DIP_iR,DIP_v,x_space_obs
 use X_m,                 ONLY:l_X_terminator
 use LIVE_t,              ONLY:live_timing
 use timing_m,            ONLY:timing
 use parallel_m,          ONLY:PAR_IND_DIPk_ibz,PAR_DIPk_ibz_index,PAR_IND_DIPk_ibz_ID,&
&                              PAR_IND_VAL_BANDS_DIP,PAR_IND_VAL_BANDS_DIP_ID,&
&                              PAR_IND_CON_BANDS_DIP,PAR_IND_CON_BANDS_DIP_ID
 !
#include<y_memory.h>
 !
 type(bz_samp), intent(in)    :: Xk
 type(levels),  intent(in)    :: Xen
 type(DIPOLE_t),intent(inout) :: Dip
 logical,       intent(in)    :: idir_todo(3)
 !
 ! Work space
 !
 integer                      :: ik,ik_mem,i_sp_pol,i_spinor,ic,iv,ic_min,&
&                                ic_wf,iv_wf,idir,n_LT_steps
 real(SP)                     :: Ev_m_Ec_KS,x_center(3)
 complex(SP)                  :: DIP_tmp(3,Dip%ib_lim(2):Dip%ib(2),Dip%ib(1):Dip%ib_lim(1))
 logical                      :: LT_msg=.FALSE.
 !
 call timing('DIPOLE_x_space',OPR='start')
 !
 if(.not.Dip%force_v_g_space) Dip%Vnl_included=.true.
 ! 
 call fft_setup(0,1,.true.)      ! FFT size
 !
 YAMBO_ALLOC(x_cc,(fft_size,3,1))
 !
 call load_cc(.false.,molecule_position,1)
 !
 n_LT_steps=n_sp_pol*PAR_IND_DIPk_ibz%n_of_elements(PAR_IND_DIPk_ibz_ID+1)*&
&                    PAR_IND_VAL_BANDS_DIP%n_of_elements(PAR_IND_VAL_BANDS_DIP_ID+1)*&
&                    PAR_IND_CON_BANDS_DIP%n_of_elements(PAR_IND_CON_BANDS_DIP_ID+1)
 !
 do i_sp_pol=1,n_sp_pol
   !
   do ik=1,Xk%nibz 
     !
     if (.not.PAR_IND_DIPk_ibz%element_1D(ik)) cycle
     !
     call WF_load(WF,0,1,(/1,Dip%ib(2)/),(/ik,ik/),sp_pol_to_load=(/i_sp_pol,i_sp_pol/),space='R',&
&             title='-Oscillators/R space/Real Space',force_WFo=l_sc_run,keep_states_to_load=.TRUE.)
     !
     ik_mem=PAR_DIPk_ibz_index(ik)
     !
     DIP_tmp=cZERO
     !
     if (n_LT_steps>0.and..not.LT_msg) call live_timing(trim(x_space_obs)//" [x-space]",n_LT_steps)
     LT_msg=.TRUE.
     !
     do iv=Dip%ib(1),Dip%ib_lim(1) 
       !
       if(.not.PAR_IND_VAL_BANDS_DIP%element_1D(iv)) cycle
       !
       ic_min=max(iv,Dip%ib_lim(2))
       if (l_X_terminator) ic_min=iv
       !
       do ic=ic_min,Dip%ib(2)
         !
         if(.not.PAR_IND_CON_BANDS_DIP%element_1D(ic)) cycle
         !
         iv_wf=WF%index(iv,ik,i_sp_pol)
         ic_wf=WF%index(ic,ik,i_sp_pol)
         !
         if (iv_wf==0) call error(" Error in parallel wf distribution (Dipole iv)")
         if (ic_wf==0) call error(" Error in parallel wf distribution (Dipole ic)")
         !
         do i_spinor=1,n_spinor
           do idir=1,3
             !
             ! Here I assume that for the non periodic direction id_red==id_cartesian
             !
             if(.not.idir_todo(idir)) cycle
             !
             DIP_tmp(idir,ic,iv) = DIP_tmp(idir,ic,iv) + &
&                       Vstar_dot_V_omp(fft_size,WF%c(1:fft_size,i_spinor,iv_wf), &
&                             x_cc(1:fft_size,idir,1)*WF%c(1:fft_size,i_spinor,ic_wf) )
             !
           enddo
         enddo
         !
         call live_timing(steps=1)
         !
       enddo 
       !
     enddo
     !
     !======
     ! <iR>
     !======
     do idir=1,3
       if(.not.idir_todo(idir)) cycle
       DIP_iR(idir,:,:,ik_mem,i_sp_pol)=cI*DIP_tmp(idir,:,:)
     enddo
     !
     !=====
     ! <v>
     !=====
     if(.not.Dip%force_v_g_space) then
       do iv=Dip%ib(1),Dip%ib_lim(1) 
         if(.not.PAR_IND_VAL_BANDS_DIP%element_1D(iv)) cycle
         do ic=max(iv,Dip%ib_lim(2)),Dip%ib(2)
           if(.not.PAR_IND_CON_BANDS_DIP%element_1D(ic)) cycle
           Ev_m_Ec_KS=Xen%E(iv,ik,i_sp_pol)-Xen%E(ic,ik,i_sp_pol)
           if (allocated(Xen%Eo)) Ev_m_Ec_KS=Xen%Eo(iv,ik,i_sp_pol)-Xen%Eo(ic,ik,i_sp_pol)
           do idir=1,3
             if(.not.idir_todo(idir)) cycle
             DIP_v(idir,ic,iv,ik_mem,i_sp_pol)=cI*DIP_tmp(idir,ic,iv)*Ev_m_Ec_KS
           enddo
         enddo
       enddo
     endif
     !
     call WF_free(WF)
     !
   enddo 
   !
 enddo
 !
 call live_timing()
 !
 call timing('DIPOLE_x_space',OPR='stop')
 !
 YAMBO_FREE(x_cc)
 !
end subroutine Dipole_x_real_space
