!
! License-Identifier: GPL
!
! Copyright (C) 2019 The Yambo Team
!
! Authors (see AUTHORS file for details): DS
!
subroutine Dipole_spin_magnetization(Xen,Xk,dip)
 !
 ! This routine returns <sigma>
 !
 use pars,                ONLY:cI
 use drivers,             ONLY:l_sc_run
 use interfaces,          ONLY:WF_load,WF_free
 use LIVE_t,              ONLY:live_timing
 use timing_m,            ONLY:timing
 use electrons,           ONLY:levels,n_spinor,n_sp_pol,n_spin
 use R_lattice,           ONLY:bz_samp
 use DIPOLES,             ONLY:DIP_spin,DIPOLE_t
 use fft_m,               ONLY:fft_size
 use wave_func,           ONLY:WF
 use wrapper,             ONLY:Vstar_dot_V
 use parallel_m,          ONLY:PAR_IND_DIPk_ibz,PAR_DIPk_ibz_index,PAR_IND_DIPk_ibz_ID,&
&                              PAR_IND_VAL_BANDS_DIP,PAR_IND_VAL_BANDS_DIP_ID,&
&                              PAR_IND_CON_BANDS_DIP,PAR_IND_CON_BANDS_DIP_ID
 !
 implicit none
 !
 type(bz_samp), intent(in)    :: Xk
 type(levels),  intent(in)    :: Xen
 type(DIPOLE_t),intent(inout) :: Dip
 !
 ! Work Space
 !
 integer                      :: ik,ic_wf(n_sp_pol),iv_wf(n_sp_pol),ic,iv,n_LT_steps,ik_mem
 logical                      :: LT_msg=.FALSE.
 !
 if (n_spin==1) return
 !
 call timing("DIPOLE_spin_M",OPR="start")
 !
 ! Parallelization Mask 
 !======================
 ! Depending on the environment DIPOLES are distributed in different ways
 !
 n_LT_steps=PAR_IND_DIPk_ibz%n_of_elements(PAR_IND_DIPk_ibz_ID+1)*&
&           PAR_IND_VAL_BANDS_DIP%n_of_elements(PAR_IND_VAL_BANDS_DIP_ID+1)*&
&           PAR_IND_CON_BANDS_DIP%n_of_elements(PAR_IND_CON_BANDS_DIP_ID+1)
 !
 ! Main loop over k in IBZ
 !
 do ik=1,Xk%nibz
   !
   if (n_LT_steps==0) cycle
   ! 
   if (.not.PAR_IND_DIPk_ibz%element_1D(ik)) cycle 
   !
   ik_mem=PAR_DIPk_ibz_index(ik)
   !
   call WF_load(WF,0,1,Dip%ib,(/ik,ik/),space='R',title='-Spin oscillators/R space',&
&               force_WFo=l_sc_run,keep_states_to_load=.TRUE.)
   !
   if (.not.LT_msg) then
     if (n_LT_steps>0) call live_timing('Spin magnetiziation',n_LT_steps)
     LT_msg=.TRUE.
   endif
   !
   do iv=Dip%ib(1),Dip%ib_lim(1) 
     !
     if(.not.PAR_IND_VAL_BANDS_DIP%element_1D(iv)) cycle
     !
     do ic=max(iv,Dip%ib_lim(2)),Dip%ib(2)
       !
       if(.not.PAR_IND_CON_BANDS_DIP%element_1D(ic)) cycle
       !
       iv_wf=WF%index(iv,ik,:)
       ic_wf=WF%index(ic,ik,:)
       !
       if (any((/iv_wf==0/))) call error(" Error in parallel wf distribution (Dipole iv)")
       if (any((/ic_wf==0/))) call error(" Error in parallel wf distribution (Dipole ic)")
       !
       ! Evaluate <iv_wf| sigma_i |ic_wf>
       !=================================
       !
       if (n_sp_pol==2) then
         !
         ! mx
         !
         DIP_spin(1,ic,iv,ik_mem,2,1)=    Vstar_dot_V(fft_size,WF%c(:,1,iv_wf(1)),WF%c(:,1,ic_wf(2)))
         DIP_spin(1,ic,iv,ik_mem,1,1)=    Vstar_dot_V(fft_size,WF%c(:,1,iv_wf(2)),WF%c(:,1,ic_wf(1)))
         !
         ! my
         !
         DIP_spin(2,ic,iv,ik_mem,2,1)=-cI*Vstar_dot_V(fft_size,WF%c(:,1,iv_wf(1)),WF%c(:,1,ic_wf(2)))
         DIP_spin(2,ic,iv,ik_mem,1,1)=+cI*Vstar_dot_V(fft_size,WF%c(:,1,iv_wf(2)),WF%c(:,1,ic_wf(1)))
         !
         ! mz
         !
         DIP_spin(3,ic,iv,ik_mem,1,1)=    Vstar_dot_V(fft_size,WF%c(:,1,iv_wf(1)),WF%c(:,1,ic_wf(1)))
         DIP_spin(3,ic,iv,ik_mem,2,1)=   -Vstar_dot_V(fft_size,WF%c(:,1,iv_wf(2)),WF%c(:,1,ic_wf(2)))
         !
       elseif (n_spinor==2) then
         !
         ! mx
         !
         DIP_spin(1,ic,iv,ik_mem,1,1)=    Vstar_dot_V(fft_size,WF%c(:,1,iv_wf(1)),WF%c(:,2,ic_wf(1))) &
         &                             +Vstar_dot_V(fft_size,WF%c(:,2,iv_wf(1)),WF%c(:,1,ic_wf(1)))
         !
         ! my
         !
         DIP_spin(2,ic,iv,ik_mem,1,1)=-cI*Vstar_dot_V(fft_size,WF%c(:,1,iv_wf(1)),WF%c(:,2,ic_wf(1))) &
         &                          +cI*Vstar_dot_V(fft_size,WF%c(:,2,iv_wf(1)),WF%c(:,1,ic_wf(1)))
         !
         ! mz
         !
         DIP_spin(3,ic,iv,ik_mem,1,1)=    Vstar_dot_V(fft_size,WF%c(:,1,iv_wf(1)),WF%c(:,1,ic_wf(1))) &
         &                             -Vstar_dot_V(fft_size,WF%c(:,2,iv_wf(1)),WF%c(:,2,ic_wf(1)))
         !
       endif
       !
       call live_timing(steps=1)
       !
     enddo    ! conduction band loop
   enddo      ! valence    band loop
   !
   call WF_free(WF)
   ! 
 enddo        ! k-points loop
 !
 call live_timing()
 !
 call timing("DIPOLE_spin_M",OPR="stop")
 !
end subroutine Dipole_spin_magnetization
