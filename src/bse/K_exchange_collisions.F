!
! License-Identifier: GPL
!
! Copyright (C) 2015 The Yambo Team
!
! Authors (see AUTHORS file for details): AM
!
subroutine K_exchange_collisions(iq,Xk,i_T_grp,NG,WF_k,WF_p,i_k_ref,i_p_ref,ipath,l_bs_exch_wf_in_loop)
 !
 use pars,           ONLY:SP,cZERO
 use D_lattice,      ONLY:nsym,i_time_rev,sop_inv,sop_inv_d
 use BS,             ONLY:BS_T_grp,K_EXCH_collision,BS_n_g_exch,&
&                         BSE_L_kind,l_BS_magnons,l_BSE_minimize_memory,BS_bands
 use wave_func,      ONLY:WAVEs,WF
 use R_lattice,      ONLY:g_rot,g_rot_d,qindx_X,bz_samp
 use stderr,         ONLY:STRING_match
 use electrons,      ONLY:n_sp_pol
 use interfaces,     ONLY:scatter_Bamp
 use timing_m,       ONLY:timing
 use openmp,         ONLY:OPENMP_update,master_thread
 use collision_el,   ONLY:elemental_collision_free,elemental_collision_alloc
 use devxlib,        ONLY:devxlib_memcpy_d2d,devxlib_memcpy_h2d,devxlib_memcpy_d2h,&
&                         devxlib_conjg_d,devxlib_memset_d
 !
#include<dev_defs.h>
 !
 implicit none
 !
 type(bz_samp),intent(in) :: Xk
 integer,      intent(in) :: iq,i_T_grp,NG(2)
 !
 ! New implementation. WFs are passed explicitly
 !
 integer,              intent(inout) :: i_k_ref,i_p_ref
 type (WAVEs), target, intent(inout) :: WF_k,WF_p
 integer,              intent(in)    :: ipath
 logical,              intent(in)    :: l_bs_exch_wf_in_loop
 !
 ! Work Space
 !
 logical :: l_load_WFs_k,l_load_WFs_p
 integer :: i_T_el,i_T_el_p,N_T_el_p,i_c,i_v,i_sp_c,i_sp_v,i_k_bz,i_k,i_s,i_g0,i_g1,i_g2,i_p_bz,i_g_p,i_p,&
&           i_sp,i_T_el_last_with_identity_sym,i_T_grp_last_with_identity_sym,i_T_grp_p
 complex(SP), pointer DEV_ATTR :: O_x_p(:,:),O_x_sym_p(:,:)
 !
 call timing('T_space EXX Osc.',OPR='start')
 !
 ! Allocate
 !
 call elemental_collision_free(K_EXCH_collision)
 !
#if defined(__NOTNOW)
 !DEV_OMP parallel default(shared), private( K_EXCH_collision, &
 !DEV_OMP &        i_T_el,i_k_bz,i_k,i_s, i_p_bz,i_p,i_sp, i_v,i_c,i_sp_c,i_sp_v, &
 !DEV_OMP &        i_T_grp_last_with_identity_sym, i_T_el_last_with_identity_sym, &
 !DEV_OMP &        i_T_grp_p, N_T_el_p, i_T_el_p, i_g1, i_g2)
 !DEV_OMP &        firstprivate(i_k_ref, i_p_ref), lastprivate(i_k_ref,i_p_ref)
#endif
 !
 call OPENMP_update(master_thread)
 !
 call elemental_collision_alloc(K_EXCH_collision,NG=BS_n_g_exch,TITLE="EXCH_COLL")
 !
#if defined(__NOTNOW)
 !DEV_OMP do schedule(dynamic)
#endif
 T_loop: do i_T_el=1,BS_T_grp(i_T_grp)%size
   !
   i_k_bz=BS_T_grp(i_T_grp)%table(i_T_el,1)
   i_k   =Xk%sstar(i_k_bz,1)
   i_s   =Xk%sstar(i_k_bz,2)
   !
   i_p_bz=qindx_X(iq,i_k_bz,1)
   i_g_p =qindx_X(iq,i_k_bz,2)
   !
   i_p   =Xk%sstar(i_p_bz,1)
   i_sp  =Xk%sstar(i_p_bz,2)
   !
   i_v   =BS_T_grp(i_T_grp)%table(i_T_el,2)
   i_c   =BS_T_grp(i_T_grp)%table(i_T_el,3)
   i_sp_c=BS_T_grp(i_T_grp)%table(i_T_el,4)
   i_sp_v=BS_T_grp(i_T_grp)%table(i_T_el,5)
   !
   K_EXCH_collision%is=(/i_c,i_k,i_s,i_sp_c/)
   K_EXCH_collision%os=(/i_v,i_p,i_sp,i_sp_v/)
   K_EXCH_collision%qs=(/i_g_p,iq,1/)
   !
   l_load_WFs_k= l_bs_exch_wf_in_loop .and. i_k_ref/=i_k
   l_load_WFs_p= l_bs_exch_wf_in_loop .and. i_p_ref/=i_p
   !
   if (l_load_WFs_k) then
     i_k_ref=i_k
#if defined(__NOTNOW) && ! defined(_CUDA)
     !$omp critical 
#endif
     call WF_load_fast(WF_k,NG(1),NG(2),BS_bands,i_k,(/1,n_sp_pol/),1,1)
#if defined(__NOTNOW) && ! defined(_CUDA)
     !$omp end critical
#endif
   endif
   !
   if (l_load_WFs_p) then
     i_p_ref=i_p
#if defined(__NOTNOW) && ! defined(_CUDA)
     !$omp critical 
#endif
     call WF_load_fast(WF_p,NG(1),NG(2),BS_bands,i_p,(/1,n_sp_pol/),1,ipath)
#if defined(__NOTNOW) && ! defined(_CUDA)
     !$omp end critical
#endif
   endif
   !
   if (iq==1) then
     !
     ! G==0 term
     !
     ! set terms on CPU and update GPU afterwards
     if (STRING_match(BSE_L_kind,"full")) BS_T_grp(i_T_grp)%O_x(1,i_T_el)=BS_T_grp(i_T_grp)%dipoles_opt(1,i_T_el,1)
     if (STRING_match(BSE_L_kind,"bar"))  BS_T_grp(i_T_grp)%O_x(1,i_T_el)=cZERO
     !
     call devxlib_memcpy_h2d( BS_T_grp(i_T_grp)%DEV_VAR(O_x), BS_T_grp(i_T_grp)%O_x, &
&                             range1=[1,1],range2=[i_T_el,i_T_el])
     !
     ! G/=0 terms
     !
     i_T_grp_last_with_identity_sym=-1
     i_T_el_last_with_identity_sym=-1
     if (i_s>1) then
       do i_T_grp_p=i_T_grp,1,-1
         if ( .not.allocated(BS_T_grp(i_T_grp_p)%O_x)      ) cycle
         if ( BS_T_grp(i_T_grp_p)%ik/=BS_T_grp(i_T_grp)%ik ) exit
         if ( i_T_grp_p==i_T_grp ) N_T_el_p=i_T_el-1
         if ( i_T_grp_p/=i_T_grp ) N_T_el_p=BS_T_grp(i_T_grp_p)%size
         do i_T_el_p=N_T_el_p,1,-1
           if ( all( BS_T_grp(i_T_grp_p)%table(i_T_el_p,:)==(/Xk%k_table(i_k,1),i_v,i_c,i_sp_c,i_sp_v/)) ) then
             i_T_grp_last_with_identity_sym=i_T_grp_p
             i_T_el_last_with_identity_sym =i_T_el_p
             exit
           endif
         enddo
       enddo
     endif
     !
     if  (i_s==1.or.i_T_grp_last_with_identity_sym<0) then
       !
       if (     l_bs_exch_wf_in_loop) call scatter_Bamp_gpu(K_EXCH_collision,WF_k,WF_p,.false.)
       if (.not.l_bs_exch_wf_in_loop) call scatter_Bamp_gpu(K_EXCH_collision,WF  ,WF  ,.false.)
       !
       call devxlib_memcpy_d2d( DEV_VAR(BS_T_grp(i_T_grp)%O_x)(:,i_T_el),&
                                DEV_VAR(K_EXCH_collision%rhotw),range1=(/2,BS_n_g_exch/) )
       !
     else
       !
       O_x_p     => DEV_VAR(BS_T_grp(i_T_grp)%O_x)
       O_x_sym_p => DEV_VAR(BS_T_grp(i_T_grp_last_with_identity_sym)%O_x)
       !
       !DEV_ACC_DEBUG data present(O_x_p,O_x_sym_p,g_rot,sop_inv)
       !DEV_ACC parallel loop private(i_g1,i_g2)
       !DEV_CUF kernel do (1) <<<*,*>>>
       !DEV_OMP parallel do default(shared), private(i_g1,i_g2)
       do i_g1=2,BS_n_g_exch
         i_g2=DEV_VAR(g_rot)(i_g1,DEV_VAR(sop_inv)(i_s))
         O_x_p(i_g1,i_T_el) = O_x_sym_p(i_g2,i_T_el_last_with_identity_sym)
       enddo
       !DEV_ACC_DEBUG end data
       !
       if (i_s>nsym/(i_time_rev+1)) call devxlib_conjg_d(BS_T_grp(i_T_grp)%DEV_VAR(O_x)(:,i_T_el),range1=(/2,BS_n_g_exch/) )
       !
     endif
     !
     call devxlib_memcpy_d2h( BS_T_grp(i_T_grp)%O_x, BS_T_grp(i_T_grp)%DEV_VAR(O_x), &
     &                             range1=[1,BS_n_g_exch],range2=[i_T_el,i_T_el])   
     !
   else
     !
     if (     l_bs_exch_wf_in_loop) call scatter_Bamp_gpu(K_EXCH_collision,WF_k,WF_p,.false.)
     if (.not.l_bs_exch_wf_in_loop) call scatter_Bamp_gpu(K_EXCH_collision,WF  ,WF  ,.false.)
     !
     ! BS_T_grp(i_T_grp)%DEV_VAR(O_x)(1,i_T_el)=cZERO
     call devxlib_memset_d(BS_T_grp(i_T_grp)%DEV_VAR(O_x)(:,i_T_el),val=cZERO,range1=[1,1])
     !
     if (STRING_match(BSE_L_kind,"full")) i_g0=1
     if (STRING_match(BSE_L_kind,"bar"))  i_g0=2
     !
     call devxlib_memcpy_d2d( DEV_VAR(BS_T_grp(i_T_grp)%O_x)(:,i_T_el), &
&                             DEV_VAR(K_EXCH_collision%rhotw), range1=(/i_g0,BS_n_g_exch/))
     !
   endif
   !
 enddo T_loop
#if defined(__NOTNOW)
 !DEV_OMP end do
#endif
 !
 call elemental_collision_free(K_EXCH_collision) 
 !
#if defined(__NOTNOW)
 !DEV_OMP end parallel
#endif
 !
 call timing('T_space EXX Osc.',OPR='stop')
 !
end subroutine K_exchange_collisions
