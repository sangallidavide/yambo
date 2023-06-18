!
!        Copyright (C) 2000-2023 the YAMBO team
!              http://www.yambo-code.org
!
! Authors (see AUTHORS file for details): AM DS AF
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
subroutine K_kernel(iq,Ken,Xk,q,X,Xw,W_bss)
 !
 ! K = <2V-W> for n_spin=1
 !
 use pars,          ONLY:SP,schlen,pi,cZERO,cI
 use LOGO,          ONLY:pickup_a_random
 use drivers,       ONLY:l_bs_fxc,l_tddft,l_rim_w
 use frequency,     ONLY:w_samp
 use electrons,     ONLY:levels,spin_occ,spin,n_spin,n_spinor,n_sp_pol
 use FFT_m,         ONLY:fft_size
 use interfaces,    ONLY:WF_load,WF_free,eval_G_minus_G,TDDFT_ALDA_R_space
 use D_lattice,     ONLY:nsym,DL_vol,i_time_rev,i_space_inv,sop_inv,sop_tab
 use R_lattice,     ONLY:G_m_G_maxval,qindx_X,bz_samp,ng_closed,&
&                        qindx_B,qindx_B_max,qindx_B_load,RIM_W,&
&                        RIM_ng,RIM_W_is_diagonal,RIM_W_ng
 use com,           ONLY:msg
 use stderr,        ONLY:intc
 use LIVE_t,        ONLY:live_timing
 use X_m,           ONLY:X_t
 use interfaces,    ONLY:el_density_and_current
 use QP_m,          ONLY:QP_ng_Vxc
 use BS,            ONLY:BS_bands,BS_res_K_corr,BS_W,BS_W_is_diagonal,&
&                        BS_res_K_exchange,BS_Block_size,&
&                        O_ng,BS_n_g_exch,BS_n_g_fxc,BS_identifier,BS_LiveTiming_steps,&
&                        BS_K_dim,BS_K_is_ALDA,BS_cpl_K_exchange,&
&                        BS_cpl_K_corr,K_EXCH_collision,K_CORR_collision,&
&                        BS_oscillators_free,WF_phase,n_BS_blks,BS_blk,BS_T_grp,&
&                        BS_nT_grps,BS_blks_free,l_BS_ares_from_res,&
&                        l_BSE_minimize_memory,l_BSE_restart,l_BSE_kernel_complete,&
&                        BS_perturbative_SOC,BS_K_cutoff,BS_max_val,l_BS_magnons
 use collision_el,  ONLY:elemental_collision_free
 use IO_int,        ONLY:io_control
 use IO_m,          ONLY:REP,OP_WR,RD,WR_CL,OP_APP,deliver_IO_error_message,WR
 use TDDFT,         ONLY:FXC_K_diagonal,F_xc_gspace,FXC_n_g_corr,io_BS_Fxc, &
&                        ALDA_cut_scheme,tddft_wf_t
 use xc_functionals,ONLY:V_xc,F_xc,F_xc_mat,magn,XC_potential_driver
 use global_XC,     ONLY:WF_xc_functional,WF_kind,WF_exx_fraction
 use wave_func,     ONLY:WF,rho_map,rho_map_thresh
 use openmp,        ONLY:OPENMP_update,master_thread
 use timing_m,      ONLY:timing
 use parallel_m,    ONLY:myid,ncpu,PAR_K_scheme,master_cpu
 use parallel_int,  ONLY:PP_redux_wait
 use MAGNONS,       ONLY:BSS_MAGN_free
 use PHOTOLUM,      ONLY:BSS_PL_free
 use parser_m,      ONLY:parser
 !
#include<memory.h>
 !
 type(levels)  ::Ken 
 type(bz_samp) ::Xk,q
 type(X_t)     ::X
 type(w_samp)  ::Xw,W_bss
 integer       ::iq
 !
 ! Work Space ...
 !
 type(tddft_wf_t), target  :: tddft_wf
 character(schlen) :: section_title
 character(1)      :: mode_now
 !
 ! ... Kernel loop
 !
 ! DEFINE BSK COMMON INDEXES
 !
 integer ::     i_k_s,i_k_s_m1,i_p_s,i_k_bz,i_p_bz,i_k,i_p,i_kp_s,& 
&               i_kmq_s,i_kmq_s_m1,i_pmq_s,i_kmq_bz,i_pmq_bz,i_kmq,i_pmq,i_kp_mq_s,& 
&               i_Tk,i_Tp,i_Tgrp_k,i_Tgrp_p,I_Tgrp_k_st,I_Tgrp_p_st,H_pos(2),& 
&               i_v_k,i_v_p,i_c_k,i_c_p,i_k_sp_pol_c,i_p_sp_pol_c,i_k_sp_pol_v,i_p_sp_pol_v,iq_W,iq_W_bz,iq_W_s,ig_W,& 
&               i_kmq_t,i_pmq_t,i_k_bz_mem,i_kmq_bz_mem,qindx_tmp(2) 
 ! OMP LAST DEF
 integer ::     i_k_bz_last,i_p_bz_last,i_pmq_bz_last,i_kmq_bz_last, & 
&               iq_W_bz_last,ig_W_last,iq_W_bz_mq_last,ig_W_mq_last 
 !
 integer    :: i_block,O_ng_shift,O_ng_exch,iHxc,NG(2),&
&              ig_kmq,ig_pmq,ig_W_mq,iq_W_bz_mq,iq_W_mq,iq_W_s_mq,&
&              bands_to_load(2),is_k(4),os_k(4),is_p(4),os_p(4)
 complex(SP):: Co,H_x,H_c,BS_mat_tmp
 complex(SP):: sqrt_f_itk,sqrt_f_itp
 real(SP), allocatable :: BS_max_val_CPU(:),rho(:)
 !
 complex(SP),   external ::K_correlation_kernel
 complex(SP),   external ::K_exchange_kernel_resonant,K_exchange_kernel_coupling
 !
 logical    :: l_bs_exch,l_bs_corr,l_matrix_init,l_load_kernel,l_write_kernel
 !
 ! I/O
 !
 integer           ::io_X_err,io_BS_Fxc_err,io_BS_err,qindx_ID,qindx_ID_frag,ID,ID_head,ID_compr,io_QINDX_err
 integer, external ::io_BS,io_BS_header
 integer, external ::qindx_B_init,qindx_B_close
 integer, external ::io_BS_PAR_init
 !
 ! Timing and restart
 !
 logical :: BS_blk_done(n_BS_blks),l_partial_kernel_loaded
 integer :: n_steps,m_steps,i_steps
 !
 ! Constants & Setup
 !==================
 !
 if(.not.BS_perturbative_SOC) Co=4._SP*real(spin_occ)*pi/DL_vol/real(q%nbz,SP)
 if(     BS_perturbative_SOC) Co=4._SP               *pi/DL_vol/real(q%nbz,SP)
 io_BS_err    =-1
 io_BS_Fxc_err=-1
 !
 call elemental_collision_free(K_CORR_collision)
 call elemental_collision_free(K_EXCH_collision)
 !
 ! K/S Table
 !
 call k_build_up_BZ_tables(Xk)
 !
 ! Titles
 !========
 !
 if (l_bs_fxc) then
   section_title="BSE Kernel -> xc-kernel"
 else
   if (BS_K_is_ALDA) then
     section_title="TDDFT Kernel"
   else
     section_title="BSE Kernel"
   endif
   section_title=trim(section_title)//" @q"//trim(intc(iq))
   !
   if (BS_res_K_corr.or.BS_res_K_exchange)     section_title=trim(section_title)//" (Resonant"
   if (BS_K_is_ALDA     )                      section_title=trim(section_title)//" ALDA"
   if (BS_res_K_corr    )                      section_title=trim(section_title)//" CORRRELATION"
   if (BS_res_K_exchange)                      section_title=trim(section_title)//" EXCHANGE"
   if (BS_res_K_corr.or.BS_res_K_exchange)     section_title=trim(section_title)//")"
   !
   if (BS_cpl_K_corr.or.BS_cpl_K_exchange)     section_title=trim(section_title)//" (Coupling"
   !if (BS_K_is_ALDA.and.BS_cpl_K_exchange)     section_title=trim(section_title)//" ALDA"
   !if (BS_cpl_K_corr    )                      section_title=trim(section_title)//" CORRRELATION"
   !if (BS_cpl_K_exchange)                      section_title=trim(section_title)//" EXCHANGE"
   if (BS_cpl_K_corr.or.BS_cpl_K_exchange)     section_title=trim(section_title)//")"
   !
 endif
 !
 call section('=',trim(section_title))
 !
 ! Oscillators Setup
 !===================
 O_ng=1
 O_ng_exch =maxval(qindx_X(iq,:,2))
 O_ng_shift=O_ng_exch
 !
 l_bs_exch=any((/BS_res_K_exchange,BS_cpl_K_exchange/)).or.l_bs_fxc.and..not.(l_BS_magnons.and.n_sp_pol==2)
 l_bs_corr=any((/BS_res_K_corr,BS_cpl_K_corr/)) .and. .not.l_bs_fxc
 !
 if (l_bs_exch) call msg('r','[BSE] Exchange components ',BS_n_g_exch)
 !
 if (l_bs_corr) then
   O_ng=G_m_G_maxval
   ! < WARNIGN: 2021/03/19
   ! The following can be safely removed with the new implementation 
   ! However it maybe needed for internal consistency. Otherwise the output of many BSE tests could change
   !! The cutoff on the non-hybrid part, BS_n_g_fxc, need to be lowered
   !call fft_check_size(1,O_ng,BS_n_g_exch,"Kx")
   !! Hybrid functional.
   !if (BS_K_is_ALDA) call fft_check_size(1,O_ng,BS_n_g_fxc,"fxc")
   ! WARNING >
   if(iq/=1) O_ng_shift=max(qindx_B_max,O_ng_shift)
 endif
 !
 if (BS_K_is_ALDA.and.trim(ALDA_cut_scheme)=="lower_GmGp") then
   FXC_n_g_corr=BS_n_g_fxc
   if (BS_n_g_fxc/= eval_G_minus_G(BS_n_g_fxc,0) ) then
     BS_n_g_fxc   = eval_G_minus_G(BS_n_g_fxc,0)
     call msg('s','[TDDFT] ALDA Fxc cutoff is forced by symmetry to be ',BS_n_g_fxc)
   endif
 endif
 !
 if (BS_K_is_ALDA) call msg('r','[BSE] ALDA     components ',BS_n_g_fxc)
 !
 call K_restart(iq,X,ID_head,ID,ID_compr,l_partial_kernel_loaded)
 !
 if( l_BSE_kernel_complete ) return
 !
 io_QINDX_err=0
 if (l_bs_corr) io_QINDX_err=qindx_B_init(qindx_ID,qindx_ID_frag)
 if(io_QINDX_err/=0) call error('Error reading qindx_B database ')
 !
 ! Screened interaction
 !======================
 !
 !Check if RIM_W is diagonal
 call parser('RIM_W_diag',RIM_W_is_diagonal)
 if(l_rim_w .and. RIM_W_ng == 0) RIM_W_ng = RIM_ng
 if (l_rim_w) call interpolate_W(X,Xw,q,'BSE')
 !
 io_X_err=0
 if (l_bs_corr) call K_screened_interaction(X,Xw,q,io_X_err)
 !
 if (io_X_err<0) then
   call deliver_IO_error_message(io_X_err,'PP/Em1s')
   call warning('Bethe Salpter section skipped. Impossible to build the kernel.')
   n_BS_blks=0
   call BS_blks_free(1)
   call BSS_MAGN_free( )
   call BSS_PL_free( )
   !
   YAMBO_FREE(BS_W)
   YAMBO_FREE(RIM_W)
   !
   return
   !
 endif
 !
 call timing('X T_space Kernel',OPR='start')
 !
 ! TDDFT xc-kernel Setup
 !======================
 if (l_bs_fxc) then
   !  
   !if (iq==1) here put dipoles I/O
   if (iq==1) stop
   !
   ! Allocation
   !-----------
   YAMBO_ALLOC(FXC_K_diagonal,(BS_K_dim(1)))
   YAMBO_ALLOC(F_xc_gspace,(FXC_n_g_corr,FXC_n_g_corr,W_bss%n_freqs))
   FXC_K_diagonal=0._SP
   F_xc_gspace=cZERO
   !
 endif
 !
 ! DB identifier 
 !===============
 BS_identifier=pickup_a_random(10000._SP)
 !
 ! BS DB description section I/O
 !===============================
 !
 call section('=','Main kernel loop')
 !
 if(BS_K_cutoff>0._SP) then
   YAMBO_ALLOC(BS_max_val_CPU,(ncpu))
   BS_max_val_CPU=0._SP
 endif
 !
#if defined _PAR_IO
 !
 call io_control(ACTION=OP_APP,COM=REP,ID=ID_head,DO_IT=.TRUE.)
 io_BS_err=io_BS_header(iq,X,ID_head,"full")
 !
 call io_control(ACTION=OP_APP,COM=REP,ID=ID,DO_IT=.TRUE.)
 io_BS_err=io_BS_PAR_init(iq,ID,"full")
 !
#else
 call io_control(ACTION=OP_WR,COM=REP,SEC=(/1/),ID=ID)
 io_BS_err=io_BS(iq,X,ID)
#endif
 !
 ! This loop is on the exchange and correlation of the BSE
 ! iHxc=1 --> K_H=v
 ! iHxc=2 --> K_TDDFT=fxc
 ! iHxc=3 --> K_c=W
 !
 l_matrix_init=.true.
 l_load_kernel=.true.
 l_write_kernel=.false.
 !
 do iHxc=1,3
 !
 if (iHxc==1) BS_blk_done=.false.
 !
 if ((.not.l_bs_exch).and.iHxc==1) cycle
 !
 if ((.not.BS_K_is_ALDA).and.iHxc==2) cycle
 !
 if ((.not.l_bs_corr).and.iHxc==3) cycle
 !
 section_title='-BSK-'//trim(intc(iHxc))
 if (l_bs_fxc) section_title='-BSK->Fxc-'//trim(intc(iHxc))
 !
 if (iHxc==1) NG=(/BS_n_g_exch,O_ng_exch/)
 if (iHxc==2) NG=(/BS_n_g_fxc, 1/)
 ! < WARNIGN: 2021/03/19
 ! The following has been changed to the commente line 
 ! However the commented line maybe needed for internal consistency. Otherwise the output of many BSE tests could change
 if (iHxc==3) NG=(/O_ng,O_ng_shift/)
 !if (iHxc==3) NG(2)=O_ng_shift
 ! WARNING >
 !
 ! ALDA
 !======
 if (iHxc==2) then
   !
   ! This part could be done just by the master bypassing the PAR_IND_WF_linear ...
   !
   if(rho_map_thresh>0._SP) then
     call WF_load(WF,NG(1),NG(2),(/1,Ken%nbm/),(/1,Xk%nibz/),space='R',title="-TDDFT",keep_states_to_load=.true.)
     YAMBO_ALLOC(rho,(fft_size))
     call el_density_and_current(Ken,Xk,rho=rho,compute_rho_map=.true.)
     YAMBO_FREE(rho)
     call WF_free(WF,keep_states_to_load=.true.)
   endif
   !
   call WF_load(WF,NG(1),NG(2),(/1,Ken%nbm/),(/1,Xk%nibz/),space='R',title="-TDDFT",keep_states_to_load=.true.)
   !
   call TDDFT_ALDA_lower_cutoff(ALDA_cut_scheme)
   !
   ! In principle here there is a factor 4 between the energy cutoff on the potential and the one of the kernel
   ! The following line, or something similar, should be used here and before
   !QP_nG=G_index_energy_factor(BS_n_g_fxc,0.25_SP)
   ! However in practice the one on the kernel is never set 4 the potential, i.e. 16 times the one onf the WFs
   QP_ng_Vxc=BS_n_g_fxc
   !
   YAMBO_ALLOC(tddft_wf%rhotwr1,(fft_size*n_spinor*n_spinor))
   YAMBO_ALLOC(tddft_wf%rhotwr2,(fft_size*n_spinor*n_spinor))
   YAMBO_ALLOC(tddft_wf%WF_symm1,(fft_size,n_spinor))
   YAMBO_ALLOC(tddft_wf%WF_symm2,(fft_size,n_spinor))
   !
   YAMBO_ALLOC(F_xc,(fft_size,n_spin,n_spin))
   YAMBO_ALLOC(V_xc,(fft_size,n_spin))
   V_xc=0._SP
   !
   if(n_spin>1) then
     YAMBO_ALLOC(magn,(fft_size,3))
   endif
   call XC_potential_driver(Ken,Xk,WF_kind,WF_xc_functional,2)
   YAMBO_ALLOC(F_xc_mat,(fft_size,n_spin,n_spin,n_spin,n_spin))
   call Build_F_xc_mat(V_xc,F_xc,F_xc_mat)
   YAMBO_FREE(F_xc)
   YAMBO_FREE(V_xc)
   if(n_spin>1) then
     YAMBO_FREE(magn)
   endif
   !
   call WF_free(WF,keep_states_to_load=.true.)
   !
 endif
 !
 call WF_load(WF,NG(1),NG(2),BS_bands,(/1,Xk%nibz/),space='R',title=trim(section_title),keep_states_to_load=.true.)
 !
 ! Wave Functions Phases
 !=======================
 if (iHxc==3) call K_WF_phases(Xk)
 ! 
 ! Spatial Inversion Test
 !========================
 if (iHxc==1) call WF_spatial_inversion(Ken,Xk)
 !
 ! Timing
 !========
 !
 if (BS_LiveTiming_steps>= 0) n_steps=BS_LiveTiming_steps
 if (BS_LiveTiming_steps==-1) n_steps=n_BS_blks
 !
 m_steps=n_steps
 !
 ! First load previously computed matrix elements.
 !=================================================
 !
#ifdef _PAR_IO
 if (l_BSE_restart.and.l_load_kernel) then
   !
   if (m_steps>0) call live_timing("Loading partial kernel",m_steps)
   !
   do i_block=1,n_BS_blks
     !
     ! Read BS_blk to check if the transition was already computed
     !
     if (.not.l_partial_kernel_loaded) then
       call io_control(ACTION=RD,ID=ID)
       call io_BS_PAR_block(iq,i_block,ID,"full")
     endif
     !
     if (BS_LiveTiming_steps>= 0) i_steps=BS_Block_size(i_block)
     if (BS_LiveTiming_steps==-1) i_steps=1
     !
     i_Tgrp_k=BS_blk(i_block)%iT_k
     i_Tgrp_p=BS_blk(i_block)%iT_p
     !
     ! check if this block need to be done
     ! skip this block if already done
     !
     BS_blk_done(i_block)=.true.
     skip_check: do i_Tp=1,BS_T_grp(i_Tgrp_p)%size
       do i_Tk=1,BS_T_grp(i_Tgrp_k)%size
         if (BS_blk(i_block)%done(i_Tk,i_Tp)=="t") cycle
         BS_blk_done(i_block)=.false.
         exit skip_check
       enddo
     enddo skip_check
     !
     if (BS_blk_done(i_block)) m_steps=m_steps-i_steps
     !
     call live_timing(steps=i_steps)
     !
   enddo
   !
   call live_timing()
   !
   call msg("s","Kernel loaded percentual ",real(n_steps-m_steps,SP)/real(n_steps,SP)*100._SP,'[%]')
   !
   l_load_kernel=.false.
   !
 endif
#endif
 !
 if (.not.l_bs_fxc) section_title='Kernel'
 if (     l_bs_fxc) section_title='Kernel->Fxc'
 !
 if (iHxc==1) section_title=trim(section_title)//' exch'
 if (iHxc==2) section_title=trim(section_title)//' tddft'
 if (iHxc==3) section_title=trim(section_title)//' corr'
 !
 if (m_steps>0) call live_timing(trim(section_title),m_steps)
 !
 !------------------
 ! LOOP over blocks
 !------------------
 !
 block_loop:&
 do i_block=1,n_BS_blks
   !
   mode_now=BS_blk(i_block)%mode
   if (.not.l_BS_ares_from_res) mode_now="F"
   !
   if (.not.trim(BS_blk(i_block)%mode)=="C") then
     if (iHxc==1) l_write_kernel=.not.(BS_K_is_ALDA.or.l_bs_corr)
     if (iHxc==2) l_write_kernel=.not. l_bs_corr
     if (iHxc==3) l_write_kernel=.true.
   else
     if (iHxc==1) l_write_kernel=.not.(BS_K_is_ALDA.or.BS_cpl_K_corr)
     if (iHxc==2) l_write_kernel=.not. BS_cpl_K_corr
     if (iHxc==3) l_write_kernel=.true.
     !
     if (iHxc==3.and..not.BS_cpl_K_corr) cycle
     !
   endif
   !
   if (BS_blk_done(i_block)) then
#if defined _PAR_IO
     if(l_partial_kernel_loaded.and.l_write_kernel) then
       call io_control(ACTION=WR,ID=ID)
       call io_BS_PAR_block(iq,i_block,ID,"full")
     endif
#endif
     cycle
   endif
   !
   ! Read BS_blk to check if the transition was already computed
   !
   if (BS_LiveTiming_steps>= 0) i_steps=BS_Block_size(i_block)
   if (BS_LiveTiming_steps==-1) i_steps=1
   !
   i_Tgrp_k   =BS_blk(i_block)%iT_k
   i_Tgrp_p   =BS_blk(i_block)%iT_p
   i_Tgrp_k_st=BS_blk(i_block)%iT_k_st
   i_Tgrp_p_st=BS_blk(i_block)%iT_p_st
   i_k        =BS_blk(i_block)%ik
   i_p        =BS_blk(i_block)%ip
   !
   ! Exchange oscillators
   !----------------------
   !
   if(iHxc==1) then
     call K_exchange_collisions(iq,Xk,i_Tgrp_k)
     if (i_Tgrp_k/=i_Tgrp_p) call K_exchange_collisions(iq,Xk,i_Tgrp_p)
   endif
   !
   ! Exchange oscillators
   !----------------------
   ! 
   if(iHxc==3) call K_correlation_collisions(iq,i_block,qindx_ID_frag,Xk,q)
   !
   ! Transition Loops
   !------------------
   !
#if defined _OPENMP && ! defined _CUDA
   !$omp parallel default(shared), &
   !$omp private(i_Tk,i_Tp, &
   !$omp &       sqrt_f_itk,i_k_bz,i_v_k,i_c_k,i_k_sp_pol_c,i_k_sp_pol_v,i_k_s,is_k,os_k,i_k_s_m1,H_pos, &
   !$omp &       sqrt_f_itp,i_p_bz,i_v_p,i_c_p,i_p_sp_pol_c,i_p_sp_pol_v,i_p_s,is_p,os_p,i_kp_s, &
   !$omp &       i_kmq_bz,ig_kmq,i_kmq,i_kmq_s,i_kmq_t,i_kmq_s_m1,i_kp_mq_s, &
   !$omp &       i_pmq_bz,ig_pmq,i_pmq,i_pmq_s,i_pmq_t, &
   !$omp &       ig_W,iq_W,iq_W_s,iq_W_bz,iq_W_bz_mq,ig_W_mq,iq_W_mq,iq_W_s_mq,H_x,H_c, &
   !$omp &       i_p_bz_last,i_k_bz_last,i_pmq_bz_last,i_kmq_bz_last, &
   !$omp &       iq_W_bz_last,ig_W_last,iq_W_bz_mq_last,ig_W_mq_last, &
   !$omp &       i_k_bz_mem,i_kmq_bz_mem,qindx_tmp,BS_mat_tmp,tddft_wf )
#endif
   !
   call OPENMP_update(master_thread)   
   !
   i_p_bz_last=0   ; i_k_bz_last=0
   i_pmq_bz_last=0 ; i_kmq_bz_last=0
   !
#if defined _OPENMP && ! defined _CUDA
   !$omp do collapse(2), schedule(dynamic)
#endif
   do i_Tp=1,BS_T_grp(i_Tgrp_p)%size
     do i_Tk=1,BS_T_grp(i_Tgrp_k)%size
       !
       ! skip this block if already done
       !
       if (BS_blk(i_block)%done(i_Tk,i_Tp)=="t".and.l_BSE_restart) cycle
       !
       H_pos(1) = sum(BS_T_grp(I_Tgrp_k_st:I_Tgrp_k-1)%size)+i_Tk+& 
       &              (BS_T_grp(I_Tgrp_k)%i_res_ares-1)*BS_K_dim(1) 
       H_pos(2) = sum(BS_T_grp(I_Tgrp_p_st:I_Tgrp_p-1)%size)+i_Tp+& 
       &              (BS_T_grp(I_Tgrp_p)%i_res_ares-1)*BS_K_dim(1) 
       if (H_pos(1)>H_pos(2)) cycle 
       !
       i_k_bz  = BS_T_grp(i_Tgrp_k)%table(i_Tk,1) 
       i_p_bz  = BS_T_grp(i_Tgrp_p)%table(i_Tp,1) 
       !
       i_k_s =Xk%sstar(i_k_bz,2) 
       i_p_s =Xk%sstar(i_p_bz,2) 
       !
       i_k_s_m1= sop_inv(i_k_s) 
       i_kp_s  = sop_tab(i_k_s_m1,i_p_s) 
       !
       i_kmq_bz=qindx_X(iq,i_k_bz,1) 
       i_pmq_bz=qindx_X(iq,i_p_bz,1) 
       !
       i_kmq   =Xk%sstar(i_kmq_bz,1) 
       i_pmq   =Xk%sstar(i_pmq_bz,1) 
       !
       i_kmq_s =Xk%sstar(i_kmq_bz,2) 
       i_pmq_s =Xk%sstar(i_pmq_bz,2) 
       !
       if (iHxc==3) then 
         i_kmq_t=BS_blk(i_block)%kp_table(1,i_kmq) 
         i_pmq_t=BS_blk(i_block)%kp_table(2,i_pmq) 
       endif 
       !
       i_kmq_s_m1 = sop_inv(i_kmq_s) 
       i_kp_mq_s  = sop_tab(i_kmq_s_m1,i_pmq_s) 
       !
       if (iHxc==3) then 
         i_k_bz_mem=PAR_K_scheme%bz_index(i_k_bz) 
         if (i_p_bz_last/=i_p_bz.or.i_k_bz_last/=i_k_bz) then 
           i_p_bz_last=i_p_bz 
           i_k_bz_last=i_k_bz 
           if (i_k_bz_mem==0) then 
             !$omp critical 
             qindx_tmp=qindx_B_load(i_p_bz,i_k_bz,qindx_ID_frag) 
             iq_W_bz=qindx_tmp(1) 
             ig_W   =qindx_tmp(2) 
             !$omp end critical 
           else 
             iq_W_bz=qindx_B(i_p_bz,i_k_bz_mem,1) 
             ig_W   =qindx_B(i_p_bz,i_k_bz_mem,2) 
           endif 
           iq_W_bz_last=iq_W_bz 
           ig_W_last   =ig_W    
         else 
           iq_W_bz=iq_W_bz_last 
           ig_W   =ig_W_last    
         endif 
         iq_W   =q%sstar( iq_W_bz ,1) 
         iq_W_s =q%sstar( iq_W_bz ,2) 
       endif 
       !
       i_v_k=BS_T_grp(i_Tgrp_k)%table(i_Tk,2) 
       i_c_k=BS_T_grp(i_Tgrp_k)%table(i_Tk,3) 
       !
       i_v_p=BS_T_grp(i_Tgrp_p)%table(i_Tp,2) 
       i_c_p=BS_T_grp(i_Tgrp_p)%table(i_Tp,3) 
       !
       i_k_sp_pol_c=BS_T_grp(i_Tgrp_k)%table(i_Tk,4) 
       i_p_sp_pol_c=BS_T_grp(i_Tgrp_p)%table(i_Tp,4) 
       i_k_sp_pol_v=BS_T_grp(i_Tgrp_k)%table(i_Tk,5) 
       i_p_sp_pol_v=BS_T_grp(i_Tgrp_p)%table(i_Tp,5) 
       !
       ! 
       ig_kmq    = qindx_X(iq,i_k_bz,2) 
       ig_pmq    = qindx_X(iq,i_p_bz,2) 
       !        
       if (iHxc==3) then 
         i_kmq_bz_mem=PAR_K_scheme%bz_index(i_kmq_bz) 
         if (i_pmq_bz_last/=i_pmq_bz.or.i_kmq_bz_last/=i_kmq_bz) then 
           i_pmq_bz_last=i_pmq_bz 
           i_kmq_bz_last=i_kmq_bz 
           if (i_kmq_bz_mem==0) then 
             !$omp critical 
             qindx_tmp=qindx_B_load(i_pmq_bz,i_kmq_bz,qindx_ID_frag) 
             iq_W_bz_mq=qindx_tmp(1) 
             ig_W_mq   =qindx_tmp(2) 
             !$omp end critical 
           else 
             iq_W_bz_mq=qindx_B(i_pmq_bz,i_kmq_bz_mem,1) 
             ig_W_mq   =qindx_B(i_pmq_bz,i_kmq_bz_mem,2) 
           endif 
           iq_W_bz_mq_last=iq_W_bz_mq 
           ig_W_mq_last   =ig_W_mq    
         else 
           iq_W_bz_mq=iq_W_bz_mq_last 
           ig_W_mq   =ig_W_mq_last    
         endif 
         iq_W_mq   =q%sstar( iq_W_bz_mq,1) 
         iq_W_s_mq =q%sstar( iq_W_bz_mq,2) 
       endif 
       ! 
       if (iHxc==2) then 
         is_k = (/i_c_k,i_k,i_k_s,i_k_sp_pol_c/) 
         os_k = (/i_v_k,i_kmq,i_kmq_s,i_k_sp_pol_v/) 
         is_p = (/i_c_p,i_p,i_p_s,i_p_sp_pol_c/) 
         os_p = (/i_v_p,i_pmq,i_pmq_s,i_p_sp_pol_v/) 
       endif 
       !
       if (iHxc==3) then 
         !
         if (iq_W_bz/=iq_W_bz_mq.or.iq_W/=iq_W_mq.or.iq_W_s/=iq_W_s_mq) call error("Wrong transferred momentum")
         !
         if ( (.not.BS_W_is_diagonal) .and. iq_W_s>nsym/(i_time_rev+1) .and. i_space_inv == 0 ) iq_W=q%nibz+iq_W
         !
       endif
       !
       H_x=cZERO
       H_c=cZERO
       !
       select case (mode_now)
       !
       case ("R","A","F")
         !
         ! :::RESONANT/ANTI-RESONANT TERM:::
         !===================================
         !
         ! Exchange
         if (iHxc==1) H_x=K_exchange_kernel_resonant(iq, BS_n_g_exch,               &
&                                   BS_T_grp(i_Tgrp_p),i_Tp, BS_T_grp(i_Tgrp_k),i_Tk)
         !
         ! ALDA
         if (iHxc==2) H_x=H_x+TDDFT_ALDA_R_space(is_k,os_k,is_p,os_p,tddft_wf,'RES')
         !
         ! Correlations
         if (iHxc==3) H_c=K_correlation_kernel(i_block,i_p,i_pmq, &
&                                   i_k_s,i_kp_s,i_c_k,i_c_p,i_kmq_s,i_kp_mq_s,i_v_k,i_v_p,     &
&                                   i_kmq_t,i_pmq_t,i_k_sp_pol_c,i_p_sp_pol_c,i_k_sp_pol_v,i_p_sp_pol_v,  &
&                                   iq_W,iq_W_s,ig_W,i_k_s_m1,iq_W_s_mq,ig_W_mq,i_kmq_s_m1)
         !         
       case ("C","Q")
         !
         ! :::COUPLING TERM:::
         !=============================
         !
         ! Exchange
         if (iHxc==1) H_x=K_exchange_kernel_coupling(iq, BS_n_g_exch,              &
&                                 BS_T_grp(i_Tgrp_p),i_Tp, BS_T_grp(i_Tgrp_k),i_Tk) 
         !
         ! ALDA
         if (iHxc==2) H_x=H_x+TDDFT_ALDA_R_space(is_k,os_k,is_p,os_p,tddft_wf,'CPL')
         !
         ! Correlations
         if (BS_cpl_K_corr.and.iHxc==3) H_c=K_correlation_kernel(i_block,i_p,i_pmq, &
&                                 i_k_s,i_kp_s,i_c_k,i_v_p,i_kmq_s,i_kp_mq_s,i_v_k,i_c_p,    &
&                                 i_kmq_t,i_pmq_t,i_k_sp_pol_c,i_p_sp_pol_v,i_k_sp_pol_v,i_p_sp_pol_c,  &
&                                 iq_W,iq_W_s,ig_W,i_k_s_m1,iq_W_s_mq,ig_W_mq,i_kmq_s_m1)
         !
       end select
       !
       if (l_tddft.and.WF_exx_fraction>0._SP.and.iHxc==3) H_c=H_c*WF_exx_fraction
       !
       sqrt_f_itk=sqrt(cmplx(BS_T_grp(i_Tgrp_k)%f(i_Tk),kind=SP))
       sqrt_f_itp=sqrt(cmplx(BS_T_grp(i_Tgrp_p)%f(i_Tp),kind=SP))
       !
       BS_mat_tmp=sqrt_f_itk*(H_x*Co-H_c)*sqrt_f_itp
       if (mode_now=="C") BS_mat_tmp= cI*BS_mat_tmp
       !
       if (     l_matrix_init) BS_blk(i_block)%mat(i_Tk,i_Tp) = BS_mat_tmp
       if (.not.l_matrix_init) BS_blk(i_block)%mat(i_Tk,i_Tp) = BS_blk(i_block)%mat(i_Tk,i_Tp)+BS_mat_tmp
       !
       if (l_write_kernel) BS_blk(i_block)%done(i_Tk,i_Tp)="t"
       !
     enddo   ! i_Tk
   enddo     ! i_Tp
#if defined _OPENMP && ! defined _CUDA
   !$omp end do
   !$omp end parallel
#endif
   !
   if(BS_K_cutoff>0._SP.and.l_write_kernel) &
&     BS_max_val_CPU(myid+1)=maxval( (/BS_max_val_CPU(myid+1),abs(BS_blk(i_block)%mat(:,:))/) )
   !
   if(l_BSE_minimize_memory.and.iHxc==1) then
     call BS_oscillators_free(i_Tgrp_p,0)
     call BS_oscillators_free(i_Tgrp_k,0)
   endif
   !
   if(iHxc==3) call BS_oscillators_free(0,i_block)
   !
   if(l_write_kernel) then
#if defined _PAR_IO
     call io_control(ACTION=WR,ID=ID)
     call io_BS_PAR_block(iq,i_block,ID,"full")
#else
     call io_control(ACTION=WR,COM=REP,SEC=(/i_block+1/),ID=ID)
     io_BS_err=io_BS(iq,X,ID)
#endif
   endif
   !
   call live_timing(steps=i_steps)
   !
 enddo block_loop
 !
 if (n_steps>0) call live_timing( )
 !
 ! GLOBAL CLEANING 
 !=================
 !
 if(iHxc==1) call BS_oscillators_free(BS_nT_grps,0)
 !
 if(iHxc==2) then
   !
   YAMBO_FREE(tddft_wf%rhotwr1)
   YAMBO_FREE(tddft_wf%rhotwr2)
   YAMBO_FREE(tddft_wf%WF_symm1)
   YAMBO_FREE(tddft_wf%WF_symm2)
   !
   if(rho_map_thresh>0._SP) deallocate(rho_map)
   !
   if (l_bs_fxc) then
     YAMBO_FREE(F_xc_gspace)
     YAMBO_FREE(FXC_K_diagonal)
   else if (BS_K_is_ALDA) then
     YAMBO_FREE(F_xc_mat)
   endif
   !
 endif
 !
 if(iHxc==3) then
   YAMBO_FREE(WF_phase)
   qindx_tmp=qindx_B_close(qindx_ID,qindx_ID_frag)
   YAMBO_FREE(BS_W)
   YAMBO_FREE(RIM_W)
 endif
 !
 call WF_free(WF,keep_states_to_load=iHxc<3)
 !
 l_matrix_init=.false.
 !
 enddo ! iHxc loop
 !
#if defined _PAR_IO
 if(BS_K_cutoff>0._SP) then
   call PP_redux_wait(BS_max_val_CPU)
   BS_max_val=maxval(BS_max_val_CPU)
   call K_compress(iq,ID_compr)
   call K_compressed_IO(iq,io_BS_err,"write")
 endif
 !
 call io_control(ACTION=WR_CL,COM=REP,ID=ID_head)
 call io_control(ACTION=WR_CL,COM=REP,ID=ID)
 call io_BS_PAR_free(ID_head,ID,"full",.true.)
#else
 call io_control(ACTION=WR_CL,COM=REP,SEC=(/n_BS_blks+2/),ID=ID)
 io_BS_err=io_BS(iq,X,ID)
#endif
 !
 call timing('X T_space Kernel',OPR='stop')
 !
end subroutine K_kernel