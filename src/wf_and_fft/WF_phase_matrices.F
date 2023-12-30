!
! License-Identifier: GPL
!
! Copyright (C) 2015 The Yambo Team
!
! Authors (see AUTHORS file for details): DS
!      
!> @brief Calculate overlaps between wave-function
!   <\psi_{Ra(kibz_i)} | \psi_{kbz_j} >
!  with for Ra(kibz_i)=kbz_j
!
!! @param[out]      WF_phases_S            Phase matrices
!
! #include<dev_defs.h>
!
subroutine WF_phase_matrices(En,Xk,blim)
 !    
 use pars,              ONLY:SP,cZERO
 use units,             ONLY:HA2EV
 use LIVE_t,            ONLY:live_timing
 use R_lattice,         ONLY:bz_samp,G_m_G,minus_G,g_vec,rl_sop,g_rot
 use D_lattice,         ONLY:nsym,sop_tab
 use electrons,         ONLY:levels
 use wave_func,         ONLY:WAVEs,io_WF_phases,&
 &                           WF_phases_b_map,WF_phases_s_map
 use wave_func,         ONLY:wf_ng_overlaps,wf_ng_1st_BZ,WF_copy
 use electrons,         ONLY:n_spinor,n_sp_pol
 use wrapper,           ONLY:Vstar_dot_V
 use hamiltonian,       ONLY:B_mat_index
 use matrix_operate,    ONLY:m3det!,UP_matrix_index
 use vec_operate,       ONLY:degeneration_finder
 use parallel_m,        ONLY:PAR_IND_WF_k,PAR_IND_WF_b,PAR_IND_B_mat
 use parallel_int,      ONLY:PP_redux_wait,PARALLEL_WF_index,PARALLEL_WF_distribute
 use interfaces,        ONLY:WF_symm_kpoint_g,WF_load,WF_free,eval_G_minus_G
 use IO_int,            ONLY:io_control
 use IO_m,              ONLY:VERIFY,DUMP,REP,OP_WR,WR,WR_CL,OP_RD_CL
 use timing_m,          ONLY:timing
 !
#include<memory.h>
 !
 integer,       intent(in) :: blim(2)
 type(levels),  intent(in) :: En
 type(bz_samp), intent(in) :: Xk 
 !
 ! Work Space
 !
 integer           :: first_el(blim(2)-blim(1)+1),n_of_el(blim(2)-blim(1)+1),&
 &                    blim_n(2),n_deg_grp,i_grp
 integer           :: ikbz,ipbz,isk,isp,istark,ibk,ibp,ibk_phase,ibp_phase,i_sp_pol,ik,ip,&
 &                    ifrag,s_count,i_wf_k,i_wf_p,iGp,ng_tmp,nb
 complex(SP)       :: det,charac(blim(2)-blim(1)+1,nsym)
 complex(SP), allocatable :: WF_phases(:,:,:,:,:,:)
 !
 type(WAVEs) :: WF_k,WF_p,WF_p_tmp
 !
 integer              :: io_err,ID
 !
 if (nsym==1) return
 !
 call io_control(ACTION=OP_RD_CL,COM=REP,SEC=(/1/),MODE=VERIFY,ID=ID)
 io_err=io_WF_phases(blim,0,0,ID)
 !
 if (io_err==0) return
 !
 call timing("WF_phases",OPR="start")
 !!
 !call k_build_up_BZ_tables(Xk)
 !
 iGp=15
 ng_tmp=eval_G_minus_G(wf_ng_1st_BZ,iGp)
 !
 nb=blim(2)-blim(1)+1
 !
 !call PARALLEL_WF_distribute(K_index=PAR_IND_WF_k,B_index=PAR_IND_WF_b,CLEAN_UP=.TRUE.)
 call PARALLEL_WF_distribute(CLEAN_UP=.TRUE.)
 call PARALLEL_WF_index( )
 !
 YAMBO_ALLOC(WF_phases_s_map,(Xk%nbz,nsym))
 YAMBO_ALLOC(WF_phases_b_map,(blim(1):blim(2),Xk%nbz,nsym,2))
 !
 call io_control(ACTION=OP_WR,COM=REP,SEC=(/1/),MODE=DUMP,ID=ID)
 io_err=io_WF_phases(blim,0,0,ID)
 !
 WF_phases_s_map=0
 WF_phases_b_map=0
 !
 !call live_timing("WF_phases:",PAR_IND_WF_k%n_of_elements(PAR_IND_DIPk_bz_ID+1)*n_sp_pol)
 call live_timing("WF_phases:",Xk%nbz*n_sp_pol)
 !
 ifrag=0
 !
 do i_sp_pol=1,n_sp_pol
   !
   do ik=1,Xk%nibz
     !
     call degeneration_finder(En%E(blim(1):blim(2),ik,i_sp_pol),nb,&
     &                           first_el,n_of_el,n_deg_grp,0.0001/HA2EV,Include_single_values=.true.)
     !
     do istark=1,Xk%nstar(ik)
       !
       !if ( .not.PAR_IND_WF_k%element_1D(ik) ) cycle
       !
       isk =Xk%star(ik,istark)
       ikbz=Xk%k_table(ik,isk)
       !
       !AMBO_ALLOC(WF_phases,(maxval(n_of_el),maxval(n_of_el),n_deg_grp,nsym,1,1))
       YAMBO_ALLOC(WF_phases,(3,3,nb,nsym,1,1))
       WF_phases=cZERO
       !
       ifrag=ifrag+1
       !
       call WF_load(WF_k,0,1,blim,(/ik,ik/),(/i_sp_pol,i_sp_pol/),&
       & space='G',title='-WF_phases',quiet=ifrag>1,keep_states_to_load=.true.)
       !
       ip=ik
       !
       call WF_load(WF_p,0,1,blim,(/ik,ik/),(/i_sp_pol,i_sp_pol/),&
       & space='G',title='-WF_phases',quiet=.true.,   keep_states_to_load=.true.)
       !
       call WF_symm_kpoint_g(blim,ik,i_sp_pol,isk,1,WF_k)
       !
       s_count=0
       !
       charac=0._SP
       !
       do isp=1,nsym
         !
         ipbz=Xk%k_table(ip,isp)
         iGp =Xk%g_table(ip,isp)
         !
         iGp=minus_G(iGp)
         !
         if (ipbz/=ikbz) cycle
         !
         s_count=s_count+1
         !
         WF_phases_s_map(ikbz,isp)=s_count
         !
         call WF_copy(WF_p,WF_p_tmp)
         !
         call WF_symm_kpoint_g(blim,ip,i_sp_pol,isp,iGp,WF_p_tmp)
         !
         do i_grp=1,n_deg_grp
           blim_n=(/1,n_of_el(i_grp)/)
           blim_n=blim_n+first_el(i_grp)-1+blim(1)-1
           !
           do ibk_phase=1,n_of_el(i_grp)
             ibk=ibk_phase+first_el(i_grp)-1+blim(1)-1
             WF_phases_b_map(ibk,ikbz,isp,:)=(/ibk_phase,i_grp/)
             i_wf_k=WF_k%index(ibk,ik,i_sp_pol)
             do ibp_phase=1,n_of_el(i_grp)
               ibp=ibp_phase+first_el(i_grp)-1+blim(1)-1
               i_wf_p=WF_p%index(ibp,ip,i_sp_pol)
               !
               WF_phases(ibk_phase,ibp_phase,i_grp,s_count,1,1)= &
               & Vstar_dot_V(wf_ng_1st_BZ*n_spinor,WF_k%c(:,:,i_wf_k),WF_p_tmp%c(:,:,i_wf_p))
               !
             enddo
             charac(i_grp,s_count)=charac(i_grp,s_count)+WF_phases(ibk_phase,ibk_phase,i_grp,s_count,1,1)
           enddo
           !
         enddo
         !
         !do ibk=blim(1),blim(2)
         !  i_wf_k=WF_k%index(ibk,ik,i_sp_pol)
         !  do ibp=blim(1),blim(2)
         !    i_wf_p=WF_p%index(ibp,ip,i_sp_pol)
         !    !if (.not.PAR_IND_B_mat_ordered%element_1D( UP_matrix_index(ib,ibp)-1 )) cycle
         !    !if (.not.PAR_IND_B_mat%element_1D( B_mat_index(ib,ibp,(/1,En%blim/))-1 )) cycle
         !    WF_phases(ibk,ibp,s_count,k_count,i_sp_pol)= &
         !    & Vstar_dot_V(wf_ng_1st_BZ*n_spinor,WF_k%c(:,:,i_wf_k),WF_p_tmp%c(:,:,i_wf_p))
         !  enddo ! ibp
         !enddo ! ibk
         !
         call WF_free(WF_p_tmp)
         !
         ! DEBUG <
         !call SERIAL_inversion(nb,WF_phases(:,:,s_count,k_count,i_sp_pol),det,.true.)
         !!
         !write(*,*) i_sp_pol,ikbz,isp,iGp,abs(det)
         ! DEBUG >
         !
       enddo ! isp
       !
       call WF_free(WF_p)
       !
       call WF_free(WF_k)
       !
       !call PP_redux_wait(WF_phases(:,:,s_count,k_count,i_sp_pol))
       !
       ! DEBUG <
       !write(*,*) "Character table at k= ",ikbz
       !write(*,*) "Number of symmetry elements= ",s_count
       !do i_grp=1,n_deg_grp
       !  write(*,*) En%E(first_el(i_grp),ik,i_sp_pol)*HA2EV,n_of_el(i_grp),real(charac(i_grp,:s_count))
       !enddo
       !write(*,*) ""
       ! DEBUG >
       !
       call live_timing(steps=1)
       !
       call io_control(ACTION=WR,COM=REP,SEC=(/2/),MODE=DUMP,ID=ID)
       io_err=io_WF_phases(blim,ikbz,i_sp_pol,ID,WF_phases)
       !
       YAMBO_FREE(WF_phases)
       !
     enddo ! istar
   enddo   ! ik
   !
 enddo ! i_sp_pol
 !
 call live_timing()
 !
 call io_control(ACTION=WR_CL,COM=REP,SEC=(/3/),MODE=DUMP,ID=ID)
 io_err=io_WF_phases(blim,ikbz,i_sp_pol,ID)
 !
 YAMBO_FREE(WF_phases_s_map)
 YAMBO_FREE(WF_phases_b_map)
 !
 call timing("WF_phases",OPR="stop")
 !
end subroutine WF_phase_matrices
