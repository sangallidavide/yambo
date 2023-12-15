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
subroutine WF_phase_matrices(En,Xk,nb)
 !    
 use pars,              ONLY:SP,cZERO
 use LIVE_t,            ONLY:live_timing
 use R_lattice,         ONLY:bz_samp,G_m_G,minus_G,g_vec,rl_sop,g_rot
 use D_lattice,         ONLY:nsym
 use electrons,         ONLY:levels
 use wave_func,         ONLY:WAVEs,WF_phases,WF_phases_k_map,WF_phases_s_map
 use wave_func,         ONLY:wf_ng_overlaps,wf_ng_1st_BZ,WF_copy
 use electrons,         ONLY:n_spinor,n_sp_pol
 use wrapper,           ONLY:Vstar_dot_V
 use hamiltonian,       ONLY:B_mat_index
 use matrix_operate,    ONLY:m3det!,UP_matrix_index
 use parallel_m,        ONLY:PAR_IND_WF_k,PAR_IND_WF_b,PAR_IND_B_mat
 use parallel_int,      ONLY:PP_redux_wait,PARALLEL_WF_index,PARALLEL_WF_distribute
 use interfaces,        ONLY:WF_symm_kpoint_g,WF_load,WF_free,eval_G_minus_G
 use IO_int,            ONLY:io_control
 use IO_m,              ONLY:VERIFY,DUMP,REP,OP_WR_CL,OP_RD_CL,OP_APP_CL
 use timing_m,          ONLY:timing
 !
#include<memory.h>
 !
 integer,       intent(in) :: nb(2)
 type(levels),  intent(in) :: En
 type(bz_samp), intent(in) :: Xk 
 !
 ! Work Space
 !
 integer           :: ikbz,ipbz,isk,isp,ibk,ibp,i_sp_pol,ik,ip,&
 &                    k_count,s_count,i_wf_k,i_wf_p,iGp,ng_tmp
 complex(SP)       :: det
 !
 type(WAVEs) :: WF_k,WF_p,WF_p_tmp
 !
 integer              :: io_err,ID
 integer, external    :: io_WF_phases
 !
 if (nsym==1) return
 !
 call timing("WF_phases",OPR="start")
 !
 call k_build_up_BZ_tables(Xk)
 !
 iGp=15
 ng_tmp=eval_G_minus_G(wf_ng_1st_BZ,iGp)
 !
 !call PARALLEL_WF_distribute(K_index=PAR_IND_WF_k,B_index=PAR_IND_WF_b,CLEAN_UP=.TRUE.)
 call PARALLEL_WF_distribute(CLEAN_UP=.TRUE.)
 call PARALLEL_WF_index( )
 !
 YAMBO_ALLOC(WF_phases,(nb(1):nb(2),nb(1):nb(2),nsym,Xk%nbz,n_sp_pol))
 YAMBO_ALLOC(WF_phases_k_map,(Xk%nbz))
 YAMBO_ALLOC(WF_phases_s_map,(nsym))
 !
 call io_control(ACTION=OP_RD_CL,COM=REP,SEC=(/1/),MODE=VERIFY,ID=ID)
 io_err=io_WF_phases(nb,ID)
 !
 if (io_err==0) then
   call io_control(ACTION=OP_RD_CL,COM=REP,SEC=(/2/),MODE=DUMP,ID=ID)
   io_err=io_WF_phases(nb,ID)
   return
 endif
 !
 call io_control(ACTION=OP_WR_CL,COM=REP,SEC=(/1/),MODE=DUMP,ID=ID)
 io_err=io_WF_phases(nb,ID)
 !
 WF_phases=cZERO
 WF_phases_k_map=0
 WF_phases_s_map=0
 !
 !call live_timing("WF_phases:",PAR_IND_WF_k%n_of_elements(PAR_IND_DIPk_bz_ID+1)*n_sp_pol)
 call live_timing("WF_phases:",Xk%nbz*n_sp_pol)
 !
 do i_sp_pol=1,n_sp_pol
   !
   k_count=0
   !
   do ikbz=1,Xk%nbz
     !
     !if ( .not.PAR_IND_WF_k%element_1D(ik) ) cycle
     !
     ik = Xk%sstar(ikbz,1)
     isk= Xk%sstar(ikbz,2)
     !
     k_count=k_count+1
     !
     WF_phases_k_map(ikbz)=k_count
     !
     call WF_load(WF_k,0,1,nb,(/ik,ik/),(/i_sp_pol,i_sp_pol/),&
     & space='G',title='-WF_phases',quiet=k_count>1,keep_states_to_load=.true.)
     !
     ip=ik
     !
     call WF_load(WF_p,0,1,nb,(/ik,ik/),(/i_sp_pol,i_sp_pol/),&
     & space='G',title='-WF_phases',quiet=.true.,   keep_states_to_load=.true.)
     !
     call WF_symm_kpoint_g(nb,ik,i_sp_pol,isk,1,WF_k)
     !
     s_count=0
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
       WF_phases_s_map(isp)=s_count
       !
       call WF_copy(WF_p,WF_p_tmp)
       !
       call WF_symm_kpoint_g(nb,ip,i_sp_pol,isp,iGp,WF_p_tmp)
       !
       do ibk=nb(1),nb(2)
         i_wf_k=WF_k%index(ibk,ik,i_sp_pol)
         do ibp=nb(1),nb(2)
           i_wf_p=WF_p%index(ibp,ip,i_sp_pol)
           !if (.not.PAR_IND_B_mat_ordered%element_1D( UP_matrix_index(ib,ibp)-1 )) cycle
           !if (.not.PAR_IND_B_mat%element_1D( B_mat_index(ib,ibp,(/1,En%nb/))-1 )) cycle
           WF_phases(ibk,ibp,s_count,k_count,i_sp_pol)= &
           & Vstar_dot_V(wf_ng_1st_BZ*n_spinor,WF_k%c(:,:,i_wf_k),WF_p_tmp%c(:,:,i_wf_p))
         enddo ! ibp
       enddo ! ibk
       !
       call WF_free(WF_p_tmp)
       !
       ! DEBUG <
       !call SERIAL_inversion(nb(2)-nb(1)+1,WF_phases(:,:,s_count,k_count,i_sp_pol),det,.true.)
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
     call live_timing(steps=1)
     !
   enddo ! ikbz
   !
 enddo ! i_sp_pol
 !
 call io_control(ACTION=OP_APP_CL,COM=REP,SEC=(/2/),MODE=DUMP,ID=ID)
 io_err=io_WF_phases(nb,ID)
 !
 call live_timing()
 !
 call timing("WF_phases",OPR="stop")
 !
end subroutine WF_phase_matrices
