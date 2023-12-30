! 
! License-Identifier: GPL
!
! Copyright (C) 2023 The Yambo Team
!
! Authors (see AUTHORS file for details): DS
!
subroutine excitons_symmetries(iq,Xen,Xk,qpt,BS_E)
 !
 use pars,           ONLY:SP,DP,cZERO,schlen
 use units,          ONLY:HA2EV
 use R_lattice,      ONLY:bz_samp
 use electrons,      ONLY:levels
 use BS_solvers,     ONLY:BS_mat,BSS_n_eig
 use BS,             ONLY:BS_bands,BS_H_dim
 use wave_func,      ONLY:io_WF_phases,WF_phases_s_map,WF_phases_b_map
 use stderr,         ONLY:intc
 use D_lattice,      ONLY:which_irr,char_mat,nclass,elem,nelem,&
 &                        name_rap,i_time_rev,nsym,nrot,dl_sop
 use vec_operate,    ONLY:sort,degeneration_finder
 use wrapper,        ONLY:Vstar_dot_V
 use IO_int,         ONLY:io_control
 use IO_m,           ONLY:DUMP,REP,OP_RD,RD_CL
 use com,            ONLY:msg,of_open_close
 !
#include<memory.h>
 !
 integer, intent(in) :: iq
 type(bz_samp)       :: Xk,qpt
 type(levels)        :: Xen
 complex(SP)         :: BS_E(BSS_n_eig)
 !
 !real(SP) :: BS_E_real(EXCITONS_n_user_states)
 !integer  :: i_c,i_l,first_exc(EXCITONS_n_user_states),n_of_exc(EXCITONS_n_user_states),n_of_exc_max
 real(SP) :: BS_E_real(BSS_n_eig)
 integer  :: S_index(BSS_n_eig)
 integer  :: first_exc(BSS_n_eig),n_of_exc(BSS_n_eig),n_deg_grp,n_of_exc_max
 integer  :: iqbz_in,iqbz,isq,isq_small_grp,isqy,iGq,nrotq,irap,iclass,ielem,il,i_grp,&
 &           io_err,ID
 !
 complex(SP), allocatable :: EXC_WF_rot1(:,:),EXC_WF_rot2(:,:)
 complex(SP) :: trace(12,BSS_n_eig)
 !
 character(schlen)    :: filename,titles(3),tmp_ch
 !
 call k_build_up_BZ_tables(Xk)
 call k_ibz2bz(Xk,'iku',.true.)
 call k_small_group(Xk)
 !
 call WF_phase_matrices(Xen,Xk,BS_bands)
 !
 call k_small_group(qpt)
 !
 ! The q point for which I want to compute the irrep.
 !isq =q%star(iq,1)
 iqbz_in=qpt%k_table(iq,1)
 !iGq =qpt%g_table(iq,isq)
 !
 nrotq=qpt%grp_nsym(iqbz)
 if ( any( qpt%grp_table(:nrotq,iqbz)>nrot ) ) nrotq=nrotq/(1+i_time_rev)
 call detect_point_group(nrotq,real(dl_sop(:,:,qpt%grp_table(:nrotq,iqbz)),DP),.true.)
 !
 YAMBO_ALLOC(WF_phases_s_map,(Xk%nbz,nsym))
 YAMBO_ALLOC(WF_phases_b_map,(BS_bands(1):BS_bands(2),Xk%nbz,nsym,2))
 !
 call io_control(ACTION=OP_RD,COM=REP,SEC=(/1,3/),MODE=DUMP,ID=ID)
 io_err=io_WF_phases(BS_bands,0,0,ID)
 !
 !do i_c=1,EXCITONS_n_user_states
 !  i_l=EXCITONS_user_indexes(i_c)
 !  BS_E_real(i_c)=real(BS_E(i_l))
 !enddo
 !
 BS_E_real=real(BS_E)
 call sort(BS_E_real,indx=S_index)
 !
 call degeneration_finder(BS_E_real,BSS_n_eig,first_exc,n_of_exc,n_deg_grp,&
 &                           0.0001/HA2EV,Include_single_values=.true.)
 !
 !n_of_exc_max=maxval(n_of_exc)
 n_of_exc_max=3
 !
 YAMBO_ALLOC(EXC_WF_rot1,(BS_H_dim,n_of_exc_max))
 YAMBO_ALLOC(EXC_WF_rot2,(BS_H_dim,n_of_exc_max))
 !
 filename="exc_syms_q"//trim(intc(iq))
 !
 call of_open_close(filename,'ot')
 !call msg("o "//trim(filename),"# Symmetries at k ",Xk%ptbz(ikbz,:))
 call msg("o "//trim(filename),"#"," Symmetries at iq "//trim(intc(iq)),INDENT=0)
 call msg("o "//trim(filename),"#")
 !
 nrotq=qpt%grp_nsym(iqbz)
 if ( any( qpt%grp_table(:nrotq,iqbz)>nrot ) ) nrotq=nrotq/(1+i_time_rev)
 call detect_point_group(nrotq,real(dl_sop(:,:,qpt%grp_table(:nrotq,iqbz)),DP),.true.)
 !
 call print_point_group("o "//trim(filename))
 !
 titles=(/" E[eV]  "," symm.  ","i_lambda"/)
 call msg("o "//trim(filename),"#")
 call msg("o "//trim(filename),'#',titles,INDENT=0,USE_TABS=.true.)
 call msg("o "//trim(filename),"#")
 !
 trace=cZERO
 !
 do iclass=1,nclass
   do ielem=1,nelem(iclass)
     !
     isq_small_grp=elem(ielem,iclass)
     !
     ! The rotated q-point
     isq =qpt%grp_table(isq_small_grp,iqbz_in) ! the symmetry of the loop
     isqy=qpt%s_table(iq,isq)                    ! the symmetry used by yambo
     !
     iqbz=qpt%k_table(iq,isq)                    ! the rotated q-point
     iGq =qpt%k_table(iq,isq)                    ! the G shift
     !
     do i_grp=1,n_deg_grp
       !
       ! Get A(l,R*q)= Ry * A(l,q)
       !
       call EXC_WF_symm_qpoint(Xk,first_exc(i_grp),n_of_exc(i_grp),&
       & iq,isqy,iGq,ID,S_index,EXC_WF_rot1(:,:n_of_exc(i_grp)))
       !
       ! Get R * A(l,q)
       !
       call EXC_WF_symm_qpoint(Xk,first_exc(i_grp),n_of_exc(i_grp),&
       & iq,isq,iGq,ID,S_index,EXC_WF_rot2(:,:n_of_exc(i_grp)))
       !
       ! Compute the trace of the Irrep
       !
       do il=1,n_of_exc(i_grp)
         trace(iclass,i_grp) = trace(iclass,i_grp) +  &
         & Vstar_dot_V(BS_H_dim,EXC_WF_rot1(:,il),EXC_WF_rot2(:,il))
       enddo
       trace(iclass,i_grp)=trace(iclass,i_grp)/real(nelem(iclass),SP)
       !
     enddo
   enddo
 enddo
 !
 call io_control(ACTION=RD_CL,COM=REP,SEC=(/0/),MODE=DUMP,ID=ID)
 io_err=io_WF_phases(BS_bands,0,0,ID)
 !
 ! Finally I compute the exciton symmetry
 !
 do i_grp=1,n_deg_grp
   do irap=1,nclass
     if ( .not. all(abs(trace(:nclass,i_grp)-char_mat(irap,which_irr(:nclass)))<0.0001_SP) ) cycle
     write(tmp_ch,"(F10.4,12x,a,12x,i4)") &
     & BS_E_real(first_exc(i_grp))*HA2EV,trim(name_rap(irap)),first_exc(i_grp)
     call msg("o "//trim(filename),'',tmp_ch,INDENT=-2)
   enddo
 enddo
 !
end subroutine excitons_symmetries
