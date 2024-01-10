! 
! License-Identifier: GPL
!
! Copyright (C) 2023 The Yambo Team
!
! Authors (see AUTHORS file for details): DS
!
subroutine excitons_symmetries(iq,Xen,Xk,qpt,BS_E)
 !
 use pars,           ONLY:SP,DP,cZERO,cI,schlen
 use units,          ONLY:HA2EV
 use YPPm,           ONLY:deg_energy,EXCITONS_n_user_states,&
 &                        EXCITONS_user_indexes,EXC_symm_thresh
 use R_lattice,      ONLY:bz_samp,coll_scattering
 use electrons,      ONLY:levels,n_sp_pol
 use BS_solvers,     ONLY:BS_mat,BSS_n_eig
 use BS,             ONLY:BS_bands,BS_H_dim,BS_K_dim
 use wave_func,      ONLY:io_WF_phases,WF_phases_b_map
 use stderr,         ONLY:intc,real2ch
 use D_lattice,      ONLY:which_irr,char_mat,nclass,elem,nelem,&
 &                        name_rap,i_time_rev,nsym,nrot,dl_sop
 use vec_operate,    ONLY:sort,degeneration_finder
 use wrapper,        ONLY:Vstar_dot_V
 use IO_int,         ONLY:io_control
 use IO_m,           ONLY:DUMP,REP,OP_RD_CL,OP_RD,RD_CL
 use com,            ONLY:msg,of_open_close
 use LIVE_t,         ONLY:live_timing
 !
#include<memory.h>
 !
 integer, intent(in) :: iq
 type(bz_samp)       :: Xk,qpt
 type(levels)        :: Xen
 complex(SP)         :: BS_E(BSS_n_eig)
 !
 real(SP) :: BS_E_real(EXCITONS_n_user_states),tmp_real,ref_val
 logical  :: l_mixed_states,l_partial_space,l_char_table_imag,l_char_not_found
 integer  :: S_index(EXCITONS_n_user_states)
 integer  :: first_exc(EXCITONS_n_user_states),n_of_exc(EXCITONS_n_user_states),n_deg_grp,n_of_exc_max
 integer  :: iqbz,iRqbz,isqc,isqy,isq_small_grp,iGqc,iGqy,nrotq,&
 &           irap,iclass,ielem,iexc,il,i_grp,io_err,ID,n_elem_used,max_deg
 !
 complex(SP), allocatable :: EXC_WF_rot1(:,:),EXC_WF_rot2(:,:)
 complex(SP) :: trace(12,EXCITONS_n_user_states)
 !
 character(schlen)    :: filename,titles(4),tmp_ch,loc_msg,loc_fmt
 !
 integer, external :: io_QINDX
 !
 ! Electrons part (pre-requisite)
 !
 call k_build_up_BZ_tables(Xk)
 call k_ibz2bz(Xk,'iku',.true.)
 call k_small_group(Xk)
 !
 call WF_phase_matrices(Xen,Xk,BS_bands,(/1,Xk%nibz/))
 !
 YAMBO_ALLOC(WF_phases_b_map,(BS_bands(1):BS_bands(2),Xk%nibz,n_sp_pol,2))
 !
 call io_control(ACTION=OP_RD,COM=REP,SEC=(/1,2/),MODE=DUMP,ID=ID)
 io_err=io_WF_phases(BS_bands,0,0,0,ID,(/0,0,0,0,0/))
 !
 call build_inverse_eh_table(Xk%nbz)
 !
 ! Exciton part
 !
 call k_build_up_BZ_tables(qpt)
 call k_ibz2bz(qpt,'iku',.true.)
 call k_small_group(qpt)
 !
 !coll_scattering=.true.
 !call io_control(ACTION=OP_RD_CL,COM=REP,SEC=(/6/),ID=ID)
 !io_err=io_QINDX(Xk,qpt,ID,'minus_q')
 !if (io_err/=0) then
 !  call msg('s',"Missing k/q scattering database")
 !  call error('Please (re)run the setup K_grids="C"')
 !endif
 !
 ! The q point for which I want to compute the irrep.
 !isq =q%star(iq,1)
 iqbz=qpt%k_table(iq,1)
 !iGq =qpt%g_table(iq,isq)
 !
 do iexc=1,EXCITONS_n_user_states
   il=EXCITONS_user_indexes(iexc)
   BS_E_real(iexc)=real(BS_E(il))
 enddo
 !
 call sort(BS_E_real,indx=S_index)
 !
 call degeneration_finder(BS_E_real,EXCITONS_n_user_states,first_exc,n_of_exc,n_deg_grp,&
 &                           deg_energy,Include_single_values=.true.)
 !
 n_of_exc_max=maxval(n_of_exc)
 !n_of_exc_max=3
 !
 YAMBO_ALLOC(EXC_WF_rot1,(BS_H_dim,n_of_exc_max))
 YAMBO_ALLOC(EXC_WF_rot2,(BS_H_dim,n_of_exc_max))
 !
 filename="exc_syms_q"//trim(intc(iq))
 !
 call of_open_close(filename,'ot')
 call msg("o "//trim(filename),"# Symmetries at qpt ",qpt%pt(iq,:),INDENT=0)
 call msg("o "//trim(filename),"#")
 !
 nrotq=qpt%grp_nsym(iqbz)
 if ( any( qpt%grp_table(:nrotq,iqbz)>nrot ) ) nrotq=nrotq/(1+i_time_rev)
 call detect_point_group(nrotq,real(dl_sop(:,:,qpt%grp_table(:nrotq,iqbz)),DP),.true.)
 !
 call print_point_group("o "//trim(filename),nrotq,dl_sop(:,:,qpt%grp_table(:nrotq,iqbz)))
 !
 titles=(/" E[eV]  "," n_deg  "," symm.  ","i_lambda"/)
 call msg("o "//trim(filename),"#")
 call msg("o "//trim(filename),'#',titles,INDENT=0,USE_TABS=.true.)
 call msg("o "//trim(filename),"#")
 !
 trace=cZERO
 !
 !call live_timing("Exciton symmetries:",sum(nelem(:))*n_deg_grp)
 call live_timing("Exciton symmetries:",nclass*n_deg_grp)
 !
 do iclass=1,nclass
   ! The following loop is useful for debugging only
   ! Any element of the class should be ok 
   !n_elem_used=nelem(iclass)
   n_elem_used=1
   do ielem=1,n_elem_used
     !
     isq_small_grp=elem(ielem,iclass)
     !
     ! The rotated q-point
     isqc=qpt%grp_table(isq_small_grp,iqbz)       ! the symmetry of the loop
     isqy=qpt%s_table(iq,isqc)                    ! the symmetry used by yambo
     !
     iRqbz=qpt%k_table(iq,isqc)                   ! the rotated q-point
     !
     iGqc=qpt%k_table(iq,isqc)                    ! the G shift
     iGqy=qpt%k_table(iq,isqy)                    ! the G shift
     !
     do i_grp=1,n_deg_grp
       !
       ! Get A(l,R*q)= Ry * A(l,q)
       !
       call EXC_WF_symm_qpoint(Xk,qpt,first_exc(i_grp),n_of_exc(i_grp),&
       & iq,isqy,iGqy,ID,S_index,EXC_WF_rot1(:,:n_of_exc(i_grp)))
       !
       ! Get R * A(l,q)
       !
       call EXC_WF_symm_qpoint(Xk,qpt,first_exc(i_grp),n_of_exc(i_grp),&
       & iq,isqc,iGqc,ID,S_index,EXC_WF_rot2(:,:n_of_exc(i_grp)))
       !
       ! Compute the trace of the Irrep
       !
       do il=1,n_of_exc(i_grp)
         trace(iclass,i_grp) = trace(iclass,i_grp) +  &
         & Vstar_dot_V(BS_H_dim,EXC_WF_rot1(:,il),EXC_WF_rot2(:,il))
       enddo
       !
       call live_timing(steps=1)
       !
     enddo
   enddo
   trace(iclass,:)=trace(iclass,:)/real(n_elem_used,SP)
 enddo
 !
 call live_timing()
 !
 call io_control(ACTION=RD_CL,COM=REP,SEC=(/0/),MODE=DUMP,ID=ID)
 io_err=io_WF_phases(BS_bands,0,0,0,ID,(/0,0,0,0,0/))
 !
 ! Finally I compute the exciton symmetry
 !
 loc_fmt="(F10.4,12x,F4.2,12x,a,12x,i4)"
 l_char_table_imag=maxval(abs(real(cI*char_mat)))>0.01_SP
 !
 max_deg=nint(maxval(real(abs(char_mat(:,1)))))
 if (l_char_table_imag) max_deg=max_deg*2
 !
 do i_grp=1,n_deg_grp
   !
   il=first_exc(i_grp)
   ! DEBUG <
   !write(*,*) BS_E_real(il)*HA2EV,i_grp,n_of_exc(i_grp),real(trace(:nclass,i_grp))
   ! DEBUG >
   !
   l_mixed_states = nint(abs(trace(1,i_grp))) > max_deg 
   if (l_mixed_states) then
     tmp_ch="mixed states"
     l_char_not_found=.true.
   endif
   !
   l_partial_space=.false.
   do iclass=1,nclass
     tmp_real=abs(trace(iclass,i_grp))
     tmp_real=abs(tmp_real-real(nint(tmp_real),SP))
     l_partial_space = tmp_real > EXC_symm_thresh
     if(l_partial_space) exit
   enddo
   if (l_partial_space) then
     tmp_ch="partial space "//trim(real2ch(tmp_real))
     l_char_not_found=.true.
   endif
   !
   if ( .not. (l_mixed_states.or.l_partial_space) ) then
     do irap=1,nclass
       !if ( .not. all(abs(trace(:nclass,i_grp)-char_mat(irap,which_irr(:nclass)))<=EXC_symm_thresh) ) cycle
       if ( .not.( all(abs(real(trace(:nclass,i_grp))      -real(char_mat(irap,which_irr(:nclass))))<=EXC_symm_thresh) .or. &
       &           all(abs(real(trace(:nclass,i_grp))/2._SP-real(char_mat(irap,which_irr(:nclass))))<=EXC_symm_thresh) )) cycle
       tmp_ch=trim(name_rap(irap))
       l_char_not_found=.false.
       exit
     enddo
     if(irap==nclass+1) then
       tmp_ch="char not found"
       l_char_not_found=.true.
     endif
   endif
   !
   write(loc_msg,loc_fmt) BS_E_real(il)*HA2EV,real(trace(1,i_grp)),trim(tmp_ch),il
   call msg("o "//trim(filename),'',loc_msg,INDENT=-2)
   if (l_char_not_found) then
     do ielem=1,2
       loc_msg=""
       do iclass=1,nclass
         if(ielem==1) write(tmp_ch,"(2x,F5.2,2x)")  real(trace(iclass,i_grp))
         if(ielem==2) write(tmp_ch,"(2x,F5.2,2x)") aimag(trace(iclass,i_grp))
         loc_msg=trim(loc_msg)//trim(tmp_ch)
       enddo
       call msg("o "//trim(filename),'',loc_msg,INDENT=-2)
     enddo
   endif
   !
 enddo
 !
 ! DEBUG <
 !write(*,*) ""
 !do irap=1,nclass
 !  write(*,*) 0.0_SP,nelem(irap),real(char_mat(irap,which_irr(:nclass)),SP)
 !enddo
 ! DEBUG >
 !
 call of_open_close(filename)
 !
end subroutine excitons_symmetries
