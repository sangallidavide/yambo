!
! License-Identifier: GPL
!
! Copyright (C) 2006 The Yambo Team
!
! Authors (see AUTHORS file for details): AM DS
!
subroutine QP_apply_interpolate_QP(i_ACTION,bands,en,k,qp)
 !
 ! Extend the QP correction from the kpt in qp to the kpt in k
 ! Requires qp%nk, qp%k, qp%table, qp%n_states, qp%E, qp%E_bare, qp%Z
 !
 ! QP_action(i_ACTION)%plan=( E W Z Interp/k Interp/E Interp/BOLTZ )
 !
 ! Note that as output of this routine W is positive(negative) above(below) the Fermi level. 
 ! This is true both at zero and finite temperature. Indeed also at finite temperature
 ! when the theory is causal the GF must be rotated to the real-energy axis.
 !
 use pars,          ONLY:SP,schlen,cZERO
 use units,         ONLY:HA2EV
 use IO_int,        ONLY:io_control
 use IO_m,          ONLY:OP_WR_CL,LOG
 use electrons,     ONLY:n_sp_pol,spin
 use com,           ONLY:msg
 use electrons,     ONLY:levels,spin
 use R_lattice,     ONLY:bz_samp,bz_samp_reset
 use stderr,        ONLY:real2ch
 use QP_m,          ONLY:QP_t,QP_alloc,QP_reset
 use QP_CTL_m,      ONLY:QP_action,QP_ctl_msg_fmt,QP_ctl_xc_STRING,At_least_one_band_corrected
 use interpolate,   ONLY:NN_n_of_nearest,INTERP_obj,INTERPOLATE_msg_fmt,&
&                        INTERPOLATION_driver_seed,INTERPOLATION_driver_do,INTERP_shell_factor
 !
#include<y_memory.h>
 !
 type(levels)  ::en
 type(QP_t)    ::qp
 type(bz_samp) ::k
 integer       ::i_ACTION,bands(2)
 !
 ! Work Space
 !
 real(SP)         :: E_diff(n_sp_pol),E_min(n_sp_pol)
 integer          :: ik,ik_fg,iband,i_spin,ID,i_qp,N_exact,N_close,N_far,DB_bands(2),ID_IO,i_err,ic
 logical          :: wrong_width_sign_warn
 character(schlen):: what,QP_interpolated_name
 type(bz_samp)    :: QP_k
 type(QP_t)       :: qp_interpolated
 !
 integer,    external    :: io_QP_and_GF
 !
 ! IBZ -> BZ (iku)
 !
 call k_ibz2bz(k,'i',.true.)
 !
 ! Transfer & interpolation
 !
 wrong_width_sign_warn=.false.
 !
 DB_bands=(/minval(qp%table(:,1)),maxval(qp%table(:,1))/)
 !
 call QP_reset(qp_interpolated)
 qp_interpolated%nk=k%nibz
 qp_interpolated%n_states=(DB_bands(2)-DB_bands(1)+1)*n_sp_pol*k%nibz
 call QP_alloc(qp_interpolated)
 ! 
 ! Setup
 !
 NN_n_of_nearest=QP_action(i_ACTION)%CTL(1)%interp_neigh
 INTERP_shell_factor=QP_action(i_ACTION)%CTL(1)%interp_shell_factor
 INTERPOLATE_msg_fmt=QP_ctl_msg_fmt
 !
 if (QP_action(i_ACTION)%plan(4).or.QP_action(i_ACTION)%plan(6)) then
   !
   ! Nearest K-point interpolation
   !===============================
   !
   call QP_transfer('init')
   !
   call k_expand(QP_k)
   !
   call INTERPOLATION_driver_seed(trim(what)//"_from_DB",ID,bands,QP_k,SEED_QP=qp,NEAREST=QP_action(i_ACTION)%plan(4))
   call INTERPOLATION_driver_do(trim(what)//"_from_DB",ID,k,FineGrid="BARE")
   !
   call QP_transfer('appl')
   !
   call INTERPOLATION_driver_end(ID)
   !
 else if (QP_action(i_ACTION)%plan(5)) then
   !
   ! Nearest Level interpolation
   !===============================
   !
   E_diff=1.E5_SP
   do i_spin=1,n_sp_pol
     E_min(i_spin) =minval(En%Eo(DB_bands(1),:,i_spin))
     do i_qp=1,qp%n_states
       if (spin(qp%table(i_qp,:))/=i_spin) cycle
       if (qp%E_bare(i_qp)<E_diff(i_spin)) E_diff(i_spin)=qp%E_bare(i_qp)
     enddo
     E_diff(i_spin)=E_diff(i_spin)-E_min(i_spin)
     if (E_diff(i_spin)>1.E-5_SP) then
       call warning(' Energy not aligned. Difference is '//trim(real2ch(E_diff(i_spin)*HA2EV))//' eV. QP interpolation maybe wrong')
     endif
   enddo
   !
   N_exact=0
   N_close=0
   N_far  =0
   !
   ic=0
   do ik=1,k%nibz
    qp_interpolated%k(ik,:)=k%pt(ik,:)
    do iband=bands(1),bands(2)
       do i_spin=1,n_sp_pol
	 !
         if (DB_bands(1)<=iband .and. iband<=DB_bands(2)) ic=ic+1
         !
         call INTERPOLATION_nearest_E(En%Eo(iband,ik,i_spin),qp,E_diff(i_spin),i_qp,i_spin)
         call QP_transfer_nearest('Stnd_Gd')
         !
         if(allocated(en%FG%E)) then
           do ik_fg=k%FGbare%k_range(ik,1),k%FGbare%k_range(ik,2)
             call INTERPOLATION_nearest_E(En%FG%E(iband,ik_fg,i_spin),qp,E_diff(i_spin),i_qp,i_spin)
             call QP_transfer_nearest('Fine_Gd')
           enddo
         endif
         !
       enddo
     enddo
   enddo
   !
   call msg("n"//trim(QP_ctl_msg_fmt),'[QP_ctl (Nearest Level)] Levels Mapped exactly ',N_exact)
   call msg(trim(QP_ctl_msg_fmt),     '[QP_ctl (Nearest Level)]               closely ',N_close)
   call msg(trim(QP_ctl_msg_fmt),     '[QP_ctl (Nearest Level)] Levels too far        ',N_far)
   !
 endif
 !
 if (index(QP_ctl_msg_fmt,"r")>0) call msg("rn","")
 !
 ! XC kind update
 !
 if (any((/QP_action(i_ACTION)%STATES_done==1/))) then
   if (len_trim(QP_ctl_xc_STRING)>0) then
     QP_ctl_xc_STRING=trim(QP_ctl_xc_STRING)//" + "//trim(QP_action(i_ACTION)%what)
   else
     QP_ctl_xc_STRING=trim(QP_action(i_ACTION)%what)
   endif
 endif
 !
 ! Output
 !--------
 !
 QP_interpolated_name='QP_interpolated'
 call io_control(ACTION=OP_WR_CL,COM=LOG,SEC=(/1,2,3/),ID=ID_IO)
 i_err=io_QP_and_GF(trim(QP_interpolated_name),qp_interpolated,ID_IO)
 call QP_reset(qp_interpolated)
 !
 ! Clean
 !
 call k_ibz2bz(k,'d',.false.)
 !
 contains
   !
   subroutine QP_transfer(action_kind)
     !
     character(4) :: action_kind
     !
     real(SP)     :: wsign
     integer      :: i_ref(3),ik1,ik2
     !
     select case(action_kind)
       !
       case('init')
         !
         call bz_samp_reset(QP_k)
         QP_k%nibz=qp%nk
         YAMBO_ALLOC(QP_k%pt,(qp%nk,3))
         YAMBO_ALLOC(QP_k%weights,(qp%nk))
         QP_k%pt=qp%k
         QP_k%weights=qp%wk
         !
         what=" "
         if (QP_action(i_ACTION)%plan(1)) what="dE"//trim(what)
         if (QP_action(i_ACTION)%plan(2)) what="W"//trim(what)
         if (QP_action(i_ACTION)%plan(3)) what="Z"//trim(what)
         !
       case('appl')
         !
         ic=0
         do ik=1,k%nibz
          qp_interpolated%k(ik,:)=k%pt(ik,:)
          do iband=bands(1),bands(2)
             !
             if (iband<DB_bands(1)) cycle
             if (iband>DB_bands(2)) cycle
             !
             do i_spin=1,n_sp_pol
               ic=ic+1
               !
               i_ref=0
               if (QP_action(i_ACTION)%plan(1)) i_ref(2)=n_sp_pol
               if (QP_action(i_ACTION)%plan(2)) i_ref(3)=i_ref(2)+n_sp_pol
               !
               if(QP_action(i_ACTION)%plan(1)) then
                 en%E(iband,ik,i_spin)=en%Eo(iband,ik,i_spin)+INTERP_obj(ID)%INTERP(iband,ik,i_ref(1)+i_spin)
                 !
                 qp_interpolated%E_bare(ic)=en%Eo(iband,ik,i_spin)
                 qp_interpolated%E(ic)=en%E(iband,ik,i_spin)
                 qp_interpolated%Z(ic)=cZERO
                 qp_interpolated%table(ic,1)=iband
                 qp_interpolated%table(ic,2)=iband
                 qp_interpolated%table(ic,3)=ik
                 if(n_sp_pol==2) qp_interpolated%table(ic,4)=i_spin
                 !
                 if(allocated(en%FG%E)) then
                   ik1=k%FGbare%k_range(ik,1)
                   ik2=k%FGbare%k_range(ik,2)
                   en%FG%E(iband,ik1:ik2,i_spin)=en%FG%E(iband,ik1:ik2,i_spin)+&
&                                                    INTERP_obj(ID)%INTERP_FineGd(iband,ik1:ik2,i_ref(1)+i_spin)
                 endif
                 At_least_one_band_corrected=.TRUE.
                 QP_action(i_ACTION)%STATES_done(iband,ik,i_spin)=1
               endif
               if(QP_action(i_ACTION)%plan(2)) then
                 if (en%Eo(iband,ik,i_spin)>0. .and. INTERP_obj(ID)%INTERP(iband,ik,i_ref(1)+i_spin)>0) then
                   if (.not.wrong_width_sign_warn) call warning(' Wrong QP width sign fixed')
                   wsign=-1._SP
                   wrong_width_sign_warn=.true.
                 else
                   wsign=+1._SP
                 endif
                 en%W(iband,ik,i_spin)=wsign*INTERP_obj(ID)%INTERP(iband,ik,i_ref(2)+i_spin)
                 if(allocated(en%FG%W)) then
                   ik1=k%FGbare%k_range(ik,1)
                   ik2=k%FGbare%k_range(ik,2)
                   en%FG%W(iband,ik1:ik2,i_spin)=wsign*INTERP_obj(ID)%INTERP_FineGd(iband,ik1:ik2,i_ref(1)+i_spin)
                 endif
                 QP_action(i_ACTION)%STATES_done(iband,ik,i_spin)=1
               endif
               if(QP_action(i_ACTION)%plan(3)) then
                 en%Z(iband,ik,i_spin)=cmplx(INTERP_obj(ID)%INTERP(iband,ik,i_ref(3)+i_spin),&
&                                            INTERP_obj(ID)%INTERP(iband,ik,i_ref(3)+n_sp_pol+i_spin),kind=SP)
                 QP_action(i_ACTION)%STATES_done(iband,ik,i_spin)=1
               endif
             enddo
           enddo
         enddo
         !
     end select
     !
   end subroutine
   !
   subroutine QP_transfer_nearest(action_kind)
     !
     character(7) :: action_kind
     real(SP)     :: wsign
     !
     if (                                E_diff(i_spin)<1.E-5_SP    ) N_exact=N_exact+1
     if (E_diff(i_spin)>1.E-5_SP    .and.E_diff(i_spin)<0.1_SP/HA2EV) N_close=N_close+1
     if (E_diff(i_spin)>0.1_SP/HA2EV                                ) N_far  =N_far  +1
     !
     if (QP_action(i_ACTION)%plan(1)) then
       if(action_kind=="Stnd_Gd") en%E(iband,ik,i_spin)=en%Eo(iband,ik,i_spin)+real(qp%E(i_qp),SP)-qp%E_bare(i_qp)
       if(action_kind=="Fine_Gd") en%FG%E(iband,ik_fg,i_spin)=en%FG%E(iband,ik_fg,i_spin)+ &
&                                                             real(qp%E(i_qp),SP)-qp%E_bare(i_qp)
       At_least_one_band_corrected=.TRUE.
       if(action_kind=="Stnd_Gd") QP_action(i_ACTION)%STATES_done(iband,ik,i_spin)=1
       !
       if(action_kind=="Stnd_Gd" .and. DB_bands(1)<=iband .and. iband<=DB_bands(2)) then
         qp_interpolated%E_bare(ic)=en%Eo(iband,ik,i_spin)
         qp_interpolated%E(ic)=en%E(iband,ik,i_spin)
         qp_interpolated%Z(ic)=cZERO
         qp_interpolated%table(ic,1)=iband
         qp_interpolated%table(ic,2)=iband
         qp_interpolated%table(ic,3)=ik
         if(n_sp_pol==2) qp_interpolated%table(ic,4)=i_spin
       endif
       !
     endif
     !
     if (QP_action(i_ACTION)%plan(2)) then
       if (real(qp%E_bare(i_qp))>0..and.aimag(qp%E(i_qp))>0) then
         if (.not.wrong_width_sign_warn) call warning(' Wrong QP width sign fixed')
         wsign=-1._SP
         wrong_width_sign_warn=.true.
       else
         wsign=+1._SP
       endif
       if(action_kind=="Stnd_Gd") en%W(iband,ik,i_spin)=wsign*aimag(qp%E(i_qp))
       if(action_kind=="Fine_Gd") en%FG%W(iband,ik_fg,i_spin)=wsign*aimag(qp%E(i_qp))
       if(action_kind=="Stnd_Gd") QP_action(i_ACTION)%STATES_done(iband,ik,i_spin)=1
     endif
     !
     if (QP_action(i_ACTION)%plan(3)) then
       if(action_kind=="Stnd_Gd") en%Z(iband,ik,i_spin)=qp%Z(i_qp)
       !DEBUG>
       !if(action_kind=="Fine_Gd") en%Z_FineGd(iband,ik_fg,i_spin)=qp%Z(i_qp)
       !DEBUG<
       if(action_kind=="Stnd_Gd") QP_action(i_ACTION)%STATES_done(iband,ik,i_spin)=1
     endif
     !
   end subroutine
   !
end subroutine
