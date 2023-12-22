! 
! License-Identifier: GPL
!
! Copyright (C) 2023 The Yambo Team
!
! Authors (see AUTHORS file for details): DS
!
subroutine electrons_symmetries(Xk,Xen)
 !
 use pars,           ONLY:DP,SP,cZERO
 use units,          ONLY:HA2EV
 use electrons,      ONLY:levels,n_sp_pol
 use R_lattice,      ONLY:bz_samp
 use D_lattice,      ONLY:elem,which_irr,char_mat,nclass,elem,nelem,&
&                         i_space_inv,i_time_rev,nsym,nrot,dl_sop,sop_tab
 use vec_operate,    ONLY:degeneration_finder
 use wave_func,      ONLY:WF_phases,WF_phases_s_map
 use IO_int,         ONLY:io_control
 use IO_m,           ONLY:VERIFY,DUMP,REP,OP_RD,RD,RD_CL
 !
#include<memory.h>
 !
 type(levels) :: Xen
 type(bz_samp):: Xk
 !
 integer     :: blim(2)
 complex(SP) :: trace(12)
 !
 integer     :: iclass,irap,i_sp_pol,ikbz,ipbz,ik,isk,ibk,istark,blim_n(2), &
 &              ielem,isp,k_count,s_count,n_deg_grp,i_grp,nb,nrotk
 integer, allocatable :: first_el(:),n_of_el(:)
 !
 integer              :: io_err,ID
 integer, external    :: io_WF_phases
 !
 blim=(/1,Xen%nb/)
 !blim=(/1,5/)
 nb=blim(2)-blim(1)+1
 !
 allocate(first_el(nb))
 allocate( n_of_el(nb))
 !
 call k_build_up_BZ_tables(Xk)
 call k_ibz2bz(Xk,'iku',.true.)
 call k_small_group(Xk)
 !
 call WF_phase_matrices(Xen,Xk,blim)
 !
 YAMBO_ALLOC(WF_phases_s_map,(Xk%nbz,nsym))
 !AMBO_ALLOC(WF_phases_b_map,(blim(1):blim(2),Xk%nbz,nsym,2))
 !
 call io_control(ACTION=OP_RD,COM=REP,SEC=(/1,3/),MODE=DUMP,ID=ID)
 io_err=io_WF_phases(blim,0,0,ID)
 !
 do i_sp_pol=1,n_sp_pol
   !
   do ik=1,Xk%nibz
     !
     call degeneration_finder(Xen%E(blim(1):blim(2),ik,i_sp_pol),nb,&
     &                           first_el,n_of_el,n_deg_grp,0.0001/HA2EV,Include_single_values=.true.)
     !
     do istark=1,Xk%nstar(ik)
       !
       !if ( .not.PAR_IND_WF_k%element_1D(ik) ) cycle
       !
       isk =Xk%star(ik,istark)
       ikbz=Xk%k_table(ik,isk)
       !
       nrotk=Xk%grp_nsym(ikbz)
       if ( any( Xk%grp_table(:nrotk,ikbz)>nrot ) ) nrotk=nrotk/(1+i_time_rev)
       call detect_point_group(nrotk,real(dl_sop(:,:,Xk%grp_table(:nrotk,ikbz)),DP),.true.)
       !
       YAMBO_ALLOC(WF_phases,(maxval(n_of_el),maxval(n_of_el),n_deg_grp,nsym,1,1))
       !
       call io_control(ACTION=RD,COM=REP,SEC=(/2/),MODE=DUMP,ID=ID)
       io_err=io_WF_phases(blim,ikbz,i_sp_pol,ID)
       !
       do i_grp=1,n_deg_grp
         !
         trace=cZERO
         do iclass=1,nclass
           do ielem=1,nelem(iclass)
             isp=Xk%grp_table(elem(ielem,iclass),ikbz)
             ipbz=Xk%k_table(ik,sop_tab(isp,isk))
             s_count=WF_phases_s_map(ikbz,sop_tab(isp,isk))
             do ibk=1,n_of_el(i_grp)
               trace(iclass)=trace(iclass)+WF_phases(ibk,ibk,i_grp,s_count,1,1)
             enddo
	   enddo
	   !
         enddo
         !write(*,*) "checking ",ikbz,i_grp,n_of_el(i_grp)
         !write(*,*) Xen%E(first_el(i_grp),ik,i_sp_pol)*HA2EV,real(trace(:nclass),SP)
         !write(*,*) "checking iraps "
         do irap=1,nclass
           if ( .not. all(abs(trace(:nclass)-char_mat(irap,which_irr(:nclass)))<0.0001_SP) ) cycle
           write(*,*) "trace= ",real(char_mat(irap,which_irr(:nclass)),SP)
           write(*,*) i_grp," corresponds to ",irap
         enddo
       enddo
       !
       YAMBO_FREE(WF_phases)
       !
     enddo
   enddo 
 enddo
 !
 call io_control(ACTION=RD_CL,COM=REP,SEC=(/0/),MODE=DUMP,ID=ID)
 io_err=io_WF_phases(blim,0,0,ID)
 !
end subroutine electrons_symmetries
