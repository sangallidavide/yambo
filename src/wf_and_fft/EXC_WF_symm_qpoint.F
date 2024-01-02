!
! License-Identifier: GPL
!
! Copyright (C) 2015 The Yambo Team
!
! Authors (see AUTHORS file for details): DS
!
subroutine EXC_WF_symm_qpoint(Xk,qpt,first_state,nstates,iq,isq,igq,ID,S_index,EXC_WFs)
 !
 use pars,          ONLY:SP,cZERO,cONE
 use electrons,     ONLY:n_sp_pol
 use D_lattice,     ONLY:idt_index,i_time_rev,nsym,sop_inv,sop_tab
 use R_lattice,     ONLY:bz_samp,qindx_X,qindx_C
 use BS_solvers,    ONLY:BS_mat,BSS_n_eig,BSS_eh_table_m1
 use BS,            ONLY:BS_bands,BS_H_dim
 use wave_func,     ONLY:io_WF_phases,WF_phases_b_map
 use IO_int,        ONLY:io_control
 use IO_m,          ONLY:DUMP,REP,RD
 !
 implicit none
 !
 type(bz_samp), intent(in)   :: Xk,qpt
 integer,       intent(in)   :: first_state,nstates
 integer,       intent(in)   :: iq,isq,igq
 integer,       intent(inout):: ID
 complex(SP),   intent(out)  :: EXC_WFs(BS_H_dim,3)
 integer,       intent(in)   :: S_index(BSS_n_eig)
 !
 ! Work space
 !
 integer :: iexc,il,ik,istark,isk,ikbz,ipbz,ipmqbz,ipmqg0,iv1,ic1,iv2,ic2,&
 &          i_v1_bnd,i_v1_grp,i_c1_bnd,i_c1_grp,i_v2_bnd,i_v2_grp,i_c2_bnd,i_c2_grp,&
 &          ikmqbz,ikmqg0,iqbz,i_sp_pol,neh1,neh2,io_err,V_bands(2),C_bands(2)
 complex(SP) :: phase_vv,phase_cc
 complex(SP) :: WF_phases_cc(3,3,BS_bands(2)-BS_bands(1)+1,nsym,1,1),&
 &              WF_phases_vv(3,3,BS_bands(2)-BS_bands(1)+1,nsym,1,1)
 !
 if(isq==idt_index) then
   do iexc=1,nstates
     il=S_index(iexc+first_state-1)
     EXC_WFs(:,iexc)=BS_mat(:,il)
   enddo
   return
 endif
 !
 iqbz=qpt%k_table(iq,isq)
 !
 !V_bands=(/BS_bands(1),maxval(Ken%nbm)/)
 !C_bands=(/minval(Ken%nbf)+1,BS_bands(2)/)
 !
 EXC_WFs=cZERO
 !
 ! First get A(l,R*q)= Ry * A(l,q)
 !
 do iexc=1,nstates
   il=S_index(iexc+first_state-1)
   !
   do i_sp_pol=1,n_sp_pol
     !
     do ik=1,Xk%nibz
       !
       do istark=1,Xk%nstar(ik)
         !
         isk =Xk%star(ik,istark)
         ikbz=Xk%k_table(ik,isk)
         !
         ! p_BZ = (R)^-1 k_BZ = (R)^-1 (Rk) k_IBZ
         ipbz=Xk%k_table(ik,sop_tab(sop_inv(isq),isk))
         !
         ! For the phases I need:
         ! - R p_BZ = k_BZ
         ! - the symmetry used to reach k_BZ, e.g. isk
         !ik_ref=ikbz
         !is_ref=isk
         !
         call io_control(ACTION=RD,COM=REP,SEC=(/2/),MODE=DUMP,ID=ID)
         io_err=io_WF_phases(BS_bands,ipbz,i_sp_pol,ID,WF_phases_cc)
         !
         ! (R)^-1 k_BZ - q_IBZ = p_BZ - q_IBZ
         ipmqbz=qindx_X(iq,ipbz,1)
         ipmqg0=qindx_X(iq,ipbz,2)
         !
         ! For the phases I need R pmq_BZ = k_BZ - R q_IBZ = k_BZ - q_BZ
         !ikmqbz=qindx_C(ikbz,iqbz,1)
         !ikmqg0=qindx_C(ikbz,iqbz,2)
         !
         call io_control(ACTION=RD,COM=REP,SEC=(/2/),MODE=DUMP,ID=ID)
         io_err=io_WF_phases(BS_bands,ipmqbz,i_sp_pol,ID,WF_phases_vv)
         !
         do iv1 = BS_bands(1),BS_bands(2) 
           !
           i_v1_bnd=WF_phases_b_map(iv1,ipmqbz,isq,1)
           i_v1_grp=WF_phases_b_map(iv1,ipmqbz,isq,2)
           !
           do iv2 = BS_bands(1),BS_bands(2) 
             !
             i_v2_bnd=WF_phases_b_map(iv2,ipmqbz,isq,1)
             i_v2_grp=WF_phases_b_map(iv2,ipmqbz,isq,2)
             !
             if (i_v1_grp/=i_v2_grp) cycle
             !
             do ic1 = BS_bands(1),BS_bands(2)
               !
               i_c1_bnd=WF_phases_b_map(ic1,ipbz,isq,1)
               i_c1_grp=WF_phases_b_map(ic1,ipbz,isq,2)
               !
               neh1=BSS_eh_table_m1(ikbz,iv1,ic1,i_sp_pol,1)
               !
               if (neh1==0) cycle
               !
               do ic2 = BS_bands(1),BS_bands(2) 
                 !
                 i_c2_bnd=WF_phases_b_map(ic2,ipbz,isq,1)
                 i_c2_grp=WF_phases_b_map(ic2,ipbz,isq,2)
                 !
                 if (i_c1_grp/=i_c2_grp) cycle
                 !
                 neh2=BSS_eh_table_m1(ipbz,iv2,ic2,i_sp_pol,1)
                 !
                 if (neh2==0) cycle
                 !
                 phase_vv=conjg(WF_phases_vv(i_v1_bnd,i_v2_bnd,i_v1_grp,isq,1,1))
                 !phase_vv=      WF_phases_vv(i_v2_bnd,i_v1_bnd,i_v1_grp,isq,1,1)
                 !
                 phase_cc=      WF_phases_cc(i_c1_bnd,i_c2_bnd,i_c1_grp,isq,1,1)
                 !
                 EXC_WFs(neh1,iexc)=EXC_WFs(neh1,iexc)+BS_mat(neh2,il)*phase_vv*phase_cc
                 !
               enddo
             enddo
             !
           enddo
         enddo
         !
       enddo
     enddo
   enddo
   !
   !write(*,*) il,"Al= ",BS_mat(1:10,il)
   !write(*,*) iexc,"Arot= ",EXC_WFs(1:10,iexc)
   !
 enddo
 !
 !! take care of time reversal if needed
 !!
 !if(isymm>nsym/(1+i_time_rev)) call dev_conjg(WF_p)
 !
end subroutine EXC_WF_symm_qpoint
