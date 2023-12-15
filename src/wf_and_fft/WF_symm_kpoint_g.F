!
! License-Identifier: GPL
!
! Copyright (C) 2015 The Yambo Team
!
! Authors (see AUTHORS file for details): MG CA DS AF
!
#include<dev_defs.h>                                                    
!
subroutine WF_symm_kpoint_g(b_to_load,ik,i_sp_pol,isymm,ig0,WF)
 !
 use pars,          ONLY:SP,cZERO
 use electrons,     ONLY:n_spinor
 use wave_func,     ONLY:WAVEs,wf_nc_k,wf_igk
 use wave_func,     ONLY:wf_ng,wf_ng_1st_BZ
 use R_lattice,     ONLY:DEV_VAR(g_rot),DEV_VAR(G_m_G)
 use deviceXlib_m,  ONLY:dev_conjg,dev_memcpy
 use D_lattice,     ONLY:idt_index,DEV_VAR(spin_sop),i_time_rev,nsym
 !
 implicit none
 !
 integer,       intent(in) :: b_to_load(2)
 integer,       intent(in) :: ik,i_sp_pol,isymm,ig0
 type (WAVEs), target, intent(inout)  :: WF
 !
 ! Work space
 !
 complex(SP)              DEV_ATTR :: WF_k_tmp(wf_ng_1st_BZ,n_spinor)
 complex(SP), pointer     DEV_ATTR :: WF_p(:,:,:)
 !
 integer :: ib,i_wf
 integer :: ic,ig,igp,i_spinor
 !
 if(isymm==idt_index) return
 !
 WF_p => DEV_VAR(WF%c)
 !
 do ib=b_to_load(1),b_to_load(2)
   !
   WF_k_tmp=cZERO
   !
   i_wf=WF%index(ib,ik,i_sp_pol)
   !
   if (n_spinor==1) then
     !
#ifdef _CUDA
     !$cuf kernel do(1) <<<*,*>>>
#endif
     do ic=1,wf_nc_k(ik)
       ig=wf_igk(ic,ik)
       igp=DEV_VAR(G_m_G)(DEV_VAR(g_rot)(ig,isymm),ig0)
       !if ( ib==b_to_load(1) ) write(*,*) ic,ig,igp
       WF_k_tmp(igp,1)=WF_p(ig,1,i_wf)
     enddo
     !
   else
     !
#ifdef _CUDA
     !$cuf kernel do(1) <<<*,*>>>
#endif
     do ic=1,wf_nc_k(ik)
       ig=wf_igk(ic,ik)
       igp=DEV_VAR(G_m_G)(DEV_VAR(g_rot)(ig,isymm),ig0)
       WF_k_tmp(igp,1)=DEV_VAR(spin_sop)(1,1,isymm)*WF_p(ig,1,i_wf)+ &
&                         DEV_VAR(spin_sop)(1,2,isymm)*WF_p(ig,2,i_wf)
       WF_k_tmp(igp,2)=DEV_VAR(spin_sop)(2,1,isymm)*WF_p(ig,1,i_wf)+ &
&                         DEV_VAR(spin_sop)(2,2,isymm)*WF_p(ig,2,i_wf)
     enddo
   endif
   !
   call dev_memcpy(WF_p(:,:,i_wf),WF_k_tmp)
   !
 enddo
 !
 ! take care of time reversal if needed
 !
 if(isymm>nsym/(1+i_time_rev)) call dev_conjg(WF_p)
 !
end subroutine WF_symm_kpoint_g
