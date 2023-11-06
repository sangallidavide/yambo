!
! License-Identifier: GPL
!
! Copyright (C) 2015 The Yambo Team
!
! Authors (see AUTHORS file for details): MG CA DS AF
!
#include<dev_defs.h>                                                    
!
subroutine WF_symm_kpoint_g(b_to_load,ik,isymm,WF)
 !
 use pars,          ONLY:SP,cZERO
 use electrons,     ONLY:n_sp_pol,n_spinor
 use wave_func,     ONLY:WAVEs
 use wave_func,     ONLY:wf_ng,wf_ng_1st_BZ
 use R_lattice,     ONLY:bz_samp,DEV_VAR(g_rot)
 use D_lattice,     ONLY:idt_index,DEV_VAR(spin_sop),i_time_rev,nsym
 use deviceXlib_m,  ONLY:dev_conjg,dev_memcpy
 !
 implicit none
 !
 integer,       intent(in)  :: ik,isymm,b_to_load(2)
 type (WAVEs), target, intent(inout)  :: WF
 !
 complex(SP)              DEV_ATTR :: WF_k_tmp(wf_ng_1st_BZ,n_spinor,b_to_load(1):b_to_load(2),n_sp_pol)
 complex(SP), pointer     DEV_ATTR :: WF_p(:,:,:)
 integer :: i_sp_pol,ib,i_wf
 integer :: i_g,i_spinor
 !
 WF_p => DEV_VAR(WF%c)
 !
 WF_k_tmp=cZERO
 !
 if(isymm==idt_index) return
   !
 do i_sp_pol=1,n_sp_pol
   do ib=b_to_load(1),b_to_load(2)
     !
     i_wf=WF%index(ib,ik,i_sp_pol)
     !
     if (n_spinor==1) then
       !
#ifdef _CUDA
       !$cuf kernel do(1) <<<*,*>>>
#endif
       do i_g=1,wf_ng
         WF_k_tmp(DEV_VAR(g_rot)(i_g,isymm),1,ib,i_sp_pol)=WF_p(i_g,1,i_wf)
       enddo
       !
     else
       !
#ifdef _CUDA
       !$cuf kernel do(1) <<<*,*>>>
#endif
       do i_g=1,wf_ng
         WF_k_tmp(DEV_VAR(g_rot)(i_g,isymm),1,ib,1)=DEV_VAR(spin_sop)(1,1,isymm)*WF_p(i_g,1,i_wf)+ &
&                                                       DEV_VAR(spin_sop)(1,2,isymm)*WF_p(i_g,2,i_wf)
         WF_k_tmp(DEV_VAR(g_rot)(i_g,isymm),2,ib,1)=DEV_VAR(spin_sop)(2,1,isymm)*WF_p(i_g,1,i_wf)+ &
&                                                       DEV_VAR(spin_sop)(2,2,isymm)*WF_p(i_g,2,i_wf)
       enddo
     endif
     !
     call dev_memcpy(WF_p(:,:,i_wf),WF_k_tmp(:,:,ib,i_sp_pol))
     !
   enddo
 enddo
 !
 ! take care of time reversal if needed
 !
 if(isymm>nsym/(1+i_time_rev)) call dev_conjg(WF_p)
 !
end subroutine WF_symm_kpoint_g
