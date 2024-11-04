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
 use devxlib,       ONLY:devxlib_conjg,devxlib_memcpy_d2d,devxlib_memset_d
 !
 implicit none
 !
 integer,       intent(in)  :: ik,isymm,b_to_load(2)
 type (WAVEs), target, intent(inout)  :: WF
 !
 complex(SP)              DEV_ATTR :: WF_k_tmp(wf_ng_1st_BZ,n_spinor,b_to_load(2)-b_to_load(1)+1,n_sp_pol)
 complex(SP), pointer     DEV_ATTR :: WF_p(:,:,:)
 integer :: i_sp_pol,ibl,ib,i_wf
 integer :: i_g,i_spinor
 !
 WF_p => DEV_VAR(WF%c)
 !
#if defined _GPU_LOC
 call devxlib_memset_d(WF_k_out,cZERO)
#else
 WF_k_tmp=cZERO
#endif
 !
 if(isymm==idt_index) return
   !
 do i_sp_pol=1,n_sp_pol
   do ib=b_to_load(1),b_to_load(2)
     !
     i_wf=WF%index(ib,ik,i_sp_pol)
     ibl=ib-b_to_load(1)+1
     !
     if (n_spinor==1) then
       !
       !DEV_ACC_DEBUG data present(WF_k_tmp,g_rot,WF_p)
       !DEV_ACC parallel loop collapse(1)
       !DEV_CUF kernel do(1) <<<*,*>>>
       !DEV_OMPGPU target map(present,alloc:WF_k_tmp,g_rot,WF_p)
       !DEV_OMPGPU teams loop collapse(1)
       do i_g=1,wf_ng
         WF_k_tmp(DEV_VAR(g_rot)(i_g,isymm),1,ibl,i_sp_pol)=WF_p(i_g,1,i_wf)
       enddo
       !
     else
       !
       !DEV_ACC_DEBUG data present(WF_k_out,g_rot,WF_p,spin_sop)
       !DEV_ACC parallel loop
       !DEV_CUF kernel do(1) <<<*,*>>>
       !DEV_OMPGPU target map(present,alloc:WF_k_out,g_rot,WF_p,spin_sop)
       !DEV_OMPGPU teams loop
       do i_g=1,wf_ng
         WF_k_tmp(DEV_VAR(g_rot)(i_g,isymm),1,ibl,1)=DEV_VAR(spin_sop)(1,1,isymm)*WF_p(i_g,1,i_wf)+ &
&                                                       DEV_VAR(spin_sop)(1,2,isymm)*WF_p(i_g,2,i_wf)
         WF_k_tmp(DEV_VAR(g_rot)(i_g,isymm),2,ibl,1)=DEV_VAR(spin_sop)(2,1,isymm)*WF_p(i_g,1,i_wf)+ &
&                                                       DEV_VAR(spin_sop)(2,2,isymm)*WF_p(i_g,2,i_wf)
       enddo
     endif
     !
#if defined _GPU
     call devxlib_memcpy_d2d(WF_p(:,:,iwf),WF_k_tmp(:,:,ib,i_sp_pol))
#else
     WF_p(:,:,iwf)=WF_k_tmp(:,:,ib,i_sp_pol)
#endif
     !
   enddo
 enddo
 !
 ! take care of time reversal if needed
 !
#ifdef _GPU_LOC
 if(is>nsym/(1+i_time_rev)) call devxlib_conjg(WF_p)
#else
 if(is>nsym/(1+i_time_rev)) WF_p=conjg(WF_p)
#endif
 !
end subroutine WF_symm_kpoint_g
