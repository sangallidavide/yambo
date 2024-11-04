!
! License-Identifier: GPL
!
! Copyright (C) 2015 The Yambo Team
!
! Authors (see AUTHORS file for details): MG CA
!
#include<dev_defs.h>
!
subroutine WF_shift_kpoint_gpu(b_to_shift,nb_to_shift,ikbz,i_sp_pol,wf_shift,Xk,WF_k_out)
 !
 use pars,          ONLY: SP,cZERO
 use electrons,     ONLY: n_spinor,n_sp_pol
 use R_lattice,     ONLY: k_map,bz_samp
 use R_lattice,     ONLY: DEV_VAR(G_m_G)
 use wave_func,     ONLY: wf_ng_1st_BZ,wf_ng_overlaps,WF
 use devxlib,       ONLY: devxlib_memcpy_d2d, devxlib_memset_d
 use gpu_m,         ONLY: have_gpu
 !
#include<y_memory.h>
 !
 integer,       intent(in) :: wf_shift(3),ikbz,i_sp_pol,b_to_shift(2),nb_to_shift
 type(bz_samp), intent(in) :: Xk
 complex(SP), intent(out) DEV_ATTR :: WF_k_out(wf_ng_overlaps,n_spinor,nb_to_shift)
 !
 ! Work space
 !
 complex(SP), pointer     DEV_ATTR :: WF_k_in(:,:,:)
 !
 integer :: id,ik,is
 complex(SP), allocatable DEV_ATTR :: WF_tmp(:,:,:)
 integer :: ig,igp,i_b1,i_b2,i_bl,i_b_ref,i_wf
 integer :: g0_idx(3,2),g0_idx_val
 integer :: i_spinor,i_b
 !
 WF_k_in => DEV_VAR(WF%c)
 !
 ik = Xk%sstar(ikbz,1)
 is = Xk%sstar(ikbz,2)
 !
 i_b1= 1
 i_b2= b_to_shift(2)-b_to_shift(1)+1
 !
 if (i_b2 /= nb_to_shift) call error("[WF_shift_kpoint] Unexpected error")
 !
 YAMBO_ALLOC_GPU(WF_tmp,(wf_ng_overlaps,n_spinor,nb_to_shift))
 call devxlib_memset_d(WF_tmp,cZERO)
 !
 do i_b=b_to_shift(1),b_to_shift(2)
   i_wf=WF%index(i_b,ik,i_sp_pol)
   i_bl=i_b-b_to_shift(1)+1
   ! dev2dev
   call devxlib_memcpy_d2d(WF_tmp(:,:,i_b),WF_k_in(:,:,i_wf), range1=(/1,wf_ng_1st_BZ/))
 enddo
 !

 !
 ! main loop
 !
 do id=1,3
   !
   WF_k_out=cZERO
   !
   if(WF_shift(id)/=0) then
     !
     g0_idx_val=k_map%g0_idx(id,WF_shift(id))
     !
     !DEV_ACC_DEBUG data present(WF_tmp,G_m_G,WF_k_out)
     !DEV_ACC parallel loop collapse(3)
     !DEV_CUF kernel do(3) <<<*,*>>>
     !DEV_OMPGPU target map(present,alloc:WF_tmp,G_m_G,WF_k_out)
     !DEV_OMPGPU teams loop collapse(3)
     !DEV_OMP parallel do default(shared), private(i_sp_pol,i_b,i_spinor,ig,igp),collapse(3)
     !
     do i_b=i_b1,i_b2
     do i_spinor=1,n_spinor
     do ig=1,wf_ng_1st_BZ
       !
       if(WF_tmp(ig,1,i_b1)==cZERO) cycle
       !
       igp=DEV_VAR(G_m_G)(ig,g0_idx_val)
       WF_k_out(igp,i_spinor,i_b)=WF_tmp(ig,i_spinor,i_b)
       !
     enddo
     enddo
     enddo
     !
     !DEV_OMP end parallel do
     !DEV_OMPGPU end target
     !DEV_ACC_DEBUG end data
     !
   else
     ! dev2dev
     call devxlib_memcpy_d2d(WF_k_out,WF_tmp)
   endif
   ! dev2dev
   call devxlib_memset_d(WF_tmp,cZERO)
   call devxlib_memcpy_d2d(WF_tmp,WF_k_out)
   !
 enddo
 !
 YAMBO_FREE_GPU(WF_tmp)
 YAMBO_FREE(WF_tmp)
 !
end subroutine WF_shift_kpoint_gpu

