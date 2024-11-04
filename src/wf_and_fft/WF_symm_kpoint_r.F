!
! License-Identifier: GPL
!
! Copyright (C) 2023 The Yambo Team
!
! Authors (see AUTHORS file for details): DS AF
!
#include<dev_defs.h>                                                    
!
subroutine WF_symm_kpoint_r(b_to_load,ik,isymm,WF)
 !
 use pars,          ONLY:SP,cZERO
 use electrons,     ONLY:n_sp_pol,n_spinor
 use wave_func,     ONLY:WAVEs
 use FFT_m,         ONLY:DEV_VAR(fft_rot_r),fft_size
 use R_lattice,     ONLY:bz_samp
 use D_lattice,     ONLY:idt_index,DEV_VAR(spin_sop),i_time_rev,nsym
 use devxlib,       ONLY:devxlib_conjg_d,devxlib_memcpy_d2d
 !
 implicit none
 !
 integer,       intent(in)           :: ik,isymm,b_to_load(2)
 type (WAVEs), target, intent(inout) :: WF
 !
 complex(SP)              DEV_ATTR :: WF_k_tmp(fft_size,n_spinor,b_to_load(1):b_to_load(2),n_sp_pol)
 complex(SP), pointer     DEV_ATTR :: WF_p(:,:,:)
 integer :: i_sp_pol,ib,i_wf
 integer :: ifft,i_spinor
 !
 if(isymm==idt_index) return
 !
 WF_p => DEV_VAR(WF%r)
 !
#if defined _GPU
 call devxlib_memset_d(WF_k_tmp,cZERO)
#else
 WF_k_tmp=cZERO
#endif
 !
 do i_sp_pol=1,n_sp_pol
   do ib=b_to_load(1),b_to_load(2)
     !
     i_wf=WF%index(ib,ik,i_sp_pol)
     !
     if (n_spinor==1) then
       !
#ifdef _GPU
       !$cuf kernel do(1) <<<*,*>>>
#else
       !$omp parallel do default(shared), private(ifft)
#endif
       do ifft = 1, fft_size
         WF_k_tmp(ifft,1,ib,i_sp_pol)=WF_p(DEV_VAR(fft_rot_r)(ifft,isymm),1,i_wf)
       enddo
       !
     else
       !
#ifdef _GPU
       !$cuf kernel do(1) <<<*,*>>>
#else
       !$omp parallel do default(shared), private(ifft)
#endif
       do ifft = 1, fft_size
         WF_k_tmp(ifft,1,ib,1)=DEV_VAR(spin_sop)(1,1,isymm)*WF_p(DEV_VAR(fft_rot_r)(ifft,isymm),1,i_wf)+ &
                 DEV_VAR(spin_sop)(1,2,isymm)*WF_p(DEV_VAR(fft_rot_r)(ifft,isymm),2,i_wf)
         WF_k_tmp(ifft,2,ib,1)=DEV_VAR(spin_sop)(2,1,isymm)*WF_p(DEV_VAR(fft_rot_r)(ifft,isymm),1,i_wf)+ &
                 DEV_VAR(spin_sop)(2,2,isymm)*WF_p(DEV_VAR(fft_rot_r)(ifft,isymm),2,i_wf)
       enddo
       !
     endif
     !
     call devxlib_memcpy_d2d(WF_p(:,:,i_wf),WF_k_tmp(:,:,ib,i_sp_pol))
     !
   enddo
 enddo
 !
 ! take care of time reversal if needed
 !
#if defined _GPU
 if(isymm>nsym/(1+i_time_rev)) call devxlib_conjg_d(WF_p)
#else
 if(isymm>nsym/(1+i_time_rev)) WF_p=conjg(WF_p)
#endif
 !
end subroutine WF_symm_kpoint_r
