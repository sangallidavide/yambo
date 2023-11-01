!
! License-Identifier: GPL
!
! Copyright (C) 2013 The Yambo Team
!
! Authors (see AUTHORS file for details): DS AF
!
subroutine DEV_SUB_ALT(WF_apply_symm)(isc,WF,WF_symm)
 !
 ! INCLUDED in: WF_apply_symm_incl.F 
 !
 use pars,          ONLY:SP
 use wave_func,     ONLY:WAVEs,rho_map_size,rho_map,WF_shifted_grids
 use FFT_m,         ONLY:DEV_VAR(fft_rot_r),fft_size
 use stderr,        ONLY:intc
 use parallel_m,    ONLY:myid
 use electrons,     ONLY:n_spinor
 use D_lattice,     ONLY:nsym,DEV_VAR(spin_sop),i_time_rev,idt_index
 use deviceXlib_m
 !
#include<yambo_wfs.h>
 !
 implicit none
 !
 integer,     intent(in)  :: isc(4)
 type (WAVEs), target, intent(in)  :: WF
 complex(SP), intent(out) DEV_ATTR :: WF_symm(rho_map_size,n_spinor)
 !
 WF_RSPACE(SP), pointer     DEV_ATTR :: WF_p(:,:,:)
 integer  :: i_wf
 integer  :: ifft,ifftp,isymm
 !
 i_wf=WF%index(isc(1),isc(2),isc(4))
 !
#if defined _MPI
 if (i_wf==0) call error('EMPTY WF index @ b'//trim(intc(isc(1)))//' k'//&
&                         trim(intc(isc(2)))//' spin'//trim(intc(isc(4)))//&
&                         ' ID'//trim(intc(myid)))
#endif
 !
 WF_p => DEV_VAR(WF%r)
 isymm = isc(3)
 !
 if(isymm==idt_index .or. WF_shifted_grids) then
   !
   call dev_memcpy(WF_symm,cmplx(DEV_VAR(WF%r)(:,:,i_wf),kind=SP))
   return
   !
 endif
 !
 if (allocated(rho_map)) call error(" rho map procedure not implemented with symmetries")
 !
 if(n_spinor==1) then
   !
#ifdef _CUDA_LOC   
    !$cuf kernel do(1) <<<*,*>>>
#else
    !$omp parallel do default(shared), private(ifft,ifftp)
#endif
   do ifft = 1, fft_size
     ifftp=DEV_VAR(fft_rot_r)(ifft,isymm)
     WF_symm(ifft,1)=cmplx(WF_p(ifftp,1,i_wf),kind=SP)
   enddo
   !
 endif
 !
 if (n_spinor==2) then
   !
#ifdef _CUDA_LOC
   !$cuf kernel do(1) <<<*,*>>>
#else
   !$omp parallel do default(shared), private(ifft,ifftp)
#endif
   do ifft = 1, fft_size
     ifftp=DEV_VAR(fft_rot_r)(ifft,isymm)
     WF_symm(ifft,1)=DEV_VAR(spin_sop)(1,1,isymm)*WF_p(ifftp,1,i_wf)+ &
                     DEV_VAR(spin_sop)(1,2,isymm)*WF_p(ifftp,2,i_wf)
     WF_symm(ifft,2)=DEV_VAR(spin_sop)(2,1,isymm)*WF_p(ifftp,1,i_wf)+ &
                     DEV_VAR(spin_sop)(2,2,isymm)*WF_p(ifftp,2,i_wf)
   enddo
   !
 endif
 !
 ! take care of time reversal if needed
 !
 if(isc(3)>nsym/(i_time_rev+1)) call dev_conjg(WF_symm)
 !
end subroutine DEV_SUB_ALT(WF_apply_symm)
