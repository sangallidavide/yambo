!
! License-Identifier: GPL
!
! Copyright (C) 2013 The Yambo Team
!
! Authors (see AUTHORS file for details): DS AF
!
subroutine DEV_SUB_ALT(WF_apply_symm)(isc,WF_symm)
 !
 ! INCLUDED in: WF_apply_symm_incl.F 
 !
 use pars,          ONLY:SP
 use wave_func,     ONLY:WF,rho_map_size,rho_map,rho_map_inv
 use FFT_m,         ONLY:DEV_VAR(fft_rot_r),fft_size
 use stderr,        ONLY:intc
 use parallel_m,    ONLY:myid
 use electrons,     ONLY:n_spinor
 use D_lattice,     ONLY:nsym,DEV_VAR(spin_sop),i_time_rev,idt_index
 use cuda_m,        ONLY:have_cuda
 use deviceXlib_m
 !
#include<yambo_wfs.h>
 !
 implicit none
 !
 integer,     intent(in)  :: isc(4)
 complex(SP), intent(out) DEV_ATTR :: WF_symm(rho_map_size,n_spinor)
 !
 WF_RSPACE(SP), pointer     DEV_ATTR :: WF_p(:,:)
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
 WF_p => DEV_VAR(WF%r)(:,:,i_wf)
 isymm = isc(3)
 !
 if(isc(3)==idt_index) then
   !
   call dev_memcpy(WF_symm,cmplx(WF_p,kind=SP))
   return
   !
 endif
 !
#ifdef _CUDA_LOC
 if (allocated(rho_map)) call error(" rho map procedure not implemented with symmetries and cuda")
#endif
 !
 if(n_spinor==1) then
   !
   if (allocated(rho_map)) then
     !$omp parallel do default(shared), private(ifft,ifftp)
     do ifft = 1, rho_map_size
       ifftp=rho_map_inv(DEV_VAR(fft_rot_r)(rho_map(ifft),isymm))
       WF_symm(ifft,1)=cmplx(WF_p(ifftp,1),kind=SP)
     enddo
   else 
#ifdef _CUDA_LOC   
     !$cuf kernel do(1) <<<*,*>>>
#else
     !$omp parallel do default(shared), private(ifft,ifftp)
#endif
     do ifft = 1, fft_size
       ifftp=DEV_VAR(fft_rot_r)(ifft,isymm)
       WF_symm(ifft,1)=cmplx(WF_p(ifftp,1),kind=SP)
     enddo
   endif
   !
 endif
 !
 if (n_spinor==2) then
   !
   if (allocated(rho_map)) then
     !$omp parallel do default(shared), private(ifft,ifftp)
     do ifft = 1, rho_map_size
       ifftp=rho_map_inv(DEV_VAR(fft_rot_r)(rho_map(ifft),isymm))
       WF_symm(ifft,1)=DEV_VAR(spin_sop)(1,1,isymm)*cmplx(WF_p(ifftp,1),kind=SP)+ &
                       DEV_VAR(spin_sop)(1,2,isymm)*cmplx(WF_p(ifftp,2),kind=SP)
       WF_symm(ifft,2)=DEV_VAR(spin_sop)(2,1,isymm)*cmplx(WF_p(ifftp,1),kind=SP)+ &
                       DEV_VAR(spin_sop)(2,2,isymm)*cmplx(WF_p(ifftp,2),kind=SP)
     enddo
   else
#ifdef _CUDA_LOC
     !$cuf kernel do(1) <<<*,*>>>
#else
     !$omp parallel do default(shared), private(ifft,ifftp)
#endif
     do ifft = 1, fft_size
       ifftp=DEV_VAR(fft_rot_r)(ifft,isymm)
       WF_symm(ifft,1)=DEV_VAR(spin_sop)(1,1,isymm)*cmplx(WF_p(ifftp,1),kind=SP)+ &
                       DEV_VAR(spin_sop)(1,2,isymm)*cmplx(WF_p(ifftp,2),kind=SP)
       WF_symm(ifft,2)=DEV_VAR(spin_sop)(2,1,isymm)*cmplx(WF_p(ifftp,1),kind=SP)+ &
                       DEV_VAR(spin_sop)(2,2,isymm)*cmplx(WF_p(ifftp,2),kind=SP)
     enddo
   endif
   !
 endif
 !
 ! take care of time reversal if needed
 !
 if(isc(3)>nsym/(i_time_rev+1)) call dev_conjg(WF_symm)
 !
end subroutine DEV_SUB_ALT(WF_apply_symm)
