!
! License-Identifier: GPL
!
! Copyright (C) 2015 The Yambo Team
!
! Authors (see AUTHORS file for details): AM
!
subroutine TDDFT_ALDA_G_space(E,k)
 !
 use pars,          ONLY:DP,SP
 use X_m,           ONLY:l_X_magn,l_X_abs
 use TDDFT,         ONLY:FXC_n_g_corr,F_xc_gspace,l_Fxc_in_Xo,l_Fxc_from_Vxc
 use QP_m,          ONLY:QP_ng_Vxc
 use com,           ONLY:msg
 use FFT_m,         ONLY:fft_size,fftw_plan,fft_dim,fft_g_table
 use electrons,     ONLY:levels,n_spin,n_spinor
 use R_lattice,     ONLY:bz_samp,G_m_G
 use wave_func,     ONLY:WF
 use xc_functionals,ONLY:V_xc,F_xc,F_xc_mat,magn,XC_potential_driver
 use global_XC,     ONLY:WF_xc_functional,WF_kind
 use interfaces,    ONLY:WF_load,WF_free,eval_G_minus_G
 !
#include<memory.h>
 !
 type(levels) ::E
 type(bz_samp)::k
 ! 
 ! Work Space
 !
 integer                  :: i1,i2,nG
 complex(DP), allocatable :: F_xc_DP(:)
 !
 integer, external        :: G_index_energy_factor
 !
 if (allocated(F_xc_gspace).and..not.l_Fxc_in_Xo) return
 if (allocated(F_xc_mat)   .and.     l_Fxc_in_Xo) return
 !
 ! Fxc RL size setup looking at the max index of G-G'
 !
 if (.not.l_Fxc_in_Xo) then
   if (FXC_n_g_corr/= eval_G_minus_G(FXC_n_g_corr,0) ) then
     FXC_n_g_corr   = eval_G_minus_G(FXC_n_g_corr,0)
     call msg('s','[TDDFT] ALDA Fxc size forced by symmetry to be ',FXC_n_g_corr)
   endif
   !
   nG=G_index_energy_factor(FXC_n_g_corr,4._SP)
 else
   !
   nG=G_index_energy_factor(FXC_n_g_corr,1._SP)
   !
 endif
 !
 ! WF distributed & load
 !=======================
 !
 call WF_load(WF,nG,1,(/1,E%nbm/),(/1,k%nibz/),title='-Rho')
 !
 ! ALDA functional
 !
 YAMBO_ALLOC(F_xc,(fft_size,n_spin,n_spin))
 if (l_X_magn.or.n_spinor==2.or.l_Fxc_from_Vxc) then
   YAMBO_ALLOC(V_xc,(fft_size,n_spin))
   V_xc=0._SP
   QP_ng_Vxc=nG
 endif
 if (n_spin>1) then
   YAMBO_ALLOC(magn,(fft_size,3))
 endif
 !
 call XC_potential_driver(E,k,WF_kind,WF_xc_functional,2)
 !
 if (l_Fxc_in_Xo) then
   !
   YAMBO_ALLOC(F_xc_mat,(fft_size,n_spin,n_spin,n_spin,n_spin))
   call Build_F_xc_mat(V_xc,F_xc,F_xc_mat)
   goto 100
   !
 endif
 !
 YAMBO_ALLOC(F_xc_DP,(fft_size))
 if(l_X_abs) then
   do i1=1,fft_size
     F_xc_DP(i1)=sum(cmplx(F_xc(i1,:,:),kind=DP))
   enddo
 else if(l_X_magn) then
   ! DS warning. This is exact for the collinear case only.
   ! The non collinear case implementation should be extended, while here the 
   ! mixing between the +-, --, -+ and ++ components is neglected in the G-space implementation
   do i1=1,fft_size
     F_xc_DP(i1)=    cmplx(F_xc(i1,1,2),kind=DP)
   enddo
 endif
 !
 ! need to reset the plan to do a BACKWARD FFT (R->G)
 !
 fftw_plan=0
#if defined _FFTW
 call fft_3d(F_xc_DP,fft_dim,-1,fftw_plan)
#else
 call fft_3d(F_xc_DP,fft_dim,-1)
#endif
 !
 ! F_xc in gspace
 !
 YAMBO_ALLOC(F_xc_gspace,(FXC_n_g_corr,FXC_n_g_corr,1))
 do i1=1,FXC_n_g_corr
   do i2=1,FXC_n_g_corr
     F_xc_gspace(i1,i2,1)=cmplx(F_xc_DP( fft_g_table( G_m_G(i1,i2),1) ),kind=SP)/real(fft_size,SP)
   enddo
 enddo
 !
 ! Clean
 !
100 continue
 !
 YAMBO_FREE(V_xc)
 YAMBO_FREE(F_xc)
 YAMBO_FREE(F_xc_DP)
 YAMBO_FREE(V_xc)
 YAMBO_FREE(magn)
 call WF_free(WF)
 !
end subroutine
