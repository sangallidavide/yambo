!
! License-Identifier: GPL
!
! Copyright (C) 2007 The Yambo Team
!
! Authors (see AUTHORS file for details): MG
!
subroutine V_Hartree(rho,vhr)
 !
 !  Hartree potential (vhr)
 !  sum_G 4pi*rho(G)/|G|**2.  with G!=0, for G=0 vhr(0)=0 
 !
 ! V_H(r) is define as the "physical hartree potential".
 ! The density in yambo rho_y is instead the "computationally more convenient" (See comment in src/xc_functionals/el_density_and_current.F)
 ! This means, V_H(r) should be constructed from rho_y(r)*fft_size/V
 ! This explains the 1/V factor (DL_vol) in V_Hartree and the missing 1/fft_size after the Fourier transform compared to other subroutines
 !
 use pars,          ONLY:pi,SP,DP
 use electrons,     ONLY:n_spin
 use R_lattice,     ONLY:bare_qpg
 use QP_m,          ONLY:QP_ng_SH  
 use D_lattice,     ONLY:DL_vol
 use fft_m,         ONLY:fft_size,fft_dim,fftw_plan,fft_g_table
 use timing_m,      ONLY:timing
 ! 
 implicit none
 !
 real(SP), intent(in)  :: rho(fft_size)
 real(SP), intent(out) :: vhr(fft_size,n_spin)
 !
 ! Work Space
 !
 complex(DP) :: rhog(fft_size),vhtmp(fft_size),vhg(QP_ng_SH)
 integer     :: ig
 !
 call timing('V_Hartree',OPR='start')
 !
!$OMP WORKSHARE  
 vhg   = (0._DP,0._DP)
 vhtmp = (0._DP,0._DP)
 rhog  = rho
!$OMP END WORKSHARE  
 !
 ! FFT call: rho in reciprocal space
 !
#if defined _FFTW
 call dfftw_destroy_plan(fftw_plan)
 fftw_plan = 0
 call fft_3d(rhog,fft_dim,-1,fftw_plan)
#else
 call fft_3d(rhog,fft_dim,-1)
#endif
 !
 !$omp parallel do default(shared), private(ig)
 do ig=2,QP_ng_SH
   vhg(ig)=((4._DP*pi/DL_vol)*(rhog(fft_g_table(ig,1))))/bare_qpg(1,ig)**2
 enddo
 !$omp end parallel do
 !
 !$omp parallel do default(shared), private(ig)
 do ig=1,QP_ng_SH
   vhtmp(fft_g_table(ig,1))=vhg(ig)
 enddo
 !$omp end parallel do
 !
 ! FFT call: vhr in direct space
 !
#if defined _FFTW
 call dfftw_destroy_plan(fftw_plan)
 fftw_plan = 0
 call fft_3d(vhtmp,fft_dim,1,fftw_plan)
#else
 call fft_3d(vhtmp,fft_dim,1)
#endif
 !
!$OMP WORKSHARE  
 vhr(:,1) = real(vhtmp,SP)
 !$OMP END WORKSHARE  
 if(n_spin==2) vhr(:,2)=vhr(:,1)
 !
 call timing('V_Hartree',OPR='stop')
 !
end subroutine V_hartree
