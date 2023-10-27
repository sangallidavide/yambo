!
! License-Identifier: GPL
!
! Copyright (C) 2015 The Yambo Team
!
! Authors (see AUTHORS file for details): DS
!
subroutine TDDFT_ALDA_eh_space_G_collisions_L(iq,Xk,i_T_grp,NG,l_bs_exch_wf_in_loop,tddft_wf)
 !
 ! Calculates the F_xc scattering 
 !
 !  (ic(1),ik(1),is(1)) --<--:...:--<-- (ic(2),ik(2),is(2))
 !                       :sqrt(Fxc):
 !
 use pars,           ONLY:SP,DP,pi,cZERO
 use FFT_m,          ONLY:fft_size
 use wave_func,      ONLY:rho_map_size,rho_map,WF
 use interfaces,     ONLY:WF_load,WF_free
 use R_lattice,      ONLY:qindx_X,bz_samp
 use BS,             ONLY:l_BS_magnons,TDDFT_oscillators_alloc,BS_T_grp,BS_bands
 use electrons,      ONLY:n_spinor,n_sp_pol,n_spin
 use timing_m,       ONLY:timing
 use openmp,         ONLY:OPENMP_update,master_thread
 !use interfaces,     ONLY:WF_apply_symm
 use TDDFT,          ONLY:tddft_wf_t
 !
 implicit none
 !
 type(bz_samp),intent(in) :: Xk
 integer,      intent(in) :: iq,i_T_grp,NG(2)
 logical,      intent(in) :: l_bs_exch_wf_in_loop
 type(tddft_wf_t), target, intent(inout) :: tddft_wf 
 ! 
 ! Work Space
 !
 logical     :: l_use_rho_map,l_load_WFs
 integer     :: i_spinor,j_spinor,ip_spinor,jp_spinor,ifft,irho,irhotw,is_yambo,os_yambo,i_spinor_y,j_spinor_y
 integer     :: i_T_el,i_T_el_p,N_T_el_p,i_c,i_v,i_sp_c,i_sp_v,i_k_bz,i_k,i_s,i_g0,i_g1,i_g2,i_p_bz,i_g_p,i_p,&
&               i_sp,i_T_grp_p,is(4),os(4),NK(2)
 !
 complex(DP), pointer :: rhotwr_p(:)
 complex(SP), pointer :: WF_symm1_p(:,:)
 complex(SP), pointer :: WF_symm2_p(:,:)
 !
 call timing('T_space ALDA Osc.',OPR='start')
 !
 ! Allocate
 !
 call TDDFT_oscillators_alloc(i_T_grp)
 !
 WF_symm1_p => tddft_wf%WF_symm1
 WF_symm2_p => tddft_wf%WF_symm2
 rhotwr_p   => tddft_wf%rhotwr_DP
 !
 ! Use pointers both for CUDA and to avoid continuous allocation and de-allocation
 l_use_rho_map=allocated(rho_map)
 !
#if defined(__NOTNOW) && ! defined(_CUDA)
 !$omp parallel default(shared), private( K_EXCH_collision, &
 !$omp &        i_T_el,i_k_bz,i_k,i_s, i_p_bz,i_p,i_sp, i_v,i_c,i_sp_c,i_sp_v, &
 !$omp &        i_T_grp_last_with_identity_sym, i_T_el_last_with_identity_sym, &
 !$omp &        i_T_grp_p, N_T_el_p, i_T_el_p, i_g1, i_g2, Nk, is, os)
#endif
 !
 call OPENMP_update(master_thread)
 !
#if defined(__NOTNOW) && ! defined(_CUDA)
 !$omp do schedule(dynamic)
#endif
 NK=-1
 do i_T_el=1,BS_T_grp(i_T_grp)%size
   !
   i_k_bz=BS_T_grp(i_T_grp)%table(i_T_el,1)
   i_k   =Xk%sstar(i_k_bz,1)
   i_s   =Xk%sstar(i_k_bz,2)
   !
   i_p_bz=qindx_X(iq,i_k_bz,1)
   i_g_p =qindx_X(iq,i_k_bz,2)
   !
   i_p   =Xk%sstar(i_p_bz,1)
   i_sp  =Xk%sstar(i_p_bz,2)
   !
   i_v   =BS_T_grp(i_T_grp)%table(i_T_el,2)
   i_c   =BS_T_grp(i_T_grp)%table(i_T_el,3)
   i_sp_c=BS_T_grp(i_T_grp)%table(i_T_el,4)
   i_sp_v=BS_T_grp(i_T_grp)%table(i_T_el,5)
   !
   is=(/i_c,i_k,i_s,i_sp_c/)
   os=(/i_v,i_p,i_sp,i_sp_v/)
   !
#if defined(__NOTNOW) && ! defined(_CUDA)
   !$omp critical
#endif
   l_load_WFs= l_bs_exch_wf_in_loop .and. (NK(1)/=min(i_k,i_p).or.NK(2)/=max(i_k,i_p))
   if (l_load_WFs) then
     if (NK(2)/=-1) call WF_free(WF,keep_fft=.true.,keep_states_to_load=.true.)
     NK=(/min(i_k,i_p),max(i_k,i_p)/)
     call WF_load(WF,NG(1),NG(2),BS_bands,NK,k_extrema_only=.true.,quiet=.true.,&
     &            space='R',title="Kernel exch",keep_states_to_load=.true.)
   endif
#if defined(__NOTNOW) && ! defined(_CUDA)
   !$omp end critical
#endif
   !
   if ( (is(4)/=os(4)) .and. ((.not.l_BS_magnons).or. n_sp_pol==1) ) call error(" TDDFT_ALDA_R_space: wrong spin index")
   if ( (is(4)==os(4)) .and. ((     l_BS_magnons).and.n_sp_pol==2) ) call error(" TDDFT_ALDA_R_space: wrong spin index")
   !
   call WF_apply_symm_cpu(is,WF_symm1_p)
   call WF_apply_symm_cpu(os,WF_symm2_p)
   !
   if( n_spinor==1 ) then
     if(     l_use_rho_map) rhotwr_p(rho_map(:))=cmplx(WF_symm1_p(:,1)*conjg(WF_symm2_p(:,1)),kind=DP)
     if(.not.l_use_rho_map) rhotwr_p(        : )=cmplx(WF_symm1_p(:,1)*conjg(WF_symm2_p(:,1)),kind=DP)
     call perform_fft_3d(rhotwr_p,BS_T_grp(i_T_grp)%O_tddft(:,i_T_el,is(4),os(4),1))
   endif
   !
   if(n_spinor==2) then
     !
     do i_spinor=1,n_spinor
       do j_spinor=1,n_spinor
         !
         do irho=1,rho_map_size
           if(     l_use_rho_map) ifft=rho_map(irho)
           if(.not.l_use_rho_map) ifft=        irho
           rhotwr_p(ifft)=cmplx(WF_symm1_p(irho,i_spinor)*conjg(WF_symm2_p(irho,j_spinor)),kind=DP)
         enddo
         call perform_fft_3d(rhotwr_p,BS_T_grp(i_T_grp)%O_tddft(:,i_T_el,i_spinor,j_spinor,1))
         !
       enddo
     enddo
   endif
   !
 enddo
 !
 nullify(WF_symm1_p)
 nullify(WF_symm2_p)
 nullify(rhotwr_p)
 !
 call timing('T_space ALDA Osc.',OPR='stop')
 !
end subroutine TDDFT_ALDA_eh_space_G_collisions_L
!
!
subroutine TDDFT_ALDA_eh_space_G_collisions_R(iq,Xk,i_T_grp,NG,l_bs_exch_wf_in_loop,tddft_wf,mode)
 !
 ! Calculates the F_xc scattering 
 !
 !  (ic(1),ik(1),is(1)) --<--:...:--<-- (ic(2),ik(2),is(2))
 !                       :sqrt(Fxc):
 !
 use pars,           ONLY:SP,DP,pi,cZERO
 use FFT_m,          ONLY:fft_size
 use wave_func,      ONLY:rho_map_size,rho_map,WF
 use xc_functionals, ONLY:F_xc_mat
 use interfaces,     ONLY:WF_load,WF_free
 use R_lattice,      ONLY:qindx_X,bz_samp
 use BS,             ONLY:l_BS_magnons,TDDFT_oscillators_alloc,BS_T_grp,BS_bands
 use electrons,      ONLY:n_spinor,n_sp_pol,n_spin
 use timing_m,       ONLY:timing
 use openmp,         ONLY:OPENMP_update,master_thread
 !use interfaces,     ONLY:WF_apply_symm
 use TDDFT,          ONLY:tddft_wf_t
 !
 implicit none
 !
 type(bz_samp),intent(in) :: Xk
 integer,      intent(in) :: iq,i_T_grp,NG(2)
 logical,      intent(in) :: l_bs_exch_wf_in_loop
 character(3), intent(in) :: mode
 type(tddft_wf_t), target, intent(inout) :: tddft_wf 
 ! 
 ! Work Space
 !
 logical     :: l_use_rho_map,l_load_WFs
 integer     :: i_spinor,j_spinor,ip_spinor,jp_spinor,ifft,irho,irhotw,is_yambo,os_yambo,i_spinor_y,j_spinor_y,&
&               i_sp_loop,o_sp_loop
 integer     :: i_T_el,i_T_el_p,N_T_el_p,i_c,i_v,i_sp_c,i_sp_v,i_k_bz,i_k,i_s,i_g0,i_g1,i_g2,i_p_bz,i_g_p,i_p,&
&               i_sp,i_T_grp_p,is(4),os(4),NK(2)
 !
 complex(DP), pointer :: rhotwr_p(:)
 complex(SP), pointer :: WF_symm1_p(:,:)
 complex(SP), pointer :: WF_symm2_p(:,:)
 !
 call timing('T_space ALDA Osc.',OPR='start')
 !
 ! Allocate
 !
 if (.not.allocated( BS_T_grp(i_T_grp)%O_tddft )) call TDDFT_oscillators_alloc(i_T_grp)
 !
 WF_symm1_p => tddft_wf%WF_symm1
 WF_symm2_p => tddft_wf%WF_symm2
 rhotwr_p   => tddft_wf%rhotwr_DP
 !
 ! Use pointers both for CUDA and to avoid continuous allocation and de-allocation
 l_use_rho_map=allocated(rho_map)
 !
#if defined(__NOTNOW) && ! defined(_CUDA)
 !$omp parallel default(shared), private( K_EXCH_collision, &
 !$omp &        i_T_el,i_k_bz,i_k,i_s, i_p_bz,i_p,i_sp, i_v,i_c,i_sp_c,i_sp_v, &
 !$omp &        i_T_grp_last_with_identity_sym, i_T_el_last_with_identity_sym, &
 !$omp &        i_T_grp_p, N_T_el_p, i_T_el_p, i_g1, i_g2s, Nk, is, os)
#endif
 !
 call OPENMP_update(master_thread)
 !
#if defined(__NOTNOW) && ! defined(_CUDA)
 !$omp do schedule(dynamic)
#endif
 NK=-1
 do i_T_el=1,BS_T_grp(i_T_grp)%size
   !
   i_k_bz=BS_T_grp(i_T_grp)%table(i_T_el,1)
   i_k   =Xk%sstar(i_k_bz,1)
   i_s   =Xk%sstar(i_k_bz,2)
   !
   i_p_bz=qindx_X(iq,i_k_bz,1)
   i_g_p =qindx_X(iq,i_k_bz,2)
   !
   i_p   =Xk%sstar(i_p_bz,1)
   i_sp  =Xk%sstar(i_p_bz,2)
   !
   i_v   =BS_T_grp(i_T_grp)%table(i_T_el,2)
   i_c   =BS_T_grp(i_T_grp)%table(i_T_el,3)
   i_sp_c=BS_T_grp(i_T_grp)%table(i_T_el,4)
   i_sp_v=BS_T_grp(i_T_grp)%table(i_T_el,5)
   !
   is=(/i_c,i_k,i_s,i_sp_c/)
   os=(/i_v,i_p,i_sp,i_sp_v/)
   !
#if defined(__NOTNOW) && ! defined(_CUDA)
   !$omp critical
#endif
   l_load_WFs= l_bs_exch_wf_in_loop .and. (NK(1)/=min(i_k,i_p).or.NK(2)/=max(i_k,i_p))
   if (l_load_WFs) then
     if (NK(2)/=-1) call WF_free(WF,keep_fft=.true.,keep_states_to_load=.true.)
     NK=(/min(i_k,i_p),max(i_k,i_p)/)
     call WF_load(WF,NG(1),NG(2),BS_bands,NK,k_extrema_only=.true.,quiet=.true.,&
     &            space='R',title="Kernel exch",keep_states_to_load=.true.)
   endif
#if defined(__NOTNOW) && ! defined(_CUDA)
   !$omp end critical
#endif
   !
   if ((is(4)/=os(4)).and.((.not.l_BS_magnons).or. n_sp_pol==1)) call error(" TDDFT_ALDA_R_space: wrong spin index")
   if ((is(4)==os(4)).and.((     l_BS_magnons).and.n_sp_pol==2)) call error(" TDDFT_ALDA_R_space: wrong spin index")
   !
   if (mode=="RES") then
     call WF_apply_symm_cpu(is,WF_symm1_p) 
     call WF_apply_symm_cpu(os,WF_symm2_p)
   else if (mode=="CPL") then
     call WF_apply_symm_cpu(is,WF_symm2_p)
     call WF_apply_symm_cpu(os,WF_symm1_p)
   endif
   !
   if( n_spin==1 ) then
     !
     if (     l_use_rho_map) rhotwr_p(rho_map(:))=cmplx(F_xc_mat(rho_map(:),1,1,1,1)*&
     &                                             WF_symm1_p(:,1)*conjg(WF_symm2_p(:,1)),kind=DP)
     if (.not.l_use_rho_map) rhotwr_p(        : )=cmplx(F_xc_mat(        : ,1,1,1,1)*&
     &                                             WF_symm1_p(:,1)*conjg(WF_symm2_p(:,1)),kind=DP)
     call perform_fft_3d(rhotwr_p,BS_T_grp(i_T_grp)%O_tddft(:,i_T_el,1,1,2))
     !
   endif
   !
   if(n_sp_pol==2) then
     !
     do i_sp_loop=1,n_sp_pol
       do o_sp_loop=1,n_sp_pol
         ! This takes into account the different conventions between the definition of the spin indexes
         ! between yambo and the derivation of fxc in spin space for magnons
         ! This happens because in the yambo convention the left indexes "v sigma_v, c sigma_c, k" are inverted
         ! with respect to the standard convention, i.e. "c sigma_c, v sigma_v, k"
         ! As a result f_\up\dn,\dn\up is in yambo f_\dn\up,\dn\up and so on
         is_yambo=i_sp_loop
         os_yambo=o_sp_loop
         ! To check/fix the spinorial version 
         if (     l_use_rho_map) rhotwr_p(rho_map(:))=cmplx(F_xc_mat(rho_map(:),is_yambo,os_yambo,is(4),os(4))* &
         &                                             WF_symm1_p(:,1)*conjg(WF_symm2_p(:,1)),kind=DP)
         if (.not.l_use_rho_map) rhotwr_p(rho_map(:))=cmplx(F_xc_mat(        : ,is_yambo,os_yambo,is(4),os(4))* &
         &                                             WF_symm1_p(:,1)*conjg(WF_symm2_p(:,1)),kind=DP)
         call perform_fft_3d(rhotwr_p,BS_T_grp(i_T_grp)%O_tddft(:,i_T_el,i_sp_loop,o_sp_loop,2))
         !
       enddo
     enddo
   endif
   !
   if(n_spinor==2) then
     !
     rhotwr_p=cZERO
     do i_spinor=1,n_spinor
       do j_spinor=1,n_spinor
         i_spinor_y=j_spinor
         j_spinor_y=i_spinor
         do ip_spinor=1,n_spinor
           do jp_spinor=1,n_spinor
             do irho=1,rho_map_size
               if(     l_use_rho_map) ifft=rho_map(irho)
               if(.not.l_use_rho_map) ifft=        irho
               rhotwr_p(ifft)=rhotwr_p(ifft)+                        &
               & F_xc_mat(ifft,i_spinor_y,j_spinor_y,ip_spinor,jp_spinor)* &
               & WF_symm1_p(irho,ip_spinor)*conjg(WF_symm2_p(irho,jp_spinor))
             enddo
           enddo
         enddo
         call perform_fft_3d(rhotwr_p,BS_T_grp(i_T_grp)%O_tddft(:,i_T_el,i_spinor,j_spinor,2))
         !
       enddo
     enddo
   endif
   !
 enddo
 !
 nullify(WF_symm1_p)
 nullify(WF_symm2_p)
 nullify(rhotwr_p)
 !
 call timing('T_space ALDA Osc.',OPR='stop')
 !
end subroutine TDDFT_ALDA_eh_space_G_collisions_R
!
!
#include<dev_defs.h>                                                    
!
!
subroutine perform_fft_3d(rhotw,O_x)
 !
 use BS,             ONLY:BS_n_g_fxc
 use pars,           ONLY:cZERO,cONE,DP,SP,DP
 use FFT_m,          ONLY:fftw_plan,cufft_plan,fft_size,fft_dim,DEV_VAR(fft_g_table)
 !
 implicit none
 !
 integer :: ig,iGo
 !
 complex(DP), intent(in)  :: rhotw(fft_size)
 complex(SP), intent(out) :: O_x(BS_n_g_fxc)
 !
 ! perform the actual FFT
 !
#if defined _CUDA
 !
 call fft_3d_cuda(rhotw,fft_dim,+1,cufft_plan)
 !
#else 
 !
#  if defined _FFTW
 call fft_3d(rhotw,fft_dim,+1,fftw_plan)
#  elif defined _FFTSG || (defined _USE_3D_FFT && defined _FFTQE)
 call fft_3d(rhotw,fft_dim,+1)
#  elif defined _FFTQE
 !
 ! the QE fft must be performed in the right
 ! direction. conjg is taken care separately
 ! AF: please note the different definition of isc%rho_tw_rs above
 !     (one extra conjg has been added)
 !
 call fft_3d(rhotw,fft_dim,-1)
 !
#  else
#    error "[CPP] Inconsistent FFT environment"
#  endif
 !
#endif
 !
 iGo=1
 !
#if defined _FFTQE && !defined _USE_3D_FFT && !defined _CUDA_LOC
 !
 ! one extra conjg is performed on rho_tw_rs
 !
 !$omp parallel do default(shared), private(ig)
 do ig = 1, BS_n_g_fxc
   O_x(ig)= cmplx(conjg(rhotw(fft_g_table)(ig,iGo)),kind=SP)
 enddo
 !
#else
 !
#ifdef _CUDA_LOC
 !$cuf kernel do(1) <<<*,*>>>
#else
 !$omp parallel do default(shared), private(ig)
#endif
 do ig = 1, BS_n_g_fxc
   O_x(ig)= cmplx(rhotw(DEV_VAR(fft_g_table)(ig,iGo)),kind=SP)
 enddo
 !
#endif
 !
end subroutine perform_fft_3d
