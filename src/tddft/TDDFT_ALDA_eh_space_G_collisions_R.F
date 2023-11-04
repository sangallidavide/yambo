!
! License-Identifier: GPL
!
! Copyright (C) 2015 The Yambo Team
!
! Authors (see AUTHORS file for details): DS
!
subroutine TDDFT_ALDA_eh_space_G_collisions_R(iq,Xk,i_T_grp,NG,l_bs_tddft_wf_in_loop,tddft_wf,mode)
 !
 ! Calculates the F_xc scattering 
 !
 !  (ic(1),ik(1),is(1)) --<--:...:--<-- (ic(2),ik(2),is(2))
 !                       :sqrt(Fxc):
 !
 use pars,           ONLY:SP,DP,pi,cZERO
 use FFT_m,          ONLY:fft_size
 use wave_func,      ONLY:rho_map_size,rho_map,WF,WAVEs
 use xc_functionals, ONLY:F_xc_mat
 use interfaces,     ONLY:WF_load,WF_free
 use R_lattice,      ONLY:qindx_X,bz_samp,minus_G
 use BS,             ONLY:l_BS_magnons,BS_T_grp,BS_bands
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
 logical,      intent(in) :: l_bs_tddft_wf_in_loop
 character(1), intent(in) :: mode
 type(tddft_wf_t), target, intent(inout) :: tddft_wf 
 !
 ! New implementation. WFs are passed explicitly
 !
 type (WAVEs), save, target :: WF_k,WF_p
 ! 
 ! Work Space
 !
 logical     :: l_use_rho_map,l_load_WFs_k,l_load_WFs_p
 integer     :: i_spinor,j_spinor,ip_spinor,jp_spinor,ifft,irho,irhotw,is_yambo,os_yambo,i_spinor_y,j_spinor_y,&
&               i_sp_loop,o_sp_loop
 integer     :: i_T_el,i_T_el_p,N_T_el_p,i_c,i_v,i_sp_c,i_sp_v,i_k_bz,i_k,i_s,i_g0,i_g1,i_g2,i_p_bz,i_g_p,i_p,&
&               i_k_ref,i_p_ref,i_sp,i_T_grp_p,is(4),os(4),qs(3)
 !
 complex(DP), pointer :: rhotwr_p(:)
 complex(SP), pointer :: WF_symm1_p(:,:)
 complex(SP), pointer :: WF_symm2_p(:,:)
 !
 call timing('T_space ALDA Osc.',OPR='start')
 !
 ! Allocate
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
 !$omp &        i_T_grp_p, N_T_el_p, i_T_el_p, i_g1, i_g2, i_k_ref, i_p_ref, is, os)
#endif
 !
 call OPENMP_update(master_thread)
 !
#if defined(__NOTNOW) && ! defined(_CUDA)
 !$omp do schedule(dynamic)
#endif
 i_k_ref=-1
 i_p_ref=-1
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
   if (mode=="R".or.mode=="Q".or.mode=="F") then
     is=(/i_c,i_k,i_s ,i_sp_c/)
     os=(/i_v,i_p,i_sp,i_sp_v/)
     qs=(/i_g_p,iq,1/)
   else if (mode=="A".or.mode=="C") then
     ! Same momenta and spin indexes, inverted bands
     is=(/i_v,i_k,i_s ,i_sp_c/)
     os=(/i_c,i_p,i_sp,i_sp_v/)
     qs=(/i_g_p,iq,1/)
   endif
   !
   l_load_WFs_k= l_bs_tddft_wf_in_loop .and. i_k_ref/=i_k
   l_load_WFs_p= l_bs_tddft_wf_in_loop .and. i_p_ref/=i_p
   !
   if (l_load_WFs_k) then
#if defined(__NOTNOW) && ! defined(_CUDA)
     !$omp critical
#endif
     if (i_k_ref/=-1) call WF_free(WF_k,keep_fft=.true.,keep_states_to_load=.true.)
     i_k_ref=i_k
     call WF_load(WF_k,NG(1),NG(2),BS_bands,(/i_k,i_k/),quiet=.true.,&
     &            space='R',title="Kernel tddft (k)",keep_states_to_load=.true.)
#if defined(__NOTNOW) && ! defined(_CUDA)
     !$omp end critical
#endif
   endif
   !
   if (l_load_WFs_p) then
#if defined(__NOTNOW) && ! defined(_CUDA)
     !$omp critical
#endif
     if (i_p_ref/=-1) call WF_free(WF_p,keep_fft=.true.,keep_states_to_load=.true.)
     i_p_ref=i_p
     call WF_load(WF_p,NG(1),NG(2),BS_bands,(/i_p,i_p/),quiet=.true.,&
     &            space='R',title="Kernel tddft (p)",keep_states_to_load=.true.)
#if defined(__NOTNOW) && ! defined(_CUDA)
     !$omp end critical
#endif
   endif
   !
   if ( (is(4)/=os(4)) .and. ((.not.l_BS_magnons).or. n_sp_pol==1) ) call error(" TDDFT_ALDA_R_space: wrong spin index")
   if ( (is(4)==os(4)) .and. ((     l_BS_magnons).and.n_sp_pol==2) ) call error(" TDDFT_ALDA_R_space: wrong spin index")
   !
   if (l_bs_tddft_wf_in_loop) then
     call WF_apply_symm_cpu(is,WF_k,WF_symm1_p)
     call WF_apply_symm_cpu(os,WF_p,WF_symm2_p)
   else
     call WF_apply_symm_cpu(is,WF  ,WF_symm1_p)
     call WF_apply_symm_cpu(os,WF  ,WF_symm2_p)
   endif
   !
   if( n_spin==1 ) then
     !
     if (     l_use_rho_map) rhotwr_p(rho_map(:))=cmplx(F_xc_mat(rho_map(:),1,1,1,1)*&
     &                                             WF_symm1_p(:,1)*conjg(WF_symm2_p(:,1)),kind=DP)
     if (.not.l_use_rho_map) rhotwr_p(        : )=cmplx(F_xc_mat(        : ,1,1,1,1)*&
     &                                             WF_symm1_p(:,1)*conjg(WF_symm2_p(:,1)),kind=DP)
     call perform_fft_3d(qs,rhotwr_p,BS_T_grp(i_T_grp)%O_tddft_R(:,i_T_el,1,1))
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
         is_yambo=o_sp_loop
         os_yambo=i_sp_loop
         ! To check/fix the spinorial version 
         if (     l_use_rho_map) rhotwr_p(rho_map(:))=cmplx(F_xc_mat(rho_map(:),is_yambo,os_yambo,is(4),os(4))* &
         &                                             WF_symm1_p(:,1)*conjg(WF_symm2_p(:,1)),kind=DP)
         if (.not.l_use_rho_map) rhotwr_p(        : )=cmplx(F_xc_mat(        : ,is_yambo,os_yambo,is(4),os(4))* &
         &                                             WF_symm1_p(:,1)*conjg(WF_symm2_p(:,1)),kind=DP)
         call perform_fft_3d(qs,rhotwr_p,BS_T_grp(i_T_grp)%O_tddft_R(:,i_T_el,i_sp_loop,o_sp_loop))
         !
       enddo
     enddo
   endif
   !
   if(n_spinor==2) then
     !
     do i_spinor=1,n_spinor
       do j_spinor=1,n_spinor
         i_spinor_y=j_spinor
         j_spinor_y=i_spinor
         rhotwr_p=cZERO
         do ip_spinor=1,n_spinor
           do jp_spinor=1,n_spinor
             do irho=1,rho_map_size
               if(     l_use_rho_map) ifft=rho_map(irho)
               if(.not.l_use_rho_map) ifft=        irho
               rhotwr_p(ifft)=rhotwr_p(ifft)+                        &
               & cmplx(F_xc_mat(ifft,i_spinor_y,j_spinor_y,ip_spinor,jp_spinor)* &
               & WF_symm1_p(irho,ip_spinor)*conjg(WF_symm2_p(irho,jp_spinor)),kind=DP)
             enddo
           enddo
         enddo
         call perform_fft_3d(qs,rhotwr_p,BS_T_grp(i_T_grp)%O_tddft_R(:,i_T_el,i_spinor,j_spinor))
         !
       enddo
     enddo
   endif
   !
 enddo
 !
 if(l_bs_tddft_wf_in_loop) then
   call WF_free(WF_k,keep_fft=.true.,keep_states_to_load=.true.)
   call WF_free(WF_p,keep_fft=.true.,keep_states_to_load=.true.)
 endif
 !
 nullify(WF_symm1_p)
 nullify(WF_symm2_p)
 nullify(rhotwr_p)
 !
 call timing('T_space ALDA Osc.',OPR='stop')
 !
end subroutine TDDFT_ALDA_eh_space_G_collisions_R
