!
! License-Identifier: GPL
!
! Copyright (C) 2022 The Yambo Team
!
! Authors (see AUTHORS file for details): DS
! 
function K_correlation_kernel_dir(i_block,i_p,i_pmq,                               &
&                                 i_k_s,i_p_s,i_n_k,i_n_p,i_kmq_s,i_pmq_s,i_m_k,i_m_p,  &
&                                 i_kmq_t,i_pmq_t,i_k_sp_pol_n,i_p_sp_pol_n,i_k_sp_pol_m,i_p_sp_pol_m, &
&                                 iq_W,iq_W_s,ig_W,i_k_s_m1,iq_W_s_mq,ig_W_mq,i_kmq_s_m1)
 !
 use pars,         ONLY:SP,cZERO,pi
 use BS,           ONLY:BS_blk,BS_n_g_W,BS_W,BS_W_is_diagonal
 use D_lattice,    ONLY:sop_inv
 use R_lattice,    ONLY:G_m_G,g_rot,g_phs
#ifdef _CUDA
 use R_lattice,    ONLY:G_m_G_d,g_rot_d,g_phs_d
#endif
 use wrapper,      ONLY:V_by_V_pwise,M_by_V,Vstar_dot_V
 !
 implicit none
 !
 integer,     intent(in) :: i_block,i_p,i_pmq,                                            &
&                                 i_k_s,i_p_s,i_n_k,i_n_p,i_kmq_s,i_pmq_s,i_m_k,i_m_p, &
&                                 i_kmq_t,i_pmq_t,i_k_sp_pol_n,i_p_sp_pol_n,i_k_sp_pol_m,i_p_sp_pol_m, &
&                                 iq_W,iq_W_s,ig_W,i_k_s_m1,iq_W_s_mq,ig_W_mq,i_kmq_s_m1
 !
 complex(SP)  :: K_correlation_kernel_dir
 !
 integer      :: iO1,iO2,i_g,i_gr,i_gs,i_gsr,iq_W_s_m1,iq_W_s_mq_m1
 complex(SP)  :: O1(BS_n_g_W),O2(BS_n_g_W),O_times_W(BS_n_g_W),gphase
#ifdef _CUDA
 attributes(device)   :: O1,O2,O_times_W
 complex(SP), pointer, device :: O_c_d_iO1(:)
 complex(SP), pointer, device :: O_c_d_iO2(:)
#endif
 !
 if (i_k_sp_pol_n/=i_p_sp_pol_n .or. i_k_sp_pol_m/=i_p_sp_pol_m) then
   K_correlation_kernel_dir=cZERO
   return
 endif
 !
 iq_W_s_m1   =sop_inv(iq_W_s)
 iq_W_s_mq_m1=sop_inv(iq_W_s_mq)
 !
 iO1 = BS_blk(i_block)%O_table(i_k_s,i_p_s,1,1,i_n_k,i_n_p,i_k_sp_pol_n)
 iO2 = BS_blk(i_block)%O_table(i_kmq_s,i_pmq_s,i_kmq_t,i_pmq_t,i_m_k,i_m_p,i_k_sp_pol_m)
 !
#ifdef _CUDA
 !
 O_c_d_iO1 => BS_blk(i_block)%O_c_d(:,iO1)
 O_c_d_iO2 => BS_blk(i_block)%O_c_d(:,iO2)
 !$cuf kernel do(1) <<<*,*>>>
 do i_g=1,BS_n_g_W
   i_gs =G_m_G_d(i_g ,ig_W  )
   i_gr =g_rot_d(i_g ,iq_W_s)
   i_gsr=G_m_G_d(i_gr,ig_W  )
   gphase=cmplx(cos(g_phs_d(i_gs,iq_W_s_m1)   ),sin(g_phs(i_gs,iq_W_s_m1)   ))
   O1(i_g)=O_c_d_iO1(i_gsr)*conjg(gphase)
 enddo
 !
 !$cuf kernel do(1) <<<*,*>>>
 do i_g=1,BS_n_g_W
   i_gs =G_m_G_d(i_g ,ig_W_mq  )
   i_gr =g_rot_d(i_g ,iq_W_s_mq)
   i_gsr=G_m_G_d(i_gr,ig_W_mq  )
   gphase=cmplx(cos(g_phs_d(i_gs,iq_W_s_mq_m1)),sin(g_phs(i_gs,iq_W_s_mq_m1)))
   O2(i_g)=O_c_d_iO2(i_gsr)*conjg(gphase)
 enddo
 !
#else
 !
 !$omp parallel do default(shared), private(i_g,i_gs,i_gr,i_gsr,gphase)
 do i_g=1,BS_n_g_W
   i_gs =G_m_G(i_g ,ig_W  )
   i_gr =g_rot(i_g ,iq_W_s)
   i_gsr=G_m_G(i_gr,ig_W  )
   gphase=cmplx(cos(g_phs(i_g,iq_W_s_m1)   ),sin(g_phs(i_g,iq_W_s_m1)   ))
   O1(i_g)=BS_blk(i_block)%O_c(i_gsr,iO1)*conjg(gphase)
 enddo
 !$omp end parallel do
 !
 !$omp parallel do default(shared), private(i_g,i_gs,i_gr,i_gsr,gphase)
 do i_g=1,BS_n_g_W
   i_gs =G_m_G(i_g ,ig_W_mq  )
   i_gr =g_rot(i_g ,iq_W_s_mq)
   i_gsr=G_m_G(i_gr,ig_W_mq  )
   gphase=cmplx(cos(g_phs(i_g,iq_W_s_mq_m1)),sin(g_phs(i_g,iq_W_s_mq_m1)))
   O2(i_g)=BS_blk(i_block)%O_c(i_gsr,iO2)*conjg(gphase)
 enddo
 !$omp end parallel do
 !
#endif
 !
 if (     BS_W_is_diagonal) call V_by_V_pwise(BS_n_g_W,O_times_W,O1,BS_W(:,1,iq_W))
 if (.not.BS_W_is_diagonal) call M_by_V('T',BS_n_g_W,BS_W(:,:,iq_W),O1,O_times_W)
 !
 K_correlation_kernel_dir=Vstar_dot_V(BS_n_g_W,O2,O_times_W)*4._SP*pi
 !
end function K_correlation_kernel_dir
