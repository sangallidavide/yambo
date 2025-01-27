!
! License-Identifier: GPL
!
! Copyright (C) 2025 The Yambo Team
!
! Authors (see AUTHORS file for details): DS
!
subroutine RT_gradk_rho_times_E(ik,i_sp_pol,H_nl_sc,A_input)
 !
 ! This subroutine computes the k-gradient of the density matrix,
 ! and multiplies it by the electric field.
 ! See Phys. Rev. B 76, 035213 (2007)
 !
 use pars,           ONLY:SP
 use real_time,      ONLY:l_length_grad_k,RT_bands
 use fields,         ONLY:gauge_field
 !
 implicit none
 !
 integer,           intent(in) :: ik,i_sp_pol
 complex(SP),    intent(inout) :: H_nl_sc(RT_bands(1):RT_bands(2),RT_bands(1):RT_bands(2))
 type(gauge_field), intent(in) :: A_input
 !
 if (.not.l_length_grad_k) return
 !
 ikbz=Xk%k_table(ik,1)
 !
 do id_red=1,3  ! loop on reciprocal lattice
   !
   if(k_map%max_kdir(id_red)==1) cycle
   !
   max_step=1
   !
   ! I need at leat 5 k-points in the line to
   ! use the second order formula
   !
   if(k_map%max_kdir(id_red)>=6.and..not.l_force_SndOrd) max_step=2
   !
   do istep=1,max_step 
     !
     idx=k_map%k_map_inv(ikbz,:)
     idx(id_red)=idx(id_red)-istep
     idx=k_periodic_idx(idx,k_map)
     ikm1=k_map%k_map_dir(idx(1),idx(2),idx(3))
     !
     idx=k_map%k_map_inv(ikbz,:)
     idx(id_red)=idx(id_red)+istep
     idx=k_periodic_idx(idx,k_map)
     ikp1=k_map%k_map_dir(idx(1),idx(2),idx(3))
     !
     !S_tmp(:,:,1)=          DIP_S(:,:,id_red+(istep-1)*3,ikbz,i_sp_pol)    ! g+(k)
     !S_tmp(:,:,2)=hermitian(DIP_S(:,:,id_red+(istep-1)*3,ikm1,i_sp_pol))   ! g-(k)=g+(k-q)
     !
     !S_kpq(:,:,1)=          DIP_S(:,:,id_red+(istep-1)*3,ikp1,i_sp_pol)    ! g+(k+q)
     Sm_kpq(:,:)=hermitian(DIP_S(:,:,id_red+(istep-1)*3,ikbz,i_sp_pol))   ! g-(k+q)=g+(k)
     !
     Sp_kmq(:,:)=          DIP_S(:,:,id_red+(istep-1)*3,ikm1,i_sp_pol)    ! g+(k-q)
     !S_kmq(:,:,2)=hermitian(DIP_S(:,:,id_red+(istep-1)*3,ikm2,i_sp_pol))   ! g-(k-q)=g+(k-2*q)
     !
     ! Option 1: Do SVD and later dagger.
     ! Coding to be imported from DIPOLE_build_derk subroutine in branch phys-dipoles-k-derivative
     !
     ! Option 2: Directly do matrix inversion in place of dagger.
     !
     ! To fix 4 indexes on the left, 2x2 on the right
     M_p_kmq(:,:,:,:)=hermitian(Sp_kmq),Sp_kmq
     M_m_kpq(:,:,:,:)=hermitian(Sm_kpq),Sm_kmq
     !
     ! To fix dipole direction
     A1p_kmq=matmul(-cI*dG(:,:,ikm1),hermitina(DIP_iR(i_dir,:,:,ikm1,i_sp_pol)))
     A1m_kpq=matmul(-cI*dG(:,:,ikp1),hermitina(DIP_iR(i_dir,:,:,ikp1,i_sp_pol)))
     !
     A2p_kmq=matmul(-cI*dG(:,:,ikm1),          DIP_iR(i_dir,:,:,ikm1,i_sp_pol)))
     A2m_kpq=matmul(-cI*dG(:,:,ikp1),          DIP_iR(i_dir,:,:,ikp1,i_sp_pol)))
     !
     do ib=RT_bands(1),RT_bands(2)
       do jb=RT_bands(1),RT_bands(2)
         do ibp=RT_bands(1),RT_bands(2)
           do jbp=RT_bands(1),RT_bands(2)
             gradk_rho(ib,jb,id_red,i_step)=gradk_rho(ib,jb,id_red,i_step)+ &
             &  -cI*M_p_kmq(ibp,ib,jbp,jb)*dG(ibp,jbp,ikpq)  &
             &  +cI*M_m_kpq(ibp,ib,jbp,jb)*dG(ibp,jbp,ikmq)  &
             &  -cI*q*M_p_kmq(ibp,ib,jbp,jb)*(A1p_kmq(ibp,ib)-A2p_kmq(ibp,ib)) &
             &  +cI*q*M_m_kpq(ibp,ib,jbp,jb)*(A1m_kpq(ibp,ib)-A2m_kpq(ibp,ib))
           enddo
         enddo
       enddo
     enddo
     !
     ! todo: rotate in cartesian coordinates
     !
   enddo ! istep
 enddo   ! id_red
 ! 
 !
 ! Length gauge
 !===============
 E_vec_pot=-A_input%vecpot_vel/SPEED_OF_LIGHT ! Gaussian units
 !
 do i_dir=1,3
   H_nl_sc(:,:)=H_nl_sc(:,:)-gradk_rho(:,:,i_dir)*E_vec_pot(i_dir)
 enddo
 !
end subroutine RT_grad_k_rho
 