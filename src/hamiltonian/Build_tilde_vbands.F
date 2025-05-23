!
! License-Identifier: GPL
!
! Copyright (C) 2018 The Yambo Team
!
! Authors (see AUTHORS file for details): MG CA
!
subroutine Build_tilde_vbands(en,Xk,ikbz,i_sp,istep,id_red,V_bands,V_tilde)
 !
 ! Build the \tilde Valence bands as defined in Eq. 53 of
 ! PRB 69, 085106 (2004), in terms of Kohn-Sham orbitals
 !
 use pars,           ONLY:SP,cZERO,cONE
 use electrons,      ONLY:levels
 use electric,       ONLY:Sm1_plus,Sm1_minus
 use QP_m,           ONLY:QP_nk
 use hamiltonian,    ONLY:H_ref_bands,H_ref_nbands
 use R_lattice,      ONLY:bz_samp,k_map
 use vec_operate,    ONLY:k_periodic_idx
 use DIPOLES,        ONLY:DIP_S
 use wrapper_omp,    ONLY:M_by_M_omp
 !
 implicit none
 !
 type(levels), intent(in)      :: en
 type(bz_samp),intent(in)      :: Xk
 integer,      intent(in)      :: ikbz
 integer,      intent(in)      :: i_sp
 integer,      intent(in)      :: istep
 integer,      intent(in)      :: id_red
 complex(SP),  intent(in)      :: V_bands(H_ref_bands(2),en%nbf(i_sp),QP_nk)
 complex(SP),  intent(out)     :: V_tilde(en%nbf(i_sp),H_ref_bands(2),2)
 !
 ! Work Space
 !
 integer :: i2,iv1,im,idx(3),id,ikp,ikm
 integer :: ikm_ibz,ikp_ibz
 complex(SP) :: tmp_Vp(H_ref_bands(2),en%nbf(i_sp)),tmp_Vm(H_ref_bands(2),en%nbf(i_sp))
 !
 V_tilde=cZERO
 !
 if(k_map%max_kdir(id_red)==1) return
 !
 idx=k_map%k_map_inv(ikbz,:)
 idx(id_red)=idx(id_red)-istep
 idx=k_periodic_idx(idx,k_map)
 ikm=k_map%k_map_dir(idx(1),idx(2),idx(3))
 !
 idx=k_map%k_map_inv(ikbz,:)
 idx(id_red)=idx(id_red)+istep
 idx=k_periodic_idx(idx,k_map)
 ikp=k_map%k_map_dir(idx(1),idx(2),idx(3))
 !
 ! Indexes in the IBZ
 ! 
 ikm_ibz=Xk%sstar(ikm,1)
 ikp_ibz=Xk%sstar(ikp,1)
 !
 ! Construct \tilde V_bands(k,k+delta) and \tilde V_bands(k,k-delta)
 !
 ! Build Vp and Vm (stored in tmp_Vp and tmp_Vm)
 !
 ! m,n indexes on valence bands
 !
 ! i,j indexes on all bands CMR basis
 !
 ! Eqs. C11 and C13 in PRB 69, 085106(2004)
 !     
 ! Vp(m,i) = \sum_{j} S(k,K+delta)_{i,j}   V(m,j)
 !
 ! Vm(m,i) = \sum_{j} S^H(k-delta,K)_{j,i} V(m,j)
 !
 ! Build \tilde V(k+delta) - \tilde V(k-delta)
 ! ===========================================
 ! where \tilde V(k+delta,n,:) = \sum_m \tilde S^-1(k,k+delta)_{m,n} * Vp(m,:)
 ! 
 call M_by_M_omp('N','N',H_ref_nbands,en%nbf(i_sp),H_ref_bands(2),cONE,  &
&  DIP_S(H_ref_bands(1):,:,id_red+(istep-1)*3,ikbz,i_sp),H_ref_nbands, &
&  V_bands(:,:en%nbf(i_sp),ikp_ibz),H_ref_bands(2),cZERO,tmp_Vp(H_ref_bands(1):,:en%nbf(i_sp)),H_ref_nbands)
 !
 call M_by_M_omp('C','N',H_ref_nbands,en%nbf(i_sp),H_ref_bands(2),cONE,  &
&  DIP_S(:,H_ref_bands(1):,id_red+(istep-1)*3,ikm ,i_sp),H_ref_bands(2), &
&  V_bands(:,:en%nbf(i_sp),ikm_ibz),H_ref_bands(2),cZERO,tmp_Vm(H_ref_bands(1):,:en%nbf(i_sp)),H_ref_nbands)
 !
 ! Build \tilde V(k+delta) - \tilde V(k-delta)
 ! ===========================================
 ! where \tilde V(k+delta,n,:) = \sum_m \tilde S^-1(k,k+delta)_{m,n} * Vp(m,:)
 ! 
!$omp parallel do default(shared), private(i2,iv1), collapse(2)
 do i2=H_ref_bands(1),H_ref_bands(2)
   do iv1=1,en%nbf(i_sp)
     V_tilde(iv1,i2,1)=sum( Sm1_plus(:en%nbf(i_sp),iv1,id_red+(istep-1)*3,ikbz,i_sp)*tmp_Vp(i2,:en%nbf(i_sp)))
     V_tilde(iv1,i2,2)=sum(Sm1_minus(:en%nbf(i_sp),iv1,id_red+(istep-1)*3,ikbz,i_sp)*tmp_Vm(i2,:en%nbf(i_sp)))
   enddo
 enddo
!$omp end parallel do
 !
end subroutine Build_tilde_vbands
