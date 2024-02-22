!
! License-Identifier: GPL
!
! Copyright (C) 2020 The Yambo Team
!
! Authors (see AUTHORS file for details): DS
!
subroutine scatter_Bamp_spin_tddft(isc,which)
 !
 ! It computes
 !
 ! rhotw(G,s,s')=<ib ik_is i_sp| e^[I(RsG-G0).r] S_i F_xc(s,s') | ob ok_os o_sp>    ! which ="i"
 !
 !
 !  is(3)=(ib,ik,is,i_sp) --<--:--<-- os(3)=(ob,ok,os,o_sp)
 !                             :
 !                            /:\ iqs=(ig0,iq,qs)
 !                             :
 !
 !  iq_is = ik_is-ok_os-Go
 !
 use pars,           ONLY:SP,cI,cZERO
 use electrons,      ONLY:n_spinor,n_sp_pol
 use collision_el,   ONLY:elemental_collision
 use interfaces,     ONLY:WF_apply_symm
 use D_lattice,      ONLY:idt_index
 use R_lattice,      ONLY:g_rot
 use FFT_m,          ONLY:fft_size,fft_dim,fftw_plan,fft_g_table
 use xc_functionals, ONLY:F_xc_mat
 !
 implicit none
 !
 type(elemental_collision)::isc
 character(1)             ::which
 !
 ! Work Space
 !
 complex(SP) :: WF_spin(fft_size,n_spinor)
 !
 integer     :: qs,iGo,ig,ir,ic,is(4),os(4),i_sp_loop,o_sp_loop,is_yambo,os_yambo
 !
 if (.not. allocated(isc%rhotw_tddft) ) call error("rhotw_tdddft not allocated")
 !
 ! Here I invert the indexes instead of calling the function with opposite indexes.
 is=isc%os
 os=isc%is
 !
 isc%rhotw_tddft=cZERO
 !
 ! This takes into account the different conventions between the definition of the spin indexes
 ! between yambo and the derivation of fxc in spin space for magnons
 ! This happens because in the yambo convention the left indexes "v sigma_v, c sigma_c, k" are inverted
 ! with respect to the standard convention, i.e. "c sigma_c, v sigma_v, k"
 ! As a result f_\up\dn,\dn\up is in yambo f_\dn\up,\dn\up and so on
 is_yambo=os(4)
 os_yambo=is(4)
 !
 ! |ib ik_is i_sp>
 call WF_apply_symm(is,isc%WF_symm_i)
 !
 ! | ob ok_os o_sp>   
 call WF_apply_symm(os,isc%WF_symm_o)
 !
 if (n_sp_pol==2) then
   if (which=="x") WF_spin=isc%WF_symm_i
   if (which=="y") WF_spin=isc%WF_symm_i*cI*(-1)**(is(4)-1)
   if (which=="z") WF_spin=isc%WF_symm_i   *(-1)**(is(4)-1)
   if (which=="+") WF_spin=isc%WF_symm_i*0.5_SP*real(1+(-1)**is(4),SP)
   if (which=="-") WF_spin=isc%WF_symm_i*0.5_SP*real(1-(-1)**is(4),SP)
 else
   if (which=="x") then
     WF_spin(:,1)=    isc%WF_symm_i(:,2)
     WF_spin(:,2)=    isc%WF_symm_i(:,1)
   endif
   if (which=="y") then
     WF_spin(:,1)= cI*isc%WF_symm_i(:,2)
     WF_spin(:,2)=-cI*isc%WF_symm_i(:,1)
   endif
   if (which=="z") then
     WF_spin(:,1)=    isc%WF_symm_i(:,2)
     WF_spin(:,2)=-   isc%WF_symm_i(:,1)
   endif
   if (which=="+") then
     WF_spin(:,1)=    0._SP
     WF_spin(:,2)=    isc%WF_symm_i(:,1)
   endif
   if (which=="-") then
     WF_spin(:,1)=    isc%WF_symm_i(:,2)
     WF_spin(:,2)=    0._SP
   endif
 endif
 !
 !$omp parallel default(shared), private(ir)
 !$omp do
 do ir = 1, fft_size
   isc%rho_tw_rs(ir) = F_xc_mat(ir,is(4),os(4),is_yambo,os_yambo)*&
   &                             conjg(WF_spin(ir,1))*isc%WF_symm_o(ir,1)
 enddo
 !$omp end do
 !
 if (n_spinor==2) then
   call error("New TDDFT with n_spinor=2 not implemented")
 !$omp do
   do ir = 1, fft_size
     isc%rho_tw_rs(ir) = isc%rho_tw_rs(ir) + conjg(WF_spin(ir,2))*isc%WF_symm_o(ir,2)
   enddo
 !$omp end do
 endif
 !$omp end parallel
 !
 ! perform the actual FFT
 ! N.B. Here there is the "-1" Fourier trasnform to compensate
 !      For the conjg operation done in X_irredux_residauls
 !
#if defined _FFTW
 call fft_3d(isc%rho_tw_rs,fft_dim,-1,fftw_plan)
#else
 call fft_3d(isc%rho_tw_rs,fft_dim,-1)
#endif
 !
 ! e^[I(Rs(G-G0)).r]
 !
 iGo=isc%qs(1)
 qs =isc%qs(3)
 !
#if defined _OPENMP
 if (qs==idt_index) then 
   !$omp parallel do default(shared), private(ig)
   do ig = 1, isc%ngrho
     isc%rhotw_tddft(ig)= isc%rho_tw_rs(fft_g_table(ig,iGo))
   enddo
   !$omp end parallel do
 else
   !$omp parallel do default(shared), private(ig)
   do ig = 1, isc%ngrho
     isc%rhotw_tddft(ig)= isc%rho_tw_rs(fft_g_table(g_rot(qs,ig),iGo))
   enddo
   !$omp end parallel do
 endif
#else
 do ig = 1,isc%ngrho
   if (qs==idt_index) isc%rhotw_tddft(ig)= isc%rho_tw_rs(fft_g_table(ig,iGo))
   if (qs/=idt_index) isc%rhotw_tddft(ig)= isc%rho_tw_rs(fft_g_table(g_rot(qs,ig),iGo))
 enddo
#endif
 !    
end subroutine scatter_Bamp_spin_tddft
