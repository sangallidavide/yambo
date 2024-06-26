!
! License-Identifier: GPL
!
! Copyright (C) 2013 The Yambo Team 
! Copyright (C) 2006 SaX developers team
!
! Authors (see AUTHORS file for details): AF
!
subroutine PP_PWscf_comp(basis,atoms)
  use pw_basis_module
  use pw_wfc_module
  use pw_atoms_module
  use numerical_module
  use num_interpolation_module
  use pw_pseudo_module
  !
  use pars,          ONLY:SP
  use pseudo,        ONLY:pp_kb,pp_kbd,pp_kbs
  !
  implicit none
  !
  type(pw_basis), target, intent(in) :: basis
  type(pw_atoms), intent(in) :: atoms
  type(pw_pseudo), pointer :: psp
  !
  integer :: itype,ibeta,nbeta,ipw,npw,l
  !
  real(SP) :: b(3,3),KpG(3),KpG_cc(3),kb_E
  real(SP) :: KpG_norm(basis%npw)
  real(SP) :: fl(basis%npw)
  real(SP) :: dfl(basis%npw)
  !
  b = basis%struct%b
  npw = basis%npw
  do ipw=1,npw
    KpG = basis%k + real(basis%g(:,ipw),SP)
    KpG_cc = matmul(b,KpG)
    KpG_norm(ipw) = sqrt(sum(KpG_cc(:)**2))
  enddo
  !
  do itype=1,atoms%ntypes
    psp => atoms%pseudo(itype)
    nbeta = psp%nbeta
    do ibeta=1,nbeta
      l=psp%lbeta(ibeta)
      kb_E = psp%d(ibeta)/basis%struct%a_omega*0.5_SP  ! 0.5 = Ry2Ha
      if(abs(kb_E)>0._SP) pp_kbs(itype,ibeta) = kb_E/abs(kb_E)
      !
      fl=0._SP
      dfl=0._SP
      do ipw=1,npw
        if(KpG_norm(ipw) < 1.E-5 .AND. l>0 ) cycle
        fl(ipw) = num_interpolation_calc(psp%interp_beta(ibeta),KpG_norm(ipw),3)
        if(KpG_norm(ipw) < 1.E-5) then
          ! for fl(x->0) ~ \delta_{l,0} * \int_{0}^\infty r beta_pp(r) dr
          pp_kb(ipw,itype,ibeta,1) = sqrt(abs(kb_E))*fl(ipw)
          pp_kbd(ipw,itype,ibeta,1) = 0._SP ! limit goes to zero
        else
          dfl(ipw) = num_interpolation_calc_der(psp%interp_beta(ibeta),KpG_norm(ipw),3,ider=1)
          pp_kb(ipw,itype,ibeta,1) = sqrt(abs(kb_E))*fl(ipw)*KpG_norm(ipw)**l
          pp_kbd(ipw,itype,ibeta,1) = sqrt(abs(kb_E))*(dfl(ipw)*KpG_norm(ipw)**l+ &
&                                         fl(ipw)*real(l,SP)*KpG_norm(ipw)**(l-1) )
        endif
      enddo
    enddo
  enddo
  !
end subroutine PP_PWscf_comp
