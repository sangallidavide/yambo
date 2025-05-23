!
! License-Identifier: GPL
!
! Copyright (C) 2015 The Yambo Team
!
! Authors (see AUTHORS file for details): DS AM MG
!
subroutine K_diago_response_functions(iq,i_BS_mat,W,BS_E,BS_R_left_abs,BS_R_right_abs,BS_E_SOC_corr,&
 &                    BS_R_left_magn,BS_R_right_magn,BS_R_left_kerr,BS_R_right_kerr,BS_R_right_dich,BS_R_PL)
 !
 ! eps2(iw)= 1 - Sum  [BS_R_left(i)*BS_R_right(i)] / [w+i*eta - E(i)]  
 !
 use pars,           ONLY:cZERO,cONE,cI,SP,pi
 use units,          ONLY:SPEED_OF_LIGHT
 use X_m,            ONLY:Resp_ii,Joint_DOS,BS_E_sorted,&
&                         global_gauge,i_G_shift
 use electrons,      ONLY:n_sp_pol,n_spinor
 use frequency,      ONLY:w_samp
 use parallel_int,   ONLY:PP_redux_wait,PARALLEL_index
 use parallel_m,     ONLY:myid,PP_indexes,PP_indexes_reset
 use LIVE_t,         ONLY:live_timing
 use BS,             ONLY:BS_perturbative_SOC,BS_res_ares_n_mat,l_BS_abs,l_BS_jdos,l_BS_esort,&
&                         BS_K_anti_resonant,l_BS_ares_from_res,BS_K_coupling
 use BS_solvers,     ONLY:BSS_n_eig,Co_factor
 use R_lattice,      ONLY:bare_qpg,q0_def_norm
 use units,          ONLY:HA2EV
 use vec_operate,    ONLY:degeneration_finder,sort
 use BS_solvers,     ONLY:B_Hall,para_term_w0
 use X_m,            ONLY:Resp_ij,X_magnons,X_dichroism
 use PHOTOLUM,       ONLY:PL,PL_prefactor
 use BS,             ONLY:l_BS_kerr,l_BS_kerr_asymm,l_BS_magnons,l_BS_photolum,l_BS_dichroism
 !
 implicit none
 !
 type(w_samp),intent(in)  :: W
 integer,     intent(in)  :: iq,i_BS_mat
 complex(SP), intent(in)  :: BS_E(BSS_n_eig)
 complex(SP), pointer, intent(in) :: BS_R_left_abs(:),BS_R_right_abs(:)
 real(SP),    pointer, intent(in) :: BS_E_SOC_corr(:,:)
 complex(SP), pointer, intent(in) :: BS_R_left_magn(:,:),BS_R_right_magn(:,:)
 complex(SP), pointer, intent(in) :: BS_R_left_kerr(:),BS_R_right_kerr(:),BS_R_right_dich(:,:)
 real(SP),    pointer, intent(in) :: BS_R_PL(:,:)
 !
 ! Work space
 !
 type(PP_indexes)          :: px
 integer                   :: i1,i_pert_SOC,n_SOC,ares_fac,i_VAR,i_grp,i_lambda,nVAR,EDIM,n_deg_grp 
 complex(SP)               :: EXC_E(2/n_sp_pol),res_abs
 real(SP), allocatable     :: vtmpR(:)
 complex(SP), allocatable  :: BS_E_sorted_tmp(:,:)
 integer, allocatable      :: sort_indx(:),first_exc(:),n_of_exc(:)
 complex(SP)               :: g_fac,res_kerr,res_dich,res_magn(2)
 logical                   :: l_ADD_the_ARES_using_the_RES
 !
 n_SOC=1
 ares_fac=(-1)**(i_BS_mat-1)
 !
 if(l_BS_abs) para_term_w0(2)=cZERO
 !
 if (trim(global_gauge)=='length') g_fac=cONE*Co_factor/bare_qpg(iq,i_G_shift)**2
 !
 if( BS_perturbative_SOC.and.n_sp_pol==1) n_SOC=2
 !
 if(l_BS_esort.and.i_BS_mat==1) then
   EDIM=BSS_n_eig
   nVAR=2
   if(l_BS_abs)       nVAR=nVAR+1
   if(l_BS_kerr)      nVAR=nVAR+1
   if(l_BS_dichroism) nVAR=nVAR+1
   if(l_BS_magnons)   nVAR=nVAR+n_spinor
   allocate(BS_E_sorted_tmp(EDIM,nVAR))
   BS_E_sorted_tmp=cZERO
 endif
 !
 l_ADD_the_ARES_using_the_RES=BS_K_anti_resonant.and.(.not.BS_K_coupling).and.l_BS_ares_from_res
 !
 ! If perturbative SOC energies are loaded compute the SOC energies
 !
 call PP_indexes_reset(px)
 call PARALLEL_index(px,(/BSS_n_eig/))
 !
 call live_timing('BSK resp. funct',px%n_of_elements(myid+1))
 !
 do i1=1,BSS_n_eig
   !
   if (.not.px%element_1D(i1)) cycle
   !
   if (l_BS_abs) then
     res_abs=BS_R_left_abs(i1)*BS_R_right_abs(i1)*ares_fac
   endif
   if (l_BS_kerr) then
     res_kerr=BS_R_left_abs(i1)*BS_R_right_kerr(i1)*ares_fac
     if(l_BS_kerr_asymm) then
       res_kerr=(BS_R_left_abs(i1)*BS_R_right_kerr(i1)-BS_R_left_kerr(i1)*BS_R_right_abs(i1))*ares_fac/2._SP
     endif
   endif
   if(l_BS_dichroism) then
     res_dich=(BS_R_left_abs(i1)*BS_R_right_dich(i1,1)+BS_R_left_kerr(i1)*BS_R_right_dich(i1,2))
     res_dich=res_dich*ares_fac/q0_def_norm/2._SP
   endif
   if(l_BS_magnons) then
     res_magn(1)=BS_R_left_magn(i1,1)*BS_R_right_magn(i1,1)*ares_fac
     res_magn(2)=BS_R_left_magn(i1,2)*BS_R_right_magn(i1,2)*ares_fac
   endif
   !
   if (     BS_perturbative_SOC) EXC_E(:)=BS_E(i1)+BS_E_SOC_corr(:,i1)
   !
   if (.not.BS_perturbative_SOC) EXC_E(1)=BS_E(i1)
   !
   do i_pert_SOC=1,n_SOC
     !
     if(l_BS_esort.and.i_BS_mat==1) then
       BS_E_sorted_tmp(i1,1)=EXC_E(1)
       nVAR=2
     endif
     !
     if (trim(global_gauge)=='velocity') g_fac=Co_factor/EXC_E(i_pert_SOC)**2
     !
     if(l_BS_jdos) then
       ! DS. Note. Jdos with coupling seems to have the wrong ordering, i.e.
       !           Re[Jdos(w)]=-Re[Jdos(-w)] while it should be Re[Jdos(w)]=Re[Jdos(-w)]
       Joint_DOS(:,2)=Joint_DOS(:,2)-ares_fac/(W%p(:)-EXC_E(i_pert_SOC))
       if (l_ADD_the_ARES_using_the_RES) then
         Joint_DOS(:,2)=Joint_DOS(:,2)+1._SP/(W%p(:)+conjg(EXC_E(i_pert_SOC)))
       endif
     endif
     if(l_BS_abs) then
       if(l_BS_esort.and.i_BS_mat==1) then
         BS_E_sorted_tmp(i1,nVAR+1:nVAR+1)=res_abs
         nVAR=nVAR+1
       endif
       Resp_ii(:,2)=Resp_ii(:,2)-res_abs/(W%p(:)-EXC_E(i_pert_SOC))
       if (l_ADD_the_ARES_using_the_RES) then
         Resp_ii(:,2)=Resp_ii(:,2)+conjg(res_abs)/(W%p(:)+conjg(EXC_E(i_pert_SOC)))
       endif
       para_term_w0(2)=para_term_w0(2)+res_abs*g_fac*EXC_E(i_pert_SOC)
     endif
     if(l_BS_kerr) then
       if(l_BS_esort.and.i_BS_mat==1) then
         BS_E_sorted_tmp(i1,nVAR+1:nVAR+1)=res_kerr
         nVAR=nVAR+1
       endif
       Resp_ij(:,2)=Resp_ij(:,2)-res_kerr/(W%p(:)-EXC_E(i_pert_SOC))
       if (l_ADD_the_ARES_using_the_RES) then
         Resp_ij(:,2)=Resp_ij(:,2)+conjg(res_kerr)/(W%p(:)+conjg(EXC_E(i_pert_SOC)))
       endif
       B_Hall(2)=B_Hall(2)+res_kerr*g_fac
     endif
     !
     if(l_BS_dichroism) then
       if(l_BS_esort.and.i_BS_mat==1) then
         BS_E_sorted_tmp(i1,nVAR+1:nVAR+1)=res_dich
         nVAR=nVAR+1
       endif
       X_dichroism(:,2)=X_dichroism(:,2)-res_dich/(W%p(:)-EXC_E(i_pert_SOC))
       if (l_ADD_the_ARES_using_the_RES) then
         X_dichroism(:,2)=X_dichroism(:,2)+conjg(res_dich)/(W%p(:)+conjg(EXC_E(i_pert_SOC)))
       endif
     endif
     !
     if(l_BS_magnons) then
       if(l_BS_esort.and.i_BS_mat==1) then
         BS_E_sorted_tmp(i1,nVAR+1:nVAR+n_spinor)=res_magn(:n_spinor)
         nVAR=nVAR+n_spinor
       endif
       X_magnons(:,1,2)=X_magnons(:,1,2)-res_magn(1)/(W%p(:)-EXC_E(i_pert_SOC))
       X_magnons(:,2,2)=X_magnons(:,2,2)-res_magn(2)/(W%p(:)-EXC_E(i_pert_SOC))
       if (l_ADD_the_ARES_using_the_RES) then
         X_magnons(:,1,2)=X_magnons(:,1,2)+conjg(res_magn(2))/(W%p(:)+conjg(EXC_E(i_pert_SOC)))
         X_magnons(:,2,2)=X_magnons(:,2,2)+conjg(res_magn(1))/(W%p(:)+conjg(EXC_E(i_pert_SOC)))
       endif
     endif
     !
     if(l_BS_photolum) then
       PL(:,2)=PL(:,2)+PL_prefactor*ares_fac*BS_R_PL(1,i1)*aimag(-1._SP/(W%p(:)-EXC_E(i_pert_SOC)))/pi
       if (l_ADD_the_ARES_using_the_RES) then
         PL(:,2)=PL(:,2)+PL_prefactor*BS_R_PL(2,i1)*aimag(-1._SP/(W%p(:)+conjg(EXC_E(i_pert_SOC))))/pi
       endif
     endif
     !
   enddo
   !
   call live_timing(steps=1)
 enddo
 call live_timing
 !
 call PP_indexes_reset(px)
 !
 if(l_BS_esort.and.i_BS_mat==1) then
   call PP_redux_wait(BS_E_sorted_tmp )
   allocate(vtmpR(EDIM)) 
   allocate(sort_indx(EDIM))
   vtmpR=real(BS_E_sorted_tmp(:,1),SP)
   call sort(arrin=vtmpR,indx=sort_indx)
   do i1=1,EDIM
     if(vtmpR(i1)>0._SP) exit
   enddo
   i1=max(1,i1-100)
   EDIM=min(BSS_n_eig-i1+1,10000)
   allocate(first_exc(EDIM)) 
   allocate(n_of_exc(EDIM))
   ! Threshold is 1meV
   call degeneration_finder(vtmpR(i1:i1-1+EDIM),EDIM,first_exc,n_of_exc,n_deg_grp,0.001_SP/HA2EV,Include_single_values=.TRUE.)
   EDIM=min(BSS_n_eig-i1+1,1000)
   do i_VAR=1,nVAR
     BS_E_sorted(:,i_VAR,2)=cZERO
     do i_grp=1,n_deg_grp
       if(i_grp>EDIM) cycle
       if(i_VAR==2) then
         BS_E_sorted(i_grp,i_VAR,2)=real(n_of_exc(i_grp),SP)
         cycle
       endif
       do i_lambda=first_exc(i_grp)+i1-1,first_exc(i_grp)+n_of_exc(i_grp)-1+i1-1
         BS_E_sorted(i_grp,i_VAR,2)=BS_E_sorted(i_grp,i_VAR,2)+BS_E_sorted_tmp(sort_indx(i_lambda),i_VAR)
       enddo
       BS_E_sorted(i_grp,i_VAR,2)=BS_E_sorted(i_grp,i_VAR,2)/real(n_of_exc(i_grp),SP)
     enddo
   enddo
   deallocate(BS_E_sorted_tmp)
 endif
 !
 if (BS_res_ares_n_mat/=i_BS_mat) return
 !
 if(l_BS_jdos)  call PP_redux_wait(Joint_DOS(:,2) )
 !
 if(l_BS_abs) then
   call PP_redux_wait(Resp_ii(:,2))
   call PP_redux_wait(para_term_w0(2))
 endif
 !
 if(l_BS_photolum) then
   call PP_redux_wait(PL(:,2))
 endif
 !
 if(l_BS_kerr) then
   call PP_redux_wait(Resp_ij(:,2))
   call PP_redux_wait(B_Hall(2))
 endif
 !
 if(l_BS_dichroism) then
   call PP_redux_wait(X_dichroism(:,2))
 endif
 !
 if(l_BS_magnons) then
   call PP_redux_wait(X_magnons(:,:,2))
 endif
 !
 if(l_BS_jdos) then
   Joint_DOS(:,2)=Joint_DOS(:,2)*Co_factor !/(4._SP*pi) is needed for the correct Jdos
 !                                         ! However it is commented for now since all test-suite
 !                                         ! references needs to be updated before insterting it
 endif
 !
 if (l_BS_abs) then
   Resp_ii(:,2)=Resp_ii(:,2)*Co_factor
   if (BS_K_anti_resonant.and..not.BS_K_coupling) para_term_w0(2)=para_term_w0(2)+conjg(para_term_w0(2))
 endif
 !
 if (l_BS_kerr) then
   Resp_ij(:,2)=Resp_ij(:,2)*Co_factor
   if (BS_K_anti_resonant.and..not.BS_K_coupling) B_Hall(2)=B_Hall(2)-conjg(B_Hall(2))
 endif
 !
 ! Factors to be fixed
 if (l_BS_dichroism) X_dichroism(:,2)=X_dichroism(:,2)*(cI/SPEED_OF_LIGHT)*Co_factor/(4._SP*pi)
 !
 if (l_BS_magnons) X_magnons(:,:,2)=X_magnons(:,:,2)*Co_factor/(4._SP*pi)
 !
end subroutine K_diago_response_functions
