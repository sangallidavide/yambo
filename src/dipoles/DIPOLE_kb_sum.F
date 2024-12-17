!
! License-Identifier: GPL
!
! Copyright (C) 2015 The Yambo Team
!
! Authors (see AUTHORS file for details): AM DS AF
!
#include<dev_defs.h>
!
subroutine DIPOLE_kb_project(pp_range,ib_tot,ib_loc,i_wf,wf_ncx,n_wf,npp,&
&                                     pp_dim_atom,pp_dim_yambo,WF,kbv,XX)
 !
 use iso_c_binding, ONLY:c_f_pointer,c_loc
 use pars,          ONLY:SP,cONE,cZERO
 use electrons,     ONLY:n_spinor
#if defined _CUDAF
 use cudafor,       ONLY:c_devloc
#endif
 use devxlib,       ONLY:devxlib_xGEMM_gpu
 !
 implicit none
 !
 integer,                      intent(in)  :: pp_range(2),ib_tot(2),ib_loc(2),i_wf(2),wf_ncx,n_wf
 integer,                      intent(in)  :: npp,pp_dim_atom,pp_dim_yambo
 complex(SP) DEV_ATTR, target, intent(in)    :: WF(wf_ncx,n_spinor,n_wf)
 complex(SP) DEV_ATTR, target, intent(in)    :: kbv(wf_ncx,n_spinor,4,pp_dim_atom)
 complex(SP) DEV_ATTR, target, intent(inout) :: XX(n_spinor,4,pp_dim_yambo,ib_tot(1):ib_tot(2))
 ! 
 ! Work Space
 !
 integer     :: dim_flat1,dim_flat2,n_bands,ib,i_pp,i1,i_spinor
 complex(SP) DEV_ATTR, target  :: XX_loc(n_spinor,4,pp_range(1):pp_range(2),n_spinor,ib_loc(1):ib_loc(2))
 complex(SP) DEV_ATTR, pointer :: kbv_p(:,:),XX_p(:,:)
 complex(SP) DEV_ATTR, pointer :: WF_p(:,:)

 !
 n_bands=ib_loc(2)-ib_loc(1)+1
 !
 dim_flat1=n_spinor*4*npp
 dim_flat2=n_bands*n_spinor
 !
 ! GPU treatment
 !
#if defined _CUDAF
 !
 call c_f_pointer(c_devloc(WF(1,1,i_wf(1))), WF_p, (/wf_ncx,dim_flat2/) )
 call c_f_pointer(c_devloc(kbv(1,1,1,1)), kbv_p, (/wf_ncx,dim_flat1/) )
 call c_f_pointer(c_devloc(XX_loc(1,1,pp_range(1),1,ib_loc(1))),  XX_p, (/dim_flat1,dim_flat2/) )
 !
#elif defined _OPENACC || defined _OPENMP_GPU
 !
 call c_f_pointer(c_loc(WF(1,1,i_wf(1))), WF_p, (/wf_ncx,dim_flat2/) )
 call c_f_pointer(c_loc(kbv(1,1,1,1)), kbv_p, (/wf_ncx,dim_flat1/) )
 call c_f_pointer(c_loc(XX_loc(1,1,pp_range(1),1,ib_loc(1))),  XX_p, (/dim_flat1,dim_flat2/) )
 !
#else
 !
 ! CPU-only
 !
 call c_f_pointer(c_loc(WF(1,1,i_wf(1))), WF_p, (/wf_ncx,dim_flat2/) )
 call c_f_pointer(c_loc(kbv(1,1,1,1)), kbv_p, (/wf_ncx,dim_flat1/) )
 call c_f_pointer(c_loc(XX_loc(1,1,pp_range(1),1,ib_loc(1))),  XX_p, (/dim_flat1,dim_flat2/) )
 ! 
#endif
 !
 ! AF: the following call falls back on CPU
 !     when GPUs are not available at compile time
 !
 call devxlib_xGEMM_gpu('T','N', dim_flat1, dim_flat2, wf_ncx, cONE, kbv_p, wf_ncx, &
&                                WF_p, wf_ncx, cZERO, XX_p, dim_flat1)
 !
 !DEV_ACC_DEBUG data present(XX,XX_loc)
 !DEV_ACC parallel loop collapse(4)
 !DEV_CUF kernel do(4)
 !DEV_OMPGPU target map(present,alloc:XX,XX_loc)
 !DEV_OMPGPU teams loop collapse(4)
 !DEV_OMP parallel do default(shared), private(ib,i_pp,i1,i_spinor), collapse(4)
 do ib=ib_loc(1),ib_loc(2)
   do i_pp=pp_range(1),pp_range(2)
     do i1=1,4
       do i_spinor=1,n_spinor
         XX(i_spinor,i1,i_pp,ib)=XX_loc(i_spinor,i1,i_pp,i_spinor,ib)
       enddo
     enddo
   enddo
 enddo
 !DEV_OMPGPU end target
 !DEV_ACC_DEBUG end data
 !
end subroutine DIPOLE_kb_project
!
!
subroutine DIPOLE_kb_sum(rho,pp_dim_yambo,ib_range,iv,ic,XX)
 !
 use pars,          ONLY:SP
 use electrons,     ONLY:n_spinor
 use gpu_m,         ONLY:have_gpu
 use devxlib,       ONLY:devxlib_memcpy_d2h
 use DIPOLES,       ONLY:drho=>DEV_VAR(DIP_work_drho)
 !
#include<dev_defs.h>
 !
 implicit none
 !
 complex(SP),          intent(inout) :: rho(3)
 integer,              intent(in)    :: pp_dim_yambo,ib_range(2),iv,ic
 complex(SP) DEV_ATTR, intent(in)    :: XX(n_spinor,4,pp_dim_yambo,ib_range(1):ib_range(2))
 ! 
 ! Work Space
 !
 integer     :: i1,i2,i_spinor
 complex(SP) :: drho_loc(n_spinor,3,pp_dim_yambo)
 !
 if (pp_dim_yambo==0) return
 !
 !DEV_ACC_DEBUG data present(drho,XX)
 !DEV_ACC parallel loop collapse(3)
 !DEV_CUF kernel do(3)
 !DEV_OMPGPU target map(present,alloc:drho,XX)
 !DEV_OMPGPU teams loop collapse(3)
 !DEV_OMP parallel do default(shared), private(i1,i2,i_spinor), collapse(3)
 do i1=1,pp_dim_yambo
   do i2=1,3
     do i_spinor=1,n_spinor
        drho(i_spinor,i2,i1) = conjg(XX(i_spinor,1,i1,iv))*XX(i_spinor,i2+1,i1,ic) &
        &             +conjg(XX(i_spinor,i2+1,i1,iv))*XX(i_spinor,1,i1,ic)
     enddo
   enddo
 enddo
 !DEV_OMPGPU end target
 !DEV_ACC_DEBUG end data
 !
 call devxlib_memcpy_d2h(drho_loc,drho)
 forall(i1=1:3) rho(i1)=rho(i1)+sum(drho_loc(:,i1,:))
 !
end subroutine DIPOLE_kb_sum

