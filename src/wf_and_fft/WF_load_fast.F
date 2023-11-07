!
! License-Identifier: GPL
!
! Copyright (C) 2006 The Yambo Team
!
! Authors (see AUTHORS file for details): DS
!
!> @brief Load wave-function in different forms using the already allocated memory
!!
!! @param[in]       bands_to_load               bands range 
!! @param[in]       kpt_to_load                 k-point to load
!! @param[in]       sp_pol_to_load              spin-polarization range 
!! @param[in]       space                       R = real-space,  G = G-space, C = ?, B = ?
!!
!! @param[out]      WF                  Wave-function object
!
subroutine WF_load_fast(WF,iG_in,iGo_max_in,bands_to_load,kpt_to_load,sp_pol_to_load,isymm)
 !
 ! Load and (eventually) FFTs the wavefunctions
 !
 use pars,          ONLY:SP,DP,schlen,cZERO
 use com,           ONLY:msg,core_io_path,fat_log
 use vec_operate,   ONLY:c2a
 use stderr,        ONLY:intc
 use electrons,     ONLY:nel,n_spinor,n_sp_pol
 use R_lattice,     ONLY:nkibz,k_pt,b
 use interfaces,    ONLY:WF_free,WF_symm_kpoint_r,WF_symm_kpoint_g
 use FFT_m,         ONLY:fft_dim_loaded,fft_dim,fft_norm,fft_g_table,fft_size,&
&                        fftw_plan,cufft_plan,fft_g_table_d
 use wave_func,     ONLY:wf_nc_k,wf_igk,wf_igk_d,rho_map,&
&                        wf_ng,wf_norm_test,wf_ncx,io_WF,wf_b_indx,&
&                        wf_nb_io,wf_nb_to_load,wf_nb_io_groups,WAVEs,states_to_load
 use IO_int,        ONLY:io_control,IO_and_Messaging_switch
 use IO_m,          ONLY:OP_RD,NONE,VERIFY,RD,RD_CL,DUMP
 use wrapper,       ONLY:Vstar_dot_V
#if defined _SC
 use SC,            ONLY:compatible_SC_DB,load_SC_components,SC_bands,found_SC_DB,SC_bands
 use global_XC,     ONLY:WF_kind,WF_xc_functional,WF_perturbation,global_XC_string
#endif
 use pseudo,        ONLY:pp_is_uspp,becp
 use qe_pseudo_m,   ONLY:vkb,tpiba
 use timing_m,      ONLY:timing
 use cuda_m,        ONLY:have_cuda
 !
#include<dev_defs.h>
#include<memory.h>
#include<yambo_wfs.h>
 !
 integer               :: iG_in,iGo_max_in,bands_to_load(2),kpt_to_load,sp_pol_to_load(2),isymm
 type(WAVEs), target   :: WF
 !
 ! Work Space
 !
 integer          ::i_sp_pol,ikibz,ifrag,ib_grp,ib_to_load,ib,ib1,ib2,i_spinor,ic,ig,igfft,&
&                   iG_max,iG_bounds_loaded(2),iGo_max,wf_grid_size,npwk,N
 real(SP)         ::mndp,mxdp,xk(3)
 complex(SP)      ::c
 complex(SP), allocatable :: wf_disk(:,:,:)
 complex(DP), allocatable :: wf_DP(:)
 !
 complex(DP), allocatable DEV_ATTR :: wf_DP_d(:)
 complex(SP), allocatable DEV_ATTR :: wf_disk_d(:,:,:)
 WF_RSPACE(SP),   pointer DEV_ATTR :: wfr_d(:,:,:)
 !
 integer :: WFN
 !
 ! I/O
 !
 integer ::io_err,ID
 !
 ! Close iG/iGo_max to the nearest shell
 !
 iG_max=iG_in
 if (iG_max==0) iG_max=wf_ng
 !
 iGo_max=iGo_max_in
 !
 call PARSER_close_G(iG_max,'tRL')
 call PARSER_close_G(iGo_max,'tRL')
 !
 ! Save dimensions
 !
 WF%b=bands_to_load
 WF%k=(/kpt_to_load,kpt_to_load/)
 WF%sp_pol=sp_pol_to_load
 !
 if(.not.allocated(states_to_load)) &
 & call error("Fast loading of WFs allowed only if states to load was computed")
 !
 ! If the cpu has no WF to load %state must be set to .FALSE. before returning.
 ! This is needed, for example, when in a group of cpus some of them do not 
 ! read the DB
 !
 WF%to_load=count(states_to_load(WF%b(1):WF%b(2),WF%k(1):WF%k(2),WF%sp_pol(1):WF%sp_pol(2)) .eqv. .true.)>0
 if (.not.WF%to_load) then
   YAMBO_FREE(WF%state)
   YAMBO_ALLOC(WF%state,(WF%b(1):WF%b(2),WF%k(1):WF%k(2),WF%sp_pol(1):WF%sp_pol(2)))
   WF%state(WF%b(1):WF%b(2),WF%k(1):WF%k(2),WF%sp_pol(1):WF%sp_pol(2))=.FALSE.
   return
 endif
 !
 if (.not.WF%to_load) return
 !
 ! Memory allocation
 !
 if (WF%space=='R')  wf_grid_size=fft_size
 if (WF%space=='G')  wf_grid_size=wf_ng
 if (WF%space=='C')  wf_grid_size=wf_ncx
 if (WF%space=='B')  wf_grid_size=wf_ncx
 !
 YAMBO_ALLOC(wf_DP,(wf_grid_size))
 !
 if (have_cuda) then
   YAMBO_ALLOC(wf_DP_d,(wf_grid_size))
   wfr_d => WF%r_d
   WFN = WF%N
 endif
 !
 YAMBO_FREE(WF%index)
 YAMBO_ALLOC(WF%index,(WF%b(1):WF%b(2),WF%k(1):WF%k(2),WF%sp_pol(1):WF%sp_pol(2)))
 WF%index=0
 !
 ! actual IO
 !
 YAMBO_ALLOC(wf_b_indx,(wf_nb_io))
 wf_b_indx=0
 !
 call io_control(ACTION=OP_RD,COM=NONE,MODE=VERIFY,SEC=(/1/),ID=ID)
 io_err=io_WF(ID)
 !
 call timing('WF_load_FFT',"start")
 !
 WF%N=0
 do i_sp_pol=1,n_sp_pol
   !
   if (any( (/i_sp_pol<WF%sp_pol(1), i_sp_pol>WF%sp_pol(2)/) )) cycle
   !
   do ikibz=1,nkibz
     !
     if (any( (/ikibz<WF%k(1), ikibz>WF%k(2)/) )) cycle
     !
     ifrag=ikibz+(i_sp_pol-1)*nkibz
     !
     do ib_grp=1,wf_nb_io_groups
       !
       ib1=wf_nb_io*(ib_grp-1)+1
       ib2=wf_nb_io*ib_grp
       !
       if (ib1<WF%b(1)) ib1=WF%b(1)
       if (ib2>WF%b(2)) ib2=WF%b(2)
       !
       wf_nb_to_load=0
       !
       do ib=ib1,ib2
         !
         if ( .not.states_to_load(ib,ikibz,i_sp_pol) ) cycle
         !
         wf_nb_to_load=wf_nb_to_load+1
         wf_b_indx(wf_nb_to_load)=ib
         !
       enddo
       !
       if (wf_nb_to_load==0) cycle
       !
       YAMBO_ALLOC(wf_disk,(wf_ncx,n_spinor,wf_nb_to_load))
       if (have_cuda) then
         YAMBO_ALLOC(wf_disk_d,(wf_ncx,n_spinor,wf_nb_to_load))
       endif
       !
       call io_control(ACTION=RD,COM=NONE,SEC=(/ifrag+1,ib_grp/),ID=ID)
       if (have_cuda) then
         io_err=io_WF(ID,wf_disk,wf_disk_d)
       else
         io_err=io_WF(ID,wf_disk)
       endif
       !
       do ib_to_load=1,wf_nb_to_load
         !
         ib=wf_b_indx(ib_to_load)
         !
         WF%N=WF%N+1
         WF%index(ib,ikibz,i_sp_pol)=WF%N
         !
#ifdef _CUDA
         ! This is a specialized CUDA version performing FFT on the device
         !
         if(WF%space=='R') then
            !
            do i_spinor=1,n_spinor
              !
              wf_DP_d=(0._DP,0._DP)
              npwk=wf_nc_k(ikibz)
              !
              !$cuf kernel do(1) <<<*,*>>>
              do ic=1,npwk
                !
                ig=wf_igk_d(ic,ikibz)
                if (ig>wf_ng) cycle
                igfft=fft_g_table_d(ig,1)
                !
                wf_DP_d(igfft)=cmplx(wf_disk_d(ic,i_spinor,ib_to_load),kind=DP)
                !
              enddo
              !
              call fft_3d_cuda(wf_DP_d,fft_dim,+1,cufft_plan)
              !
              WFN = WF%N
              if(allocated(rho_map)) then
                !$cuf kernel do(1)
                do ic=1,rho_map_size
                  wfr_d(ic,i_spinor,WFN)=wfcmplx(wf_DP_d(rho_map(ic)),kind=SP)*fft_norm
                enddo
              else
                !$cuf kernel do(1)
                do ic=1,fft_size
                  wfr_d(ic,i_spinor,WFN)=wfcmplx(wf_DP_d(ic),kind=SP)*fft_norm
                enddo
              endif
              !
            enddo  ! i_spinor
            !
         else
#endif
            !
            do i_spinor=1,n_spinor
              !
              wf_DP=(0._DP,0._DP)
              do ic=1,wf_nc_k(ikibz)
                ig=wf_igk(ic,ikibz)
                if (ig>wf_ng) cycle
                igfft=ig
                if (WF%space=='R') igfft=fft_g_table(ig,1)
                if (WF%space=='C'.or.WF%space=='B') igfft=ic
                wf_DP(igfft)=cmplx(wf_disk(ic,i_spinor,ib_to_load),kind=DP)
              enddo
              !
              if (WF%space=='G'.or.WF%space=='C'.or.WF%space=='B') then
                WF%c(:,i_spinor,WF%N)=cmplx(wf_DP(:),kind=SP)
                cycle
              endif
              !
#if defined _FFTW
              call fft_3d(wf_DP,fft_dim,+1,fftw_plan)
#else
              !
              ! Here we explot (+2 instead of +1) the wfc cutoff when using FFTQE
              !
              call fft_3d(wf_DP,fft_dim,+2)
#endif
              if(     allocated(rho_map)) WF%r(:,i_spinor,WF%N)=wfcmplx(wf_DP(rho_map(:)),kind=SP)*fft_norm
              if(.not.allocated(rho_map)) WF%r(:,i_spinor,WF%N)=wfcmplx(wf_DP(:)         ,kind=SP)*fft_norm
              !
            enddo  ! i_spinor
#ifdef _CUDA
         endif
#endif
         !
       enddo    ! i_b
       !
       YAMBO_FREE(wf_disk)
       YAMBO_FREE(wf_disk_d)
       !
     enddo      ! i_b_grp
   enddo        ! i_k
 enddo          ! i_sp_pol
 !
 if (have_cuda .and. WF%space/='R') WF%c_d = WF%c
 !
 if (isymm>1) then
   if (WF%space=='R') call WF_symm_kpoint_r(WF%b,WF%k(1),isymm,WF)
   if (WF%space=='G') call WF_symm_kpoint_g(WF%b,WF%k(1),isymm,WF)
 endif
 !
 if (have_cuda .and. WF%space=='R') WF%r = WF%r_d
 !
 call timing('WF_load_FFT',"stop")
 !
 ! CLEAN
 !
 call io_control(ACTION=RD_CL,COM=NONE,MODE=VERIFY,SEC=(/0/),ID=ID)
 io_err=io_WF(ID)
 !
 YAMBO_FREE(wf_b_indx)
 !
 YAMBO_FREE(wf_DP)
 YAMBO_FREE(wf_DP_d)
 !
 YAMBO_FREE(WF%state)
 !
 YAMBO_ALLOC(WF%state,(WF%b(1):WF%b(2),WF%k(1):WF%k(2),WF%sp_pol(1):WF%sp_pol(2)))
 !
 WF%state(WF%b(1):WF%b(2),WF%k(1):WF%k(2),WF%sp_pol(1):WF%sp_pol(2))= &
&         states_to_load(WF%b(1):WF%b(2),WF%k(1):WF%k(2),WF%sp_pol(1):WF%sp_pol(2))
 !
 return
 !
end subroutine WF_load_fast
