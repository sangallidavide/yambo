!
! License-Identifier: GPL
!
! Copyright (C) 2006 The Yambo Team
!
! Authors (see AUTHORS file for details): AM AF 
!
#include<dev_defs.h>
!
subroutine X_redux(iq,what,X_par,X_low,X_alda_up,X_alda_dn,Xw,X,call_kind)
 !
 ! Here we solve the Dyson equation using a combination of Serial(SER)/Parallel(PAR)/SLK environments.
 !
 ! More precisely the routine handles:
 !
 !  1. Xo(XUP/XDN/SER): INPUT
 !  2. Dyson(SER/SLK) : SOLVER
 !  3. X(PAR/SER)     : OUTPUT
 !
 ! Therefore the loop is dived in a PAR/SER -> SER/SLK interface, a SER/SLK solver and, finally, a
 ! SER/SLK -> PAR/SER interface
 !
 use pars,          ONLY:SP,pi,cZERO,cONE
 use com,           ONLY:msg
 use drivers,       ONLY:l_bs_fxc,l_alda_fxc,l_lrc_fxc
 use LIVE_t,        ONLY:live_timing
 use parallel_int,  ONLY:PP_wait,PP_redux_wait
 use parallel_m,    ONLY:PP_indexes,PP_indexes_reset,PAR_COM_X_WORLD_RL_resolved,PAR_COM_HOST,&
&                        PAR_COM_X_WORLD,PAR_COM_SLK_INDEX_local,PAR_COM_SLK,PAR_n_freqs,PAR_FREQS_index
 use R_lattice,     ONLY:bare_qpg, bare_qpg_d
 use frequency,     ONLY:w_samp
 use linear_algebra,ONLY:INV,SVD_INV,LIN_SYS,MAT_MUL
 use matrix,        ONLY:PAR_matrix,MATRIX_reset
 use wrapper,       ONLY:M_by_M
 use interfaces,    ONLY:LINEAR_ALGEBRA_driver,MATRIX_transfer
 use stderr,        ONLY:intc
 use X_m,           ONLY:X_t,X_use_lin_sys,X_use_gpu,X_par_lower_triangle,X_redux_build_kernel,&
&                        i_G_shift,l_X_magn
 use TDDFT,         ONLY:FXC_n_g_corr,F_xc_gspace,FXC_LRC_alpha,FXC_LRC_beta,FXC_SVD_digits,&
&                        l_Fxc_in_Xo
 use timing_m,      ONLY:timing
 use devxlib,       ONLY:devxlib_memset_d,devxlib_memcpy_h2d,devxlib_memcpy_d2d,&
                         devxlib_memcpy_d2h
 use devxlib_async, ONLY:devxlib_async_synchronize
 use gpu_m
 !
#include<y_memory.h>
 !
 type(PAR_matrix) :: X_par,X_low,X_alda_up,X_alda_dn
 type(X_t)        :: X
 type(w_samp)     :: Xw
 integer          :: iq
 character(*)     :: call_kind
 character(*)     :: what
 !
 ! Work Space
 !
 type(PP_indexes) :: PAR_IND_freqs
 integer          :: ig1,ig2,iw_par,iw_par_buffer,buffer_nfreq,iw,iw_X(1),INV_MODE,i_LA_pool,&
&                    PAR_n_freqs_global(PAR_COM_SLK_INDEX_local%n_CPU)
 integer          :: Xo_rows(2), Xo_cols(2)
 integer(IPL)     :: HOST_SIZE(3),TMP_SIZE(3)
 logical          :: compute_on_gpu
 type(PAR_matrix) :: K_FXC
 type(PAR_matrix), target      :: BUFFER,BUFFER2,KERNEL,Xo
 complex(SP), allocatable      :: BUFFER_SAVE(:)
 complex(SP), pointer DEV_ATTR :: BUFFER_p(:,:,:)
 complex(SP), pointer DEV_ATTR :: KERNEL_p(:,:,:)
 complex(SP), pointer DEV_ATTR :: Xo_p(:,:,:)
 !
 ! Setup
 !
 call timing(what//' (procedure)',OPR='start')
 !
 ! Barrier...
 !
 call PP_wait(COMM=PAR_COM_X_WORLD%COMM)
 !
 ! compute X on GPU if compiled with GPU support,
 ! unless GPU-usage is switch-off from input
 !
 ! Main logic:
 ! if SLK are present .and. ncpu_SLK >1 .and. CUDA => SLK on CPU
 ! otherwise a problem is obtained since SLK doesn't work on GPU yet
 !
 compute_on_gpu=have_gpu.and.X_use_gpu
 !
 if (PAR_COM_SLK%n_CPU>=4.and.have_gpu) compute_on_gpu=.false.
 !
 ! AF_FIX: openacc + GPU linalg should be enabled
 if (have_openacc) compute_on_gpu=.false.
 !
 ! AMD-GPUs + GPU linalg not yet implemented
 if (have_hip) compute_on_gpu=.false.
 !
 ! Intel-GPUs + GPU linalg not yet implemented
 if (have_mklgpu) compute_on_gpu=.false.
 !
 if (have_gpu) call msg('s','[LA] '//what//' compute_on_gpu ',compute_on_gpu )
 !
 if (.not.X_use_lin_sys) call warning("X inversion implementaiton to be finlized")
 !
 ! ... and frequencies distribution
 !         ========================
 call PARALLEL_FREQS_setup(Xw,PAR_IND_freqs,PAR_COM_SLK_INDEX_local)
 PAR_n_freqs_global =0
 PAR_n_freqs_global(PAR_COM_SLK_INDEX_local%CPU_id+1)=PAR_n_freqs
 call PP_redux_wait(PAR_n_freqs_global,COMM=PAR_COM_SLK_INDEX_local%COMM)
 !
 if (PAR_n_freqs>0) call live_timing(what//'@q['//trim(intc(iq))//'] ',PAR_n_freqs)
 !
 ! ... and local BLACS structures for the response function components
 !         ======================
 call MATRIX_init( "SLK", KERNEL , X%ng, 1 )
 !TMP_SIZE=(/KERNEL%BLCrows(2)-KERNEL%BLCrows(1)+1,KERNEL%BLCcols(2)-KERNEL%BLCcols(1)+1,2/)
 !AMBO_PAR_ALLOC_CHECK3(KERNEL%blc,TMP_SIZE)
 YAMBO_ALLOC(KERNEL%blc,(KERNEL%BLCrows(1):KERNEL%BLCrows(2),KERNEL%BLCcols(1):KERNEL%BLCcols(2),1))
 KERNEL%blc=cZERO
 !
 if(call_kind=="screening")  buffer_nfreq=max(PAR_n_freqs,1)
 if(call_kind=="response" )  buffer_nfreq=1
 !
 call MATRIX_init( "SLK", BUFFER , X%ng, buffer_nfreq )
 !TMP_SIZE=(/BUFFER%BLCrows(2)-BUFFER%BLCrows(1)+1,BUFFER%BLCcols(2)-BUFFER%BLCcols(1)+1,buffer_nfreq*2/)
 !AMBO_PAR_ALLOC_CHECK3(BUFFER%blc,TMP_SIZE)
 YAMBO_ALLOC(BUFFER%blc,(BUFFER%BLCrows(1):BUFFER%BLCrows(2),BUFFER%BLCcols(1):BUFFER%BLCcols(2),buffer_nfreq))
 BUFFER%blc=cZERO
 !
 if (.not.X_use_lin_sys) then
   call MATRIX_init( "SLK", BUFFER2 , X%ng, buffer_nfreq )
   YAMBO_ALLOC(BUFFER2%blc,(BUFFER%BLCrows(1):BUFFER%BLCrows(2),BUFFER%BLCcols(1):BUFFER%BLCcols(2),buffer_nfreq))
   BUFFER2%blc=cZERO
 endif
 !
#if !defined _SCALAPACK
 if (BUFFER%kind=="SLK") call error("[LA] SLK library required but not available")
#endif
 !
 if (compute_on_gpu) then
   YAMBO_ALLOC_GPU(DEV_VAR(KERNEL%blc),(KERNEL%BLCrows(1):KERNEL%BLCrows(2),KERNEL%BLCcols(1):KERNEL%BLCcols(2),1))
   YAMBO_ALLOC_GPU(DEV_VAR(BUFFER%blc),(BUFFER%BLCrows(1):BUFFER%BLCrows(2),BUFFER%BLCcols(1):BUFFER%BLCcols(2),buffer_nfreq))
   call devxlib_memset_d(DEV_VAR(KERNEL%blc),cZERO)
   call devxlib_memset_d(DEV_VAR(BUFFER%blc),cZERO)
 endif
 !
 if(call_kind=="response") then
   buffer_nfreq=max(PAR_n_freqs,1)
   YAMBO_ALLOC(BUFFER_SAVE,(buffer_nfreq))
   BUFFER_SAVE=cZERO
 endif
 !
 if (l_alda_fxc.and.BUFFER%kind=="SLK".and..not.l_Fxc_in_Xo) then
   call MATRIX_init( "SLK", K_FXC , FXC_n_g_corr, 1 )
   YAMBO_ALLOC(K_FXC%blc,(K_FXC%BLCrows(1):K_FXC%BLCrows(2),K_FXC%BLCcols(1):K_FXC%BLCcols(2),1))
   K_FXC%blc=cZERO
   call MATRIX_transfer(M=-F_xc_gspace(:,:,1),M_out=K_FXC)
 endif
 !
 do iw_par=1,maxval(PAR_n_freqs_global)
   !
   if(call_kind=="response" ) iw_par_buffer=1
   if(call_kind=="screening") iw_par_buffer=iw_par
   !
   BUFFER%I   = iw_par_buffer
   BUFFER2%I   = iw_par_buffer
   !
   ! [1] Interface PAR/SER(Xo) => SLK/SER(X_redux)
   ! =============================================
   ! In this delicate part I need to transfer the PAR structure of Xo in the BLACS/serial
   ! by filling the appropriate frequency element. This is why I need to loop on the frequency pools.
   !
   do i_LA_pool=1,PAR_COM_SLK_INDEX_local%n_CPU
     !
     ! Empty pool at this frequency?
     !
     if (iw_par>PAR_n_freqs_global(i_LA_pool)) cycle
     !
     ! The PP_redux forces me to define the iw_X only using the first CPU of each freq pool
     !
     iw_X    = 0
     if (i_LA_pool==PAR_COM_SLK_INDEX_local%CPU_id+1.and.PAR_COM_SLK%CPU_id==0) then
       iw_X  = PAR_FREQS_index(iw_par)
     endif
     call PP_redux_wait(iw_X,COMM=PAR_COM_X_WORLD%COMM)
     !
     ! Only the CPU of this POOL has to store the X matrix (RCV logical).
     !
     call MATRIX_transfer(M_in=X_par,M_out=BUFFER,INDEX_in=iw_X(1),&
&                                    SND=X_par%INTRA_comm%CPU_id==0,&
&                                    RCV=i_LA_pool==PAR_COM_SLK_INDEX_local%CPU_id+1,&
&                                    COMM=PAR_COM_X_WORLD,COMM_name="X_redux_1_"//trim(intc(i_LA_pool)))
     if(allocated(X_low%blc)) then
       call MATRIX_transfer(M_in=X_low,M_out=BUFFER,INDEX_in=iw_X(1),&
&                                      SND=X_par%INTRA_comm%CPU_id==0,&
&                                      RCV=i_LA_pool==PAR_COM_SLK_INDEX_local%CPU_id+1,&
&                                      COMM=PAR_COM_X_WORLD,COMM_name="X_redux_1_"//trim(intc(i_LA_pool)))
     endif
     !
     if (l_alda_fxc.and.l_Fxc_in_Xo) then
       KERNEL%I = iw_par_buffer
       call MATRIX_transfer(M_in=X_alda_up,M_out=KERNEL,INDEX_in=iw_X(1),&
       &                             SND=X_par%INTRA_comm%CPU_id==0,&
       &                             RCV=i_LA_pool==PAR_COM_SLK_INDEX_local%CPU_id+1,&
       &                             COMM=PAR_COM_X_WORLD,COMM_name="X_tddft_redux_1_"//trim(intc(i_LA_pool)))
       if(allocated(X_alda_dn%blc)) then
         call MATRIX_transfer(M_in=X_alda_dn,M_out=KERNEL,INDEX_in=iw_X(1),&
         &                             SND=X_par%INTRA_comm%CPU_id==0,&
         &                             RCV=i_LA_pool==PAR_COM_SLK_INDEX_local%CPU_id+1,&
         &                             COMM=PAR_COM_X_WORLD,COMM_name="X_tddft_redux_1_"//trim(intc(i_LA_pool)))
       endif
     endif
     !
   enddo
   !
   if (iw_par>PAR_n_freqs) cycle
   !
   ! [2] SER/SLK Solver
   ! ==================
   !
   iw = PAR_FREQS_index(iw_par)
   BUFFER_p => DEV_VAR(BUFFER%blc)
   !
   if (.not.(l_alda_fxc.and.l_Fxc_in_Xo)) then
     KERNEL_p => DEV_VAR(KERNEL%blc)
     if (compute_on_gpu) call devxlib_memset_d(KERNEL_p,cZERO)
   endif
   !
   Xo_cols(:)=BUFFER%cols(:)
   Xo_rows(:)=BUFFER%rows(:)
   !
   if (FXC_SVD_digits>0) INV_MODE=SVD_INV
   !
   ! TDDFT Kernel. Different procedure depending on the kernel
   !
   ! Kind: BS,ALDA,LRC.
   !
   if (l_bs_fxc) then
     !
   else if (l_lrc_fxc) then
     !
     ! LRC Fxc
     !
     if (BUFFER%cols(1)==1) then
       !$omp parallel do default(shared), private(ig1)
       do ig1=BUFFER%rows(1),BUFFER%rows(2)
         KERNEL%blc(ig1,1,1)=-BUFFER%blc(ig1,1,BUFFER%I)*(FXC_LRC_alpha + FXC_LRC_beta*abs(Xw%p(iw))**2)/bare_qpg(iq,1)**2
       enddo
       !$omp end parallel do
     endif
     !
   else if (l_alda_fxc.and..not.l_FXC_in_Xo) then
     !
     ! ALDA Fxc
     !
     if (BUFFER%kind=="SLK") then
       call PARALLEL_M_by_M(BUFFER,K_FXC,KERNEL,"N","N",X%ng,FXC_n_g_corr,FXC_n_g_corr)
     else
       call M_by_M('N','N',X%ng,FXC_n_g_corr,FXC_n_g_corr,-cONE,BUFFER%blc(:FXC_n_g_corr,:FXC_n_g_corr,BUFFER%I),X%ng,&
       &                  F_xc_gspace(:,:,1),FXC_n_g_corr,cZERO,KERNEL%blc(:FXC_n_g_corr,:FXC_n_g_corr,1),       X%ng )
     endif
     !
   endif
   !
   if (compute_on_gpu) then
     if (l_bs_fxc.or.l_lrc_fxc.or.l_alda_fxc) call devxlib_memcpy_h2d(KERNEL_p,KERNEL%blc)
     ! This may not be needed since BUFFER is already sincronized later
     ! It was needed previously when using Xo
     call devxlib_memcpy_h2d(BUFFER_p(:,:,BUFFER%I),BUFFER%blc(:,:,BUFFER%I))
   endif
   !
   ! A this point
   ! - BUFFER contains Xo
   ! - KERNEL is either 0 [ no tddft] or -[Xo*F_xc](g1,g2) [tddft]
   !
   ! The following call adds to KERNEL (delta_(g1,g2)-Xo(g1,g2)*v(g2))
   call X_redux_build_kernel(KERNEL,BUFFER,Xo_rows,Xo_cols,iw_par_buffer,compute_on_gpu,iq)
   !
   ! X(g,gp)=Sum_gpp[KERNEL]^-1_(g,gpp)*Xo(gpp,gp)
   !
   ! on input:  Xo irreducible polarizability
   ! on output: Xo is overwritten
   !
   if (X_use_lin_sys) then
     !
     ! Linear System
     !
     call LINEAR_ALGEBRA_driver(LIN_SYS,M_slk=KERNEL,B_slk=BUFFER)
     !
   else
     !
     ! Matrix Inversion + Matmul
     !
     if (compute_on_gpu) call error("[LA] inversion+matmul not available on GPUs")
     !
     call LINEAR_ALGEBRA_driver(INV,M_slk=KERNEL)
     !
     BUFFER2%blc=cZERO
     call LINEAR_ALGEBRA_driver(MAT_MUL,M_slk=KERNEL,B_slk=BUFFER,C_slk=BUFFER2)
     do ig2=BUFFER%cols(1),BUFFER%cols(2)
       do ig1=BUFFER%rows(1),BUFFER%rows(2)
         BUFFER%blc(ig1,ig2,BUFFER%I)=BUFFER2%blc(ig1,ig2,BUFFER%I)
       enddo
     enddo
     !
   endif
   !
   ! When calculating EPS1 multiply X left and right by v_coul^1/2 (diagonal in reciprocal space)
   !
   if (X%whoami/=1 .and. .not.l_X_magn) then
     !
     if (compute_on_gpu) then
        !
        !DEV_ACC_DEBUG data present(Xo_p,bare_qpg)
        !DEV_ACC parallel loop collapse(2)
        !DEV_CUF kernel do(2)
        !DEV_OMPGPU target map(present,alloc:Xo_p,bare_qpg)
        !DEV_OMPGPU teams loop collapse(2)
        do ig2=Xo_cols(1),Xo_cols(2)
          do ig1=Xo_rows(1),Xo_rows(2)
            BUFFER_p(ig1,ig2,BUFFER%I)=BUFFER_p(ig1,ig2,BUFFER%I)*4._SP*pi/DEV_VAR(bare_qpg)(iq,ig1)/DEV_VAR(bare_qpg)(iq,ig2)
          enddo
        enddo
        !DEV_OMPGPU end target
        !DEV_ACC_DEBUG end data
        !
        !call devxlib_async_synchronize(async_d2h)
        call devxlib_memcpy_d2h(BUFFER%blc(:,:,BUFFER%I),BUFFER_p(:,:,BUFFER%I), async_d2h)
        !
     else
        !
        !$omp parallel do default(shared), private(ig1,ig2), collapse(2)
        do ig2=BUFFER%cols(1),BUFFER%cols(2)
          do ig1=BUFFER%rows(1),BUFFER%rows(2)
            BUFFER%blc(ig1,ig2,BUFFER%I)=BUFFER%blc(ig1,ig2,BUFFER%I)*4._SP*pi/bare_qpg(iq,ig1)/bare_qpg(iq,ig2)
          enddo
        enddo
        !$omp end parallel do
        !
     endif
     !
   endif
   !
   if (PAR_n_freqs>0) call live_timing(steps=1)
   !
   if(call_kind=="response") BUFFER_SAVE(iw_par)=BUFFER%blc(i_G_shift,i_G_shift,BUFFER%I)
   !
 enddo
 !
 !call devxlib_async_synchronize(async_d2h)
 !
 if (PAR_n_freqs>0) call live_timing()
 !
 ! Clean-up of local BLACS structures
 ! ========
 call MATRIX_reset(KERNEL)
 if (l_alda_fxc.and.BUFFER%kind=="SLK") call MATRIX_reset(K_FXC)
 !
 call timing(what//' (procedure)',OPR='stop')
 !
 if(call_kind=="screening") call X_redux_X_par_micro(X,Xw,what,X_par,X_low,BUFFER,PAR_n_freqs_global)
 !
 if(call_kind=="response")  call X_redux_X_par_macro(Xw,X_par,BUFFER_SAVE,buffer_nfreq,PAR_n_freqs_global)
 !
 ! Clean-up of FREQ structure
 !==========
 call MATRIX_reset(BUFFER)
 YAMBO_FREE(PAR_FREQS_index)
 call PP_indexes_reset(PAR_IND_freqs)
 if(call_kind=="response") then
   YAMBO_FREE(BUFFER_SAVE)
 endif
 !
end subroutine
   !
   !
   subroutine X_redux_build_kernel(KERNEL,Xo,Xo_rows,Xo_cols,iw_par_buffer,compute_on_gpu,iq)
     use R_lattice,     ONLY:bare_qpg, DEV_VAR(bare_qpg)
     use pars,          ONLY:pi,cONE,SP
     use matrix,        ONLY:PAR_matrix
     use X_m,           ONLY:l_X_magn
     implicit none
     type(PAR_matrix), target :: KERNEL,Xo
     integer, intent(in)      :: Xo_rows(2),Xo_cols(2),iq,iw_par_buffer
     logical, intent(in)      :: compute_on_gpu
     !
     integer :: ig1,ig2
     !
     complex(SP), pointer DEV_ATTR :: KERNEL_p(:,:,:)
     complex(SP), pointer DEV_ATTR :: Xo_p(:,:,:)
     !
     if (compute_on_gpu) then
       !
       KERNEL_p => DEV_VAR(KERNEL%blc)(:,:,:)
       Xo_p     => DEV_VAR(Xo%blc)(:,:,:)
       !
       !DEV_ACC_DEBUG data present(KERNEL_p,Xo_p,bare_qpg)
       !DEV_ACC parallel loop collapse(2)
       !DEV_CUF kernel do(2) <<<*,*>>>
       !DEV_OMPGPU target map(present,alloc:KERNEL_p,Xo_p,bare_qpg)
       !DEV_OMPGPU teams loop collapse(2)
       do ig2=Xo_cols(1),Xo_cols(2)
         do ig1=Xo_rows(1),Xo_rows(2)
           if (.not.l_X_magn) KERNEL_p(ig1,ig2,1)=KERNEL_p(ig1,ig2,1) &
           &                  -Xo_p(ig1,ig2,iw_par_buffer)*4._SP*pi/DEV_VAR(bare_qpg)(iq,ig2)**2
           if (ig1==ig2) KERNEL_p(ig1,ig1,1)=KERNEL_p(ig1,ig1,1)+cONE
         enddo
       enddo
       !DEV_OMPGPU end target
       !DEV_ACC_DEBUG end data
       !
     else
       !
       !$omp parallel do default(shared), private(ig1,ig2), collapse(2)
       do ig2=Xo%cols(1),Xo%cols(2)
         do ig1=Xo%rows(1),Xo%rows(2)
           if (.not.l_X_magn) KERNEL%blc(ig1,ig2,1) =KERNEL%blc(ig1,ig2,1) &
           &                      -Xo%blc(ig1,ig2,iw_par_buffer)*4._SP*pi/bare_qpg(iq,ig2)**2
           if (ig1==ig2) KERNEL%blc(ig1,ig1,1)=KERNEL%blc(ig1,ig1,1)+cONE
         enddo
       enddo
       !
     endif
     !
   end subroutine X_redux_build_kernel
   !
   !
   subroutine X_redux_X_par_micro(X,Xw,what,X_par,X_low,BUFFER,PAR_n_freqs_global)
     use matrix,        ONLY:PAR_matrix,MATRIX_reset
     use frequency,     ONLY:w_samp
     use X_m,           ONLY:X_t
     use stderr,        ONLY:intc
     use timing_m,      ONLY:timing
     use parallel_int,  ONLY:PP_redux_wait,PP_bcast
     use interfaces,    ONLY:MATRIX_transfer
     use parallel_m,    ONLY:PAR_COM_X_WORLD_RL_resolved,PAR_COM_X_WORLD,&
&                            PAR_COM_SLK_INDEX_local,PAR_COM_SLK,PAR_n_freqs,PAR_FREQS_index
     !
     implicit none
     !
     type(PAR_matrix) :: X_par,X_low,BUFFER
     type(X_t)        :: X
     type(w_samp)     :: Xw
     integer          :: PAR_n_freqs_global(PAR_COM_SLK_INDEX_local%n_CPU)
     character(*)     :: what
     !
     ! Local variables
     !
     integer :: i_LA_pool,iw,iw_par,iw_X(1)
     !
     call timing(what//' (REDUX)',OPR='start')
     !
     ! [3] Interface SER/SLK(X_redux) => SER/PAR(X)
     ! ============================================
     !
     ! The old X_par structure need to be cleaned and replaced by a non symmetric structure
     !
     if (X_par%INTER_comm%n_CPU>1) call X_ALLOC_parallel(X_par,X_low,X%ng,Xw%n_freqs,"X,CPU_ONLY")
     !
     ! Now the interfaces. 
     !
     ! From all CPU of the SLK group to one of the RL group.
     !
     do i_LA_pool=1,PAR_COM_SLK_INDEX_local%n_CPU
       do iw_par=1,PAR_n_freqs_global(i_LA_pool)
         !
         ! The PP_redux forces me to define the iw_X only using the first CPU of each freq pool
         !
         iw_X    = 0
         if (i_LA_pool==PAR_COM_SLK_INDEX_local%CPU_id+1.and.PAR_COM_SLK%CPU_id==0) then
           iw_X  = PAR_FREQS_index(iw_par)
         endif
         call PP_redux_wait(iw_X,COMM=PAR_COM_X_WORLD%COMM)
         !
         ! iw_par is used ONLY by the CPU's that SEND the matrix (i_LA_pool==PAR_COM_SLK_INDEX_local%CPU_id+1)
         ! All other CPU's do not use the BUFFER and can set INDEX_in as a dummy 1. This is needed as 
         ! BUFFER is dimensioned with the PAR #of freqs.
         !
         iw=iw_par
         if (iw_par>PAR_n_freqs) iw=1
         !
         call MATRIX_transfer(M_in=BUFFER,M_out=X_par,INDEX_in=iw,INDEX_out=iw_X(1),&
&                             SND=i_LA_pool==PAR_COM_SLK_INDEX_local%CPU_id+1,&
&                             RCV=PAR_COM_X_WORLD_RL_resolved%CPU_id==0,&
&                             COMM=PAR_COM_X_WORLD,COMM_name="X_redux_2_"//trim(intc(i_LA_pool)))
         !
       enddo
     enddo
     !
     ! complete the transfer with a collective communication
     !
     call PP_bcast(X_par%blc,node=0,COMM=PAR_COM_X_WORLD_RL_resolved%COMM)
     !
     call timing(what//' (REDUX)',OPR='stop')
     !
   end subroutine X_redux_X_par_micro
   !
   !
   subroutine X_redux_X_par_macro(Xw,X_par,BUFFER_SAVE,buffer_nfreq,PAR_n_freqs_global)
     use pars,          ONLY:cZERO,SP
     use matrix,        ONLY:PAR_matrix,MATRIX_reset
     use frequency,     ONLY:w_samp
     use parallel_int,  ONLY:PP_redux_wait
     use parallel_m,    ONLY:PAR_COM_X_WORLD_RL_resolved,&
&                        PAR_COM_SLK_INDEX_local,PAR_COM_SLK,PAR_n_freqs,PAR_FREQS_index
     !
#include<y_memory.h>
     !
     type(PAR_matrix) :: X_par
     type(w_samp)     :: Xw
     integer          :: buffer_nfreq
     complex(SP)      :: BUFFER_SAVE(buffer_nfreq)
     integer          :: PAR_n_freqs_global(PAR_COM_SLK_INDEX_local%n_CPU)
     !
     ! Local variables
     !
     integer :: i_LA_pool,iw,iw_par,iw_X(1)
     !
     call MATRIX_reset(X_par)
     YAMBO_ALLOC(X_par%blc,(1,1,Xw%n_freqs))
     X_par%blc=cZERO
     !
     do i_LA_pool=1,PAR_COM_SLK_INDEX_local%n_CPU
       do iw_par=1,PAR_n_freqs_global(i_LA_pool)
         if (iw_par>PAR_n_freqs) cycle
         !
         ! The PP_redux forces me to define the iw_X only using the first CPU of each freq pool
         !
         iw_X    = 0
         if (i_LA_pool==PAR_COM_SLK_INDEX_local%CPU_id+1.and.PAR_COM_SLK%CPU_id==0) then
           iw_X  = PAR_FREQS_index(iw_par)
           iw=iw_par
           X_par%blc(1,1,iw_X(1))=BUFFER_SAVE(iw)
         endif
         !
       enddo
     enddo
     !
     call PP_redux_wait(X_par%blc,COMM=PAR_COM_X_WORLD_RL_resolved%COMM)
     !
   end subroutine X_redux_X_par_macro
