!
! License-Identifier: GPL
!
! Copyright (C) 2006 The Yambo Team
!
! Authors (see AUTHORS file for details): AM AF
!
!> @brief Load wave-function in different forms
!!
!! @param[in]       iG_in                       max-number of G-vectors
!! @param[in]       iGo_max_in                  ???
!! @param[in]       bands_to_load               bands range 
!! @param[in]       kpts_to_load                k-points range 
!! @param[in]       k_extrema_only              k-points, load only !first and last in range
!! @param[in]       sp_pol_to_load              spin-polarization range 
!! @param[in]       space                       R = real-space,  G = G-space, C = ?, B = ?
!! @param[in]       title                       purpose of the wave-function loading
!! @param[in]       impose_free_and_alloc       force deallocation of already loaded WF 
!! @param[in]       force_WFo                   do not rotate the WF (for yambo_sc)
!! @param[in]       keep_states_to_load         keep (.TRUE.) or not (.FALSE.) states already loaded 
!! @param[in]       quiet                       turn off report and log for the WF_laod 
!!
!! @param[out]      WF                  Wave-function object
!
subroutine WF_load(WF,iG_in,iGo_max_in,bands_to_load,kpts_to_load,sp_pol_to_load,k_extrema_only,&
&                  space,title,impose_free_and_alloc,force_WFo,keep_states_to_load,quiet)
 !
 ! Load and (eventually) FFTs the wavefunctions
 !
 use pars,          ONLY:SP,DP,schlen,cZERO
 use com,           ONLY:msg,core_io_path,fat_log
 use vec_operate,   ONLY:c2a
 use stderr,        ONLY:intc
 use electrons,     ONLY:nel,n_spinor,n_sp_pol
 use R_lattice,     ONLY:nkibz,k_pt,b
 use interfaces,    ONLY:WF_free
 use FFT_m,         ONLY:fft_dim_loaded,fft_dim,fft_norm,fft_g_table,fft_size,&
&                        fftw_plan,cufft_plan,fft_g_table_d
 use wave_func,     ONLY:wf_nc_k,wf_igk,wf_igk_d,WF_buffer,&
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
 !
 integer               :: iG_in,iGo_max_in,bands_to_load(2),kpts_to_load(2)
 type(WAVEs), target   :: WF
 integer,     optional :: sp_pol_to_load(2)
 character(*),optional :: space
 character(*),optional :: title
 logical,     optional :: k_extrema_only
 logical,     optional :: impose_free_and_alloc
 logical,     optional :: force_WFo
 logical,     optional :: keep_states_to_load
 logical,     optional, intent(in) :: quiet
 !
 ! Work Space
 !
 character(1)     ::wf_space_here
 character(schlen)::wf_title
 integer          ::i_sp_pol,ikibz,ifrag,ib_grp,ib_to_load,ib,ib1,ib2,i_spinor,ic,ig,igfft,&
&                   iG_max,iG_bounds_loaded(2),iGo_max,sp_pol_to_load_(2),wf_grid_size,&
&                   npwk,N
 real(SP)         ::mndp,mxdp,xk(3)
 complex(SP)      ::c
 logical          ::loaded_bounds_ok,free_and_alloc,buffer_is_ok,clean_up_states,&
                    k_extrema_only_,force_WFo_,QUIET_wf,SIZE_msg
 complex(SP), allocatable :: wf_disk(:,:,:)
 complex(DP), allocatable :: wf_DP(:)
 !
 complex(DP), allocatable DEV_ATTR :: wf_DP_d(:)
 complex(SP), allocatable DEV_ATTR :: wf_disk_d(:,:,:)
 complex(SP),     pointer DEV_ATTR :: wfc_d(:,:,:)
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
 wf_space_here='R'
 if (present(space)) wf_space_here=space
 !
 wf_title=""
 if (present(title)) wf_title=title
 !
 force_WFo_=.false.
 if(present(force_WFo)) force_WFo_=force_WFo
 !
 k_extrema_only_=.false.
 if(present(k_extrema_only)) k_extrema_only_=k_extrema_only
 !
 sp_pol_to_load_=(/1,n_sp_pol/)
 if(present(sp_pol_to_load)) sp_pol_to_load_=sp_pol_to_load
 !
 ! Imposed Freeing and Allocating of WF's...
 !
 free_and_alloc=.FALSE.
 !... when a single k-point is loaded each time... (no _SC)
#if !defined _SC
 if (kpts_to_load(1)==kpts_to_load(2)) free_and_alloc =kpts_to_load(1)==1
#endif
 !... or when it is externally imposed
 if (present(impose_free_and_alloc)) then
   free_and_alloc=impose_free_and_alloc
 endif
 !
 ! [1]: check that loaded bounds are larger(equal) then bounds 
 !      asked now
 !
 call check_bounds(WF,loaded_bounds_ok)
 if (loaded_bounds_ok) loaded_bounds_ok= WF%space==wf_space_here
 !
 call check_bounds(WF_buffer,buffer_is_ok)
 !
 buffer_is_ok=.FALSE.
 !
 if (allocated(states_to_load)) then
   !
   if (loaded_bounds_ok) then
     !
     ! See if the %c in WF is already ok
     !
     if (allocated(WF%state)) call check_states(WF,loaded_bounds_ok)
     !
   endif
   !
   if (allocated(WF_buffer%state).and.wf_space_here/='B') then
     !
     ! See if the %c Buffered are ok
     !
     buffer_is_ok=.TRUE.
     !
     call check_states(WF_buffer,buffer_is_ok)
     !
     ! if not clean the buffer
     !
     if (.not.buffer_is_ok) call WF_free(WF_buffer)
     !
   endif
   !
 endif
 !
 ! [2]: Check FFT size
 !
 if (loaded_bounds_ok.and..not.free_and_alloc) then
   !
   if (WF%space=="G".or.WF%space=="C".or.WF%space=="B") return
   !
   call fft_setup(iG_max,iGo_max,.true.)
   !
   if (all(fft_dim<=fft_dim_loaded)) then
     iG_bounds_loaded=shape(fft_g_table)
     if (iG_bounds_loaded(1)>=iG_max.and.&
&        iG_bounds_loaded(2)>=iGo_max) then
       !
       ! Reset dimensions to fft dim loaded 
       !
       fft_dim=fft_dim_loaded
       fft_size=product(fft_dim)
       return
     endif
   endif
   !
 endif
 !
 ! In case of k by k I/O prevent multiple free/alloc
 !
 QUIET_wf = .false.
 if (kpts_to_load(1)==kpts_to_load(2).and.kpts_to_load(2)>1) QUIET_wf = .true.
 if (k_extrema_only_.and.kpts_to_load(1)>1)                  QUIET_wf = .true.
 if (present(quiet))                                         QUIET_wf = quiet
 if (QUIET_wf) call IO_and_Messaging_switch("-report -log")
 SIZE_msg=.not.QUIET_wf.and.(any((/fft_dim_loaded/=fft_dim/)).or.any((/fft_dim==0/)))
 !
 ! First WF_free then fft_setup.
 !
 call WF_free(WF)
 !
 ! Save dimensions
 !
 WF%b=bands_to_load
 WF%k=kpts_to_load
 WF%sp_pol=sp_pol_to_load_
 WF%space=wf_space_here
 !
 ! If the cpu has no WF to load %state must be set to .FALSE. before returning.
 ! This is needed, for example, when in a group of cpus some of them do not 
 ! read the DB
 !
 if (allocated(states_to_load)) then
   WF%to_load=count(states_to_load(WF%b(1):WF%b(2),WF%k(1):WF%k(2),WF%sp_pol(1):WF%sp_pol(2)) .eqv. .true.)>0
   if (.not.WF%to_load) then
     YAMBO_FREE(WF%state)
     YAMBO_ALLOC(WF%state,(WF%b(1):WF%b(2),WF%k(1):WF%k(2),WF%sp_pol(1):WF%sp_pol(2)))
     WF%state(WF%b(1):WF%b(2),WF%k(1):WF%k(2),WF%sp_pol(1):WF%sp_pol(2))=.FALSE.
     YAMBO_FREE(states_to_load)
   endif
 endif
 !
 ! AM 28 Sept 2016. When no states need to be loaded still the calculation
 ! is not stopped. If we are in an HF calculation also the XC potential is calculated.
 ! In that case the fft_size must be defined for all cpu's of the density COMM otherwise
 ! the PP_redux in el_density hungs up forever
 !
 if (WF%space=='R') call fft_setup(iG_max,iGo_max,.false.)
 if (WF%space=="R") fft_dim_loaded=fft_dim
 !
 if (.not.WF%to_load) return
 !
#if defined _SC
 !
 ! First check and fix the number of bands 
 !
 call load_SC_components('check_WF_DB',n_bands=bands_to_load(2))
 !
 if (found_SC_DB.and.compatible_SC_DB) then
   !
   if(allocated(states_to_load).and..not.force_WFo_) then
     YAMBO_FREE(states_to_load)
     call warning(' [WF'//trim(wf_title)//'] WF-distribution turned off, not compatible with SC wave functions')
   endif
   !
   if (any ( WF%b/=(/1,SC_bands(2)/) )) then 
     call warning(' [WF'//trim(wf_title)//'] Forced '//trim(intc(SC_bands(2)))//' bands to be loaded')
     if(allocated(states_to_load)) states_to_load(WF%b(2)+1:SC_bands(2),:,:)=.TRUE.
   endif
   WF%b=(/1,SC_bands(2)/)
   !
 else if (found_SC_DB.and..not.compatible_SC_DB) then
   !
   call error(' [WF'//trim(wf_title)//'] Too many bands required ('//&
&             trim(intc(bands_to_load(2)))//'). Only '//trim(intc(SC_bands(2)))//' available.')
   !
 endif
 !
#endif
 !
 ! Init USPP if needed
 !
 call PP_uspp_init()
 !
 ! Memory allocation
 !
 call WF_alloc(WF,k_extrema_only_)
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
   wfc_d => WF%c_d
   WFN = WF%N
 endif
 !
 ! LOG the states loaded
 !
 if ((.not.QUIET_wf).and.fat_log) then
   if (buffer_is_ok.and.WF%space/="B") then
     call msg('rns','[WF'//trim(wf_title)//'] Loading Wave-Functions from buffer')
   else
     call msg('rns','[WF'//trim(wf_title)//'] Performing Wave-Functions I/O from '//trim(core_io_path)//'/SAVE')
   endif
 endif
 !
 if (WF%space=='R'.and.SIZE_msg) call msg('ns','[FFT'//trim(wf_title)//'] Mesh size',fft_dim)
 !
 ! actual IO
 !
 if (.not.buffer_is_ok) then
   !
   YAMBO_ALLOC(wf_b_indx,(wf_nb_io))
   wf_b_indx=0
   !
   call io_control(ACTION=OP_RD,COM=NONE,MODE=VERIFY,SEC=(/1/),ID=ID)
   io_err=io_WF(ID)
   !
 endif
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
     if (k_extrema_only_ .and. all( (/ikibz>WF%k(1), ikibz<WF%k(2)/) )) cycle
     !
     ! USPP init
     !
     if (pp_is_uspp) then
       !
       call c2a(b,k_pt(ikibz,:),xk,mode="ki2c")
       xk=xk/real(tpiba,SP)
       !
       npwk=wf_nc_k(ikibz)
       call init_us_2(npwk,wf_igk(1:npwk,ikibz),xk,vkb)
       !
     endif
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
         if ( allocated(states_to_load) ) then
           if ( .not.states_to_load(ib,ikibz,i_sp_pol) ) cycle
         endif
         !
         wf_nb_to_load=wf_nb_to_load+1
         wf_b_indx(wf_nb_to_load)=ib
         !
       enddo
       !
       if (wf_nb_to_load==0) cycle
       !
       if (.not.buffer_is_ok) then
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
         if (pp_is_uspp) then
            call PP_compute_becp( becp(ikibz,i_sp_pol), wf_nc_k(ikibz), wf_nb_to_load, wf_disk, wf_b_indx)
         endif
         !
       else
         !
         !if (pp_is_uspp) call PP_compute_becp(becp(ikibz,i_sp_pol),WF_buffer)
         if (pp_is_uspp) call error(" [PP] buffering not yet implemented")
         !
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
         if(buffer_is_ok) call error(" [CUDA WF_LOAD] buffering not yet implemented")
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
              !$cuf kernel do(1)
              do ic=1,fft_size
                wfc_d(ic,i_spinor,WFN)=cmplx(wf_DP_d(ic),kind=SP)*fft_norm
              enddo
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
                if (buffer_is_ok) then
                  wf_DP(igfft)=WF_buffer%c(ic,i_spinor,WF_buffer%index(ib,ikibz,i_sp_pol))
                else
                  wf_DP(igfft)=cmplx(wf_disk(ic,i_spinor,ib_to_load),kind=DP)
                endif
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
              WF%c(:,i_spinor,WF%N)=cmplx(wf_DP(:),kind=SP)*fft_norm
              !
            enddo  ! i_spinor
#ifdef _CUDA
         endif
#endif
         !
       enddo    ! i_b
       !
       if (.not.buffer_is_ok) then
         YAMBO_FREE(wf_disk)
         YAMBO_FREE(wf_disk_d)
       endif
       !
     enddo      ! i_b_grp
   enddo        ! i_k
 enddo          ! i_sp_pol

 if (have_cuda) then
   !
   call msg('rns', '[WF] Copying WF data from GPU device')
   if(WF%space=='R') then
     WF%c = WF%c_d
   else
     WF%c_d = WF%c
   endif
   !
 endif
 !
 call timing('WF_load_FFT',"stop")
 !
 ! CLEAN
 !
 if (.not.buffer_is_ok) then
   !
   call io_control(ACTION=RD_CL,COM=NONE,MODE=VERIFY,SEC=(/0/),ID=ID)
   io_err=io_WF(ID)
   !
   YAMBO_FREE(wf_b_indx)
   !
 endif
 !
 YAMBO_FREE(wf_DP)
 YAMBO_FREE(wf_DP_d)

 !
 ! The states_to_load is a temporary array that points to the states to be loaded
 ! here. It is defined in PARALLEL_WF_distribute and, therefore, only by calling 
 ! WF_load after a PARALLEL_WF_distribute OR by using the buffered components 
 ! the states_to_load are used.
 !
 ! With the next lines I assume that whatever next call of WF_load will 
 ! occur AFTER a PARALLEL_WF_distribute. But this is not true when I call
 ! WF_load with the same distribution of states but different spaces ("G" or "R").
 ! 
 ! This is the reason of the KEEP_STATES_TO_LOAD logical
 !
 if(allocated(states_to_load)) then
   !
   YAMBO_FREE(WF%state)
   !
   YAMBO_ALLOC(WF%state,(WF%b(1):WF%b(2),WF%k(1):WF%k(2),WF%sp_pol(1):WF%sp_pol(2)))
   !
   WF%state(WF%b(1):WF%b(2),WF%k(1):WF%k(2),WF%sp_pol(1):WF%sp_pol(2))= &
&           states_to_load(WF%b(1):WF%b(2),WF%k(1):WF%k(2),WF%sp_pol(1):WF%sp_pol(2))
   !
   ! When using the buffer the states_to_load must not be cleaned
   !
   clean_up_states=.not.allocated(WF_buffer%c)
   if (present(keep_states_to_load)) clean_up_states=.not.keep_states_to_load
   !
   if (clean_up_states) then
     YAMBO_FREE(states_to_load)
   endif
   !
 endif
 !
#if defined _SC
 !
 ! Rotate in the SC basis
 !========================
 !
 if (found_SC_DB.and..not.force_WFo_) then
   !
   call load_SC_components('WF',COM_=NONE,MODE_=DUMP,kind=WF_kind,xc_functional=WF_xc_functional,&
&                          perturbation=WF_perturbation)
   !
   call msg('rns','[WF'//trim(wf_title)//' loader] '//trim(global_XC_string(WF_kind,WF_xc_functional,WF_perturbation))//&
&                             ' wavefunctions loaded')
   !
   if (WF%space=='R') wf_grid_size=fft_size
   if (WF%space=='G') wf_grid_size=wf_ng
   do i_sp_pol=WF%sp_pol(1),WF%sp_pol(2)
     do ikibz=WF%k(1),WF%k(2)
       call WF_rotate(ikibz,i_sp_pol,WF,wf_grid_size)
     enddo
   enddo
   !
   call load_SC_components('WF_clean')
   !
 endif
 !
#endif
 !
 if (.not.wf_norm_test) then
   if (QUIET_wf) call IO_and_Messaging_switch("+report +log ")
   return
 endif
 !
 ! Check normalization @ 1st k-point, spin up only.
 !
 mndp=10._SP
 mxdp=-1._SP
 do ib1=1,min(int(nel)+5,WF%b(2))
   do ib2=1,min(int(nel)+5,WF%b(2))
     !
     if ( WF%index(ib1,1,1)==0 .or. WF%index(ib2,1,1)==0 ) cycle
     !
     N=size(WF%c(:,:,WF%index(ib1,1,1)))
     c=Vstar_dot_V( N, WF%c(:,:,WF%index(ib1,1,1) ) , WF%c(:,:,WF%index(ib2,1,1)) )
     !
     if (abs(c)>mxdp) mxdp=abs(c)
     if (abs(c)<mndp) mndp=abs(c)
     !
   enddo
 enddo
 !
 wf_norm_test=.false.
 call msg('rn','[WF'//trim(wf_title)//' loader] Normalization (few states)  min/max  ',(/mndp,mxdp/))
 !
 if (QUIET_wf) call IO_and_Messaging_switch("+report +log ")
 !
 contains
   !
   subroutine check_bounds(WF_,flag_)
     type(WAVEs)  :: WF_
     logical      :: flag_
     flag_=all((/bands_to_load(1)>=WF_%b(1),bands_to_load(2)<=WF_%b(2),&
&                kpts_to_load(1)>=WF_%k(1),kpts_to_load(2)<=WF_%k(2),&
&                sp_pol_to_load_(1)>=WF_%sp_pol(1),sp_pol_to_load_(2)<=WF_%sp_pol(2)/))
   end subroutine
   !
   subroutine check_states(WF_,flag_)
     type(WAVEs)  :: WF_
     logical      :: flag_
     !
     i_sp_pol=WF_%sp_pol(1)
     do while(i_sp_pol<=WF_%sp_pol(2).and.flag_)
       ikibz=WF_%k(1)
       do while(ikibz<=WF_%k(2).and.flag_) 
         ib   =WF_%b(1)
         do while(ib<=WF_%b(2).and.flag_) 
           flag_=states_to_load(ib,ikibz,i_sp_pol).eqv.WF_%state(ib,ikibz,i_sp_pol)
           ib=ib+1
         enddo
         ikibz=ikibz+1
       enddo
       i_sp_pol=i_sp_pol+1
     enddo
   end subroutine
   !
end subroutine
