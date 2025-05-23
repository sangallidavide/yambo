!
! License-Identifier: GPL
!
! Copyright (C) 2014 The Yambo Team
!
! Authors (see AUTHORS file for details): AM DS
!
#include<dev_defs.h>
!
subroutine DIPOLE_g_space(Xen,Xk,Dip)
 !
 ! This routine returns <iR> and <P^2>.
 !
 use pars,                ONLY:SP,cZERO,schlen
 use drivers,             ONLY:l_sc_run
 use interfaces,          ONLY:WF_load,WF_free
 use LIVE_t,              ONLY:live_timing
 use electrons,           ONLY:levels,n_spinor,n_sp_pol,spin_string
 use parallel_m,          ONLY:PAR_IND_DIPk_ibz,PAR_DIPk_ibz_index,PAR_IND_DIPk_ibz_ID,&
&                              PAR_IND_VAL_BANDS_DIP,PAR_IND_CON_BANDS_DIP
 use pseudo,              ONLY:PP_free,pp_table,pp_n_l_comp,pp_kbv_dim_yambo,pp_kbv_table,pp_factor,pp_kbv_dim_atom,pp_kbv_dim
 use D_lattice,           ONLY:n_atomic_species,n_atoms_species
 use R_lattice,           ONLY:bz_samp,g_vec_d,g_vec
 use DIPOLES,             ONLY:DIPOLE_t,DIP_iR,P_square,DIP_P,DIP_v,DIP_P_spinor,g_space_obs
 use ALLOC,               ONLY:DIPOLE_alloc_elemental 
 use X_m,                 ONLY:l_X_terminator
 use IO_int,              ONLY:io_control
 use IO_m,                ONLY:manage_action,RD_CL_IF_END,OP_RD,OP_RD_CL
 use wave_func,           ONLY:wf_ncx,WF
 use com,                 ONLY:msg
 use timing_m,            ONLY:timing
 use gpu_m,               ONLY:have_gpu,gpu_devsync
 use devxlib,             ONLY:devxlib_memcpy_h2d,devxlib_memset_d
 use stderr,              ONLY:STRING_match
 !
#include<y_memory.h>
 !
 type(bz_samp), intent(in)    :: Xk
 type(levels),  intent(in)    :: Xen
 type(DIPOLE_t),intent(inout) :: Dip
 !
 ! Work Space
 !
 logical                      :: l_use_comm,l_check
 logical, allocatable         :: element_1D_spin(:)
 integer                      :: ik,ic_wf,iv_wf,i_wf,ic,iv,ib,ifrag,pp_dim_yambo,pp_dim_atom,&
&                                i_sp_pol,N_elements_todo,ik_mem,ic_min,is,ia,pp_range(2),npp,npp_tab
 complex(SP), allocatable          :: kbv(:,:,:,:)
 complex(SP), allocatable DEV_ATTR :: kbv_d(:,:,:,:)
 complex(SP), allocatable DEV_ATTR :: XX(:,:,:,:,:)
 complex(SP)                  :: rho_spinor(3,n_spinor,n_spinor),rho(3),P2
 real(SP)                     :: Ev_m_Ev_KS,Dipole_Energy_treshold
 logical                      :: LT_msg=.FALSE.
 !
 !I/O and external functions
 !
 integer                      :: ID,IO_ACT_NOW,io_err(3),io_KB_abinit_err,io_KB_pwscf_err
 integer, external            :: io_KB_abinit,io_KB_pwscf
 !
 call timing('DIPOLE_g_space',OPR='start')
 !
 ! CPU-depedent work to be done 
 !==============================
 N_elements_todo=0
 do iv=Dip%ib(1),Dip%ib_lim(1) 
  if(.not.PAR_IND_VAL_BANDS_DIP%element_1D(iv)) cycle
  ic_min=max(iv,Dip%ib_lim(2))
  if (l_X_terminator) ic_min=iv
  do ic=ic_min,Dip%ib(2)
    if(.not.PAR_IND_CON_BANDS_DIP%element_1D(ic)) cycle
    N_elements_todo=N_elements_todo+1
  enddo
 enddo
 N_elements_todo=N_elements_todo*n_sp_pol*PAR_IND_DIPk_ibz%n_of_elements(PAR_IND_DIPk_ibz_ID+1)
 !
 ! Synchronize some memory
#ifdef _GPU
 YAMBO_FREE_GPU(DEV_VAR(g_vec))
 YAMBO_ALLOC_GPU_SOURCE(DEV_VAR(g_vec),g_vec)
 if(STRING_match(g_space_obs,"P2")) call warning("P2 not coded with GPUs, CPU only algorithm for it")
#endif
 !
 l_use_comm=STRING_match(g_space_obs,"V").or.STRING_match(g_space_obs,"R")
 !
 if(l_use_comm) then
   !
   call DIPOLE_kb_init(Xen,Dip,io_err,ID)
   !
   io_KB_abinit_err= io_err(1)
   io_KB_pwscf_err = io_err(2)
   !
   IO_ACT_NOW=OP_RD
   if (N_elements_todo==0) IO_ACT_NOW=OP_RD_CL
   !
   ! Define pp_kbv_dim_yambo and pp_kbv_table
   !
   if (io_KB_pwscf_err==0) then
     call io_control(ACTION=IO_ACT_NOW,SEC=(/2/),ID=ID)
     io_KB_pwscf_err=io_KB_pwscf(ID)
     call DIPOLE_kb_pwscf_def_dim()
   endif
   !
   if (io_KB_abinit_err==0) then
     call io_control(ACTION=IO_ACT_NOW,SEC=(/2/),ID=ID)
     io_KB_abinit_err=io_KB_abinit(ID)
     call DIPOLE_kb_abinit_def_dim()
   endif
   !
   l_use_comm=l_use_comm.and.(io_KB_abinit_err==0.or.io_KB_pwscf_err==0).and.pp_kbv_dim_yambo>0
   !
 endif
 !
 if(l_use_comm) then
   !
   ! allocation of DEV variabless with zero dims is a problem (PGI, CUDA For)
   pp_dim_yambo=max(pp_kbv_dim_yambo,1)
   pp_dim_atom =max(pp_kbv_dim_atom, 1)
   YAMBO_ALLOC(kbv,(wf_ncx,n_spinor,4,pp_dim_atom))
   !
   if (have_gpu) then
     YAMBO_ALLOC_GPU(XX,(n_spinor,n_spinor,4,pp_dim_yambo,Dip%ib(1):Dip%ib(2)))
     YAMBO_ALLOC_GPU(DEV_VAR(kbv),(wf_ncx,n_spinor,4,pp_dim_atom))
     call devxlib_memset_d(DEV_VAR(kbv),cZERO)
   else
     YAMBO_ALLOC(XX,(n_spinor,n_spinor,4,pp_dim_yambo,Dip%ib(1):Dip%ib(2)))
   endif
   !
   ! Eo    are always the eigenvalues consistent with the WFs
   ! E     do also include the QP corrections when Eo are allocated
   !
   if (allocated(Xen%Eo).and.STRING_match(g_space_obs,"R")) then
     call msg('rns','[X] Using energies without QP corrections for the dipoles')
   endif
   !
   allocate(element_1D_spin(Xk%nibz*n_sp_pol))
   if (n_sp_pol==1) element_1D_spin=PAR_IND_DIPk_ibz%element_1D
   if (n_sp_pol==2) element_1D_spin=(/PAR_IND_DIPk_ibz%element_1D(:),PAR_IND_DIPk_ibz%element_1D(:)/)
   !
   ! sum_kb workspace
   !
   call DIPOLE_alloc_elemental("DIP_work",[pp_dim_yambo])
   !
 endif
 !
 ! Main loop over k in IBZ
 !
 if (N_elements_todo>0) call live_timing(trim(g_space_obs)//" [g-space]",N_elements_todo)
 !
 do i_sp_pol=1,n_sp_pol
   !
   if (N_elements_todo==0) cycle
   !
   do ik=1,Xk%nibz
     ! 
     if (.not.PAR_IND_DIPk_ibz%element_1D(ik)) cycle
     !
     ifrag=ik+(i_sp_pol-1)*Xk%nibz
     ! 
     ! [Vnl,r]
     !
     if(l_use_comm) then
       !
       IO_ACT_NOW=manage_action(RD_CL_IF_END,ifrag,1,Xk%nibz*n_sp_pol,element_1D=element_1D_spin)
       call io_control(ACTION=IO_ACT_NOW,SEC=(/ifrag+1/),ID=ID)
       !
       if (io_KB_pwscf_err==0 .and.ifrag>1)  io_KB_pwscf_err=io_KB_pwscf(ID)
       if (io_KB_abinit_err==0.and.ifrag>1) io_KB_abinit_err=io_KB_abinit(ID)
       !
     endif
     !
     ik_mem=PAR_DIPk_ibz_index(ik)
     !
     ! WF load
     !
     call WF_load(WF,0,1,Dip%ib,(/ik,ik/),sp_pol_to_load=(/i_sp_pol,i_sp_pol/),space='C',&
&                 title='-Oscillators/G space/Transverse '//spin_string(i_sp_pol),force_WFo=l_sc_run,keep_states_to_load=.TRUE.)
     !
     if (l_use_comm) then
       !
       npp=0
       npp_tab=0
       pp_range(2)=0
       call devxlib_memset_d(XX,cZERO)
       !
       do is = 1,n_atomic_species
         do ia = 1,n_atoms_species(is)
           !
           pp_range(1)=pp_range(2)+1
           !
           if (io_KB_abinit_err==0) call DIPOLE_kb_abinit_comp(npp,npp_tab,ia,is,ik,i_sp_pol,Xk,kbv)
           if (io_KB_pwscf_err ==0) call DIPOLE_kb_pwscf_comp(npp,npp_tab,ia,is,ik,Xk,kbv)
           !
           if (npp==0) cycle
           !
           pp_range(2)=pp_range(1)+npp-1
           !
           if (have_gpu) call devxlib_memcpy_h2d(DEV_VAR(kbv),kbv)
           !
           do ib=Dip%ib(1),Dip%ib(2)
             !
             if( ib<=Dip%ib_lim(1) ) then
               l_check= .not.(PAR_IND_VAL_BANDS_DIP%element_1D(ib).or.PAR_IND_CON_BANDS_DIP%element_1D(ib))
             else
               l_check= .not.PAR_IND_CON_BANDS_DIP%element_1D(ib)
             endif
             if ( l_check ) cycle
             !
             i_wf=WF%index(ib,ik,i_sp_pol)
             !
             call DIPOLE_kb_project(pp_range,Dip%ib,ib,i_wf,wf_ncx,size(WF%c,3),npp,&
             &                      pp_dim_atom,pp_dim_yambo,DEV_VAR(WF%c),DEV_VAR(kbv),XX)
             !
           enddo
           !
         enddo
       enddo
       !
     endif
     !
     do iv=Dip%ib(1),Dip%ib_lim(1) 
       !
       if(.not.PAR_IND_VAL_BANDS_DIP%element_1D(iv)) cycle
       !
       ic_min=max(iv,Dip%ib_lim(2))
       if (l_X_terminator) ic_min=iv
       !
       do ic=ic_min,Dip%ib(2)
         !
         if(.not.PAR_IND_CON_BANDS_DIP%element_1D(ic)) cycle
         !
         iv_wf=WF%index(iv,ik,i_sp_pol)
         ic_wf=WF%index(ic,ik,i_sp_pol)
         !
         if (iv_wf==0) call error(" Error in parallel wf distribution (Dipole iv)")
         if (ic_wf==0) call error(" Error in parallel wf distribution (Dipole ic)")
         !
         ! Evaluate <iv_wf|   p   |ic_wf>
         !================================
         call DIPOLE_p_matrix_elements(Xk,iv_wf,ic_wf,ik,rho,rho_spinor,P2,&
&                                      STRING_match(g_space_obs,"P2"),STRING_match(g_space_obs,"M_spin"))
         !
         !=====
         ! <P>
         !=====
         ! P_vc=rho
         !
         DIP_P(:,ic,iv,ik_mem,i_sp_pol)=rho
         !
         !======
         ! <P^2>
         !======
         !
         if (STRING_match(g_space_obs,"P2")) P_square(ic,iv,ik_mem,i_sp_pol)= P2
         !
         !============
         ! <P_spinor>
         !============
         !
#if defined _RT
         if (STRING_match(g_space_obs,"P_spinor")) DIP_P_spinor(:,:,:,ic,iv,ik_mem)=rho_spinor
#endif
         !
         if(.not.(STRING_match(g_space_obs,"V").or.STRING_match(g_space_obs,"R"))) then
           call live_timing(steps=1)
           cycle
         endif
         !
         ! Evaluate <iv_wf|[x,Vnl]|ic_wf>
         !================================
         if (l_use_comm) call DIPOLE_kb_sum(rho,pp_dim_yambo,Dip%ib,iv,ic,XX)
         !
         ! define a different thr for val-val transitions
         ! used for XTerm (here Ev_m_Ev_KS may become very large)
         !
         Dipole_Energy_treshold=Dip%Energy_treshold
         ! DEBUG <
         ! if(l_X_terminator.and.ic<=Dip%ib_lim(1)) Dipole_Energy_treshold=Dip%Energy_treshold_vv
         ! DEBUG >
         !
         !=====
         ! <v>
         !=====
         ! v_vc=rho
         !
         if (STRING_match(g_space_obs,"V")) DIP_v(:,ic,iv,ik_mem,i_sp_pol)=rho
         !
         !======
         ! <iR>
         !======
         !
         if (STRING_match(g_space_obs,"R")) then
           ! 
           !  [x,p_x]=i we get [x,H] = [x,p^2/2]+[x,Vnl]
           ! 
           ! we gave that 
           ! 
           ! [x,p^2/2]= i p_x 
           !
           ! So from the local part of rho (that is P_vc) we can extract the corresponding local part of <x>
           !
           ! DIP_iR(c,v) = i <v|r|c> = i <v|[r,H]|c>/(Ec-Ev) = 
           !             = i <v|i p|c>/(Ec-Ev) + i<v|[x,Vnl]|c>/(Ec-Ev) =
           !             = - <v|-i grad|c>/(Ec-Ev) +i<v|[x,Vnl]|c>/(Ec-Ev) =
           !             =   <v|-i grad|c>/(Ev-Ec) -i<v|[x,Vnl]|c>/(Ev-Ec) = 
           !             =   P_vc/(Ev-Ec) -i<v|[x,Vnl]|c>/(Ev-Ec)
           !
           if (.not.allocated(Xen%Eo)) Ev_m_Ev_KS=Xen%E (iv,ik,i_sp_pol)-Xen%E (ic,ik,i_sp_pol)
           if (     allocated(Xen%Eo)) Ev_m_Ev_KS=Xen%Eo(iv,ik,i_sp_pol)-Xen%Eo(ic,ik,i_sp_pol)
           !
           if (abs(Ev_m_Ev_KS)> Dipole_Energy_treshold) DIP_iR(:,ic,iv,ik_mem,i_sp_pol)=rho/Ev_m_Ev_KS
           if (abs(Ev_m_Ev_KS)<=Dipole_Energy_treshold) DIP_iR(:,ic,iv,ik_mem,i_sp_pol)=cZERO
           !
         endif
         !
         call live_timing(steps=1)
         !
       enddo    ! conduction band loop
     enddo      ! valence    band loop
     !
     call WF_free(WF)
     ! 
   enddo        ! k-points loop
 enddo          ! sp_pol   loop
 !
 call live_timing()
 !
 ! CLEAN
 !
 if (l_use_comm) then
   call DIPOLE_alloc_elemental("DIP_work")
   if (have_gpu) then
     YAMBO_FREE_GPU(DEV_VAR(kbv))
     YAMBO_FREE_GPU(XX)
   else
     YAMBO_FREE(XX)
   endif
   YAMBO_FREE(DEV_VAR(kbv))
   YAMBO_FREE(pp_n_l_comp)
   YAMBO_FREE(pp_table)
   YAMBO_FREE(pp_kbv_table)
   YAMBO_FREE(pp_factor)
   call PP_free()
 endif
 if (have_gpu) then
   YAMBO_FREE_GPU(DEV_VAR(g_vec))
 endif
 !
 call timing('DIPOLE_g_space',OPR='stop')
 !
end subroutine DIPOLE_g_space
