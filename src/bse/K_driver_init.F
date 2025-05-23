!
! License-Identifier: GPL
!
! Copyright (C) 2021 The Yambo Team
!
! Authors (see AUTHORS file for details): AM
!
subroutine K_driver_init(what,iq,Ken,Xk)
 !
 use pars,         ONLY:SP,schlen
 use R_lattice,    ONLY:bz_samp
 use electrons,    ONLY:levels,n_sp_pol,n_spin
 use parser_m,     ONLY:parser
 use parallel_m,   ONLY:PARALLEL_default_mode
 use stderr,       ONLY:STRING_match,STRING_split
 use BS_solvers,   ONLY:BSS_mode,BSS_slepc_matrix_format,BSS_uses_DbGd,BSS_uses_FKE,BSS_ydiago_solver,&
&                       BSS_slepc_pseudo_herm,BSS_slepc_double_grp,run_Slepc,BSS_mode,BSS_first_eig,&
&                       l_abs_prop_chi_bse,l_eels_can_be_computed,l_eels_from_inversion,&
&                       BSS_kernel_IO_on_the_fly,BSS_kernel_stored_size
 use BS,           ONLY:BSE_L_kind,BSE_mode,BS_K_is_ALDA,BS_dip_size,l_BSE_minimize_memory,BS_perturbative_SOC,&
&                       BS_perturbative_SOC,l_BS_magnons,l_BS_photolum,BSK_IO_mode,&
&                       BS_cpl_K_exchange,BS_n_g_exch,BS_res_K_exchange,BS_K_coupling,BS_res_ares_n_mat,&
&                       BS_n_eh_spaces,l_BS_ares_from_res,BS_bands,BS_K_anti_resonant
 use TDDFT,        ONLY:FXC_is_retarded,FXC_n_g_corr,l_Fxc_from_Vxc
 use drivers,      ONLY:l_col_cut,l_rpa_IP,l_alda_fxc,l_bs_fxc
 use D_lattice,    ONLY:i_time_rev,i_space_inv,l_3D
#if defined _RT
 use drivers,      ONLY:l_rt_carriers_in_use
 use RT_control,   ONLY:NEQ_Kernel,EQ_Transitions,EQ_NoOcc,NEQ_Residuals,RT_BSE_Occ_Mode
#endif
#if defined _CUDA
 use TDDFT,        ONLY:FXC_mode
#endif
 !
 implicit none
 !
 character(*)  ::what
 integer       ::iq
 type(levels)  ::Ken
 type(bz_samp) ::Xk
 !
 ! Work Space
 !
 integer :: i1,BSK_IO_n_strings
 logical :: l_flag,l_slepc_special,l_Fxc_Libxc
 character(schlen) ::BSK_IO_strings(10)
 !
 ! Defaults
 !==========
 BS_n_eh_spaces=1
 BS_res_ares_n_mat=1
 !
 ! default mode
 if(STRING_match(BSE_L_kind,"default")) BSE_L_kind="default-bar"
 !
 if (l_col_cut.and.what=="init") then
   ! We are here in 0D, 1D or 2D
   ! The G=0 term included in L-full is representative
   ! of the portion of the line (2D) / plane (1D) / volume (0D)  around G=0
   ! which contributes to the L-T splitting.
   ! Such term goes to zero as the volume goes to infinity.
   if (iq==1) then
     ! At q=0 there is no L-T spliting. So Lbar and Lfull give the same result
     ! It is useless to compute the non analytical Lfull and better to
     ! directly remove the G=0 term
     call warning('Coulomb cutoff at q=0. Default/Suggested Lkind= bar')
     if (STRING_match(BSE_L_kind,"default")) BSE_L_kind="default-bar"
   endif
   if (iq/=1) then
     ! We are here in 1D or 2D.
     ! At q>0 there is L-T spliting in 1D or 2D
     ! Lbar and Lfull give the same result only in the limit V->infinity
     ! The above mentioned G=0 term is part of this effect.
     ! Here we assume it would be better to include it, to be double checked
     call warning('Coulomb cutoff at q>0. Default/Suggested Lkind= full')
     if (STRING_match(BSE_L_kind,"default")) BSE_L_kind="default-full"
     ! N.B.: One could compute the calculation without L-T slitting,
     ! by removing such full line (2D) / plane (1D). Not coded at the moment.
     ! I would define it as L="transverse"
   endif
 endif
 !
 ! Is absorption (alpha or epsilon) proportional to chi
 !===================================
 !
 ! * In 3D Alpha and epsilon are proporitonal to Lbar
 !
 ! * In 2D, 1D, and 0D Alpha is proportional to L(low_dim), where
 !   L(low_dim)=Lfull=Lbar in the V-> infinity volume
 !
 !   - If the cut-off is not used Lbar is a better approximation to L(low_dim)
 !     e.g. same situation of the 3D case
 !
 !   - If the cut-off is used Lfull is probably a better approximation to L(low_dim),
 !     e.g. opposite situation of the 3D case
 !
 ! Thus, overall, only without the coulomb cutoff Alpha il proporional to Lbar,
 ! e.g. it is proportional to epsilon. This is important without cut-off only
 !
 if (       l_col_cut .and. iq> 1) l_abs_prop_chi_bse=STRING_match(BSE_L_kind,"full")
 if ( (.not.l_col_cut) .or. iq==1) l_abs_prop_chi_bse=STRING_match(BSE_L_kind,"bar")
 !
 ! However in general the difference with the coulomb cutoff is expected to be tiny.
 ! Advanced user can impose to change behaviour.
 ! With Slepc this must be imposed by the code
 call parser('ImposeAbsPropChiBSE',l_flag)
 if (.not.(l_3D.or.l_abs_prop_chi_bse).and.l_flag) then
   l_abs_prop_chi_bse=.true.
   call warning("Imposting absorption ot be proportional to computed L")
 endif
 if (run_Slepc .and. .not.l_abs_prop_chi_bse) then
   if (l_col_cut .and. iq>1) then
     call warning("Slepc with coulomb cutoff at finite-q, and Lbar. abs not printed")
     call warning("You can print it by adding in input the flag ImposeAbsPropChiBse")
   else if (.not.l_3D) then
     call warning("Slepc without coulomb cutoff or q=0, and Lfull in low dim. abs not printed")
     call warning("You can print it by adding in input the flag ImposeAbsPropChiBse")
   endif
 endif
 !
 ! Either I'm forced to get the epsilon from extension + inversion
 ! or I want to get the epsilon from extension + inversion
 ! to account of G=0 term in  coupling mode
 ! Not needed if already in coupling mode
 ! Not possible with slepc
 l_eels_from_inversion= (l_abs_prop_chi_bse .or. .not.trim(BSE_mode)=="coupling")
 !
 ! If I cannot get the eels from the inversion and I'm in Lbar case I'm in troubles
 l_eels_can_be_computed = l_eels_from_inversion .or. .not.l_abs_prop_chi_bse
 !
 if (.not.l_eels_can_be_computed) call warning("EELs evaluation will be switched-off")
 !
 ! RES/ANTI-RES (q-dependent)
 !==============
 !
 ! Do I need the anti-res part?
 if (BS_K_coupling) BS_K_anti_resonant=.TRUE.
 if (STRING_match(BSE_mode,'causal').or.STRING_match(BSE_mode,'retarded')) BS_K_anti_resonant=.TRUE.
#if defined _RT 
 if (l_rt_carriers_in_use) BS_K_anti_resonant=.TRUE.
#endif
 !
 ! Is ARES derivable from RES?
 !-----------------------------
 ! Developer options. Without the code works in default mode
 ! With the user can switch to experimental coding
 call parser('ImposeAsym',l_flag)
 if (l_flag)        l_BS_ares_from_res=.FALSE.
 !
 if (iq/=1.and.i_time_rev==0.and.i_space_inv==0 )  l_BS_ares_from_res=.FALSE.
 if (l_BS_photolum)                                l_BS_ares_from_res=.FALSE.
 if (l_BS_magnons.and.n_sp_pol==2)                 l_BS_ares_from_res=.FALSE.
 ! The next line is to switch off sections of the code due to the calculation
 ! of the ARES part without checking the BSE_mode (see io_BS_PAR_init for example)
 if (trim(BSE_mode)=="resonant")                   l_BS_ares_from_res=.TRUE.
 !
 ! I need ARES and I cannot get it from RES
 ! 
 if (BS_K_anti_resonant.and..not.l_BS_ares_from_res) then
   BS_n_eh_spaces=2
   if (.not.BS_K_coupling) BS_res_ares_n_mat=2
 endif
 !
 if (what=="loop-init") return
 !
 ! BSE_L_kind check
 !==============
 if (.not.STRING_match(BSE_L_kind,"bar").and..not.STRING_match(BSE_L_kind,"full")) &
 & call error('set Lkind = DEFAULT or BAR or FULL')
 !
 BS_K_is_ALDA=l_alda_fxc
 if (l_alda_fxc) then
   call parser('FxcLibxc',l_Fxc_Libxc)
   l_Fxc_from_Vxc=(n_spin==1).and..not.l_Fxc_Libxc
 endif
#if defined _CUDA
 if ( BS_K_is_alda.and. index(FXC_mode,"G-")>0 ) then
   call warning(" Tddft with G-integrals is not GPU ported. Fallback to R-integrals")
   FXC_mode="R-def"
 endif
#endif
 BS_dip_size=1
 !
 ! Memory saving
 !===============
 l_BSE_minimize_memory=((Xk%nbz==1).or.index(PARALLEL_default_mode,"memory")>0)
 !
 if (l_BSE_minimize_memory) call warning(" BSE calculation with minimal memory utilization")
 !
 ! SOC
 !=====
 if (Ken%nb_SOC/=0) BS_perturbative_SOC=.TRUE.
 !
 ! Copatibility with old inputs
 !==============================
 if (STRING_match(BSE_mode,'causal')) BSE_mode='retarded'
 !
 ! BS-Fxc Setup
 !==============
 !
 ! Carriers and bs-fxc non compatible!
 !-------------------------------------
#if defined _RT 
 if(l_rt_carriers_in_use .and. STRING_match(BSS_mode,'t')) call error('RT occupations and BSE via f_xc not compatible')
#endif
 !
 ! If the em1s DB has been calculated in the same run
 ! of a BS-Fxc kernel construction that l_bs_fxc is set .FALSE.
 !      =============
 ! in TDDFT_do_X_W_typs. So here we need to check again if
 ! l_bs_fxc has been provided
 !
 if (.not.l_bs_fxc) call parser('bs_fxc',l_bs_fxc)
 !
 if (l_bs_fxc) then
   BS_n_g_exch=FXC_n_g_corr
   if (STRING_match(BSS_mode,'t')) then
     ! Remove the exchange if TDDFT Fxc is obtained from K
     BS_res_K_exchange=.FALSE.
     BS_cpl_K_exchange=.FALSE.
     BSE_mode='retarded'
     call parser('FxcRetarded',FXC_is_retarded)
   endif
 endif
 !
 ! Remove exchange if computing the irreducible response
 if (STRING_match(BSE_L_kind,"tilde") ) then
   BS_res_K_exchange=.FALSE.
   BS_cpl_K_exchange=.FALSE.
 endif
 !
 ! CORE Logicals
 !===============
 if (l_rpa_IP.and.STRING_match(BSE_mode,'coupling')) BSE_mode='retarded'
 !
 ! TR-ABS logicals
 !-----------------
#if defined _RT
 NEQ_Kernel   =STRING_match(RT_BSE_Occ_Mode,"k").and.allocated(Ken%fo)
 NEQ_Residuals=STRING_match(RT_BSE_Occ_Mode,"r").and.allocated(Ken%fo)
 call parser('ForceEqTrans',EQ_Transitions)
 call parser('ForceEqNoOcc',EQ_NoOcc)
#endif
 !
 ! SLEPC
 !=======
 !
#if defined _SLEPC && !defined _NL
 !
 if(index(BSS_slepc_matrix_format,"shell")==0.and.index(BSS_slepc_matrix_format,"explicit")==0) &
 & call error("Wrong value of slepc matrix format: "//trim(BSS_slepc_matrix_format))
 !
 ! Special case, slepc with distributed matrix need to explicitly compute the anti-resonant block
 ! This was needed in the past. Now the older implementation can be imposed
 ! chosing shell-debug matric format
 l_slepc_special = index(BSS_mode,'s')/=0 .and. BS_K_coupling .and. (BSS_slepc_matrix_format=="shell-debug")
 !
 BSS_slepc_double_grp=l_BS_ares_from_res.and.l_slepc_special
 if (l_slepc_special) l_BS_ares_from_res=.FALSE.
 !
 BSS_slepc_pseudo_herm=BS_K_coupling .and. l_BS_ares_from_res .and. &
 &                     (.not.BSS_slepc_double_grp) .and. index(BSS_slepc_matrix_format,"debug")==0
 !
#endif
 !
 ! Double Grid support
 !=====================
 BSS_uses_DbGd = ((Ken%FG%nb>=BS_bands(2)) .and. (index(BSS_mode,'i')/=0)) &
 &          .or. ((Ken%FG%nb>=BS_bands(2)) .and. (index(BSS_mode,'h')/=0))
 !
 call parser('BSS_FKE',BSS_uses_FKE)
 !
 BSS_kernel_IO_on_the_fly=index(BSK_IO_mode,"read_on_the_fly")>0 .and. (index(BSS_mode,'h')/=0 .or. index(BSS_mode,'s')/=0)
 if (BSS_kernel_IO_on_the_fly) then
   call STRING_split(BSK_IO_mode,BSK_IO_strings,n_non_empty_strings=BSK_IO_n_strings)
   BSS_kernel_stored_size=0._SP
   do i1=1,BSK_IO_n_strings
     if (.not.STRING_match(BSK_IO_strings(i1),"GB")) cycle
     read (BSK_IO_strings(i1-1),*) BSS_kernel_stored_size
   enddo
 endif
 !
#if defined _ELPA
 if ( index(BSS_mode,'d')/=0 .and. BSS_ydiago_solver=='e' .and. BSS_first_eig>1) then
   call warning('ELPA does not support BSS_First_eig>1')
   BSS_first_eig=1
 endif
#endif

 !
end subroutine
