!
! License-Identifier: GPL
!
! Copyright (C) 2007 The Yambo Team
!
! Authors (see AUTHORS file for details): AM DS
!
integer function variables_BS(ID,iq,desc,CLOSE_the_menu,X)
 !
 use pars,          ONLY:SP,IP,LP,IP_YIO,LP_YIO,CR_YIO
 use units,         ONLY:HA2EV
 use drivers,       ONLY:l_td_hf
 use X_m,           ONLY:X_t
 use R_lattice,     ONLY:RIM_ng,RIM_n_rand_pts,cut_geometry
 use BS,            ONLY:BS_identifier,BS_K_dim,BS_bands,BSE_L_kind,&
&                        BS_K_coupling,BS_res_K_exchange,BS_res_K_corr,&
&                        BS_n_g_exch,BS_n_g_W,BS_eh_en,BS_eh_win,BS_W_is_diagonal,&
&                        BS_K_is_ALDA,BS_cpl_K_exchange,&
&                        BS_cpl_K_corr,BSqpts,BS_perturbative_SOC
 use IO_m,          ONLY:io_status,variable_exist,write_is_on
 use IO_int,        ONLY:def_variable_elemental,io_variable_elemental,ver_is_gt_or_eq
 use descriptors,   ONLY:IO_desc,IO_desc_reset
#if defined _PAR_IO
 use BS,            ONLY:BSK_IO_mode
#endif
 !
 implicit none
 !
 integer            ::ID,iq
 type(IO_desc)      ::desc
 logical            ::CLOSE_the_menu
 type(X_t), optional::X
 !
 ! Work Space
 !
 logical             ::res_corr_disk,res_exch_disk,cpl_disk,&
&                      cpl_corr_disk,cpl_exch_disk,td_hf_disk
 !
 variables_BS=0
 !
 call IO_desc_reset(desc)
 !
 call def_variable_elemental(ID,"BS_identifier",1,IP_YIO,0)
 call io_variable_elemental(ID,VAR="[K] Identifier",I0=BS_identifier,DESCRIPTOR=desc)
 !
 if(iq>1) then
   call def_variable_elemental(ID,"Q-point",3,SP,0)
   call io_variable_elemental(ID,VAR="[K] Q-point",&
&                             R1=BSqpts(:,iq),CHECK=.true.,OP=(/"==","==","=="/),DESCRIPTOR=desc)
 endif
 !
 call def_variable_elemental(ID,"Dimension",1,IP_YIO,0)
 call io_variable_elemental(ID,VAR="[K] Dimension",I0=BS_K_dim(1),CHECK=.true.,OP=(/"=="/),DESCRIPTOR=desc)
 !
 call def_variable_elemental(ID,"Bands",2,IP_YIO,0)
 call io_variable_elemental(ID,VAR="[K] Bands",I1=BS_bands,CHECK=.true.,OP=(/"==","=="/),DESCRIPTOR=desc)
 !
 call def_variable_elemental(ID,"RES_EXCH",1,LP_YIO,0)
 call io_variable_elemental(ID,DB_L0=res_exch_disk,VAR="[K] Exchange res",&
&                           L0=BS_res_K_exchange,CHECK=.true.,OP=(/"=="/),DESCRIPTOR=desc)
 !
 call def_variable_elemental(ID,"RES_CORR",1,LP_YIO,0)
 call io_variable_elemental(ID,DB_L0=res_corr_disk,&
&     VAR="[K] Correlation res",L0=BS_res_K_corr,CHECK=.true.,OP=(/"=="/),DESCRIPTOR=desc)
 !
 call def_variable_elemental(ID,"COUPLING",1,LP_YIO,0)
 call io_variable_elemental(ID,DB_L0=cpl_disk,&
&     VAR="[K] Kernel`s coupling        ",L0=BS_K_coupling,CHECK=.true.,OP=(/"=="/),DESCRIPTOR=desc)
 !
 call def_variable_elemental(ID,"CPL_EXCH",1,LP_YIO,0)
 call io_variable_elemental(ID,DB_L0=cpl_exch_disk,&
&     VAR="[K] Exchange cpl",L0=BS_cpl_K_exchange,CHECK=.true.,OP=(/"=="/),DESCRIPTOR=desc)
 !
 call def_variable_elemental(ID,"CPL_CORR",1,LP_YIO,0)
 call io_variable_elemental(ID,DB_L0=cpl_corr_disk,&
&     VAR="[K] Correlation cpl",L0=BS_cpl_K_corr,CHECK=.true.,OP=(/"=="/),DESCRIPTOR=desc)
 !
 call def_variable_elemental(ID,"W_BARE",1,LP_YIO,0)
 call io_variable_elemental(ID,DB_L0=td_hf_disk,&
&       VAR="[K] W interaction is bare",L0=l_td_hf,CHECK=.true.,OP=(/"=="/),DESCRIPTOR=desc)
 !
 call def_variable_elemental(ID,"ALDA_KERNEL",1,LP_YIO,0)
 call io_variable_elemental(ID,&
&     VAR="[K] ALDA kernel in R-space",L0=BS_K_is_ALDA,CHECK=.true.,OP=(/"=="/),DESCRIPTOR=desc)
 !
 if (res_exch_disk.or.cpl_exch_disk) then
   call def_variable_elemental(ID,"RL_EXCH",1,IP_YIO,0)
   call io_variable_elemental(ID,&
&      VAR="[K] RL vectors exchange",I0=BS_n_g_exch,CHECK=.true.,OP=(/"=="/),DESCRIPTOR=desc)
 endif
 if (res_corr_disk.or.cpl_corr_disk) then
   call def_variable_elemental(ID,"RL_CORR",1,IP_YIO,0)
   call io_variable_elemental(ID,&
&     VAR="[K] RL vectors correlation",I0=BS_n_g_W,CHECK=.true.,OP=(/"=="/),DESCRIPTOR=desc)
 endif
 !
 call def_variable_elemental(ID,"EH_E_RANGE",2,SP,0)
 call io_variable_elemental(ID,UNIT=HA2EV,&
&     VAR="[K] E/h energy range",TERMINATOR="[eV]",R1=BS_eh_en,CHECK=.true.,OP=(/"==","=="/),DESCRIPTOR=desc)
 !
 call def_variable_elemental(ID,"CPL_RANGE",2,SP,0)
 call io_variable_elemental(ID,&
&     VAR="[K] Coupling range",TERMINATOR="[o/o]",R1=BS_eh_win,CHECK=.true.,OP=(/"==","=="/),DESCRIPTOR=desc)
 !
 ! Interaction
 !
 if (res_corr_disk.or.cpl_corr_disk) then
   call def_variable_elemental(ID,"BS_W_DIAGONAL",1,LP_YIO,0)
   call io_variable_elemental(ID,&
&       VAR="[W] Interaction is diagonal",L0=BS_W_is_diagonal,CHECK=.true.,OP=(/"=="/),DESCRIPTOR=desc)
   !
   call def_variable_elemental(ID,"BS_W_NG",1,IP_YIO,0)
   call io_variable_elemental(ID,VAR="[W] Matrix size",I0=X%ng,WARN=.true.,OP=(/"=="/),DESCRIPTOR=desc)
   !
   ! Only if full screened interaction is used
   !
   if (.not.td_hf_disk) then
     !
     call def_variable_elemental(ID,"BS_W_BANDS",2,IP_YIO,0)
     call io_variable_elemental(ID,VAR="[W] Bands",I1=X%ib,CHECK=.true.,OP=(/"==","=="/),DESCRIPTOR=desc)
     !
     call def_variable_elemental(ID,"BS_W_EH_RANGE",2,SP,0)
     call io_variable_elemental(ID,UNIT=HA2EV,&
&         VAR="[W] e/h energy range",TERMINATOR="[eV]",R1=X%ehe,CHECK=.true.,OP=(/">=","<="/),DESCRIPTOR=desc)
     !
     call def_variable_elemental(ID,"BS_W_POLES",1,SP,0)
     call io_variable_elemental(ID,VAR="[W] Poles",TERMINATOR="[o/o]",R0=X%cg_percentual,CHECK=.true.,OP=(/"=="/),DESCRIPTOR=desc)
     !
     call def_variable_elemental(ID,"BS_W_DIP_NG",1,IP_YIO,0)
     call io_variable_elemental(ID,&
&         VAR="[W] RL vectors in the sum",I0=X%ngostnts,WARN=.true.,OP=(/"<="/),DESCRIPTOR=desc)
     !
     call def_variable_elemental(ID,"BS_W_DIP_rVnl",1,LP_YIO,0)
     call io_variable_elemental(ID,&
&         VAR="[W] [r,Vnl] included",L0=X%Vnl_included,WARN=.true.,OP=(/"=="/),DESCRIPTOR=desc)
     !
     call def_variable_elemental(ID,"BS_W_q_dir",3,SP,0)
     call io_variable_elemental(ID,VAR="[W] Field direction",R1=X%q0,CHECK=.true.,OP=(/"==","==","=="/),DESCRIPTOR=desc)
     call def_variable_elemental(ID,"",0,0,0)
     !
     call def_variable_elemental(ID,"W_Cutoff",1,CR_YIO,0)
     call io_variable_elemental(ID,VAR="[W] Coulomb Cutoff",CH0=cut_geometry,CHECK=.true.,OP=(/"=="/),DESCRIPTOR=desc)
     !
     ! W_TDDFT_KERNEL
     !
     call def_variable_elemental(ID,"W_TDDFT_KERNEL",1,CR_YIO,0)
     call io_variable_elemental(ID,VAR='[W] Hxc Kernel',CH0=X%KERNEL,CHECK=.true.,DESCRIPTOR=desc)
     !
   endif
   !
   ! W_RIM
   !
   ! Here the name W_RIM is correct. 
   ! However the ndb.BS_PAR_Q1 used for restart in the test-suite should
   ! be updated, otherwise the hBN/GW-OPTICS/ corresponding test fails
   !
   call def_variable_elemental(ID,"W_DbGd",2,IP_YIO,0)
     call io_variable_elemental(ID,I0=RIM_ng,&
&     VAR='[RIM] RL components col',CHECK=.true.,OP=(/"=="/),DESCRIPTOR=desc,DO_NOT_DUMP=.TRUE.)
     call io_variable_elemental(ID,I0=RIM_n_rand_pts,&
&     VAR='[RIM] Random points col',CHECK=.true.,OP=(/"=="/),DESCRIPTOR=desc,DO_NOT_DUMP=.TRUE.)
   !
   ! Perturbative SOC
   !
   if (write_is_on(ID).or.variable_exist(ID,"PERT_SOC")) then
     !
     ! SOC
     !
     call def_variable_elemental(ID,"PERT_SOC",1,LP_YIO,MENU=0)
     call io_variable_elemental(ID,L0=BS_perturbative_SOC,&
&       VAR='[SOC] Perturbative SOC',CHECK=.true.,OP=(/"=="/),DESCRIPTOR=desc,DO_NOT_DUMP=.TRUE.)
   endif
   !
 endif
 !
#if defined _PAR_IO
 !
 ! Check if BSE kernel is written in triangular form
 !
 if (write_is_on(ID).or.variable_exist(ID,"BSE_KERNEL_IO_MODE")) then
   call def_variable_elemental(ID,"BSE_KERNEL_IO_MODE",1,CR_YIO,0)
   call  io_variable_elemental(ID,VAR="[IO] Kernel IO mode",CH0=BSK_IO_mode)
 else
   BSK_IO_mode="2D_standard"
 endif
 !
#endif
 !
 if (write_is_on(ID).or.variable_exist(ID,"X_kind")) then
   call def_variable_elemental(ID,"X_kind",1,CR_YIO,0)
   call io_variable_elemental(ID,VAR="[X] Response Kind",CH0=BSE_L_kind,CHECK=.true.,OP=(/"=="/),DESCRIPTOR=desc)
 endif
 !
 if (     CLOSE_the_menu) call def_variable_elemental(ID,"",0,0,1)
 if (.not.CLOSE_the_menu) call def_variable_elemental(ID,"",0,0,0)
 !
 variables_BS=io_status(ID)
 !
end function
