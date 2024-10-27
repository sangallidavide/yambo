!
! License-Identifier: GPL
!
! Copyright (C) 2013 The Yambo Team
!
! Authors (see AUTHORS file for details): AM DS
!
subroutine PARALLEL_global_indexes(E,Xk,q,ENVIRONMENT,X,Dip,RESET,Dip_limits_pre_defined,X_all_bands)
 !
 ! LOGICALS
 use parallel_m,    ONLY:l_par_X_T,l_par_RT,l_par_SE,l_par_X_G,l_par_X_G_finite_q,&
&                        l_par_DIP,l_par_BZINDX,l_par_NL,l_par_SC,l_par_OSCLL
 ! COMMUNICATORS
 use parallel_m,    ONLY:PAR_COM_CON_INDEX_X,PAR_COM_VAL_INDEX_X
 ! IND
 use parallel_m,    ONLY:PAR_IND_CON_BANDS_X,PAR_IND_VAL_BANDS_X
 ! ID's
 use parallel_m,    ONLY:PAR_IND_CON_BANDS_X_ID,PAR_IND_VAL_BANDS_X_ID
 ! DIMENSIONS 
 use parallel_m,    ONLY:PAR_Dip_ib,PAR_Dip_ib_lim,PAR_X_ib,PAR_X_iq,PAR_n_c_bands,PAR_n_v_bands
 ! WORLDS
 use parallel_m,    ONLY:PAR_COM_WORLD,PAR_COM_X_WORLD,PAR_COM_X_WORLD_RL_resolved,COMM_copy
 use parallel_int,  ONLY:PARALLEL_index,PARALLEL_live_message
 use DIPOLES,       ONLY:DIPOLE_t,DIPOLES_reset,DIPOLES_duplicate
 use X_m,           ONLY:X_t,l_X_terminator,X_reset,X_duplicate
 use electrons,     ONLY:levels
 use R_lattice,     ONLY:bz_samp
 !
#include<y_memory.h>
 !
 type(levels)               ::E
 type(bz_samp)              ::Xk,q
 character(*)               ::ENVIRONMENT
 type(X_t),     optional    ::X
 type(DIPOLE_t),optional    ::Dip
 logical,       optional    ::RESET,Dip_limits_pre_defined,X_all_bands
 !
 ! Work Space
 !
 logical              :: Do_dip_limits,X_all_bands_
 integer              :: X_type,bands_loc(2)
 character(8)         :: MORE
 ! 
 type(X_t)      :: X_loc
 type(DIPOLE_t) :: Dip_loc
 !
 !===========================
 ! Response function pre-sets
 !===========================
 X_all_bands_=.false.
 if (present(X_all_bands)) X_all_bands_=X_all_bands
 !
 call X_reset(X_loc,1)
 if (present(X)) then
   call X_duplicate(X,X_loc)
   X_loc%whoami=X%whoami
 endif
 !
 X_type=X_loc%whoami
 if (ENVIRONMENT=="Response_T_space") X_type=5
 !
 !========
 ! Resets
 !========
 if (present(RESET)) then
   if (RESET) then
     call PARALLEL_global_reset(ENVIRONMENT,X_type)
     return
   endif
 else
   call PARALLEL_global_reset(ENVIRONMENT,X_type)
 endif
 !
 !===========================
 ! Dipoles function pre-sets
 !===========================
 call DIPOLES_reset(Dip_loc)
 if (present(Dip)) call DIPOLES_duplicate(Dip,Dip_loc)
 !
 MORE=""
 !
 if (present(Dip)) then
   Do_dip_limits=.TRUE.
   if (present( Dip_limits_pre_defined )) then
     Do_dip_limits=.not.Dip_limits_pre_defined
   endif
   if (trim(Dip%approach)=="Shifted grids") MORE="_SHIFTED"
 endif
 !
 !========
 ! WORLDS
 !========
 call COMM_copy(PAR_COM_WORLD,PAR_COM_X_WORLD)
 call COMM_copy(PAR_COM_WORLD,PAR_COM_X_WORLD_RL_resolved)
 !
 !==========
 ! Logicals
 !==========
 !
 l_par_BZINDX       =ENVIRONMENT=="BZ_Indexes"
 l_par_DIP          =ENVIRONMENT=="DIPOLES"
 l_par_X_T          =ENVIRONMENT=="Response_T_space"
 l_par_X_G_finite_q =ENVIRONMENT=="Response_G_space_and_IO"
 l_par_X_G          =ENVIRONMENT=="Response_G_space"
 l_par_SE           =ENVIRONMENT=="Self_Energy"
 l_par_RT           =ENVIRONMENT=="Real_Time"
 l_par_NL           =ENVIRONMENT=="Non_Linear"
 l_par_SC           =ENVIRONMENT=="Self_Energy"
 l_par_OSCLL        =ENVIRONMENT=="Oscillators"
 !
 !==================================
 ! USER provided PARALLEL structure
 !==================================
 !
 !... via the ENVIRONMENT variable
 !
 call PARALLEL_get_ENVIRONMENT_structure(ENVIRONMENT)
 !
 !...via the input file
 !
 call PARALLEL_get_user_structure(ENVIRONMENT,.TRUE.)
 !
 !============
 ! Dimensions
 !============
 !
 call PARALLEL_global_dimensions(E,Xk,q,X_loc,Dip_loc,Do_dip_limits,ENVIRONMENT)
 !
 !==========
 ! DEFAULTS
 !==========
 !
 CALL PARALLEL_global_defaults(ENVIRONMENT//trim(MORE))
 !
 !==============================
 ! ENVIRONMENT DEPENDENT SCHEMES
 !==============================
 !
 if (index(ENVIRONMENT, "Response_G_space")>0) call PARALLEL_global_Response_G(E,Xk,q,ENVIRONMENT,X_type)
 !
 if (      ENVIRONMENT=="BZ_Indexes")          call PARALLEL_global_BZINDX(E,Xk,q,ENVIRONMENT)
 !
 if (      ENVIRONMENT=="DIPOLES")             call PARALLEL_global_DIPOLES(E,Xk,q,ENVIRONMENT)
 !
 if (      ENVIRONMENT=="Self_Energy")         call PARALLEL_global_Self_Energy(E,Xk,q)
 !
 if (      ENVIRONMENT=="Response_T_space")    call PARALLEL_global_Response_T_base(Xk)
 !
#if defined _RT
 if (      ENVIRONMENT=="Real_Time")           call PARALLEL_global_Real_Time(E,Xk,q,X_type)
#endif
 !
#if defined _NL
 if (      ENVIRONMENT=="Non_Linear")          call PARALLEL_global_Non_Linear(E,Xk,q)
 if (      ENVIRONMENT=="Oscillators")         call PARALLEL_global_Oscillators(E,Xk,q)
#endif
 !
#if defined _SCALAPACK
 if (      ENVIRONMENT=="ScaLapacK")           call PARALLEL_global_ScaLapack( )
#endif
 !
 !==============================================================================
 if ( ENVIRONMENT=="Response_G_space".or.ENVIRONMENT=="Response_G_space_and_IO") then
   !============================================================================
   !
   ! Response functions conduction bands
   !
   bands_loc=PAR_n_c_bands
   if (l_X_terminator.or.X_all_bands_) bands_loc=PAR_X_ib
   !
   call PARALLEL_index(PAR_IND_CON_BANDS_X(X_type),(/bands_loc(2)/),low_range=(/bands_loc(1)/),&
&                        COMM=PAR_COM_CON_INDEX_X(X_type),CONSECUTIVE=.TRUE.,NO_EMPTIES=.TRUE.)
   !
   PAR_IND_CON_BANDS_X_ID(X_type)=PAR_COM_CON_INDEX_X(X_type)%CPU_id
   !
   call PARALLEL_live_message("CON bands",ENVIRONMENT=ENVIRONMENT,&
&                             LOADED=PAR_IND_CON_BANDS_X(X_type)%n_of_elements(PAR_COM_CON_INDEX_X(X_type)%CPU_id+1),&
&                             TOTAL=bands_loc(2)-bands_loc(1)+1,&
&                             NCPU=PAR_COM_CON_INDEX_X(X_type)%n_CPU)
   !
   ! Response functions valence bands
   !
   bands_loc=PAR_n_v_bands
   if (l_X_terminator.or.X_all_bands_) bands_loc=PAR_X_ib
   !
   call PARALLEL_index(PAR_IND_VAL_BANDS_X(X_type),(/bands_loc(2)/),low_range=(/bands_loc(1)/),&
&                      COMM=PAR_COM_VAL_INDEX_X(X_type),CONSECUTIVE=.TRUE.,NO_EMPTIES=.TRUE.)
   !
   PAR_IND_VAL_BANDS_X_ID(X_type)=PAR_COM_VAL_INDEX_X(X_type)%CPU_id
   !
   call PARALLEL_live_message("VAL bands",ENVIRONMENT=ENVIRONMENT,&
&                             LOADED=PAR_IND_VAL_BANDS_X(X_type)%n_of_elements(PAR_COM_VAL_INDEX_X(X_type)%CPU_id+1),&
&                             TOTAL=bands_loc(2)-bands_loc(1)+1,&
&                             NCPU=PAR_COM_VAL_INDEX_X(X_type)%n_CPU)
   !
 endif
 !
end subroutine
