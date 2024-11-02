!
! License-Identifier: GPL
!
! Copyright (C) 2006 The Yambo Team
!
! Authors (see AUTHORS file for details): AM
!
subroutine RT_ELPH_initialize(E,k,q)
 ! 
 use parser_m,       ONLY:parser
 use real_time,      ONLY:RT_bands
 use electrons,      ONLY:levels
 use drivers,        ONLY:l_elph_scatt,l_phel_scatt,l_rt_coh_phon
 use R_lattice,      ONLY:bz_samp
 use IO_int,         ONLY:io_control
 use IO_m,           ONLY:OP_RD_CL,DUMP,REP
 use ELPH_intfcs,    ONLY:ELPH_alloc
 use ELPH,           ONLY:PH_freqs_sq,PH_W_debye,GKKP_ep_sq,GKKP_pe_sq,l_GKKP_expanded_DB_exists,elph_nQ,&
 &                        GKKP
 !
#include<y_memory.h>
 !
 type(levels)             :: E
 type(bz_samp),intent(in) :: k,q
 !
 ! Work Space
 !
 logical             :: l_use_Debye_E
 integer             :: ID_gkkp,io_err
 integer, external   :: io_ELPH
 !
 call section('=','El-Ph matrix elements')
 !                 ---------------------
 ! Logicals
 !==========
 if (.not.(l_elph_scatt.or.l_phel_scatt.or.l_rt_coh_phon)) return
 !
 call parser('UseDebyeE'   , l_use_Debye_E)
 !
 ! El-Ph databases check
 !=======================
 call ELPH_databases_check(RT_bands(2))
 if (.not.l_GKKP_expanded_DB_exists) then
   call warning('ELPH (expanded) database not found and/or incompatible OR QP widths not provided. E-P scattering switched off')
   l_elph_scatt   =.FALSE.
   return
 endif
 !
 ! Q/K grid remapping
 !====================
 call ELPH_databases_grids_map(k,q)
 !
 ! El-Ph databases I/O
 !=====================
 if (l_elph_scatt) call ELPH_databases_load(GKKP_ep_sq,E,k,q,RT_bands,(/1,elph_nQ/),"Q-BZ K-IBZ K-sym","dressed-dressed")
 !
 ! Pj-El databases I/O
 !=====================
 if (l_phel_scatt) call ELPH_databases_load(GKKP_pe_sq,E,k,q,RT_bands,(/1,elph_nQ/),"Q-IBZ K-BZ Q-sym","dressed-dressed")
 !
 ! Q=0 databases for coherent phonon dynamics
 !=====================
 if (l_rt_coh_phon) then
  !
  call ELPH_alloc("dV",GKKP=GKKP)
  call io_control(ACTION=OP_RD_CL,COM=REP,MODE=DUMP,SEC=(/1,2/),ID=ID_gkkp)
  io_err=io_ELPH(ID_gkkp,'gkkp_expanded')                                
  !
 endif
 !
 ! Test of the effect of using a single optical mode
 !                               ===================
 !
 if (l_use_Debye_E) PH_freqs_sq=PH_W_debye**2
 !
 ! Q weights 
 !===========
 !
 call section('=','El-Ph BZ integrals')
 !                 -------------------
 call ELPH_acoustic_phonon_properties(q)
 !
#if defined _ELPH_SCATT
 call RT_ELPH_Q_weights(q)
#endif
 !
 ! THETA matrix indexes
 !======================
#if defined _ELPH_ITERATIVE
 call RT_ELPH_iterative_init( )
#endif
 !
end subroutine RT_ELPH_initialize
