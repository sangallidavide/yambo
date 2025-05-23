!
! License-Identifier: GPL
!
! Copyright (C) 2017 The Yambo Team
!
! Authors (see AUTHORS file for details): AM
!
subroutine PARALLEL_global_Self_Energy(E,Xk,q)
 !
 use electrons,     ONLY:levels
 use R_lattice,     ONLY:bz_samp,nqbz,nXkibz,nXkbz
 use drivers,       ONLY:l_eval_collisions,l_sc_electric,l_sc_run,l_real_time
 use parallel_int,  ONLY:PARALLEL_index,PARALLEL_assign_chains_and_COMMs,&
&                        PARALLEL_live_message,PARALLEL_MATRIX_distribute
 use IO_m,          ONLY:io_COLLs
 use IO_int,        ONLY:IO_and_Messaging_switch
 use QP_m,          ONLY:QP_n_states,QP_nk
 use openmp,        ONLY:n_threads_SE,OPENMP_set_threads
 use parallel_m,    ONLY:HEAD_QP_cpu,HEAD_b_cpu,master_cpu,n_WF_bands_to_load,COMM_copy,&
&                        PAR_Q_bz_range,PAR_Q_ibz_range,PAR_n_bands,PAR_build_index,ncpu,PAR_INDEX_copy
 ! COMMUNICATORS
 use parallel_m,    ONLY:PAR_COM_QP_INDEX,PAR_COM_Q_INDEX,PAR_COM_RL_INDEX,PAR_COM_G_b_INDEX,&
&                        PAR_COM_QP_A2A,PAR_COM_Q_A2A,PAR_COM_RL_A2A,PAR_COM_SE_WORLD,&
&                        PAR_COM_SE_WORLD_RL_resolved
 ! IND
 use parallel_m,    ONLY:PAR_IND_QP,PAR_IND_Q_ibz,PAR_IND_G_b,PAR_IND_Q_bz,PAR_IND_Xk_ibz,PAR_IND_G_k,&
&                        PAR_IND_DIPk_bz,PAR_IND_DIPk_ibz,PAR_IND_OVLPk_ibz
 ! INDEX
 use parallel_m,    ONLY:PAR_QP_index,PAR_Q_ibz_index,PAR_Q_bz_index,PAR_G_bands_index,PAR_IND_Bp_mat
 ! DIMENSIONS
 use parallel_m,    ONLY:PAR_nG_bands,PAR_nQ_ibz,PAR_nQ_bz,PAR_nQP,PAR_n_Bp_mat_elements
 ! ID's
 use parallel_m,    ONLY:PAR_IND_Xk_ibz_ID,PAR_IND_G_b_ID,PAR_IND_QP_ID,PAR_IND_Q_ibz_ID,PAR_IND_Q_bz_ID
#if defined _SC
 use parallel_m,    ONLY:PAR_COM_Xk_ibz_INDEX,PAR_IND_Xk_bz
#endif
#if defined _SC || defined _RT || defined _QED
 use drivers,       ONLY:l_elphoton_corr,l_life,l_real_time,l_sc_run
 use collision_ext, ONLY:COLL_bands,COH_collisions,HXC_collisions,P_collisions
 use parallel_m,    ONLY:COMM_copy
 ! COMMUNICATORS
 use parallel_m,    ONLY:PAR_COM_QP_INDEX,PAR_COM_G_b_INDEX
 ! IND
 use parallel_m,    ONLY:PAR_IND_QP,PAR_IND_Q_ibz,PAR_IND_Q_bz
 ! INDEX
 use parallel_m,    ONLY:PAR_DIPk_ibz_index,PAR_DIPk_bz_index
 ! DIMENSIONS
 use parallel_m,    ONLY:PAR_DIPk_nibz
#endif
#if defined _ELPH
 use drivers,       ONLY:l_elph_corr
 use ELPH,          ONLY:elph_use_q_grid
#endif
 !
#include<y_memory.h>
 !
 type(levels)         :: E
 type(bz_samp)        :: Xk,q
 !
 ! Work space
 !
 character(10)        :: WHAT
#if defined _SC || defined _RT
 integer              :: i_k,nb_mat
#endif
 !
 ! AF: here we disable the RL parallelism (g) for SC and RT, where
 !     the parallel structure needs to be adapted to the SE g parallelism
 !
 if (l_sc_run.or.l_real_time) then
   !
   call PARALLEL_structure(3,(/"q ","qp","b "/))
   call PARALLEL_assign_chains_and_COMMs(3,COMM_index_1=PAR_COM_Q_INDEX,&
&                                          COMM_index_2=PAR_COM_QP_INDEX,&
&                                          COMM_index_3=PAR_COM_G_b_INDEX,&
&                                          COMM_A2A_1=PAR_COM_Q_A2A,&
&                                          COMM_A2A_2=PAR_COM_QP_A2A)
   !
 else
   !
   call PARALLEL_structure(4,(/"q ","g ","qp","b "/))
   call PARALLEL_assign_chains_and_COMMs(4,COMM_index_1=PAR_COM_Q_INDEX,&
&                                          COMM_index_2=PAR_COM_RL_INDEX,&
&                                          COMM_index_3=PAR_COM_QP_INDEX,&
&                                          COMM_index_4=PAR_COM_G_b_INDEX,&
&                                          COMM_A2A_1=PAR_COM_Q_A2A,&
&                                          COMM_A2A_2=PAR_COM_RL_A2A,&
&                                          COMM_A2A_3=PAR_COM_QP_A2A)
 endif
 !
 ! COMMs setup
 !
 call COMM_copy(PAR_COM_Q_A2A,PAR_COM_SE_WORLD)
 call COMM_copy(PAR_COM_Q_A2A,PAR_COM_SE_WORLD_RL_resolved)
 !
 !
 if (PAR_COM_RL_INDEX%n_CPU>1) call COMM_copy(PAR_COM_RL_A2A,PAR_COM_SE_WORLD_RL_resolved)
 !
 !
 ! The routine PARALLEL_assign_chains_and_COMMs cannot define COMMUNICATORS for internal
 ! A2A when there is no internal distribution
 !
 if (PAR_COM_QP_INDEX%n_CPU==1) then
   call COMM_copy(PAR_COM_Q_A2A,PAR_COM_QP_A2A)
 endif
 !
 ! QP states
 !
 call PARALLEL_index(PAR_IND_QP,(/QP_n_states/),COMM=PAR_COM_QP_INDEX,NO_EMPTIES=.TRUE.)
 PAR_IND_QP_ID=PAR_COM_QP_INDEX%CPU_id
 PAR_nQP=PAR_IND_QP%n_of_elements(PAR_IND_QP_ID+1)
 YAMBO_ALLOC(PAR_QP_index,(QP_n_states))
 PAR_QP_index=0
 call PAR_build_index(PAR_IND_QP,QP_n_states,PAR_QP_index,PAR_nQP)
 call PARALLEL_live_message("QPs",ENVIRONMENT="Self_Energy",&
&                           LOADED=PAR_IND_QP%n_of_elements(PAR_IND_QP_ID+1),TOTAL=QP_n_states,&
&                           NCPU=PAR_COM_QP_INDEX%n_CPU)
 !
 ! Q-points
 !
 WHAT="ibz"
 if (l_eval_collisions)                       WHAT="bz"
#if defined _ELPH
 if (l_elph_corr.and.elph_use_q_grid)         WHAT="bz"
 if (l_elph_corr.and..not.elph_use_q_grid)    WHAT="RIM"
#endif
 !
 if (trim(WHAT)=="ibz") then
   !
   call PARALLEL_index(PAR_IND_Q_ibz,(/PAR_Q_ibz_range(2)/),COMM=PAR_COM_Q_INDEX,CONSECUTIVE=.TRUE.,NO_EMPTIES=.TRUE.)
   PAR_IND_Q_ibz_ID=PAR_COM_Q_INDEX%CPU_id
   PAR_nQ_ibz=PAR_IND_Q_ibz%n_of_elements(PAR_IND_Q_ibz_ID+1)
   call PARALLEL_live_message("Q("//trim(WHAT)//")",ENVIRONMENT="Self_Energy",&
&                             LOADED=PAR_IND_Q_ibz%n_of_elements(PAR_IND_Q_ibz_ID+1),TOTAL=PAR_Q_ibz_range(2),&
&                             NCPU=PAR_COM_Q_INDEX%n_CPU)
   !
   YAMBO_ALLOC(PAR_Q_ibz_index,(PAR_Q_ibz_range(2)))
   call PAR_build_index(PAR_IND_Q_ibz,PAR_Q_ibz_range(2),PAR_Q_ibz_index,PAR_nQ_ibz)
   !
 else if (trim(WHAT)=="bz".or.trim(WHAT)=="RIM") then
   !
   call PARALLEL_index(PAR_IND_Q_bz,(/PAR_Q_bz_range(2)/),COMM=PAR_COM_Q_INDEX,CONSECUTIVE=.TRUE.,NO_EMPTIES=.TRUE.)
   PAR_IND_Q_bz_ID=PAR_COM_Q_INDEX%CPU_id
   PAR_nQ_bz=PAR_IND_Q_bz%n_of_elements(PAR_IND_Q_bz_ID+1)
   call PARALLEL_live_message("Q("//trim(WHAT)//")",ENVIRONMENT="Self_Energy",&
&                              LOADED=PAR_IND_Q_bz%n_of_elements(PAR_IND_Q_bz_ID+1),TOTAL=PAR_Q_bz_range(2),&
&                              NCPU=PAR_COM_Q_INDEX%n_CPU)
   !
   YAMBO_ALLOC(PAR_Q_bz_index,(PAR_Q_bz_range(2)))
   call PAR_build_index(PAR_IND_Q_bz,PAR_Q_bz_range(2),PAR_Q_bz_index,PAR_nQ_bz)
   !
 endif
 !
 if (trim(WHAT)=="ibz") then
   YAMBO_ALLOC(PAR_Q_bz_index,(nqbz))
   call PARALLEL_distribute_BZk_using_IBZk(PAR_COM_Q_INDEX,q,PAR_IND_Q_ibz,PAR_IND_Q_ibz_ID,&
&                                                   PAR_IND_Q_bz, PAR_IND_Q_ibz_ID,&
&                                                   PAR_Q_bz_index,PAR_nQ_bz)
   call PAR_build_index(PAR_IND_Q_bz,nqbz,PAR_Q_bz_index,PAR_nQ_bz)
 endif
 !
#if defined _SC || defined _RT
 !
 if (l_sc_run.or.l_real_time.or.l_eval_collisions) then
   !
   ! 00000             <- q
   ! x0 x0    00 00    <- qp
   ! 0 0 0 0  0 0 0 0  <- b
   !
   ! QP_cpu corresponds to x marked CPU's. This flag is used when isolated QP loops are performed.
   !
   HEAD_QP_cpu=PAR_COM_Q_INDEX%CPU_id==0.and.PAR_COM_QP_A2A%CPU_id==0
   if (PAR_COM_Q_INDEX%n_CPU==1.and.PAR_COM_QP_A2A%n_CPU==1) then
     HEAD_QP_cpu=PAR_COM_G_b_index%CPU_id==0
   endif
   !
 else
   !   
#endif
   !
   ! 0000     0000     <- q
   ! 0000     0000     <- g
   ! x0 x0    x0 x0    <- qp
   ! 0 0 0 0  0 0 0 0  <- b
   !
   ! QP_cpu corresponds to x marked CPU's. This flag is used when no b loops are done
   !
   HEAD_QP_cpu=PAR_COM_QP_A2A%CPU_id==0
   !
#if defined _SC || defined _RT
 endif
#endif
 !
 ! K-points
 !
 if (trim(WHAT)=="ibz".or.trim(WHAT)=="RIM") then
   call PARALLEL_add_Q_to_K_list('k_qp_q_'//trim(WHAT),PAR_IND_QP,PAR_IND_QP_ID,PAR_IND_Xk_ibz,PAR_IND_Xk_ibz_ID,&
&                                PAR_IND_Q_ibz,PAR_COM_QP_INDEX,(/0,0/),Xk,q)
 else if (trim(WHAT)=="bz") then
   call PARALLEL_add_Q_to_K_list('k_qp_q_'//trim(WHAT),PAR_IND_QP,PAR_IND_QP_ID,PAR_IND_Xk_ibz,PAR_IND_Xk_ibz_ID,&
&                                PAR_IND_Q_bz,PAR_COM_QP_INDEX,(/0,0/),Xk,q)
 endif
 !
 ! G bands
 !=========
 !
 ! WF bands to load
 !
 n_WF_bands_to_load=PAR_n_bands(2)
 !
 call PARALLEL_index(PAR_IND_G_b,(/PAR_n_bands(2)/),low_range=(/PAR_n_bands(1)/),&
&                    COMM=PAR_COM_G_b_INDEX,CONSECUTIVE=.TRUE.,NO_EMPTIES=.TRUE.)
 PAR_IND_G_b_ID=PAR_COM_G_b_INDEX%CPU_id
 PAR_nG_bands=PAR_IND_G_b%n_of_elements(PAR_IND_G_b_ID+1)
 YAMBO_ALLOC(PAR_G_bands_index,(PAR_n_bands(2)))
 call PAR_build_index(PAR_IND_G_b,PAR_n_bands(2),PAR_G_bands_index,PAR_nG_bands)
 !
 HEAD_b_cpu=PAR_COM_G_b_index%CPU_id==0
 !
 call PARALLEL_live_message("G bands",ENVIRONMENT="Self_Energy",&
&                           LOADED=PAR_IND_G_b%n_of_elements(PAR_IND_G_b_ID+1),&
&                           TOTAL=PAR_n_bands(2)-PAR_n_bands(1)+1,NCPU=PAR_COM_G_b_INDEX%n_CPU)
 !
#if defined _SC || defined _RT || defined _QED
 !
 ! Lamb is a special case. Collisions can be used also in simple Self-Energy mode.
 !
 if (l_elphoton_corr.or.l_eval_collisions) then
   call PARALLEL_collisions( Xk,   P_collisions )
 endif
 !
 ! To eval/use the COLLISIONS I need the PAR_IND_Bp_mat, built here.
 !
 if (l_sc_run.or.l_real_time.or.l_eval_collisions) then
   !
   !.........................................................................
   !   "COLLISIONS"
   !.........................................................................
   !
   call PARALLEL_collisions( Xk, COH_collisions )
   call PARALLEL_collisions( Xk, HXC_collisions )
   !
   nb_mat=(COLL_bands(2)-COLL_bands(1)+1)**2
   call PARALLEL_MATRIX_distribute(PAR_COM_G_b_INDEX,PAR_IND_Bp_mat,COLL_bands,PAR_n_elements=PAR_n_Bp_mat_elements)
   !
   call PARALLEL_live_message("Bands Matrix (prime)",ENVIRONMENT="Self_Energy",&
&                             LOADED=PAR_n_Bp_mat_elements,TOTAL=nb_mat,NCPU=PAR_COM_G_b_INDEX%n_CPU)
   !
   call PARALLEL_index(PAR_IND_G_k,(/QP_nk/),COMM=PAR_COM_Q_INDEX,CONSECUTIVE=.TRUE.,NO_EMPTIES=.TRUE.)
   !
   call PARALLEL_live_message("k-q",ENVIRONMENT="Self_Energy",LOADED=count( PAR_IND_G_k%element_1D ),&
&                             TOTAL=QP_nk)
   !
 endif
 !
 if (l_elphoton_corr.and.l_life) then
   !
   ! Oscillators needed to evaluate the Radiative Lifetimes
   !
   YAMBO_FREE(PAR_DIPk_ibz_index)
   YAMBO_FREE(PAR_DIPk_bz_index)
   YAMBO_ALLOC(PAR_DIPk_ibz_index,(nXkibz))
   do i_k=1,nXkibz
     PAR_DIPk_ibz_index(i_k)=i_k
   enddo
   YAMBO_ALLOC(PAR_DIPk_bz_index,(nXkbz))
   do i_k=1,nXkbz
     PAR_DIPk_bz_index(i_k)=i_k
   enddo
   PAR_DIPk_nibz=nXkibz
   PAR_DIPk_nibz=nXkibz
   ! 
 endif
 !
#endif
 !
#if defined _SC
 !
 if (l_sc_electric) then
   !
   ! Parallelization indexes for the Berry Polarization
   !
   call PARALLEL_index(PAR_IND_Xk_bz,(/nXkbz/),COMM=PAR_COM_Q_INDEX,CONSECUTIVE=.TRUE.)
   call COMM_copy(PAR_COM_Q_INDEX,PAR_COM_Xk_ibz_INDEX)
   !
   call PAR_INDEX_copy(PAR_IND_Xk_bz,PAR_IND_DIPk_bz)
   call PAR_INDEX_copy(PAR_IND_Xk_ibz,PAR_IND_DIPk_ibz)
   !
   allocate(PAR_DIPk_ibz_index(Xk%nibz))
   call PAR_build_index(PAR_IND_DIPk_ibz,Xk%nibz,PAR_DIPk_ibz_index,PAR_DIPk_nibz)
   !
   ! Overlap indeces for WF distribution
   !
   call PAR_INDEX_copy(PAR_IND_DIPk_ibz,PAR_IND_OVLPk_ibz)
   !
 endif
 !
#endif
 !
 ! I/O privileges
 !
 if (ncpu>1) then
   if (l_eval_collisions.and.io_COLLs) then
     call IO_and_Messaging_switch("+io_out",CONDITION=PAR_COM_QP_A2A%CPU_id==0.and.PAR_COM_Q_INDEX%CPU_id==0)
   else 
     call IO_and_Messaging_switch("+io_out",CONDITION=master_cpu)
   endif
 else
   call IO_and_Messaging_switch("+io_out",CONDITION=.TRUE.)
 endif
 !
 call IO_and_Messaging_switch("+output",CONDITION=master_cpu)
 !
 call OPENMP_set_threads(n_threads_in=n_threads_SE)
 !
end subroutine PARALLEL_global_Self_Energy
