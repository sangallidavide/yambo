!
! License-Identifier: GPL
!
! Copyright (C) 2017 The Yambo Team
!
! Authors (see AUTHORS file for details): CA MG DS
!
subroutine PARALLEL_global_Non_Linear(E,Xk,q)
 !
 use electrons,     ONLY:levels
 use R_lattice,     ONLY:bz_samp
 use collision_ext, ONLY:HXC_collisions
 use nl_optics,     ONLY:NL_bands
 use IO_int,        ONLY:IO_and_Messaging_switch
 use parallel_int,  ONLY:PARALLEL_index,PARALLEL_assign_chains_and_COMMs,PARALLEL_live_message
 use parallel_m,    ONLY:ncpu,COMM_copy,PAR_build_index,PAR_INDEX_copy
 ! COMMUNICATORS
 use parallel_m,    ONLY:PAR_COM_freqs_INDEX,PAR_COM_freqs_A2A,PAR_COM_freqs,&
&                        PAR_COM_Xk_ibz_INDEX,PAR_COM_WF_k_INDEX,PAR_COM_Xk_ibz_A2A
 ! IND
 use parallel_m,    ONLY:PAR_IND_freqs,PAR_IND_Xk_ibz, &
&                        PAR_IND_QP,PAR_IND_Q_ibz,PAR_IND_WF_k,PAR_IND_Xk_bz
 ! INDEX
 use parallel_m,    ONLY:PAR_freqs_index,PAR_Xk_ibz_index,PAR_Xk_bz_index,PAR_QP_index
 ! DIMENSIONS
 use parallel_m,    ONLY:PAR_n_freqs,PAR_nQP,PAR_Xk_nbz,PAR_Xk_nibz
 ! ID's
 use parallel_m,    ONLY:PAR_IND_freqs_ID,PAR_IND_Q_ibz_ID,PAR_IND_WF_k_ID,&
&                        PAR_IND_Xk_bz_ID,PAR_IND_Xk_ibz_ID
 use nl_optics,     ONLY:n_runs
 use openmp,        ONLY:n_threads_NL,OPENMP_set_threads
 use QP_m,          ONLY:QP_n_states,QP_table
 !
 implicit none
 !
 type(levels)         :: E
 type(bz_samp)        :: Xk,q
 !
 ! Work space
 !
 integer      ::  i_qp,ib1,ib2,i_k
 !
 CALL PARALLEL_structure(2,(/"w ","k "/))
 !
 call PARALLEL_assign_chains_and_COMMs(2,COMM_index_1=PAR_COM_freqs_INDEX, &
&                                        COMM_index_2=PAR_COM_Xk_ibz_INDEX,&
&                                        COMM_A2A_1=PAR_COM_freqs_A2A)
 !
 ! The routine PARALLEL_assign_chains_and_COMMs cannot define COMMUNICATORS for internal
 ! A2A when there is no internal distribution
 !
 if (PAR_COM_Xk_ibz_INDEX%n_CPU==1) then
   call COMM_copy(PAR_COM_freqs_A2A,PAR_COM_Xk_ibz_A2A)
 endif
 !
 ! Freq-points
 !
 call PARALLEL_index(PAR_IND_freqs,(/n_runs/),COMM=PAR_COM_freqs_INDEX)
 PAR_IND_freqs_ID=PAR_COM_freqs_INDEX%CPU_id
 PAR_n_freqs     =PAR_IND_freqs%n_of_elements(PAR_IND_Q_ibz_ID+1)
 call PARALLEL_live_message("Freqs",ENVIRONMENT="Non_Linear",         &
&                           LOADED=PAR_IND_freqs%n_of_elements(PAR_IND_freqs_ID+1), &
&                           TOTAL=n_runs,NCPU=PAR_COM_freqs_INDEX%n_CPU)
 !
 allocate(PAR_freqs_index(n_runs)) 
 call PAR_build_index(PAR_IND_freqs,n_runs,PAR_freqs_index,PAR_n_freqs)
 !
 ! K-points
 !
 call PARALLEL_index(PAR_IND_Xk_ibz,(/Xk%nibz/),COMM=PAR_COM_Xk_ibz_INDEX)
 PAR_IND_Xk_ibz_ID=PAR_COM_Xk_ibz_INDEX%CPU_id
 PAR_Xk_nibz=PAR_IND_Xk_ibz%n_of_elements(PAR_IND_Xk_ibz_ID+1)
 !
 ! Dipoles
 ! 
 allocate(PAR_Xk_ibz_index(Xk%nibz))
 call PAR_build_index(PAR_IND_Xk_ibz,Xk%nibz,PAR_Xk_ibz_index,PAR_Xk_nibz)
 !
 ! Define the Xk_bz indeces for BZ parallelization
 ! 
 allocate(PAR_Xk_bz_index(Xk%nbz))
 call PARALLEL_distribute_BZk_using_IBZk(PAR_COM_Xk_ibz_INDEX,Xk,PAR_IND_Xk_ibz,PAR_IND_Xk_ibz_ID,&
&                                                       PAR_IND_Xk_bz, PAR_IND_Xk_bz_ID,&
&                                                       PAR_Xk_bz_index,PAR_Xk_nbz)

 !
 ! Overlap indeces for WF distribution
 !
 call PARALLEL_live_message("H(ibz)",ENVIRONMENT="Non_Linear",&
          LOADED=PAR_IND_Xk_ibz%n_of_elements(PAR_IND_Xk_ibz_ID+1),TOTAL=Xk%nibz,&
          NCPU=PAR_COM_Xk_ibz_INDEX%n_CPU)
 !
 ! Copy to WF
 !
 call PAR_INDEX_copy(PAR_IND_Xk_ibz,PAR_IND_WF_k)
 call COMM_copy(PAR_COM_Xk_ibz_INDEX,PAR_COM_WF_k_INDEX)
 PAR_IND_WF_k_ID=PAR_COM_WF_k_INDEX%CPU_id
 !
 ! Build QP indexes for the Collisions
 !
 allocate(PAR_QP_index(QP_n_states),PAR_IND_QP%n_of_elements(ncpu),PAR_IND_QP%element_1D(QP_n_states)) 
 PAR_IND_QP%element_1D=.FALSE.
 PAR_IND_QP%n_of_elements=0
 PAR_nQP=0
 PAR_QP_index=0
 do i_qp=1,QP_n_states
   ib1 = QP_table(i_qp,1)
   ib2 = QP_table(i_qp,2)
   i_k = QP_table(i_qp,3)
   if (ib1<NL_bands(1)) cycle
   if (ib2<NL_bands(1)) cycle
   if (.not.PAR_IND_Xk_ibz%element_1D(i_k)) cycle
   PAR_nQP=PAR_nQP+1
   PAR_QP_index(i_qp)=PAR_nQP
   PAR_IND_QP%element_1D(i_qp)=.TRUE.
   PAR_IND_QP%n_of_elements=PAR_nQP
 enddo
 call PARALLEL_live_message("QPs",ENVIRONMENT="Non_Linear",LOADED=PAR_nQP,TOTAL=QP_n_states)

 !.........................................................................
 !   "COLLISIONS"
 !.........................................................................
 !
 call PARALLEL_collisions( Xk,    HXC_collisions )
 ! 
 !
 ! Io and messaging
 !
 call IO_and_Messaging_switch("+io_out",CONDITION=PAR_COM_Xk_ibz_INDEX%CPU_id==0)
 !
 call IO_and_Messaging_switch("+output",CONDITION=PAR_COM_Xk_ibz_INDEX%CPU_id==0)
 !
 call OPENMP_set_threads(n_threads_in=n_threads_NL)
 !
end subroutine PARALLEL_global_Non_Linear
