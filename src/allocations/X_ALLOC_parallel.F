!
! License-Identifier: GPL
!
! Copyright (C) 2017 The Yambo Team
!
! Authors (see AUTHORS file for details): AM
!
subroutine X_ALLOC_parallel(X_par,NG,NW,mode)
 !
 use pars,          ONLY:cZERO
 use parallel_m,    ONLY:COMM_copy,PAR_COM_RL_INDEX,PAR_COM_RL_A2A,PAR_COM_HOST
 use matrix,        ONLY:MATRIX_reset,PAR_matrix
 use X_m,           ONLY:X_par_lower_triangle,X_rows,X_cols,&
&                        X_FILL_UP_matrix_only,X_lower_triangle_matrix_in_use
 use parallel_int,  ONLY:PARALLEL_live_message,PP_redux_wait
 use cuda_m,        ONLY:have_cuda
 !
#include<memory.h>
 !
 type(PAR_matrix) :: X_par
 integer          :: NG,NW
 integer(IPL)     :: HOST_SIZE(3),TMP_SIZE(3)
 character(*)     :: mode
 !
 ! Work Space
 !
 logical :: l_XUP,l_PAR,l_CPU_ONLY
 !
 ! Reset...
 !
 call MATRIX_reset(X_par)
 !
 ! Clean-up
 !
 if (mode=="CLEAN") return
 !
 ! Type of allocation....
 !
 l_XUP      = index(mode,"XUP")>0
 l_PAR      = index(mode,"X")>0.and.index(mode,"XUP")==0
 l_CPU_ONLY = index(mode,"CPU_ONLY")>0
 !
 ! Initialization
 !
 if (l_XUP) then
   call MATRIX_init( "XUP", X_par , NG, NW )
 else if (l_PAR) then
   call MATRIX_init( "X"  , X_par , NG, NW )
 endif
 !
 X_rows = X_par%rows
 X_cols = X_par%cols
 !
 if (l_XUP.and..not.X_FILL_UP_matrix_only.and.X_par%INTER_comm%n_CPU>1) then
   !
   call MATRIX_reset(X_par_lower_triangle)
   X_lower_triangle_matrix_in_use=.TRUE.
   X_par_lower_triangle%kind ="XDN"
   X_par_lower_triangle%rows =X_par%cols
   X_par_lower_triangle%cols =X_par%rows
   X_par_lower_triangle%nrows=X_par%ncols
   X_par_lower_triangle%ncols=X_par%nrows
   !
   !TMP_SIZE=(/X_par%cols(2)-X_par%cols(1)+1,X_par%rows(2)-X_par%rows(1)+1,NW*2/)
   !AMBO_PAR_ALLOC_CHECK3(X_par_lower_triangle%blc,TMP_SIZE)
   YAMBO_ALLOC(X_par_lower_triangle%blc,(X_par%cols(1):X_par%cols(2),X_par%rows(1):X_par%rows(2),NW))
   X_par_lower_triangle%blc=cZERO
   if (have_cuda.and..not.l_CPU_ONLY) then
     YAMBO_ALLOC(X_par_lower_triangle%blc_d,(X_par%cols(1):X_par%cols(2),X_par%rows(1):X_par%rows(2),NW))
     X_par_lower_triangle%blc_d=cZERO
   endif
   !
   call COMM_copy(PAR_COM_RL_A2A,   X_par_lower_triangle%INTRA_comm)
   call COMM_copy(PAR_COM_RL_INDEX, X_par_lower_triangle%INTER_comm)
 endif
 !
 !TMP_SIZE=(/X_par%rows(2)-X_par%rows(1)+1,X_par%cols(2)-X_par%cols(1)+1,NW*2/)
 !AMBO_PAR_ALLOC_CHECK3(X_par%blc,TMP_SIZE)
 YAMBO_ALLOC(X_par%blc,(X_par%rows(1):X_par%rows(2),X_par%cols(1):X_par%cols(2),NW))
 X_par%blc=cZERO
 if (have_cuda.and..not.l_CPU_ONLY) then
   YAMBO_ALLOC(X_par%blc_d,(X_par%rows(1):X_par%rows(2),X_par%cols(1):X_par%cols(2),NW))
   X_par%blc_d=cZERO
 endif
 !
 call PARALLEL_live_message("RL vectors(X)",LOADED=X_par%ncols*X_par%nrows,TOTAL=NG**2,NCPU=X_par%INTER_comm%n_CPU)
 !
end subroutine
