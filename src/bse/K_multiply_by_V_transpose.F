!
! License-Identifier: GPL
!
! Copyright (C) 2021 The Yambo Team
!
! Authors (see AUTHORS file for details): AM MG DS IMA
!
subroutine K_multiply_by_V_transpose(Vi,Vo,iter_par,iq)
 !
 ! Given K distributed on the processors, and Vi distributed on the processors
 ! multiply Tranpose(K) by Vi.
 !
 ! See mod_BS for details on the K structure
 !
 ! IA (2022/01/19):
 ! Extra dimension of %fragment and %E hardcoded to 1
 ! To be modified if DbGd applied to Slepc
 ! 
 use pars,           ONLY:SP,cONE,cZERO,cI
 use parallel_int,   ONLY:PP_redux_wait
 use parallel_m,     ONLY:PAR_IND_T_Haydock,PAR_COM_T_Haydock
 use BS_solvers,     ONLY:Haydock_vector,BSS_Slepc_double_grp,BSS_kernel_last_stored_block
 use BS,             ONLY:BS_nT_grps,BS_blk,n_BS_blks,BS_T_grp,l_BS_ares_from_res,BS_ID
 use IO_int,         ONLY:io_control
 use IO_m,           ONLY:RD
 use wrapper,        ONLY:M_by_V
 !
 implicit none
 !
 integer,     intent(in)::iter_par,iq
 type(Haydock_vector)   ::Vi(BS_nT_grps),Vo(BS_nT_grps)
 !
 ! Work Space
 !
 integer         ::i_B,i_Tk,i_Tp,i_g
 complex(SP)     ::alpha
 !
 ! The factor cI is needed for the coupling when l_BS_ares_from_res=.TRUE.
 !
 if(.not.BSS_Slepc_double_grp) alpha=-real(iter_par,SP)*cI
 if(     BSS_Slepc_double_grp) alpha=-real(iter_par,SP)
 !
 forall(i_g=1:BS_nT_grps,PAR_IND_T_Haydock%element_1D(i_g)) Vo(i_g)%fragment = cZERO
 !
 do i_B=1,n_BS_blks
   !
   i_Tk  =BS_blk(i_B)%iT_k
   i_Tp  =BS_blk(i_B)%iT_p
   !
   if (i_B>BSS_kernel_last_stored_block) then
     call io_control(ACTION=RD,ID=BS_ID(2))
     call io_BS_PAR_block(iq,i_B,BS_ID(2),"matonly")
     call K_block_symmetrize(i_B)
   endif
   !
   select case(BS_blk(i_B)%mode)
   case("R")
     !
     ! Resonant part 
     !
     if (i_Tk==i_Tp) then
       !
       call M_by_V('t',BS_blk(i_B)%size(1),BS_blk(i_B)%size(2),cONE,&
&                      BS_blk(i_B)%mat,BS_blk(i_B)%size(1),Vi(i_Tp)%fragment(:,1),1,cONE,&
&                      Vo(i_Tk)%fragment(:,1),1)
       !
       ! Add the energy diffs on the diagonal
       Vo(i_Tk)%fragment(:,1) = Vo(i_Tk)%fragment(:,1) + BS_T_grp(i_Tk)%E(:,1)*Vi(i_Tk)%fragment(:,1)
       !
     else
       !
       call M_by_V('n',BS_blk(i_B)%size(1),BS_blk(i_B)%size(2),cONE,&
&                      conjg(BS_blk(i_B)%mat),BS_blk(i_B)%size(1),Vi(i_Tp)%fragment(:,1),1,cONE,&
&                      Vo(i_Tk)%fragment(:,1),1)
       !
       ! Resonant part is hermitian
       !
       call M_by_V('t',BS_blk(i_B)%size(1),BS_blk(i_B)%size(2),cONE,&
&                      BS_blk(i_B)%mat,BS_blk(i_B)%size(1),Vi(i_Tk)%fragment(:,1),1,cONE,&
&                      Vo(i_Tp)%fragment(:,1),1)
       !
     endif
     !
     !
     !
   case("C")
     !
     if (l_BS_ares_from_res) then
       !
       ! Coupling part
       !
       call M_by_V('n',BS_blk(i_B)%size(1),BS_blk(i_B)%size(2),-alpha,&
&                      conjg(BS_blk(i_B)%mat),BS_blk(i_B)%size(1),conjg(Vi(i_Tp)%fragment(:,1)),1,cONE,&
&                      Vo(i_Tk)%fragment(:,1),1)
       !
       ! Coupling part can be symmetric
       !
       if ( i_Tk.ne.i_Tp ) then
         call M_by_V('t',BS_blk(i_B)%size(1),BS_blk(i_B)%size(2),-alpha,&
&                        conjg(BS_blk(i_B)%mat),BS_blk(i_B)%size(1),conjg(Vi(i_Tk)%fragment(:,1)),1,cONE,&
&                        Vo(i_Tp)%fragment(:,1),1)
       endif
       !
     else
       !
       ! Coupling part
       !
       call M_by_V('n',BS_blk(i_B)%size(1),BS_blk(i_B)%size(2),-cONE,&
&                      conjg(BS_blk(i_B)%mat),BS_blk(i_B)%size(1),Vi(i_Tp)%fragment(:,1),1,cONE,&
&                      Vo(i_Tk)%fragment(:,1),1)
       !
       ! Anti-coupling Q = -(-(C)*)^dag=C^t
       !
       call M_by_V('t',BS_blk(i_B)%size(1),BS_blk(i_B)%size(2),cONE,&
&                      BS_blk(i_B)%mat,BS_blk(i_B)%size(1),Vi(i_Tk)%fragment(:,1),1,cONE,&
&                      Vo(i_Tp)%fragment(:,1),1)
       !
     endif
     !
   case("A")
     !
     ! Anti-Resonant part 
     !
     if (i_Tk==i_Tp) then
       !
       call M_by_V('t',BS_blk(i_B)%size(1),BS_blk(i_B)%size(2),cONE,&
&                      BS_blk(i_B)%mat,BS_blk(i_B)%size(1),Vi(i_Tp)%fragment(:,1),1,cONE,&
&                      Vo(i_Tk)%fragment(:,1),1)
       !
       ! Add the energy diffs on the diagonal
       !
       Vo(i_Tk)%fragment(:,1) = Vo(i_Tk)%fragment(:,1) + BS_T_grp(i_Tk)%E(:,1)*Vi(i_Tk)%fragment(:,1)
       !
     else
       !
       call M_by_V('n',BS_blk(i_B)%size(1),BS_blk(i_B)%size(2),cONE,&
&                      conjg(BS_blk(i_B)%mat),BS_blk(i_B)%size(1),Vi(i_Tp)%fragment(:,1),1,cONE,&
&                      Vo(i_Tk)%fragment(:,1),1)
       !
       ! Anti-Resonant part is hermitian
       !
       call M_by_V('t',BS_blk(i_B)%size(1),BS_blk(i_B)%size(2),cONE,&
&                      BS_blk(i_B)%mat,BS_blk(i_B)%size(1),Vi(i_Tk)%fragment(:,1),1,cONE,&
&                      Vo(i_Tp)%fragment(:,1),1)
       !
     endif
     !
   end select
   !
   if (i_B>BSS_kernel_last_stored_block) call K_block_free(i_B,"matonly")
   !
 end do
 !
#if defined _MPI
 do i_g=1,BS_nT_grps
   if (.not.PAR_IND_T_Haydock%element_1D(i_g)) cycle
   call PP_redux_wait(Vo(i_g)%fragment,COMM=PAR_COM_T_Haydock(i_g)%COMM)
 enddo
#endif
 !
end subroutine
