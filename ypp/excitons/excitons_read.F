!
! License-Identifier: GPL
!
! Copyright (C) 2019 The Yambo Team
!
! Authors (see AUTHORS file for details): AM DV DS
!
subroutine excitons_read(k,Xk,en,Xen,iq,what,io_err)
 !
 use pars,          ONLY:pi
 use R_lattice,     ONLY:bz_samp
 use stderr,        ONLY:intc
 use X_m,           ONLY:X_reset
 use QP_CTL_m,      ONLY:QP_apply
 use electrons,     ONLY:levels,n_spin,n_sp_pol,spin
 use IO_int,        ONLY:io_control,io_disconnect
 use IO_m,          ONLY:OP_RD_CL,DUMP,NONE,REP,VERIFY
 use BS_solvers,    ONLY:BSS_eh_table,BSS_write_eig_2_db,io_BSS_diago,BS_mat,BSS_n_eig
 use BS,            ONLY:BS_H_dim,BS_K_dim,BS_bands,BS_K_coupling,BSqpts
 use YPPm,          ONLY:BS_R_left,BS_R_right,BS_E,BS_E_SOC_corr,BS_R_PL,BS_R_kerr,Xbsk   
 !
#include<y_memory.h>
 !
 type(bz_samp) ::Xk,k
 type(levels)  ::Xen,en
 integer, intent(in) :: iq
 integer, intent(out):: io_err
 character(*)        :: what
 !
 ! Work Space
 !
 !
 !... I/0
 !
 integer           ::io_ID,i_sec
 integer, external ::io_BS,io_E_SOC_map
 !
 integer                  :: neh
 logical                  :: l_check,l_eigenval,l_eigenvec,l_clean
 !
 l_check   =index(what,"check")        >0
 l_eigenval=index(what,"eigenvalues")  >0
 l_eigenvec=index(what,"eigenvectors") >0
 l_clean   =index(what,"clean")        >0
 !
 if(l_check) then
   !
   ! E SOC perturbative
   ! 
   if(n_spin==1) then
     call io_control(ACTION=OP_RD_CL,COM=NONE,SEC=(/1/),MODE=DUMP,ID=io_ID)
     io_err=io_E_SOC_map(en,k,io_ID)
   endif
   !
   YAMBO_ALLOC(BSqpts,(3,iq))
   !
   ! X type. Also used in PARALLEL setup
   !
   Xbsk%whoami=2
   !
   ! Here I read information from the BSS_diago database
   !
   call io_control(ACTION=OP_RD_CL,COM=NONE,MODE=DUMP,SEC=(/1/),ID=io_ID)
   io_err=io_BSS_diago(iq,1,io_ID,Xbsk)
   !
   if (io_err/=0) then
     call warning('ndb.BS_diago not found')
     return
   endif
   !
   ! Apply quasi-particle correction if presents
   ! 
   call QP_apply(BS_bands,Xen,Xk,'G',msg_fmt='rs')
   !
   BS_H_dim=BS_K_dim(1)
   if (BS_K_coupling) BS_H_dim=2*BS_K_dim(1)
   !
 endif
 !
 if(l_eigenval) then
   !
   ! Allocation of energies and residuals, then also eigenvector
   !
   YAMBO_ALLOC_P(BS_R_right,(BSS_n_eig))
   YAMBO_ALLOC_P(BS_R_left,(BSS_n_eig))
   YAMBO_ALLOC(BS_E,(BSS_n_eig))
   if(en%nb_SOC/=0) then
     YAMBO_ALLOC_P(BS_E_SOC_corr,(2/n_sp_pol,BSS_n_eig))
   endif
   !
   do i_sec=1,2
     call io_control(ACTION=OP_RD_CL,COM=REP,MODE=VERIFY,SEC=(/i_sec/),ID=io_ID)
     io_err=io_BSS_diago(iq,1,io_ID,Xbsk,BS_E,BS_R_left,BS_R_right,BsE_corr=BS_E_SOC_corr,bsR_kerr=BS_R_kerr,bsR_pl=BS_R_PL)
     if (io_err/=0) then
       if(i_sec==1) then
         call warning('Error reading SEC1 in ndb.BSS_diago, this can be due to the missing ndb.dipoles or ndb.ems1')
       else
         call warning('ndb.BSS_diago, SEC '//trim(intc(i_sec))//', not compatible with ypp input')
       endif
       return
     endif
   enddo
   !
   call io_disconnect(ID=io_ID)
   !
 endif
 !
 if(l_eigenvec) then
   ! 
   if(.not.BSS_write_eig_2_db) then
     call error(' Diagonalization database does not contain excitonic eigenstates.'//&
&               ' Re-run yambo using the WRbsWF flag in the input file.')
   endif
   !
   ! Loading tables and eigenvectors
   !
   YAMBO_ALLOC(BSS_eh_table,(BS_H_dim,3+n_sp_pol-1))
   YAMBO_ALLOC(BS_mat,(BS_H_dim,BSS_n_eig))
   !
   call io_control(ACTION=OP_RD_CL,COM=REP,MODE=VERIFY,SEC=(/3/),ID=io_ID)
   io_err=io_BSS_diago(iq,1,io_ID,Xbsk,BS_E,BS_R_left,BS_R_right,&
   &                   BsE_corr=BS_E_SOC_corr,bsR_kerr=BS_R_kerr,bsR_pl=BS_R_PL,BS_mat=BS_mat)
   if (io_err/=0) call error(' ndb.BSS_diago, SEC 3, not compatible with ypp input')
   !
   ! Create the anti-resonant part of the eh_table
   !
   if (BS_K_coupling) then
     do neh=BS_K_dim(1)+1,BS_H_dim
       !
       BSS_eh_table(neh,:)=BSS_eh_table(neh-BS_K_dim(1),:)
       !
       ! invert conduction <-> valence order
       !
       BSS_eh_table(neh,2)=BSS_eh_table(neh-BS_K_dim(1),3)
       BSS_eh_table(neh,3)=BSS_eh_table(neh-BS_K_dim(1),2)
       !
     enddo
   endif
   !
 endif
 !
 if(l_clean) then
   !
   ! CLEAN
   !
   YAMBO_FREE_P(BS_R_left)
   YAMBO_FREE_P(BS_R_right)
   YAMBO_FREE(BS_E)
   YAMBO_FREE(BSS_eh_table)
   YAMBO_FREE(BS_mat)
   YAMBO_FREE(BSqpts)
   if(associated(BS_E_SOC_corr)) then
     YAMBO_FREE_P(BS_E_SOC_corr)
   endif
   call X_reset(Xbsk)
   !
 endif
 !
end subroutine
