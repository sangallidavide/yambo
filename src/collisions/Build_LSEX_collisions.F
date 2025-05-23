!
! License-Identifier: GPL
!
! Copyright (C) 2020 the YAMBO team
!
! Authors (see AUTHORS file for details): CA
! 
subroutine Build_LSEX_collisions(k,q,i_qp,COLL) 
 !
 ! This subroutine calculate the long-range part of the Screened exchange
 !
 use electrons,      ONLY:n_sp_pol
 use pars,           ONLY:SP,cI,cZERO
 use collision_ext,  ONLY:OSCLL,W,ng_oscll,COLL_bands
 use QP_m,           ONLY:QP_table,QP_nk,QP_n_states
 use R_lattice,      ONLY:bz_samp,qindx_S
 ! 
#include<y_memory.h>
 !
 type(bz_samp), intent(in)  :: k,q
 integer,       intent(in)  :: i_qp
 !
 ! Work Space
 !
 integer :: i_mp  ! self-energy external indexes
 integer :: i_np,iqbz,i_kmq ! self-energy internal indexes
 integer :: ig1,ig2
 integer :: i_n,i_m
 integer :: i_k
 !
 ! Collisions output
 !
 complex(SP), intent(out) :: COLL(COLL_bands(1):COLL_bands(2),COLL_bands(1):COLL_bands(2),QP_nk)
 !
 ! Coll indexes
 !
 i_n=QP_table(i_qp,1)
 i_m=QP_table(i_qp,2)
 i_k=QP_table(i_qp,3)
 !
 ! Build Collisions from the Oscillators
 !
 ! COLL(n,m,n',m',k,k-q) = \sum_{G,G'} \rho(n,n',k,q,G')*\rho^*(m,m',k,q,G)*W(q,G',G)
 !
 do i_mp=COLL_bands(1),COLL_bands(2)
    do i_np=COLL_bands(1),COLL_bands(2)
      !
      COLL(i_np,i_mp,:)=cZERO
      ! 
      do iqbz=1,q%nbz
         i_kmq   = k%sstar(qindx_S(i_k,iqbz,1),1)
         !
         ! Build Collision on the fly
         !
         do ig1=1,ng_oscll
           do ig2=1,ng_oscll
             COLL(i_np,i_mp,i_kmq)=COLL(i_np,i_mp,i_kmq)+OSCLL(i_n,i_np,i_k,iqbz,ig2) &
&                                                 *conjg(OSCLL(i_m,i_mp,i_k,iqbz,ig1))*W(iqbz,ig2,ig1) 
            enddo
         enddo
         !
       enddo
   enddo
 enddo
 !
end subroutine Build_LSEX_collisions
