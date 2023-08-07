!
! Copyright (C) 2000-2008 D. Sangalli and the YAMBO team 
!              http://www.yambo-code.org
! 
! This file is distributed under the terms of the GNU 
! General Public License. You can redistribute it and/or 
! modify it under the terms of the GNU General Public 
! License as published by the Free Software Foundation; 
! either version 2, or (at your option) any later version.
!
! This program is distributed in the hope that it will 
! be useful, but WITHOUT ANY WARRANTY; without even the 
! implied warranty of MERCHANTABILITY or FITNESS FOR A 
! PARTICULAR PURPOSE.  See the GNU General Public License 
! for more details.
!
! You should have received a copy of the GNU General Public 
! License along with this program; if not, write to the Free 
! Software Foundation, Inc., 59 Temple Place - Suite 330,Boston, 
! MA 02111-1307, USA or visit http://www.gnu.org/copyleft/gpl.txt.
!
subroutine K_doubles(Ken)
 !
 ! This subroutine computes and write to disk in ndb.BS_C_doubles
 ! the coupling part of the second order kernel:
 ! C^{i}_{sd} with s={c,v} ; d={c1,c2,v1,v2} ; i=1,8
 ! see comments around lines 380 for a more detailed description of C
 !
 ! The subroutine K_doubles_solvers.F will read this file and pass the
 ! matrix BS_mat_doubles to the subroutine BS_doubles_mount.F which will
 ! compute the second order kernel for a specific frequency
 !
 ! WARNING: at present coded only for the k=q=0 case
 !
 use pars,           ONLY:cZERO,SP,pi
 use drivers,        ONLY:l_2nd_born
 use electrons,      ONLY:levels,n_sp_pol,spin_occ,spin
 use wave_func,      ONLY:WF
 use interfaces,     ONLY:WF_load,WF_free
 use parallel_int,   ONLY:PP_wait,PP_redux_wait,PARALLEL_index
 use parallel_m,     ONLY:pp_indexes,myid,pp_indexes_reset
 use com,            ONLY:msg
 use LIVE_t,         ONLY:live_timing
 use wrapper,        ONLY:V_dot_V,Vstar_dot_V
 use BS,             ONLY:BS_bands,BS_W,BS_W_is_diagonal,BS_n_g_W,BS_K_dim
 use BS_solvers,     ONLY:BSS_eh_table
 use BS_doubles,     ONLY:BS_mat_doubles,BS_res_K_doubles,BS_K_doubles_dim,  &
&                         BS_doubles_table,BS_bands_doubles,BS_doubles_E,    &
&                         doubles_corr,static_approx
 !use X_m,           ONLY:DIP_iq_dot_r
 !use R_lattice,      ONLY:bare_qpg
 use collision_el,   ONLY:elemental_collision,elemental_collision_free,elemental_collision_alloc
 use IO_int,         ONLY:io_control
 use IO_m,           ONLY:OP_RD_CL,REP,VERIFY,OP_WR_CL
 !
 implicit none
 !
 type(levels)     :: Ken  
 type(pp_indexes) :: px
 !
 ! Main loop variables
 !
 integer    :: iv1,ic1,iv2,ic2,iv3,ic3,icv1,icv2,icv3,icin,&
&              iOvv,iOcc,i1,i2,i3,ncin,&
&              icv,ic,iv,nc,nv,nv_doubles,nc_doubles,i_doubles
 complex(SP),allocatable :: Interaction_cccv(:,:,:,:),Interaction_vvcv(:,:,:,:)
 !
 ! Elements for the construction of the Kernel:
 ! 1 - interaction
 type(elemental_collision) ::isc
 ! 2 - oscillators
 complex(SP),allocatable::Ovv(:,:,:),Occ(:,:,:),Ocv(:,:,:)
 !
 ! I/O
 integer           ::io_err,ID     
 integer, external ::io_BS_doubles
 !
 call K_doubles_init() 
 !
 if(.not.BS_res_K_doubles) return
 !
#if defined _DEBUG
 open(unit=3456, file='debug.txt', status='unknown')
 open(unit=3457, file='interaction_cccv.txt', status='unknown')
 open(unit=3458, file='interaction_vvcv.txt', status='unknown')
 open(unit=3459, file='BS_mat_doubles.txt', status='unknown')
#endif
 !
 call section('=','Bethe-Salpeter Kernel + DOUBLES')
 !
 ! Check for the doubles DB, if already present exit from the subroutine
 call io_control(ACTION=OP_RD_CL,COM=REP,MODE=VERIFY,SEC=(/1/),ID=ID)
 io_err=io_BS_doubles(ID)
 if (io_err==0) then
   return
 endif
 !
 ! *************************
 ! ****** Allocations ******
 ! *************************
 call K_doubles_alloc(1) 
 !
 ! ********************************
 ! ***** Before the main loop *****
 ! ********************************
 ! Doubles tables for indexes and energies
 i_doubles=0
 icv1=0
 do iv1=BS_bands_doubles(1),Ken%nbf(1)
   do ic1=Ken%nbf(1)+1,BS_bands_doubles(2)
     icv1=icv1+1
     do iv2=BS_bands_doubles(1),Ken%nbf(1)
       do ic2=Ken%nbf(1)+1,BS_bands_doubles(2)
         i_doubles=i_doubles+1
         BS_doubles_table(i_doubles,:)=(/iv1,iv2,ic1,ic2/)
	 !BS_doubles_E(i_doubles)=Ken%E(ic1,1,1)+Ken%E(ic2,1,1)-Ken%E(iv1,1,1)-Ken%E(iv2,1,1)
         BS_doubles_E(i_doubles)=Ken%Eo(ic1,1,1)+Ken%Eo(ic2,1,1)-Ken%Eo(iv1,1,1)-Ken%Eo(iv2,1,1)
         ! STATIC APPROX <
         ! Check for BS_doubles_E. If greater than ?? eV then the corrisponding 
         ! index is recorded so that will be possible to do a static approximation
         ! for the kernel for that index
         !if ( BS_doubles_E(i_doubles)*27.2116> ?? ) then
         !  static_approx(i_doubles)=1
#if defined _DEBUG
         !  write(3456,*) i_doubles
#endif
         !endif
         ! STATIC APPROX >
#if defined _DEBUG
         !write(3456,*) ic1,ic2,iv1,iv2,BS_doubles_E(i_doubles)*27.2116
#endif
       enddo
     enddo
   enddo
 enddo
 !
 ! Load wave funtions
 call WF_load(WF,BS_n_g_W,1,(/BS_bands_doubles(1),BS_bands_doubles(2)/),(/1,1/),&
&             space='R',title='-BSK_doubles')
 !
 ! Oscillators construction
 ! 1 - preparations lines
 call elemental_collision_free(isc)
 isc%is(:)=1
 isc%os(:)=1
 isc%qs(:)=1
 isc%ngrho=BS_n_g_W
 isc%iqref=1
 !
 call elemental_collision_alloc(isc,NG=BS_n_g_W,NG_GAMP=(/BS_n_g_W,1/),TITLE='2ndBorn')
 !
 ! 2 - Ovv construction
 ! Paral <
 call pp_indexes_reset(px)
 call PARALLEL_index(px,(/Ken%nbf(1),Ken%nbf(1)/),(/BS_bands(1),BS_bands_doubles(1)/))
 call PP_wait()
 call live_timing('Ovv',px%n_of_elements(myid+1))
 ! Paral >
 !
 do iv=BS_bands(1),Ken%nbf(1)
   do iv1=BS_bands_doubles(1),Ken%nbf(1)
     ! Paral <
     if (.not.px%element_2D(iv,iv1)) cycle
     call live_timing(steps=1)
     ! Paral >
     isc%is(1)=iv1
     isc%os(1)=iv
     call scatter_Bamp(isc)
     forall(i1=1:BS_n_g_W) Ovv(i1,iv1-BS_bands_doubles(1)+1,iv-BS_bands(1)+1)=isc%rhotw(i1)
#if defined _DEBUG
     write(3456,*) iv1,iv,Ovv(:,iv1,iv-BS_bands(1)+1)
#endif
   enddo
 enddo
 !
 ! Paral <
 call live_timing()
 call PP_redux_wait(Ovv)
 call PP_indexes_reset(px)
 !
 ! 3 -Occ construction
 call PARALLEL_index(px,(/BS_bands(2),BS_bands_doubles(2)/),(/1+Ken%nbf(1),1+Ken%nbf(1)/))
 call PP_wait()
 call live_timing('Occ',px%n_of_elements(myid+1))
 ! Paral >
 !
 do ic=1+Ken%nbf(1),BS_bands(2)
   do ic1=1+Ken%nbf(1),BS_bands_doubles(2)
     ! Paral <
     if (.not.px%element_2D(ic,ic1)) cycle
     call live_timing(steps=1)
     ! Paral >
     isc%is(1)=ic
     isc%os(1)=ic1
     call scatter_Bamp(isc)
     forall(i1=1:BS_n_g_W) Occ(i1,ic-Ken%nbf(1),ic1-Ken%nbf(1))=isc%rhotw(i1)
#if defined _DEBUG
     write(3456,*) ic,ic1,Occ(:,ic-Ken%nbf(1),ic1-Ken%nbf(1))
#endif
   enddo
 enddo
 !
 ! Paral <
 call live_timing()
 call PP_redux_wait(Occ)
 call PP_indexes_reset(px)
 !
 ! 4 -Ocv construction
 call PARALLEL_index(px,(/Ken%nbf(1),BS_bands_doubles(2)/),(/BS_bands_doubles(1),1+Ken%nbf(1)/))
 call PP_wait()
 call live_timing('Ocv',px%n_of_elements(myid+1))
 ! Paral >
 !
 do iv2=BS_bands_doubles(1),Ken%nbf(1)
   do ic2=1+Ken%nbf(1),BS_bands_doubles(2)
     ! Paral <
     if (.not.px%element_2D(iv2,ic2)) cycle
     call live_timing(steps=1)
     ! Paral >
     isc%is(1)=ic2
     isc%os(1)=iv2
     call scatter_Bamp(isc)
     forall(i1=1:BS_n_g_W) Ocv(i1,ic2-Ken%nbf(1),iv2-BS_bands_doubles(1)+1)=isc%rhotw(i1)
     ! Warning, should evaluate the term  G=0 of Ocv but use this only in the 
     ! bubble diagrams when c=c' ; v=v'. Presently Ocv=0 for G=0
!     Ocv(1,ic2-Ken%nbf(1),iv2-BS_bands_doubles(1)+1)=-conjg(DIP_iq_dot_r(ic2,iv2,1,1))/bare_qpg(1,1)
#if defined _DEBUG
     !write(3456,*) ic2,iv2,Ocv(:,ic2-Ken%nbf(1),iv2)
#endif
   enddo
 enddo
 !
 ! Paral <
 call live_timing()
 call PP_redux_wait(Ocv)
 call PP_indexes_reset(px)
 ! Paral >
 ! 
 call elemental_collision_free(isc)
 !
 call WF_free(WF)
 !
 !  ********************************************
 !  *******       Main loop: start       *******
 !  ********************************************
 if (BS_res_K_doubles) then
   !
   !  Step 1: Here the code constructs U_{cc1,c2v2} and U_{v1v,c2v2} 
   !          starting Ocv, Occ and Ovv, and BS_W
   !
   !          1.a: U_{c c1 c2 v2} = Occ * BS_W * Ocv
   !
   ! Paral <
   call pp_indexes_reset(px)
   call PARALLEL_index(px,(/nc,nc_doubles/),(/1,1/))
   call PP_wait()
   call live_timing('Interaction_cccv',px%n_of_elements(myid+1))
   ! Paral >
   !
   do ic=1,nc
     do ic1=1,nc_doubles
       !
       ! Paral <
       if (.not.px%element_2D(ic,ic1)) cycle
       call live_timing(steps=1)
       ! Paral >
       if (BS_W_is_diagonal) then
         forall(i1=1:BS_n_g_W) Occ(i1,ic,ic1)=Occ(i1,ic,ic1)*BS_W(i1,1,1)          
       else
         !
         do i1=1,BS_n_g_W 
           Occ(i1,ic,ic1)=V_dot_V(BS_n_g_W,Occ(:,ic,ic1),BS_W(:,i1,1))  
         enddo
       endif
       do iv2=1,nv_doubles
         do ic2=1,nc_doubles
           !
           Interaction_cccv(ic,ic1,ic2,iv2)=Vstar_dot_V(BS_n_g_W,Ocv(:,ic2,iv2),Occ(:,ic,ic1))*4._SP*pi
	   ! 
#if defined _DEBUG
           write(3457,*) ic+Ken%nbf(1),ic1+Ken%nbf(1),ic2+Ken%nbf(1),iv2+BS_bands_doubles(1)-1,&
&                        Interaction_cccv(ic,ic1,ic2,iv2)
#endif
	   !
         enddo
       enddo
     enddo
   enddo
   !
   ! Paral <
   call live_timing()
   call pp_redux_wait(Interaction_cccv)
   call pp_indexes_reset(px)
   !
   call PARALLEL_index(px,(/nv,nv_doubles/),(/1,1/))
   call PP_wait
   call live_timing('Interaction_vvcv',px%n_of_elements(myid+1))
   ! Paral >
   !
   !          1.b: U_{v v1 c2 v2} = Ovv * BS_W * Ocv
   !
   do iv=1,nv
     do iv1=1,nv_doubles
       !
       ! Paral <
       if (.not.px%element_2D(iv,iv1)) cycle
       call live_timing(steps=1)
       ! Paral >
       if (BS_W_is_diagonal) then
         forall(i1=1:BS_n_g_W) Ovv(i1,iv1,iv)=Ovv(i1,iv1,iv)*BS_W(i1,1,1) 
       else
         !
         do i1=1,BS_n_g_W 
           Ovv(i1,iv1,iv)=V_dot_V(BS_n_g_W,Ovv(:,iv1,iv),BS_W(:,i1,1))
         enddo
       endif
       do iv2=1,nv_doubles
         do ic2=1,nc_doubles
           Interaction_vvcv(iv1,iv,ic2,iv2)=Vstar_dot_V(BS_n_g_W,Ocv(:,ic2,iv2),Ovv(:,iv1,iv))*4._SP*pi
	   ! 
#if defined _DEBUG
           write(3458,*) iv+BS_bands(1)-1,iv1+BS_bands_doubles(1)-1,ic2+Ken%nbf(1),iv2+BS_bands_doubles(1)-1,&
&                        Interaction_vvcv(iv1,iv,ic2,iv2)
#endif
           !
         enddo
       enddo
     enddo
   enddo 
   !
   ! Paral <
   call live_timing()
   call pp_redux_wait(Interaction_vvcv)
   call pp_indexes_reset(px)
   !
 endif
 !
 call K_doubles_free(1) 
 !
 if (BS_res_K_doubles) then
   !
   call K_doubles_alloc(2)
   !
   ! TMP coding, needed for BSS_eh_table
   call K_components_folded_in_serial_arrays(1)
   !
   !  Step 2: Here the code constructs the half-diagrams
   !          D1(icv,i_doubles) and D2(icv,i_doubles),
   !          starting from U_{cc1,c2v2} and U_{v1v,c2v2}
   !
   !  The code is just for zero momentum. However here I draw the
   !  general case at finite momentum. There are 6 momenta index,
   !  e.g. 2 more compared to SBSE, but also two extra conservation,
   !  which brings all downd to the same number of momenta
   !
   !                             p2=k2-q_W
   !
   !               D1                             D2
   !    (v1=v; k1=k; p1=k-q_W)          (c1=c; p1=k; k1=k-q_W)
   !
   !   c,k ===<===:===<=== c1,p1       c,k ===<===x===<===  c1,p1
   !              :
   !         q_W /:\                                ==<===  c2,k2
   !              :                                /
   !              :===<=== c2,k2                  :===>===  v2,p2
   !               \                              : 
   !                ==>=== v2,p2               q \:/
   !                                              :
   ! v,k-q ===>===x===>=== v1,k1-q   v,k-q ===>===:===>===  v1,k1-q
   !
   !          We have, for D1 v=v1, k=k1 ; 
   !
   !          and the final simmetrized expression:
   !          C_{sd}=0.5*[D1+D2 -(c1<->c2) -(v1<->v2) +(c1,v1<->c2,v2)]
   !          Here the prefactor 0.5 to each half-diagram is needed as we have 8 half-diagrams
   !          which, composed, will give 8*8=64 diagrams, but equal at groups of 4
   !          (i.e. D1*D2 = D1(c1<->c2)*D2(c1<->c2) = ...) so we have to divide by 4(=0.5*0.5)
   !
   !  In the code:
   !  BS_mat_doubles(:,:,1)=sqrt(0.5)*[D1+D2 +(c1,v1<->c2,v2)]*sqrt(2)                 spin:(up,down)
   !  BS_mat_doubles(:,:,2)=sqrt(0.5)*[D1+D2 +(c1,v1<->c2,v2) -(c1<->c2) -(v1<->v2) ]  spin:(up,up)
   !
#if defined _DEBUG
   write(3459,*) '                                                                       ',&
&                '       kin1            kin2            kin3            kin4     ',&
&                '       kin5            kin6            kin7            kin8     ',&
&                '       kin1            kin2            kin3            kin4     ',&
&                '       kin5            kin6            kin7            kin8     '
#endif
   !
   ! Paral <
   call PARALLEL_index(px,(/4*ncin,BS_K_doubles_dim/))
   call PP_wait()
   call live_timing('BS_mat_doubles',px%n_of_elements(myid+1))
   ! Paral >
   do icin=1,4*ncin
     do i_doubles=1,BS_K_doubles_dim
       ! Paral <
       if (.not.px%element_2D(icin,i_doubles)) cycle
       call live_timing(steps=1)
       ! Paral >
       iv1=BS_doubles_table(i_doubles,1)-BS_bands_doubles(1)+1
       iv2=BS_doubles_table(i_doubles,2)-BS_bands_doubles(1)+1
       ic1=BS_doubles_table(i_doubles,3)-Ken%nbf(1)
       ic2=BS_doubles_table(i_doubles,4)-Ken%nbf(1)
       !
       do icv=1,BS_K_dim(1)
         iv=BSS_eh_table(icv,2)-BS_bands(1)+1
         ic=BSS_eh_table(icv,3)-Ken%nbf(1)
         select case (icin)
           case(1)
             if ( (iv+BS_bands(1)) /= (iv1+BS_bands_doubles(1)) ) cycle
             BS_mat_doubles(icv,i_doubles,ncin+1)=&
&                   Interaction_cccv(ic,ic1,ic2,iv2)                  !D1
           case(2)
             if ( (iv+BS_bands(1)) /= (iv2+BS_bands_doubles(1)) ) cycle 
             BS_mat_doubles(icv,i_doubles,ncin+1)=&
&                   Interaction_cccv(ic,ic2,ic1,iv1)                  !D1cv
           case(3)
             if (ic/=ic1) cycle
             BS_mat_doubles(icv,i_doubles,ncin+1)=&
&                   Interaction_vvcv(iv1,iv,ic2,iv2)                  !D2
           case(4)
             if (ic/=ic2) cycle 
             BS_mat_doubles(icv,i_doubles,ncin+1)=&
&                   Interaction_vvcv(iv2,iv,ic1,iv1)                  !D2cv
           case(5)
             if ( (iv+BS_bands(1)) /= (iv1+BS_bands_doubles(1)) ) cycle
             BS_mat_doubles(icv,i_doubles,ncin+1)=&
&                  -Interaction_cccv(ic,ic2,ic1,iv2)                  !D1c
           case(6)
             if ( (iv+BS_bands(1)) /= (iv2+BS_bands_doubles(1)) ) cycle 
             BS_mat_doubles(icv,i_doubles,ncin+1)=&
&                  -Interaction_cccv(ic,ic1,ic2,iv1)                  !D1v
           case(7)
             if (ic/=ic1) cycle
             BS_mat_doubles(icv,i_doubles,ncin+1)=&
&                  -Interaction_vvcv(iv2,iv,ic2,iv1)                  !D2v
           case(8)
             if (ic/=ic2) cycle
             BS_mat_doubles(icv,i_doubles,ncin+1)=&
&                  -Interaction_vvcv(iv1,iv,ic1,iv2)                  !D2c
         end select
#if defined _DEBUG
       write(3459,*) ic+Ken%nbf(1),iv+BS_bands(1)-1,&
&                    ic1+Ken%nbf(1),ic2+Ken%nbf(1),iv1+BS_bands_doubles(1)-1,iv2+BS_bands_doubles(1)-1,&
&                    (real(BS_mat_doubles(icv,i_doubles,i1)),i1=1,8),&
&                    (aimag(BS_mat_doubles(icv,i_doubles,i1)),i1=1,8)
#endif
       enddo
     enddo
     !
     ! sqrt(0.5) is due to second order
     ! For BS_mat_doubles_1 we obtain an overall factor 0.5
     ! For BS_mat_doubles_2 we obtain an overall factor 0.5
     if(icin<=4) BS_mat_doubles(:,:,1)=BS_mat_doubles(:,:,1)+0.5*BS_mat_doubles(:,:,ncin+1)
     if(.not.doubles_corr) BS_mat_doubles(:,:,2)=&
&                    BS_mat_doubles(:,:,2)+0.5*sqrt(0.5)*BS_mat_doubles(:,:,ncin+1)
     BS_mat_doubles(:,:,ncin+1)=(0.,0.)
   enddo
   !
   ! Paral <
   call pp_redux_wait(BS_mat_doubles)
   call live_timing()
   call pp_indexes_reset(px)
   ! Paral >
   !
   !if (doubles_corr) call K_doubles_correlated()
   !
   ! DS 2023/08/07
   ! Instead of introducing a new ndb, if would be better to directly do
   ! the I/O of the coulomb integrals, called Interaction_ccvc, etc, here
   !
   ! Write BS_mat_doubles to disk
   call io_control(ACTION=OP_WR_CL,COM=REP,SEC=(/1,2/),ID=ID)
   io_err=io_BS_doubles(ID)
   !
   call K_doubles_free(2)
   !
 endif
 !
 !  ********************************************
 !  *******       Main loop: end         *******
 !  ********************************************
 !
#if defined _DEBUG
 close(unit=3456)
 close(unit=3457)
 close(unit=3458)
 close(unit=3459)
#endif
!
contains
 !
 subroutine K_doubles_init()
  !
  implicit none
  !
  ! Davide doubles
  ! preparation lines:
  BS_res_K_doubles=l_2nd_born
  ! For now only standard 2nd Born
  doubles_corr=.true.
  !doubles_corr=.false.
  BS_bands_doubles(1)=BS_bands(1)
  BS_bands_doubles(2)=BS_bands(2)
  !
  nc=BS_bands(2)-Ken%nbf(1)
  nv=Ken%nbf(1)-BS_bands(1)+1
  nc_doubles=BS_bands_doubles(2)-Ken%nbf(1)
  nv_doubles=Ken%nbf(1)-BS_bands_doubles(1)+1
  !
  BS_K_doubles_dim=(nc_doubles*nv_doubles)**2
  !
  ncin=2
  if(doubles_corr) ncin=1
  !
  ! Check fo the dimension of the kernel
  if (BS_K_doubles_dim==0) call error(' Null BSE doubles kernel block dimension(s) found. Increase doubles range')
  !
  call msg('rsn','[BSE] Doubles dimension      ',BS_K_doubles_dim)
  !
 end subroutine K_doubles_init
 !
 !
 subroutine K_doubles_alloc(isec)
  !
#include<memory.h>
  !
  integer, intent(in) :: isec
  !
  select case(isec)
  case(1)
    ! Allocation: Table and energies
    YAMBO_ALLOC(BS_doubles_table,(BS_K_doubles_dim,4))
    BS_doubles_table=0
    !
    allocate(BS_doubles_E(BS_K_doubles_dim))
    BS_doubles_E=cZERO
    !
    !allocate(static_approx(BS_K_doubles_dim))
    !static_approx=0
    !
    ! Allocation: Oscillators
    YAMBO_ALLOC(Ovv,(BS_n_g_W,nv_doubles,nv))
    YAMBO_ALLOC(Occ,(BS_n_g_W,nc,nc_doubles))
    YAMBO_ALLOC(Ocv,(BS_n_g_W,nc_doubles,nv_doubles))
    Occ=CZERO
    Ovv=cZERO
    Ocv=cZERO
    !
    ! Allocation: Interaction and BS_mat_doubles
    YAMBO_ALLOC(Interaction_cccv,(nc,nc_doubles,nc_doubles,nv_doubles))
    YAMBO_ALLOC(Interaction_vvcv,(nv_doubles,nv,nc_doubles,nv_doubles))
    Interaction_cccv=cZERO
    Interaction_vvcv=cZERO
  case(2)
    YAMBO_ALLOC(BS_mat_doubles,(BS_K_dim(1),BS_K_doubles_dim,ncin+1))
    BS_mat_doubles=cZERO
  end select
  !
 end subroutine K_doubles_alloc
 !
 !
 subroutine K_doubles_free(isec)
  !
#include<memory.h>
  !
  integer, intent(in) :: isec
  !
  select case(isec)
  case(1)
    YAMBO_FREE(Occ)
    YAMBO_FREE(Ovv)
    YAMBO_FREE(Ocv)
    !
  case(2)
    ! Final deallocation
    YAMBO_FREE(Interaction_cccv)
    YAMBO_FREE(Interaction_vvcv)
    YAMBO_FREE(BS_mat_doubles)
    YAMBO_FREE(BS_doubles_table)
    deallocate(BS_doubles_E)
    !deallocate(static_approx)
    !
  end select
  !
 end subroutine K_doubles_free
 !
end subroutine K_doubles
