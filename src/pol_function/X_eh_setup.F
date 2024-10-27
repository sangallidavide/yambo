!
! License-Identifier: GPL
!
! Copyright (C) 2007 The Yambo Team
!
! Authors (see AUTHORS file for details): AM DS
!
integer function X_eh_setup(iq,X,Xen,Xk,minmax_ehe)
 !
 use pars,          ONLY:SP
 use X_m,           ONLY:X_t,X_poles,X_Ein_poles,X_poles_tab,l_X_terminator,X_term_E,&
 &                       l_X_abs,l_X_magn
 use R_lattice,     ONLY:qindx_X,bz_samp
 use electrons,     ONLY:levels,n_sp_pol,spin_occ,i_spin_majority,get_spin_majority
 use frequency,     ONLY:ordered_grid_index
 use parallel_m,    ONLY:PAR_IND_CON_BANDS_X,PAR_IND_Xk_bz,PAR_IND_VAL_BANDS_X
 !
#include<y_memory.h>
 !
 type(levels)  ::Xen 
 type(bz_samp) ::Xk
 type(X_t)     ::X
 integer       ::iq
 real(SP)      ::minmax_ehe(2)
 !
 ! Work Space
 !
 integer :: ik_bz,i_sym,ik,ik_m_q,iv,ic,i_sp,c_sp_pol,v_sp_pol,i_pole,N_tot, &
 &          Nv(n_sp_pol),Nc(n_sp_pol),VB_lims(2,n_sp_pol),CB_lims(2,n_sp_pol)
 logical :: MPI_logicals_on,l_XTerm_vv
 logical :: l_res,l_ares_T,l_ares_R
 real(SP):: E_eh,E_in,f_eh,local_tresh
 real(SP), allocatable :: poles(:),Ein_poles(:)
 !
 l_res   = X%ordering=="R".or.X%ordering=="T".or.X%ordering=="r"
 l_ares_T= X%ordering=="T".or.X%ordering=="Ta"
 l_ares_R= X%ordering=="R".or.X%ordering=="Ra"
 !
 VB_lims(1,:)=1
 CB_lims(1,:)=1
 VB_lims(2,:)=0
 CB_lims(2,:)=0
 !
 if (l_X_magn) then
   !
   call get_spin_majority(Xen)
   !
   ! 2024/02/28 DS
   ! In the magnons case I loop on both the resonant and the anti-resonant transitions.
   ! In the tag phys-magn-v0 I had coded an implementation where
   ! the valence and conduction limits are selected as usually, doing a loop only on the resonant part.
   ! In such implementation, the anti-resonant term is included via the opposite spin transitions.
   ! Such implementation is significantly more complex, but also more efficient in distributing
   ! WFs in memory over bands (see new changes in PARALLEL_global_indexes in the present branch)
   ! The structure of such implementation is also more similar to the one in the BSE case
   !
   if (l_ares_T.or.l_ares_R) then
     do i_sp=1,n_sp_pol
       if(i_sp==i_spin_majority) cycle
       VB_lims(:,i_sp)=(/X%ib(1),X%ib(2)/)
       CB_lims(:,i_sp)=(/X%ib(1),X%ib(2)/)
     enddo
   else
     do i_sp=1,n_sp_pol
       if(i_sp==i_spin_majority) cycle
       c_sp_pol=i_sp
       v_sp_pol=mod(c_sp_pol,n_sp_pol)+1
       VB_lims(:,i_sp)=(/X%ib(1),Xen%nbm(v_sp_pol)/)
       CB_lims(:,i_sp)=(/Xen%nbf(c_sp_pol)+1,X%ib(2)/)
     enddo
   endif
   !
 else
   !
   do i_sp=1,n_sp_pol
     VB_lims(:,i_sp)=(/X%ib(1),Xen%nbm(i_sp)/)
     CB_lims(:,i_sp)=(/Xen%nbf(i_sp)+1,X%ib(2)/)
   enddo
   if (l_X_terminator) then
     CB_lims(1,:)=X%ib(1)
   endif
   !
 endif
 !
 N_tot=0
 do i_sp=1,n_sp_pol
   Nv(i_sp)=VB_lims(2,i_sp)-VB_lims(1,i_sp)+1
   Nc(i_sp)=CB_lims(2,i_sp)-CB_lims(1,i_sp)+1
   N_tot=N_tot+Nc(i_sp)*Nv(i_sp)
 enddo
 !
 YAMBO_ALLOC(poles,(Xk%nbz*N_tot))
 YAMBO_ALLOC(Ein_poles,(Xk%nbz*N_tot))
 !
 MPI_logicals_on=allocated(PAR_IND_Xk_bz%element_1D)
 !
 ! Note that in metals not all the elements are defined
 !
 poles=0._SP 
 Ein_poles=0._SP
 !
 X_eh_setup=0
 local_tresh=epsilon(1._SP)
 !
 do i_sp=1,n_sp_pol
   !
   c_sp_pol = i_sp
   v_sp_pol = i_sp
   !
   if (l_X_magn.and.n_sp_pol==2) then
     if(c_sp_pol==i_spin_majority) cycle
     v_sp_pol=mod(c_sp_pol,n_sp_pol)+1
   endif
   !
   do ik_bz=1,Xk%nbz
     !
     if (MPI_logicals_on) then
       if (.not.PAR_IND_Xk_bz%element_1D(ik_bz)) cycle 
     endif
     !
     do iv=VB_lims(1,i_sp),VB_lims(2,i_sp)
       !
       if (MPI_logicals_on) then
         if (.not.PAR_IND_VAL_BANDS_X(X%whoami)%element_1D(iv)) cycle 
       endif
       !
       do ic=CB_lims(1,i_sp),CB_lims(2,i_sp)
         !
         if (MPI_logicals_on) then
           if (.not.PAR_IND_CON_BANDS_X(X%whoami)%element_1D(ic)) cycle 
         endif
         !
         i_sym  =qindx_X(iabs(iq),ik_bz,1)
         ik     =Xk%sstar(ik_bz,1) 
         ik_m_q =Xk%sstar(i_sym,1)
         !
         ! terminator handling
         if (ic>=X%ib(1).and.ic<=Xen%nbm(c_sp_pol).and.l_X_terminator) then
           l_XTerm_vv=.TRUE.
           E_eh=X_term_E-Xen%E(iv,ik_m_q,v_sp_pol)       
         else
           l_XTerm_vv=.FALSE.
           E_eh=Xen%E(ic,ik,c_sp_pol)-Xen%E(iv,ik_m_q,v_sp_pol)
         endif 
         !
         E_in=Xen%E(iv,ik_m_q,v_sp_pol)
         !
         ! Note that all possible E_eh signs are accepted. Negative
         ! transitions energies appear at finite temperature.
         !
         ! The way to distinguish between resonant and anti-reonant transitions
         ! is to check fv(1-fc) factor that comes from the t>0 ordering
         ! of the G's function. In this way, however, E_eh can be negative as
         ! shown below
         !
         ! n(E) ----     
         !          *    
         !          4*  3
         !            -.   
         !              *
         !          1   2-._______
         !          -Eeh-
         !
         !   1,2 = v , 3,4 = c
         !
         !   both transitions 1->3 , 2->4 are resonant but 2->4 has negative energy
         !
         ! When iv--> iv transitions are counted, f_eh=0. 
         ! Therefore  if (abs(f_eh)<local_tresh) cycle hold only 
         ! when the terminator  is off 
         !
         if (.not.l_XTerm_vv) f_eh=Xen%f(iv,ik_m_q,v_sp_pol)*(spin_occ-Xen%f(ic,ik,c_sp_pol))/spin_occ
         if (     l_XTerm_vv) f_eh= Xen%f(iv,ik_m_q,v_sp_pol)*Xen%f(ic,ik,c_sp_pol) /spin_occ
         if (      l_X_magn ) f_eh=(Xen%f(iv,ik_m_q,v_sp_pol)-Xen%f(ic,ik,c_sp_pol))/spin_occ
         !
         if (abs(f_eh)<local_tresh)  cycle
         !                                                              
         if (any((/abs(E_eh)<X%ehe(1),abs(E_eh)>X%ehe(2).and.X%ehe(2)>0._SP/))) cycle
         !
         X_eh_setup=X_eh_setup+1
         poles(X_eh_setup)=E_eh
         Ein_poles(X_eh_setup)=E_in
         !
         i_pole=X_eh_setup
         if (allocated(ordered_grid_index)) i_pole=ordered_grid_index(X_eh_setup)
         if (iq>0) X_poles_tab(i_pole,:)=(/ik_bz,iv,ic,c_sp_pol,v_sp_pol/)
         !
       enddo
     enddo
   enddo
 enddo
 !
 minmax_ehe=(/max(minval(poles(:X_eh_setup))-0.1_SP,0._SP),maxval(poles(:X_eh_setup))+0.1_SP/)
 !
 if (iq<0) then
   if (.not.allocated(X_poles)) then
     YAMBO_ALLOC(X_poles,(X_eh_setup))
     X_poles=0._SP
   endif
   X_poles=X_poles+poles(:X_eh_setup)
   !
   if (.not.allocated(X_Ein_poles)) then
     YAMBO_ALLOC(X_Ein_poles,(X_eh_setup))
     X_Ein_poles=0._SP
   endif
   X_Ein_poles=X_Ein_poles+Ein_poles(:X_eh_setup)
 endif
 !
 YAMBO_FREE(poles)
 YAMBO_FREE(Ein_poles)
 !
end function
