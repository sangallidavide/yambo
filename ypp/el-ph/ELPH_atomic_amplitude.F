!
! License-Identifier: GPL
!
! Copyright (C) 2020 The Yambo Team
!
! Authors (see AUTHORS file for details): AM
!
subroutine ELPH_atomic_amplitude(q)
 !
 use pars,                ONLY:SP
 use units,               ONLY:BO2ANG,pi
 use com,                 ONLY:msg
 use R_lattice,           ONLY:bz_samp,RL_vol,RIM_n_rand_pts
 use D_lattice,           ONLY:n_atoms
 use ELPH,                ONLY:ph_modes,PH_freqs_sq,elph_nQ,elph_use_q_grid,&
&                              PH_pol_vector,PH_qpt
 use IO_int,              ONLY:io_control
 use IO_m,                ONLY:manage_action,OP_RD,RD_CL_IF_END,DUMP,RD_CL
 use stderr,              ONLY:intc
 use functions,           ONLY:bose_f
 use interfaces,          ONLY:LINEAR_ALGEBRA_driver
 use ALLOC,               ONLY:ELPH_alloc
 use linear_algebra,      ONLY:DIAGO,USE_LK
 use stderr,              ONLY:intc
 !
#include<y_memory.h>
 type(bz_samp)::q
 !
 ! Work Space 
 !
 integer           ::nq_to_sum,iq,im,is,i,j
 complex(SP)       ::mean_ui_uj(n_atoms,3,3)
 real(SP)          ::omega,U_eigenvalues(3)
 real(SP),allocatable :: q_weight(:)
 !
 !I/O
 !
 integer           ::io_err,ID,IO_ACT
 integer, external ::io_ELPH
 !
 call section("*","== Electron-Phonon Interface: mean atomic amplitudes ==")
 !---------------------------------------------------------------------------
 !
 call io_control(ACTION=OP_RD,SEC=(/1/),MODE=DUMP,ID=ID)
 io_err=io_ELPH(ID,'no_matrix_elements')
 !
 call msg('s',' Checking database ...')
 if (io_err<0) then
   call msg('s','not found')
   return
 endif
 call msg('s','sane')
 !
 nq_to_sum=elph_nQ
 if (elph_use_q_grid) nq_to_sum=q%nibz
 !
 ! calculate Q weight due to symmetry
 !
 call k_expand(q) 
 !
 ! ... and to 1/q divergence using the RIM integration
 !
 YAMBO_ALLOC(q_weight,(nq_to_sum))
 !
 RIM_n_rand_pts=100000
 if (elph_use_q_grid) then
   q_weight=q%weights
 else
   if (RIM_n_rand_pts>0) call msg('s','Using RIM with '//trim(intc(RIM_n_rand_pts))//' point')
   call rim_spherical(nq_to_sum,PH_qpt(1:nq_to_sum,1:3),q_weight,(3.*RL_vol/nq_to_sum/4./pi)**(1./3.),2,.TRUE.)
 endif
 !
 mean_ui_uj=(0.,0.)
 !
 do iq=1,nq_to_sum 
   !
   ! I/O
   !
   IO_ACT=manage_action(RD_CL_IF_END,iq,1,nq_to_sum)
   call io_control(ACTION=IO_ACT,SEC=(/iq+1/),ID=ID)
   io_err=io_ELPH(ID,'no_matrix_elements')  
   !
   if (.not.allocated(PH_pol_vector)) then
     call warning('Polarization vectors not found!')
     call ELPH_alloc('FREE')
     return
   endif
   !
   do im=1,ph_modes
     !
     omega = sqrt( MAX( PH_freqs_sq(iq,im),0. ) )
     !
     if (omega<1.E-10) cycle
     !
     do is=1,n_atoms
       do i=1,3
         do j=1,3
           mean_ui_uj(is,i,j)=mean_ui_uj(is,i,j)+q_weight(iq)*1./2./omega*&
&                             conjg(PH_pol_vector(im,is,i,iq))*PH_pol_vector(im,is,j,iq)*&
&                             (bose_f(omega)+1.)
           mean_ui_uj(is,i,j)=real ( mean_ui_uj(is,i,j) )
         enddo
       enddo
     enddo
     !
   enddo
   !
 enddo
 !
 do is=1,n_atoms
   !
   call LINEAR_ALGEBRA_driver(DIAGO,M=mean_ui_uj(is,:,:),lib_in=USE_LK,E_real=U_eigenvalues)
   !
   call msg('s','Atom '//intc(is))
   call msg('s','Displacement  [a.u.]/A #1 ',sqrt(U_eigenvalues(1))*(/1._SP,BO2ANG/) )
   call msg('s','Normalized direction   #1 ',real(mean_ui_uj(is,:,1)))
   call msg('s','                       #2 ', sqrt(U_eigenvalues(2))*(/1._SP,BO2ANG/) )
   call msg('s','                       #2 ',real(mean_ui_uj(is,:,2)))
   call msg('s','                       #3 ', sqrt(U_eigenvalues(3))*(/1._SP,BO2ANG/) )
   call msg('s','                       #3 ',real(mean_ui_uj(is,:,3)))
   !
 enddo
 !
 !CLEAN
 !
 YAMBO_FREE(q_weight)
 call ELPH_alloc('FREE')
 !
end subroutine
