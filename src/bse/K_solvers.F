!
! License-Identifier: GPL
!
! Copyright (C) 2006 The Yambo Team
!
! Authors (see AUTHORS file for details): AM IMA
!
subroutine K_solvers(iq,Ken,Xk,q,X_static,W_bss)
 !
 use stderr,         ONLY:intc
 use frequency,      ONLY:w_samp
 use electrons,      ONLY:levels
 use R_lattice,      ONLY:bz_samp
 use X_m,            ONLY:X_t
 use BS_solvers,     ONLY:BSS_desc,BSS_free,BS_mat,run_Haydock,&
&                         run_inversion,run_Diago,run_Slepc,BSS_kernel_IO_on_the_fly
 use BS,             ONLY:BS_K_has_been_calculated_loaded,BS_ID,&
&                         l_BS_abs,l_BS_kerr,l_BS_magnons,l_BS_dichroism
 use parallel_m,     ONLY:master_cpu,CPU_structure
 use IO_int,         ONLY:IO_and_Messaging_switch
 use linear_algebra, ONLY:INV,DIAGO
 use PHOTOLUM,       ONLY:BSS_PL_free
 use MAGNONS,        ONLY:BSS_MAGN_free
 use DICHROISM,      ONLY:BSS_DICH_free
 !
#include<y_memory.h>
 !
 type(levels) ::Ken
 type(bz_samp)::Xk,q
 type(X_t)    ::X_static
 type(w_samp) ::W_bss
 integer      ::iq
 !
 logical :: l_partial_kernel_loaded
 !
 ! If the BS database was not loaded in memory, call K_restart
 !
 if (.not. BS_K_has_been_calculated_loaded ) call K_restart(iq,X_static,BS_ID(1),BS_ID(2),BS_ID(3),l_partial_kernel_loaded)
 !
 if (.not. BS_K_has_been_calculated_loaded ) call error("Something went wrong in K_restart, cannot proceed with solver")
 !
 call section('-','BSE solver(s) @q'//trim(intc(iq)))
 !
 ! The TDDFT way 
 ! DEBUG_AM
 ! if (index(BSS_mode,'t')/=0) then
 !   call tddft_BSK_disk_2_FXC(iq,Xk,W_bss,X_static) 
 !   l_bs_fxc=.true.
 ! endif
 ! DEBUG_AM
 ! Solvers
 !
 ! I/O. All read. Only master writes.
 ! 
 call IO_and_Messaging_switch("+io_out",CONDITION=master_cpu)
 !
 ! Diagonalization
 !
 if (run_Diago) then 
   !
   call PARALLEL_assign_LA_COMMs("Response_T_space",DIAGO,CPU_structure(4)%nCPU_lin_algebra_DIAGO)
   !
   call K_diago_driver(iq,W_bss,X_static)
   !
 endif
 !
 if (run_Slepc) call K_diago_driver(iq,W_bss,X_static)
 !
 ! Inversion
 !
 if (run_inversion) then 
   !
   call PARALLEL_assign_LA_COMMs("Response_T_space",INV,CPU_structure(4)%nCPU_lin_algebra_INV)
   !
   call K_inversion_driver(iq,W_bss,Ken,Xk,q)
   !
 endif
 !
 ! CLEAN: diago/invertion workspace
 !
 call BSS_free( )
 call BSS_MAGN_free( )
 call BSS_DICH_free( )
 call BSS_PL_free( )
 !
 ! CLEAN: BS_mat
 !
 YAMBO_FREE(BS_mat)
 !
 ! Haydock
 !
 if (run_Haydock) then
   !
   ! Remove Inv specific descriptors
   !
   if (run_inversion.and.master_cpu) BSS_desc%n=BSS_desc%n-3
   !
   ! Note: for developers (DS 2021/01/03)
   ! The jdos (l_BS_jdos) cannot be computed easely via the Haydock solver.
   ! it could be obtained as the average of N calls to K_Haydock
   ! with random intial vectors. See https://doi.org/10.1002/nla.2170
   !
   ! A possible alternative could be to use the algorithm with Cf
   ! from both the left and the right hand side to get from a single Haydock loop
   !
   if (l_BS_abs)     call K_Haydock(iq,W_bss,"optics",Xk)
   if (l_BS_kerr)    call K_Haydock(iq,W_bss,"kerr",Xk)
   if (l_BS_magnons) call K_Haydock(iq,W_bss,"magnons",Xk)
   if (l_BS_dichroism) call K_Haydock(iq,W_bss,"dichroism",Xk)
   !
 endif
 !
 ! CLEAN (final)
 !
 call BS_blks_free(1)
 !
end subroutine
