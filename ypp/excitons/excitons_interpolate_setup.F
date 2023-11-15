!
! License-Identifier: GPL
!
! Copyright (C) 2020 The Yambo Team
!
! Authors (see AUTHORS file for details): CA
!
subroutine excitons_interpolate_setup(k,Xk,en,Xen,q,ID_INTERP_EXC,report)
 !
 ! This subroutine read and interpolate exciton 
 !
 use pars,          ONLY:SP,schlen
 use units,         ONLY:HA2EV
 use R_lattice,     ONLY:bz_samp
 use electrons,     ONLY:levels
 use interpolate,   ONLY:INTERPOLATION_BZ,INTERP_obj,INTERPOLATION_coefficients
 use YPPm,          ONLY:BS_E,EXCITONS_user_indexes,EXCITONS_n_user_states,coo_out,K_transform,&
 &                       l_magnons,BS_R_left,BS_R_right,BS_R_left_magn,BS_R_right_magn
 use BS_solvers,    ONLY:BSS_n_eig
 use com,           ONLY:msg
 use stderr,        ONLY:intc
 !
#include<memory.h> 
 !
 type(bz_samp) ::Xk,k,q
 type(levels)  ::Xen,en
 logical, intent(in) :: report
 integer, intent(inout) :: ID_INTERP_EXC
 !
 ! Work space
 !
 integer :: iq,i_c,i_l,io_err,EXC_n_states_ref,BSS_n_eig_ref,n_eig
 real(SP):: BS_all_E(EXCITONS_n_user_states,q%nibz)
 real(SP), allocatable    :: values(:)
 real(SP) :: Res(BSS_n_eig),BS_E_degs(BSS_n_eig),tmp_q(3)
 character(schlen), allocatable :: headings(:)
 !
 EXC_n_states_ref=EXCITONS_n_user_states
 BSS_n_eig_ref=BSS_n_eig
 ! 
 do i_c=1,EXC_n_states_ref
   i_l=EXCITONS_user_indexes(i_c)
   BS_all_E(i_c,1)=BS_E(i_l)
 enddo
 !
 call section('*','Excitons Interpolation')
 call excitons_read(k,Xk,en,Xen,1,"clean",io_err)
 !
 ! Read all eigenvalues and eigenvectors
 !
 n_eig=BSS_n_eig
 do iq=2,q%nibz
   call msg('s','Reading excitons at @ Q-index #',iq)
   call excitons_read(k,Xk,en,Xen,iq,"check",io_err)
   call excitons_read(k,Xk,en,Xen,iq,"eigenvalues",io_err)
   write(*,*) BSS_n_eig,n_eig
   n_eig=min(n_eig,BSS_n_eig)
   if (     l_magnons)  Res(:n_eig)=real(BS_R_left_magn(1,:n_eig)*BS_R_right_magn(1,:n_eig),SP)
   if (.not.l_magnons)  Res(:n_eig)=real(BS_R_left(:n_eig)*BS_R_right(:n_eig),SP)
   write(*,*) "Residuals filled"
   call excitons_get_user_states(BS_E_degs(:n_eig),Res(:n_eig),n_eig,.true.)
   do i_c=1,EXC_n_states_ref
     i_l=EXCITONS_user_indexes(i_c)
     BS_all_E(i_c,iq)=BS_E(i_l)
   enddo
   call excitons_read(k,Xk,en,Xen,iq,"clean",io_err)
   !
 enddo
 !
 EXCITONS_n_user_states=EXC_n_states_ref
 BSS_n_eig=BSS_n_eig_ref
 !
 if(report.eqv..true.) then
   !
   ! Write in the report the exciton energies
   !
   YAMBO_ALLOC(values,(EXCITONS_n_user_states+4))
   allocate(headings(EXCITONS_n_user_states+4))
   !
   call msg('nr','Excitons Energies:')
   headings(1:3)=(/"q_x ("//trim(coo_out)//")","q_y ("//trim(coo_out)//")","q_z ("//trim(coo_out)//")"/)
   !
   do i_c=1,EXCITONS_n_user_states
      i_l=EXCITONS_user_indexes(i_c)
      headings(i_c+3)=' e'//trim(intc(i_l))//' [eV]'
   enddo
   call msg('r',"#",headings(1:EXCITONS_n_user_states+3))
   do iq=1,q%nibz
     tmp_q=q%pt(iq,:3)
     call K_transform(tmp_q,'iku')
     values(1:3)=tmp_q
     do i_c=1,EXCITONS_n_user_states
        i_l=EXCITONS_user_indexes(i_c)
        values(i_c+3)=BS_all_E(i_c,iq)*HA2EV
     enddo
     call msg('r',' ',values(1:EXCITONS_n_user_states+3))
   enddo
   !
   YAMBO_FREE(values)
   deallocate(headings)
   !
 endif
 !
 ! Fourier interpolation of exciton dispersion
 !
 call INTERPOLATION_BZ_setup(q)
 INTERP_obj(ID_INTERP_EXC)%what="excitons"
 call INTERPOLATION_coefficients(R1D=BS_all_E,k=q,NK=q%nibz,ID=ID_INTERP_EXC,ID_obj=ID_INTERP_EXC)
 !
end subroutine
