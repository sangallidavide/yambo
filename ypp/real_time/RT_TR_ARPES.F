!
!        Copyright (C) 2000-2021 the YAMBO team
!              http://www.yambo-code.org
!
! Authors (see AUTHORS file for details): DS
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
subroutine RT_transient_ARPES(Xen,Xk,Xq)
 !
 use pars,            ONLY:SP,cONE,cI,cZERO,schlen
 use units,           ONLY:HA2EV
 use R_lattice,       ONLY:bz_samp,bz_samp_duplicate_nosymm
 use BS,              ONLY:BS_K_dim,BS_bands
 use BS_solvers,      ONLY:BSS_eh_table,BS_mat !,BSS_n_eig
 use YPPm,            ONLY:BS_E,DIPs,EXCITONS_user_indexes !,EXCITONS_n_user_states
 use YPP_real_time,   ONLY:RT_conf,GreenF_T_and_W
 use parallel_m,      ONLY:PAR_IND_DIPk_ibz,PAR_IND_VAL_BANDS_DIP,PAR_IND_CON_BANDS_DIP
 use parallel_int,    ONLY:PP_redux_wait,PARALLEL_global_indexes
 use electrons,       ONLY:levels,spin
 use com,             ONLY:msg,of_open_close
 use interpolate,     ONLY:INTERPOLATE_use_symm
 !
#include<memory.h>
 !
 type(levels),  intent(in)    :: Xen
 type(bz_samp), intent(inout) :: Xk,Xq
 !
 type(bz_samp)         :: k_nosymm
 integer               :: i_l_ref,i_c,i_v,i_cvk,i_kibz,i_kbz,i_sp_pol,i_w,io_err
 real(SP)              :: F_w_real(RT_conf%W_n_steps,Xk%nbz,2),dw
 complex(SP)           :: F_tmp(RT_conf%W_n_steps)
 !
 character(schlen)     :: headings(3),file_name
 real(SP)              :: odata(3)
 !
 complex(SP),external      ::Lorentzian_FT
 !
 !
 ! Par Setup
 !============
 DIPs%ib = BS_bands
 call PARALLEL_global_indexes(Xen,Xk,Xq," ",RESET=.TRUE.)
 call PARALLEL_global_indexes(Xen,Xk,Xq,"DIPOLES",Dip=DIPs)
 !
 ! ... I/O
 call DIPOLE_IO(Xk,Xen,DIPs,'read',io_err,'DIP')
 if (io_err/=0) return
 !
 ! Frequency setup
 !============
 !call RT_1D_Fourier_setup(GreenF_T_and_W,RT_conf,1)
 !AMBO_FREE(GreenF_T_and_W%F_w)
 !AMBO_FREE(GreenF_T_and_W%F_t)
 !
 YAMBO_ALLOC(GreenF_T_and_W%W,(RT_conf%W_n_steps))
 GreenF_T_and_W%W=cZERO
 !
 dw=(RT_conf%X_W_range(2) - RT_conf%X_W_range(1))/(RT_conf%W_n_steps-1)
 do i_w = 1, RT_conf%W_n_steps
   GreenF_T_and_W%W(i_w)=cmplx((i_w-1)*dw+RT_conf%X_W_range(1),RT_conf%damp_factor)
 enddo
 !
 ! Set the exciton
 !============
 call RT_TRabs_rotate_exciton("DIP_iR",cONE,Xen,Xk)
 !
 i_l_ref=EXCITONS_user_indexes(1)
 !
 ! Set the kpts
 !============
 !
 call k_ibz2bz(Xk,'i',.true.)
 !
 call bz_samp_duplicate_nosymm(Xk,k_nosymm)
 ! 
 F_w_real=cZERO
 !
 do i_cvk = 1,BS_K_dim(1)
   !
   i_kbz    = BSS_eh_table(i_cvk,1)
   i_v      = BSS_eh_table(i_cvk,2)
   i_c      = BSS_eh_table(i_cvk,3)
   i_sp_pol = spin(BSS_eh_table(i_cvk,:))
   !
   i_kibz   = Xk%sstar(i_kbz,1)
   !
   if (.not. PAR_IND_CON_BANDS_DIP%element_1D(i_c)) cycle
   if (.not. PAR_IND_VAL_BANDS_DIP%element_1D(i_v)) cycle
   if (.not. PAR_IND_DIPk_ibz%element_1D(i_kibz)) cycle
   !
   do i_w=1,RT_conf%W_n_steps
     F_tmp(i_w)=Lorentzian_FT(GreenF_T_and_W%W(i_w),(Xen%E(i_v,i_kibz,i_sp_pol)+BS_E(i_l_ref)),"r")
   enddo
   !
   F_tmp=F_tmp*abs(BS_mat(i_cvk,i_l_ref))**2
   !
   F_w_real(:,i_kbz,1)=F_w_real(:,i_kbz,1)+ real(F_tmp(:),SP)
   F_w_real(:,i_kbz,2)=F_w_real(:,i_kbz,2)+aimag(F_tmp(:))
   !
 enddo
 !
 call PP_redux_wait(F_w_real)
 !
 ! Print G(w) integrated over k
 !
 file_name="G_w_integrated"
 call of_open_close(trim(file_name),'ot')
 headings=(/" E[eV]  ","Im[G(w)]","Re[G(w)]"/)
 call msg('o '//trim(file_name),"#",headings,INDENT=0,USE_TABS=.true.)
 !
 F_tmp=cZERO
 do i_kbz=1,k_nosymm%nibz
   F_tmp(:)=F_tmp(:)+cmplx(F_w_real(:,i_kbz,1),F_w_real(:,i_kbz,2),kind=SP) 
 enddo
 !
 do i_w=1,RT_conf%W_n_steps
   odata=(/real(GreenF_T_and_W%W(i_w),SP)*HA2EV,real(-cI*F_tmp(i_w),SP),real(F_tmp(i_w),SP)/)
   call msg('o '//trim(file_name),'',odata,INDENT=-2,USE_TABS=.true.)
 enddo
 !
 call of_open_close(trim(file_name))
 !
 INTERPOLATE_use_symm=.false.
 call RT_G_two_times_interp_and_plot(Xen,k_nosymm,k_nosymm%nibz,F_w_real)
 !
end subroutine RT_transient_ARPES
