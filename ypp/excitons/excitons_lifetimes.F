!
!        Copyright (C) 2000-2017 the YAMBO team
!              http://www.yambo-code.org
!
! Authors (see AUTHORS file for details): AM DV
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
subroutine excitons_lifetimes(Xbsk,Xen,Xk,Xq,BS_H_dim,BS_R,BS_E,BS_E_degs,lambda_s,n_lambda)
 !
 use pars,          ONLY:SP,schlen,pi,cZERO
 use units,         ONLY:HA2EV,AU2SEC,SPEED_OF_LIGHT
 use electrons,     ONLY:spin_occ,levels,spin
 use stderr,        ONLY:intc,real2ch
 use BS_solvers,    ONLY:BSS_eh_table,BS_mat
 use BS,            ONLY:BS_bands
 use com,           ONLY:msg,of_open_close
 use R_lattice,     ONLY:q0_def_norm,nkbz,bz_samp,q_norm,bare_qpg
 use D_lattice,     ONLY:DL_vol,non_periodic_directions,l_1D,l_2D,l_0D,l_3D
 use X_m,           ONLY:X_t,DIP_iR,eps_2_alpha
 use YPP,           ONLY:EXC_temp,EXC_mass,EXC_Phi_Grid
 use interfaces,    ONLY:PARALLEL_global_indexes
 use vec_operate,   ONLY:v_norm
 use IO_m,          ONLY:io_DIP
 !
 ! implicit none 
 !
#include<memory.h>
 !
 integer            :: BS_H_dim
 type(levels)       :: Xen
 type(X_t)          :: Xbsk
 type(bz_samp)      :: Xk,Xq
 !
 complex(SP)        :: BS_E(BS_H_dim), BS_R(BS_H_dim) 
 !
 character(schlen)  :: ch_dummy
 !
 ! Work Space
 !
 integer            :: i_lambda,i_grid
 integer            :: neh,i_l,i_l_grp,n_lambda_deg,i_v,i_c,i_kbz,i_kibz,i_spin 
 integer            :: Phi_grid_size
 !
 complex(SP)        :: Ir_EXC(3),rotated_iR(3)
 !
 real(SP)           :: P_EXC_square
 real(SP)           :: gamma0(BS_H_dim)
 real(SP)           :: gamma_exp, gamma0_grp,f_cv,alpha,tau_exp,tau0_grp
 integer            :: BS_E_degs(BS_H_dim),lambda_s(n_lambda),n_lambda
 !
 !
 real(SP)           :: phi(EXC_Phi_Grid),exc_energy
 real(SP)           :: phi0(BS_H_dim),vdummy(3),tau0
 real(SP)           :: gamma_thermal(EXC_Phi_Grid),tau_thermal(EXC_Phi_Grid)
 !
 ! Excitonc times 
 !
 character(3), parameter:: R_normalize="yes"
 !
 !
 !
 !
  Phi_Grid_size = EXC_Phi_Grid
 !
 call section('=','Computing Excitons radiative lifetimes')
 !
 if (trim(non_periodic_directions)/="none")  call check_periodic_directions()
 !
 do i_grid=1,Phi_grid_size
   phi(i_grid)=2._SP*pi*real(i_grid-1,SP)/real(Phi_grid_size-1,SP)
 enddo
 !
 Xbsk%ib=BS_bands
 ! 
 YAMBO_ALLOC(q_norm,(1))
 YAMBO_ALLOC(bare_qpg,(1,1))
 !
 call PARALLEL_global_indexes(Xen,Xk,Xq," ",RESET=.TRUE.)
 call PARALLEL_global_indexes(Xen,Xk,Xq,"Response_G_space_Zero_Momentum",X=Xbsk)
 !
 io_DIP=.FALSE.
 call DIPOLE_driver(Xen,Xk,Xbsk,Xbsk%q0)
 !
 YAMBO_FREE(q_norm)
 YAMBO_FREE(bare_qpg)
 !
 do i_lambda=1,n_lambda
!
   iR_EXC = cZERO
   gamma_thermal(:)=0._SP
   gamma_exp =0._SP
   gamma0_grp=0._SP
   tau_thermal(:)=0._SP
   tau_exp =0._SP
   tau0_grp=0._SP
   !
   i_l_grp=lambda_s(i_lambda)
   !
   n_lambda_deg=count(BS_E_degs==BS_E_degs(i_l_grp))
   !
   if (n_lambda_deg>1) call msg('s',':: State '//trim(intc(i_l_grp))//' Merged with states '//&
&                                     trim(intc(BS_E_degs(i_l_grp)))//' -> '//&
&                                     trim(intc(BS_E_degs(i_l_grp)+n_lambda_deg-1)))
   !
   ch_dummy='exc_lifetimes_at_'//trim(intc(i_l_grp))
   call of_open_close(ch_dummy,'ot')
   !
   ! First compute the excitonic dipoles
   !==========================================
   !
   !
   do i_l=BS_E_degs(i_l_grp),BS_E_degs(i_l_grp)+n_lambda_deg-1
     do neh = 1,BS_H_dim !(size BSE matrix, cvk)
       !
       i_kbz   = BSS_eh_table(neh,1)
       i_kibz  = Xk%sstar(i_kbz,1)
       i_v     = BSS_eh_table(neh,2)
       i_c     = BSS_eh_table(neh,3)
       i_spin  = spin(BSS_eh_table(neh,:))
       !
       f_cv=      (Xen%f(i_v,i_kibz,i_spin)- Xen%f(i_c,i_kibz,i_spin))/spin_occ
       !
       call GET_the_rotated_iR()
       !
       iR_EXC(:)=iR_EXC(:) + rotated_iR(:)*BS_mat(neh,i_l)*sqrt(f_cv)
       !
     enddo
      write(*,*) 'ir_EXC*q0_def_norm,BS_R', iR_EXC(1)*q0_def_norm, BS_R(i_l) 
     !
     ! Angle of the excitonic dipoles against the x direction
     !
     ! Warning define N degenerate angles phi(n)
     vdummy(:) = real(iR_exc(:)) / v_norm(real(iR_exc(:)))
     phi0(i_l)=acos(vdummy(1))
     !
     P_EXC_square =(real(iR_EXC(1)*conjg(iR_EXC(1)) + &
&                 iR_EXC(2)*conjg(iR_EXC(2)) + &
&                 iR_EXC(3)*conjg(iR_EXC(3)) )*real(BS_E(i_l))**2)
     !
     ! Second: compute the lifetime at q=0
     !==========================================
     if (l_2D) gamma0(i_l) = (eps_2_alpha/DL_vol)*P_EXC_square*4.*pi &
&                             /(real(BS_E(i_l))*SPEED_OF_LIGHT*real(nkbz))
     if (l_1D) gamma0(i_l) = (eps_2_alpha/DL_vol)*P_EXC_square*2.*pi*real(nkbz)*SPEED_OF_LIGHT**2
     if (l_3D) gamma0(i_l) =                      P_EXC_square*4.*0.333/(real(nkbz)*SPEED_OF_LIGHT**3)
     if (l_0D) gamma0(i_l) =                      P_EXC_square*4.*0.333/(real(nkbz)*SPEED_OF_LIGHT**3)
     !
     !
     ! Third: define temperature and angle dependent lifetime
     !=======================================================
     !
     alpha=real(BS_E(i_l))**2/(2._SP*EXC_mass*EXC_temp*SPEED_OF_LIGHT**2)
     !
     gamma0_grp=gamma0_grp+gamma0(i_l)
     !
     !
     gamma_thermal(:)=gamma_thermal(:)+&
&       gamma0(i_l)*alpha*(2._SP*sin(phi(:)+phi0(i_l))**2+2._SP/3._SP*cos(phi(:)+phi0(i_l))**2)
     !
     gamma_exp=gamma_exp+gamma0(i_l)*alpha*(4._SP/3._SP)
     !
   exc_energy = real(BS_E(i_l)*HA2EV)
   write(10,*) i_l, exc_energy, AU2SEC/gamma0(i_l), DL_vol,eps_2_alpha
   enddo
   tau0_grp= AU2SEC/gamma0_grp
   tau_exp = AU2SEC/gamma_exp
   !
   call msg('o exc_life',"# ","exc. E = "//trim(real2ch(exc_energy)),INDENT=0)
   call msg('o exc_life',"# ","gamma zero = "//trim(real2ch(gamma0_grp)),INDENT=0)
   call msg('o exc_life',"# ","gamma exp  = "//trim(real2ch(gamma_exp)),INDENT=0)
   !
   call msg('o exc_life',"# ","tau zero = "//trim(real2ch(tau0_grp)),INDENT=0)
   call msg('o exc_life',"# ","tau exp  = "//trim(real2ch(tau_exp)),INDENT=0)
   do i_grid=1,Phi_grid_size
    tau_thermal(i_grid)= AU2SEC/gamma_thermal(i_grid)
     call msg('o exc_life',' ',(/phi(i_grid),gamma_thermal(i_grid),tau_thermal(i_grid)/),INDENT=0,USE_TABS=.true.)
   enddo
   !
 enddo
 !
 contains
   !
   subroutine GET_the_rotated_iR( )
     !
     use D_lattice,      ONLY:i_time_rev,dl_sop,sop_inv,nsym
     use parallel_m,     ONLY:PAR_DIPk_ibz_index
     !
     integer :: ik,is,ik_mem
     logical :: t_rev
     !
     ik=Xk%sstar(i_kbz,1)
     is=Xk%sstar(i_kbz,2)
     !
     if ( is<= nsym/(i_time_rev+1) ) t_rev=.false.
     if ( is > nsym/(i_time_rev+1) ) t_rev=.true.
     !
     is = sop_inv( is )
     !      
     ik_mem    = PAR_DIPk_ibz_index(ik)
     rotated_iR = matmul( dl_sop(:,:,is) , DIP_iR(:,i_c,i_v,ik_mem,i_spin) )
     !
     if (t_rev) rotated_iR = conjg( rotated_iR )
     !
   end subroutine
 !
end subroutine
