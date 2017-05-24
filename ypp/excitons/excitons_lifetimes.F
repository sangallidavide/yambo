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
subroutine excitons_sort_and_report(Xen,Xk,Xq,BS_E,BS_E_SOC_corr)
 !
 use pars,          ONLY:SP,schlen,pi
!cmau use units,         ONLY:HA2EV
!cmau
 use units,         ONLY:HA2EV,AU2SEC,SPEED_OF_LIGHT
!cmau
 use electrons,     ONLY:spin_occ,levels,spin,E_SOC_nbands
 use stderr,        ONLY:intc
 use BS,            ONLY:BS_H_dim
 use BS_solvers,    ONLY:BSS_description,BSS_n_descs
 use com,           ONLY:msg,of_open_close
 use R_lattice,     ONLY:d3k_factor,q0_def_norm,nkbz
!cmau
 use D_lattice,     ONLY:n_atomic_species,n_atoms_species,atom_pos,Z_species &
&                   ,DL_vol,non_periodic_directions
!cmau
 use X_m,           ONLY:eps_2_alpha,l_1D,l_2D,l_0D,l_3D
 use YPP,           ONLY:l_Rtimes
!cmau
 use vec_operate,   ONLY:sort
 !
 implicit none 
 !
 type(levels)       :: E
 complex(SP)        :: BS_R(BS_H_dim),BS_E(BS_H_dim)
 real(SP),pointer   :: BS_E_SOC_corr(:,:)
 real(SP),optional  :: S_z(BS_H_dim)
 real(SP),optional  :: S_sq(BS_H_dim)
 !
 ! Work Space
 !
 integer            :: j1,j2,ia,is,i_mode,ic,n_elements
 integer            :: S_indx(BS_H_dim)
 real(SP)           :: Residuals(BS_H_dim),v2sort(BS_H_dim),Rmax,rv(10)
 !
 real(SP)           :: phi(Phi_grid_size)
 !
 ! Excitonc amplitude
 !
 character(3), parameter:: R_normalize="yes"
 !
 write_spin    = present(S_sq)
 write_widths  = any(abs(aimag(BS_E))>1.E-5)
 write_soc_corr= (E_SOC_nbands/=0)
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
 YAMBO_ALLOC(q_norm,(1))
 YAMBO_ALLOC(bare_qpg,(1,1))
 !
 call PARALLEL_global_indexes(Xen,Xk,Xq," ",RESET=.TRUE.)
 call PARALLEL_global_indexes(Xen,Xk,Xq,"Response_G_space_Zero_Momentum",X=Xbsk)
 !
 io_DIP=.FALSE.
 call DIPOLE_driver(Xen,Xk,Xbsk,Xbsk%q0)
 !
 iR_EXC=cZERO
 !
 do i_lambda=1,n_lambda
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
   ! First compute the excionic dipoles
   !==========================================
   !
   P_EXC_group(i_l_grp)=0._SP
   gamma_thermal(:)=0._SP
   gamma_exp =0._SP
   gamma0_grp=0._SP
   !
   do i_l=BS_E_degs(i_l_grp),BS_E_degs(i_l_grp)+n_lambda_deg-1
     !
     do neh = 1,BS_H_dim !(size BSE matrix, cvk)
       !
       i_kbz   = BSS_eh_table(neh,1)
       i_kibz  = Xk%sstar(i_kbz,1)
       i_v     = BSS_eh_table(neh,2)
       i_c     = BSS_eh_table(neh,3)
       i_spin  = spin(BSS_eh_table(neh,:))
       !
       f_cv=      ( Xen%f(i_v,i_kibz,i_spin)- Xen%f(i_c,i_kibz,i_spin))/spin_occ
       !
       iR_EXC(:,i_l)=iR_EXC(i_l)-DIP_iR(:,ic,iv,ik_mem,i_sp_pol)*BS_mat(neh,i_l)*sqrt(f_cv)*q0_def_norm
       !
     enddo
     !
     ! Angle of the excitonic dipoles against the x direction
     !
     ! Warning define N degenerate angles phi(n)
     phi0(i_l)=acos(iR_EXC(1,i_l))
     !
     P_EXC_square(i_l) =( iR_EXC(1,i_l)*conjg(iR_EXC(1,i_l))**2+ &
&                         iR_EXC(2,i_l)*conjg(iR_EXC(2,i_l))**2+ &
&                         iR_EXC(3,i_l)*conjg(iR_EXC(3,i_l))**2 ) * real(BS_E(i_l))**2
     !
     ! Second: compute the lifetime at q=0
     !==========================================
     if (l_2D) gamma0(i_l) = (eps_2_alpha/DL_vol)*P_EXC_square(i_l)*8.*pi/(real(nkbz)*SPEED_OF_LIGHT)
     if (l_1D) gamma0(i_l) = (eps_2_alpha/DL_vol)*P_EXC_square(i_l)*2.*pi*real(nkbz)*SPEED_OF_LIGHT**2
     if (l_3D) gamma0(i_l) =                      P_EXC_square(i_l)*4.*0.333/(real(nkbz)*SPEED_OF_LIGHT**3)
     if (l_0D) gamma0(i_l) =                      P_EXC_square(i_l)*4.*0.333/(real(nkbz)*SPEED_OF_LIGHT**3)
     !
     !
     ! Third: define temperature and angle dependent lifetime
     !=======================================================
     !
     alpha=BS_E(i_l)**2/(2._SP*EXC_mass*EXC_temp*SPEED_OF_LIGHT**2)
     !
     gamma0_grp=gamma0_grp+gamma0(i_l)
     !
     gamma_thermal(:)=gamma_thermal(:)+&
&       gamma0(i_l)*alpha*(2._SP*sin(phi(:)+phi0(i_l))**2+2._SP/3._SP*cos(phi(:)+phi0(i_l))**2)
     !
     gamma_exp=gamma_exp+gamma0(i_l)*alpha*(4._SP/3._SP)
     !
   enddo
   !
   call msg('o exc_life',"# ","gamma zero = "//trim(real2ch(gamma0_grp)),INDENT=0)
   call msg('o exc_life',"# ","gamma exp  = "//trim(real2ch(gamma_exp)),INDENT=0)

   do i_grid=1,Phi_grid_size
     call msg('o exc_life',' ',(/phi(i_grid),gamma_thermal(i_grid)/),INDENT=0,USE_TABS=.true.)
   enddo
   !
 enddo
 !
end subroutine
