!
!        Copyright (C) 2000-2019 the YAMBO team
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
subroutine RT_G_two_times_driver(en,kpt,qpt)
 !
 use pars,           ONLY:SP
 use R_lattice,      ONLY:bz_samp
 use electrons,      ONLY:levels
 use YPP_real_time,  ONLY:GreenF_T_and_W,RT_conf
 use RT_control,     ONLY:RTibz_to_RT_k
 use real_time,      ONLY:RT_nk
 !
#include<memory.h>
 !
 type(levels)  :: en
 type(bz_samp) :: kpt,qpt
 !
 ! Workspace
 !
 type(bz_samp) :: RT_k
 real(SP) :: G_w_t_real(RT_conf%W_n_steps,RT_nk,2)
 !real(SP),pointer, contiguous :: G_w_t_pt(:,:,:) => null()
 !
 call RT_G_two_times_build(en,kpt,qpt)
 !
 !!GF_size=(/2,size(GreenF_T_and_W%F_W,1),size(GreenF_T_and_W%F_W,2)/)
 !GF_size=(/2,RT_conf%W_n_steps,RT_nk/)
 !G_w_t_pt => linkC2R(GreenF_T_and_W%F_W,isize_loc)
 !
 G_w_t_real(:,:,1)=real(      GreenF_T_and_W%F_W, SP)
 G_w_t_real(:,:,2)=real(aimag(GreenF_T_and_W%F_W),SP)
 !
 YAMBO_FREE(GreenF_T_and_W%F_W)
 !
 call RTibz_to_RT_k(RT_k)
 !
 call k_ibz2bz(RT_k,'iku',.TRUE.)
 !
 call RT_G_two_times_interp_and_plot(en,RT_k,RT_nk,G_w_t_real)
 !
end subroutine RT_G_two_times_driver
