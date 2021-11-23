!
!        Copyright (C) 2000-2021 the YAMBO team
!              http://www.yambo-code.org
!
! Authors (see AUTHORS file for details): MG, CA, DS
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
subroutine DIPOLE_build_derk(Xen,Xk,Dip)
 !
 ! Build the dipole as k-derivative using
 ! the overlaps  to contruct local rotation matrices 
 ! see Phys. Rev. B 76, 035213 (2007) and  Phys. Rev. Research 2, 013357 (2020)
 !
 use pars,           ONLY:SP,cI,cZERO,pi
 use LIVE_t,         ONLY:live_timing
 use electrons,      ONLY:levels,n_sp_pol
 use DIPOLES,        ONLY:DIPOLE_t,DIP_S,DIP_iR,DIP_v
 use R_lattice,      ONLY:bz_samp,k_map
 use vec_operate,    ONLY:k_periodic_idx
 use matrix_operate, ONLY:hermitian
 use parallel_m,     ONLY:PAR_IND_DIPk_ibz,PAR_DIPk_ibz_index,PAR_IND_DIPk_ibz_ID,&
&                         PAR_DIPk_nbz,PAR_IND_CON_BANDS_DIP, PAR_IND_VAL_BANDS_DIP
 use com,            ONLY:msg
 use D_lattice,      ONLY:a
 use parser_m,       ONLY:parser
 use timing_m,       ONLY:timing

 !
 implicit none
 !
 type(levels),           intent(in)     :: Xen
 type(bz_samp),          intent(inout)  :: Xk
 type(DIPOLE_t),         intent(inout)  :: Dip
 !
 ! Work Space
 !
 logical     :: idir_not_done(3)
 logical     :: l_eval_polarization,l_force_SndOrd
 integer     :: ik,iv,ic,ic_min,iv_max,i_sp_pol,ikm1,id_red,idir,idx(3),ikbz,ik_mem
 integer     :: istep,max_step
 !
 call timing("DIPOLE_buil_derk",OPR="start")
 !
 if(.not.Dip%force_v_g_space) Dip%Vnl_included=.true.
 !
 idir_not_done(:)= (k_map%max_kdir(:)==1)
 !
 if(any(idir_not_done)) call warning(' Der-k Dipoles not in 3d. Assuming non periodic system.')
 !
 if(any(k_map%max_kdir(:)>=6.and..not.l_force_SndOrd)) &
 &   call msg('rs','Using forth order approximation for der-k dipoles')
 !
 call live_timing('Der-k Dipoles',PAR_IND_DIPk_ibz%n_of_elements(PAR_IND_DIPk_ibz_ID+1)*n_sp_pol)
 !
 ! Notice that DIP_iR/DIP_v are defined only in the IBZ
 ! while the overlaps are in the BZ 
 !
 call live_timing()
 !
 if(any(idir_not_done)) call DIPOLE_x_real_space(Xen,Xk,Dip,idir_not_done)
 !
 call timing("DIPOLE_buil_derk",OPR="stop")
 return
 !
end subroutine DIPOLE_build_derk
