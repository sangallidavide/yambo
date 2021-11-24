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
subroutine DIPOLE_build_der_k(Xen,Xk,Dip)
 !
 ! Build the dipole as k-derivative using
 ! the overlaps  to contruct local rotation matrices 
 ! see Phys. Rev. B 76, 035213 (2007) and  Phys. Rev. Research 2, 013357 (2020)
 !
 use pars,           ONLY:SP,cI,cZERO,pi
 use LIVE_t,         ONLY:live_timing
 use electrons,      ONLY:levels,n_sp_pol
 use DIPOLES,        ONLY:DIPOLE_t,DIP_S,DIP_iR,DIP_v,DIP_alloc,l_force_SndOrd
 use R_lattice,      ONLY:bz_samp,k_map
 use vec_operate,    ONLY:k_periodic_idx
 use matrix_operate, ONLY:hermitian
 use parallel_m,     ONLY:PAR_IND_DIPk_ibz,PAR_DIPk_ibz_index,PAR_IND_DIPk_ibz_ID,&
&                         PAR_DIPk_nbz,PAR_IND_CON_BANDS_DIP, PAR_IND_VAL_BANDS_DIP, &
&                         PAR_IND_DIPk_bz,PAR_DIPk_bz_index,PAR_IND_DIPk_bz_ID,PAR_COM_DIPk_ibz_INDEX
 use parallel_int,   ONLY:PP_redux_wait
 use com,            ONLY:msg
 use D_lattice,      ONLY:a
 use parser_m,       ONLY:parser
 use timing_m,       ONLY:timing
 !
#include <memory.h>

 type(levels),           intent(in)     :: Xen
 type(bz_samp),          intent(inout)  :: Xk
 type(DIPOLE_t),         intent(inout)  :: Dip
 !
 ! Work Space
 !
 logical     :: idir_not_done(3)
 integer     :: ikbz,i_sp_pol,id
 integer     :: istep,max_step
 complex(SP), allocatable :: U_l(:,:,:),U_r(:,:,:)
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
 ! Notice that DIP_iR/DIP_v are defined only in the IBZ
 ! while the overlaps are in the BZ 
 !
 call live_timing("Build Rotation Matrices:",PAR_IND_DIPk_bz%n_of_elements(PAR_IND_DIPk_bz_ID+1)*n_sp_pol)
 !
 YAMBO_ALLOC(U_l,(Dip%ib(2),Dip%ib(2),2))
 YAMBO_ALLOC(U_r,(Dip%ib(2),Dip%ib(2),2))
 !
 do i_sp_pol=1,n_sp_pol
   !
   do ikbz=1,Xk%nbz
     !
     if(.not.PAR_IND_DIPk_bz%element_1D(ikbz)) cycle
     !
     do id=1,3  ! Loop on the three crystal directions in the BZ
       !
       if(k_map%max_kdir(id)==1) cycle ! Non-periodic direction, skip calculation
       !
       max_step=1 ! Only the first neighbor
       if(k_map%max_kdir(id)>=6.and..not.l_force_SndOrd) max_step=2
       !
       do istep=1,max_step
!         call SERIAL_SVD(Dip%ib(2),DIP_U(:,:,id+(istep-1)*3,ikbz,i_sp_pol),'uni',0)
       enddo
       !
     enddo
     !
     call live_timing(steps=1)
     !
   enddo ! BZ loop
   !
   call PP_redux_wait(DIP_S(:,:,:,:,i_sp_pol),COMM=PAR_COM_DIPk_ibz_INDEX%COMM)
   !
 enddo ! i_sp_pol
 !
 call live_timing()
 !
 YAMBO_FREE(U_l)
 YAMBO_FREE(U_r)
 !
 ! Build non-periodic directions if necessary
 ! 
 if(any(idir_not_done)) call DIPOLE_x_real_space(Xen,Xk,Dip,idir_not_done)
 !
 call timing("DIPOLE_buil_derk",OPR="stop")
 return
 !
end subroutine DIPOLE_build_der_k
