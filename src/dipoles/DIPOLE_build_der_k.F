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
 use pars,           ONLY:SP,cI,cZERO,pi,cONE
 use units,          ONLY:HA2EV
 use LIVE_t,         ONLY:live_timing
 use electrons,      ONLY:levels,n_sp_pol
 use DIPOLES,        ONLY:DIPOLE_t,DIP_S,DIP_iR,DIP_v,DIP_alloc,l_force_SndOrd
 use R_lattice,      ONLY:bz_samp,k_map
 use vec_operate,    ONLY:k_periodic_idx,degeneration_finder
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
 integer     :: ik,i_sp_pol,id_red,idx(3),idir,ic,ik_mem,ikbz,iv
 integer     :: istep,max_step,ikm1,ikp1
 complex(SP), allocatable :: U_l(:,:,:),U_r(:,:,:),Der_k_red(:,:,:),Der_k_cart(:,:,:),delta(:,:)
 !
 logical     :: deg_table_k(Dip%ib(2),Dip%ib(2)),deg_table_kp(Dip%ib(2),Dip%ib(2)),deg_table_km(Dip%ib(2),Dip%ib(2))
 integer     :: ib,ibp,i_grp,first_el_k(Dip%ib(2)),n_of_el_k(Dip%ib(2)),n_deg_grp_k
 real(SP)    :: Ev_m_Ec_KS,deg_thresh
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
 YAMBO_ALLOC(Der_k_cart,(Dip%ib(2),Dip%ib(2),3))
 YAMBO_ALLOC(Der_k_red,(Dip%ib(2),Dip%ib(2),2))
 YAMBO_ALLOC(delta,(Dip%ib(2),Dip%ib(2)))
 !
 delta=cZERO
 do ib=1,Dip%ib(2)
   delta(ib,ib)=cONE
 enddo
 !
 Dip%Energy_threshold=max(Dip%Energy_threshold,0.001/HA2EV)
 deg_thresh=Dip%Energy_threshold
 !
 do i_sp_pol=1,n_sp_pol
   !
   do ik=1,Xk%nibz 
     !
     if (.not.PAR_IND_DIPk_ibz%element_1D(ik)) cycle
     !
     ik_mem=PAR_DIPk_ibz_index(ik)
     !
     ikbz=Xk%k_table(ik,1)
     !
     call degeneration_finder(Xen%E(1:Dip%ib(2),ik,i_sp_pol),Dip%ib(2),&
     &                           first_el_k,n_of_el_k,n_deg_grp_k,deg_thresh,Include_single_values=.true.)
     !
     deg_table_k=.false.
     do i_grp=1,n_deg_grp_k
       do ib=first_el_k(i_grp),first_el_k(i_grp)+n_of_el_k(i_grp)-1
         do ibp=first_el_k(i_grp),first_el_k(i_grp)+n_of_el_k(i_grp)-1
           deg_table_k(ib,ibp)=.true.
           deg_table_k(ibp,ib)=.true.
         enddo
       enddo
     enddo
     !
     U_l       =cZERO
     U_r       =cZERO
     !
     Der_k_cart=cZERO
     Der_k_red =cZERO
     !
     do id_red=1,3  ! Loop on the three crystal directions in the BZ
       !
       if(k_map%max_kdir(id_red)==1) cycle ! Non-periodic direction, skip calculation
       !
       max_step=1 ! Only the first neighbor
       if(k_map%max_kdir(id_red)>=6.and..not.l_force_SndOrd) max_step=2
       !
       do istep=1,max_step
         ! 
         idx=k_map%k_map_inv(ikbz,:)
         idx(id_red)=idx(id_red)+istep
         idx=k_periodic_idx(idx,k_map)
         ikp1=k_map%k_map_dir(idx(1),idx(2),idx(3))
         !
         call degeneration_finder(Xen%E(1:Dip%ib(2),Xk%sstar(ikp1,1),i_sp_pol),Dip%ib(2),&
         &                           first_el_k,n_of_el_k,n_deg_grp_k,deg_thresh,Include_single_values=.true.)
         !
         deg_table_kp=.false.
         do i_grp=1,n_deg_grp_k
           do ib=first_el_k(i_grp),first_el_k(i_grp)+n_of_el_k(i_grp)-1
             do ibp=first_el_k(i_grp),first_el_k(i_grp)+n_of_el_k(i_grp)-1
               deg_table_kp(ib,ibp)=.true.
               deg_table_kp(ibp,ib)=.true.
             enddo
           enddo
         enddo
         ! 
         idx=k_map%k_map_inv(ikbz,:)
         idx(id_red)=idx(id_red)-istep
         idx=k_periodic_idx(idx,k_map)
         ikm1=k_map%k_map_dir(idx(1),idx(2),idx(3))
         !
         call degeneration_finder(Xen%E(1:Dip%ib(2),Xk%sstar(ikm1,1),i_sp_pol),Dip%ib(2),&
         &                           first_el_k,n_of_el_k,n_deg_grp_k,deg_thresh,Include_single_values=.true.)
         !
         deg_table_km=.false.
         do i_grp=1,n_deg_grp_k
           do ib=first_el_k(i_grp),first_el_k(i_grp)+n_of_el_k(i_grp)-1
             do ibp=first_el_k(i_grp),first_el_k(i_grp)+n_of_el_k(i_grp)-1
               deg_table_km(ib,ibp)=.true.
               deg_table_km(ibp,ib)=.true.
             enddo
           enddo
         enddo
         !
         U_r(:,:,istep)=          DIP_S(:,:,id_red+(istep-1)*3,ikbz,i_sp_pol)
         U_l(:,:,istep)=hermitian(DIP_S(:,:,id_red+(istep-1)*3,ikm1,i_sp_pol))
         !
         do ib=1,Dip%ib(2)
           do ibp=1,Dip%ib(2)
             if(deg_table_k(ibp,ib).or.deg_table_kp(ibp,ib)) cycle
             U_r(ib,ibp,istep)=cZERO
           enddo
         enddo
         !
         do ib=1,Dip%ib(2)
           do ibp=1,Dip%ib(2)
             if(deg_table_k(ibp,ib).or.deg_table_km(ibp,ib)) cycle
             U_l(ib,ibp,istep)=cZERO
           enddo
         enddo
         ! 
         call SERIAL_SVD(Dip%ib(2),U_r(:,:,istep),'uni',0)
         call SERIAL_SVD(Dip%ib(2),U_l(:,:,istep),'uni',0)
         !
         Der_k_red(:,:,istep)=+matmul(          DIP_S(:,:,id_red+(istep-1)*3,ikbz,i_sp_pol) ,hermitian(U_r(:,:,istep)) )  &
&                             -matmul(hermitian(DIP_S(:,:,id_red+(istep-1)*3,ikm1,i_sp_pol)),hermitian(U_l(:,:,istep)) )
         !
         Der_k_red(:,:,istep)=Der_k_red(:,:,istep)*cI*k_map%max_kdir(id_red)/(4._SP*pi*real(istep,SP))/2._SP
         !
       enddo
       !
       ! Transform in cartesian coordinates 
       !
       if(max_step==1) then
         !      
         ! First order formula D(dk)      
         !
         do idir=1,3
           Der_k_cart(:,:,idir)=Der_k_cart(:,:,idir)+transpose(Der_k_red(:,:,1))*a(id_red,idir)
         enddo
         !
       else
         !      
         ! Second order formula   ( 4*D(dk) - D(2*dk) ) / 3
         !
         do idir=1,3
           Der_k_cart(:,:,idir)=Der_k_cart(:,:,idir)+transpose(4._SP*Der_k_red(:,:,1)-Der_k_red(:,:,2))/3._SP*a(id_red,idir)
         enddo
         !
       endif
       !
     enddo  ! id_red loop on crystal directions 
     !
     !======
     ! <iR>
     !======
     do iv=Dip%ib(1),Dip%ib_lim(1)
       if(.not.PAR_IND_VAL_BANDS_DIP%element_1D(iv)) cycle
       do ic=max(iv,Dip%ib_lim(2)),Dip%ib(2)
         if(.not.PAR_IND_CON_BANDS_DIP%element_1D(ic)) cycle
         do idir=1,3
           DIP_iR(idir,ic,iv,ik_mem,i_sp_pol)=(Der_k_cart(ic,iv,idir)+conjg(Der_k_cart(iv,ic,idir)))/cI
         enddo
       enddo
     enddo
     !
     !=====
     ! <v>
     !=====
     if(.not.Dip%force_v_g_space) then
       do iv=Dip%ib(1),Dip%ib_lim(1) 
         if(.not.PAR_IND_VAL_BANDS_DIP%element_1D(iv)) cycle
         do ic=max(iv,Dip%ib_lim(2)),Dip%ib(2)
           if(.not.PAR_IND_CON_BANDS_DIP%element_1D(ic)) cycle
           Ev_m_Ec_KS   =Xen%E (iv,ik,i_sp_pol)-Xen%E (ic,ik,i_sp_pol)
           if (allocated(Xen%Eo)) Ev_m_Ec_KS   =Xen%Eo(iv,ik,i_sp_pol)-Xen%Eo(ic,ik,i_sp_pol)
           do idir=1,3
             DIP_v(idir,ic,iv,ik_mem,i_sp_pol)=DIP_iR(idir,ic,iv,ik_mem,i_sp_pol)*Ev_m_Ec_KS
           enddo
           !
         enddo
       enddo
       !
     endif
     !
     call live_timing(steps=1)
     !
   enddo ! ik iiBZ loop
   !
 enddo ! i_sp_pol
 !
 call live_timing()
 !
 YAMBO_FREE(U_l)
 YAMBO_FREE(U_r)
 YAMBO_FREE(Der_k_red)
 YAMBO_FREE(Der_k_cart)
 YAMBO_FREE(delta)
 !
 ! Build non-periodic directions if necessary
 ! 
 if(any(idir_not_done)) call DIPOLE_x_real_space(Xen,Xk,Dip,idir_not_done)
 !
 call timing("DIPOLE_buil_derk",OPR="stop")
 return
 !
end subroutine DIPOLE_build_der_k
