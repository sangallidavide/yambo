!
!        Copyright (C) 2000-2021 the YAMBO team
!              http://www.yambo-code.org
!
! Authors (see AUTHORS file for details): DS AM
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
subroutine RT_TRabs_residuals(what,N_trans,N_dirs,DIP_dir,Xen,Xk,nb,l_EQ)
 !
 use pars,            ONLY:SP,pi,cI
 use drivers,         ONLY:l_rt_carriers_in_use
 use R_lattice,       ONLY:bz_samp
 use electrons,       ONLY:levels,n_sp_pol,spin,spin_occ
 use X_m,             ONLY:global_gauge
 use DIPOLES,         ONLY:DIP_rotated
 use YPP_real_time,   ONLY:TRabs_RES_left,TRabs_RES_right,TRabs_E,TRabs_Eeh_treshold,&
&                          TRabs_Exc_Ref
 use BS,              ONLY:BS_K_dim
 use BS_solvers,      ONLY:BSS_eh_table,BS_mat
 use YPPm,            ONLY:BS_E
 use real_time,       ONLY:RT_bands
 use parallel_m,      ONLY:PAR_IND_DIPk_ibz,PAR_IND_VAL_BANDS_DIP,PAR_IND_CON_BANDS_DIP
 !
#include<memory.h>
 !
 character(*)      :: what
 integer           :: N_trans,nb(2),N_dirs
 type(levels)      :: Xen
 type(bz_samp)     :: Xk
 real(SP)          :: DIP_dir(3,3)
 logical           :: l_EQ
 !
 ! Work space
 !
 integer           :: i_cv,i_v,i_c,i_kbz,i_kibz,i_spin,i_l,i_l_ref,bands(2,2),i_n,i_np,i_dir, &
&                     i_cv_p,i_v_p,i_c_p,i_kbz_p,i_spin_p,i_trans
 real(SP)          :: E_eh
 complex(SP)       :: dipole_rotated(3),f_eh(2),DIP_dot
 !
 ! BSE residuals
 ! -------------
 if (index(what,"BSE")>0) then
   !
   do i_l=1,BS_K_dim(1)
     !
     TRabs_E(i_l)   = real(BS_E(i_l),SP)
     !
     do i_cv = 1,BS_K_dim(1)
       !
       i_kbz   = BSS_eh_table(i_cv,1)
       i_v     = BSS_eh_table(i_cv,2)
       i_c     = BSS_eh_table(i_cv,3)
       i_spin  = spin(BSS_eh_table(i_cv,:))
       !
       i_kibz  = Xk%sstar(i_kbz,1)
       !
       if (i_v<RT_bands(1) .or. i_c>RT_bands(2) ) cycle
       !
       if (.not. PAR_IND_CON_BANDS_DIP%element_1D(i_c)) cycle
       if (.not. PAR_IND_VAL_BANDS_DIP%element_1D(i_v)) cycle
       if (.not. PAR_IND_DIPk_ibz%element_1D(i_kibz)) cycle
       !
       f_eh=gimme_f_and_delta_f(i_v,i_c,i_kibz,i_spin)
       dipole_rotated=conjg(-cI*DIP_rotated(i_c,i_v,i_kbz,i_spin,"DIP_iR",Xk))
       !
       do i_dir=1,N_dirs
         !
         DIP_dot=dot_product(DIP_dir(i_dir,:),dipole_rotated)
         !
         TRabs_RES_right(i_l,i_dir,1)= TRabs_RES_right(i_l,i_dir,1)+&
&                                      conjg(DIP_dot)*BS_mat(i_cv,i_l)*sqrt(f_eh(1))
         TRabs_RES_left(i_l,i_dir,1) = TRabs_RES_left(i_l,i_dir,1)+&
&                                      DIP_dot*conjg(BS_mat(i_cv,i_l))*sqrt(f_eh(1))
         !
         if (l_rt_carriers_in_use) then
           TRabs_RES_right(i_l,i_dir,2)= TRabs_RES_right(i_l,i_dir,2)+&
&                                        conjg(DIP_dot)*BS_mat(i_cv,i_l)*f_eh(2)/sqrt(f_eh(1))/sqrt(2._SP)
           TRabs_RES_left(i_l,i_dir,2) = TRabs_RES_left(i_l,i_dir,2)+&
&                                        DIP_dot*conjg(BS_mat(i_cv,i_l))*f_eh(2)/sqrt(f_eh(1))/sqrt(2._SP)
         endif
         !
       enddo
       !
     enddo
     !
   enddo
   !
   N_trans=BS_K_dim(1)
   !
   return
   !
 endif
 !
 ! EXC-EXC residuals
 ! -------------
 if (index(what,"EXC")>0) then
   !
   i_l_ref=TRabs_Exc_Ref
   !
   do i_l=1,BS_K_dim(1)
     !
     i_trans=i_l+N_trans
     TRabs_E(i_trans) = real(BS_E(i_l)-BS_E(i_l_ref),SP)
     !
     if (i_l==i_l_ref) cycle
     !
     do i_cv = 1,BS_K_dim(1)
       !
       i_kbz   = BSS_eh_table(i_cv,1)
       i_v     = BSS_eh_table(i_cv,2)
       i_c     = BSS_eh_table(i_cv,3)
       i_spin  = spin(BSS_eh_table(i_cv,:))
       !
       i_kibz  = Xk%sstar(i_kbz,1)
       !
       if (.not. PAR_IND_CON_BANDS_DIP%element_1D(i_c)) cycle
       if (.not. PAR_IND_VAL_BANDS_DIP%element_1D(i_v)) cycle
       if (.not. PAR_IND_DIPk_ibz%element_1D(i_kibz)) cycle
       !
       ! c--> c' term
       !
       do i_cv_p = 1,BS_K_dim(1)
         !
         i_kbz_p   = BSS_eh_table(i_cv_p,1)
         i_v_p     = BSS_eh_table(i_cv_p,2)
         i_c_p     = BSS_eh_table(i_cv_p,3)
         i_spin_p  = spin(BSS_eh_table(i_cv_p,:))
         !
         i_kibz  = Xk%sstar(i_kbz,1)
         !
         if(i_spin_p/=i_spin ) cycle
         if(i_kbz_p /=i_kbz_p) cycle
         if (i_v    /=i_v_p  ) cycle
         !
         E_eh=Xen%E(i_c,i_kibz,i_spin)-Xen%E(i_c_p,i_kibz,i_spin)
         if (abs(E_eh) < TRabs_Eeh_treshold.and. trim(global_gauge)=="length") cycle
         !
         dipole_rotated=conjg(-cI*DIP_rotated(i_c_p,i_c,i_kbz,i_spin,"DIP_iR",Xk))
         !
         do i_dir=1,N_dirs
           !
           DIP_dot=dot_product(DIP_dir(i_dir,:),dipole_rotated)
           !
           TRabs_RES_right(i_trans,i_dir,1)= TRabs_RES_right(i_l,i_dir,1)+&
&                                            conjg(BS_mat(i_cv,i_l_ref))*conjg(DIP_dot)*BS_mat(i_cv_p,i_l)
           TRabs_RES_left(i_trans,i_dir,1) = TRabs_RES_left(i_l,i_dir,1)+&
&                                            BS_mat(i_cv,i_l_ref)*DIP_dot*conjg(BS_mat(i_cv_p,i_l))
           !
         enddo
         !
       enddo
       !
       ! v --> v' term
       !
       do i_cv_p = 1,BS_K_dim(1)
         !
         i_kbz_p   = BSS_eh_table(i_cv_p,1)
         i_v_p     = BSS_eh_table(i_cv_p,2)
         i_c_p     = BSS_eh_table(i_cv_p,3)
         i_spin_p  = spin(BSS_eh_table(i_cv_p,:))
         !
         i_kibz  = Xk%sstar(i_kbz,1)
         !
         if(i_spin_p/=i_spin ) cycle
         if(i_kbz_p /=i_kbz_p) cycle
         if (i_c    /=i_c_p  ) cycle
         !
         E_eh=Xen%E(i_v,i_kibz,i_spin)-Xen%E(i_v_p,i_kibz,i_spin)
         if (abs(E_eh) < TRabs_Eeh_treshold.and. trim(global_gauge)=="length") cycle
         !
         dipole_rotated=conjg(-cI*DIP_rotated(i_v_p,i_v,i_kbz,i_spin,"DIP_iR",Xk))
         !
         do i_dir=1,N_dirs
           !
           DIP_dot=dot_product(DIP_dir(i_dir,:),dipole_rotated)
           !
           TRabs_RES_right(i_trans,i_dir,1)= TRabs_RES_right(i_l,i_dir,1)-&
&                                            conjg(BS_mat(i_cv,i_l_ref))*conjg(DIP_dot)*BS_mat(i_cv_p,i_l)
           TRabs_RES_left(i_trans,i_dir,1) = TRabs_RES_left(i_l,i_dir,1)-&
&                                            BS_mat(i_cv,i_l_ref)*DIP_dot*conjg(BS_mat(i_cv_p,i_l))
           !
         enddo
         !
       enddo
       !
     enddo
     !
   enddo
   !
   N_trans=N_trans+BS_K_dim(1)
   !
 endif
 !
 ! IP residuals
 ! -------------
 if (index(what,"IP")>0) then
   !
   if (index(what,"cv")>0) then
     bands(1,:)=(/nb(1)+1,RT_bands(2)/)
     bands(2,:)=(/RT_bands(1),nb(2)/)
   else if (index(what,"vv")>0) then
     bands(1,:)=(/RT_bands(1),nb(1)/)
     bands(2,:)=bands(1,:)
   else if (index(what,"cc")>0) then
     bands(1,:)=(/nb(2)+1,RT_bands(2)/)
     bands(2,:)=bands(1,:)
   endif
   !
   do i_kbz=1,Xk%nbz
     i_kibz=Xk%sstar(i_kbz,1)
     do i_spin=1,n_sp_pol
       do i_n=bands(1,1),bands(1,2)
         do i_np=bands(2,1),bands(2,2)
           !
           E_eh=Xen%E(i_n,i_kibz,i_spin)-Xen%E(i_np,i_kibz,i_spin)
           f_eh=gimme_f_and_delta_f(i_np,i_n,i_kibz,i_spin)
           !
           ! Select just "resonant" (in the causal sense) transitions
           if (real(f_eh(1)+f_eh(2))<0._SP) cycle
           !
           if (l_EQ) f_eh(2)=f_eh(2)+f_eh(1)
           !
           if (abs(E_eh) < TRabs_Eeh_treshold.and. trim(global_gauge)=="length") cycle
           !
           N_trans=N_trans+1
           !
           if (.not. PAR_IND_DIPk_ibz%element_1D(i_kibz)) cycle
           if (.not. PAR_IND_CON_BANDS_DIP%element_1D(i_np)) cycle
           if (.not. PAR_IND_VAL_BANDS_DIP%element_1D(i_n)) cycle
           !
           TRabs_E(N_trans)=E_eh
           dipole_rotated=conjg(-cI*DIP_rotated(i_np,i_n,i_kbz,i_spin,"DIP_iR",Xk))
           !
           do i_dir=1,N_dirs
             TRabs_RES_right(N_trans,i_dir,1)= conjg(dot_product(DIP_dir(i_dir,:),dipole_rotated))*sqrt(f_eh(2))
             TRabs_RES_left(N_trans,i_dir,1) =       dot_product(DIP_dir(i_dir,:),dipole_rotated) *sqrt(f_eh(2))
           enddo
           !
         enddo
       enddo
     enddo
   enddo
   !
   return
   !
 endif
 !
 contains
   !
   function gimme_f_and_delta_f(ib,ibp,ik,i_sp)
     complex(SP) :: gimme_f_and_delta_f(2)
     real(SP)    :: Delta_f,Delta_df
     integer     :: ib,ibp,i_sp,ik
     if(l_rt_carriers_in_use) then
       Delta_f=(Xen%fo(ib,ik,i_sp)-Xen%fo(ibp,ik,i_sp))/spin_occ
       !Delta_df=(Xen%df(ib,ik,i_sp)-Xen%df(ibp,ik,i_sp))/spin_occ
       Delta_df=(Xen%f(ib,ik,i_sp)-Xen%f(ibp,ik,i_sp))/spin_occ-Delta_f
     else
       Delta_f=(Xen%f(ib,ik,i_sp)-Xen%f(ibp,ik,i_sp))/spin_occ
       Delta_df=0._SP
     endif
     gimme_f_and_delta_f(1)=Delta_f
     gimme_f_and_delta_f(2)=Delta_df
   end function 
   !
end subroutine RT_TRabs_residuals
