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
 use parser_m,        ONLY:parser
 use units,           ONLY:HA2EV
 use pars,            ONLY:SP,cONE,cI,cZERO,schlen
 use R_lattice,       ONLY:bz_samp,bz_samp_duplicate_nosymm
 use D_lattice,       ONLY:lattice,nsym
 use BS,              ONLY:BS_K_dim,BS_bands
 use BS_solvers,      ONLY:BSS_eh_table,BS_mat               !,BSS_n_eig
 use YPPm,            ONLY:BS_E,DIPs,EXCITONS_user_indexes   !,EXCITONS_n_user_states
 use YPPm,            ONLY:BANDS_path,BANDS_steps,CIRCUIT_k  !,USER_k
 use YPP_interfaces,  ONLY:k_special
 use YPP_real_time,   ONLY:RT_conf,GreenF_T_and_W
 use parallel_m,      ONLY:PAR_IND_DIPk_ibz,PAR_IND_VAL_BANDS_DIP,PAR_IND_CON_BANDS_DIP,myid
 use parallel_int,    ONLY:PP_redux_wait,PARALLEL_global_indexes
 use electrons,       ONLY:levels,spin,n_sp_pol
 use com,             ONLY:msg,of_open_close
 use interpolate,     ONLY:INTERPOLATE_use_symm,INTERP_obj, &
 &                         INTERPOLATION_driver_seed,INTERPOLATION_driver_do
 !
#include<memory.h>
 !
 type(levels),  intent(in)    :: Xen
 type(bz_samp), intent(inout) :: Xk,Xq
 !
 type(bz_samp)         :: k_nosymm
 logical               :: l_full_expansion
 integer               :: i_l,i_l_ref,i_c,i_v,i_cvk,i_kibz,i_kbz,i_karpes,i_ks,i_sp_pol,i_w,i_l_v, &
 &                        io_err,ID_ARPES,ID_ARPES_int,ID_bands,NV,NK
 real(SP)              :: dw,Ev,N_l
 complex(SP)           :: F_tmp(RT_conf%W_n_steps),PES(RT_conf%W_n_steps)
 real(SP),allocatable  :: ARPES_int(:,:,:)
 !
 character(schlen)     :: headings(3),file_name
 real(SP)              :: odata(3)
 !
 complex(SP),external      ::Lorentzian_FT
 !
 call parser("FullExpansion",l_full_expansion)
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
 NV=Xen%nbf(1)-BS_bands(1)+1
 NK=Xk%nibz
 if (l_full_expansion) NK=Xk%nbz
 !
 YAMBO_ALLOC(ARPES_int,(BS_bands(1):Xen%nbf(1),NK,i_l_ref:i_l_ref))   ! last dimension is for number of excitons
 ! 
 PES       = cZERO
 ARPES_int = cZERO
 !
 do i_l=i_l_ref,i_l_ref
   N_l=1._SP
   do i_cvk = 1,BS_K_dim(1)
     !
     i_kbz    = BSS_eh_table(i_cvk,1)
     i_v      = BSS_eh_table(i_cvk,2)
     i_c      = BSS_eh_table(i_cvk,3)
     i_sp_pol = spin(BSS_eh_table(i_cvk,:))
     !
     i_kibz   = Xk%sstar(i_kbz,1)
     i_ks     = Xk%sstar(i_kbz,2)
     !
     if ( (.not.l_full_expansion) .and. i_ks>1 ) cycle
     if (.not. PAR_IND_CON_BANDS_DIP%element_1D(i_c)) cycle
     if (.not. PAR_IND_VAL_BANDS_DIP%element_1D(i_v)) cycle
     if (.not. PAR_IND_DIPk_ibz%element_1D(i_kibz)) cycle
     !
     if(     l_full_expansion) i_karpes=i_kbz
     if(.not.l_full_expansion) i_karpes=i_kibz
     !
     ! Quantity to be interpolated
     !
     ARPES_int(i_v,i_karpes,i_l)=ARPES_int(i_v,i_karpes,i_l)+abs(BS_mat(i_cvk,i_l))**2
     !
     ! Constract the PES integrating over k
     !
     do i_w=1,RT_conf%W_n_steps
       F_tmp(i_w)=Lorentzian_FT(GreenF_T_and_W%W(i_w),(Xen%E(i_v,i_kibz,i_sp_pol)+BS_E(i_l)),"r")
     enddo
     !
     F_tmp=F_tmp*abs(BS_mat(i_cvk,i_l))**2*N_l
     !
     PES(:)=PES(:)+F_tmp(:)
     !
   enddo
 enddo
 !
 call PP_redux_wait(PES)
 call PP_redux_wait(ARPES_int)
 !
 ! Print the integrated PES
 !
 file_name="G_w_integrated"
 call of_open_close(trim(file_name),'ot')
 headings=(/" E[eV]  ","Im[G(w)]","Re[G(w)]"/)
 call msg('o '//trim(file_name),"#",headings,INDENT=0,USE_TABS=.true.)
 !
 do i_w=1,RT_conf%W_n_steps
   odata=(/real(GreenF_T_and_W%W(i_w),SP)*HA2EV,real(cI*PES(i_w),SP),real(PES(i_w),SP)/)
   call msg('o '//trim(file_name),'',odata,INDENT=-2,USE_TABS=.true.)
 enddo
 !
 call of_open_close(trim(file_name))
 !
 call section('=',"Loading special Points for the "//trim(lattice)//" lattice")
 !============================================================================= 
 !
 call k_special()
 !
 call k_circuit_driver( BANDS_path, BANDS_steps )
 !
 call section('=','Interpolation@work: Circuit')
 !============================================================================= 
 !
 ! First energies
 !
 call INTERPOLATION_driver_seed("E_K_PATH",ID_bands,BS_bands,Xk,SEED_E=Xen,NEAREST=.false.)
 call INTERPOLATION_driver_do("E_K_PATH",ID_bands,CIRCUIT_k)
 !
 call plot_interpolated_values &
 &   (Xk,Xen,0._SP,CIRCUIT_k,ID_bands,BS_bands,1,n_sp_pol,"interpolated","bands",0)
 !
 ! Then ARPES intensities
 !
 if (l_full_expansion) then
   !
   INTERPOLATE_use_symm=.false.
   call k_ibz2bz(Xk,'i',.true.)
   call bz_samp_duplicate_nosymm(Xk,k_nosymm)
   !
   call INTERPOLATION_driver_seed("ARPES",ID_ARPES_int,(/BS_bands(1),Xen%nbf(1)/),k_nosymm,SEED_R3=ARPES_int)
   !
 else
   !
   call INTERPOLATION_driver_seed("ARPES",ID_ARPES_int,(/BS_bands(1),Xen%nbf(1)/),Xk,SEED_R3=ARPES_int)
   !
 endif
 call INTERPOLATION_driver_do("ARPES",ID_ARPES_int,CIRCUIT_k)
 !
 YAMBO_FREE(ARPES_int)
 !
 ID_ARPES=3
 YAMBO_ALLOC(INTERP_obj(ID_ARPES)%INTERP,(RT_conf%W_n_steps,CIRCUIT_K%nibz,2))
 INTERP_obj(ID_ARPES)%INTERP=0._SP
 !
 i_karpes=0
 !
 do i_l=i_l_ref,i_l_ref
   !
   N_l=1._SP
   !
   do i_kibz=1,CIRCUIT_K%nibz
     !
     if (all(abs(CIRCUIT_K%pt(i_kibz,:))<1.E-5_SP)) i_karpes=i_kibz
     !
     do i_v=BS_bands(1),Xen%nbf(1)
       !
       Ev=INTERP_obj(ID_bands)%INTERP(i_v,i_kibz,1)
       !
       do i_w=1,RT_conf%W_n_steps
         F_tmp(i_w)=Lorentzian_FT(GreenF_T_and_W%W(i_w),(Ev+BS_E(i_l)),"r")
       enddo
       !
       F_tmp=F_tmp*INTERP_obj(ID_ARPES_int)%INTERP(i_v,i_kibz,i_l)*N_l
       !
       INTERP_obj(ID_ARPES)%INTERP(:,i_kibz,1)=INTERP_obj(ID_ARPES)%INTERP(:,i_kibz,1)+ real(F_tmp,SP)
       INTERP_obj(ID_ARPES)%INTERP(:,i_kibz,2)=INTERP_obj(ID_ARPES)%INTERP(:,i_kibz,2)-aimag(F_tmp)
       !
     enddo
   enddo
   !
 enddo
 !
 ! Print G(w,k) for specific k  (only gamma for now)
 !
 file_name="G_k1_w"
 call of_open_close(trim(file_name),'ot')
 headings=(/" E[eV]  ","Im[G(w)]","Re[G(w)]"/)
 call msg('o '//trim(file_name),"#",headings,INDENT=0,USE_TABS=.true.)
 !
 do i_w=1,RT_conf%W_n_steps
   odata=(/real(GreenF_T_and_W%W(i_w),SP)*HA2EV,      &
   &      INTERP_obj(ID_ARPES)%INTERP(i_w,i_karpes,1), &
   &      INTERP_obj(ID_ARPES)%INTERP(i_w,i_karpes,2) /)
   call msg('o '//trim(file_name),'',odata,INDENT=-2,USE_TABS=.true.)
 enddo
 !
 call of_open_close(trim(file_name))
 !
 call section('=','Plotting')
 !============================================================================= 
 !
 if(CIRCUIT_k%nbz>0) call plot_interpolated_values(Xk,Xen,0._SP,CIRCUIT_k,&
 &                                                  (/ID_ARPES/),(/1,RT_conf%W_n_steps/),1,(/2/),"interpolated","G_k_w",0)
 !if(   USER_k%nbz>0) call plot_interpolated_values(Xk,Xen,0._SP,USER_k,   &
 !&                                                  (/ID_ARPES/),(/1,NBW/),1,(/2/),"built_in",    "G_k_w",0)
 !
 call INTERPOLATION_driver_end(0)
 !
end subroutine RT_transient_ARPES
