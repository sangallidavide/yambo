!
!        Copyright (C) 2000-2022 the YAMBO team
!              http://www.yambo-code.org
!
! Authors (see AUTHORS file for details): AM, DS
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
subroutine io_RT_components_V_bands(ID,T_size)
 !
 use pars,           ONLY:SP,schlen
 use drivers,        ONLY:l_nl_p_and_p
 use R_lattice,      ONLY:nkibz
 use stderr,         ONLY:intc
 use iso_c_binding,  ONLY:C_LOC,C_F_POINTER
 use IO_int,         ONLY:def_variable_bulk,io_variable_bulk
 use IO_m,           ONLY:io_sec,io_unit,netcdf_call,nf90_sync
 use electrons,      ONLY:n_sp_pol,n_full_bands
 use RT_control,     ONLY:Vbands_RT_IO_t,NE_i_start_time,RT_do_IO
 use real_time,      ONLY:RT_bands,RT_nbands,NE_i_time!,RT_dyn_step,RT_step
 use nl_optics,      ONLY:V_bands
 use netcdf,         ONLY:NF90_UNLIMITED
 use QP_m,           ONLY:QP_nk
 !
 implicit none
 !
 integer, intent(in)  :: ID,T_size
 !
 logical           :: l_H_rot,l_H_sig,l_H_EQ
 integer           :: ib,i_sp_pol,i_mem,i_mem_sorted,VAR_SZ(6),IPOS(6),&
&                     PAR_G_k_range_ibz(2),k_range(2),ip,TMP_VAR_SZ
 character(schlen) :: dim_names(6)
 complex(SP),pointer :: tmp_link(:,:,:,:,:)
 !
 ! Setup
 !=======
 !
 ! Variables SIZE and strings...
 !-------------------------------
 VAR_SZ(1:5)=(/2,RT_bands(2),maxval(n_full_bands),QP_nk,n_sp_pol/)
 dim_names(1:5)  =[character(schlen) :: 'complex', 'NL_bands','nb_full','QP_nk','n_sp_pol']
 IPOS(1:5)=1
 !
 VAR_SZ(6)    = 1
 if (RT_do_IO("V_bands_K_section")) VAR_SZ(6)    = NF90_UNLIMITED
 dim_names(6) ='V_bands_IO_time_step'
 IPOS(6)=Vbands_RT_IO_t%N
 ! 
 ! Sizes...
 !----------
 !ik = maxval(io_sec(ID,:))-2
 !k_range      =(/1,QP_nk/)
 !G_IO_steps   =Vbands_RT_IO_t%INTERVAL_steps/nint(RT_dyn_step/RT_step)
 ! 
 ! Variable Definitions...
 !-------------------------
 if( NE_i_time==NE_i_start_time.or..not.l_nl_p_and_p) then
   !call def_variable_bulk(ID,"V_bands",4,VAR_SZ,SP,dim_names=dim_names,par_io_kind='collective',silent=.true.)
   call def_variable_bulk(ID,"V_bands",4,VAR_SZ,SP,dim_names=dim_names,silent=.true.)
 endif
 !
 ! Variable I/O...
 !-----------------
 call C_F_POINTER(C_LOC(V_bands),tmp_link,VAR_SZ(2:6))
 call io_variable_bulk(ID,4,C5=tmp_link,IPOS=IPOS)
 nullify(tmp_link)
 !
 ! In case of frequent I/O this gives a significant slowdown of the simulation
 ! if (G_IO_steps>=10) call netcdf_call(nf90_sync(io_unit(ID)),ID)
 !
end subroutine io_RT_components_V_bands
