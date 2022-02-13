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
subroutine RT_G_two_times_interp_and_plot(en,kpt,NK,F_w_real)
 !
 use pars,           ONLY:SP
 use YPP_interfaces, ONLY:k_special
 use YPPm,           ONLY:BANDS_path,BANDS_steps,CIRCUIT_k,USER_k
 use YPP_real_time,  ONLY:RT_conf,l_RT_G_two_times
 use electrons,      ONLY:levels
 use R_lattice,      ONLY:bz_samp
 use D_lattice,      ONLY:lattice
 use interpolate,    ONLY:INTERPOLATION_driver_seed,INTERPOLATION_driver_do
 !
 implicit none
 !
 type(levels)  :: en
 type(bz_samp) :: kpt
 integer,  intent(in) :: NK
 real(SP), intent(in) :: F_w_real(RT_conf%W_n_steps,NK,2)
 !
 ! Work space
 !
 integer  :: ID_RT
 !
 if (l_RT_G_two_times) return
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
 call INTERPOLATION_driver_seed("RT_G_w",ID_RT,(/1,RT_conf%W_n_steps/),kpt,SEED_R3=F_w_real)
 call INTERPOLATION_driver_do("RT_G_w",ID_RT,CIRCUIT_k)
 !
 call section('=','Plotting')
 !============================================================================= 
 !
 if(CIRCUIT_k%nbz>0) call plot_interpolated_values(kpt,en,0._SP,CIRCUIT_k,&
 &                                                  (/ID_RT/),(/1,RT_conf%W_n_steps/),1,(/2/),"interpolated","G_k_w",0)
 if(   USER_k%nbz>0) call plot_interpolated_values(kpt,en,0._SP,USER_k,   &
 &                                                  (/ID_RT/),(/1,RT_conf%W_n_steps/),1,(/2/),"built_in",    "G_k_w",0)
 !
 call INTERPOLATION_driver_end(0)
 !
end subroutine RT_G_two_times_interp_and_plot
