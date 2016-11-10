!
!        Copyright (C) 2000-2016 the YAMBO team
!              http://www.yambo-code.org
!
! Authors (see AUTHORS file for details): AM CA DS
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
subroutine RT_propagate_fields(E,A,T,T_step)
 !
 ! The external and total fields in Yambo are defined and organized as follow:
 !
 ! The general vector potential is defined as
 !
 !   A (t)=-c Eo    a(t) theta(t-to)                                     [1]
 !   A'(t)=-c Eo ( a'(t) theta(t-to)+ a(0) delta(t-to))                  [2]
 !   A"(t)=-c Eo ( a"(t) theta(t-to)+ a'(0) delta(t-to) - a"(0) sign(t)) [3]
 !
 ! the function a(t) and theta(t) are defined in mod_fields.
 !
 ! Using [1-3] yambo :
 !
 ! a. calculates the analytical form of A/A'/A" in RT_Ext_fields 
 ! b. uses A" (the ext current/polarization) to evolve the Maxwell equations when
 !    the induced current is included in the dynamics ( RT_propagate_fields )
 ! 
 use pars,           ONLY:SP,pi
 use units,          ONLY:SPEED_OF_LIGHT
 use electrons,      ONLY:levels
 use X_m,            ONLY:global_gauge
 use real_time,      ONLY:RT_ind_J_prev,RT_ind_J,RT_P_prev,RT_P,&
&                         l_RT_induced_field,integrator_step,&
&                         RT_dyn_step,l_RT_EULER
 use fields,         ONLY:gauge_field,A_ind,A_ind_previous,A_ext
 !
 implicit none
 !
 type(levels),         intent(in) :: E
 real(SP),             intent(in) :: T,T_step
 type(gauge_field),   intent(out) :: A
 !
 ! Get the analytic form of Ext field @ T
 !----------------------------------------
 !
 call RT_Ext_fields(A_ext,T)
 !
 ! (See also by Yaban and co. PRB 2012)
 ! The approach which          includes the induced field is also called longitudinal gauge.
 ! The approach which does not include  the induced field is also called transverse   gauge.
 ! This has nothing to do with the length and the velocity gauge.
 ! It is instead related to distinction between \chi and overline{\chi} or
 ! between transverse and longitudinal excitons/phonons in the literature
 !
 if (l_RT_induced_field) then
   !
   ! OPTION 1: "velocity gauge"
   !
   ! Induced potential A(t) from j(t)
   !
   ! Integrate the eq. for the vector potential
   !
   ! d^2 A_ind / d t^2 = 4 pi * c * J_ind(t)
   !
   ! A''(t+dt) =         4 pi * j(t+dt) * c           [this is defined but it is useless]
   ! A' (t+dt) = A'(t) + 4 pi * j(t)    * c * dt
   ! A  (t+dt) = A (t) + 4 pi * P(t)    * c * dt      [this is defined but it is useless as well]
   !
   ! Results are almost ok, but, due to numerical errors, derived epsilon and eels are very bad
   !
   if(trim(global_gauge)=="velocity") then
     A_ind%vecpot_acc =                             real(RT_ind_J)     *4._SP*pi*SPEED_OF_LIGHT
     A_ind%vecpot_vel = A_ind_previous%vecpot_vel + real(RT_ind_J_prev)*4._SP*pi*SPEED_OF_LIGHT*T_step
     A_ind%vecpot     = A_ind_previous%vecpot     + real(RT_P_prev)    *4._SP*pi*SPEED_OF_LIGHT*T_step
   endif
   !
   !
   ! Option 2: length gauge
   !
   !  Induced potential A(t) from P(t)
   !
   !
   ! Integrate the eq. for the vector potential
   !
   ! d A_ind / d t = 4 pi * c * P(t)
   !
   ! A''(t+dt) =         4 pi * j(t+dt) * c        [this is defined but it is useless]
   ! A' (t+dt) =         4 pi * P(t+dt) * c
   ! A  (t+dt) = A (t) + 4 pi * P(t)    * c * dt   [this is defined but it is useless as well]
   !
   ! It works fine
   !
   if(trim(global_gauge)=="length") then
     A_ind%vecpot_acc =                         real(RT_ind_J) *4._SP*pi*SPEED_OF_LIGHT 
     A_ind%vecpot_vel =                         real(RT_P)     *4._SP*pi*SPEED_OF_LIGHT
     A_ind%vecpot     = A_ind_previous%vecpot + real(RT_P_prev)*4._SP*pi*SPEED_OF_LIGHT*T_step
   endif
   !
   if (integrator_step==1) A_ind_previous  =A_ind
   !
   A%vecpot_acc=A_ext%vecpot_acc+A_ind%vecpot_acc
   A%vecpot_vel=A_ext%vecpot_vel+A_ind%vecpot_vel
   A%vecpot    =A_ext%vecpot    +A_ind%vecpot
   !
 else
   !
   A=A_ext
   !
 endif
 !
 ! Fluence accumulation
 !======================
 !
 ! I(t) = \frac{c}{8\pi} E^2(t) and E(t)=A_vel(t)/c
 ! F(t) = \int I(t) dt
 !
 if(integrator_step==3.or.l_RT_EULER) then
   A%intensity=(abs(A%vecpot_vel(1))**2+abs(A%vecpot_vel(2))**2+abs(A%vecpot_vel(3))**2)/SPEED_OF_LIGHT/4._SP/pi
   A%fluence=A%fluence+A%intensity*RT_dyn_step
   A_ext%intensity=(abs(A_ext%vecpot_vel(1))**2+abs(A_ext%vecpot_vel(2))**2+abs(A_ext%vecpot_vel(3))**2)/SPEED_OF_LIGHT/4._SP/pi
   A_ext%fluence=A_ext%fluence+A_ext%intensity*RT_dyn_step
 endif
 !
end subroutine RT_propagate_fields