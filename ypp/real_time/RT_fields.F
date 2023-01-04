!
!        Copyright (C) 2000-2022 the YAMBO team
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
subroutine RT_Fields()
 !
 use pars,           ONLY:schlen,n_fields_defs_max
 use RT_control,     ONLY:OBS_RT_IO_t
 use fields,         ONLY:Efield
 use stderr,         ONLY:STRING_split
 !
 implicit none
 !
 ! Work Space
 !
 logical           :: l_load_from_file,l_from_time_fft 
 character(schlen) :: field_defs(n_fields_defs_max)
 !
 ! Time Initialization from input
 call RT_time_configuration_setup('FIELDS',OBS_RT_IO_t)
 !
 ! Analize the input variables for the external field
 call STRING_split(Efield(1)%ef_name,field_defs)
 !
 ! Field in frequency space from external file
 if ( field_defs(1)=='FROM_W_FILE') call RT_Field_w_load()
 ! 
 ! Field in frequency space from FT of real time expression
 if ( field_defs(1)/='FROM_W_FILE') call RT_Field_t_to_w()
 !
 ! Chirp the field and convert back into real time
 call RT_field_chirp()
 !
end subroutine RT_fields
!
!
subroutine RT_Field_w_load()
 !
 ! TODO
 !
end subroutine RT_Field_w_load
!
!
subroutine RT_Field_t_to_w()
 !
 ! TODO
 !
end subroutine RT_Field_t_to_w
!
!
subroutine RT_Field_chirp()
 !
 ! TODO
 !
end subroutine RT_Field_chirp

