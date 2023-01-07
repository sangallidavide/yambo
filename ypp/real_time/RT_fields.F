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
subroutine RT_fields()
 !
 use units,          ONLY:FS2AUT,AUT2FS
 use stderr,         ONLY:intc
 use com,            ONLY:msg
 use pars,           ONLY:SP,schlen,pi,n_fields_defs_max
 use com,            ONLY:com_path
 use YPP_real_time,  ONLY:RT_conf,ElField,RT_time
 use fields,         ONLY:Efield
 use stderr,         ONLY:STRING_split
 !
#include<memory.h>
 !
 ! Work Space
 !
 integer           :: i_field,i_t
 logical           :: l_load_from_file,l_from_time_fft 
 character(schlen) :: field_defs(n_fields_defs_max),filename
 real(SP)          :: field_data(7),Field_time
 !
 RT_conf%T_n_steps=nint( (RT_conf%T_range(2)-RT_conf%T_range(1))/RT_conf%delta_T )+1
 !
 ! Time Initialization from input
 !call RT_time_configuration_setup('FIELDS',OBS_RT_IO_t)
 YAMBO_ALLOC(RT_time,(RT_conf%T_n_steps))
 do i_t=1,RT_conf%T_n_steps
   RT_time(i_t)=(i_t-1)*RT_conf%delta_T+RT_conf%T_range(1)
 enddo
 !
 call msg('s',' Time range         :',(/minval(RT_time),maxval(RT_time)/)*AUT2FS,"[fs]")
 call msg('s',' Time step          :',RT_conf%delta_T*AUT2FS,"[fs]")
 call msg('s',' Time steps         :',RT_conf%T_n_steps)
 !
 ! Analize the input variables for the external field
 call STRING_split(Efield(1)%ef_name,field_defs)
 !
 i_field=1
 !
 ! Field in frequency space from external file
 if ( field_defs(1)=='FROM_W_FILE') call RT_Field_w_load(i_field)
 ! 
 ! Field in frequency space from FT of real time expression
 if ( field_defs(1)/='FROM_W_FILE') call RT_Field_t_to_w(i_field)
 !
 ! Chirp the field and convert back into real time
 ! Warning a pi factor is missing in the FT
 call RT_1D_Fourier_Transform('W2T_LIN_CHIRP',ElField,RT_conf,1)
 !
 filename="YPP_EXTFIELD_LINCHIRP_"//trim(intc(i_field))//".time"
 !
 open(unit=99+i_field,file=trim(com_path)//"/"//trim(filename))
 write(99+i_field,*) RT_conf%T_n_steps
 !
 do i_t=1,RT_conf%T_n_steps
   !
   Field_time=RT_time(i_t)-Efield(i_field)%t_initial
   !
   field_data(1)=Field_time/FS2AUT
   field_data(2)=real(ElField%F_t(i_t,1))/pi
   field_data(3)=real(ElField%F_t(i_t,3))/pi
   field_data(4)=real(ElField%F_t(i_t,5))/pi
   field_data(5)=real(ElField%F_t(i_t,2))/pi
   field_data(6)=real(ElField%F_t(i_t,4))/pi
   field_data(6)=real(ElField%F_t(i_t,6))/pi
   !
   write(99+i_field,*) field_data
   !
 enddo
 !
 close(99+i_field)
 !
end subroutine RT_fields
!
!
subroutine RT_Field_w_load(i_field)
 !
 ! TODO
 !
 implicit none
 !
 integer, intent(in) :: i_field
 !
 call error(' Loading field frequency from file not yet coded')
 !
end subroutine RT_Field_w_load
!
!
subroutine RT_Field_t_to_w(i_field)
 !
 use units,          ONLY:FS2AUT
 use stderr,         ONLY:intc
 use pars,           ONLY:SP,schlen
 use fields,         ONLY:Efield,small_a
 use com,            ONLY:com_path
 use YPP_real_time,  ONLY:RT_conf,ElField,RT_time
 !
 implicit none
 !
 integer, intent(in) :: i_field
 !
 integer   :: i_t
 character(schlen) :: filename
 real(SP)  :: field_data(7),Field_time
 !
 ElField%title="Field"
 ElField%N    =6
 call RT_1D_Fourier_setup(ElField,RT_conf,6) 
 !
 do i_t=1,RT_conf%T_n_steps
   !
   Field_time=RT_time(i_t)-Efield(i_field)%t_initial
   !
   ElField%F_t(i_t,1:2) =small_a(Field_time,Efield(i_field),0)
   ElField%F_t(i_t,3:4) =small_a(Field_time,Efield(i_field),1)
   ElField%F_t(i_t,5:6) =small_a(Field_time,Efield(i_field),2)
   !
 enddo
 !
 filename="YPP_EXTFIELD_"//trim(intc(i_field))//".time"
 !
 open(unit=99+i_field,file=trim(com_path)//"/"//trim(filename))
 write(99+i_field,*) RT_conf%T_n_steps
 !
 do i_t=1,RT_conf%T_n_steps
   !
   Field_time=RT_time(i_t)-Efield(i_field)%t_initial
   !
   field_data(1)=Field_time/FS2AUT
   field_data(2)=real(ElField%F_t(i_t,1))
   field_data(3)=real(ElField%F_t(i_t,3))
   field_data(4)=real(ElField%F_t(i_t,5))
   field_data(5)=real(ElField%F_t(i_t,2))
   field_data(6)=real(ElField%F_t(i_t,4))
   field_data(6)=real(ElField%F_t(i_t,6))
   !
   write(99+i_field,*) field_data
   !
 enddo
 !
 close(99+i_field)
 !
 call RT_1D_Fourier_Transform('T2W',ElField,RT_conf,1)
 !
end subroutine RT_Field_t_to_w
