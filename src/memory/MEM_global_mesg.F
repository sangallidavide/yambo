!
!        Copyright (C) 2000-2023 the YAMBO team
!              http://www.yambo-code.org
!
! Authors (see AUTHORS file for details): AM DS
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
subroutine MEM_global_mesg(what,xkind,HOST_SIZE)
 !
 use pars,         ONLY:SP,DP,IPL
 use parallel_m,   ONLY:host_name
 use stderr,       ONLY:real2ch
 use units,        ONLY:Kilobyte
 use LIVE_t,       ONLY:LIVE_message
 use memory,       ONLY:mem_string
 !
 implicit none
 !
 character(*), intent(in) :: what
 integer,      intent(in) :: xkind
 integer(IPL), intent(in) :: HOST_SIZE(:)
 !
 integer   :: i1
 real(SP)  :: MEM_size
 !
 MEM_size=real(HOST_SIZE(1),SP)
 do i1=2,size(HOST_SIZE)
   MEM_size=MEM_size*real(HOST_SIZE(i1),SP)
 enddo
 MEM_size=MEM_size*real(xkind,SP)/Kilobyte
 !
 ! Error Message
 !===============
 call LIVE_message("[PAR-MEMORY] Parallel distribution of "//what//" on HOST "//trim(host_name)//&
&                  " with total size "//trim(mem_string(MEM_size)))
 !
 call MEM_error("GLOBAL_"//what)
 !
end subroutine MEM_global_mesg
