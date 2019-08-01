/*
        Copyright (C) 2000-2019 the YAMBO team
              http://www.yambo-code.org

 Authors (see AUTHORS file for details): HM AM
 
 This file is distributed under the terms of the GNU 
 General Public License. You can redistribute it and/or 
 modify it under the terms of the GNU General Public 
 License as published by the Free Software Foundation; 
 either version 2, or (at your option) any later version.

 This program is distributed in the hope that it will 
 be useful, but WITHOUT ANY WARRANTY; without even the 
 implied warranty of MERCHANTABILITY or FITNESS FOR A 
 PARTICULAR PURPOSE.  See the GNU General Public License 
 for more details.

 You should have received a copy of the GNU General Public 
 License along with this program; if not, write to the Free 
 Software Foundation, Inc., 59 Temple Place - Suite 330,Boston, 
 MA 02111-1307, USA or visit http://www.gnu.org/copyleft/gpl.txt.
 
*/
 use memory, ONLY:MEM_err,MEM_dri
 implicit none
#define YAMBO_ALLOC_P(x,SIZE) \
  allocate(x SIZE,stat=MEM_err)NEWLINE \
  if (     associated(x)) &NEWLINE& call MEM_dri(QUOTES x QUOTES,x)NEWLINE \
  if (.not.associated(x)) &NEWLINE& call MEM_dri(QUOTES x QUOTES)
#define YAMBO_ALLOC(x,SIZE) \
  allocate(x SIZE,stat=MEM_err)NEWLINE \
  if (     allocated(x)) &NEWLINE& call MEM_dri(QUOTES x QUOTES,x)NEWLINE \
  if (.not.allocated(x)) &NEWLINE& call MEM_dri(QUOTES x QUOTES)
#define YAMBO_FREE(x) \
  if (.not.allocated(x)) &NEWLINE& call MEM_dri(QUOTES x QUOTES,int(-1,KIND=8))NEWLINE \
  if (     allocated(x)) &NEWLINE& call MEM_dri(QUOTES x QUOTES,size(x,KIND=8))NEWLINE \
  if (     allocated(x)) &NEWLINE& deallocate(x)
#define YAMBO_FREE_P(x) \
  if (.not.associated(x)) &NEWLINE& call MEM_dri(QUOTES x QUOTES,int(-1,KIND=8))NEWLINE \
  if (     associated(x)) &NEWLINE& call MEM_dri(QUOTES x QUOTES,size(x,KIND=8))NEWLINE \
  if (     associated(x)) &NEWLINE& deallocate(x);nullify(x)
