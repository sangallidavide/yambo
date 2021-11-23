!
!        Copyright (C) 2000-2021 the YAMBO team
!              http://www.yambo-code.org
!
! Authors (see AUTHORS file for details): AM
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
subroutine SERIAL_SVD(n,M,mode,SVD_digits)
 !
 use pars,           ONLY:SP,cZERO,cONE
 use wrapper,        ONLY:M_by_M
 use linear_algebra, ONLY:LINEAR_ALGEBRA_error,SV_decomposition,&
&                         LINEAR_ALGEBRA_WS_reset,LALGEBRA_WS
 use timing_m,       ONLY:timing
#include<memory.h>
 !
 integer,     intent(in)    :: n
 complex(SP), intent(inout) :: M(n,n)
 integer,     intent(in)    :: SVD_digits 
 character(3),intent(in)    :: mode
 !
 !ws
 integer           :: i_loop 
 real(SP)          :: SV_min
 type(LALGEBRA_WS) :: WS
 !
 call timing('SERIAL_SVD',OPR='start')
 !
 YAMBO_ALLOC(WS%v_real,(5*n))
 YAMBO_ALLOC(WS%vp_real,(n))
 YAMBO_ALLOC(WS%m2_cmplx,(n,n))
 YAMBO_ALLOC(WS%m1_cmplx,(n,n))
 allocate(WS%v_cmplx(1))
 !
 call SV_decomposition(n,M,WS%vp_real,WS%m2_cmplx,WS%m1_cmplx,WS%v_cmplx,-1,WS%v_real,WS%i_fail)
 !
 if(WS%i_fail/=0) call LINEAR_ALGEBRA_error('GESVD (SVD)','WS%v_cmplxspace failed')
 !
 WS%dim=int(real(WS%v_cmplx(1)))
 deallocate(WS%v_cmplx)
 YAMBO_ALLOC(WS%v_cmplx,(WS%dim))
 call SV_decomposition(n,M,WS%vp_real,WS%m2_cmplx,WS%m1_cmplx,WS%v_cmplx,WS%dim,WS%v_real,WS%i_fail)
 !
 if(WS%i_fail/=0) call LINEAR_ALGEBRA_error('GESVD (SVD)','failed')
 !
 if(index(mode,'inv')/=0) then   ! Matrix inversion
   !      
   SV_min=10._SP**(-SVD_digits)
   !
   do i_loop=1,n
     if (abs(WS%vp_real(i_loop))<SV_min) then
       WS%m1_cmplx(i_loop,:)=cZERO
     else
       WS%m1_cmplx(i_loop,:)=WS%m1_cmplx(i_loop,:)/WS%vp_real(i_loop)
     endif
   enddo
   !
 elseif(index(mode,'uni')/=0)  then  
   !      
   ! Transform in a unitary matrix
   ! U= M1*I*M2
   !
 else
   call error('Error SVD mode ')
 endif
 !
 call M_by_M('c','c',n,WS%m1_cmplx,WS%m2_cmplx,M)
 !
 call LINEAR_ALGEBRA_WS_reset(WS)
 !
 call timing('SERIAL_SVD',OPR='stop')
 !
end subroutine
