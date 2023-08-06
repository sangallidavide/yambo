!
! Copyright (C) 2000-2008 D. Sangalli and the YAMBO team 
!              http://www.yambo-code.org
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
subroutine K_doubles_mount(omega)
 !
 ! In this subrioutine the complete second order kernel is
 ! constructed for a specific frequency, starting from the BS_mat_doubles
 ! coupling term which has been printed by the subroutine K_doubles.F
 ! throught ioBS_doubles.F into the file
 ! ndb.BS_C_doubles = C_{sd}^{i}    with i=1,2
 !
 !  K_{ss'}= 1/2 \sum_{ij} \sum_d
 !       \frac{C^{i}_{sd}C^{j}*_{ds'}\}{\omega-E_d} 
 !
 !  s={c,v}    d={c1,c2,v1,v2}   E_d={e_{c1}+e_{c2}-e_{v1}+e_{v2}}
 !
 use pars,          ONLY:SP,cZERO,cONE
 use BS,            ONLY:BS_K_dim
 use BS_solvers,    ONLY:BS_mat
 use wrapper,       ONLY:M_by_M
 use BS_doubles,    ONLY:BS_mat_doubles,BS_doubles_E,BS_res_K_doubles, &
&                        BS_doubles_table,BS_K_doubles_dim,&
&                        static_approx,BS_mat_tmp,&
! DOUBLES CORR PROJECT
&                        doubles_corr
 !
 implicit none
 ! Kernel construction
 integer             :: i_doubles,icv1,icv2,icin,ncin,BS_mat_dim
 complex             :: omega
 !
 if(.not.BS_res_K_doubles) return
 !
#if defined _DEBUG
 open(unit=3455, file='Kernel_double.txt', status='unknown')
 write(3455,*) BS_doubles_E(:)*27.2116
#endif
 !
 ncin=2
 if(doubles_corr) ncin=1
 !
 BS_mat_dim=BS_K_dim(1)
 !
 do icin=1,ncin
   BS_mat_tmp(:,:,2)=cZERO
   forall(i_doubles=1:BS_K_doubles_dim) BS_mat_doubles(:,i_doubles,ncin+1)=&
&     BS_mat_doubles(:,i_doubles,icin)/(real(omega)-BS_doubles_E(i_doubles))
   call M_by_M('n','c',BS_mat_dim,BS_mat_dim,BS_K_doubles_dim,cONE,&
&    BS_mat_doubles(:,:,ncin+1),BS_mat_dim,BS_mat_doubles(:,:,icin),BS_mat_dim,&
&    cZERO,BS_mat_tmp(:,:,2),BS_mat_dim)    
   BS_mat(:,:)=BS_mat(:,:)+BS_mat_tmp(:,:,2)
 enddo
 !
 forall(icv1=1:BS_mat_dim) BS_mat(icv1,icv1)=real(BS_mat(icv1,icv1))
 !
#if defined _DEBUG
 write(3455,*) BS_mat(:,:)
 close(unit=3455)
#endif
 !
end subroutine
!
!
subroutine K_doubles_unmount()
 !
 use BS_solvers,    ONLY:BS_mat
 use BS_doubles,    ONLY:BS_mat_tmp,BS_res_K_doubles
 !
 implicit none
 !
 if(.not.BS_res_K_doubles) return
 !
 BS_mat(:,:)=BS_mat_tmp(:,:,1)
 !
end subroutine K_doubles_unmount
!
!
subroutine K_doubles_mount_init()
 !
 use IO_int,        ONLY:io_control
 use IO_m,          ONLY:OP_RD_CL,REP,VERIFY
 use BS,            ONLY:BS_K_dim
 use BS_solvers,    ONLY:BS_mat
 use BS_doubles,    ONLY:BS_mat_doubles,BS_mat_tmp,BS_K_doubles_dim,    &
&                        BS_doubles_table,BS_bands_doubles,BS_doubles_E,&
&                        doubles_corr,BS_res_K_doubles
 !
#include<memory.h>
 !
 integer :: BS_mat_dim,ncin,i_err,ID
 integer, external :: io_BS_doubles
 !
 if(.not.BS_res_K_doubles) return
 !
 ncin=2                                                                 
 if (doubles_corr) ncin=1
 !
 BS_mat_dim=BS_K_dim(1)
 !
 YAMBO_ALLOC(BS_mat_tmp,(BS_mat_dim,BS_mat_dim,2))                              
 YAMBO_ALLOC(BS_mat_doubles,(BS_mat_dim,BS_K_doubles_dim,ncin+1))             
 YAMBO_ALLOC(BS_doubles_E,(BS_K_doubles_dim))                               
 YAMBO_ALLOC(BS_doubles_table,(BS_K_doubles_dim,4))                         
 !
 BS_mat_tmp(:,:,1)=BS_mat
 !                                                                      
 call io_control(action=OP_RD_CL,COM=REP,MODE=VERIFY,SEC=(/1,2/),ID=ID) 
 i_err=io_BS_doubles(ID)                                                 
 !
end subroutine K_doubles_mount_init
!
!
subroutine K_doubles_mount_reset()
 !
 use BS_doubles,    ONLY:BS_mat_doubles,BS_mat_tmp,  &
&                        BS_doubles_table,BS_doubles_E, &
&                        BS_res_K_doubles
 !
#include<memory.h>
 !
 if(.not.BS_res_K_doubles) return
 !
 YAMBO_FREE(BS_mat_tmp)                              
 YAMBO_FREE(BS_mat_doubles)             
 YAMBO_FREE(BS_doubles_E)                               
 YAMBO_FREE(BS_doubles_table)                         
 !
end subroutine K_doubles_mount_reset
