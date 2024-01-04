!
! License-Identifier: GPL
!
! Copyright (C) 2023 The Yambo Team
!
! Authors (see AUTHORS file for details): DS
!
subroutine build_symmetries()
 !
 use pars,           ONLY:SP
 use stderr,         ONLY:intc
 use com,            ONLY:msg
 use matrix_operate, ONLY:m3det
 use D_lattice,      ONLY:alat,sp_nsym,nsym,nrot,i_time_rev,dl_sop,&
&                         atoms_spatial_inversion,i_space_inv,inv_index,&
&                         symmetry_group_table,mag_syms,idt_index,      &
&                         symmetry_equiv_atoms,pure_rotation
 use R_lattice,      ONLY:rl_sop
 !
#include <memory.h>
 !
 integer           :: i1,i2,i3,is,iclass,jclass,ich,len_class(12),max_len
 character(5)      :: vec_ch(12)
 !
 ! Spin symmetries & xc spin for libxc
 !
 call build_spin_sop()
 !
 YAMBO_ALLOC(rl_sop,(3,3,nsym))
 !
 do is=1,nsym
   forall (i2=1:3,i3=1:3) rl_sop(i2,i3,is)=dl_sop(i2,i3,is)*alat(i2)/alat(i3)
 enddo
 !
 ! Indentity index
 !
 idt_index=-1
 do is=1,nsym
   if (all(nint(reshape(dl_sop(:,:,is),(/9/)))==(/1,0,0,0,1,0,0,0,1/)) ) idt_index=is
 enddo
 !
 if (idt_index==-1) then
   call warning('Identity not found among the given symmetry list')
 endif
 !
 ! Inversion index
 !
 inv_index=-1
 do is=1,nsym
   if ( all(nint(reshape(dl_sop(:,:,is),(/9/)))==(/-1,0,0,0,-1,0,0,0,-1/)) ) inv_index=is
 enddo
 !
 ! Defines i_space_inv value
 call atoms_spatial_inversion()
 !
 ! Symmetries Multiplication Table
 !
 !  R_i*R_j=R_stab(i,j)
 !
 call symmetry_group_table('r')
 call symmetry_equiv_atoms()
 !
 sp_nsym=nsym/(1+i_time_rev)
 !
 YAMBO_ALLOC(pure_rotation,(nsym))
 pure_rotation=.false.
 !
 nrot=0
 do is=1,sp_nsym
   if( abs(m3det(dl_sop(:,:,is))-1._SP) > 1.E-5 ) cycle
   nrot=nrot+1
   pure_rotation(is)=.true.
 enddo
 !
 call msg('r','Inversion symmetry    ',inv_index>0)
 call msg('r','Spatial inversion     ',i_space_inv==1)
 if (inv_index>0) call msg('r','Inversion index       ',inv_index)
 !
 ! Time reversal
 !
 call msg('r','K-space Time-reversal ',(i_time_rev==1))
 call msg('r','Magnetic symmetries   ',mag_syms)
 if (i_time_rev==1) call msg('r','Time-reversal derived K-space symmetries',(/nsym/2+1,nsym/))
 !
end subroutine build_symmetries
!
!
subroutine detect_point_group(sp_nsym,dl_sop)
 !
 use pars,           ONLY:DP
 use electrons,      ONLY:n_spinor
 use D_lattice,      ONLY:group_code,group_name,nclass,nelem,elem,which_irr,&
&                         nclass_ref,char_mat,name_rap,name_class,ir_ram
 !
 implicit none
 !
 integer,  intent(in)  :: sp_nsym
 real(DP), intent(in)  :: dl_sop(3,3,sp_nsym)
 !
 ! Double point group detection not coded
 if (n_spinor==2) call warning("Detecting no SOC point group")
 !
 ! Symmetries and point groups
 !
 call find_group(sp_nsym,dl_sop,group_name,group_code)
 !
 call divide_class(group_code,sp_nsym,dl_sop,nclass,nelem,elem,which_irr)
 !
 ! DS: to check if nclass and nclass_ref can be different in any case
 call set_irr_rap( group_code, nclass_ref, char_mat, name_rap, name_class, ir_ram )
 !
end subroutine detect_point_group
!
!
subroutine print_point_group(fname,sp_nsym,dl_sop)
 !
 use com,            ONLY:msg
 use pars,           ONLY:SP,schlen
 use electrons,      ONLY:n_spinor
 use stderr,         ONLY:intc
 use D_lattice,      ONLY:group_code,group_name,nclass,nelem,elem,which_irr,&
&                         nclass_ref,char_mat,name_rap,name_class,ir_ram,&
&                         i_time_rev,nsym
 !
 implicit none
 !
 character(*), intent(in) :: fname
 integer,  intent(in)  :: sp_nsym
 real(SP), intent(in)  :: dl_sop(3,3,sp_nsym)
 !
 integer           :: i1,i2,i3,is,iclass,jclass,ich,len_class(12),max_len
 character(5)      :: vec_ch(12)
 logical           :: trev
 character(2)      :: tmp_ch
 character(schlen) :: dumb_ch
 !
 tmp_ch=""
 if (fname(1:1)=="o") tmp_ch="# "
 !
 ! Symmetries list
 !
 if (fname(1:1)=="r") call msg('rn','Symmetries units      ',' [cc]')
 !
 do is=1,sp_nsym
   dumb_ch=trim(intc(is))//'] '
   i1=9
   do while (i1<nsym)
     if(is<=i1) dumb_ch=' '//trim(dumb_ch)
     i1=i1*10+9
   enddo
   dumb_ch=tmp_ch//' [S '//trim(dumb_ch)
   call msg(fname,trim(dumb_ch),reshape(dl_sop(:,:,is),(/9/)),INDENT=0)  
 enddo
 !
 if (fname(1:1)=="r") then
   do is=sp_nsym+1,nsym
     dumb_ch=trim(intc(is))//'] '
     i1=9
     do while (i1<nsym)
       if(is<=i1) dumb_ch=' '//trim(dumb_ch)
       i1=i1*10+9
     enddo
     dumb_ch=tmp_ch//' [S*'//trim(dumb_ch)
     call msg(fname,trim(dumb_ch),reshape(dl_sop(:,:,is),(/9/)),INDENT=0)  
   enddo
 endif
 !
 call msg(fname,tmp_ch//'')
 call msg(fname,tmp_ch//'Point group is        ',group_name,INDENT=0)
 call msg(fname,tmp_ch//'')
 call msg(fname,tmp_ch//'Number of classes        ',nclass,INDENT=0)
 call msg(fname,tmp_ch//'Symmetries per class     ',nelem(:nclass),INDENT=0)
 call msg(fname,tmp_ch//'Name of classes          ',name_class(:nclass),INDENT=0)
 !
 len_class=-1
 do iclass=1,nclass
   len_class(iclass)=len_trim(name_class(iclass))
   call msg(fname,tmp_ch//'Symmetries in class '//name_class(iclass),elem(:nelem(iclass),iclass),INDENT=0)
 enddo
 max_len=max(maxval(len_class),5)
 !
 call msg(fname,tmp_ch//'')
 call msg(fname,tmp_ch//'Number of rapresentations',nclass_ref,INDENT=0)
 call msg(fname,tmp_ch//'')
 call msg(fname,tmp_ch//'Character table')
 !
 vec_ch=""
 do iclass=1,nclass_ref
   vec_ch(iclass)=trim(name_class(iclass))
   do ich=len_trim(vec_ch(iclass)),max_len-1
     vec_ch(iclass)=" "//trim(vec_ch(iclass))
   enddo
 enddo
 !
 call msg(fname,tmp_ch//'Symm./classes  ',vec_ch(:nclass),INDENT=0)
 !
 do iclass=1,nclass_ref
   do jclass=1,nclass_ref
     write (vec_ch(jclass),'(F5.2)') real(char_mat(iclass,jclass),SP)
     do ich=len_trim(vec_ch(jclass)),max_len-1
       vec_ch(jclass)=" "//vec_ch(jclass)
     enddo
   enddo
   call msg(fname,tmp_ch//name_rap(iclass),vec_ch(:nclass),INDENT=0)
 enddo
 !
 call msg(fname,tmp_ch//'')
 !
end subroutine print_point_group
