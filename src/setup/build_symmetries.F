!
! License-Identifier: GPL
!
! Copyright (C) 2023 The Yambo Team
!
! Authors (see AUTHORS file for details): DS
!
subroutine build_symmetries()
 !
 use pars,           ONLY:DP,schlen
 use electrons,      ONLY:n_spinor
 use stderr,         ONLY:intc
 use com,            ONLY:msg
 use D_lattice,      ONLY:alat,nsym,i_time_rev,dl_sop,&
&                         atoms_spatial_inversion,i_space_inv,inv_index,&
&                         symmetry_group_table,mag_syms,idt_index,symmetry_equiv_atoms,&
&                         group_code,group_name,nrot,nclass,nelem,elem,which_irr,&
&                         nclass_ref,char_mat,name_rap,name_class,ir_ram
 use R_lattice,      ONLY:rl_sop
 !
#include <memory.h>
 !
 integer           :: i1,i2,i3,is,iclass
 character(schlen) :: dumb_ch
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
 ! Symmetries Multiplication Table
 !
 call symmetry_group_table('r')
 call symmetry_equiv_atoms()
 !
 ! Symmetries list
 !
 call msg('rn','Symmetries units      ',' [cc]')
 do is=1,nsym
   dumb_ch=trim(intc(is))//'] '
   i1=9
   do while (i1<nsym)
     if(is<=i1) dumb_ch=' '//trim(dumb_ch)
     i1=i1*10+9
   enddo
   if (is<=nsym/(1+i_time_rev)) dumb_ch=' [S '//trim(dumb_ch)
   if (is> nsym/(1+i_time_rev)) dumb_ch=' [S*'//trim(dumb_ch)
   call msg('r',trim(dumb_ch),reshape(dl_sop(:,:,is),(/9/)))  
 enddo
 !
 ! Symmetries and point groups
 !
 nrot=nsym/(1+max(i_time_rev,i_space_inv))
 call find_group(nrot,real(dl_sop(:,:,:nrot),DP),group_name,group_code)
 call msg('nr','Point group is        ',group_name)
 !
 allocate(nelem(12*n_spinor))
 allocate(elem(8+4*(n_spinor-1),12*n_spinor))
 allocate(which_irr(12*n_spinor))
 if (n_spinor==1) then
   call divide_class(group_code,nrot,real(dl_sop(:,:,:nrot),DP),nclass,nelem,elem,which_irr)
 else
   ! TO DO
   !allocate(has_e())
   !call divide_class_so(group_code,nrot,dl_sop(:,:,:nrot),nclass,nelem,elem,which_irr)
 endif
 call msg('nr','Number of classes     ',nclass)
 call msg('r' ,'Symmetries per class  ',nelem(:nclass))
 !
 ! DS: to check if nclass and nclass_ref can be different in any case
 call set_irr_rap( group_code, nclass_ref, char_mat, name_rap, name_class, ir_ram )
 !
 call msg('r' ,'Name of classes       ',name_class(:nclass))
 do iclass=1,nclass
   call msg('r','Symmetries in class '//name_class(iclass),elem(:nelem(iclass),iclass))
 enddo
 !
 call msg('nr' ,'Number of rapresentations',nclass_ref)
 do iclass=1,nclass_ref
   dumb_ch=trim(intc(iclass))//'] '
   i1=9
   do while (i1<nclass_ref)
     if(iclass<=i1) dumb_ch=' '//trim(dumb_ch)
     i1=i1*10+9
   enddo
   dumb_ch=trim(intc(iclass))//'] '
   call msg('r','[R'//trim(dumb_ch),name_rap(iclass))
 enddo
 !
end subroutine build_symmetries
