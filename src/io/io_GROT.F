!
! License-Identifier: GPL
!
! Copyright (C) 2015 The Yambo Team
!
! Authors (see AUTHORS file for details): AM
!
#include<dev_defs.h>
!
integer function io_GROT(io_db)
 !
 use D_lattice,      ONLY:nsym,inv_index
 use R_lattice,      ONLY:n_g_shells,ng_closed,ng_in_shell,E_of_shell,g_rot,minus_G
 use R_lattice,      ONLY:DEV_VAR(g_rot),DEV_VAR(minus_G)
 use IO_m,           ONLY:io_status,io_sec,read_is_on
 use IO_int,         ONLY:io_connect,io_disconnect,io_elemental,&
&                         io_bulk,io_header,ver_is_gt_or_eq
 use timing_m,       ONLY:timing
 !
#include<y_memory.h>
 !
 integer :: io_db,nsym_db
 !
 call timing('io_GROT',OPR='start')
 !
 io_GROT=io_connect(desc='gops',type=1,ID=io_db)
 if (io_GROT/=0) goto 1
 !
 io_GROT=io_header(io_db,IMPOSE_SN=.true.)
 if (io_GROT/=0) goto 1
 !
 if (.not.ver_is_gt_or_eq(io_db,revision=15192)) then
   io_GROT=-1 
   goto 1
 endif
 !
 if (any((/io_sec(io_db,:)==1/))) then
   call io_elemental(io_db,VAR="PARS",VAR_SZ=3,MENU=0)
   call io_elemental(io_db,I0=nsym,DB_I0=nsym_db)
   call io_elemental(io_db,I0=ng_closed, CHECK=.true.,OP=(/"<="/))
   call io_elemental(io_db,I0=n_g_shells,CHECK=.true.,OP=(/"<="/))
   call io_elemental(io_db,VAR="",VAR_SZ=0,MENU=0)
   io_GROT=io_status(io_db)
   if (io_GROT/=0) goto 1
 endif
 !
 if (.not.any((/io_sec(io_db,:)==2/))) goto 1
 !
 if (read_is_on(io_db)) then
   YAMBO_ALLOC(ng_in_shell,(n_g_shells))
   YAMBO_ALLOC(g_rot,(ng_closed,nsym_db))
   YAMBO_ALLOC(E_of_shell,(n_g_shells))
   YAMBO_ALLOC(minus_G,(ng_closed))
 endif
 !
 call io_bulk(io_db,VAR="ng_in_shell",VAR_SZ=(/n_g_shells/))
 call io_bulk(io_db,I1=ng_in_shell(:n_g_shells))
 !
 call io_bulk(io_db,VAR="E_of_shell",VAR_SZ=(/n_g_shells/))
 call io_bulk(io_db,R1=E_of_shell(:n_g_shells))
 !
 call io_bulk(io_db,VAR="GROT",VAR_SZ=(/ng_closed,nsym_db/))
 call io_bulk(io_db,I2=g_rot(:ng_closed,:))
 !
 if (ver_is_gt_or_eq(io_db,revision=955).and.inv_index<0) then
   call io_bulk(io_db,VAR="minusG",VAR_SZ=(/ng_closed/))
   call io_bulk(io_db,I1=minus_G(:ng_closed))
 endif
 !
#ifdef _GPU
 if (read_is_on(io_db)) then
   YAMBO_ALLOC_GPU_SOURCE(DEV_VAR(g_rot),g_rot)
   if (inv_index<0) then
     YAMBO_ALLOC_GPU_SOURCE(DEV_VAR(minus_G),minus_G)
   endif
 endif
#endif
 !
 call initactivate(1,'Gmode')
 !
1 call io_disconnect(ID=io_db)
 !
 call timing('io_GROT',OPR='stop')
 !
end function
