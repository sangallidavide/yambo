!
! License-Identifier: GPL
!
! Copyright (C) 2011 The Yambo Team
!
! Authors (see AUTHORS file for details): AM DS
!
integer function variables_X(ID,local_desc,X,EM1S_NOSOC)
 !
 use units,         ONLY:HA2EV
 use descriptors,   ONLY:IO_desc,IO_desc_reset
 use pars,          ONLY:SP,IP,schlen,CR_YIO
 use X_m,           ONLY:X_t,use_X_DbGd,X_DbGd_nkpts,X_DbGd_percentual,OPTICAL_averaged_dirs_string
 use IO_m,          ONLY:io_status,read_is_on,variable_exist
 use IO_int,        ONLY:def_variable_elemental,io_variable_elemental,ver_is_gt_or_eq
 !
 implicit none
 !
 integer            ::ID
 type(IO_desc)      ::local_desc
 type(X_t)          ::X
 logical            ::EM1S_NOSOC
 !
 ! Work Space
 !
 integer             ::sec_size,tmp_ib(2)
 character(schlen)   ::X_label,VAR_NAME
 !
 call IO_desc_reset(local_desc)
 X_label=' '
 if (ver_is_gt_or_eq(ID,revision=955)) X_label="X_"
 !
 ! PARS_1
 !
 call def_variable_elemental(ID,trim(X_label)//"PARS_1",5,SP,0)
 !
 call io_variable_elemental(ID,DB_I0=X%ng_db,&
&     VAR="X matrix size",I0=X%ng,CHECK=.true.,OP=(/"<="/),DESCRIPTOR=local_desc)
 !
 tmp_ib=X%ib
 if (read_is_on(ID).and.EM1S_NOSOC) then
   ! Here I'm in the case where I'm using a DB computed without SOC
   ! The bands in X%ib are with SOC, while on disk I have half of the bands
   tmp_ib(1)=(X%ib(1)-1)/2+1
   tmp_ib(2)= X%ib(2)   /2
 endif
 !
 call io_variable_elemental(ID,&
&     VAR="X band range",I1=tmp_ib,CHECK=.true.,OP=(/"==","=="/),DESCRIPTOR=local_desc)
 !
 if (read_is_on(ID).and..not.EM1S_NOSOC) X%ib=tmp_ib
 !
 call io_variable_elemental(ID,UNIT=HA2EV,&
&     VAR="X e/h energy range ",TERMINATOR="[eV]",R1=X%ehe,CHECK=.true.,OP=(/">=","<="/),DESCRIPTOR=local_desc)
 !
 ! PARS_2
 !
 VAR_NAME=trim(X_label)//"Time_ordering"
 if (read_is_on(ID) .and. .not. variable_exist(ID,trim(VAR_NAME))) VAR_NAME=trim(X_label)//"PARS_2"
 call def_variable_elemental(ID,trim(VAR_NAME),1,CR_YIO,0)
 call io_variable_elemental(ID,VAR="X Time ordering",CH0=X%ordering,WARN=.true.,DESCRIPTOR=local_desc)
 !
 if (ver_is_gt_or_eq(ID,(/3,0,5/))) then
   !
   ! TDDFT_KERNEL
   !
   call def_variable_elemental(ID,trim(X_label)//"TDDFT_KERNEL",1,CR_YIO,0)
   call io_variable_elemental(ID,VAR="X Hxc Kernel",CH0=X%KERNEL,CHECK=.true.,DESCRIPTOR=local_desc)
   !
   ! DRUDE
   !
   call def_variable_elemental(ID,trim(X_label)//"DRUDE",2,SP,0)
   call io_variable_elemental(ID,UNIT=HA2EV,&
&     VAR="X Drude frequency",C0=X%Wd,CHECK=.true.,OP=(/"=="/),DESCRIPTOR=local_desc)
   !
 endif
 !
 ! PARS_3
 !
 sec_size=6
 if (ver_is_gt_or_eq(ID,(/3,0,1/)))      sec_size=8
 if (ver_is_gt_or_eq(ID,(/3,0,9/)))      sec_size=9
 if (ver_is_gt_or_eq(ID,revision=15206)) sec_size=10
 call def_variable_elemental(ID,trim(X_label)//"PARS_3",sec_size,SP,0)
 call io_variable_elemental(ID,&
&     VAR="X poles",TERMINATOR="[o/o]",R0=X%cg_percentual,CHECK=.true.,OP=(/"=="/),DESCRIPTOR=local_desc)
 call io_variable_elemental(ID,&
&     VAR="RL vectors in the sum",I0=X%ngostnts,WARN=.true.,OP=(/"<="/),DESCRIPTOR=local_desc)
 call io_variable_elemental(ID,&
&     VAR="[r,Vnl] included",L0=X%Vnl_included,WARN=.true.,OP=(/"=="/),DESCRIPTOR=local_desc)
 !
 call io_variable_elemental(ID,&
&     VAR="Field direction",R1=X%q0,CHECK=.true.,OP=(/"==","==","=="/),DESCRIPTOR=local_desc)
 !
 if (ver_is_gt_or_eq(ID,(/3,0,1/))) then
   call io_variable_elemental(ID,&
&       VAR="BZ energy Double Grid",L0=use_X_DbGd,WARN=.true.,OP=(/"=="/),DESCRIPTOR=local_desc)
   call io_variable_elemental(ID,&
&       VAR="BZ energy DbGd points",I0=X_DbGd_nkpts,WARN=.true.,OP=(/"=="/),DESCRIPTOR=local_desc)
   if (ver_is_gt_or_eq(ID,revision=15206)) then
     call io_variable_elemental(ID,&
&         VAR="BZ Q point size factor  ",R0=X_DbGd_percentual,WARN=.true.,OP=(/"=="/),DESCRIPTOR=local_desc)
   endif
 endif
 !
 if (ver_is_gt_or_eq(ID,revision=21155)) then
   call def_variable_elemental(ID,trim(X_label)//"OPTICAL_AVERAGE",1,CR_YIO,0)
   call io_variable_elemental(ID,VAR="X Optical Averaged directions",CH0=OPTICAL_averaged_dirs_string,&
&                                CHECK=.true.,DESCRIPTOR=local_desc)
 endif
 !
 variables_X=io_status(ID)
 !
end function
