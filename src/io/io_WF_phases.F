!
! License-Identifier: GPL
!
! Copyright (C) 2011 The Yambo Team
!
! Authors (see AUTHORS file for details): DS
!
integer function io_WF_phases(nb,ikbz,i_sp_pol,ID,WF_phases)
 !
 use R_lattice,     ONLY:nXkbz
 use D_lattice,     ONLY:nsym
 use electrons,     ONLY:n_sp_pol
 use pars,          ONLY:schlen,SP,SP_YIO,IP_YIO
 use IO_m,          ONLY:io_sec,io_extension
 use IO_int,        ONLY:io_connect,io_disconnect,def_variable_bulk,&
&                        io_variable_bulk,io_header,def_variable_elemental,io_variable_elemental
 use wave_func,     ONLY:WF_phases_b_map,WF_phases_s_map
 use timing_m,      ONLY:timing
 !
 implicit none
 !
 integer                :: nb(2),ikbz,i_sp_pol,ID
 complex(SP), optional  :: WF_phases(3,3,nb(2)-nb(1)+1,nsym,1,1)
 !
 ! Work Space
 !
 integer :: VAR_ID,VAR_SIZE(7),VAR_POS(7)
 character(schlen)    :: VAR_name
 !
 call timing('io_WF_phases',OPR='start')
 !
 io_extension(ID)='wf_phases'
 !
 io_WF_phases=io_connect(desc=trim(io_extension(ID)),type=2,ID=ID)
 if (io_WF_phases/=0) goto 1
 !
 if (io_sec(ID,1)==0) goto 1
 !
 if (any((/io_sec(ID,:)==1/))) then
   !
   io_WF_phases=io_header(ID,R_LATT=.true.,WF=.true.,IMPOSE_SN=.true.,TEMP=.true.)
   if (io_WF_phases/=0) goto 1
   !
   call def_variable_elemental(ID,"nbands",2,IP_YIO,0)
   call io_variable_elemental(ID,VAR="Number of bands",I1=nb,CHECK=.true.,OP=(/"==","=="/))
   !call io_variable_elemental(ID,VAR="Number of bands",I1=nb,CHECK=.true.,OP=(/">=","<="/))
   !
   VAR_ID=1
   VAR_name='WF_phases'
   VAR_SIZE=(/2,3,3,nb(2)-nb(1)+1,nsym,nXkbz,n_sp_pol/)
   call def_variable_bulk(ID,trim(VAR_name),VAR_ID,VAR_SIZE,SP_YIO)
   !
 endif
 !
 if (any((/io_sec(ID,:)==2/))) then
   !
   ! WF_phases
   !
   VAR_ID=1
   VAR_name='WF_phases'
   VAR_POS=(/1,1,1,1,1,ikbz,i_sp_pol/)
   call io_variable_bulk(ID,VAR_ID,C6=WF_phases,IPOS=VAR_POS)
   !
 endif
 !
 if (any((/io_sec(ID,:)==2/))) then
   !
   ! WF_phases_k_map
   !
   VAR_ID=2
   VAR_name='WF_phases_s_map'
   call def_variable_bulk(ID,trim(VAR_name),VAR_ID,shape(WF_phases_s_map),IP_YIO)
   call  io_variable_bulk(ID,VAR_ID,I2=WF_phases_s_map)
   !
   ! WF_phases_b_map
   !
   VAR_ID=3
   VAR_name='WF_phases_b_map'
   call def_variable_bulk(ID,trim(VAR_name),VAR_ID,shape(WF_phases_b_map),IP_YIO)
   call  io_variable_bulk(ID,VAR_ID,I4=WF_phases_b_map)
   !
 endif
 !
1 call io_disconnect(ID=ID)
 !
 call timing('io_WF_phases',OPR='stop')
 !
end function
