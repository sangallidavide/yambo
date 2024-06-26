!
! License-Identifier: GPL
!
! Copyright (C) 2022 The Yambo Team
!
! Authors (see AUTHORS file for details): HM DS
!
integer function io_BS_header(iq,X,ID,mode)
 !
 ! Intialize variables on BSE netcdf parallel file
 !
 use parallel_int,  ONLY:PP_wait
 use com,           ONLY:msg
 use pars,          ONLY:schlen,lchlen
 use stderr,        ONLY:intc
 use X_m,           ONLY:X_t
 use parallel_m,    ONLY:master_cpu
 use BS,            ONLY:BS_K_cutoff_done,l_BSE_restart,l_BSE_kernel_complete
 use BS_solvers,    ONLY:variables_BS,BSS_desc
 use IO_m,          ONLY:frag_BS_K,write_is_on,read_is_on,io_mode,io_BS_K,DUMP
 use IO_int,        ONLY:get_name,io_connect,io_disconnect,io_elemental,io_header
 use descriptors,   ONLY:IO_desc,IO_desc_duplicate
 !
 implicit none
 !
 type(X_t) :: X 
 integer   :: iq, ID
 character(*), intent(in) :: mode
 !
 logical, external :: file_exists
 integer, external :: variables_BS_more
 !
 type(IO_desc)     :: local_desc
 character(schlen) :: db_name,tmp_ch
 character(lchlen) :: full_db_name
 !
 ! Connect to parallel netCDF file.
 !
 if(.not.trim(mode)=="connect") then
   !
   db_name='BS_PAR_Q'//trim(intc(iq))
   full_db_name=get_name(trim(db_name)," ",2,CORE_IO=.FALSE.,MORE_IO=.TRUE.,COM_IO=.FALSE.)
   !
   l_BSE_kernel_complete=.false.
   !
   ! If l_BSE_restart true at this point, I'm already re-starting from a compressed kernel
   !
   if (write_is_on(ID).and.io_BS_K.and..not.l_BSE_restart.and.io_mode(ID)/=DUMP) then
     if (.not.file_exists(trim(full_db_name))) then
       l_BSE_restart=.false.
     else if(trim(mode)=="full") then
       call msg("s","Partial BSE file found in "//trim(full_db_name)//" used as restart.")
       l_BSE_restart = .true.
     endif
   endif
   !
   call PP_wait()
   !
 endif
 !
 if(.not. (read_is_on(ID) .or. (write_is_on(ID).and.master_cpu)) ) return
 !
 db_name='BS_head_Q'//trim(intc(iq))
 !
 io_BS_header=io_connect(desc=trim(db_name),type=2,ID=ID)
 if (io_BS_header/=0)       return
 if (trim(mode)=="connect") return
 !
 ! Check if the file already contains data and it is consistent
 !
 io_BS_header=io_header(ID,QPTS=.true.,R_LATT=.true.,WF=.true.,&
&                       IMPOSE_SN=.true.,XC_KIND="K_WF force Xs",CUTOFF=.true.,FRAG=frag_BS_K)
 !
 io_BS_header=variables_BS(ID,iq,local_desc,.TRUE.,X=X)
 !
 if (io_BS_header/=0.and.io_mode(ID)/=DUMP) then
   call io_elemental(ID,VAR="",VAR_SZ=0,MENU=1)
   call warning('Invalid dimensions found. BSE restart file '//trim(full_db_name)//' not accepted')
   return
 endif
 !
 io_BS_header=variables_BS_more(ID,'start')
 !
 if(read_is_on(ID).and.io_BS_K.and.io_mode(ID)/=DUMP) then
   if(l_BSE_kernel_complete) then
     call msg("s","Trying to load complete BSE file  ",trim(full_db_name))
   else
     call msg("s","Possible partial BSE file not accepted in reading mode ",trim(full_db_name))
   endif
 endif
 !
 if(read_is_on(ID).and.allocated(BS_K_cutoff_done).and.l_BSE_kernel_complete.and.io_mode(ID)/=DUMP) then
   if (     io_BS_K) tmp_ch="Backup option will be compressed BSE file"
   if (.not.io_BS_K) tmp_ch="Trying to load compressed BSE file"
   call msg("s",trim(tmp_ch))
   call msg("s","Previously computed compressed files ",BS_K_cutoff_done)
 endif
 !
 call io_disconnect(ID)
 !
 if (read_is_on(ID)) call PP_wait()
 !
 ! Descriptions transfer
 !   
 call IO_desc_duplicate(local_desc,BSS_desc)
 !
end function io_BS_header
