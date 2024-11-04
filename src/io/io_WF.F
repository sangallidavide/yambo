!
! License-Identifier: GPL
!
! Copyright (C) 2015 The Yambo Team
!
! Authors (see AUTHORS file for details): AM DS
!
integer function io_WF(ID,wf,ipath,keep_fragments_open)
 !
 use pars,             ONLY:SP,schlen
 use electrons,        ONLY:n_spinor
 use stderr,           ONLY:intc
 use R_lattice,        ONLY:nkibz
 use wave_func,        ONLY:wf_ncx,wf_nb_io,wf_nb_io_groups,wf_frag_id
 use IO_m,             ONLY:io_sec,io_resetable,frag_WF,NETCDF_file_is_open
 use IO_int,           ONLY:io_connect,io_disconnect,io_header,&
&                           io_fragment,io_elemental,io_variable_bulk,def_variable_bulk
 use timing_m,         ONLY:timing
 !use devxlib
 !
 implicit none
 integer,       intent(in) :: ID
 integer,       optional, intent(in) :: ipath
 logical,       optional, intent(in) :: keep_fragments_open
 complex(SP),   optional   :: wf(:,:,:)
 !
 logical                   :: keep_frag_open,skip_var_def
 integer                   :: i_sp_pol,ik,ifrag,ib_grp,ib_block,n_blocks,ib_lim(3,wf_nb_io),head_io,ID_frag,ipath_
 character(schlen)         :: VAR_NAME
 !
 call timing('io_WF',OPR='start')
 !
 keep_frag_open = .false.
 if (present(keep_fragments_open)) keep_frag_open = keep_fragments_open
 !
 ipath_=1
 if (present(ipath)) ipath_=ipath
 !
 io_WF=0
 head_io=0
 select case ( count((/io_sec(ID,:)/=0/)) )
   case(0)
     io_resetable(ID)=.true.
     goto 1
   case(1)
     ifrag=0
     ib_grp=0
     head_io=1
   case(2)
     ifrag=io_sec(ID,1)-1
     ib_grp=io_sec(ID,2)
     if (.not. present(wf)) call error(" [WF] unexpected: wf not present (1)")
   case(3)
     ifrag=io_sec(ID,2)-1
     ib_grp=io_sec(ID,3)
     head_io=1
     if (.not. present(wf)) call error(" [WF] unexpected: wf not present (2)")
 end select
 !
 if(ifrag<=nkibz) then
   ik=ifrag
   i_sp_pol=1
 else
   ik=ifrag-nkibz
   i_sp_pol=2
 endif
 !
 io_WF=io_connect(desc="wf",type=0,ID=ID,FORCE_READ_MODE=(head_io/=1).and.frag_WF)
 if (io_WF/=0) goto 1
 !
 ! S/N
 !
 if (head_io==1) then
   !
   io_WF=io_header(ID,IMPOSE_SN=.true.,FRAG=frag_WF)
   !
   call io_elemental(ID, VAR="BAND_GROUPS",VAR_SZ=2,MENU=0)
   call io_elemental(ID, VAR=" Bands in each block              ",I0=wf_nb_io)
   call io_elemental(ID, VAR=" Blocks                           ",I0=wf_nb_io_groups)
   call io_elemental(ID, VAR="",VAR_SZ=0)
   !
   if (io_WF/=0) goto 1
   !
 endif
 !
 if (ifrag==0) goto 1
 !
 ! Define number of bocks to read inside single fragment
 !
 call compute_ib_lims(ib_lim,n_blocks)
 if (n_blocks==0) goto 1
 !
 ! Fragmentation
 !
 skip_var_def=.false.
 if (keep_frag_open) then 
   if (wf_frag_id(ifrag,ib_grp,ipath)/=0) then 
     skip_var_def=.true.
     ID_frag=wf_frag_id(ifrag,ib_grp,ipath)
    endif
 endif
 !
 if (.not.skip_var_def) then
   !
   call io_fragment(ID,ID_frag,i_fragment=ifrag,j_fragment=ib_grp,keep_open=keep_frag_open)
   if (keep_frag_open) wf_frag_id(ifrag,ib_grp,ipath)=ID_frag
   !
   ! Define variables
   !
   VAR_NAME="WF_COMPONENTS_@_SP_POL"//trim(intc(i_sp_pol))//"_K"//trim(intc(ik))//"_BAND_GRP_"//trim(intc(ib_grp))
   call def_variable_bulk(ID_frag,VAR=trim(VAR_NAME),VAR_ID=1,VAR_SZ=(/2,wf_ncx,n_spinor,wf_nb_io/),VAR_KIND=SP)
   !
 endif
 !
 do ib_block=1,n_blocks
   if(ib_grp>1) ib_lim(1,ib_block)=mod(ib_lim(1,ib_block)-1,wf_nb_io)+1
   call io_variable_bulk(ID_frag,1,C3=wf(:,:,ib_lim(2,ib_block):ib_lim(3,ib_block)),IPOS=(/1,1,1,ib_lim(1,ib_block)/))
 enddo
 !
 call io_fragment_disconnect(ID,ID_frag)
 !
1 call io_disconnect(ID=ID)
 !
 call timing('io_WF',OPR='stop')
 !
end function io_WF
!
!
integer function io_WF_frag(ID_frag,wf)
 !
 use pars,             ONLY:SP,schlen
 use electrons,        ONLY:n_spinor
 use stderr,           ONLY:intc
 use R_lattice,        ONLY:nkibz
 use wave_func,        ONLY:wf_ncx,wf_nb_io,wf_nb_io_groups
 use IO_m,             ONLY:io_sec,NETCDF_file_is_open
 use IO_int,           ONLY:def_variable_bulk,io_variable_bulk,io_connect
 use timing_m,         ONLY:timing
 implicit none
 integer,       intent(in) :: ID_frag
 complex(SP),   optional   :: wf(:,:,:)
 !
 integer                   :: i_sp_pol,ik,ifrag,ib_grp,ib_block,n_blocks,ib_lim(3,wf_nb_io)
 character(schlen)         :: VAR_NAME
 !
 call timing('io_WF_frag',OPR='start')
 !
 io_WF_frag=0
 ifrag=io_sec(ID_frag,1)-1
 ib_grp=io_sec(ID_frag,2)
 !
 if(ifrag<=nkibz) then
   ik=ifrag
   i_sp_pol=1
 else
   ik=ifrag-nkibz
   i_sp_pol=2
 endif
 !
 ! Define number of bocks to read inside single fragment
 !
 call compute_ib_lims(ib_lim,n_blocks)
 if (n_blocks==0) return
 !
 if(.not.NETCDF_file_is_open(ID_frag)) then
   !
   ! File definition
   !
   io_WF_frag=io_connect(desc="wf_fragments_"//trim(intc(ifrag))//"_"//trim(intc(ib_grp)),type=0,ID=ID_frag)
   if (io_WF_frag/=0) return
   !
   VAR_NAME="WF_COMPONENTS_@_SP_POL"//trim(intc(i_sp_pol))//"_K"//trim(intc(ik))//"_BAND_GRP_"//trim(intc(ib_grp))
   call def_variable_bulk(ID_frag,VAR=trim(VAR_NAME),VAR_ID=1,VAR_SZ=(/2,wf_ncx,n_spinor,wf_nb_io/),VAR_KIND=SP)
   !
 endif
 !
 do ib_block=1,n_blocks
   if(ib_grp>1) ib_lim(1,ib_block)=mod(ib_lim(1,ib_block)-1,wf_nb_io)+1
   call io_variable_bulk(ID_frag,1,C3=wf(:,:,ib_lim(2,ib_block):ib_lim(3,ib_block)),IPOS=(/1,1,1,ib_lim(1,ib_block)/))
 enddo
 !
 call timing('io_WF_frag',OPR='stop')
 !
end function io_WF_frag
!
!
subroutine compute_ib_lims(ib_lim,n_blocks)
 !
 use wave_func,          ONLY:wf_nb_io,wf_b_indx,wf_nb_to_load
 !
 implicit none
 !
 integer, intent(inout) :: n_blocks
 integer, intent(out)   :: ib_lim(3,wf_nb_io)
 !
 integer :: ib_block,nb_io,ib_prev,ib,ib_io
 !
 if(.not.allocated(wf_b_indx)) then
  ib=1
  nb_io=wf_nb_io
 else
  ib=wf_b_indx(1)
  nb_io=wf_nb_to_load
 endif
 !
 if(nb_io==0) then
  n_blocks=0
  return
 endif
 !
 n_blocks=1
 ib_prev=ib
 !
 ib_lim(1,n_blocks)=ib
 ib_lim(2,n_blocks)=1
 !
 if(nb_io==1)        ib_lim(3,n_blocks)=1
 !
 if(nb_io==wf_nb_io) ib_lim(3,n_blocks)=nb_io
 !
 if(.not.(nb_io/=1 .and. nb_io/=wf_nb_io)) return
 !
 do ib_io=2,nb_io
   !
   if(.not.allocated(wf_b_indx)) ib=ib_io
   if(     allocated(wf_b_indx)) ib=wf_b_indx(ib_io)
   !
   if(ib_prev+1==ib) then
     ib_prev=ib
     cycle
   endif
   !
   ib_lim(3,n_blocks)=ib_io-1
   !
   n_blocks=n_blocks+1
   !
   ib_prev=ib
   ib_lim(1,n_blocks)=ib
   ib_lim(2,n_blocks)=ib_io
   !
 enddo
 !
 ib_lim(3,n_blocks)=ib_io-1
 !
end subroutine compute_ib_lims
 
