!
! License-Identifier: GPL
!
! Copyright (C) 2020 The Yambo Team
!
! Authors (see AUTHORS file for details): AM DS AF
!
subroutine K_restart(iq,X,ID_head,ID,ID_compr,l_partial_kernel_loaded)
 !
 use pars,         ONLY:SP,schlen
 use LIVE_t,       ONLY:live_timing
 use com,          ONLY:msg
 use X_m,          ONLY:X_t
 use stderr,       ONLY:real2ch
 use units,        ONLY:Gigabyte
 use drivers,      ONLY:l_bs_fxc,l_bss
 use BS,           ONLY:l_BSE_restart,l_BSE_kernel_complete,&
&                       BS_K_cutoff,BS_K_cutoff_done,n_BS_blks,BS_blk
 use BS_solvers,   ONLY:BSS_kernel_IO_on_the_fly,BSS_kernel_last_stored_block,BSS_kernel_stored_size
 use IO_int,       ONLY:io_control
 use IO_m,         ONLY:OP_RD_CL,REP,VERIFY,NONE,OP_RD,manage_action,&
&                       RD,WR,RD_CL,RD_CL_IF_END,WR_CL,OP_APP,DUMP
 !
 implicit none
 !
 type(X_t), intent(inout) :: X
 integer, intent(in)      :: iq
 integer, intent(out)     :: ID_head,ID,ID_compr
 logical, intent(out)     :: l_partial_kernel_loaded
 !
 ! Work space
 !
 integer           :: i_file,i_block,io_BS_err,IO_ACT,ID_now
 real(SP)          :: BS_K_cutoff_inp,size_now,Sz
 character(schlen) :: timing_string,IO_mode
 !
 integer, external :: io_BS,io_BS_PAR_init,io_BS_header
 real(SP),external :: best_cutoff_allocated
 !
 ! BS/FXC DB I/O
 !===============
 if (l_bs_fxc) then
   !
   !call io_control(ACTION=OP_RD_CL,COM=REP,MODE=VERIFY,SEC=(/1/),ID=ID)
   !
   !io_BS_Fxc_err=io_BS_Fxc(iq,W_bss,ID,X=X)
   !
   return
   !
 endif
 !
#if defined _PAR_IO
 !
 BS_K_cutoff_inp=BS_K_cutoff
 l_BSE_restart=.false.
 l_partial_kernel_loaded=.false.
 IO_mode="full"   ! Full loads BS_blk%done from header & BS_blk%mat from BS_PAR 
 IO_ACT=OP_RD_CL
 !
 call io_control(ACTION=OP_RD_CL,COM=REP,MODE=VERIFY,ID=ID_head)
 io_BS_err=io_BS_header(iq,X,ID_head,IO_mode)
 if(io_BS_err/=0) return
 !
 call io_control(ACTION=OP_RD_CL,COM=REP,MODE=DUMP,ID=ID_head)
 io_BS_err=io_BS_header(iq,X,ID_head,IO_mode)
 !
 BS_K_cutoff=best_cutoff_allocated(BS_K_cutoff_inp)
 !
 ! Check if ndb.BS_PAR exist.
 ! In case use it
 call io_control(ACTION=OP_RD,COM=REP,MODE=DUMP,ID=ID)
 io_BS_err=io_BS_PAR_init(iq,ID,IO_mode)
 !
 if(io_BS_err/=0.and.allocated(BS_K_cutoff_done)) then
   IO_mode="compressed_head"
   ! Check for ndb.BS_PAR compressed from the list of previous calculations
   do i_file=1,size(BS_K_cutoff_done)
     call io_control(ACTION=OP_RD,COM=REP,MODE=DUMP,ID=ID_compr)
     io_BS_err=io_BS_PAR_init(iq,ID_compr,IO_mode)
     if(io_BS_err==0) exit
     BS_K_cutoff=best_cutoff_allocated(BS_K_cutoff_inp)
   enddo
 endif
 !
 if(l_BSE_kernel_complete.and.l_bss.and.io_BS_err==0) then
   if (BSS_kernel_stored_size/=0._SP) then
   if(trim(IO_mode)=="full") then
     timing_string="Loading full BSE kernel"
     ID_now=ID
   endif
   if(trim(IO_mode)=="compressed_head") then
     timing_string="Loading BSE_done from compressed file"
     ID_now=ID_compr
   endif
   if(n_BS_blks>0) call live_timing(trim(timing_string),n_BS_blks)
   if(.not.BSS_kernel_IO_on_the_fly) BSS_kernel_last_stored_block=n_BS_blks
   size_now=0
   do i_block=1,n_BS_blks
     call io_control(ACTION=RD,ID=ID_now)
     call io_BS_PAR_block(iq,i_block,ID_now,IO_MODE)
     call live_timing(steps=1)
     Sz=BS_blk(i_block)%size(1)*BS_blk(i_block)%size(2)
     size_now=size_now+Sz*kind(size_now)*2/Gigabyte
     if (BSS_kernel_IO_on_the_fly.and.size_now>=BSS_kernel_stored_size) then
       BSS_kernel_last_stored_block=i_block
       exit
     endif
   enddo
   if(n_BS_blks>0) call live_timing()
   if (i_block<n_BS_blks) call msg("s","BSE kernel stored in memory for "//trim(real2ch(size_now))//" GB")
   if (i_block>n_BS_blks) BSS_kernel_last_stored_block=n_BS_blks
   endif
   if (BSS_kernel_last_stored_block==n_BS_blks) BSS_kernel_IO_on_the_fly=.false.
   if (BSS_kernel_IO_on_the_fly.and.trim(IO_mode)=="full") then
     call msg("s","BSE kernel will be loaded on the fly")
   endif
   !
   ! Only BS_blk%done was loaded
   ! Now I need to estimate compressed size and load
   ! also BS_blk%mat from compressed file
   !
   if(trim(IO_mode)=="compressed_head") then
     call K_compressed_IO(iq,io_BS_err,"read")
     if (BS_K_cutoff_inp<BS_K_cutoff.and.io_BS_err==0) then
       call msg("s","Found kernel more compressed than needed. Switching to restart mode")
       l_BSE_kernel_complete=.false.
       l_BSE_restart=.true.
       l_partial_kernel_loaded=.true.
     endif
     if (io_BS_err/=0) then
       call msg("s","Compressed kernel file missing.")
       l_BSE_kernel_complete=.false.
       l_BSE_restart=.false.
     endif
   endif
   !
   ! Finally I check if I need to write the compressed file
   !
   if(BS_K_cutoff_inp>0._SP.and.l_BSE_kernel_complete) then
     ! If my compression level now is different from existing one
     if(BS_K_cutoff_inp/=BS_K_cutoff) then
       call msg("s","Performing extra compression requested from input.")
       BS_K_cutoff=BS_K_cutoff_inp
       ! Perform compression
       call K_compress(iq,ID_compr)
       call K_compressed_IO(iq,io_BS_err,"write")
       ! Open HEAD in write mode
       call io_control(ACTION=OP_APP,COM=REP,ID=ID_head)
       io_BS_err=io_BS_header(iq,X,ID_head,"connect")
       IO_ACT=WR_CL
     endif
   endif
   !
 endif
 !
 call io_control(ACTION=IO_ACT,COM=REP,ID=ID_head)
 if(.not.BSS_kernel_IO_on_the_fly) call io_control(ACTION=RD_CL,COM=REP,ID=ID)
 if(     l_BSE_kernel_complete) call io_BS_PAR_free(ID_head,ID,IO_mode,.true.)
 if(.not.l_BSE_kernel_complete) call io_BS_PAR_free(ID_head,ID,IO_mode,.false.)
 !
 BS_K_cutoff=BS_K_cutoff_inp
 !
#else
 !
 ID_head=0
 ID_compr=0
 l_BSE_kernel_complete=.false.
 l_partial_kernel_loaded=.false.
 l_BSE_kernel_complete=.false.
 !
 call io_control(ACTION=OP_RD_CL,COM=REP,MODE=VERIFY,SEC=(/1/),ID=ID)
 io_BS_err=io_BS(iq,X,ID)
 !
 if (io_BS_err==0.and.l_bss) then
   if(n_BS_blks>0) call live_timing("Loading BSE kernel",n_BS_blks)
   call io_control(ACTION=OP_RD,COM=NONE,MODE=DUMP,SEC=(/1/),ID=ID)
   io_BS_err=io_BS(iq,X,ID)
   do i_block=1,n_BS_blks
     IO_ACT=manage_action(RD_CL_IF_END,i_block,FIRST=1,LAST=n_BS_blks)
     call io_control(ACTION=IO_ACT,COM=NONE,MODE=DUMP,SEC=(/i_block+1/),ID=ID)
     io_BS_err=io_BS(iq,X,ID)
     if(n_BS_blks>0) call live_timing(steps=1)
   enddo
   if(n_BS_blks>0) call live_timing()
 endif
 !
 if (io_BS_err==0) l_BSE_kernel_complete=.true.
 !
#endif
 !
end subroutine K_restart
!
!
function best_cutoff_allocated(BS_K_cutoff_inp)
 !
 use vec_operate,  ONLY:sort
 use pars,         ONLY:SP
 use BS,           ONLY:BS_K_cutoff_done,BS_K_cutoff
 !
 implicit none
 !
 real(SP)             :: best_cutoff_allocated
 real(SP), intent(in) :: BS_K_cutoff_inp
 !
 integer :: i1
 !
 best_cutoff_allocated=BS_K_cutoff_inp
 !
 if(.not.allocated(BS_K_cutoff_done)) return
 !
 call sort(BS_K_cutoff_done)
 !
 do i1=size(BS_K_cutoff_done),1,-1
   if(BS_K_cutoff_done(i1)>BS_K_cutoff) cycle
   exit
 enddo
 !
 if(i1==0) i1=minloc(BS_K_cutoff_done,1)
 !
 best_cutoff_allocated=BS_K_cutoff_done(i1)
 BS_K_cutoff_done(i1)=1
 !
end function best_cutoff_allocated
