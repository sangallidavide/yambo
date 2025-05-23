!
! License-Identifier: GPL
!
! Copyright (C) 2021 The Yambo Team
!
! Authors (see AUTHORS file for details): AM
!
subroutine X_ALLOC_global(WHAT,W)
 !
 use pars,          ONLY:cZERO
 use ALLOC,         ONLY:X_ALLOC_elemental
 use X_m,           ONLY:Resp_ii,Resp_ij,Epsilon_ii,Epsilon_ij,&
&                        Alpha,Joint_DOS,X_magnons,X_dichroism,EEL,X_do_obs
 use frequency,     ONLY:w_samp
 use drivers,       ONLY:l_optics
 use BS_solvers,    ONLY:B_Hall
 use PHOTOLUM,      ONLY:PL,PL_init
 !
#include<y_memory.h>
 !
 character(*)   :: WHAT
 type(w_samp)   :: W
 !
 if (WHAT=="ALLOC") then
   !
   if (l_optics) then
     call X_ALLOC_elemental("RESP",(/W%n_freqs,4/))
     Resp_ii(:,1)=W%p(:)
   endif
   !
   if (X_do_obs("eps")) then
     call X_ALLOC_elemental("EPS",(/W%n_freqs,4/))
     Epsilon_ii(:,1)=W%p(:)
   endif
   !
   if (X_do_obs("alpha")) then
     call X_ALLOC_elemental("ALPHA",(/W%n_freqs,2/))
     Alpha(:,1)=W%p(:)
   endif
   !
   if (X_do_obs("kerr")) then
     call X_ALLOC_elemental("OFF",(/W%n_freqs,4/))
     Resp_ij(:,1)=W%p(:)
     B_Hall=cZERO
     if (X_do_obs("eps")) then
       Epsilon_ij(:,1)=W%p(:)
     endif
   endif
   !
   if (X_do_obs("magnons")) then
     call X_ALLOC_elemental("MAGNONS",(/W%n_freqs,2,4/))
     X_magnons=cZERO
     X_magnons(:,1,1)=W%p(:)
     X_magnons(:,2,1)=W%p(:)
   endif
   !
   if (X_do_obs("CD")) then
     call X_ALLOC_elemental("DICH",(/W%n_freqs,4/))
     X_dichroism(:,1)=W%p(:)
   endif
   !
   if (X_do_obs("pl")) then
     call PL_init( )
     call X_ALLOC_elemental("PL",(/W%n_freqs,3/))
     PL(:,1) = W%p(:)  
   endif
   !
   if (X_do_obs("eel")) call X_ALLOC_elemental("EEL",(/W%n_freqs,2/))
   !
   if (X_do_obs("jdos")) then
     call X_ALLOC_elemental("JDOS",(/W%n_freqs,4/))
     Joint_DOS(:,1)=W%p(:)
   endif
   !
   if (X_do_obs("Esort")) call X_ALLOC_elemental("BS_E_SORTED",(/1000,9,2/))
   if (X_do_obs("fxc")) call X_ALLOC_elemental("FXC",(/W%n_freqs/))
   !
 else if (WHAT=="FREE") then
   call X_ALLOC_elemental("FREE ALL")
 endif
 !
end subroutine
