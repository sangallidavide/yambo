! 
! License-Identifier: GPL
!
! Copyright (C) 2023 The Yambo Team
!
! Authors (see AUTHORS file for details): DS
!
subroutine excitons_symmetries(Xk,BS_R_right,BS_E,BS_E_degs)
 !
 use pars,          ONLY:SP
 use R_lattice,     ONLY:bz_samp
 use BS_solvers,    ONLY:BS_mat,BSS_n_eig
 use BS,            ONLY:BS_bands
 !
 implicit none
 !
 type(bz_samp) :: Xk
 integer       :: BS_E_degs(BSS_n_eig)
 complex(SP)   :: BS_R_right(BSS_n_eig),BS_E(BSS_n_eig)
 !
 !
end subroutine excitons_symmetries
