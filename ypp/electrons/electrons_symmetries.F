! 
! License-Identifier: GPL
!
! Copyright (C) 2023 The Yambo Team
!
! Authors (see AUTHORS file for details): DS
!
subroutine electrons_symmetries(Xk,Xen)
 !
 use electrons,      ONLY:levels
 use R_lattice,      ONLY:bz_samp
 !
 implicit none
 !
 type(levels) :: Xen
 type(bz_samp):: Xk
 !
 call WF_phase_matrices(Xen,Xk,(/1,Xen%nb/))
 !
end subroutine electrons_symmetries
