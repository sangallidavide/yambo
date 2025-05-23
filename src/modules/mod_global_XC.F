!
! License-Identifier: GPL
!
! Copyright (C) 2008 The Yambo Team
!
! Authors (see AUTHORS file for details): AM
!
module global_XC
 !
 use pars,      ONLY:lchlen,SP,RZERO
 !
 implicit none
 !
 ! general notice:
 !
 ! Yambo uses integers to store locally the xc specs. These integers
 ! are transformed in strings when setup_global_XC is called.
 ! setup_global_XC is called in
 !
 !  mod_QP_ctl: when using external QP corrections
 !  setup: when using SC energies/WF
 !
 ! The strings obtained are then written in the corresponding DB or veryfied
 ! against the value in the DB. Note that in the case of WFs only the loaded_WF_xc_string
 ! is set in SC_initial_WF_setup. All other strings are eventually read from the DB's
 !
 ! Energies (linked to the three energy types) ...
 !
 !  ... Integers ...
 !
 integer     , public :: G_kind,G_xc_functional,G_perturbation
 integer     , public :: X_kind,X_xc_functional,X_perturbation
 integer     , public :: K_kind,K_xc_functional,K_perturbation
 !
 ! ... and Strings... 
 !
 character(2*lchlen) , public :: G_E_xc_string
 character(2*lchlen) , public :: X_E_xc_string(4)
 character(2*lchlen) , public :: K_E_xc_string
 !
 ! Wavefunctions ...
 !
 ! ... integers ...
 !
 integer     , public :: WF_kind,WF_xc_functional,WF_perturbation
 real(SP)    , public :: WF_exx_FRACTION, WF_exx_SCREENING
 !
 character(lchlen) , public :: loaded_WF_xc_string  
 !
 ! ... and Strings
 ! ... these are read from the DB (when writing, instead, all these string are set equal to  loaded_WF_xc_string)
 !
 character(lchlen) , public :: G_WF_xc_string
 character(lchlen) , public :: X_WF_xc_string(4)
 character(lchlen) , public :: K_WF_xc_string
 character(lchlen) , public :: Dipole_WF_xc_string ! Dipole strengths DB uses its own field (so no integer associated)
 !
 ! QP database
 !
 integer     , public :: QP_DB_kind
 !
 ! Self Energies KINDS
 !
 ! Only diagonal term in the hamiltonian basis-set
 integer, public, parameter :: &
   QP_SE_NONE          = 400, &
   QP_SE_COHSEX        = 401, &
   QP_SE_GoWo_PPA      = 402, &
   QP_SE_GoWo          = 403, &
   QP_SE_GWo_PPA       = 404, &
   QP_SE_GWo           = 405, &
   QP_SE_POLARON       = 406, &
   QP_SE_COH           = 409, &
   QP_SE_SEX           = 410, & 
   QP_SE_GW_PPA        = 411
 !
 ! Perturbations
 integer, public, parameter :: &
   EXT_NONE             = 601, &
   EXT_B_LANDAU         = 602, &
   EXT_B_PAULI          = 603, &
   EXT_B_MAGNETISM      = 604
 !
 ! Scattering Self-Energies (to be renamed)
 integer, public, parameter :: &
   QP_SE_GW_NEQ        = 407, &
   QP_SE_GW_QED        = 408
 !
 ! Full matrix in the hamiltonian basis-set (SC or RT)
 integer, public, parameter :: &
   H_SE_Hartree       =1000, &  ! Hartree
   H_SE_EXX           = 301, &  ! OEP-EXX
   H_SE_Fock          = 302, &  ! Fock
   H_SE_EXXC          = 304, &  ! OEP-EXX + LDA correlation
   H_SE_SRPA          = 305, &  ! OEP-COHSEX (static RPA)
   H_SE_EXX_KLI       = 311, &  ! OEP-EXX (KLI apprx)
   H_SE_EXXC_KLI      = 314, &  ! OEP-EXX (KLI apprx) + LDA correlation
   H_SE_SRPA_KLI      = 315, &  ! OEP-COHSEX (KLI apprx)
   H_SE_EXX_SLT       = 321, &  ! OEP-EXX (SLT apprx)
   H_SE_EXXC_SLT      = 324, &  ! OEP-EXX (SLT apprx) + LDA correlation
   H_SE_SRPA_SLT      = 325, &  ! OEP-COHSEX (SLT apprx)
   H_SE_COH           = 326, &  ! Coulomb-Hole 
   H_SE_SEX           = 327, &  ! Screened-eXchange
   H_SE_LSEX          = 328, &  ! Long range-screened-eXchange
   H_SE_LHF           = 329, &  ! Long range-eXchange
   H_SE_COHSEX        = H_SE_SEX+H_SE_COH  ! Coulomb-Hole + Screened-eXchange
 !
 ! Phonons
 integer, public, parameter :: &
   QPH_X_RPA         = 501    ! Simple RPA PH self-energy
 !
 contains
   !
   subroutine  setup_global_XC(what,kind,xc_functional,perturbation,exx_fraction,exx_screening)
     character(*)           :: what          
     integer                :: kind          
     integer,      optional :: xc_functional 
     integer,      optional :: perturbation 
     real(SP),     optional :: exx_fraction
     real(SP),     optional :: exx_screening
     !
     ! Work Space
     !
     integer :: local_xc_functional,local_perturbation
     real(SP):: local_exx_fraction,local_exx_screening
     !
     local_xc_functional=0
     local_perturbation=EXT_NONE
     local_exx_fraction = RZERO
     local_exx_screening = RZERO
     if (present(xc_functional)) local_xc_functional=xc_functional
     if (present(perturbation))  local_perturbation=perturbation
     if (present(exx_fraction))  local_exx_fraction=exx_fraction
     if (present(exx_screening)) local_exx_screening=exx_screening
     !
     if (what=='E'.or.what=='preset'.or.what=='G') then
       G_kind=kind
       G_xc_functional=local_xc_functional
       G_perturbation=local_perturbation
       G_E_xc_string=global_XC_string(kind,local_xc_functional,local_perturbation)
     endif
     !
     if (what=='X'.or.what=='preset') then
       X_kind=kind
       X_xc_functional=local_xc_functional
       X_perturbation=local_perturbation
       X_E_xc_string=global_XC_string(kind,local_xc_functional,local_perturbation)
     endif
     !
     if (what=='K'.or.what=='preset') then
       K_kind=kind
       K_xc_functional=local_xc_functional
       K_perturbation=local_perturbation
       K_E_xc_string=global_XC_string(kind,local_xc_functional,local_perturbation)
     endif
     !
     if (what=='WF'.or.what=='preset') then
       WF_kind=kind
       WF_xc_functional=local_xc_functional
       WF_exx_FRACTION =local_exx_fraction
       WF_exx_SCREENING =local_exx_screening
       !
       ! loaded_WF_xc_string represents the xc kind of the WF that are loaded trough wf_load.
       ! G/X/K/Dipole_WF_xc_string, instead, may be overwritten by the value written in the databases.
       !
       loaded_WF_xc_string=global_XC_string(kind,local_xc_functional,local_perturbation)
       !
       if (what=='preset') G_WF_xc_string     =loaded_WF_xc_string
       if (what=='preset') X_WF_xc_string     =loaded_WF_xc_string
       if (what=='preset') K_WF_xc_string     =loaded_WF_xc_string
       if (what=='preset') Dipole_WF_xc_string=loaded_WF_xc_string
       !
     endif
     !
     if (what=='preset') QP_DB_kind=QP_SE_none
     !
   end subroutine
   !
   subroutine  correct_global_XC(what,STRING)
    character(*), intent(in) :: what
    character(*), intent(in) :: STRING
    !
    ! Work Space
    !
    integer :: i_X
    !
    if (what=='E'.or.what=='G')         G_E_xc_string=      trim(G_E_xc_string)//"+"//STRING
    if (what=='X') then
      forall(i_X=1:size(X_E_xc_string)) X_E_xc_string(i_X)= trim(X_E_xc_string(i_X))//"+"//STRING
    endif
    if (what=='K')                      K_E_xc_string=      trim(K_E_xc_string)//"+"//STRING
    if (what=='WF')               loaded_WF_xc_string=trim(loaded_WF_xc_string)//"+"//STRING
    !
   end subroutine
   !
   character(lchlen) function global_XC_string(kind,xc_functional,perturbation)
    !
    use xc_functionals,  ONLY:xc_string,XC_CORRELATION
    !
    integer           :: kind          ! <300 (LDA's GGA's ...) >300 See above
    integer           :: xc_functional ! In case kind < 300 gives DFT functional
    integer           :: perturbation
    ! 
    integer           :: local_kind
    character(1)      :: ch_plus
    character(lchlen) :: ch
    global_XC_string=' '
    do local_kind=kind-1000,kind,1000
      select case (local_kind)
        case(-1)
          global_XC_string='IP'
        case(H_SE_Hartree) !=1000
          global_XC_string='Hartree'
        case(1:300)
          global_XC_string=trim(xc_string(xc_functional))
        case(H_SE_EXX)
          global_XC_string='SC Exact Exchange'
        case(H_SE_EXXC)
          global_XC_string='SC Exact Exchange +'//trim(xc_string(xc_functional))
        case(H_SE_FOCK)
          global_XC_string='SC Fock'
        case(H_SE_COH)
          global_XC_string='SC COH'
        case(H_SE_SEX)
          global_XC_string='SC SEX'
        case(H_SE_COHSEX)
          global_XC_string='SC COHSEX'
        case(QP_SE_COH)
          global_XC_string='COH'
        case(QP_SE_SEX)
          global_XC_string='SEX'
        case(QP_SE_COHSEX)
          global_XC_string='COHSEX'
        case(QP_SE_GoWo_PPA)
          global_XC_string='GoWo Plasmon Pole approx.'
        case(QP_SE_GoWo)
          global_XC_string='GoWo on real axis'
        case(QP_SE_GWo_PPA)
          global_XC_string='GWo Plasmon Pole approx.'
        case(QP_SE_GWo)
          global_XC_string='GWo on real axis'
        case(QP_SE_POLARON)
          global_XC_string='Polarons'
        case(QP_SE_GW_NEQ)
          global_XC_string='OUT-of-EQ GW'
        case(QP_SE_GW_QED)
          global_XC_string='QED-based GW'
        case(QPH_X_RPA)
          global_XC_string='QPH RPA approx'
        case default
          if(trim(global_XC_string)/="".and.local_kind>1000) global_XC_string='SC HARTREE +'//trim(global_XC_string)
      end select
    enddo
    !
    ch_plus='+'
    if (kind==0.or.len_trim(global_XC_string)==0) ch_plus=' '
    !
    if (perturbation/=EXT_NONE) then
      ch=global_XC_string
      select case (perturbation)
        case(EXT_B_PAULI)
          global_XC_string=trim(ch)//trim(ch_plus)//'Pauli Magnetism'
        case(EXT_B_LANDAU)
          global_XC_string=trim(ch)//trim(ch_plus)//'Landau Magnetism'
        case(EXT_B_MAGNETISM)
          global_XC_string=trim(ch)//trim(ch_plus)//'Magnetism'
      end select
    endif
    !
   end function
   !
end module 
