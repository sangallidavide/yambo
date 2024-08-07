!
! License-Identifier: GPL
!
! Copyright (C) 2006 The Yambo Team
!
! Authors (see AUTHORS file for details): AM CA
!
module interpolate
 !
 use pars,       ONLY:SP,DP,schlen
 use R_lattice,  ONLY:bz_samp
 !
 ! Interpolation according to PRB 38 p2721
 ! Code copied from BolzTraP 
 ! Notice that this interpolation is the same interpolation
 ! implemented in AbiPy (or YamboPy)
 ! where is called Shankland-Koelling-Wood(SKW)
 !
 implicit none
 !
 integer               :: nshells
 real(SP), allocatable :: int_sop(:,:,:)
 real(SP), pointer     :: lattice_vectors(:,:)
 real(SP) :: INTERP_shell_factor
 real(SP) :: metric(3,3)
 !
 ! Roughness function parameters
 !
 real(DP), parameter :: R_C1=0.75_DP
 real(DP), parameter :: R_C2=0.75_DP
 !  
 character(schlen) :: INTERP_mode = "NN" ! NN/BOLTZ
 character(schlen) :: INTERPOLATION_name = " "
 logical           :: INTERPOLATE_initialized=.FALSE.
 logical           :: INTERPOLATE_is_quiet=.FALSE.
 logical           :: INTERPOLATE_is_serial=.FALSE.
 logical           :: INTERPOLATE_use_symm=.TRUE.
 character(schlen) :: INTERPOLATE_msg_fmt="rs"
 !
 integer, parameter :: max_interpls=40
 !
 ! NN method
 !
 integer :: OUT_k_nk=0
 integer :: NN_nk_exact=0
 integer :: NN_nk_border=0
 integer :: NN_n_of_nearest=1
 integer :: NN_MAX_real_nearest=0
 !
 ! What
 !
 logical :: l_INTERP_f
 logical :: l_INTERP_E
 logical :: l_INTERP_dE
 logical :: l_INTERP_Z
 logical :: l_INTERP_W
 logical :: l_integral_respect=.TRUE. ! In the case of f/df keep the same integral of the SEED and of the interpolated
 character(schlen) :: INTERP_FineGd   = "NONE"
 !
 ! ID's
 !
 integer :: INTERP_E_id =0
 integer :: INTERP_dE_id=0
 integer :: INTERP_f_id=0
 integer :: INTERP_RT_id=0
 !
 ! K-grids
 !
 type(bz_samp):: GRID_k
 !
 ! Interpolation type
 !
 character(schlen), parameter :: REAL1D   = "REAL1D"
 character(schlen), parameter :: REAL2D   = "REAL2D"
 character(schlen), parameter :: CMPLX1D  = "CMPLX1D"
 character(schlen), parameter :: CMPLX2D  = "CMPLX2D"
 ! 
 type INTERP_coeff
   character(schlen)        :: interp_type = " "
   complex(DP), allocatable :: engre(:,:)
   integer                  :: ndim = 0
 end type INTERP_coeff
 ! 
 type INTERP_element
   integer                  :: nb(2)=0
   integer                  :: nk=0
   integer                  :: D=0  ! 3rd dimension of SEED
   integer                  :: ID=0
   logical                  :: in_use=.FALSE.
   logical                  :: use_NN=.FALSE.
   logical                  :: has_fineGd=.FALSE.
   character(schlen)        :: seed_kind=" "
   character(schlen)        :: what=" "
   !
   real(SP),    allocatable :: SEED(:,:,:)
   real(SP),    allocatable :: k(:,:)
   real(SP),    allocatable :: weights(:)
   !
   real(SP),    allocatable :: INTERP(:,:,:)
   real(SP),    allocatable :: INTERP_FineGd(:,:,:)
 end type INTERP_element
 !
 type(INTERP_coeff)   :: interpls(max_interpls)
 type(INTERP_element), target :: INTERP_obj(max_interpls)
 !
 interface
   !
   subroutine INTERPOLATION_nearest_k(K_in,V_in,K_out,V_out)
     use pars,           ONLY:SP
     real(SP)                :: K_in(:,:)
     real(SP)                :: K_out(:,:)
     real(SP)                :: V_in(:,:,:)
     real(SP)                :: V_out(:,:,:)
   end subroutine
   !
   subroutine INTERPOLATION_coefficients(R1D,R2D,C1D,C2D,E,k,Nk,ID,ID_obj)
     use pars,       ONLY:SP
     use R_lattice,  ONLY:bz_samp
     use electrons,  ONLY:levels
     implicit none
     type(bz_samp), intent(in)           :: k
     real(SP),      intent(in), optional :: R1D(:,:),R2D(:,:,:)
     complex(SP),   intent(in), optional :: C1D(:,:),C2D(:,:,:)
     type(levels),  intent(in), optional :: E
     integer,       intent(in), optional :: Nk
     integer,       intent(in)           :: ID_obj
     integer,       intent(out)          :: ID
   end subroutine  INTERPOLATION_coefficients
   !
   subroutine INTERPOLATION_BZ(K,NK,ID,R1D,R2D,C1D,C2D,E)
     use pars,       ONLY:SP
     use R_lattice,  ONLY:bz_samp
     use electrons,  ONLY:levels
     implicit none
     type(bz_samp), intent(in)            :: K
     integer,       intent(in)            :: ID,NK
     real(SP),      intent(out), optional :: R1D(:,:),R2D(:,:,:)
     complex(SP),   intent(out), optional :: C1D(:,:),C2D(:,:,:)
     type(levels),  intent(out), optional :: E
   end subroutine INTERPOLATION_BZ
   !
   subroutine electrons_bands_interpolate(Xen,Xk,BANDS_range,ID_bands,CIRCUIT_E,GRID_E,CIRCUIT_dE)
     use pars,       ONLY:SP
     use electrons,  ONLY:levels
     use R_lattice,  ONLY:bz_samp
     implicit none
     type(bz_samp), intent(in)  :: Xk
     type(levels),  intent(in)  :: Xen
     integer     ,  intent(in)  :: BANDS_range(2)
     integer     ,  intent(out) :: ID_bands(2)
     type(levels),  intent(out) :: CIRCUIT_E
     type(levels),  intent(out), optional  :: GRID_E
     type(levels),  intent(out), optional  :: CIRCUIT_dE
   end subroutine electrons_bands_interpolate
   !
#if defined _RT
   subroutine INTERPOLATION_driver_seed(operation,ID,dims,SEED_k,SEED_E,SEED_RT,SEED_QP,SEED_R3,&
&                                       NEAREST,QUIET,SERIAL)
     use real_time,   ONLY:RT_t
#else
   subroutine INTERPOLATION_driver_seed(operation,ID,dims,SEED_k,SEED_E,        SEED_QP,SEED_R3,&
&                                       NEAREST,QUIET,SERIAL)
#endif
     use pars,        ONLY:SP
     use R_lattice,   ONLY:bz_samp
     use electrons,   ONLY:levels
     use QP_m,        ONLY:QP_t
     implicit none
     character(*),  intent(in)              :: operation
     integer,       intent(out)             :: ID
     integer,       intent(in)              :: dims(2)
     type(bz_samp), intent(in)              :: SEED_k
     real(SP),      intent(in),    optional :: SEED_R3(:,:,:)
#if defined _RT
     type(RT_t),    intent(in),    optional :: SEED_RT
#endif
     type(QP_t),    intent(in),    optional :: SEED_QP
     type(levels),  intent(in),    optional :: SEED_E
     logical,       intent(in),    optional :: NEAREST
     logical,       intent(in),    optional :: QUIET
     logical,       intent(in),    optional :: SERIAL
   end subroutine INTERPOLATION_driver_seed
   !
   subroutine INTERPOLATION_driver_do(operation,ID,OUT_k,OUT_E,FineGrid)
     use R_lattice,   ONLY:bz_samp
     use electrons,   ONLY:levels
     implicit none
     character(*),  intent(in)              :: operation
     integer,       intent(in)              :: ID
     type(bz_samp), intent(inout)           :: OUT_k
     type(levels),  intent(inout), optional :: OUT_E
     character(*),  intent(in)   , optional :: FineGrid
   end subroutine INTERPOLATION_driver_do
   !
 end interface  
 !
end module interpolate
