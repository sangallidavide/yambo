!
! License-Identifier: GPL
!
! Copyright (C) 2006 The Yambo Team
!
! Authors (see AUTHORS file for details): AM
!
module collision_ext
 !
 !====================
 ! EXTENDED COLLISION
 !====================
 !
 use pars,       ONLY:SP,SP6,LP,schlen
 use parallel_m, ONLY:PP_indexes
 !
 implicit none
 !
 integer :: COLL_bands(2)
 integer :: COLL_ID
 integer :: LAST_COLL_sync
 integer :: PAR_COLL_min
 !
 ! Cut off on exchange matrix elements
 ! Sx_cutoff=0 all, (Sx_cutoff>1 or Sx_cutoff<0)=none
 !
 real(SP)                 :: COLLISIONS_cutoff
 real(SP)                 :: COLLISIONS_max
 !
 logical :: COLLISIONS_parallel_SC_structure=.FALSE.
 logical :: COLLISIONS_have_HARTREE=.FALSE.
 logical :: COLLISIONS_HXC_local=.FALSE.
 logical :: COLLISIONS_HXC_MB=.FALSE.
 logical :: COLLISIONS_HXC_use_TDDFT=.FALSE.
 logical :: COLLISIONS_CV_only=.FALSE.
 logical :: COLLISIONS_load_SP=.FALSE.  ! Force loading of collision in single precision
 logical :: COLLISIONS_compr  =.FALSE. ! If parallel-io is used collisions in memory are compressed
 !
 ! Elements used in the code
 !
 ! Evaluation (COLLISION_sex ...)
 !
 ! SEX, HF (eval): v3_c
 !          (I/O): v_c
 !      (compose): v_c
 !
 ! GW(NEQ),(QED) (eval): v3_r
 !                (I/O): v_r
 !                (USE): v3_r
 !
 type COLLISIONS_element
   integer               :: I=0         ! State Index 
   integer               :: N=0         ! Linear Dimension
   character,      allocatable  :: table(:,:,:)
   complex(SP),  allocatable  :: v_c(:)      
   complex(SP6), allocatable  :: v_c_SP(:)    ! SEX and HF only in single precision
   complex(SP), allocatable  :: v3_c(:,:,:)
   real(SP)   , allocatable  :: v_r(:)       ! GW NEQ collisions
   real(SP)   , allocatable  :: v3_r(:,:,:)
 end type
 !
 type COLLISIONS_group
   integer               :: D(3,2)=0    ! Dimension (v3 of corresponding COLLISIONS_element)
   character(schlen)     :: name="none"
   integer               :: kind=0
   integer               :: N=0         ! # of collisions (dimension of linked COLLISIONS_element)
   integer,     allocatable  :: state(:,:)
   !
   ! I/O
   !
   integer,     allocatable  :: IO_status(:)
   !
   ! Parallel Substructure
   !
   type(PP_indexes)      :: PAR_IND
   integer               :: PAR_N=0
   integer, allocatable  :: PAR_map(:)
   !
 end type
 !
 ! COLLISIONS groups...
 !======================
 type(COLLISIONS_group), save :: GW_NEQ_collisions
 type(COLLISIONS_group), save :: COH_collisions
 type(COLLISIONS_group), save :: HXC_collisions
 type(COLLISIONS_group), save :: P_collisions
 !
 ! .. and ELEMENTS
 !=================
 type(COLLISIONS_element), allocatable :: GW_NEQ_COLL_element(:)
 type(COLLISIONS_element), allocatable :: COH_COLL_element(:)
 type(COLLISIONS_element), allocatable :: HXC_COLL_element(:)
 type(COLLISIONS_element), allocatable :: P_COLL_element(:)
 !
 ! On the fly long-range SEX (LSEX) collisions
 ! ====================================
 integer :: ng_oscll                  ! number of G-vectors in the Oscillators
 complex(SP), allocatable :: W(:,:,:) ! screened interaction
 complex(SP), allocatable :: OSCLL(:,:,:,:,:) ! Oscillators
 complex(SP), allocatable :: OSCLL_k(:,:,:,:) ! Oscillators in k, just for the evaluation
 !
 contains
   !
   subroutine COLLISIONS_naming( H_SE_kind , what )
     !
     ! In this routine I define the %name/%kind of the COLLISIONS.
     ! If no H_SE_kind is given I set defaults, otherwise I define the proper
     ! componentes of the HXC *ONLY* part.
     !
     use global_XC,      ONLY:H_SE_COH,H_SE_HARTREE,H_SE_SEX,QP_SE_GW_QED,QP_SE_GW_NEQ,H_SE_FOCK
     use xc_functionals, ONLY:GS_xc_KIND
     !
     integer            :: H_SE_kind
     character(*)       :: what
     !
     integer, parameter :: nc=7
     integer            :: i_c,local_kind(nc)
     character(schlen)  :: local_name(nc)
     !
     if (H_SE_kind==0) then
       !
       COH_collisions%kind=H_SE_COH
       COH_collisions%name="COH"
       HXC_collisions%kind=H_SE_HARTREE+H_SE_SEX
       HXC_collisions%name="HARTREE+SEX"
       P_collisions%kind=QP_SE_GW_QED
       P_collisions%name="Momentum"
       GW_NEQ_collisions%kind=QP_SE_GW_NEQ
       GW_NEQ_collisions%name="GW_NEQ"
       !
     else if (what=="HXC") then
       !
       HXC_collisions%kind=0
       HXC_collisions%name=" "
       !
       local_kind=(/H_SE_FOCK,H_SE_SEX,H_SE_Hartree,H_SE_FOCK+H_SE_Hartree,H_SE_SEX+H_SE_Hartree,&
&                   GS_xc_KIND,GS_xc_KIND+H_SE_Hartree/)
       local_name(1)="FOCK"
       local_name(2)="SEX"
       local_name(3)="HARTREE"
       local_name(4)="HARTREE+FOCK"
       local_name(5)="HARTREE+SEX"
       local_name(6)="GS_xc"
       local_name(7)="HARTREE+GS_XC"
       !
       do i_c=1,nc
         if (H_SE_kind==local_kind(i_c).or.H_SE_kind==local_kind(i_c)+H_SE_COH) then
           HXC_collisions%kind=local_kind(i_c)
           HXC_collisions%name=local_name(i_c)
           return
         endif
       enddo
       !
     endif
   end subroutine
   !
end module collision_ext
