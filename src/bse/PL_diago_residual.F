!
!        Copyright (C) 2000-2022 the YAMBO team
!              http://www.yambo-code.org
!
! Authors (see AUTHORS file for details): AM
! 
! This file is distributed under the terms of the GNU 
! General Public License. You can redistribute it and/or 
! modify it under the terms of the GNU General Public 
! License as published by the Free Software Foundation; 
! either version 2, or (at your option) any later version.
!
! This program is distributed in the hope that it will 
! be useful, but WITHOUT ANY WARRANTY; without even the 
! implied warranty of MERCHANTABILITY or FITNESS FOR A 
! PARTICULAR PURPOSE.  See the GNU General Public License 
! for more details.
!
! You should have received a copy of the GNU General Public 
! License along with this program; if not, write to the Free 
! Software Foundation, Inc., 59 Temple Place - Suite 330,Boston, 
! MA 02111-1307, USA or visit http://www.gnu.org/copyleft/gpl.txt.
!
subroutine PL_diago_herm_residual(BS_mat,BS_R_PL) 
 !
 ! Resonant PL Residuals
 !
  ! Resonant Residuals
 !====================
 ! BS_R_right(i) = \sum_k <q0|k><k|i> =
 !                 \sum_k BS_mat(k,i) x conjg( BSS_dipoles(k) )
 ! BS_R_left(i)  = conjg( BS_R_right(i) )
 !
 use pars,           ONLY:SP,rZERO,cZERO,cI
 use wrapper_omp,    ONLY:V_dot_V_omp,Vstar_dot_V_omp
 use parallel_m,     ONLY:PP_indexes,myid,PP_indexes_reset
 use parallel_int,   ONLY:PP_redux_wait,PARALLEL_index 
 use LIVE_t,         ONLY:live_timing
 use BS,             ONLY:BS_H_dim,BS_K_dim,BS_dip_size
 use BS_solvers,     ONLY:BSS_eh_f_RES,BSS_n_eig, BS_E
 use PHOTOLUM,       ONLY:BSS_dipoles_PL,BSS_PL_f
 use functions,      ONLY:boltzman_f
 use RT_control,     ONLY:NEQ_Residuals
 use D_lattice,      ONLY:Boltz_Temp
 !
 implicit none
 !
 complex(SP),intent(in)  :: BS_mat(BS_K_dim(1),BSS_n_eig)
 real(SP)   ,intent(out) :: BS_R_PL(2,BSS_n_eig)
 !
 ! Workspace
 !
 type(PP_indexes)     ::px
 logical              ::compute_ares_explicitly
 integer              ::i_l,i_c,i_K,i_lr,i_ll, ik
 complex(SP)          ::P_x_fZ(BS_H_dim,3),Pi_x_R_x_adA_left(3),Pi_x_R_x_adA_right(3),adA_x_R_over_R_x_A,&
&                       R_over_R_x_A(BS_K_dim(1))
!Write to a file
 ! Parallel indexes
 !
 call PP_indexes_reset(px)
 call PARALLEL_index(px,(/BSS_n_eig/))
 !
 ! sqrt(R_i) Re[ (Pi_i)* A^l_i (A^l_j)* R<_j/R_j A^l_j (A^l_k)* Pi_k ] sqrt(R_k)
 !
 ! Zeroing
 !
 BS_R_PL       =rZERO
 P_x_fZ        =cZERO
 R_over_R_x_A  =cZERO
 !
 ! Live-Timing
 !
 call live_timing('PL    residuals',px%n_of_elements(myid+1))
 !
 ! Hermitian residual
 !====================
 !
 ! 1=> left
 ! 2=> right
 !
 if (NEQ_Residuals .and. Boltz_Temp .gt. 0) then
    call warning ('The user can specify either boltzmann or external residual occupations, not both')
 endif 
 do i_c=1,BS_dip_size
   do i_K=1,BS_K_dim(1)
     P_x_fZ(i_K,i_c)=conjg(BSS_dipoles_PL(i_c,i_K))*sqrt(cmplx(BSS_eh_f_RES(i_K))) !Pi^* left sqrt(R)
   enddo
 enddo
 !
 do i_l=1,BSS_n_eig
  !
  if (.not.px%element_1D(i_l)) cycle !parallelization
  !                       
  ! R_over_R_x_A = R<_j/R_j A^r_ij
  ! Here I create an If statement if I want to populate L lesser with electornic or excitonic populations
  if (NEQ_Residuals) then
    do ik = 1, BS_K_dim(1)
      call msg('s','Running with NEQ Residuals')
      R_over_R_x_A(ik) =  BSS_PL_f(ik)*BS_mat(ik,i_l)/BSS_eh_f_RES(ik) !resonant R</R*Ar
    enddo
    !
    ! adA_x_R_over_R_x_A = (A^l_j)* R<_j/R_j A^l_j
    !
    adA_x_R_over_R_x_A=V_dot_V_omp(BS_K_dim(1),conjg(BS_mat(1:BS_K_dim(1),i_l)),& !BS_mat should be conjg but maybe V_dot_V
            !does it
            R_over_R_x_A(1:BS_K_dim(1))) !resonant this is dipoles_PL*BS_mat*fermifactors= sum_k <q0|k><k|i>
    !antiresonant, no error if two arrays have different dimensions, We don't know if it's conjg or not
    !Al*R</R*Ar
  else if (Boltz_Temp .gt. 0) then 
    adA_x_R_over_R_x_A=boltzman_f(real(BS_E(i_l)-BS_E(1)))
  endif
  do i_c=1,BS_dip_size
    !dim(P_x_fz) = BS_H_dim
    ! we have a problem in this multiplication
    Pi_x_R_x_adA_left(i_c) =V_dot_V_omp(BS_K_dim(1),P_x_fZ(:BS_K_dim(1),i_c),BS_mat(:,i_l)) !resonan Pi^* Ar*sqrt(R)                                   
    ! Pi_x_R_x_adA_left =   Pi sqrt(R) dipoles_PL*BS_mat/R
  enddo
  !
  !Check dimensions for these vectors to be consistent that we have only resonant and antiresonant
  BS_R_PL(1,i_l)=real(dot_product(Pi_x_R_x_adA_left,Pi_x_R_x_adA_left)*adA_x_R_over_R_x_A) !resonant
  !Pi*Al*Pi^**Ar*Al*R</R*Ar
  !
 enddo
 !
 BS_R_PL(2,:)=-BS_R_PL(1,:)
 !
 call live_timing()
 !
 call PP_redux_wait(BS_R_PL)
 !
 call PP_indexes_reset(px)
 !
end subroutine PL_diago_herm_residual

subroutine PL_diago_non_herm_residual(BS_V_left,BS_V_right,BS_R_PL,BS_overlap) 
 !
 ! Non-Hermitian residual
 !========================
 !
 ! Resonant PL Residuals
 !
 use pars,           ONLY:SP,rZERO,cZERO,cI
 use wrapper_omp,    ONLY:V_dot_V_omp,Vstar_dot_V_omp
 use parallel_m,     ONLY:PP_indexes,myid,PP_indexes_reset
 use parallel_int,   ONLY:PP_redux_wait,PARALLEL_index 
 use LIVE_t,         ONLY:live_timing
 use BS,             ONLY:BS_H_dim,BS_K_dim,BS_dip_size
 use BS_solvers,     ONLY:BSS_eh_f_RES,BSS_n_eig
 use PHOTOLUM,       ONLY:BSS_dipoles_PL,BSS_PL_f
 !
 implicit none
 !
 complex(SP),intent(in)  :: BS_V_left(BS_H_dim,BSS_n_eig),BS_V_right(BS_H_dim,BSS_n_eig)
 real(SP)   ,intent(out) :: BS_R_PL(2,BSS_n_eig)
 complex(SP),optional,intent(in)  :: BS_overlap(BSS_n_eig,BSS_n_eig)
 !
 ! Workspace
 !
 type(PP_indexes)     ::px
 integer              ::i_l,i_c,i_K,i_lr,i_ll, ik
 complex(SP)          ::P_x_fZ(BS_H_dim,3,2),Pi_x_R_x_adA_left(3),Pi_x_R_x_adA_right(3),adA_x_R_over_R_x_A(2),&
&                       R_over_R_x_A(2,BS_H_dim),B_var(BS_H_dim)
 !
 ! Parallel indexes
 !
 call PP_indexes_reset(px)
 call PARALLEL_index(px,(/BSS_n_eig/))
 !
 ! sqrt(R_i) Re[ (Pi_i)* A^l_i (A^l_j)* R<_j/R_j A^l_j (A^l_k)* Pi_k ] sqrt(R_k)
 !
 ! Zeroing
 !
 BS_R_PL       =rZERO
 P_x_fZ        =cZERO
 R_over_R_x_A  =cZERO
 !
 ! Live-Timing
 !
 call live_timing('PL    residuals',px%n_of_elements(myid+1))
 !
 ! 1=> left
 ! 2=> right
 !
 do i_c=1,BS_dip_size
    do i_K=1,BS_H_dim
      P_x_fZ(i_K,i_c,1)=conjg(BSS_dipoles_PL(i_c,i_K))*sqrt(cmplx(BSS_eh_f_RES(i_K)))
      ! Here there is a conjg becaus we are getting a -1 somewhere from
      ! the occupation factors. This seems to fix it. Not sure it is the
      ! proper solution
      P_x_fZ(i_K,i_c,2)=      BSS_dipoles_PL(i_c,i_K)*conjg(sqrt(cmplx(BSS_eh_f_RES(i_K))))
    enddo
  enddo
  !
  do i_lr=1,BS_H_dim
    !
    if (.not.px%element_1D(i_lr)) cycle
    !
    ! B_i = \sum_ll Overlap^*(ll,lr) A^ll_i
    !
    B_var=cZERO
    do i_ll=1,BS_H_dim
      B_var(:)=B_var(:)+conjg(BS_overlap(i_ll,i_lr))*BS_V_left(:,i_ll)
    enddo
    !
    ! R_over_R_x_A = R<_j/R_j B^l_j
    !
    do i_ll=1,BS_H_dim
      R_over_R_x_A(1,i_ll) = BSS_PL_f(i_ll)*B_var(i_ll)/BSS_eh_f_RES(i_ll)
    enddo
    !
    ! adA_x_R_over_R_x_A = B^l*_j    R<_j/R_j B^l_j
    !
    adA_x_R_over_R_x_A(1)=Vstar_dot_V_omp(BS_H_dim,B_var,R_over_R_x_A(1,:))
    !
    do i_c=1,BS_dip_size
      !
      ! Pi_x_R_x_adA_left =   sqrt(R) Pi_i* A^lr_i
      !
      Pi_x_R_x_adA_left(i_c) =V_dot_V_omp(BS_H_dim,P_x_fZ(:,i_c,1),BS_V_right(:,i_lr))
      !
      ! Pi_x_R_x_adA_right=   A^lr*_i Pi_i sqrt(R)  
      !
      Pi_x_R_x_adA_right(i_c)=Vstar_dot_V_omp(BS_H_dim,BS_V_right(:,i_lr),P_x_fZ(:,i_c,2))
      !
      BS_R_PL(1,i_lr)=BS_R_PL(1,i_lr)+real(Pi_x_R_x_adA_left(i_c)*Pi_x_R_x_adA_right(i_c))
      !
    enddo
    !
    BS_R_PL(1,i_lr)=BS_R_PL(1,i_lr)*adA_x_R_over_R_x_A(1)
    !
    call live_timing(steps=1)
    !
  enddo
 !
 call live_timing()
 !
 call PP_redux_wait(BS_R_PL)
 !
 BS_R_PL(2,:)=-BS_R_PL(1,:)
 !
 call PP_indexes_reset(px)
 !
end subroutine PL_diago_non_herm_residual
