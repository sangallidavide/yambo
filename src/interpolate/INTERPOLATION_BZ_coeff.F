!
! License-Identifier: GPL
!
! Copyright (C) 2011 The Yambo Team
!
! Authors (see AUTHORS file for details): CA DS AM
!
subroutine INTERPOLATION_BZ_coeff(Xk,VALUEs,engre,nvalues,Nk,ID)
 !
 ! Interpolation scheme (PRB 38 p2721)
 ! Aug2002 Georg Madsen : First version based on subroutine from D.Singh
 ! Apr2004 Georg Madsen : blocked dgemm version
 ! Code take from BoltzTrap
 !
 ! https://www.imc.tuwien.ac.at//forschungsbereich_theoretische_chemie/forschungsgruppen/prof_dr_gkh_madsen_theoretical_materials_chemistry/boltztrap/
 !
 use pars,              ONLY:SP,DP,pi,cI,cZERO
 use interpolate,       ONLY:nshells,lattice_vectors,metric,int_sop,INTERP_obj,INTERPOLATE_is_quiet,&
&                            INTERPOLATE_is_serial,INTERP_FineGd,INTERPOLATE_use_symm,R_C1,R_C2
 use interpolate_tools, ONLY:make_star
 use D_lattice,         ONLY:nsym
 use R_lattice,         ONLY:bz_samp
 use vec_operate,       ONLY:c2a
 use LIVE_t,            ONLY:live_timing
 use parallel_m,        ONLY:PP_indexes,PP_indexes_reset,myid,PAR_build_index,ncpu
 use parallel_int,      ONLY:PP_redux_wait,PARALLEL_index 
 !
#include<y_memory.h>
 !
 integer,       intent(in)  :: nvalues
 integer,       intent(in)  :: Nk,ID
 type(bz_samp), intent(in)  :: Xk
 real(SP),      intent(in)  :: VALUEs(nvalues,Nk)
 complex(DP),   intent(out) :: engre(nshells,nvalues)
 !
 ! Work Space
 !
 type(PP_indexes)  :: px
 real(DP)          :: R2_min,R2,pre_factor
 real(SP)          :: v(3)
 integer           :: ik,ikp,iv,i_shell,i_star,nstar,Nkm1,info,&
&                     PAR_nshells,i_par_shell
 integer           :: nsym_interp
 !
 ! Parameters
 !
 real(DP), parameter :: twopi=2._SP*pi
 !
 ! Allocatable
 !
 real(DP),    allocatable :: rho(:)
 complex(SP), allocatable :: h_mat(:,:),delta(:,:)
 complex(DP), allocatable :: BZ_star(:,:),PAR_engre(:,:)
 integer,     allocatable :: ipiv(:),SHELL_index(:)
 real(SP),    allocatable :: star_vec(:,:)
 !
 ! External
 !
 Nkm1=Nk-1           ! number of k-points minus 1
 !
 YAMBO_ALLOC(delta,(Nkm1,nvalues))
 YAMBO_ALLOC(ipiv,(Nkm1))
 !
 ! Construct delta_epsilon eq. 10
 !
 delta=cZERO
 !
 nsym_interp=nsym
 if(.not.INTERPOLATE_use_symm) nsym_interp=1
 YAMBO_ALLOC(star_vec,(3,nsym_interp))
 !
 do ik=1,Nkm1
   delta(ik,:)=VALUEs(:,ik)-VALUEs(:,Nk)
 enddo
 !
 PAR_nshells=nshells
 call PP_indexes_reset(px)
 YAMBO_ALLOC(SHELL_index,(nshells))
 if (INTERPOLATE_is_serial) then
   YAMBO_ALLOC(px%element_1D,(nshells))
   YAMBO_ALLOC(px%n_of_elements,(ncpu))
   px%n_of_elements=nshells
   px%element_1D=.TRUE.
   do i_shell=1,nshells
     SHELL_index(i_shell)=i_shell
   enddo
 else
   call PARALLEL_index(px,(/nshells/))
   call PAR_build_index(px,nshells,SHELL_index,px%n_of_elements(myid+1))
   PAR_nshells=px%n_of_elements(myid+1)
 endif
 !
 ! Rho
 !
 YAMBO_ALLOC(rho,(PAR_nshells))
 rho=0._DP
 R2_min=dot_product(lattice_vectors(:,2),matmul(metric,lattice_vectors(:,2)))
 do i_shell=1,nshells
   if (.not.px%element_1D(i_shell)) cycle
   i_par_shell=SHELL_index(i_shell)
   if (i_shell==1) cycle
   R2=dot_product(lattice_vectors(:,i_shell),matmul(metric,lattice_vectors(:,i_shell)))
   rho(i_par_shell)=roughness_func(R2,R2_min)
 enddo
 !
 YAMBO_ALLOC(BZ_star,(PAR_nshells,Nk))
 BZ_star=(0._DP,0._DP)
 !
 if (.not.INTERPOLATE_is_quiet) call live_timing("Interpolation@"//trim(INTERP_obj(ID)%what)//" (star)",PAR_nshells)
 !==================================================================================================================
 do i_shell=1,nshells
   !
   if (i_shell==1) then
     if (.not.INTERPOLATE_is_quiet) call live_timing(steps=1)
     cycle
   endif
   !
   if (.not.px%element_1D(i_shell)) cycle
   i_par_shell=SHELL_index(i_shell)
   !
   call make_star(lattice_vectors(:,i_shell),nsym_interp,int_sop,nstar,star_vec) 
   do ik=1,Nk
     if (trim(INTERP_FineGd)=="NONE") call c2a(v_in=Xk%pt(ik,:),v_out=v,mode='ki2a')
     if (trim(INTERP_FineGd)=="IBZ")  call c2a(v_in=Xk%FGibz%pt(ik,:),v_out=v,mode='ki2a')
     if (trim(INTERP_FineGd)=="BZ")   call c2a(v_in=Xk%FGbz%pt(ik,:),v_out=v,mode='ki2a')
     if (trim(INTERP_FineGd)=="BARE") call c2a(v_in=Xk%FGbare%pt(ik,:),v_out=v,mode='ki2a')
     do i_star=1,nstar
       BZ_star(i_par_shell,ik)=BZ_star(i_par_shell,ik) + exp(cI*twopi*dot_product(v(:),star_vec(:,i_star)))
     enddo
   enddo
   !
   BZ_star(i_par_shell,:)=BZ_star(i_par_shell,:)/real(nstar,DP)
   !
   do ik=1,Nkm1
     BZ_star(i_par_shell,ik) = BZ_star(i_par_shell,ik) - BZ_star(i_par_shell,Nk)
     BZ_star(i_par_shell,ik) = BZ_star(i_par_shell,ik) / sqrt(rho(i_par_shell))
   enddo
   !
   if (.not.INTERPOLATE_is_quiet) call live_timing(steps=1)
   !
 enddo
 !
 if (.not.INTERPOLATE_is_quiet)  call live_timing()
 !
 if (.not.INTERPOLATE_is_quiet) call live_timing("Interpolation@"//trim(INTERP_obj(ID)%what)//" (H fill)",Nkm1*(Nkm1-1)/2)
 !=======================================================================================================================
 YAMBO_ALLOC(h_mat,(Nkm1,Nkm1))
 h_mat  =cZERO
 !
 do ik=1,Nkm1
   h_mat(ik,ik)=real(dot_product(BZ_star(:,ik),BZ_star(:,ik)),kind=SP)
   do ikp=ik+1,Nkm1
     h_mat(ik,ikp)=cmplx(dot_product(BZ_star(:,ikp),BZ_star(:,ik)),kind=SP)
     h_mat(ikp,ik)=conjg(h_mat(ik,ikp))
     if (.not.INTERPOLATE_is_quiet) call live_timing(steps=1)
   enddo
 enddo
 !
 if (.not.INTERPOLATE_is_serial) call PP_redux_wait(h_mat)
 if (.not.INTERPOLATE_is_quiet)  call live_timing()
 !
 ! Lapack CALL
 !============
 !
#if defined _DOUBLE
 call ZGETRF(Nkm1,Nkm1,h_mat,Nkm1,ipiv,info)
#else
 call CGETRF(Nkm1,Nkm1,h_mat,Nkm1,ipiv,info)
#endif
 !
 if(info/=0) call error("[INTERPOLATION] Error in factorization ")
 !
#if defined _DOUBLE
 call  ZGETRS('N',Nkm1,nvalues,h_mat,Nkm1,ipiv,delta,Nkm1,info)
#else
 call  CGETRS('N',Nkm1,nvalues,h_mat,Nkm1,ipiv,delta,Nkm1,info)
#endif
 !
 if(info/=0) call error("[INTERPOLATION] Error in GETRS")
 !
 if (.not.INTERPOLATE_is_quiet) call live_timing("Interpolation@"//trim(INTERP_obj(ID)%what)//" (computing engre)",PAR_nshells)
 !=======================================================================================================================
 YAMBO_ALLOC(PAR_engre,(PAR_nshells,nvalues))
 PAR_engre=cZERO
 engre    =cZERO
 !
 do ik=1,Nkm1
   BZ_star(:,ik) = BZ_star(:,ik) * sqrt(rho(:))
 enddo
 do i_shell=1,nshells
   if (i_shell==1) then
     if (.not.INTERPOLATE_is_quiet) call live_timing(steps=1)
     cycle
   endif
   if (.not.px%element_1D(i_shell)) cycle
   i_par_shell=SHELL_index(i_shell)
   do iv=1,nvalues
     PAR_engre(i_par_shell,iv)= dot_product(BZ_star(i_par_shell,1:Nkm1)/rho(i_par_shell),cmplx(delta(1:Nkm1,iv),kind=DP))
   enddo
   engre(i_shell,:)=PAR_engre(i_par_shell,:)
   if (.not.INTERPOLATE_is_quiet) call live_timing(steps=1)
 enddo
 !
 if (.not.INTERPOLATE_is_quiet)  call live_timing()
 !
 pre_factor=1._DP/real(ncpu,DP)
 if (INTERPOLATE_is_serial) pre_factor=1._DP
 do iv=1,nvalues
 engre(1,iv)=VALUEs(iv,Nk)*pre_factor - dot_product(conjg(PAR_engre(:,iv)),BZ_star(:,Nk))
 enddo
 !
 if (.not.INTERPOLATE_is_serial) call PP_redux_wait( engre )
 !
 ! Deallocation
 !
 YAMBO_FREE(h_mat)
 YAMBO_FREE(BZ_star)
 YAMBO_FREE(ipiv)
 YAMBO_FREE(delta)
 YAMBO_FREE(rho)
 YAMBO_FREE(SHELL_index)
 YAMBO_FREE(PAR_engre)
 YAMBO_FREE(star_vec)
 !
 call k_compare(Xk,VALUEs,engre,nshells,nvalues,Nk)
 !
 contains
  !
  pure real(DP) function roughness_func(R2,R2_min)
   real(DP),  intent(in) :: R2,R2_min
   !
   roughness_func=(1._DP-R_C1*R2/R2_min)**2 + R_C2*(R2/R2_min)**3
   !
  end function roughness_func
  !
  subroutine k_compare(Xk,VALUEs,engre,nshells,nvalues,Nk)
   use units,          ONLY:HA2EV
   use pars,           ONLY:SP,lchlen
   use R_lattice,      ONLY:bz_samp
   use com,            ONLY:msg
   !
   type(bz_samp), intent(in) :: Xk
   integer,       intent(in) :: nvalues,nshells,Nk
   real(SP),      intent(in) :: VALUEs(nvalues,Nk)
   complex(DP),   intent(in) :: engre(nshells,nvalues)
   !
   ! Work Space
   !
   integer :: ik,iv
   character(lchlen)     :: dump_ch
   real(SP), allocatable :: new_VALUEs(:,:)
   real(SP) :: ave_err,max_err,rel_ave_err
   !
   YAMBO_ALLOC(new_VALUEs,(nvalues,Nk))
   !
   call INTERPOLATION_fourier(Xk,new_VALUEs,engre,nvalues,Nk,ID)
   !
   max_err=0._SP
   ave_err=0._SP
   rel_ave_err=0._SP
   !
   do ik=1,Nk
     do iv=1,nvalues
       ave_err=ave_err+abs(new_VALUEs(iv,ik)-VALUEs(iv,ik))
       if(abs(new_VALUEs(iv,ik))>0._SP) rel_ave_err=rel_ave_err+abs(new_VALUEs(iv,ik)-VALUEs(iv,ik))/abs(new_VALUEs(iv,ik))
       if(abs(new_VALUEs(iv,ik)-VALUEs(iv,ik))>max_err) max_err=abs(new_VALUEs(iv,ik)-VALUEs(iv,ik))
     enddo
   enddo
   !
   ave_err=ave_err/real(nvalues*Nk)
   rel_ave_err=rel_ave_err/real(nvalues*Nk)
   !
   write(dump_ch,'(a,3es16.8)') '[INTERPOLATION@'//trim(INTERP_obj(ID)%what)//&
&                               '] Errors: Max, Avg, Avg_rel = ',max_err,ave_err,rel_ave_err
   call msg('s',dump_ch)
   !
   YAMBO_FREE(new_VALUES)
   !
  end subroutine k_compare 
  !
end subroutine INTERPOLATION_BZ_coeff
