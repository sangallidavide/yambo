!
! License-Identifier: GPL
!
! Copyright (C) 2016 The Yambo Team
!
! Authors (see AUTHORS file for details): AM
!
subroutine k_grids(en,k,Xk,q)
 !
 use pars,                ONLY:SP,schlen,zero_dfl
 use electrons,           ONLY:levels
 use R_lattice,           ONLY:bz_samp,bz_samp_reset,nXkibz,k_the_nearest
 use D_lattice,           ONLY:nsym,lattice
 use stderr,              ONLY:intc
 use com,                 ONLY:msg
 use parser_m,            ONLY:parser
 use vec_operate,         ONLY:c2a,rlu_v_is_zero,v_is_zero,v_norm
 use YPPm,                ONLY:coo_in,coo_out,l_k_grid,l_q_grid,l_regular_grid,&
&                              l_shifted_grid,l_high_sym_pts,K_transform,&
&                              n_u_k_pts,n_u_q_pts,l_random_grid,PtsPath,CIRCUIT_k,INTERP_grid,&
&                              U_k_pts,U_q_pts,U_q_shifts,USER_k,k_grid_printer,N_path_pts,PtsPath
 use LIVE_t,              ONLY:live_timing_is_on
 use zeros,               ONLY:k_iku_zero
 use YPP_interfaces,      ONLY:k_special
 !
#include<y_memory.h>
 !
 type(levels) ::en       
 type(bz_samp)::k,Xk,q
 !
 ! Work Space
 !
 type(bz_samp)::FINAL_k,Q_user
 real(SP)     ::v(3),ka(3),q0(3),kout(Xk%nibz,3),k_plus_q(Xk%nibz,3)
 integer      ::i1,i2,is,ik,iq,ik_near,i_err
 integer,  allocatable :: int_ctl(:)
 real(SP), allocatable :: BZ_weights(:),K_dist(:)
 character(schlen)     :: ch
 logical               :: No_Weight,Force_User_points,List_internal_points,Expand_internal_points,l_new_gw_point,l_shift01
 !
 call bz_samp_reset(USER_k)
 call bz_samp_reset(Q_user)
 call bz_samp_reset(FINAL_k)
 !
 call section('*',"BZ grid analyzer/generator")  
 !=============================================
 !
 ! Setup logicals 
 !
 if (l_high_sym_pts) l_high_sym_pts=trim(lattice)/='none'
 !
 ! Q/Xk (IBZ->BZ)
 !
 call k_ibz2bz( q,'i',.false.)
 call k_ibz2bz(Xk,'i',l_random_grid)
 !
 call msg("s",'Q-points  (BZ)', q%nbz)
 call msg("s",'K-points  (BZ)',Xk%nbz)
 !
 call parser('NoWeights'   ,No_Weight)
 call parser('ForceUserPts',Force_User_points)
 call parser('ListPts'     ,List_internal_points)
 call parser('ExpandPts'   ,Expand_internal_points)
 call parser('Shift01'     ,l_shift01)
 !
 if (l_random_grid) then
   call k_random(Xk,No_Weight)
   goto 1
 endif
 !
 if (len_trim(coo_in)==0) return
 !
 if (l_k_grid) then
   !
   call k_special( )
   !
   call k_circuit_made_of_special_points(PtsPath,0)
   !
   USER_k%nbz=n_u_k_pts+CIRCUIT_k%nbz
   YAMBO_ALLOC(USER_k%pt,(USER_k%nbz,3))
   YAMBO_ALLOC(K_dist,(Xk%nbz))
   !
   call section('=',"== GW K-grid analyzer ==")
   !===========================================
   !
   USER_k%nbz=0
   do ik=1,n_u_k_pts+CIRCUIT_k%nbz
     !
     l_new_gw_point=.TRUE.
     !
     if (ik<=n_u_k_pts) then
       v=U_k_pts(ik,:)
       if (trim(coo_in)=="cc")  call c2a(v_in=v,mode="kc2i")
       if (trim(coo_in)=="rlu") call c2a(v_in=v,mode="ka2i")
     else
       v=CIRCUIT_k%pt(ik-n_u_k_pts,:)
     endif
     !
     ik_near=k_the_nearest(v,Xk%ptbz,Xk%nbz,k_dist=K_dist,i_err=i_err)
     !
     call msg('sr','User point         ',v,"[iku]")
     if (i_err==0) then
       call msg('sr','Grid correspondance',Xk%sstar(ik_near,1))
       l_new_gw_point=.FALSE.
     else
       call msg('sr','Grid nearest       ',Xk%sstar(ik_near,1))
       call msg('sr','Distance           ',K_dist(ik_near),"[iku]")
       l_new_gw_point=.TRUE.
     endif
     !
     if (l_new_gw_point) then
       USER_k%nbz=USER_k%nbz+1
       USER_k%pt(USER_k%nbz,:)=v
     endif
     !
   enddo
   !
   call msg("sr",'User K-points ',USER_k%nbz)
   !
   ! Extended grid
   !
   YAMBO_ALLOC(FINAL_k%ptbz,(Xk%nibz+USER_k%nbz*q%nbz,3))
   FINAL_k%ptbz(:Xk%nibz,:)=k%pt(:Xk%nibz,:)
   FINAL_k%nbz=Xk%nibz
   do i1=1,USER_k%nbz
     do i2=1,q%nbz
       FINAL_k%nbz=FINAL_k%nbz+1
       FINAL_k%ptbz(FINAL_k%nbz,:)=USER_k%pt(i1,:)-q%ptbz(i2,:)
     enddo
   enddo
   !
   if (USER_k%nbz>0) call reduce_and_print(FINAL_k,.TRUE.)
   !
   call bz_samp_reset(FINAL_k)
   YAMBO_FREE(K_dist)
   !
 endif
 !
 if (l_q_grid) then
   !
   call section('=',"== Q-grid analyzer ==")
   !========================================
   !
   ! Input file parsing
   !
   Q_user%nibz=1
   q_main_loop: do iq=1,n_u_q_pts
     !
     YAMBO_FREE(Q_user%pt)
     YAMBO_ALLOC(Q_user%pt,(Q_user%nibz,3))
     Q_user%pt(Q_user%nibz,:)=U_q_pts(iq,:)
     !
     do i1=1,q%nbz
       v=q%ptbz(i1,:)
       if (trim(coo_in)=="cc")  call c2a(v_in=q%ptbz(i1,:),v_out=v,mode="ki2c")
       if (trim(coo_in)=="rlu") call c2a(v_in=q%ptbz(i1,:),v_out=v,mode="ki2a")
       if (v_is_zero(v(:)-Q_user%pt(Q_user%nibz,:),zero_=k_iku_zero).and..not.Force_User_points) then
         call warning(' Point already in the grid. Listing internal Q-points.')
         cycle q_main_loop
       endif
     enddo
     !
     Q_user%nibz=Q_user%nibz+1
     !
   enddo q_main_loop
   Q_user%nibz=Q_user%nibz-1
   !
   if (Q_user%nibz>=1) call msg("s",'User Q-points ',Q_user%nibz)
   !
   if (Force_User_points) then
     call k_grid_printer(Q_user%nibz,Q_user%pt,trim(coo_in),"Y")
     goto 1
   endif
   !
   ! Extended grid
   !
   YAMBO_ALLOC(FINAL_k%ptbz,(Xk%nibz+Xk%nbz*Q_user%nibz,3))
   FINAL_k%ptbz(:Xk%nibz,:)=k%pt(:Xk%nibz,:)
   FINAL_k%nbz=Xk%nibz
   do i1=1,Q_user%nibz
     v=Q_user%pt(i1,:)
     if (trim(coo_in)=="cc")  call c2a(v_in=Q_user%pt(i1,:),v_out=v,mode="kc2i")
     if (trim(coo_in)=="rlu") call c2a(v_in=Q_user%pt(i1,:),v_out=v,mode="ka2i")
     do i2=1,Xk%nbz
       FINAL_k%nbz=FINAL_k%nbz+1
       FINAL_k%ptbz(FINAL_k%nbz,:)=v(:)-Xk%ptbz(i2,:)
     enddo
   enddo
   !
   if (Q_user%nibz>0) call reduce_and_print(FINAL_k,.FALSE.)
   !
   call bz_samp_reset(FINAL_k)
   !
 endif
 !
 if (l_shifted_grid) then
   !
   call section('=',"== Built-in IBZ K-grid ["//trim(coo_out)//"] ==")
   YAMBO_ALLOC(BZ_weights,(Xk%nibz))
   BZ_weights(:)=Xk%weights(:)*real(Xk%nbz,SP)
   call k_grid_printer(Xk%nibz,Xk%pt,'iku',"PW",wk=BZ_weights)
   YAMBO_FREE(BZ_weights)
   !
   ! Generate shifted k-points set for dipole calculation
   !
   call section('=',"== Shifted grids generator ==")
   !================================================
   !
   if(any(abs(U_q_shifts)/=0._SP)) then
     !
     kout(:,:) = Xk%pt(:,:)
     !
     ! Convert IBZ mesh to coo_out
     !
     do ik=1,Xk%nibz
       call K_transform(kout(ik,:),'iku')
     enddo
     !
     i2=0
     do i1=1,3
       !
       if (v_norm(U_q_shifts(i1,:))>1.E-4.or.v_norm(U_q_shifts(i1,:))<1.E-6) then
         q0 = U_q_shifts(i1,:)*1.E-4/v_norm(U_q_shifts(i1,:))
         call msg("s","Renormalizing shift to :",v_norm(q0))
       else 
         q0 = U_q_shifts(i1,:)
       endif
       !
       i2 = i2 + 1
       !
       ! Convert U_q_shifts[iku] to coo_out
       !
       call K_transform(q0,trim(coo_in))
       !
       ! Print shift vector in new basis
       !
       write(ch,'(a,i1,a,3(f13.10,a),2a)') "Shifted K-grid[",i2,"]: {K} + (",q0(1),",",q0(2),",",q0(3),") [",trim(coo_out),"]"
       call section("=",trim(ch))
       !
       ! Apply the shift
       !
       do ik=1,Xk%nibz
         !
         k_plus_q(ik,:)=kout(ik,:)+q0(:)
         !
         if (v_norm(k_plus_q(ik,:)).lt.1.E-7) k_plus_q(ik,:) = 0.0_SP ! remove SP error k points
         !
       enddo
       !
       call k_grid_printer(Xk%nibz,k_plus_q,trim(coo_out),"PW")
       !
     enddo
      !
   endif
   !
 end if
 !
 if (l_high_sym_pts) then
   !
   call k_special( )
   !
   call k_circuit_driver(PtsPath,N_path_pts)
   !
 endif
 !
 if(l_regular_grid) then
   if(product(INTERP_grid)<=0) call error("Nk_x Nk_y Nk_z must be positive ")
   call INTERPOLATION_grid(USER_k,'i',INTERP_grid)
   call k_ibz2bz(USER_k,'a',.false.)
   USER_k%ptbz(:,:)=USER_k%ptbz(:,:)-nint(USER_k%ptbz(:,:)) ! Bring in the BZ
   USER_k%pt(:,:)  =USER_k%pt(:,:)  -nint(USER_k%pt(:,:)) ! Bring in the BZ
   if(l_shift01) then
     where(USER_k%ptbz(:,:)<-zero_dfl)
       USER_k%ptbz(:,:)=1._SP+USER_k%ptbz(:,:)
     end where
     where(USER_k%pt(:,:)<-zero_dfl)
       USER_k%pt(:,:)=1._SP+USER_k%pt(:,:)
     end where
   endif
   !
   call msg("s",'K-points (BZ) PW-formatted')
   call k_grid_printer(USER_k%nbz,USER_k%ptbz,'rlu',"PW")
   call msg("s",'K-points (IBZ) PW-formatted')
   call k_grid_printer(USER_k%nibz,USER_k%pt,'rlu',"PW")
   !
   goto 1
   !
 endif
 !
 if (List_internal_points) then
   !
   write(ch,'(a,a,a)') "== Built-in grids [",trim(coo_out),"] =="
   !==================================================================
   call section('+',trim(ch))
   !
   if (l_k_grid) then
     if (Expand_internal_points) then
       call k_expand(Xk)
       call msg("s",'K-points (BZ) formatted')
       call k_grid_printer(Xk%nbz,Xk%ptbz,'iku',"Y")
       call msg("s",'K-points (BZ) un-formatted')
       call k_grid_printer(Xk%nbz,Xk%ptbz,'iku',"none")
       call msg("s",'K-points (BZ) PW-formatted')
       call k_grid_printer(Xk%nbz,Xk%ptbz,'iku',"PW")
     else
       YAMBO_ALLOC(BZ_weights,(Xk%nibz))
       BZ_weights(:)=Xk%weights(:)*real(Xk%nbz,SP)
       call msg("s",'K-points (IBZ) formatted')
       call k_grid_printer(Xk%nibz,Xk%pt,'iku',"Y",wk=BZ_weights)
       call msg("s",'K-points (IBZ) un-formatted')
       call k_grid_printer(Xk%nibz,Xk%pt,'iku',"none",wk=BZ_weights)
       call msg("s",'K-points (BZ) PW-formatted')
       call k_grid_printer(Xk%nibz,Xk%pt,'iku',"PW",wk=BZ_weights)
       YAMBO_FREE(BZ_weights)
     endif
   endif
   !
   if (l_q_grid) then
     if (Expand_internal_points) then
       call k_expand(q)
       call msg("s",'Q-points (BZ) formatted')
       call k_grid_printer(q%nbz,q%ptbz,'iku',"Y")
       call msg("s",'Q-points (BZ) un-formatted')
       call k_grid_printer(q%nbz,q%ptbz,'iku',"none")
       call msg("s",'Q-points (BZ) PW-formatted')
       call k_grid_printer(q%nbz,q%ptbz,'iku',"PW")
     else
       YAMBO_ALLOC(BZ_weights,(q%nibz))
       BZ_weights(:)=q%weights(:)*real(q%nbz,SP)
       call msg("s",'Q-points (IBZ) formatted')
       call k_grid_printer(q%nibz,q%pt,'iku',"Y",wk=BZ_weights)
       call msg("s",'Q-points (IBZ) un-formatted')
       call k_grid_printer(q%nibz,q%pt,'iku',"none",wk=BZ_weights)
       call msg("s",'Q-points (IBZ) PW-formatted')
       call k_grid_printer(q%nibz,q%pt,'iku',"PW",wk=BZ_weights)
       YAMBO_FREE(BZ_weights)
     endif
   endif
   !
   goto 1
   !
 endif
 !
 live_timing_is_on=.true.
 !
1 continue
 !
 ! CLEAN
 !
 call bz_samp_reset(USER_k)
 call bz_samp_reset(Q_user)
 call k_ibz2bz(q,'d',.false.)
 call k_ibz2bz(Xk,'d',.false.)
 !
 contains
   !
   subroutine reduce_and_print(K_grid,GW_grid)
     !
     use zeros,   ONLY:k_rlu_zero,define_zeros
     use stderr,  ONLY:intc
     type(bz_samp)   ::K_grid
     logical         ::GW_grid
     !
     ! Work Space
     !
     real(SP), allocatable :: GWK_table(:)
     !
     ! Before doing any operation I need to redefine the zeros module
     ! components. This is beacuse K_grid contains additional points (the GW ones, for example)
     !
     call define_zeros(vector_=K_grid%ptbz,zero_=k_rlu_zero,RLU=.TRUE.)
     !
     call msg("s","Reducing & Expanding the "//trim(intc(K_grid%nbz))//" k-points ...")
     call k_reduce(K_grid,.true.)
     call k_expand(K_grid)
     !
     call msg("l","done")
     call msg("s","Reduced K-grid points:",K_grid%nibz)
     !
     ! the K_grid contains the final grid.
     ! When this contains the {k}+k_gw-{q} grids
     ! before reporting the points I want to sign the 
     ! position in the final grid of the given QP k-points 
     ! (read from the input file)
     !
     YAMBO_ALLOC(int_ctl,(K_grid%nibz))
     int_ctl=0
     !
     ! int_ctl =0 -> nothing
     ! int_ctl/=0 -> GW (index)
     !
     if (GW_grid) then
       do i1=1,K_grid%nibz
         call c2a(v_in=K_grid%pt(i1,:),v_out=ka,mode='ki2a')
         do i2=1,USER_k%nbz
           call c2a(v_in=USER_k%pt(i2,:),v_out=v,mode="ki2a")
           do is=1,nsym
             if (rlu_v_is_zero(v(:)-ka(:))) then
               if (int_ctl(i1)/=0) cycle
               int_ctl(i1)=i2
               exit
             endif
           enddo
         enddo
       enddo
     else
       int_ctl=0
     endif
     !
     YAMBO_ALLOC(GWK_table,(USER_k%nbz))
     !
     live_timing_is_on=.false.
     call msg("s","       ---------- Reduced K-grid ----------")
     GWK_table=0
     do i2=0,USER_k%nbz
       do i1=1,K_grid%nibz
         if (int_ctl(i1)/=i2.and..not.( i2==0.and.int_ctl(i1)/=0.and.i1<nXkibz )) cycle
         if (int_ctl(i1)/=0) then
            if(GWK_table(int_ctl(i1))>0) cycle
            GWK_table(int_ctl(i1))=1
         endif
         !
         call K_transform(K_grid%pt(i1,:),'iku')
         !
         if(No_Weight) then
           if (int_ctl(i1)==0) write (ch,'(3f12.7)') K_grid%pt(i1,:)
           if (int_ctl(i1)/=0) write (ch,'(3f12.7,i3)') K_grid%pt(i1,:),int_ctl(i1)
         else
           if (int_ctl(i1)==0) write (ch,'(4f12.7)') K_grid%pt(i1,:),K_grid%weights(i1)
           if (int_ctl(i1)/=0) write (ch,'(4f12.7,i3)') K_grid%pt(i1,:),&
&                                                     K_grid%weights(i1),int_ctl(i1)
         endif
         call msg("s",trim(ch))
       enddo
     enddo
     live_timing_is_on=.true.
     !
     YAMBO_FREE(GWK_table)
     YAMBO_FREE(int_ctl)
     !
   end subroutine
   !
end subroutine
