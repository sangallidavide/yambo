!
! License-Identifier: GPL
!
! Copyright (C) 2025 The Yambo Team
!
! Authors (see AUTHORS file for details): DS
!
subroutine RT_gradk_rho_times_E(ik,i_sp_pol,H_nl_sc,dG_in,A_input,Xk)
 !
 ! This subroutine computes the k-gradient of the density matrix,
 ! and multiplies it by the electric field.
 ! See Phys. Rev. B 76, 035213 (2007)
 !
 use pars,           ONLY:SP,cI,cZERO,pi
 use units,          ONLY:SPEED_OF_LIGHT
 use real_time,      ONLY:l_length_grad_k,RT_bands
 use fields,         ONLY:gauge_field
 use real_time,      ONLY:RT_bands,NE_i_time,NE_i_last_field
 use vec_operate,    ONLY:k_periodic_idx!,degeneration_finder
 use matrix_operate, ONLY:hermitian
 use DIPOLES,        ONLY:DIP_S,DIP_iR,l_force_SndOrd
 use R_lattice,      ONLY:bz_samp,k_map
 use D_lattice,      ONLY:a
 use parallel_m,     ONLY:PAR_G_k_range
 !
 implicit none
 !
 integer,           intent(in) :: ik,i_sp_pol
 complex(SP),    intent(inout) :: H_nl_sc(RT_bands(1):RT_bands(2),RT_bands(1):RT_bands(2))
 complex(SP),       intent(in) :: dG_in(RT_bands(1):RT_bands(2),RT_bands(1):RT_bands(2),PAR_G_k_range(1):PAR_G_k_range(2))
 type(gauge_field), intent(in) :: A_input
 type(bz_samp),     intent(in) :: Xk
 !
 ! Workspace
 !
 integer     :: id_red,max_step,istep,idx(3),idir,ib1,ib2,jb1,jb2,nb(2),&
                ikm1,ikp1,ikbz,ikm1_ibz,ikm1_is,ikp1_ibz,ikp1_is
 real(SP)    :: q_fac,E_field(3)
 complex(SP) :: Dip_kp_c(3),Dip_km_c(3)
 !
 !complex(SP) :: Mp_kmq(RT_bands(2),RT_bands(2),RT_bands(2),RT_bands(2)),&
 !&              Mm_kpq(RT_bands(2),RT_bands(2),RT_bands(2),RT_bands(2)),&
 !&              Sp_kmq(RT_bands(2),RT_bands(2)),Sm_kpq(RT_bands(2),RT_bands(2)),&
 !&              A1p_kmq(RT_bands(2),RT_bands(2)),A1m_kpq(RT_bands(2),RT_bands(2)),&
 !&              A2p_kmq(RT_bands(2),RT_bands(2)),A2m_kpq(RT_bands(2),RT_bands(2)),&
 !&              DIP_kp(RT_bands(2),RT_bands(2)),DIP_km(RT_bands(2),RT_bands(2))
 !
 complex(SP), allocatable :: Mp_kmq(:,:,:,:),Mm_kpq(:,:,:,:),&
 &              Sp_kmq(:,:),Sm_kpq(:,:),&
 &              A1p_kmq(:,:),A1m_kpq(:,:),&
 &              A2p_kmq(:,:),A2m_kpq(:,:),&
 &              DIP_kp(:,:),DIP_km(:,:)
 !
 complex(SP) :: gradk_rho(RT_Bands(1):RT_bands(2),RT_bands(1):RT_bands(2),3),&
 &          gradk_rho_red(RT_Bands(1):RT_bands(2),RT_Bands(1):RT_bands(2),2)
 !
 ! Warning: kpts parallelization to be fixed, since rho is distributed
 !
 if (.not.l_length_grad_k) return
 !
 ! Length gauge
 !===============
 E_field=-A_input%vecpot_vel/SPEED_OF_LIGHT ! Gaussian units
 !
 if (NE_i_time>NE_i_last_field) return
 !
 nb=RT_bands
 !
 allocate(Mp_kmq(nb(1):nb(2),nb(1):nb(2),nb(1):nb(2),nb(1):nb(2)))
 allocate(Mm_kpq(nb(1):nb(2),nb(1):nb(2),nb(1):nb(2),nb(1):nb(2)))
 !
 allocate(Sp_kmq(nb(1):nb(2),nb(1):nb(2)))
 allocate(Sm_kpq(nb(1):nb(2),nb(1):nb(2)))
 !
 allocate(A1p_kmq(nb(1):nb(2),nb(1):nb(2)))
 allocate(A1m_kpq(nb(1):nb(2),nb(1):nb(2)))
 !
 allocate(A2p_kmq(nb(1):nb(2),nb(1):nb(2)))
 allocate(A2m_kpq(nb(1):nb(2),nb(1):nb(2)))

 allocate(DIP_kp(nb(1):nb(2),nb(1):nb(2)))
 allocate(DIP_km(nb(1):nb(2),nb(1):nb(2)))
 !
 ikbz=Xk%k_table(ik,1)
 !
 gradk_rho=cZERO
 !
 do id_red=1,3  ! loop on reciprocal lattice
   !
   if(k_map%max_kdir(id_red)==1) cycle
   !
   max_step=1
   !
   ! I need at leat 5 k-points in the line to
   ! use the second order formula
   !
   if(k_map%max_kdir(id_red)>=6.and..not.l_force_SndOrd) max_step=2
   !
   gradk_rho_red=cZERO
   !
   do istep=1,max_step 
     !
     idx=k_map%k_map_inv(ikbz,:)
     idx(id_red)=idx(id_red)-istep
     idx=k_periodic_idx(idx,k_map)
     ikm1=k_map%k_map_dir(idx(1),idx(2),idx(3))
     !
     ikm1_ibz=Xk%sstar(ikm1,1)
     ikm1_is =Xk%sstar(ikm1,2)
     !
     idx=k_map%k_map_inv(ikbz,:)
     idx(id_red)=idx(id_red)+istep
     idx=k_periodic_idx(idx,k_map)
     ikp1=k_map%k_map_dir(idx(1),idx(2),idx(3))
     !
     ikp1_ibz=Xk%sstar(ikp1,1)
     ikp1_is =Xk%sstar(ikp1,2)
     !
     !S_tmp(:,:,1)=          DIP_S(:,:,id_red+(istep-1)*3,ikbz,i_sp_pol)    ! g+(k)
     !S_tmp(:,:,2)=hermitian(DIP_S(:,:,id_red+(istep-1)*3,ikm1,i_sp_pol))   ! g-(k)=g+(k-q)
     !
     !S_kpq(:,:,1)=          DIP_S(:,:,id_red+(istep-1)*3,ikp1,i_sp_pol)    ! g+(k+q)
     Sm_kpq(:,:)=hermitian(DIP_S(nb(1):,nb(1):,id_red+(istep-1)*3,ikbz,i_sp_pol))   ! g-(k+q)=g+(k)
     !
     Sp_kmq(:,:)=          DIP_S(nb(1):,nb(1):,id_red+(istep-1)*3,ikm1,i_sp_pol)    ! g+(k-q)
     !S_kmq(:,:,2)=hermitian(DIP_S(:,:,id_red+(istep-1)*3,ikm2,i_sp_pol))   ! g-(k-q)=g+(k-2*q)
     !
     ! Option 1: Do SVD and later dagger.
     ! Coding to be imported from DIPOLE_build_derk subroutine in branch phys-dipoles-k-derivative
     !
     ! Option 2: Directly do matrix inversion in place of dagger.
     !
     ! Option 3: Use the dagger without SVD.
     do ib1=RT_bands(1),RT_bands(2)
       do ib2=RT_bands(1),RT_bands(2)
         Mp_kmq(ib1,ib2,:,:)=conjg(Sp_kmq(ib2,ib1))*Sp_kmq(:,:)
         Mm_kpq(ib1,ib2,:,:)=conjg(Sm_kpq(ib2,ib1))*Sm_kpq(:,:)
       enddo
     enddo
     !
     ! I need to rotate the dipoles in k-space,
     ! and then project them along the id_red direction
     ! Note: this term is quitre slow.
     ! Pre-expanding the dipoles with the old
     ! DIPOLE expand wuold speed-up this part
     !
     DIP_km=cZERO
     DIP_kp=cZERO
     do ib1=RT_bands(1),RT_bands(2)
       do jb1=RT_bands(1),RT_bands(2)
         call DIPOLE_rotate(ib1,jb1,ikm1,i_sp_pol,"DIP_iR",Xk,DIP_km_c)
         call DIPOLE_rotate(ib1,jb1,ikp1,i_sp_pol,"DIP_iR",Xk,DIP_kp_c)
	 ! DEBUG lines to verify the time gain when non-expanding on the fly
         !DIP_km_c=DIP_iR(:,ib1,jb1,ikm1_ibz,i_sp_pol)
         !DIP_kp_c=DIP_iR(:,ib1,jb1,ikp1_ibz,i_sp_pol)
         do idir=1,3
           DIP_km(ib1,jb1)=DIP_km(ib1,jb1)+cI*DIP_km_c(idir)*a(id_red,idir)
           DIP_kp(ib1,jb1)=DIP_kp(ib1,jb1)+cI*DIP_kp_c(idir)*a(id_red,idir)
         enddo
       enddo
     enddo
     !
     !
     A1p_kmq=matmul(-cI*dG_in(:,:,ikm1_ibz),hermitian(DIP_km))
     A1m_kpq=matmul(-cI*dG_in(:,:,ikp1_ibz),hermitian(DIP_kp))
     !
     A2p_kmq=matmul(-cI*dG_in(:,:,ikm1_ibz),          DIP_km )
     A2m_kpq=matmul(-cI*dG_in(:,:,ikp1_ibz),          DIP_kp )
     !
     ! 1/q factor = k_map%max_kdir(id_red)/(4._SP*pi*real(istep,SP))
     q_fac=k_map%max_kdir(id_red)/(4._SP*pi*real(istep,SP))/2._SP
     !
     do ib1=RT_bands(1),RT_bands(2)
       do jb1=RT_bands(1),RT_bands(2)
         do ib2=RT_bands(1),RT_bands(2)
           do jb2=RT_bands(1),RT_bands(2)
             gradk_rho_red(ib1,jb1,istep)=gradk_rho_red(ib1,jb1,istep)+ &
             &  -cI*Mp_kmq(ib2,ib1,jb2,jb1)*dG_in(ib2,jb2,ikp1_ibz)  &
             &  +cI*Mm_kpq(ib2,ib1,jb2,jb1)*dG_in(ib2,jb2,ikm1_ibz)  &
             &  -cI/q_fac*Mp_kmq(ib2,ib1,jb2,jb1)*(A1p_kmq(ib2,ib1)-A2p_kmq(ib2,ib1)) &
             &  +cI/q_fac*Mm_kpq(ib2,ib1,jb2,jb1)*(A1m_kpq(ib2,ib1)-A2m_kpq(ib2,ib1))
           enddo
         enddo
       enddo
     enddo
     !
     gradk_rho_red(:,:,istep)=gradk_rho_red(:,:,istep)*q_fac/2._SP
     !
   enddo ! istep
   !
   ! Rotate in cartesian coordinates
   !
   if(max_step==1) then
     !      
     ! First order formula D(dk)      
     !
     do idir=1,3
       gradk_rho(:,:,idir)=gradk_rho(:,:,idir)+transpose(gradk_rho_red(:,:,1))*a(id_red,idir)
     enddo
     !
   else
     !      
     ! Second order formula   ( 4*D(dk) - D(2*dk) ) / 3
     !
     do idir=1,3
       gradk_rho(:,:,idir)=gradk_rho(:,:,idir)+transpose(4._SP*gradk_rho_red(:,:,1)-gradk_rho_red(:,:,2))/3._SP*a(id_red,idir)
     enddo
     !
   endif
   !
 enddo   ! id_red
 !
 do idir=1,3
   H_nl_sc(:,:)=H_nl_sc(:,:)-gradk_rho(:,:,idir)*E_field(idir)
 enddo
 !
end subroutine RT_gradk_rho_times_E
 
