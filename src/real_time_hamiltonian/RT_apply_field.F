!
! License-Identifier: GPL
!
! Copyright (C) 2010 The Yambo Team
!
! Authors (see AUTHORS file for details): AM CA DS
!
subroutine RT_apply_field(ik,i_sp_pol,H_nl_sc,A_input)
 !
 ! This subroutine add the vector potential to the Hamiltonian
 !
 use pars,           ONLY:SP,cI
 use drivers,        ONLY:l_nl_optics
 use electrons,      ONLY:n_spinor
 use units,          ONLY:SPEED_OF_LIGHT
 use fields,         ONLY:gauge_field
 use DIPOLES,        ONLY:DIP_iR,DIP_P,DIP_v,DIP_P_spinor
 use X_m,            ONLY:global_gauge
 use real_time,      ONLY:RT_bands,l_update_SOC,l_velocity_gauge_corr
 use vec_operate,    ONLY:cross_product
 use D_lattice,      ONLY:sigma_x,sigma_y,sigma_z
 use parallel_m,     ONLY:PAR_Xk_ibz_index
 !
 implicit none
 !
 integer,           intent(in) :: ik,i_sp_pol
 complex(SP),    intent(inout) :: H_nl_sc(RT_bands(1):RT_bands(2),RT_bands(1):RT_bands(2))
 type(gauge_field), intent(in) :: A_input
 ! 
 ! Work Space
 !
 integer     :: i1,i2,i_spinor,j_spinor
 complex(SP) :: A_dot_J(2),E_dot_P(2),&
&               ExJs(3,n_spinor,n_spinor),A_vec_pot(3),E_vec_pot(3),&
&               POL(3),Jtot(3),Jspinor(3,n_spinor,n_spinor),sigma_dot_ExJs
 !
 select case(trim(global_gauge))
 !
 case('velocity')
   !
   ! Velocity gauge
   !================
   !
   A_vec_pot=A_input%vecpot/SPEED_OF_LIGHT ! Gaussian units
   !
   ! Add the vector potential to the singular part of the Self-Energy
   ! in case of complex Efield puts the correct Hermitian matrix
   ! and the intraband optical matrix elements have been set to zero
   !
!$omp parallel do default(shared), private(i1,Jtot,A_dot_J)
   do i1=RT_bands(1),RT_bands(2)
     !
     ! Case A: the pseudo-potential is correctly gauged, so A*p enters the hamiltonian
     if (     l_velocity_gauge_corr) Jtot(1:3)= -DIP_P(1:3,i1,i1,PAR_Xk_ibz_index(ik),i_sp_pol)
     ! Case b: the pseudo-potential is   *not*   gauged, so A*v enters the hamiltonian.
     !         A*v accounts for the linear expansion of the gauge correction to the non local potential
     if (.not.l_velocity_gauge_corr) Jtot(1:3)= -DIP_v(1:3,i1,i1,PAR_Xk_ibz_index(ik),i_sp_pol)
     !	   
     if(.not.l_nl_optics) then
       !	   
       ! This term can be set to zero with the freedom of the velocity gauge
       ! see discussion in sec. 4.A of  J. of Optical Society B, 39(2), 535 (2022)
       ! we do not include it in yambo_nl
       !	   
       Jtot(1:3)=Jtot(1:3)+A_vec_pot(1:3)/2._SP
       !	   
     endif
     !
     A_dot_J(1)=sum(A_vec_pot(1:3)*      Jtot(1:3) )
     A_dot_J(2)=sum(A_vec_pot(1:3)*conjg(Jtot(1:3)))
     H_nl_sc(i1,i1)=H_nl_sc(i1,i1)-real(A_dot_J(1),SP)
   enddo
!$omp end parallel do
   !
!$omp parallel do default(shared), private(i1,i2,Jtot,A_dot_J)
   do i1=RT_bands(1),RT_bands(2)
     do i2=i1+1,RT_bands(2)
       !
       ! Note that (P)_cv = DIP_P (:,iv,ic,:,:)
       ! the off-diagonal diamagnetic term is zero
       !     
       !Jtot(:)=-DIP_P(1:3,i2,i1,PAR_Xk_ibz_index(ik),i_sp_pol) 
       ! Case A: the pseudo-potential is correctly gauged, so A*p enters the hamiltonian
       if (     l_velocity_gauge_corr) Jtot(:)= -DIP_P(1:3,i2,i1,PAR_Xk_ibz_index(ik),i_sp_pol)
       ! Case b: the pseudo-potential is   *not*   gauged, so A*v enters the hamiltonian.
       !         A*v accounts for the linear expansion of the gauge correction to the non local potential
       if (.not.l_velocity_gauge_corr) Jtot(:)= -DIP_v(1:3,i2,i1,PAR_Xk_ibz_index(ik),i_sp_pol)
       A_dot_J(1)=sum(A_vec_pot(1:3)*      Jtot(1:3) )
       A_dot_J(2)=sum(A_vec_pot(1:3)*conjg(Jtot(1:3)))
       H_nl_sc(i1,i2)=H_nl_sc(i1,i2)-A_dot_J(1)
       H_nl_sc(i2,i1)=H_nl_sc(i2,i1)-A_dot_J(2)
       !
     enddo
   enddo
!$omp end parallel do
   !
 case('length')
   !
   ! Length gauge
   !===============
   !
   E_vec_pot=-A_input%vecpot_vel/SPEED_OF_LIGHT ! Gaussian units
   !
   ! Intraband terms
   !
   do i1=RT_bands(1),RT_bands(2)
     POL(:)=cI*DIP_iR(:,i1,i1,PAR_Xk_ibz_index(ik),i_sp_pol)
     E_dot_P(1)=sum(E_vec_pot(1:3)*      POL(1:3) )
     E_dot_P(2)=sum(E_vec_pot(1:3)*conjg(POL(1:3)))
     H_nl_sc(i1,i1)=H_nl_sc(i1,i1)-real(E_dot_P(1),SP)
     !
     if(l_update_SOC.and.n_spinor==2) then
       !
       ! Updte of SOC interaction, i.e. last term in 
       ! Eq.(1) JCTC 11, 4870 (2015)
       !
       Jspinor(:,:,:)=DIP_P_spinor(:,:,:,i1,i1,PAR_Xk_ibz_index(ik))
       do i_spinor=1,n_spinor
         do j_spinor=1,n_spinor
           ExJs(:,i_spinor,j_spinor)=cross_product(conjg(E_vec_pot),Jspinor(:,i_spinor,j_spinor))
         enddo
       enddo
       sigma_dot_ExJs=sum(sigma_x*ExJS(1,:,:))+ &
&                     sum(sigma_y*ExJS(2,:,:))+ &
&                     sum(sigma_z*ExJS(3,:,:))
       H_nl_sc(i1,i1)=H_nl_sc(i1,i1)-0.25_SP*real(sigma_dot_ExJs)
     endif
     !
   enddo
   !
   ! Inter-band terms
   !
   do i1=RT_bands(1),RT_bands(2)
     do i2=i1+1,RT_bands(2)
       !
       ! Note that (iR)_cv = DIP_iR (:,iv,ic,:,:)
       !
       POL(:)=cI*DIP_iR(:,i2,i1,PAR_Xk_ibz_index(ik),i_sp_pol)
       E_dot_P(1)=sum(E_vec_pot(1:3)*      POL(1:3) )
       E_dot_P(2)=sum(E_vec_pot(1:3)*conjg(POL(1:3)))
       H_nl_sc(i1,i2)=H_nl_sc(i1,i2)-E_dot_P(1)
       H_nl_sc(i2,i1)=H_nl_sc(i2,i1)-E_dot_P(2)
       !
       if(l_update_SOC.and.n_spinor==2) then
         Jspinor(:,:,:)=DIP_P_spinor(:,:,:,i2,i1,PAR_Xk_ibz_index(ik))
         do i_spinor=1,n_spinor
           do j_spinor=1,n_spinor
             ExJs(:,i_spinor,j_spinor)=cross_product(E_vec_pot,Jspinor(:,i_spinor,j_spinor))
           enddo
         enddo
         sigma_dot_ExJs=sum(sigma_x(:,:)*ExJs(1,:,:))+ &
&                       sum(sigma_y(:,:)*ExJs(2,:,:))+ &
&                       sum(sigma_z(:,:)*ExJs(3,:,:))
         H_nl_sc(i1,i2)=H_nl_sc(i1,i2)-      sigma_dot_ExJs
         H_nl_sc(i2,i1)=H_nl_sc(i2,i1)-conjg(sigma_dot_ExJs)
       endif
       !
     enddo
   enddo
   !
 end select
 !
end subroutine RT_apply_field
