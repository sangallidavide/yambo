!
!        Copyright (C) 2000-2017 the YAMBO team
!              http://www.yambo-code.org
!
! Authors (see AUTHORS file for details): MP DS
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
subroutine excitons_lifetimes(Dip,Xen,Xk,Xq,BS_R,BS_E,BS_E_degs)
 !
 use pars,          ONLY:SP,DP,schlen,pi,cZERO
 use units,         ONLY:HA2EV,AU2SEC,SPEED_OF_LIGHT
 use electrons,     ONLY:spin_occ,levels,spin
 use stderr,        ONLY:intc,real2ch
 use BS_solvers,    ONLY:BSS_eh_table,BS_mat,BSS_n_eig
 use BS,            ONLY:BS_bands,BS_H_dim
 use com,           ONLY:msg,of_open_close
 use R_lattice,     ONLY:q0_def_norm,nkbz,bz_samp,q_norm
 use D_lattice,     ONLY:DL_vol,non_periodic_directions,l_1D,l_2D,l_0D,l_3D
 use X_m,           ONLY:eps_2_alpha
 use DIPOLES,       ONLY:DIPOLE_t,DIP_iR
 use YPPm,          ONLY:EXC_temp,EXC_mass,EXC_Phi_Grid,EXC_shift,epsi0,&
&                        EXCITONS_user_indexes,EXCITONS_n_user_states
 use interfaces,    ONLY:PARALLEL_global_indexes
 use vec_operate,   ONLY:v_norm
 use IO_m,          ONLY:io_DIP
 !HYC  add
 use parallel_m,     ONLY:PAR_DIPk_ibz_index
 !
 ! implicit none 
 !
#include<memory.h>
 !
 type(levels)       :: Xen
 type(DIPOLE_t)     :: Dip
 type(bz_samp)      :: Xk,Xq
 !
 complex(SP)        :: BS_E(BSS_n_eig), BS_R(BSS_n_eig) 
 integer            :: BS_E_degs(BSS_n_eig)
 !
 ! Work Space
 !
 character(schlen)  :: titles(7), ch_dummy
 !
 integer            :: i_lambda,i_grid,io_err
 integer            :: neh,i_l,i_l_grp,n_lambda_deg,i_v,i_c,i_kbz,i_kibz,i_spin 
 integer            :: Phi_grid_size
 !
 complex(SP)        :: Ir_EXC(3),rotated_iR(3),P_EXC(3)
 !
 real(SP)           :: P_EXC_square,Px_EXC_square,Py_EXC_square,Pz_EXC_square
 real(SP)           :: gamma0(BSS_n_eig),gamma0_proj(BSS_n_eig)
 real(SP)           :: gamma_exp, gamma0_grp,f_cv,alpha,tau_exp,tau0_grp
 !
 !
 real(SP)           :: phi(EXC_Phi_Grid),exc_energy
 real(SP)           :: phi0(BSS_n_eig),vdummy(3),tau0
 real(SP)           :: gamma_thermal(EXC_Phi_Grid),tau_thermal(EXC_Phi_Grid)
 real(SP)           :: alpha1(EXC_Phi_Grid)
 real(SP)           :: num,den,tau_eff,gamma_eff
 real(DP)           :: arg
 !
 !HYC add 
 integer :: ik_mem
 complex(SP):: iR(3) , iR_EXC_TRY(3)
 !
 !VAJ add - trapezoidal integration
 integer, parameter :: n_trap=1000
 real(SP)           :: u_trap,du_trap,end_trap,start_trap,integrand,A,B,C,D,E,F
 real(SP)           :: integral
 integer            :: i_trap
 !
 call DIPOLE_dimensions(Xen,Dip,BS_bands,(/0.,0.,0./))
 !
 call PARALLEL_global_indexes(Xen,Xk,Xq," ",RESET=.TRUE.)
 call PARALLEL_global_indexes(Xen,Xk,Xq,"DIPOLES",Dip=Dip)
 !
 call DIPOLE_IO(Xk,Xen,Dip,'read ',io_err,'DIP')
 !
 titles(1)='phi [deg]'
 titles(2)='Gamma_Therm'
 titles(3)='Tau_Therm'
  Phi_Grid_size = EXC_Phi_Grid
 !
 call section('=','Computing Excitons radiative lifetimes')
 !
 if (trim(non_periodic_directions)/="none")  call check_periodic_directions()
 !
 do i_grid=1,Phi_grid_size
   phi(i_grid)=2._SP*pi*real(i_grid-1,SP)/real(Phi_grid_size-1,SP)
 enddo
 !
 do i_lambda=1,EXCITONS_n_user_states
   !
   gamma_thermal(:)=0._SP
   gamma_exp =0._SP
   gamma0_grp=0._SP
   tau_thermal(:)=0._SP
   tau_exp =0._SP
   tau0_grp=0._SP
   num = 0._SP
   den = 0._SP
   alpha = 0._SP
   alpha1(:) = 0._SP
   !
   i_l_grp=EXCITONS_user_indexes(i_lambda)
   !
   if(i_lambda/=BS_E_degs(i_l_grp)) cycle
   !
   n_lambda_deg=count(BS_E_degs==BS_E_degs(i_l_grp))
   !
   if (n_lambda_deg>1) call msg('s',':: State '//trim(intc(i_l_grp))//' Merged with states '//&
&                                     trim(intc(BS_E_degs(i_l_grp)))//' -> '//&
&                                     trim(intc(BS_E_degs(i_l_grp)+n_lambda_deg-1)))
   !
   !
   ch_dummy='exc_lifetimes_at_'//trim(intc(i_l_grp))
   call of_open_close(ch_dummy,'ot')
   !
   ! First compute the excitonic dipoles <mu>  starting from indipendent-particle
   ! ones
   !==========================================
   !
   do i_l=BS_E_degs(i_l_grp),BS_E_degs(i_l_grp)+n_lambda_deg-1
     iR_EXC = cZERO
     iR_EXC_TRY= cZERO
     do neh = 1,BS_H_dim !(size BSE matrix, cvk)
       !
       i_kbz   = BSS_eh_table(neh,1)
       i_kibz  = Xk%sstar(i_kbz,1)
       i_v     = BSS_eh_table(neh,2)
       i_c     = BSS_eh_table(neh,3)
       i_spin  = spin(BSS_eh_table(neh,:))
       !
       f_cv=      (Xen%f(i_v,i_kibz,i_spin)- Xen%f(i_c,i_kibz,i_spin))/spin_occ
       !
       call GET_the_rotated_iR()
       !
       iR_EXC(:)=iR_EXC(:) - rotated_iR(:)*BS_mat(neh,i_l)*sqrt(f_cv)
       !
     enddo
 !    do neh = 1,BS_H_dim !(size BSE matrix, cvk)
       !
 !      i_kbz   = BSS_eh_table(neh,1)
 !      i_kibz  = Xk%sstar(i_kbz,1)
 !      i_v     = BSS_eh_table(neh,2)
 !      i_c     = BSS_eh_table(neh,3)
 !      i_spin  = spin(BSS_eh_table(neh,:))
       !
 !      f_cv=      (Xen%f(i_v,i_kibz,i_spin)- Xen%f(i_c,i_kibz,i_spin))/spin_occ
       !
       !call GET_the_rotated_iR()
 !      ik_mem    = PAR_DIPk_ibz_index(i_kibz)
 !      iR=DIP_iR(:,i_c,i_v,ik_mem,i_spin)
       !
 !      iR_EXC_TRY(:)=iR_EXC_TRY(:) - iR(:)*BS_mat(neh,i_l)*sqrt(f_cv)

       !

 !    enddo
     !OPEN(UNIT=10, FILE='dipole', STATUS='UNKNOWN', ACCESS='APPEND') 
     !WRITE(10,*) "i_l_grp=",i_l_grp, "BS_E_degs(i_l_grp)", BS_E_degs(i_l_grp)
     !WRITE(10,*) iR_EXC(1)
     !WRITE(10,*) iR_EXC(2)
     !WRITE(10,*) iR_EXC(3)
     !WRITE(10,*) iR_EXC_TRY(1)
     !WRITE(10,*) iR_EXC_TRY(2)
     !WRITE(10,*) iR_EXC_TRY(3)
     !close(10)
     
    
     
     !
     ! Angle of the excitonic dipoles against the x direction
     !
     ! Warning define N degenerate angles phi(n)
     vdummy(:) = real(iR_exc(:)) / v_norm(real(iR_exc(:)))
!     vdummy(:) = aimag(iR_exc(:)) / v_norm(aimag(iR_exc(:)))
     phi0(i_l)=acos(vdummy(1))
!     phi0(i_l) = 0._SP
     !
!  exciton dipole modulus square |P_S|^2 
!        P_EXC_square =real(BS_E(i_l))**2*(real(iR_EXC(1)*conjg(iR_EXC(1)) + &
!&                 iR_EXC(2)*conjg(iR_EXC(2)) + &
!&                 0 ))
!VAJ add
        P_EXC_square =real(BS_E(i_l))**2*(real(iR_EXC(1)*conjg(iR_EXC(1)) + &
&                 iR_EXC(2)*conjg(iR_EXC(2)) + &
                  iR_EXC(3)*conjg(iR_EXC(3)) ))


!     P_EXC_square =real(BS_E(i_l))**2*(real(iR_EXC(1)*conjg(iR_EXC(1))  &
!&                 ))
!  exciton dipole P_Sx, P_Sy,P_Sz components

      P_EXC(:) =BS_E(i_l)*iR_EXC(:)
      Px_EXC_square = real(P_EXC(1)*conjg(P_EXC(1)))
      Py_EXC_square = real(P_EXC(2)*conjg(P_EXC(2)))
      Pz_EXC_square = real(P_EXC(3)*conjg(P_EXC(3))) !VAJ add
     
     OPEN(UNIT=10, FILE='momentum', STATUS='UNKNOWN', ACCESS='APPEND') 
     WRITE(10,*) "i_l_grp=",i_l_grp, "BS_E_degs(i_l_grp)", BS_E_degs(i_l_grp)
     WRITE(10,*) P_EXC(1), Px_EXC_square
     WRITE(10,*) P_EXC(2), Py_EXC_square
     WRITE(10,*) P_EXC(3), Pz_EXC_square
     WRITE(10,*) 'epsi0', epsi0
     WRITE(10,*) 'EXC_mass', EXC_mass

     close(10)
      
     ! compute the lifetime at q=0  gamma(0) and the prefactors for finite T
     !==========================================
     if (l_2D) then
      P_EXC_square = Px_EXC_square + Py_EXC_square 
      gamma0(i_l) = (eps_2_alpha/DL_vol)*P_EXC_square*4.*pi/((real(BS_E(i_l))-EXC_shift)*SPEED_OF_LIGHT*real(nkbz))
      gamma0_proj(i_l) = (eps_2_alpha/DL_vol)*real(BS_R(i_l)*conjg(BS_R(i_l)))* &
&       8._SP*pi*real(BS_E(i_l))/(q0_def_norm**2*SPEED_OF_LIGHT*real(nkbz))
      alpha=(real(BS_E(i_l))-EXC_shift)**2/(2._SP*(EXC_mass(1))*EXC_temp*SPEED_OF_LIGHT**2) 
      alpha1(:)=2._SP*sin(phi(:))**2+2._SP/3._SP*cos(phi(:))**2 *  &
&       (Px_EXC_square/P_EXC_square) + &
&         2._SP*cos(phi(:))**2+2._SP/3._SP*sin(phi(:))**2 *& 
&       (Py_EXC_square/P_EXC_square) + &
&       (-4._SP)/3._SP*cos(phi(:))*sin(phi(:))*  &
&        (conjg(P_EXC(1))*P_EXC(2) + conjg(P_EXC(2))*P_EXC(1))/P_EXC_SQUARE
      OPEN(UNIT=10, FILE='dipole', STATUS='UNKNOWN', ACCESS='APPEND')
      WRITE(10,*) 'gamma0', gamma0(i_l)
      WRITE(10,*) 'E^2/2McKbT', alpha
      WRITE(10,*) 'unit AU to Sec', AU2SEC
      close(10)
      
      OPEN(UNIT=7777, FILE='thermal_time', STATUS='UNKNOWN', ACCESS='APPEND')
      WRITE(7777,*) 'exciton number', i_l_grp
      WRITE(7777,*) 'energy(shifted)', (real(BS_E(i_l))-EXC_shift)*HA2EV
      WRITE(7777,*) 'energy(not shifted)', (real(BS_E(i_l)))*HA2EV
      WRITE(7777,*) 'temperature', EXC_temp
      WRITE(7777,*) 'gamma0', gamma0(i_l)/AU2SEC
      WRITE(7777,*) 'gamma_therm', gamma0(i_l)*alpha*(4._SP/3._SP)/AU2SEC
      WRITE(7777,*) 'tau0', (gamma0(i_l)/AU2SEC)**(-1)
      WRITE(7777,*) 'tau_therm', (gamma0(i_l)*alpha*(4._SP/3._SP)/AU2SEC)**(-1)
      WRITE(7777,*) '  '
      close(7777) 
     endif 
     if (l_1D) then 
      write(6,*) 'not yet implemented' 
      stop
      gamma0(i_l) = (eps_2_alpha/DL_vol)*P_EXC_square*2.*pi*real(nkbz)*SPEED_OF_LIGHT**2
      alpha= 4.*SQRT((real(BS_E(i_l))-EXC_shift)**2/(2._SP*EXC_mass(1)*EXC_temp*SPEED_OF_LIGHT**2))/3./SQRT(PI)
      alpha1(:)= 1.
     endif
     if (l_3D) then
      ! VAJ: check dielectric symmetry for isotropic system
      !if ((epsi0(1).eq.epsi0(2)).and.(epsi0(1).eq.epsi0(3))) then
      !   gamma0(i_l) =  (sqrt(epsi0(1))/DL_vol)*P_EXC_square*(8./3.)*sqrt(pi)/((real(BS_E(i_l))-EXC_shift)**2*real(nkbz))
      !   alpha=((real(BS_E(i_l))-EXC_shift)**2/(2._SP*EXC_mass(1)*EXC_temp*SPEED_OF_LIGHT**2))**1.5
      !   alpha1(:)=sin(phi(:)+phi0(i_l))**2
      !endif
      ! VAJ: check dielectric symmetry for unaxial system
      !if ((epsi0(1).eq.epsi0(2)).and.(epsi0(1).ne.epsi0(3))) then
         ! VAJ: Trapezoidal integration loop for P bar squared
         !end_trap = 1.
         !start_trap = -1.
         !A = 0.5*Px_EXC_square*epsi0(3)**2 + &
         !    0.5*Py_EXC_square*epsi0(3)**2 - &
         !        Pz_EXC_square*epsi0(1)**2
         !B = Pz_EXC_square*epsi0(1)**2
         !C = epsi0(1)
         !D = epsi0(3) - epsi0(1)
         !E = epsi0(1)**2
         !F = epsi0(3)**2 - epsi0(1)**2
         !integral = 0.0
         !do i_trap=0,n_trap
         !   u_trap = (i_trap*(end_trap-start_trap))/n_trap + start_trap
         !   integrand = (A*u_trap**2+B)/(sqrt(C+D*u_trap**2)*(E+F*u_trap**2))
         !   if ((i_trap.eq.0).or.(i_trap.eq.n_trap)) then
         !      integral = integral + integrand
         !   else
         !      integral = integral+2.*integrand
         !   endif
         !enddo
         !du_trap = (end_trap-start_trap)/n_trap
         !integral = sqrt(epsi0(3))*integral*du_trap/2.
         !gamma0(i_l) =  2.*(sqrt(pi*epsi0(1))/DL_vol)/((real(BS_E(i_l))-EXC_shift)**2*real(nkbz))*&
         !               (Px_EXC_square+Py_EXC_square+integral)
         gamma0(i_l) =  4.*pi*(sqrt(pi*epsi0(1))/DL_vol)/((real(BS_E(i_l))-EXC_shift)**2*real(nkbz))*&
                        ((2./3.*epsi0(3)/epsi0(1)+2)*(Px_EXC_square+Py_EXC_square)+8./3.*Pz_EXC_square)
         alpha=((real(BS_E(i_l))-EXC_shift)**2/(2._SP*(EXC_mass(1)**2*EXC_mass(2))**(1./3.)*EXC_temp*SPEED_OF_LIGHT**2))**1.5

         OPEN(UNIT=10, FILE='momentum', STATUS='UNKNOWN', ACCESS='APPEND')
         WRITE(10,*) 'gamma0',gamma0(i_l)
         WRITE(10,*) 'BvK SC volume',DL_vol*real(nkbz)
         WRITE(10,*) '(E^2/2Mc^2KbT)^(3/2)', alpha
         WRITE(10,*) 'unit AU to Sec', AU2SEC
         WRITE(10,*) 'eps2alpha', eps_2_alpha
         WRITE(10,*) 'epsi0', epsi0
         close(10)
     ! endif
      OPEN(UNIT=7777, FILE='thermal_time', STATUS='UNKNOWN', ACCESS='APPEND')
      WRITE(7777,*) 'exciton number', i_l_grp
      WRITE(7777,*) 'energy(shifted)', (real(BS_E(i_l))-EXC_shift)*HA2EV, 'eV'
      WRITE(7777,*) 'energy(not shifted)', (real(BS_E(i_l)))*HA2EV, 'eV'
      WRITE(7777,*) 'temperature', EXC_temp*HA2EV, 'eV'
      WRITE(7777,*) 'gamma0', gamma0(i_l)/AU2SEC, 's^-1'
      WRITE(7777,*) 'gamma_therm', gamma0(i_l)*alpha/AU2SEC, 's^-1'
      WRITE(7777,*) 'tau0', (gamma0(i_l)/AU2SEC)**(-1), 's'
      WRITE(7777,*) 'tau_therm', (gamma0(i_l)*alpha/AU2SEC)**(-1), 's'
      WRITE(7777,*) '  '
      close(7777)
     endif
     if (l_0D) gamma0(i_l) =   P_EXC_square*4.*0.333/(real(nkbz)*SPEED_OF_LIGHT**3)
     !
     gamma0_grp=gamma0_grp+gamma0(i_l)
     gamma_thermal(:)=gamma_thermal(:)+ gamma0(i_l)*alpha*alpha1(:)
!     gamma_thermal(:)=gamma_thermal(:)+ alpha1(:)
     if (l_1D) stop
     if (l_0D) stop 
     if (l_2D) gamma_exp=gamma_exp+gamma0(i_l)*alpha*(4._SP/3._SP)
     if (l_3D) gamma_exp=gamma_exp+gamma0(i_l)*alpha ! *(2._SP/3._SP)
     !
     exc_energy = (real(BS_E(i_l)-EXC_shift)*HA2EV)
     tau0_grp= AU2SEC/gamma0_grp
     tau_exp = AU2SEC/gamma_exp
!
     arg =real(BS_E(i_l)-BS_E(1))/EXC_temp
     num = num + dexp(-arg)*gamma_exp
     den = den + dexp(-arg)
     gamma_eff = num/den
     tau_eff = AU2SEC/gamma_eff
     write (9,*) i_l, BS_E(1)*HA2EV,BS_E(i_l)*HA2EV,gamma_eff,tau_eff
     write(10,111) i_l,phi0(i_l),exc_energy, alpha1(1),alpha1(5),alpha1(10),alpha1(20)
!     write(11,112) i_l,i_spin,exc_energy,real(BS_R(i_l)*conjg(BS_R(i_l)))/q0_def_norm**2, &
!&    real(iR_EXC(1)*conjg(iR_EXC(1)))
   enddo
   111 format (i4, 6d15.7)
! 112 format (2i4,3d15.7)
   !
   call msg('o exc_life',"# ","exc. E = "//trim(real2ch(exc_energy)),INDENT=0)
   call msg('o exc_life',"# ","gamma zero = "//trim(real2ch(gamma0_grp)),INDENT=0)
   call msg('o exc_life',"# ","gamma exp  = "//trim(real2ch(gamma_exp)),INDENT=0)
   !
   call msg('o exc_life',"# ","tau zero = "//trim(real2ch(tau0_grp)),INDENT=0)
   call msg('o exc_life',"# ","tau exp  = "//trim(real2ch(tau_exp)),INDENT=0)
   call msg('o exc_life',"# ","tau eff  = "//trim(real2ch(tau_eff)),INDENT=0)
!
   call msg('o exc_life','#',titles(:3),INDENT=0,USE_TABS=.true.)
!
   do i_grid=1,Phi_grid_size
    tau_thermal(i_grid)= AU2SEC/gamma_thermal(i_grid)
     call msg('o exc_life',' ',(/phi(i_grid)*180._SP/PI,gamma_thermal(i_grid),tau_thermal(i_grid)/),INDENT=0,USE_TABS=.true.)
   enddo
   !
  call of_open_close(ch_dummy)
   !
 enddo
 !
 contains
   !
   subroutine GET_the_rotated_iR( )
     !
     use D_lattice,      ONLY:i_time_rev,dl_sop,sop_inv,nsym
     use parallel_m,     ONLY:PAR_DIPk_ibz_index
     !
     implicit none
     !
     integer :: ik,is,ik_mem,i_dir
     logical :: t_rev
     real(SP):: field_dir(3,3),field_dir_rot(3)
     complex(SP):: iR(3)
     !
     ik=Xk%sstar(i_kbz,1)
     is=sop_inv(Xk%sstar(i_kbz,2))
     !
     if ( is<= nsym/(i_time_rev+1) ) t_rev=.false.
     if ( is > nsym/(i_time_rev+1) ) t_rev=.true.
     !
     ik_mem    = PAR_DIPk_ibz_index(ik)
     iR(:)=DIP_iR(:,i_c,i_v,ik_mem,i_spin)
     !
     ! OPTION 1
     !
     field_dir(:,1)=(/1._SP,0._SP,0._SP/)
     field_dir(:,2)=(/0._SP,1._SP,0._SP/)
     field_dir(:,3)=(/0._SP,0._SP,1._SP/)
     !
     do i_dir=1,3
        !
        field_dir_rot = matmul( dl_sop(:,:,is), field_dir(:,i_dir) )
        !
        if (.not.t_rev) rotated_iR(i_dir) = dot_product( field_dir_rot, iR )
        if (     t_rev) rotated_iR(i_dir) = dot_product( iR, field_dir_rot ) 
        !
     enddo
     !
     ! OPTION 2
     !rotated_iR = matmul( dl_sop(:,:,is) , DIP_iR(:,i_c,i_v,ik_mem,i_spin) )
     !if (t_rev) rotated_iR = -conjg( rotated_iR )
     !
     !
   end subroutine GET_the_rotated_iR
   !
end subroutine excitons_lifetimes
