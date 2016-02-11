!
!        Copyright (C) 2000-2016 the YAMBO team
!              http://www.yambo-code.org
!
! Authors (see AUTHORS file for details): DS
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
subroutine WF_derivative(derivative_bands,derivative_kpt,do_derivative)
 !
 !  Remember wf_t = -i \der_t wf   with t=x,y,z
 !  is the periodic part of the derivative
 !
 use pars,           ONLY:SP,DP,pi,cZERO
 use interfaces,     ONLY:WF_load
 use wave_func,      ONLY:WF,wf_ng,wf_x,wf_y,wf_z,WF_derivative_alloc
 use electrons,      ONLY:n_sp_pol
 use FFT_m,          ONLY:fft_size,fft_dim,fft_g_table,fft_norm,fftw_plan
 use R_lattice,      ONLY:g_vec,nkibz
 use stderr,         ONLY:intc
 use LIVE_t,         ONLY:live_timing
 use memory_m,       ONLY:mem_est
 use R_lattice,      ONLY:k_pt
 use D_lattice,      ONLY:alat
#if defined _MAGNETIC
 use SC,             ONLY:SC_bands,found_SC_DB
#endif
 use com,            ONLY:warning
 !
 implicit none
 ! Input variables
 integer :: derivative_bands(2)
 integer :: derivative_kpt(2)
 logical :: do_derivative(3)
 !
 ! Internal variables
 integer :: bands_to_load(2)
 logical :: use_live_timing
 complex(SP), allocatable  :: wf_der(:)
 complex(DP), allocatable  :: wf_DP(:)
 ! Dummys
 integer :: i1,i_spin,ib,ik,ifft
 !
 ! wf_load could change bands_to_load(:) if it founds an SC-database 
 bands_to_load(:)=derivative_bands(:)
 !
 ! Wavefunctions in G space 
 !==========================
 call WF_load(WF,0,1,bands_to_load,derivative_kpt,space='G',title='-DERIV',impose_free_and_alloc=.TRUE.)
 !
 ! Need to call fft_setup to know the size of the FFT mesh
 ! 
 call fft_setup(0,1,.false.)
 !
#if defined _MAGNETIC
 if (found_SC_DB.and.any ( bands_to_load/=(/1,SC_bands(2)/) )) then
   bands_to_load=(/1,SC_bands(2)/)
   call warning(' [WF-DERIV] Forced '//trim(intc(SC_bands(2)))//' bands to be derived')
 endif
#endif
 !
 call WF_derivative_alloc(do_derivative)
 !
 allocate(wf_der(wf_ng))                              ! G space
 call mem_est("wf_der",(/wf_ng/))
 allocate(wf_DP(fft_size))                            ! real space
 call mem_est("wf_DP",(/fft_size/),(/DP/))
 !
 use_live_timing=.false.
 if (derivative_kpt(1)==1 .and. derivative_kpt(2)==nkibz) use_live_timing=.true.
 if (use_live_timing) call live_timing('[WF-DERIV] Evaluation loop',WF%N)
#if defined _FFTW
 fftw_plan = 0
#endif
 !
 do ik=derivative_kpt(1),derivative_kpt(2)
   do ib=bands_to_load(1),bands_to_load(2)
     do i_spin=1,n_sp_pol   
       ifft=WF%index(ib,ik,i_spin)
       !
       ! Evaluate the wf_derivative in G space
       ! wf_derivative(:,i)=-i grad_i wf(:)=p_i wf(:)    i=x,y
       ! wf_derivative(G,i)= (k_i+G_i) wf(:)
       ! wf_z is not needed
       !
       do i1=1,3
         if (.not.do_derivative(i1)) cycle
         wf_der=cZERO
         wf_der(:)=WF%c(:,ifft)*(k_pt(ik,i1)+g_vec(: wf_ng ,i1))*2.0_SP*pi/alat(i1)
         !
         ! Here we FFT wf_derivative to real space
         wf_DP=(0._DP,0._DP)
         wf_DP(fft_g_table(1:wf_ng,1))= wf_der(1:wf_ng)*fft_norm
#if defined _FFTW
         call fft_3d(wf_DP,fft_dim,1,fftw_plan)
#else
         !call fft_3d(wf_DP,fft_dim,1)
         !
         ! here we use +2 in order to exploit the wfc cutoff 
         ! when using FFTQE
         call fft_3d(wf_DP,fft_dim,2)
#endif
         if(i1==1) wf_x(:,ifft)=cmplx(wf_DP(:),kind=SP)
         if(i1==2) wf_y(:,ifft)=cmplx(wf_DP(:),kind=SP)
         if(i1==3) wf_z(:,ifft)=cmplx(wf_DP(:),kind=SP)
         !
       enddo
       if(use_live_timing) call live_timing(steps=1)
       !
     enddo
   enddo
 enddo
 !
 if(use_live_timing) call live_timing()
 call WF_free(WF)
 !
 deallocate(wf_DP)
 call mem_est("wf_DP")
 deallocate(wf_der)
 call mem_est("wf_der")
 !
end subroutine