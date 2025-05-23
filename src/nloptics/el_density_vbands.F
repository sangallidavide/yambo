!
! License-Identifier: GPL
!
! Copyright (C) 2016 The Yambo Team
!
! Authors (see AUTHORS file for details): MG CA   
!
subroutine el_density_vbands(en,Xk,rho,VAL_BANDS)
 !
 !Electronic density from the time-dependent valence bands 
 !(much faster than the standard el_density_vbands() subsoutine)
 !
 !Note that here the identity
 !
 ! \sum_{S_k} F(S_k r) = 1/R_k \sum_{S} F(S r) (1)
 ! 
 !where S_k is a symm. op. in the star of k, S is a generic symm. op.
 !and R_k is the rank of the small group at k.
 !
 !Using (1) the density is calculated in two steps
 !
 ! rho(r) = 1/nkibz \sum_{ n k S_k } f_{n k} |wf_{n k}(S_k^-1 r)|^2=
 !        = 1/nkibz \sum_{S} \sum_{n k} f_{n k}/R_k |wf_{n k}(S^-1 r)|^2 =
 !        = 1/nsym \sum_{S} ( \sum_{n k} f_{n k} w_k |wf_{n k}(S^-1 r)|^2 )
 !
 ! where we have used the relation
 !
 ! w_k = nsym/(nkibz * rank(k))
 !
 !
 use pars,          ONLY:SP,rZERO
 use electrons,     ONLY:levels,n_sp_pol,spin_occ
 use R_lattice,     ONLY:bz_samp
 use D_lattice,     ONLY:nsym,i_time_rev
 use FFT_m,         ONLY:fft_size,fft_rot_r
 use parallel_m,    ONLY:PAR_IND_Xk_ibz,PAR_COM_Xk_ibz_INDEX,PAR_Xk_ibz_index,PAR_Xk_nibz
 use parallel_int,  ONLY:PP_redux_wait
 !
 implicit none
 type(bz_samp),intent(in)  :: Xk
 type(levels), intent(in)  :: en       
 real(SP),     intent(out) :: rho(fft_size)
 complex(SP),  intent(in)  :: VAL_BANDS(fft_size,maxval(en%nbf),PAR_Xk_nibz,n_sp_pol)
 !
 ! Work Space 
 !
 integer     :: i1,ik,is,rho_syms,ir,ik_mem
 real(SP)    :: rho_no_sym(fft_size)
 !
 rho_no_sym=rZERO
 !
 do ik=1,Xk%nibz
   !  
   if (.not.PAR_IND_Xk_ibz%element_1D(ik)) cycle
   ik_mem=PAR_Xk_ibz_index(ik)
   !
   !$omp parallel do default(shared), private(is,ir),reduction(+:rho_no_sym),collapse(2)
   do is=1,n_sp_pol
     do ir=1,fft_size
       rho_no_sym(ir)=rho_no_sym(ir)+real(spin_occ,SP)*Xk%weights(ik)*&
&                         sum(abs(VAL_BANDS(ir,:en%nbf(is),ik_mem,is))**2._SP)
     enddo
   enddo
   !$omp end parallel do
   !
 enddo
 !
 call PP_redux_wait(rho_no_sym,COMM=PAR_COM_Xk_ibz_INDEX%COMM)
 !
 ! Simmetrization
 !
 rho_syms=nsym/(i_time_rev+1)
 rho(:)=rZERO
 ! 
!$omp parallel do default(shared), private(ir), reduction(+:rho)
 do ir=1,fft_size
   rho(ir)=rho(ir)+sum(rho_no_sym(fft_rot_r(ir,1:rho_syms)))/real(nsym,SP)
 enddo
!$omp end parallel do
 !
!$omp parallel do
 do ir=1,fft_size
   rho(ir)=real(1._SP+i_time_rev,SP)*rho(ir)
 enddo
!$omp end parallel do
 !
end subroutine

