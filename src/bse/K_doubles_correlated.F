!
! Copyright (C) 2000-2008 D. Sangalli and the YAMBO team 
!              http://www.yambo-code.org
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
   !
   ! DS: 2023/08/07
   ! This should be redone and coded in the same spirit of the exc-ph kernel
   ! where the static BSE is already available and loaded from disk
!
subroutine K_doubles_correlated(Ken)
   !
   ! DUOBLES CORRELATION PROJECT
   ! *******************************
   ! We do a static diagonalization to get the RPA eigenvalues and
   ! eigenvectors. Then construct from these the RPA doubles eigenvectors.
   ! (Write a new database for the doubles with CU instead of C (this for restart) 
   ! and BS_doubles_E(i_doubles)=E_1^RPA+E_2^RPA
 !
 use interfaces,     ONLY:LINEAR_algebra_driver
 use BS_solvers,     ONLY:BSS_eh_table,BS_mat
 use BS_doubles,     ONLY:BS_mat_tmp,Eigen_doubles,doubles_corr
 use linear_algebra, ONLY:DIAGO,USE_LK
 !
 implicit none
 !
 integer    ,allocatable :: BSS_eh_table_m1(:,:)
 real(SP)   ,allocatable :: BS_E_tmp(:)
 !
      allocate(BSS_eh_table_m1(BS_bands_doubles(2),BS_bands_doubles(2)))
      BSS_eh_table_m1=0
      allocate(BS_E_tmp(BS_K_dim(1)))
      allocate(Eigen_doubles(BS_K_dim(1)**2,1))
   !
   if (doubles_corr) then
     do iv1=BS_bands_doubles(1),Ken%nbf(1)
       do ic1=Ken%nbf(1)+1,BS_bands_doubles(2)
         icv1=icv1+1
         BSS_eh_table_m1(iv1,ic1)=icv1
       enddo
     enddo
     BS_mat_doubles(:,:,2)=(0.,0.)
     BS_doubles_E(:)=0.
     do icv1=2,BS_K_dim(1)
       forall(icv2=1:icv1-1) BS_mat_tmp(icv1,icv2,1)=conjg(BS_mat_tmp(icv2,icv1,1))
     enddo
     do icv1=1,BS_K_dim(1)
       iv1=BSS_eh_table(icv1,2)
       ic1=BSS_eh_table(icv1,3)
       BS_mat_tmp(icv1,icv1,1)=BS_mat_tmp(icv1,icv1,1)+Ken%E(ic1,1,1)-Ken%E(iv1,1,1)
     enddo
     call LINEAR_algebra_driver(DIAGO,M=BS_mat_tmp(:,:,1),E_real=BS_E_tmp)
     !
     ! Paral <
     call PARALLEL_index(px,(/BS_K_doubles_dim/),(/1/))
     call PP_wait()
     call live_timing('Doubles correlation',px%n_of_elements(myid+1))
     ! Paral >
     !
     do i_doubles=1,BS_K_doubles_dim
       ! Paral <
       if (.not.px%element_1D(i_doubles)) cycle
       call live_timing(steps=1)
       ! Paral >
       iv1=BS_doubles_table(i_doubles,1)
       iv2=BS_doubles_table(i_doubles,2)
       ic1=BS_doubles_table(i_doubles,3)
       ic2=BS_doubles_table(i_doubles,4)
       icv1=BSS_eh_table_m1(iv1,ic1)
       icv2=BSS_eh_table_m1(iv2,ic2)
       BS_doubles_E(i_doubles)=BS_E_tmp(icv1)+BS_E_tmp(icv2)
       ! Important: here the indexes loops in the same way in i_doubles and icv1,icv2 in the
       ! tables, i_doubles=icv1,ivc2=iv1,ic1,iv2,ic2 ---> the last index goes faster
       do icv=1,BS_K_dim(1)
         Eigen_doubles( 1+(icv-1)*BS_K_dim(1) : icv*BS_K_dim(1) ,1)=&
&                  BS_mat_tmp(:,icv1,1)*BS_mat_tmp(icv,icv2,1)
       enddo
       call cgemv('n',BS_K_dim(1),BS_K_doubles_dim,cONE,&
&        BS_mat_doubles(:,:,1),BS_K_dim(1),Eigen_doubles(:,1),1,&
&        cZERO,BS_mat_doubles(:,i_doubles,2),1)
     enddo
     !
     ! Paral <
     call pp_redux_wait(BS_mat_doubles(:,:,2))
     call pp_redux_wait(BS_doubles_E)
     call live_timing()
     call pp_indexes_reset(px)
     ! Paral >
     BS_mat_doubles(:,:,1)=BS_mat_doubles(:,:,2)
     !
#if defined _DEBUG
     do icv=1,BS_K_dim(1)
       write(*,*) BS_E_tmp(icv)*27.2116
     enddo
#endif
     !
     deallocate(BSS_eh_table_m1)
     deallocate(BS_E_tmp)
     deallocate(Eigen_doubles)
     !
   endif
   !
end subroutine K_doubles_correlated
