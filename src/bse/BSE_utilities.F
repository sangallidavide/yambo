!
! License-Identifier: GPL
!
! Copyright (C) 2023 The Yambo Team
!
! Authors (see AUTHORS file for details): AM MG DS
!
subroutine build_inverse_BS_eh_table(Xk)                                
 !                                                                      
 use electrons,       ONLY:n_sp_pol                                     
 use R_lattice,       ONLY:bz_samp                                      
 use BS,              ONLY:BS_bands,BS_K_dim,BS_blk                     
 use BS_solvers,      ONLY:BSS_eh_table_m1                              
 !                                                                      
 implicit none                                                          
 !                                                                      
 type(bz_samp), intent(in) :: Xk                                        
 !                                                                      
 integer :: i_H,ik_bz,iv,ic,i_sp_c,i_sp_v                               
 !                                                                      
 allocate(BSS_eh_table_m1(Xk%nbz,BS_bands(1):BS_bands(2),BS_bands(1):BS_bands(2),n_sp_pol,n_sp_pol))
 BSS_eh_table_m1=0                                                      
 !call live_timing('BSE inverse eh-table',BS_H_dim)                     
 !                                                                      
 do i_H=1,2*BS_K_dim(1)                                                 
   !                                                                    
   ik_bz    = BS_blk(i_H)%table(1,1)                                    
   iv       = BS_blk(i_H)%table(2,1)                                    
   ic       = BS_blk(i_H)%table(3,1)                                    
   i_sp_c   = BS_blk(i_H)%table(4,1)                                    
   i_sp_v   = BS_blk(i_H)%table(5,1)                                    
   !                                                                    
   BSS_eh_table_m1(ik_bz,iv,ic,i_sp_c,i_sp_v)=i_H                       
   !                                                                    
   !call live_timing(steps=1)                                           
   !                                                                    
 enddo                                                                  
 !call live_timing( )                                                   
 !                                                                      
end subroutine build_inverse_BS_eh_table                                

