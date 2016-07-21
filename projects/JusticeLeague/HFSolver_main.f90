program HFSolver
  use :: types
  use :: variables
  use :: HartreeFock
  use :: Minnesota
  implicit none 
  integer :: i,j,k
  real(dp) :: EHF,tr_rho
  Nparticles = 2
  do Nsize = 6,6,2
     call initialize_HF
     call Initialize_Minnseota
     do i = 1,20
        call Construct_rho
        call Construct_gamma
        h_mat = t_mat + gamma_mat
        call Diagonalize_h 
        write(*,'(4f15.8)') (E_values(j),j=1,4)
        delta_E = sum(abs(E_values - E_prev))/real(Nsize,kind=dp)
        if(delta_E.lt.small) then
           write(*,*) 'Eigenvalues converged'
           exit 
        endif
        E_prev = E_values
     enddo
     call Construct_rho
     call Construct_gamma
     h_mat = t_mat + gamma_mat
     write(*,'(4f15.8)') (E_values(i),i=1,4)
     EHF = Trace_product(t_mat,rho_mat) + Trace_product(h_mat,rho_mat)
     write(*,*) EHF
  enddo
  
end program HFSolver
