program HFSolver
  use :: types
  use :: variables
  use :: HartreeFock
  use :: Minnesota
  implicit none 
  integer :: i
  Nsize = 10
  call initialize_HF
  call Initialize_Minnseota
  do i = 1,10
     call Construct_rho
     call Construct_gamma
     h_mat = t_mat + gamma_mat
     call Diagonalize_h 
     delta_E = sum(abs(E_values - E_prev))/real(Nsize,kind=dp)
     E_prev = E_values
     write(*,'(10f10.4)') E_values
     if(delta_E.lt.small) exit 
  enddo

  
end program HFSolver
