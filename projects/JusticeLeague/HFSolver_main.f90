program HFSolver
  use :: types
  use :: variables
  use :: HartreeFock
  use :: Minnesota
  implicit none 
  integer :: i,j,k
  real(dp) :: EHF
  Nparticles = 2
  do Nsize = 6,6,2
     call initialize_HF
     call Initialize_Minnseota
     do i = 1,10
        call Construct_rho
        call Construct_gamma
        h_mat = t_mat + gamma_mat
        call Diagonalize_h 
        delta_E = sum(abs(E_values - E_prev))/real(Nsize,kind=dp)
        E_prev = E_values
        if(delta_E.lt.small) exit 
     enddo
     write(*,'(4f15.8)') (E_values(i),i=1,4)
     EHF = 0
     do i = 1,1!Nsize
        do j = 1,1!Nsize
           EHF = EHF + 2*t_mat(i,j)*rho_mat(j,i)+(Gamma_mat(i,j)*rho_mat(j,i))*0.5_dp*4
        enddo
     enddo

     write(*,*) EHF
  enddo
  
end program HFSolver
