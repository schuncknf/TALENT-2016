program HFSolver
  use :: types
  use :: variables
  use :: HartreeFock
  use :: Minnesota
  implicit none 
  integer :: i,j,k
  real(dp) :: EHF,tr_rho
  Nparticles = 2
  n_orbitals = 64
!  Nsize = n_orbitals/2!
  call read_orbitals
!  write(*,*) nsize
!  write(*,*) n_ho(1:2*nsize)
!  write(*,*) ho_flag
!  write(*,*) l_ho(1:nsize)
!  write(*,*) j_ho(1:nsize)
  !write(*,*) ho_flag
!  stop
! write(*,*) ho_index(1:nsize)
  call read_TBME
!  stop
! write(*,*) v_mat
!  stop
!  do Nsize = 4,4,2
     call initialize_HF
     call Initialize_Minnseota
!     stop
!     v_mat = 0
!     write(*,*) t_mat 
!     stop
     do i = 1,20
        call Construct_rho
        call Construct_gamma
        h_mat = t_mat + gamma_mat
!        write(*,*) h_mat
!        stop
        call Diagonalize_h 
        write(*,*) (E_values(j),j=1,4)
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
     write(*,*) (E_values(i),i=1,12)
     EHF = (Trace_product(t_mat,rho_mat) + Trace_product(h_mat,rho_mat))
     write(*,*) EHF
!  enddo
  
end program HFSolver
