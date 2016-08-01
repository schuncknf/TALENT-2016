program HFSolver
  use :: types
  use :: variables
  use :: HartreeFock
  use :: Minnesota
  implicit none 
!  integer :: i,j,k
  integer :: i,j,k!,Noccupied
!  real(dp) :: EHF,tr_rho
!  integer :: i,j,k
  real(dp) :: EHF,tr_rho,r,dr
  Nparticles = 8 !< Set the number of neutrons in the system.
  n_orbitals = 216 !< Set the number of available orbitals in your basis space.
  call read_orbitals
  Noccupied=fermi_level()
  write(*,*) Noccupied
!  write(*,*) nsize
!  write(*,*) n_ho(1:2*nsize)
!  write(*,*) ho_flag
!  write(*,*) l_ho(1:nsize)
!  write(*,*) j_ho(1:nsize)

! write(*,*) ho_index(1:nsize)
  call initialize_HF
  call Initialize_Minnseota
!  call read_TBME
  do i = 1,100 !< Set the number of Hartree-Fock iterations.
!     call Construct_rho
!     call Construct_gamma
!     write(*,*) i, gamma_mat(1,1)
!     call  plot_rho_LDA(i)
!     write(*,*) Trace_rho_LDA()
     call calculate_gamma_LDA 
     h_mat = t_mat + gamma_mat
     D_prev = D_mat
     call Diagonalize_h 
     write(*,'(10f11.4)') (E_values(j),j=1,10)
     delta_E = sum(abs(E_values - E_prev))/real(Nsize,kind=dp)
     if(delta_E.lt.small) then
        write(*,*) '  Eigenvalues converged'
        exit 
     endif
     E_prev = E_values
  enddo
  call Construct_rho
!  call Construct_gamma
  call calculate_gamma_LDA 
  h_mat = t_mat + gamma_mat
  EHF = (Trace_product(t_mat,rho_mat)+Trace_product(h_mat,rho_mat))/2._dp !< $E_{HF}=Tr(t\rho)+Tr(h\rho)$
  write(*,*) '  Hartree-Fock energy in MeV'
  write(*,*) EHF
  
end program HFSolver
