!> Solves the problem of \f$n\f$ neutrons in a harmonic oscillator trap
!! by a self-consistent, Hartree-Fock calculation.
!! $E_{HF}=Tr(t\rho)+Tr(h\rho)$
!! Set the number of neutrons in the trap by changing \f$Nparticles\f$,
!! and change the number of orbitals in your basis by changing
!! \f$n_orbitals\f$.
program HFSolver
  use :: types
  use :: variables
  use :: HartreeFock
  use :: Minnesota
  use :: LDA
  implicit none 
  integer :: i,j,k
  real(dp) :: EHF,tr_rho,r,dr
  Nparticles = 8 
  n_orbitals = 216 
  call read_orbitals
  Noccupied=fermi_level()
  call initialize_HF
  call Initialize_Minnseota
!  call read_TBME
  k_Fermi = 3.7_dp
!  do 
!     if (k_Fermi.gt.7.5_dp) exit
!     calc_ivc = .true.
     calc_couplings = .true.
     do i = 1,100
        ! call Construct_rho
        ! write(*,*) Trace(rho_mat)
        ! call Construct_gamma
        ! call  plot_rho_LDA(i)
        ! write(*,*) Trace_rho_LDA()
        ! call sample_rho_LDA
        call sample_DME_fields
        call calculate_gamma_LDA 
        h_mat = t_mat + gamma_mat
        D_prev = D_mat
        call Diagonalize_h 
!        write(*,'(10f11.4)') (E_values(j),j=1,10)
        delta_E = sum(abs(E_values - E_prev))/real(Nsize,kind=dp)
        if(delta_E.lt.small) then
!           write(*,*) '  Eigenvalues converged'
           exit 
        endif
        E_prev = E_values
     enddo
     call Construct_rho
     ! call Construct_gamma
     ! write(*,*)  Trace_rho_LDA()
     ! call plot_rho_lda
     ! call sample_rho_LDA
     call sample_DME_fields
!     call plot_DME_fields
     call calculate_gamma_LDA 
     h_mat = t_mat + gamma_mat
     EHF = (Trace_product(t_mat,rho_mat)+Trace_product(h_mat,rho_mat))/2._dp
!     write(*,*) '  Hartree-Fock energy in MeV'
!     write(*,*) "C_rhorho, C_rhotau = ", C_rhorho, C_rhotau
     write(*,*) k_Fermi, EHF,rms_DME(),i
!     k_Fermi = K_Fermi + 0.01_dp
!  enddo
  
end program HFSolver
