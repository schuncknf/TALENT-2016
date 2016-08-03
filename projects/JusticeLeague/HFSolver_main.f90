!> \mainpage Installation and Instructions For Use
!! 
!!  \section intro_sec Introduction
!! 
!! Solves the problem of \f$n\f$ neutrons in a harmonic oscillator trap
!! by a self-consistent, Hartree-Fock calculation.
!! \f$E_{HF}=Tr(t\rho)+Tr(h\rho)\f$
!! 
!!  \section install_sec Installation
!! 
!!  The source code for \f$\texttt{HFSolver}\f$ includes the following
!!  Fortran modules:
!!
!!  -\f$\texttt{HFSolver\_main.f90}\f$
!!  -\f$\texttt{HF\_HartreeFock.f90}\f$
!!  -\f$\texttt{HF\_HOQuadrature.f90}\f$
!!  -\f$\texttt{HF\_LDA.f90}\f$
!!  -\f$\texttt{HF\_Minnesota.f90}\f$
!!  -\f$\texttt{HF\_variables.f90}\f$
!!  -\f$\texttt{types.f90}\f$
!!
!!  as well as one input file, \f$\texttt{HF\_input.dat}\f$ (to be
!!  described in the following section).
!!
!!  To compile, simply type \f$\texttt{make}\f$. The
!!  executable is called \f$\texttt{HFSolver}\f$.
!! 
!! \section instructions Instructions For Use
!! 
!! \subsection input_file HF_input.dat
!!
!! A template input file is included. No default arguments are provided
!! by the code, so you MUST specify every item in the output file and
!! in the order they are given in the template; otherwise the
!! code will not run. Furthermore, the input file is case-
!! sensitive. The individual entries are explained below:
!! 
!! \subsubsection Nparticles Number of particles 
!! Set the number of neutrons in the trap
!!
!! \subsubsection max_iter Maximum number of iterations
!! Set the maximum number of iterations
!! 
!! \subsubsection k_fermi Fermi Momentum
!! Specify the value of the Fermi momentum \f$k_F\f$ (in
!! \f$fm^{-1}\f$) to be used in LDA and DME calculations. Typical
!! values are between 1 and 2.
!!
!! \subsubsection orbitals_file Orbitals file
!! List the name of the input file in which you index your
!! m-scheme orbitals. If they were generated using Morten's
!! \f$\texttt{vNN}\f$ code, the file will probably be called
!! \f$\texttt{spM.dat}\f$. The file should have the following format:
!!
!! \f$\texttt{Legend:         n     l     2j   2mj   2tz  }\f$\n
!! \f$\texttt{Orbit number:  1     0   0   1  -1  -1}\f$\n
!! \f$\texttt{Orbit number:  2     0   0   1   1  -1}\f$\n
!! \f$\texttt{Orbit number:  3     0   0   1  -1   1}\f$\n
!! \f$\texttt{Orbit number:  4     0   0   1   1   1}\f$\n
!! \n\f$\texttt{\dots}\f$\n
!!
!! \subsubsection matrix_file Matrix elements file
!! List the name of the input file in which you list your two-body
!! m-scheme matrix elements. If they were generated using Morten's
!! \f$\texttt{vNN}\f$ code, the file will probably be called
!! \f$\texttt{VM-scheme.dat}\f$. The file should have the following
!! format:
!!
!! \f$\texttt{Legend:         n     l     2j   2mj   2tz  }\f$\n
!! \f$\texttt{  Matrix elements < ab | V | cd> in mscheme}\f$\n
!! \f$\texttt{ Legend: a  b  c  d and  < ab | V | cd>}\f$\n
!! \f$\texttt{  5  6  5 13     0.113470E+01}\f$\n
!! \f$\texttt{  5  6  5 22     -.287678E+00}\f$\n
!! \f$\texttt{  5  6  5 29     0.406838E+00}\f$\n
!! \f$\texttt{  5  6  5 38     -.659420E-01}\f$\n
!! \n\f$\texttt{\dots}\f$\n
!!
!! where \f$a,b,c,d\f$ refer to the orbitals indexed in
!! \f$\mathbf{Orbitals\ file}\f$
!!
!! \subsubsection out_file Output file
!! Specify the name of your output file.
!!
!! \subsubsection dens_file
!!
!! Specify the name of the density file, which will contain the density
!! distribution on a format suitable for plotting utilities such as
!! gnuplot.
!!
!! \subsubsection potential_type truncated or spherical
!!
!! If the option \f$\texttt{truncated}\f$ is selected, the program
!! assumes \f$l=0\f$ and calculates the two-body matrix elements
!! exactly using the Minnesota potential (using the parameters given in
!! HF_truncated_v2.pdf).
!!
!! If the option \f$\texttt{spherical}\f$ is selected, the program will
!! read the matrix elements from a file,
!! \f$\mathbf{Matrix\ elements\ file}\f$.
!! There is (in principle) no limitation on \f$l\f$ in this case.
!!
!! If the option \f$\texttt{LDA}\f$ is selected, the program will
!! use the Local Density Approximation to the Hamiltonian.
!!
!! If the option \f$\texttt{spherical}\f$ is selected, the program will
!! use the Density Matrix Expansion to approximate the Hamiltonian.
!!
!! \subsubsection n_max N_max
!! Specify the maximum value of the orbital quantum number \f$n\f$.
!! If using the option \f$\texttt{spherical}\f$, \f$\mathbf{N\_max}\f$
!! can be no
!! larger than the maximum value of \f$n\f$ specified in
!! \f$\mathbf{Orbitals\ file}\f$.
!!
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
