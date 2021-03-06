module variables
  use :: types
  implicit none
  integer :: Nsize, Nparticles,n_orbitals,N_n_orbitals,Noccupied,Maxit
  integer, parameter :: N_quad = 95
  real(dp), dimension(1:N_quad) :: w_quad, x_quad, rho_quad, tau_quad, delrho_quad, C_rhorho_quad, C_rhotau_quad, C_rhodelrho_quad
  real(dp), dimension(1:N_quad) :: dC_rhorho_quad, dC_rhotau_quad, dC_rhodelrho_quad
  real(dp), allocatable, dimension(:,:) :: D_mat,rho_mat,h_mat,t_mat,Gamma_mat!, D_prev
  real(dp), allocatable, dimension(:,:,:,:) :: v_mat
  real(dp), allocatable, dimension(:) :: E_values, E_prev
  real(dp) ::  small = 1.e-4_dp
  real(dp) ::  k_Fermi = 3.01_dp
  real(dp) :: Delta_E, C_hartree, C_rhorho, C_rhotau, C_rhodelrho, Trrho
  integer, allocatable, dimension(:) :: n_ho, l_ho,HO_index,m_ho,j_ho,ho_flag,HO_inverse,tz_ho, n_hf,l_hf,j_hf
  logical :: calc_Ivc, calc_couplings, truncated, approximated_rho
  character(30) :: orbitals_file, elements_file, type_of_calculation, output_file, density_file


end module variables
