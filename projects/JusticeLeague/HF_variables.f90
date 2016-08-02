module variables
  use :: types
  implicit none
  integer :: Nsize, Nparticles,n_orbitals,N_n_orbitals,Noccupied
  integer, parameter :: N_quad = 95
  real(dp), dimension(1:N_quad) :: w_quad, x_quad, rho_quad, tau_quad, delrho_quad
  real(dp), allocatable, dimension(:,:) :: D_mat,rho_mat,h_mat,t_mat,Gamma_mat, D_prev
  real(dp), allocatable, dimension(:,:,:,:) :: v_mat
  real(dp), allocatable, dimension(:) :: E_values, E_prev
  real(dp), parameter :: small = 1.e-4_dp
  real(dp) ::  k_Fermi = 3.01_dp 
  real(dp) :: Delta_E, C_hartree, C_rhorho, C_rhotau, C_rhodelrho
  integer, allocatable, dimension(:) :: n_ho, l_ho,HO_index,m_ho,j_ho,ho_flag,HO_inverse,tz_ho, n_hf,l_hf,j_hf
  logical :: calc_Ivc, calc_couplings


end module variables
