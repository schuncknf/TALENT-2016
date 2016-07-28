module variables
  use :: types
  implicit none
  integer :: Nsize, Nparticles,n_orbitals,N_n_orbitals,Noccupied
  real(dp), allocatable, dimension(:,:) :: D_mat,rho_mat,h_mat,t_mat,Gamma_mat, D_prev
  real(dp), allocatable, dimension(:,:,:,:) :: v_mat
  real(dp), allocatable, dimension(:) :: E_values, E_prev
  real(dp), parameter :: small = 1.e-6_dp, k_Fermi = 5._dp 
  real(dp) :: Delta_E
  integer, allocatable, dimension(:) :: n_ho, l_ho,HO_index,m_ho,j_ho,ho_flag,HO_inverse,tz_ho, n_hf,l_hf,j_hf


end module variables
