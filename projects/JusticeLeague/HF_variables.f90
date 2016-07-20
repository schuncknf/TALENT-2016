module variables
  use :: types
  implicit none
  integer :: Nsize, Nparticles
  real(dp), allocatable, dimension(:,:) :: D_mat,rho_mat,h_mat,t_mat,Gamma_mat
  real(dp), allocatable, dimension(:,:,:,:) :: v_mat
  real(dp), allocatable, dimension(:) :: E_values, E_prev
  real(dp), parameter :: small = 1.e-5_dp
  real(dp) :: Delta_E
  integer, allocatable, dimension(:) :: n_ho, l_ho

end module variables
