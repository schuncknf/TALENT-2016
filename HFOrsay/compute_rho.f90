subroutine compute_rho(rho,D,dim)

implicit none
! use constants

integer, intent(in):: dim
double precision, dimension (dim,dim), intent(out) :: rho
double precision, dimension (dim,dim), intent(in) :: D


! subroutine local variables
integer :: i,j,k
rho= 0.d0

! compute rho from D and D_star
do i=1,dim
  do j=1,dim
    do k=1,dim
      rho(i,j) = rho(i,j) + D(i,k)*(D(j,k))
    enddo
  enddo
enddo

end subroutine compute_rho
