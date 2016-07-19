subroutine compute_rho(D,rho,dim)

implicit none
! use constants

integer :: dim
complex, dimension (dim,dim), intent(in) :: D
double precision, dimension (dim,dim), intent(out) :: rho


! subroutine local variables
integer :: i,j,k
rho(:,:) = 0.d0

! compute rho from D and D_star
do i=1,dim
  do j=1,dim
    do k=1,dim
      rho(i,j) = rho(i,j) + D(i,k)*conjg(D(j,k))
    enddo
  enddo
enddo

end subroutine compute_rho
