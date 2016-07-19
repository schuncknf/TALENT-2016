subroutine compute_h(h,t,gamma_matrix,dim)

implicit none
!use constants

integer, intent(in) :: dim
double precision, dimension (dim,dim), intent(in) :: t
double precision, dimension (dim,dim), intent(in) :: gamma_matrix
double precision, dimension (dim,dim), intent(out) :: h

!subroutine local variables
integer :: i,j
h(:,:) = 0.d0

do i=1,dim
  do j=1,dim
    h(i,j) = t(i,j) + gamma_matrix(i,j)
  enddo
enddo

end subroutine
