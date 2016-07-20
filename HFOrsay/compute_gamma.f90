subroutine compute_gamma(gamma_matrix,TBME,rho,dim)

implicit none
!use constants


integer, intent(in) :: dim
double precision, dimension (dim,dim,dim,dim), intent(in) :: TBME
double precision, dimension (dim,dim), intent(in) :: rho
double precision, dimension (dim,dim), intent(out) :: gamma_matrix

!subroutine local variables
integer ::  i,j,k,l
gamma_matrix(:,:) = 0.d0

!compute the gamma matrix out of the TBMEs and the rho matrix
do i=1,dim
  do j=1,dim
    do k=1,dim
      do l=1,dim
        gamma_matrix(i,j) = gamma_matrix(i,j) + TBME(i,j,k,l)*rho(k,l)
      enddo
    enddo
  enddo
enddo

end subroutine compute_gamma
