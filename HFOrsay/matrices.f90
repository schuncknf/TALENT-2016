subroutine compute_rho(rho,D,dim)
implicit none
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


subroutine compute_gamma(gamma_matrix,TBME,rho,dim)
implicit none
integer, intent(in) :: dim
double precision, dimension (dim,dim,dim,dim), intent(in) :: TBME
double precision, dimension (dim,dim), intent(in) :: rho
double precision, dimension (dim,dim), intent(out) :: gamma_matrix

!subroutine local variables
integer :: n1,n2,n3,n4
gamma_matrix = 0.d0

!Compute the gamma matrix out of the TBMEs and the rho matrix
do n2=1,dim
 do n4=1,dim
  do n1=1,dim
   do n3=1,dim
        gamma_matrix(n2,n4) = gamma_matrix(n2,n4) + TBME(n1,n2,n3,n4)*rho(n3,n4)
      enddo
    enddo
  enddo
enddo
!Petar idea change order and introduce gamma_temp

end subroutine compute_gamma


subroutine compute_h(h,t,gamma_matrix,dim)
implicit none
integer, intent(in) :: dim
double precision, dimension (dim,dim), intent(in) :: t
double precision, dimension (dim,dim), intent(in) :: gamma_matrix
double precision, dimension (dim,dim), intent(out) :: h
!subroutine local variables
integer :: i,j
h = 0.d0

do i=1,dim
  do j=1,dim
    h(i,j) = t(i,j) + gamma_matrix(i,j)
  enddo
enddo

end subroutine


