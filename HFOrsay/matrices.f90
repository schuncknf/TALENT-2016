subroutine compute_rho(rho,Dt,di)
use constants
implicit none
integer, intent(in):: di
double precision, dimension (di,di), intent(out) :: rho
double precision, dimension (di,di), intent(in) :: Dt
double precision,dimension(di,di)::D
double precision::rho_temp
! subroutine local variables
integer :: i,j,k
rho= 0.d0
rho_temp =0.d0
D=Dt
! compute rho from D and D_star
do i=1,di
  do j=1,di
    do k=1,npart/2
      rho(i,j) = rho(i,j)+ 2.d0*D(i,k)*D(j,k)
    enddo
!    rho(i,j) = rho_temp
!    rho(j,i) = rho_temp
  enddo
enddo
end subroutine compute_rho


subroutine compute_gamma(gamma_matrix,mtrxel,rho,di)
implicit none
integer, intent(in) :: di
double precision, dimension (di,di,di,di), intent(in) :: mtrxel
double precision, dimension (di,di), intent(in) :: rho
double precision, dimension (di,di), intent(out) :: gamma_matrix
double precision::gammatemp
integer :: n1,n2,n3,n4

!subroutine local variables
gamma_matrix = 0.d0
gammatemp = 0.d0

!Compute the gamma matrix out of the TBMEs and the rho matrix
do n1=1,di
 do n2=1,di
  do n3=1,di
   do n4=1,di
        gammatemp = gammatemp + mtrxel(n1,n4,n2,n3)*rho(n3,n4)
      enddo
    enddo
    gamma_matrix(n1,n2) = gammatemp
  enddo
enddo
!Petar idea change order and introduce gamma_temp

end subroutine compute_gamma


subroutine compute_h(h,t,gamma_matrix,di)
implicit none
integer, intent(in) :: di
double precision, dimension (di,di), intent(in) :: t
double precision, dimension (di,di), intent(in) :: gamma_matrix
double precision, dimension (di,di), intent(out) :: h
!subroutine local variables
integer :: i,j
h = 0.d0

do i=1,di
  do j=1,di
    h(i,j) = t(i,j) + gamma_matrix(i,j)
  enddo
enddo

end subroutine


