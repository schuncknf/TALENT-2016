subroutine compute_gamma(gamma_matrix,TBME,rho,dim)

implicit none
!use constants


integer, intent(in) :: dim
double precision, dimension (dim,dim,dim,dim), intent(in) :: TBME
double precision, dimension (dim,dim), intent(in) :: rho
double precision, dimension (dim,dim), intent(out) :: gamma_matrix

!subroutine local variables
integer :: n1,n2,n3,n4
gamma_matrix(:,:) = 0.d0

!:qcompute the gamma matrix out of the TBMEs and the rho matrix
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
