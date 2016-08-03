subroutine compute_rho(rho,Dt,di,nl,nj,iteration)
use constants
use basis
implicit none
integer, intent(in):: di
double precision, dimension (di,di), intent(out) :: rho
double precision, dimension (di,di), intent(in) :: Dt
double precision,dimension(di,di)::D
double precision::rho_temp
integer::nl,nj,iteration
! subroutine local variables
integer :: i,j,k,ibl
rho= 0.d0
D = 0.d0
if (iteration .eq. 1) then
 do i=1,di
  if (nocc(tag_hf(i-1,nl,nj)) .gt. 0) D(i,i) = 1.d0
 enddo
else
D=Dt
endif
! compute rho from D and D_star
do i=1,di
  do j=1,di
rho_temp =0.d0
    do k=1,di
      ibl=tag_hf(k-1,nl,nj)
      rho_temp = rho_temp+ nocc(ibl)*D(i,k)*D(j,k)
    enddo
    rho(i,j) = rho_temp
  enddo
enddo
end subroutine compute_rho


subroutine compute_gamma(gamma_matrix,rho,di,nl,nj)
use basis
implicit none
integer, intent(in) :: di,nl,nj
double precision, dimension (lmin:lmax,jmin:jmax,di,di), intent(in) :: rho
double precision, dimension (di,di), intent(out) :: gamma_matrix
double precision::gammatemp
integer :: i1,i2,i3,i4
integer :: ib1,ib2,ib3,ib4
integer::n3,n4,l3,l4,j3,j4,nocc3,nocc4
integer::n1,n2,l1,l2,j1,j2,nocc2,nocc1

!subroutine local variables
gamma_matrix = 0.d0

!Compute the gamma matrix out of the TBMEs and the rho matrix
do i1=1,di 
  ib1 = tag_hf(i1-1,nl,nj) ! ib1 = (n1,l1,j1)
 do i3=1,di
   ib3 = tag_hf(i3-1,nl,nj) ! ib3 = (n3,l3,j3)
     gammatemp = 0.d0
  do i2=1,red_size
   l2=l_red(i2);j2=j_red(i2);n2=n_red(i2)
   do i4=1,red_size
   l4=l_red(i4);j4=j_red(i4);n4=n_red(i4)
      if (l2 .eq. l4 .and. j2 .eq. j4) then
        gammatemp = gammatemp + 1.d0*tbme_ext(ib1,i2,ib3,i4)*rho(l2,j2,n2+1,n4+1)
        endif
      enddo
    enddo
    gamma_matrix(i1,i3) = gammatemp
  enddo
enddo
end subroutine compute_gamma

