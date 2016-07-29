subroutine compute_rho(rho,Dt,di,nl,nj)
use constants
use basis
implicit none
integer, intent(in):: di
double precision, dimension (di,di), intent(out) :: rho
double precision, dimension (di,di), intent(in) :: Dt
double precision,dimension(di,di)::D
double precision::rho_temp
integer::nl,nj
! subroutine local variables
integer :: i,j,k,ibl
rho= 0.d0
D=Dt
! compute rho from D and D_star
do i=1,di
  do j=1,di
rho_temp =0.d0
    do k=1,di+1
      ibl=tag_hf(k-1,nl,nj)
      rho_temp = rho_temp+ nocc(ibl)*D(i,k)*D(j,k)
!      if (nocc(ibl) .ne. 0) write(*,*) "In rho, nocc",nocc(ibl)
      !rho_temp = rho_temp+ occ_num(k)*D(i,k)*D(j,k)
    enddo
    rho(i,j) = rho_temp
  enddo
enddo
end subroutine compute_rho


subroutine compute_gamma(gamma_matrix,mtrxel,rho,di,nl,nj)
use basis
implicit none
integer, intent(in) :: di,nl,nj
double precision, dimension (di,di,di,di), intent(in) :: mtrxel
double precision, dimension (di,di), intent(in) :: rho
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
  ib1 = tag_hf(i1-1,nl,nj)
 do i2=1,di
   ib2 = tag_hf(i2-1,nl,nj)
     gammatemp = 0.d0
 !    if (l1 .eq. l2 .and. j1 .eq. j2) then
    ! if (nocc2 .ne. 0 .and. nocc1 .ne. 0) then
  do i3=1,di
   ib3 = tag_hf(i3-1,nl,nj)
   do i4=1,di
   ib4 = tag_hf(i4-1,nl,nj)
     !if (nocc4 .ne. 0 .and. nocc3 .ne. 0) then
!       if (j3 .eq. j4 .and. l3 .eq. l4) then
     !  if (n3 .ne. 0) then
 ! write(6,'(6i3)') n3,n4,l3,l4,j3,j4
  !endif
        gammatemp = gammatemp + mtrxel(ib1,ib4,ib2,ib3)*rho(i3,i4)
!        endif ! J3=J4;L3=L4
    !endif !Occupied states only
      enddo
    enddo
    gamma_matrix(i1,i2) = gammatemp
!    endif
   ! endif
  enddo
enddo
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


