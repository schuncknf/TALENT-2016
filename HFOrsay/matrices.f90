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
D=Dt
! compute rho from D and D_star
do i=1,di
  do j=1,di
rho_temp =0.d0
    do k=1,npart/2
      rho_temp = rho_temp+ 1.d0*D(i,k)*D(j,k)
    enddo
    rho(i,j) = rho_temp
  enddo
enddo


! --- Petar

!   do i = 1, nt
!      do j = 1, nt
!	 rho_temp=0.d0
!	 do k = 1, nt (sum over all states, when the state is not occupied it will not contribute, occupation numbers calculated in basis)
!	      rho_temp = rho_temp + occ(i)*D(i,k)*D(j,k)
!	 enddo
!         rho(i,j)=rho_temp
!       enddo
!   enddo
! --- Petar


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

!Compute the gamma matrix out of the TBMEs and the rho matrix
do n1=1,di 
 do n2=1,di
gammatemp = 0.d0
  do n3=1,di
   do n4=1,di
        gammatemp = gammatemp + mtrxel(n1,n4,n2,n3)*rho(n3,n4)
      enddo
    enddo
    gamma_matrix(n1,n2) = gammatemp
  !  gamma_matrix(n2,n1) = gammatemp
  enddo
enddo


!---- Petar 

!do i=1,nt
! do k=1,nt
!gammatemp = 0.d0
!  do j=1,nt ! sum over all states, the unoccupied states will not contribute
! read quantum numbers n2,l2,j2 that corresponds to j state
!   do l=1,nt ! sum over all states, the unoccupied states will not contribute
! read quantum numbers n4,l4,j4 that corresponds to l state
!	 if(occj .eq. 0 .or. occl .eq. 0) then	 
!            matel=0
!	 elseif(l2 .ne. l4 .or. j2 .ne. j4) then
!	    matel = 0
!	 else
!	    matel = mtrxel(n1,n4,n2,n3)*rho(n3,n4)
!	 endif
!         gammatemp = gammatemp + mtrxel(n1,n2,n3,n4)*rho(n4,n2)
!      enddo
!    enddo
!    gamma_matrix(i,k) = gammatemp
!  enddo
! enddo

!---- Petar 


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


!--- Petar

!do i=1,nt
!  do j=1,nt
!    h(i,j) = t(i,j) + gamma_matrix(i,j)
!  enddo
!enddo

!--- Petar



end subroutine


