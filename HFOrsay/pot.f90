!> This module contains the required potentials
module pot
  contains
!> This function compute in r-space the value of the Minessota potential
!> \param r1, point 1 in r-space
!> \param r2, point 2 in r-space
function minessota(r1,r2) result(v)
use constants
implicit none
double precision::v1,v2,r1,r2,v,ri,rj
ri=dsqrt(r1)*bosc
rj=dsqrt(r2)*bosc
v1 =v0r/kr*(exp(-kr*(ri**2+rj**2-2*ri*rj))-exp(-kr*(ri**2+rj**2+2*ri*rj)))/(4.d0*ri*rj) 
v2=-v0s/ks*(exp(-ks*(ri**2+rj**2-2*ri*rj))-exp(-ks*(ri**2+rj**2+2*ri*rj)))/(4.d0*ri*rj)
v=(v1+v2)*half
end function
!> This subroutine construct by block the Kinetic Matrix for an harmonic trap
!>\param n, principal quantum number
!>\param l, orbital quantum number
!>\param t_mat, the block diagonal kinetic matrix
subroutine t_bloc(n,l,t_mat)
 use constants
 implicit none
 double precision::t_mat(n,n)
 integer::l,n,i
 do i=1,n
   t_mat(i,i) = (2.d0*(i-1)+l+1.5d0)*ama*2.d0/(bosc**2)
 enddo
 end subroutine

!
end module
