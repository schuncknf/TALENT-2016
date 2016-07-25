module pot
  contains
  function potential(r1,r2,v0,mu) result(v)
  use constants
  implicit none
  double precision::r1,r2,v,mu,v0
  r1=dsqrt(r1)*bosc
  r2=dsqrt(r2)*bosc
  v=(v0*(exp(-mu*(r1**2+r2**2))*sinh(2.d0*mu*r1*r2)/(4.d0*r1*r2*mu)))
  end function

function mine(r1,r2) result(v)
use constants
implicit none
double precision::v1,v2,r1,r2,v,ri,rj
ri=dsqrt(r1)*bosc
rj=dsqrt(r2)*bosc
v1 =v0r/kr*(exp(-kr*(ri**2+rj**2-2*ri*rj))-exp(-kr*(ri**2+rj**2+2*ri*rj)))/(4.d0*ri*rj) 
v2=-v0s/ks*(exp(-ks*(ri**2+rj**2-2*ri*rj))-exp(-ks*(ri**2+rj**2+2*ri*rj)))/(4.d0*ri*rj)
!v1 =v0r/kr*(exp(-kr*(ri-rj)**2)-exp(-kr*(ri+rj)**2))/(4.d0*ri*rj) 
!v2 =v0s/ks*(exp(-ks*(ri-rj)**2)-exp(-ks*(ri+rj)**2))/(4.d0*ri*rj) 
v=(v1+v2)*half
end function

subroutine kinetic(kin,n)
 use constants
 implicit none
 double precision::kin(n,n)
 integer::n,i
 kin=0.d0
 do i=1,n
   kin(i,i) = (2.d0*(i-1)+1.5d0)*ama*2.d0/(bosc**2)
 enddo
 end subroutine
end module
