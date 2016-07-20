module ho
contains
function ho_rad_wf(n,l,r) result(wf)
use maths
use lag
use constants
use lag_pol
implicit none
integer::n,l
double precision::r,nnorm,norm,xi,wf,expf
xi = r/bosc
norm=dsqrt(2**(n+l+2)*fac(n)/(dsqrt(pi)*ffac(2*n+2*l+1)))
nnorm=norm/bosc**(3.d0/2.d0)
expf = exp(-(xi**2)*0.5d0)
wf = nnorm*xi**l*expf*laguerre(n,dble(l+1.d0/2.d0),xi**2)
!wf = nnorm*xi**l*laguerre(n,dble(l+1.d0/2.d0),xi**2)
end function ho_rad_wf

function prod_2rad_wf(n1,n2,l1,l2,r1,r2) result(wf2)
implicit none
integer::n1,n2,l1,l2
double precision::r1,r2,nnorm,norm,wf2
wf2 = ho_rad_wf(n1,l1,r1)*ho_rad_wf(n2,l2,r2)
end function prod_2rad_wf
end module ho
