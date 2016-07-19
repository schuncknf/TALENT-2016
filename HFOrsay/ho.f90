module ho
contains
function ho_rad_wf(n,l,r) result(wf)
use maths
use lag
use constants
use lag_pol
implicit none
integer::n,l
double precision::r,nnorm,norm,xi,wf
xi = r/bosc
norm=dsqrt(2**(n+l+2)*fac(n)/(dsqrt(pi)*ffac(2*n+2*l+1)))
nnorm=norm/bosc**(3.d0/2.d0)
wf = nnorm*xi**l*exp(-xi**2/2.d0)*laguerre(n,dble(l+1.d0/2.d0),xi**2)

end function ho_rad_wf
end module ho
