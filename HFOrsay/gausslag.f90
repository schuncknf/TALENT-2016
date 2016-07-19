subroutine gausslag(n,l)
use maths
use lag
use lag_pol
use ho
use constants
implicit none
double precision::inte,wi,xxi,xi,lag_norm,inte2,har_norm
integer::i,n1,n2,n,l
n1=3
n2=3
inte=0.d0
inte2=0.d0
do i=1,n
wi=lag_w(i)
xxi=(lag_zeros(i))
xi=dsqrt(xxi)
lag_norm = gamma(dble(n1+1)+1/2.d0)/fac(n1)
har_norm = half
!inte = inte + wi*laguerre(n1,1/2.d0,xi)*laguerre(n2,1/2.d0,xi)
inte = inte + wi*ho_rad_wf(n1,l,xi)*ho_rad_wf(n2,l,xi)*exp(xxi)*xxi**(l)
enddo 

write(*,*) "Integral value= ",inte*har_norm

end 
