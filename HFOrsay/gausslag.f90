subroutine gausslag(n,n1,n2,inte) 
use ho
use lag_pol
use constants
use lag
implicit none
double precision::inte,wi,xxi,xi,lag_norm,inte2,har_norm
double precision::wj,xxj,xj,phi1,phi2,phi3,phi4,poten
double precision::coeffi,coeffj
integer::i,n1,n2,n3,n4,n,l,j
inte=0.d0
do i=1,n
wi=lag_w(i)
xxi=(lag_zeros(i))
xi = dsqrt(xxi)
coeffi = exp(xxi)*xxi**(0.d0)
!inte = inte + wi*prod_2rad_wf(n1,n2,0,0,xi,xi)*coeffi
!inte = inte + wi*laguerre(n1,0.5d0,xxi)*laguerre(n2,0.5d0,xxi)
inte = inte + wi*ho_rad_wf(n1,0,xi)*ho_rad_wf(n2,0,xi)*coeffi
!prod_2rad_wf(n1,n2,0,0,xi,xi)*coeffi


!inte = inte + wi*coeffi*coeffj*phi1*phi2*poten*phi3*phi4
enddo 
write(*,*) "Integral result",inte
end 
