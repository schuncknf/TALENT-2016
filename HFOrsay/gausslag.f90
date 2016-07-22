subroutine gausslag(n,n1,n2,inte) 
use ho
use lag_pol
use constants
use maths
implicit none
double precision::inte,wi,xxi,xi,lag_norm,inte2,har_norm
double precision::wj,xxj,xj,phi1,phi2,phi3,phi4,poten
double precision::coeffi,coeffj,norm_lag,lag10,testwf
integer::i,n1,n2,n3,n4,n,l,j
inte=0.d0
testwf =0.d0
coeffi = 0.d0
!nosc1=ho_norm(n1,0.d0)
!nosc2=ho_norm(n2,0.d0)
do i=1,n
wi=lag_w(i)
testwf = testwf + wi
xxi=(lag_zeros(i))!*bosc**2
norm_lag=Gamma(dble(n1)+1.d0+0.5d0)/fac(n1)
xi = dsqrt(xxi)
coeffi = exp(xi)*xi**(-0.5d0)
lag10=laguerre(n1,0.5d0,xxi)*laguerre(n2,0.5d0,xxi)
inte = inte + wi*lag10
!laguerre(n1,0.5d0,xxi)*laguerre(n2,0.5d0,xxi)
!inte = inte + wi*xxi*ho_rad_wf(n1,0,xi)*ho_rad_wf(n2,0,xi)!*coeffi
enddo 
write(*,*) "n1,n2,res",n1,n2,inte/norm_lag
write(*,*) "Test weigth",testwf
end 
