subroutine tbme(n,n1,n2,n3,n4,resu) 
use ho
use lag_pol
use constants
use lag
use pot
implicit none
integer::i,n1,n2,n3,n4,n,j
double precision::inte1,inte2,resu
double precision::osc1,osc2
double precision::coeffi,coeffj
double precision::xxi,xi,xxj,xj
double precision::wi,wj
inte1=0.d0
inte2=0.d0
resu=0.d0
do i=1,n
    wi=lag_w(i)
    xxi=(lag_zeros(i))
    xi = dsqrt(xxi)
    coeffi = exp(xxi)*xxi**(0.d0)
    osc1 = coeffi*ho_rad_wf(n1,0,xi)*ho_rad_wf(n3,0,xi)
 do j=1,n
    wj = lag_w(j)
    xxj=(lag_zeros(j))
    xj = dsqrt(xxj)
    coeffj = exp(xxj)*xxj**(0.d0)
    osc2 = coeffj*ho_rad_wf(n2,0,xj)*ho_rad_wf(n4,0,xj)
    inte1 = inte1 + osc1*osc2*minnesota(xxi,xxj,v0r,kr)*coeffi*coeffj
    inte2 = inte2 + osc1*osc2*minnesota(xxi,xxj,v0s,ks)*coeffi*coeffj
  enddo !j
 enddo !i
resu = -inte1 + inte2
write(*,*) "Integral result",resu
end 
