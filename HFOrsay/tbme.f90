subroutine tbme(n,n1,n2,n3,n4,resu)
use ho
use lag_pol
use constants
use lag
use pot
implicit none
integer::i,n1,n2,n3,n4,n,j,l1,l2,l3,l4
double precision::inte1,inte2,resu
double precision::osc,osc1,osc2
double precision::coeffi,coeffj
double precision::xxi,xi,xxj,xj
double precision::wi,wj
l1=0;l2=0;l3=0;l4=0
inte1=0.d0
coeffi =0.d0
coeffj=0.d0
osc2=0.d0
osc=0.d0
inte2=0.d0
resu=0.d0
xxi = 0.d0
xxj = 0.d0
xi = 0.d0
xj = 0.d0
do i=1,n
    osc = 0.d0
    osc1 = 0.d0
    osc2 = 0.d0
    wi=lag_w(i)
    xxi=(lag_zeros(i))
    xi = dsqrt(xxi)
    coeffi =1.d0!exp(xxi)*xxi**(0.d0)
    osc1 = coeffi*ho_rad_wf(n1,l1,xi)*ho_rad_wf(n3,l3,xi)
 do j=1,n
    wj = lag_w(j)
    xxj=(lag_zeros(j))
    xj = dsqrt(xxj)
    coeffj = 1.d0!exp(xxj)*xxj**(0.d0)
    osc2 = coeffj*ho_rad_wf(n2,l2,xj)*ho_rad_wf(n4,l4,xj)
    inte1 = inte1 + xxi**2*xxj**2*osc1*osc2*potential(xxi,xxj,v0r,kr)*coeffi*coeffj
    inte2 = inte2 + xxi**2*xxj**2*osc1*osc2*potential(xxi,xxj,v0s,ks)*coeffi*coeffj
  enddo !j
 enddo !i
resu = -inte1 + inte2
!write(*,*) "Integral result",resu
end
